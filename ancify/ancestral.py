"""Phase 2: Infer ancestral alleles from projected outgroup sequences.

Uses a two-tier outgroup voting scheme:

* **Inner outgroup** -- closely related species; the most frequent
  nucleotide among them forms the inner consensus.
* **Outer outgroup** -- more distantly related species; serves as an
  independent confirmation.

Confidence is encoded via letter case in the output FASTA:

========  ===========  ==========================================
Char      Confidence   Condition
========  ===========  ==========================================
``ACGT``  High         Inner and outer outgroups agree
``acgt``  Low          Only one tier has data
``n``     Unresolved   Inner and outer disagree
``N``     Missing      Both tiers lack data
========  ===========  ==========================================
"""

import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from .utils import read_fasta, write_fasta, majority_vote, VALID_ALLELES
from .backend import vectorized_ancestral_call, get_available_gpus
from .parsimony import (
    VALID_BASES,
    fitch_ancestral,
    get_leaf_names,
    parse_newick,
)

logger = logging.getLogger(__name__)


def call_ancestral_base(inner_bases, outer_bases,
                        min_inner_freq=1, min_outer_freq=1):
    """Infer the ancestral allele at a single position.

    Parameters
    ----------
    inner_bases : list of str
        Nucleotides from the inner (closely related) outgroup species.
    outer_bases : list of str
        Nucleotides from the outer (distantly related) outgroup species.
    min_inner_freq, min_outer_freq : int
        Minimum allele count to accept a majority-vote consensus.

    Returns
    -------
    str
        Single character with case-encoded confidence (see module docstring).
    """
    inner = majority_vote(inner_bases, min_inner_freq)
    outer = majority_vote(outer_bases, min_outer_freq)

    if inner != "N" and inner == outer:
        return inner
    if inner == "N" and outer != "N":
        return outer.lower()
    if inner != "N" and outer == "N":
        return inner.lower()
    if inner == "N" and outer == "N":
        return "N"
    return "n"


def _call_chromosome(args):
    """Worker: call ancestral states for one chromosome."""
    chrom, inner_paths, outer_paths, out_path, min_inner, min_outer = args

    inner_seqs = [read_fasta(p)[1] for p in inner_paths]
    outer_seqs = [read_fasta(p)[1] for p in outer_paths]

    length = len(inner_seqs[0])
    for s in inner_seqs + outer_seqs:
        if len(s) != length:
            raise ValueError(
                f"Length mismatch on {chrom}: expected {length}, got {len(s)}"
            )

    anc = []
    for i in range(length):
        ib = [s[i].upper() for s in inner_seqs]
        ob = [s[i].upper() for s in outer_seqs]
        anc.append(call_ancestral_base(ib, ob, min_inner, min_outer))

    write_fasta(out_path, f">{chrom}", "".join(anc))
    return chrom


def _call_chromosome_vectorized(args):
    """Vectorized worker: call ancestral states for one chromosome.

    Uses bulk NumPy (CPU) or PyTorch (GPU) operations instead of a
    per-position Python loop.  *device_str* is ``None`` for the CPU path
    or a CUDA device string like ``"cuda:0"`` for the GPU path.
    """
    chrom, inner_paths, outer_paths, out_path, min_inner, min_outer, device_str = args

    device = None
    if device_str is not None:
        import torch
        device = torch.device(device_str)

    inner_seqs = [read_fasta(p)[1] for p in inner_paths]
    outer_seqs = [read_fasta(p)[1] for p in outer_paths]

    length = len(inner_seqs[0])
    for s in inner_seqs + outer_seqs:
        if len(s) != length:
            raise ValueError(
                f"Length mismatch on {chrom}: expected {length}, got {len(s)}"
            )

    anc = vectorized_ancestral_call(
        inner_seqs, outer_seqs, min_inner, min_outer, device,
    )
    write_fasta(out_path, f">{chrom}", anc)
    return chrom


def call_ancestral_base_parsimony(tree, species_bases):
    """Infer the ancestral allele at a single position using Fitch parsimony.

    Parameters
    ----------
    tree : TreeNode
        Phylogenetic tree of outgroup species.
    species_bases : dict
        Mapping of species name → observed nucleotide (single character).

    Returns
    -------
    str
        Single character with case-encoded confidence:
        uppercase = unambiguous, lowercase = ambiguous, ``N`` = all missing.
    """
    leaf_names = get_leaf_names(tree)
    has_data = any(
        species_bases.get(name, "N").upper() in VALID_BASES
        for name in leaf_names
    )
    if not has_data:
        return "N"

    allele, is_ambiguous = fitch_ancestral(tree, species_bases)
    if is_ambiguous:
        return allele.lower()
    return allele.upper()


def _call_chromosome_parsimony(args):
    """Worker: call ancestral states for one chromosome using Fitch parsimony."""
    chrom, species_paths, tree_text, out_path = args

    tree = parse_newick(tree_text)
    species_seqs = {}
    length = None
    for name, path in species_paths.items():
        _, seq = read_fasta(path)
        species_seqs[name] = seq
        if length is None:
            length = len(seq)
        elif len(seq) != length:
            raise ValueError(
                f"Length mismatch on {chrom}: expected {length}, got {len(seq)}"
            )

    anc = []
    for i in range(length):
        bases = {name: seq[i].upper() for name, seq in species_seqs.items()}
        anc.append(call_ancestral_base_parsimony(tree, bases))

    write_fasta(out_path, f">{chrom}", "".join(anc))
    return chrom


def run_ancestral_calling(config):
    """Execute Phase 2: call ancestral alleles for every chromosome.

    Reads projected FASTA files from ``<work_dir>/projected/<species>/``
    and writes ancestral FASTA files to ``<output_dir>/<chrom>.fa``.

    Supported methods: ``"voting"`` (default), ``"parsimony"``, ``"ml"``.
    """
    method = getattr(config, "method", "voting")

    if method == "parsimony":
        return _run_parsimony(config)
    if method == "ml":
        return _run_ml(config)

    _run_voting(config)


def _run_parsimony(config):
    """Parsimony-based ancestral calling for all chromosomes."""
    chromosomes = config.resolve_chromosomes()
    work = Path(config.work_dir)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_outgroups = config.outgroups_inner + config.outgroups_outer

    tasks = []
    for chrom in chromosomes:
        species_paths = {
            og.name: str(work / "projected" / og.name / f"{chrom}.fa")
            for og in all_outgroups
        }
        out_path = str(out_dir / f"{chrom}.fa")
        tasks.append((chrom, species_paths, config.tree, out_path))

    logger.info(
        "Phase 2: calling ancestral states for %d chromosomes [parsimony]",
        len(tasks),
    )

    with ProcessPoolExecutor(max_workers=config.num_cpus) as pool:
        futures = {
            pool.submit(_call_chromosome_parsimony, t): t for t in tasks
        }
        for future in as_completed(futures):
            chrom = future.result()
            logger.info("  Completed %s", chrom)

    logger.info("Phase 2 complete.")


def _run_ml(config):
    """ML-based ancestral calling for all chromosomes."""
    from .ml import _call_chromosome_ml

    chromosomes = config.resolve_chromosomes()
    work = Path(config.work_dir)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    model_path = config.ml_model_path
    if not model_path:
        raise ValueError(
            "method 'ml' requires 'ml_model_path' pointing to a trained model. "
            "Run 'ancify train' first."
        )

    tasks = []
    for chrom in chromosomes:
        inner_paths = [
            str(work / "projected" / og.name / f"{chrom}.fa")
            for og in config.outgroups_inner
        ]
        outer_paths = [
            str(work / "projected" / og.name / f"{chrom}.fa")
            for og in config.outgroups_outer
        ]
        out_path = str(out_dir / f"{chrom}.fa")
        tasks.append((
            chrom, inner_paths, outer_paths, out_path,
            model_path, config.ml_high_threshold, config.ml_low_threshold,
        ))

    logger.info(
        "Phase 2: calling ancestral states for %d chromosomes [ml]",
        len(tasks),
    )

    with ProcessPoolExecutor(max_workers=config.num_cpus) as pool:
        futures = {
            pool.submit(_call_chromosome_ml, t): t for t in tasks
        }
        for future in as_completed(futures):
            chrom = future.result()
            logger.info("  Completed %s", chrom)

    logger.info("Phase 2 complete.")


def _run_voting(config):
    """Two-tier voting ancestral calling (original algorithm)."""
    chromosomes = config.resolve_chromosomes()
    work = Path(config.work_dir)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── resolve backend and GPU devices ──────────────────────────────
    backend_mode = getattr(config, "backend", "auto")
    gpu_device_ids = getattr(config, "gpu_devices", None)

    use_gpu = False
    devices = []

    if backend_mode in ("auto", "gpu"):
        gpus = get_available_gpus()
        if gpu_device_ids is not None:
            gpus = [g for g in gpus if g in gpu_device_ids]
        if gpus:
            import torch  # noqa: F811 – guarded import
            devices = [torch.device(f"cuda:{i}") for i in gpus]
            use_gpu = True

    if backend_mode == "gpu" and not use_gpu:
        logger.warning(
            "GPU backend requested but no CUDA devices available; "
            "falling back to CPU vectorised path"
        )

    # ── build task list ──────────────────────────────────────────────
    tasks = []
    for i, chrom in enumerate(chromosomes):
        inner_paths = [
            str(work / "projected" / og.name / f"{chrom}.fa")
            for og in config.outgroups_inner
        ]
        outer_paths = [
            str(work / "projected" / og.name / f"{chrom}.fa")
            for og in config.outgroups_outer
        ]
        out_path = str(out_dir / f"{chrom}.fa")
        device_str = str(devices[i % len(devices)]) if use_gpu else None
        tasks.append((
            chrom, inner_paths, outer_paths, out_path,
            config.min_inner_freq, config.min_outer_freq, device_str,
        ))

    label = f"gpu ({len(devices)} device(s))" if use_gpu else "cpu"
    logger.info(
        "Phase 2: calling ancestral states for %d chromosomes [%s]",
        len(tasks), label,
    )

    # ── dispatch ─────────────────────────────────────────────────────
    if use_gpu:
        # GPU ops release the GIL; threads overlap I/O and compute.
        max_workers = min(len(tasks), len(devices) * 2) or 1
        pool_cls = ThreadPoolExecutor
    else:
        max_workers = config.num_cpus
        pool_cls = ProcessPoolExecutor

    with pool_cls(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_call_chromosome_vectorized, t): t for t in tasks
        }
        for future in as_completed(futures):
            chrom = future.result()
            logger.info("  Completed %s", chrom)

    logger.info("Phase 2 complete.")
