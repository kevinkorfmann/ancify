"""Phase 3: Evaluate ancestral calls against a reference and/or VCF variants.

All evaluation steps are optional and driven by the ``evaluation`` block
in the pipeline configuration file:

* **Coverage statistics** are always produced (no extra data needed).
* **Reference comparison** requires a directory of reference ancestral
  FASTA files (e.g. Ensembl EPO) and a filename pattern.
* **VCF comparison** requires a directory of VCF files and a filename
  pattern, plus the ``scikit-allel`` package.
"""

import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from .utils import read_fasta, chrom_id, VALID_ALLELES

logger = logging.getLogger(__name__)


def compute_coverage_stats(sequence):
    """Compute coverage and confidence statistics for an ancestral sequence."""
    arr = np.array(list(sequence), dtype="U1")
    upper = np.char.upper(arr)

    valid = np.isin(upper, list(VALID_ALLELES))
    high_conf = np.isin(arr, list(VALID_ALLELES))
    low_conf = valid & ~high_conf
    total = len(arr)

    return {
        "total_positions": total,
        "high_confidence": int(np.sum(high_conf)),
        "low_confidence": int(np.sum(low_conf)),
        "missing": total - int(np.sum(valid)),
        "prop_nonmissing": float(np.sum(valid) / total),
        "prop_high_confidence": float(np.sum(high_conf) / total),
    }


def compare_to_reference(test_seq, ref_seq, positions=None):
    """Compare two ancestral sequences, optionally at specific positions."""
    test = np.array(list(test_seq), dtype="U1")
    ref = np.array(list(ref_seq), dtype="U1")

    if positions is not None:
        test = test[positions]
        ref = ref[positions]

    alleles = list(VALID_ALLELES)
    tu = np.char.upper(test)
    ru = np.char.upper(ref)

    tv = np.isin(tu, alleles)
    rv = np.isin(ru, alleles)
    both = tv & rv

    n_both = int(np.sum(both))
    agree = tu[both] == ru[both]

    return {
        "test_nonmissing": float(np.mean(tv)),
        "ref_nonmissing": float(np.mean(rv)),
        "both_nonmissing": n_both,
        "agreement_rate": float(np.mean(agree)) if n_both > 0 else float("nan"),
        "disagreement_rate": float(1 - np.mean(agree)) if n_both > 0 else float("nan"),
    }


def compare_to_vcf(ancestral_seq, vcf_path):
    """Compare ancestral calls against VCF REF/ALT alleles.

    Requires ``scikit-allel`` (install with ``pip install scikit-allel``).
    """
    try:
        import allel
    except ImportError:
        raise ImportError(
            "scikit-allel is required for VCF evaluation.  "
            "Install with:  pip install 'ancify[evaluate]'"
        )

    vcf = allel.read_vcf(
        str(vcf_path),
        fields=["variants/POS", "variants/REF", "variants/ALT"],
    )
    if vcf is None:
        return None

    pos = vcf["variants/POS"]
    ref = vcf["variants/REF"]
    alt = vcf["variants/ALT"]
    if alt.ndim == 2:
        alt = alt[:, 0]

    anc = np.array(list(ancestral_seq), dtype="U1")[pos - 1]
    anc_upper = np.char.upper(anc)

    alleles = list(VALID_ALLELES)
    valid = np.isin(anc_upper, alleles)

    matches_ref = anc_upper[valid] == ref[valid]
    matches_alt = anc_upper[valid] == alt[valid]
    matches_either = matches_ref | matches_alt

    return {
        "num_variants": len(pos),
        "prop_nonmissing": float(np.mean(valid)),
        "matches_ref": float(np.mean(matches_ref)),
        "matches_alt": float(np.mean(matches_alt)),
        "matches_either": float(np.mean(matches_either)),
    }


def _format_pattern(pattern, chrom):
    """Substitute ``{chrom}`` and ``{chrom_id}`` placeholders."""
    return pattern.format(chrom=chrom, chrom_id=chrom_id(chrom))


def _evaluate_chromosome(args):
    """Worker: evaluate one chromosome."""
    chrom, anc_path, ref_path, vcf_path, summary_path = args

    _, anc_seq = read_fasta(anc_path)
    results = {"chromosome": chrom}

    results["coverage"] = compute_coverage_stats(anc_seq)

    if ref_path and Path(ref_path).exists():
        _, ref_seq = read_fasta(ref_path)
        results["reference_comparison"] = compare_to_reference(anc_seq, ref_seq)

    if vcf_path and Path(vcf_path).exists():
        vcf_stats = compare_to_vcf(anc_seq, vcf_path)
        if vcf_stats:
            results["vcf_comparison"] = vcf_stats

    if summary_path:
        Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            for section, stats in results.items():
                if isinstance(stats, dict):
                    f.write(f"[{section}]\n")
                    for k, v in stats.items():
                        f.write(f"  {k}: {v}\n")
                    f.write("\n")

    return results


def run_evaluation(config):
    """Execute Phase 3: evaluate ancestral calls for every chromosome.

    Writes per-chromosome summary files to ``<output_dir>/evaluation/``.
    """
    chromosomes = config.resolve_chromosomes()
    out_dir = Path(config.output_dir) / "evaluation"
    out_dir.mkdir(parents=True, exist_ok=True)

    eval_cfg = config.evaluation
    tasks = []

    for chrom in chromosomes:
        anc_path = str(Path(config.output_dir) / f"{chrom}.fa")

        ref_path = None
        if eval_cfg and eval_cfg.reference_dir:
            ref_path = str(
                Path(eval_cfg.reference_dir)
                / _format_pattern(eval_cfg.reference_pattern, chrom)
            )

        vcf_path = None
        if eval_cfg and eval_cfg.vcf_dir:
            vcf_path = str(
                Path(eval_cfg.vcf_dir)
                / _format_pattern(eval_cfg.vcf_pattern, chrom)
            )

        summary_path = str(out_dir / f"{chrom}.evaluation.txt")
        tasks.append((chrom, anc_path, ref_path, vcf_path, summary_path))

    logger.info("Phase 3: evaluating %d chromosomes", len(tasks))

    all_results = []
    with ProcessPoolExecutor(max_workers=config.num_cpus) as pool:
        futures = {pool.submit(_evaluate_chromosome, t): t for t in tasks}
        for future in as_completed(futures):
            result = future.result()
            all_results.append(result)
            logger.info("  Completed %s", result["chromosome"])

    logger.info("Phase 3 complete.")
    return all_results
