"""Configuration loading and validation for ancify."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import yaml


@dataclass
class OutgroupSpec:
    """Specification for a single outgroup species."""
    name: str
    alignment: str


@dataclass
class EvaluationConfig:
    """Optional evaluation settings."""
    reference_dir: Optional[str] = None
    reference_pattern: str = "{chrom}.fa"
    vcf_dir: Optional[str] = None
    vcf_pattern: str = "{chrom}.vcf.gz"


@dataclass
class PipelineConfig:
    """Complete pipeline configuration."""
    focal_species: str
    chromosome_lengths: str
    outgroups_inner: List[OutgroupSpec]
    outgroups_outer: List[OutgroupSpec]
    chromosomes: Optional[List[str]] = None
    work_dir: str = "."
    output_dir: str = "./ancestral_calls"
    min_inner_freq: int = 1
    min_outer_freq: int = 1
    num_cpus: int = 4
    backend: str = "auto"
    gpu_devices: Optional[List[int]] = None
    method: str = "voting"
    tree: Optional[str] = None
    ml_model_path: Optional[str] = None
    ml_training_reference: Optional[str] = None
    ml_high_threshold: float = 0.8
    ml_low_threshold: float = 0.5
    substitution_model: str = "JC69"
    model_kappa: float = 2.0
    model_base_freqs: Optional[List[float]] = None
    model_rates: Optional[List[float]] = None
    likelihood_high_threshold: float = 0.8
    likelihood_low_threshold: float = 0.5
    evaluation: Optional[EvaluationConfig] = None

    def resolve_chromosomes(self):
        """Return the list of chromosomes to process.

        If *chromosomes* was not set explicitly, reads all chromosome
        names from the chromosome-lengths file.
        """
        if self.chromosomes is None:
            from .utils import read_chromosome_lengths
            lengths = read_chromosome_lengths(self.chromosome_lengths)
            self.chromosomes = list(lengths.keys())
        return self.chromosomes

    @property
    def all_outgroups(self) -> "List[OutgroupSpec]":
        """All outgroup species (inner + outer), deduplicated by name."""
        seen: set = set()
        result: list = []
        for og in self.outgroups_inner + self.outgroups_outer:
            if og.name not in seen:
                seen.add(og.name)
                result.append(og)
        return result


def load_config(path):
    """Load a YAML configuration file and return a PipelineConfig."""
    with open(path) as f:
        raw = yaml.safe_load(f)

    inner = [OutgroupSpec(**s) for s in raw["outgroups"]["inner"]]
    outer = [OutgroupSpec(**s) for s in raw["outgroups"]["outer"]]

    eval_cfg = None
    if raw.get("evaluation"):
        eval_cfg = EvaluationConfig(**raw["evaluation"])

    tree_raw = raw.get("tree")
    if tree_raw is not None:
        tree_path = Path(tree_raw)
        if tree_path.suffix in (".nwk", ".tree", ".newick") or (
            "(" not in tree_raw and tree_path.is_file()
        ):
            tree_raw = tree_path.read_text().strip()

    cfg = PipelineConfig(
        focal_species=raw["focal_species"],
        chromosome_lengths=raw["chromosome_lengths"],
        outgroups_inner=inner,
        outgroups_outer=outer,
        chromosomes=raw.get("chromosomes"),
        work_dir=raw.get("work_dir", "."),
        output_dir=raw.get("output_dir", "./ancestral_calls"),
        min_inner_freq=raw.get("min_inner_freq", 1),
        min_outer_freq=raw.get("min_outer_freq", 1),
        num_cpus=raw.get("num_cpus", 4),
        backend=raw.get("backend", "auto"),
        gpu_devices=raw.get("gpu_devices"),
        method=raw.get("method", "voting"),
        tree=tree_raw,
        ml_model_path=raw.get("ml_model_path"),
        ml_training_reference=raw.get("ml_training_reference"),
        ml_high_threshold=float(raw.get("ml_high_threshold", 0.8)),
        ml_low_threshold=float(raw.get("ml_low_threshold", 0.5)),
        substitution_model=raw.get("substitution_model", "JC69"),
        model_kappa=float(raw.get("model_kappa", 2.0)),
        model_base_freqs=raw.get("model_base_freqs"),
        model_rates=raw.get("model_rates"),
        likelihood_high_threshold=float(raw.get("likelihood_high_threshold", 0.8)),
        likelihood_low_threshold=float(raw.get("likelihood_low_threshold", 0.5)),
        evaluation=eval_cfg,
    )

    validate_config(cfg)
    return cfg


def validate_config(cfg):
    """Check that a PipelineConfig is self-consistent.

    Raises ValueError with a descriptive message on any problem.
    """
    if cfg.method not in ("voting", "parsimony", "ml", "likelihood"):
        raise ValueError(
            f"Unknown method: {cfg.method!r} "
            "(must be 'voting', 'parsimony', 'ml', or 'likelihood')"
        )

    if cfg.method == "voting":
        if not cfg.outgroups_inner:
            raise ValueError("At least one inner outgroup species is required.")
        if not cfg.outgroups_outer:
            raise ValueError("At least one outer outgroup species is required.")
    elif cfg.method == "ml":
        if not cfg.outgroups_inner:
            raise ValueError("At least one inner outgroup species is required.")
        if not cfg.outgroups_outer:
            raise ValueError("At least one outer outgroup species is required.")
        if cfg.ml_model_path and not Path(cfg.ml_model_path).exists():
            raise ValueError(
                f"ML model file not found: {cfg.ml_model_path}"
            )
        if not (0.0 <= cfg.ml_low_threshold <= cfg.ml_high_threshold <= 1.0):
            raise ValueError(
                "ML thresholds must satisfy 0 <= low <= high <= 1; "
                f"got low={cfg.ml_low_threshold}, high={cfg.ml_high_threshold}"
            )
    elif cfg.method == "likelihood":
        all_outgroups = cfg.outgroups_inner + cfg.outgroups_outer
        if not all_outgroups:
            raise ValueError(
                "At least one outgroup species is required for likelihood."
            )
        valid_models = ("JC69", "K80", "HKY85", "GTR")
        if cfg.substitution_model.upper() not in valid_models:
            raise ValueError(
                f"Unknown substitution_model: {cfg.substitution_model!r} "
                f"(must be one of {valid_models})"
            )
        if cfg.substitution_model.upper() == "GTR":
            if cfg.model_rates is None or len(cfg.model_rates) != 6:
                raise ValueError(
                    "substitution_model 'GTR' requires 'model_rates' "
                    "with exactly 6 exchangeability rates [AC, AG, AT, CG, CT, GT]"
                )
        if cfg.model_base_freqs is not None:
            if len(cfg.model_base_freqs) != 4:
                raise ValueError(
                    "'model_base_freqs' must have exactly 4 values [A, C, G, T]"
                )
            total = sum(cfg.model_base_freqs)
            if abs(total - 1.0) > 0.01:
                raise ValueError(
                    f"'model_base_freqs' must sum to ~1.0 (got {total:.4f})"
                )
        if not (0.0 <= cfg.likelihood_low_threshold <= cfg.likelihood_high_threshold <= 1.0):
            raise ValueError(
                "Likelihood thresholds must satisfy 0 <= low <= high <= 1; "
                f"got low={cfg.likelihood_low_threshold}, high={cfg.likelihood_high_threshold}"
            )
    else:
        all_outgroups = cfg.outgroups_inner + cfg.outgroups_outer
        if not all_outgroups:
            raise ValueError(
                "At least one outgroup species is required for parsimony."
            )

    if not Path(cfg.chromosome_lengths).exists():
        raise ValueError(
            f"Chromosome lengths file not found: {cfg.chromosome_lengths}"
        )
    for og in cfg.outgroups_inner + cfg.outgroups_outer:
        if not Path(og.alignment).exists():
            raise ValueError(
                f"Alignment file not found for {og.name}: {og.alignment}"
            )

    if cfg.method in ("parsimony", "likelihood"):
        if not cfg.tree:
            raise ValueError(
                f"method {cfg.method!r} requires a 'tree' field "
                "(Newick string or path to .nwk file)"
            )
        from .parsimony import get_leaf_names, parse_newick

        tree = parse_newick(cfg.tree)
        leaf_names = set(get_leaf_names(tree))
        outgroup_names = {
            og.name for og in cfg.outgroups_inner + cfg.outgroups_outer
        }
        missing = leaf_names - outgroup_names
        if missing:
            raise ValueError(
                f"Tree leaf names not found in outgroups: {sorted(missing)}"
            )
        extra = outgroup_names - leaf_names
        if extra:
            raise ValueError(
                f"Outgroup names not found as tree leaves: {sorted(extra)}"
            )
