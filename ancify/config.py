"""Configuration loading and validation for the polarization pipeline."""

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


def load_config(path):
    """Load a YAML configuration file and return a PipelineConfig."""
    with open(path) as f:
        raw = yaml.safe_load(f)

    inner = [OutgroupSpec(**s) for s in raw["outgroups"]["inner"]]
    outer = [OutgroupSpec(**s) for s in raw["outgroups"]["outer"]]

    eval_cfg = None
    if raw.get("evaluation"):
        eval_cfg = EvaluationConfig(**raw["evaluation"])

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
        evaluation=eval_cfg,
    )

    validate_config(cfg)
    return cfg


def validate_config(cfg):
    """Check that a PipelineConfig is self-consistent.

    Raises ValueError with a descriptive message on any problem.
    """
    if not cfg.outgroups_inner:
        raise ValueError("At least one inner outgroup species is required.")
    if not cfg.outgroups_outer:
        raise ValueError("At least one outer outgroup species is required.")
    if not Path(cfg.chromosome_lengths).exists():
        raise ValueError(
            f"Chromosome lengths file not found: {cfg.chromosome_lengths}"
        )
    for og in cfg.outgroups_inner + cfg.outgroups_outer:
        if not Path(og.alignment).exists():
            raise ValueError(
                f"Alignment file not found for {og.name}: {og.alignment}"
            )
