"""Tests for ancify.config."""

import textwrap

import pytest

from ancify.config import (
    EvaluationConfig,
    OutgroupSpec,
    PipelineConfig,
    load_config,
    validate_config,
)


def _write_config(tmp_path, *, extra=""):
    """Helper to write a valid config file with optional extra lines."""
    inner_axt = tmp_path / "inner.axt"
    outer_axt = tmp_path / "outer.axt"
    inner_axt.write_text("")
    outer_axt.write_text("")

    lengths = tmp_path / "chromoLens.txt"
    lengths.write_text("chr1\t100\nchr2\t200\n")

    lines = [
        f"focal_species: test",
        f"chromosome_lengths: {lengths}",
        f"outgroups:",
        f"  inner:",
        f"    - name: sp_inner",
        f"      alignment: {inner_axt}",
        f"  outer:",
        f"    - name: sp_outer",
        f"      alignment: {outer_axt}",
        f"work_dir: {tmp_path}",
        f"output_dir: {tmp_path / 'output'}",
        f"num_cpus: 2",
    ]
    if extra:
        lines.append(extra)

    p = tmp_path / "config.yaml"
    p.write_text("\n".join(lines) + "\n")
    return p


class TestLoadConfig:
    def test_minimal_config(self, tmp_path):
        cfg_path = _write_config(tmp_path)
        cfg = load_config(cfg_path)
        assert cfg.focal_species == "test"
        assert len(cfg.outgroups_inner) == 1
        assert len(cfg.outgroups_outer) == 1
        assert cfg.num_cpus == 2

    def test_defaults(self, tmp_path):
        cfg_path = _write_config(tmp_path)
        cfg = load_config(cfg_path)
        assert cfg.min_inner_freq == 1
        assert cfg.min_outer_freq == 1

    def test_resolve_chromosomes_from_file(self, tmp_path):
        cfg_path = _write_config(tmp_path)
        cfg = load_config(cfg_path)
        chroms = cfg.resolve_chromosomes()
        assert chroms == ["chr1", "chr2"]

    def test_explicit_chromosomes(self, tmp_path):
        extra = "chromosomes:\n  - chr2"
        cfg_path = _write_config(tmp_path, extra=extra)
        cfg = load_config(cfg_path)
        assert cfg.resolve_chromosomes() == ["chr2"]

    def test_evaluation_config(self, tmp_path):
        ref_dir = tmp_path / "ref"
        ref_dir.mkdir()
        extra = (
            f"evaluation:\n"
            f"  reference_dir: {ref_dir}\n"
            f'  reference_pattern: "anc_{{chrom_id}}.fa"'
        )
        cfg_path = _write_config(tmp_path, extra=extra)
        cfg = load_config(cfg_path)
        assert cfg.evaluation is not None
        assert cfg.evaluation.reference_pattern == "anc_{chrom_id}.fa"

    def test_no_evaluation_is_none(self, tmp_path):
        cfg_path = _write_config(tmp_path)
        cfg = load_config(cfg_path)
        assert cfg.evaluation is None


class TestValidateConfig:
    def test_missing_inner_outgroup(self, tmp_path):
        cfg = PipelineConfig(
            focal_species="test",
            chromosome_lengths=str(tmp_path / "chromoLens.txt"),
            outgroups_inner=[],
            outgroups_outer=[OutgroupSpec("sp", str(tmp_path / "x.axt"))],
        )
        (tmp_path / "chromoLens.txt").write_text("chr1\t100\n")
        (tmp_path / "x.axt").write_text("")
        with pytest.raises(ValueError, match="inner"):
            validate_config(cfg)

    def test_missing_outer_outgroup(self, tmp_path):
        cfg = PipelineConfig(
            focal_species="test",
            chromosome_lengths=str(tmp_path / "chromoLens.txt"),
            outgroups_inner=[OutgroupSpec("sp", str(tmp_path / "x.axt"))],
            outgroups_outer=[],
        )
        (tmp_path / "chromoLens.txt").write_text("chr1\t100\n")
        (tmp_path / "x.axt").write_text("")
        with pytest.raises(ValueError, match="outer"):
            validate_config(cfg)

    def test_missing_chrom_lengths_file(self, tmp_path):
        cfg = PipelineConfig(
            focal_species="test",
            chromosome_lengths=str(tmp_path / "nonexistent.txt"),
            outgroups_inner=[OutgroupSpec("sp1", str(tmp_path / "a.axt"))],
            outgroups_outer=[OutgroupSpec("sp2", str(tmp_path / "b.axt"))],
        )
        with pytest.raises(ValueError, match="Chromosome lengths"):
            validate_config(cfg)

    def test_missing_alignment_file(self, tmp_path):
        lengths = tmp_path / "chromoLens.txt"
        lengths.write_text("chr1\t100\n")
        cfg = PipelineConfig(
            focal_species="test",
            chromosome_lengths=str(lengths),
            outgroups_inner=[OutgroupSpec("sp1", str(tmp_path / "missing.axt"))],
            outgroups_outer=[OutgroupSpec("sp2", str(tmp_path / "b.axt"))],
        )
        with pytest.raises(ValueError, match="Alignment file not found"):
            validate_config(cfg)
