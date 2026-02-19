"""Tests for ancify.cli."""

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest


def run_ancify(*args):
    """Run ancify as a subprocess and return the CompletedProcess."""
    return subprocess.run(
        [sys.executable, "-m", "ancify", *args],
        capture_output=True,
        text=True,
    )


class TestCLIHelp:
    def test_help(self):
        r = run_ancify("--help")
        assert r.returncode == 0
        assert "Ancestral allele polarization" in r.stdout

    def test_no_args_shows_help(self):
        r = run_ancify()
        assert r.returncode == 1
        assert "usage:" in r.stdout.lower() or "usage:" in r.stderr.lower()

    def test_project_help(self):
        r = run_ancify("project", "--help")
        assert r.returncode == 0
        assert "--config" in r.stdout

    def test_call_help(self):
        r = run_ancify("call", "--help")
        assert r.returncode == 0

    def test_evaluate_help(self):
        r = run_ancify("evaluate", "--help")
        assert r.returncode == 0

    def test_run_help(self):
        r = run_ancify("run", "--help")
        assert r.returncode == 0

    def test_init_help(self):
        r = run_ancify("init", "--help")
        assert r.returncode == 0


class TestCLIInit:
    def test_init_creates_file(self, tmp_path):
        out = tmp_path / "my_config.yaml"
        r = run_ancify("init", "-o", str(out))
        assert r.returncode == 0
        assert out.exists()
        content = out.read_text()
        assert "focal_species" in content
        assert "outgroups" in content
        assert "inner" in content
        assert "outer" in content

    def test_init_default_name(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        r = run_ancify("init")
        assert r.returncode == 0
        assert (tmp_path / "config.yaml").exists()


class TestCLIIntegration:
    """End-to-end integration tests using tiny synthetic data."""

    def test_full_run(self, full_pipeline_dir):
        """Run the full pipeline on synthetic data."""
        tmp_path, config_path = full_pipeline_dir
        r = run_ancify("run", "-c", str(config_path))
        assert r.returncode == 0, f"STDOUT: {r.stdout}\nSTDERR: {r.stderr}"

        out_dir = tmp_path / "ancestral_calls"
        assert out_dir.exists()
        assert (out_dir / "chr1.fa").exists()

        content = (out_dir / "chr1.fa").read_text()
        assert content.startswith(">chr1")

    def test_project_only(self, full_pipeline_dir):
        """Run only Phase 1."""
        tmp_path, config_path = full_pipeline_dir
        r = run_ancify("project", "-c", str(config_path))
        assert r.returncode == 0, f"STDERR: {r.stderr}"

        for name in ["inner1", "inner2", "outer1"]:
            assert (tmp_path / "projected" / name / "chr1.fa").exists()

    def test_project_then_call(self, full_pipeline_dir):
        """Run Phase 1 then Phase 2 separately."""
        tmp_path, config_path = full_pipeline_dir

        r1 = run_ancify("project", "-c", str(config_path))
        assert r1.returncode == 0

        r2 = run_ancify("call", "-c", str(config_path))
        assert r2.returncode == 0

        assert (tmp_path / "ancestral_calls" / "chr1.fa").exists()

    def test_missing_config_fails(self, tmp_path):
        r = run_ancify("run", "-c", str(tmp_path / "nonexistent.yaml"))
        assert r.returncode != 0

    def test_num_cpus_override(self, full_pipeline_dir):
        """The -n flag overrides num_cpus."""
        _, config_path = full_pipeline_dir
        r = run_ancify("run", "-c", str(config_path), "-n", "1")
        assert r.returncode == 0

    def test_verbose_flag(self, full_pipeline_dir):
        _, config_path = full_pipeline_dir
        r = run_ancify("-v", "run", "-c", str(config_path))
        assert r.returncode == 0

    def test_ancestral_output_has_correct_length(self, full_pipeline_dir):
        """The ancestral FASTA should match the chromosome length."""
        tmp_path, config_path = full_pipeline_dir
        run_ancify("run", "-c", str(config_path))

        from ancify.utils import read_fasta
        _, seq = read_fasta(tmp_path / "ancestral_calls" / "chr1.fa")
        assert len(seq) == 30  # matches chromoLens.txt

    def test_ancestral_output_correct_bases(self, full_pipeline_dir):
        """Where all outgroups agree, the call should be high-confidence."""
        tmp_path, config_path = full_pipeline_dir
        run_ancify("run", "-c", str(config_path))

        from ancify.utils import read_fasta
        _, seq = read_fasta(tmp_path / "ancestral_calls" / "chr1.fa")
        # Positions 0-9 had all species agreeing on ACGTACGTAC
        assert seq[0] == "A"
        assert seq[1] == "C"
        assert seq[2] == "G"
        assert seq[3] == "T"
        for i in range(10):
            assert seq[i].isupper()
