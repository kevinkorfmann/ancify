#!/usr/bin/env python3
"""Plot voting run-time: CPU vs GPU per chromosome.

Reads a CSV with columns: chromosome, backend, time_sec
or (stratified): chromosome, backend, phase1_sec, phase2_sec, time_sec.
Writes grouped bar chart(s) to the given output path.
If phase1_sec/phase2_sec are present, draws two panels: Phase 1 (project) and Phase 2 (voting).

Usage:
    python scripts/plot_voting_benchmark.py <timings.csv> <output.png>
"""

import csv
import sys
from pathlib import Path


def _read_float(row, key, default=-1):
    try:
        return float(row.get(key, default))
    except (ValueError, TypeError):
        return default


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    if not csv_path.exists():
        print(f"File not found: {csv_path}")
        sys.exit(1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib is required: pip install 'ancify[evaluate]'")
        sys.exit(1)

    # Load data
    rows = []
    has_phase = False
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        has_phase = "phase1_sec" in fieldnames and "phase2_sec" in fieldnames
        for row in reader:
            chrom = row.get("chromosome", "").strip()
            backend = row.get("backend", "").strip().lower()
            t = _read_float(row, "time_sec")
            if chrom and backend and t >= 0:
                p1 = _read_float(row, "phase1_sec") if has_phase else -1
                p2 = _read_float(row, "phase2_sec") if has_phase else -1
                rows.append((chrom, backend, p1, p2, t))

    if not rows:
        print("No valid rows in CSV (need chromosome, backend, time_sec with time_sec >= 0).")
        sys.exit(1)

    chroms = []
    seen = set()
    for r in rows:
        c = r[0]
        if c not in seen:
            seen.add(c)
            chroms.append(c)

    n = len(chroms)
    x = np.arange(n)
    width = 0.35

    def series(key_idx):
        # key_idx: 2=phase1, 3=phase2, 4=total
        data = {(r[0], r[1]): r[key_idx] for r in rows}
        cpu = [data.get((c, "cpu")) or 0 for c in chroms]
        gpu = [data.get((c, "gpu")) or 0 for c in chroms]
        return cpu, gpu

    def draw_bars(ax, cpu_vals, gpu_vals, title, ylabel="Run-time (seconds)"):
        bars_cpu = ax.bar(x - width / 2, cpu_vals, width, label="CPU", color="#2563eb", edgecolor="white", linewidth=0.8)
        bars_gpu = ax.bar(x + width / 2, gpu_vals, width, label="GPU", color="#059669", edgecolor="white", linewidth=0.8)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlabel("Chromosome", fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(chroms, fontsize=10)
        ax.legend(loc="upper right", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for b in bars_cpu:
            h = b.get_height()
            if h > 0:
                ax.annotate(f"{h:.0f}s", xy=(b.get_x() + b.get_width() / 2, h), xytext=(0, 4),
                            textcoords="offset points", ha="center", va="bottom", fontsize=8, fontweight="bold")
        for b in bars_gpu:
            h = b.get_height()
            if h > 0:
                ax.annotate(f"{h:.0f}s", xy=(b.get_x() + b.get_width() / 2, h), xytext=(0, 4),
                            textcoords="offset points", ha="center", va="bottom", fontsize=8, fontweight="bold")

    if has_phase and any(r[2] >= 0 or r[3] >= 0 for r in rows):
        # Stratified: two panels (Phase 1 and Phase 2)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(10, n * 2.2), 5))
        cpu_p1, gpu_p1 = series(2)
        cpu_p2, gpu_p2 = series(3)
        draw_bars(ax1, cpu_p1, gpu_p1, "Phase 1: Project (AXT → FASTA)")
        draw_bars(ax2, cpu_p2, gpu_p2, "Phase 2: Voting (ancestral calling)")
        fig.suptitle("Voting run-time: CPU vs GPU (human hg38) — by phase", fontsize=13, fontweight="bold", y=1.02)
    else:
        # Total only
        fig, ax = plt.subplots(figsize=(max(6, n * 1.2), 5))
        cpu_times, gpu_times = series(4)
        draw_bars(ax, cpu_times, gpu_times, "Voting run-time: CPU vs GPU (human hg38)")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
