#!/usr/bin/env python3
"""Plot benchmark: CPU vs GPU per chromosome.

Reads a CSV with columns:
    chromosome, backend, phase1_sec, phase2_sec, time_sec

Produces a multi-panel figure:
  - Top:    stacked bar chart (Phase 1 + Phase 2) per chromosome
  - Bottom: Phase 2 speedup (CPU / GPU) with summary statistics

Infers the method name from the CSV filename (e.g. ``voting_timings.csv``
→ "voting", ``likelihood_timings.csv`` → "likelihood").

Usage:
    python scripts/plot_voting_benchmark.py <timings.csv> <output.png>
"""

import csv
import sys
from pathlib import Path


def _float(row, key, default=-1.0):
    try:
        return float(row.get(key, default))
    except (ValueError, TypeError):
        return default


def _chrom_sort_key(name):
    """Sort chr1..chr22 numerically, then chrX, chrY, others."""
    s = name.replace("chr", "")
    if s.isdigit():
        return (0, int(s))
    return (1, ord(s[0]) if s else 0)


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
        from matplotlib.patches import Patch
        import numpy as np
    except ImportError:
        print("matplotlib is required: pip install 'ancify[evaluate]'")
        sys.exit(1)

    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row.get("chromosome", "").strip()
            backend = row.get("backend", "").strip().lower()
            p1 = _float(row, "phase1_sec")
            p2 = _float(row, "phase2_sec")
            total = _float(row, "time_sec")
            if chrom and backend and total >= 0:
                rows.append({"chrom": chrom, "backend": backend,
                             "p1": max(p1, 0), "p2": max(p2, 0),
                             "total": total})

    if not rows:
        print("No valid rows in CSV.")
        sys.exit(1)

    chroms = sorted({r["chrom"] for r in rows}, key=_chrom_sort_key)
    n = len(chroms)

    lookup = {(r["chrom"], r["backend"]): r for r in rows}

    cpu_p1 = np.array([lookup.get((c, "cpu"), {}).get("p1", 0) for c in chroms])
    cpu_p2 = np.array([lookup.get((c, "cpu"), {}).get("p2", 0) for c in chroms])
    gpu_p1 = np.array([lookup.get((c, "gpu"), {}).get("p1", 0) for c in chroms])
    gpu_p2 = np.array([lookup.get((c, "gpu"), {}).get("p2", 0) for c in chroms])

    # Speedup (CPU / GPU) for Phase 2
    with np.errstate(divide="ignore", invalid="ignore"):
        speedup_p2 = np.where(gpu_p2 > 0, cpu_p2 / gpu_p2, np.nan)

    # ── Palette ──────────────────────────────────────────────
    C_CPU_P1 = "#93c5fd"  # light blue
    C_CPU_P2 = "#2563eb"  # blue
    C_GPU_P1 = "#6ee7b7"  # light green
    C_GPU_P2 = "#059669"  # teal
    C_SPEEDUP = "#7c3aed"  # purple
    C_GRID = "#e5e7eb"

    # ── Infer method name from CSV filename ────────────────
    method = csv_path.stem.replace("_timings", "").replace("_", " ").title()

    # ── Figure layout ────────────────────────────────────────
    fig_w = max(12, n * 0.6)
    fig, (ax_bars, ax_speed) = plt.subplots(
        2, 1, figsize=(fig_w, 8),
        gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.35},
    )

    x = np.arange(n)
    w = 0.35

    # ── Top panel: stacked Phase 1 + Phase 2 bars ───────────
    ax_bars.bar(x - w / 2, cpu_p1, w, label="CPU Phase 1 (project)",
                color=C_CPU_P1, edgecolor="white", linewidth=0.5)
    ax_bars.bar(x - w / 2, cpu_p2, w, bottom=cpu_p1,
                label="CPU Phase 2 (voting)",
                color=C_CPU_P2, edgecolor="white", linewidth=0.5)
    ax_bars.bar(x + w / 2, gpu_p1, w, label="GPU Phase 1 (project)",
                color=C_GPU_P1, edgecolor="white", linewidth=0.5)
    ax_bars.bar(x + w / 2, gpu_p2, w, bottom=gpu_p1,
                label="GPU Phase 2 (voting)",
                color=C_GPU_P2, edgecolor="white", linewidth=0.5)

    cpu_total = cpu_p1 + cpu_p2
    gpu_total = gpu_p1 + gpu_p2

    label_step = 1 if n <= 8 else (2 if n <= 16 else 3)
    for i in range(n):
        if i % label_step == 0:
            if cpu_total[i] > 0:
                ax_bars.text(x[i] - w / 2, cpu_total[i] + 2, f"{cpu_total[i]:.0f}",
                             ha="center", va="bottom", fontsize=6.5,
                             color=C_CPU_P2, fontweight="bold")
            if gpu_total[i] > 0:
                ax_bars.text(x[i] + w / 2, gpu_total[i] + 2, f"{gpu_total[i]:.0f}",
                             ha="center", va="bottom", fontsize=6.5,
                             color=C_GPU_P2, fontweight="bold")

    ax_bars.set_ylabel("Run-time (seconds)", fontsize=11, fontweight="medium")
    ax_bars.set_xticks(x)
    ax_bars.set_xticklabels([c.replace("chr", "") for c in chroms],
                            fontsize=9)
    ax_bars.set_xlabel("")
    ax_bars.set_xlim(-0.6, n - 0.4)
    ax_bars.set_title(f"{method} — CPU vs GPU  (4 outgroups, Phase 1 + Phase 2 stacked)",
                      fontsize=13, fontweight="bold", pad=12)

    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)
    ax_bars.yaxis.grid(True, color=C_GRID, linewidth=0.6, zorder=0)
    ax_bars.set_axisbelow(True)

    legend_handles = [
        Patch(facecolor=C_CPU_P1, edgecolor="white", label="CPU Phase 1"),
        Patch(facecolor=C_CPU_P2, edgecolor="white", label="CPU Phase 2"),
        Patch(facecolor=C_GPU_P1, edgecolor="white", label="GPU Phase 1"),
        Patch(facecolor=C_GPU_P2, edgecolor="white", label="GPU Phase 2"),
    ]
    ax_bars.legend(handles=legend_handles, loc="upper right", fontsize=8,
                   ncol=2, framealpha=0.9, edgecolor=C_GRID)

    # ── Bottom panel: Phase 2 speedup ────────────────────────
    valid = ~np.isnan(speedup_p2)
    colors = np.where(speedup_p2[valid] >= 1, C_SPEEDUP, "#dc2626")

    ax_speed.vlines(x[valid], 1, speedup_p2[valid], color=colors,
                    linewidth=2.5, zorder=3)
    ax_speed.scatter(x[valid], speedup_p2[valid], color=colors,
                     s=40, zorder=4, edgecolors="white", linewidths=0.5)

    ax_speed.axhline(y=1.0, color="#94a3b8", linewidth=1, linestyle="--",
                     zorder=2)

    mean_speedup = np.nanmean(speedup_p2)
    if not np.isnan(mean_speedup):
        ax_speed.axhline(y=mean_speedup, color=C_SPEEDUP, linewidth=1,
                         linestyle=":", alpha=0.6, zorder=2)
        ax_speed.text(n - 0.5, mean_speedup, f" avg {mean_speedup:.1f}×",
                      va="center", ha="left", fontsize=8.5,
                      color=C_SPEEDUP, fontweight="bold")

    if n <= 12:
        for i in range(n):
            if valid[i]:
                ax_speed.text(x[i], speedup_p2[i] + 0.08,
                              f"{speedup_p2[i]:.1f}×",
                              ha="center", va="bottom", fontsize=6.5,
                              color=C_SPEEDUP, fontweight="bold")

    ax_speed.set_ylabel("Phase 2 speedup\n(CPU / GPU)", fontsize=10,
                        fontweight="medium")
    ax_speed.set_xlabel("Chromosome", fontsize=11, fontweight="medium")
    ax_speed.set_xticks(x)
    ax_speed.set_xticklabels([c.replace("chr", "") for c in chroms],
                             fontsize=9)
    ax_speed.set_xlim(-0.6, n - 0.4)
    ymin = min(0.5, np.nanmin(speedup_p2) - 0.3) if valid.any() else 0.5
    ymax = max(2.0, np.nanmax(speedup_p2) + 0.5) if valid.any() else 2.0
    ax_speed.set_ylim(ymin, ymax)
    ax_speed.set_title("Phase 2 GPU speedup (higher = GPU faster)",
                       fontsize=11, fontweight="bold", pad=8)

    ax_speed.spines["top"].set_visible(False)
    ax_speed.spines["right"].set_visible(False)
    ax_speed.yaxis.grid(True, color=C_GRID, linewidth=0.6, zorder=0)
    ax_speed.set_axisbelow(True)

    # ── Summary text ─────────────────────────────────────────
    cpu_genome = cpu_total.sum()
    gpu_genome = gpu_total.sum()
    overall_speedup = cpu_genome / gpu_genome if gpu_genome > 0 else float("nan")
    cpu_p2_total = cpu_p2.sum()
    gpu_p2_total = gpu_p2.sum()
    p2_speedup = cpu_p2_total / gpu_p2_total if gpu_p2_total > 0 else float("nan")

    summary = (
        f"Genome-wide totals:  "
        f"CPU {cpu_genome/60:.0f} min  |  GPU {gpu_genome/60:.0f} min  |  "
        f"overall {overall_speedup:.2f}×    "
        f"Phase 2 only:  "
        f"CPU {cpu_p2_total:.0f}s  |  GPU {gpu_p2_total:.0f}s  |  "
        f"{p2_speedup:.1f}× speedup"
    )
    fig.text(0.5, -0.01, summary, ha="center", va="top", fontsize=9,
             fontstyle="italic", color="#64748b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
