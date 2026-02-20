#!/usr/bin/env python3
"""Visualize and compare ancestral calling results across methods.

Reads per-method FASTA output and generates:
  1. Stacked bar chart — confidence breakdown per method
  2. Pairwise agreement heatmap
  3. Sliding-window disagreement plot along the chromosome

Usage:
    python scripts/visualize_comparison.py <output_dir> <chrom> [methods...]

Example:
    python scripts/visualize_comparison.py ancify_compare/output_chr22 chr22 voting parsimony likelihood
"""

import sys
from pathlib import Path

import numpy as np


def load_fasta_seq(path):
    lines = Path(path).read_text().splitlines()
    return "".join(l for l in lines if not l.startswith(">"))


def classify(seq):
    arr = np.array(list(seq), dtype="U1")
    upper = np.char.upper(arr)
    high = np.isin(arr, list("ACGT"))
    low = np.isin(arr, list("acgt"))
    unresolved = arr == "n"
    missing = arr == "N"
    return arr, upper, high, low, unresolved, missing


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    out_base = Path(sys.argv[1])
    chrom = sys.argv[2]
    methods = sys.argv[3:] if len(sys.argv) > 3 else ["voting", "parsimony", "likelihood"]

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError:
        print("matplotlib is required: pip install 'ancify[evaluate]'")
        sys.exit(1)

    seqs = {}
    uppers = {}
    stats = {}

    for m in methods:
        fa = out_base / m / f"{chrom}.fa"
        if not fa.exists():
            print(f"  skipping {m}: {fa} not found")
            continue
        seq = load_fasta_seq(fa)
        arr, upper, high, low, unresolved, missing = classify(seq)
        seqs[m] = arr
        uppers[m] = upper
        total = len(arr)
        stats[m] = {
            "high": high.sum() / total,
            "low": low.sum() / total,
            "unresolved": unresolved.sum() / total,
            "missing": missing.sum() / total,
        }

    methods = list(seqs.keys())
    if len(methods) < 2:
        print("Need at least 2 methods with output to compare.")
        sys.exit(1)

    fig_dir = out_base / "figures"
    fig_dir.mkdir(exist_ok=True)

    colors = {
        "high": "#2563eb",
        "low": "#60a5fa",
        "unresolved": "#f59e0b",
        "missing": "#d1d5db",
    }

    # ── Figure 1: Stacked bar chart ──────────────────────────
    fig, ax = plt.subplots(figsize=(max(4, len(methods) * 1.5), 5))
    x = np.arange(len(methods))
    bottoms = np.zeros(len(methods))

    for cat in ["high", "low", "unresolved", "missing"]:
        vals = [stats[m][cat] for m in methods]
        ax.bar(x, vals, bottom=bottoms, label=cat.capitalize(),
               color=colors[cat], width=0.6, edgecolor="white", linewidth=0.5)
        for i, v in enumerate(vals):
            if v > 0.03:
                ax.text(i, bottoms[i] + v / 2, f"{v:.1%}",
                        ha="center", va="center", fontsize=8,
                        color="white" if cat in ("high", "low") else "black")
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels([m.capitalize() for m in methods], fontsize=11)
    ax.set_ylabel("Fraction of positions", fontsize=11)
    ax.set_title(f"Confidence breakdown — {chrom}", fontsize=13, fontweight="bold")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path1 = fig_dir / "confidence_breakdown.png"
    fig.savefig(path1, dpi=150)
    plt.close(fig)
    print(f"  saved {path1}")

    # ── Figure 2: Pairwise agreement heatmap ─────────────────
    n = len(methods)
    agree_matrix = np.ones((n, n))
    shared_matrix = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(n):
            a, b = uppers[methods[i]], uppers[methods[j]]
            both_called = np.isin(a, list("ACGT")) & np.isin(b, list("ACGT"))
            shared = both_called.sum()
            shared_matrix[i, j] = shared
            if shared > 0:
                agree_matrix[i, j] = (a[both_called] == b[both_called]).sum() / shared

    fig, ax = plt.subplots(figsize=(max(4, n * 1.2 + 1), max(4, n * 1.2)))
    cmap = LinearSegmentedColormap.from_list("agree", ["#ef4444", "#fbbf24", "#22c55e"], N=256)
    vmin = max(0.95, agree_matrix.min() - 0.005)
    im = ax.imshow(agree_matrix, cmap=cmap, vmin=vmin, vmax=1.0, aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    labels = [m.capitalize() for m in methods]
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_yticklabels(labels, fontsize=11)

    for i in range(n):
        for j in range(n):
            val = agree_matrix[i, j]
            txt = f"{val:.3%}" if i != j else "—"
            ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                    color="white" if val < (vmin + 1) / 2 else "black")

    ax.set_title(f"Pairwise agreement — {chrom}", fontsize=13, fontweight="bold")
    fig.colorbar(im, ax=ax, label="Agreement rate", shrink=0.8)
    fig.tight_layout()
    path2 = fig_dir / "pairwise_agreement.png"
    fig.savefig(path2, dpi=150)
    plt.close(fig)
    print(f"  saved {path2}")

    # ── Figure 3: Sliding-window disagreement ────────────────
    length = min(len(seqs[m]) for m in methods)
    window = max(10_000, length // 500)

    fig, axes = plt.subplots(n - 1, 1, figsize=(12, max(3, 2.5 * (n - 1))),
                             sharex=True, squeeze=False)

    positions = np.arange(length)
    for idx in range(1, n):
        ax = axes[idx - 1, 0]
        a = uppers[methods[0]][:length]
        b = uppers[methods[idx]][:length]
        both = np.isin(a, list("ACGT")) & np.isin(b, list("ACGT"))
        disagree = both & (a != b)

        kernel = np.ones(window) / window
        rate = np.convolve(disagree.astype(float), kernel, mode="same")
        coverage = np.convolve(both.astype(float), kernel, mode="same")
        rate_safe = np.divide(rate, coverage, out=np.zeros_like(rate),
                              where=coverage > 0.01)

        pos_mb = positions / 1e6
        ax.fill_between(pos_mb, rate_safe, alpha=0.5, color="#ef4444")
        ax.plot(pos_mb, rate_safe, linewidth=0.5, color="#b91c1c")
        ax.set_ylabel("Disagree rate", fontsize=9)
        ax.set_title(f"{methods[0].capitalize()} vs {methods[idx].capitalize()}",
                     fontsize=10, loc="left")
        ax.set_ylim(0, max(0.01, rate_safe.max() * 1.3))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    axes[-1, 0].set_xlabel(f"Position on {chrom} (Mb)", fontsize=11)
    fig.suptitle(f"Sliding-window disagreement ({window // 1000}kb windows) — {chrom}",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    path3 = fig_dir / "disagreement_windows.png"
    fig.savefig(path3, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path3}")

    # ── Figure 4: Venn-style unique calls ────────────────────
    if n == 3:
        fig, ax = plt.subplots(figsize=(7, 5))

        a_arr = uppers[methods[0]][:length]
        b_arr = uppers[methods[1]][:length]
        c_arr = uppers[methods[2]][:length]
        a_called = np.isin(a_arr, list("ACGT"))
        b_called = np.isin(b_arr, list("ACGT"))
        c_called = np.isin(c_arr, list("ACGT"))

        categories = {
            f"Only {methods[0]}": a_called & ~b_called & ~c_called,
            f"Only {methods[1]}": ~a_called & b_called & ~c_called,
            f"Only {methods[2]}": ~a_called & ~b_called & c_called,
            f"{methods[0]}+{methods[1]}": a_called & b_called & ~c_called,
            f"{methods[0]}+{methods[2]}": a_called & ~b_called & c_called,
            f"{methods[1]}+{methods[2]}": ~a_called & b_called & c_called,
            "All three": a_called & b_called & c_called,
            "None": ~a_called & ~b_called & ~c_called,
        }

        labels = list(categories.keys())
        counts = [int(v.sum()) for v in categories.values()]
        total = sum(counts)

        bar_colors = ["#2563eb", "#22c55e", "#f59e0b", "#8b5cf6",
                      "#06b6d4", "#ec4899", "#374151", "#d1d5db"]
        non_zero = [(l, c, col) for l, c, col in zip(labels, counts, bar_colors) if c > 0]
        if non_zero:
            labs, cnts, cols = zip(*non_zero)
            y = np.arange(len(labs))
            ax.barh(y, cnts, color=cols, edgecolor="white", height=0.6)
            for i, c in enumerate(cnts):
                ax.text(c + total * 0.005, i, f"{c:,} ({c/total:.1%})",
                        va="center", fontsize=9)
            ax.set_yticks(y)
            ax.set_yticklabels(labs, fontsize=10)
            ax.set_xlabel("Number of positions", fontsize=11)
            ax.set_title(f"Called-site overlap — {chrom}", fontsize=13, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.invert_yaxis()

        fig.tight_layout()
        path4 = fig_dir / "called_site_overlap.png"
        fig.savefig(path4, dpi=150)
        plt.close(fig)
        print(f"  saved {path4}")

    print(f"\nAll figures saved to {fig_dir}/")


if __name__ == "__main__":
    main()
