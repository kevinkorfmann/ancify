"""Command-line interface for ancify.

Usage::

    ancify init [-o config.yaml]         # generate a template config
    ancify project -c config.yaml        # Phase 1: project alignments
    ancify call -c config.yaml           # Phase 2: call ancestral states
    ancify evaluate -c config.yaml       # Phase 3: evaluate calls
    ancify run -c config.yaml            # run all phases
    ancify train -c config.yaml [-o model.lgb]  # train ML model

The package can also be invoked as ``python -m ancify <cmd> ...``.
"""

import argparse
import logging
import sys
import textwrap

from .config import load_config

EXAMPLE_CONFIG = textwrap.dedent("""\
    # ── Ancestral Allele Polarization Pipeline ──────────────────────
    # Adapt this template for your focal species and outgroups.

    # Informal label for the focal species (used only in log messages).
    focal_species: my_species

    # Tab-separated file with at least two columns:
    #   chromosome_name <TAB> length
    # Additional columns are ignored.
    chromosome_lengths: chromLens.txt

    # Optional: restrict to specific chromosomes.
    # If omitted, every chromosome in the lengths file is processed.
    # chromosomes:
    #   - chr1
    #   - chr2
    #   - chrX

    outgroups:
      # Inner outgroup: one or more closely related species.
      # A majority vote determines the inner consensus.
      inner:
        - name: close_species_1
          alignment: focal.closeSpecies1.net.axt.gz
        - name: close_species_2
          alignment: focal.closeSpecies2.net.axt.gz

      # Outer outgroup: one or more distantly related species.
      # Serves as an independent confirmation of the inner call.
      outer:
        - name: distant_species
          alignment: focal.distantSpecies.net.axt.gz

    # Working directory for intermediate projected sequences.
    work_dir: .

    # Directory for final ancestral FASTA files.
    output_dir: ./ancestral_calls

    # Minimum allele count to accept a majority-vote consensus.
    min_inner_freq: 1
    min_outer_freq: 1

    # Number of parallel worker processes.
    num_cpus: 4

    # ── Ancestral inference method ────────────────────────────────────
    # "voting" (default): two-tier inner/outer outgroup voting.
    # "parsimony": Fitch parsimony on a phylogenetic tree.
    #
    # method: parsimony
    #
    # Newick tree topology (required when method is "parsimony").
    # Can be an inline Newick string or a path to a .nwk file.
    # Leaf names must match outgroup 'name' fields.
    # tree: "(((close_species_1,close_species_2),distant_species),outgroup2)"

    # Optional evaluation block.  Both sub-sections are independent;
    # include either, both, or neither.
    #
    # Patterns support two placeholders:
    #   {chrom}    -- full chromosome name (e.g. "chr1")
    #   {chrom_id} -- name without leading "chr" (e.g. "1")
    #
    # evaluation:
    #   reference_dir: ./reference_ancestral/
    #   reference_pattern: "ancestor_{chrom_id}.fa"
    #   vcf_dir: ./vcf/
    #   vcf_pattern: "variants.{chrom}.vcf.gz"
""")


def _setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def cmd_init(args):
    out = args.output or "config.yaml"
    with open(out, "w") as f:
        f.write(EXAMPLE_CONFIG)
    print(f"Template configuration written to {out}")


def cmd_project(args):
    cfg = load_config(args.config)
    if args.num_cpus:
        cfg.num_cpus = args.num_cpus
    from .project import run_projection
    run_projection(cfg)


def cmd_call(args):
    cfg = load_config(args.config)
    if args.num_cpus:
        cfg.num_cpus = args.num_cpus
    from .ancestral import run_ancestral_calling
    run_ancestral_calling(cfg)


def cmd_evaluate(args):
    cfg = load_config(args.config)
    if args.num_cpus:
        cfg.num_cpus = args.num_cpus
    from .evaluate import run_evaluation
    run_evaluation(cfg)


def cmd_train(args):
    cfg = load_config(args.config)
    if args.num_cpus:
        cfg.num_cpus = args.num_cpus
    if args.output:
        cfg.ml_model_path = args.output
    from .ml import train_from_config
    out_path = train_from_config(cfg)
    print(f"Trained ML model saved to {out_path}")


def cmd_run(args):
    cfg = load_config(args.config)
    if args.num_cpus:
        cfg.num_cpus = args.num_cpus

    from .project import run_projection
    from .ancestral import run_ancestral_calling

    run_projection(cfg)
    run_ancestral_calling(cfg)

    if cfg.evaluation:
        from .evaluate import run_evaluation
        run_evaluation(cfg)


def main():
    parser = argparse.ArgumentParser(
        prog="ancify",
        description="Ancestral allele polarization using outgroup species.",
    )
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="enable debug logging")

    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("init", help="generate a template config file")
    p.add_argument("-o", "--output", help="output path (default: config.yaml)")
    p.set_defaults(func=cmd_init)

    for name, helptext, func in [
        ("project",  "Phase 1: project outgroup alignments", cmd_project),
        ("call",     "Phase 2: call ancestral states",       cmd_call),
        ("evaluate", "Phase 3: evaluate ancestral calls",    cmd_evaluate),
        ("run",      "run all phases end-to-end",            cmd_run),
    ]:
        p = sub.add_parser(name, help=helptext)
        p.add_argument("-c", "--config", required=True,
                       help="path to YAML configuration file")
        p.add_argument("-n", "--num-cpus", type=int,
                       help="override num_cpus from config")
        p.set_defaults(func=func)

    p = sub.add_parser("train", help="train an ML model for ancestral calling")
    p.add_argument("-c", "--config", required=True,
                   help="path to YAML configuration file")
    p.add_argument("-o", "--output",
                   help="output path for trained model (default: from config or work_dir)")
    p.add_argument("-n", "--num-cpus", type=int,
                   help="override num_cpus from config")
    p.set_defaults(func=cmd_train)

    args = parser.parse_args()
    _setup_logging(getattr(args, "verbose", False))

    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)
