#!/usr/bin/env python3
"""Generate sbatch commands for all batch TSVs.

Produces a shell script where each line is an ``sbatch --wrap`` command
that runs the full pipeline on one batch TSV.

Usage:
    python scripts/generate_sbatch.py \
        --batch-dir data/batches \
        --root-outdir /groups/rubin/projects/kuang/out/IS_cycle \
        -o submit_all.sh

    # Then review and submit:
    bash submit_all.sh
"""

import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Generate sbatch --wrap commands for batch pipeline runs.",
    )
    parser.add_argument(
        "--batch-dir", default="data/batches",
        help="Directory containing batch_NNN.tsv files.",
    )
    parser.add_argument(
        "--root-outdir", default="/groups/rubin/projects/kuang/out/IS_cycle",
        help="Root output directory on shared storage.",
    )
    parser.add_argument(
        "--pipeline-script", default=None,
        help="Path to run_pipeline.py. Default: auto-detect from this script's location.",
    )
    parser.add_argument(
        "--conda-env", default="opfi",
        help="Conda environment name. Default: opfi",
    )
    parser.add_argument(
        "--partition", default="standard",
        help="SLURM partition. Default: standard",
    )
    parser.add_argument(
        "--qos", default="standard",
        help="SLURM QOS. Default: standard",
    )
    parser.add_argument(
        "--cpus", type=int, default=48,
        help="CPUs per job. Default: 48 (full node)",
    )
    parser.add_argument(
        "--mem", default="192G",
        help="Memory per job. Default: 192G (full node)",
    )
    parser.add_argument(
        "--time", default="2-00:00:00",
        help="Wall time limit. Default: 2-00:00:00",
    )
    parser.add_argument(
        "--threads", type=int, default=44,
        help="Threads for minimap2/samtools inside pipeline. Default: 44",
    )
    parser.add_argument(
        "--sort-memory", default="4G",
        help="Memory per samtools sort thread. Default: 4G",
    )
    parser.add_argument(
        "--steps", default=None,
        help="Pipeline steps to run (space-separated). Default: all",
    )
    parser.add_argument(
        "-o", "--output", default="submit_all.sh",
        help="Output shell script. Default: submit_all.sh",
    )
    args = parser.parse_args()

    # Locate pipeline script
    if args.pipeline_script:
        pipeline = Path(args.pipeline_script).resolve()
    else:
        pipeline = (Path(__file__).resolve().parent / "run_pipeline.py")
    if not pipeline.exists():
        print(f"ERROR: pipeline script not found: {pipeline}", file=sys.stderr)
        sys.exit(1)

    # Find batch TSVs
    batch_dir = Path(args.batch_dir).resolve()
    batch_files = sorted(batch_dir.glob("batch_*.tsv"))
    if not batch_files:
        print(f"ERROR: no batch_*.tsv found in {batch_dir}", file=sys.stderr)
        sys.exit(1)

    root_outdir = Path(args.root_outdir)
    log_dir = root_outdir / "logs"

    lines = [
        "#!/bin/bash",
        f"# Auto-generated sbatch commands for {len(batch_files)} batches",
        f"# Root output: {root_outdir}",
        f"# Pipeline: {pipeline}",
        "",
        f"mkdir -p {log_dir}",
        "",
    ]

    for batch_tsv in batch_files:
        batch_name = batch_tsv.stem  # e.g. batch_000
        outdir = root_outdir / batch_name
        log_file = log_dir / f"{batch_name}_%j.log"

        # Build pipeline command
        pipeline_cmd = (
            f"source activate {args.conda_env} && "
            f"python {pipeline}"
            f" --metadata {batch_tsv}"
            f" --outdir {outdir}"
            f" --threads {args.threads}"
            f" --sort-memory {args.sort_memory}"
        )
        if args.steps:
            pipeline_cmd += f" --steps {args.steps}"

        sbatch_cmd = (
            f"sbatch"
            f" --job-name={batch_name}"
            f" --partition={args.partition}"
            f" --qos={args.qos}"
            f" --cpus-per-task={args.cpus}"
            f" --mem={args.mem}"
            f" --time={args.time}"
            f" --output={log_file}"
            f" --error={log_file}"
            f" --wrap '{pipeline_cmd}'"
        )
        lines.append(sbatch_cmd)

    lines.append("")

    output_path = Path(args.output)
    output_path.write_text("\n".join(lines))
    print(f"Wrote {len(batch_files)} sbatch commands to {output_path}")
    print(f"Review, then run: bash {output_path}")


if __name__ == "__main__":
    main()
