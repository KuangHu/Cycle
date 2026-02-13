#!/usr/bin/env python3
"""Install tldr and exonerate into the current conda environment.

tldr requires Python >= 3.6, plus external tools: minimap2, samtools,
mafft, and exonerate.  This script installs exonerate via conda and
tldr via pip (from GitHub).

Usage:
    python scripts/install_tldr.py
    python scripts/install_tldr.py --env-name opfi   # explicit env name

After installation, tldr is available directly in the active env:
    tldr -b sample.sorted.bam -e is_reference.fa -r reference.fa
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys

logger = logging.getLogger(__name__)

TLDR_REPO = "https://github.com/adamewing/tldr.git"
REQUIRED_TOOLS = ["minimap2", "samtools", "mafft"]


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a command, logging it and checking for errors."""
    logger.info(f"  $ {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, text=True, **kwargs)


def conda_exe() -> str:
    """Find conda or mamba executable."""
    for name in ("mamba", "conda"):
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError("Neither conda nor mamba found in PATH")


def current_env_name() -> str:
    """Get the name of the currently active conda env."""
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix:
        return os.path.basename(prefix)
    return "base"


def main():
    parser = argparse.ArgumentParser(
        description="Install tldr + exonerate into the current conda environment.",
    )
    parser.add_argument(
        "--env-name", default=None,
        help="Conda env to install into. Default: current active env.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Reinstall even if tldr is already present.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    conda = conda_exe()
    env_name = args.env_name or current_env_name()
    logger.info(f"Target env: {env_name}")

    # ── Check prerequisites ───────────────────────────────────────────
    missing = [t for t in REQUIRED_TOOLS if not shutil.which(t)]
    if missing:
        logger.warning(
            f"Missing tools (expected in PATH): {missing}. "
            f"tldr will need these at runtime."
        )

    # ── Step 1: Install exonerate via conda ───────────────────────────
    if shutil.which("exonerate") and not args.force:
        logger.info("exonerate already installed, skipping")
    else:
        logger.info("Installing exonerate")
        run([
            conda, "install", "-n", env_name,
            "exonerate",
            "-c", "bioconda", "-c", "conda-forge",
            "-y",
        ])

    # ── Step 2: Install tldr via pip ──────────────────────────────────
    if shutil.which("tldr") and not args.force:
        logger.info("tldr already installed, skipping")
    else:
        logger.info("Installing tldr from GitHub")
        run([
            "pip", "install",
            f"git+{TLDR_REPO}",
        ])

    # ── Step 3: Verify ────────────────────────────────────────────────
    logger.info("Verifying installation")
    all_ok = True

    for tool, cmd in [("tldr", ["tldr", "--version"]), ("exonerate", ["exonerate", "--version"])]:
        try:
            ret = subprocess.run(cmd, capture_output=True, text=True)
            if ret.returncode == 0:
                version = (ret.stdout or ret.stderr).strip().splitlines()[0]
                logger.info(f"  {tool}: {version}")
            else:
                logger.error(f"  {tool} check failed")
                all_ok = False
        except Exception as e:
            logger.error(f"  {tool} verification failed: {e}")
            all_ok = False

    if not all_ok:
        logger.error("Installation incomplete — see errors above")
        sys.exit(1)

    logger.info(
        f"\nInstallation complete in '{env_name}' env. Run tldr directly:\n"
        f"  tldr -b sample.sorted.bam -e is_reference.fa -r reference.fa"
    )


if __name__ == "__main__":
    main()
