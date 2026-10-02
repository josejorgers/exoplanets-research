"""One entry point for the whole pipeline.

Usage:
    uv run python -m src.cli run --sector 104 --n-jobs 8            # every step, in order
    uv run python -m src.cli run --sector 104 --from vet             # resume from a step
    uv run python -m src.cli run --sector 104 --from search --to vet
    uv run python -m src.cli <step> [that step's own options]        # one step, e.g.:
    uv run python -m src.cli vet --rethreshold --snr-min 6
    uv run python -m src.cli inject --n 1000 --n-jobs 10

Steps, in pipeline order: targets, fetch, detrend, search, vet, known, report
(plus inject, which is not part of `run`: it is a calibration of steps 3-5).

Each step skips stars whose output already exists unless it is given --force (`run --force`
passes it to targets, fetch, detrend and search). Everything logged by the main process also
goes to data/logs/sector_{N}.log; `run` records package versions in
data/logs/sector_{N}_versions.json. (Lines logged inside joblib worker processes only go to
the console.)
"""

from __future__ import annotations

import argparse
import importlib
import json
import logging
import os
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import dotenv

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
STEPS = ["targets", "fetch", "detrend", "search", "vet", "known", "report"]
EXTRA = ["inject"]
PACKAGES = ["numpy", "pandas", "astropy", "astroquery", "lightkurve", "wotan", "transitleastsquares",
            "tess-point", "batman-package", "scipy", "matplotlib"]

log = logging.getLogger("cli")


def setup_logging(sector: int | None) -> Path | None:
    """Console + file logging on the root logger (the modules' own basicConfig then does nothing)."""
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s [%(name)s] %(message)s")
    if not root.handlers:
        h = logging.StreamHandler()
        h.setFormatter(fmt)
        root.addHandler(h)
    path = None
    if sector:
        path = DATA_DIR / "logs" / f"sector_{sector}.log"
        path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(path)
        fh.setFormatter(fmt)
        root.addHandler(fh)
    logging.getLogger("lightkurve").setLevel(logging.ERROR)
    return path


def record_versions(sector: int) -> None:
    v = {}
    for pkg in PACKAGES:
        try:
            v[pkg] = version(pkg)
        except PackageNotFoundError:
            v[pkg] = None
    v["python"] = sys.version.split()[0]
    v["run_started"] = time.strftime("%Y-%m-%d %H:%M:%S")
    path = DATA_DIR / "logs" / f"sector_{sector}_versions.json"
    path.write_text(json.dumps(v, indent=2))
    log.info("package versions: %s", v)


def run_step(step: str, argv: list[str]) -> None:
    """Run src.<step>.main() as if called from the command line with `argv`."""
    module = importlib.import_module(f"src.{step}")
    old = sys.argv
    sys.argv = [f"src.{step}"] + argv
    t0 = time.time()
    log.info("=== %s %s", step, " ".join(argv))
    try:
        module.main()
    finally:
        sys.argv = old
    log.info("=== %s done in %.0f s", step, time.time() - t0)


def step_args(step: str, sector: int, n_jobs: int, force: bool) -> list[str]:
    """The options `run` passes to each step."""
    args = ["--sector", str(sector)]
    if step in ("fetch",):
        args += ["--n-jobs", str(min(n_jobs, 4))]  # MAST rate-limits
    elif step in ("detrend", "search", "vet"):
        args += ["--n-jobs", str(n_jobs)]
    if force and step in ("targets", "fetch", "detrend", "search"):
        args.append("--force")
    return args


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] in STEPS + EXTRA:
        step, rest = sys.argv[1], sys.argv[2:]
        sector = None
        if "--sector" in rest:
            sector = int(rest[rest.index("--sector") + 1])
        setup_logging(sector or int(os.environ.get("SECTOR", 0)) or None)
        run_step(step, rest)
        return

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the pipeline steps in order")
    r.add_argument("--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None)
    r.add_argument("--n-jobs", type=int, default=8)
    r.add_argument("--from", dest="start", choices=STEPS, default=STEPS[0])
    r.add_argument("--to", dest="stop", choices=STEPS, default=STEPS[-1])
    r.add_argument("--force", action="store_true", help="redo targets/fetch/detrend/search for every star")
    for s in STEPS + EXTRA:
        sub.add_parser(s, help=f"run src.{s} (all further options are passed to it)", add_help=False)
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    log_path = setup_logging(a.sector)
    record_versions(a.sector)
    steps = STEPS[STEPS.index(a.start): STEPS.index(a.stop) + 1]
    log.info("sector %d: running %s (log: %s)", a.sector, " -> ".join(steps), log_path)
    t0 = time.time()
    for s in steps:
        run_step(s, step_args(s, a.sector, a.n_jobs, a.force))
    log.info("pipeline finished in %.1f min", (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
