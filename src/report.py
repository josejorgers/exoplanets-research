"""Step 7 — the report: one page per candidate, plus a run summary, in one PDF.

Usage:
    uv run python -m src.report                          # sector from .env
    uv run python -m src.report --tics 87329149,23817315 # add specific stars (e.g. for review)
    uv run python -m src.report --all-known              # also pages for known/blend signals (default: yes)

Input : data/search/sector_{N}_known.parquet  (from src.known: vet metrics + known_class)
        data/targets/sector_{N}.parquet         (stellar parameters)
        data/targets/sector_{N}_manifest.parquet (light-curve author)
        data/lc_flat/{tic}.parquet               (raw + flattened light curves)
        data/search/sector_{N}_vetted_meta.json  (momentum dumps)
Output: data/reports/sector_{N}.pdf             summary page + one page per signal, new candidates first
        data/reports/sector_{N}/{rank:02d}_{tic}.png  the same pages as PNG (for quick viewing)
        data/reports/sector_{N}_candidates.csv  one row per page: ephemeris, derived planet properties, flags
        data/reports/sector_{N}_summary.json    counts at every step of the pipeline

TLS is re-run for each reported star (periodogram and model are not stored by src.search);
results are cached in data/reports/cache/tls_{tic}.npz.

How to read a page: notebooks/pipeline-7-report-walkthrough.ipynb
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import warnings
from pathlib import Path

import dotenv
import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from transitleastsquares import cleaned_array, transitleastsquares  # noqa: E402

from src.known import CLASS_ORDER  # noqa: E402
from src.search import tls_params  # noqa: E402
from src.vet import HARD, THRESH, _phase_time  # noqa: E402

log = logging.getLogger("report")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
REPORT_DIR = DATA_DIR / "reports"

G = 6.674e-11
R_SUN = 6.957e8
M_SUN = 1.989e30
AU = 1.496e11
DAY = 86400.0
T_SUN = 5772.0


# ------------------------------------------------------------------------------------------------
# derived planet properties
# ------------------------------------------------------------------------------------------------


def derived(period_d: float, rad: float, mass: float, teff: float) -> dict:
    """Orbit and irradiation of a planet on a circular orbit around a star of rad/mass (solar), teff (K).

    a from Kepler's third law; S = L/(4 pi a^2) in Earth units, L/L_sun = R^2 (T/T_sun)^4;
    T_eq = T_eff sqrt(R_star / 2a) for zero albedo and full heat redistribution.
    """
    if not all(np.isfinite([period_d, rad, mass, teff])):
        return {"a_au": np.nan, "a_rs": np.nan, "insolation_earth": np.nan, "teq_k": np.nan}
    a = (G * mass * M_SUN * (period_d * DAY) ** 2 / (4 * np.pi**2)) ** (1 / 3)
    lum = rad**2 * (teff / T_SUN) ** 4
    return {
        "a_au": a / AU,
        "a_rs": a / (rad * R_SUN),
        "insolation_earth": lum / (a / AU) ** 2,
        "teq_k": teff * np.sqrt(rad * R_SUN / (2 * a)),
    }


# ------------------------------------------------------------------------------------------------
# TLS periodogram and model (cached)
# ------------------------------------------------------------------------------------------------


def tls_cached(tic_id: int, flat_file: str, rad: float, mass: float, cache_dir: Path) -> dict:
    """Periodogram + model of the best TLS signal, same parameters as src.search, cached as npz."""
    path = cache_dir / f"tls_{tic_id}.npz"
    if path.exists():
        return dict(np.load(path))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lc = pd.read_parquet(flat_file, columns=["time", "flux", "flux_err"])
        t, y, dy = cleaned_array(lc["time"].values, lc["flux"].values, lc["flux_err"].values)
        r = transitleastsquares(t, y, dy).power(**tls_params(rad, mass, t.max() - t.min()))
    out = {
        "periods": np.asarray(r.periods), "power": np.asarray(r.power), "period": float(r.period),
        "model_time": np.asarray(r.model_lightcurve_time), "model_flux": np.asarray(r.model_lightcurve_model),
        "fold_phase": np.asarray(r.model_folded_phase), "fold_model": np.asarray(r.model_folded_model),
    }
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **out)
    return out


# ------------------------------------------------------------------------------------------------
# panels
# ------------------------------------------------------------------------------------------------


def _binned(x, y, edges):
    idx = np.digitize(x, edges)
    vals = np.array([np.mean(y[idx == i]) if np.any(idx == i) else np.nan for i in range(1, len(edges))])
    return 0.5 * (edges[1:] + edges[:-1]), vals


def _ylim_for(y_in, depth):
    lo = min(np.nanpercentile(y_in, 2) if len(y_in) else 1, 1 - depth)
    return lo - 0.4 * (1 - lo), 1 + 0.5 * (1 - lo)


def panel_raw(ax, lc):
    ax.plot(lc["time"], lc["flux_raw"], ",", color="0.55")
    ax.plot(lc["time"], lc["trend"], "-", color="C3", lw=0.8)
    ax.set_ylabel("raw flux")
    ax.set_title("1. raw light curve (grey) and the biweight trend (red)", fontsize=8, loc="left")


def panel_flat(ax, lc, s, dumps):
    ax.plot(lc["time"], lc["flux"], ",k")
    n0 = int(np.floor((lc["time"].min() - s.T0) / s.period))
    n1 = int(np.ceil((lc["time"].max() - s.T0) / s.period))
    for n in range(n0, n1 + 1):
        ax.axvline(s.T0 + n * s.period, ymin=0.9, ymax=1.0, color="C3", lw=0.8)
    for td in dumps:
        ax.axvline(td, color="C0", ls="--", lw=0.8)
    ax.set_ylabel("flattened")
    ax.set_xlabel("time [BTJD]", fontsize=7)
    ax.set_title("2. flattened light curve; red ticks = predicted transits, blue dashed = momentum dumps",
                 fontsize=8, loc="left")


def panel_fold(ax, lc, s, tls):
    x = _phase_time(lc["time"].values, s.period, s.T0) * 24
    w = max(3 * s.duration_h, 1.5)
    sel = np.abs(x) < w
    ax.plot(x[sel], lc["flux"].values[sel], ".", ms=1, color="0.6")
    bx, by = _binned(x[sel], lc["flux"].values[sel], np.linspace(-w, w, 40))
    ax.plot(bx, by, "o", ms=2.5, color="C0")
    if tls is not None:
        mx = (tls["fold_phase"] - 0.5) * tls["period"] * 24
        m = np.abs(mx) < w
        ax.plot(mx[m], tls["fold_model"][m], "-", color="C3", lw=1)
    ax.set_xlim(-w, w)
    ax.set_ylim(*_ylim_for(lc["flux"].values[np.abs(x) < s.duration_h / 2], s.depth_ppm / 1e6))
    ax.set_xlabel("hours from mid-transit", fontsize=7)
    ax.set_title("3. fold: binned (blue),\nTLS model (red)", fontsize=8, loc="left")


def panel_odd_even(ax, lc, s):
    t, y = lc["time"].values, lc["flux"].values
    x = _phase_time(t, s.period, s.T0) * 24
    epoch = np.round((t - s.T0) / s.period).astype(int)
    w = max(3 * s.duration_h, 1.5)
    for parity, c, lab in [(1, "C0", "odd"), (0, "C3", "even")]:
        m = (epoch % 2 == parity) & (np.abs(x) < w)
        bx, by = _binned(x[m], y[m], np.linspace(-w, w, 30))
        ax.plot(bx, by, "o-", ms=2, lw=0.8, color=c, label=lab)
    ax.axhline(1 - s.depth_ppm / 1e6, color="0.5", ls=":", lw=0.8)
    ax.set_xlim(-w, w)
    ax.set_ylim(*_ylim_for(y[np.abs(x) < s.duration_h / 2], s.depth_ppm / 1e6))
    ax.legend(fontsize=6, loc="lower right")
    ax.set_xlabel("hours from mid-transit", fontsize=7)
    ax.set_title(f"4. odd (blue) vs even (red):\n{s.odd_even_sigma:.1f}σ apart", fontsize=8, loc="left")


def panel_secondary(ax, lc, s):
    """Zoom on phase 0.5 (the secondary eclipse of a circular orbit), same scale as the transit."""
    t, y = lc["time"].values, lc["flux"].values
    x = _phase_time(t - s.period / 2, s.period, s.T0) * 24  # hours from phase 0.5
    w = max(3 * s.duration_h, 1.5)
    sel = np.abs(x) < w
    ax.plot(x[sel], y[sel], ".", ms=1, color="0.6")
    bx, by = _binned(x[sel], y[sel], np.linspace(-w, w, 40))
    ax.plot(bx, by, "o", ms=2.5, color="C0")
    xt = _phase_time(t, s.period, s.T0) * 24
    ax.set_ylim(*_ylim_for(y[np.abs(xt) < s.duration_h / 2], s.depth_ppm / 1e6))
    ax.set_xlim(-w, w)
    ax.set_xlabel("hours from phase 0.5", fontsize=7)
    ax.set_title(f"5. phase 0.5 (worst dip: {s.secondary_snr:.1f}σ\nat phase {s.secondary_phase:.2f})",
                 fontsize=8, loc="left")


def panel_periodogram(ax, tls, s):
    if tls is None:
        ax.text(0.5, 0.5, "no TLS cache", ha="center", transform=ax.transAxes)
        return
    ax.plot(tls["periods"], tls["power"], "k-", lw=0.6)
    ax.axvline(s.period, color="C3", alpha=0.35, lw=3)
    for n in (2, 3):
        for p in (s.period * n, s.period / n):
            ax.axvline(p, color="C0", alpha=0.35, ls="--", lw=0.8)
    ax.axhline(THRESH["sde_min"], color="0.5", ls=":", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("period [d]", fontsize=7)
    ax.set_ylabel("SDE", fontsize=7)
    ax.set_title(f"6. TLS periodogram (red: best, blue: ×2, ×3, ÷2, ÷3; dotted: SDE {THRESH['sde_min']:.0f})",
                 fontsize=8, loc="left")


def panel_table(ax, s, d):
    ax.axis("off")
    rows = [
        ("odd/even [σ]", s.odd_even_sigma, "odd_even"), ("secondary [σ]", s.secondary_snr, "secondary"),
        ("dur / expected", s.duration_ratio, "duration"), ("Rp [R⊕]", s.rp_earth, "radius"),
        ("var amp / depth", s.var_amp_frac, "variability"), ("1-event survival", s.single_event_survival, "single_event"),
        ("dump survival", s.dump_survival, "momentum_dump"), ("edge survival", s.edge_survival, "data_edges"),
        ("SNR (TLS)", s.snr_tls, "snr"), ("SDE", s.SDE, "sde"),
    ]
    lines = ["hard tests", "-" * 28]
    for label, val, key in rows:
        mark = "✓" if s[f"pass_{key}"] else "✗"
        lines.append(f"{mark} {label:17s}{val:9.2f}")
    lines += ["", "soft flags", "-" * 28]
    lines.append(f"{'⚑' if s.flag_v_shaped else ' '} T23/T14       {s.trap_flat_frac:9.2f}")
    lines.append(f"{'⚑' if s.flag_centroid else ' '} centroid [σ]  {s.centroid_sigma:9.1f}")
    lines.append(f"{'⚑' if s.flag_contamination else ' '} contratio     {s.contratio:9.2f}")
    lines.append(f"{'⚑' if s.flag_variable_star else ' '} var FAP       {s.var_fap:9.1e}")
    ax.text(0, 1.05, "\n".join(lines), va="top", ha="left", family="DejaVu Sans Mono", fontsize=6.6,
            transform=ax.transAxes)


def candidate_page(s, star, lc: pd.DataFrame, tls: dict | None, dumps, rank: int, sector: int, author: str):
    """One figure (US letter, portrait) for one signal. s = row of the known table, star = target row."""
    d = derived(s.period, star.rad, star.mass, star.teff)
    fig = plt.figure(figsize=(8.5, 11))
    gs = fig.add_gridspec(5, 3, height_ratios=[0.7, 1, 1, 1.25, 1.25], hspace=0.6, wspace=0.32,
                          left=0.07, right=0.97, top=0.97, bottom=0.04)
    head = fig.add_subplot(gs[0, :])
    head.axis("off")
    verdict = "PASSES all hard tests" if s.passes_all else f"FAILS: {s.failed}"
    known = s.known_class + (f" — {s.known_ids}" if s.known_ids else "") + (
        f" — neighbour: {s.nb_id} at {s.nb_sep_arcsec:.0f}\" (likely source: {s.likely_source})" if s.nb_id else "")
    head.text(0, 1, f"#{rank}  TIC {s.tic_id}   sector {sector}   {author}", fontsize=13, weight="bold", va="top")
    head.text(0, 0.62,
              f"P = {s.period:.5f} d   T0 = {s.T0:.4f} BTJD   depth = {s.depth_ppm:.0f} ppm   "
              f"duration = {s.duration_h:.2f} h   Rp = {s.rp_earth:.2f} R⊕   SNR = {s.snr:.1f}   SDE = {s.SDE:.1f}",
              fontsize=8.5, va="top")
    head.text(0, 0.36,
              f"star: Tmag {star.tmag:.2f}   Teff {star.teff:.0f} K   R {star.rad:.3f} R☉   M {star.mass:.3f} M☉   "
              f"d {star.dist_pc:.1f} pc   contratio {star.contratio:.2f}", fontsize=8.5, va="top")
    color = "C2" if s.passes_all and s.known_class in ("new", "new_in_known_system") else "C3" if not s.passes_all else "C1"
    head.text(0, 0.1, f"{verdict}   |   flags: {s['flags'] or '-'}   |   {known}", fontsize=8, va="top", color=color)
    head.text(0, -0.16,
              f"planet if circular: a = {d['a_au']:.4f} au = {d['a_rs']:.1f} R*   S = {d['insolation_earth']:.1f} S⊕   "
              f"Teq = {d['teq_k']:.0f} K   (β = {s.beta:.2f}, {s.n_transits_used} transits with data)",
              fontsize=8, va="top")
    head.text(1, 1, f"exofop.ipac.caltech.edu/tess/target.php?id={s.tic_id}", fontsize=7, ha="right", va="top",
              color="C0")

    panel_raw(fig.add_subplot(gs[1, :]), lc)
    panel_flat(fig.add_subplot(gs[2, :]), lc, s, dumps)
    panel_fold(fig.add_subplot(gs[3, 0]), lc, s, tls)
    panel_odd_even(fig.add_subplot(gs[3, 1]), lc, s)
    panel_secondary(fig.add_subplot(gs[3, 2]), lc, s)
    panel_periodogram(fig.add_subplot(gs[4, 0:2]), tls, s)
    panel_table(fig.add_subplot(gs[4, 2]), s, d)
    for ax in fig.axes:
        ax.tick_params(labelsize=6.5)
    return fig, d


# ------------------------------------------------------------------------------------------------
# run summary
# ------------------------------------------------------------------------------------------------


def run_summary(sector: int, known_tab: pd.DataFrame) -> dict:
    tg = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet", columns=["tic_id"])
    mf = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_manifest.parquet")
    dt = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_detrend.parquet")
    se = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}.parquet", columns=["status"])
    funnel, alive = {}, known_tab["status"] != "error"
    for k in HARD:
        alive &= known_tab[f"pass_{k}"]
        funnel[k] = int(alive.sum())
    cand = known_tab[known_tab["passes_all"]]
    return {
        "sector": sector,
        "targets": len(tg),
        "light_curves_downloaded": int((mf["status"] == "ok").sum()),
        "light_curves_by_author": {k: int(v) for k, v in mf[mf["status"] == "ok"]["author"].value_counts().items()},
        "detrended": int((dt["status"] == "ok").sum()),
        "searched": int((se["status"] == "ok").sum()),
        "vetting_funnel": funnel,
        "pass_vetting": len(cand),
        "classes_of_candidates": {k: int(v) for k, v in cand["known_class"].value_counts().items()},
        "new_candidates": int(cand["known_class"].isin(["new", "new_in_known_system"]).sum()),
        "new_without_flags": int((cand["known_class"].isin(["new", "new_in_known_system"]) & (cand["n_flags"] == 0)).sum()),
        "thresholds": THRESH,
    }


def summary_page(summary: dict, cand: pd.DataFrame):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.07, 0.96, f"TESS sector {summary['sector']} — M-dwarf transit search", fontsize=16, weight="bold")
    steps = [
        ("targets (TIC M dwarfs on silicon)", summary["targets"]),
        ("light curves downloaded", summary["light_curves_downloaded"]),
        ("detrended", summary["detrended"]),
        ("searched with TLS", summary["searched"]),
    ] + [(f"vetting: + {k}", v) for k, v in summary["vetting_funnel"].items()]
    ax = fig.add_axes([0.32, 0.52, 0.6, 0.4])
    labels, vals = zip(*steps)
    y = np.arange(len(steps))[::-1]
    ax.barh(y, vals, color=["0.55"] * 4 + ["C0"] * (len(steps) - 4))
    ax.set_xscale("log")
    ax.set_yticks(y, labels, fontsize=8)
    for yi, v in zip(y, vals):
        ax.text(v * 1.08, yi, f"{v:,}", va="center", fontsize=7.5)
    ax.set_xlabel("stars / signals left (log scale)", fontsize=8)
    ax.tick_params(labelsize=7.5)
    cls = summary["classes_of_candidates"]
    fig.text(0.07, 0.47, "Signals passing all hard vetting tests, by known-object class:", fontsize=10, weight="bold")
    fig.text(0.07, 0.445, "   ".join(f"{k}: {cls.get(k, 0)}" for k in CLASS_ORDER if cls.get(k, 0)), fontsize=9)
    rows = ["rank  TIC          class          P [d]    depth[ppm]  Rp[R⊕]   SNR   flags"]
    for r in cand.itertuples():
        rows.append(f"{r.rank:4d}  {r.tic_id:<11d}  {r.known_class:<13s} {r.period:8.4f} {r.depth_ppm:10.0f}  "
                    f"{r.rp_earth:6.2f} {r.snr:6.1f}   {r.flags or '-'}")
    fig.text(0.07, 0.42, "\n".join(rows[:60]), va="top", family="DejaVu Sans Mono", fontsize=5.6)
    return fig


# ------------------------------------------------------------------------------------------------
# sector
# ------------------------------------------------------------------------------------------------


def report_sector(sector: int, extra_tics: list[int] | None = None, n_jobs: int = 4) -> pd.DataFrame:
    tab = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}_known.parquet")
    targets = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet").set_index("tic_id")
    manifest = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_manifest.parquet").set_index("tic_id")
    detrend = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_detrend.parquet").set_index("tic_id")
    meta_path = DATA_DIR / "search" / f"sector_{sector}_vetted_meta.json"
    dumps = json.loads(meta_path.read_text())["momentum_dumps_btjd"] if meta_path.exists() else []

    sel = tab["passes_all"] | tab["tic_id"].isin(extra_tics or [])
    pages = tab[sel].copy()
    pages["_order"] = pages["known_class"].map(CLASS_ORDER.index)
    # new candidates first, then by class; within a class passing ones first, fewest flags, highest SNR
    pages = pages.sort_values(["_order", "passes_all", "n_flags", "snr"], ascending=[True, False, True, False])
    pages = pages.drop(columns="_order").reset_index(drop=True)
    pages.insert(0, "rank", np.arange(1, len(pages) + 1))
    log.info("sector %d: %d pages (%d passing vetting, %d extra)", sector, len(pages), int(pages["passes_all"].sum()),
             int((~pages["passes_all"]).sum()))

    cache = REPORT_DIR / "cache"
    tls = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(tls_cached)(int(r.tic_id), detrend.loc[r.tic_id, "file_path"],
                                   float(targets.loc[r.tic_id, "rad"]), float(targets.loc[r.tic_id, "mass"]), cache)
        for r in pages.itertuples()
    )

    png_dir = REPORT_DIR / f"sector_{sector}"
    png_dir.mkdir(parents=True, exist_ok=True)
    summary = run_summary(sector, tab)
    pdf_path = REPORT_DIR / f"sector_{sector}.pdf"
    rows = []
    with PdfPages(pdf_path) as pdf:
        fig = summary_page(summary, pages)
        pdf.savefig(fig)
        plt.close(fig)
        for (_, s), t in zip(pages.iterrows(), tls):
            star = targets.loc[s.tic_id]
            lc = pd.read_parquet(detrend.loc[s.tic_id, "file_path"])
            fig, d = candidate_page(s, star, lc, t, dumps, int(s["rank"]), sector, manifest.loc[s.tic_id, "author"])
            pdf.savefig(fig)
            fig.savefig(png_dir / f"{int(s['rank']):02d}_{s.tic_id}.png", dpi=110)
            plt.close(fig)
            rows.append({**{k: s[k] for k in ["rank", "tic_id", "known_class", "passes_all", "failed", "flags", "period",
                                               "T0", "duration_h", "depth_ppm", "rp_earth", "snr", "SDE", "known_ids",
                                               "nb_id", "likely_source"]},
                         "tmag": star.tmag, "teff": star.teff, "rad": star.rad, "mass": star.mass,
                         "dist_pc": star.dist_pc, **d,
                         "exofop": f"https://exofop.ipac.caltech.edu/tess/target.php?id={s.tic_id}"})
    out = pd.DataFrame(rows)
    out.to_csv(REPORT_DIR / f"sector_{sector}_candidates.csv", index=False)
    (REPORT_DIR / f"sector_{sector}_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("wrote %s (%d pages), %s, PNGs in %s", pdf_path, len(pages) + 1,
             REPORT_DIR / f"sector_{sector}_candidates.csv", png_dir)
    log.info("summary: %s", {k: v for k, v in summary.items() if k != "thresholds"})
    return out


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument("--tics", type=str, default="", help="comma-separated extra TIC IDs to include")
    p.add_argument("--n-jobs", type=int, default=4, help="parallel TLS re-runs for uncached stars")
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    extra = [int(x) for x in a.tics.split(",") if x.strip()]
    out = report_sector(a.sector, extra, a.n_jobs)
    print(out[["rank", "tic_id", "known_class", "period", "depth_ppm", "rp_earth", "teq_k", "snr", "flags"]].to_string())


if __name__ == "__main__":
    main()
