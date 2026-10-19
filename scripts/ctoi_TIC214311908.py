"""Summary plots and CTOI parameters for TIC 214311908.

Reads the follow-up products in data/followup/214311908/ (transit fit, multi-sector light curve,
Gaia neighbours, pixel exclusions) and writes to docs/ctoi/TIC214311908/:

    fig1_lightcurve.png   fold of all six sectors, each sector alone, odd/even, phase 0.5, depth per sector
    fig2_pixels.png       PRF scene and difference image (sector 104), delta-chi2 of each star as the source
    fig3_neighbours.png   every Gaia DR3 star within 60": how deep it would have to eclipse, and why it is excluded
    pixel_dchi2.csv       PRF-fit delta-chi2 of each Gaia source against the target
    summary.pdf           the three figures, one per page
    ctoi_parameters.json  the numbers for the CTOI form

Usage: uv run --extra followup python scripts/ctoi_TIC214311908.py [--skip-pixels]
(the pixel figure downloads six TESScut cutouts unless they are in the lightkurve cache).
"""

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

import batman
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.followup import masked_flatten  # noqa: E402

warnings.filterwarnings("ignore")
TIC = 214311908
GAIA_TARGET = 5925373272850784768
FU = ROOT / "data" / "followup" / str(TIC)
OUT = ROOT / "docs" / "ctoi" / f"TIC{TIC}"
BTJD0 = 2457000.0
SECTOR_YEAR = {12: 2019, 39: 2021, 66: 2023, 93: 2025, 103: 2026, 104: 2026}
# Gaia DR3 eclipsing binaries listed by SIMBAD near the target (vari_eclipsing_binary, checked 2026-10-19)
KNOWN_EB = {5925373272827278208: 0.31608, 5925373204107739648: 0.26093}
NAMES = {GAIA_TARGET: "target", 5925373272850782208: "co-moving companion", 5925373272827253632: "companion's 0.4\" source",
         5925373268528684288: "red giant (LPV)", 5925373272827260544: "3.1\" star", 5925373272827262208: "5.7\" star",
         5925373272827256576: "12.0\" star"}


def load():
    s = json.loads((FU / "summary.json").read_text())
    star = pd.read_parquet(ROOT / "data" / "targets" / "sector_104.parquet").set_index("tic_id").loc[TIC]
    chain = np.load(FU / "fit_samples.npz")["chain"]  # T0, P, k, b, log10(a/R*), q1, q2
    return s, star, chain


def params(fit_row):
    t0, per, k, b, la, q1, q2 = fit_row
    pm = batman.TransitParams()
    a = 10**la
    sq = np.sqrt(q1)
    pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = t0, per, k, a, np.degrees(np.arccos(b / a)), 0.0, 90.0
    pm.limb_dark, pm.u = "quadratic", [2 * sq * q2, sq * (1 - 2 * q2)]
    return pm


def contact_times(chain):
    P, k, b, a = chain[:, 1], chain[:, 2], chain[:, 3], 10 ** chain[:, 4]
    inc = np.arccos(b / a)

    def t(x):
        return P / np.pi * np.arcsin(np.clip(np.sqrt(np.clip(x**2 - b**2, 0, None)) / (a * np.sin(inc)), 0, 1))
    return t(1 + k), t(1 - k)


def ctoi_parameters(s, star, chain):
    fit = s["fit"]
    med = np.array([fit[n]["median"] for n in ["T0", "P", "k", "b", "log10_a_rs", "q1", "q2"]])
    T14, _ = contact_times(chain)
    rng = np.random.default_rng(1)
    tt = np.linspace(-0.05, 0.05, 2001)
    depths = []
    for row in chain[rng.choice(len(chain), 400, replace=False)]:
        pm = params(row)
        pm.t0 = 0.0
        depths.append(1 - batman.TransitModel(pm, tt).light_curve(pm).min())
    depths = np.array(depths) * 1e6
    a_rs = 10 ** chain[:, 4]
    a_au = a_rs * star.rad * 6.957e8 / 1.496e11
    teq = star.teff * np.sqrt(1 / (2 * a_rs))
    lum = star.rad**2 * (star.teff / 5772) ** 4
    s_earth = lum / a_au**2

    def pct(x):
        lo, m, hi = np.percentile(x, [16, 50, 84])
        return dict(value=float(m), err=float((hi - lo) / 2))

    out = {
        "tic_id": TIC,
        "period_d": dict(value=fit["P"]["median"], err=max(fit["P"]["plus"], fit["P"]["minus"])),
        "epoch_bjd_tdb": dict(value=fit["T0"]["median"] + BTJD0, err=max(fit["T0"]["plus"], fit["T0"]["minus"])),
        "depth_ppm": pct(depths),
        "duration_h": pct(T14 * 24),
        "rp_earth": dict(value=fit["rp_earth"]["median"], err=max(fit["rp_earth"]["plus"], fit["rp_earth"]["minus"])),
        "rp_rs": dict(value=fit["k"]["median"], err=max(fit["k"]["plus"], fit["k"]["minus"])),
        "impact_parameter": dict(value=fit["b"]["median"], err=max(fit["b"]["plus"], fit["b"]["minus"])),
        "a_rs": pct(a_rs),
        "a_au": pct(a_au),
        "teq_k_albedo0": pct(teq),
        "insolation_earth": pct(s_earth),
        "stellar_density_ratio_fit_over_tic": s["verdict"]["density_ratio"],
        "snr_all_sectors": None,
        "fpp": dict(value=s["fpp"]["FPP"], err=s["fpp"]["FPP_std"], n_runs=s["fpp"]["n_runs"],
                    without_pixel_exclusions=s["verdict"]["FPP_without_exclusions"]),
        "nfpp": dict(value=s["fpp"]["NFPP"], err=s["fpp"].get("NFPP_std")),
        "sectors": s["sectors"],
        "star": dict(tmag=float(star.tmag), teff=float(star.teff), rad=float(star.rad), mass=float(star.mass),
                     dist_pc=float(star.dist_pc), ra=float(star.ra), dec=float(star.dec)),
    }
    sec = pd.read_csv(FU / "sectors.csv").set_index("data")
    out["snr_all_sectors"] = float(sec.loc["all sectors", "snr_expected"])
    return out, med


def binned(x, y, edges, n_min=3):
    idx = np.digitize(x, edges) - 1
    xb, yb, eb = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= n_min:
            xb.append(x[m].mean()); yb.append(np.mean(y[m])); eb.append(np.std(y[m]) / np.sqrt(m.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


def fig_lightcurve(s, chain, med, par):
    P, T0 = med[1], med[0]
    T14, T23 = contact_times(chain)
    dur = float(np.median(T14))
    ms = pd.read_parquet(FU / "multisector.parquet")
    d = masked_flatten(ms, P, T0, dur)
    d = d[np.isfinite(d.flux_fit)]
    epoch = np.round((d.time - T0) / P)
    d["dt_h"] = (d.time - T0 - epoch * P) * 24
    d["epoch"] = epoch
    pm = params(med)
    pm.t0 = 0.0
    tm = np.linspace(-4, 4, 1000)
    model = batman.TransitModel(pm, tm / 24, supersample_factor=5, exp_time=200 / 86400).light_curve(pm)
    depth = par["depth_ppm"]["value"]

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 6, height_ratios=[1.3, 1, 1], hspace=0.45, wspace=0.35)
    ax = fig.add_subplot(gs[0, :4])
    w = d[np.abs(d.dt_h) < 4]
    ax.plot(w.dt_h, (w.flux_fit - 1) * 1e3, ".", color="0.75", ms=2, zorder=0)
    xb, yb, eb = binned(w.dt_h.to_numpy(), (w.flux_fit.to_numpy() - 1) * 1e3, np.arange(-4, 4.01, 1 / 6))
    ax.errorbar(xb, yb, eb, fmt="o", color="k", ms=4, label="10-min bins")
    ax.plot(tm, (model - 1) * 1e3, "C3-", lw=2, label="fit")
    ax.set_xlim(-4, 4); ax.set_ylim(-30, 15)
    ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("relative flux [ppt]")
    ax.set_title(f"TIC {TIC}: all six sectors (2019-2026) folded at P = {P:.6f} d "
                 f"(depth {depth:,.0f} ppm, T14 {par['duration_h']['value']:.2f} h)", fontsize=10, loc="left")
    ax.legend(fontsize=8, loc="lower left")

    # depth per sector: mean of the full-depth part against the out-of-transit level
    axd = fig.add_subplot(gs[0, 4:])
    t23 = float(np.median(T23)) * 24
    rows = []
    for sec, g in d.groupby("sector"):
        intr = g[np.abs(g.dt_h) < t23 / 2].flux_fit
        oot = g[(np.abs(g.dt_h) > dur * 24) & (np.abs(g.dt_h) < 6)].flux_fit
        if len(intr) >= 2:
            rows.append((sec, (oot.mean() - intr.mean()) * 1e6,
                         np.hypot(intr.std() / np.sqrt(len(intr)), oot.std() / np.sqrt(len(oot))) * 1e6))
    rows = np.array(rows)
    full_depth = (1 - model.min()) * 1e6
    axd.errorbar(range(len(rows)), rows[:, 1], rows[:, 2], fmt="o", color="k")
    axd.axhline(full_depth, color="C3", lw=1.5, label="fit (full depth)")
    axd.axhline(0, color="0.6", lw=0.8)
    axd.set_xticks(range(len(rows)), [f"s{int(r)}\n{SECTOR_YEAR[int(r)]}" for r in rows[:, 0]], fontsize=8)
    axd.set_ylabel("depth [ppm]"); axd.legend(fontsize=8)
    axd.set_title("depth in each sector (full-depth part)", fontsize=10, loc="left")

    for i, (sec, g) in enumerate(d.groupby("sector")):
        a = fig.add_subplot(gs[1, i])
        w = g[np.abs(g.dt_h) < 4]
        cad = int(g.cadence_s.iloc[0])
        a.plot(w.dt_h, (w.flux_fit - 1) * 1e3, ".", color="k" if cad >= 1800 else "0.75", ms=4 if cad >= 1800 else 2)
        bw = max(1 / 6, cad / 3600)
        xb, yb, eb = binned(w.dt_h.to_numpy(), (w.flux_fit.to_numpy() - 1) * 1e3, np.arange(-4, 4.01, bw))
        a.errorbar(xb, yb, eb, fmt="o", color="k", ms=3)
        a.plot(tm, (model - 1) * 1e3, "C3-", lw=1)
        a.set_ylim(-35, 20); a.set_xlim(-4, 4)
        n_tr = int(g[np.abs(g.dt_h) < dur * 12].epoch.nunique())
        a.set_title(f"s{sec} ({SECTOR_YEAR[int(sec)]}, {cad} s)\n{n_tr} transits observed", fontsize=8)
        a.tick_params(labelsize=7)
        if i:
            a.set_yticklabels([])
        else:
            a.set_ylabel("flux [ppt]")
        a.set_xlabel("hours", fontsize=8)

    w = d[np.abs(d.dt_h) < 4]
    ao = fig.add_subplot(gs[2, :3])
    for par_, c, lab in [(0, "C0", "even transits"), (1, "C1", "odd transits")]:
        m = (w.epoch % 2) == par_
        xb, yb, eb = binned(w.dt_h[m].to_numpy(), (w.flux_fit[m].to_numpy() - 1) * 1e3, np.arange(-4, 4.01, 1 / 3))
        ao.errorbar(xb + (0.03 if par_ else -0.03), yb, eb, fmt="o", color=c, ms=4, label=lab)
    ao.plot(tm, (model - 1) * 1e3, "k-", lw=1)
    ao.set_xlim(-4, 4); ao.set_ylim(-25, 10); ao.legend(fontsize=8)
    ao.set_xlabel("hours from mid-transit"); ao.set_ylabel("flux [ppt]")
    ao.set_title("odd and even transits (20-min bins)", fontsize=10, loc="left")

    a2 = fig.add_subplot(gs[2, 3:])
    phase = ((d.time - T0) / P) % 1
    dt2 = (phase - 0.5) * P * 24
    m = np.abs(dt2) < 4
    xb, yb, eb = binned(dt2[m].to_numpy(), (d.flux_fit[m].to_numpy() - 1) * 1e3, np.arange(-4, 4.01, 1 / 3))
    a2.errorbar(xb, yb, eb, fmt="o", color="k", ms=4)
    a2.axhline(0, color="0.6", lw=0.8)
    a2.axvspan(-dur * 12, dur * 12, color="C3", alpha=0.1, label="where a secondary eclipse would be")
    a2.set_xlim(-4, 4); a2.set_ylim(-25, 10); a2.legend(fontsize=8)
    a2.set_xlabel("hours from phase 0.5"); a2.set_title("phase 0.5 (20-min bins): no secondary eclipse", fontsize=10,
                                                          loc="left")
    return fig


def fig_pixels():
    import lightkurve as lk
    from src import pixelfit as pf

    fit = json.loads((FU / "summary.json").read_text())["fit"]
    chain = np.load(FU / "fit_samples.npz")["chain"]
    _, T23 = contact_times(chain)
    P, T0, dur_in = fit["P"]["median"], fit["T0"]["median"], float(np.median(T23))
    src = pd.read_csv(FU / "gaia_dr3_240arcsec_G19.5.csv")
    cache = ROOT / "notebooks" / ".lightkurve" / "cache"
    sectors = [12, 39, 66, 93, 103, 104]
    res = {}
    for sec in sectors:
        tpf = lk.search_tesscut(f"TIC {TIC}", sector=sec)[0].download(cutout_size=15, download_dir=str(cache))
        prf = pf.load_prf(tpf.camera, tpf.ccd, tpf.column + 7, tpf.row + 7, sec, pf.cache_dir())
        res[sec] = pf.source_test(tpf, src, P, T0, dur_in, prf, GAIA_TARGET, epoch=2016.0,
                                  n_in_min=1 if sec == 12 else 3, tmag_max_cand=17.5)

    r = res[104]
    ids = list(r["sources"].ID)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw=dict(width_ratios=[1, 1, 1.2]))
    vmax = np.nanmax(np.abs(r["diff"]))
    for a, img, ttl, cmap, lim in [(ax[0], r["scene"]["model"], "sector 104: PRF scene model (out of transit)", "viridis",
                                    None),
                                   (ax[1], r["diff"], "sector 104: difference image (out − in transit)", "RdBu_r", vmax)]:
        im = a.imshow(img, origin="lower", cmap=cmap, vmin=None if lim is None else -lim, vmax=lim)
        plt.colorbar(im, ax=a, fraction=0.046)
        a.contour(r["mask"], levels=[0.5], colors="k", linewidths=0.6)
        for sid, mk, cl in [(GAIA_TARGET, "x", "k"), (5925373272850782208, "+", "m"), (5925373268528684288, "^", "orange")]:
            if sid in ids:
                i = ids.index(sid)
                a.plot(r["sources"].x[i], r["sources"].y[i], mk, color=cl, ms=9, mew=2)
        a.set_title(ttl, fontsize=9)
        a.set_xlabel("column [px, 21\"]"); a.set_ylabel("row [px]")
    ax[1].plot(r["free"]["x"], r["free"]["y"], "o", mfc="none", mec="lime", ms=12, mew=2)
    ax[0].text(0.02, 0.02, "x target   + co-moving companion   ^ red giant", transform=ax[0].transAxes, fontsize=8,
               color="w")
    ax[1].text(0.02, 0.02, "o best-fit free position of the lost flux", transform=ax[1].transAxes, fontsize=8)

    # sector 12 (30-min cadence, one point per transit) is left out of the sums
    tab = pd.concat([rr["table"].assign(sector=sec) for sec, rr in res.items() if sec != 12])
    dchi2 = tab.groupby("ID").agg(sep_px=("sep_px", "first"), depth_needed=("depth_needed", "median"),
                                  dchi2=("dchi2_vs_target", "sum"))
    dchi2.to_csv(OUT / "pixel_dchi2.csv")
    dc = pd.Series({name: dchi2.dchi2.get(sid, np.nan) for sid, name in NAMES.items()}).dropna()
    order = [n for n in NAMES.values() if n in dc.index and n != "target"]
    vals = dc[order]
    colors = ["C3" if v > 25 else "0.6" for v in vals]
    ax[2].barh(range(len(order)), vals, color=colors)
    ax[2].set_yticks(range(len(order)), order, fontsize=9)
    ax[2].axvline(25, color="k", ls=":", lw=1)
    ax[2].text(26, len(order) - 0.6, "Δχ² = 25 (5σ)", fontsize=8)
    ax[2].set_xlabel("Δχ² of 'this star loses the light' vs 'the target does'\n(sum over sectors 39–104)")
    ax[2].set_title("which star loses the light? (red: Δχ² > 25)", fontsize=9, loc="left")
    for i, v in enumerate(vals):
        ax[2].text(max(v, 0) + 2, i, f"{v:.0f}", va="center", fontsize=8)
    plt.tight_layout()
    return fig


def fig_neighbours(par):
    g = pd.read_csv(FU / "gaia_dr3_60arcsec.csv")
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    tgt = g[g.Source == GAIA_TARGET].iloc[0]
    c0 = SkyCoord(tgt.RA_ICRS * u.deg, tgt.DE_ICRS * u.deg)
    g["sep"] = c0.separation(SkyCoord(g.RA_ICRS.values * u.deg, g.DE_ICRS.values * u.deg)).arcsec
    g = g[g.Source != GAIA_TARGET].copy()
    depth = par["depth_ppm"]["value"] / 1e6
    g["needed"] = depth * 10 ** (0.4 * (g.Gmag - tgt.Gmag))  # eclipse depth the neighbour would need
    excl = json.loads((FU / "pixel_excluded.json").read_text())
    shape_lim = excl["shape_limit_dm_99.7"]
    dchi2 = pd.read_csv(OUT / "pixel_dchi2.csv").set_index("ID").dchi2  # written by fig_pixels
    g["dchi2"] = g.Source.map(dchi2)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    too_faint = g.needed > 1
    ax.scatter(g.sep[too_faint], g.needed[too_faint] * 100, s=10, color="0.75",
               label="too faint: would need > 100 % eclipse")
    pix = (~too_faint) & (g.dchi2 > 25)
    ax.scatter(g.sep[pix], g.needed[pix] * 100, s=50, color="C3", marker="D",
               label="PRF pixel fit: Δχ² > 25 against the target (sectors 39–104)")
    shape = (~too_faint) & ~pix & (g.Gmag - tgt.Gmag > shape_lim)
    ax.scatter(g.sep[shape], g.needed[shape] * 100, s=40, color="C0", marker="s",
               label=f"transit shape only: host ≤ {shape_lim:.1f} mag fainter (99.7 %)")
    rest = (~too_faint) & ~shape & ~pix
    if rest.any():
        ax.scatter(g.sep[rest], g.needed[rest] * 100, s=60, color="k", marker="o", label="not excluded")
    for sid, per in KNOWN_EB.items():
        m = g.Source == sid
        if m.any():
            ax.scatter(g.sep[m], g.needed[m] * 100, s=200, facecolor="none", edgecolor="C2", lw=2)
            ax.annotate(f"Gaia EB, P = {per:.3f} d", (g.sep[m].iloc[0], g.needed[m].iloc[0] * 100),
                        xytext=(8, 8), textcoords="offset points", fontsize=8, color="C2")
    for sid, name in NAMES.items():
        m = g.Source == sid
        if m.any() and sid != GAIA_TARGET:
            ax.annotate(name, (g.sep[m].iloc[0], g.needed[m].iloc[0] * 100), xytext=(6, -12),
                        textcoords="offset points", fontsize=8)
    ax.axhline(100, color="k", lw=1)
    ax.axvline(21, color="0.5", ls=":", lw=1)
    ax.text(21.5, 0.25, "1 TESS pixel", fontsize=8, color="0.4")
    ax.set_yscale("log"); ax.set_ylim(0.1, 1e4); ax.set_xlim(0, 66)
    ax.set_xlabel("separation from the target [arcsec]")
    ax.set_ylabel("eclipse depth the star would need [%]")
    ax.set_title(f"Gaia DR3 sources within 60\" of TIC {TIC}: which could produce a {depth * 1e6:,.0f}-ppm dip?",
                 fontsize=10, loc="left")
    ax.legend(fontsize=8, loc="upper left")
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--skip-pixels", action="store_true", help="skip the pixel figure (needs TESScut downloads)")
    a = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    s, star, chain = load()
    par, med = ctoi_parameters(s, star, chain)
    (OUT / "ctoi_parameters.json").write_text(json.dumps(par, indent=2))
    figs = [fig_lightcurve(s, chain, med, par)]
    if not a.skip_pixels or not (OUT / "pixel_dchi2.csv").exists():
        figs.append(fig_pixels())
    figs.append(fig_neighbours(par))
    names = ["fig1_lightcurve", "fig2_pixels", "fig3_neighbours"] if len(figs) == 3 else ["fig1_lightcurve",
                                                                                         "fig3_neighbours"]
    with PdfPages(OUT / "summary.pdf") as pdf:
        for f, n in zip(figs, names):
            f.savefig(OUT / f"{n}.png", dpi=130, bbox_inches="tight")
            pdf.savefig(f, bbox_inches="tight")
    print(json.dumps({k: v for k, v in par.items() if k != "star"}, indent=1))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    os.chdir(ROOT)
    main()
