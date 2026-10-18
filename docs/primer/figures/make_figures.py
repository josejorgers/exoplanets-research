"""Figures for the primer articles (docs/primer/). Synthetic data only, so they need no download.

    uv run python docs/primer/figures/make_figures.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

OUT = Path(__file__).parent
RNG = np.random.default_rng(42)

# reference palette (dataviz skill), light mode
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
POINTS = "#b5b4ae"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
STAR = "#f2c14e"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#e6e5e0",
        "grid.linewidth": 0.8,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "font.size": 10.5,
        "lines.linewidth": 2,
        "legend.frameon": False,
    }
)


def transit(t, t0, depth, dur, ingress_frac=0.15):
    """Flux of a trapezoid transit (1 = no dip), with a slightly rounded bottom like a real star.

    dur is the total duration; ingress_frac the fraction of it spent on each of ingress/egress.
    """
    x = np.abs(t - t0)
    half = dur / 2
    ing = ingress_frac * dur
    f = np.ones_like(t)
    full = x <= half - ing
    # rounded bottom: deepest at the centre (stars are darker at the edge: "limb darkening")
    f[full] = 1 - depth * (1.0 - 0.12 * (x[full] / (half - ing)) ** 2)
    edge = (x > half - ing) & (x < half)
    f[edge] = 1 - depth * 0.88 * (half - x[edge]) / ing
    return f


def periodic(t, period, t0, depth, dur, ingress_frac=0.15):
    ph = (t - t0 + period / 2) % period - period / 2
    return transit(ph, 0.0, depth, dur, ingress_frac)


def save(fig, name):
    fig.savefig(OUT / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ------------------------------------------------------------------------------------------------
# 1. one transit, with the geometry above it
# ------------------------------------------------------------------------------------------------
def fig_single_transit():
    fig = plt.figure(figsize=(9, 5.2))
    top = fig.add_axes([0.08, 0.68, 0.88, 0.28])
    ax = fig.add_axes([0.08, 0.12, 0.88, 0.5])

    t = np.linspace(-3, 3, 2000)  # hours
    f = transit(t, 0, 0.01, 2.0, 0.18)
    ax.plot(t, (f - 1) * 100, color=BLUE)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-1.25, 0.25)
    ax.set_xlabel("time (hours)")
    ax.set_ylabel("change in brightness (%)")
    ax.annotate(
        "dip depth: 1 %",
        xy=(0, -1.0),
        xytext=(1.05, -0.75),
        color=INK,
        arrowprops=dict(arrowstyle="-", color=INK2, lw=1),
    )
    ax.annotate("", xy=(-1, 0.12), xytext=(1, 0.12), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1))
    ax.text(0, 0.16, "duration: 2 hours", ha="center", va="bottom", color=INK)

    # cartoon: the star as seen from Earth at five moments (not to scale in time)
    top.set_xlim(-3.3, 3.3)
    top.set_ylim(-1, 1)
    top.set_aspect("equal")
    top.axis("off")
    moments = [(-2.6, -0.62, "before"), (-1.3, -0.42, "entering"), (0, 0, "in front"), (1.3, 0.42, "leaving"),
               (2.6, 0.62, "after")]
    for x, dx, label in moments:
        top.add_patch(Circle((x, 0.1), 0.42, color=STAR, zorder=1))
        top.add_patch(Circle((x + dx, 0.1), 0.1, color=INK, zorder=2))
        top.text(x, -0.55, label, ha="center", va="top", color=INK2, fontsize=9.5)
    save(fig, "single_transit.png")


# ------------------------------------------------------------------------------------------------
# 2. a month of data: transits hidden in noise, then folded
# ------------------------------------------------------------------------------------------------
MONTH = dict(period=3.7, t0=1.3, depth=0.006, dur=2.0 / 24, noise=0.003)


def month_times():
    t = np.arange(0, 27, 10 / 60 / 24)  # one point every 10 minutes for 27 days
    return t[(t < 13.2) | (t > 14.2)]  # the gap when the spacecraft sends its data home


def month_data(planet=True, seed=7):
    t = month_times()
    rng = np.random.default_rng(seed)
    f = 1 + rng.normal(0, MONTH["noise"], t.size)
    if planet:
        f *= periodic(t, MONTH["period"], MONTH["t0"], MONTH["depth"], MONTH["dur"])
    return t, f


def fig_month_and_fold():
    t, f = month_data()
    period, t0 = MONTH["period"], MONTH["t0"]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6.6), gridspec_kw=dict(hspace=0.55))
    a1.plot(t, (f - 1) * 100, ".", ms=2, color=POINTS)
    times = t0 + period * np.arange(8)
    for tt in times[times < 27]:
        a1.annotate("", xy=(tt, -1.45), xytext=(tt, -2.0), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.2))
    a1.set_xlim(0, 27)
    a1.set_ylim(-2.05, 1.5)
    a1.set_xlabel("time (days)")
    a1.set_ylabel("brightness change (%)")
    a1.set_title("A. 27 days of measurements: the dips (arrows) are hidden in the noise")

    ph = ((t - t0 + period / 2) % period - period / 2) * 24
    a2.plot(ph, (f - 1) * 100, ".", ms=2, color=POINTS, label="every measurement, stacked")
    bins = np.arange(-6, 6.01, 0.25)
    idx = np.digitize(ph, bins)
    centres = 0.5 * (bins[1:] + bins[:-1])
    means = np.array([np.mean(f[idx == i]) if np.any(idx == i) else np.nan for i in range(1, len(bins))])
    a2.plot(centres, (means - 1) * 100, "o-", ms=5, lw=1.5, color=BLUE, label="average in 15-minute bins")
    a2.set_xlim(-6, 6)
    a2.set_ylim(-1.5, 1.5)
    a2.set_xlabel("hours from the predicted middle of the dip")
    a2.set_ylabel("brightness change (%)")
    a2.set_title(f"B. Same data cut every {period} days and stacked: the dip appears")
    a2.legend(loc="lower right")
    save(fig, "month_and_fold.png")


# ------------------------------------------------------------------------------------------------
# 3. depth -> size
# ------------------------------------------------------------------------------------------------
def fig_depth_vs_size():
    r_earth_km, r_sun_km = 6371.0, 695700.0
    r = np.logspace(np.log10(0.5), np.log10(30), 300)  # Earth radii
    fig, ax = plt.subplots(figsize=(9, 4.8))
    for rstar, color, name in [(1.0, ORANGE, "a Sun-like star"), (0.4, BLUE, "a small red dwarf (40 % of the Sun)")]:
        depth = (r * r_earth_km / (rstar * r_sun_km)) ** 2 * 100
        ax.plot(r, depth, color=color, label=f"passing in front of {name}")
    for name, rr in [("Earth", 1.0), ("Neptune", 3.9), ("Jupiter", 11.2), ("a small star", 21.8)]:
        d = (rr * r_earth_km / (0.4 * r_sun_km)) ** 2 * 100
        ax.plot(rr, d, "o", ms=8, color=BLUE, mec=SURFACE, mew=2, zorder=3)
        ax.annotate(f"{name}: {d:.2g} %", (rr, d), xytext=(-8, 8), textcoords="offset points", ha="right", color=INK)
    d_earth_sun = (r_earth_km / r_sun_km) ** 2 * 100
    ax.plot(1.0, d_earth_sun, "o", ms=8, color=ORANGE, mec=SURFACE, mew=2, zorder=3)
    ax.annotate(
        f"Earth: {d_earth_sun:.3f} %", (1.0, d_earth_sun), xytext=(8, -14), textcoords="offset points", color=INK
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.5, 1, 2, 5, 10, 20])
    ax.set_xticklabels(["0.5", "1", "2", "5", "10", "20"])
    ax.set_yticks([0.001, 0.01, 0.1, 1, 10, 100])
    ax.set_yticklabels(["0.001 %", "0.01 %", "0.1 %", "1 %", "10 %", "100 %"])
    ax.set_xlabel("size of the passing object (Earth = 1)")
    ax.set_ylabel("depth of the dip")
    ax.legend(loc="upper left")
    save(fig, "depth_vs_size.png")


# ------------------------------------------------------------------------------------------------
# 4. impostors
# ------------------------------------------------------------------------------------------------
def fig_odd_even():
    """Eclipsing binary of two unequal stars: alternating dips; the search finds half the true period."""
    p_true = 4.0
    t = np.arange(0, 16, 5 / 60 / 24)
    f = periodic(t, p_true, 1.0, 0.06, 0.12, 0.45) * periodic(t, p_true, 3.0, 0.035, 0.12, 0.45)
    f = f + RNG.normal(0, 0.002, t.size)
    fig, ax = plt.subplots(figsize=(9, 3.4))
    ax.plot(t, (f - 1) * 100, ".", ms=2, color=POINTS)
    for i, tt in enumerate(np.arange(1.0, 16, 2.0)):
        ax.text(tt, 1.0, "odd" if i % 2 == 0 else "even", ha="center", color=ORANGE if i % 2 == 0 else BLUE)
    ax.axhline(-6, color=ORANGE, lw=1, ls="--")
    ax.axhline(-3.5, color=BLUE, lw=1, ls="--")
    ax.set_xlim(0, 16)
    ax.set_ylim(-7, 2)
    ax.set_xlabel("time (days)")
    ax.set_ylabel("brightness change (%)")
    ax.set_title("Dips every 2 days, but every other one is deeper: two stars eclipsing each other")
    save(fig, "eb_odd_even.png")


def fig_secondary():
    """Folded binary: big dip at phase 0, small one at phase 0.5 (the companion hidden behind the star)."""
    period = 2.0
    t = np.arange(0, 27, 10 / 60 / 24)
    f = periodic(t, period, 0.5, 0.012, 0.09, 0.4) * periodic(t, period, 1.5, 0.002, 0.09, 0.4)
    f = f + RNG.normal(0, 0.002, t.size)
    phase = ((t - 0.5) / period + 0.25) % 1 - 0.25
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.plot(phase, (f - 1) * 100, ".", ms=2, color=POINTS)
    bins = np.linspace(-0.25, 0.75, 161)
    idx = np.digitize(phase, bins)
    c = 0.5 * (bins[1:] + bins[:-1])
    m = np.array([np.mean(f[idx == i]) for i in range(1, len(bins))])
    ax.plot(c, (m - 1) * 100, color=BLUE, lw=1.8)
    ax.annotate("main dip", (0, -1.2), xytext=(0.08, -1.1), color=INK, arrowprops=dict(arrowstyle="-", color=INK2))
    ax.annotate(
        "second, smaller dip exactly half-way",
        (0.5, -0.22),
        xytext=(0.47, -1.32),
        ha="right",
        color=INK,
        arrowprops=dict(arrowstyle="-", color=INK2),
    )
    ax.set_xlim(-0.25, 0.75)
    ax.set_ylim(-1.5, 0.8)
    ax.set_xticks([-0.25, 0, 0.25, 0.5, 0.75])
    ax.set_xticklabels(["-¼", "0", "¼", "½", "¾"])
    ax.set_xlabel("position in the cycle (0 = main dip, ½ = half-way to the next one)")
    ax.set_ylabel("brightness change (%)")
    save(fig, "eb_secondary.png")


def fig_shapes():
    t = np.linspace(-2.5, 2.5, 1000)
    fig, ax = plt.subplots(figsize=(9, 3.8))
    ax.plot(t, (transit(t, 0, 0.01, 2.0, 0.12) - 1) * 100, color=BLUE, label="planet: steep sides, flat-ish bottom (U)")
    ax.plot(t, (transit(t, 0, 0.01, 2.0, 0.5) - 1) * 100, color=ORANGE, label="star grazing a star: no flat part (V)")
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-1.2, 0.3)
    ax.set_xlabel("time (hours)")
    ax.set_ylabel("brightness change (%)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    save(fig, "shapes.png")


def fig_blend():
    """Pixel cartoon: a bright target and a faint neighbour in the same camera pixels."""
    n = 9
    y, x = np.mgrid[0:n, 0:n]

    def psf(x0, y0, amp, s=0.9):
        return amp * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * s**2))

    img = psf(3.6, 4.2, 1.0) + psf(5.6, 3.6, 0.1)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4.2), gridspec_kw=dict(width_ratios=[1, 1.25], wspace=0.45))
    a1.imshow(img ** 0.5, cmap="Greys", origin="lower")
    a1.add_patch(Rectangle((1.5, 2.5), 5, 4, fill=False, ec=ORANGE, lw=2))
    a1.annotate("our star", (3.6, 4.2), xytext=(0.2, 7.8), color=INK, arrowprops=dict(arrowstyle="-", color=INK2))
    a1.annotate(
        "faint neighbour", (5.6, 3.6), xytext=(4.0, 0.3), color=INK, arrowprops=dict(arrowstyle="-", color=INK2)
    )
    a1.text(1.5, 6.65, "light we add up", color=ORANGE, fontsize=9.5)
    a1.set_xticks(np.arange(-0.5, n, 1))
    a1.set_yticks(np.arange(-0.5, n, 1))
    a1.set_xticklabels([])
    a1.set_yticklabels([])
    a1.tick_params(length=0)
    a1.grid(color="#d0cfc8", lw=0.6)
    a1.set_title("Camera view (1 square = 1 pixel)", fontsize=10.5)

    labels = ["our star", "neighbour", "total we measure"]
    before = np.array([100, 10, 110])
    during = np.array([100, 5, 105])  # the neighbour loses half its light in an eclipse
    xs = np.arange(3)
    a2.bar(xs - 0.2, before, 0.38, color=BLUE, label="normally")
    a2.bar(xs + 0.2, during, 0.38, color=ORANGE, label="while the neighbour is eclipsed")
    a2.set_xticks(xs)
    a2.set_xticklabels(labels)
    a2.set_ylabel("light (arbitrary units)")
    a2.set_ylim(0, 125)
    a2.text(2, 116, "−4.5 %", ha="center", color=INK, fontweight="bold")
    a2.text(1, 15, "−50 %", ha="center", color=INK, fontweight="bold")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    a2.grid(axis="x", visible=False)
    a2.set_title("Neighbour loses 50 %: we see a 4.5 % dip", fontsize=10.5)
    save(fig, "blend.png")


def fig_star_and_instrument():
    t = np.arange(0, 27, 10 / 60 / 24)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6), gridspec_kw=dict(hspace=0.55))

    # spotted, rotating star with flares
    f = 1 + 0.012 * np.sin(2 * np.pi * t / 2.3) + 0.004 * np.sin(2 * np.pi * t / 1.15 + 1) + RNG.normal(0, 0.0015, t.size)
    for tf in [3.1, 9.7, 17.4, 22.8]:
        k = t >= tf
        f[k] += 0.03 * np.exp(-(t[k] - tf) / 0.02)
    a1.plot(t, (f - 1) * 100, ".", ms=2, color=POINTS)
    a1.set_xlim(0, 27)
    a1.set_ylim(-2.2, 3.5)
    a1.set_ylabel("brightness change (%)")
    a1.set_xlabel("time (days)")
    a1.set_title("A. The star itself: dark spots rotating in and out of view, plus flares")
    a1.annotate("flare", (9.72, 2.9), xytext=(11.0, 2.7), color=INK, arrowprops=dict(arrowstyle="-", color=INK2))

    # instrument: ramp at the start of each orbit, a step, a glitch
    g = np.ones_like(t)
    for start in [0.0, 14.2]:
        k = t >= start
        g[k] -= 0.02 * np.exp(-(t[k] - start) / 0.4)
    g[(t > 6.0)] += 0.0
    g[(t > 19.5) & (t < 21.0)] -= 0.006
    glitch = (t > 9.0) & (t < 9.05)
    g = g + RNG.normal(0, 0.0015, t.size)
    g[glitch] -= 0.012 + RNG.normal(0, 0.006, glitch.sum())
    gap = (t > 13.2) & (t < 14.2)
    g[gap] = np.nan
    a2.plot(t, (g - 1) * 100, ".", ms=2, color=POINTS)
    a2.set_xlim(0, 27)
    a2.set_ylim(-2.5, 1)
    a2.set_ylabel("brightness change (%)")
    a2.set_xlabel("time (days)")
    a2.set_title("B. The camera: ramps after each data download, a step, a short glitch")
    for xy, txt, xt in [((0.25, -1.4), "ramp", (1.5, -2.2)), ((14.45, -1.4), "ramp", (15.7, -2.2)),
                        ((9.02, -1.6), "glitch", (10.2, -2.2)), ((20.2, -0.6), "step", (21.6, -1.6))]:
        a2.annotate(txt, xy, xytext=xt, color=INK, arrowprops=dict(arrowstyle="-", color=INK2))
    a2.text(13.7, 0.6, "gap", ha="center", color=INK2)
    save(fig, "star_and_instrument.png")


# ------------------------------------------------------------------------------------------------
# geometry: which orbits produce a transit
# ------------------------------------------------------------------------------------------------
def fig_geometry():
    a = 6.0  # orbit radius in star radii (real close-in planets: 10-30; Earth-Sun: 215)
    cases = [(0.0, BLUE, "edge-on: crosses the middle"), (8.0, AQUA, "tilted 8°: crosses near the edge"),
             (20.0, ORANGE, "tilted 20°: passes above, no dip")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.4), gridspec_kw=dict(width_ratios=[1.7, 1], wspace=0.15))

    a1.set_aspect("equal")
    a1.axis("off")
    a1.axhspan(-1, 1, xmin=0.5, color="#efeee9", zorder=0)
    a1.set_xlim(-7, 9.5)
    a1.set_ylim(-3, 3)
    a1.axhline(1, color=INK2, lw=0.8, ls="--", xmin=0.43)
    a1.axhline(-1, color=INK2, lw=0.8, ls="--", xmin=0.43)
    a1.add_patch(Circle((0, 0), 1, color=STAR, zorder=2))
    for th, c, lab in cases:
        r = np.radians(th)
        a1.plot([-a * np.cos(r), a * np.cos(r)], [-a * np.sin(r), a * np.sin(r)], color=c, lw=1.5, zorder=1)
        a1.add_patch(Circle((a * np.cos(r), a * np.sin(r)), 0.22, color=c, zorder=3))
    a1.annotate("", xy=(9.3, -1.8), xytext=(6.9, -1.8), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    a1.text(6.7, -1.8, "to Earth", ha="right", va="center", color=INK)
    a1.text(7.2, 1.25, "the star's shadow band", color=INK2, fontsize=9.5)
    a1.text(-6.9, 2.75, "Side view (not to scale)", color=INK, fontweight="bold", fontsize=11)

    a2.set_aspect("equal")
    a2.axis("off")
    a2.set_xlim(-2.2, 2.2)
    a2.set_ylim(-1.4, 2.6)
    a2.add_patch(Circle((0, 0), 1, color=STAR, zorder=1))
    for th, c, lab in cases:
        h = a * np.sin(np.radians(th))
        a2.plot([-2.1, 2.1], [h, h], color=c, lw=1.5, ls=":", zorder=2)
        a2.add_patch(Circle((0.6, h), 0.12, color=c, zorder=3))
    a2.text(-2.1, 2.45, "What we see from Earth", color=INK, fontweight="bold", fontsize=11)
    handles = [plt.Line2D([], [], color=c, lw=2) for _, c, _ in cases]
    fig.legend(handles, [lab for _, _, lab in cases], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04))
    save(fig, "transit_geometry.png")


# ------------------------------------------------------------------------------------------------
# noise
# ------------------------------------------------------------------------------------------------
def gauss(x, s):
    return np.exp(-0.5 * (x / s) ** 2)


def fig_bell():
    x = np.linspace(-4, 4, 800)
    fig, ax = plt.subplots(figsize=(9, 3.6))
    y = gauss(x, 1) / np.sqrt(2 * np.pi)
    for k, c in [(3, "#e6eefa"), (2, "#bcd2f1"), (1, "#8db3e8")]:
        m = np.abs(x) <= k
        ax.fill_between(x[m], 0, y[m], color=c, lw=0)
    ax.plot(x, y, color=BLUE)
    ax.text(0, 0.17, "68 %", ha="center", color=INK, fontweight="bold")
    ax.text(1.5, 0.045, "95 %\n(both shades)", ha="center", color=INK, fontsize=9)
    ax.text(3.0, 0.045, "99.7 %", ha="center", color=INK, fontsize=9)
    ax.set_xticks(range(-4, 5))
    ax.set_xticklabels(["−4σ", "−3σ", "−2σ", "−σ", "0", "+σ", "+2σ", "+3σ", "+4σ"])
    ax.set_yticks([])
    ax.set_xlabel("how far a measurement falls from the true value")
    ax.set_ylabel("how often")
    ax.set_xlim(-4, 4)
    ax.set_ylim(0, 0.43)
    save(fig, "bell_curve.png")


def fig_averaging():
    sigma, depth = 0.3, 0.6
    rows = [(1, "one measurement"), (12, "average of 12 measurements (one 2-hour dip)"),
            (84, "average of 84 measurements (seven dips stacked)")]
    x = np.linspace(-1.1, 0.9, 1000)
    fig, axes = plt.subplots(3, 1, figsize=(9, 6.4), sharex=True, gridspec_kw=dict(hspace=0.6))
    for ax, (n, lab) in zip(axes, rows):
        s = sigma / np.sqrt(n)
        ax.fill_between(x, 0, gauss(x, s), color="#bcd2f1", lw=0)
        ax.plot(x, gauss(x, s), color=BLUE, lw=1.5)
        ax.axvline(-depth, color=ORANGE, lw=2)
        ax.set_yticks([])
        ax.set_ylim(0, 1.15)
        ax.set_title(f"{lab}: spread {s:.3f} %, dip = {depth / s:.1f} spreads", fontsize=10.5)
    axes[0].text(-depth - 0.02, 0.85, "the planet's dip (−0.6 %)", color=INK, ha="right")
    axes[-1].set_xlabel("average brightness change (%) if there were NO planet: noise only")
    save(fig, "averaging.png")


def fig_white_red():
    t = np.arange(0, 3, 10 / 60 / 24)
    rng = np.random.default_rng(3)
    white = rng.normal(0, 0.3, t.size)
    slow = np.convolve(rng.normal(0, 1, t.size + 60), np.ones(60) / 60, mode="same")[30 : 30 + t.size]
    red = rng.normal(0, 0.24, t.size) + slow / slow.std() * 0.18
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.4), sharex=True, gridspec_kw=dict(hspace=0.45))
    for ax, y, title in [(axes[0], white, "A. Independent ('white') noise: 2-hour averages stay close to zero"),
                         (axes[1], red, "B. Noise with slow wiggles ('red'): 2-hour averages wander and can fake a dip")]:
        ax.plot(t, y, ".", ms=2.5, color=POINTS)
        nb = t.size // 12
        tb = t[: nb * 12].reshape(nb, 12).mean(1)
        yb = y[: nb * 12].reshape(nb, 12).mean(1)
        ax.plot(tb, yb, "o-", ms=5, lw=1.5, color=BLUE, label="2-hour averages")
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_ylim(-1.2, 1.2)
        ax.set_ylabel("change (%)")
        ax.set_title(title, fontsize=10.5)
    axes[0].legend(loc="upper right")
    axes[1].set_xlabel("time (days)")
    save(fig, "white_vs_red.png")


# ------------------------------------------------------------------------------------------------
# folding: a toy example small enough to do by hand
# ------------------------------------------------------------------------------------------------
TOY_DIP_HOURS = [2, 8, 14, 20]


def toy_data():
    rng = np.random.default_rng(11)
    t = np.arange(24)
    f = 100 + np.round(rng.normal(0, 0.8, 24))
    f[TOY_DIP_HOURS] -= 3
    return t, f


def fig_toy_fold():
    t, f = toy_data()
    fig, axes = plt.subplots(3, 1, figsize=(9, 8.2), gridspec_kw=dict(hspace=0.65, height_ratios=[1, 1, 1]))
    a0 = axes[0]
    a0.plot(t, f, "-", color="#d8d7d0", lw=1, zorder=1)
    a0.plot(t, f, "o", ms=7, color=POINTS, zorder=2)
    a0.plot(TOY_DIP_HOURS, f[TOY_DIP_HOURS], "o", ms=8, color=ORANGE, zorder=3, label="the planet's dips (every 6 h)")
    for k in range(1, 4):
        a0.axvline(6 * k - 0.5, color=INK2, lw=0.8, ls="--")
    a0.set_xticks(range(0, 24, 2))
    a0.set_xlim(-0.6, 23.6)
    a0.set_ylim(94.5, 102.5)
    a0.set_xlabel("hour")
    a0.set_ylabel("brightness")
    a0.set_title("A. 24 hourly measurements; dashed lines cut them into slices of 6 hours", fontsize=10.5)
    a0.legend(loc="lower left")

    for ax, period, title in [(axes[1], 6, "B. Slices of 6 hours stacked: the 4 dips land in the same slot"),
                              (axes[2], 5, "C. Slices of 5 hours (wrong period): the dips land in different slots")]:
        slot = t % period
        dip = np.isin(t, TOY_DIP_HOURS)
        jitter = (t // period - 1.5) * 0.06
        ax.plot(slot[~dip] + jitter[~dip], f[~dip], "o", ms=7, color=POINTS)
        ax.plot(slot[dip] + jitter[dip], f[dip], "o", ms=8, color=ORANGE)
        means = [f[slot == k].mean() for k in range(period)]
        ax.plot(range(period), means, "s-", ms=8, color=BLUE, lw=1.5, label="average of each slot")
        for k, m in enumerate(means):
            ax.text(k + 0.15, m - 0.05, f"{m:.2f}", color=INK, fontsize=9.5, va="top")
        ax.set_xticks(range(period))
        ax.set_xticklabels([f"slot {k}" for k in range(period)])
        ax.set_xlim(-0.5, 5.6)
        ax.set_ylim(95.5, 102.5)
        ax.set_ylabel("brightness")
        ax.set_title(title, fontsize=10.5)
        ax.legend(loc="lower right")
    save(fig, "toy_fold.png")


# ------------------------------------------------------------------------------------------------
# the search: try every period, every start, a few durations; score = SNR of the box
# ------------------------------------------------------------------------------------------------
def box_search(t, f, periods, durs_h=(1, 2, 3), bin_min=20):
    """Best box SNR per trial period. Same idea as BLS/TLS, written for clarity, not speed."""
    sigma = 1.4826 * np.median(np.abs(f - np.median(f)))
    y = f - np.median(f)
    binw = bin_min / 1440
    ks = [round(d * 60 / bin_min) for d in durs_h]
    out = np.zeros(len(periods))
    for i, P in enumerate(periods):
        nb = int(np.ceil(P / binw))
        idx = np.minimum(((t % P) / binw).astype(int), nb - 1)
        S = np.bincount(idx, weights=y, minlength=nb)
        N = np.bincount(idx, minlength=nb).astype(float)
        kmax = max(ks)
        cs = np.concatenate([[0], np.cumsum(np.concatenate([S, S[:kmax]]))])
        cn = np.concatenate([[0], np.cumsum(np.concatenate([N, N[:kmax]]))])
        best = 0.0
        for k in ks:
            s = cs[k : k + nb] - cs[:nb]
            n = cn[k : k + nb] - cn[:nb]
            ok = n >= 6
            if ok.any():
                best = max(best, np.max(-s[ok] / (sigma * np.sqrt(n[ok]))))
        out[i] = best
    return out


def period_grid(pmin=0.5, pmax=13.0, baseline=27.0, dur=2 / 24):
    df = dur / (3 * baseline * pmax)
    return 1 / np.arange(1 / pmax, 1 / pmin, df)[::-1]


def _best_noise(seed):
    t, f = month_data(planet=False, seed=1000 + seed)
    return box_search(t, f, period_grid()).max()


def fig_search(n_noise=100):
    from joblib import Parallel, delayed

    periods = period_grid()
    t, f = month_data()
    p_planet = box_search(t, f, periods)
    t, g = month_data(planet=False, seed=1000)
    p_noise = box_search(t, g, periods)
    best_noise = np.array(Parallel(n_jobs=-1)(delayed(_best_noise)(s) for s in range(n_noise)))

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, sharey=True, gridspec_kw=dict(hspace=0.4))
    for ax, p, title in [(axes[0], p_planet, "A. Star with a planet (the data of the figure above)"),
                         (axes[1], p_noise, "B. Star without a planet: noise only")]:
        ax.plot(periods, p, color=BLUE, lw=0.8)
        j = np.argmax(p)
        ax.annotate(f"best: {periods[j]:.2f} days, SNR {p[j]:.1f}", (periods[j], p[j]), xytext=(-8, 2),
                    textcoords="offset points", color=INK, va="bottom", ha="right")
        ax.set_ylabel("dip strength (SNR)")
        ax.set_title(title, fontsize=10.5)
    for mult, lab in [(2, "2 × period"), (0.5, "½ × period")]:
        pp = MONTH["period"] * mult
        j = np.argmin(np.abs(periods - pp))
        axes[0].annotate(lab, (pp, p_planet[max(0, j - 30) : j + 30].max()), xytext=(8, 6),
                         textcoords="offset points", color=INK2, fontsize=9)
    axes[1].set_xlabel("trial period (days)")
    axes[1].set_xscale("log")
    axes[1].set_xticks([0.5, 1, 2, 3, 5, 10])
    axes[1].set_xticklabels(["0.5", "1", "2", "3", "5", "10"])
    axes[0].set_ylim(0, 22)
    save(fig, "periodogram.png")

    fig, ax = plt.subplots(figsize=(9, 3.4))
    ax.hist(best_noise, bins=np.arange(3.5, 7.01, 0.1), color=BLUE, edgecolor=SURFACE, lw=1.5)
    ax.axvline(7, color=ORANGE, lw=2)
    ax.text(6.95, ax.get_ylim()[1] * 0.85, "threshold 7", ha="right", color=INK)
    ax.set_xlabel(f"best SNR found by the search in a star with no planet ({n_noise} simulated stars)")
    ax.set_ylabel("number of stars")
    ax.grid(axis="x", visible=False)
    save(fig, "noise_best_snr.png")

    print(f"period grid: {len(periods)} trial periods")
    print(f"planet star: best P = {periods[np.argmax(p_planet)]:.3f} d, SNR = {p_planet.max():.1f}")
    print(f"noise star : best P = {periods[np.argmax(p_noise)]:.3f} d, SNR = {p_noise.max():.1f}")
    print("noise stars best SNR: min %.2f median %.2f max %.2f" % (best_noise.min(), np.median(best_noise), best_noise.max()))


# ------------------------------------------------------------------------------------------------
# eclipsing binary over one orbit
# ------------------------------------------------------------------------------------------------
EB = dict(light_a=1000, light_b=100, area_a=100, area_b=25)  # B: half the width, 4x dimmer per unit area


def fig_eb_orbit():
    la, lb, aa, ab = EB["light_a"], EB["light_b"], EB["area_a"], EB["area_b"]
    lost_primary = la * ab / aa
    fig = plt.figure(figsize=(9, 5.6))
    top = fig.add_axes([0.1, 0.66, 0.86, 0.3])
    ax = fig.add_axes([0.1, 0.1, 0.86, 0.5])
    top.set_xlim(-0.08, 1.08)
    top.set_ylim(-0.14, 0.112)
    top.set_aspect("equal")
    top.axis("off")
    ra, rb = 0.07, 0.035
    for ph, label in [(0, "B in front of A"), (0.25, "side by side"), (0.5, "B behind A"), (0.75, "side by side"),
                      (1.0, "B in front again")]:
        dx = 0.11 * np.sin(2 * np.pi * ph)
        front_b = np.cos(2 * np.pi * ph) > 0
        za, zb = (1, 2) if front_b else (2, 1)
        top.add_patch(Circle((ph, 0.03), ra, color=STAR, zorder=za))
        top.add_patch(Circle((ph + dx, 0.03), rb, color=ORANGE, zorder=zb))
        top.text(ph, -0.06, label, ha="center", va="top", color=INK2, fontsize=9)
    fig.text(0.1, 0.98, "A: big, hot star (yellow)     B: smaller, cooler star (orange)", color=INK, fontsize=10)

    ph = np.linspace(0, 1, 2000)
    total = np.full_like(ph, la + lb, dtype=float)
    total -= lost_primary * (1 - transit(ph, 0, 1, 0.08, 0.3)) + lost_primary * (1 - transit(ph, 1, 1, 0.08, 0.3))
    total -= lb * (1 - transit(ph, 0.5, 1, 0.08, 0.3))
    ax.plot(ph, total, color=BLUE)
    ax.axhline(la + lb, color=INK2, lw=0.8, ls="--")
    ax.text(0.25, la + lb + 15, f"A + B = {la + lb}: what we normally receive", ha="center", color=INK)
    ax.annotate(f"B hides part of A: −{lost_primary:.0f} ({lost_primary / (la + lb) * 100:.0f} %)",
                (0.03, la + lb - lost_primary), xytext=(0.09, 870), color=INK,
                arrowprops=dict(arrowstyle="-", color=INK2))
    ax.annotate(f"A hides all of B: −{lb} ({lb / (la + lb) * 100:.0f} %)", (0.5, la), xytext=(0.56, 930),
                color=INK, arrowprops=dict(arrowstyle="-", color=INK2))
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(820, 1150)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "¼", "½", "¾", "1 orbit"])
    ax.set_xlabel("time, as a fraction of one orbit")
    ax.set_ylabel("total light received")
    save(fig, "eb_orbit.png")


if __name__ == "__main__":
    fig_single_transit()
    fig_month_and_fold()
    fig_depth_vs_size()
    fig_odd_even()
    fig_secondary()
    fig_shapes()
    fig_blend()
    fig_star_and_instrument()
    fig_geometry()
    fig_bell()
    fig_averaging()
    fig_white_red()
    fig_toy_fold()
    fig_eb_orbit()
    fig_search()

    from math import erfc, log, sqrt

    for k in (1, 2, 3, 5, 7):
        print(f"P(noise below -{k} sigma) = {0.5 * erfc(k / sqrt(2)):.3g}")
    for m in (1e3, 1e5, 1e7, 1e9):
        print(f"typical max of {m:.0e} independent tries ~ {sqrt(2 * log(m)):.1f} sigma")
    t, f = toy_data()
    print("toy data:", f.astype(int).tolist())
    for period in (6, 5):
        print(f"toy fold {period}:", [round(f[t % period == k].mean(), 2) for k in range(period)])

    # numbers quoted in the article
    r_sun_au = 695700 / 1.496e8
    m, rs, p = 0.4, 0.4, 3.0  # solar masses, solar radii, days
    a_au = (m * (p / 365.25) ** 2) ** (1 / 3)
    dur_h = p / np.pi * (rs * r_sun_au / a_au) * 24
    print(f"M dwarf 0.4/0.4, P=3 d: a = {a_au:.4f} AU, central duration = {dur_h:.2f} h")
    print(f"Earth around the Sun, P=365 d: duration = {365.25 / np.pi * r_sun_au * 24:.1f} h")
    print(f"M dwarf P=3 d: a/R* = {a_au / (rs * r_sun_au):.1f}, transit probability = {rs * r_sun_au / a_au:.3f}")
    for p in (1, 3, 10, 30):
        a = (m * (p / 365.25) ** 2) ** (1 / 3)
        print(f"  P={p} d: a = {a:.4f} AU, probability R*/a = {rs * r_sun_au / a * 100:.1f} %")
    print(f"Earth-Sun probability = {r_sun_au * 100:.2f} %")
