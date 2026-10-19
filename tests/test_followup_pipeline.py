"""Automated follow-up stage: verdict rules, survivor selection, resume and failure handling, CLI wiring."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from src import cli
from src import followup as fu

P, DEPTH, SECTOR = 2.0, 3000.0, 104


def _sectors(repeat: bool = True, alias: bool = False):
    """A sector table: discovery + two informative sectors that do (or do not) return P."""
    other_p = P * 1.001 if repeat else 1.37
    rows = [("s104", P, 9.0, DEPTH, 9.0), ("s101", other_p, 8.0, 2900, 9.0), ("s102", other_p, 7.0, 3100, 9.0),
            ("all sectors", P, 12.0, DEPTH, np.nan)]
    if alias:
        rows += [("s101@2P", 2 * P, 10.0, 7000, np.nan), ("s102@2P", 2 * P, 9.0, 7200, np.nan)]
    return [dict(zip(["data", "period", "snr", "depth_ppm", "snr_expected"], r)) for r in rows]


def _di(off_px=0.05, err=0.1, cat=None, strongest=12.0, at_target=10.0, dist=0.3, off_target=None):
    cat = off_px if cat is None else cat
    sig = off_px / err
    if off_target is None:
        off_target = sig > 3 and min(off_px, cat) > 0.5
    return {"offset_px": off_px, "offset_err_px": err, "offset_sigma": sig, "offset_from_catalog_px": cat,
            "strongest_loss_snr": strongest, "loss_snr_at_target": at_target, "strongest_loss_dist_px": dist,
            "off_target": off_target}


def _summary(sec=None, di=None, rho=None, fpp=None):
    s = {"tic_id": 1, "P": P, "T0": 4180.0, "depth_ppm": DEPTH, "discovery_sector": SECTOR,
         "sector_search": sec if sec is not None else _sectors(), "difference_image": di if di is not None else _di()}
    if rho is not None:
        s["fit"] = {"rho_ratio_median": rho, "rp_earth": {"median": 2.0, "plus": 0.1, "minus": 0.1}}
    if fpp is not None:
        s["fpp"] = {"FPP": fpp[0], "NFPP": fpp[1], "fpp_background": True}
    return s


@pytest.mark.parametrize("kwargs, category", [
    ({"rho": 0.9, "fpp": (0.0148, 0.0)}, "validated"),        # LHS 3844 b: 0.0148 must pass "< 0.015"
    ({"rho": 0.9, "fpp": (0.2, 0.01)}, "candidate"),
    ({"rho": 0.84, "fpp": (0.64, 0.64)}, "ambiguous"),        # TIC 214311908
    ({"rho": 0.9}, "incomplete"),                             # no FPP yet
    ({"rho": 0.1, "fpp": (0.45, 0.44)}, "blend"),             # TIC 64336445: density says blend despite FPP < 0.5
    ({"sec": _sectors(repeat=False)}, "not_repeating"),
    ({"sec": _sectors(repeat=False, alias=False)[:1] + _sectors()[3:]}, "untestable"),  # no other sector
    ({"di": _di(off_px=2.2, err=0.02)}, "blend"),             # off target wins even though it repeats
    ({"di": _di(off_px=2.2, err=0.02), "sec": _sectors(repeat=False)}, "blend"),  # ... and when it does not
])
def test_verdict_rules(kwargs, category):
    assert fu.verdict(_summary(**kwargs))["category"] == category


def test_eb_at_half_period():
    s = _summary(sec=_sectors(repeat=True, alias=True))
    s["sector_search"][1]["depth_ppm"] = s["sector_search"][2]["depth_ppm"] = 3000
    v = fu.verdict({**s, "difference_image": _di(off_px=2.3, err=0.01)})
    assert v["category"] == "eb" and "neighbour" in v["verdict"]


def test_pixel_verdicts():
    assert fu.pixel_verdict({})["pixels"] == "missing"
    assert fu.pixel_verdict({"error": "no cutout"})["pixels"] == "error"
    # no significant loss anywhere: the pixels cannot say (TIC 64336445: centroid 6.5 sigma off, but S/N < 5)
    assert fu.pixel_verdict(_di(off_px=2.4, err=0.37, cat=0.76, strongest=3.9, at_target=2.5))["pixels"] == "inconclusive"
    # off target only through the strongest pixel, centroid on target: weak (TIC 372868160)
    weak = _di(off_px=0.4, err=1.5, cat=2.1, strongest=12.0, at_target=2.4, dist=9.8, off_target=True)
    assert fu.pixel_verdict(weak)["pixels"] == "off_weak"
    assert fu.verdict(_summary(sec=_sectors(repeat=False), di=weak))["category"] == "not_repeating"
    assert fu.verdict(_summary(di=weak, rho=1.0, fpp=(0.1, 0.0)))["category"] == "blend"
    # crowded field: far from the light centroid, close to the catalog position -> on target (TIC 214311908)
    assert fu.pixel_verdict(_di(off_px=0.76, err=0.15, cat=0.30, dist=0.7))["pixels"] == "on"


def test_survivors():
    known = pd.DataFrame({"tic_id": [1, 2, 3, 4, 5],
                          "passes_all": [True, True, True, False, True],
                          "known_class": ["new", "new_in_known_system", "known_planet", "new", "neighbour"]})
    assert fu.survivors(known)["tic_id"].tolist() == [1, 2]


@pytest.fixture
def fake_data(tmp_path, monkeypatch):
    """A known table with two survivors and a FOLLOWUP_DIR in tmp_path."""
    (tmp_path / "search").mkdir()
    known = pd.DataFrame({"tic_id": [11, 22, 33], "passes_all": [True, True, True],
                          "known_class": ["new", "new", "known_planet"], "period": [P, 3.0, 1.0],
                          "T0": [4180.0, 4181.0, 4182.0], "depth_ppm": [DEPTH, DEPTH, DEPTH]})
    known.to_parquet(tmp_path / "search" / f"sector_{SECTOR}_known.parquet", index=False)
    monkeypatch.setattr(fu, "DATA_DIR", tmp_path)
    monkeypatch.setattr(fu, "FOLLOWUP_DIR", tmp_path / "followup")
    return tmp_path


def test_followup_sector_survives_failures_and_resumes(fake_data, monkeypatch):
    calls = []

    def fake_one(tic_id, sector, auto=False, **_):
        calls.append(tic_id)
        if tic_id == 22:
            raise RuntimeError("MAST down")
        s = {**_summary(rho=0.9, fpp=(0.2, 0.0)), "tic_id": tic_id, "T0": 4180.0,
             "followup_version": fu.FOLLOWUP_VERSION, "sectors": [101, 102, 104]}
        s["verdict"] = fu.verdict(s)
        d = fu.FOLLOWUP_DIR / str(tic_id)
        d.mkdir(parents=True, exist_ok=True)
        (d / "summary.json").write_text(json.dumps(s, default=float))
        return s

    monkeypatch.setattr(fu, "followup_one", fake_one)
    tab = fu.followup_sector(SECTOR, n_jobs=1).set_index("tic_id")
    assert sorted(calls) == [11, 22]                       # survivors only, not the known planet
    assert tab.loc[11, "category"] == "candidate"
    assert tab.loc[22, "category"] == "error" and "MAST down" in tab.loc[22, "error"]
    assert (fake_data / "followup" / f"sector_{SECTOR}_followup.parquet").exists()
    calls.clear()
    fu.followup_sector(SECTOR, n_jobs=1)                   # resume: 11 is up to date, the failed 22 is retried
    assert calls == [22]
    calls.clear()
    fu.followup_sector(SECTOR, n_jobs=1, force=True)
    assert sorted(calls) == [11, 22]


def test_up_to_date_needs_same_ephemeris_and_version(fake_data):
    d = fu.FOLLOWUP_DIR / "11"
    d.mkdir(parents=True)
    row = pd.Series({"period": P, "T0": 4180.0})
    s = {"P": P, "T0": 4180.0, "followup_version": fu.FOLLOWUP_VERSION, "verdict": {}}
    (d / "summary.json").write_text(json.dumps(s))
    assert fu.up_to_date(11, row)
    (d / "summary.json").write_text(json.dumps({**s, "P": P * 1.01}))
    assert not fu.up_to_date(11, row)
    (d / "summary.json").write_text(json.dumps({**s, "followup_version": fu.FOLLOWUP_VERSION - 1}))
    assert not fu.up_to_date(11, row)


def test_cli_runs_followup_before_report():
    assert cli.STEPS.index("known") < cli.STEPS.index("followup") < cli.STEPS.index("report")
    assert cli.selected_steps("vet", "report") == ["vet", "known", "followup", "report"]
    assert cli.selected_steps("vet", "report", skip_followup=True) == ["vet", "known", "report"]
    args = cli.step_args("followup", 104, 64, force=False)
    assert args[:2] == ["--sector", "104"]
    assert int(args[args.index("--n-jobs") + 1]) <= max(1, (__import__("os").cpu_count() or 2) // 2)


def test_network_errors_are_retried(fake_data, monkeypatch):
    class RemoteServiceError(Exception):
        pass

    n = {"calls": 0}

    def flaky(tic_id, sector, auto=False, **_):
        n["calls"] += 1
        if n["calls"] < 2:
            raise RemoteServiceError("Request failed.")
        return {"tic_id": tic_id}

    monkeypatch.setattr(fu, "followup_one", flaky)
    assert fu._safe_one(11, SECTOR, True, wait_s=0) == {"tic_id": 11} and n["calls"] == 2

    def bug(tic_id, sector, auto=False, **_):  # not a network error: recorded at once, no retry
        n["calls"] += 1
        raise KeyError("rad")

    n["calls"] = 0
    monkeypatch.setattr(fu, "followup_one", bug)
    assert "KeyError" in fu._safe_one(11, SECTOR, True, wait_s=0)["error"] and n["calls"] == 1


def test_pixel_status_from_null_images_wins():
    """Summaries with difference_image.status (pseudo-transit noise) use it; fixed thresholds are the fallback."""
    # TIC 372868160: old fields say off target (noise pixel), the null-calibrated status says inconclusive
    old = _di(off_px=0.4, err=1.5, cat=2.1, strongest=12.0, at_target=2.4, dist=9.8, off_target=True)
    assert fu.pixel_verdict(old)["pixels"] == "off_weak"
    assert fu.pixel_verdict({**old, "off_target": False, "status": "inconclusive"})["pixels"] == "inconclusive"
    # TIC 229804549: off target only through a significant peak 5.9 px away (centroid 2.7 sigma): conclusive
    peak = {**_di(off_px=2.52, err=0.93, strongest=19.2, at_target=7.2, dist=5.9, off_target=True),
            "status": "off_target"}
    assert fu.pixel_verdict(peak)["pixels"] == "off"
    assert fu.verdict(_summary(sec=_sectors(repeat=False), di=peak))["category"] == "blend"
    assert fu.pixel_verdict({**_di(), "status": "on_target"})["pixels"] == "on"


def test_fit_and_fpp_reuse():
    """An earlier fit / FPP for the same ephemeris is kept; an FPP without background scenarios is redone."""
    old = {"P": P, "T0": 4180.0, "fit": {"rho_ratio_median": 0.9}, "fpp": {"FPP": 0.6, "NFPP": 0.6}}
    assert fu._reusable(old, "fit", P, 4180.0)
    assert not fu._reusable(old, "fpp", P, 4180.0)                       # pre-TRILEGAL FPP: optimistic
    assert fu._reusable({**old, "fpp": {"FPP": 0.85, "fpp_background": True}}, "fpp", P, 4180.0)
    assert not fu._reusable(old, "fit", P * 1.01, 4180.0)                # the signal changed
    assert not fu._reusable({**old, "fit": {"error": "x"}}, "fit", P, 4180.0)


def test_refresh_verdicts(fake_data):
    d = fu.FOLLOWUP_DIR / "11"
    d.mkdir(parents=True)
    s = {**_summary(rho=0.9, fpp=(0.2, 0.0)), "tic_id": 11, "verdict": {"category": "stale"}}
    (d / "summary.json").write_text(json.dumps(s, default=float))
    tab = fu.refresh_verdicts(SECTOR).set_index("tic_id")
    assert tab.loc[11, "category"] == "candidate" and tab.loc[22, "category"] == "error"
    assert json.loads((d / "summary.json").read_text())["verdict"]["category"] == "candidate"


def test_fpp_without_background_is_never_a_candidate():
    s = _summary(rho=0.9, fpp=(0.0148, 0.0))
    s["fpp"]["fpp_background"] = False
    v = fu.verdict(s)
    assert v["category"] == "incomplete" and "optimistic" in v["verdict"] and v["fpp_background"] is False
    s["fpp"].update(fpp_background=True, fpp_excluded=["1", "2"], without_exclusions={"FPP": 0.45}, FPP_std=0.05)
    v = fu.verdict(s)
    assert v["category"] == "validated" and v["n_fpp_excluded"] == 2 and v["FPP_without_exclusions"] == 0.45
