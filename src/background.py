"""Background-star population for TRICERATOPS: a TRILEGAL Galaxy-model simulation per target.

TRICERATOPS' background scenarios (BTP, BEB and their x2P variants, plus DTP/DEB for diluted
light curves) need a list of the faint stars expected along the line of sight, simulated by the
TRILEGAL web service (Girardi et al. 2005; http://stev.oapd.inaf.it/cgi-bin/trilegal).

Why this module exists: triceratops' own query (`funcs.query_TRILEGAL`) always failed for us,
silently, and the background scenarios were dropped. The cause is not the server. It sends only
its leaf certificate and leaves out the ZeroSSL intermediate, so a strict TLS client (requests +
certifi) cannot build the chain and the query dies with "unable to get local issuer
certificate". Browsers fill in the missing certificate from the AIA URL; we do the same: download
the intermediate once, append it to certifi's bundle, and verify against that. Verification
stays on.

The simulated population is cached as `data/followup/<tic>/trilegal.csv`, in the format
`triceratops.funcs.save_trilegal` writes (so it can be passed as `trilegal_fname`).
"""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

TRILEGAL_URL = "https://stev.oapd.inaf.it/cgi-bin/trilegal_1.6"
TRILEGAL_BASE = "https://stev.oapd.inaf.it/"
# the intermediate the server omits (Authority Information Access "CA Issuers" URI of its leaf)
INTERMEDIATE_URL = "http://crt.sectigo.com/ZeroSSLRSADVSSLCA2.crt"
DATA_DIR = Path(os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data"))
CERT_DIR = DATA_DIR / "followup" / "_certs"
DONE_LINE = "#TRILEGAL normally terminated"
# columns trilegal_results() reads
TRILEGAL_COLUMNS = ["Mact", "logg", "logTe", "[M/H]", "TESS", "J", "H", "Ks"]


def ca_bundle(cert_dir: Path = CERT_DIR) -> str:
    """certifi's CA bundle plus the ZeroSSL intermediate, built once and cached."""
    import certifi
    import requests
    from cryptography import x509
    from cryptography.hazmat.primitives.serialization import Encoding

    bundle = cert_dir / "ca_bundle_with_zerossl.pem"
    if bundle.exists():
        return str(bundle)
    cert_dir.mkdir(parents=True, exist_ok=True)
    der = requests.get(INTERMEDIATE_URL, timeout=30).content
    # the intermediate is fetched over plain http (as browsers do); it is only accepted if it
    # chains to a root already in certifi, which TLS checks on every connection
    pem = x509.load_der_x509_certificate(der).public_bytes(Encoding.PEM).decode()
    bundle.write_text(Path(certifi.where()).read_text() + "\n" + pem)
    return str(bundle)


def parse_trilegal(text: str) -> pd.DataFrame:
    """TRILEGAL's whitespace table (header line starts with '#') -> DataFrame.

    Raises if the run did not terminate normally or lacks the columns triceratops needs.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[-1].startswith(DONE_LINE):
        raise RuntimeError("TRILEGAL output incomplete (no termination line)")
    header = lines[0].lstrip("#").split()
    rows = [ln.split() for ln in lines[1:] if not ln.startswith("#")]
    df = pd.DataFrame(rows, columns=header).apply(pd.to_numeric, errors="coerce")
    missing = [c for c in TRILEGAL_COLUMNS if c not in df.columns]
    if missing:
        raise RuntimeError(f"TRILEGAL output lacks columns {missing}")
    return df


def write_triceratops_csv(df: pd.DataFrame, path: Path) -> None:
    """Write the population as triceratops expects: it reads `read_csv(f)[:-2]`, because its own
    writer keeps two trailing non-star rows (the termination comment). Add two padding rows."""
    pad = pd.DataFrame([[float("nan")] * df.shape[1]] * 2, columns=df.columns)
    pd.concat([df, pad], ignore_index=True).to_csv(path)


def query_trilegal(ra: float, dec: float, field_deg2: float = 0.1, mag_lim: float = 21.0,
                   timeout_s: float = 600.0) -> pd.DataFrame:
    """Submit the TRILEGAL 1.6 form (same settings as triceratops) and wait for the result."""
    import requests
    from mechanicalsoup import StatefulBrowser

    verify = ca_bundle()
    br = StatefulBrowser(soup_config={"features": "html.parser"})
    br.session.verify = verify
    br.open(TRILEGAL_URL, timeout=60)
    br.select_form(nr=0)
    br["gal_coord"] = "2"
    br["eq_alpha"] = str(ra)
    br["eq_delta"] = str(dec)
    br["field"] = str(field_deg2)
    br["photsys_file"] = "tab_mag_odfnew/tab_mag_TESS_2mass.dat"
    br["icm_lim"] = "1"          # magnitude limit applies to the 1st filter, TESS
    br["mag_lim"] = str(mag_lim)
    br["binary_kind"] = "0"      # no unresolved binaries in the population (triceratops' choice)
    br.submit_selected()
    links = [a.get("href") for a in br.get_current_page().select("a") if a.get("href")]
    out = [h for h in links if h.endswith(".dat")] or links[:1]
    if not out:
        raise RuntimeError("TRILEGAL returned no result link")
    url = TRILEGAL_BASE + out[0].lstrip("./")
    t0 = time.time()
    while True:
        r = requests.get(url, verify=verify, timeout=60)
        if r.ok and DONE_LINE in r.text[-200:]:
            return parse_trilegal(r.text)
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"TRILEGAL result not ready after {timeout_s:.0f} s: {url}")
        time.sleep(15)


def trilegal_file(tic_id: int, ra: float, dec: float, out_dir: Path) -> str:
    """Path to the cached TRILEGAL population for this target, querying it if needed."""
    path = out_dir / "trilegal.csv"
    if not path.exists():
        log.info("TIC %d: querying TRILEGAL at (%.4f, %.4f)", tic_id, ra, dec)
        df = query_trilegal(ra, dec)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_triceratops_csv(df, path)
        log.info("TIC %d: TRILEGAL population of %d stars saved", tic_id, len(df))
    return str(path)
