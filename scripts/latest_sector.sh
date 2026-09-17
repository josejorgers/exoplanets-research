#!/usr/bin/env bash
# Find the newest TESS sector available for each light-curve product on MAST,
# and estimate whether the newest folder is still being uploaded.
#
# Usage:  scripts/latest_sector.sh [-p qlp|tess-spoc|tglc] [-n N] [-j JOBS]
#   -p  product to inspect in depth (default: qlp)
#   -n  how many previous sectors to compare against (default: 3)
#   -j  parallel curl jobs for the folder census (default: 8)
#
# How "upload state" is judged (heuristics, no official API exists):
#   1. DRN     – does https://archive.stsci.edu/tess/tess_drn.html list the sector?
#                MAST publishes Data Release Notes when a sector is released.
#   2. Census  – number of TIC-prefix folders (s0NNN/0000/XXXX/) compared with the
#                median of the previous N sectors. A sector still uploading has
#                noticeably fewer. Ratio >= 0.85 -> "complete", 0.5-0.85 -> "partial",
#                < 0.5 -> "uploading".
#   3. Bulk    – for tess-spoc, a *-bulk-dl.tar.gz in the sector folder means the
#                sector was packaged after upload finished.
set -euo pipefail
export LC_ALL=C   # "." as decimal separator in awk regardless of locale

PRODUCT=qlp; NPREV=3; JOBS=8
while getopts "p:n:j:h" o; do
  case $o in
    p) PRODUCT=$OPTARG ;; n) NPREV=$OPTARG ;; j) JOBS=$OPTARG ;;
    *) sed -n '2,20p' "$0"; exit 0 ;;
  esac
done

BASE="https://archive.stsci.edu/hlsps"
fetch() { curl -sSL --retry 3 --retry-delay 2 "$1"; }

# --- sector folders for a product (4-digit names, real sectors only: <= 999) ---
sectors_of() {
  fetch "$BASE/$1/" | grep -oE 'href="s[0-9]{4}/"' | grep -oE '[0-9]{4}' \
    | awk '$1+0 <= 999' | sort -un
}

# --- count TIC-prefix folders under s0NNN/0000/ (parallel) ---
census() {  # $1=product $2=sector(4 digits)
  local url="$BASE/$1/s$2/0000/"
  [[ $1 == tess-spoc ]] && url="$BASE/$1/s$2/target/0000/"
  fetch "$url" | grep -oE 'href="[0-9]{4}/"' | grep -oE '[0-9]{4}' \
    | xargs -P "$JOBS" -I{} sh -c "curl -sSL --retry 3 '$url{}/' | grep -cE 'href=\"[0-9]{4}/\"' || true" \
    | awk '{s+=$1} END{print s+0}'
}

echo "== Newest sector per product (MAST HLSP folder listings)"
for p in qlp tess-spoc tglc; do
  printf "  %-10s %s\n" "$p" "$(sectors_of "$p" | tail -1 | sed 's/^0*//')"
done
spoc=$(fetch "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html" \
       | grep -oE 'tesscurl_sector_[0-9]+_lc\.sh' | grep -oE '[0-9]+' | awk '$1<=999' | sort -n | tail -1)
printf "  %-10s %s\n" "spoc-2min" "$spoc"
drn_max=$(fetch "https://archive.stsci.edu/tess/tess_drn.html" | grep -oiE 'sector[ _]?[0-9]+' \
          | grep -oE '[0-9]+' | awk '$1<=999' | sort -n | tail -1)
printf "  %-10s %s\n" "DRN" "$drn_max"

echo
echo "== Upload state of newest $PRODUCT sector"
mapfile -t secs < <(sectors_of "$PRODUCT")
newest=${secs[-1]}
prev=("${secs[@]: -$((NPREV+1)):$NPREV}")
echo "  newest: s$newest   comparing with: ${prev[*]/#/s}"

n_new=$(census "$PRODUCT" "$newest")
declare -a n_prev=()
for s in "${prev[@]}"; do n_prev+=("$(census "$PRODUCT" "$s")"); done
median=$(printf "%s\n" "${n_prev[@]}" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
printf "  TIC-prefix folders: newest=%s  previous=%s  median=%s\n" "$n_new" "${n_prev[*]}" "$median"

ratio=$(awk -v a="$n_new" -v b="$median" 'BEGIN{ if(b==0) print 0; else printf "%.2f", a/b }')
state=$(awk -v r="$ratio" 'BEGIN{ if(r>=0.85) print "complete"; else if(r>=0.5) print "partial"; else print "uploading" }')
if (( 10#$newest <= drn_max )); then drn_txt="yes"; else drn_txt="no"; fi
printf "  census ratio: %s -> %s\n" "$ratio" "$state"
printf "  DRN published for sector %s: %s\n" "$((10#$newest))" "$drn_txt"
if [[ $PRODUCT == tess-spoc ]]; then
  bulk=$(fetch "$BASE/$PRODUCT/s$newest/" | grep -c 'bulk-dl.tar.gz' || true)
  printf "  bulk tarball present: %s\n" "$([[ $bulk -gt 0 ]] && echo yes || echo no)"
fi

echo
if [[ $state == complete && $drn_txt == yes ]]; then
  echo "RECOMMENDED: sector $((10#$newest)) ($PRODUCT)"
else
  fallback=${prev[-1]}
  echo "Newest sector looks $state / DRN=$drn_txt. RECOMMENDED: sector $((10#$fallback)) ($PRODUCT); re-check later for $((10#$newest))."
fi
