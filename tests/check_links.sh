#!/bin/bash
# Check that every URL in the download_links files is still alive.
#
# usage: tests/check_links.sh [download_links.txt ...]
#
# Without arguments, every download_links/*.download_links.txt is checked.
# Only headers are requested, so no data file is downloaded. Exits with 1 if
# any URL is dead.
#
# Environment variables:
#   JOBS     number of URLs checked in parallel (default: 8)
#   TIMEOUT  seconds to wait for one request (default: 60)
set -u

REPO=$(cd "$(dirname "$0")/.." && pwd)
JOBS=${JOBS:-8}
TIMEOUT=${TIMEOUT:-60}
export TIMEOUT

# Check one URL and print "<OK|DEAD>\t<size in bytes or ->\t<detail>\t<url>".
check_one() {
  local url=$1 headers size rc detail
  headers=$(curl --silent --show-error --head --location --fail \
                 --connect-timeout 20 --max-time "$TIMEOUT" --retry 2 \
                 "$url" 2>&1)
  rc=$?
  if [ $rc -ne 0 ] ; then
    # Some servers refuse HEAD requests. Ask for the first byte instead.
    headers=$(curl --silent --show-error --location --fail --range 0-0 \
                   --connect-timeout 20 --max-time "$TIMEOUT" --retry 2 \
                   --dump-header - --output /dev/null "$url" 2>&1)
    rc=$?
  fi
  if [ $rc -eq 0 ] ; then
    # With --range the total size is the part after "/" in Content-Range.
    size=$(printf '%s\n' "$headers" | tr -d '\r' \
      | awk 'tolower($1) == "content-range:" { sub(/.*\//, "", $0) ; r = $0 }
             tolower($1) == "content-length:" { l = $2 }
             END { if (r != "") print r ; else if (l != "") print l ; else print "-" }')
    printf 'OK\t%s\t-\t%s\n' "$size" "$url"
  else
    detail=$(printf '%s\n' "$headers" | grep '^curl:' | tail -n 1 | tr '\t' ' ')
    printf 'DEAD\t-\t%s\t%s\n' "${detail:-curl exit code $rc}" "$url"
  fi
}

if [ "${1:-}" = "--check-one" ] ; then
  check_one "$2"
  exit 0
fi

if ! command -v curl > /dev/null ; then
  echo "ERROR: curl is required" 1>&2
  exit 2
fi

if [ $# -gt 0 ] ; then
  FILES=("$@")
else
  FILES=("$REPO"/download_links/*.download_links.txt)
fi
for f in "${FILES[@]}" ; do
  if [ ! -f "$f" ] ; then
    echo "ERROR: $f is not found" 1>&2
    exit 2
  fi
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Each distinct URL is requested once, even if several files list it.
cat "${FILES[@]}" | grep . | sort -u > "$TMP/urls.txt"
echo "Checking $(wc -l < "$TMP/urls.txt" | tr -d ' ') URLs from ${#FILES[@]} files ..." 1>&2
xargs -P "$JOBS" -n 1 "$0" --check-one < "$TMP/urls.txt" > "$TMP/first.tsv"

# Servers that limit simultaneous connections can time out under parallel
# requests. Check the failed URLs once more, one at a time.
grep '^OK' "$TMP/first.tsv" > "$TMP/results.tsv"
grep '^DEAD' "$TMP/first.tsv" | cut -f 4 > "$TMP/retry.txt"
if [ -s "$TMP/retry.txt" ] ; then
  echo "Checking $(wc -l < "$TMP/retry.txt" | tr -d ' ') failed URLs again ..." 1>&2
  while read -r url ; do
    check_one "$url" >> "$TMP/results.tsv"
  done < "$TMP/retry.txt"
fi

n_ok=0
n_dead=0
for f in "${FILES[@]}" ; do
  echo ""
  echo "# $(basename "$f")"
  while read -r url ; do
    [ -n "$url" ] || continue
    line=$(awk -F '\t' -v u="$url" '$4 == u' "$TMP/results.tsv" | head -n 1)
    status=$(printf '%s' "$line" | cut -f 1)
    size=$(printf '%s' "$line" | cut -f 2)
    detail=$(printf '%s' "$line" | cut -f 3)
    if [ "$status" = "OK" ] ; then
      n_ok=$((n_ok + 1))
      printf 'ok   %14s  %s\n' "$size" "$url"
    else
      n_dead=$((n_dead + 1))
      printf 'DEAD %14s  %s\n' "-" "$url"
      echo "     ${detail:-no result}"
    fi
  done < "$f"
  # The same file name twice in one list would be downloaded only once.
  dup=$(grep . "$f" | sed 's|.*/||' | sort | uniq -d)
  if [ -n "$dup" ] ; then
    echo "WARN duplicated file name: $(echo $dup)"
  fi
done

echo ""
echo "alive: $n_ok, dead: $n_dead"
[ "$n_dead" -eq 0 ]
