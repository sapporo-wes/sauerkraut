#!/bin/bash
# Tests for scripts/download.sh and the download commands in README.md.
#
# usage: tests/test_download.sh
#
# Nothing is fetched from the internet. Most tests replace wget with a stub;
# the last one runs the real wget against a local HTTP server and is skipped
# when wget or python3 is not installed.
set -u

REPO=$(cd "$(dirname "$0")/.." && pwd)
DOWNLOAD="$REPO/scripts/download.sh"
REAL_PATH=$PATH
TMP=$(mktemp -d)
SERVER_PID=""
PASS=0
FAIL=0

cleanup() {
  if [ -n "$SERVER_PID" ] ; then
    kill "$SERVER_PID" 2> /dev/null
    wait "$SERVER_PID" 2> /dev/null
  fi
  rm -rf "$TMP"
}
trap cleanup EXIT

ok()   { PASS=$((PASS + 1)) ; echo "ok   - $1" ; }
ng()   { FAIL=$((FAIL + 1)) ; echo "FAIL - $1" ; }
check() { # check <description> <command...>
  local desc=$1 ; shift
  if "$@" ; then ok "$desc" ; else ng "$desc" ; fi
}

# A stub wget that understands "wget <url> -O <path>". It logs every URL to
# $WGET_LOG. URLs containing "FAIL" leave a truncated file behind and return
# an error, as an interrupted download would, unless $WGET_NO_FAIL is set.
mkdir -p "$TMP/stub"
cat > "$TMP/stub/wget" <<'EOF'
#!/bin/sh
url=$1 ; out=$3
echo "$url" >> "$WGET_LOG"
case "$url" in
  *FAIL*) if [ -z "${WGET_NO_FAIL:-}" ] ; then echo "trunc" > "$out" ; exit 1 ; fi ;;
esac
echo "content of $url" > "$out"
EOF
chmod +x "$TMP/stub/wget"

# Start each test in a fresh working directory with the stub wget on PATH.
setup() {
  WORK="$TMP/work.$1"
  mkdir -p "$WORK"
  cd "$WORK"
  export WGET_LOG="$WORK/wget.log"
  : > "$WGET_LOG"
  unset WGET_NO_FAIL
  PATH="$TMP/stub:$REAL_PATH"
}

count_files() { find "$1" -type f | wc -l | tr -d ' ' ; }
count_lines() { grep -c . "$1" ; }

# ---------------------------------------------------------------------------
echo "# argument handling"
setup args

"$DOWNLOAD" > out.txt 2> err.txt ; rc=$?
check "no argument exits with 1" [ $rc -eq 1 ]
check "no argument prints usage to stderr" grep -q "usage:" err.txt

"$DOWNLOAD" does_not_exist.txt > out.txt 2> err.txt ; rc=$?
check "nonexistent links file exits with 1" [ $rc -eq 1 ]
check "nonexistent links file downloads nothing" [ "$(count_lines "$WGET_LOG")" -eq 0 ]

# ---------------------------------------------------------------------------
echo "# output directory"
setup outdir
printf '%s\n' "https://example.org/a/one.fasta" "https://example.org/b/two.vcf.gz" \
  > mydata_hg38.download_links.txt

"$DOWNLOAD" mydata_hg38.download_links.txt > /dev/null 2>&1 ; rc=$?
check "default run exits with 0" [ $rc -eq 0 ]
check "default output directory is the dataset name" [ -f mydata_hg38/one.fasta ]
check "file name is the basename of the URL" [ -f mydata_hg38/two.vcf.gz ]
check "file holds what wget fetched for its URL" \
  grep -q "content of https://example.org/b/two.vcf.gz" mydata_hg38/two.vcf.gz

"$DOWNLOAD" mydata_hg38.download_links.txt nested/out/dir > /dev/null 2>&1 ; rc=$?
check "explicit nested output directory is created" [ -f nested/out/dir/one.fasta ]
check "explicit output directory gets every file" [ "$(count_files nested/out/dir)" -eq 2 ]

printf '%s\n' "https://example.org/one.txt" "" "https://example.org/two.txt" "" > blank.download_links.txt
"$DOWNLOAD" blank.download_links.txt > /dev/null 2>&1 ; rc=$?
check "blank lines in the links file are ignored" [ $rc -eq 0 -a "$(count_files blank)" -eq 2 ]

# ---------------------------------------------------------------------------
echo "# re-run"
setup rerun
printf '%s\n' "https://example.org/one.txt" "https://example.org/two.txt" > d.download_links.txt
mkdir d
echo "already here" > d/one.txt

"$DOWNLOAD" d.download_links.txt > /dev/null 2>&1
check "existing file is not downloaded again" [ "$(count_lines "$WGET_LOG")" -eq 1 ]
check "existing file is left untouched" grep -q "already here" d/one.txt
check "missing file is downloaded" [ -f d/two.txt ]

: > "$WGET_LOG"
"$DOWNLOAD" d.download_links.txt > /dev/null 2>&1
check "second run downloads nothing" [ "$(count_lines "$WGET_LOG")" -eq 0 ]

# ---------------------------------------------------------------------------
echo "# failed download"
setup failure
printf '%s\n' "https://example.org/one.txt" "https://example.org/FAIL.bam" "https://example.org/three.txt" \
  > f.download_links.txt

"$DOWNLOAD" f.download_links.txt > /dev/null 2>&1 ; rc=$?
check "failed download exits with non-zero" [ $rc -ne 0 ]
check "failed download leaves no file under the final name" [ ! -e f/FAIL.bam ]

: > "$WGET_LOG"
export WGET_NO_FAIL=1
"$DOWNLOAD" f.download_links.txt > /dev/null 2>&1 ; rc=$?
check "re-run after a failure exits with 0" [ $rc -eq 0 ]
check "re-run retries the failed file" grep -q "FAIL.bam" "$WGET_LOG"
check "re-run completes the failed file" grep -q "content of" f/FAIL.bam
check "re-run fetches the files after the failed one" [ -f f/three.txt ]
check "re-run skips the file that had succeeded" sh -c "! grep -q one.txt '$WGET_LOG'"
check "no temporary .part file is left after success" [ -z "$(find f -name '*.part')" ]

# ---------------------------------------------------------------------------
echo "# README commands"
setup readme
# The README tells the user to clone the repository into the working directory.
ln -s "$REPO" sauerkraut
grep -E '^\$ sauerkraut/scripts/download\.sh ' "$REPO/README.md" | sed 's/^\$ //' > commands.txt

n_commands=$(count_lines commands.txt)
n_links=$(ls "$REPO"/download_links/*.download_links.txt | wc -l | tr -d ' ')
check "README has download commands" [ "$n_commands" -gt 0 ]
check "README has one command per links file ($n_links)" [ "$n_commands" -eq "$n_links" ]

while read -r cmd links outdir ; do
  name=$(basename "$links" .download_links.txt)
  check "README [$name]: links file exists" [ -f "$links" ]
  $cmd "$links" "$outdir" > /dev/null 2>&1 ; rc=$?
  check "README [$name]: command exits with 0" [ $rc -eq 0 ]
  while read -r url ; do
    [ -n "$url" ] || continue
    [ -f "$outdir/$(basename "$url")" ] || { ng "README [$name]: $outdir/$(basename "$url") is missing" ; continue 2 ; }
  done < "$links"
  ok "README [$name]: every listed file is in $outdir"
done < commands.txt

# Inputs the create_test_data scripts look for, at the paths the README passes.
for f in reference_hg38/Homo_sapiens_assembly38.fasta \
         reference_hg19/Homo_sapiens_assembly19.fasta ; do
  check "README layout provides $f" [ -f "$f" ]
done

check "README has no inline wget loop left" sh -c "! grep -q 'wget \$url' '$REPO/README.md'"

# ---------------------------------------------------------------------------
echo "# real wget against a local HTTP server"
setup real
PATH=$REAL_PATH
if ! command -v wget > /dev/null || ! command -v python3 > /dev/null ; then
  echo "skip - wget or python3 is not installed"
else
  mkdir -p served/sub
  echo "alpha" > served/alpha.txt
  head -c 200000 /dev/urandom > served/sub/beta.bin
  port=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1])')
  python3 -m http.server "$port" --bind 127.0.0.1 --directory served > server.log 2>&1 &
  SERVER_PID=$!
  for _ in $(seq 50) ; do
    wget -q -O /dev/null "http://127.0.0.1:$port/alpha.txt" 2> /dev/null && break
    sleep 0.1
  done

  printf '%s\n' "http://127.0.0.1:$port/alpha.txt" "http://127.0.0.1:$port/sub/beta.bin" \
    > real.download_links.txt
  "$DOWNLOAD" real.download_links.txt > /dev/null 2>&1 ; rc=$?
  check "real wget: exits with 0" [ $rc -eq 0 ]
  check "real wget: text file is identical" cmp -s served/alpha.txt real/alpha.txt
  check "real wget: binary file is identical" cmp -s served/sub/beta.bin real/beta.bin

  printf '%s\n' "http://127.0.0.1:$port/alpha.txt" "http://127.0.0.1:$port/no_such_file.txt" \
    > broken.download_links.txt
  "$DOWNLOAD" broken.download_links.txt > /dev/null 2>&1 ; rc=$?
  check "real wget: HTTP 404 exits with non-zero" [ $rc -ne 0 ]
  check "real wget: HTTP 404 leaves no file under the final name" [ ! -e broken/no_such_file.txt ]
fi

# ---------------------------------------------------------------------------
echo ""
echo "passed: $PASS, failed: $FAIL"
[ "$FAIL" -eq 0 ]
