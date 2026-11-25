#!/bin/bash
# Download all files listed in a download_links text file.
# Files that already exist in the output directory are skipped.
set -eu

usage() {
  echo "usage: download.sh <path to download_links.txt> [output directory]" 1>&2
  echo "  The output directory defaults to the dataset name of the links file," 1>&2
  echo "  e.g. reference_hg38 for reference_hg38.download_links.txt" 1>&2
}

PATH_TO_DOWNLOAD_LINKS_TXT=${1:-}
if [ -z "$PATH_TO_DOWNLOAD_LINKS_TXT" ] || [ ! -f "$PATH_TO_DOWNLOAD_LINKS_TXT" ] ; then
  echo "ERROR: download links text is missing or not specified" 1>&2
  usage
  exit 1
fi

OUTDIR=${2:-$(basename "$PATH_TO_DOWNLOAD_LINKS_TXT" | cut -d. -f 1)}
mkdir -p "$OUTDIR"

for url in $(cat "$PATH_TO_DOWNLOAD_LINKS_TXT") ; do
  echo "$url"
  file=$(basename "$url")
  if [ ! -f "${OUTDIR}/$file" ] ; then
    # Download under a temporary name so that an interrupted download is
    # never mistaken for a complete file on the next run.
    wget "$url" -O "${OUTDIR}/$file.part"
    mv "${OUTDIR}/$file.part" "${OUTDIR}/$file"
  fi
done
