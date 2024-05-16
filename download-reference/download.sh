#!/bin/bash
set -e
PATH_TO_DOWNLOAD_LINKS_TXT=$1
if [ -z $PATH_TO_DOWNLOAD_LINKS_TXT ] || [ ! -e $PATH_TO_DONWLOAD_LINKS_TXT ] ; then
  echo "ERROR: download links text is missing or not specified" 1>&2
  echo "usage: ./download.sh <path to download_links.txt>" 1>&2
  exit 1
fi

set -u
OUTDIR=$(basename $PATH_TO_DOWNLOAD_LINKS_TXT | cut -d. -f 1)
mkdir -p $OUTDIR
for url in `cat $PATH_TO_DOWNLOAD_LINKS_TXT` ; do
  echo $url
  file=`basename $url`
  if [ ! -f ${OUTDIR}/$file ] ; then
    wget $url -O ${OUTDIR}/$file
  fi
done
