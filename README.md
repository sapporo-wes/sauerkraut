# sauerkraut

## List of test data
- **reference_hg38:** Reference fasta file, index files, and resource files according to human genome assembly hg38.
- **reference_hg19:** Reference fasta file, index files, and resource files according to human genome assembly hg19.
- **germlineWGS_hg38:** Germline whole genome sequence (WGS) data mapped on hg38.
- **somaticWGS_hg38:** Somatic WGS data mapped on hg38.
- **somaticCNV_hg19:** Somatic copy number variation (CNV) data mapped on hg19.
- **germlineRNA_hg19:** Germline RNA-Seq data mapped on hg19.

## How to create test data

### reference_hg38
Download reference file, index files, and resource files from the URLs listed in **[reference_hg38.download_links.txt](./download_links/reference_hg38.download_links.txt)**. 
You can download the input files by executing the following commnd:
```
$ OUTDIR=reference_hg38 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg38.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```


### reference_hg19


### germlineWGS_hg38



### somaticWGS_hg38


### somaticSNV_hg19


### germlineRNA_hg19





