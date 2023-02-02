# sauerkraut

## List of test data
- **`reference_hg38`:** Reference fasta file, index files, and resource files according to human genome assembly hg38.
- **`reference_hg19`:** Reference fasta file, index files, and resource files according to human genome assembly hg19.
- **`germlineWGS_hg38`:** Germline whole genome sequence (WGS) data mapped on hg38.
- **`somaticWGS_hg38`:** Somatic WGS data mapped on hg38.
- **`somaticCNV_hg19`:** Somatic copy number variation (CNV) data mapped on hg19.
- **`germlineRNA_hg19`:** Germline RNA-Seq data mapped on hg19.

## Steps to create test data

### Step 1. Install `cwltool`

You can install `cwltool` by executing the following commands: 
```
$ cd /path/to/working/directory/
$ python -m venv cwlenv
$ . cwlenv/bin/activate
$ python -m pip install --upgrade pip
```

Please confirm the `cwltool` version:
```
$ cwltool --version
/path/to/working/directory/cwlenv/bin/cwltool 3.1.20230201224320
```

### Step 2. Clone this repository

Clone this repository by executing the following commands:
```
$ cd /path/to/working/directory/
$ git clone https://github.com/sapporo-wes/sauerkraut.git
```


### Step 3. Download and create test data related to `reference_hg38`
Download reference file, index files, and resource files from the URLs listed in **[reference_hg38.download_links.txt](./download_links/reference_hg38.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg38 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg38.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

Download structural variation (SV)-related resource files from the URLs listed in **[reference_hg38_sv.download_links.txt](./download_links/reference_hg38_sv.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg38/sv ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg38_sv.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

Download mitochondria (MT)-related resource files from the URLs listed in **[reference_hg38_mt.download_links.txt](./download_links/reference_hg38_mt.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg38/mt ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg38_mt.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

Calculate bwa index files and create small test data by executing the following commands:
```
$ cd /path/to/working/directory/
$ WORKDIR=$(pwd)
$ /bin/sh sauerkraut/scripts/reference_hg38.create_test_data.sh \
     $WORKDIR/sauerkraut \
     $WORKDIR/cwlenv \
     $WORKDIR/reference_hg38
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**


### reference_hg19


### germlineWGS_hg38



### somaticWGS_hg38


### somaticSNV_hg19


### germlineRNA_hg19





