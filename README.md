# sauerkraut

## About this project

This repository contains source codes to download and create the following test datasets:
- **`reference_hg38`:** Reference fasta file, index files, and resource files according to human genome assembly hg38.
     - Our program will download sequence, index, and resource files for a full reference dataset from the URLs listed in **[reference_hg38.download_links.txt](./download_links/reference_hg38.download_links.txt)**. Then, our program will calculate additional index files for the full reference as well as squence, index, and resource files for a small reference dataset. 
     - **`Homo_sapiens_assembly38.fasta`:** A sequence file for the full reference dataset. The file size in a Linux environment was `3249912778`. Index files for the full reference (`.fasta.fai`, `.fasta.alt`, and `.dict`) will be downloaded. Additional index files for the full reference (`.fasta.amb`, `.fasta.ann`, `.fasta.bwt`, `.fasta.pac`, and `.fasta.sa`) will be calculated. Resource files for the full reference dataset such as `Homo_sapiens_assembly38.dbsnp138.vcf` and `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` will be downloaded.
     - **`Homo_sapiens_assembly38.small.fasta`:** A sequence file for the small reference dataset. The file size in a Linux environment was `202524`. This small file contains reference sequence of the region `chr1:1000001-1200000`. Index and resource files for the small reference will be calculated. 
- **`reference_hg19`:** Reference fasta file, index files, and resource files according to human genome assembly hg19.
     - Our program will download sequence, index, and resource files for a full reference dataset from the URLs listed in **[reference_hg19.download_links.txt](./download_links/reference_hg19.download_links.txt)**. 
     - **`Homo_sapiens_assembly19.fasta`:** A sequence file for the full reference dataset. The file size in a Linux environment was `3140756381`. Index files for the full reference (`.fasta.fai` and `.dict`) will be downloaded. Resource files for the full reference dataset such as `dbsnp_138.b37.vcf.gz` and `Mills_and_1000G_gold_standard.indels.b37.vcf.gz` will be downloaded. 
- **`germlineWGS_hg38`:** Germline whole genome sequence (WGS) data mapped on hg38.
     - Our program will download germline WGS data in the `cram` format (10 samples: `HG00446`, `NA12718`, `NA12775`, `NA12842`, `NA18536`, `NA18549`, `NA20752`, `NA20757`, `NA20764`, and `NA20769`) from the URLs listed in **[germlineWGS_hg38.download_links.txt](./download_links/germlineWGS_hg38.download_links.txt)**. Then our program will calculate `bam`, `unmap.bam`, and `fastq.gz` (paired and interleaved) files from the `cram` file. In addition, our program will calculate small `cram`, `bam`, `unmap.bam`, and `fastq.gz` (paired and interleaved) files by extracting sequence reads mapped onto the region `chr1:1000001-1200000`. Mitochondria-specific small `cram`, `bam`, `unmap.bam`, and `fastq.gz` (paired and interleaved) files will be calculated by extracting sequence reads mapped onto the `chrM`. Index files (`.crai` and `.bai`) will be calculated for all `cram` and `bam` files.  
     - **`HG00446_full.cram`:** A `cram` file for the full dataset. The file size in a Linux environment was `15575753770`.
     - **`HG00446_full.bam`:** A `bam` file for the full dataset. The file size in a Linux environment was `41284415391`.
     - **`HG00446_full.unmap.bam`:** A `unmap.bam` file for the full dataset. The file size in a Linux environment was `46725838640`.
     - **`HG00446_full_1.fastq.gz` and `HG00446_full_2.fastq.gz`:** Paired `fastq.gz` files for the full dataset. The file sizes in a Linux environment were `15602892502` and `16931512190`, respectively.
     - **`HG00446_full_interleaved.fastq.gz`:** A interleaved `fastq.gz` file for the full dataset. The file size in a Linux environment was `22051523852`.
     - **`HG00446_small.cram`:** A `cram` file for the small dataset. The file size in a Linux environment was `1268473`.
     - **`HG00446_small.bam`:** A `bam` file for the small dataset. The file size in a Linux environment was `3393270`.
     - **`HG00446_small.unmap.bam`:** A `unmap.bam` file for the small dataset. The file size in a Linux environment was `3613702`.
     - **`HG00446_small_1.fastq.gz` and `HG00446_small_2.fastq.gz`:** Paired `fastq.gz` files for the small dataset. The file sizes in a Linux environment were `1285315` and `1404875`, respectively.
     - **`HG00446_small_interleaved.fastq.gz`:** A interleaved `fastq.gz` file for the small dataset. The file size in a Linux environment was `1874817`.
     - **`HG00446_mtsmall.cram`:** A `cram` file for the mitochondria-specific small dataset. The file size in a Linux environment was `28956734`.
     - **`HG00446_mtsmall.bam`:** A `bam` file for the mitochondria-specific small dataset. The file size in a Linux environment was `58621592`.
- **`somaticWGS_hg38`:** Somatic WGS data mapped on hg38.
     - Our program will download somatic WGS data in the `bam` format (a tumor samples `hcc1143_T` and a matched control sample `hcc1143_N`) from the URLs listed in **[somaticWGS_hg38.download_links.txt](./download_links/somaticWGS_hg38.download_links.txt)**. Then our program will calculate `cram`, `unmap.bam`, and `fastq.gz` (paired) files from the `bam` file. In addition, our program will calculate small `cram`, `bam`, `unmap.bam`, and `fastq.gz` (paired) files by extracting sequence reads mapped onto the region `chr1:1000001-1200000`. Index files (`.crai` and `.bai`) will be calculated for all `cram` and `bam` files.  
     - **`hcc1143_T_full.cram`:** A `cram` file for the full dataset. The file size in a Linux environment was `8488487437`.
     - **`hcc1143_T_full.bam`:** A `bam` file for the full dataset. The file size in a Linux environment was `13478097244`.
     - **`hcc1143_T_full.unmap.bam`:** A `unmap.bam` file for the full dataset. The file size in a Linux environment was `13863994168`.
     - **`hcc1143_T_full_1.fastq.gz` and `hcc1143_T_full_2.fastq.gz`:** Paired `fastq.gz` files for the full dataset. The file sizes in a Linux environment were `3038008245` and `2526825944`, respectively.
     - **`hcc1143_T_small.cram`:** A `cram` file for the small dataset. The file size in a Linux environment was `1068517`.
     - **`hcc1143_T_small.bam`:** A `bam` file for the small dataset. The file size in a Linux environment was `1698508`.
     - **`hcc1143_T_small.unmap.bam`:** A `unmap.bam` small for the full dataset. The file size in a Linux environment was `1829678`.
     - **`hcc1143_T_small_1.fastq.gz` and `hcc1143_T_small_2.fastq.gz`:** Paired `fastq.gz` files for the small dataset. The file sizes in a Linux environment were `416961` and `347292`, respectively.
- **`somaticCNV_hg19`:** Somatic copy number variation (CNV) data mapped on hg19.
     - Our program will download somatic WGS data used for CNV analysis in the `bam` format (a tumor samples `SM-74P4M` and a matched control sample `SM-74NEG`) from the URLs listed in **[somaticCNV_hg19.download_links.txt](./download_links/somaticCNV_hg19.download_links.txt)**. Then our program will calculate `cram`, `unmap.bam`, and `fastq.gz` (paired and interleaved) files from the `bam` file. In addition, our program will calculate small `cram`, `bam`, `unmap.bam`, and `fastq.gz` (paired and interleaved) files by extracting sequence reads mapped onto the region `1:1000001-1200000`. Index files (`.crai` and `.bai`) will be calculated for all `cram` and `bam` files.  
     - **`SM-74P4M_full.cram`:** A `cram` file for the full dataset. The file size in a Linux environment was `18274132829`.
     - **`SM-74P4M_full.bam`:** A `bam` file for the full dataset. The file size in a Linux environment was `30766633211`.
     - **`SM-74P4M_full.unmap.bam`:** A `unmap.bam` file for the full dataset. The file size in a Linux environment was `21238770786`.
     - **`SM-74P4M_full_1.fastq.gz` and `SM-74P4M_full_2.fastq.gz`:** Paired `fastq.gz` files for the full dataset. The file sizes in a Linux environment were `8966831711` and `8984898864`, respectively.
     - **`SM-74P4M_full_interleaved.fastq.gz`:** A interleaved `fastq.gz` file for the full dataset. The file size in a Linux environment was `14653053842`.
     - **`SM-74P4M_small.cram`:** A `cram` file for the small dataset. The file size in a Linux environment was `3288611`.
     - **`SM-74P4M_small.bam`:** A `bam` file for the small dataset. The file size in a Linux environment was `5206596`.
     - **`SM-74P4M_small.unmap.bam`:** A `unmap.bam` small for the full dataset. The file size in a Linux environment was `3544499`.
     - **`SM-74P4M_small_1.fastq.gz` and `SM-74P4M_small_2.fastq.gz`:** Paired `fastq.gz` files for the small dataset. The file sizes in a Linux environment were `1456295` and `1460187`, respectively.
     - **`SM-74P4M_small_interleaved.fastq.gz`:** A interleaved `fastq.gz` file for the small dataset. The file size in a Linux environment was `2381216`.
- **`germlineRNA_hg19`:** Germline RNA-Seq data mapped on hg19.
     - Our program will download germline RNA-Seq data in the `unmap.bam` format (a sample `NA12878`) from the URLs listed in **[germlineRNA_hg19.download_links.txt](./download_links/germlineRNA_hg19.download_links.txt)**. Then our program will calculate small `unmap.bam` file from by down-sampling sequence. Index files (`.bai`) will be calculated for the `unmap.bam` files.
     - **`NA12878_full.bam`:** A `unmap.bam` file for the full dataset. The file size in a Linux environment was `3321881612`.
     - **`NA12878_small.bam`:** A `unmap.bam` file for the small dataset. The file size in a Linux environment was `5313719`.

The above datasets are designed to use as test datasets for the following workflows:
- **`Data pre-processing (hg38)`:** Workflows to perform the data pre-processing for variant discovery from short-read WGS data.
     - Best practice for this step is described in [GATK web site (Data-pre-processing-for-variant-discovery)](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery). The raw sequence data (either in `fastq.gz (paired)` or `unmap.map` format) will be mapped to a reference genome, and subsequently, the alignment data will be cleaned up by removing PCR duplicates and calibrating base quality scores. 
- **`Germline short variant call from bam/cram (hg38)`:**
- **`Somatic short variant call from bm/cram (hg38)`:**
- **`Germline structural variation call from bam/cram (hg38)`:**
- **`Mitochondria short variant call from bam/cram (hg38)`:**
- **`Somatic CNV call from bam/cram (hg19)`:**
- **`Germline short variant call from RNA-Seq data (hg19)`:**


## Steps to create test data

### Requirements
To download and create test data, the following software is needed:
- python (version 3+)
- singularity
- wget

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

Download somatic-related resource files from the URLs listed in **[reference_hg38_somatic.download_links.txt](./download_links/reference_hg38_somatic.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg38/somatic ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg38_somatic.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
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


The following files will be downloaded or created:
```{reference_hg38}
/path/to/working/directory/reference_hg38
|--1000G_omni2.5.hg38.vcf.gz
|--1000G_omni2.5.hg38.vcf.gz.tbi
|--1000G_omni2.5.hg38.vcf.small.vcf.gz
|--1000G_omni2.5.hg38.vcf.small.vcf.gz.tbi
|--1000G_phase1.snps.high_confidence.hg38.vcf.gz
|--1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
|--1000G_phase1.snps.high_confidence.hg38.vcf.small.vcf.gz
|--1000G_phase1.snps.high_confidence.hg38.vcf.small.vcf.gz.tbi
|--Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
|--Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
|--Homo_sapiens_assembly38.contam.UD
|--Homo_sapiens_assembly38.contam.bed
|--Homo_sapiens_assembly38.contam.mu
|--Homo_sapiens_assembly38.dbsnp138.vcf
|--Homo_sapiens_assembly38.dbsnp138.vcf.gz
|--Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
|--Homo_sapiens_assembly38.dbsnp138.vcf.idx
|--Homo_sapiens_assembly38.dbsnp138.vcf.small.vcf
|--Homo_sapiens_assembly38.dbsnp138.vcf.small.vcf.idx
|--Homo_sapiens_assembly38.dict
|--Homo_sapiens_assembly38.fasta
|--Homo_sapiens_assembly38.fasta.alt
|--Homo_sapiens_assembly38.fasta.amb
|--Homo_sapiens_assembly38.fasta.ann
|--Homo_sapiens_assembly38.fasta.bwt
|--Homo_sapiens_assembly38.fasta.fai
|--Homo_sapiens_assembly38.fasta.pac
|--Homo_sapiens_assembly38.fasta.sa
|--Homo_sapiens_assembly38.haplotype_database.txt
|--Homo_sapiens_assembly38.known_indels.vcf.gz
|--Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
|--Homo_sapiens_assembly38.known_indels.vcf.small.vcf.gz
|--Homo_sapiens_assembly38.known_indels.vcf.small.vcf.gz.tbi
|--Homo_sapiens_assembly38.small.dict
|--Homo_sapiens_assembly38.small.fasta
|--Homo_sapiens_assembly38.small.fasta.alt
|--Homo_sapiens_assembly38.small.fasta.amb
|--Homo_sapiens_assembly38.small.fasta.ann
|--Homo_sapiens_assembly38.small.fasta.bwt
|--Homo_sapiens_assembly38.small.fasta.fai
|--Homo_sapiens_assembly38.small.fasta.pac
|--Homo_sapiens_assembly38.small.fasta.sa
|--Homo_sapiens_assembly38.small.interval_list
|--Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
|--Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
|--Mills_and_1000G_gold_standard.indels.hg38.vcf.small.vcf.gz
|--Mills_and_1000G_gold_standard.indels.hg38.vcf.small.vcf.gz.tbi
|--abralab_igvtools:v2.8.2.sif
|--biocontainers_bwa:v0.7.17_cv1.sif
|--broadinstitute_gatk:4.3.0.0.sif
|--hapmap_3.3.hg38.vcf.gz
|--hapmap_3.3.hg38.vcf.gz.tbi
|--hapmap_3.3.hg38.vcf.small.vcf.gz
|--hapmap_3.3.hg38.vcf.small.vcf.gz.tbi
|--hg38.even.handcurated.20k.intervals
|--mgibio_samtools:1.16.1.sif
|--mt
|  |--Homo_sapiens_assembly38.chrM.dict
|  |--Homo_sapiens_assembly38.chrM.fasta
|  |--Homo_sapiens_assembly38.chrM.fasta.amb
|  |--Homo_sapiens_assembly38.chrM.fasta.ann
|  |--Homo_sapiens_assembly38.chrM.fasta.bwt
|  |--Homo_sapiens_assembly38.chrM.fasta.fai
|  |--Homo_sapiens_assembly38.chrM.fasta.pac
|  |--Homo_sapiens_assembly38.chrM.fasta.sa
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac
|  |--Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa
|  |--ShiftBack.chain
|  |--blacklist_sites.hg38.chrM.bed
|  |--blacklist_sites.hg38.chrM.bed.idx
|  |--blacklist_sites.hg38.chrM.shifted_by_8000_bases.fixed.bed
|  |--blacklist_sites.hg38.chrM.shifted_by_8000_bases.fixed.bed.idx
|  |--chrMWithFinalNuMTs.hg38.interval_list
|  |--control_region_shifted.chrM.interval_list
|  |--non_control_region.interval_list
|--quay.io_biocontainers_bcftools:1.16--hfe4b78e_1.sif
|--somatic
|  |--1000g_pon.hg38.vcf.gz
|  |--1000g_pon.hg38.vcf.gz.tbi
|  |--Homo_sapiens_assembly38.index_bundle
|  |--af-only-gnomad.hg38.vcf.gz
|  |--af-only-gnomad.hg38.vcf.gz.tbi
|  |--funcotator_dataSources.v1.6.20190124s.tar.gz
|  |--small_exac_common_3.hg38.vcf.gz
|  |--small_exac_common_3.hg38.vcf.gz.tbi
|  |--transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt
|--sv
|  |--GRCh38_Nmask.bed
|  |--MANE.GRCh38.v0.95.select_ensembl_genomic.gtf
|  |--PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz
|  |--allosome.fai
|  |--autosome.fai
|  |--bin_exclude.hg38.gatkcov.bed.gz
|  |--contig.fai
|  |--cytobands_hg38.bed.gz
|  |--delly_human.hg38.excl.tsv
|  |--depth_blacklist.sorted.bed.gz
|  |--empty.file
|  |--gencode.lincRNA.gtf.gz
|  |--gnomad_v2.1_sv.sites.GRCh38.bed.gz
|  |--hg38.SD_gaps_Cen_Tel_Heter_Satellite_lumpy.blacklist.sorted.merged.bed.gz
|  |--hg38.contig_ploidy_priors_homo_sapiens.tsv
|  |--hg38.genome
|  |--hg38.randomForest_blacklist.withRepMask.bed.gz
|  |--hg38.wgs.blacklist.wPAR.bed
|  |--hg38_primary_contigs.bed
|  |--mei_hg38.bed.gz
|  |--melt_standard_vcf_header.txt
|  |--noncoding.sort.hg38.bed
|  |--preprocessed_intervals.interval_list
|  |--primary_contigs.list
|  |--primary_contigs_plus_mito.bed.gz
|  |--promoter.bed
|  |--seed_cutoff.txt
|  |--wgd_scoring_mask.hg38.gnomad_v3.bed
|  |--wham_whitelist.bed
|--wgs_calling_regions.hg38.interval_list
|--wgs_evaluation_regions.hg38.interval_list
```


### Step 4. Download and create test data related to `reference_hg19`
Download reference file, index files, and resource files from the URLs listed in **[reference_hg19.download_links.txt](./download_links/reference_hg19.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg19 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg19.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

Download somatic-related resource files from the URLs listed in **[reference_hg19_somatic.download_links.txt](./download_links/reference_hg19_somatic.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg19/somatic ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg19_somatic.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

Download RNA-seq related resource files from the URLs listed in **[reference_hg19_rna.download_links.txt](./download_links/reference_hg19_rna.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg19/rna ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/reference_hg19_rna.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

The following files will be downloaded or created:
```{reference_hg19}
/path/to/working/directory/reference_hg19
|--Homo_sapiens_assembly19.dict
|--Homo_sapiens_assembly19.fasta
|--Homo_sapiens_assembly19.fasta.fai
|--Homo_sapiens_assembly19.known_indels_20120518.vcf
|--Homo_sapiens_assembly19.known_indels_20120518.vcf.idx
|--Mills_and_1000G_gold_standard.indels.b37.vcf.gz
|--Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi
|--dbsnp_138.b37.vcf.gz
|--dbsnp_138.b37.vcf.gz.tbi
|--rna
|  |--star.gencode.v19.transcripts.patched_contigs.gtf
|--somatic
|  |--common_snps.interval_list
|  |--funcotator_dataSources.v1.6.20190124s.tar.gz
|  |--ice_targets.tsv.interval_list
|  |--transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt
|  |--wes-do-gc.pon.hdf5
```


### Step 5. Download and create test data related to `germlineWGS_hg38`
Download cram files from the URLs listed in **[germlineWGS_hg38.download_links.txt](./download_links/germlineWGS_hg38.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=germlineWGS_hg38 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/germlineWGS_hg38.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**

Create test data by executing the following commands:
```
$ cd /path/to/working/directory/
$ WORKDIR=$(pwd)
$ /bin/sh sauerkraut/scripts/germlineWGS_hg38.create_test_data.sh \
     $WORKDIR/sauerkraut \
     $WORKDIR/cwlenv \
     $WORKDIR/reference_hg38 \
     $WORKDIR/germlineWGS_hg38 \
     HG00446
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**
Repeat the above commands by replacing sample ID (`HG00446`) to other sample IDs:
- `NA12718`
- `NA12775`
- `NA12842`
- `NA18536`
- `NA18549`
- `NA20752`
- `NA20757`
- `NA20764`
- `NA20769`

The following files will be downloaded or created:
```{germlineWGS_hg38}
/path/to/working/directory/germlineWGS_hg38
|--HG00446.final.cram
|--HG00446_full.bam
|--HG00446_full.bam.bai
|--HG00446_full.cram
|--HG00446_full.cram.crai
|--HG00446_full.unmap.bam
|--HG00446_full.unmap.bam.bai
|--HG00446_full_1.fastq.gz
|--HG00446_full_2.fastq.gz
|--HG00446_full_interleaved.fastq.gz
|--HG00446_mtsmall.bam
|--HG00446_mtsmall.bam.bai
|--HG00446_mtsmall.cram
|--HG00446_mtsmall.cram.crai
|--HG00446_small.bam
|--HG00446_small.bam.bai
|--HG00446_small.cram
|--HG00446_small.cram.crai
|--HG00446_small.unmap.bam
|--HG00446_small.unmap.bam.bai
|--HG00446_small_1.fastq.gz
|--HG00446_small_2.fastq.gz
|--HG00446_small_interleaved.fastq.gz
|--NA12718.final.cram
|--NA12718_full.bam
|--NA12718_full.bam.bai
|--NA12718_full.cram
|--NA12718_full.cram.crai
|--NA12718_full.unmap.bam
|--NA12718_full.unmap.bam.bai
|--NA12718_full_1.fastq.gz
|--NA12718_full_2.fastq.gz
|--NA12718_full_interleaved.fastq.gz
|--NA12718_mtsmall.bam
|--NA12718_mtsmall.bam.bai
|--NA12718_mtsmall.cram
|--NA12718_mtsmall.cram.crai
|--NA12718_small.bam
|--NA12718_small.bam.bai
|--NA12718_small.cram
|--NA12718_small.cram.crai
|--NA12718_small.unmap.bam
|--NA12718_small.unmap.bam.bai
|--NA12718_small_1.fastq.gz
|--NA12718_small_2.fastq.gz
|--NA12718_small_interleaved.fastq.gz
|--NA12775.final.cram
|--NA12775_full.bam
|--NA12775_full.bam.bai
|--NA12775_full.cram
|--NA12775_full.cram.crai
|--NA12775_full.unmap.bam
|--NA12775_full.unmap.bam.bai
|--NA12775_full_1.fastq.gz
|--NA12775_full_2.fastq.gz
|--NA12775_full_interleaved.fastq.gz
|--NA12775_mtsmall.bam
|--NA12775_mtsmall.bam.bai
|--NA12775_mtsmall.cram
|--NA12775_mtsmall.cram.crai
|--NA12775_small.bam
|--NA12775_small.bam.bai
|--NA12775_small.cram
|--NA12775_small.cram.crai
|--NA12775_small.unmap.bam
|--NA12775_small.unmap.bam.bai
|--NA12775_small_1.fastq.gz
|--NA12775_small_2.fastq.gz
|--NA12775_small_interleaved.fastq.gz
|--NA12842.final.cram
|--NA12842_full.bam
|--NA12842_full.bam.bai
|--NA12842_full.cram
|--NA12842_full.cram.crai
|--NA12842_full.unmap.bam
|--NA12842_full.unmap.bam.bai
|--NA12842_full_1.fastq.gz
|--NA12842_full_2.fastq.gz
|--NA12842_full_interleaved.fastq.gz
|--NA12842_mtsmall.bam
|--NA12842_mtsmall.bam.bai
|--NA12842_mtsmall.cram
|--NA12842_mtsmall.cram.crai
|--NA12842_small.bam
|--NA12842_small.bam.bai
|--NA12842_small.cram
|--NA12842_small.cram.crai
|--NA12842_small.unmap.bam
|--NA12842_small.unmap.bam.bai
|--NA12842_small_1.fastq.gz
|--NA12842_small_2.fastq.gz
|--NA12842_small_interleaved.fastq.gz
|--NA18536.final.cram
|--NA18536_full.bam
|--NA18536_full.bam.bai
|--NA18536_full.cram
|--NA18536_full.cram.crai
|--NA18536_full.unmap.bam
|--NA18536_full.unmap.bam.bai
|--NA18536_full_1.fastq.gz
|--NA18536_full_2.fastq.gz
|--NA18536_full_interleaved.fastq.gz
|--NA18536_mtsmall.bam
|--NA18536_mtsmall.bam.bai
|--NA18536_mtsmall.cram
|--NA18536_mtsmall.cram.crai
|--NA18536_small.bam
|--NA18536_small.bam.bai
|--NA18536_small.cram
|--NA18536_small.cram.crai
|--NA18536_small.unmap.bam
|--NA18536_small.unmap.bam.bai
|--NA18536_small_1.fastq.gz
|--NA18536_small_2.fastq.gz
|--NA18536_small_interleaved.fastq.gz
|--NA18549.final.cram
|--NA18549_full.bam
|--NA18549_full.bam.bai
|--NA18549_full.cram
|--NA18549_full.cram.crai
|--NA18549_full.unmap.bam
|--NA18549_full.unmap.bam.bai
|--NA18549_full_1.fastq.gz
|--NA18549_full_2.fastq.gz
|--NA18549_full_interleaved.fastq.gz
|--NA18549_mtsmall.bam
|--NA18549_mtsmall.bam.bai
|--NA18549_mtsmall.cram
|--NA18549_mtsmall.cram.crai
|--NA18549_small.bam
|--NA18549_small.bam.bai
|--NA18549_small.cram
|--NA18549_small.cram.crai
|--NA18549_small.unmap.bam
|--NA18549_small.unmap.bam.bai
|--NA18549_small_1.fastq.gz
|--NA18549_small_2.fastq.gz
|--NA18549_small_interleaved.fastq.gz
|--NA20752.final.cram
|--NA20752_full.bam
|--NA20752_full.bam.bai
|--NA20752_full.cram
|--NA20752_full.cram.crai
|--NA20752_full.unmap.bam
|--NA20752_full.unmap.bam.bai
|--NA20752_full_1.fastq.gz
|--NA20752_full_2.fastq.gz
|--NA20752_full_interleaved.fastq.gz
|--NA20752_mtsmall.bam
|--NA20752_mtsmall.bam.bai
|--NA20752_mtsmall.cram
|--NA20752_mtsmall.cram.crai
|--NA20752_small.bam
|--NA20752_small.bam.bai
|--NA20752_small.cram
|--NA20752_small.cram.crai
|--NA20752_small.unmap.bam
|--NA20752_small.unmap.bam.bai
|--NA20752_small_1.fastq.gz
|--NA20752_small_2.fastq.gz
|--NA20752_small_interleaved.fastq.gz
|--NA20757.final.cram
|--NA20757_full.bam
|--NA20757_full.bam.bai
|--NA20757_full.cram
|--NA20757_full.cram.crai
|--NA20757_full.unmap.bam
|--NA20757_full.unmap.bam.bai
|--NA20757_full_1.fastq.gz
|--NA20757_full_2.fastq.gz
|--NA20757_full_interleaved.fastq.gz
|--NA20757_mtsmall.bam
|--NA20757_mtsmall.bam.bai
|--NA20757_mtsmall.cram
|--NA20757_mtsmall.cram.crai
|--NA20757_small.bam
|--NA20757_small.bam.bai
|--NA20757_small.cram
|--NA20757_small.cram.crai
|--NA20757_small.unmap.bam
|--NA20757_small.unmap.bam.bai
|--NA20757_small_1.fastq.gz
|--NA20757_small_2.fastq.gz
|--NA20757_small_interleaved.fastq.gz
|--NA20764.final.cram
|--NA20764_full.bam
|--NA20764_full.bam.bai
|--NA20764_full.cram
|--NA20764_full.cram.crai
|--NA20764_full.unmap.bam
|--NA20764_full.unmap.bam.bai
|--NA20764_full_1.fastq.gz
|--NA20764_full_2.fastq.gz
|--NA20764_full_interleaved.fastq.gz
|--NA20764_mtsmall.bam
|--NA20764_mtsmall.bam.bai
|--NA20764_mtsmall.cram
|--NA20764_mtsmall.cram.crai
|--NA20764_small.bam
|--NA20764_small.bam.bai
|--NA20764_small.cram
|--NA20764_small.cram.crai
|--NA20764_small.unmap.bam
|--NA20764_small.unmap.bam.bai
|--NA20764_small_1.fastq.gz
|--NA20764_small_2.fastq.gz
|--NA20764_small_interleaved.fastq.gz
|--NA20769.final.cram
|--NA20769_full.bam
|--NA20769_full.bam.bai
|--NA20769_full.cram
|--NA20769_full.cram.crai
|--NA20769_full.unmap.bam
|--NA20769_full.unmap.bam.bai
|--NA20769_full_1.fastq.gz
|--NA20769_full_2.fastq.gz
|--NA20769_full_interleaved.fastq.gz
|--NA20769_mtsmall.bam
|--NA20769_mtsmall.bam.bai
|--NA20769_mtsmall.cram
|--NA20769_mtsmall.cram.crai
|--NA20769_small.bam
|--NA20769_small.bam.bai
|--NA20769_small.cram
|--NA20769_small.cram.crai
|--NA20769_small.unmap.bam
|--NA20769_small.unmap.bam.bai
|--NA20769_small_1.fastq.gz
|--NA20769_small_2.fastq.gz
|--NA20769_small_interleaved.fastq.gz
|--broadinstitute_gatk:4.3.0.0.sif
|--mgibio_samtools:1.16.1.sif
|--quay.io_biocontainers_bbmap:39.01--h5c4e2a8_0.sif
```


### Step 6. Download and create test data related to `somaticWGS_hg38`
Download bam files from the URLs listed in **[somaticWGS_hg38.download_links.txt](./download_links/somaticWGS_hg38.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=somaticWGS_hg38 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/somaticWGS_hg38.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**

Create test data by executing the following commands:
```
$ cd /path/to/working/directory/
$ WORKDIR=$(pwd)
$ /bin/sh sauerkraut/scripts/somaticWGS_hg38.create_test_data.sh \
     $WORKDIR/sauerkraut \
     $WORKDIR/cwlenv \
     $WORKDIR/reference_hg38 \
     $WORKDIR/somaticWGS_hg38 \
     hcc1143_N
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**
Repeat the above commands by replacing sample ID (`hcc1143_N`) to other sample ID:
- `hcc1143_T`

The following files will be downloaded or created:
```{somaticWGS_hg38}
/path/to/working/directory/somaticWGS_hg38
|--broadinstitute_gatk:4.3.0.0.sif
|--hcc1143_N_clean.bam
|--hcc1143_N_full.bam
|--hcc1143_N_full.bam.bai
|--hcc1143_N_full.cram
|--hcc1143_N_full.cram.crai
|--hcc1143_N_full.unmap.bam
|--hcc1143_N_full.unmap.bam.bai
|--hcc1143_N_full_1.fastq.gz
|--hcc1143_N_full_2.fastq.gz
|--hcc1143_N_small.bam
|--hcc1143_N_small.bam.bai
|--hcc1143_N_small.cram
|--hcc1143_N_small.cram.crai
|--hcc1143_N_small.unmap.bam
|--hcc1143_N_small.unmap.bam.bai
|--hcc1143_N_small_1.fastq.gz
|--hcc1143_N_small_2.fastq.gz
|--hcc1143_T_clean.bam
|--hcc1143_T_full.bam
|--hcc1143_T_full.bam.bai
|--hcc1143_T_full.cram
|--hcc1143_T_full.cram.crai
|--hcc1143_T_full.unmap.bam
|--hcc1143_T_full.unmap.bam.bai
|--hcc1143_T_full_1.fastq.gz
|--hcc1143_T_full_2.fastq.gz
|--hcc1143_T_small.bam
|--hcc1143_T_small.bam.bai
|--hcc1143_T_small.cram
|--hcc1143_T_small.cram.crai
|--hcc1143_T_small.unmap.bam
|--hcc1143_T_small.unmap.bam.bai
|--hcc1143_T_small_1.fastq.gz
|--hcc1143_T_small_2.fastq.gz
|--mgibio_samtools:1.16.1.sif
|--quay.io_biocontainers_bbmap:39.01--h5c4e2a8_0.sif
```



### Step 7. Download and create test data related to `somaticCNV_hg19`
Download bam files from the URLs listed in **[somaticCNV_hg19.download_links.txt](./download_links/somaticCNV_hg19.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=somaticCNV_hg19 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/somaticCNV_hg19.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**

Create test data by executing the following commands:
```
$ cd /path/to/working/directory/
$ WORKDIR=$(pwd)
$ /bin/sh sauerkraut/scripts/somaticCNV_hg19.create_test_data.sh \
     $WORKDIR/sauerkraut \
     $WORKDIR/cwlenv \
     $WORKDIR/reference_hg19 \
     $WORKDIR/somaticCNV_hg19 \
     SM-74P4M
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**
Repeat the above commands by replacing sample ID (`SM-74P4M`) to other sample ID:
- `SM-74NEG`

The following files will be downloaded or created:
```{somaticCNV_hg19}
/path/to/working/directory/somaticCNV_hg19
|--SM-74NEG.bam
|--SM-74NEG_full.bam
|--SM-74NEG_full.bam.bai
|--SM-74NEG_full.cram
|--SM-74NEG_full.cram.crai
|--SM-74NEG_full.unmap.bam
|--SM-74NEG_full.unmap.bam.bai
|--SM-74NEG_full_1.fastq.gz
|--SM-74NEG_full_2.fastq.gz
|--SM-74NEG_full_interleaved.fastq.gz
|--SM-74NEG_small.bam
|--SM-74NEG_small.bam.bai
|--SM-74NEG_small.cram
|--SM-74NEG_small.cram.crai
|--SM-74NEG_small.unmap.bam
|--SM-74NEG_small.unmap.bam.bai
|--SM-74NEG_small_1.fastq.gz
|--SM-74NEG_small_2.fastq.gz
|--SM-74NEG_small_interleaved.fastq.gz
|--SM-74P4M.bam
|--SM-74P4M_full.bam
|--SM-74P4M_full.bam.bai
|--SM-74P4M_full.cram
|--SM-74P4M_full.cram.crai
|--SM-74P4M_full.unmap.bam
|--SM-74P4M_full.unmap.bam.bai
|--SM-74P4M_full_1.fastq.gz
|--SM-74P4M_full_2.fastq.gz
|--SM-74P4M_full_interleaved.fastq.gz
|--SM-74P4M_small.bam
|--SM-74P4M_small.bam.bai
|--SM-74P4M_small.cram
|--SM-74P4M_small.cram.crai
|--SM-74P4M_small.unmap.bam
|--SM-74P4M_small.unmap.bam.bai
|--SM-74P4M_small_1.fastq.gz
|--SM-74P4M_small_2.fastq.gz
|--SM-74P4M_small_interleaved.fastq.gz
|--broadinstitute_gatk:4.3.0.0.sif
|--mgibio_samtools:1.16.1.sif
|--quay.io_biocontainers_bbmap:39.01--h5c4e2a8_0.sif
```

### Step 8. Download and create test data related to `germlineRNA_hg19`
Download bam files from the URLs listed in **[germlineRNA_hg19.download_links.txt](./download_links/germlineRNA_hg19.download_links.txt)** by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=germlineRNA_hg19 ; mkdir -p $OUTDIR ; for url in `cat sauerkraut/download_links/germlineRNA_hg19.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**

Create test data by executing the following commands:
```
$ cd /path/to/working/directory/
$ WORKDIR=$(pwd)
$ /bin/sh sauerkraut/scripts/germlineRNA_hg19.create_test_data.sh \
     $WORKDIR/sauerkraut \
     $WORKDIR/cwlenv \
     $WORKDIR/reference_hg19 \
     $WORKDIR/germlineRNA_hg19 \
     NA12878
```
Note that this step may take several hours, and therefore, **the use of `nohup` or job scheduler such as `slurm` is recommended.**
There is no need to repeat the above commands. 

The following files will be downloaded or created:
```{germlineRNA_hg19}
/path/to/working/directory/germlineRNA_hg19
|--NA12878.unmapped.bam
|--NA12878_full.bam
|--NA12878_full.bam.bai
|--NA12878_small.bam
|--NA12878_small.bam.bai
|--broadinstitute_gatk:4.3.0.0.sif
|--mgibio_samtools:1.16.1.sif
```

