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

