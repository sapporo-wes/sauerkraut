#!/bin/sh

# Usage
usage(){
    echo "Usage: "
    echo "CMD <SRC_DIR> <ENV_DIR> <OUT_DIR>"
}

if [ $# -ne 3 ] ; then
    usage
    exit
fi

SRC_DIR=$1
ENV_DIR=$2
REF_DIR=$3

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.fasta ] ; then
    echo "ERROR: ${REF_DIR}/Homo_sapiens_assembly38.fasta is required, but not found."
    exit
fi

. ${ENV_DIR}/bin/activate
if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.dict ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-DICT"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-CreateSequenceDictionary.cwl \
	--fasta ${REF_DIR}/Homo_sapiens_assembly38.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.interval_list ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-INTERVAL_LIST"
    cat ${REF_DIR}/Homo_sapiens_assembly38.dict > ${REF_DIR}/Homo_sapiens_assembly38.small.interval_list 
    echo -e "chr1\t1000000\t1200000\t+\tchr1:1000000-1200000" >> ${REF_DIR}/Homo_sapiens_assembly38.small.interval_list
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-FASTA"
    cd ${REF_DIR}
    cwltool --singularity \
	--leave-tmpdir \
	${SRC_DIR}/Tools/gatk4-ExtractSequences.cwl \
	--fasta ${REF_DIR}/Homo_sapiens_assembly38.fasta \
	--interval_list ${REF_DIR}/Homo_sapiens_assembly38.small.interval_list  \
	--interval_name small
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.fai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-FAI"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-fasta.cwl \
	--fasta ${REF_DIR}/Homo_sapiens_assembly38.small.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.dict ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-DICT"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-CreateSequenceDictionary.cwl \
	--fasta ${REF_DIR}/Homo_sapiens_assembly38.small.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.amb ] || 
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.ann ] || 
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.bwt ] || 
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.pac ] ||  
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.sa ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-BWA-INDEX"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bwa-index.cwl \
	--fasta ${REF_DIR}/Homo_sapiens_assembly38.small.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.fasta.amb ] || 
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.fasta.ann ] || 
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.fasta.bwt ] || 
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.fasta.pac ] ||  
    [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.fasta.sa ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-BWA-INDEX"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bwa-index.cwl \
	--fasta ${REF_DIR}/Homo_sapiens_assembly38.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.small.fasta.alt ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-ALT"
    cd ${REF_DIR}
    cp Homo_sapiens_assembly38.fasta.alt Homo_sapiens_assembly38.small.fasta.alt 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-DBSNP-VCFGZ"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bgzip.cwl \
	--vcf ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-DBSNP-VCFGZ-TBI"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bcftools-index-for-vcfgz.cwl \
	--vcf ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.small.vcf ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-DBSNP"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bcftools-extract-by-region-from-vcfgz-to-vcf.cwl \
	--vcf ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
	--region chr1:1000000-1200000 \
	--region_name small
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.small.vcf.idx ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-DBSNP-IDX"
    cd ${REF_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/igvtools-index-for-vcf.cwl \
	--vcf ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.small.vcf 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

for file in Mills_and_1000G_gold_standard.indels.hg38 \
    Homo_sapiens_assembly38.known_indels \
    hapmap_3.3.hg38 \
    1000G_omni2.5.hg38 \
    1000G_phase1.snps.high_confidence.hg38 \
    ; do 
    if [ ! -f ${REF_DIR}/$file.vcf.small.vcf.gz ] ; then
	echo ""
	echo "----------------------------------------------"
	echo "Create $file.vcf.small.vcf.gz"
	cd ${REF_DIR}
	cwltool --singularity \
	    ${SRC_DIR}/Tools/bcftools-extract-by-region-from-vcfgz-to-vcfgz.cwl \
	    --vcf ${REF_DIR}/$file.vcf.gz \
	    --region chr1:1000000-1200000 \
	    --region_name small
	echo "...done."
	echo "----------------------------------------------"
	echo ""
    fi

    if [ ! -f ${REF_DIR}/$file.vcf.small.vcf.gz.tbi ] ; then
	echo ""
	echo "----------------------------------------------"
	echo "Create $file.vcf.small.vcf.gz.tbi"
	cd ${REF_DIR}
	cwltool --singularity \
	    ${SRC_DIR}/Tools/bcftools-index-for-vcfgz.cwl \
	    --vcf ${REF_DIR}/$file.vcf.small.vcf.gz 
	echo "...done."
	echo "----------------------------------------------"
	echo ""
    fi
done

