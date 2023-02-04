#!/bin/sh

# Usage
usage(){
    echo "Usage: "
    echo "CMD <SRC_DIR> <ENV_DIR> <REF_DIR> <OUT_DIR> <SAMPLE_ID>"
}

if [ $# -ne 5 ] ; then
    usage
    exit
fi

SRC_DIR=$1
ENV_DIR=$2
REF_DIR=$3
DATA_DIR=$4
sample=$5

if [ ! -f ${DATA_DIR}/${sample}.bam ] ; then
    echo "ERROR: ${DATA_DIR}/${sample}_clean.bam is required, but not found."
    exit
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly19.fasta ] ; then
    echo "ERROR: ${REF_DIR}/Homo_sapiens_assembly19.fasta is required, but not found."
    exit
fi

if [ ! -f ${DATA_DIR}/${sample}_full.bam ] ; then
    ln -s ${DATA_DIR}/${sample}.bam ${DATA_DIR}/${sample}_full.bam
fi

. ${ENV_DIR}/bin/activate
if [ ! -f ${DATA_DIR}/${sample}_full.bam.bai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create index for full-bam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-bam.cwl \
	--bam ${DATA_DIR}/${sample}_full.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_full.cram ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-cram"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-convert-bam-to-cram.cwl \
	--bam ${DATA_DIR}/${sample}_full.bam \
	--reference ${REF_DIR}/Homo_sapiens_assembly19.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_full.cram.crai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create index for full-cram"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-cram.cwl \
	--cram ${DATA_DIR}/${sample}_full.cram
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small.cram ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-cram"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-extract-by-region-from-cram-to-cram.cwl \
	--cram ${DATA_DIR}/${sample}_full.cram \
	--reference ${REF_DIR}/Homo_sapiens_assembly19.fasta \
	--region chr1:1000001-1200000 \
	--prefix ${sample}_small
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small.cram.crai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create index for small-cram"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-cram.cwl \
	--cram ${DATA_DIR}/${sample}_small.cram
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small.bam ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-bam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-convert-cram-to-bam.cwl \
	--cram ${DATA_DIR}/${sample}_small.cram \
	--reference ${REF_DIR}/Homo_sapiens_assembly19.fasta 
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small.bam.bai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create index for small-bam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-bam.cwl \
	--bam ${DATA_DIR}/${sample}_small.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small.unmap.bam ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-ubam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-RevertSam \
	--bam ${DATA_DIR}/${sample}_small.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small.unmap.bam.bai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create index for small-ubam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-bam.cwl \
	--bam ${DATA_DIR}/${sample}_small.unmap.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_full.unmap.bam ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-ubam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-RevertSam \
	--bam ${DATA_DIR}/${sample}_full.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_full.unmap.bam.bai ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create index for small-ubam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/samtools-index-for-bam.cwl \
	--bam ${DATA_DIR}/${sample}_full.unmap.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small_1.fastq.gz ] || [ ! -f ${DATA_DIR}/${sample}_small_2.fastq.gz ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-fastq"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-BamToFastq.cwl \
	--bam ${DATA_DIR}/${sample}_small.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_full_1.fastq.gz ] || [ ! -f ${DATA_DIR}/${sample}_full_2.fastq.gz ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-fastq"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-BamToFastq.cwl \
	--bam ${DATA_DIR}/${sample}_full.bam
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_small_interleaved.fastq.gz ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-interleaved-fastq"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bbmap-reformat-from-pairedfastq-to-interleavedfastq.cwl \
	--fastq1 ${DATA_DIR}/${sample}_small_1.fastq.gz \
	--fastq2 ${DATA_DIR}/${sample}_small_2.fastq.gz \
	--outprefix ${sample}_small_interleaved
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

if [ ! -f ${DATA_DIR}/${sample}_full_interleaved.fastq.gz ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create full-interleaved-fastq"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/bbmap-reformat-from-pairedfastq-to-interleavedfastq.cwl \
	--fastq1 ${DATA_DIR}/${sample}_full_1.fastq.gz \
	--fastq2 ${DATA_DIR}/${sample}_full_2.fastq.gz \
	--outprefix ${sample}_full_interleaved
    echo "...done."
    echo "----------------------------------------------"
    echo ""
fi

