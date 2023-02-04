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

if [ ! -f ${DATA_DIR}/$sample.unmapped.bam ] ; then
    echo "ERROR: ${DATA_DIR}/$sample.unmapped.bam is required, but not found."
    exit
fi

if [ ! -f ${REF_DIR}/Homo_sapiens_assembly19.fasta ] ; then
    echo "ERROR: ${REF_DIR}/Homo_sapiens_assembly19.fasta is required, but not found."
    exit
fi

if [ ! -f ${DATA_DIR}/${sample}_full.bam ] ; then
    ln -s ${DATA_DIR}/${sample}.unmapped.bam ${DATA_DIR}/${sample}_full.bam
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

if [ ! -f ${DATA_DIR}/${sample}_small.bam ] ; then
    echo ""
    echo "----------------------------------------------"
    echo "Create small-bam"
    cd ${DATA_DIR}
    cwltool --singularity \
	${SRC_DIR}/Tools/gatk4-downsample-for-bam.cwl \
	--bam ${DATA_DIR}/${sample}_full.bam \
	--probability 0.001 \
	--prefix ${sample}_small
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
