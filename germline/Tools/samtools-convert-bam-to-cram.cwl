#!/usr/bin/env cwl-runner

class: CommandLineTool
id: samtools-convert-bam-to-cram
label: samtools-convert-bam-to-cram
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: mgibio/samtools:1.16.1
    
baseCommand: [ samtools, view ]

inputs:
  bam:
    type: File
    doc: "A BAM file containing alignments"
    secondaryFiles:
      - .bai
    inputBinding:
      position: 1

  reference:
    type: File
    doc: "A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed by samtools faidx. If an index is not present one will be generated for you, if the reference file is local"
    secondaryFiles: 
      - .fai
    inputBinding: 
      prefix: --reference
      position: 2

outputs:
  - id: cram
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot).cram

arguments:
  - position: 3
    prefix: -o 
    valueFrom: $(inputs.bam.nameroot).cram
  - position: 5
    valueFrom: --cram
