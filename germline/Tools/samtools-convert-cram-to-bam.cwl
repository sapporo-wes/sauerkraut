#!/usr/bin/env cwl-runner

class: CommandLineTool
id: samtools-convert-cram-to-bam
label: samtools-convert-cram-to-bam
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: mgibio/samtools:1.16.1
    
baseCommand: [ samtools, view ]

inputs:
  cram:
    type: File
    doc: "A CRAM file containing variants"
    secondaryFiles:
      - .crai
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
  - id: bam
    type: File
    outputBinding:
      glob: $(inputs.cram.nameroot).bam

arguments:
  - position: 3
    prefix: -o 
    valueFrom: $(inputs.cram.nameroot).bam
  - position: 5
    valueFrom: --bam
