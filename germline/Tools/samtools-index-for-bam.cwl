#!/usr/bin/env cwl-runner

class: CommandLineTool
id: samtools-index-for-bam
label: samtools-index-for-bam
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: mgibio/samtools:1.16.1
    
baseCommand: [ samtools, index ]

inputs:
  bam:
    type: File
    doc: A BAM file containing variants
    inputBinding:
      position: 1

outputs:
  - id: bai
    type: File
    outputBinding:
      glob: $(inputs.bam.basename).bai

arguments:
  - position: 2
    prefix: -o 
    valueFrom: $(inputs.bam.basename).bai
