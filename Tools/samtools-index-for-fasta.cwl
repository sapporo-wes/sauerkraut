#!/usr/bin/env cwl-runner

class: CommandLineTool
id: samtools-index-for-fasta
label: samtools-index-for-fasta
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: mgibio/samtools:1.16.1
    
baseCommand: [ samtools, faidx ]

inputs:
  fasta:
    type: File
    doc: A FASTA file
    inputBinding:
      position: 1

outputs:
  - id: fai
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename).fai

arguments:
  - position: 2
    prefix: -o 
    valueFrom: $(inputs.fasta.basename).fai
