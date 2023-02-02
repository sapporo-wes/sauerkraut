#!/usr/bin/env cwl-runner

class: CommandLineTool
id: samtools-index-for-cram
label: samtools-index-for-cram
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: mgibio/samtools:1.16.1
    
baseCommand: [ samtools, index ]

inputs:
  cram:
    type: File
    doc: A CRAM file containing variants
    inputBinding:
      position: 1

outputs:
  - id: crai
    type: File
    outputBinding:
      glob: $(inputs.cram.basename).crai

arguments:
  - position: 2
    prefix: -o 
    valueFrom: $(inputs.cram.basename).crai
