#!/usr/bin/env cwl-runner

class: CommandLineTool
id: samtools-extract-by-region-from-cram-to-cram
label: samtools-extract-by-region-from-cram-to-cram
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
      position: 3

  region: 
    type: string
    doc: "Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and all position coordinates are 1-based"
    inputBinding:
      position: 2

  region_name:
    type: string
    doc: "Region name"

outputs:
  - id: cram
    type: File
    outputBinding:
      glob: $(inputs.cram.nameroot).$(inputs.region_name).cram

arguments:
  - position: 4
    prefix: -o 
    valueFrom: $(inputs.cram.nameroot).$(inputs.region_name).cram
  - position: 5
    valueFrom: --cram
