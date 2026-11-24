#!/usr/bin/env cwl-runner

class: CommandLineTool
id: bwa-index
label: bwa-index
cwlVersion: v1.1

$namespaces:
  edam: 'http://edamontology.org/'

requirements:
  DockerRequirement:
    dockerPull: biocontainers/bwa:v0.7.17_cv1
  ShellCommandRequirement: {}

baseCommand: [ bwa, index ]

inputs:
  fasta:
    type: File
    doc: A reference FASTA file
    inputBinding:
      position: 2

outputs:
  amb:
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename).amb
  ann:
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename).ann
  bwt:
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename).bwt
  pac:
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename).pac
  sa:
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename).sa

arguments:
  - position: 1
    prefix: -p
    valueFrom: $(inputs.fasta.basename)
