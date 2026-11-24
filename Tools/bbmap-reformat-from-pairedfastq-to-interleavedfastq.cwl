#!/usr/bin/env cwl-runner

class: CommandLineTool
id: bbmap-reformat-from-pairedfastq-to-interleavedfastq
label: bbmap-reformat-from-pairedfastq-to-interleavedfastq
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0
  ShellCommandRequirement: {}

baseCommand: reformat.sh

inputs:
  fastq1:
    type: File
    doc: Paired FASTQ file 1

  fastq2:
    type: File
    doc: Paired FASTQ file 2

  outprefix: 
    type: string
    doc: Output prefix

outputs:
  out_fastq:
    type: File
    outputBinding:
      glob: $(inputs.outprefix).fastq.gz

arguments:
  - position: 1
    valueFrom: "in1=$(inputs.fastq1.path)"
  - position: 2
    valueFrom: "in2=$(inputs.fastq2.path)"
  - position: 3
    valueFrom: "out=$(inputs.outprefix).fastq.gz"
