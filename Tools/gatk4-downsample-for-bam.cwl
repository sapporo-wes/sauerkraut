#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gatk4-downsample-for-bam
label: gatk4-downsample-for-bam
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: broadinstitute/gatk:4.3.0.0
  ShellCommandRequirement: {}
    
baseCommand: /usr/bin/java

inputs:
  bam:
    type: File
    doc: A BAM file containing variants
    secondaryFiles:
      - .bai
    inputBinding:
      prefix: --INPUT
      position: 4

  java_options:
    type: string?
    default: -XX:-UseContainerSupport -Xmx50g -Xms50g
    inputBinding:
      position: 1
      shellQuote: false

  probability:
    type: float
    inputBinding:
      position: 5
      prefix: --PROBABILITY

  prefix: 
    type: string
    doc: Output prefix

outputs:
  out_bam:
    type: File
    outputBinding:
      glob: $(inputs.prefix).bam

arguments:
  - position: 2
    prefix: -jar
    valueFrom: /gatk/gatk-package-4.3.0.0-local.jar
  - position: 3
    valueFrom: DownsampleSam
  - position: 6
    prefix: -O
    valueFrom: $(inputs.prefix).bam
