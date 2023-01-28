#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gatk4-RevertSam
label: gatk4-RevertSam
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
    default: -XX:-UseContainerSupport -Xmx5g -Xms5g
    inputBinding:
      position: 1
      shellQuote: false

outputs:
  ubam:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot).unmap.bam

arguments:
  - position: 2
    prefix: -jar
    valueFrom: /gatk/gatk-package-4.3.0.0-local.jar
  - position: 3
    valueFrom: RevertSam
  - position: 5
    prefix: --OUTPUT_BY_READGROUP
    valueFrom: "false"
  - position: 6
    prefix: --VALIDATION_STRINGENCY
    valueFrom: "LENIENT"
  - position: 7
    prefix: --ATTRIBUTE_TO_CLEAR
    valueFrom: "FT"
  - position: 8
    prefix: --ATTRIBUTE_TO_CLEAR
    valueFrom: "CO"
  - position: 9
    prefix: --SORT_ORDER
    valueFrom: "coordinate"
  - position: 10
    prefix: -O
    valueFrom: $(inputs.bam.nameroot).unmap.bam
