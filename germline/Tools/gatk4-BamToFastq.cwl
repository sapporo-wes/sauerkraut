#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gatk4-BamToFastq
label: gatk4-BamToFastq
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

outputs:
  fastq1:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot)_1.fastq.gz
  fastq2:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot)_2.fastq.gz

arguments:
  - position: 2
    prefix: -jar
    valueFrom: /gatk/gatk-package-4.3.0.0-local.jar
  - position: 3
    valueFrom: SamToFastq
  - position: 5
    prefix: --INCLUDE_NON_PF_READS
    valueFrom: "true"
  - position: 6
    prefix: --VALIDATION_STRINGENCY
    valueFrom: "SILENT"
  - position: 6
    prefix: --FASTQ
    valueFrom: $(inputs.bam.nameroot)_1.fastq.gz
  - position: 7
    prefix: --SECOND_END_FASTQ
    valueFrom: $(inputs.bam.nameroot)_2.fastq.gz
