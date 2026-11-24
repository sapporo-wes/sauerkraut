#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gatk4-ExtractSequences
label: gatk4-ExtractSequences
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: broadinstitute/gatk:4.3.0.0
  ShellCommandRequirement: {}
    
baseCommand: /usr/bin/java

inputs:
  fasta:
    type: File
    doc: A reference FASTA file
    secondaryFiles: 
      - .fai
    inputBinding:
      position: 4
      prefix: --REFERENCE_SEQUENCE

  interval_list:
    type: File
    doc: Interval list describing intervals to be extracted from the reference sequence
    inputBinding: 
      position: 5
      prefix: --INTERVAL_LIST

  interval_name:
    type: string
    doc: interval name

  java_options:
    type: string?
    default: -XX:-UseContainerSupport -Xmx32g -Xms32g
    inputBinding:
      position: 1
      shellQuote: false

outputs:
  out_fasta:
    type: File
    outputBinding:
      glob: $(inputs.fasta.nameroot).$(inputs.interval_name).fasta

arguments:
  - position: 2
    prefix: -jar
    valueFrom: /gatk/gatk-package-4.3.0.0-local.jar
  - position: 3
    valueFrom: ExtractSequences
  - position: 6
    prefix: -O
    valueFrom: $(inputs.fasta.nameroot).$(inputs.interval_name).fasta
