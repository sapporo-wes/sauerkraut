#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gatk4-CreateSequenceDictionary
label: gatk4-CreateSequenceDictionary
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
    inputBinding:
      position: 4
      prefix: --REFERENCE

  java_options:
    type: string?
    default: -XX:-UseContainerSupport -Xmx64g -Xms64g
    inputBinding:
      position: 1
      shellQuote: false

outputs:
  dict:
    type: File
    outputBinding:
      glob: $(inputs.fasta.nameroot).dict

arguments:
  - position: 2
    prefix: -jar
    valueFrom: /gatk/gatk-package-4.3.0.0-local.jar
  - position: 3
    valueFrom: CreateSequenceDictionary
  - position: 5
    prefix: -O
    valueFrom: $(inputs.fasta.nameroot).dict
