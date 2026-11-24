#!/usr/bin/env cwl-runner

class: CommandLineTool
id: bcftools-index-for-vcfgz
label: bcftools-index-for-vcfgz
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bcftools:1.16--hfe4b78e_1
    
baseCommand: [ bcftools, index ]

inputs:
  vcf:
    type: File
    doc: A VCF file containing variants
    inputBinding:
      position: 1

outputs:
  - id: tbi
    type: File
    outputBinding:
      glob: $(inputs.vcf.basename).tbi

arguments:
  - position: 2
    prefix: -o 
    valueFrom: $(inputs.vcf.basename).tbi
