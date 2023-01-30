#!/usr/bin/env cwl-runner

class: CommandLineTool
id: bcftools-extract-by-region-from-vcfgz-to-vcf
label: bcftools-extract-by-region-from-vcfgz-to-vcf
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bcftools:1.16--hfe4b78e_1
    
baseCommand: [ bcftools, view ]

inputs:
  vcf:
    type: File
    doc: A VCF file containing variants
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 1

  region: 
    type: string
    doc: "Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and all position coordinates are 1-based"
    inputBinding:
      position: 2

  region_name:
    type: string
    doc: "Region name"

outputs:
  - id: out_vcf
    type: File
    outputBinding:
      glob: $(inputs.vcf.nameroot).$(inputs.region_name).vcf

arguments:
  - position: 3
    prefix: -o 
    valueFrom: $(inputs.vcf.nameroot).$(inputs.region_name).vcf
  - position: 4
    prefix: -O
    valueFrom: v
