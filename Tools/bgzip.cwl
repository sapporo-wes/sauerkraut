#!/usr/bin/env cwl-runner

class: CommandLineTool
id: bgzip
label: bgzip
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.vcf)
        writable: true
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bcftools:1.16--hfe4b78e_1

baseCommand: bgzip

inputs:
  vcf:
    type: File
    inputBinding:
      position: 1

outputs:
  vcf_gz:
    type: File
    outputBinding:
      glob: $(inputs.vcf.basename).gz

