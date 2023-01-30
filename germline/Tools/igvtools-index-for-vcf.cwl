#!/usr/bin/env cwl-runner

class: CommandLineTool
id: igvtools-index-for-vcf
label: igvtools-index-for-vcf
cwlVersion: v1.1

$namespaces:
  edam: http://edamontology.org/

requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.vcf)
        writable: true
  DockerRequirement:
    dockerPull: abralab/igvtools:v2.8.2
    
baseCommand: [ igvtools, index ]

inputs:
  vcf:
    type: File
    doc: A VCF file containing variants
    inputBinding:
      position: 1

outputs:
  - id: idx
    type: File
    outputBinding:
      glob: $(inputs.vcf.basename).idx

