# BV-BRC outbreak-tree-automation

Proposed Protocol:
  *   pull metatdata list of genome
  *   pull latest sequence dumps from sftp [done]
  *   identify genomes missing from sftp-seq-dumps
  *   pull missing genomes via CLI
  *   Build fasta
  *   Align fasta -> MSA
  *   Build tree

To Run:
```
snakemake --cores 1 --printshellcmd download
```
