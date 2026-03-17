# BV-BRC outbreak-tree-automation

## Proposed Protocol:
  *   pull metatdata list of genome [done]
  *   pull latest sequence dumps from sftp [done]
  *   identify genomes missing from sftp-seq-dumps [done]
  *   parallel pull missing genomes via CLI [done]
  *   Build per-segment fastas based on metadat:genome.segment 
  *   Align fasta -> MSA
  *   Update MSA.fasta seq headers from metadata for tree builder
  *   Build tree(s)

## Quick Start

To Run on a CELS machine at anl.gov:
```
# set up environment
source /vol/patric3/cli/user-env.sh

# clean caches and outputs
snakemake --cores 1 clean

# downlaod and build
#   ~30 minute runtime for config/flu-h5n1.yaml
#   ./config/flu-h5n1.yaml
#       family=Orthomyxoviridae
#       subtype=H5N1
#       collection_year=[2024,2025,2026]
#       genome_status=[Complete,Partial]
#       contigs=1
snakemake --rerun-incomplete --cores 10 --printshellcmd all
```
## Configuration

Configuration files are in [./config](./config]
  * cels.anl.gov.yaml
    * Purpose: defines the execution environment.
    * Sections:
       * ftp.*: params for bulk FTP downloads
       * sif: path to Singularity image for standard tools; comment out to run on MacOS
  * flu-h5n1.yaml
    * Purpose: defines the taxononmy and metadata fileters for the dataset
    * Multiplicity: There should be one of these files for each outbreak; it can be overridden on the snakemake CLI
    * Sections:
        * family: controls bulk ftp download and p3-genome metadata query
        * p3_genome_filters: additional genome selection filters
	* TBA: controls for fasta sequence names, tree decorations, MSA parameters, etc.

## Snakemake Targets

```
# go ask snakemake
snakemake --cores 1 help
```
gives
```
Loaded configs:
    # Snakefile settings
    DATA_CACHE_DIR= bv-brc-cache
    # RUNTIME ENV ------
    FTP.server= anonymous:guest@ftp://ftp.bv-brc.org:/viruses
    sif= singularity exec /vol/patric3/production/containers/ubuntu-dev-160-build12-2.sif
    # TAXONOMY ---------
    Family= Orthomyxoviridae
    p3_genome_filters= --eq 'subtype,H5N1' --in 'collection_year,2024,2025,2026' --in 'genome_status,Complete,Partial'
Targets:
    help - list targets
    clean - delete all downloads and outputs
    metadata - CLI download genome list/metadata
    bulk - ftp bulk sequence download
    missing - build list of missing sequences
    delta - CLI download of missing sequences
    all
```
   
## Limitations

* only valid for clean runs, for 3 reasons:
 * metadata timestamp
   * there is nothing in place to check if CLI metadata query returns exactly what it returned last time. The query is always run, if the results file doesn't exist, and not run if the results file exists. There is no timestamp to check on the server
     * we could implement logic to only update the local file if the new download actually had new data, but the the text file processing is fast, so this code complexity is not worth it. 
 * bulk timestamp
   * Snakemake can NOT read the modification timestamp from an ftp*s* server. Thus, when it's building it's execution plan, it can't make a download decision about the bulk sequence dump. While the curl command commend will make a decision about whether or not to skip the download, and will preserve the server timestamp, by the time the curl command runs, the Snakemake excution plan is already fixed
     * this is a limitation of Snakemake's ```storage.ftp``` extension
     * workaround: have a driver script run snakemake twice:
       * ```snakemake metadata``` - download latest metadata (always downloads)
       * ```snakemake bulk``` - download bulk (curl - only downloads if updated)
       * ```snakemake all``` - rebuild based on bulk's timestamp
  * ```delta``` CLI sequence download
    * currently this is not able to re-use the CLI download sequence. It always re-downloads all the sequences missing from bulk. 
    * this could be fixed by having the rules that makes the "sequence_missing" list take into account andy CLI downloaded sequence already in the cache
    * this would be a great optimization and would make the pipeline run MUCH faster.


    
       
   