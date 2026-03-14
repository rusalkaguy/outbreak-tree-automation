#!/usr/bin/env snakemake --cores 1 help
#
# Framework for automated outbreak tree building
#
# Run as:
#
#   snakemake --cores 1 --printshellcmd download
#
# Full Protocol:
#   pull metatdata list of genome
#   pull latest sequence dumps from sftp
#   identify genomes missing from sftp-seq-dumps
#   pull missing genomes via CLI
#   Build fasta
#   Align fasta -> MSA
#   Build tree
#

# ----------------------------------------------------------------------
#
# Config
#
# ----------------------------------------------------------------------

# BV-BRC ftp setup
configfile: "config/bv-brc.org.yaml"

# FLU outbreak config
configfile: "config/flu-h5n1.yaml"

print("Loaded configs:")
print("\tSFTP.server=",config['sftp']['server'])
print("\tFamily=",config['family'])

# download cache
DATA_CACHE_DIR="bv-brc-cache"
METADATA_FILE=f"{DATA_CACHE_DIR}/BVBRC_genome.txt"
BULK_FNA_FILE=f"{DATA_CACHE_DIR}/{config['family']}.fna"

localrules: help
rule help:
    run:
        print("Targets:")
        print("    help")
        print("    metadata")
        print("    download")

rule all:
    input:
        metadata=METADATA_FILE,
        bulk_fna=BULK_FNA_FILE
        
# ----------------------------------------------------------------------
#
# try storage.ftp to get remote mtime integrated in the DAG
#
# FAILS: storage.ftp has bad ftps support
# ----------------------------------------------------------------------

# comment out, otherwise requires snakemake-storage-plugin-ftp to be installed
#include: "rules/test_storage_ftps.smk"

# ----------------------------------------------------------------------
#
# CLEAN
#
# ----------------------------------------------------------------------
rule clean:
    shell: """
       rm -rf {DATA_CACHE_DIR}
    """
    
# ----------------------------------------------------------------------
# Bulk FNA download
#
# @ANL: ~ 0m8sec for 3.3G Orthomyxoviridae.fna
# @UAB: ~ 1m     for 3.3G Orthomyxoviridae.fna
# ----------------------------------------------------------------------
localrules: download
rule download:
    input:
        bulk_fna=BULK_FNA_FILE

rule bulk_fna_download:
    output:
        fna_path=BULK_FNA_FILE
    params:
        user = config["sftp"]["user"],
        server = config["sftp"]["server"],
        virus_dir = config["sftp"]["virus_dir"],
        fna_name = f"{config['family']}.fna"
    # force Snakemake to always try and pull
    #always_run: True
    shell:
        "curl --ssl-reqd "
        # complain on failure
        "  --fail --silent --show-error "
        # conserve server mod date/time
        "  -R "
        # mod date/time check between server and this file
        "  -z {output.fna_path} "
        # force output filename and path
        "  -o {output.fna_path} "
        # connection params
        " --user '{params.user}' "
        # server file
        " '{params.server}/{params.virus_dir}/{params.fna_name}' "

rule bulk_fna_faidx:
    output:
        fai_path=f"{BULK_FNA_FILE}.fai",
        err_path=f"{BULK_FNA_FILE}.errors"
    input:
        fna_path=BULK_FNA_FILE
    params:
        sif=config['sif']
    shell:
        "singularity exec {params.sif} samtools faidx {input.fna_path} 2>{output.err_path} "
        " && "
        "wc -l {output.err_path} "
        
# ----------------------------------------------------------------------
#
# Download ALL metadata (CLI)
#
# https://www.bv-brc.org/view/Taxonomy/11320#view_tab=genomes&filter=eq(subtype,%22H5N1%22)
#
# expect 264,483
# ----------------------------------------------------------------------
rule metadata:
    input: METADATA_FILE
    
rule download_family_metadata:
    output: METADATA_FILE
    params:
        family=config['family'],
        subtype=config['subtype']
    shell: """
        p3-all-genomes \
            --eq 'family,{params.family}' \
            --eq 'subtype,{params.subtype}' \
            --eq 'contigs,1' \
            --attr genbank_accessions,segment \
        > {output}
           """
    
