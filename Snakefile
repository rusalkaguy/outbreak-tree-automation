#!/usr/bin/env snakemake
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

# download cache
SEQ_CACHE="sequence_cache"

# BV-BRC ftp setup
configfile: "config/ftp.bv-brc.org.yaml"

# FLU outbreak config
configfile: "config/flu-h5n1.yaml"

print("Loaded configs:")
print("\tSFTP.server=",config['sftp']['server'])
print("\tFamily=",config['family'])


localrules: help
rule help:
    run:
        print("Targets:")
        print("    test")
        print("    download")

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
       rm -rf {SEQ_CACHE}
    """
    
# ----------------------------------------------------------------------
# Bulk FNA download
#
# @UAB: ~ 1 minu for 3.3G Orthomyxoviridae.fna
# ----------------------------------------------------------------------
localrules: download
rule download:
    input:
        bulk_fna=f"{SEQ_CACHE}/{config['family']}.fna"

rule bulk_fna_download:
    output:
        fna_path=f"{SEQ_CACHE}/{config['family']}.fna"
    params:
        user = config["sftp"]["user"],
        server = config["sftp"]["server"],
        virus_dir = config["sftp"]["virus_dir"],
        fna_name = f"{config['family']}.fna"
    shell:
        "curl --ssl-reqd "
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
