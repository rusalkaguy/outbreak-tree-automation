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
#   pull latest sequence dumps from ftp
#   identify genomes missing from ftp-seq-dumps
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

# environment setup
# run at cels.anl.gov
ENV_CONFIG="config/cels.anl.gov.yaml" 
configfile: ENV_CONFIG

# taxonomy & tree parameters
# Flu (Orthomyxoviridae) [segmented virus!]
TAX_CONFIG="config/flu-h5n1.yaml"
configfile: TAX_CONFIG

print("Loaded configs:")
print("    # RUNTIME ENV ------")
print("    FTP.server=",f"{config['ftp']['user']}@{config['ftp']['server']}:/{config['ftp']['virus_dir']}")
print("    sif=",config['sif'])
print("    # TAXONOMY ---------")
print("    Family=",config['family'])
print("    p3_genome_filters=",config['p3_genome_filters'])

# download cache
DATA_CACHE_DIR="bv-brc-cache"

METADATA_FILE=f"{DATA_CACHE_DIR}/BVBRC_genome.txt"
METADATA_COLLAPSED_FILE=f"{DATA_CACHE_DIR}/BVBRC_genome.collapsed.txt"
METADATA_ERROR_FILE=f"{DATA_CACHE_DIR}/BVBRC_genome.removed.txt"

BULK_FNA_FILE=f"{DATA_CACHE_DIR}/{config['family']}.fna"
BULK_ID_FNA_FILE=f"{DATA_CACHE_DIR}/{config['family']}.genome_id.fna"

GENOMES_MISSING_HIST=f"{DATA_CACHE_DIR}/BVBRC_genome_id.date_modified.hist.txt"
GENOMES_MISSING_SEQ=f"{DATA_CACHE_DIR}/BVBRC_genome_id.missing_seq.txt"
MISSING_CONTIG_TSV=f"{DATA_CACHE_DIR}/missing_contigs_download.tsv"

localrules: help
rule help:
    run:
        print("Targets:")
        print("    help - list targets")
        print("    metadata - CLI download genome list/metadata")
        print("    bulk - ftp bulk sequence download")
        print("    missing - build list of missing sequences")
        print("    delta - CLI download of missing sequences")
        print("    all")

rule all:
    input:
        metadata=METADATA_FILE,
        bulk_fna=BULK_FNA_FILE,
        missing_list=GENOMES_MISSING_SEQ,
        delta=MISSING_CONTIG_TSV
        
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
localrules: bulk
rule bulk:
    input:
        bulk_fna=BULK_ID_FNA_FILE

rule bulk_fna_download:
    output:
        fna_path=BULK_FNA_FILE
    params:
        user = config["ftp"]["user"],
        server = config["ftp"]["server"],
        virus_dir = config["ftp"]["virus_dir"],
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

rule re_id_bulk_fna:
    output:
        fna_path=BULK_ID_FNA_FILE
    input:
        fna_path=BULK_FNA_FILE,
        script="scripts/fasta_convert_to_genome_id.gawk"
    shell:
        "gawk -f {input.script} {input.fna_path} > {output.fna_path}"
        
rule bulk_fna_faidx:
    output:
        fai_path=f"{BULK_ID_FNA_FILE}.fai",
        err_path=f"{BULK_ID_FNA_FILE}.errors"
    input:
        fna_path=BULK_ID_FNA_FILE
    params:
        sif=config['sif']
    shell:
        "{params.sif} samtools faidx {input.fna_path} 2>{output.err_path} "
        " && "
        "wc -l {output.err_path} "
        
# ----------------------------------------------------------------------
#
# Download ALL metadata (CLI)
#
# https://www.bv-brc.org/view/Taxonomy/11320#view_tab=genomes&filter=eq(subtype,%22H5N1%22)
# expect 264,483
#
# https://www.bv-brc.org/view/Taxonomy/11320#view_tab=genomes&filter=and(eq(subtype,%22H5N1%22),or(eq(collection_year,%222024%22),eq(collection_year,%222025%22),eq(collection_year,%222026%22)),or(eq(genome_status,%22Complete%22),eq(genome_status,%22Partial%22)))&defaultColumns=-cds,h5_clade,segment,sra_accession,collection_date&defaultSort=genome_name,segment
# expect 161,532
# ----------------------------------------------------------------------

rule metadata:
    input:
        raw=METADATA_FILE,
        dedup=METADATA_COLLAPSED_FILE,
        error_report=METADATA_ERROR_FILE

    
rule download_family_metadata:
    output: METADATA_FILE
    params:
        family=config['family'],
        p3_genome_filters=config['p3_genome_filters']
    shell: """
        p3-all-genomes \
            --eq 'family,{params.family}' \
            --eq 'contigs,1' \
            {params.p3_genome_filters} \
            --attr date_modified,genbank_accessions \
            --attr segment \
            --attr genome_length,genome_quality,genome_status \
            --attr genome_name \
        > {output}
           """
    
#
# remove duplicate genomes
#
# If there are multiple genome_id's for the same genbank_accessions,
# then remove all but the oldest,
# and add a column, genome.collapsed_ct, with the original count
#
rule dedup_metadata:
    output:
        metadata=METADATA_COLLAPSED_FILE
    input:
        metadata=METADATA_FILE,
        script="scripts/collapse_bvbrc_genome_by_accession.py"
    shell:
        "{input.script} {input.metadata} {output.metadata}"

#
# dump a report of removed (duplicate) genome_id's
#
rule dedup_metadata_qc_report:
    output:
        report=METADATA_ERROR_FILE
    input:
        orig=METADATA_FILE,
        fixed=METADATA_COLLAPSED_FILE
    shell:
        # headers
        "head -1 {input.orig} > {output}"
        " && " 
        # data
        "join -v1 -t $'\t' "
        "  <(tail -n +2 {input.orig} | sort -k1,1)"
        "  <(cut -f 1 {input.fixed} | tail -n +2 | sort -k1,1) "
        ">> {output.report}"
        " && "
        # stats during snake run
        "wc -l {input.orig} "
        " && "
        "wc -l {input.fixed} {output}"


# ----------------------------------------------------------------------
#
# Download DELTA with CLI
#
# list all genomes with metadata and no sequence.
#
# download those genomes with the CLI
#
# ----------------------------------------------------------------------
rule missing:
    input:
        list=GENOMES_MISSING_SEQ,
        hist=GENOMES_MISSING_HIST

rule list_genomes_missing_sequence:
    output:
        list=GENOMES_MISSING_SEQ
    input:
        metadata=METADATA_COLLAPSED_FILE,
        bulk_fna_fai_path=f"{BULK_ID_FNA_FILE}.fai"
    shell:
        "echo genome_id > {output.list} "
        " && "
        "join -v1 "
        "  <(p3-extract -i {input.metadata} genome_id | tail -n +2 | sort) "
        "  <(cut -f 1 {input.bulk_fna_fai_path} | sort) "
        ">> {output.list}"

rule missing_sequence_date_histogram:
    input:
        download_id_list=GENOMES_MISSING_SEQ,
        metadata=METADATA_FILE
    output:
        hist=GENOMES_MISSING_HIST
    shell:
        # header
        "echo '  count date_modified/day' > {output.hist}"
        " && "
        # get genome_id,date_modified for missing seq genome"
        "grep -w -f <(tail -n +2 {input.download_id_list}) <(cut -f 1-2 {input.metadata} ) "
        # extact just date_modified
        " | cut -f 2 "
        # extract just YYYY-MM-DD
        " | cut -c 1-10 "
        # histogram the dates
        " | sort | uniq -c "
        " >> {output.hist} "
        
rule delta:
    input:
        delta=MISSING_CONTIG_TSV

rule pull_missing_contigs_via_cli:
    output:
        tsv=MISSING_CONTIG_TSV
    input:
        ids=GENOMES_MISSING_SEQ,
        p3_echo_stdin="scripts/p3-echo-stdin.awk"
    threads: 5
    shell: r"""
        set -euo pipefail

        tail -n +2 {input.ids} \
          | parallel -j {threads} -k --pipe -N 500 \
            '{input.p3_echo_stdin} -v title=genome_id | p3-get-genome-contigs --col genome_id --attr sequence' \
          > {output.tsv}
    """    
