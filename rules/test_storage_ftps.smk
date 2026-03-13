# ----------------------------------------------------------------------
#
# try storage.ftp to get remote mtime integrated in the DAG
#
# FAILS: bad ftps: support
# ----------------------------------------------------------------------

storage:
    provider="ftp",
    username="anonymous",
    password="guest"

localrules: test_storage_ftps
rule test_storage_ftps:
    input:
        storage.ftp("ftps://ftp.bv-brc.org/viruses/Orthomyxoviridae.fna")
    output:
        "sequence_cache/Orthomyxoviridae.via_storage.fna"
    shell:
        "cp {input} {output}"

