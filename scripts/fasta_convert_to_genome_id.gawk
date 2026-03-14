#!/usr/bin/env gawk -f
#
# rename_fasta_ids.awk
# For each FASTA header, extract the trailing digits.digits]$
# and replace the full header with just >digits.digits
#
# Example:
#  >4BBL_Y   Cryo-electron microscopy reconstruction of the helical part of influenza A virus ribonucleoprotein isolated from virions.   [Influenza A virus | 11320.6579]
# becomes
#  >11320.6579

#
# detect and re-write DEF_LINE
#
/^>/ {
    if (match($0, /[0-9]+\.[0-9]+\]$/)) {
        id = substr($0, RSTART, RLENGTH-1)
        print ">" id " " substr($0,2)
    } else {
        printf "ERROR: no trailing /digits.digits]$/ (aka genome_id) pattern in header:\n%s\n", $0 > "/dev/stderr"
        exit 1
    }
    next
}

#
# pass through sequence unmodified
#
{ print }
