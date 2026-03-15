#!/usr/bin/awk -f
#
# like p3-echo but adds a header to stdin instead of to args
#

# add header
BEGIN {
    print(title)
}

# copy all other lines
{print}
