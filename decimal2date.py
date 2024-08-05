#!/usr/bin/env python
# Sam Shepard - 4.26.2011 - decimalDate.py
# Converts fasta file with header format: >date_virus
# To the decimal date format used by beast.
# Ported from the decimalDate.php by Justin Bahl & Yongmei Liu.

import sys
import re

if len(sys.argv) < 3:
    print(
        "Usage:\n\tpython "
        + sys.argv[0]
        + "\t<output_file> <input.fasta> [additional_files.fasta...]"
    )
    exit(2)

try:
    outfile = open(sys.argv[1], "w")
except IOError:
    print("Error with filename.")
    exit(2)

days1 = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
days2 = (0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
dds = 0
ddfrac = 0

for a in range(2, len(sys.argv)):
    try:
        f = open(sys.argv[a], "r")
    except IOError:
        print("Error with filename.")
        exit(2)
    for line in f:
        if line.startswith(">"):
            pieces = ()
            headerFields = line.split("_")
            headerFields[0] = headerFields[0].lstrip(">")
            headerFields[0] = headerFields[0].rstrip("|")

            if headerFields[0].find("/") != -1:
                pieces = map(int, headerFields[0].split("/"))
            elif headerFields[0].find("-") != -1:
                pieces = map(int, headerFields[0].split("-"))
            else:
                pieces = (int(headerFields[0]),)

            if len(pieces) == 3:
                (yy, mm, dd) = pieces
            elif len(pieces) == 2:
                (yy, mm) = pieces
                dd = 15
            else:
                yy = pieces[0]
                (mm, dd) = (6, 15)

            if yy % 4 == 0:
                dds = 0
                for i in range(1, mm):
                    dds = dds + days2[i]
                dds = dds + dd
                ddfrac = float(dds) / 366.0
            else:
                dds = 0
                for i in range(1, mm):
                    dds = dds + days1[i]
                dds = dds + dd
                ddfrac = float(dds) / 365.0
            if dd == 31 and mm == 12:
                ddfrac = 0.999
            headerFields[0] = (">%s" % yy) + ("%.3f" % ddfrac).lstrip("0")

            line = "_".join(headerFields)
        outfile.write(line)
    f.close()
outfile.close()
