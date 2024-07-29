#!/usr/bin/env python3
"""
This script takes a GenBank file with and converts it to a standard fasta file.
The genbank format is structured like so:
LOCUS               defline           length
DEFINITION          defline           length
TITLE               defline
ORIGIN

position    seq beginning
...
position    seq end
//

This is done for each fasta entry. The sequence is also broken up into segments
of 10 bases delimited by spaces. I will set "DEFINITION" as the default field for
capturing the defline while allowing the user to specify a different source.

Generally the order of sequences in the input is maintained in the output. When
defline duplicates are found they are appended to the end of the fasta in no
intentional order and a warning is printed to stderr.

To use this script simply type "python genbank2Fasta.py fileName.gb" or 
"genbank2Fasta.py fileName.gb" if genbank2Fasta.py is in your PATH

Created by  David E. Hufnagel on Feb 7, 2024
last update: April 26, 2024
"""

import sys, os


# Add items into a dictionary. When keys are duplicated append to a list in the value
def SaveIntoDict(key, val, dictx):
    if key in dictx:
        dictx[key].append(val)
    else:
        dictx[key] = [val]


# If an alternative defline source has been requested collect it in alt_def_source
alt_def_source = "DEFINITION"
if len(sys.argv) == 3:
    alt_def_source = sys.argv[2]
    if alt_def_source not in ["LOCUS", "DEFINITION", "TITLE"]:
        print(
            "ERROR: Incorrect defline source or too many arguments provided",
            file=sys.stderr,
        )
        sys.exit()
# no arguments provided
elif len(sys.argv) < 2 or sys.argv[1] in [
    "-h",
    "--help",
    "-help",
    "--version",
]:
    print(
        "Too few arguments provided.\n\n\
Usage instructions for genbank2Fasta.py:\n\
python genbank2Fasta.py <input_file> [defline indicator, default is 'TITLE'] > <output_file>\n",
        file=sys.stderr,
    )
    sys.exit()
# too many arguments provided
elif len(sys.argv) > 3:
    print(
        "Too many arguments provided.\n\n\
Usage instructions for genbank2Fasta.py:\n\
python genbank2Fasta.py <input_file> [defline indicator, default is 'TITLE'] > <output_file>\n",
        file=sys.stderr,
    )
    sys.exit()


# Parse the input file and, for every entry, capture the defline using TITLE or
#  a specified alternate source, also capture the sequence ignoring all other
#  data and formatting; Save the result in a dict of key: defline val: seq.
#  Also save all deflines in a list to maintain the original order.
inp = open(sys.argv[1])
def_list = []
seq_dict = {}
seq = ""
def_str = ""
is_seq = False
is_def = False
for line in inp:
    if is_seq:  # capture sequence data between "ORIGIN" and "//"
        seq_line = "".join(line.strip().split()[1:])  # the sequence data from this line
        seq += seq_line

    if line.startswith("ORIGIN"):  # turn on sequence capture at "ORIGIN"
        is_seq = True
    elif line.startswith("//"):  # reset seq when at the end of the entry
        SaveIntoDict(def_str, seq, seq_dict)
        seq = ""
        def_str = ""
        is_seq = False

    if line.startswith(alt_def_source):  # capture defline line 1
        is_def = True
    # capture subsequent defline lines
    elif not line.startswith(" ") and is_def == True:
        is_def = False
        def_list.append(def_str)

    # capture defline starting with the defline indicator and running through the lines that don't start with fields
    if is_def:
        if line.startswith(" "):
            def_str = def_str + "_" + "_".join(line.split())
        else:
            def_str = "_".join(line.split()[1:])


# Go through the defline list, grab sequences from the dictionary, and output
#  the results in fasta format. Defline duplicates are placed at the end and
#  come with a warning.
to_warn = False  # a boolean for when to warn that deflines are duplicated
def_dups = []
for defline in def_list:
    seqs = seq_dict[defline]
    if (
        len(seqs) > 1
    ):  # When the defline is duplicated add a number to the defline and print a warning
        to_warn = True
        dup_num = 1
        # iterate through each duplicate
        for seq in seqs:
            newlines = f">{defline}_dup{dup_num}\n{seq}\n"
            def_dups.append(newlines)
            dup_num += 1
    else:  # simple non-duplicated cases
        seq = seqs[0]
        newlines = f">{defline}\n{seq}\n"
        print(newlines, end="")


# Print warning for duplicates and write them to the output after converting to a
#  set so as to not duplicate the duplicates in the output.
def_dups = set(def_dups)
if to_warn:
    print("WARNING: deflines are duplicated in the input file", file=sys.stderr)

    for dup in def_dups:
        print(dup, end="")


inp.close()
