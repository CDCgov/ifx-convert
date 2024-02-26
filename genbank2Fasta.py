#!/usr/bin/env python3
"""
This script takes a GenBank file with and converts it to a standard fasta file.
The genbank (in this case .fas) format is structured like so:
LOCUS               defline           length
DEFINITION          defline           length
TITLE               defline
ORIGIN

position    seq beginning
...
position    seq end
//

This is done for each fasta entry. The sequence is also broken up into segments
of 10 bases delimited by spaces. I will set "TITLE" as the default field for
capturing the defline while allowing the user to specify a different source.

Generally the order of sequences in the input is maintained in the output. When
defline duplicates are found they are appended to the end of the fasta in no
intentional order and a warning is printed to stderr.

To use this script simply type "python genbank2Fasta.py fileName.fas" or 
"genbank2Fasta.py fileName.fas" if genbank2Fasta.py is in your PATH

Created by  David E. Hufnagel on Feb 7, 2024
last update: Feb 22, 2024
"""

import sys, os


# Add items into a dictionary. When keys are duplicated append to a list in the value
def SaveIntoDict(key, val, dictx):
    if key in dictx:
        dictx[key].append(val)
    else:
        dictx[key] = [val]


# If an alternative defline source has been requested collect it in alt_def_source
alt_def_source = "TITLE"
if len(sys.argv) > 2:
    alt_def_source = sys.argv[2]
    if alt_def_source not in ["LOCUS", "DEFINITION", "TITLE"]:
        print("ERROR: Incorrect defline source or too many arguments provided", file=sys.stderr)
        sys.exit()


# Parse the input file and, for every entry, capture the defline using TITLE or
#  a specified alternate source, also capture the sequence ignoring all other
#  data and formatting; Save the result in a dict of key: defline val: seq.
#  Also save all deflines in a list to maintain the original order.
inp = open(sys.argv[1])
def_list = []
seq_dict = {}
seq = ""
for line in inp:

    if line.startswith(alt_def_source):
        defline = line.split()[1]
        def_list.append(defline)
    elif line.startswith(" "):
        seq_line = "".join(line.strip().split()[1:])  # the sequence data from this line
        seq += seq_line
    elif line.startswith("//"):  # reset seq when at the end of the entry
        SaveIntoDict(defline, seq, seq_dict)
        seq = ""


# Go through the defline list, grab sequences from the dictionary, and output
#  the results in fasta format. Defline duplicates are placed at the end and
#  come with a warning.
to_warn = False  # a boolean for when to warn that deflines are duplicated
def_dups = []
for defline in def_list:
    seqs = seq_dict[defline]
    if (len(seqs) > 1):  # When the defline is duplicated add a number to the defline and print a warning
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
