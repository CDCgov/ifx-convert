#!/usr/bin/env python3
# Sam Shepard - 2020.01
# See defaults and options: https://arrow.apache.org/docs/python/index.html

import sys
if len(sys.argv) != 3:
	print("Usage:\n\tpython "+sys.argv[0]+"\t<input.txt> <output.parquet>")
	exit(1)


import pyarrow.parquet as pq
from pyarrow import csv

input_table 	= sys.argv[1]
output_parquet	= sys.argv[2]

try:
	table = csv.read_csv(input_table,parse_options = csv.ParseOptions(delimiter='\t'))
except IOError:
	print("Error reading/parsing: " + input_table)
	exit(2)


try:
	pq.write_table(table, output_parquet)
except IOError:
	print("Error writing parquet: " + output_parquet)
	exit(2)
