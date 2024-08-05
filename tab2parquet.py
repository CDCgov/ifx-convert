#!/usr/bin/env python3
# Sam Shepard - 2020.01
# See defaults and options: https://arrow.apache.org/docs/python/index.html

import sys

if len(sys.argv) != 3:
    print("Usage:\n\tpython " + sys.argv[0] + "\t<input.txt> <output.parquet>")
    exit(1)


import pyarrow.parquet as pq
from pyarrow import csv

input_table = sys.argv[1]
output_parquet = sys.argv[2]

# Attribution: https://github.com/apache/arrow/blob/f7ef65e5fc367f1f5649dfcea0754e413fcca394/cpp/src/arrow/csv/options.cc#L28-L30
# Added Impala-style null strings
# Have to remove "NA" if strings can be NULL
null_values = [
    "",
    "#N/A",
    "#N/AN/A",
    "-1.#IND",
    "-1.#QNAN",
    "-NaN",
    "-nan",
    "1.#IND",
    "1.#QNAN",
    "N/A",
    "NULL",
    "NaN",
    "n/a",
    "nan",
    "null",
    "\\N",
]

try:
    table = csv.read_csv(
        input_table,
        parse_options=csv.ParseOptions(delimiter="\t"),
        read_options=csv.ReadOptions(autogenerate_column_names=True),
        convert_options=csv.ConvertOptions(
            null_values=null_values, strings_can_be_null=True
        ),
    )

except IOError:
    print("Error reading/parsing: " + input_table)
    exit(2)


try:
    pq.write_table(table, output_parquet)
except IOError:
    print("Error writing parquet: " + output_parquet)
    exit(2)
