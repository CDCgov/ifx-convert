# `CONVERT` collection

A collection of miscellaneous conversion scripts. Most involve file format conversions though some may involve other types of conversion such as time format.

## Brief descriptions

| Author(s)                                        | Script                | Convert From         | Convert To            | Notes                                               |
| :----------------------------------------------- | :-------------------- | :------------------- | :-------------------- | :-------------------------------------------------- |
| S. S. Shepard                                    | `afa2phy.pl`          | Aligned FASTA        | Phylip (sequential)   |                                                     |
| S. S. Shepard                                    | `collectXSV.pl`       | Delimited (CSV, TSV) | concatenated CSV, TSV |                                                     |
| Unknown                                          | `convertHMMformat.pl` | SAM                  | HMMER                 | Deprecated: find replacement                        |
| S. S. Shepard                                    | `csv2sqm.pl`          | CSV                  | Square Matrix         | Esoteric: used for [LABEL training](https://git.biotech.cdc.gov/vfn4/label)|
| S. S. Shepard                                    | `date2epiweek.pl`     | Date                 | Epi Week              | Deprecated: see [UDF bioutils](https://git.biotech.cdc.gov/flu-informatics/udf-bioutils) or [Zoe](https://git.biotech.cdc.gov/vfn4/zoe)|
| S. S. Shepard, port of a port from J. Bahl (PHP) | `decimal2date.pl`     | Date                 | Decimal date          |                                                     |
| S. S. Shepard, ported by work from J. Bahl (PHP) | `decimal2date.py`     | Date                 | Decimal date          |                                                     |
| S. S. Shepard                                    | `delim2fasta.pl`      | TSV                  | FASTA                 | Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)|
| S. S. Shepard                                    | `fa2delim.pl`         | FASTA                | TSV                   | Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)|
| S. S. Shepard                                    | `fa2fastq.pl`         | FASTA                | FASTQ                 | Minimal. For reverse, see also fastQ-converter in [IRMA](https://git.biotech.cdc.gov/vfn4/irma) and [IRMA-core](https://git.biotech.cdc.gov/vfn4/irma-core)|
| D. E. Hufnagel    |   `genbank2Fasta.py`    |   GenBank |   FASTA |
| S. S. Shepard                                    | `nt2aa.pl`            | NT FASTA             | AA FASTA              | Amino acid translation. Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)|
| S. S. Shepard                                    | `sam2fasta.pl`        | SAM                  | FASTA                 | Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)|
| S. S. Shepard                                    | `sam2qtbl.pl`         | SAM                  | TSV                   | Quality table for diagnostics                       |
| S. S. Shepard                                    | `sqm2csv.pl`          | Square Matrix        | CSV                   | Esoteric: used for [LABEL training](https://git.biotech.cdc.gov/vfn4/label)|
| S. S. Shepard                                    | `tab2parquet.py`      | TSV                  | Parquet               |                                                     |

## License

Works currently do not have a specified license, but may generally be considered public domain where written exclusively by S. S. Shepard (@vfn4).
