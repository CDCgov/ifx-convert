# ifx-convert: a small software collection of bioinformatics conversion utilities

`ifx-convert`is a small software collection of bioinformatics conversion utilities used by CDC/NCIRD/ID Informatics. Most involve file format conversions though some may involve other types of conversion such as time format. This repository is maintained by Sam S. Shepard and David E. Hufnagel. For direct correspondence, feel free to contact: [Samuel S. Shepard](mailto:vfn4@cdc.gov), Centers for Disease Control and Prevention

## Brief descriptions

| Author(s)                                        | Script                | Convert From         | Convert To            | Notes                                                                                                                                                       |
| :----------------------------------------------- | :-------------------- | :------------------- | :-------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| S. S. Shepard                                    | `afa2phy.pl`          | Aligned FASTA        | Phylip (sequential)   |                                                                                                                                                             |
| S. S. Shepard                                    | `collectXSV.pl`       | Delimited (CSV, TSV) | concatenated CSV, TSV |                                                                                                                                                             |
| Unknown                                          | `convertHMMformat.pl` | SAM                  | HMMER                 | Deprecated: find replacement                                                                                                                                |
| S. S. Shepard                                    | `csv2sqm.pl`          | CSV                  | Square Matrix         | Esoteric: used for [LABEL training](https://git.biotech.cdc.gov/vfn4/label)                                                                                 |
| S. S. Shepard                                    | `date2epiweek.pl`     | Date                 | Epi Week              | Deprecated: see [UDF bioutils](https://git.biotech.cdc.gov/flu-informatics/udf-bioutils) or [Zoe](https://git.biotech.cdc.gov/vfn4/zoe)                     |
| S. S. Shepard, port of a port from J. Bahl (PHP) | `decimal2date.pl`     | Date                 | Decimal date          |                                                                                                                                                             |
| S. S. Shepard, ported by work from J. Bahl (PHP) | `decimal2date.py`     | Date                 | Decimal date          |                                                                                                                                                             |
| S. S. Shepard                                    | `delim2fasta.pl`      | TSV                  | FASTA                 | Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)                                                                          |
| S. S. Shepard                                    | `fa2delim.pl`         | FASTA                | TSV                   | Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)                                                                          |
| S. S. Shepard                                    | `fa2fastq.pl`         | FASTA                | FASTQ                 | Minimal. For reverse, see also fastQ-converter in [IRMA](https://git.biotech.cdc.gov/vfn4/irma) and [IRMA-core](https://git.biotech.cdc.gov/vfn4/irma-core) |
| D. E. Hufnagel                                   | `genbank2Fasta.py`    | GenBank              | FASTA                 |                                                                                                                                                             |
| S. S. Shepard                                    | `nt2aa.pl`            | NT FASTA             | AA FASTA              | Amino acid translation. Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)                                                  |
| S. S. Shepard                                    | `sam2fasta.pl`        | SAM                  | FASTA                 | Used by [DAIS-ribosome](https://git.biotech.cdc.gov/flu-informatics/dais-ribosome)                                                                          |
| S. S. Shepard                                    | `sam2qtbl.pl`         | SAM                  | TSV                   | Quality table for diagnostics                                                                                                                               |
| S. S. Shepard                                    | `sqm2csv.pl`          | Square Matrix        | CSV                   | Esoteric: used for [LABEL training](https://git.biotech.cdc.gov/vfn4/label)                                                                                 |
| S. S. Shepard                                    | `tab2parquet.py`      | TSV                  | Parquet               |                                                                                                                                                             |

## Notices

Copyright is public domain but [attributions to the original authors are welcomed](CITATION.bib). The software license uses [ASL2](http://www.apache.org/licenses/LICENSE-2.0.html).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal public domain dedication. All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see <http://www.apache.org/licenses/LICENSE-2.0.html>

The source code forked from other open source projects will inherit its license.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).
