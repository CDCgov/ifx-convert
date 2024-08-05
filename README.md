# ifx-convert: a small software collection of bioinformatics conversion utilities

## Overview

This is the public repository to host CDC/NCIRD/ID Informatics conversion scripts. Most involve file format conversions though some may involve other types of conversion such as time format.

This repository is maintained by Sam S. Shepard and David E. Hufnagel. For direct correspondence, feel free to contact: [Samuel S. Shepard](mailto:vfn4@cdc.gov), Centers for Disease Control and Prevention

## License

This work is licensed under Apache 2.0 Licence. see LICENCE. Overall this means that these scripts can be used for any legal purpose including copying, modifying, or selling this code by internal or external parties so long as credit is given to the original authors. To meet this requirement please cite this github page when referring to code hosted here or when buliding code from this code.

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

General disclaimer: This repository was created for use by CDC programs to collaborate on public health related projects in support of [the CDC mission](https://www.cdc.gov/about/cdc/index.html). GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise.

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through [the CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md) and [Code of Conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md). For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://docs.github.com/en/get-started/getting-started-with-git/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the XXXXXXX License vXX or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through [the CDC web site](https://www.cdc.gov/).

### Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
