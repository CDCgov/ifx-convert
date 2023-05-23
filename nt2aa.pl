#!/usr/bin/env perl
# Samuel Shepard - 2.2015
# Version 1.0
# Converts codons to amino acids.

use English qw(-no_match_vars);
use Getopt::Long;
GetOptions(
            'strip-gapped|S'       => \$stripGapped,
            'no-end|N'             => \$noEndCodon,
            'warning-skip|W'       => \$warningSkip,
            'reading-frame-mode|R' => \$selectFrame,
            'partial-codon|P'      => \$usePartialCodons,
            'adjust-gaps|A'        => \$adjustGaps,
            'missing-data|M'       => \$missingData,
            'stop-translation|T'   => \$stopTranslation,
            'output-codons|C'      => \$codonsOnly,
            'end-3p-missing|E'     => \$chopDownstreamMissing,
            'write-updated-nt|U=s' => \$writeUpdated,
            'no-updated-nt-chop'   => \$noUpdatedNTChop
);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <nts.fasta> [options]\n"
         . "\t\t-S|--strip-gapped\t\tStrip gaps before translating.\n"
         . "\t\t-W|--warning-skip\t\tSkip and print a warning rather than throwing an error.\n"
         . "\t\t-N|--no-end\t\t\tSkip the last codon (useful when it is the stop site).\n"
         . "\t\t-R|--read-frame-mode\t\tFind longest ORF in sequence for translation.\n"
         . "\t\t-A|--adjust-gaps\t\tTry to adjust partial codons to be full codons.\n"
         . "\t\t-P|--partial-codon\t\tUse '~' to express partial codons.\n"
         . "\t\t-M|--missing-data\t\tUse '.' to express missing data.\n"
         . "\t\t-T|--stop-translation\t\tEnd sequence after first stop codon.\n"
         . "\t\t-U|--write-updated <FILE>\tWrite updated nucleotide FASTA (as with '-T').\n"
         . "\t\t-E|--end-3p-missing\t\tChop downstream missing.\n"
         . "\t\t--no-updated-nt-chop\t\tPreserve any right padding for 3' missing data (assumes -U).\n"
         . "\n" );
}

# ambiguation mappings of nucleotides
%M = (
       'AT'   => 'W',
       'CG'   => 'S',
       'AC'   => 'M',
       'GT'   => 'K',
       'AG'   => 'R',
       'CT'   => 'Y',
       'CGT'  => 'B',
       'AGT'  => 'D',
       'ACT'  => 'H',
       'ACG'  => 'V',
       'ACGT' => 'N'
);

# reverse mappings for non-canonical nucleotides
%DM = (
        'W' => 'AT',
        'S' => 'CG',
        'M' => 'AC',
        'K' => 'GT',
        'R' => 'AG',
        'Y' => 'CT',
        'B' => 'CGT',
        'D' => 'AGT',
        'H' => 'ACT',
        'V' => 'ACG',
        'N' => 'ACGT'
);

# augmented translation table
%gc = (
        'TAA' => '*',
        'TAG' => '*',
        'TAR' => '*',
        'TGA' => '*',
        'TRA' => '*',
        'GCA' => 'A',
        'GCB' => 'A',
        'GCC' => 'A',
        'GCD' => 'A',
        'GCG' => 'A',
        'GCH' => 'A',
        'GCK' => 'A',
        'GCM' => 'A',
        'GCN' => 'A',
        'GCR' => 'A',
        'GCS' => 'A',
        'GCT' => 'A',
        'GCV' => 'A',
        'GCW' => 'A',
        'GCY' => 'A',
        'TGC' => 'C',
        'TGT' => 'C',
        'TGY' => 'C',
        'GAC' => 'D',
        'GAT' => 'D',
        'GAY' => 'D',
        'GAA' => 'E',
        'GAG' => 'E',
        'GAR' => 'E',
        'TTC' => 'F',
        'TTT' => 'F',
        'TTY' => 'F',
        'GGA' => 'G',
        'GGB' => 'G',
        'GGC' => 'G',
        'GGD' => 'G',
        'GGG' => 'G',
        'GGH' => 'G',
        'GGK' => 'G',
        'GGM' => 'G',
        'GGN' => 'G',
        'GGR' => 'G',
        'GGS' => 'G',
        'GGT' => 'G',
        'GGV' => 'G',
        'GGW' => 'G',
        'GGY' => 'G',
        'CAC' => 'H',
        'CAT' => 'H',
        'CAY' => 'H',
        'ATA' => 'I',
        'ATC' => 'I',
        'ATH' => 'I',
        'ATM' => 'I',
        'ATT' => 'I',
        'ATW' => 'I',
        'ATY' => 'I',
        'AAA' => 'K',
        'AAG' => 'K',
        'AAR' => 'K',
        'CTA' => 'L',
        'CTB' => 'L',
        'CTC' => 'L',
        'CTD' => 'L',
        'CTG' => 'L',
        'CTH' => 'L',
        'CTK' => 'L',
        'CTM' => 'L',
        'CTN' => 'L',
        'CTR' => 'L',
        'CTS' => 'L',
        'CTT' => 'L',
        'CTV' => 'L',
        'CTW' => 'L',
        'CTY' => 'L',
        'TTA' => 'L',
        'TTG' => 'L',
        'TTR' => 'L',
        'YTA' => 'L',
        'YTG' => 'L',
        'YTR' => 'L',
        'ATG' => 'M',
        'AAC' => 'N',
        'AAT' => 'N',
        'AAY' => 'N',
        'CCA' => 'P',
        'CCB' => 'P',
        'CCC' => 'P',
        'CCD' => 'P',
        'CCG' => 'P',
        'CCH' => 'P',
        'CCK' => 'P',
        'CCM' => 'P',
        'CCN' => 'P',
        'CCR' => 'P',
        'CCS' => 'P',
        'CCT' => 'P',
        'CCV' => 'P',
        'CCW' => 'P',
        'CCY' => 'P',
        'CAA' => 'Q',
        'CAG' => 'Q',
        'CAR' => 'Q',
        'AGA' => 'R',
        'AGG' => 'R',
        'AGR' => 'R',
        'CGA' => 'R',
        'CGB' => 'R',
        'CGC' => 'R',
        'CGD' => 'R',
        'CGG' => 'R',
        'CGH' => 'R',
        'CGK' => 'R',
        'CGM' => 'R',
        'CGN' => 'R',
        'CGR' => 'R',
        'CGS' => 'R',
        'CGT' => 'R',
        'CGV' => 'R',
        'CGW' => 'R',
        'CGY' => 'R',
        'MGA' => 'R',
        'MGG' => 'R',
        'MGR' => 'R',
        'AGC' => 'S',
        'AGT' => 'S',
        'AGY' => 'S',
        'TCA' => 'S',
        'TCB' => 'S',
        'TCC' => 'S',
        'TCD' => 'S',
        'TCG' => 'S',
        'TCH' => 'S',
        'TCK' => 'S',
        'TCM' => 'S',
        'TCN' => 'S',
        'TCR' => 'S',
        'TCS' => 'S',
        'TCT' => 'S',
        'TCV' => 'S',
        'TCW' => 'S',
        'TCY' => 'S',
        'ACA' => 'T',
        'ACB' => 'T',
        'ACC' => 'T',
        'ACD' => 'T',
        'ACG' => 'T',
        'ACH' => 'T',
        'ACK' => 'T',
        'ACM' => 'T',
        'ACN' => 'T',
        'ACR' => 'T',
        'ACS' => 'T',
        'ACT' => 'T',
        'ACV' => 'T',
        'ACW' => 'T',
        'ACY' => 'T',
        'GTA' => 'V',
        'GTB' => 'V',
        'GTC' => 'V',
        'GTD' => 'V',
        'GTG' => 'V',
        'GTH' => 'V',
        'GTK' => 'V',
        'GTM' => 'V',
        'GTN' => 'V',
        'GTR' => 'V',
        'GTS' => 'V',
        'GTT' => 'V',
        'GTV' => 'V',
        'GTW' => 'V',
        'GTY' => 'V',
        'TGG' => 'W',
        'TAC' => 'Y',
        'TAT' => 'Y',
        'TAY' => 'Y'
);

my $UPDATED;
if ( defined $writeUpdated ) {
    open( $UPDATED, '>', $writeUpdated ) or die("Cannot open '$writeUpdated' for writing: $OS_ERROR\n");
}

# Process records.
local $RS = ">";
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = join( q{}, @lines );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    } elsif ( $length % 3 != 0 ) {
        if ( defined($selectFrame) ) {
            $selectedFrame = 0;
            $maxCoding     = 0;
            $minStop       = 0;
            for ( $frame = 0; $frame < 3; $frame++ ) {
                $numCoding = 0;
                $numStop   = 0;
                for ( $i = $frame; $i < $length; $i += 3 ) {
                    $codon = uc( substr( $sequence, $i, 3 ) );
                    if ( $codon !~ /[-.~]/ && defined( $gc{$codon} ) && $gc{$codon} ne '*' ) {
                        $numCoding++;
                    } elsif ( $gc{$codon} eq '*' ) {
                        $numStop++;
                    }
                }
                if ( $numCoding > $maxCoding ) {
                    $maxCoding     = $numCoding;
                    $selectedFrame = $frame;
                    $minStop       = $numStop;
                }
            }
            $length2  = int( $length / 3 );
            $sequence = substr( $sequence, $frame, $length2 );
            print STDERR "min $minStop num $maxCoding \n";
            $length = $length2;
            if ( $minStop > 1 ) {
                print STDERR "WARNING, $header not in triplets ($length).\n";
                next;
            }

        } elsif ( defined($warningSkip) ) {
            print STDERR "WARNING, $header not in triplets ($length).\n";
            next;
        } else {
            die("$header not in triplets ($length).\n");
        }
    }

    if ($adjustGaps) {
        while ( $sequence =~ /([A-Za-z]{3})((---)+)([A-Za-z]{3})/g ) {
            ( $left, $gaps, $gapT, $right ) = ( $1, $2, $3, $4 );
            $frame = $-[1] % 3;
            if ( $frame == 1 ) {

                # RIGHT SHIFT
                $replacement = substr( $left, 0, -1 ) . $gaps . substr( $left, -1 );
                substr( $sequence, $-[1], length($replacement) ) = $replacement;
            } elsif ( $frame == 2 ) {

                # LEFT SHIFT
                $replacement = $left . substr( $right, 0, 1 ) . $gaps;
                substr( $sequence, $-[1], length($replacement) ) = $replacement;
            }
        }
    }

    if ($codonsOnly) {
        $codons = q{};
        for ( $i = 0; $i < $length; $i += 3 ) {
            $codon = uc( substr( $sequence, $i, 3 ) );
            $codons .= $codon;
            if ( $stopTranslation && $gc{$codon} eq '*' ) {
                last;
            }
        }

        if ($noEndCodon) {
            $codons = substr( $codons, 0, -3 );

        }

        if ($chopDownstreamMissing) {
            $codons =~ s/(\.{3})+$//;
        }

        if ( length($codons) % 3 != 0 ) {
            die("Unexpected error: $header not in codon triplets.\n");
        }

        print '>', $header, "\n", $codons, "\n";

    } else {
        $aa = q{};
        my $new_length = $length;
        for ( $i = 0; $i < $length; $i += 3 ) {
            $codon = uc( substr( $sequence, $i, 3 ) );
            if ( $codon =~ /[-.~]/ ) {
                if ($stripGapped) {
                    next;
                } else {
                    if ( $missingData && $codon eq '...' ) {
                        $aa .= '.';
                    } elsif ( $usePartialCodons && $codon =~ /[^.~-]/ ) {
                        $aa .= '~';
                    } else {
                        $aa .= '-';
                    }
                }
            } elsif ( !defined $gc{$codon} ) {
                $aa .= 'X';
            } else {
                $aa .= $gc{$codon};
                if ( $stopTranslation && $gc{$codon} eq '*' ) {
                    if ( defined $writeUpdated ) {
                        $new_length = $i + 3;
                    }
                    last;
                }
            }
        }

        if ($noEndCodon) {
            chop($aa);
            if ( defined $writeUpdated ) {
                $new_length -= 3;
            }
        }

        if ($chopDownstreamMissing) {
            $aa =~ s/([^.])(\.+)$/$1/smx;
            if ( defined $writeUpdated && !defined $noUpdatedNTChop ) {
                $new_length -= length($2) * 3;
            }
        }

        print '>', $header, "\n", $aa, "\n";

        if ( defined $writeUpdated ) {
            print $UPDATED '>', $header, "\n", substr( $sequence, 0, $new_length ), "\n";
        }
    }
}

if ( defined $writeUpdated ) {
    close $UPDATED or croak("Cannot close file: $OS_ERROR\n");
}
