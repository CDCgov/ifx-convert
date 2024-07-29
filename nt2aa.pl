#!/usr/bin/env perl
#
# Filename:         nt2aa.pl
#
# Description:      Converts codons to amino acids with extra options for
#                   padding, missing data, and termination.
#
# Date dedicated:   2023-06-23
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Unpublished
#
# =============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#  This source code file or script constitutes a work of the United States
#  Government and is not subject to domestic copyright protection under 17 USC §
#  105. This file is in the public domain within the United States, and
#  copyright and related rights in the work worldwide are waived through the CC0
#  1.0 Universal public domain dedication:
#  https://creativecommons.org/publicdomain/zero/1.0/
#
#  The material embodied in this software is provided to you "as-is" and without
#  warranty of any kind, express, implied or otherwise, including without
#  limitation, any warranty of fitness for a particular purpose. In no event
#  shall the Centers for Disease Control and Prevention (CDC) or the United
#  States (U.S.) government be liable to you or anyone else for any direct,
#  special, incidental, indirect or consequential damages of any kind, or any
#  damages whatsoever, including without limitation, loss of profit, loss of
#  use, savings or revenue, or the claims of third parties, whether or not CDC
#  or the U.S. government has been advised of the possibility of such loss,
#  however caused and on any theory of liability, arising out of or in
#  connection with the possession, use or performance of this software.
#
#  Please provide appropriate attribution in any work or product based on this
#  material.
#
## no critic (ControlStructures::ProhibitCStyleForLoops)

use English qw(-no_match_vars);
use Getopt::Long;
use warnings;

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
            'right-pad-cds'        => \$rightPadCDS
);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
    die(  "Usage:\n\tperl $PROGRAM_NAME <nts.fasta> [options]\n"
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
        . "\t\t--right-pad-cds\t\t\tPreserve any right padding for 3' missing data (masking CDS for early termination). Assumes -U.\n"
        . "\n" );
}

#<<< augmented translation table
%gc = (
	'TAA'=>'*','TAG'=>'*','TAR'=>'*','TGA'=>'*','TRA'=>'*','GCA'=>'A','GCB'=>'A','GCC'=>'A','GCD'=>'A','GCG'=>'A','GCH'=>'A',
	'GCK'=>'A','GCM'=>'A','GCN'=>'A','GCR'=>'A','GCS'=>'A','GCT'=>'A','GCV'=>'A','GCW'=>'A','GCY'=>'A','TGC'=>'C','TGT'=>'C',
	'TGY'=>'C','GAC'=>'D','GAT'=>'D','GAY'=>'D','GAA'=>'E','GAG'=>'E','GAR'=>'E','TTC'=>'F','TTT'=>'F','TTY'=>'F','GGA'=>'G',
	'GGB'=>'G','GGC'=>'G','GGD'=>'G','GGG'=>'G','GGH'=>'G','GGK'=>'G','GGM'=>'G','GGN'=>'G','GGR'=>'G','GGS'=>'G','GGT'=>'G',
	'GGV'=>'G','GGW'=>'G','GGY'=>'G','CAC'=>'H','CAT'=>'H','CAY'=>'H','ATA'=>'I','ATC'=>'I','ATH'=>'I','ATM'=>'I','ATT'=>'I',
	'ATW'=>'I','ATY'=>'I','AAA'=>'K','AAG'=>'K','AAR'=>'K','CTA'=>'L','CTB'=>'L','CTC'=>'L','CTD'=>'L','CTG'=>'L','CTH'=>'L',
	'CTK'=>'L','CTM'=>'L','CTN'=>'L','CTR'=>'L','CTS'=>'L','CTT'=>'L','CTV'=>'L','CTW'=>'L','CTY'=>'L','TTA'=>'L','TTG'=>'L',
	'TTR'=>'L','YTA'=>'L','YTG'=>'L','YTR'=>'L','ATG'=>'M','AAC'=>'N','AAT'=>'N','AAY'=>'N','CCA'=>'P','CCB'=>'P','CCC'=>'P',
	'CCD'=>'P','CCG'=>'P','CCH'=>'P','CCK'=>'P','CCM'=>'P','CCN'=>'P','CCR'=>'P','CCS'=>'P','CCT'=>'P','CCV'=>'P','CCW'=>'P',
	'CCY'=>'P','CAA'=>'Q','CAG'=>'Q','CAR'=>'Q','AGA'=>'R','AGG'=>'R','AGR'=>'R','CGA'=>'R','CGB'=>'R','CGC'=>'R','CGD'=>'R',
	'CGG'=>'R','CGH'=>'R','CGK'=>'R','CGM'=>'R','CGN'=>'R','CGR'=>'R','CGS'=>'R','CGT'=>'R','CGV'=>'R','CGW'=>'R','CGY'=>'R',
	'MGA'=>'R','MGG'=>'R','MGR'=>'R','AGC'=>'S','AGT'=>'S','AGY'=>'S','TCA'=>'S','TCB'=>'S','TCC'=>'S','TCD'=>'S','TCG'=>'S',
	'TCH'=>'S','TCK'=>'S','TCM'=>'S','TCN'=>'S','TCR'=>'S','TCS'=>'S','TCT'=>'S','TCV'=>'S','TCW'=>'S','TCY'=>'S','ACA'=>'T',
	'ACB'=>'T','ACC'=>'T','ACD'=>'T','ACG'=>'T','ACH'=>'T','ACK'=>'T','ACM'=>'T','ACN'=>'T','ACR'=>'T','ACS'=>'T','ACT'=>'T',
	'ACV'=>'T','ACW'=>'T','ACY'=>'T','GTA'=>'V','GTB'=>'V','GTC'=>'V','GTD'=>'V','GTG'=>'V','GTH'=>'V','GTK'=>'V','GTM'=>'V',
	'GTN'=>'V','GTR'=>'V','GTS'=>'V','GTT'=>'V','GTV'=>'V','GTW'=>'V','GTY'=>'V','TGG'=>'W','TAC'=>'Y','TAT'=>'Y','TAY'=>'Y'
);
#>>>

my $UPDATED;
if ( defined $writeUpdated ) {
    open( $UPDATED, '>', $writeUpdated ) or die("Cannot open '$writeUpdated' for writing: $OS_ERROR\n");
}

# Process records.
local $RS = ">";
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/smx, $record );
    $header   = shift(@lines);
    $sequence = join( q{}, @lines );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    } elsif ( $length % 3 != 0 ) {
        if ( defined $selectFrame ) {
            $selectedFrame = 0;
            $maxCoding     = 0;
            $minStop       = 0;
            foreach my $frame ( 0 .. 2 ) {
                $numCoding = 0;
                $numStop   = 0;
                for ( my $i = $frame; $i < $length; $i += 3 ) {
                    $codon = uc( substr( $sequence, $i, 3 ) );
                    if ( $codon !~ /[-.~]/smx && defined $gc{$codon} && $gc{$codon} ne '*' ) {
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

            $length = $length2;
            if ( $minStop > 1 ) {
                print STDERR "WARNING, $header not in triplets ($length).\n";
                next;
            }

        } elsif ( defined $warningSkip ) {
            print STDERR "WARNING, $header not in triplets ($length).\n";
            next;
        } else {
            die("$header not in triplets ($length).\n");
        }
    }

    if ($adjustGaps) {
        while ( $sequence =~ /([A-Za-z]{3})((---)+)([A-Za-z]{3})/gsmx ) {
            my ( $left_side, $gaps, $right_side ) = ( $1, $2, $4 );
            $frame = $LAST_MATCH_START[1] % 3;
            if ( $frame == 1 ) {

                # RIGHT SHIFT
                $replacement = substr( $left_side, 0, -1 ) . $gaps . substr( $left_side, -1 );
                substr( $sequence, $LAST_MATCH_START[1], length($replacement), $replacement );
            } elsif ( $frame == 2 ) {

                # LEFT SHIFT
                $replacement = $left_side . substr( $right_side, 0, 1 ) . $gaps;
                substr( $sequence, $LAST_MATCH_START[1], length($replacement), $replacement );
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
            $codons =~ s/(\.{3})+$//smx;
        }

        if ( length($codons) % 3 != 0 ) {
            die("Unexpected error: $header not in codon triplets.\n");
        }

        print STDOUT '>', $header, "\n", $codons, "\n";

    } else {
        $aa = q{};
        my $new_length = $length;
        for ( $i = 0; $i < $length; $i += 3 ) {
            $codon = uc( substr( $sequence, $i, 3 ) );
            if ( $codon =~ /[-.~]/smx ) {
                if ($stripGapped) {
                    next;
                } else {
                    if ( $missingData && $codon eq '...' ) {
                        $aa .= '.';
                    } elsif ( $usePartialCodons && $codon =~ /[^.~-]/smx ) {
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
                    last;
                }
            }
        }

        if ($noEndCodon) {
            chop($aa);
        }

        if ($chopDownstreamMissing) {
            $aa =~ s/([^.])(\.+)$/$1/smx;
        }

        print STDOUT '>', $header, "\n", $aa, "\n";

        if ( defined $writeUpdated ) {
            my $codon_length = length($aa) * 3;
            my $remainder    = $length - $codon_length;

            if ($rightPadCDS) {
                substr( $sequence, $codon_length, $remainder, '.' x $remainder );
                print $UPDATED '>', $header, "\n", $sequence, "\n";
            } else {
                print $UPDATED '>', $header, "\n", substr( $sequence, 0, $codon_length ), "\n";
            }
        }
    }
}

if ( defined $writeUpdated ) {
    close $UPDATED or croak("Cannot close file: $OS_ERROR\n");
}
