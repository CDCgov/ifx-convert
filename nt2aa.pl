#!/usr/bin/env perl
# Samuel Shepard - 2.2015
# Version 1.0
# Converts codons to amino acids.

use Getopt::Long;
GetOptions(	'strip-gapped|S'=>\$stripGapped, 
		'no-end|N' => \$noEndCodon,
		'warning-skip|W' => \$warningSkip,
		'reading-frame-mode|R' => \$selectFrame,
		'partial-codon|P' => \$usePartialCodons,
		'adjust-gaps|A' => \$adjustGaps,
		'missing-data|M' => \$missingData
	);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <nts.fasta> [options]\n";
	$message .= "\t\t--strip-gapped|-S\tStrip gaps before translating.\n";
	$message .= "\t\t--warning-skip|-W\tSkip and print a warning rather than throwing an error.\n";
	$message .= "\t\t--no-end|-N\t\tSkip the last codon (useful when it is the stop site).\n";
	$message .= "\t\t--read-frame-mode|-R\tFind longest ORF in sequence for translation.\n";
	$message .= "\t\t--adjust-gaps|-A\tTry to adjust partial codons to be full codons.\n";
	$message .= "\t\t--partial-codon|-P\tUse '~' to express partial codons.\n";
	$message .= "\t\t--missing-data|-M\tUse '.' to express missing data.\n";
	die($message."\n");
}

# mappings of nucleotides
%M = (
	'AT' => 'W', 'CG' => 'S', 'AC' => 'M', 'GT' => 'K', 'AG' => 'R', 'CT' => 'Y',
	'CGT' => 'B', 'AGT' => 'D', 'ACT' => 'H', 'ACG' => 'V', 'ACGT' => 'N'
);

%DM = (
        'W' => 'AT', 'S' => 'CG', 'M' => 'AC', 'K' => 'GT', 'R' => 'AG', 'Y' => 'CT',
	'B' => 'CGT', 'D' => 'AGT', 'H' => 'ACT', 'V' => 'ACG', 'N' => 'ACGT'
);


# Simple translation table.
%gc = (
	'TCA' => 'S', # Serine
	'TCC' => 'S', # Serine
	'TCG' => 'S', # Serine
	'TCT' => 'S', # Serine
	'TTC' => 'F', # Phenylalanine
	'TTT' => 'F', # Phenylalanine
	'TTA' => 'L', # Leucine
	'TTG' => 'L', # Leucine
	'TAC' => 'Y', # Tyrosine
	'TAT' => 'Y', # Tyrosine
	'TAA' => '*', # Stop
	'TAG' => '*', # Stop
	'TGC' => 'C', # Cysteine
	'TGT' => 'C', # Cysteine
	'TGA' => '*', # Stop
	'TGG' => 'W', # Tryptophan
	'CTA' => 'L', # Leucine
	'CTC' => 'L', # Leucine
	'CTG' => 'L', # Leucine
	'CTT' => 'L', # Leucine
	'CCA' => 'P', # Proline
	'CAT' => 'H', # Histidine
	'CAA' => 'Q', # Glutamine
	'CAG' => 'Q', # Glutamine
	'CGA' => 'R', # Arginine
	'CGC' => 'R', # Arginine
	'CGG' => 'R', # Arginine
	'CGT' => 'R', # Arginine
	'ATA' => 'I', # Isoleucine
	'ATC' => 'I', # Isoleucine
	'ATT' => 'I', # Isoleucine
	'ATG' => 'M', # Methionine
	'ACA' => 'T', # Threonine
	'ACC' => 'T', # Threonine
	'ACG' => 'T', # Threonine
	'ACT' => 'T', # Threonine
	'AAC' => 'N', # Asparagine
	'AAT' => 'N', # Asparagine
	'AAA' => 'K', # Lysine
	'AAG' => 'K', # Lysine
	'AGC' => 'S', # Serine
	'AGT' => 'S', # Serine
	'AGA' => 'R', # Arginine
	'AGG' => 'R', # Arginine
	'CCC' => 'P', # Proline
	'CCG' => 'P', # Proline
	'CCT' => 'P', # Proline
	'CAC' => 'H', # Histidine
	'GTA' => 'V', # Valine
	'GTC' => 'V', # Valine
	'GTG' => 'V', # Valine
	'GTT' => 'V', # Valine
	'GCA' => 'A', # Alanine
	'GCC' => 'A', # Alanine
	'GCG' => 'A', # Alanine
	'GCT' => 'A', # Alanine
	'GAC' => 'D', # Aspartic Acid
	'GAT' => 'D', # Aspartic Acid
	'GAA' => 'E', # Glutamic Acid
	'GAG' => 'E', # Glutamic Acid
	'GGA' => 'G', # Glycine
	'GGC' => 'G', # Glycine
	'GGG' => 'G', # Glycine
	'GGT' => 'G'  # Glycine
);

# Process records.
$/ = ">";
while($record = <> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	} elsif ( $length % 3 != 0 ) {
		if ( defined($selectFrame) ) {
			$selectedFrame = 0;
			$maxCoding = 0; $minStop = 0;
			for($frame=0;$frame<3;$frame++) {
				$numCoding = 0; $numStop = 0;
				for($i=$frame;$i<$length;$i+=3) {
					$codon = uc(substr($sequence, $i, 3));
					if ( $codon !~ /[-.~]/ && defined($gc{$codon}) && $gc{$codon} ne '*' ) {
						$numCoding++;
					} elsif ( $gc{$codon} eq '*' ) {
						$numStop++;
					}
				}
				if ( $numCoding > $maxCoding ) {
					$maxCoding = $numCoding;
					$selectedFrame = $frame;
					$minStop = $numStop;
				}
			}
			$length2 = int($length/3);
			$sequence = substr($sequence,$frame,$length2);
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


	if ( $adjustGaps ) {
		while ( $sequence =~ /([A-Za-z]{3})((---)+)([A-Za-z]{3})/g ) {
			($left,$gaps,$gapT,$right) = ($1,$2,$3,$4);
			$frame = $-[1]%3;
			if ( $frame == 1 ) {
			# RIGHT SHIFT
				$replacement = substr($left,0,-1).$gaps.substr($left,-1);
				substr($sequence,$-[1],length($replacement)) =  $replacement;
			} elsif ( $frame == 2 ) {
			# LEFT SHIFT
				$replacement = $left.substr($right,0,1).$gaps;
				substr($sequence,$-[1],length($replacement)) =  $replacement;
			}
		}
	}

	$aa = '';
	for( $i = 0; $i < $length; $i +=3 ) {
		$codon = uc(substr($sequence, $i, 3));
		if ( $codon =~ /[-.~]/ ) {
			if ( $stripGapped ) {
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
		} elsif ( !defined($gc{$codon}) ) {
			$aa .= 'X';
		} else {
			$aa .= $gc{$codon};
		}
	}

	if ( $noEndCodon ) {
		chop($aa);
	}
	print '>',$header,"\n",$aa,"\n";
}
close(IN);
