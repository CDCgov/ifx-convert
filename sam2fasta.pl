#!/usr/bin/env perl
# Sam Shepard -- print sam as fasta format
# 3.2014

use Storable;
use Getopt::Long;
GetOptions(
		'use-storable|S' => \$useStorable,
		'typical-alignment|T' => \$typicalFormat,
		'print-inserts|P' => \$printInserts,
		'print-inline-inserts|I' => \$inlineInserts,
		'amend-missing|M=s' => \$amendMissingFile,
		'extend-to-stop|E' => \$extendToStop,
		'reverse-complement|R' => \$rev,
		'add-header-field|H=s' => \$headerField
	);

if ( scalar(@ARGV) != 3 && scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> [prefix]\n";
	$message .= "\t\t-O|--use-storable\t\tUse storable objects rather than FASTA.\n";
	$message .= "\t\t-T|--typical-alignment\t\tTypical alignment format.\n";
	$message .= "\t\t-P|--print-inserts\t\tPrint insertion table (STDERR or file based on prefix).\n";
	$message .= "\t\t-I|--print-inline-inserts\tPrint inline insertions as lowercase.\n";
	$message .= "\t\t-M|--amend-missing <FILE>\tAmend missing 5' and 3' from alignment using the original file.\n";
	$message .= "\t\t-E|--extend-to-stop\t\tAssumes -M is used. Allows post-sequence insertion extension.\n";
	$message .= "\t\t-R|--reverse-complement\t\tReverse complement the read if original was on complementary strand.\n";
	$message .= "\t\t-H|--add-header-field <STR>\tAdds constant header field (pipe-delimited) to ID using specified argument.\n";
	die($message."\n");
}

if ( scalar(@ARGV) == 2 ) {
	if ( $useStorable ) { die("Need the [prefix] argument for storable results.\n"); }
	$stdout = 1;
} else {
	$stdout = 0;
}

$REisBase = qr/[ATCG]/;
$REgetMolID = qr/(.+?)[_ ]([12]):.+/;
$headerField = defined($headerField) ? '|'.$headerField : '';

if ( defined($amendMissingFile) ) {
	open(MISS,'<',$amendMissingFile) or die("Cannot open $amendMissingFile for reading.\n");
	$/ = ">"; %seqByID = ();
	while($record = <MISS>) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$id = shift(@lines);
		$seq = uc(join('',@lines));
		if ( length($seq) < 1 ) {
			next;
		}
		$seqByID{$id} = $seq;
	}
	close(MISS);
	$amendMissing = 1;
} else {
	$amendMissing = 0;
}

$/ = '>';
my %references = ();
my $REF_NAME = '';
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$seq = join('',@lines);
	$length = length($seq);
	if ( $length < 1 ) {
		next;
	} else {
		$references{$REF_NAME} = $length;
	}
}
close(REF);

$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
@sam = <SAM>; chomp(@sam);
close(SAM);

# baseline gap & spacer format
if ( $inlineInserts ) {
	$D = ''; $N = '.';
} elsif ( $typicalFormat ) {
	$D = '-'; $N = '-';	
#SAM-like format
} else {
	$D = '.'; $N = 'N';	
}

if ( $useStorable ) {
	%pairs = %insByIndex = ();
	for($K=0;$K<scalar(@sam);$K++) {
		if ( substr($sam[$K],0,1) eq '@' ) {
			next;
		}

		($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
		if ( $qname =~ $REgetMolID ) {
			$qMolID = $1;
			$qSide = $2;
		}

		if ( defined($references{$rn}) ) {
			$REF_N = $references{$rn};
			@NTs = split('',uc($seq));
			@QCs = split('',$qual);
			@Qint = unpack("c* i*",$qual);
			@cigars = split('',$cigar);
			$rpos=$pos-1;
			$qpos=0;
			
			if ( $rpos > 0 ) {
				$aln = $D x $rpos; 
				$qAln = ' ' x $rpos;
			} else {
				$aln = '';
				$qAln = '';
			}
			
			while($cigar =~ /(\d+)([MIDNSHP])/g ) {
				$inc=$1;
				$op=$2;
				if ( $op eq 'M' ) {
					for(1..$inc) {
						$qAln .= $QCs[$qpos];
						$aln .= $NTs[$qpos];
						$qpos++; $rpos++;
					}
				} elsif ( $op eq 'D' ) {
					$qAln .= ' ' x $inc;
					$aln .= '-' x $inc;
					for(1..$inc) {
						$rpos++;
					}
				} elsif( $op eq 'I' ) {
					$insByIndex{$K}{$rpos} = [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
					$qpos += $inc;
				} elsif( $op eq 'S' ) {
					$qpos += $inc;
				} elsif( $op eq 'N' ) {
					$aln .= $N x $inc;
					$qAln .= ' ' x $inc;
					$rpos += $inc;
				} elsif ( $op eq 'H' ) {
					next;
				} else {
					die("Extended CIGAR ($op) not yet supported.\n");
				}
			}
			$aln .= $D x (($REF_N)-$rpos);
			$qAln .= ' ' x (($REF_N)-$rpos);
			$pairs{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
		}
	}
	store(\%pairs, $ARGV[2].'.aln');
	store(\%insByIndex, $ARGV[2].'.ins');
} else {
	if ( $stdout ) {
		*FASTA = *STDOUT;
		if ( $printInserts ) {
			*INSRT = *STDERR;
		}
	} else {
		open(FASTA,'>',$ARGV[2].'.fasta') or die("Cannot open $ARGV[2].fasta for writing.\n");
		if ( $printInserts ) {
			open(INSRT,'>',$ARGV[2].'.ins.txt') or die("Cannot open $ARGV[2].ins.txt for writing.\n");
		}
	}

	for($K=0;$K<scalar(@sam);$K++) {
		if ( substr($sam[$K],0,1) eq '@' ) {
			next;
		}

		($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
		if ( $qname =~ $REgetMolID ) {
			$qMolID = $1;
			$qSide = $2;
		}

		if ( defined($references{$rn}) ) {
			$REF_N = $references{$rn};
			@NTs = split('',uc($seq));
			@cigars = split('',$cigar);
			$rpos=$pos-1;
			$qpos=0;
			
			if ( $rpos > 0 ) {
				$preAln = $D x $rpos; 
			} else {
				$preAln = ''; 
			}
			
			$aln = '';
			while($cigar =~ /(\d+)([MIDNSHP])/g ) {
				$inc=$1;
				$op=$2;
				if ( $op eq 'M' ) {
					for(1..$inc) {
						$aln .= $NTs[$qpos];
						$qpos++; $rpos++;
					}
				} elsif ( $op eq 'D' ) {
					$aln .= '-' x $inc;
					for(1..$inc) {
						$rpos++;
					}
				} elsif( $op eq 'I' ) {
					if ( $inlineInserts ) { $aln .= lc(substr($seq,$qpos,$inc)); }
					if ( $printInserts ) { print INSRT $qname,$headerField,"\t",($rpos),"\t",substr($seq,$qpos,$inc),"\n"; }
					if ( $amendMissing ) { $inserts{$qname}{$rpos} = substr($seq,$qpos,$inc); }
					$qpos += $inc;
				} elsif( $op eq 'S' ) {
					$qpos += $inc;
				} elsif( $op eq 'N' ) {
					$aln .= $N x $inc;
					$rpos += $inc;
				} elsif ( $op eq 'H' ) {
					next;
				} else {
					die("Extended CIGAR ($op) not yet supported.\n");
				}
			}

			$postAln = '';
			$postAln = $D x (($REF_N)-$rpos);

			if ( $amendMissing ) {
				if ( defined($seqByID{$qname}) ) {
					$unaligned = $preAln.$aln.$postAln; $offset = 0;
					foreach $pos ( sort { $a <=> $b } keys(%{$inserts{$qname}}) ) {
						$insert = $inserts{$qname}{$pos};

						substr($unaligned,$pos+$offset,0) = $insert;
						$offset += length($insert);
					}
					$unaligned =~ tr/.-//d;		# unaligned sequence, removing deletions and adding back insertions

					$original = $seqByID{$qname};
					$O = length($original);		# length of the original query sequence
					$L = length($aln);		# length of alignable sequence to the reference

					if ( $original =~ /\Q$unaligned\E/ ) {
						# Alignable sequence to Start & Stop in Original Sequence Coords
						($start,$stop) = ($-[0],$+[0]);
						# Nucleotide count between the alignable and original sequence on 5' & 3' ends (pre/post)
						($preO,$postO) = ( $start, ($O - $stop) );
						# Nucleotide count between the alignable and reference sequence on 5' & 3' ends (pre/post)
						($preR,$postR) = (length($preAln),length($postAln));
					
						# The nucleotides to Move is the lesser of the available nucleotides from the original,
						# and the nucleotides needed to complete the reference alignment.
						if ( $preR > $preO ) { $preM = $preO; } else { $preM = $preR; }
						if ( $postR > $postO ) { $postM = $postO; } else { $postM = $postR; }

						# Only move over nucleotides if there are nucleotides that are movable.
						# AND, the unalignable flanking region to the reference shall not exceed 9 nucleotides.
						if ( $preM > 0 && $preR < 10 ) {
							$preAln = ($D x ($preR-$preM)) . substr($original,$start - $preM,$preM);
						}

						if ( $postM > 0 && $postR < 10 ) {
							$postAln = substr($original,$stop,$postM) . ($D x ($postR-$postM));
							# if insertions are to be added, the alignable position has been increased
							$L += $postM;	
						}
					}

					# We may also extend our alignment to the first stop, and write out 3' insertions to reference.
					if ( defined($extendToStop) ) {
						$suffix = $aln.$postAln;
						$last3 = substr($suffix,-3);
						if ( $last3 =~ /[A-Za-z]{3}/ && $last3 !~ /(TGA|TAA|TAG|TAR|TRA)/i ) {
							if ( $original =~ /\Q$suffix\E/ ) {
								$start2 = $+[0];
								$tail  = substr($original,$start2,$O-$start2);
								$T = length($tail);
								if ( $T >= 3 ) {
									$i = 0; $codons = '';
									for($i=0;$i<$T;$i+=3) {
										$codon = substr($tail,$i,3);
										$codons .= $codon;

										if ( $codon =~ m/(TGA|TAA|TAG|TAR|TRA)/i ) {
											# Given that end of suffix is not empty, reference length (O) could be used.
											# However, there may have been a reason in the original code not to do this.
											# For the most part: L+len(preAln)+postM == O
											# I theorize that for (suffix = aln) and postAln empty,  
											# non-shunt extension may take place. If shunting occurs aln != suffix.
											# Under that scenario it may be best to leave as-is.
											print INSRT $qname,$headerField,"\t",($L+length($preAln)),"\t",$codons,"\n";
											last;
										}
									}
								}
							}
						}
					}
				}
			}

			print FASTA '>',$qname,$headerField,"\n";
			if ( defined($rev) && $flag == 16 ) {
				$tmp = reverse($preAln.$aln.$postAln);
				$tmp =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
				print $tmp,"\n";
			} else {
				print FASTA $preAln,$aln,$postAln,"\n";
			}
		}
	}

	if ( !$stdout ) {
		if ( $printInserts ) {
			close(INSRT);
		}
		close(FASTA);
	}
}
