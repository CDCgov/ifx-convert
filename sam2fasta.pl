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
		'reverse-complement|R' => \$rev
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


$/ = ">"; $REF_NAME = $REF_SEQ = $REF_N = '';
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	$REF_SEQ = [split('',uc($seq))];
	$REF_N = length($seq);
	last;
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

		if ( $REF_NAME eq $rn ) {
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

		if ( $REF_NAME eq $rn ) {
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
					if ( $printInserts ) { print INSRT $qname,"\t",($rpos),"\t",substr($seq,$qpos,$inc),"\n"; }
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

					$aln2 = $preAln.$aln.$postAln; $offset = 0;
					foreach $pos ( sort { $a <=> $b } keys(%{$inserts{$qname}}) ) {
						$insert = $inserts{$qname}{$pos};

						substr($aln2,$pos+$offset,0) = $insert;
						$offset += length($insert);
					}
					$aln2 =~ tr/.-//d;
					$original = $seqByID{$qname};
					$O = length($original);
					$L = length($aln);
					if ( $original =~ /\Q$aln2\E/ ) {
						($start,$stop) = ($-[0],$+[0]);
						($preO,$postO) = ( $start, ($O - $stop) );
						($preL,$postL) = (length($preAln),length($postAln));
						
						if ( $preL > $preO ) { $preM = $preO; } else { $preM = $preL; }
						if ( $postL > $postO ) { $postM = $postO; } else { $postM = $postL; }

						if ( $preM > 0 && $preL < 10 ) {
							$preAln = ($D x ($preL-$preM)) . substr($original,$start - $preM,$preM);
						}

						if ( $postM > 0 && $postL < 10 ) {
							$postAln = substr($original,$stop,$postM) . ($D x ($postL-$postM));
						}
					}

					if ( defined($extendToStop) ) {
						$suffix = $aln.$postAln;
						$last3 = substr($suffix,-3);
						if ( $last3 =~ /[A-Za-z]{3}/ && $last3 !~ /(TGA|TAA|TAG)/i ) {
							if ( $original =~ /\Q$suffix\E/ ) {
								$start2 = $+[0];
								$tail  = substr($original,$start2,$O-$start2);
								$T = length($tail);
								if ( $T >= 3 ) {
									$i = 0; $codons = '';
									for($i=0;$i<$T;$i+=3) {
										$codon = substr($tail,$i,3);
										$codons .= $codon;
										if ( $codon =~ m/(TGA|TAA|TAG)/i ) {
											print INSRT $qname,"\t",$L,"\t",$codons,"\n";
											last;
										}
									}
								}
							}
						}
					}
				}
			}

			print FASTA '>',$qname,"\n";
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
