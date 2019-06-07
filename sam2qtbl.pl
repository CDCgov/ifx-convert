#!/usr/bin/env perl
# Sam Shepard -- print sam as quality table for diagnostics
# 6.2019

use Storable;
use Getopt::Long;
GetOptions(
		'typical-alignment|T' => \$typicalFormat,
		'reverse-complement|R' => \$rev
	);

if ( scalar(@ARGV) != 3 && scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> [prefix]\n";
	$message .= "\t\t-T|--typical-alignment\t\tTypical alignment format.\n";
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
$REgetMolID = qr/(.+?)[_ ]([123]):.+/;

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
if ( $typicalFormat ) {
	$D = '-'; $N = '-';	
#SAM-like format
} else {
	$D = '.'; $N = 'N';	
}


for($K=0;$K<scalar(@sam);$K++) {
	if ( substr($sam[$K],0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
	if ( $qname =~ $REgetMolID ) {
		$qMolID = $1;
		$qSide = $2;
	} else {
		$qMolID = $qname;
		$qSide = '';
	}

	$strand = (($flag & 16) == 16) ? '-' : '+';

	if ( defined($references{$rn}) ) {
		$REF_N = $references{$rn};
		@NTs = split('',uc($seq));
		@cigars = split('',$cigar);
		@QS = unpack("c* i*",$qual);
		$rpos=$pos-1;
		$qpos=0;
		
		$aln = '';
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				for(1..$inc) {
					print STDOUT $qMolID,"\t",$qSide,"\t",$op,"\t",($rpos+1),"\t",($qpos+1),"\t",$NTs[$qpos],"\t",($QS[$qpos]-33),"\t",$strand,"\n";
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				$aln .= '-' x $inc;
				for(1..$inc) {
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
			    $q_insert = 0; 
			    foreach my $x ( @QS[($qpos+1)..($qpos+$inc)] ) {
		        	$q_insert += $x;
		        }
		        $q_insert = ($q_insert-$inc*33)/$inc;
				print STDOUT $qMolID,"\t",$qSide,"\t",$op,"\t",($rpos),"\t",($qpos+1),"\t",substr($seq,$qpos,$inc),"\t",$q_insert,"\t",$strand,"\n";
				$qpos += $inc;
			} elsif( $op eq 'S' ) {
				$qpos += $inc;
			} elsif( $op eq 'N' ) {
				$rpos += $inc;
			} elsif ( $op eq 'H' ) {
				next;
			} else {
				die("Extended CIGAR ($op) not yet supported.\n");
			}
		}

#		if ( defined($rev) && $flag == 16 ) {
#			$tmp = reverse($preAln.$aln.$postAln);
#			$tmp =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
#			print $tmp,"\n";
#		}
	}
}
