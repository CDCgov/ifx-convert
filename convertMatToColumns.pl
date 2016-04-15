#!/usr/bin/env perl
# Sam Shepard - 2013

use Getopt::Long;
GetOptions( 	'csv|C=s'=> \$csv,
		'use-index|I' => \$useIndex,
		'dist-name|N=s'=>\$distName, 
		'sample-set-name|S=s' => \$sampleSet,
		'vsX|V=s' => \$vsX,
		'full-matrix|F' => \$fullMatrix
		);

# VERSUS X
# Annotation field
# Strip annots from ID
if( scalar(@ARGV) != 1 ) {
	die("Usage:\n\tperl $0 <square_matrix.tab> [-C|--csv <file.csv>] [-I|--use-index] [-N|--dist-name <NAME>]\n");
}


if ( !defined($sampleSet) ) {
	$sampleSet = '';
} else {
	$sampleSet = ",". $sampleSet;
}

open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0].\n");
$/ = "\n"; $i = 0;
while($line = <IN>) {
	chomp($line);
	@values = split("\t",$line);
	$id = shift(@values);
	$distHash{$id} = [@values];
	$indexByID{$id} = $i;
	@IDs[$i] = $id;
	$i++;
}
close(IN);
$N = $i;
if ( !defined($distName) ) {
	$distName = 'dist';
}

if ( defined($vsX) ) {
	print 'vsX,',$distName,"\n";
} else {
	print 'id1,id2,',$distName,"\n";
}
if ( defined($csv) ) {
	open(IN,'<', $csv) or die("Cannot open $csv.\n");
	$header = <IN>;
	while($line = <IN>) {
		chomp($line);
		@values = split(',',$line);
		$id1 = $values[0];
		$id2 = $values[1];

		if ( !defined($vsX) ) {
			print $id1,',',$id2,',',$distHash{$id1}[$indexByID{$id2}],$sampleSet,"\n";
		} elsif( $id1 eq $vsX ) {
			print $id2,',',$distHash{$id1}[$indexByID{$id2}],$sampleSet,"\n";
		} elsif( $id2 eq $vsX ) {
			print $id1,',',$distHash{$id1}[$indexByID{$id2}],$sampleSet,"\n";
		}
	}
	close(IN);
} else {
	if ( $fullMatrix ) {
		for($i=0;$i<$N;$i++) {
			for($j=0;$j<$N;$j++) {
				if ( !defined($vsX ) ) {
					if ( $useIndex ) {
						print $i,',',$j,',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],',',$IDs[$j],',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$i] eq $vsX ) {
					if ( $useIndex ) {
						print $j,',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$j],',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$j] eq $vsX ) {
					if ( $useIndex ) {
						print $i,',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				}
			}
		}
	} else {
		for($i=1;$i<$N;$i++) {
			for($j=0;$j<$i;$j++) {
				if ( !defined($vsX ) ) {
					if ( $useIndex ) {
						print $i,',',$j,',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],',',$IDs[$j],',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$i] eq $vsX ) {
					if ( $useIndex ) {
						print $j,',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$j],',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$j] eq $vsX ) {
					if ( $useIndex ) {
						print $i,',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],',',$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				}
			}
		}
	}
}
