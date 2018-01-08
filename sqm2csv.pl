#!/usr/bin/env perl
# Sam Shepard - 2013

use Getopt::Long;
GetOptions( 	'csv|C=s'=> \$csv,
		'use-index|I' => \$useIndex,
		'dist-name|N=s'=>\$distName, 
		'set-name|S=s' => \$sampleSet,
		'vsX|V=s' => \$vsX,
		'full-matrix|F' => \$fullMatrix,
		'vs-list|L=s' => \$vsList,
		'remove-refs' => \$removeRefs,
		'delim|D=s' => \$delim
		);

# VERSUS X
# Annotation field
# Strip annots from ID
if( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <square_matrix.tab>\n";
	$message .= "\t\t-D|--delim\t\t\tDelimiter, default is TAB.\n";
	$message .= "\t\t-C|--csv <file.csv>\t\tA different CSV to use for ID-pair order.\n";
	$message .= "\t\t-I|--use-index\t\t\tUse numerical indices rather than the IDs.\n";
	$message .= "\t\t-N|--dist-name <NAME>\t\tRename the distance header.\n";
	$message .= "\t\t-S|--set-name\t\t\tAdd a column for the set name.\n";
	$message .= "\t\t-V|--vsX <ID>\t\t\tOnly print versus ID.\n";
	$message .= "\t\t-L|--vs-list <ID1[,ID2,..]>\tOne column is used per list reference (sqm subset).\n";
	$message .= "\t\t-R|--remove-refs\t\tGiven above, remove reference rows.\n";
	die($message."\n");
}


if ( !defined($sampleSet) ) {
	$sampleSet = '';
} else {
	$sampleSet = ",". $sampleSet;
}

if ( defined($delim) ) {
	$d = $delim;
} else {
	$d = "\t";
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
close(IN); $N = $i;
if ( !defined($distName) ) { $distName = 'dist'; }

# process vs List
if ( defined($vsList) ) {
	%vsID = ();
	foreach $id ( split(',',$vsList) ) {
		if ( defined($indexByID{$id}) ) {
			$vsID{$id} = $indexByID{$id};
		} else {
			print STDERR "Invalid reference, $id. Skipping.\n";
		}
	}
	$vsList = 1;
} else {
	$vsList = 0;
}

if ( !$vsList ) {	
	if ( defined($vsX) ) {
		print "vs $vsX,$distName\n";
	} else {
		print "id1,id2,$distName\n";
	}
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
			print $id1,$d,$id2,$d,$distHash{$id1}[$indexByID{$id2}],$sampleSet,"\n";
		} elsif( $id1 eq $vsX ) {
			print $id2,$d,$distHash{$id1}[$indexByID{$id2}],$sampleSet,"\n";
		} elsif( $id2 eq $vsX ) {
			print $id1,$d,$distHash{$id1}[$indexByID{$id2}],$sampleSet,"\n";
		}
	}
	close(IN);
} else {
	if ( $vsList ) {
		@columns = sort { $a <=> $b } values(%vsID);
		print 'ID'; 
		if ( $distName ) { $suffix = '_'.$distName; } else { $suffix = ''; }
		foreach $col ( @columns ) { print $d,$IDs[$col],$suffix; } 
		print "\n";
		for($i=0;$i<$N;$i++) {
			if ( $removeRefs && defined($vsID{$IDs[$i]}) ) { next; }

			print $IDs[$i];
			foreach $col ( @columns ) {
				print $d,$distHash{$IDs[$i]}[$col];
			}
			print "\n";
		}
	} elsif ( $fullMatrix ) {
		for($i=0;$i<$N;$i++) {
			for($j=0;$j<$N;$j++) {
				if ( !defined($vsX ) ) {
					if ( $useIndex ) {
						print $i,$d,$j,$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],$d,$IDs[$j],$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$i] eq $vsX ) {
					if ( $useIndex ) {
						print $j,$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$j],$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$j] eq $vsX ) {
					if ( $useIndex ) {
						print $i,$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				}
			}
		}
	} else {
		for($i=1;$i<$N;$i++) {
			for($j=0;$j<$i;$j++) {
				if ( !defined($vsX ) ) {
					if ( $useIndex ) {
						print $i,$d,$j,$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],$d,$IDs[$j],$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$i] eq $vsX ) {
					if ( $useIndex ) {
						print $j,$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$j],$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				} elsif( $IDs[$j] eq $vsX ) {
					if ( $useIndex ) {
						print $i,$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					} else {
						print $IDs[$i],$d,$distHash{$IDs[$i]}[$j],$sampleSet,"\n";
					}
				}
			}
		}
	}
}
