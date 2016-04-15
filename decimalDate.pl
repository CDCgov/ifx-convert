#!/usr/bin/perl
# decimalDate.pl - Sam Shepard - 2011.26.08
# version 1.1
# Converts normal dates into the decimal format.

use Getopt::Long;
GetOptions(
		'delimiter|D:s'=> \$delim,
		'field|F:i' => \$field,
		'ascending|A' => \$ascending,
		'us-format|U' => \$usFormat
	);


if (  scalar(@ARGV) != 2  ) {
	$message = "\n$0 <input.fasta> <output.fasta> [OPTIONS]\n";
	$message .= "\t-A|--ascending\t\t\t\tAscending date: DD/MM/YYYY. Default is YYYY/MM/DD.\n";
	$message .= "\t-U|--us-format\t\t\t\tUS format date: MM/DD/YYYY. Default is YYYY/MM/DD.\n";
	$message .= "\t-D|--delimiter CHARACTER\t\tSingle character delimiter within FASTA header. Default is '|'.\n";
	$message .= "\t-F|--field POSITIVE_NUMBER\t\tHeader field to use for dates. Default is 1.\n\n";
	die($message);
}


if ( !defined($field) ) {
	$field = 0;
} elsif($field == 0 ) {
	die("ERROR: field must be specified.\n");
} elsif($field < 0) {
	die("ERROR: field must be a positive number.\n");
} else {
	$field -= 1;
}

if ( !defined($delim) ) {
	$delim = '|';
} elsif( $delim eq '' ) {
	die("ERROR: No delimiter argument detected.\n");
} elsif( length($delim) > 1 ) {
	die("ERROR: single character delimiter expected instead of '$delim'.\n");
}

open( IN, '<', $ARGV[0] ) or die("ERROR: Could not open $ARGV[0].\n");
open( OUT, '>', $ARGV[1] ) or die("ERROR: Could not open $ARGV[1].\n");
$/ = ">"; $i = 0;
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));
	$length = length($sequence);

	@values = split(/\Q$delim\E/,$header);
	if ( $length == 0 ) {
		next;	
	} elsif ( $header !~ /\Q$delim\E/ && $delim ne '' ) {
		die("ERROR: delimiter '$delim' not found in header:\n\t$header\n");
	} elsif ( $#values < $field && 0 ) {
		$x = $#values + 1; $field++;
		die("ERROR: Only $x fields found in header but field $field selected:\n\t$header\n");
	}

	if ( $values[$field] eq '' ) {
		print "WARNING: No date information found. Skipping:\n\t$header\n\n";
		$numberSkipped++;
		next;
	} else {
		$values[$field] = decimalDate($values[$field]);
	}
	print OUT '>',join($delim,@values),"\n",$sequence,"\n";
}

close(IN);
close(OUT);
sub decimalDate($) {
	# Sam Shepard - 6.2011 - decimalDate
	# Converts fasta file with header format: >date_virus
	# To the decimal date format used by beast.
	# Ported from the decimalDate.php by Justin Bahl & Yongmei Liu.
	# Updated to handle empty values found within a database. -SS
	my $date = $_[0];
	my @days1 = (0,31,28,31,30,31,30,31,31,30,31,30,31);
	my @days2 = (0,31,29,31,30,31,30,31,31,30,31,30,31);
	my $dds = 0; 
	my $ddfrac = 0;
	my $i = 0;
	my ($yy, $mm, $dd);

	my @pieces = split(/\/|-|\./, $date);
	my $numPieces = scalar(@pieces);
	if ( $numPieces == 3 ) {
		if ( $ascending ) {
			($dd, $mm, $yy ) = @pieces;
		} elsif ( $usFormat ) {
			($mm, $dd, $yy ) = @pieces;
		} else {
			($yy, $mm, $dd ) = @pieces;
		}

		if ( $dd == 0 ) {
			$dd = 15;
		}
		if ( $mm == 0 ) {
			$mm = 6;
		}
	} elsif ( $numPieces == 2 ) {
		if ( $ascending || $usFormat ) {
			($mm, $yy) = @pieces;
		} else {
			($yy, $mm) = @pieces;
		}
		$dd = 15;
	} else {
		$yy = $pieces[0];
		($mm, $dd) = (7, 1);
	}

	if ( $yy % 4 == 0 ) {
		$dds = 0;
		for ($i = 1; $i < $mm; $i++) {
			$dds = $dds + $days2[$i];
		}
		$dds = $dds + $dd;
		$ddfrac = $dds / 366.0;
	} else {
		$dds = 0;
		for ($i = 1; $i < $mm; $i++ ) { 
			$dds = $dds + $days1[$i];
		}
		$dds = $dds + $dd;
		$ddfrac = $dds / 365.0;
	}

	if ($dd == 31 && $mm == 12) {
		$ddfrac = 0.999;
	}
	
	$ddfrac = sprintf("%.3f", $ddfrac);
	return sprintf('%s.%s', $yy, substr($ddfrac,2) );
}
