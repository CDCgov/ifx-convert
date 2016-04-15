#!/usr/bin/env perl

use Getopt::Long;
GetOptions( 'input-delimiter|I=s' => \$inD,
	    'output-delimiter|O=s' => \$outD );

if ( -t STDIN && ! scalar(@ARGV) ) {
	$message = "Usage:\n\tperl $0 <delimited.txt> [options]\n";
	$message .= "\t\t-I|--input-delimiter <CHAR>\tParsing delimiter for input, default <tab>.\n";
	$message .= "\t\t-O|--output-delimiter <CHAR>\tParsing delimiter for output header, default is pipe '|'.\n";
	die($message."\n");
}

if ( !defined($inD) ) {
	$inD = "\t";
}

if ( !defined($outD) ) {
	$outD = '|';
}

# discover sequence field
$/ = undef;
@lines = split(/\r\n|\n|\r/, <>);
chomp(@lines);
$firstLine = shift(@lines);
@fields = split($inD,$firstLine);
$maxLen = 0; $seqField  = 0;
foreach $i ( 0 .. $#fields ) {
	$tmp = $fields[$i];
	$tmp =~ tr/ //d;
	if ( $tmp =~ /^[A-Za-z.~-]+$/ ) {
		if ( length($tmp) > $maxLen ) {
			$seqField = $i;
			$maxLen = length($tmp);
		}
	}
}


# output first line
$seq = $fields[$seqField]; $header = '';
foreach $i ( 0 .. $#fields ) {
	if ( $i != $seqField ) {
		if ( $header ne '' ) {
			$header .= $outD.$fields[$i];
		} else {
			$header = $fields[$i];
		}
	}
}
$seq =~ tr/ //d;
if ( substr($header,0,1) ne '>' ) { print '>'; }
print $header,"\n",$seq,"\n";

# loop over the rest
foreach $line ( @lines ) {
	@fields = split($inD,$line);
	$seq = $fields[$seqField]; $header = '';
	$seq =~ tr/ //d;
	foreach $i ( 0 .. $#fields ) {
		if ( $i != $seqField ) {
			if ( $header ne '' ) {
				$header .= $outD.$fields[$i];
			} else {
				$header = $fields[$i];
			}
		}
	}

	if ( substr($header,0,1) ne '>' ) { print '>'; }
	print $header,"\n",$seq,"\n";
}
