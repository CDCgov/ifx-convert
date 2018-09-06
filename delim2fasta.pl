#!/usr/bin/env perl

use Getopt::Long;
GetOptions( 'input-delimiter|I=s' => \$inD,
	    'output-delimiter|O=s' => \$outD,
	    'strip-gap-chars|G' => \$stripGaps,
	    'jump-header-line|J' => \$jumpHeader 
	    );

if ( -t STDIN && ! scalar(@ARGV) ) {
	$message = "Usage:\n\tperl $0 <delimited.txt> [options]\n";
	$message .= "\t\t-I|--input-delimiter <CHAR>\tParsing delimiter for input, default <tab>.\n";
	$message .= "\t\t-J|--jump-header-line\t\tSkip first line when processing.\n";
	$message .= "\t\t-O|--output-delimiter <CHAR>\tParsing delimiter for output header, default is pipe '|'.\n";
	$message .= "\t\t-G|--strip-gap-chars\t\tStrip gap characters found in sequences: '.-~:'\n";
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
@lines = split(/\r\n|\n|\r/, <>); chomp(@lines);
if ( scalar(@lines) == 0 ) { exit(0); }
if ( $jumpHeader ) { $junk = shift(@lines); }
$firstLine = $lines[0];
@fields = split($inD,$firstLine);
$maxLen = 0; $seqField  = 0;
foreach $i ( 0 .. $#fields ) {
	$tmp = $fields[$i];
	$tmp =~ tr/ //d;
	if ( $tmp =~ /^[A-Za-z*.~:-]+$/ ) {
		if ( length($tmp) > $maxLen ) {
			$seqField = $i;
			$maxLen = length($tmp);
		}
	}
}

# loop over the rest
foreach $line ( @lines ) {
	@fields = split($inD,$line);
	$seq = $fields[$seqField]; $header = '';
	$seq =~ tr/ //d;
	if ( $stripGaps ) {
		$seq =~ tr/.:~-//d;
	}
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
