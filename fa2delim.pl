#!/usr/bin/env perl
# Converts fasta to a delimited format


use Getopt::Long;
GetOptions( 
		'csv-format|C' => \$CSV,
		'per-site|P' => \$perSite,
		'enclose|E' => \$enclose,
		'delimiter|D=s' => \$delim,
		'extra-name|N=s' => \$name,
		'single-line-delimiter|S=s' => \$sDelim
	);

if ( -t STDIN && ! scalar(@ARGV) ) {
	$message = "Usage:\n\tperl $0 <annotated.fasta> [options]\n";
	$message .= "\t\t-P|--per-site\t\tOne row per site.\n";
	$message .= "\t\t-E|--enclose\t\tEnclose fields in quotes.\n";
	$message .= "\t\t-C|--csv-format\t\tUse csv delimiter (tab default).\n";
	$message .= "\t\t-D|--delimiter <CHAR>\tDelimiter for output.\n";
	$message .= "\t\t-N|--extra-name <STR>\tExtra name field to add to every row.\n";
	die($message."\n");
}

if ( defined($sDelim) ) {
	$singleLine = 1;
} else {
	$singleLine = 0;
}

if ( defined($CSV) ) {
	$delim = ',';
	$enclose = 1;
} elsif ( !defined($delim) ) {
	$delim = "\t";
	$enclose = 0;
}

if ( $enclose ) {
	$q = "'";
} else {
	$q = '';
}

if ( !defined($name) ) {
	$extraField = "";
} else {
	$extraField = $delim.$q.$name.$q;
}

$/ = ">";
while( $record = <> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = trim(shift(@lines));
	$sequence = uc(join('',@lines));

	$length = length($sequence);
	if ( $length == 0 ) {
		next;	
	}

	if ( $enclose ) {
		$header =~ tr/'/\'/;
		$sequence =~ tr/'//d;
	}

	$length = length($sequence);

	if ( $perSite ) {
		for( $pos=0;$pos<$length;$pos++ ) {
			print $q,$header,$q,$extraField,$delim,$q,($pos+1),$q,$delim,$q,substr($sequence,$pos,1),$q,"\n";
		}
	} elsif ( $singleLine ) {
		print $q,$header,$extraField,$q,$delim,$q,join($sDelim,split('',$sequence)),$q,"\n";
	} else {
		print $q,$header,$extraField,$q,$delim,$q,$sequence,$q,"\n";
	}
}

# Trim function.
# # Removes whitespace from the start and end of the string
 sub trim($) {
 	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
 	return $1;
}
