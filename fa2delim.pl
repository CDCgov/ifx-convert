#!/usr/bin/env perl
# Converts fasta to a delimited format


use Getopt::Long;
GetOptions( 
		'csv-format|C' => \$CSV
		);

if ( -t STDIN && ! scalar(@ARGV) ) {
	$message = "Usage:\n\tperl $0 <annotated.fasta> [options]\n";
	$message = "\t\t-C|--csv-format\tUse csv delimiter (tab default).\n";
	die($message."\n");
}


if ( defined($CSV) ) {
	$delim = ',';
	$enclose = 1;
}

if ( !defined($delim) ) {
	$delim = "\t";
	$enclose = 0;
}

if ( $enclose ) {
	$q = "'";
} else {
	$q = '';
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

	print $q,$header,$q,$delim,$q,$sequence,$q,"\n";
}

# Trim function.
# # Removes whitespace from the start and end of the string
 sub trim($) {
 	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
 	return $1;
}

