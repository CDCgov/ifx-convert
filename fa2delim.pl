#!/usr/bin/env perl
# Converts fasta to a delimited format


use Getopt::Long;
GetOptions( 
		'csv-format|C' => \$CSV,
		'per-site|P' => \$perSite,
		'enclose|E' => \$enclose,
		'delimiter|D=s' => \$delim,
		'extra-name|N=s' => \$name,
		'seq-delim|S=s' => \$sDelim,
		'codon-triplets|T' => \$triplets,
		'header-delim|H=s' => \$headerDelim
	);

if ( -t STDIN && ! scalar(@ARGV) ) {
	$message = "Usage:\n\tperl $0 <annotated.fasta> [options]\n";
	$message .= "\t\t-P|--per-site\t\tOne row per site.\n";
	$message .= "\t\t-E|--enclose\t\tEnclose fields in quotes.\n";
	$message .= "\t\t-C|--csv-format\t\tUse csv delimiter (tab default).\n";
	$message .= "\t\t-D|--delimiter <CHAR>\tDelimiter for output.\n";
	$message .= "\t\t-N|--extra-name <STR>\tExtra name field to add to every row.\n";
	$message .= "\t\t-S|--seq-delim <CHAR>\tCharacter delimiter for the string.\n";
	$message .= "\t\t-T|--codon-triplets\tAssumes data is in triplets for (-P and -S options).\n";
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

	if ( $headerDelim ) {
		@fields = split(/\Q$headerDelim\E/,$header);
		for($i=0;$i<scalar(@fields);$i++ ) {
			$fields[$i] = $q.$fields[$i].$q;
		}
		$header = join($delim,@fields);
	}

	if ( $perSite ) {
		if ( $triplets ) {
			for( $pos=0;$pos<$length;$pos += 3 ) {
				print $header,$extraField,$delim,$q,(int($pos/3)+1),$q,$delim,$q,substr($sequence,$pos,3),$q,"\n";
			}
		} else {
			for( $pos=0;$pos<$length;$pos++ ) {
				print $header,$extraField,$delim,$q,($pos+1),$q,$delim,$q,substr($sequence,$pos,1),$q,"\n";
			}
		}
	} elsif ( $singleLine ) {
		if ( $triplets ) {
			$sequence = join($sDelim, ($sequence =~ /.{3}/g) );
		} else{
			$sequence = join($sDelim,split('',$sequence));
		}
		print $header,$extraField,$delim,$q,$sequence,$q,"\n";
	} else {
		print $header,$extraField,$delim,$q,$sequence,$q,"\n";
	}
}

# Trim function.
# # Removes whitespace from the start and end of the string
 sub trim($) {
 	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
 	return $1;
}
