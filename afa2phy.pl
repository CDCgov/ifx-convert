#!/usr/bin/env perl
# afa2phy.pl - Samuel S. Shepad - 1.2013
# Convert aligned fasta to phylip sequential format.
# Version 0.1.0

if (  scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\t$0 <aligned.fasta>\n";
	die($message);
}

# process parameters
open( IN, '<', $ARGV[0] ) or die("ERROR: Could not open $ARGV[0].\n");

@headers = (); @sequences = ();
$/ = ">"; $i = 0;
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$headers[$i] = shift(@lines);
	$sequences[$i] = uc(join('',@lines));

	$headers[$i] =~ s/^\s*(.*?)\s*$/\1/;
	$headers[$i] =~ s/[\s:]/_/g;
	$headers[$i] =~ s/[\(;]/-/g;
	$headers[$i] =~ tr/;',\)\[\]//d;

	$length = length($sequences[$i]);
	if ( $length == 0 ) {
		next;	
	}

	$i++;
}
close(IN);
$N = $i;
$L = length($sequences[0]);

print $N,' ',$L,"\n";
for($i = 0; $i < $N; $i++) {
	print $headers[$i],' ',$sequences[$i],"\n";
}
