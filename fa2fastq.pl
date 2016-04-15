#!/usr/bin/env perl

$/ = ">";
while($record = <>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$header =~ tr/ /_/;
	$sequence = join('',@lines);

	$len = length($sequence);
	if ( $len > 0 ) {
		print '@',$header,"\n";
		print lc($sequence),"\n";
		print '+',"\n";
		print '~'x$len,"\n";
	}
}
