#!/usr/bin/env perl
# Sam Shepard - 2014
# Converts FASTA to FASTQ trivially.

$/ = ">";
while ( $record = <> ) {
    chomp($record);
    @lines  = split( /\r\n|\n|\r/, $record );
    $header = shift(@lines);
    $header =~ tr/ /_/;
    $sequence = lc( join( '', @lines ) );
    $sequence =~ tr/wsmkrybdhv/n/;

    $len = length($sequence);
    if ( $len > 0 ) {
        print '@', $header, "\n";
        print lc($sequence), "\n";

        print '+',        "\n";
        print '~' x $len, "\n";
    }
}
