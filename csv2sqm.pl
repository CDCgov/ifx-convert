#!/usr/bin/env perl

use Getopt::Long;
GetOptions( 'annotation-field|A' => \$annotField,
            'tab-delim|T'        => \$tabDelim );

if ( scalar(@ARGV) != 1 ) {
    die("Usage:\n\tperl $0 <2-col-matrix> [-A|--annotation-field] [-T|--tab-delim]\n");
}

if ($tabDelim) {
    $d = "\t";
} else {
    $d = ',';
}

%annots = %distMatrix = ();
open( IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[1].\n");
$header = <IN>;
while ( $line = <IN> ) {
    chomp($line);
    @values                                 = split( $d, $line );
    $distMatrix{ $values[0] }{ $values[1] } = $values[2];
    $distMatrix{ $values[1] }{ $values[0] } = $values[2];
    if ( $annotField && $values[0] =~ /{(.+?)}/ ) {
        $annots{ $values[0] } = $1;
    } else {
        $annots{ $values[0] } = 'unknown';
    }

    if ( $annotField && $values[1] =~ /{(.+?)}/ ) {
        $annots{ $values[1] } = $1;
    } else {
        $annots{ $values[1] } = 'unknown';
    }
}
close(IN);
$N   = scalar( keys(%annots) );
@IDs = sort( keys(%annots) );

for ( $i = 0; $i < $N; $i++ ) {
    $distMatrix{ $IDs[$i] }{ $IDs[$i] } = 0;
    print $IDs[$i];
    if ($annotField) {
        print "\t", $annots{ $IDs[$i] };
    }
    for ( $j = 0; $j < $N; $j++ ) {
        print "\t", $distMatrix{ $IDs[$i] }{ $IDs[$j] };
    }
    print "\n";
}

