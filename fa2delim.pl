#!/usr/bin/env perl
# Sam Shepard - 2016
# Converts fasta to a delimited format

use Digest::SHA qw(sha1_hex);
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use English qw(-no_match_vars);
use Carp qw(croak);
use warnings;

GetOptions(
            'csv-format|C'      => \$CSV,
            'per-site|P'        => \$perSite,
            'enclose|E'         => \$enclose,
            'delimiter|D=s'     => \$delim,
            'extra-name|N=s'    => \$name,
            'seq-delim|S=s'     => \$sDelim,
            'codon-triplets|T'  => \$triplets,
            'header-delim|H=s'  => \$headerDelim,
            'add-length|L'      => \$addLength,
            'add-size|Z'        => \$addSize,
            'add-hash|A'        => \$addHash,
            'inserts|I=s'       => \$insertsFile,
            'use-original|O'    => \$useOriginal,
            'use-unaligned|U'   => \$useUnaligned,
            'use-both|B'        => \$useBoth,
            'show-if-insertion' => \$insertionApplied,
            'show-shift-indel'  => \$frameShifted,
            'nt-id'             => \$doNtID,
            'fna-id'            => \$fnaID
);

if ( -t STDIN && !scalar(@ARGV) ) {
    $message = "Usage:\n\tperl $PROGRAM_NAME <annotated.fasta> [options]\n";
    $message .= "\t\t-P|--per-site\t\tOne row per site.\n";
    $message .= "\t\t-E|--enclose\t\tEnclose fields in quotes.\n";
    $message .= "\t\t-C|--csv-format\t\tUse csv delimiter (tab default).\n";
    $message .= "\t\t-D|--delimiter <CHAR>\tDelimiter for output.\n";
    $message .= "\t\t-N|--extra-name <STR>\tExtra name field to add to every row.\n";
    $message .= "\t\t-Z|--add-size\t\tAdd field for ungapped size (unaligned + ins).\n";
    $message .= "\t\t-A|--add-hash\t\tAdd field for sequence hash.\n";
    $message .= "\t\t-L|--add-length\t\tAdd seq length field (aligned + ins).\n";
    $message .= "\t\t-S|--seq-delim <CHAR>\tCharacter delimiter for the string.\n";
    $message .= "\t\t-I|--inserts <FILE>\tTab delimited insertion file: ID<t>NTS<T>AA\n";
    $message .= "\t\t-U|--use-unaligned\tPrint unaligned sequence/sites: unaligned + ins. Def: aligned + ins\n";
    $message .= "\t\t-O|--use-original\tPrint aligned sequence/sites: aligned - ins. Def: aligned + ins\n";
    $message .= "\t\t-B|--use-both\t\tPrints unaligned and the original (aligned) sequences.\n";
    $message .= "\t\t-T|--codon-triplets\tAssumes data is in triplets for (-P and -S options).\n";
    $message .= "\t\t   --show-if-insertion\tBoolean (true/false) added if insertion was put into the sequence.\n";
    $message .= "\t\t   --show-shift-indel\tBoolean (true/false) added if an indel not divisible by 3.\n";
    $message .= "\t\t   --nt-id\t\tPerform the nt_id instead of the variant_hash, given --add-hash option.\n";
    $message .= "\t\t   --fna-id\t\tGets the ID portion of the fasta header for NCBI fna files.\n";
    die( $message . "\n" );
}

# Take the nucleotide ID as reckoned by PubSeq
sub nt_id2 {
    my $seq = uc( shift // q{} );
    $seq =~ tr/\n\r\t :.~-//d;
    if ( $seq ne q{} ) {
        return ( sha1_hex($seq), $seq );
    } else {
        return ( '\N', '\N' );
    }
}

# Two argument version
sub variant_hash2 {
    my $seq = uc( shift // q{} );
    $seq =~ tr/\n\r\t :.-//d;
    if ( $seq ne q{} ) {
        return ( md5_hex($seq), $seq );
    } else {
        return ( '\N', '\N' );
    }
}

# Trim function.
# Removes whitespace from the start and end of the string
sub trim {
    my $string = shift // q{};
    if ( $string =~ /^\s*(.*?)\s*$/smx ) {
        return $1;
    } else {
        return $string;
    }
}

if ( defined $doNtID ) {
    *doHash = \&nt_id2;
} else {
    *doHash = \&variant_hash2;
}

if ( defined $sDelim ) {
    $singleLine = 1;
} else {
    $singleLine = 0;
}

if ( defined $CSV ) {
    $delim   = ',';
    $enclose = 1;
} elsif ( !defined $delim ) {
    $delim = "\t";
}

if ($enclose) {
    $q = "'";
} else {
    $q = q{};
}

if ( defined $useBoth ) {
    $useBoth      = 1;
    $useUnaligned = 1;
    $useOriginal  = 0;
} else {
    $useBoth = 0;
}

if ( !defined $name ) {
    $extraField = q{};
} else {
    $extraField = $delim . $q . $name . $q;
}

if ( defined $insertsFile ) {
    %inserts = ();
    local $RS = "\n";

    open( my $INS, '<', $insertsFile ) or die("Cannot open $insertsFile for reading.\n");
    @lines = <$INS>;
    chomp(@lines);
    foreach $line (@lines) {
        ( $id, $pos, $insert ) = split( "\t", $line );

        # position is upstream base for index-1
        # ergo, no adjustment is needed for index-0
        $inserts{$id}{$pos} = $insert;
    }
    close($INS) or croak("Cannot close file: $OS_ERROR");
    $tryInsertions = 1;
} else {
    $tryInsertions = 0;
}

local $RS = ">";
my ( $seqExtra, $lengthField, $hashField, $hasInsertion, $hasShifted ) = ( q{}, q{}, q{}, q{}, q{} );
while ( $record = <> ) {
    chomp($record);
    @lines            = split( /\r\n|\n|\r/smx, $record );
    $id               = $header   = trim( shift(@lines) );
    $sequenceOriginal = $sequence = uc( join( q{}, @lines ) );

    if ( length($sequence) == 0 )             { next; }
    if ($enclose)                             { $header =~ tr/'/\'/; $sequence =~ tr/'//d; }
    if ( $fnaID && $header =~ /^(\S+)\s/smx ) { $header = $1; }

    if ($insertionApplied) { $hasInsertion = $delim . 'false'; }
    if ($frameShifted)     { $hasShifted   = $delim . 'false'; }
    if ( $tryInsertions && defined $inserts{$id} ) {
        $sequence =~ s/[.]+$//smx;
        $offset = 0;
        foreach my $pos ( sort { $a <=> $b } keys( %{ $inserts{$id} } ) ) {
            if ( length($sequence) >= ( int($pos) + $offset ) ) {
                $insert = $inserts{$id}{$pos};
                substr( $sequence, int($pos) + $offset, 0, $insert );
                $offset += length($insert);
                if ( defined $frameShifted && length($insert) % 3 != 0 ) { $hasShifted   = $delim . 'true'; }
                if ($insertionApplied)                                   { $hasInsertion = $delim . 'true'; }
            }
        }
    }

    if ($frameShifted) {
        while ( $sequenceOriginal =~ m/[^-](-+)[^-]/gsmx ) {
            if ( length($1) % 3 != 0 ) { $hasShifted = $delim . 'true'; }
        }
    }

    if ($headerDelim) {
        @fields = split( /\Q$headerDelim\E/smx, $header );
        foreach my $i ( 0 .. $#fields ) {
            $fields[$i] = $q . $fields[$i] . $q;
        }
        $header = join( $delim, @fields );
    }

    # length of the original sequence, plus any insertion
    if ($addLength) { $lengthField = $delim . length($sequence); }

    if ( $addSize || $addHash || $useUnaligned ) {
        ( $hash, $sequenceForHash ) = doHash($sequence);

        # the protein size does not include any alignment characters
        if ($addSize) { $sizeField = $delim . length($sequenceForHash); }
        if ($addHash) { $hashField = $delim . $q . $hash . $q; }
    }

    if ($useOriginal) {
        $sequence = $sequenceOriginal;
    } elsif ($useUnaligned) {
        $sequence = $sequenceForHash;
    }

    # the length for each codon and/or site can be unaligned (the one used in the Hash) or the original
    # the original sequence skips insertions, so they will not be printed herein
    $length = length($sequence);
    if ($perSite) {

        # one site/codon per record
        if ($triplets) {

            ## no critic (ControlStructures::ProhibitCStyleForLoops)
            for ( my $pos = 0; $pos < $length; $pos += 3 ) {
                print STDOUT $header, $extraField, $hashField // q{}, $lengthField, $sizeField // q{}, $hasInsertion,
                  $hasShifted, $delim, $q, ( int( $pos / 3 ) + 1 ), $q, $delim, $q, substr( $sequence, $pos, 3 ), $q, "\n";
            }
        } else {
            foreach my $pos ( 0 .. ( $length - 1 ) ) {
                print STDOUT $header, $extraField, $hashField // q{}, $lengthField, $sizeField // q{}, $hasInsertion,
                  $hasShifted, $delim, $q, ( $pos + 1 ), $q, $delim, $q, substr( $sequence, $pos, 1 ), $q, "\n";
            }
        }
    } elsif ($singleLine) {

        # one sequence per record
        if ($triplets) {
            $sequence = join( $sDelim, ( $sequence =~ /.{3}/gsmx ) );
        } else {
            $sequence = join( $sDelim, split( q{}, $sequence ) );
        }
        print STDOUT $header, $extraField, $hashField // q{}, $lengthField, $sizeField // q{}, $hasInsertion, $hasShifted,
          $delim, $q, $sequence, $q, "\n";
    } else {
        $extraSeq = $useBoth ? $delim . $q . $sequenceOriginal . $q : q{};
        print STDOUT $header, $extraField, $hashField // q{}, $lengthField, $sizeField // q{}, $hasInsertion, $hasShifted,
          $delim, $q, $sequence, $q, $extraSeq, "\n";
    }
}

