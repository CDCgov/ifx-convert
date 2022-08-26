#!/usr/bin/env perl
## no critic (ValuesAndExpressions::ProhibitConstantPragma,ControlStructures::ProhibitPostfixControls)
# Sam Shepard

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use English qw(-no_match_vars);
use Carp qw(croak);

use constant { true => 1, false => 0 };

binmode( STDOUT, ":utf8" );

my $NULL           = '\N';
my $hasInputHeader = false;
my $noOutputHeader = false;
my $convertIfCSV   = false;
my $printIfGlobSet = false;
my $getFileField   = false;
my $noEmptyEnd     = false;
my $nullifyEnd     = false;
my $fileField      = 0;
my $fileDelim      = '_';
my $PROG           = "CollectXSV";

my @ext = ( '.txt', '.tsv', '.csv', '.TXT', '.TSV', '.CSV' );

sub uniq(@) {
    my @a = @_;
    my %h = ();
    return ( grep { !$h{$_}++ } @a );
}

GetOptions(
            'has-input-header|H'  => \$hasInputHeader,
            'no-output-header|N'  => \$noOutputHeader,
            'convert-if-csv|C'    => \$convertIfCSV,
            'file-field|F=i'      => \$getFileField,
            'file-delim|D=s'      => \$fileDelim,
            'print-if-glob-set|G' => \$printIfGlobSet,
            'no-empty-end|E'      => \$noEmptyEnd,
            'nullify-ends|Z'      => \$nullifyEnd
);

if ( scalar(@ARGV) < 1 ) {
    die(   "\nUsage:\n\t$PROGRAM_NAME\t<path> [OPTS]\n"
         . "\t\t-H|--has-input-header\n"
         . "\t\t-N|=-no-output-header\n"
         . "\t\t-C|--convert-if-csv\n"
         . "\t\t-G|--print-if-glob-set\n"
         . "\t\t-E|--no-empty-end\n"
         . "\t\t-Z|--nullify-ends\n"
         . "\n" );
}

if ( $getFileField > 0 ) {
    $fileField    = $getFileField - 1;
    $getFileField = true;
} else {
    $getFileField = false;
}

my $first = true;
if ($noOutputHeader) {
    $first = false;
}

if ( $nullifyEnd && $noEmptyEnd ) {
    die("$PROGRAM_NAME ERROR: -Z and -E options cannot both be used together\n");
}

my $printGlobSet = false;
my $globField    = 0;
if ( $printIfGlobSet && $ARGV[0] =~ /\{.+?\}/smx ) {
    my @parts = split( '/', $ARGV[0] );
    for my $i ( 0 .. $#parts ) {
        if ( $parts[$i] =~ /\{.+?\}/smx ) {
            $globField    = $i;
            $printGlobSet = true;
            last;
        }
    }
}

my @files = glob( $ARGV[0] );
foreach my $i ( 1 .. $#ARGV ) {
    push( @files, glob( $ARGV[$i] ) );
}
@files = uniq( sort(@files) );

my $d             = "\t";     # output delim
my $files_written = false;    # check if output was actually written
if ( scalar @files > 0 ) {
    for my $file (@files) {
        if ( !-R $file ) {
            print STDERR "$PROG WARNING: $file' not readable, skipping.\n";
            next;
        } elsif ( !-s $file ) {
            print STDERR "$PROG WARNING: '$file' was empty, skipping.\n";
            next;
        } else {
            print STDERR "$PROG INFO:    '$file' found.\n";
            $files_written = true;
            open( my $handle, '<:encoding(utf-8)', $file ) or die("Cannot open file $file.\n");
            my $convertCSV = $convertIfCSV && $file =~ /csv$/ismx ? true : false;

            my $globSetValue = $NULL;
            if ($printGlobSet) {
                my @parts = split( '/', $file );
                $globSetValue = $globField > $#parts ? $NULL : $parts[$globField];
            }

            my $fileFieldValue = $NULL;
            if ($getFileField) {
                my $basefile = basename( $file, @ext );
                my @parts    = split( /$fileDelim/smx, $basefile );
                $fileFieldValue = $fileField > $#parts ? $NULL : $parts[$fileField];
            }

            if ($hasInputHeader) {
                my $header = <$handle>;
                chomp($header);
                $header =~ tr/\r//d;
                $header =~ tr/,/\t/   if $convertCSV;
                $header =~ s/\t$//smx if $noEmptyEnd;

                if ($first) {
                    print STDOUT $header;
                    print STDOUT $d, 'Set'                   if $printGlobSet;
                    print STDOUT $d, 'F', ( $fileField + 1 ) if $getFileField;
                    print STDOUT "\n";
                    $first = false;
                }
            }
            while ( my $line = <$handle> ) {
                chomp($line);
                $line =~ tr/\r//d;
                $line =~ tr/,/\t/        if $convertCSV;
                $line =~ s/\t$//smx      if $noEmptyEnd;
                $line =~ s/\t$/$NULL/smx if $nullifyEnd;

                print STDOUT $line;
                print STDOUT $d, $globSetValue   if $printGlobSet;
                print STDOUT $d, $fileFieldValue if $getFileField;
                print STDOUT "\n";
            }
            close $handle or croak("Cannot close file: $OS_ERROR\n");
        }
    }
}

if ( !$files_written ) {
    die("$PROG ERROR:   NO data written!\n");
}
