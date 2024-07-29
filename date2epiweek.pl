#!/usr/bin/perl
# Samuel Shepard - 7.2011
# date2Epiweek.pl
# Converts dates to Epiweek.

# Version 1.

use Getopt::Long;
GetOptions(
            'no-newline|N'   => \$noNewline,
            'year-format|Y'  => \$byYear,
            'pretty-print|P' => \$pretty
);

if ( scalar(@ARGV) < 1 ) {
    $message = "\n$0 <input_date> [OPTIONS]\n";
    $message .= "\t-N|--no-newline\t\t\tDo not print newline.\n";
    $message .= "\t-Y|--year-format\t\tOutput Epiweek in yearweek format.\n";
    $message .= "\t-P|--pretty-print\t\tOutputs a prettier format.\n";
    die($message);
}

$newline = "\n";
if ( defined($noNewline) ) {
    $newline = "";
}

%startDates = ();
$date       = $ARGV[0];
$date       = weekDate( $date, \%startDates );

if ( defined($pretty) ) {
    $week = substr( $date, -2 );
    $year = substr( $date, 0, 4 );
    print sprintf( "%4d Wk. %2d", $year, $week ), $newline;
} elsif ( defined($byYear) ) {
    print $date, $newline;
} else {
    print substr( $date, -2 ), $newline;
}

sub startDate($) {
    my $qYear = $_[0];
    my ( $eDay, $eMonth, $eYear, $sDay, $sMonth, $sYear, $y );

    $eDay = 9;
    for ( $y = 1972; $y <= $qYear; $y++ ) {
        if ( ( $y - 1 ) % 4 == 0 ) {
            $eDay -= 2;
        } else {
            $eDay--;
        }

        if ( $eDay < 4 ) {
            $eDay += 7;
        }

        $eMonth = 1;
        $eYear  = $y;
        $sDay   = $eDay - 6;
        if ( $sDay < 1 ) {
            $sDay   = 31 + $sDay;
            $sMonth = 12;
            $sYear  = $y - 1;
        } else {
            $sMonth = 1;
            $sYear  = $y;
        }
    }

    return ( $sYear, $sMonth, $sDay, $eYear, $eMonth, $eDay );
}

sub weekDate($$) {
    my @days1 = ( 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );
    my @days2 = ( 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );
    my $dds   = 0;

    # relies on startDate to compute the weekDate for a given date field
    my $date       = $_[0];
    my $startDates = $_[1];
    my ( $yy, $mm, $dd, $i, $tmp );
    my ( $dYear, $dWeek );

    my @pieces    = split( /\/|-/, $date );
    my $numPieces = scalar(@pieces);
    if ( $numPieces == 3 ) {
        ( $yy, $mm, $dd ) = @pieces;
    } elsif ( $numPieces == 2 ) {
        ( $yy, $mm ) = @pieces;
        $dd = 15;
    } else {
        $yy = $pieces[0];
        ( $mm, $dd ) = ( 6, 15 );
    }

    if ( !exists( $startDates->{$yy} ) ) {
        $startDates->{$yy} = [startDate($yy)];
    }

    $tmp = $yy - 1;
    if ( !exists( $startDates->{$tmp} ) && $yy > 1972 ) {
        $startDates->{$tmp} = [startDate($tmp)];

    }

    $tmp = $yy + 1;
    if ( !exists( $startDates->{$tmp} ) ) {
        $startDates->{$tmp} = [startDate($tmp)];
    }

    if ( $mm == 1 ) {
        if ( $yy == $startDates->{$yy}[0] && $mm == $startDates->{$yy}[1] && $dd < $startDates->{$yy}[2] ) {
            $yy--;
            if ( $yy % 4 == 0 ) {
                $dds = 366;
            } else {
                $dds = 365;
            }
        }
    } elsif ( $mm == 12 ) {

        # will create an array if doesn't exist by accessing element
        if ( $yy == $startDates->{ $yy + 1 }[0] && $mm == $startDates->{ $yy + 1 }[1] && $dd >= $startDates->{ $yy + 1 }[2] )
        {
            if ( $yy % 4 == 0 ) {
                $dds = -366;
            } else {
                $dds = -365;
            }
            $yy++;
        }
    }

    if ( $yy % 4 == 0 ) {
        for ( $i = 1; $i < $mm; $i++ ) {
            $dds = $dds + $days2[$i];
        }
        $dds = $dds + $dd;
    } else {
        for ( $i = 1; $i < $mm; $i++ ) {
            $dds = $dds + $days1[$i];
        }
        $dds = $dds + $dd;
    }

    if ( $startDates->{$yy}[0] < $yy ) {
        $dds = $dds - ( $startDates->{$yy}[2] - 32 ) - 1;
    } else {
        $dds = $dds - ( $startDates->{$yy}[2] - 1 ) - 1;
    }
    $dWeek = int( $dds / 7 ) + 1;
    $dYear = $yy;

    return sprintf( "%4d%02d", $dYear, $dWeek );
}
