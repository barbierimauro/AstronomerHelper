#!/usr/bin/perl

# Produce ut1.dat for Tempo from the new IERS Bulletin B format

# Process many bulletin B input files
$nfile = scalar(@ARGV);

if ($nfile==0) {
    print STDERR "Usage: do.iers.ut1.new bulletinb.XXX (bulletinb.YYY ...)\n";
    print STDERR "\n";
    print STDERR "  Many bulletin B files can be processed at once, however\n";
    print STDERR "  they must be a continuous set.  This program only deals\n";
    print STDERR "  with \"new-format\" bulletin B files\n";
    print STDERR "\n";
    exit(1);
}


$mjd0 = 0;
$mjd1 = 0;
%tai_ut1 = ();
for $i (0..$nfile-1) {
    $file = $ARGV[$i];
    # Assume last file given is the one we need extension values from
    if ($i==$nfile-1) { $do_ext=1; } else { $do_ext=0; }
    # print STDERR "Processing $file\n";
    open (BULB, $file) or die "Error opening $file\n";
    $read=0;
    %ut1_utc = ();
    while (<BULB>) {
        # This loop gets the UT1-UTC values, including preliminary
        # extensions for the last file in the set.
        chomp;
        # Detect old-format files
        if (/1 - EARTH ORIENTATION PARAMETERS/) { 
            die "Error: $file appears to be an \"old-format\" bulletin B\n";
        }
        if (/Mean formal error/) { $read=1; next; }
        if ($read==0) { next; }
        if (/^\s+$/) { next; }
        if ($do_ext==0 && /Preliminary extension/) { last; }
        if ($do_ext==1 && /DAILY/) { last; }
        @line = split;
        if (scalar(@line)!=14) { next; }
        $mjd = $line[3];
        $dut1 = $line[6] / 1000.0; # Convert units to seconds
        if ($mjd0==0 || $mjd<$mjd0) { $mjd0=$mjd; }
        if ($mjd1==0 || $mjd>$mjd1) { $mjd1=$mjd; }
        #print "$mjd $dut1\n"; # DEBUG
        $ut1_utc{$mjd} = $dut1;
    }
    $tai_utc = 0;
    while (<BULB>) {
        # This loop looks for the TAI-UTC difference.  We'll need
        # to do more parsing if there has been a leap second in
        # the relevant date range.  So far this has not yet happened
        # during the existence of "new-format" files, so I don't
        # know the correct syntax.
        chomp;
        if (/TAI\s+-\s+UTC\s+=\s+(\d+)\s+s/) { $tai_utc = $1; }
    }
    #print "TAI - UTC = $tai_utc\n"; # DEBUG
    close BULB;

    # Update the final TAI - UT1 table
    foreach $mjd (keys %ut1_utc) {
        $tai_ut1{$mjd} = $tai_utc - $ut1_utc{$mjd};
    }
}
#print "$mjd0 $mjd1\n";

# Now the values need to be reformatted into ut1.dat format
$ref_mjd = 40219; # First MJD in the ut1.dat file - maybe read it from the file?
%lines = ();
for ($mjd=$mjd0; $mjd<$mjd1; $mjd++) {

    # We only need values every 5 days
    unless ((($mjd - $ref_mjd) % 5) == 0) { next; }

    unless (exists $tai_ut1{$mjd}) { die "Error: missing data for $mjd\n"; }
    #print "$mjd $tai_ut1{$mjd}\n";

    # Start MJD for the line that this date belongs in
    $line_mjd = int(($mjd - $ref_mjd)/30.0)*30.0 + $ref_mjd;
    # Position in line
    $line_pos = int(($mjd - $line_mjd)/5);

    #print "$line_mjd $line_pos\n";

    # Skip any incomplete first lines
    if ($line_pos!=0 and !exists $lines{$line_mjd}) { next; }

    # Round value to int
    $ut1val = int($tai_ut1{$mjd}*1e4 + 0.5);

    # Create line and first entry
    if ($line_pos==0) {
        $lines{$line_mjd} = sprintf("%10d%7d", $line_mjd, $ut1val);
    } else {
        $lines{$line_mjd} .= sprintf(" %7d", $ut1val);
    }
}

# Add final line notation
for ($i=$line_pos+1; $i<6; $i++) {
    $lines{$line_mjd} .= sprintf("        ");
}
$lines{$line_mjd} .= sprintf("  %02d", $line_pos + 1);

# Print output
# TODO: make this automatically update ut1.dat like David's script?
foreach $mjd (sort keys %lines) { print $lines{$mjd} . "\n"; }

