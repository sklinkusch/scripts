#!/usr/bin/perl

use strict;
use warnings;

my $file = join('', $ARGV[0]);
my $debtoau = 2.5417464;

open my $log, $file or die "Count not open $file: $!";

while (my $line = <$log>) {
    if ($line =~ /GROUND STATE \(SCF\) DIPOLE/) {
	my @gsdipxyz = split(/[ ]+/, parseFloat($line));
	foreach my $x (@gsdipxyz) { $x = $x / $debtoau; }
	print $gsdipxyz[0], " ",  $gsdipxyz[1], " ",  $gsdipxyz[2], "\n";
    }
    if ($line =~ /TRANSITION DIPOLE.*BOHR/) {
	my @trdipxyz = split(/[ ]+/, parseFloat($line));
	print $trdipxyz[0], " ", $trdipxyz[1], " ", $trdipxyz[2], "\n";
    }
    if ($line =~ /STATE DIPOLE.*BOHR/) {
	my @stdipxyz = split(/[ ]+/, parseFloat($line));
        print "\n", $stdipxyz[0], " ", $stdipxyz[1], " ", $stdipxyz[2], "\n ";
    }
}

close $log;

sub parseFloat {
    my $str = shift;
    $str =~ /([-]*[\d\.]+[ ]*[-]*[\d\.]+[ ]*[-]*[\d\.]+)/;
    return $1;
}
