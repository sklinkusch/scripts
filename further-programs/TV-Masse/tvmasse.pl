#!/usr/bin/perl -w

use strict;
use warnings;

print_error() if ($#ARGV != 1);
my $diag = join('',$ARGV[0]);
my $ratio = join('',$ARGV[1]);
my $inchtocm = 2.54;
my $wpl;
if ( $ratio eq 'A' or $ratio eq 'a'){
  $wpl = 9/16;
}elsif( $ratio eq 'B' or $ratio eq 'b'){
  $wpl = 3/4;
}else{
  print_error();
}
if($ratio eq 'a' or $ratio eq 'b'){
 $diag *= $inchtocm;
}
my $wplsq = $wpl**2;
my $corsq = 1 + $wplsq;

my $diagsq = $diag**2;
my $lengsq = $diagsq / $corsq;
my $length = sqrt($lengsq);
my $width  = $wpl * $length;
my $lengthin = $length/$inchtocm;
my $widthin  = $width/$inchtocm;
my $diagin   = $diag/$inchtocm;
printf("Width:   %5.1f cm  %5.1f inch\n",$length, $lengthin);
printf("Height:  %5.1f cm  %5.1f inch\n",$width, $widthin);
printf("Diagonal:%5.1f cm  %5.1f inch\n",$diag, $diagin);

sub print_error {
 print "Usage: ./tvmasse.pl <diagonal> <aspect ratio>\n";
 print "<aspect ratio> = A, a, B, or b\n";
 print "A: 16:9, diagonal in cm\n";
 print "a: 16:9, diagonal in inch\n";
 print "B:  4:3, diagonal in cm\n";
 print "B:  4:3, diagonal in inch\n";
 exit;
}
