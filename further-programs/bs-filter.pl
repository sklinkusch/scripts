#!/usr/bin/perl -w

# a perl script reading the output from browser-sync and swallowing it after the first rows

use strict;
use warnings;

my $row = 0;
while (my $line = <>){
 print $line if ($row < 8);
 $row++;
}
