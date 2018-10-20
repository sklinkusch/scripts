#!/usr/bin/perl -w

use strict;
use warnings;

my $datfile = join('',$ARGV[0]);
my $popfile = join('',$ARGV[1]);
my $number = join('',$ARGV[2]);
my $targeten;
my %dips;
my %ens;
my %pops;
my $probability;
my %prob;
my %newprob;
my $totprob;

open my $dat, $datfile or die "Could not open $datfile: $!\n";
while(my $line = <$dat>){
    $line =~ s/->/ /;
    $line =~ s/^[ \t]+//;
    $line =~ s/\n//;
    my @string = split(/[ ]+/,$line);
    if($string[1] == $number){
	my $totdip = sqrt(($string[2]*$string[2])+($string[3]*$string[3])+($string[4]*$string[4]));
	$dips{$string[0]} = $totdip;
    }
}
close $dat;

open my $pop, $popfile or die "Could not open $popfile: $!\n";
while(my $line = <$pop>){
    $line =~ s/^[ \t]+//;
    $line =~ s/\n//;
    my @string = split(/[ \t]+/,$line);
    $ens{$string[0]} = $string[1];
    $pops{$string[0]} = $string[2];
}
close $pop;

$targeten = $ens{$number};
# remove states with negative populations
foreach my $key (keys %pops){
    if($pops{$key} <= 0){
	delete $pops{$key};
	delete $ens{$key};
	delete $dips{$key};
    }
}

# remove states with higher energies than the target state
foreach my $key (keys %ens){
    if($ens{$key} >= $targeten or $ens{$key} == 0.){
	delete $ens{$key};
	delete $pops{$key};
	delete $dips{$key};
    }
}

# remove states with zero dipole moment to target state
foreach my $key (keys %ens){
    if(!defined $dips{$key}){
	delete $pops{$key};
	delete $ens{$key};
    }
}

# compute probability to excite system from state |i> to target state
$totprob = 0.;
foreach my $key (keys %dips){
    $probability = $dips{$key} * $pops{$key};
    $prob{$key} = $probability;
    $totprob += $probability;
}

# divide probability by total probability
foreach my $key (keys %prob){
    $newprob{$key} = $prob{$key}/$totprob;
}

foreach my $key (sort {$newprob{$b} <=> $newprob{$a}} keys %newprob){
    printf("%3u   %20.18f   %20.18f\n",$key, $prob{$key}, $newprob{$key});
}
