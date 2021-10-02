#!/usr/bin/perl -w

use strict;
use warnings;
no warnings 'experimental';
use List::Util qw(reduce any all none notall first max maxstr min minstr product sum sum0 pairs unpairs pairkeys pairvalues pairfirst pairgrep pairmap shuffle uniq uniqnum uniqstr);
use POSIX qw/ceil/;
use Algorithm::Combinatorics qw(combinations);
use FindBin;

# Überprüfen, ob die Eingabe korrekt ist
print_info() if($#ARGV < 3 or $#ARGV % 2 == 0);

# Übergabe der Parameter
my $modus_a = $ARGV[0];
shift(@ARGV);
my $modus_b = $ARGV[0];
shift(@ARGV);
my @arguments = @ARGV;

# Einsortieren als Parteinamen und Sitzen
my @unsorted_keys = pairkeys @arguments;
my @unsorted_values = pairvalues @arguments;

# Bildung eines Hashes
my %results;
foreach my $x (0..$#unsorted_keys) {
	if($unsorted_values[$x] != 0){
		$results{$unsorted_keys[$x]} = $unsorted_values[$x];
	}
}

# Sortieren der Parteien nach Größe (absteigend)
my @keys;
my @values;
foreach my $key (sort { $results{$b} <=> $results{$a} } (keys(%results))) {
	push(@keys,$key);
	push(@values,$results{$key});
}

# Summe und Mehrheit
my $summe = sum @values;
print "Gesamt: $summe Mandate\n";
my $half  = $summe / 2;
my $absmh;
my $mdreg = ceil(0.4 * $summe);
my $percent;
if($summe % 2 == 0){
	$absmh = $half + 1;
} else {
	$absmh = ceil($half);
}
print "Absolute Mehrheit ab $absmh Sitzen\n";

# Ermitteln aller Koalitionen (mind. eine Partei, maximal alle Parteien)
my @combkeys;
my @combvals;
my @combnames;
foreach my $n (1..(scalar @keys)){
	my $iter = combinations(\@keys,$n);
	while (my $c = $iter->next) {
		my $comb = join("+",@$c);
		my $sum = summe(@$c);
		push(@combkeys,$comb);
		push(@combvals,$sum);
		push(@combnames,naming(@$c));
	}
}

# Ermitteln der erfolgreichen Koalitionen (d.h. mit absoluter Mehrheit)
my @succkeys;
my @succvals;
my @succnames;
if ($modus_a eq "A" or $modus_a eq "B"){
	foreach my $s (0..$#combkeys){
		push(@succkeys,$combkeys[$s]) if ($combvals[$s] >= $absmh);
		push(@succvals,$combvals[$s]) if ($combvals[$s] >= $absmh);
		push(@succnames,$combnames[$s]) if ($combvals[$s] >= $absmh);
	}
} elsif ($modus_a eq "C" or $modus_a eq "D"){
	foreach my $t (0..$#combkeys){
		push(@succkeys,$combkeys[$t]) if ($combvals[$t] >= $mdreg);
		push(@succvals,$combvals[$t]) if ($combvals[$t] >= $mdreg);
		push(@succnames,$combnames[$t]) if ($combvals[$t] >= $mdreg);
	}
} else {
	@succkeys = @combkeys;
	@succvals = @combvals;
	@succnames = @combnames;
}

# Ausschließen unwahrscheinlicher/ungewollter Koalitionen
my @exclkeys;
my @exclvals;
my @exclnames;
foreach my $z (0..$#succkeys){
	push(@exclkeys,$succkeys[$z]) if(exclusion($succkeys[$z]) == 1);
	push(@exclvals,$succvals[$z]) if(exclusion($succkeys[$z]) == 1);
	push(@exclnames,$succnames[$z]) if(exclusion($succkeys[$z]) == 1);
}

# Ausschließen übergroßer Koalitionen (d.h. mit überflüssigen Koalitionspartnern)
my @coalkeys;
my @coalvals;
my @coalnames;
if ($modus_b eq "A"){
	if ($modus_a eq "A" or $modus_a eq "B"){
		foreach my $y (0..$#exclkeys){
			my $index = 1;
			foreach my $ya (0..($y-1)){
				if (compare($exclkeys[$y],$exclkeys[$ya]) == 0){
					$index = 0;
					last;
				}
			}
			push(@coalkeys,$exclkeys[$y]) if($index == 1);
			push(@coalvals,$exclvals[$y]) if($index == 1);
			push(@coalnames,$exclnames[$y]) if($index == 1);
		}
	} else {
		@coalkeys = @exclkeys;
		@coalvals = @exclvals;
		@coalnames = @exclnames;
	}
} else {
	@coalkeys = @exclkeys;
	@coalvals = @exclvals;
	@coalnames = @exclnames;
}

# Überführen in ein Hash
my %coalitions;
my %coalitionnames;
foreach my $cc (0..$#coalkeys){
	$coalitions{$coalkeys[$cc]} = $coalvals[$cc];
	$coalitionnames{$coalkeys[$cc]} = $coalnames[$cc];
}

# Sortieren der Parteien nach Größe (absteigend)
my @coalparties;
my @coalvotes;
my @coalalliances;
foreach my $key (sort { $coalitions{$b} <=> $coalitions{$a} } (keys(%coalitions))) {
	push(@coalparties,$key);
	push(@coalvotes,$coalitions{$key});
	push(@coalalliances,$coalitionnames{$key});
}

# Ausgabe aller übriggebliebenen Koalitionen
foreach my $m (0..$#coalparties) {
	my $num_specchars = count_specchars($coalparties[$m]);
	my $num_width = 40 + $num_specchars;
	printf("%-${num_width}s: %3d (%5.1f %%) %-40s\n",$coalparties[$m],$coalvotes[$m],100*$coalvotes[$m]/$summe, $coalalliances[$m]);
}

# Ausgabetext falls keine Koalition übrig
print "Keine möglichen Koalitionen\n" if (scalar @coalparties == 0);

## Subroutinen
# Ausschluss ungewollter Bündnisse
sub exclusion {
	my $coalition = shift;
	my @exclusion = readExcluded();
	my @exclusion_a = pairkeys @exclusion;
	my @exclusion_b = pairvalues @exclusion;
	if ($modus_a eq "A" or $modus_a eq "C" or $modus_a eq "E"){
		foreach my $x (0..$#exclusion_a) {
			return (0) if($coalition =~ /$exclusion_a[$x]/ and $coalition =~ /$exclusion_b[$x]/);
		}
	}
	return (1);
}

# Lesen der ausgeschlossenen Bündnisse
sub readExcluded {
	my @exclusion;
	my $exclusionFile = "$FindBin::RealBin/exclusion.list";
	open(EXCL, $exclusionFile) || die "Kann Ausschlussdatei nicht öffnen\n";
	while (my $line = <EXCL>){
		chomp($line);
		if($line =~ "#") {
			next;
		} else {
			my @vals = split(/[ \t]+/,$line);
			push(@exclusion,$vals[0],$vals[1]);
		}
	}
	close(EXCL);
	return @exclusion;
}

# Summierung der Mandate einer Koalition
sub summe {
	my @arr = @_;
	my $summe = 0;
	foreach my $x (0..$#arr){
		$summe += $results{$arr[$x]};
	}
	return $summe;

}

# Vergleich zum Ausschluss übergroßer Koalitionen
sub compare {
	my ($a, $b) = @_;
	if ($modus_b eq "A"){
		if(index($a,$b) != -1){
			return 0;
		}
	}
	return 1;
}

# Zähle Sonderzeichen (Umlaute, Ligaturen)
sub count_specchars {
	my $str = shift;
	my $number = $str =~ tr/ÄäÖöÜüß/ÄäÖöÜüß/;
	$number /= 2;
	return $number;
}

# Infotext bei Abbruch wegen unzureichender Angaben
sub print_info {
	print "Argumente fehlen\n";
	print "Gebrauch: mehrheit.pl <modus_a> <modus_b> <Partei> <Sitze> ...\n";
	print "+---------+-----------------------------------------------+\n";
	print "| Modus A |                                               |\n";
	print "‡=========‡===============================================‡\n";
	print "| A       | nur über 50%, ohne ausgeschlossene Bündnisse  |\n";
	print "| B       | nur über 50%, mit ausgeschlossenen Bündnissen |\n";
	print "| C       | nur über 40%, ohne ausgeschlossene Bündnisse  |\n";
	print "| D       | nur über 40%, mit ausgeschlossenen Bündnissen |\n";
	print "| E       | alles, ohne ausgeschlossene Bündnisse         |\n";
	print "| F       | alles, mit ausgeschlossenen Bündnissen        |\n";
	print "+---------+-----------------------------------------------+\n";
	print "\n";
	print "+---------+-----------------------------------------------+\n";
	print "| Modus B |                                               |\n";
	print "‡=========‡===============================================‡\n";
	print "| A       | keine übergroßen Koalitionen                  |\n";
	print "| B       | auch übergroße Koalitionen                    |\n";
	print "+---------+-----------------------------------------------+\n";
	exit;
}


sub naming {
	my @parties = @_;
	if (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "SPD" ~~ @parties and "Grüne" ~~ @parties and "FDP" ~~ @parties and $#parties == 3) {
		return "Simbabwe-Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "SPD" ~~ @parties and "Grüne" ~~ @parties and $#parties == 2) {
		return "Kenia-Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "SPD" ~~ @parties and "FDP" ~~ @parties and $#parties == 2) {
		return "Deutschland-Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "SPD" ~~ @parties and $#parties == 1) {
		return "Große Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "Grüne" ~~ @parties and "FDP" ~~ @parties and $#parties == 2) {
		return "Jamaika-Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "FDP" ~~ @parties and "AfD" ~~ @parties and $#parties == 2) {
		return "Bahamas-Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "Grüne" ~~ @parties and $#parties == 1) {
		return "Kiwi-Koalition";
	} elsif (("CDU" ~~ @parties or "CSU" ~~ @parties or "CDU/CSU" ~~ @parties) and "FDP" ~~ @parties and $#parties == 1) {
		return "Schwarz-gelbe Koalition";
	} elsif ("SPD" ~~ @parties and "Grüne" ~~ @parties and "FDP" ~~ @parties and "Linke" ~~ @parties and $#parties == 3) {
		return "Rot-rot-grün-gelbe Koalition (R2G2)";
	} elsif ("SPD" ~~ @parties and "Grüne" ~~ @parties and "Linke" ~~ @parties and "AfD" ~~ @parties and $#parties == 3) {
		return "Eritrea-Koalition";
	} elsif ("SPD" ~~ @parties and "Grüne" ~~ @parties and "FDP" ~~ @parties and $#parties == 2) {
		return "Ampelkoalition";
	} elsif ("SPD" ~~ @parties and "Grüne" ~~ @parties and "Linke" ~~ @parties and $#parties == 2) {
		return "Rot-rot-grüne Koalition (R2G)";
	} elsif ("SPD" ~~ @parties and "Grüne" ~~ @parties and "SSW" ~~ @parties and $#parties == 2) {
		return "Dänen-Ampel oder Gambia-Koalition";
	} elsif ("SPD" ~~ @parties and "Grüne" ~~ @parties and $#parties == 1) {
		return "Rot-grüne Koalition";
	} elsif ("SPD" ~~ @parties and "Linke" ~~ @parties and $#parties == 1) {
		return "Rot-rote Koalition";
	} elsif ("SPD" ~~ @parties and "FDP" ~~ @parties and $#parties == 1) {
		return "Sozialliberale Koalition";
	} elsif ("Grüne" ~~ @parties and "FDP" ~~ @parties and $#parties == 1) {
		return "Limetten- oder Zitruskoalition";
	} else {
		return "";
	}
}