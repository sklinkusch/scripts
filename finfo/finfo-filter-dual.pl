#!/usr/bin/perl 

use strict;
use warnings;
use utf8;
use Encode;
use open ':encoding(utf8)';
use open ':std';
use FindBin;
use lib $FindBin::RealBin;
use Fahrinfo_ubuntu;

### set variables
print_exit() if ($#ARGV < 4);

my $haltestelle = join('',$ARGV[0]);
my $xnumm = join('',$ARGV[2]);
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellenum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo_ubuntu::get_number($haltestellenum);
my $station = $haltestellenarr[0];

my $filtre = join('',$ARGV[1]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo_ubuntu::calc_filter($filtre);

my $num;

my $fliniespec = join('',$ARGV[3]);
print_exit() if ($fliniespec ne 'f' and $fliniespec ne 'n');

my @flinie;
for (my $xln = 4; $xln <= $#ARGV; $xln++) {
  push(@flinie,$ARGV[$xln]);
}

my $arrtypus = 'arr';
my $deptypus = 'dep';
my $arrurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$arrtypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $depurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$deptypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $arrcommand = "elinks -dump -dump-width 200 \"$arrurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $depcommand = "elinks -dump -dump-width 200 \"$depurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $spaced = ' ';

# call the power checker now
checkNet ($num);

############
# end main #
############

#
# subroutine to run the checking command, parse its output and 
# call the update subroutine
#
sub checkNet {
  my ($num) = @_;
  my (@text);
  my (@ftext);
  my (@deptext); my (@fdeptext); 
  my (@arrtext); my (@farrtext); 
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$arrcommand |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_finfo(\*ACPI,$station,\@arrtext);
  close ACPI;
  Fahrinfo_ubuntu::ft_sgl($fliniespec,\@flinie,\@arrtext,\@farrtext);
  open (ACPJ, "$depcommand |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_finfo(\*ACPJ,$station,\@deptext);
  close ACPJ;
  Fahrinfo_ubuntu::ft_sgl($fliniespec,\@flinie,\@deptext,\@fdeptext);
  Fahrinfo_ubuntu::pmax_dual(\@farrtext,\@fdeptext);
  Fahrinfo_ubuntu::cmp_dual(\@farrtext,\@fdeptext,\@text);
  my $nfetext = $#text;
  if($nfetext > -1){
    print "@text";
  }
}

sub haltnummer {
  my $halt = shift;
  open(DATA, "$FindBin::RealBin/../data/fahrinfo-elinks2.dat") || die "can't open 'fahrinfo-elinks2.dat'";
  my $numma = -1;
  my $str;
  my $haltu;
  my @haltarr;
  my @haltnrarr;
  while (my $line = <DATA>) {
    $haltu = decode_utf8($halt);
    $str = decode_utf8($line);
    if(index($line, $halt) != -1) {
      $numma++;
      $str =~ /(\"[a-zA-Z0-9\.\,\-\/\+\[\]\(\) ]+\")([\t])([0-9]+)/;
      push(@haltarr,$1);
      push(@haltnrarr,$3);
    }
  }
  my @testarr;
  if($numma == 0){
    @testarr = ($haltarr[$numma],$haltnrarr[$numma]);
    return @testarr;
  }elsif($numma > 0 and defined $xnumm){
    @testarr = ($haltarr[$xnumm],$haltnrarr[$xnumm]);
    return @testarr;
  }else{
    my $lnum = $numma + 1;
    print "$lnum Haltestellen gefunden: \n";
    foreach my $axnum (0..$#haltarr) {
      print "$axnum : $haltarr[$axnum]\n";
    }
    exit;
  }
}

sub print_exit {
  print "Gebrauch: finfo-filter-dual.pl \"Haltestelle\" <code> <Zeile (optional)> <filter-spec> <filter vals>\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  print "filter-spec: f: nur Einträge mit <filter-vals> berücksichtigen, n: Einträge mit <filter-vals> nicht berücksichtigen\n";
  exit;
}
