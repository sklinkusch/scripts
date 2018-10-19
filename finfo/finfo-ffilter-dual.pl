#!/usr/bin/perl 

use strict;
use warnings;
use utf8;
use Encode;
use open ':encoding(utf8)';
use open ':std';
use FindBin;
use lib "/home/stefan/bin/";
use Fahrinfo;


### set variables
print_exit() if ($#ARGV < 5);

my $haltestelle = join('',$ARGV[0]);
my $xnumm = join('',$ARGV[2]);
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellenum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo::get_number($haltestellenum);
my $station = $haltestellenarr[0];

my $filtre = join('',$ARGV[1]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo::calc_filter($filtre);

my $num;

my $farrliniespec = join('',$ARGV[3]);
print_exit() if ($farrliniespec ne 'f' and $farrliniespec ne 'n');

my $fdepliniespec = join('',$ARGV[4]);
print_exit() if ($fdepliniespec ne 'f' and $fdepliniespec ne 'n');

my @farrlinie;
my @fdeplinie;
my $read_spec = 'a';
for (my $xln = 5; $xln <= $#ARGV; $xln++) {
  if($ARGV[$xln] eq '/'){
    $read_spec = 'd';
    next;
  }
  push(@farrlinie,$ARGV[$xln]) if($read_spec eq 'a');
  push(@fdeplinie,$ARGV[$xln]) if($read_spec eq 'd');
}

my $arrtypus = 'arr';
my $deptypus = 'dep';
my $arrtipp = 'Ankunft';
my $deptipp = 'Abfahrt';
my $arrurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$arrtypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&";
my $depurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$deptypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&";
my $arrcommand = "elinks -dump -dump-width 200 \"$arrurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $depcommand = "elinks -dump -dump-width 200 \"$depurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

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
  my $num = shift;
  my (@text);
  my (@ftext);
  my (@deptext); my (@fdeptext); 
  my (@arrtext); my (@farrtext); 
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$arrcommand |") || die "can't open pipe!";
  Fahrinfo::read_finfo(\*ACPI,$station,\@arrtext);
  close ACPI;
  Fahrinfo::ft_sgl($farrliniespec,\@farrlinie,\@arrtext,\@farrtext);
  open (ACPJ, "$depcommand |") || die "can't open pipe!";
  Fahrinfo::read_finfo(\*ACPJ,$station,\@deptext);
  close ACPJ;
  Fahrinfo::ft_sgl($fdepliniespec,\@fdeplinie,\@deptext,\@fdeptext);
  Fahrinfo::pmax_dual(\@farrtext,\@fdeptext);
  Fahrinfo::cmp_dual(\@farrtext,\@fdeptext,\@text);
  my $nfetext = $#text;
  if($nfetext > -1){
    print "@text";
  }
}

sub haltnummer {
  my $halt = shift;
  open(DATA, '/home/stefan/bin/fahrinfo-elinks2.dat') || die "can't open 'fahrinfo-elinks2.dat'";
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
  print "Gebrauch: fahrinfo-ffilter-dual.pl \"Haltestelle\" <code> <Zeile> <filter-spec-arr> <filter-spec-dep> <filter-vals-arr> / <filter-vals-dep> \n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  print "filter-spec: f: nur Einträge mit <filter-vals> berücksichtigen, n: Einträge mit <filter-vals> nicht berücksichtigen\n";
  exit;
}
