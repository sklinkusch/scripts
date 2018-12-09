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
print_exit() if ($#ARGV < 6);

my $haltestelle = join('',$ARGV[0]);
my $xnumm = join('',$ARGV[3]);
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellenum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo_ubuntu::get_number($haltestellenum);
my $station = $haltestellenarr[0];

my $typus = join('',$ARGV[1]);
print_exit() if ($typus ne 'dep' and $typus ne 'arr');

my $filtre = join('',$ARGV[2]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo_ubuntu::calc_filter($filtre);

my $num;
my $debug = 0;

my $faliniespec = join('',$ARGV[4]);
print_exit() if ($faliniespec ne 'f' and $faliniespec ne 'n');

my $fbliniespec = join('',$ARGV[5]);
print_exit() if ($fbliniespec ne 'f' and $fbliniespec ne 'n');

my @falinie;
my @fblinie;
my $read_spec = 'a';
for (my $xln = 6; $xln <= $#ARGV; $xln++) {
  if($ARGV[$xln] eq '/'){
    $read_spec = 'b';
    next;
  }
  push(@falinie,$ARGV[$xln]) if($read_spec eq 'a');
  push(@fblinie,$ARGV[$xln]) if($read_spec eq 'b');
}

my $url = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$typus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&";
my $command = "elinks -dump -dump-width 200 \"$url\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

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
  my (@btext); my (@fbtext); 
  my (@atext); my (@fatext); 
  my (@xtext); 
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$command |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_finfo(\*ACPI,$station,\@xtext);
  close ACPI;
  Fahrinfo_ubuntu::ft_sgl($faliniespec,\@falinie,\@xtext,\@fatext);
  Fahrinfo_ubuntu::ft_sgl($fbliniespec,\@fblinie,\@xtext,\@fbtext);
  Fahrinfo_ubuntu::pmax_dual(\@fatext,\@fbtext);
  Fahrinfo_ubuntu::cmp_dual(\@fatext,\@fbtext,\@text);
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
  print "Gebrauch: finfo-ffilter-ddual.pl \"Haltestelle\" <typus> <code> <Zeile> <filter-spec-arr> <filter-spec-dep> <filter-vals-arr> / <filter-vals-dep> \n";
  print "typus: dep: Abfahrtplan, arr: Ankunftplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  print "filter-spec: f: nur Einträge mit <filter-vals> berücksichtigen, n: Einträge mit <filter-vals> nicht berücksichtigen\n";
  exit;
}
