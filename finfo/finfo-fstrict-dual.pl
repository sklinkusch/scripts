#!/usr/bin/perl 

use strict;
use warnings;
use utf8;
use Encode;
use open ':encoding(utf8)';
use open ':std';
use FindBin;
use lib $FindBin::Bin;
use Fahrinfo_ubuntu;

### set variables
print_exit() unless ($#ARGV == 1 or $#ARGV == 2);

my $haltestelle = join('',$ARGV[0]);
my $xnumm = join('',$ARGV[2]) if ($#ARGV == 2);
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellenum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo_ubuntu::get_number($haltestellenum);
my $station = $haltestellenarr[0];

my $filtre = join('',$ARGV[1]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo_ubuntu::calc_filter($filtre);

my $num;

my $arrtypus = 'arr';
my $deptypus = 'dep';
my $arrurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$arrtypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $depurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$deptypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
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
  my @deptext; my @arrtext;
  my @depptext; my @arrptext;
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$arrcommand |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_str(\*ACPI,$station,\@arrptext);
  close ACPI;
  Fahrinfo_ubuntu::make_fstr(\@arrptext,\@arrtext);

  open (ACPJ, "$depcommand |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_str(\*ACPJ,$station,\@depptext);
  close ACPJ;
  Fahrinfo_ubuntu::make_fstr(\@depptext,\@deptext);
 
  Fahrinfo_ubuntu::pmax_dual(\@arrtext,\@deptext);
  Fahrinfo_ubuntu::cmp_fstr_dual(\@arrtext,\@deptext,\@text);
  my $nrtext = $#text;
  if($nrtext > -1){
    print "@text";
  }
}

sub haltnummer {
  my $halt = shift;
  open(DATA, "$FindBin::Bin/../data/fahrinfo-elinks2.dat") || die "can't open 'fahrinfo-elinks2.dat'";
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
  print "Gebrauch: finfo-fstrict-dual.pl \"Haltestelle\" <code> <Zeile (optional)>\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}
