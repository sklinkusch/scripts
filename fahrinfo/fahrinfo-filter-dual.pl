#!/usr/bin/perl 

#
# A Perl/Tk wrapper to graphically show the arrival and departure table for a certain station/tramway stop/bus stop
#
# It calls 'elinks', takes the content of the page, and writes it to # the graphical window. It is further updated every 10 seconds.
#
# Written by SK
#

use Tk;
use strict;
use warnings;
use utf8;
use Encode;
use open ':encoding(utf8)';
use open ':std';
use FindBin;
use lib $FindBin::RealBin;
use Fahrinfo;

### set variables
print_exit() if ($#ARGV < 4);

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
my $debug = 0;

my $fliniespec = join('',$ARGV[3]);
print_exit() if ($fliniespec ne 'f' and $fliniespec ne 'n');

my @flinie;
for (my $xln = 4; $xln <= $#ARGV; $xln++) {
  push(@flinie,$ARGV[$xln]);
}

my $del;                                        # delay time (seconds)
$del = Fahrinfo::read_params(\$del);            # read parameter from external file
my $DELAY=$del * 1000;                          # delay time (milliseconds)

my $arrtypus = 'arr';
my $deptypus = 'dep';
my $arrurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$arrtypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $depurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$deptypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $arrcommand = "elinks -dump -dump-width 200 \"$arrurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $depcommand = "elinks -dump -dump-width 200 \"$depurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $spaced = ' ';

my $rheight=10;
my $rwidth=5;

# create main window
my $mw = MainWindow->new();
$mw->geometry('1250x550');
$mw->resizable(1, 1);

# set up quit key bindings
$mw->bind('<q>', sub {exit});
$mw->bind('<Escape>', sub {exit;});


my $label= $mw->Label(-font=>"Courier -12 bold",-justify=>'left')->pack;


# call the power checker now
checkNet ($label, $num);
$mw->repeat($DELAY, \&checkNet, $label, $num);

MainLoop;

############
# end main #
############

#
# subroutine to run the checking command, parse its output and 
# call the update subroutine
#
sub checkNet {
  my ($label, $num) = @_;
  my (@text);
  my (@enctext);
  my ($encline);
  my (@ftext);
  my (@fenctext);
  my (@deptext); my (@depenctext); my (@fdeptext); my (@fdepenctext);
  my (@arrtext); my (@arrenctext); my (@farrtext); my (@farrenctext);
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$arrcommand |") || die "can't open pipe!";
  Fahrinfo::read_fahrinfo(\*ACPI,$station,\@arrtext,\@arrenctext);
  close ACPI;
  Fahrinfo::filter_sgl($fliniespec,\@flinie,\@arrtext,\@arrenctext,\@farrtext,\@farrenctext);
  open (ACPJ, "$depcommand |") || die "can't open pipe!";
  Fahrinfo::read_fahrinfo(\*ACPJ,$station,\@deptext,\@depenctext);
  close ACPJ;
  Fahrinfo::filter_sgl($fliniespec,\@flinie,\@deptext,\@depenctext,\@fdeptext,\@fdepenctext);
  Fahrinfo::popmax_dual(\@farrtext,\@farrenctext,\@fdeptext,\@fdepenctext);
  Fahrinfo::compose_dual(\@farrtext,\@farrenctext,\@fdeptext,\@fdepenctext,\@text,\@enctext);
  my $nfetext = $#enctext;
  if($nfetext > -1){
    print "@text" if ($debug);
  # set the labels and battery indicator
  if ($debug == 0){
    if ($#enctext > 0 and defined $enctext[0]){
      $label->configure(-text=>"@enctext");
    }
  }
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
    print "$haltarr[$numma]\n";
    @testarr = ($haltarr[$numma],$haltnrarr[$numma]);
    return @testarr;
  }elsif($numma > 0 and defined $xnumm){
    print "$xnumm : $haltarr[$xnumm]\n";
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
  print "Gebrauch: fahrinfo-filter-dual.pl \"Haltestelle\" <code> <Zeile (optional)> <filter-spec> <filter vals>\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  print "filter-spec: f: nur Einträge mit <filter-vals> berücksichtigen, n: Einträge mit <filter-vals> nicht berücksichtigen\n";
  exit;
}
