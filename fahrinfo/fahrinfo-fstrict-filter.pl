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
use lib "/home/stefan/bin/";
use Fahrinfo;

### set variables
print_exit() if ($#ARGV < 5);

my $haltestelle = join('',$ARGV[0]);
my $xnumm = join('',$ARGV[3]);
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellenum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo::get_number($haltestellenum);
my $station = $haltestellenarr[0];

my $typus = join('',$ARGV[1]);
print_exit() if ($typus ne 'dep' and $typus ne 'arr');

my $filtre = join('',$ARGV[2]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo::calc_filter($filtre);

my $num;
my $debug = 0;

my $fliniespec = join('',$ARGV[4]);
print_exit() if ($fliniespec ne 'f' and $fliniespec ne 'n');

my @flinie;
for (my $xln = 5; $xln <= $#ARGV; $xln++) {
  push(@flinie,$ARGV[$xln]);
}

my $del;                                        # delay time (seconds)
$del = Fahrinfo::read_params(\$del);            # read parameter from external file
my $DELAY=$del*1000;                            # delay time (milliseconds)

my $url = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$typus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&";
my $command = "elinks -dump -dump-width 200 \"$url\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

my $rheight=10;
my $rwidth=5;

# create main window
my $mw = MainWindow->new();
$mw->geometry('1000x550');
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
  my (@ftext);
  my (@fenctext);
  my @pretext; my @preenctext;
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$command |") || die "can't open pipe!";
  Fahrinfo::read_strict(\*ACPI,$station,\@pretext,\@preenctext);
  close ACPI;
  Fahrinfo::make_fstrict(\@pretext,\@preenctext,\@text,\@enctext);
  Fahrinfo::filter_sgl($fliniespec,\@flinie,\@text,\@enctext,\@ftext,\@fenctext);
  Fahrinfo::popmax_sgl(\@ftext,\@fenctext);
  Fahrinfo::add_linebreaks(\@ftext,\@fenctext);
  my $nrtext = $#ftext;
  if($nrtext > -1){
    print "@ftext" if ($debug);
  }
  # set the labels and battery indicator
  if($#fenctext > -1 and defined $fenctext[0]){
    $label->configure(-text=>"@fenctext");
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
  print "Gebrauch: fahrinfo-fstrict-filter.pl \"Haltestelle\" <dep/arr> <code> <Zeile> <filter-spec> <filter-vals>\n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  print "filter-spec: f/n, f: nur <vals> beruecksichtigen, n: <vals> nicht beruecksichtigen\n";
  exit;
}
