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
print_exit() unless ($#ARGV >= 0);

my $nrhaltestellen = join('',$ARGV[0]);
my $nroarg = 3 * $nrhaltestellen + 3;
print_exit() unless ($#ARGV >= $nroarg);
my @haltestellen;
my @xnumm;
my $maxj = 500;
my @haltestellennr;
my @station;
my @temparr;
my @url;
my @command;
foreach my $xh (0..($nrhaltestellen - 1)){
 $haltestellen[$xh] = join('',$ARGV[($xh + 1)]);
 $xnumm[$xh] = join('',$ARGV[($xh + $nrhaltestellen + 3)]);
 @temparr = haltnummer($haltestellen[$xh],$xnumm[$xh]);
 push(@station,$temparr[0]);
 push(@haltestellennr,Fahrinfo::get_number($temparr[1]));
}
my $typus = join('',$ARGV[($nrhaltestellen + 1)]);
print_exit() if ($typus ne 'dep' and $typus ne 'arr');

my $filtre = join('',$ARGV[($nrhaltestellen + 2)]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo::calc_filter($filtre);

my @fls;
foreach my $xh (0..$#haltestellen){
  push(@fls,$ARGV[($xh + 2 * $nrhaltestellen + 3)]);
  print_exit() if ($fls[$#fls] ne 'f' and $fls[$#fls] ne 'n');
}
my @fli;
my @xfli;
my $ref_xfli;
my $temparg;
my $nrslash = 0;
my @slash;
my @tmparr;
foreach my $xf ((3 * $nrhaltestellen + 3)..$#ARGV){
  $temparg = join('',$ARGV[$xf]);
  push(@tmparr,$temparg);
}
foreach my $xj (0..$#tmparr){
    push (@fli,$tmparr[$xj]) if ($tmparr[$xj] ne '/') ;
    push (@slash,$#fli) if ($tmparr[$xj] eq '/');
    $nrslash++ if ($tmparr[$xj] eq '/');
}
push(@slash,$#fli);
print_exit() unless ($nrslash == $#haltestellen);

my $num;
my $debug = 0;

my $del;                                        # delay time (seconds)
$del = Fahrinfo::read_params(\$del);            # read parameter from external file
my $DELAY=$del*1000;                            # delay time (milliseconds)

foreach my $xh (0..$#haltestellennr) {
 $url[$xh] = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellennr[$xh]&boardType=$typus&maxJourneys=$maxj&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&";
 $command[$xh] = "elinks -dump -dump-width 200 \"$url[$xh]\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g;s/è/e/g;s/ê/e/g;s/à/a/g;s/á/a/g;s/â/a/g\'";
}

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
  my @pretext; my @preenctext;
  my @preftext; my @prefenctext;
  my $xn = 0;
  my @fxs;
  my @prenr;
  foreach my $xh (0..$#haltestellennr) {
   open(ACPI, "$command[$xh] |") || die "can't open pipe $xh!";
   Fahrinfo::read_strict_pn(\*ACPI,$xh,$station[$xh],\@prenr,\@pretext,\@preenctext);
   close ACPI;
  }
  foreach my $xh (0..$#haltestellennr) {
    Fahrinfo::filter_sgl_pn($xh,\@fls,\@fli,\@slash,\@prenr,\@pretext,\@preenctext,\@preftext,\@prefenctext,\@fxs);
  }
  Fahrinfo::sort_entries_n($#haltestellennr,\@fxs,\@preftext,\@prefenctext,\@text,\@enctext);
  Fahrinfo::popmax_sgl(\@text,\@enctext);
  Fahrinfo::add_linebreaks(\@text,\@enctext);
  my $nrtext = $#text;
  if($nrtext > -1){
    print "@text" if ($debug);
  # set the labels and battery indicator
    $label->configure(-text=>"@enctext");
  }
}

sub haltnummer {
  my ($halt,$xnumm) = @_;
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
  print "Gebrauch: fahrinfo-multi-sort-filter.pl <#Haltestellen> \"Haltestelle1\" ... <dep/arr> <code> <Zeile1> ... <filter-spec> <filter-vals>\n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  print "filter-spec: f: nur Einträge mit <filter-vals> berücksichtigen, n: Einträge mit <filter-vals> nicht berücksichtigen\n";
  exit;
}
