package Fahrinfo_ubuntu;
use strict;
use List::Util qw(reduce any all none notall first max maxstr min minstr product sum sum0 pairs unpairs pairkeys pairvalues pairgrep pairfirst pairmap shuffle uniq uniqnum uniqstr);

sub add_lbr {
  my $in  = $_[0];
  foreach my $x (0..$#$in){
    $$in[$x] = sprintf("%s\n",$$in[$x]) unless ($x == 0);
    $$in[$x] = sprintf(" %s\n",$$in[$x]) if ($x == 0);
  }
}

sub calculate_differenz {
   my ($plantime, $realtime) = (@_);
   $plantime =~ /([012][0-9])\:([012345][0-9])/;
   my $planstd = $1;
   my $planmte = $2;
   $realtime =~ /([012][0-9])\:([012345][0-9])/;
   my $realstd = $1;
   my $realmte = $2;
   my $differenz;
   $differenz = 60*($realstd - $planstd) + ($realmte - $planmte);
   if ($differenz < -120){
     $differenz += 1440;
   }
   return $differenz;
}

sub calc_filter {
  my $filtra = shift;
  my $fbahn = 0; my $rbahn = 0; my $sbahn = 0; my $ubahn = 0; my $tram = 0; my $bus = 0; my $ferry = 0;
  RETURN:
    if ($filtra >= 64 and $filtra < 128){
      $fbahn = 1;
      $filtra -= 64;
      goto RETURN;
    }
    if ($filtra >= 32 and $filtra < 64){
      $rbahn = 1;
      $filtra -= 32;
      goto RETURN;
    }
    if ($filtra >= 16 and $filtra < 32){
      $sbahn = 1;
      $filtra -= 16;
      goto RETURN;
    }
    if ($filtra >= 8 and $filtra < 16){
      $ubahn = 1;
      $filtra -= 8;
      goto RETURN;
    }
    if ($filtra >= 4 and $filtra < 8){
      $tram = 1;
      $filtra -= 4;
      goto RETURN;
    }
    if ($filtra == 2 or $filtra == 3){
      $bus = 1;
      $filtra -= 2;
      goto RETURN;
    }
    if ($filtra == 1){
      $ferry = 1;
    }
  my $filtnum = sprintf("%1d%1d%1d%1d%1d%1d%1d%1d%1d", $sbahn, $ubahn, $tram, $bus, $ferry, $fbahn, $rbahn, 0, 1);
  return $filtnum;
}

sub calc_hour {
  my $ctime = shift;
  $ctime = /([\d]{2})\:([\d]{2})/;
  my $hour = $1;
  return $hour;
}

sub cmp_dual {
  my $ina = $_[0];
  my @inax = @{$_[0]};
  my $ind = $_[1];
  my @indx = @{$_[1]};
  my $t = $_[2];
  my $et = $_[3];
  my $nat = $#inax;
  my $ndt = $#indx;
  my $nz = min($nat, $ndt);
  foreach my $xa (0..$nz){
    $$t[$xa] = $inax[$xa] unless ($xa == 0);
    $$t[$xa] = sprintf(" %s",$inax[$xa]) if ($xa == 0);
  }
  if($ndt > $nat){
    foreach my $xz (($nz+1)..$ndt){
      $$t[$xz] = sprintf("%-87s",$main::spaced);
    }
  }elsif($ndt < $nat){
    foreach my $xz (($nz+1)..$nat){
      push(@indx,"\n");
      $$t[$xz] = $$t[$xz] . "\n";
    }
  }
  foreach my $xd (0..$ndt){
    $$t[$xd] = $$t[$xd] . $indx[$xd] . "\n";
  }
}

sub cmp_fstr_dual {
  my $ina = $_[0];
  my @inax = @{$_[0]};
  my $ind = $_[1];
  my @indx = @{$_[1]};
  my $t = $_[2];
  my $nat = $#inax;
  my $ndt = $#indx;
  my $spaced = ' ';
  foreach my $xa (0..$nat){
    if($xa == 0){
      $$t[$xa] = sprintf(" %-87s",$inax[$xa]);
    }else{
      $$t[$xa] = $inax[$xa];
    }
  }
  if($ndt > $nat){
    foreach my $xz (($nat+1)..$ndt){
      $$t[$xz] = sprintf("%-87s",$spaced);
    }
  }elsif($ndt < $nat){
    foreach my $xz (($ndt+1)..$nat){
      push(@indx,"\n");
      $$t[$xz] = $$t[$xz] . "\n";
    }
  }
  foreach my $xd (0..$ndt){
    $$t[$xd] = $$t[$xd] . $indx[$xd] . "\n";
  }
}

sub def_spec {
 my ($dtime,$A,$d,$D,$f,$F,$s,$S,$N) = @_;
 my $spec;
 $spec = 'A' if ($dtime <= $A);
 $spec = 'a' if ($dtime > $A and $dtime < 0);
 $spec = 'P' if ($dtime == 0);
 $spec = 'd' if ($dtime > 0 and $dtime <= $d);
 $spec = 'D' if ($dtime > $d and $dtime <= $D);
 $spec = 'f' if ($dtime > $D and $dtime <= $f);
 $spec = 'F' if ($dtime > $f and $dtime <= $F);
 $spec = 's' if ($dtime > $F and $dtime <= $s);
 $spec = 'S' if ($dtime > $s and $dtime <= $S);
 $spec = 'N' if ($dtime > $S and $dtime <= $N);
 $spec = 'x' if ($dtime > $N);
 return $spec;
}

sub define_product {
   my @arr = @_;
   my $linie = $arr[0];
   my $fahrtziel = $arr[1];
   my $spec = '0';
   $spec = 'F' if ($linie =~ /ICE|IC|EC|RJ|D|EN|NJ|TGV|THA|HBX|FLX/);
   $spec = 'R' if ($linie =~ /[I]?[R]{1}[BE]{1}/);
   $spec = 'S' if ($linie =~ /[S]{1}(1|2|25|26|3|41|42|45|46|47|5|7|75|8|85|9)$/);
   $spec = 'S' if ($linie =~ /^[0-9]{4,6}$/ and $fahrtziel =~ /^S /);
   $spec = 'S' if ($linie =~ /^[0-9]{4,6}$/ and $fahrtziel =~ /^Ringbahn S4[12]/);
   $spec = 'U' if ($linie =~ /[U]{1}(1|2|3|4|5|55|6|7|8|9)$/);
   $spec = 'M' if ($linie =~ /^[M]{1}(1|2|4|5|6|8|10|13|17)$/);
   $spec = 'm' if ($linie =~ /^[M]{1}(11|19|21|27|32|37|41|44|45|46|48|49|76|77|82|85)$/);
   $spec = 'T' if ($linie =~ /^[1-9]{1}[0-9]{1}$/);
   $spec = 'x' if ($linie =~ /^[X]{1}[1-9][0-9]?$/);
   $spec = 'x' if ($linie =~ /TXL/);
   $spec = 'b' if ($linie =~ /^[1-9]{1}[0-9]{2}$/);
   $spec = 'n' if ($linie =~ /^[N]{1}[1-9][0-9]?$/);
   $spec = 'f' if ($linie =~ /^[F]{1}[0-9]{1,2}$/);
   $spec = 'e' if ($spec eq '0');
   return $spec;
}

sub define_specifier {
   my @inp = @_;
   my $dtime = $inp[0];
   my $tspec = $inp[1];
   my $ctime = $inp[2];
   my $hour = Fahrinfo_ubuntu::calc_hour($ctime);
   my $spec;
   goto BUS   if ($tspec eq 'b' or $tspec eq 'x');
   goto ELSE  if ($tspec eq 'e');
   goto FERN  if ($tspec eq 'F');
   goto FERRY if ($tspec eq 'f');
   goto METRO if (($tspec eq 'M' or $tspec eq 'm') and $hour > 4);
   goto NACHT if (($tspec eq 'M' or $tspec eq 'm') and $hour < 5);
   goto NACHT if ($tspec eq 'n');
   goto REGIO if ($tspec eq 'R');
   goto SBAHN if ($tspec eq 'S');
   goto TRAM  if ($tspec eq 'T');
   goto UBAHN if ($tspec eq 'U');
   FERN:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-10,9,29,59,89,119,179,209);
   goto AUSW;
   REGIO:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-5,4,9,19,29,59,89,119);
   goto AUSW;
   SBAHN:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,4,9,14,19,39);
   goto AUSW;
   UBAHN:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,3,5,7,9,19);
   goto AUSW;
   METRO:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,3,5,7,9,29);
   goto AUSW;
   TRAM:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,4,9,14,19,39);
   goto AUSW;
   BUS:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,4,9,14,19,39);
   goto AUSW;
   NACHT:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,4,9,19,29,59);
   goto AUSW;
   FERRY:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-2,1,2,3,4,5,6,9);
   goto AUSW;
   ELSE:
   $spec = Fahrinfo_ubuntu::def_spec($dtime,-3,2,3,7,9,14,19,39);
   goto AUSW;
   AUSW:
   return $spec;
}

sub ft_sgl {
  my $fls = $_[0];
  my $fl  = $_[1];
  my $ot  = $_[2];
  my $ft  = $_[3];
  my $and = 0;
  my @open;
  my @close;
  my @ort = @$ot;
  my @prt;
  push (@$ft,$$ot[0]);
  $and = 1 if (join(' ',@$fl) =~ /[\{\}]/);
  if ($and == 1){
    foreach my $xz (0..$#$fl){
      if ($$fl[$xz] =~ /^\{.+$/){
        push(@open,$xz);
        $$fl[$xz] = substr($$fl[$xz],1);
      }
      if ($$fl[$xz] =~ /^.+\}$/){
        push(@close,$xz);
        $/ = "\}";
        chomp($$fl[$xz]);
      }
    }
    main::print_exit() if(($#open - $#close) != 0);
    foreach my $xz (0..$#open){
      main::print_exit() if(($close[$xz] - $open[$xz]) != 1);
    }
    if ($fls eq 'f'){
      foreach my $xz (1..$#$ot){
        foreach my $xy (0..$#open){
          if ($$ot[$xz] =~ /$$fl[$open[$xy]]/ and $$ot[$xz] =~ /$$fl[$close[$xy]]/){
            push(@prt,$$ot[$xz]);
            last;
          }
        }
      }
    }else{
      foreach my $xz (1..$#$ot){
        foreach my $xy (0..$#open){
          if ($$ot[$xz] =~ /$$fl[$open[$xy]]/ and $$ot[$xz] =~ /$$fl[$close[$xy]]/){
            $$ot[$xz]  = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (0..$#open){
      $$fl[$open[$xy]]   = undef;
      $$fl[$close[$xy]] = undef;
    }
  }
  if ($fls eq 'f'){
    foreach my $xz (1..$#$ot) {
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($$ot[$xz] =~ /$$fl[$xy]/){
            push(@prt,$$ot[$xz]);
            last;
          }
        }
      }
    }
  }else{
    foreach my $xz (1..$#$ot) {
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($$ot[$xz] =~ /$$fl[$xy]/){
            $$ot[$xz]  = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (1..$#$ot) {
      if (defined $$ot[$xy]) {
        push(@prt,$$ot[$xy]);
      }
    }
  }
  foreach my $xy (1..$#ort){
    foreach my $xz (0..$#prt){
      if ($ort[$xy] eq $prt[$xz]){
        push(@$ft,$ort[$xy]);
        last;
      }
    }
  }
}

sub ft_sgl_p {
  my $fls = $_[0];
  my $fl  = $_[1];
  my $ot  = $_[2];
  my $ft  = $_[3];
  my $and = 0;
  my @open;
  my @close;
  my @ort = @$ot;
  my @prt;
  $and = 1 if (join(' ',@$fl) =~ /[\{\}]/);
  if ($and == 1){
    foreach my $xz (0..$#$fl){
      if ($$fl[$xz] =~ /^\{.+$/){
        push(@open,$xz);
        $$fl[$xz] = substr($$fl[$xz],1);
      }
      if ($$fl[$xz] =~ /^.+\}$/){
        push(@close,$xz);
        $/ = "\}";
        chomp($$fl[$xz]);
      }
    }
    main::print_exit() if(($#open - $#close) != 0);
    foreach my $xz (0..$#open){
      main::print_exit() if(($close[$xz] - $open[$xz]) != 1);
    }
    if ($fls eq 'f'){
      foreach my $xz (0..$#$ot){
        foreach my $xy (0..$#open){
          if ($$ot[$xz] =~ /$$fl[$open[$xy]]/ and $$ot[$xz] =~ /$$fl[$close[$xy]]/){
            push(@prt,$$ot[$xz]);
            last;
          }
        }
      }
    }else{
      foreach my $xz (0..$#$ot){
        foreach my $xy (0..$#open){
          if ($$ot[$xz] =~ /$$fl[$open[$xy]]/ and $$ot[$xz] =~ /$$fl[$close[$xy]]/){
            $$ot[$xz]  = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (0..$#open){
      $$fl[$open[$xy]]   = undef;
      $$fl[$close[$xy]] = undef;
    }
  }
  if ($fls eq 'f'){
    foreach my $xz (0..$#$ot) {
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($$ot[$xz] =~ /$$fl[$xy]/){
            push(@prt,$$ot[$xz]);
            last;
          }
        }
      }
    }
  }else{
    foreach my $xz (0..$#$ot) {
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($$ot[$xz] =~ /$$fl[$xy]/){
            $$ot[$xz]  = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (0..$#$ot) {
      if (defined $$ot[$xy]) {
        push(@prt,$$ot[$xy]);
      }
    }
  }
  foreach my $xy (0..$#ort){
    foreach my $xz (0..$#prt){
      if ($ort[$xy] eq $prt[$xz]){
        push(@$ft,$ort[$xy]);
        last;
      }
    }
  }
}

sub get_number {
 my $str = shift;
 $str =~ /([0-9]{1})([0-9]{6})/;
 my $add = '00';
 my $res = $1 . $add . $2;
 return $res;
}

sub make_fstr {
  my $pt  = $_[0];
  my $t   = $_[1];
  my $exs;
  my $ini;
  foreach my $x (0..$#$pt){
    if($x != 0){
      my $exs = substr($$pt[$x],87,55);
      if($exs =~ /^[ ]*$/ or $exs =~ /Gleis/){
        $ini = substr($$pt[$x],0,87);
        push(@$t,$ini);
      }else{
        next;
      }
    }else{
      $ini = substr($$pt[$x],0,87);
      push(@$t,$ini);
    }
  }
}

sub pmax_dual {
  my $ina = $_[0];
  my $ind = $_[1];
  my $nat = $#$ina;
  my $ndt = $#$ind;
  my $max = 35;
  my $diff;
  my $tempvar;
  if($nat > $max){
    $diff = $nat - $max;
    while ($diff > 0){
      $tempvar = pop(@$ina);
      $diff--;
    }
  }
  if($ndt > $max){
    $diff = $ndt - $max;
    while ($diff > 0){
      $tempvar = pop(@$ind);
      $diff--;
    }
  }
}

sub pmax_sgl {
  my $in  = $_[0];
  my $nt  = $#$in;
  my $max = 35;
  my $diff;
  my $tempvar;
  if($nt > $max){
    $diff = $nt - $max;
    while ($diff > 0){
      $tempvar = pop(@$in);
      $diff--;
    }
  }
}

sub read_finfo {
  my $fs = $_[0];
  my $station = $_[1];
  my $t = $_[2];
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $product;
  my $specifier;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $spaced = ' ';
  my $fahrt;
  while (my $line = <$fs>){
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf("%-86s",$begline);
      push(@$t,$begaline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo_ubuntu::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            if(defined $realzeit){
              $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $realzeit);
              $specifier = Fahrinfo_ubuntu::define_specifier($difftime,$product,$currtime);
            }else{
              $difftime = 0;
              $specifier = '?';
            }
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z0-9]{2,3})[ ]{1,}/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([ ]{1})([A-Z0-9]{2,3})([ ]{1,})/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $skip_request == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = 'X';
            my $resttime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $currtime);
            $difftime = undef;
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $skip_request == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and defined $planzeit and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ \t]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            }
            if(defined $realzeit){
              $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit,$realzeit);
            }else{
              $difftime = 0;
            }
            $specifier = 'Z';
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            if (defined $difftime){
              my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
              if (defined $realzeit){
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel);
              }else{
                if (defined $linie){
                  $fahrt = sprintf("%-5s ----- %-1s %+4d $-1s %-9s %-55s",$planzeit,$specifier,$difftime,$product,$linie,$fziel);
                }else{
                  $fahrt = sprintf("%-5s ----- %-1s %+4d %-1s %-9s %-55s",$planzeit,$specifier,$difftime,$product,$spaced,$fziel);
                }
              }
            }else{
              my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
              if (defined $realzeit){
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel);
              }else{
                if (defined $linie){
                  $fahrt = sprintf("%-5s ----- %-1s %-4s %-1s %-9s %-55s",$planzeit,$specifier,$spaced,$product,$linie,$fziel);
                }else{
                  $fahrt = sprintf("%-5s ----- %-1s %-4s %-1s %-9s %-55s",$planzeit,$specifier,$spaced,$product,$spaced,$fziel);
                }
              }
            }
            push(@$t,$fahrt);
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
}

sub read_str {
  my $fs = $_[0];
  my $station = $_[1];
  my $t = $_[2];
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $def_extra = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $specifier;
  my $extra_stop;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $product;
  my $spaced = ' ';
  while (my $line = <$fs>) {
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf ("%-86s",$begline);
      push(@$t,$begaline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo_ubuntu::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9][ ]{1,}[A-Z]{0,3}[ ]?[0-9]{0,5}[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis [A-Za-z0-9\,\.\-\+\/\(\)\[\]]*/){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $extra_stop = $6;
            $def_extra = 1;
            $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo_ubuntu::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo_ubuntu::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $1;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $2;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $difftime = undef;
            $specifier = '?';
            $extra_stop = $3;
            $def_extra = 1;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $specifier = 'X';
            $difftime = undef;
            my $resttime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $currtime);
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            }
            $specifier = 'Z';
            $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit,$realzeit);
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $curline =~ /Haltestelle\/Gleis/) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
            $curline =~ /Haltestelle\/Gleis ([A-Za-z0-9\+\-\.\,\/\(\)\[\] ]*)/;
            $extra_stop = $1;
            $def_extra = 1;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            my $fahrt;
            if ($def_extra == 1 and defined $extra_stop){
              if (defined $difftime){
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }
            }else{
              if (defined $difftime){
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$spaced) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$spaced) unless ($specifier eq 'x');
              }
            }
            push(@$t,$fahrt) unless ($specifier eq 'x');
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            $def_extra = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
}

sub read_str_p {
  my $fs = $_[0];
  my $station = $_[1];
  my $t = $_[2];
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $def_extra = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $specifier;
  my $extra_stop;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $product;
  my $spaced = ' ';
  while (my $line = <$fs>) {
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf ("%-86s",$begline);
#      push(@$t,$begaline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo_ubuntu::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9][ ]{1,}[A-Z]{0,3}[ ]?[0-9]{0,5}[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis [A-Za-z0-9\,\.\-\+\/\(\)\[\]]*/){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $extra_stop = $6;
            $def_extra = 1;
            $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo_ubuntu::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo_ubuntu::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $1;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $2;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $difftime = undef;
            $specifier = '?';
            $extra_stop = $3;
            $def_extra = 1;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $specifier = 'X';
            $difftime = undef;
            my $resttime = Fahrinfo_ubuntu::calculate_differenz($planzeit, $currtime);
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            }
            $specifier = 'Z';
            $difftime = Fahrinfo_ubuntu::calculate_differenz($planzeit,$realzeit);
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $curline =~ /Haltestelle\/Gleis/) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
            $curline =~ /Haltestelle\/Gleis ([A-Za-z0-9\+\-\.\,\/\(\)\[\] ]*)/;
            $extra_stop = $1;
            $def_extra = 1;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo_ubuntu::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            my $fahrt;
            if ($def_extra == 1 and defined $extra_stop){
              if (defined $difftime){
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }
            }else{
              if (defined $difftime){
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$station) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo_ubuntu::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$station) unless ($specifier eq 'x');
              }
            }
            if($fahrt =~ /$station/){
              push(@$t,$fahrt) unless ($specifier eq 'x');
            }
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            $def_extra = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
}

sub remove_blanks {
 my $str = shift;
 $str =~ s/^\s*//;
 return $str;
}

sub search_linie {
  my $curline = shift;
  my $linie;
  if($curline =~ / [A-Z]{2,3} [0-9]{1,5} /){
    $curline =~ / ([A-Z]{2,3} [0-9]{1,5}) /;
    $linie = $1;
  }elsif($curline =~ /TXL/){
    $linie = 'TXL';
  }elsif($curline =~ / [MNSUX]{1}[0-9]{1,2} /){
    $curline =~ / ([MNSUX]{1}[0-9]{1,2}) /;
    $linie = $1;
  }elsif($curline =~ / [0-9]{3} /){
    $curline =~ / ([0-9]{3}) /;
    $linie = $1;
  }
  return $linie;
}

sub st_entr {
  my $at  = @_[0];
  my $bt  = @_[1];
  my $t   = @_[2];
  for (my $x = 0; $x <= ($#$at + $#$bt); $x++){
    if(defined $$at[0] and defined $$bt[0]){
      if($$at[0] lt $$bt[0]){
        push(@$t,$$at[0]);
        shift(@$at);
      }else{
        push(@$t,$$bt[0]);
        shift(@$bt);
      }
    }elsif(defined $$at[0]){
      push(@$t,$$at[0]);
      shift(@$at);
    }else{
      push(@$t,$$bt[0]);
      shift(@$bt);
    }
  }
}

sub subst_text {
      my $lina = shift;
      $lina =~ s/\[[1-9][0-9]\]/    /g;
      $lina =~ s/\[[1-9][0-9]{2}\]/     /g;
      $lina =~ s/\[[1-9][0-9]{3}\]/      /g;
      $lina =~ s/\[IMG\]/     /g;
      $lina =~ s/^[ \t]*$//g;
      return $lina;
}

1;
