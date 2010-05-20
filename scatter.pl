#! /usr/bin/perl -w

use strict;

main: {
  my @raw = (1,2,3,4,5,6,
        40,43,44,47,48,51,
        111,112,113,117,118,119);

  my @diag = ();
  my @diag2 = ();
  my @diag40 = ();
  my @cume = ();
  my @other = ();
  
  for my $point (@raw) {
    push @cume, "$point, -3";
    push @cume, "-3, $point";
    for my $flip (@raw) {
      if ($flip == $point) {
        push @diag, "$point, $flip";
      } elsif ($flip == $point+2) {
        push @diag2, "$point, $flip";
      } elsif ($flip == $point+40) {
        push @diag40, "$point, $flip";
      } else {
        push @other, "$point, $flip";
      }
    }
  }

  for (my $ii = 0; $ii <= $#other; $ii++){
    print $other[$ii];
    if (0 && $ii <= $#cume) {
      print ", $cume[$ii]";
    }
    if ($ii <= $#diag) {
      print ", $diag[$ii]";
    }
    if ($ii <= $#diag2) {
      print ", $diag2[$ii]";
    }
    if ($ii <= $#diag40) {
      print ", $diag40[$ii]";
    }
    print "\n";
  }
}