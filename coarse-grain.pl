#! /usr/bin/perl -w

use strict;

main: {
  my $N = 29;    ## number of bins, labeled 1 through N-1
  my @_ngr = (12, 6, 3);
  my $fuzz = 1e-9;

  for my $ngr (@_ngr) {
    my $gs = $N/$ngr;
    print "$gs\n";
    for (my $ii = 0;  $ii <= $N + $fuzz; $ii += $gs) {
      print int($ii + .5), "\n";
    }
  }
}
