#!/usr/bin/perl -w
#
use Math::BigInt lib => 'GMP'; # warn if BigInt is not present
#use Math::BigInt ;
use bigint ;

#  calculate probability for being active or inactive based on the probabilities
#  bayes_mmpf.pl
#
# requires following other programs:
#
# libmath-bigint-gmp-perl
#
# (apt-get install libmath-bigint-gmp-perl)
#
#
# Usage: bayes_mmpf.pl compounds_with_fps mmprior_db.asc
#
# compounds_with_fps = list of molecule with their fingerprints
# mmpriordb.asc = probabilities as previously calculated by gen_mmprior.pl
#
# Format of compounds_with_fps:
#
# c1ccccc1 benzene 0x0ecf1s...(MACSS) 0x176000...(PubChem)
#
# Output:
# SMILES#1 mmfp-similarity
# SMILES#2 mmfp-similarity
# [...]

# masking strings derived from cleaned Dark Chemical Matter set
# Bits that are not set in these strings do not appear and are
# neglected upon calculating the similarity
#
$mk[1] = Math::BigInt->new('0xffffc70f9fffffffffffffffffffffffffffffffe') ; # MACCS
$mk[2] = Math::BigInt->new('0x1e1ffe5cccccccc88cccc863c400065e200afe301ff84aff8a57fc512dffffffffffffffffffffffffffdeeffffefeefffbfffeff6ffdeffebff77fddb7d87e7fedff5ffefffc1fc3d7ffbd7fffffff7bffffffffff004c477df1f06800ffffffffffefcf000000000000400037bc') ; # PubChem

 
$SMILES = $ARGV[0] ; 
$MMFPDB = "" ;

if (defined($ARGV[1])) {
  $MMFPDB = $ARGV[1];
}
else {
  print "missing a MMFP file!\n";
  die ;
}


if (-e "$SMILES") { print "list of SMILES $SMILES found\n" } 
  else { print "$SMILES not found\n" ; die }

if (-e "$MMFPDB") { print "file $MMFPDB found\n" } 
  else { print "$MMFPDB not found\n" ; die }

my @cols ;               # holds the columns in each line
$nfp = 0; # number of fingerprints

#
# read in MMFP file first
 
open(MMFP,$MMFPDB) or die "error opening MMFPDB : $!\n" ;

$nline = 0 ;

while(defined($i = <MMFP>)) {
  $nline++;
  my @cols = split(/\s+/,$i); # split line at spaces
# my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs
  if ($cols[1] =~ m/[A-Z]/) { # start of new fingerprint
# if (defined($cols[1])) {  # start of new fingerprint
#   print "\n" ;
    $nfp++;
    $lfp[$nfp] = $cols[0] ; # number of fingerprint bits
    $nc = 0 ; # reset counter 
    next ;
  }
  $nc++;
  $mmfp[$nfp][$nc] = sprintf '%10.6f', $cols[1] ; # log(prob(i|A)) actives
  $mmfd[$nfp][$nc] = sprintf '%10.6f', $cols[2] ; # log(prob(i|I)) decoys
# print "$mmfp[$nfp][$nc] $mmfd[$nfp][$nc]\n" ;
}
#print "\n" ;

close(MMFP);

print "number of fingerprints $nfp\n" ;

open(LIST,$SMILES) or die "error opening $SMILES : $!\n"; 

$ncomp = -1;  # number of compounds

while(defined($i = <LIST>)) {
  $ncomp++;
  if ($ncomp == 0) { # first line holds text
    next ;
  }
  #  my @cols = split(/\s+/,$i);
# c1ccccc1 benzene 0x0ecf1s...(MACSS) 0x176000...(PubChem)
  my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs
  $smile[$ncomp] = $cols[0] ;
  $name[$ncomp]  = $cols[1] ;
  for ($j = 1; $j <= $nfp; $j++) { # read fingerprints as hex numbers
     $fp[$ncomp][$j] = Math::BigInt->new($cols[$j+1]) ;
  }
}

close(LIST);

$lfpmax = 0 ; # number of bits in largest fingerprint
for ($k = 1; $k <= $nfp; $k++) { 
  if ($lfp[$k] > $lfpmax) {$lfpmax = $lfp[$k]}
}
#print "lfpmax $lfpmax \n" ;

for ($k = 0; $k < $lfpmax; $k++) { # init bits for logical AND comparison
  $bit[$k+1] = Math::BigInt->new(2**$k) ;
# print $bit[$k+1]->as_hex(), "\n" ;
}

print "SMILES name log_odds_scores\n";

for ($l = 1; $l <= $ncomp; $l++) { # loop over compounds
  print "$smile[$l] $name[$l] ";
  for ($m = 1; $m <= $nfp; $m++) { # loop over fingerprints
    $sim = sprintf '%7.5f', 0.0e0 ;
    for ($k = 1; $k <= $lfp[$m]; $k++) { # loop over bits in fingerprint
##    $bmask = $bit[$k] & $mk[$m] ; # bit is not masked
#     print "$m $bmask\n" ;
##    if ($bmask < 1) {next}
#     $tsum = $bit[$k]->band($fp[$l][$m]); # both version of AND work
      $tsum = $bit[$k] & $fp[$l][$m] ; # bit in query molecules is set
#     print "$tsum " ;
      if ($tsum > 0) { # if bit is set in molecule
        $sim = $sim + $mmfp[$m][$k] - $mmfd[$m][$k] ; # log odds ratio
#       $mmfp[$nfp][$nc]  # log(prob(i|A)) actives
#       $mmfd[$nfp][$nc]  # log(prob(i|I)) decoys
      }
    }
    printf ("%7.5f",$sim);
    print " ";
  }
  print "\n" ;
}

exit(0);

