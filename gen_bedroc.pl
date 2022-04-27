#!/usr/bin/perl -w
#
# compute the BEDROC of the MACCS and PubChem fingerprint from
# the scores of actives and inactives
# for alpha values of 20 and 100 (earlier part of accumulation curve will count more).
#
# Usage: gen_bedroc.pl <scores_actives> <scores_inactives>
#
# Format of scores_{actives|inactives}:
#
# list of SMILES ACE_actives_db.asc found
# file ACE_actives_db_mmfp.asc found
# number of fingerprints 2
# SMILES name similarities
# C[C@H](N[C@@H](CCC1=CC=CC=C1)C(O)=O)C(=O)N2CCC[C@H]2C(O)=O enalaprilat 0.34326 0.59024
#
# text will be ignored, only third and forth column are used.
#
# 

use Scalar::Util qw(looks_like_number);
 
$ACTIVES = "" ; 
if (defined($ARGV[0])) {$ACTIVES = $ARGV[0]}
$DECOYS  = "" ;
if (defined($ARGV[1])) {$DECOYS = $ARGV[1]}


# read in database(s) of actives (and decoy) molecules

# variables used for properties
my @cols ;               # holds the columns in each line

if(-e "$ACTIVES") { print "file $ACTIVES found\n" } 
  else { print "file $ACTIVES not found\n" ; die }

if(-e "$DECOYS") { print "file $DECOYS found\n" } 
  else { print "file $DECOYS not found\n" ; die }


open(ALIST,$ACTIVES) or die "error opening $ACTIVES : $!\n"; 

$nactives = 0;  # number of actives
$ndecoys  = 0;  # number of inactives 
$ncomp    = 0;  # total number of compounds


while(defined($i = <ALIST>)) {
  my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs
  $nactives++;
  $ncomp++ ;

  if (defined($cols[2])) {   
    if (looks_like_number($cols[2])) { # is a number?
      $actf1[$ncomp]=$cols[2] ;
      $act[$ncomp] = "1" ;
    }
    else {
      $nactives--;
      $ncomp--;
      next ;
    }
  }
  
  if (defined($cols[3])) {   
    if (looks_like_number($cols[3])) {
      $actf2[$ncomp]=$cols[3] ;
    }
  }

  next ;
}

close(ALIST);
 
print "number of actives $nactives\n";

#for ($i = 1; $i <= $nactives; $i++) { 
#  print "$actf1[$i] $actf2[$i] \n" ;
#}

open(DLIST,$DECOYS) or die "error opening $DECOYS : $!\n"; 

while(defined($j = <DLIST>)) {
  my @cols = split(/[\s\t]/,$j); # split line at spaces or tabs
  $ndecoys++;
  $ncomp++ ;

  if (defined($cols[2])) {   
    if (looks_like_number($cols[2])) {
      $actf1[$ncomp]=$cols[2] ;
      $act[$ncomp] = "0" ;
    }
    else {
      $ndecoys--;
      $ncomp--;
      next ;
    }
  }
  
  if (defined($cols[3])) {   
    if (looks_like_number($cols[3])) {
      $actf2[$ncomp]=$cols[3] ;
    }
  }

  next ;
}

close(DLIST);

print "number of inactives $ndecoys\n";


# alpha values for RIE and BEDROC
$alpha[1] =  20 ; 
$alpha[2] =  50 ; 
$alpha[3] = 100 ; 
$maxal = 3 ;

# compute RIE_min(alpha) and RIE_max(alpha)
#
$ttt = $ncomp / $nactives ;
for ($l = 1; $l <= $maxal; $l++) {
  $riemin[$l] = $ttt * (1.0 - exp($alpha[$l] * $nactives / $ncomp)) ;
  $riemin[$l] = $riemin[$l] / (1.0 - exp($alpha[$l])) ;
  $riemax[$l] = $ttt * (1.0 - exp($alpha[$l] * -1.0 * $nactives / $ncomp)) ;
  $riemax[$l] = $riemax[$l] / (1.0 - exp($alpha[$l] * -1.0)) ;
}

for ($l = 1; $l <= $maxal; $l++) {
  print "RIEmin($alpha[$l]) $riemin[$l]\n" ;
  print "RIEmax($alpha[$l]) $riemax[$l]\n" ;
}

# set up matrices to be sorted
# elements of the matrix start at zero
#
for ($l = 1; $l <= $ncomp; $l++) {
  push(@matf1,[$actf1[$l],$act[$l]]) ; # fingerprint 1
  push(@matf2,[$actf2[$l],$act[$l]]) ; # fingerprint 2
}

#print "$matf1[0][0] $matf1[0][1]\n" ;
#print "$matf1[1][0] $matf1[1][1]\n" ;
#print "$matf1[2][0] $matf1[2][1]\n" ;

# sort matrices according to activity (hightest first)
#
@smatf1 = sort{$b->[0] <=> $a->[0]} @matf1 ;
@smatf2 = sort{$b->[0] <=> $a->[0]} @matf2 ;


# debug print sorted matrices
#for ($l = 0; $l <= 50; $l++) {
#  $a0 = $smatf1[$l][0] ;
#  $a1 = $smatf1[$l][1] ;
#  $b0 = $smatf2[$l][0] ;
#  $b1 = $smatf2[$l][1] ;
#  print "$a0 $a1   $b0 $b1\n" ;
#}

# determine rank of each active
#
$nactf1 = 0 ;
$nactf2 = 0 ;
for ($l = 1; $l <= $ncomp; $l++) {
  $t1 = $smatf1[$l-1][1] ;
  $t2 = $smatf2[$l-1][1] ;
  if ($t1 eq "1") { 
    $nactf1++ ;
    $rif1[$nactf1] = $l ; # rank of active fingerprint 1
  }
  if ($t2 eq "1") { 
    $nactf2++ ;
    $rif2[$nactf2] = $l ; # rank of active fingerprint 1
  }
}

if ($nactf1 > $nactives) {print "FP1 more actives $nactf1 than expected $nactives\n"} 
if ($nactf2 > $nactives) {print "FP2 more actives $nactf2 than expected $nactives\n"} 
if ($nactf1 < $nactives) {print "FP1 less actives $nactf1 than expected $nactives\n"}
if ($nactf2 < $nactives) {print "FP2 less actives $nactf2 than expected $nactives\n"}

# compute RIE(alpha) and BEDROC(alpha) for each alpha (and both fingerprints)
#
for ($l = 1; $l <= $maxal; $l++) { # loop over alpha values
  $sumf1[$l] = 0 ;
  $sumf2[$l] = 0 ;
  for ($k = 1; $k <= $nactives; $k++) { # loop over actives
    $sumf1[$l] = $sumf1[$l] + exp($alpha[$l] * -1.0 * $rif1[$k] / $ncomp) ; 
    $sumf2[$l] = $sumf2[$l] + exp($alpha[$l] * -1.0 * $rif2[$k] / $ncomp) ; 
  }
  $kkk = (1.0 - exp($alpha[$l] * -1.0)) / (exp($alpha[$l] / $ncomp) - 1.0) ;
  if ($kkk < 1.0e-10) {$kkk = 1.0e-10} ;
  $rief1[$l] = $ttt * $sumf1[$l] / $kkk ;
  $rief2[$l] = $ttt * $sumf2[$l] / $kkk ;

  $bedf1[$l] = ($rief1[$l] - $riemin[$l]) / ($riemax[$l] - $riemin[$l]) ;
  $bedf2[$l] = ($rief2[$l] - $riemin[$l]) / ($riemax[$l] - $riemin[$l]) ;
}


print "\nfingerprint MACCS PubChem\n" ;
for ($l = 1; $l <= $maxal; $l++) { # loop over alpha values
  print "   RIE($alpha[$l]) $rief1[$l] $rief2[$l]\n" ;
  print "BEDROC($alpha[$l]) $bedf1[$l] $bedf2[$l]\n" ;
}

exit(0);

