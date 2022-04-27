#!/usr/bin/perl -w
#
use Math::BigInt lib => 'GMP'; # warn if BigInt is not present
#use Math::BigInt ;
use bigint ;
use List::Util qw(max);

#
#  generate probabilites (multi molecule fingerprint) for later Bayesian classification
#  gen_mmprior.pl
#
# requires following other programs:
#
# libmath-bigint-gmp-perl
#
# (apt-get install libmath-bigint-gmp-perl)
#
#
# Usage: gen_prior.pl database_actives database_inactives
#
# Format of database_{actives|inactives}
#
# 1st line: compound name MACCS Pubchem 
# lines:    C1CC1 cyclopropan 0x240000100000000000000000000000000000800
#
# Output: database_actives_mmprior.asc
#
# output format
#
# 166 MACCS     1st line holds number of bits in that fingerprint
# 0.00000       next lines contain log(probabilities) active inactives 
# -0.03030
# 0.22222
# ...
# 880 PubChem
# 0.00000
# 0.00000
# -0.09091
# ...
#
 
$lfp[1] = 166 ; # number of bits in MACCS fingerprint
$lfp[2] = 880 ; # number of bits in PubChem fingerprint

$ACTIVES = $ARGV[0] ; 
$DECOYS  = "" ;

$MMFPDB = substr($ACTIVES,0,length($ACTIVES)-4) ;
$MMFPDB = $MMFPDB."_mmprior.asc" ; 


if (defined($ARGV[1])) {
  $DECOYS = $ARGV[1] ;
}
else {
  print "no database file for inactives found\n" ; die ;
}

# read in database(s) of actives (and decoy) molecules

# variables used for properties
my @cols ;               # holds the columns in each line

if(-e "$ACTIVES") { print "database file $ACTIVES found\n" } 
  else { print "database file $ACTIVES not found\n" ; die }

if ($DECOYS ne "") {
  if(-e "$DECOYS") { print "database file $DECOYS found\n" } 
    else { print "database file $DECOYS not found\n" ; die }
}


open(ALIST,$ACTIVES) or die "error opening $ACTIVES : $!\n"; 

$nactives = -1;  # number of compounds
$ndecoys  =  0;  # number of inactives resp. decopys
$nfp      =  0;  # number of different fingerprints of actives
$nfpd     =  0;  # number of different fingerprints of decoys

while(defined($i = <ALIST>)) {
  $nactives++;
  #  my @cols = split(/\s+/,$i);
  my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs

  if($nactives == 0) { # first line: get number of different fingerprints
    if(defined($cols[2])) {
      $nfp++ ;
      $nfpnam[1] = $cols[2] ;
    }
    if(defined($cols[3])) {
      $nfp++ ;
      $nfpnam[2] = $cols[3] ;
    }
    if(defined($cols[4])) {
      $nfp++ ;
      $nfpnam[3] = $cols[4] ;
    }
    if($nfp < 1) {
      die "no fingerprints in 1st line defined!\n"; 
    }
    next ;
  }
# lines:    C1CC1 cyclopropan 0x240000100000000000000000000000000000800

# $SMILES[$nactives] = $cols[0] ; # read SMILES string
# $smname[$nactives] = $cols[1] ; # read name of compound

  for ($j = 1; $j <= $nfp; $j++) { # read fingerprints as hex numbers
    $temp = Math::BigInt->new($cols[$j + 1]) ;
    $fp[$nactives][$j] = $temp ;
  }

# print "$SMILES[$nactives] $smname[$nactives] " ;
# print $fp[$nactives][1]->as_hex(),"\n" ;
# if ($nfp > 1) {print $fp[$nactives][2]->as_hex(),"\n" }
# print "$cols[2] \n\n" ;

}

print "number of active molecules $nactives\n";
print "number of fingerprints $nfp\n";

close(ALIST);


if ($DECOYS ne "") {
  $ndecoys = -1;
  open(DLIST,$DECOYS) or die "error opening $DECOYS : $!\n"; 

  while(defined($j = <DLIST>)) {
    $ndecoys++;
    my @cols = split(/[\s\t]/,$j); # split line at spaces or tabs
    if($ndecoys == 0) { # first line: get number of different fingerprints
      if(defined($cols[2])) {
        $nfpd++ ;
        $nfpdnam[1] = $cols[2] ;
      }
      if(defined($cols[3])) {
        $nfpd++ ;
        $nfpdnam[2] = $cols[3] ;
      }
      if(defined($cols[4])) {
        $nfpd++ ;
        $nfpdnam[3] = $cols[4] ;
      }
      if($nfpd < 1) {
        die "no fingerprints in 1st line defined!\n"; 
      }
      if ($nfpd != $nfp) {
	print "number of fingerprints in actives $nfp\n";
        print "differs from that in decoys $nfpd !\n" ;
        die ;
      }
      next ;
    }

#   $DSMILES[$ndecoys] = $cols[0] ; # read SMILES string
#   $dsmname[$ndecoys] = $cols[1] ; # read name of compound

    for ($i = 1; $i <= $nfpd; $i++) { # read fingerprints as hex numbers
      $temp = Math::BigInt->new($cols[$i + 1]) ;
      $fd[$ndecoys][$i] = $temp ;
#     print $fd[$ndecoys][$i]->as_hex(), "\n" ;
    }

  }
  close(DLIST);
  print "number of decoy molecules $ndecoys\n";
}

# compute frequency of each bit in the fingerprints
#

for ($k = 0; $k <= $lfp[2]; $k++) { # init bits for logical AND comparison
  $bit[$k+1] = Math::BigInt->new(2**$k) ;
# print $bit[$k+1]->as_hex(), "\n" ;
}

for ($j = 1; $j <= $nfp; $j++) { # init priors
  for ($k = 1; $k <= $lfp[$j]; $k++) { 
     $mmfp[$j][$k] = 0; # actives
     $mmfd[$j][$k] = 0; # decoys
   }
}

for ($i = 1; $i <= $nfp; $i++) { # loop over fingerprints

  for ($j = 1; $j <= $nactives; $j++) { # loop over actives

    for ($k = 1; $k <= $lfp[$i]; $k++) { # loop over bits in fingerprint
#     $tsum = $bit[$k]->band($fp[$j][$i]); # both version of AND work
#     $tsum = $bit[$k] & $fp[$j][$i] ;
#     print "$tsum " ;
#     if ($bit[$k]->band($fp[$j][$i]) > 0) { # this bit is set
      if (($bit[$k] & $fp[$j][$i]) > 0) { # this bit is set
        $mmfp[$i][$k]++ ; 
      }
    }
#   print "fp $i $smname[$j]\n" ;
  }

# normalize mmfd
# print "\n" ;
  for ($l = 1; $l <= $lfp[$i]; $l++) {
#   print "$mmfp[$i][$l] " ;
    $t1 = sprintf '%10.1f', $nactives ;
#   print "$t1 ";
    $t2 = sprintf '%10.1f', $mmfp[$i][$l] ;
    $tt = $t2 / $t1; # divide by occurrence = P(i|A)
    $mmfp[$i][$l] = log(max($tt,"1.0e-4")) ;
#   print "$mmfp[$i][$l]\n" ;
  }
}

if ($ndecoys > 0) { 
  for ($i = 1; $i <= $nfp; $i++) { # loop over fingerprints
    for ($j = 1; $j <= $ndecoys; $j++) { # loop over decoys
      for ($k = 1; $k <= $lfp[$i]; $k++) { # loop over bits in fingerprint
#       $tsum = $bit[$k]->band($fd[$j][$i]); 
#       print "$tsum " ;
#       print $bit[$k]->as_hex(), "\n" ;
#       if ($bit[$k]->band($fd[$j][$i]) > 0) {
        if (($bit[$k] & $fd[$j][$i]) > 0) {
          $mmfd[$i][$k]++ ; 
        }
      }
    }
#   normalize mmfd
#   print "\n" ;
    for ($l = 1; $l <= $lfp[$i]; $l++) {
      $t1 = sprintf '%10.1f', $ndecoys;
#     print "$t1 ";
      $t2 = sprintf '%10.1f', $mmfd[$i][$l] ;
      $tt = $t2 / $t1; # divide by occurrence = P(i|I)
      $mmfd[$i][$l] = log(max($tt,"1.0e-4")) ;
#     print "$mmfd[$i][$l]\n" ;
    }
  }
}

open(DATA,"> $MMFPDB") or die "error opening $MMFPDB file : $!\n"; 

for ($i = 1; $i <= $nfp; $i++) { # loop over fingerprints
  print DATA "$lfp[$i] " ; # length of fingerprint
  if ($i == 1) {print DATA "MACCS\n"} ;
  if ($i == 2) {print DATA "PubChem\n"} ;
  for ($k = 1; $k <= $lfp[$i]; $k++) { # loop over bits in fingerprint
    $t1 = $mmfp[$i][$k] ; # log(p(i|active))
    $t2 = $mmfd[$i][$k] ; # log(p(i|decoy))
    printf DATA ("%10.6f%10.6f",$t1,$t2) ;
    print  DATA "\n" ;
  } 
}

close(DATA) or die "error closing $MMFPDB file : $!\n";

exit(0);

