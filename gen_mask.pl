#!/usr/bin/perl -w
#
use Math::BigInt lib => 'GMP'; # warn if BigInt is not present
#use Math::BigInt ;
use bigint ;

# generated masking strings (non-appearing bits in the molecules gives as
# input will be set to zero)
#
# requires following other programs:
#
# libmath-bigint-gmp-perl
#
# (apt-get install libmath-bigint-gmp-perl)
#
#
# Usage: gen_mask.pl database.asc
#
# database.asc = list of molecules and their fingerprints
#
# Format of database.asc
#
# c1ccccc1 benzene 0x0ecf1s...(MACSS) 0x176000...(PubChem)
#
# Output:
# masking string "large hex number" (MACCS)
# masking string "large hex number" (Pubchem)
# 
#
 
$DBASE = $ARGV[0] ; 


if (-e "$DBASE") { print "input file $DBASE found\n" } 
  else { print "$DBASE not found\n" ; die }

my @cols ;               # holds the columns in each line
$nfp = 0; # number of fingerprints

#
# read in database of molecules with their fingerprints
 
open(BAS,$DBASE) or die "error opening input file  : $!\n" ;

$ncomp = -1 ; # number of compounds
$lfpmax = 0 ; # number of bits in largest fingerprint

while(defined($i = <BAS>)) {
# c1ccccc1 benzene 0x0ecf1s...(MACSS) 0x176000...(PubChem)
# my @cols = split(/\s+/,$i);
  my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs
  $ncomp++;
  if ($ncomp == 0) { # first line holds text, determine number of fingerprints
    if (defined($cols[2])) { $nfp++ } ;
    if (defined($cols[3])) { $nfp++ }
    if (defined($cols[4])) { $nfp++ }
    if (defined($cols[5])) { $nfp++ }
    if ($nfp < 1) {
      print "no fingerprint defined (reading line 1)?\n";
      exit(1) ;
    }
#   determine number of bits in each fingerprint
#   0x176000...   = 4**(length - 2)  each hex digit corresponds to 4 bits
#
    for ($l = 1; $l <= $nfp; $l++) { 
      $jfp[$l] = Math::BigInt->new(0x0)  ; # initial bitstring mask
      $lenfp[$l] = length($cols[$l+1]) - 2 ; # length of each fingerprint string
      for ($j = 1; $j <= $nfp; $j++) { 
      $bitfp[$j] = Math::BigInt->new(4**$lenfp[$j]) ;
#     $bitfp[$j] = 16**$lenfp[$j] ;
#     print "bitfpj $bitfp[$j] \n" ;
      if ($bitfp[$j] > $lfpmax) { $lfpmax = $bitfp[$j] }
      }
#     print "lfpmax $lfpmax \n" ;
    }
  next ;
  }

# $smile[$ncomp] = $cols[0] ;
# $name[$ncomp]  = $cols[1] ;
   
# bitwise OR the fingerprints while reading them and do not store them
# at all (saves large amounts of memory)
#
  for ($j = 1; $j <= $nfp; $j++) { # read fingerprints as hex numbers
    $fp[$j] = Math::BigInt->new($cols[$j+1]) ;
    $jfp[$j] = $jfp[$j] | $fp[$j] ;
  }
}

close(BAS);

print "number of compounds $ncomp  number of fingerprints $nfp\n" ;



# generate masking strings
# report unused features in the fingerprints
print "\nunused featurs and masking strings\n\n" ;

for ($m = 1; $m <= $nfp; $m++) { # loop over fingerprints
  for ($n = 0; $n < $bitfp[$m]; $n++) { # loop over bits in fingerprint
    $bset = $jfp[$m] & 2**$n ; # respective bit is set
    if ($bset < 1) {
      $nn = $n + 1;
      print "feature $nn\n" ;
    }
  }
  print "\n" ;
  print $jfp[$m]->as_hex(),"\n\n" ;
}

print "\n" ;

exit(0);

