#!/usr/bin/perl -w
#
# compute the AUC of the MACCS and PubChem fingerprint from
# the scores of actives and inactives
# Difference to usual AUC: hrer the score covers the range from -1.0 to +1.0 !
#
# Usage: gen_auc.pl <scores_actives> <scores_inactives>
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
# Output: auc.txt
#
# output format (space separated) for use as spreadsheet
#
# 

use Scalar::Util qw(looks_like_number);

 
$ACTIVES = "" ; 
if (defined($ARGV[0])) {$ACTIVES = $ARGV[0]}
$DECOYS  = "" ;
if (defined($ARGV[1])) {$DECOYS = $ARGV[1]}

$AUCF ="auc.txt" ;


$stepw = 50 ; # step width (number of bins)
$step2 = $stepw / 2.0 ;

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

while(defined($i = <ALIST>)) {
  my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs
  $nactives++;

  if (defined($cols[2])) {   
    if (looks_like_number($cols[2])) {$actf1[$nactives]=$cols[2]} # is a number?
      else {
        $nactives--;
        next ;
      }
  }
  
  if (defined($cols[3])) {   
    if (looks_like_number($cols[3])) {$actf2[$nactives]=$cols[3]} # is a number?
      else {
        $nactives--;
	next ;
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

  if (defined($cols[2])) {   
    if (looks_like_number($cols[2])) {$inaf1[$ndecoys]=$cols[2]} # is a number?
      else {
        $ndecoys--;
        next ;
      }
  }
  
  if (defined($cols[3])) {   
    if (looks_like_number($cols[3])) {$inaf2[$ndecoys]=$cols[3]} # is a number?
      else {
        $ndecoys--;
	next ;
      }
  }

  next ;
}

close(DLIST);

print "number of inactives $ndecoys\n";


# init arrays

$sumtp1 = 0 ; # sum of true positives fingerprint 1 
$sumfp1 = 0 ; # sum of false positives fingerprint 1 
$sumtn1 = 0 ; # sum of true negatives
$sumfn1 = 0 ; # sum of false negatives
$sumtp2 = 0 ; # sum of true positives fingerprint 2 
$sumfp2 = 0 ; # sum of false positives fingerprint 2 
$sumtn2 = 0 ; # sum of true negatives
$sumfn2 = 0 ; # sum of false negatives

# Remark:
# In this context all correct cases are true positives (=true positives + true negatives)
# whereas all wrong cases are false positives (=false positives + false negatives)

for ($l = 0; $l <= $stepw; $l++) {
   $tp1[$l] = 0 ;      # true positives fingerprint 1
   $fp1[$l] = 0 ;      # false positives
   $tn1[$l] = 0 ;      # true negatives
   $fn1[$l] = 0 ;      # false negatives
   $tp2[$l] = 0 ;      # true positives fingerprint 2
   $fp2[$l] = 0 ;      # false positives
   $tn2[$l] = 0 ;      # true negatives
   $fn2[$l] = 0 ;      # false negatives
   $ttp1[$l] = 0 ;     # sum of true positives fingerprint 1    
   $tfp1[$l] = 0 ;     # sum of false positives 
   $ttp2[$l] = 0 ;     # sum of true positives fingerprint 2
   $tfp2[$l] = 0 ;     # sum of false positives 
   $tpr1[$l] = 0.0 ;   # true positive rate fingerprint 1
   $fpr1[$l] = 0.0 ;   # false positive rate
   $tpr2[$l] = 0.0 ;   # true positive rate fingerprint 2
   $fpr2[$l] = 0.0 ;   # false positive rate
   $aucf1[$l] = 0.0 ;  # AUC fraction fingerprint 1
   $aucf2[$l] = 0.0 ;  # AUC fraction fingerprint 2
}

# init true positive rate and false positive rate
# both rates start from zero and approach unity
$ttp1[-1] = 0.00 ; 
$tfp1[-1] = 0.00 ; 
$ttp2[-1] = 0.00 ; 
$tfp2[-1] = 0.00 ; 

$tpr1[-1] = 0.00 ;
$fpr1[-1] = 0.00 ; 
$tpr2[-1] = 0.00 ;
$fpr2[-1] = 0.00 ; 

$tpr1[$stepw + 1] = 1.00 ;
$fpr1[$stepw + 1] = 1.00 ;
$tpr2[$stepw + 1] = 1.00 ;
$fpr2[$stepw + 1] = 1.00 ;


# loop over actives

for ($l = 1; $l <= $nactives; $l++) {
  $thr = int($actf1[$l] * $step2) + $step2 ;
  if ($actf1[$l] < 0.0) { # false negative, here a "false positive" 
    $fp1[$thr]++ ;
    $fn1[$thr]++ ;
    }
    else {  # true positive
      $tp1[$thr]++ ;
    }

  $thr = int($actf2[$l] * $step2) + $step2 ;
  if ($actf2[$l] < 0.0) { # false negative, here a "false positive"
    $fp2[$thr]++ ;
    $fn2[$thr]++ ;
    }
    else {
      $tp2[$thr]++ ;
    }
}


# loop over inactives

for ($l = 1; $l <= $ndecoys; $l++) {
  $thr = int($inaf1[$l] * $step2) + $step2 ;
  if ($inaf1[$l] < 0.0) {  # true negative, here a "true positive" 
    $tp1[$thr]++ ;
    $tn1[$thr]++ ;
    }
    else {  # false postive
      $fp1[$thr]++ ;
    }

  $thr = int($inaf2[$l] * $step2) + $step2 ;
  if ($inaf2[$l] < 0.0) { # true negative, here a "true positive"
    $tp2[$thr]++ ;
    $tn2[$thr]++ ;
    }
    else {  # false postive
      $fp2[$thr]++ ;
    }
}

# compute sum of true positives and sum of false positives

for ($l = 0; $l <= $stepw; $l++) {
  $sumtp1 = $sumtp1 + $tp1[$l] ; # sum of true positives fingerprint 1 
  $sumfp1 = $sumfp1 + $fp1[$l] ; # sum of false positives fingerprint 1 
  $sumtn1 = $sumtn1 + $tn1[$l] ; # sum of true negatives
  $sumfn1 = $sumfn1 + $fn1[$l] ; # sum of false negatives

  $sumtp2 = $sumtp2 + $tp2[$l] ; # sum of true positives fingerprint 2 
  $sumfp2 = $sumfp2 + $fp2[$l] ; # sum of false positives fingerprint 2 
  $sumtn2 = $sumtn2 + $tn2[$l] ; # sum of true negatives
  $sumfn2 = $sumfn2 + $fn2[$l] ; # sum of false negatives
}

# now loop over threshold range to obtain sum of true and false positives 

for ($l = 0; $l <= $stepw; $l++) {
  $ttp1[$l] = $ttp1[$l-1] + $tp1[$l] ; # sum up true positives fingerprint 1
  $tfp1[$l] = $tfp1[$l-1] + $fp1[$l] ; # sum of false positives 

  $ttp2[$l] = $ttp2[$l-1] + $tp2[$l] ; # fingerprint 2
  $tfp2[$l] = $tfp2[$l-1] + $fp2[$l] ;

  $tpr1[$l] = $ttp1[$l] / ($sumtp1 + $sumfn1) ; # true positive rate
  $fpr1[$l] = $tfp1[$l] / ($sumfp1 + $sumtn1) ; # false positive rate

  $tpr2[$l] = $ttp2[$l] / ($sumtp2 + $sumfn2) ; # true positive rate
  $fpr2[$l] = $tfp2[$l] / ($sumfp2 + $sumtn2) ; # false positive rate

}

if ($tpr1[$stepw] < 1.0) {$tpr1[$stepw] = 1.00}
if ($fpr1[$stepw] < 1.0) {$fpr1[$stepw] = 1.00}
if ($tpr2[$stepw] < 1.0) {$tpr2[$stepw] = 1.00}
if ($fpr2[$stepw] < 1.0) {$fpr2[$stepw] = 1.00}


# loop over threshold range 

$auct1 = 0.00 ; # total AUC fingerprint 1
$auct2 = 0.00 ; # total AUC fingerprint 2

for ($l = 0; $l < $stepw; $l++) {
  # area of a trapez = hf * 1/2 (ta + tc)
  $tc = $tpr1[$l+1] ; 
  $ta = $tpr1[$l] ;
  $hf = $fpr1[$l+1] - $fpr1[$l] ;
  $aucf1[$l] = $hf * 0.5 * ($ta + $tc) ;
  $auct1 = $auct1 + $aucf1[$l] ;

  $tc = $tpr2[$l+1] ; 
  $ta = $tpr2[$l] ;
  $hf = $fpr2[$l+1] - $fpr2[$l] ;
  $aucf2[$l] = $hf * 0.5 * ($ta + $tc) ;
  $auct2 = $auct2 + $aucf2[$l] ;
}

# since false negative were counted as "false positives" (= fail) and
# likewise true negatives counted as "true positives" (= success), 
# the sum of TN, respectively FN have to be corrected

$sens1 = $sumtp1 / ($sumtp1 + $sumfn1) ; # sensitivity fingerprint 1
$sens2 = $sumtp2 / ($sumtp2 + $sumfn2) ; # sensitivity fingerprint 2

$spec1 = $sumtn1 / ($sumtn1 + $sumfp1)  ; # specificity fingerprint 1
$spec2 = $sumtn2 / ($sumtn2 + $sumfp2)  ; # specificity fingerprint 2

#$acc1  = ($sumtp1 + $sumtn1) / ($nactives + $ndecoys) ; # accuracy fingerprint 1
#$acc2  = ($sumtp2 + $sumtn2) / ($nactives + $ndecoys) ; # accuracy fingerprint 2

$acc1  = ($sumtp1 + $sumtn1) / ($sumtp1 + $sumtn1 + $sumfp1 + $sumfn1) ; # accuracy fingerprint 1
$acc2  = ($sumtp2 + $sumtn2) / ($sumtp2 + $sumtn2 + $sumfp2 + $sumfn2) ; # accuracy fingerprint 1


print "fingerprint sensitivity specificity accuracy AUC\n" ;
print "MACCS      ";
printf "%8.5f%8.5f%8.5f%9.5f\n", $sens1, $spec1, $acc1, $auct1 ;
print "TP TN FP FN" ;
printf "%6i%6i%6i%6i\n", $sumtp1, $sumtn1, $sumfp1, $sumfn1 ;
print "PubChem    ";
printf "%8.5f%8.5f%8.5f%9.5f\n", $sens2, $spec2, $acc2, $auct2 ;
print "TP TN FP FN" ;
printf "%6i%6i%6i%6i\n", $sumtp2, $sumtn2, $sumfp2, $sumfn2 ;


# print range, TP, FP, sum_TP, sum_FP, FPR, TPR, AUC for each range

open(DATA,"> $AUCF") or die "error opening $AUCF file : $!\n"; 

print DATA "range TP FP sumTP sumFP FPR TPR AUC" ;
print DATA " TP2 FP2 sumTP2 sumFP2 FPR2 TPR2 AUC2\n" ;

for ($l=0; $l<=$stepw; $l++) { 

  $thr = -1.0 + ($l * 2.0 / $stepw) ; 
  printf DATA ("%5.3f %i %i %i %i %7.5f %7.5f %8.5f",$thr,$tp1[$l],$fp1[$l],$ttp1[$l],$tfp1[$l],$fpr1[$l],$tpr1[$l],$aucf1[$l]) ;
  printf DATA (" %i %i %i %i %7.5f %7.5f %8.5f",$tp2[$l],$fp2[$l],$ttp2[$l],$tfp2[$l],$fpr2[$l],$tpr2[$l],$aucf2[$l]) ;
  print  DATA "\n" ;
}
printf DATA ("%s %8.5f %s %8.5f","AUC1",$auct1,"AUC2",$auct2);
print  DATA "\n" ;
print  DATA "\n" ;

printf DATA ("%s %i\n","number of actives", $nactives) ;
printf DATA ("%s %i\n","number of decoys", $ndecoys) ;
print  DATA "\n" ;

print DATA "fingerprint sensitivity specificity accuracy AUC\n" ;
print DATA "MACCS ";
printf DATA "%8.5f%8.5f%8.5f%9.5f\n", $sens1, $spec1, $acc1, $auct1 ;
print DATA "TP TN FP FN " ;
printf DATA "%i %i %i %i\n", $sumtp1, $sumtn1, $sumfp1, $sumfn1 ;
print DATA "PubChem ";
printf DATA "%8.5f%8.5f%8.5f%9.5f\n", $sens2, $spec2, $acc2, $auct2 ;
print DATA "TP TN FP FN " ;
printf DATA "%i %i %i %i\n", $sumtp2, $sumtn2, $sumfp2, $sumfn2 ;


close(DATA) or die "error closing $AUCF file : $!\n";

exit(0);

