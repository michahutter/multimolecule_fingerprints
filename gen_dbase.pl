#!/usr/bin/perl -w
#
use Math::BigInt ;
use bigint ;

#
#  generate a database containing the fingerprints of
#  the SMILES given in the input.
#
# requires following other programs:
#
# libmath-bigint-gmp-perl
#
# (apt-get install libmath-bigint-gmp-perl)
#
# Open Babel version 2.3.2 or later 
#
# (apt-get install babel)
#
# Usage: gen_dbase.pl listname.smi  
#
# Output: listname_db.asc
#
# listname may either contain SMILES or SMILES and names:
# CCCCO butyl-alcohol 
#
#
# directory where babel and obgrep are located
$BABELDIR = "/usr/bin/" ;
$obgrep = $BABELDIR."obgrep " ;
#$babel =  $BABELDIR."babel " ;

$CLIST  = $ARGV[0] ; 

# read in list of compounds

# variables used for properties
my @cols ;               # holds the columns in each line

if(-e "$CLIST") { print "list of compounds found\n" } 
  else { print "list of compounds $CLIST not found\n" ; die }

open(COMPLIST,$CLIST) or die "error opening $CLIST : $!\n"; 

$DBFILE = substr($CLIST,0,length($CLIST)-4) ;
$DBFILE = $DBFILE."_db.asc" ; 

$ncomp = 0;

while(defined($i = <COMPLIST>)) {
  $ncomp++;
  #  my @cols = split(/\s+/,$i);
  my @cols = split(/[\s\t]/,$i); # split line at spaces or tabs
  if (defined($cols[1])) {
    $names = 1;     #  list of SMILES with names in second column
    $smname[$ncomp] = $cols[1] ;
  }
  else {
    $names = 0;     #  list of SMILES only
    $cols[1] = "";
    $smname[$ncomp] = "c"."$ncomp" ; # generic compound names
  }
 #print "$cols[0] $cols[1] \n";
  $SMILES[$ncomp]=$cols[0] ;
}

print "number of compounds $ncomp \n";

close(COMPLIST);


$matchl = 1 ;
$flag = 0 ;
$jjj = "jjj" ; # temporary file


# public MACCS keys as defined in MACCS.cpp of RDKIT
#
$maxmekeys = 146 ; # maximum number of "easy" MACCS keys   

# define MACCS keys as SMARTS

$mc[1]   = "[!#6!#1]1~*~*~*~1";              # bit8
$mc[2]   = "*1~*~*~*~1";                     # bit11
$mc[3]   = "[#8]~[#7](~[#6])~[#6]";          # bit13
$mc[4]   = "[#16]-[#16]";                    # bit14
$mc[5]   = "[#8]~[#6](~[#8])~[#8]";          # bit15
$mc[6]   = "[!#6!#1]1~*~*~1";                # bit16
$mc[7]   = "[#6]#[#6]";                      # bit17
$mc[8]   = "*1~*~*~*~*~*~*~1";               # bit19
$mc[9]   = "[#14]";                          # bit20
$mc[10]  = "[#6]=[#6](~[!#6!#1])~[!#6!#1]";  # bit21
$mc[11]  = "*1~*~*~1";                       # bit22
$mc[12]  = "[#7]~[#6](~[#8])~[#8]";          # bit23
$mc[13]  = "[#7]-[#8]";                      # bit24
$mc[14]  = "[#7]~[#6](~[#7])~[#7]";          # bit25
$mc[15]  = "[#6]=@[#6](@*)@*";               # bit26
$mc[16]  = "[!#6!#1]~[CH2]~[!#6!#1]";        # bit28
$mc[17]  = "[#6]~[!#6!#1](~[#6])(~[#6])~*";  # bit30
$mc[18]  = "[!#6!#1]~[F,Cl,Br,I]";           # bit31
$mc[19]  = "[#6]~[#16]~[#7]";                # bit32
$mc[20]  = "[#7]~[#16]";                     # bit33
$mc[21]  = "[CH2]=*";                        # bit34
$mc[22]  = "[#16R]";                         # bit36
$mc[23]  = "[#7]~[#6](~[#8])~[#7]";          # bit37
$mc[24]  = "[#7]~[#6](~[#6])~[#7]";          # bit38
$mc[25]  = "[#8]~[#16](~[#8])~[#8]";         # bit39
$mc[26]  = "[#16]-[#8]";                     # bit40
$mc[27]  = "[#6]#[#7]";                      # bit41
$mc[28]  = "[!#6!#1!H0]~*~[!#6!#1!H0]";      # bit43
$mc[29]  = "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]";# bit44
$mc[30]  = "[#6]=[#6]~[#7]";                 # bit45
$mc[31]  = "[#16]~*~[#7]";                   # bit47
$mc[32]  = "[#8]~[!#6!#1](~[#8])~[#8]";      # bit48
$mc[33]  = "[!+0]";                          # bit49 charge
$mc[34]  = "[#6]=[#6](~[#6])~[#6]";          # bit50
$mc[35]  = "[#6]~[#16]~[#8]";                # bit51
$mc[36]  = "[#7]~[#7]";                      # bit52
$mc[37]  = "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]";  # bit53
$mc[38]  = "[!#6!#1!H0]~*~*~[!#6!#1!H0]";    # bit54
$mc[39]  = "[#8]~[#16]~[#8]";                # bit55
$mc[40]  = "[#8]~[#7](~[#8])~[#6]";          # bit56
$mc[41]  = "[#8R]";                          # bit57
$mc[42]  = "[!#6!#1]~[#16]~[!#6!#1]";        # bit58
$mc[43]  = "[#16]!:*:*";                     # bit59
$mc[44]  = "[#16]=[#8]";                     # bit60
$mc[45]  = "*~[#16](~*)~*";                  # bit61
$mc[46]  = "*@*!@*@*";                       # bit62
$mc[47]  = "[#7]=[#8]";                      # bit63
$mc[48]  = "*@*!@[#16]";                     # bit64
$mc[49]  = "c:n";                            # bit65
$mc[50]  = "[#6]~[#6](~[#6])(~[#6])~*";      # bit66
$mc[51]  = "[!#6!#1]~[#16]";                 # bit67
$mc[52]  = "[!#6!#1!H0]~[!#6!#1!H0]";        # bit68
$mc[53]  = "[!#6!#1]~[!#6!#1!H0]";           # bit69
$mc[54]  = "[!#6!#1]~[#7]~[!#6!#1]";         # bit70
$mc[55]  = "[#7]~[#8]";                      # bit71
$mc[56]  = "[#8]~*~*~[#8]";                  # bit72
$mc[57]  = "[#16]=*";                        # bit73
$mc[58]  = "[CH3]~*~[CH3]";                  # bit74
$mc[59]  = "*!@[#7]@*";                      # bit75
$mc[60]  = "[#6]=[#6](~*)~*";                # bit76
$mc[61]  = "[#7]~*~[#7]";                    # bit77
$mc[62]  = "[#6]=[#7]";                      # bit78
$mc[63]  = "[#7]~*~*~[#7]";                  # bit79
$mc[64]  = "[#7]~*~*~*~[#7]";                # bit80
$mc[65]  = "[#16]~*(~*)~*";                  # bit81
$mc[66]  = "*~[CH2]~[!#6!#1!H0]";            # bit82
$mc[67]  = "[!#6!#1]1~*~*~*~*~1";            # bit83
$mc[68]  = "[NH2]";                          # bit84
$mc[69]  = "[#6]~[#7](~[#6])~[#6]";          # bit85
$mc[70]  = "[C;H2,H3][!#6!#1][C;H2,H3]";     # bit86
$mc[71]  = "[F,Cl,Br,I]!@*@*";               # bit87
$mc[72]  = "[#8]~*~*~*~[#8]";                # bit89
$mc[73]  = "[#8]~[#6](~[#7])~[#6]";          # bit92
$mc[74]  = "[!#6!#1]~[CH3]";                 # bit93
$mc[75]  = "[!#6!#1]~[#7]";                  # bit94
$mc[76]  = "[#7]~*~*~[#8]";                  # bit95
$mc[77]  = "*1~*~*~*~*~1";                   # bit96
$mc[78]  = "[#7]~*~*~*~[#8]";                # bit97
$mc[79]  = "[!#6!#1]1~*~*~*~*~*~1";          # bit98
$mc[80]  = "[#6]=[#6]";                      # bit99
$mc[81]  = "*~[CH2]~[#7]";                   # bit100
$mc[82]  = "[!#6!#1]~[#8]";                  # bit102
$mc[83]  = "[!#6!#1!H0]~*~[CH2]~*";          # bit104
$mc[84]  = "*@*(@*)@*";                      # bit105
$mc[85]  = "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]"; # bit106
$mc[86]  = "[F,Cl,Br,I]~*(~*)~*";            # bit107
$mc[87]  = "[CH3]~*~*~*~[CH2]~*";            # bit108
$mc[88]  = "*~[CH2]~[#8]";                   # bit109
$mc[89]  = "[#7]~[#6]~[#8]";                 # bit110
$mc[90]  = "[#7]~*~[CH2]~*";                 # bit111
$mc[91]  = "*~*(~*)(~*)~*";                  # bit112
$mc[92]  = "[#8]!:*:*";                      # bit113
$mc[93]  = "[CH3]~[CH2]~*";                  # bit114
$mc[94]  = "[CH3]~*~[CH2]~*";                # bit115
$mc[95]  = "[#7]~*~[#8]";                    # bit117
$mc[96]  = "[#7]=*";                         # bit119
$mc[97]  = "[!#6R]";                         # bit120
$mc[98]  = "[#7R]";                          # bit121
$mc[99]  = "*~[#7](~*)~*";                   # bit122
$mc[100] = "[#8]~[#6]~[#8]";                 # bit123
$mc[101] = "[!#6!#1]~[!#6!#1]";              # bit124
$mc[102] = "*!@[#8]!@*";                     # bit126
$mc[103] = "*@*!@[#8]";                      # bit127
$mc[104] = "[!#6!#1!H0]";                    # bit131
$mc[105] = "[#8]~*~[CH2]~*";                 # bit132
$mc[106] = "*@*!@[#7]";                      # bit133
$mc[107] = "[#7]!:*:*";                      # bit135
$mc[108] = "[#8]=*";                         # bit136
$mc[109] = "[!C!cR]";                        # bit137
$mc[110] = "[!#6!#1]~[CH2]~*";               # bit138
$mc[111] = "[O!H0]";                         # bit139
$mc[112] = "[#8]";                           # bit140
$mc[113] = "[CH3]";                          # bit141
$mc[114] = "[#7]";                           # bit142
$mc[115] = "*!:*:*!:*";                      # bit144
$mc[116] = "*1~*~*~*~*~*~1";                 # bit145
$mc[117] = "*~[!#6!#1](~*)~*";               # bit148
$mc[118] = "[CH3]";                          # bit149
$mc[119] = "*!@*@*!@*";                      # bit150
$mc[120] = "[#7!H0]";                        # bit151
$mc[121] = "[#8]~[#6](~[#6])~[#6]";          # bit152
$mc[122] = "[#6]=[#8]";                      # bit154
$mc[123] = "*!@[CH2]!@*";                    # bit155
$mc[124] = "[#7]~*(~*)~*";                   # bit156
$mc[125] = "[#6]-[#8]";                      # bit157
$mc[126] = "[#6]-[#7]";                      # bit158
$mc[127] = "a";                              # bit162
$mc[128] = "[R]";                            # bit165
$mc[129] = "[#3,#11,#19,#37,#55,#87]";       # bit35
$mc[130] = "[#4,#12,#20,#38,#56,#88]";       # bit10
$mc[131] = "[#5,#13,#31,#49,#81]";           # bit18
$mc[132] = "[#9]";                           # bit42
$mc[133] = "[#9,#17,#35,#53]";               # bit134
$mc[134] = "[#16]";                          # bit88
$mc[135] = "[#17]";                          # bit103 
$mc[136] = "[#21,#22,#39,#40,#72]";          # bit5
$mc[137] = "[#23,#24,#25,#41,#42,#73,#74,#75]";         # bit7
$mc[138] = "[#26,#27,#28,#44,#45,#46,#47,#76,#77,#78]"; # bit9
$mc[139] = "[#29,#30,#47,#48,#79,#80]";      # bit12
$mc[140] = "[#32,#33,#34,#50,#51,#52,#82,#83,#84]";     # bit3
$mc[141] = "[#35]";                          # bit46
$mc[142] = "[#53]";                          # bit27
$mc[143] = "[#57,#58,#59,#60,#61,#62,#63,#64,#65,#66,#67,#68"; 
$mc[143] = $mc[143].",#69,#70,#71]" ;        # bit6 
$mc[144] = "[#89,#90,#91,#92,#93,#94,#95,#96,#97,#98,#99"; 
$mc[144] = $mc[144].",#100,#101,#102,#103]"; # bit4
$mc[145] = "[#104]";                         # bit2
$mc[146] = "[r8]";                           # bit101 8-membered ring

$nmcmult  =  7 ; # number of multiple keys
$ncmult[1] = 3 ; # 147
$ncmult[2] = 4 ; # 148
$ncmult[3] = 2 ; # 149
$ncmult[4] = 2 ; # 150
$ncmult[5] = 2 ; # 151
$ncmult[6] = 3 ; # 152
$ncmult[7] = 2 ; # 153
$md[1][1] = "[!#6!#1!H0]~*~*~[CH2]~*";        # bit90
$md[1][2] = "[!#6!#1!H0R]1@[R]@[R]@[CH2R]1";  # bit90
$md[1][3] = "[!#6!#1!H0]~[R]1@[R]@[CH2R]1";   # bit90
$md[2][1] = "[!#6!#1!H0]~*~*~*~[CH2]~*";         # bit91
$md[2][2] = "[!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1"; # bit91
$md[2][3] = "[!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1";  # bit91
$md[2][4] = "[!#6!#1!H0]~*~[R]1@[R]@[CH2R]1";    # bit91
$md[3][1] = "[CH3]~*~*~[CH2]~*";              # bit116
$md[3][2] = "[CH3]~*1~*~[CH2]1";              # bit116
$md[4][1] = "*~[CH2]~[CH2]~*";                # bit118
$md[4][2] = "*1~[CH2]~[CH2]1";                # bit118
$md[5][1] = "*~[CH2]~*~*~*~[CH2]~*";             # bit128
$md[5][2] = "[R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1";   # bit128
$md[5][3] = "*~[CH2]~[R]1@[R]@[R]@[CH2R]1";      # bit128
$md[5][4] = "*~[CH2]~*~[R]1@[R]@[CH2R]1";        # bit128
$md[6][1] = "*~[CH2]~*~*~[CH2]~*";            # bit129
$md[6][2] = "[R]1@[CH2]@[R]@[R]@[CH2R]1";     # bit129
$md[6][3] = "*~[CH2]~[R]1@[R]@[CH2R]1";       # bit129
$md[7][1] = "*~[CH2]~[CH2]~*";                # bit147
$md[7][2] = "[R]1@[CH2R]@[CH2R]1";            # bit147

#  MACCS keys that set additional bits on multiple occurence
#
#  obgrep -t2 "pattern" file  #matches 2 occurrences of pattern
#
$nmcmbits = 8; # number of pattern that set multiple bits
$mcmbit[1][1] = $mc[101] ;    # bit124
$mcmbit[1][2] = 154 ; # >1 occurrence
$mcmbit[1][3] = 0   ; # third occurrence
$mcmbit[1][4] = 0   ; # third occurrence
$mcmbit[2][1] = $mc[103] ;    # bit127
$mcmbit[2][2] = 155 ; # >1 occurrence
$mcmbit[2][3] = 0   ; 
$mcmbit[2][4] = 0   ; 
$mcmbit[3][1] = $mc[110] ;    # bit138
$mcmbit[3][2] = 156 ; # >1 occurrence
$mcmbit[3][3] = 0   ; 
$mcmbit[3][4] = 0   ; 
$mcmbit[4][1] = $mc[113] ;    # bit141
$mcmbit[4][2] = 157 ; # >1 occurrence
$mcmbit[4][3] = 0   ; 
$mcmbit[4][4] = 0   ; 
$mcmbit[5][1] = $mc[116] ;    # bit145
$mcmbit[5][2] = 158 ; # >1 occurrence
$mcmbit[5][3] = 0   ; 
$mcmbit[5][4] = 0   ; 
$mcmbit[6][1] = $mc[118] ;    # bit149
$mcmbit[6][2] = 159 ; # >1 occurrence
$mcmbit[6][3] = 0   ; 
$mcmbit[6][4] = 0   ; 
$mcmbit[7][1] = "[ar6]1~[ar6]~[ar6]~[ar6]~[ar6]~[ar6]1" ;  # bit125 >1 aromatic rings 
$mcmbit[7][2] = 160 ; # >1 occurrence
$mcmbit[7][3] = 0   ; 
$mcmbit[7][4] = 0   ; 
$mcmbit[8][1] = $mc[112] ;    # bit140
$mcmbit[8][2] = 161 ; # >1 occurrence
$mcmbit[8][3] = 162 ; # >2 
$mcmbit[8][4] = 163 ; # >3 

#
#
#  MACCS keys that use $ have to be treated differently when
#  using open babel, otherwise occurrences in rings are not detected
#  bit 90, 91, 116, 118, 128, 129, 147 
#  A = *
#  Q = [!#6!#1!H0]
#  bit 90 QHAACH2A    
#  bit 91 QHAAACH2A   
#  bit 116  CH3AACH2A         [CH3]~*~*~[CH2]~*
#  bit 118  ACH2CH2A > 1      *~[CH2]~[CH2]~*  
#  bit 128  ACH2AAACH2A
#  bit 129  ACH2AACH2A
#  bit 147  ACH2CH2A
#  
#  bit 124  set bit 124 , multiple occurrences set bit 130
#  bit 127  set bit 127 ,    # set bit 143
#  bit 131  only set if occurrences > 1
#  bit 136  only set if occurrences > 1
#  bit 138  set bit 153,  multiple occurrences set bit 138
#  bit 140  set bit 164, > 1 bit 159, > 2 bit 146, > 3 bit 140
#  bit 141  set bit 161, > 1 bit 161
#  bit 145  set bit 163, > 1 bit 145
#  bit 149  set bit 160, > 1 bit 149
#  bit 166  set if more than one fragement is present (.)
#  bit 125  number of aromatic rings > 1


#  PubChem fingerprint version 1.3 (May 1st, 2009) as given
#  on the PubChem ftp server

$maxpekeys = 753  ; # maximum number of "easy" PubChem keys

# define MACCS keys as SMARTS  (>= 1 occurrences)
#
$pc[1]   = "[#3]";     # bit4   Li  section 1 elements        
$pc[2]   = "[#5]";     # bit6   B     
$pc[3]   = "[#7]";     # bit14  N       
$pc[4]   = "[#8]";     # bit18  O      
$pc[5]   = "[#9]";     # bit23  F
$pc[6]   = "[#11]";    # bit26  Na
$pc[7]   = "[#14]";    # bit28  Si         
$pc[8]   = "[#15]";    # bit30  P         
$pc[9]   = "[#16]";    # bit33  S         
$pc[10]  = "[#17]";    # bit37  Cl         
$pc[11]  = "[#19]";    # bit41  K        
$pc[12]  = "[#35]";    # bit43  Br         
$pc[13]  = "[#53]";    # bit46  I         
$pc[14]  = "[#4]";     # bit49  Be         
$pc[15]  = "[#12]";    # bit50  Mg         
$pc[16]  = "[#13]";    # bit51  Al         
$pc[17]  = "[#20]";    # bit52  Ca         
$pc[18]  = "[#21]";    # bit53  Sc         
$pc[19]  = "[#22]";    # bit54  Ti         
$pc[20]  = "[#23]";    # bit55  V    
$pc[21]  = "[#24]";    # bit56  Cr         
$pc[22]  = "[#25]";    # bit57  Mn         
$pc[23]  = "[#26]";    # bit58  Fe         
$pc[24]  = "[#27]";    # bit59  Co         
$pc[25]  = "[#28]";    # bit60  Ni        
$pc[26]  = "[#29]";    # bit61  Cu         
$pc[27]  = "[#30]";    # bit62  Zn         
$pc[28]  = "[#31]";    # bit63  Ga         
$pc[29]  = "[#32]";    # bit64  Ge       
$pc[30]  = "[#33]";    # bit65  As         
$pc[31]  = "[#34]";    # bit66  Se
$pc[32]  = "[#36]";    # bit67  Kr         
$pc[33]  = "[#37]";    # bit68  Rb         
$pc[34]  = "[#38]";    # bit69  Sr
$pc[35]  = "[#39]";    # bit70  Y         
$pc[36]  = "[#40]";    # bit71  Zr         
$pc[37]  = "[#41]";    # bit72  Nb         
$pc[38]  = "[#42]";    # bit73  Mo         
$pc[39]  = "[#44]";    # bit74  Ru         
$pc[40]  = "[#45]";    # bit75  Rh         
$pc[41]  = "[#46]";    # bit76  Pd     
$pc[42]  = "[#47]";    # bit77  Ag     
$pc[43]  = "[#48]";    # bit78  Cd     
$pc[44]  = "[#49]";    # bit79  In     
$pc[45]  = "[#50]";    # bit80  Sn     
$pc[46]  = "[#51]";    # bit81  Sb     
$pc[47]  = "[#52]";    # bit82  Te
$pc[48]  = "[#54]";    # bit83  Xe
$pc[49]  = "[#55]";    # bit84  Cs  
$pc[50]  = "[#56]";    # bit85  Ba
$pc[51]  = "[#71]";    # bit86  Lu
$pc[52]  = "[#72]";    # bit87  Hf        
$pc[53]  = "[#73]";    # bit88  Ta         
$pc[54]  = "[#74]";    # bit89  W         
$pc[55]  = "[#75]";    # bit90  Re         
$pc[56]  = "[#76]";    # bit91  Os         
$pc[57]  = "[#77]";    # bit92  Ir         
$pc[58]  = "[#78]";    # bit93  Pt         
$pc[59]  = "[#79]";    # bit94  Au         
$pc[60]  = "[#80]";    # bit95  Hg         
$pc[61]  = "[#81]";    # bit96  Tl         
$pc[62]  = "[#82]";    # bit97  Pb     
$pc[63]  = "[#83]";    # bit98  Bi         
$pc[64]  = "[#57]";    # bit99  La         
$pc[65]  = "[#58]";    # bit100 Ce         
$pc[66]  = "[#59]";    # bit101 Pr         
$pc[67]  = "[#60]";    # bit102 Ne         
$pc[68]  = "[#61]";    # bit103 Pm         
$pc[69]  = "[#62]";    # bit104 Sm         
$pc[70]  = "[#63]";    # bit105 Eu         
$pc[71]  = "[#64]";    # bit106 Gd         
$pc[72]  = "[#65]";    # bit107 Tb         
$pc[73]  = "[#66]";    # bit108 Dy         
$pc[74]  = "[#67]";    # bit109 Ho         
$pc[75]  = "[#68]";    # bit110 Er         
$pc[76]  = "[#69]";    # bit111 Tm         
$pc[77]  = "[#70]";    # bit112 Yb
$pc[78]  = "[#43]";    # bit113 Tc         
$pc[79]  = "[#92]";    # bit114 U         
$pc[80]  = "*1~*~*1";               # bit115 3-ring   section 2 SSSM         
$pc[81]  = "[#6]1~[#6]~[#6]1";      # bit116 3-ring carbons only   
$pc[82]  = "[#6]1~[#7]~[#6,#7]1";   # bit117 3-ring N containing
$pc[83]  = "[#6]1~[!#6!#1]~[#6,#7]1";    # bit118 3-ring hetero
$pc[84]  = "C1=CC1";                # bit119 unsaturated 3-ring C only
$pc[85]  = "N1=[C,N]~C1";           # bit120 unsaturated 3 ring N         
$pc[86]  = "C1=[!C!#1]~[C,N]1";     # bit121 unsaturated 3 ring hetero         
$pc[87]  = "*1~*~*~*1";             # bit129 4-ring
$pc[88]  = "[#6]1~[#6]~[#6]~[#6]1"; # bit130 4-ring carbon only         
$pc[89]  = "[#7]1~[#6,#7]~[#6,#7]~[#6]1"; # bit131 4-ring N containing         
$pc[90]  = "[!#6!#1]1~*~*~*1";      # bit132 4-ring hetero         
$pc[91]  = "C1!-C~C~C1";            # bit133 unsaturated 4-ring C only         
$pc[92]  = "C1~N!-[C,N]~[C,N]1";    # bit134 unsaturated 4-ring N         
$pc[93]  = "*1~[!C]~*!-*1";         # bit135 unsaturated 4-ring hetero         
$pc[94]  = "*1~*~*~*~*1";           # bit143 5-ring   
$pc[95]  = "[#6]1~[#6]~[#6]~[#6]~[#6]1";  # bit144 5-ring C only
$pc[96]  = "[#7]1~[#6,#7]~[#6,#7]~[#6,#7]~[#6]1"; # bit145 5-ring with N 
$pc[97]  = "[!#6!#1]1~*~*~*~*1";    # bit146 5-ring hetero         
$pc[98]  = "C1!-C~C~C~C1";          # bit147 5-ring unsat.C only         
$pc[99]  = "C1=[C,N]~[C,N]~[C,N]~[C,N]1"; # bit148 5-ring unsat. N         
$pc[100] = "*1~[!C]~*!-*~*1";       # bit149 5-ring unsat. hetero         
$pc[101] = "*1~*~*~*~*~*1";         # bit178 6-ring
$pc[102] = "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1"; # bit179 6-ring C only         
$pc[103] = "[#7]1~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]~[#6]1"; # bit180 6-ring N 
$pc[104] = "[!#6!#1]1~*~*~*~*~*1";  # bit181 6-ring hetero
$pc[105] = "C1!-C~C~C~C~C1";        # bit182 6-ring unsat. C only
$pc[106] = "C1=[C,N]~[C,N]~[C,N]~[C,N]~[C,N]1"; # bit183 6-ring unsat. N
$pc[107] = "C1!-[C,N]~*~*~*~*1";    # bit184 6-ring unsat. hetero
$pc[108] = "*1~*~*~*~*~*~*1";       # bit213 7-ring         
$pc[109] = "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]1"; # bit214 7-ring C only
$pc[110] = "[#7]1~*~*~*~*~*~*1";    # bit215 7-ring N
$pc[111] = "[!#6!#1]1-*-*-*-*-*-*1";  # bit216 7-ring hetero
$pc[112] = "C1!-C~C~C~C~C~C1";      # bit217 7-ring unsat. C only
$pc[113] = "C1=[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]1"; # bit218 7-ring unsat. N
$pc[114] = "C1!-[C,N]~*~*~*~*~*1";  # bit219 7-ring unsat. hetero
$pc[115] = "[r8]";                  # bit227 8-ring
$pc[116] = "C1CCCCCCC1";            # bit228 8-ring C only
$pc[117] = "[#7]1-*-*-*-*-*-*-*1";  # bit229 8-ring N
$pc[118] = "[!#6!#1]1-*-*-*-*-*-*-*1"; # bit230 8-ring hetero
$pc[119] = "C1!-C~C~C~C~C~C~C1";    # bit231 8-ring unsat. C only
$pc[120] = "C1=[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]1"; # bit232 8-ring unsat. N
$pc[121] = "C1!-[C,N]~*~*~*~*~*~*1";   # bit233 8-ring unsat. hetero
$pc[122] = "[r9]";                  # bit241 9-ring
$pc[123] = "C1CCCCCCCC1";           # bit242 9-ring C only
$pc[124] = "[#7]1-*-*-*-*-*-*-*-*1";   # bit243 9-ring N
$pc[125] = "[!#6!#1]1-*-*-*-*-*-*-*-*1"; # bit244 9-ring hetero
$pc[126] = "C1!-C~C~C~C~C~C~C~C1";  # bit245 9-ring unsat. C only
$pc[127] = "C1=[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]1"; # bit246 9-ring unsat. N
$pc[128] = "C1!-[C,N]~*~*~*~*~*~*~*1"; # bit247 9-ring unsat. hetero
$pc[129] = "[r10]"; # bit248 10-ring
$pc[130] = "C1CCCCCCCC1";           # bit249 10-ring C only
$pc[131] = "[#7]1-*-*-*-*-*-*-*-*-*1"; # bit250 10-ring N
$pc[132] = "[!#6!#1]1-*-*-*-*-*-*-*-*-*1"; # bit251 10-ring hetero
$pc[133] = "C1!-C~C~C~C~C~C~C~C~C1";   # bit252 10-ring unsat. C only
$pc[134] = "C1=[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]~[C,N]1"; # bit253 10-ring unsat. N
$pc[135] = "C1!-[C,N]~*~*~*~*~*~*~*~*1"; # bit254 10-ring unsat. hetero
$pc[136] = "[#3]~[#1]";      # bit263  Li-H section 3 atom pairs
$pc[137] = "[#3]~[#3]";      # bit264  Li-Li
$pc[138] = "[#3]~[#5]";      # bit265  Li-B         
$pc[139] = "[#3]~[#6]";      # bit266  Li-C         
$pc[140] = "[#3]~[#8]";      # bit267  Li-O
$pc[141] = "[#3]~[#9]";      # bit268  Li-F         
$pc[142] = "[#3]~[#15]";     # bit269  Li-P 
$pc[143] = "[#3]~[#16]";     # bit270  Li-S         
$pc[144] = "[#3]~[#17]";     # bit271  Li-Cl         
$pc[145] = "[#5]~[#1]";      # bit272  B-H
$pc[146] = "[#5]~[#1]";      # bit273  B-B         
$pc[147] = "[#5]~[#6]";      # bit274  B-C
$pc[148] = "[#5]~[#7]";      # bit275  B-N
$pc[149] = "[#5]~[#8]";      # bit276  B-O
$pc[150] = "[#5]~[#9]";      # bit277  B-F
$pc[151] = "[#5]~[#14]";     # bit278  B-Si         
$pc[152] = "[#5]~[#15]";     # bit279  B-P
$pc[153] = "[#5]~[#16]";     # bit280  B-S
$pc[154] = "[#5]~[#17]";     # bit281  B-Cl         
$pc[155] = "[#5]~[#35]";     # bit282  B-Br         
$pc[156] = "[#6]~[#1]";      # bit283  C-H
$pc[157] = "[#6]~[#1]";      # bit284  C~C
$pc[158] = "[#6]~[#5]";      # bit285  C~N
$pc[159] = "[#6]~[#8]";      # bit286  C~O
$pc[160] = "[#6]~[#9]";      # bit287  C-F         
$pc[161] = "[#6]~[#11]";     # bit288  C-Na
$pc[162] = "[#6]~[#12]";     # bit289  C-Mg
$pc[163] = "[#6]~[#13]";     # bit290  C-Al
$pc[164] = "[#6]~[#14]";     # bit291  C~Si
$pc[165] = "[#6]~[#15]";     # bit292  C~P
$pc[166] = "[#6]~[#16]";     # bit293  C~S
$pc[167] = "[#6]~[#17]";     # bit294  C-Cl
$pc[168] = "[#6]~[#33]";     # bit295  C~As         
$pc[169] = "[#6]~[#34]";     # bit296  C~Se         
$pc[170] = "[#6]~[#35]";     # bit297  C-Br         
$pc[171] = "[#6]~[#53]";     # bit298  C-I         
$pc[172] = "[#7]~[#1]";      # bit299  N-H         
$pc[173] = "[#7]~[#7]";      # bit300  N~N 
$pc[174] = "[#7]~[#8]";      # bit301  N~O
$pc[175] = "[#7]~[#9]";      # bit302  N-F
$pc[176] = "[#7]~[#14]";     # bit303  N~Si
$pc[177] = "[#7]~[#15]";     # bit304  N~P
$pc[178] = "[#7]~[#16]";     # bit305  N~S         
$pc[179] = "[#7]~[#17]";     # bit306  N-Cl
$pc[180] = "[#7]~[#35]";     # bit307  N-Br
$pc[181] = "[#8]~[#1]";      # bit308  O-H         
$pc[182] = "[#8]~[#8]";      # bit309  O-O         
$pc[183] = "[#8]~[#12]";     # bit310  O-Mg
$pc[184] = "[#8]~[#11]";     # bit311  O-Na
$pc[185] = "[#8]~[#13]";     # bit312  O-Al
$pc[186] = "[#8]~[#14]";     # bit313  O~Si         
$pc[187] = "[#8]~[#15]";     # bit314  O~P
$pc[188] = "[#8]~[#19]";     # bit315  O-K
$pc[189] = "[#9]~[#15]";     # bit316  F-P         
$pc[190] = "[#9]~[#16]";     # bit317  F-S
$pc[191] = "[#13]~[#1]";     # bit318  Al-H
$pc[192] = "[#13]~[#17]";    # bit319  Al-Cl
$pc[193] = "[#14]~[#1]";     # bit320  Si-H
$pc[194] = "[#14]~[#14]";    # bit321  Si~Si
$pc[195] = "[#14]~[#17]";    # bit322  Si-Cl
$pc[196] = "[#15]~[#1]";     # bit323  P-H
$pc[197] = "[#15]~[#15]";    # bit324  P~P
$pc[198] = "[#33]~[#1]";     # bit325  As-H         
$pc[199] = "[#33]~[#33]";    # bit326  As~As
$pc[200] = "[#6](~[#35])(~[#6])";      # bit327  Section 4: simple atom nearest neighbor  C(~Br)(~C)
$pc[201] = "[#6](Br)(~[#6])(~[#6])";   # bit328  C(~Br)(~C)(~C)
$pc[202] = "[#6](Br)(~[#1])";          # bit329  C(~Br)(~H)
$pc[203] = "[#6](Br)c";                # bit330  C(~Br)(:C)
$pc[204] = "[#6](Br)n";                # bit331  C(~Br)(:N) 
$pc[205] = "[#6](~[#6])(~[#6])";       # bit332  C(~C)(~C)
$pc[206] = "[#6](~[#6])(~[#6])~[#6]";  # bit333  C(~C)(~C)(~C)
$pc[207] = "[#6](~[#6])(~[#6])(~[#6])~[#6]"; # bit334  C(~C)(~C)(~C)(~C)
$pc[208] = "[#6](~[#6])(~[#6])(~[#6])~[#1]"; # bit335  C(~C)(~C)(~C)(~H)
$pc[209] = "[#6](~[#6])(~[#6])(~[#6])~[#7]"; # bit336  C(~C)(~C)(~C)(~N)
$pc[210] = "[#6](~[#6])(~[#6])(~[#6])~[#8]"; # bit337  C(~C)(~C)(~C)(~O)
$pc[211] = "[#6](~[#6])(~[#6])(~[#1])~[#7]"; # bit338  C(~C)(~C)(~H)(~N)
$pc[212] = "[#6](~[#6])(~[#6])(~[#1])~[#8]"; # bit339  C(~C)(~C)(~H)(~O)
$pc[213] = "[#6](~[#6])(~[#6])~[#7]";        # bit340  C(~C)(~C)(~N) 
$pc[214] = "[#6](~[#6])(~[#6])~[#8]";        # bit341  C(~C)(~C)(~O)
$pc[215] = "[#6](~[#6])(~[#17])";            # bit342  C(~C)(~Cl)
$pc[216] = "[#6](~[#6])(~[#17])~[#1]";       # bit343  C(~C)(~Cl)(~H)
$pc[217] = "[#6](~[#6])~[#1]";               # bit344  C(~C)(~H)
$pc[218] = "[#6](~[#6])(~[#1])~[#7]";        # bit345  C(~C)(~H)(~N)
$pc[219] = "[#6](~[#6])(~[#1])~[#8]";        # bit346  C(~C)(~H)(~O)
$pc[220] = "[#6](~[#6])(~[#1])(~[#8])~[#8]"; # bit347  C(~C)(~H)(~O)(~O)
$pc[221] = "[#6](~[#6])(~[#1])~[#15]"; # bit348  C(~C)(~H)(~P)
$pc[222] = "[#6](~[#6])(~[#1])~[#16]"; # bit349  C(~C)(~H)(~S)
$pc[223] = "[#6](~[#6])~[#53]";        # bit350  C(~C)(~I)
$pc[224] = "[#6](~[#6])~[#7]";         # bit351  C(~C)(~N) 
$pc[225] = "[#6](~[#6])~[#8]";         # bit352  C(~C)(~O)
$pc[226] = "[#6](~[#6])~[#16]";        # bit353  C(~C)(~S)          
$pc[227] = "[#6](~[#6])~[#14]";        # bit354  C(~C)(~Si)
$pc[228] = "[#6](~[#6])c";             # bit355  C(~C)(:C)
$pc[229] = "[#6](~[#6])(c)c";          # bit356  C(~C)(:C)(:C)
$pc[230] = "[#6](~[#6])(c)n";          # bit357  C(~C)(:C)(:N)
$pc[231] = "[#6](~[#6])n";             # bit358  C(~C)(:N)
$pc[232] = "[#6](~[#6])(n)n";          # bit359  C(~C)(:N)(:N)
$pc[233] = "[#6](~[#17])~[#17]";       # bit360  C(~Cl)(~Cl)
$pc[234] = "[#6](~[#17])~[#1]";        # bit361  C(~Cl)(~H)
$pc[235] = "[#6](~[#17])c";            # bit362  C(~Cl)(:C)
$pc[236] = "[#6](~[#9])~[#9]";         # bit363  C(~F)(~F)
$pc[237] = "[#6](~[#9])c";             # bit364  C(~F)(:C)
$pc[238] = "[#6](~[#1])~[#7]";         # bit365  C(~H)(~N)
$pc[239] = "[#6](~[#1])~[#8]";         # bit366  C(~H)(~O) 
$pc[240] = "[#6](~[#1])(~[#8])~[#8]";  # bit367  C(~H)(~O)(~O)
$pc[241] = "[#6](~[#1])~[#16]";        # bit368  C(~H)(~S)
$pc[242] = "[#6](~[#1])~[#14]";        # bit369  C(~H)(~Si)
$pc[243] = "[#6](~[#1])c";             # bit370  C(~H)(:C)
$pc[244] = "[#6](~[#1])(c)c";          # bit371  C(~H)(:C)(:C)
$pc[245] = "[#6](~[#1])(c)n";          # bit372  C(~H)(:C)(:N)
$pc[246] = "[#6](~[#1])n";             # bit373  C(~H)(:N)
$pc[247] = "[#6](~[#1])(~[#1])~[#1]";  # bit374  C(~H)(~H)(~H)
$pc[248] = "[#6](~[#7])~[#7]";         # bit375  C(~N)(~N)
$pc[249] = "[#6](~[#7])c";             # bit376  C(~N)(:C)
$pc[250] = "[#6](~[#7])(c)c";          # bit377  C(~N)(:C)(:C)
$pc[251] = "[#6](~[#7])(c)n";          # bit378  C(~N)(:C)(:N)
$pc[252] = "[#6](~[#7])n";             # bit379  C(~N)(:N)
$pc[253] = "[#6](~[#8])~[#8]";         # bit380  C(~O)(~O)
$pc[254] = "[#6](~[#8])c";             # bit381  C(~O)(:C)         
$pc[255] = "[#6](~[#8])(c)c";          # bit382  C(~O)(:C)(:C)
$pc[256] = "[#6](~[#16])c";            # bit383  C(~S)(:C)
$pc[257] = "[#6](c)c";                 # bit384  C(:C)(:C)
$pc[258] = "[#6](c)(c)c";              # bit385  C(:C)(:C)(:C)
$pc[259] = "[#6](c)(c)n";              # bit386  C(:C)(:C)(:N)
$pc[260] = "[#6](c)n";                 # bit387  C(:C)(:N)
$pc[261] = "[#6](c)(n)n";              # bit388  C(:C)(:N)(:N)
$pc[262] = "[#6](n)n";                 # bit389  C(:N)(:N)
$pc[263] = "[#7](~[#6])~[#6]";         # bit390  N(~C)(~C)
$pc[264] = "[#7](~[#6])(~[#6])~[#6]";  # bit391  N(~C)(~C)(~C)
$pc[265] = "[#7](~[#6])(~[#6])~[#1]";  # bit392  N(~C)(~C)(~H)
$pc[266] = "[#7](~[#6])~[#1]";         # bit393  N(~C)(~H)
$pc[267] = "[#7](~[#6])(~[#1])~[#7]";  # bit394  N(~C)(~H)(~N)
$pc[268] = "[#7](~[#6])~[#8]";         # bit395  N(~C)(~O)
$pc[269] = "[#7](~[#6])c";             # bit396  N(~C)(:C) 
$pc[270] = "[#7](~[#6])(c)c";          # bit397  N(~C)(:C)(:C) 
$pc[271] = "[#7](~[#1])~[#7]";         # bit398  N(~H)(~N)
$pc[272] = "[#7](~[#1])c";             # bit399  N(~H)(:C)
$pc[273] = "[#7](~[#1])(c)c";          # bit400  N(~H)(:C)(:C)
$pc[274] = "[#7](~[#8])~[#8]";         # bit401  N(~O)(~O)
$pc[275] = "[#7](~[#8])o";             # bit402  N(~O)(:O)
$pc[276] = "[#7](c)c";                 # bit403  N(:C)(:C)
$pc[277] = "[#7](c)(c)c";              # bit404  N(:C)(:C)(:C)
$pc[278] = "[#8](~[#6])~[#6]";         # bit405  O(~C)(~C)
$pc[279] = "[#8](~[#6])~[#1]";         # bit406  O(~C)(~H)
$pc[280] = "[#8](~[#6])~[#15]";        # bit407  O(~C)(~P)
$pc[281] = "[#8](~[#1])~[#16]";        # bit408  O(~H)(~S)
$pc[282] = "[#8](c)c";                 # bit409  O(:C)(:C)
$pc[283] = "[#15](~[#6])~[#6]";        # bit410  P(~C)(~C)
$pc[284] = "[#15](~[#8])~[#8]";        # bit411  P(~O)(~O)
$pc[285] = "[#16](~[#6])~[#6]";        # bit412  S(~C)(~C)
$pc[286] = "[#16](~[#6])~[#1]";        # bit413  S(~C)(~H)
$pc[287] = "[#16](~[#6])~[#8]";        # bit414  S(~C)(~O)
$pc[288] = "[#14](~[#6])~[#6]";        # bit415  Si(~C)(~C)
$pc[289] = "[#6]=[#6]";         # bit416  Section 5: Detailed neighborhood  C=C
$pc[290] = "[#6]#[#6]";         # bit417  C#C       
$pc[291] = "[#6]=[#7]";         # bit418  C=N
$pc[292] = "[#6]#[#7]";         # bit419  C#N
$pc[293] = "[#6]=[#8]";         # bit420  C=O
$pc[294] = "[#6]=[#16]";        # bit421  C=S
$pc[295] = "[#7]=[#7]";         # bit422  N=N
$pc[296] = "[#7]=[#8]";         # bit423  N=O
$pc[297] = "[#7]=[#15]";        # bit424  N=P
$pc[298] = "[#15]=[#8]";        # bit425  P=O         
$pc[299] = "[#15]=[#15]";       # bit426  P=P
$pc[300] = "[#6](#[#6])[#6]";   # bit427  C(#C)(-C)
$pc[301] = "[#6](#[#6])[#1]";   # bit428  C(#C)(-H)
$pc[302] = "[#6](#[#7])[#6]";   # bit429  C(#N)(-C)
$pc[303] = "[#6]([#6])([#6])=[#6]";   # bit430  C(-C)(-C)(=C)
$pc[304] = "[#6]([#6])([#6])=[#7]";   # bit431  C(-C)(-C)(=N)
$pc[305] = "[#6]([#6])([#6])=[#8]";   # bit432  C(-C)(-C)(=O)
$pc[306] = "[#6]([#6])([#17])=[#8]";  # bit433  C(-C)(-Cl)(=O)
$pc[307] = "[#6]([#6])([#1])=[#6]";   # bit434  C(-C)(-H)(=C)
$pc[308] = "[#6]([#6])([#1])=[#7]";   # bit435  C(-C)(-H)(=N)
$pc[309] = "[#6]([#6])([#1])=[#8]";   # bit436  C(-C)(-H)(=O)
$pc[310] = "[#6]([#6])([#7])=[#6]";   # bit437  C(-C)(-N)(=C)
$pc[311] = "[#6]([#6])([#7])=[#7]";   # bit438  C(-C)(-N)(=N)
$pc[312] = "[#6]([#6])([#7])=[#8]";   # bit439  C(-C)(-N)(=O)
$pc[313] = "[#6]([#6])([#8])=[#8]";   # bit440  C(-C)(-O)(=O)
$pc[314] = "[#6]([#6])=[#6]";         # bit441  C(-C)(=C)          
$pc[315] = "[#6]([#6])=[#7]";         # bit442  C(-C)(=N)
$pc[316] = "[#6]([#6])=[#8]";         # bit443  C(-C)(=O)
$pc[317] = "[#6]([#17])=[#8]";        # bit444  C(-Cl)(=O)     
$pc[318] = "[#6]([#1])([#7])=[#6]";   # bit445  C(-H)(-N)(=C)
$pc[319] = "[#6]([#1])=[#6]";         # bit446  C(-H)(=C)
$pc[320] = "[#6]([#1])=[#7]";         # bit447  C(-H)(=N)
$pc[321] = "[#6]([#1])=[#8]";         # bit448  C(-H)(=O)
$pc[322] = "[#6]([#7])=[#6]";         # bit449  C(-N)(=C)
$pc[323] = "[#6]([#7])=[#7]";         # bit450  C(-N)(=N)
$pc[324] = "[#6]([#7])=[#8]";         # bit451  C(-N)(=O)
$pc[325] = "[#6]([#8])=[#8]";         # bit452  C(-O)(=O)
$pc[326] = "[#7]([#6])=[#6]";         # bit453  N(-C)(=C)
$pc[327] = "[#7]([#6])=[#8]";         # bit454  N(-C)(=O)  
$pc[328] = "[#7]([#8])=[#6]";         # bit455  N(-O)(=O)
$pc[329] = "[#15]([#8])=[#8]";        # bit456  P(-O)(=O)
$pc[330] = "[#16]([#6])=[#8]";        # bit457  S(-C)(=O)
$pc[331] = "[#16]([#8])=[#8]";        # bit458  S(-O)(=O)
$pc[332] = "[#16]([#8])=[#8]";        # bit459  S(=O)(=O)
$pc[333] = "CCC#C";        # bit460  Section 6: simple SMARTS patterns C-C-C#C
$pc[334] = "OCC=N";        # bit461  O-C-C=N           
$pc[335] = "OCC=O";        # bit462  O-C-C=O 
$pc[336] = "nc[#16][#1]";  # bit463  N:C-S-[#1]
$pc[337] = "NCC=C";        # bit464  N-C-C=C 
$pc[338] = "O=SCC";        # bit465  O=S-C-C
$pc[339] = "N#CC=C";       # bit466  N#C-C=C
$pc[340] = "C=NNC";        # bit467  C=N-N-C
$pc[341] = "O=SCN";        # bit468  O=S-C-N
$pc[342] = "SScc";         # bit469  S-S-C:C
$pc[343] = "ccC=C";        # bit470  C:C-C=C
$pc[344] = "sccc";         # bit471  S:C:C:C         
$pc[345] = "cncC";         # bit472  C:N:C-C 
$pc[346] = "Scnc";         # bit473  S-C:N:C
$pc[347] = "sccn";         # bit474  S:C:C:N
$pc[348] = "SC=NC";        # bit475  S-C=N-C         
$pc[349] = "COOC=C";       # bit476  C-O-C=C 
$pc[350] = "NNcc";         # bit477  N-N-C:C
$pc[351] = "SC=N[#1]";     # bit478  S-C=N-[#1]
$pc[352] = "SCSC";         # bit479  S-C-S-C
$pc[353] = "cscC";         # bit480  C:S:C-C
$pc[354] = "OScc";         # bit481  O-S-C:C
$pc[355] = "cn-cc";        # bit482  C:N-C:C
$pc[356] = "NScc";         # bit483  N-S-C:C
$pc[357] = "Ncnc";         # bit484  N-C:N:C
$pc[358] = "nccn";         # bit485  N:C:C:N
$pc[359] = "Ncnn";         # bit486  N-C:N:N;
$pc[360] = "NC=NC";        # bit487  N-C=N-C
$pc[361] = "NC=N[#1]";     # bit488  N-C=N-[#1]
$pc[362] = "NCSC";         # bit489  N-C-S-C
$pc[363] = "CCC=C";        # bit490  C-C-C=C
$pc[364] = "Cnc[#1]";      # bit491  C-N:C-[#1]
$pc[365] = "Ncoc";         # bit492  N-C:O:C      
$pc[366] = "O=Ccc";        # bit493  O=C-C:C
$pc[367] = "O=Ccn";        # bit494  O=C-C:N
$pc[368] = "CNcc";         # bit495  C-N-C:C
$pc[369] = "nnC[#1]";      # bit496  N:N-C-[#1]
$pc[370] = "Occn";         # bit497  O-C:C:N          
$pc[371] = "OC=CC";        # bit498  O-C=C-C
$pc[372] = "Nccn";         # bit499  N-C:C:N
$pc[373] = "CScc";         # bit500  C-S-C:C
$pc[374] = "[#17]ccC";     # bit501  Cl-C:C-C
$pc[375] = "NC=C[#1]";     # bit502  N-C=C-[#1]
$pc[376] = "[#17]cc[#1]";  # bit503  Cl-C:C-[#1]
$pc[377] = "ncnC";         # bit504  N:C:N-C
$pc[378] = "[#17]ccO";     # bit505  Cl-C:C-O
$pc[379] = "Cnnc";         # bit506  C-C:N:C
$pc[380] = "CCSC";         # bit507  C-C-S-C
$pc[381] = "S=CNC";        # bit508  S=C-N-C       
$pc[382] = "[#35]ccC";     # bit509  Br-C:C-C
$pc[383] = "[#1]NN[#1]";   # bit510  [#1]-N-N-[#1]
$pc[384] = "S=CN[#1]";     # bit511  S=C-N-[#1]
$pc[385] = "C[#33]O[#1]";  # bit512  C-[As]-O-[#1]
$pc[386] = "scc[#1]";      # bit513  S:C:C-[#1]
$pc[387] = "ONCC";         # bit514  O-N-C-C
$pc[388] = "NNCC";         # bit515  N-N-C-C
$pc[389] = "[#1]C=C[#1]";  # bit516  [#1]-C=C-[#1]
$pc[390] = "NNCN";         # bit517  N-N-C-N
$pc[391] = "O=CNN";        # bit518  O=C-N-N
$pc[392] = "N=CNC";        # bit519  N=C-N-C
$pc[393] = "C=Ccc";        # bit520  C=C-C:C
$pc[394] = "cnC[#1]";      # bit521  C:N-C-[#1]
$pc[395] = "CNN[#1]";      # bit522  C-N-N-[#1]
$pc[396] = "nccC";         # bit523  N:C:C-C
$pc[397] = "CC=CC";        # bit524  C-C=C-C
$pc[398] = "[#33]cc[#1]";  # bit525  [As]-C:C-[#1]
$pc[399] = "[#17]cc[#17]"; # bit526  Cl-C:C-Cl
$pc[400] = "ccn[#1]";      # bit527  C:C:N-[#1]
$pc[401] = "[#1]NC[#1]";   # bit528  [#1]-N-C-[#1]
$pc[402] = "[#17]CC[#17]"; # bit529  Cl-C-C-Cl
$pc[403] = "nc-cc";        # bit530  N:C-C:C
$pc[404] = "SccC";         # bit531  S-C:C-C
$pc[405] = "Scc[#1]";      # bit532  S-C:C-[#1]
$pc[406] = "SccN";         # bit533  S-C:C-N
$pc[407] = "SccO";         # bit534  S-C:C-O
$pc[408] = "O=CCC";        # bit535  O=C-C-C
$pc[409] = "O=CCN";        # bit536  O=C-C-N
$pc[410] = "O=CCO";        # bit537  O=C-C-O
$pc[411] = "N=CCC";        # bit538  N=C-C-C       
$pc[412] = "N=CC[#1]";     # bit539  N=C-C-[#1]
$pc[413] = "CNC[#1]";      # bit540  C-N-C-[#1]
$pc[414] = "OccC";         # bit541  O-C:C-C
$pc[415] = "Occ[#1]";      # bit542  O-C:C-[#1]
$pc[416] = "OccN";         # bit543  O-C:C-N
$pc[417] = "OccO";         # bit544  O-C:C-O
$pc[418] = "NccC";         # bit545  N-C:C-C
$pc[419] = "Ncc[#1]";      # bit546  N-C:C-[#1]
$pc[420] = "NccN";         # bit547  N-C:C-N
$pc[421] = "OCcc";         # bit548  O-C-C:C
$pc[422] = "NCcc";         # bit549  N-C-C:C
$pc[423] = "[#17]CCC";     # bit550  Cl-C-C-C
$pc[424] = "[#17]CCO";     # bit551  Cl-C-C-O
$pc[425] = "cc-cc";        # bit552  C:C-C:C
$pc[426] = "O=CC=C";       # bit553  O=C-C=C
$pc[427] = "[#35]CCC";     # bit554  Br-C-C-C
$pc[428] = "N=CC=C";       # bit555  N=C-C=C
$pc[429] = "C=CCC";        # bit556  C=C-C-C
$pc[430] = "ncO[#1]";      # bit557  N:C-O-[#1]
$pc[431] = "O=Ncc";        # bit558  O=N-C:C
$pc[432] = "OCN[#1]";      # bit559  O-C-N-[#1]
$pc[433] = "NCNC";         # bit560  N-C-N-C
$pc[434] = "[#17]CC=O";    # bit561  Cl-C-C=O
$pc[435] = "[#35]CC=O";    # bit562  Br-C-C=O
$pc[436] = "OCOC";         # bit563  O-C-O-C
$pc[437] = "C=CC=C";       # bit564  C=C-C=C
$pc[438] = "ccOC";         # bit565  C:C-O-C
$pc[439] = "OCCN";         # bit566  O-C-C-N
$pc[440] = "OCCO";         # bit567  O-C-C-O
$pc[441] = "N#CCC";        # bit568  N#C-C-C
$pc[442] = "NCCN";         # bit569  N-C-C-N
$pc[443] = "ccCC";         # bit570  C:C-C-C
$pc[444] = "[#1]CO[#1]";   # bit571  [#1]-C-O-[#1]
$pc[445] = "ncnc";         # bit572  N:C:N:C
$pc[446] = "OCC=C";        # bit573  O-C-C=C
$pc[447] = "OCccC";        # bit574  O-C-C:C-C
$pc[448] = "OCccO";        # bit575  O-C-C:C-O
$pc[449] = "N=Ccc[#1]";    # bit576  N=C-C:C-[#1]
$pc[450] = "ccNcc";        # bit577  C:C-N-C:C
$pc[451] = "Ccc-cc";       # bit578  C-C:C-C:C
$pc[452] = "O=CCCC";       # bit579  O=C-C-C-C
$pc[453] = "OCCCN";        # bit580  O=C-C-C-N
$pc[454] = "O=CCCO";       # bit581  O=C-C-C-O
$pc[455] = "CCCCC";        # bit582  C-C-C-C-C
$pc[456] = "[#17]ccOC";    # bit583  Cl-C:C-O-C
$pc[457] = "ccC=CC";       # bit584  C:C-C=C-C
$pc[458] = "CccNC";        # bit585  C-C:C-N-C
$pc[459] = "CSCCC";        # bit586  C-S-C-C-C
$pc[460] = "NccO[#1]";     # bit587  N-C:C-O-[#1]
$pc[461] = "O=CCC=O";      # bit588  O=C-C-C=O
$pc[462] = "CccOC";        # bit589  C-C:C-O-C
$pc[463] = "CccO[#1]";     # bit590  C-C:C-O-[#1]
$pc[464] = "[#17]CCCC";    # bit591  Cl-C-C-C-C 
$pc[465] = "NCCCC";        # bit592  N-C-C-C-C  
$pc[466] = "NCCCN";        # bit593  N-C-C-C-N
$pc[467] = "COCC=C";       # bit594  C-O-C-C=C
$pc[468] = "ccCCC";        # bit595  C:C-C-C-C
$pc[469] = "N=CNCC";       # bit596  N=C-N-C-C
$pc[470] = "O=CCcc";       # bit597  O=C-C-C:C
$pc[471] = "[#17]cccC";    # bit598  Cl-C:C:C-C
$pc[472] = "[#1]CC=C[#1]"; # bit599  [#1]-C-C=C-[#1]
$pc[473] = "NcccC";        # bit600  N-C:C:C-C
$pc[474] = "NcccN";        # bit601  N-C:C:C-N
$pc[475] = "O=CCNC";       # bit602  O=C-C-N-C
$pc[476] = "CcccC";        # bit603  C-C:C:C-C
$pc[477] = "COCcc";        # bit604  C-O-C-C:C
$pc[478] = "O=CCOC";       # bit605  O=C-C-O-C
$pc[479] = "OccCC";        # bit606  O-C:C-C-C
$pc[480] = "NCCcc";        # bit607  N-C-C-C:C
$pc[481] = "CCCcc";        # bit608  C-C-C-C:C
$pc[482] = "[#17]CCNC";    # bit609  Cl-C-C-N-C
$pc[483] = "COCOC";        # bit610  C-O-C-O-C
$pc[484] = "NCCNC";        # bit611  N-C-C-N-C
$pc[485] = "NCOCC";        # bit612  N-C-O-C-C
$pc[486] = "CNCCC";        # bit613  C-N-C-C-C
$pc[487] = "CCOCC";        # bit614  C-C-O-C-C
$pc[488] = "NCCOC";        # bit615  N-C-C-O-C
$pc[489] = "ccnnc";        # bit616  C:C:N:N:C
$pc[490] = "CCCO[#1]";     # bit617  C-C-C-O-[#1]
$pc[491] = "ccCcc";        # bit618  C:C-C-C:C
$pc[492] = "OCC=CC";       # bit619  O-C-C=C-C
$pc[493] = "ccOCC";        # bit620  C:C-O-C-C
$pc[494] = "Ncccn";        # bit621  N-C:C:C:N
$pc[495] = "O=COcc";       # bit622  O=C-O-C:C
$pc[496] = "O=CccC";       # bit623  O=C-C:C-C
$pc[497] = "O=CccN";       # bit624  O=C-C:C-N
$pc[498] = "O=CccO";       # bit625  O=C-C:C-O
$pc[499] = "COccC";        # bit626  C-O-C:C-C
$pc[500] = "O=[#33]ccc";   # bit627  O=[As]-C:C:C
$pc[501] = "CNCcc";        # bit628  C-N-C-C:C
$pc[502] = "ScccN";        # bit629  S-C:C:C-N
$pc[503] = "OccOC";        # bit630  O-C:C-O-C
$pc[504] = "OccO[#1]";     # bit631  O-C:C-O-[#1]
$pc[505] = "CCOcc";        # bit632  C-C-O-C:C
$pc[506] = "NCccC";        # bit633  N-C-C:C-C
$pc[507] = "CCccC";        # bit634  C-C-C:C-C
$pc[508] = "NNCN[#1]";     # bit635  N-N-C-N-[#1]
$pc[509] = "CNCNC";        # bit636  C-N-C-N-C
$pc[510] = "OCCCC";        # bit637  O-C-C-C-C
$pc[511] = "OCCCN";        # bit638  O-C-C-C-N
$pc[512] = "OCCCO";        # bit639  O-C-C-C-O
$pc[513] = "C=CCCC";       # bit640  C=C-C-C-C
$pc[514] = "OCCC=C";       # bit641  O-C-C-C=C
$pc[515] = "OCCC=O";       # bit642  O-C-C-C=O
$pc[516] = "[#1]CCN[#1]";  # bit643  [#1]-C-C-N-[#1]
$pc[517] = "CC=NNC";       # bit644  C-C=N-N-C
$pc[518] = "O=CNCC";       # bit645  O=C-N-C-C
$pc[519] = "O=CNC[#1]";    # bit646  O=C-N-C-[#1]
$pc[520] = "O=CNCN";       # bit647  O=C-N-C-N
$pc[521] = "O=NccN";       # bit648  O=N-C:C-N
$pc[522] = "O=NccO";       # bit649  O=N-C:C-O      
$pc[523] = "O=CNC=O";      # bit650  O=C-N-C=O
$pc[524] = "OcccC";        # bit651  O-C:C:C-C
$pc[525] = "OcccN";        # bit652  O-C:C:C-N
$pc[526] = "OcccO";        # bit653  O-C:C:C-O
$pc[527] = "NCNCC";        # bit654  N-C-N-C-C
$pc[528] = "OCCcc";        # bit655  O-C-C-C:C
$pc[529] = "CCNCC";        # bit656  C-C-N-C-C
$pc[530] = "CNccC";        # bit657  C-N-C:C-C
$pc[531] = "CCSCC";        # bit658  C-C-S-C-C
$pc[532] = "OCCNC";        # bit659  O-C-C-N-C
$pc[533] = "CC=CCC";       # bit660  C-C=C-C-C
$pc[534] = "OCOCC";        # bit661  O-C-O-C-C
$pc[535] = "OCCOC";        # bit662  O-C-C-O-C
$pc[536] = "OCCO[#1]";     # bit663  O-C-C-O-[#1]
$pc[537] = "CC=CC=C";      # bit664  C-C=C-C=C
$pc[538] = "NccCC";        # bit665  N-C:C-C-C
$pc[539] = "C=CCOC";       # bit666  C=C-C-O-C
$pc[540] = "C=CCO[#1]";    # bit667  C=C-C-O-[#1]
$pc[541] = "CccCC";        # bit668  C-C:C-C-C
$pc[542] = "[#17]ccC=O";   # bit669  Cl-C:C-C=O
$pc[543] = "[#35]cccC";    # bit670  Br-C:C:C-C
$pc[544] = "O=CC=CC";      # bit671  O=C-C=C-C
$pc[545] = "O=CC=C[#1]";   # bit672  O=C-C=C-[#1]
$pc[546] = "O=CC=CN";      # bit673  O=C-C=C-N
$pc[547] = "NCNcc";        # bit674  N-C-N-C:C
$pc[548] = "[#35]CCcc";    # bit675  Br-C-C-C:C
$pc[549] = "N#CCCC";       # bit676  N#C-C-C-C
$pc[550] = "CC=Ccc";       # bit677  C-C=C-C:C
$pc[551] = "CCC=CC";       # bit678  C-C-C=C-C
$pc[552] = "CCCCCC";       # bit679  C-C-C-C-C-C
$pc[553] = "OCCCCC";       # bit680  O-C-C-C-C-C
$pc[554] = "OCCCCO";       # bit681  O-C-C-C-C-O
$pc[555] = "OCCCCN";       # bit682  O-C-C-C-C-N
$pc[556] = "NCCCCC";       # bit683  N-C-C-C-C-C
$pc[557] = "O=CCCCC";      # bit684  O=C-C-C-C-C
$pc[558] = "O=CCCCN";      # bit685  O=C-C-C-C-N
$pc[559] = "O=CCCCO";      # bit686  O=C-C-C-C-O
$pc[560] = "O=CCCC=O";     # bit687  O=C-C-C-C=O
$pc[561] = "CCCCCCC";      # bit688  C-C-C-C-C-C-C
$pc[562] = "OCCCCCC";      # bit689  O-C-C-C-C-C-C
$pc[563] = "OCCCCCO";      # bit690  O-C-C-C-C-C-O
$pc[564] = "OCCCCCN";      # bit691  O-C-C-C-C-C-N
$pc[565] = "O=CCCCCC";     # bit692  O=C-C-C-C-C-C
$pc[566] = "O=CCCCCO";     # bit693  O=C-C-C-C-C-O
$pc[567] = "O=CCCCC=O";    # bit694  O=C-C-C-C-C=O
$pc[568] = "O=CCCCCN";     # bit695  O=C-C-C-C-C-N
$pc[569] = "CCCCCCCC";     # bit696  C-C-C-C-C-C-C-C
$pc[570] = "CCCCCC(C)C";   # bit697  C-C-C-C-C-C(C)-C
$pc[571] = "OCCCCCCC";     # bit698  O-C-C-C-C-C-C-C
$pc[572] = "OCCCCC(C)C";   # bit699  O-C-C-C-C-C(C)-C 
$pc[573] = "OCCCCCOC";     # bit700  O-C-C-C-C-C-O-C
$pc[574] = "OCCCCC(O)C";   # bit701  O-C-C-C-C-C(O)-C
$pc[575] = "OCCCCCNC";     # bit702  O-C-C-C-C-C-N-C
$pc[576] = "OCCCCC(N)C";   # bit703  O-C-C-C-C-C(N)-C
$pc[577] = "O=CCCCCCC";    # bit704  O=C-C-C-C-C-C-C
$pc[578] = "O=CCCCC(O)C";  # bit705  O=C-C-C-C-C(O)-C
$pc[579] = "O=CCCCC(=O)C"; # bit706  O=C-C-C-C-C(=O)-C
$pc[580] = "O=CCCCC(N)C";  # bit707  O=C-C-C-C-C(N)-C
$pc[581] = "CC(C)CC";      # bit708  C-C(C)-C-C
$pc[582] = "CC(C)CCC";     # bit709  C-C(C)-C-C-C
$pc[583] = "CCC(C)CC";     # bit710  C-C-C(C)-C-C
$pc[584] = "CC(C)(C)CC";   # bit711  C-C(C)(C)-C-C
$pc[585] = "CC(C)C(C)C";   # bit712  C-C(C)-C(C)-C
$pc[586] = "Cc1ccc(C)cc1";      # bit713  Section 7: Complex SMARTS patterns Cc1ccc(C)cc1
$pc[587] = "Cc1ccc(O)cc1";      # bit714  Cc1ccc(O)cc1         
$pc[588] = "Cc1ccc(S)cc1";      # bit715  Cc1ccc(S)cc1        
$pc[589] = "Cc1ccc(N)cc1";      # bit716  Cc1ccc(N)cc1
$pc[590] = "Cc1ccc(Cl)cc1";     # bit717  Cc1ccc(Cl)cc1
$pc[591] = "Cc1ccc(Br)cc1";     # bit718  Cc1ccc(Br)cc1
$pc[592] = "Oc1ccc(O)cc1";      # bit719  Oc1ccc(O)cc1
$pc[593] = "Oc1ccc(S)cc1";      # bit720  Oc1ccc(S)cc1
$pc[594] = "Oc1ccc(N)cc1";      # bit721  Oc1ccc(N)cc1
$pc[595] = "Oc1ccc(Cl)cc1";     # bit722  Oc1ccc(Cl)cc1
$pc[596] = "Oc1ccc(Br)cc1";     # bit723  Oc1ccc(Br)cc1
$pc[597] = "Sc1ccc(S)cc1";      # bit724  Sc1ccc(S)cc1
$pc[598] = "Sc1ccc(N)cc1";      # bit725  Sc1ccc(N)cc1
$pc[599] = "Sc1ccc(Cl)cc1";     # bit726  Sc1ccc(Cl)cc1
$pc[600] = "Sc1ccc(Br)cc1";     # bit727  Sc1ccc(Br)cc1
$pc[601] = "Nc1ccc(N)cc1";      # bit728  Nc1ccc(N)cc1
$pc[602] = "Nc1ccc(Cl)cc1";     # bit729  Nc1ccc(Cl)cc1
$pc[603] = "Nc1ccc(Br)cc1";     # bit730  Nc1ccc(Br)cc1
$pc[604] = "Clc1ccc(Cl)cc1";    # bit731  Clc1ccc(Cl)cc1
$pc[605] = "Clc1ccc(Br)cc1";    # bit732  Clc1ccc(Br)cc1
$pc[606] = "Brc1ccc(Br)cc1";    # bit733  Brc1ccc(Br)cc1
$pc[607] = "Cc1cc(C)ccc1";      # bit734  Cc1cc(C)ccc1
$pc[608] = "Cc1cc(O)ccc1";      # bit735  Cc1cc(O)ccc1
$pc[609] = "Cc1cc(S)ccc1";      # bit736  Cc1cc(S)ccc1
$pc[610] = "Cc1cc(N)ccc1";      # bit737  Cc1cc(N)ccc1
$pc[611] = "Cc1cc(Cl)ccc1";     # bit738  Cc1cc(Cl)ccc1
$pc[612] = "Cc1cc(Br)ccc1";     # bit739  Cc1cc(Br)ccc1
$pc[613] = "Oc1cc(O)ccc1";      # bit740  Oc1cc(O)ccc1
$pc[614] = "Oc1cc(S)ccc1";      # bit741  Oc1cc(S)ccc1
$pc[615] = "Oc1cc(N)ccc1";      # bit742  Oc1cc(N)ccc1
$pc[616] = "Oc1cc(Cl)ccc1";     # bit743  Oc1cc(Cl)ccc1
$pc[617] = "Oc1cc(Br)ccc1";     # bit744  Oc1cc(Br)ccc1
$pc[618] = "Sc1cc(S)ccc1";      # bit745  Sc1cc(S)ccc1
$pc[619] = "Sc1cc(N)ccc1";      # bit746  Sc1cc(N)ccc1
$pc[620] = "Sc1cc(Cl)ccc1";     # bit747  Sc1cc(Cl)ccc1
$pc[621] = "Sc1cc(Br)ccc1";     # bit748  Sc1cc(Br)ccc1
$pc[622] = "Nc1cc(N)ccc1";      # bit749  Nc1cc(N)ccc1
$pc[623] = "Nc1cc(Cl)ccc1";     # bit750  Nc1cc(Cl)ccc1
$pc[624] = "Nc1cc(Br)ccc1";     # bit751  Nc1cc(Br)ccc1
$pc[625] = "Clc1cc(Cl)ccc1";    # bit752  Clc1cc(Cl)ccc1
$pc[626] = "Clc1cc(Br)ccc1";    # bit753  Clc1cc(Br)ccc1
$pc[627] = "Brc1cc(Br)ccc1";    # bit754  Brc1cc(Br)ccc1
$pc[628] = "Cc1c(C)cccc1";      # bit755  Cc1c(C)cccc1
$pc[629] = "Cc1c(O)cccc1";      # bit756  Cc1c(O)cccc1
$pc[630] = "Cc1c(S)cccc1";      # bit757  Cc1c(S)cccc1
$pc[631] = "Cc1c(N)cccc1";      # bit758  Cc1c(N)cccc1
$pc[632] = "Cc1c(Cl)cccc1";     # bit759  Cc1c(Cl)cccc1
$pc[633] = "Cc1c(Br)cccc1";     # bit760  Cc1c(Br)cccc1
$pc[634] = "Oc1c(O)cccc1";      # bit761  Oc1c(O)cccc1
$pc[635] = "Oc1c(S)cccc1";      # bit762  Oc1c(S)cccc1
$pc[636] = "Oc1c(N)cccc1";      # bit763  Oc1c(N)cccc1
$pc[637] = "Oc1c(Cl)cccc1";     # bit764  Oc1c(Cl)cccc1
$pc[638] = "Oc1c(Br)cccc1";     # bit765  Oc1c(Br)cccc1
$pc[639] = "Sc1c(S)cccc1";      # bit766  Sc1c(S)cccc1
$pc[640] = "Sc1c(N)cccc1";      # bit767  Sc1c(N)cccc1
$pc[641] = "Sc1c(Cl)cccc1";     # bit768  Sc1c(Cl)cccc1
$pc[642] = "Sc1c(Br)cccc1";     # bit769  Sc1c(Br)cccc1
$pc[643] = "Nc1c(N)cccc1";      # bit770  Nc1c(N)cccc1
$pc[644] = "Nc1c(Cl)cccc1";     # bit771  Nc1c(Cl)cccc1
$pc[645] = "Nc1c(Br)cccc1";     # bit772  Nc1c(Br)cccc1       
$pc[646] = "Clc1c(Cl)cccc1";    # bit773  Clc1c(Cl)cccc1
$pc[647] = "Clc1c(Br)cccc1";    # bit774  Clc1c(Br)cccc1
$pc[648] = "Brc1c(Br)cccc1";    # bit775  Brc1c(Br)cccc1
$pc[649] = "CC1CCC(C)CC1";      # bit776  CC1CCC(C)CC1
$pc[650] = "CC1CCC(O)CC1";      # bit777  CC1CCC(O)CC1
$pc[651] = "CC1CCC(S)CC1";      # bit778  CC1CCC(S)CC1
$pc[652] = "CC1CCC(N)CC1";      # bit779  CC1CCC(N)CC1
$pc[653] = "CC1CCC(Cl)CC1";     # bit780  CC1CCC(Cl)CC1
$pc[654] = "CC1CCC(Br)CC1";     # bit781  CC1CCC(Br)CC1
$pc[655] = "OC1CCC(O)CC1";      # bit782  OC1CCC(O)CC1
$pc[656] = "OC1CCC(S)CC1";      # bit783  OC1CCC(S)CC1
$pc[657] = "OC1CCC(N)CC1";      # bit784  OC1CCC(N)CC1
$pc[658] = "OC1CCC(Cl)CC1";     # bit785  OC1CCC(Cl)CC1
$pc[659] = "OC1CCC(Br)CC1";     # bit786  OC1CCC(Br)CC1
$pc[660] = "SC1CCC(S)CC1";      # bit787  SC1CCC(S)CC1
$pc[661] = "SC1CCC(N)CC1";      # bit788  SC1CCC(N)CC1
$pc[662] = "SC1CCC(Cl)CC1";     # bit789  SC1CCC(Cl)CC1
$pc[663] = "SC1CCC(Br)CC1";     # bit790  SC1CCC(Br)CC1
$pc[664] = "NC1CCC(N)CC1";      # bit791  NC1CCC(N)CC1
$pc[665] = "NC1CCC(Cl)CC1";     # bit792  NC1CCC(Cl)CC1 
$pc[666] = "NC1CCC(Br)CC1";     # bit793  NC1CCC(Br)CC1
$pc[667] = "ClC1CCC(Cl)CC1";    # bit794  ClC1CCC(Cl)CC1
$pc[668] = "ClC1CCC(Br)CC1";    # bit795  ClC1CCC(Br)CC1
$pc[669] = "BrC1CCC(Br)CC1";    # bit796  BrC1CCC(Br)CC1
$pc[670] = "CC1CC(C)CCC1";      # bit797  CC1CC(C)CCC1
$pc[671] = "CC1CC(C)CCC1";      # bit798  CC1CC(C)CCC1
$pc[672] = "CC1CC(S)CCC1";      # bit799  CC1CC(S)CCC1
$pc[673] = "CC1CC(N)CCC1";      # bit800  CC1CC(N)CCC1
$pc[674] = "CC1CC(Cl)CCC1";     # bit801  CC1CC(Cl)CCC1
$pc[675] = "CC1CC(Br)CCC1";     # bit802  CC1CC(Br)CCC1
$pc[676] = "OC1CC(O)CCC1";      # bit803  OC1CC(O)CCC1
$pc[677] = "OC1CC(S)CCC1";      # bit804  OC1CC(S)CCC1
$pc[678] = "OC1CC(N)CCC1";      # bit805  OC1CC(N)CCC1
$pc[679] = "OC1CC(Cl)CCC1";     # bit806  OC1CC(Cl)CCC1
$pc[680] = "OC1CC(Br)CCC1";     # bit807  OC1CC(Br)CCC1
$pc[681] = "SC1CC(S)CCC1";      # bit808  SC1CC(S)CCC1
$pc[682] = "SC1CC(N)CCC1";      # bit809  SC1CC(N)CCC1
$pc[683] = "SC1CC(Cl)CCC1";     # bit810  SC1CC(Cl)CCC1
$pc[684] = "SC1CC(Br)CCC1";     # bit811  SC1CC(Br)CCC1
$pc[685] = "NC1CC(N)CCC1";      # bit812  NC1CC(N)CCC1
$pc[686] = "NC1CC(Cl)CCC1";     # bit813  NC1CC(Cl)CCC1
$pc[687] = "NC1CC(Br)CCC1";     # bit814  NC1CC(Br)CCC1
$pc[688] = "ClC1CC(Cl)CCC1";    # bit815  ClC1CC(Cl)CCC1
$pc[689] = "ClC1CC(Br)CCC1";    # bit816  ClC1CC(Br)CCC1
$pc[690] = "BrC1CC(Br)CCC1";    # bit817  BrC1CC(Br)CCC1
$pc[691] = "CC1C(C)CCCC1";      # bit818  CC1C(C)CCCC1
$pc[692] = "CC1C(O)CCCC1";      # bit819  CC1C(O)CCCC1
$pc[693] = "CC1C(S)CCCC1";      # bit820  CC1C(S)CCCC1
$pc[694] = "CC1C(N)CCCC1";      # bit821  CC1C(N)CCCC1
$pc[695] = "CC1C(Cl)CCCC1";     # bit822  CC1C(Cl)CCCC1
$pc[696] = "CC1C(Br)CCCC1";     # bit823  CC1C(Br)CCCC1
$pc[697] = "OC1C(O)CCCC1";      # bit824  OC1C(O)CCCC1
$pc[698] = "OC1C(S)CCCC1";      # bit825  OC1C(S)CCCC1
$pc[699] = "OC1C(N)CCCC1";      # bit826  OC1C(N)CCCC1
$pc[700] = "OC1C(Cl)CCCC1";     # bit827  OC1C(Cl)CCCC1
$pc[701] = "OC1C(Br)CCCC1";     # bit828  OC1C(Br)CCCC1
$pc[702] = "SC1C(S)CCCC1";      # bit829  SC1C(S)CCCC1
$pc[703] = "SC1C(N)CCCC1";      # bit830  SC1C(N)CCCC1
$pc[704] = "SC1C(Cl)CCCC1";     # bit831  SC1C(Cl)CCCC1
$pc[705] = "SC1C(Br)CCCC1";     # bit832  SC1C(Br)CCCC1
$pc[706] = "NC1C(N)CCCC1";      # bit833  NC1C(N)CCCC1
$pc[707] = "NC1C(Cl)CCCC1";     # bit834  NC1C(Cl)CCCC1
$pc[708] = "NC1C(Br)CCCC1";     # bit835  NC1C(Br)CCCC1
$pc[709] = "ClC1C(Cl)CCCC1";    # bit836  ClC1C(Cl)CCCC1
$pc[710] = "ClC1C(Br)CCCC1";    # bit837  ClC1C(Br)CCCC1
$pc[711] = "BrC1C(Br)CCCC1";    # bit838  BrC1C(Br)CCCC1
$pc[712] = "CC1CC(C)CC1";       # bit839  CC1CC(C)CC1
$pc[713] = "CC1CC(O)CC1";       # bit840  CC1CC(O)CC1
$pc[714] = "CC1CC(S)CC1";       # bit841  CC1CC(S)CC1
$pc[715] = "CC1CC(N)CC1";       # bit842  CC1CC(N)CC1
$pc[716] = "CC1CC(Cl)CC1";      # bit843  CC1CC(Cl)CC1
$pc[717] = "CC1CC(Br)CC1";      # bit844  CC1CC(Br)CC1
$pc[718] = "OC1CC(O)CC1";       # bit845  OC1CC(O)CC1
$pc[719] = "OC1CC(S)CC1";       # bit846  OC1CC(S)CC1
$pc[720] = "OC1CC(N)CC1";       # bit847  OC1CC(N)CC1
$pc[721] = "OC1CC(Cl)CC1";      # bit848  OC1CC(Cl)CC1
$pc[722] = "OC1CC(Br)CC1";      # bit849  OC1CC(Br)CC1
$pc[723] = "SC1CC(S)CC1";       # bit850  SC1CC(S)CC1
$pc[724] = "SC1CC(N)CC1";       # bit851  SC1CC(N)CC1
$pc[725] = "SC1CC(Cl)CC1";      # bit852  SC1CC(Cl)CC1
$pc[726] = "SC1CC(Br)CC1";      # bit853  SC1CC(Br)CC1
$pc[727] = "NC1CC(N)CC1";       # bit854  NC1CC(N)CC1
$pc[728] = "NC1CC(Cl)CC1";      # bit855  NC1CC(Cl)CC1
$pc[729] = "NC1CC(Br)CC1";      # bit856  NC1CC(Br)CC1
$pc[730] = "ClC1CC(Cl)CC1";     # bit857  ClC1CC(Cl)CC1
$pc[731] = "ClC1CC(Br)CC1";     # bit858  ClC1CC(Br)CC1
$pc[732] = "BrC1CC(Br)CC1";     # bit859  BrC1CC(Br)CC1
$pc[733] = "CC1C(C)CCC1";       # bit860  CC1C(C)CCC1
$pc[734] = "CC1C(O)CCC1";       # bit861  CC1C(O)CCC1
$pc[735] = "CC1C(S)CCC1";       # bit862  CC1C(S)CCC1
$pc[736] = "CC1C(N)CCC1";       # bit863  CC1C(N)CCC1
$pc[737] = "CC1C(Cl)CCC1";      # bit864  CC1C(Cl)CCC1
$pc[738] = "CC1C(Br)CCC1";      # bit865  CC1C(Br)CCC1
$pc[739] = "OC1C(O)CCC1";       # bit866  OC1C(O)CCC1
$pc[740] = "OC1C(S)CCC1";       # bit867  OC1C(S)CCC1
$pc[741] = "OC1C(N)CCC1";       # bit868  OC1C(N)CCC1
$pc[742] = "OC1C(Cl)CCC1";      # bit869  OC1C(Cl)CCC1
$pc[743] = "OC1C(Br)CCC1";      # bit870  OC1C(Br)CCC1
$pc[744] = "SC1C(S)CCC1";       # bit871  SC1C(S)CCC1
$pc[745] = "SC1C(N)CCC1";       # bit872  SC1C(N)CCC1
$pc[746] = "SC1C(Cl)CCC1";      # bit873  SC1C(Cl)CCC1
$pc[747] = "SC1C(Br)CCC1";      # bit874  SC1C(Br)CCC1
$pc[748] = "NC1C(N)CCC1";       # bit875  NC1C(N)CCC1
$pc[749] = "NC1C(Cl)CC1";       # bit876  NC1C(Cl)CC1
$pc[750] = "NC1C(Br)CCC1";      # bit877  NC1C(Br)CCC1
$pc[751] = "ClC1C(Cl)CCC1";     # bit878  ClC1C(Cl)CCC1
$pc[752] = "ClC1C(Br)CCC1";     # bit879  ClC1C(Br)CCC1
$pc[753] = "BrC1C(Br)CCC1";     # bit880  BrC1C(Br)CCC1

$npbmkeys = 42;  # number of multiple occurrences keys 
#
#  PubChem keys that set additional bits on multiple occurence
#
$mpmbit[1][1] = $pc[80] ; # 3-ring carbon only
$mpmbit[1][2] = 786 ; #  bit122  >= 2 occurrences
$mpmbit[1][3] = 0   ; # >= 3 
$mpmbit[1][4] = 0   ; # >= 4 
$mpmbit[1][5] = 0   ; # >= 5
$mpmbit[2][1] = $pc[81] ; # 3-ring carbon only
$mpmbit[2][2] = 787 ; #  bit123  >= 2 occurrences
$mpmbit[2][3] = 0   ; # >= 3 
$mpmbit[3][1] = $pc[82] ; # 3-ring with nitrogen
$mpmbit[3][2] = 788 ; #  bit124  >= 2 occurrences
$mpmbit[3][3] = 0   ; # >= 3 
$mpmbit[4][1] = $pc[83] ; # 3-ring hetero
$mpmbit[4][2] = 789 ; #  bit125  >= 2 occurrences
$mpmbit[4][3] = 0   ; # >= 3 
$mpmbit[5][1] = $pc[84] ; # 3-ring unsaturated carbon only
$mpmbit[5][2] = 790 ; #  bit126  >= 2 occurrences
$mpmbit[5][3] = 0   ; # >= 3 
$mpmbit[6][1] = $pc[85] ; # 3-ring unsaturated with nitrogen
$mpmbit[6][2] = 791 ; #  bit127  >= 2 occurrences
$mpmbit[6][3] = 0   ; # >= 3 
$mpmbit[7][1] = $pc[86] ; # 3-ring unsaturated hetero
$mpmbit[7][2] = 792 ; #  bit128  >= 2 occurrences
$mpmbit[7][3] = 0   ; # >= 3 
$mpmbit[8][1] = $pc[87] ; # 4-ring 
$mpmbit[8][2] = 793 ; #  bit136  >= 2 occurrences
$mpmbit[8][3] = 0   ; # >= 3 
$mpmbit[9][1] = $pc[88] ; # 4-ring carbon only 
$mpmbit[9][2] = 794 ; #  bit137  >= 2 occurrences
$mpmbit[9][3] = 0   ; # >= 3 
$mpmbit[10][1] = $pc[89] ; # 4-ring with N
$mpmbit[10][2] = 795 ; #  bit138  >= 2 occurrences
$mpmbit[10][3] = 0   ; # >= 3 
$mpmbit[11][1] = $pc[90] ; # 4-ring hetero
$mpmbit[11][2] = 796 ; #  bit139  >= 2 occurrences
$mpmbit[11][3] = 0   ; # >= 3 
$mpmbit[12][1] = $pc[91] ; # 4-ring unsaturated C only
$mpmbit[12][2] = 797 ; #  bit140  >= 2 occurrences
$mpmbit[12][3] = 0   ; # >= 3 
$mpmbit[13][1] = $pc[92] ; # 4-ring unsaturated N
$mpmbit[13][2] = 798 ; #  bit141  >= 2 occurrences
$mpmbit[13][3] = 0   ; # >= 3 
$mpmbit[14][1] = $pc[93] ; # 4-ring unsaturated hetero
$mpmbit[14][2] = 799 ; #  bit142  >= 2 occurrences
$mpmbit[14][3] = 0   ; # >= 3 
$mpmbit[15][1] = $pc[94] ; # 5-ring 
$mpmbit[15][2] = 799 ; #  bit150  >= 2 occurrences
$mpmbit[15][3] = 800 ; #  bit157  >= 3 
$mpmbit[15][3] = 801 ; #  bit164  >= 4 
$mpmbit[15][3] = 802 ; #  bit171  >= 5 
$mpmbit[16][1] = $pc[95] ; # 5-ring C only 
$mpmbit[16][2] = 803 ; #  bit151  >= 2 occurrences
$mpmbit[16][3] = 804 ; #  bit158  >= 3 
$mpmbit[16][3] = 805 ; #  bit165  >= 4 
$mpmbit[16][3] = 806 ; #  bit172  >= 5 
$mpmbit[17][1] = $pc[96] ; # 5-ring with N
$mpmbit[17][2] = 807 ; #  bit152  >= 2 occurrences
$mpmbit[17][3] = 808 ; #  bit159  >= 3 
$mpmbit[17][3] = 809 ; #  bit166  >= 4 
$mpmbit[17][3] = 810 ; #  bit173  >= 5 
$mpmbit[18][1] = $pc[97] ; # 5-ring hetero
$mpmbit[18][2] = 811 ; #  bit153  >= 2 occurrences
$mpmbit[18][3] = 812 ; #  bit160  >= 3 
$mpmbit[18][3] = 813 ; #  bit167  >= 4 
$mpmbit[18][3] = 814 ; #  bit174  >= 5 
$mpmbit[19][1] = $pc[98] ; # 5-ring unsat. C only
$mpmbit[19][2] = 815 ; #  bit154  >= 2 occurrences
$mpmbit[19][3] = 816 ; #  bit161  >= 3 
$mpmbit[19][3] = 817 ; #  bit168  >= 4 
$mpmbit[19][3] = 818 ; #  bit175  >= 5 
$mpmbit[20][1] = $pc[99] ; # 5-ring unsat. N
$mpmbit[20][2] = 819 ; #  bit155  >= 2 occurrences
$mpmbit[20][3] = 820 ; #  bit162  >= 3 
$mpmbit[20][3] = 821 ; #  bit169  >= 4 
$mpmbit[20][3] = 822 ; #  bit176  >= 5 
$mpmbit[21][1] = $pc[100] ; # 5-ring unsat. hetero
$mpmbit[21][2] = 823 ; #  bit156  >= 2 occurrences
$mpmbit[21][3] = 824 ; #  bit163  >= 3 
$mpmbit[21][3] = 825 ; #  bit170  >= 4 
$mpmbit[21][3] = 826 ; #  bit177  >= 5 
$mpmbit[22][1] = $pc[101] ; # 6-ring 
$mpmbit[22][2] = 827 ; #  bit185  >= 2 occurrences
$mpmbit[22][3] = 828 ; #  bit192  >= 3 
$mpmbit[22][3] = 829 ; #  bit199  >= 4 
$mpmbit[22][3] = 830 ; #  bit206  >= 5 
$mpmbit[23][1] = $pc[102] ; # 6-ring C only 
$mpmbit[23][2] = 831 ; #  bit186  >= 2 occurrences
$mpmbit[23][3] = 832 ; #  bit193  >= 3 
$mpmbit[23][3] = 833 ; #  bit200  >= 4 
$mpmbit[23][3] = 834 ; #  bit207  >= 5 
$mpmbit[24][1] = $pc[103] ; # 6-ring with N
$mpmbit[24][2] = 835 ; #  bit187  >= 2 occurrences
$mpmbit[24][3] = 836 ; #  bit194  >= 3 
$mpmbit[24][3] = 837 ; #  bit201  >= 4 
$mpmbit[24][3] = 838 ; #  bit208  >= 5 
$mpmbit[25][1] = $pc[104] ; # 6-ring hetero
$mpmbit[25][2] = 839 ; #  bit188  >= 2 occurrences
$mpmbit[25][3] = 840 ; #  bit195  >= 3 
$mpmbit[25][3] = 841 ; #  bit202  >= 4 
$mpmbit[25][3] = 842 ; #  bit209  >= 5 
$mpmbit[26][1] = $pc[105] ; # 6-ring unsat. C only
$mpmbit[26][2] = 843 ; #  bit189  >= 2 occurrences
$mpmbit[26][3] = 844 ; #  bit196  >= 3 
$mpmbit[26][3] = 845 ; #  bit203  >= 4 
$mpmbit[26][3] = 846 ; #  bit210  >= 5 
$mpmbit[27][1] = $pc[106] ; # 6-ring unsat. N
$mpmbit[27][2] = 847 ; #  bit190  >= 2 occurrences
$mpmbit[27][3] = 848 ; #  bit197  >= 3 
$mpmbit[27][3] = 849 ; #  bit204  >= 4 
$mpmbit[27][3] = 850 ; #  bit211  >= 5 
$mpmbit[28][1] = $pc[107] ; # 6-ring unsat. hetero
$mpmbit[28][2] = 851 ; #  bit191  >= 2 occurrences
$mpmbit[28][3] = 852 ; #  bit198  >= 3 
$mpmbit[28][3] = 853 ; #  bit205  >= 4 
$mpmbit[28][3] = 854 ; #  bit212  >= 5 
$mpmbit[29][1] = $pc[108] ; # 7-ring
$mpmbit[29][2] = 855 ; #  bit220  >= 2 occurrences
$mpmbit[29][3] = 0 ; 
$mpmbit[30][1] = $pc[109] ; # 7-ring C only
$mpmbit[30][2] = 856 ; #  bit221  >= 2 occurrences
$mpmbit[30][3] = 0 ; 
$mpmbit[31][1] = $pc[110] ; # 7-ring with N
$mpmbit[31][2] = 857 ; #  bit222  >= 2 occurrences
$mpmbit[31][3] = 0 ; 
$mpmbit[31][1] = $pc[111] ; # 7-ring hetero
$mpmbit[31][2] = 858 ; #  bit223  >= 2 occurrences
$mpmbit[31][3] = 0 ; 
$mpmbit[32][1] = $pc[112] ; # 7-ring unsat. C only
$mpmbit[32][2] = 859 ; #  bit224  >= 2 occurrences
$mpmbit[32][3] = 0 ; 
$mpmbit[33][1] = $pc[113] ; # 7-ring unsat. N
$mpmbit[33][2] = 860 ; #  bit225  >= 2 occurrences
$mpmbit[33][3] = 0 ; 
$mpmbit[34][1] = $pc[114] ; # 7-ring unsat. N
$mpmbit[34][2] = 861 ; #  bit226  >= 2 occurrences
$mpmbit[34][3] = 0 ; 
$mpmbit[35][1] = $pc[115] ; # 8-ring 
$mpmbit[35][2] = 862 ; #  bit234  >= 2 occurrences
$mpmbit[35][3] = 0 ; 
$mpmbit[36][1] = $pc[116] ; # 8-ring C only 
$mpmbit[36][2] = 863 ; #  bit235  >= 2 occurrences
$mpmbit[36][3] = 0 ; 
$mpmbit[37][1] = $pc[117] ; # 8-ring with N
$mpmbit[37][2] = 864 ; #  bit236  >= 2 occurrences
$mpmbit[37][3] = 0 ; 
$mpmbit[38][1] = $pc[118] ; # 8-ring hetero
$mpmbit[38][2] = 865 ; #  bit237  >= 2 occurrences
$mpmbit[38][3] = 0 ; 
$mpmbit[39][1] = $pc[119] ; # 8-ring unsat. C only
$mpmbit[39][2] = 866 ; #  bit238  >= 2 occurrences
$mpmbit[39][3] = 0 ; 
$mpmbit[40][1] = $pc[120] ; # 8-ring unsat. N
$mpmbit[40][2] = 867 ; #  bit239  >= 2 occurrences
$mpmbit[40][3] = 0 ; 
$mpmbit[41][1] = $pc[121] ; # 8-ring unsat. hetero
$mpmbit[41][2] = 868 ; #  bit240  >= 2 occurrences
$mpmbit[41][3] = 0 ; 
$mpmbit[42][1] = "c1ccccc1" ; # aromatic ring C only
$mpmbit[42][2] = 869 ; #  bit257  >= 2 occurrences
$mpmbit[42][3] = 870 ; #  bit259  >= 3
$mpmbit[42][4] = 871 ; #  bit261  >= 4
$mpmbit[42][5] = 0 ; 
$mpmbit[42][6] = 872 ; #  bit255  >= 1 occurrences 


for ($i=1 ; $i <= $ncomp; $i++) { # loop over all compounds

# print "$SMILES[$i] "; 

  open(SMIL,"> kkkk.smi") or die "error opening SMILES file : $!\n"; 
    print SMIL $SMILES[$i];
  close(SMIL);

# determine MACCS keys fingerprint for each SMILES
 
  $mcfp[$i]= 0 ;
  for ($j=1 ; $j <= $maxmekeys; $j++) {   # "easy" keys
    if ($mc[$j] ne "") {
      $SMART = " -c \"".$mc[$j]."\"" ;
#     print "$SMART\n" ;
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
# look for matching line from obgrep output
      $matchl = system("/bin/grep 1 jjj > jjjj") ;
# UNIX grep returns 0 if selected lines are found
# otherwise 1
      if($matchl < 1 ) {
        $flag = 1 ;
#       print "$j "; #   print number of key that is set on
      }
      else {
        $flag = 0 ;
      }

      $mcfp[$i] = $mcfp[$i] + $flag * 2**$j;
    }
  }

  for ($j=1 ; $j <= $nmcmult; $j++) {   #  keys with $
    $npos = $maxmekeys + $j ;
    for ($k=1 ; $k <= $ncmult[$j]; $k++) {   
      $flag = 0 ;
      $SMART = " -c \"".$md[$j][$k]."\"" ;
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
      $matchl = system("/bin/grep 1 jjj > jjjj") ;
      if($matchl < 1 ) {
        $flag = 1 ;
      }
    }
    if ($flag > 0) {
#     print "$npos ";  # print number of key being set on
      $mcfp[$i] = $mcfp[$i] +  2**$npos;
    }
  }

#  MACCS keys that set additional bits on multiple occurence
 
  for ($j=1 ; $j <= $nmcmbits; $j++) {   # multiple occurrences
    $flagm = 0 ; # number of matches
    $flag  = 0 ; # matches at all 

    $kmax = 4 ; # default: check up to 4 occurrences
    if ($j == 5) {$kmax = 10} #  6-membered rings
    if ($j == 6) {$kmax = 16} #  CH3 groups
    if ($j == 7) {$kmax = 10} #  number of aromatic 6 rings
    if ($j == 8) {$kmax = 16} #  number of oxygen atoms

    for ($k=1 ; $k <= $kmax ; $k++) {   
      $SMART = " -t$k \"".$mcmbit[$j][1]."\"" ;
      unlink("jjj");
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
      #     -s means: exists and has nonzero size
      if (-s "$jjj") {$flagm = $flagm + $k}  # number of exact matches
 #    print "$SMART $k $flagm \n" ;
    }

    $SMART = " -c \"".$mcmbit[$j][1]."\"" ;
    $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
    system($cline) ;
    $matchl = system("/bin/grep 1 jjj > jjjj") ; 
    if($matchl < 1 ) {
      $flag = 1 ;   # at least one match
    }
    $npos = 0 ;
    if ($flag > 0) {  # at least 1 match
      if ($flagm > 1) {  # multiple matches
        if ($mcmbit[$j][3] == 0 && $mcmbit[$j][4] == 0) {
          $npos = $mcmbit[$j][2] ;
#         print "$npos "; # print number of key being set on
          $mcfp[$i] = $mcfp[$i] +  2**$npos;
        }
	else {
          $npos = $mcmbit[$j][2] ; # >1 matches
#         print "$npos "; # print number of key being set on
          $mcfp[$i] = $mcfp[$i] +  2**$npos;
          if ($flagm > 2) { # >2 matches
            $npos = $mcmbit[$j][3] ;
#           print "$npos "; # print number of key being set on
            $mcfp[$i] = $mcfp[$i] +  2**$npos;
          }
          if ($flagm > 3) { #  more than 3 matches
            $npos = $mcmbit[$j][4] ;
#           print "$npos "; # print number of key being set on
            $mcfp[$i] = $mcfp[$i] +  2**$npos;
          }
        }
      }
    }
  }

# last MACCS key (fragments present, SMILES contains ".")
#
  $npos = 164 ;
  $flag = 0 ;
  if ($SMILES[$i] =~ /\./){
    $flag = 1 ;
#   print "$npos "; # print number of key being set on
  } 
  $mcfp[$i] = $mcfp[$i] + $flag * 2**$npos;


# determine PubChem fingerprint for each SMILES
 
  $pbfp[$i]= 0 ;
  for ($j=1 ; $j <= $maxpekeys; $j++) {   # "easy" keys
#   print "pc $j \n" ;
    if ($pc[$j] ne "") {
      $SMART = " -c \"".$pc[$j]."\"" ;
#     print "$SMART\n" ;
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
# look for matching line from obgrep output
      $matchl = system("/bin/grep 1 jjj > jjjj") ;
# UNIX grep returns 0 if selected lines are found
# otherwise 1
      if($matchl < 1 ) {
        $flag = 1 ; # >= 1 matches
      }
      else {
        $flag = 0 ;
      }

      $pbfp[$i] = $pbfp[$i] + $flag * 2**$j;
    }
  }

# count number of matches in a string
#my $string = 'John Doe John Done';
#my $matches = () = $string =~ /John/g;
#print $matches; #prints 2
#
  my $cam = () = $SMILES[$i] =~ /Ca/g ; # calcium
  my $cdm = () = $SMILES[$i] =~ /Cd/g ; # cadmium
  my $cem = () = $SMILES[$i] =~ /Ce/g ; # cer
  my $clm = () = $SMILES[$i] =~ /Cl/g ; # chlorine
  my $crm = () = $SMILES[$i] =~ /Cr/g ; # chromium
  my $csm = () = $SMILES[$i] =~ /Cs/g ; # caesium
  my $cum = () = $SMILES[$i] =~ /Cu/g ; # copper
  my $cc1m = () = $SMILES[$i] =~ /C/g ;  # carbon upper case
  my $cc2m = () = $SMILES[$i] =~ /c/g ;  # carbon lower case
# print "$SMILES[$i]\n" ;
# print "Ca $cam  Cl $clm  Cr $crm  Cu $cum  C $cc1m  c $cc2m\n" ;
  $ccm = $cc2m + $cc1m - $cam - $cdm - $cem - $clm - $crm - $csm - $cum ; 

  my $lim = () = $SMILES[$i] =~ /Li/g ; # lithium

  my $bam = () = $SMILES[$i] =~ /Ba/g ; # barium
  my $bem = () = $SMILES[$i] =~ /Be/g ; # beryllium
  my $bim = () = $SMILES[$i] =~ /Bi/g ; # bismuth
  my $brm = () = $SMILES[$i] =~ /Br/g ; # bromine
  my $bb1m = () = $SMILES[$i] =~ /B/g ; # bor upper case
  my $bb2m = () = $SMILES[$i] =~ /b/g ; # bor lower case
  $bbm = $bb2m + $bb1m - $bam - $bem - $bim - $brm ; # bor

  my $nam = () = $SMILES[$i] =~ /Na/g ; # sodium
  my $nbm = () = $SMILES[$i] =~ /Nb/g ; # niob
  my $ndm = () = $SMILES[$i] =~ /Nd/g ; # neodym
  my $nem = () = $SMILES[$i] =~ /Ne/g ; # neon
  my $nn1m = () = $SMILES[$i] =~ /N/g ; # nitrogen upper case
  my $nn2m = () = $SMILES[$i] =~ /n/g ; # nitrogen lower case
  $nnm = $nn2m + $nn1m - $nam - $nbm - $ndm - $nem ; # nitrogen

  my $osm = () = $SMILES[$i] =~ /Os/g ; # osmium
  my $oo1m = () = $SMILES[$i] =~ /O/g ; # oxygen upper case
  my $oo2m = () = $SMILES[$i] =~ /o/g ; # oxygen lower case
  $oom = $oo2m + $oo1m - $osm ; # oxygen

  my $frm = () = $SMILES[$i] =~ /Fr/g ; # francium
  my $ffm = () = $SMILES[$i] =~ /F/g ;  # fluorine
  $ffm = $ffm - $frm ; # fluorine

  my $sbm = () = $SMILES[$i] =~ /Sb/g ; # antimon
  my $scm = () = $SMILES[$i] =~ /Sc/g ; # scandium
  my $sem = () = $SMILES[$i] =~ /Se/g ; # selen
  my $sim = () = $SMILES[$i] =~ /Si/g ; # silicium
  my $smm = () = $SMILES[$i] =~ /Sm/g ; # samarium
  my $snm = () = $SMILES[$i] =~ /Sn/g ; # tin
  my $srm = () = $SMILES[$i] =~ /Sr/g ; # strontium
  my $ss1m = () = $SMILES[$i] =~ /S/g ; # sulfur upper case
  my $ss2m = () = $SMILES[$i] =~ /s/g ; # sulfur lower case
  $ssm = $ss2m + $ss1m - $sbm - $scm - $sem - $sim - $smm - $snm - $srm ;

  my $pam = () = $SMILES[$i] =~ /Pa/g ; # [#91]
  my $pdm = () = $SMILES[$i] =~ /Pd/g ; # palladium
  my $pbm = () = $SMILES[$i] =~ /Pb/g ; # lead
  my $pom = () = $SMILES[$i] =~ /Po/g ; # pollonium
  my $ptm = () = $SMILES[$i] =~ /Pt/g ; # platinum
  my $pp1m = () = $SMILES[$i] =~ /P/g ; # phosphorus upper case
  my $pp2m = () = $SMILES[$i] =~ /p/g ; # phosphorus lower case
  $ppm = $pp2m + $pp1m - $pam - $pdm - $pbm - $pom - $ptm ; # phosphorus

  my $krm = () = $SMILES[$i] =~ /Kr/g ; # krypton
  my $kkm = () = $SMILES[$i] =~ /K/g ;  # potassium
  $kkm = $kkm - $krm ; # potassium

  my $inm = () = $SMILES[$i] =~ /In/g ; # indium
  my $irm = () = $SMILES[$i] =~ /Ir/g ; # iridium
  my $iim = () = $SMILES[$i] =~ /I/g ;  # iodine
  $iim = $iim - $inm - $irm ; # iodine

# print "$SMILES[$i]\n" ;
# print "B $bbm  C $ccm  N $nnm  O $oom  F $ffm  P $ppm  S $ssm\n" ;


#  PubChem keys set on multiple occurence
#
  $npos = $maxpekeys + 1;

  if($lim >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit5  >= 2 Li
  $npos++;
  if($bbm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit7  >= 2 B
  $npos++;
  if($bbm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit8  >= 2 B
  $npos++;
  if($ccm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit9  >= 2 C
  $npos++;
  if($ccm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit10 >= 4 C
  $npos++;
  if($ccm >= 8){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit11 >= 8 C
  $npos++;
  if($ccm >= 16){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit12 >= 16 C
  $npos++;
  if($ccm >= 32){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit13 >= 32 C
  $npos++;
  if($nnm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit15 >= 2 N
  $npos++;
  if($nnm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit16 >= 4 N
  $npos++;
  if($nnm >= 8){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit17 >= 8 N
  $npos++;
  if($oom >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit19 >= 2 O
  $npos++;
  if($oom >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit20 >= 4 O
  $npos++;
  if($oom >= 8){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit21 >= 8 O
  $npos++;
  if($oom >= 16){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit22 >= 16 O
  $npos++;
  if($ffm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit24 >= 2 F
  $npos++;
  if($ffm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit25 >= 4 F
  $npos++;
  if($nam >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit27 >= 2 Na
  $npos++;
  if($sim >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit29 >= 2 Si
  $npos++;
  if($ppm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit31 >= 2 P
  $npos++;
  if($ppm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit32 >= 4 P
  $npos++;
  if($ssm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit34 >= 2 S
  $npos++;
  if($ssm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit35 >= 4 S
  $npos++;
  if($ssm >= 8){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit36 >= 8 S
  $npos++;
  if($clm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit38 >= 2 Cl
  $npos++;
  if($clm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit39 >= 4 Cl
  $npos++;
  if($clm >= 8){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit40 >= 8 Cl
  $npos++;
  if($kkm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit42 >= 2 K
  $npos++;
  if($brm >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit44 >= 2 Br
  $npos++;
  if($brm >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit45 >= 4 Br
  $npos++;
  if($iim >= 2){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit47 >= 2 I
  $npos++;
  if($iim >= 4){$pbfp[$i] = $pbfp[$i] + $flag * 2**$npos} # bit48 >= 4 I

# print "npos after mult pubchem keys $npos \n" ; # => 785

 
  for ($j=1 ; $j <= $npbmkeys; $j++) {   # multiple occurrences
    $flagm = 0 ; # number of matches
    $flag  = 0 ; # matches at all 

#   $ll = $npos + $j ; # index 

    $SMART = " -c \"".$mpmbit[$j][1]."\"" ;
    $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
    system($cline) ;
    $matchl = system("/bin/grep 1 jjj > jjjj") ; 
    if($matchl < 1 ) {
      $flag = 1 ;   # at least one match
    }
    $npos = 0 ;
    
    $kmax = 8 ; # maximum occurrences of pattern tested
                # (ring systems of various sizes)

    for ($k=1 ; $k <= $kmax ; $k++) {   
      $SMART = " -t$k \"".$mpmbit[$j][1]."\"" ;
      unlink("jjj");
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
      #     -s means: exists and has nonzero size
      if (-s "$jjj") {$flagm = $flagm + $k}  # number of exact matches
#     print "$SMART $k $flagm \n" ;
    }

    if ($flag > 0) {  # at least 1 match
        if (defined($mpmbit[$j][6])) {
          $npos = $mpmbit[$j][6] ; # >= 1 matches (previously not set)
#         print "$npos "; # print number of key being set on
          $pbfp[$i] = $pbfp[$i] +  2**$npos;
        }
      if ($flagm > 1) {  # multiple matches
          $npos = $mpmbit[$j][2] ; # >= 2 matches
#         print "$npos "; # print number of key being set on
          $pbfp[$i] = $pbfp[$i] +  2**$npos;
      }
      if ($flagm > 2) {  # >= 3 matches
          $npos = $mpmbit[$j][3] ;
	  if ($npos > 0) {
#           print "$npos "; # print number of key being set on
            $pbfp[$i] = $pbfp[$i] +  2**$npos;
          }
	  else {next}
      }
      if ($flagm > 3) {  # >= 4 matches
        if (defined($mpmbit[$j][4])) {
          $npos = $mpmbit[$j][4] ;
        }
        else {next}
        if ($npos > 0) {
#         print "$npos "; # print number of key being set on
          $pbfp[$i] = $pbfp[$i] +  2**$npos;
        }
        else {next}
      }
      if ($flagm > 4) {  # >= 5 matches
        if (defined($mpmbit[$j][5])) {
          $npos = $mpmbit[$j][5] ;
        }
        else {next}
        if ($npos > 0) {
#         print "$npos "; # print number of key being set on
          $pbfp[$i] = $pbfp[$i] +  2**$npos;
        }
        else {next}
      }
    }
  }

#  number of heteroaromatic rings 
  $flag  = 0 ;
  $flagm = 0 ; # number of matches
  $flag5 = 0 ; # matches at all 
  $flag6 = 0 ; # matches at all 
  $flag7 = 0 ; # matches at all 

  $SMART5 = " -c \""."[!#6]1aaaa1"."\"" ;
  $cline = $obgrep.$SMART5." kkkk.smi > jjj" ;
  system($cline) ;
  $matchl = system("/bin/grep 1 jjj > jjjj") ; 
  if($matchl < 1 ) {
    $flag5 = 1 ;   # >= 1 heteroaromatic 5-ring
  }
  $SMART6 = " -c \""."[!#6]1aaaaa1"."\"" ;
  $cline = $obgrep.$SMART6." kkkk.smi > jjj" ;
  system($cline) ;
  $matchl = system("/bin/grep 1 jjj > jjjj") ; 
  if($matchl < 1 ) {
    $flag6 = 1 ;   # >= 1 heteroaromatic 6-ring
  }
  $SMART7 = " -c \""."[!#6]1aaaaaa1"."\"" ;
  $cline = $obgrep.$SMART7." kkkk.smi > jjj" ;
  system($cline) ;
  $matchl = system("/bin/grep 1 jjj > jjjj") ; 
  if($matchl < 1 ) {
    $flag7 = 1 ;   # >= 1 heteroaromatic 7-ring
  }

  $flag = $flag5 + $flag6 + $flag7 ;

  if ($flag > 0) { # bit256  >= 1 heteroaromatic ring
    $pbfp[$i] = $pbfp[$i] +  2**877;
  }

  $kmax = 8 ; # maximum occurrences of patterns tested

  for ($k=1 ; $k <= $kmax ; $k++) {   
    if ($flag5 >0) {
      $SMART = " -t$k \""."[!#6]1aaaa1"."\"" ;
      unlink("jjj");
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
      if (-s "$jjj") {$flagm = $flagm + $k}  # number of exact matches
    }
    if ($flag6 >0) {
      $SMART = " -t$k \""."[!#6]1aaaaa1"."\"" ;
      unlink("jjj");
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
      if (-s "$jjj") {$flagm = $flagm + $k}  # number of exact matches
    }
    if ($flag7 >0) {
      $SMART = " -t$k \""."[!#6]1aaaaaa1"."\"" ;
      unlink("jjj");
      $cline = $obgrep.$SMART." kkkk.smi > jjj" ;
      system($cline) ;
      if (-s "$jjj") {$flagm = $flagm + $k}  # number of exact matches
    }
  }

  if ($flagm > 1) { # bit258 => 2 heteroaromatic rings
    $pbfp[$i] = $pbfp[$i] +  2**878;
  }
  if ($flagm > 2) { # bit260 => 3 heteroaromatic rings
    $pbfp[$i] = $pbfp[$i] +  2**879;
  }
  if ($flagm > 3) { # bit262 => 4 heteroaromatic rings
    $pbfp[$i] = $pbfp[$i] +  2**880;
  }

# end of PubChem fingerprint

# $hex = sprintf("0x%X", $mcfp[$i]);
# print "$hex " ;
# my $nhex = sprintf("0x%X", $hex); # number is to long for this representation
# print "$nhex \n" ;

  if ($i == 1) { # create database file
    open(DATA,"> $DBFILE") or die "error opening $DBFILE file : $!\n"; 
    print DATA "compound name MACCS PubChem\n";
  }
  print DATA "$SMILES[$i] $smname[$i] " ;
  print DATA $mcfp[$i]->as_hex()," ";  # MACSS fingerprint
  print DATA $pbfp[$i]->as_hex();      # PubChem fingerprint
  print DATA "\n" ; 

}

close(DATA) or die "error closing file $DBFILE : $!\n";

exit(0);

