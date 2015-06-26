#!/usr/bin/perl
use Cwd;
use Switch;
use List::Util qw(min max);
use POSIX;
use strict;
use warnings;

##################################
# read in number of atom $NA and #
# second line information $info2 #
##################################
open XYZ, "<", "$ARGV[0]" or die $!;
my $stem = $ARGV[0];
$stem =~ s/\..*//;
my $NA = <XYZ>;
$NA =~ s/\n//;
my $info2 = <XYZ>;
my $Nve = 0;
my $multi = 1;
my $charge = 0;
my @ATOMS;
my $itr;

##################################
# read atom type and coordinates #
##################################
for(my $itr=0;$itr<$NA;$itr++){
  # read and string into componets
  my @input = split(' ',<XYZ>);
  my @ATOM;
  $ATOM[0] = atom_a2n($input[0]);
  for(my $i=1;$i<4;$i++){$ATOM[$i] = $input[$i];}
  # assign multiD data array
  push (@{$ATOMS[$itr]}, @ATOM);
  $Nve += atom_ve($input[0]);
}
$multi += $Nve%2;
if($multi != 1){
  if($ARGV[0] =~ m/^catDIR/){
    $charge = 1;
    $Nve -= 1;
  }else{
    $charge = -1;
    $Nve += 1;
  }
  $multi = 1;
}
my $strNveu = sprintf("%-14d",$Nve/2);
my $strNved = sprintf("%-14d",$Nve/2);
my $strNeq  = sprintf("%-14d",$ARGV[1]/10);
my $strNst  = sprintf("%-14d",$ARGV[1]);
my $strNwt  = sprintf("%-14d",$ARGV[1]/10);
my $strBlk  = sprintf("%-14d",10);

print <<EOF;
#-------------------#
# CASINO input file #
#-------------------#

# converted from $ARGV[0]

# SYSTEM
neu               : $strNveu #*! Number of up electrons (Integer)
ned               : $strNved #*! Number of down electrons (Integer)
periodic          : F              #*! Periodic boundary conditions (Boolean)
atom_basis_type   : blip           #*! Basis set type (Text)
psi_s             : slater         #*! Type of [anti]symmetrizing wfn (Text)
complex_wf        : F              #*! Wave function real or complex (Boolean)

# RUN
runtype           : vmc_opt        #*! Type of calculation (Text)
newrun            : T              #*! New run or continue old (Boolean)
testrun           : F              #*! Test run flag (Boolean)

# VMC
vmc_equil_nstep   : $strNeq #*! Number of equilibration steps (Integer)
vmc_nstep         : $strNst #*! Number of steps (Integer)
vmc_nblock        : $strBlk #*! Number of checkpoints (Integer)
vmc_nconfig_write : $strNwt #*! Number of configs to write (Integer)

# DMC
dmc_equil_nstep   : 2000           #*! Number of steps (Integer)
dmc_equil_nblock  : 1              #*! Number of checkpoints (Integer)
dmc_stats_nstep   : 10000          #*! Number of steps (Integer)
dmc_target_weight : 1000.d0        #*! Total target weight in DMC (Real)
dtdmc             : 0.002          #*! DMC time step (Real)
use_tmove         : F              #*! Casula nl pp for DMC (Boolean)

# OPTIMIZATION
opt_method        : varmin         #*! Opt method (varmin/madmin/emin/...)
opt_cycles        : 4              #*! Number of optimization cycles (Integer)
opt_jastrow       : T              #*! Optimize Jastrow factor (Boolean)
opt_det_coeff     : F              #*! Optimize determinant coeffs (Boolean)
opt_backflow      : F              #*! Optimize backflow parameters (Boolean)
opt_orbitals      : F              #*! Optimize orbital parameters (Boolean)

# RMC

# GENERAL PARAMETERS
use_jastrow       : T              #*! Use a Jastrow function (Boolean)
backflow          : F              #*! Use backflow corrections (Boolean)
use_orbmods       : F              #*! Use orbital modifications (Boolean)
expot             : F              #*! Use external potential (Boolean)
timing_info       : F              #*! Activate subroutine timers (Boolean)
esupercell        : F              #*! Energy/supercell in output (Boolean)
relativistic      : F              #*! Relativistic effects (Boolean)
neighprint        : 0              #*! Neighbour analysis (Integer)
mpc_cutoff        : 30.d0 hartree  #*! G vector cutoff for MPC (Physical)
forces            : F              #*! Evaluate forces on atoms (Boolean)
checkpoint        : 1              #*! Checkpoint level (Integer)

# EXPECTATION VALUES
density           : F              #*! Accumulate density (Boolean)
spin_density      : F              #*! Accumulate spin densities (Boolean)
pair_corr         : F              #*! Accumulate rec. space PCF (Boolean)
pair_corr_sph     : F              #*! Accumulate sph. real space PCF (Boolean)
loc_tensor        : F              #*! Accumulate localization tensor (Boolean)
structure_factor  : F              #*! Accumulate structure factor (Boolean)
struc_factor_sph  : F              #*! Accumulate sph. struc. factor (Boolean)
onep_density_mat  : F              #*! Accumulate 1p density matrix (Boolean)
twop_density_mat  : F              #*! Accumulate 2p density matrix (Boolean)
cond_fraction     : F              #*! Accumulate cond fraction (Boolean)
dipole_moment     : F              #*! Accumulate elec. dipole moment (Boolean)
expval_cutoff     : 30.d0 hartree  #*! G vector cutoff for expval (Physical)
permit_den_symm   : F              #*! Symmetrize QMC charge data (Boolean)
qmc_density_mpc   : F              #*! Use QMC density in MPC int (Boolean)
int_sf            : F              #*! Calc ee int from strucfac (Boolean)

# BLOCK INPUT
EOF

close(XYZ);

###################
# useful routines #
##################
# count valence electrons
sub atom_a2n{
  switch($_[0]){
    case "H" {return  1}
    case "B" {return  5}
    case "C" {return  6}
    case "N" {return  7}
    case "O" {return  8}
    case "F" {return  9}
    case "Si"{return 14}
    case "P" {return 15}
    case "S" {return 16}
    case "Cl"{return 17}
    case "Ge"{return 32}
    case "As"{return 33}
    case "Se"{return 34}
    case "Br"{return 35}
    case "Sn"{return 50}
    case "Sb"{return 51}
    case "Te"{return 52}
    case "I" {return 53}
    case "1" {return "H" }
    case "5" {return "B" }
    case "6" {return "C" }
    case "7" {return "N" }
    case "8" {return "O" }
    case "9" {return "F" }
    case "14" {return "Si"}
    case "15" {return "P" }
    case "16" {return "S" }
    case "17" {return "Cl"}
    case "32" {return "Ge"}
    case "33" {return "As"}
    case "34" {return "Se"}
    case "35" {return "Br"}
    case "50" {return "Sn"}
    case "51" {return "Sb"}
    case "52" {return "Te"}
    case "53" {return "I" }
    #else    {
    #  print "element: $_[0] not found\n";
    #  die;
    #}
  }
}

sub atom_ve{
  switch($_[0]){
    case "H" {return 1}
    case "B" {return 3}
    case "C" {return 4}
    case "Si"{return 4}
    case "Ge"{return 4}
    case "Sn"{return 4}
    case "N" {return 5}
    case "P" {return 5}
    case "As"{return 5}
    case "Sb"{return 5}
    case "O" {return 6}
    case "S" {return 6}
    case "Se"{return 6}
    case "Te"{return 6}
    case "F" {return 7}
    case "Cl"{return 7}
    case "Br"{return 7}
    case "I" {return 7}
    else     {
      print "element: $_[0] not found\n";
      die;
    }
  }
}

