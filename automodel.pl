#!/usr/bin/perl

#  easyAmber: A comprehensive toolbox to automate the molecular dynamics simulation of proteins
#  https://biokinet.belozersky.msu.ru/easyAmber
#
#  Copyright @ 2020-2021 Dmitry Suplatov, Yana Sharapova, Vytas Svedas 
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or of
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use File::Basename;
use 5.010;
use POSIX ":sys_wait_h";
use File::Copy;
$| = 1;

print "\n";
print "    AUTOMODEL master-script part of the easyAMBER toolkit\n";
print "\n"; 
print "    Prepares the AMBER parameter files for full-atom molecular system\n";
print "    of a protein or a protein-ligand complex in a cubic water box\n";
print "    Script v. 2.0 2019-12-15\n";
print "\n";
print " If you find this software or its results useful please cite our work\n";
print " More information available at https://biokinet.belozersky.msu.ru/easyAMBER\n";
print "\n";
print " Don't hesitate to send your questions and feedback to d.a.suplatov\@belozersky.msu.ru\n";
print "\n";

my $FF = "FF14SB";
my $GAFF = "GAFF";
my $WATER = "TIP4PEW";
my $BOX = 12;
my $ADDNA = 50;
my $ADDCL = $ADDNA;
my $IONS_FAST = 0;
my $IONS="";
my $FFFILE="";
my $BONDS="";
my $ATOMS="";

if (@ARGV<2) {
    print "BASIC USAGE\n";
    print "Usage: $0 <protein.pdb> <outputmask>\n";
    print "\n";
    print "OPTIONS\n";
    print "              Force field control (Default=$FF and $GAFF)\n";
    print "              ff=ff14sb          #Use the FF14SB force field\n";
    print "              ff=ff15ipq         #Use the FF15IPQ force field (compatible only with the SPCE/B water model)\n";
    print "              gaff=off           #Do not load the General Amber Force-Field\n";
    print "              gaff=true          #Load the General Amber Force-Field\n";
    print "              gaff2=true         #Load the General Amber Force-Field-2\n";
    print "              ffmisc=<file>      #Read custom force field commands for leap from a file\n";
    print "\n";
    print "              Water control (Default=$WATER)\n";
    print "              water=no           #no water\n";
    print "              water=tip3p        #Solvate the protein with TIP3P waters\n";
    print "              water=tip4pew      #Solvate the protein with TIP4P-Ew waters\n";
    print "              water=spceb        #Solvate the protein with SPCE/B waters (for the FF15IPQ force field)\n";
    print "\n";
    print "              Box size control\n";
    print "              box=n              #Minimum distance between any atom of the protein and the edge of the solvation box (Default=$BOX)\n";
    print "\n";
    print "              New atom types\n";
    print "              atoms=<file>       #Read in the custom leap commands (e.g. specify new Atom Types)\n";   
    print "\n";
    print "              Bonds\n";
    print "              bonds=<file>       #Read in the parameter file to set additional covalent bonds within the protein (e.g. SS-bonds)\n";   
    print "\n";
    print "              Charge/Ionic strength control\n";
    print "              ions=type          #Load \"hfe\" set for divalent ions to reproduce Hydration Free Energy (HFE) values\n";
    print "                                  Load \"iod\" set for divalent ions to reproduce Ion-Oxygen Distance (IOD) values\n";
    print "                                  Load \"cm\" set for divalent ions to reproduce Coordination Number (CN) values of the first solvation shell\n";
    print "              na=n               #Add Na+ ions to the solvation box (Default=$ADDNA)\n";
    print "              cl=n               #Add Cl- ions to the solvation box (Default=$ADDCL)\n";
    print "              fastions=true      #Add ions using a fast addions routine instead of slow addions2\n";
    print "\n";
    print "NOTES\n";
    print "      NOTE 1. Disulfide bridges will not be created automatically (even if CYX residues have already been set in the PDB)\n";
    print "              Define S-S bonds manually in a file and provide it using the bonds=SS_BONDS.dat command-line flag\n";
    print "      NOTE 2. To set new atom types use the atoms=<file> command-line flag\n";
    print "      NOTE 3. To load auxiliary force-fields use the ffmisc=<file> command-line flag.\n";
    print "      The deals are available at https://biokinet.belozersky.msu.ru/easyAmber#techdoc\n";
    print "\n";
    exit 0;
}

my $gaff_check = 0;

foreach my $i (2 .. $#ARGV) {
    if ($ARGV[$i] =~ /^ff=ff14sb/) {$FF = "FF14SB"; next;}
    if ($ARGV[$i] =~ /^ff=ff15ipq/) {$FF = "FF15IPQ"; next;}
    if ($ARGV[$i] =~ /^gaff=off/) {$GAFF = "OFF"; $gaff_check++; next;}
    if ($ARGV[$i] =~ /^gaff=true/) {$GAFF = "GAFF"; $gaff_check++; next;}
    if ($ARGV[$i] =~ /^gaff2=true/) {$GAFF = "GAFF2"; $gaff_check++; next;}
    if ($ARGV[$i] =~ /^ffmisc=(.*)/) {$FFFILE = $1; $FFFILE =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^water=spceb$/) {$WATER = "SPCE/B"; next;}
    if ($ARGV[$i] =~ /^water=tip3p$/) {$WATER = "TIP3P"; next;}
    if ($ARGV[$i] =~ /^water=tip4pew$/) {$WATER = "TIP4PEW"; next;}
    if ($ARGV[$i] =~ /^water=no$/) {$WATER = "NULL"; next;}
    if ($ARGV[$i] =~ /^box=(.*)/) {$BOX = $1; $BOX =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^na=(.*)/) {$ADDNA = $1; $ADDNA =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^cl=(.*)/) {$ADDCL = $1; $ADDCL =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^ions=(.*)/) {$IONS = $1; $IONS =~ s/\s+//g; next;}    
    if ($ARGV[$i] =~ /^fastions=true$/) {$IONS_FAST = 1; next;}    
    if ($ARGV[$i] =~ /^bonds=(.*)$/) {$BONDS = $1; $BONDS =~ s/\s+//g; next;}        
    if ($ARGV[$i] =~ /^atoms=(.*)$/) {$ATOMS = $1; $ATOMS =~ s/\s+//g; next;}        
    
    print "Error: Option $ARGV[$i] was not recognized\n";
    exit 1;
}

if ($gaff_check>1) {
    print "Error: Use of multiple command-line flags to define the GAFF settings is invalid\n";
    exit 1;
}

my $PROTEIN = $ARGV[0];
my $OUTPUT = $ARGV[1];

my $AMBER = `echo \$AMBERHOME`; chomp($AMBER);
my $ANTECHAMBER = "$AMBER/bin/antechamber";
my $PARMCHK     = "$AMBER/bin/parmchk2"; #Amber18 - parmchk2 only !!!
my $LEAP        = "$AMBER/bin/tleap";
my $AMBPDB      = "$AMBER/bin/ambpdb";

if ($FF eq "FF15IPQ" && $WATER ne "NULL") {
    $WATER="SPCE/B";
    print "Warning: The water model has been automatically set to the $WATER\n";    
}

print "Input: Protein         $PROTEIN\n";
print "Input: Output filemask $OUTPUT\n";
print "Input: Force field is  $FF\n";
print "Input: GAFF is         $GAFF\n";
print "Input: Misc. FF file   $FFFILE\n";
print "Input: Water model is  $WATER\n";
print "Input: Atoms file is   $ATOMS\n";
print "Input: Box size is     $BOX\n";
print "Input: Bonds file is   $BONDS\n";
print "Input: Na+ ions        $ADDNA\n";
print "Input: Cl- ions        $ADDCL\n";
if ($IONS_FAST) {
    print "Input: Fast ions       true\n";
}


print "\n";

#-----------------------------------------------------------------------------------------------------

print "Info: Checking for AmberTools binaries\n";
if ($AMBER eq "" || ! -d $AMBER) {die "Error: Amber home folder not set or does not exist at \"$AMBER\"\n";} else {print "Info: Amber home folder $AMBER ... OK\n";}
if (! -x $ANTECHAMBER || ! -f $ANTECHAMBER) {die "Error: The antechamber binary does not exist at \"$ANTECHAMBER\" or is not executable\n";} else {print "Info: The antechamber binary $ANTECHAMBER ... OK\n";}
if (! -x $PARMCHK || ! -f $PARMCHK) {die "Error: The parmchk binary does not exist at \"$PARMCHK\" or is not executable\n";} else {print "Info: The parmchk binary $PARMCHK ... OK\n";}
if (! -x $LEAP || ! -f $LEAP) {die "Error: The tleap binary does not exist at \"$LEAP\" or is not executable\n";} else {print "Info: The tleap binary $LEAP ... OK\n";}
if (! -x $AMBPDB || ! -f $AMBPDB) {die "Error: The ambpdb binary does not exist at \"$AMBPDB\" or is not executable\n";} else {print "Info: The ambpdb binary $AMBPDB ... OK\n";}
print "\n";

#-----------------------------------------------------------------------------------------------------

my $output_leap = "$OUTPUT.leap";
my $output_leaplog = "$OUTPUT.leap_log";
my $output_prmtop = "$OUTPUT.prmtop";
my $output_inpcrd = "$OUTPUT.inpcrd";
my $output_pdb = "$OUTPUT.pdb";

print "Info: Checking files from previous runs\n";
my @remove = ($output_leap, $output_leaplog, $output_prmtop, $output_inpcrd, $output_pdb);
my @files = glob("$OUTPUT*");
my $last_bkp_i = 0;
foreach my $file (@remove) {
    foreach my $dirfile (@files) {
        next if ($dirfile =~ m/^\./);
        if ($dirfile =~ /$file\_bkp(\d+)/) {
            if ($last_bkp_i < $1) {$last_bkp_i = $1;}                
        }                        
    }    
}
my $i = $last_bkp_i+1;
foreach my $file (@remove) {
    if (-f $file) {
        my $file_copy = "$file\_bkp$i";
        move($file, $file_copy);
        print "Info: Creating a backup copy of the file $file as $file_copy\n";
    }
}
print "\n";
#-----------------------------------------------------------------------------------------------------

$PROTEIN =~ /^(.*)\.(pdb)$/;
my $protein_mask = $1;
my $format = $2;

my $leap_string = "";
my @complex = ("PROT");

#Custom commands 
if ($ATOMS) {
    print "\nInfo: Appending user-defined custom commands from the additional parameter file $ATOMS:\n";
    open(A, $ATOMS) || die "Error: Failed to open the additional bonds parameter file $ATOMS\n";
    while(<A>){
	chomp;
	my $line = $_;
	$leap_string .= "$line\n";
	print "Info: 	$line\n";
    }
    close(A);
    print "Info: Closing the additional bonds parameter file $ATOMS\n\n";
}

#$leap_string .= "addPdbAtomMap {\n";
#$leap_string .= "{ \"O1P\" \"OP1\" }\n";
#$leap_string .= "{ \"O2P\" \"OP2\" }\n";
#$leap_string .= "{ \"O3P\" \"OP3\" }\n";
#$leap_string .= "{ \"O4P\" \"OP4\" }\n";
#$leap_string .= "}\n";

if ($GAFF eq "OFF") {
    print "Info: GAFF will not be loaded\n";
}
elsif ($GAFF eq "GAFF") {
    print "Info: Loading GAFF\n";
    $leap_string .= "source leaprc.gaff\n";
}
elsif ($GAFF eq "GAFF2") {
    print "Info: Loading GAFF-2\n";
    $leap_string .= "source leaprc.gaff2\n";
}
else {
    print "Error: Invalid GAFF keyword\n";
    exit 1;
}

if ($FF eq "FF14SB") {    
    #my $defaultff = "leaprc.ff14SB";
    my $defaultff = "leaprc.protein.ff14SB";
    print "Info: Loading the $defaultff forcefield\n";
    $leap_string .= "source $defaultff\n";
    
    print "Info: Loading the forcefield parameters for water ($WATER) and ions\n";
    if ($WATER eq "TIP3P") {
        #$leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionsjc_tip3p\n";    
        $leap_string .= "source leaprc.water.tip3p\n";    
    }
    elsif ($WATER eq "TIP4PEW") {
        $leap_string .= "source leaprc.water.tip4pew\n";
        #$leap_string .= "loadOff tip4pbox.off\n";
        #$leap_string .= "loadOff tip4pewbox.off\n";
        #$leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.tip4pew\n";    
        #$leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionsjc_tip4pew\n";
        
        $leap_string .= "HOH=TP4\n"; #Treat all residues named HOH as TP4 waters (default is TIP3P)
        $leap_string .= "WAT=TP4\n"; #Treat all residues named WAT as TP4 waters (default is TIP3P)    
    }
    elsif ($WATER eq "SPCE/B") {
	$leap_string .= "source leaprc.water.spceb\n";
    }
    elsif ($WATER eq "NULL") {
    }
    else {
        print "Error: Failed to load parameters for water model $WATER\n";
        exit 1;
    }
    
    if ($IONS eq "iod") {
        print "Info: Loading IOD ion parameters for divalent ions (same parameters for all waters)\n";
        $leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionslrcm_iod"
    }
    elsif ($IONS eq "cm") {
        if ($WATER eq "TIP3P") {
        print "Info: Loading CM ion parameters for divalent ions for tip3p water model\n";
        $leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionslrcm_cm_tip3p\n";
    }
    elsif ($WATER eq "TIP4PEW") {
            print "Info: Loading CM ion parameters for divalent ions for tip4pew water model\n";
        $leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionslrcm_cm_tip4pew\n";
        }
    }
    elsif ($IONS eq "hfe") {
        if ($WATER eq "TIP3P") {
            print "Info: Loading HFE ion parameters for divalent ions for tip3p water model\n";
        $leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionslrcm_hfe_tip3p\n";
    }
    elsif ($WATER eq "TIP4PEW") {
            print "Info: Loading HFE ion parameters for divalent ions for tip4pew water model\n";
        $leap_string .= "loadamberparams $AMBER/dat/leap/parm/frcmod.ionslrcm_hfe_tip4pew\n";
        }
    }
}
elsif ($FF eq "FF15IPQ") {
    if ($WATER eq "NULL") {
    }
    if ($WATER ne "NULL" && $WATER ne "SPCE/B") {
        print "Error: The selected water model $WATER is not compatible with the selected force field $FF\n";
        exit 1;
    }
    else {
	$leap_string .= "source leaprc.protein.ff15ipq\n";
	$leap_string .= "source leaprc.water.spceb\n";
    }
}
else {
   print "Error: The selected forcefield not supported\n";
   exit 1;
}

#Add SS-bonds here
if ($FFFILE) {
    print "\nInfo: Appending force field parameters to the LEAP file from the additional file $FFFILE:\n";
    open(F, $FFFILE) || die "Error: Failed to open the additional force field parameter file $FFFILE\n";
    while(<F>){
	chomp;
	my $line = $_;
	$leap_string .= "$line\n";
	print "Info: 	$line\n";
    }
    close(F);
    print "Info: Closing the additional force field parameter file $BONDS\n\n";
}

print "Info: Preparing parameters for the protein / protein-ligand complex\n";
print "Info: Reading the protein file from $PROTEIN\n";
my $pdb = pdbReformat(readpdb($PROTEIN));

#-----------------------------------------------------------------------------------------------------

print "Info: Checking for ligands in the current folder and all subfolders\n";
my %ligands = ();
foreach my $file (@{rfiles("./")}) {

    if ($file !~ m/\.mol2$/ && $file !~ m/\.frcmod$/ && $file !~ m/\.lib$/) {next;}
    
    if ($file =~ /_amber\.mol2$/) {
        my $filetitle = fileparse($file, ".mol2");
        $ligands{$filetitle}{"mol2"} = $file;
    }
    elsif ($file =~ /_amber\.frcmod$/) {
        my $filetitle = fileparse($file, ".frcmod");
        $ligands{$filetitle}{"frcmod"} = $file;
    }
    elsif ($file =~ /_amber\.lib$/) {
        my $filetitle = fileparse($file, ".lib");
        $ligands{$filetitle}{"lib"} = $file;
    }
}

foreach my $ligand (keys %ligands) {    
    if ($ligands{$ligand}{"lib"} && $ligands{$ligand}{"frcmod"} && $ligands{$ligand}{"mol2"}) {
        my $ligand_name1 = libtitle($ligands{$ligand}{"lib"});
        my $ligand_name2 = mol2title($ligands{$ligand}{"mol2"});
        if ($ligand_name1 ne $ligand_name2) {
            print "Error: Ligand titles in the LIB and MOL2 files do not match ($ligand_name1 $ligand_name2)\n";
            exit 1;
        }
        
        print "Info: Adding ligand $ligand_name1\n";
        print "      Amber library:    ".$ligands{$ligand}{"lib"}."\n";
        print "      Amber parameters: ".$ligands{$ligand}{"frcmod"}."\n";
        print "      MOL2 file:        ".$ligands{$ligand}{"mol2"}."\n";
        
        my $ligand_occurences = checkLigand($pdb, $ligand_name1);
        if ($ligand_occurences == 0) {
            print "      ZERO instances of the ligand found in PDB\n";
            print "      Loading the ligand position from the MOL2 file\n";
            
            my $molid = $#complex;            
            $leap_string .= "loadoff ".$ligands{$ligand}{"lib"}."\n";
            $leap_string .= "loadamberparams ".$ligands{$ligand}{"frcmod"}."\n";
            $leap_string .= "MOL$molid = loadmol2 ".$ligands{$ligand}{"mol2"}."\n";
            push(@complex, "MOL$molid");
        }
        else {
            print "      Instances of the ligand found in PDB: ".$ligand_occurences."\n";
            $leap_string .= "loadoff ".$ligands{$ligand}{"lib"}."\n";
            $leap_string .= "loadamberparams ".$ligands{$ligand}{"frcmod"}."\n";        
        }        
    }
    
    elsif ($ligands{$ligand}{"lib"} && $ligands{$ligand}{"frcmod"}) {
        my $ligand_name = libtitle($ligands{$ligand}{"lib"});        
        print "Info: Adding ligand $ligand_name\n";
        print "      Amber library:    ".$ligands{$ligand}{"lib"}."\n";
        print "      Amber parameters: ".$ligands{$ligand}{"frcmod"}."\n";        
        my $ligand_occurences = checkLigand($pdb, $ligand_name);        
        if ($ligand_occurences == 0) {print " !!!  Warning: ZERO instances of the ligand found in PDB\n";}
        else {print "      Instances of the ligand found in PDB: ".$ligand_occurences."\n";}

        $leap_string .= "loadoff ".$ligands{$ligand}{"lib"}."\n";
        $leap_string .= "loadamberparams ".$ligands{$ligand}{"frcmod"}."\n";
    }
    
    elsif ($ligands{$ligand}{"mol2"} && $ligands{$ligand}{"frcmod"}) {
        my $ligand_name = mol2title($ligands{$ligand}{"mol2"});
        print "Info: Adding ligand $ligand_name\n";
        print "      MOL2 file:        ".$ligands{$ligand}{"mol2"}."\n";
        print "      Amber parameters: ".$ligands{$ligand}{"frcmod"}."\n";        
        my $ligand_occurences = checkLigand($pdb, $ligand_name);        
        if ($ligand_occurences == 0) {print "      ZERO instances of the ligand found in PDB\n";}
        else {print " !!!! Warning: ".$ligand_occurences." instances of the ligand found in PDB but still we will load MOL2 coordinates\n";}
        print "      Loading the ligand position from the MOL2 file\n";
        
        my $molid = $#complex;
        $leap_string .= "MOL$molid = loadmol2 ".$ligands{$ligand}{"mol2"}."\n";
        $leap_string .= "loadamberparams ".$ligands{$ligand}{"frcmod"}."\n";
        push(@complex, "MOL$molid");
    }    
}

#-----------------------------------------------------------------------------------------------------
#Load the protein now (after ligand parameters have been loaded -- in case the ligands are included in the protein structure)

my $cyx_occurences = checkLigand($pdb, "CYX");
if ($cyx_occurences) {print "\nWarning: Your PDB has $cyx_occurences CYX residues. You must set the SS-bonds manually in a separate paremeter file!\n\n";}

my $cys_occurences = checkLigand($pdb, "CYS");
if ($cys_occurences) {print "\nWarning: Your PDB has $cys_occurences CYS residues. Cysteins can be CYS, CYX or CYM. You might want to check this assignment again\n\n";}

$leap_string .= "PROT = loadpdb $PROTEIN\n";

#Add SS-bonds here
if ($BONDS) {
    print "\nInfo: Appending new lines to the LEAP file from the additional bonds parameter file $BONDS:\n";
    open(B, $BONDS) || die "Error: Failed to open the additional bonds parameter file $BONDS\n";
    while(<B>){
	chomp;
	my $line = $_;
	$leap_string .= "$line\n";
	print "Info: 	$line\n";
    }
    close(B);
    print "Info: Closing the additional bonds parameter file $BONDS\n\n";
}

if ($#complex == 0) {
    $leap_string .= "complex = $complex[0]\n";
}
else {
    $leap_string .= "complex = combine {";
    foreach my $str (@complex) {
        $leap_string .= " $str";
    }
    $leap_string .= "}\n";    
}
$leap_string .= "check complex\n";

print "Info: Preparing the water box of size $BOX A and compensating for non-zero charge of the system (if any)\n";

if ($WATER eq "TIP3P" && $BOX !=0) {
    $leap_string .= "solvateBox complex TIP3PBOX $BOX\n";
}
elsif ($WATER eq "TIP4PEW" && $BOX !=0) {
    $leap_string .= "solvateBox complex TIP4PEWBOX $BOX\n";
}
elsif ($WATER eq "SPCE/B" && $BOX !=0) {
    $leap_string .= "solvateBox complex SPCBOX $BOX\n";
}
elsif ($WATER eq "NULL" || $BOX ==0) {
    print "Info: No solvation (vacuum)\n";
}
else {
    print "Error: Failed to set box for water model $WATER\n";
    exit 1;
}

if ($ADDNA !=0) {
    print "Info: Counterbalancing the charge with positive inorganic ions (Na+=$ADDNA)\n";
    if ($IONS_FAST) {$leap_string .= "addIons complex Na+ $ADDNA\n";}
    else {$leap_string .= "addIons2 complex Na+ $ADDNA\n";}
}

if ($ADDCL !=0) {
    print "Info: Counterbalancing the charge with negative inorganic ions (Cl-=$ADDCL)\n";
    if ($IONS_FAST) {$leap_string .= "addIons complex Cl- $ADDCL\n";}
    else {$leap_string .= "addIons2 complex Cl- $ADDCL\n";}
}

print "Info: Preparing the leap run file\n";
$leap_string .= "saveamberparm complex $output_prmtop $output_inpcrd\n";
$leap_string .= "quit\n";
open (FO, ">$output_leap") || die "Error: Failed to open file $output_leap for writing\n";
print FO $leap_string;
close(FO);

print "\n";
#-----------------------------------------------------------------------------------------------------

print "Info: Running leap\n";
my $pid = fork();
if (not defined $pid) {
    print "Error: Failed to fork a thread for leap\n";
    exit 1;
}
#Execute leap in the forked thread
if (not $pid) {
    system("$LEAP -f $output_leap > $output_leaplog");
    exit;
}
#Monitor leap errors in the main process
print "------------------------leap warning/error log begins------------------------------------\n";
my %leap_error=();
while (1) {
    print "Waiting for update ... \r";
    
    sleep(1);    

    my $res = waitpid($pid, WNOHANG); #Request status code from the child thread    
    open(FI, "$output_leaplog") or die "Error: Failed to open file $output_leaplog for reading\n";
    my $lineid = 0;
    my $printnext = 0;
    while (<FI>) {
        chomp;
        $lineid++;
        if (/[Ee][Rr][Rr][Oo][Rr]/ || /[Ww][Aa][Rr][Nn][Ii][Nn][Gg]/ || /[Ff][Aa][Tt][Aa][Ll]/ || /Added missing heavy atom/ || /Leap added.*missing atoms according to residue templates/ || /bond:/ || /usage:/ ) {
            my $line = $_;
            if (!$leap_error{$lineid}) {
                print $line."\n";
                $leap_error{$lineid}=1;
                $printnext=$lineid+1;
                next;
            }            
        }
        if ($printnext == $lineid && !$leap_error{$lineid}){
            my $line = $_;
            print $line."\n";
            $leap_error{$lineid}=1;
            $printnext=0;
        }
    }
    close(FI);
    
    if ($res == -1) {
        print "Error: Some error occurred while forking leap". $? ."\n";
        exit();
    }
    if ($res) {last;} #The child has completed
    if (!$res) {next;} #The child is still running
}
print "\nSee the output file $output_leaplog for the full leap log\n";
print "------------------------leap warning/error log ends--------------------------------------\n";

if (!-f "$output_prmtop" || !-f "$output_inpcrd" || !-s "$output_prmtop" || !-s "$output_inpcrd") {
    print "Error: The required amber molecular topology, force field parameters, atom and residue names file ($output_prmtop)\n";
    print "       and/or initial coordinates file ($output_inpcrd) are invalid or do not exist\n";
    exit 1;
}

print "\n";

#-----------------------------------------------------------------------------------------------------

print "Info: Writing the molecular system to $output_pdb\n";
print "------------------------ambpdb log starts--------------------------------------\n";
system("$AMBPDB -p $output_prmtop < $output_inpcrd > $output_pdb");
print "-------------------------ambpdb log ends---------------------------------------\n";
if (!-f $output_pdb || !-s $output_pdb) {
    print "Error: Failed to write the molecular system in PDB format to $output_pdb";
}
print "\n";

#-----------------------------------------------------------------------------------------------------

print "\n";
print "Warning: Check the charge of the molecular system and re-run to add a different amount of ions if necessery\n";
print "Instruction: In order to set the number if ions to be added in order to simulate a particular ionic strenth\n";
print "             you can follow the procedure discribed below\n";
print "             1. Create a copy of the initial unsolvated PDB by removing all non-standart residues\n";
print "             2. Run the vmd_autoionize.tcl - this will create a solvated and ionized pdb box \n";
print "             (use the same box size when running both this script and vmd_autoionize.tcl)\n";
print "             3. Calculate the number of Na+ (SOD) and Cl- (CLA) ions \n";
print "             4. Re-run this script with the new setup\n";
print "\n";
print "Info: Automatic protein / protein-ligand model preparation has completed\n";
print "Info: The amber molecular topology, parameters, names file: $output_prmtop\n";
print "Info: The initial coordinates files:                        $output_inpcrd\n";
print "Info: The molecular system in PDB format:                   $output_pdb\n";
print "Info: Prefix to launch autorun.pl:                          $OUTPUT\n";
print "Done!\n";
print "\n";

#===========================================================================

#Recersively collect all file names in all subdirectories starting from $root
sub rfiles{
    my $root = shift;
    my @return = ();
    if ( -f $root ) {push(@return, $root); return \@return;}
    opendir (DIR, $root) or die "Error: Failed to open input folder $root\n";
    my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
    closedir DIR;

    foreach my $file(@files) {
      my $subreturn = rfiles("$root/$file");
      foreach my $subfile (@$subreturn) {
        push(@return, $subfile);
      }
    }
    return \@return;
}

sub mol2title {
    my $mol2file = shift;
    open(F, $mol2file) or die "Error: Failed to open file $mol2file\n";
    my $line_no = 0;
    my $title="";
    while (<F>) {
        chomp;
        $line_no++;
        if ($line_no==2) {
            $title = $_;
            last;
        }        
    }    
    close(F);
    
    if ($title eq "") {
        print "Error: Unvalid (emtpy) title in MOL2 file $mol2file\n";
        exit 1;
    }
    
    return $title;    
}

sub libtitle {
    my $libfile = shift;
    open(F, $libfile) or die "Error: Failed to open file $libfile\n";
    
    my $readline = 0;
    my $title="";
    while (<F>) {
        chomp;
        if (/^!!index array str$/) {
            $readline = 1;
            next;
        }
        if ($readline) {
            if (/"(.*)"/) {
                $title = $1;
            }
            else {
                print "Error: Invalid format of the ligand title in library file $libfile\n";
                exit 1;
            }
            $readline=0;
            last;0
        }        
    }    
    close(F);
    
    if ($title eq "") {
        print "Error: Unvalid (emtpy) title in library file $libfile\n";
        exit 1;
    }
    
    return $title;    
}


#Read in a pdb file
sub readpdb {
    my $input = shift;
    
    my %pdb = ();
    open(FI, $input) or die "ERROR: Can not open file $input for reading\n\n";
    
    #Read in the PDB file
    
    while(<FI>)
    {
        chomp;                
        
        #Remove Windows ^M sign
        s/\r//g;
        
        my $keyword = substr($_,0,6);    # look for ATOM record
        if ($keyword eq "ATOM  " || $keyword eq "HETATM") {
        
            my $length = length $_;
            if ($length < 54) {                
                print "\n$_\n";
                print "ERROR: ATOM and HETATM fields must be at least 54 characters in length (found $length))\n";
                pdbFormatError();
            }
            
            my $num = substr($_,6,5);
            $num =~ s/\s+//g;
            while ($pdb{$num}) {
                $num++;
            }

            $pdb{$num}{"head"} = substr($_,0,6);
            $pdb{$num}{"atom"} = substr($_,12,4);
            my $type = $pdb{$num}{"atom"};
            $type=substr($type,0,1);
            $pdb{$num}{"type"} = $type;
            $pdb{$num}{"res"} = substr($_,16,4);
            $pdb{$num}{"chainid"} = substr($_,21,1);
            $pdb{$num}{"resn"} = substr($_,22,4);
            $pdb{$num}{"xcor"} = substr($_,30,8);
            $pdb{$num}{"ycor"} = substr($_,38,8);
            $pdb{$num}{"zcor"} = substr($_,46,8);
            #$pdb{$num}{"occ"} = substr($_,54,6);
            #$pdb{$num}{"bfac"} = substr($_,60,6);
            #$pdb{$num}{"segid"} = substr($_,72,4);
            #$pdb{$num}{"elm"} = substr($_,76,2);
            #$pdb{$num}{"charge"} = substr($_,78,2);
            
            #Remove spaces
            foreach my $key (keys %{$pdb{$num}}) {
                $pdb{$num}{$key} =~ s/\s+//g;
            }
        }
    }
    close(FI);
    
    return \%pdb;
}

sub pdbReformat {
    my $pdb = shift;    
    
    my %hash = ();
    foreach my $num (sort {$a <=> $b} keys %{$pdb}) {
        my $chain = $pdb->{$num}->{"chainid"};
        my $resname = $pdb->{$num}->{"res"};
        my $resid = $pdb->{$num}->{"resn"};
        my $atom = $pdb->{$num}->{"atom"};
        
        my $x = $pdb->{$num}->{"xcor"};
        my $y = $pdb->{$num}->{"ycor"};
        my $z = $pdb->{$num}->{"zcor"};
                
        my @data = ($x, $y, $z);
        
        if ($hash{$chain}{$resname}{$resid}{$atom}) {
            print "Warning: Overwriting duplicate coordinates for atom $chain:$resname:$resid:$atom}\n";
        }
        
        $hash{$chain}{$resname}{$resid}{$atom} = \@data;        
    }
    
    return \%hash;
}

sub checkLigand {
    my $pdb = shift;
    my $ligand = shift;
    
    my $count = 0;
    foreach my $chain (sort keys %{$pdb}) {
        if ($pdb->{$chain}->{$ligand}) {
            my $chainCount = scalar keys %{$pdb->{$chain}->{$ligand}};
            $count += $chainCount;
        }        
    }
    return $count;
}
