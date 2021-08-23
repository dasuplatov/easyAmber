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
$| = 1;

print "\n";
print "    AUTOLIGAND master-script part of the easyAMBER toolkit\n";
print "\n";
print "    Prepare the AMBER parameter files for a custom low-molecular-weight ligand\n";
print "    Script v. 2.0.2 2021-01-09\n";
print "\n";
print " If you find this software or its results useful please cite our work\n";
print " More information available at https://biokinet.belozersky.msu.ru/easyAMBER\n";
print "\n";
print " Don't hesitate to send your questions and feedback to d.a.suplatov\@belozersky.msu.ru\n";
print "\n";


if (@ARGV!=3) {    
    print "Usage: $0 <ligand.pdb> charge=<integer> ff=<gaff/gaff2/amber>\n";
    print "Usage: $0 <ligand.mol> charge=<integer> ff=<gaff/gaff2/amber>\n";    
    print "NB!: Launch in the same directory where the ligand is located and use local path\n";
    print "\n";
    exit 0;
}

my $LIGAND = $ARGV[0];
my $CHARGE = "N/A";
my $FF = "N/A";

foreach my $i (1 .. $#ARGV) {
    if ($ARGV[$i] =~ /^charge=(.*)/) {$CHARGE = $1; $CHARGE =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^ff=(.*)/) {$FF = $1; $FF =~ s/\s+//g; next;}    
}

my $AMBER = `echo \$AMBERHOME`; chomp($AMBER);
my $ANTECHAMBER = "$AMBER/bin/antechamber";
my $PARMCHK     = "$AMBER/bin/parmchk";
my $PARMCHK2     = "$AMBER/bin/parmchk2";
my $LEAP        = "$AMBER/bin/tleap";
my $OBMINIMIZE  = `which obminimize`; chomp($OBMINIMIZE);

my $defaultff = "";
#my $defaultffparm = "";
#my $defaultffparm_path = "";
if ($FF eq "gaff") {
    $defaultff = "leaprc.gaff";
}
elsif ($FF eq "gaff2") {
    $defaultff = "leaprc.gaff2";
}
elsif ($FF eq "amber") {
    $defaultff = "leaprc.ff14SB";
    #$defaultffparm = "frcmod.ff14SB";

    #$defaultffparm_path = "";# `find $AMBER -name \"$defaultffparm\"`; chomp($defaultffparm_path);
    #if (! -f $defaultffparm_path) {
	#print "Error: Forcefield parameter file $defaultffparm was not found or multiple instances found in $AMBER\n";
	#exit 1;
    #}
}
else {
    print "Error: Force field invalid or unset ($FF)\n";
    exit 1;
}

if ($CHARGE !~ /^-?\d+$/) {
    print "Error: Charge invalid or unset ($CHARGE)\n";
    exit 1;
}

print "Input: Ligand        $LIGAND\n";
print "Input: Ligand charge $CHARGE\n";
#if ($FF eq "amber") {
    print "Input: Forcefield    $FF ($defaultff)\n"; #, $defaultffparm)\n";
#}
#elsif ($FF eq "gaff") {
#    print "Input: Forcefield    GAFF ($defaultff)\n";
#}
print "\n";

print "Info: Checking for AmberTools binaries\n";
if ($AMBER eq "" || ! -d $AMBER) {die "Error: Amber home folder not set or does not exist at \"$AMBER\"\n";} else {print "Info: Amber home folder $AMBER ... OK\n";}
if (! -x $ANTECHAMBER || ! -f $ANTECHAMBER) {die "Error: The antechamber binary does not exist at \"$ANTECHAMBER\" or is not executable\n";} else {print "Info: The antechamber binary $ANTECHAMBER ... OK\n";}

if (! -x $PARMCHK || ! -f $PARMCHK) 
{
    if (! -x $PARMCHK2 || ! -f $PARMCHK2) 
    {
        die "Error: The parmchk binary does not exist at \"$PARMCHK\" or is not executable\n";
    }
    $PARMCHK = $PARMCHK2;
} 
else {print "Info: The parmchk binary $PARMCHK ... OK\n";}

if (! -x $LEAP || ! -f $LEAP) {die "Error: The tleap binary does not exist at \"$LEAP\" or is not executable\n";} else {print "Info: The tleap binary $LEAP ... OK\n";}
print "\n";

print "Info: Checking for miscellaneous binaries\n";
if (! -x $OBMINIMIZE || ! -f $OBMINIMIZE) {die "Error: The obminimize binary does not exist at \"$OBMINIMIZE\" or is not executable\n";} else {print "Info: The obminimize binary $OBMINIMIZE ... OK\n";}
print "\n";

$LIGAND =~ /^(.*)\.(mol2|pdb)$/;
my $ligand_mask = $1;
my $format = $2;

my $OUTFOLDER = "AMBPAR_$ligand_mask";
if (-d "$OUTFOLDER") {
    print "Error: Folder $OUTFOLDER already exists and must be deleted prior to re-launching this script\n";
    exit 1;
}
else {mkdir $OUTFOLDER; system("cp $LIGAND $OUTFOLDER/$LIGAND");}
print "Info: Created folder $OUTFOLDER to store Amber files\n";

my $output = "$ligand_mask\_amber.mol2";

if (! makeParams($OUTFOLDER, $ligand_mask, $format, $CHARGE, $output)) {    
    print "Info: Attempting to perform structure optimization (you can try mopac manually instead)\n";
    
    print "Info: Running at most 2500 steps of conjugate gradients energy minimization\n";
    system("cd $OUTFOLDER; $OBMINIMIZE -cg -n 2500 $LIGAND > $ligand_mask\_aem.pdb 2> $ligand_mask\_em.log");
    if (!-f "$OUTFOLDER/$ligand_mask\_aem.pdb") {print "Error: Failed at energy minimization\n"; exit 1;}
    
    print "Info: Submitting new structure for parametrization\n";
    makeParams($OUTFOLDER, "$ligand_mask\_aem", "pdb", $CHARGE, $output);
}
if (-f "$OUTFOLDER/$output") {print "Info: Atomic charges written to $OUTFOLDER/$output\n";}
else {print "Error: Failed calculating atomic charges\n"; exit 1;}
print "\n";

print "Info: Preparing parameter file\n";
my $output_frcmod = "$ligand_mask\_amber.frcmod";

#if ($defaultffparm_path ne "") { 
#    system("cd $OUTFOLDER; $PARMCHK -p $defaultffparm_path -pf 1 -i $output -f mol2 -o $output_frcmod -a Y -w Y");
#}
#else {
    system("cd $OUTFOLDER; $PARMCHK -i $output -f mol2 -o $output_frcmod -a Y -w Y");
#}

if (-f "$OUTFOLDER/$output_frcmod") {print "Info: Amber parameters written to $OUTFOLDER/$output_frcmod\n";}
else {print "Error: Failed praparing Amber parameters"; exit 1;}
print "\n";

print "Info: Preparing library file in the $defaultff forcefield\n";
my $output_lib = "$ligand_mask\_amber.lib";
my $ligand_title = mol2title("$OUTFOLDER/$output");
my $leap = "source $defaultff\n";
$leap .= "$ligand_title = loadmol2 $output\n";
$leap .= "saveoff $ligand_title $output_lib\n";
$leap .= "quit\n";
open(FO, ">$OUTFOLDER/$ligand_mask.leap") or die "Error: Failed to open file $OUTFOLDER/$ligand_mask.leap for writing\n";
print FO $leap;
close (FO);
print "------------------------leap log begins------------------------------------\n";
system("cd $OUTFOLDER; $LEAP -f $ligand_mask.leap");
print "------------------------leap log ends------------------------------------\n";
if (-f "$OUTFOLDER/$output_lib") {print "Info: Amber library file written to $OUTFOLDER/$output_lib\n";}
else {print "Error: Failed praparing Amber library file"; exit 1;}
print "\n";

print "Info: The automatic parametrization has been completed\n";
print "\n";
print "Info: Please note that in order to complete parametrization of a non-standart residue\n";
print "Info: (i.e., introduced into the protein backbone) further manual actions are required\n";
print "Info: Zero-step, re-run with amber force-field (i.e., not gaff/gaff2)\n";
print "Info: First, delete the the library and frcmod files and rename the mol2 file\n";
print "Info:    bash>rm $OUTFOLDER/$output_frcmod $OUTFOLDER/$output_lib\n";
print "Info:    bash>mv $OUTFOLDER/$output $OUTFOLDER/$output\_bkp\n";
print "Info: Second, in case a tripeptide has been submitted as a ligand, run the following commands\n";
print "Info:    leap>xleap -s -f $defaultff #Load the xleap with amber forcefield\n";
print "Info:    leap>$ligand_title = loadmol2 $OUTFOLDER/$output\_bkp #Load the ligand\n";
print "Info:    leap>edit $ligand_title #Use the editor to remove atoms which do not belong to the query resudue\n";
print "Info: Third, set the N-terminal and C-terminal atoms to be used to connect the backbone\n";
print "Info:    leap>set $ligand_title head $ligand_title.$ligand_title.N #Set the N-terminal atom\n";
print "Info:    leap>set $ligand_title tail $ligand_title.$ligand_title.C #Set the C-terminal atom\n";
print "Info:    leap>desc $ligand_title #Check the results\n";
print "Info: Forth, rebuild lib and frcmod files\n";
print "Info:    leap>savemol2 $ligand_title $OUTFOLDER/$output 1 #Save the mol2 file (1 = AMBER atom types)\n";
print "Info:    leap>saveoff $ligand_title $OUTFOLDER/$output_lib #Save the lib file\n";
print "Info:    leap>quit #close tleap\n";
print "Info:    bash>$PARMCHK -i $OUTFOLDER/$output -f mol2 \\ \n";
print "Info:    bash>  -o $OUTFOLDER/$output_frcmod -a Y -w Y #use parmchk to rebuild the frcmod file\n";
print "\n";
print "Info: The automatic parametrization has been completed\n";
print "Info: If the query ligand is not covalently attached to the protein you can use the prepared file below straight-ahead\n";
print "Output: Ligand topology and charges: $OUTFOLDER/$output\n";
print "Output: Amber parameters:            $OUTFOLDER/$output_frcmod\n";
print "Output: Amber library:               $OUTFOLDER/$output_lib\n";
print "Done!\n";
print "\n";
exit 0;

sub makeParams {
    my $folder = shift;
    my $ligand_mask = shift;
    my $format = shift;
    my $charge = shift;
    my $output = shift; 
    
    my $ligand = $ligand_mask. "." .$format;
        
    print "Info: Calculating AM1-BCC charges for ligand $ligand with default parameters\n";
    
    print "\n";
    print "------------------------antechamber log begins------------------------------------\n";

    if ($FF eq "amber") {system("cd $folder; $ANTECHAMBER -at amber -i $ligand -fi $format -o $output -fo mol2 -c bcc -pf y -nc $charge");}
    elsif ($FF eq "gaff") {system("cd $folder; $ANTECHAMBER -at gaff -i $ligand -fi $format -o $output -fo mol2 -c bcc -pf y -nc $charge");}
    
    print "-------------------------antechamber log ends-------------------------------------\n";
    print "\n";
    if (-f "$folder/$output") {return 1;}
    print "Info: The attempt has failed\n";
    
    print "Info: Calculating AM1-BCC charges for ligand $ligand with -j 5 option\n";
    print "\n";
    print "------------------------antechamber log begins------------------------------------\n";

    if ($FF eq "amber") {system("cd $folder; $ANTECHAMBER -at amber -i $ligand -fi $format -o $output -fo mol2 -c bcc -j 5 -pf y -nc $charge");}
    elsif ($FF eq "gaff") {system("cd $folder; $ANTECHAMBER -at gaff -i $ligand -fi $format -o $output -fo mol2 -c bcc -j 5 -pf y -nc $charge");}

    print "-------------------------antechamber log ends-------------------------------------\n";
    print "\n";
    if (-f "$folder/$output") {return 1;}
    print "Info: The attempt has failed\n";
    
=cut
    print "Info: Calculating atomic charges for ligand $ligand with -c mul and -j 5 options\n";
    print "\n";
    print "------------------------antechamber log begins------------------------------------\n";

    if ($FF eq "amber") {system("cd $folder; $ANTECHAMBER -at amber -i $ligand -fi $format -o $output -fo mol2 -c mul -j 5 -pf y -nc $charge");}
    elsif ($FF eq "gaff") {system("cd $folder; $ANTECHAMBER -at gaff -i $ligand -fi $format -o $output -fo mol2 -c mul -j 5 -pf y -nc $charge");}

    print "-------------------------antechamber log ends-------------------------------------\n";
    print "\n";
    if (-f "$folder/$output") {return 1;}
=cut

    print "Info: Failed to calculate atomic charges for ligand $ligand\n";
    return 0;
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






