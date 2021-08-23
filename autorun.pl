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

my $TEMPERATURE = "300.0"; # 300 K
my $DT_STEP = "0.002"; # 2 fs
my $CUTOFF = 10; # Angstroms
my $RESTRAINTS_CONSTANT = 3; # in kcal/(mol Angstroms^2)
my $EQTIME = 2600; # default is 2.6 ns in total
my $EQSTEPS = 0; # to be calculated internally
my $EQDEC = 0.25; # in kcal/(mol Angstroms^2) 
my $RUNTIME = 500000; # 500 ns 

print "\n";
print "    AUTORUN master-script part of the easyAMBER toolkit\n";
print "\n";
print "    Operate the seven-step MD simulation pipeline\n";
print "    Script v. 2.1 2020-12-18\n";
print "\n";
print " If you find this software or its results useful please cite our work\n";
print " More information available at https://biokinet.belozersky.msu.ru/easyAMBER\n";
print "\n";
print " Don't hesitate to send your questions and feedback to d.a.suplatov\@belozersky.msu.ru\n";
print "\n";


if (@ARGV<1) {
    print "Usage: autorun.pl help #prints help\n";
    print "Usage: autorun.pl ex   #prints examples\n";
    print "\n";
    exit 0;
}

if (@ARGV==1 && $ARGV[0] eq "help") {
    
    print "Usage: autorun.pl <filemask> #Run all steps (default)\n";
    print "\n";
    print "Usage: autorun.pl <filemask> [options]  #Use options\n";
    print "Misc options: Print configuration files control (use only once to print the config files)\n";
    print "              print=step         #Print configuration for the selected step only and exit\n";
    print "              cutoff=<float>     #Set PME electrostatic cutoff (default=$CUTOFF)\n";
    print "              temp=<int>         #Set target temperature in K (default=$TEMPERATURE)\n";
    print "              rest=<float>       #Set the restraints constant in kcal/(mol Angstroms^2) (default=$RESTRAINTS_CONSTANT)\n";
    print "              eqdec=<float>      #Set the decrement to apply to the restraints constant in kcal/(mol Angstroms^2)\n";
    print "                                 # at each consequtive equilibration step to effectively regulate number of these steps (default=$EQDEC)\n";
    print "              eqtime=<int>       #Set the total time for equilibration steps (all substeps of the step #5) in ps (default=$EQTIME)\n";
    print "              runtime=<int>      #Set the total time for the production run (steps 6 or 7) in ps (default=$RUNTIME)\n";
    print "                                 #FYC: 1 ns = 1000 ps; 1300 ps = 1.3 ns; 500000 ps = 500 ns; 1000000 ps = 1 microsecond\n";
    print "\n";
    print "Misc options: Protocol control\n";
    print "              dummy=true         #Backup files from previous run, update configuration file, print the command but do not execute it.\n";
    print "                                 # Use the dummy mode to prepare the run, but execute it manually, e.g., in a special mode on a supercomputer\n";
    print "              emonly=true        #Run only the em (will autodetect configuration for water or vacuum)\n";
    print "              nofree=true        #Run em, heating, equil and prepare free but do not execute it\n";
    print "              freeonly=true      #Run only the free step using previously prepared files\n";
    print "              amd=true           #After the freeMD step run accelerated MD\n";
    print "              amdprep=true       #Prepare configuration for accelerated MD based on the canonical freeMD step and exit\n";
    print "              amdonly=true       #Run only the accelerated MD step\n";
    print "              nostep2=true       #Do not run water optimization before energy minimization with no constraints\n";
    print "\n";
    print "Misc options: Amber binary control\n";
    print "              pmemdcuda=true     #Use GPU PMEMD implementation of sander [single-node for local use, MPI at the supercomputer] (default)\n";
    print "              pmemdmpi=true      #Use MPI PMEMD implementation of sander\n";
    print "              sandermpi=true     #Use the original implementation of SANDER with MPI\n";
    #print "              blockcudampi=true  #When running on a supercomputer use single-node build of the GPU PMEMD instead of MPI build (default=false)\n";
    print "\n";
    print "Misc options: Localhost MPI control\n";
    print "              uselocalmpi=true   #Launch Amber binary locally with mpirun\n";
    print "              cpu=n              #Launch on #n cpus (default = all local cores)\n";
    print "\n";
    print "Misc options: Supercomputer control\n";
    print "              slurm=true         #Launch on a supercomputer using Slurm manager\n";
    print "              script=<string>    #Slurm startup script - impi or ompi to run in parallel on multiple nodes; run to run on a single node\n";
    print "              clustercpu=n       #Launch on #n CPUs\n";
    print "              gpunodes=n         #Launch on #n GPU nodes\n";
    print "              queue=name         #Queue\n";
    print "              maxtime=minutes    #Set time limit for the task in minutes (FYC: 3h=179m, 6h=359m, 1d=1439m, 3d=4319m)\n";    
    print "\n";
    exit 0;
}

if (@ARGV==1 && $ARGV[0] eq "ex") {
    print "Examples:\n";
    print "#Run for the first time to create configuration files for the steps 1-6 with the default settings:\n";
    print ">autorun.pl my_file \n";
    print "\n";
    print "#Run for the first time to create configuration files for the steps 1-6 with the user-defined temperature (310 K)\n";
    print "#and length of equilibrations and production runs set to 10 ns and 100 ns, respectively:\n";
    print ">autorun.pl my_file temp=310 eqtime=10000 runtime=100000\n";
    print "\n";
    print "#Run for the second time to execute the protocol using the configuration files prepared during the previous run\n";
    print "#Use my_file.pdb, my_file.prmtop and my_file.incprd to run all steps locally with pmemdcuda:\n";
    print ">autorun.pl my_file \n";
    print "\n";
    print "#Same, but stops when the freeMD step reached (to be launched separatelly on a cluster)\n";
    print ">autorun.pl my_file nofree=true\n";
    print "\n";
    print "#Same, but runs only the freeMD step (does not check files from other steps -- for the cluster)\n";
    print ">autorun.pl my_file freeonly=true\n";
    print "\n";
    print "#Run all steps with pmemd.mpi on all local cores\n";
    print ">autorun.pl my_file pmemdmpi=true uselocalmpi=true\n";
    print "\n";
    print "#Run all steps with sander.mpi on 4 local cores\n";
    print ">autorun.pl my_file sandermpi=true uselocalmpi=true cpu=4\n";
    print "\n";
    print "#Run the free step on 28 CPU-cores at the Lomonosov supercomputer\n";
    print ">autorun.pl my_file freeonly=true slurm=true clustercpu=28 queue=compute maxtime=3000 script=ompi\n";
    print "\n";
    print "#Run the free step on 4 GPU-nodes (4 GPU cards) in the COMPUTE queue at the Lomonosov supercomputer\n";
    print ">autorun.pl my_file freeonly=true slurm=true gpunodes=4 queue=compute maxtime=3000 script=ompi\n";
    print "\n";
    print "#Run the free step on 1 GPU-node (2 GPU cards) in the PASCAL queue at the Lomonosov supercomputer\n";
    print ">autorun.pl my_file freeonly=true slurm=true gpunodes=1 queue=pascal maxtime=3000 script=ompi\n";
    print "\n";
    print "#Prepare the run (dummy mode) of the free step on 1 GPU-node (2 GPU cards) in the PASCAL queue at the Lomonosov supercomputer\n";
    print ">autorun.pl my_file dummy=true freeonly=true slurm=true gpunodes=1 queue=pascal maxtime=3000 script=ompi\n";
    print "#This will update the input/output/config files and print (but not execute) the run command\n";
    print "\n";
    print "#Run the free step on a single GPU card in the COMPUTE queue at the Lomonosov supercomputer\n";
    print ">autorun.pl my_file freeonly=true slurm=true gpunodes=1 queue=compute maxtime=3000 script=run\n";
    #print ">autorun.pl my_file freeonly=true lomonosov=true gpunodes=1 blockcudampi=true queue=compute maxtime=3000 script=run\n";
    print "\n";
    print "#Create configuration file for the aMD and execute the step on a supercomputer (given that the output of all previous six steps is available)\n";
    print ">autorun.pl my_file print=step7_amd\n";
    print ">autorun.pl my_file amdprep=true\n";
    print ">autorun.pl my_file amdonly=true slurm=true gpunodes=4 queue=compute maxtime=3000 script=ompi\n";
    print "\n";    
    exit 0;
}

#Setup allowed cluster queues and number of CPUs per node
my %cpu_per_node = ("pascal" => 12,
                    "compute" => 14,
                    "regular4" => 8,
                    "regular6" => 12, 
                    "hdd4" => 8,
                    "hdd6" => 12,
                    "gpu" => 8,
                    "smp" => 128,
                    "test" => 8,
                    "gputest" => 8);

#Setup cluster queues with MULTIPLE Graphical processing units
# DO NOT APPEND QUEUE WITH SINGLE-GPU NODES !!!
my %gpu_per_node = ("pascal" => 2);

my $INPUT = $ARGV[0];
my $EMONLY = 0;
my $NOFREE = 0;
my $FREEONLY = 0;
my $AMBERBIN_FLAG = 1;
#my $BLOCKGPUMPI = 0;
my $MPI_FLAG = 0;
my $CPU = 0;
my $AMD = 0;
my $AMDPREP = 0;
my $AMDONLY = 0;
my $NOSTEP2 = 0;
my $PRINT_STEP=0;

my $DUMMY=0;

my $LOMONOSOV=0;
my $SCRIPT="ompi";
my $QUEUE = "";
my $CLUSTER_CPU=0;
my $GPUNODES=0;
my $MAXTIME=0;

foreach my $i (1 .. $#ARGV) {
    if ($ARGV[$i] =~ /^print=(.*)/) {$PRINT_STEP = $1; $PRINT_STEP =~ s/\s+//g; next;}    
    if ($ARGV[$i] =~ /^dummy=true/) {$DUMMY=1; next;}    
    if ($ARGV[$i] =~ /^cutoff=(.*)/) {$CUTOFF = $1; $CUTOFF =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^temp=(.*)/) {$TEMPERATURE = $1; $TEMPERATURE =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^rest=(.*)/) {$RESTRAINTS_CONSTANT = $1; $RESTRAINTS_CONSTANT =~ s/\s+//g; next;}    
    if ($ARGV[$i] =~ /^eqdec=(.*)/) {$EQDEC = $1; $EQDEC =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^eqtime=(.*)/) {$EQTIME = $1; $EQTIME =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^runtime=(.*)/) {$RUNTIME = $1; $RUNTIME =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^pmemdcuda=true$/) {$AMBERBIN_FLAG = 1; next;}
    if ($ARGV[$i] =~ /^pmemdmpi=true$/) {$AMBERBIN_FLAG = 2; next;}
    if ($ARGV[$i] =~ /^sandermpi=true$/) {$AMBERBIN_FLAG = 3; next;}
#   if ($ARGV[$i] =~ /^blockcudampi=true$/) {$BLOCKGPUMPI = 1; next;}        
    if ($ARGV[$i] =~ /^emonly=true$/) {$EMONLY = 1; next;}
    if ($ARGV[$i] =~ /^nofree=true$/) {$NOFREE = 1; next;}
    if ($ARGV[$i] =~ /^freeonly=true$/) {$FREEONLY = 1; next;}
    if ($ARGV[$i] =~ /^uselocalmpi=true$/) {$MPI_FLAG = 1; next;}
    if ($ARGV[$i] =~ /^cpu=(\d+)$/) {$CPU = $1; next;}
    if ($ARGV[$i] =~ /^amd=true$/) {$AMD = 1; next;}
    if ($ARGV[$i] =~ /^amdprep=true$/) {$AMDPREP = 1; next;}
    if ($ARGV[$i] =~ /^amdonly=true$/) {$AMDONLY = 2; next;}
    if ($ARGV[$i] =~ /^nostep2=true$/) {$NOSTEP2 = 1; next;}    
    if ($ARGV[$i] =~ /^slurm=true$/) {$LOMONOSOV = 1; next;}
    if ($ARGV[$i] =~ /^lomonosov=true$/) {$LOMONOSOV = 1; next;}
    if ($ARGV[$i] =~ /^script=impi$/) {$SCRIPT = "impi"; next;}
    if ($ARGV[$i] =~ /^script=ompi$/) {$SCRIPT = "ompi"; next;}
    if ($ARGV[$i] =~ /^script=run$/) {$SCRIPT = "run"; next;}
    if ($ARGV[$i] =~ /^clustercpu=(\d+)/) {$CLUSTER_CPU = $1; next;}
    if ($ARGV[$i] =~ /^gpunodes=(\d+)/) {$GPUNODES = $1; next;}
    if ($ARGV[$i] =~ /^queue=(.*)/) {$QUEUE = $1; $QUEUE =~ s/\s+//g; next;}
    if ($ARGV[$i] =~ /^maxtime=(\d+)/) {$MAXTIME = $1; next;}    
    
    print "Error: Option $ARGV[$i] was not recognized\n";
    exit 1;
}

if ($NOFREE and $FREEONLY) {
    print "Error: Conflicting options nofree and freeonly can not be used at the same time\n";
    exit 1;
}

my $AMBER        = `echo \$AMBERHOME`; chomp($AMBER);
my $PMEMDCUDA    = "$AMBER/bin/pmemd.cuda";
my $PMEMDCUDAMPI = "$AMBER/bin/pmemd.cuda.MPI";
my $PMEMDMPI     = "$AMBER/bin/pmemd.MPI";
my $SANDERMPI    = "$AMBER/bin/sander.MPI";
my $MPIRUN       = `which mpirun`; chomp($MPIRUN);
my $AMBPDB       = "$AMBER/bin/ambpdb";

#if ($BLOCKGPUMPI) {
#    print "Warning: Forcing the use of the single-node build of GPU PMEMD instead of the MPI build\n";
#    print "\n";
#    $PMEMDCUDAMPI=$PMEMDCUDA;
#}

#-----------------------------------------------------------------------------------------------------
if (!$CPU) {
    $CPU = `cat /proc/cpuinfo  | grep processor | wc | awk '{print \$1}'`; chomp($CPU);
    print "Info: Number of CPUs in this system is $CPU\n";
} else {
    print "Info: Using $CPU CPUs to execute the task\n";
}
print "\n";

#-----------------------------------------------------------------------------------------------------

print "Input: Filemask $INPUT\n";
print "\n";

#-----------------------------------------------------------------------------------------------------

if ($LOMONOSOV) {
    print "Info: Running on Lomonosov-1/2 supercomputer\n";
    print "Info: Slurm startup script=$SCRIPT\n";
    print "Info: Queue=$QUEUE\n";
    
    if (!$cpu_per_node{$QUEUE}) {
        print "Error: Invalid queue name $QUEUE\n";
        exit 1;
    }
    
    if ($CLUSTER_CPU) {
        $GPUNODES = $CLUSTER_CPU / $cpu_per_node{$QUEUE};
        print "Info: Cluster CPUs=$CLUSTER_CPU\n";
    }
    else {
        if ($gpu_per_node{$QUEUE}) {
            my $nGPU = $gpu_per_node{$QUEUE};
            print "Info: GPUs-per-node=$nGPU\n";
            print "Info: Total number of GPUs to use=$nGPU*$GPUNODES\n";
        }
    }
    
    print "Info: Nodes=$GPUNODES\n";
    print "Info: Maxtime=$MAXTIME\n";    
    
    if ($GPUNODES == 1 && !$gpu_per_node{$QUEUE}) {
        print "Info: Using single-GPU build $PMEMDCUDA instead of MPI build to run on a single GPU\n";
        $PMEMDCUDAMPI=$PMEMDCUDA;
    }
    
    print "\n";    
}
# Define rules for local execution here
else {
    #Use single-unit build for local execution
    $PMEMDCUDAMPI=$PMEMDCUDA;
}

#-----------------------------------------------------------------------------------------------------
print "Info: Checking for the required input files\n";

my $input_prmtop = "$INPUT.prmtop";
my $input_inpcrd = "$INPUT.inpcrd";
my $input_pdb = "$INPUT.pdb";

my @check = ($input_prmtop, $input_inpcrd, $input_pdb);
foreach my $file (@check) {
    if (!-f $file || !-s $file) {print "Error: The required input file $file is not available or has zero size\n"; exit 1;}
    else {print "Info: Checking for file $file ... OK\n";}
}
print "\n";

#-----------------------------------------------------------------------------------------------------

if ($EMONLY) {
    print "Input: This script will run only the Energy Minimization protocol\n";
    
    print "Info: Determining if the solvent is present in $INPUT.pdb\n";
    my $water = calcWater(readpdb("$INPUT.pdb"));
    if ($water == 0) {
        print "Info: The molecular system does not contain detectable solvent and thus EM in vacuum will be performed\n";
        $EMONLY = 2;
        
        if ($AMBERBIN_FLAG != 3) {
            print "Error: Energy minimization in vacuum must be performed by sander.MPI\n";
            exit 1;
        }
        
    }
    else {
        print "Info: The molecular system contains $water solvent molecules and thus EM in water will be performed\n";
        $EMONLY = 3;
    }
    print "\n";
}

if ($NOFREE) {print "Input: This script will stop at the free MD step\n";}
if ($FREEONLY) {print "Input: This script will run only the free MD step\n";}

if ($NOSTEP2) {
    print "Input: This script will not run water optimization (step2) before executing energy minization with no constraints\n";
}

if ($AMD) {print "Input: This script will run the accelerated MD after the free MD step\n";}
if ($AMDPREP) {print "Input: This script will prepare configuration for the accelerated MD step\n";}
if ($AMDONLY) {print "Input: This script will run only the accelerated MD step\n";}

#-----------------------------------------------------------------------------------------------------
#Choose and check the amber binary
print "Info: Checking for Amber binaries\n";

if ($AMBER eq "" || ! -d $AMBER) {die "Error: Amber home folder not set or does not exist at \"$AMBER\"\n";} else {print "Info: Amber home folder $AMBER ... OK\n";}

if ($MPI_FLAG) {
    if (! -x $MPIRUN || ! -f $MPIRUN) {die "Error: The mpirun binary does not exist at \"$MPIRUN\" or is not executable\n";} else {print "Info: The mpirun binary $MPIRUN ... OK\n";}
}

if ($AMBERBIN_FLAG == 1) {
    if ($LOMONOSOV) {
        if (! -x $PMEMDCUDAMPI || ! -f $PMEMDCUDAMPI) {die "Error: The pmemd.cuda.MPI binary does not exist at \"$PMEMDCUDAMPI\" or is not executable\n";} else {print "Info: The pmemd.cuda.MPI binary $PMEMDCUDAMPI ... OK\n";}	
    }
    else {
        if (! -x $PMEMDCUDA || ! -f $PMEMDCUDA) {die "Error: The pmemd.cuda binary does not exist at \"$PMEMDCUDA\" or is not executable\n";} else {print "Info: The pmemd.cuda binary $PMEMDCUDA ... OK\n";}	
    }
}
elsif ($AMBERBIN_FLAG == 2) {
    if (! -x $PMEMDMPI || ! -f $PMEMDMPI) {die "Error: The pmemd.MPI binary does not exist at \"$PMEMDMPI\" or is not executable\n";} else {print "Info: The pmemd.MPI binary $PMEMDMPI ... OK\n";}
}
elsif ($AMBERBIN_FLAG == 3) {
    if (! -x $SANDERMPI || ! -f $SANDERMPI) {die "Error: The sander.MPI binary does not exist at \"$SANDERMPI\" or is not executable\n";} else {print "Info: The sander.MPI binary $SANDERMPI ... OK\n";}
}
else {
    print "Error: Flag AMBERBIN_FLAG=$AMBERBIN_FLAG is not supported\n";
    exit 1;
}

if (! -x $AMBPDB || ! -f $AMBPDB) {die "Error: The ambpdb binary does not exist at \"$AMBPDB\" or is not executable\n";} else {print "Info: The ambpdb binary $AMBPDB ... OK\n";}

print "\n";

#-----------------------------------------------------------------------------------------------------

#Creating filemasks for all steps of the MD simulation
my $em1 = "step1_em1";
my $em1_vacuum = "step1_em1_vacuum";
my $water = "step2_water";
my $em2 = "step3_em2";
my $heat = "step4_heat";
my $equil = "step5_equil";
my @equil_substeps = ();
if ($RESTRAINTS_CONSTANT == 0) {
    push(@equil_substeps, $equil);
}
else {
    my $decrement = 0;
    my $current_restraint_constant = $RESTRAINTS_CONSTANT;
    while ($current_restraint_constant > 0) {
        push(@equil_substeps, "$equil\__$current_restraint_constant");
        $current_restraint_constant -= $EQDEC; #decrement the constant
        $EQSTEPS++;
    }
    push(@equil_substeps, "$equil\__0"); #Final step with zero restraints    
    $EQSTEPS++;
}
my $free = "step6_free";
my $amd = "step7_amd";

# ..............................................................
#Set the order of the steps
my @steps = ();

if ($EMONLY == 2) {@steps = ($em1_vacuum);}
elsif ($EMONLY == 3) {@steps = ($em1);}
else {
    if ($NOSTEP2) {
        @steps = ($em1, $em2, $heat, @equil_substeps, $free);    
    }
    else {
        @steps = ($em1, $water, $em2, $heat, @equil_substeps, $free);
    }
        
    if ($AMD || $AMDONLY || $AMDPREP) {push(@steps, $amd);}
}
# ..............................................................

#Set extention for the configuration files (-i)
my $config_ext = "conf";
#Set extention for the output files (-o)
my $out_ext = "out";
#Set extention for the output files (-inf)
my $mdinfo_ext = "mdinfo";
#Set extention for the final coordinates and velocities (-r)
my $rst_ext = "rst";
#Set extention for the md trajectories (-x)
my $mdcrd_ext = "nc"; #nc extention for binary NetCDF output
#Set extention for the last pdb of each step
my $pdb_ext = "last.pdb";
#List all extentions except for the config_file extention
my @extentions = ($out_ext, $mdinfo_ext, $rst_ext, $mdcrd_ext, $pdb_ext);

#-----------------------------------------------------------------------------------------------------

if ($TEMPERATURE > 300) {
    $DT_STEP = 0.001;
    print "Warning: The dt integration step has been automatically changed to $DT_STEP ns due to user-requested temperature increase\n";
}
#-----------------------------------------------------------------------------------------------------

if ($PRINT_STEP) {    
    print "Info: Printing configuration for step $PRINT_STEP using the parameters below and exiting\n";    
    print "Input: Electrostatic cutoff $CUTOFF A\n";
    print "Input: Target temperature $TEMPERATURE K\n";
    print "Input: The integration step $DT_STEP ns\n";    
    print "Input: Restraints constant $RESTRAINTS_CONSTANT kcal/mol\n";    
    print "Input: Length of the equilibration run ".sprintf("%.1f", $EQTIME/1000)." ns\n";
    print "Input: Equilibration (step5_equil) will be performed in $EQSTEPS substeps with a $EQDEC kcal/mol*A2 decrement\n";
    print "Input: Length of the production runs ".sprintf("%.1f",$RUNTIME/1000)." ns\n";
    
    printConfig($PRINT_STEP);
    exit 0;
}

print "Info: Checking for the Amber configuration files\n";

my $config_exist = 0;
foreach my $step (@steps) {
    
    if ($step ne $free && $FREEONLY) {next;}
    if ($step ne $amd && $AMDONLY) {next;} 
    
    my $file = "$INPUT.$step.$config_ext";
    if (-f $file && -s $file) {$config_exist++; print "Info: Checking for the configuration file $file (step $step) ... OK\n";}    
    else {print "Info: Checking for file $file (step $step) ... N/A\n";}            
}

if ($config_exist == 0) {
    if ($FREEONLY) {
        print "Error: Script has been launched with the freeonly option but the configuration file for the $free step is not available\n";
        exit 1;
    }
    if ($AMDONLY) {
        print "Error: Script has been launched with the aimdonly option but the configuration file for the $free step is not available\n";
        exit 1;
    }
    
    print "\n";
    print "Info: Configuration files will be printed to the current folder using the parameters below\n";
    print "Input: Electrostatic cutoff $CUTOFF A\n";
    print "Input: Target temperature $TEMPERATURE K\n";
    print "Input: The integration step $DT_STEP ns\n";
    print "Input: Restraints constant $RESTRAINTS_CONSTANT kcal/mol\n";
    print "Input: Length of the equilibration run ".sprintf("%.1f", $EQTIME/1000)." ns\n";
    print "Input: Equilibration (step5_equil) will be performed in $EQSTEPS substeps with a $EQDEC kcal/mol*A2 decrement\n";
    print "Input: Length of the production runs ".sprintf("%.1f",$RUNTIME/1000)." ns\n";
    print "\n";
    
    foreach my $step (@steps) {
        printConfig($step);    
    }
    
    print "Info: All configuration files have been written to the current folder\n";
    print "Info: You may edit the content of the files but you must not change their names\n";
    print "Info: Restart the script after you have finished editing the configuration of your MD\n";
    print "\n";
    exit 0;
}
elsif ($config_exist != @steps && !$AMDONLY && !$FREEONLY) {
    print "Error: Some configuration files are present and some aren`t - this is not allowed (all or nothing)\n";
    exit 1;
}
print "\n";

#-----------------------------------------------------------------------------------------------------

print "Info: All input files are ready\n";
print "Info: Starting the MD simulation\n";
print "\n";

#-----------------------------------------------------------------------------------------------------

foreach my $stepid (0..$#steps) {
    
    my $step = $steps[$stepid];
    
    if ($step ne $free && $FREEONLY) {next;}
    
    if ($step ne $amd && $AMDONLY) {next;}
    
    if ($step eq $free && $NOFREE) {
        print "Info: Reached the $free step and will now stop\n";
        print "\n";
        last;
    }
        
    #Read the first line of the configuration file as the step title
    my $title = "";
    open(FI, "$INPUT.$step.$config_ext") or die "Error: Failed to open configuration file $INPUT.$step.$config_ext for reading\n";
    while (<FI>) {
        chomp;
        $title = $_;
        last;
    }    
    close(FI);
    
    print "Info: $title\n";    
    
    preparePDB($step, 1);
    
    if (checkOutput($step)) {
        print "Info: All output files are already present in the current folder\n";
        print "Info: This step will be skipped\n";
    }    
    else {        
        bkpOutput($step);

        launch($step, $stepid);

        preparePDB($step, 2); #Mode 2 means the method will exit with error if the output pdb is already present
        
        #Check files
        if (!checkOutput($step)) {
            print "Error: The required output files were not created\n";
            print "Error: Check the log files for errors\n";
            print "\n";
            print "Important: If you get the 'SIGSEGV, segmentation fault' when simulating very large molecular systems\n";
            print "Important: you may need to increase the stack size\n";
            print "Important: Check the current stack size limit by running 'ulimit -s'\n";
            print "Important: Increase hard limit for the stack size in the /etc/security/limits.conf file\n";
            print "Important: Set the pam_limits.so module, i.e., add line 'session required pam_limits.so to the /etc/pam.d/login'\n";
            print "Important: Relogin to your account to activate new settings\n";
            print "\n";
            exit 1;
        }
                
    }
    print "\n";
    
    #Prepare the configuration file for the AMD
    if ($step eq $free && ($AMD || $AMDPREP)) {
        #Read the average energies from the output file of the previous step and update the config file of the current step
        updateAMD("$INPUT.$step.$out_ext", "$INPUT.$amd.$config_ext");
        
        if ($AMDPREP) {
            print "Info: AMD configuration has been prepared. Exiting.\n";
            exit 0;
        }
        
    }
}

print "Done!\n";
print "\n";

#=====================================================================================================

#The main method which executes the md
sub launch {
    
    my $step = shift;
    my $stepid = shift;
    
    if ($steps[$stepid] ne $step) {print "Error: Internal inconsistency of step`s priorities [@steps] [$step] [$stepid]\n"; exit 1;}        
    
    #-------------------------------------------------------------------------------
    #Create the command to launch amber
    my $command = "";
    
    #-------------------------------------------------------------------------------
    #Define GPUs if required
    if ($gpu_per_node{$QUEUE}) {        
        my @gpus = ();
        foreach my $i (0..$gpu_per_node{$QUEUE}-1) {push(@gpus, $i);}
        $command .= "export CUDA_VISIBLE_DEVICES=".join(",", @gpus)."; ";
    }
    
    #-------------------------------------------------------------------------------    
    #Choose if to use MPI    
    if ($MPI_FLAG == 1) {
        $command .= "$MPIRUN -np $CPU -hosts localhost"
    }    
    elsif ($LOMONOSOV) {
        $command .= "sbatch";
        if ($CLUSTER_CPU) {$command .= " -N $GPUNODES --ntasks-per-node=$cpu_per_node{$QUEUE} -p $QUEUE -t $MAXTIME $SCRIPT";}
        else {
            if ($gpu_per_node{$QUEUE}) {
                $command .= " -N $GPUNODES -n ".($GPUNODES*$gpu_per_node{$QUEUE})." -p $QUEUE -t $MAXTIME $SCRIPT";
            }
            else {
                $command .= " -N $GPUNODES -p $QUEUE -t $MAXTIME $SCRIPT";
            }
        }
    }
    
    #-------------------------------------------------------------------------------
    #Choose sander or pmemd or pmemd.cuda
    if ($AMBERBIN_FLAG == 1) {
	if ($LOMONOSOV) {
        if ($CLUSTER_CPU) {
            $command .= " $PMEMDMPI";
        }
        else {
            $command .= " $PMEMDCUDAMPI";
        }
    }
	else {$command .= " $PMEMDCUDA";}
    }
    elsif ($AMBERBIN_FLAG == 2) {$command .= " $PMEMDMPI";}
    elsif ($AMBERBIN_FLAG == 3) {$command .= " $SANDERMPI";}    
    
    #-------------------------------------------------------------------------------
    #Set input configuration file and prmtop
    
    $command .= " -i $INPUT.$step.$config_ext -p $input_prmtop";
    
    #-------------------------------------------------------------------------------
    #Choose and check for the inpcrd / rst file
    #For first step em1 use $input_inpcrd
    if ($stepid == 0) {
        if (-f $input_inpcrd && -s $input_inpcrd) {
            $command .= " -c $input_inpcrd"; # <<< COMMAND    
        }
        else {
            print "Error: File $input_inpcrd does not exist or has zero size\n";
            exit 1;
        }        
    }
    else {
        #Check if restart files for the current step are available
        my $last_rst_bkp = 0;
        opendir (DIR, "./") or die "Error: Failed to list files in the current folder\n";
        while (my $file = readdir(DIR)) {
            next if ($file =~ m/^\./);
            if ($file =~ /$INPUT.$steps[$stepid].$rst_ext\_bkp(\d+)/) {
                if ($last_rst_bkp < $1) {
                    $last_rst_bkp = $1;
                }                
            }                        
        }
        closedir(DIR);
        
        #For free-equil, freeMD and AMD runs - check the latest restart file
        if (($step eq "$equil\__0" || $step eq $free || $step eq $amd) && $last_rst_bkp) {
            print "Info: The most recent restart file for the current step $steps[$stepid] is $INPUT.$steps[$stepid].$rst_ext\_bkp$last_rst_bkp\n";
            
            #If the rst restart file is available for the current step but is invalid - terminate
            if (!-s "$INPUT.$steps[$stepid].$rst_ext\_bkp$last_rst_bkp") {
                print "Error: The available restart file is invalid (zero-size)\n";
                print "Error: This script will now terminate. Further execution requires curation.\n";
                exit 1;
            }
            
            #If the rst restart file is OK but the output file is invalid - terminate
            if (!-s "$INPUT.$steps[$stepid].$out_ext\_bkp$last_rst_bkp") {
                print "Error: The available restart file is invalid (zero-size)\n";
                print "Error: This script will now terminate. Further execution requires curation.\n";
                exit 1;
            }

            #Determine the last step            
            my $laststep = 0;
            my $cur_bkp = $last_rst_bkp;
            while ($cur_bkp != 0) {                
                if (-f "$INPUT.$steps[$stepid].$out_ext\_bkp$cur_bkp" && -f "$INPUT.$steps[$stepid].$mdcrd_ext\_bkp$cur_bkp"
                    && -s "$INPUT.$steps[$stepid].$out_ext\_bkp$cur_bkp" && -s "$INPUT.$steps[$stepid].$mdcrd_ext\_bkp$cur_bkp") {
                    my $laststep_cur = 0;
                    open(FI, "$INPUT.$steps[$stepid].$out_ext\_bkp$cur_bkp") or die "Error: Failed to open file $INPUT.$steps[$stepid].$out_ext\_bkp$cur_bkp for reading\n";
                    while (<FI>) {
                        if (/NSTEP\s+=(\s*?\d+)\s+/) {
                            $laststep_cur = $1;
                        }
                    }
                    close(FI);
                    print "Info: The last step described in $INPUT.$steps[$stepid].$out_ext\_bkp$cur_bkp is $laststep_cur\n";
                    $laststep += $laststep_cur;
                }
                else {print "Info: Files $INPUT.$steps[$stepid].$out_ext\_bkp$cur_bkp and/or $INPUT.$steps[$stepid].$mdcrd_ext\_bkp$cur_bkp are not available or emtpy\n";}
                
                $cur_bkp --;
            }                                        
            
            if (!$laststep) {
                print "Error: Failed to determine the last step in the crashed run\n";
                exit 1;
            }
            print "Info: The total number of steps described in backed-up output files is $laststep\n";
            
            #my $ig = -1;
            #my $lastout = "$INPUT.$steps[$stepid].$out_ext\_bkp$last_rst_bkp";
            #print "Info: Retrieving the random seed from the last non-null output file $lastout\n";
            #open(FI, $lastout) or die "Error: Failed to open file $lastout for reading\n";
            #my $parameter_section = 0;
            #my $langevin_parameter_section = 0;
            #
            #while (<FI>) {
            #    if (/CONTROL  DATA  FOR  THE  RUN/) {
            #        $parameter_section=1
            #    }
            #    
            #    if ($parameter_section && /Langevin dynamics temperature regulation:/) {
            #        $langevin_parameter_section = 1;
            #    }
            #                    
            #    if ($langevin_parameter_section && /ig      =(.*)/) {                
            #        if ($ig != -1) {
            #            print "Error: Multiple definitions for ig found in the same output file $lastout\n";
            #            exit 0;
            #        }                    
            #        $ig = $1;
            #        $ig =~ s/\s+//g;
            #    }                
            #}                        
            #close(FI);
            #if ($ig==-1) {
            #    print "Warning: Failed to retrieve the random seed and using $ig instead\n";
            #}
            #else {
            #    print "Info: Random seed was set to $ig in the previous non-null run (output file $lastout)\n";
            #} 
            
            #Modify the configuration file
            print "Info: Modifying the configuration file $INPUT.$step.$config_ext\n";
            my $set_nstlim=0;
            #my $set_ig=0;
            my $config_tmp = "";
            my $newnsteps = 0;
            open(FI, "$INPUT.$step.$config_ext") or die "Error: Failed to open configuration file $INPUT.$step.$config_ext for reading\n";
            while (<FI>) {
                chomp;
                my $line = $_;
                if (/total_steps=(\d+)/) {
                    my $totalsteps = $1;
                    $newnsteps = $totalsteps - $laststep;
                    print "Info: The total number of steps to calculate during the full run is $totalsteps\n";                    
                }
                if (/nstlim=(\d+)/) {
                    if ($newnsteps <= 0) {
                        print "Error: The remaining number of steps to calculate is a non-positive value ($newnsteps) or the keyword total_steps is missing from the configuration file\n";
                        exit 1;
                    }
                    else {
                        print "Info: The remaining number of steps to calculate is $newnsteps\n";
                    }
                    $line =~ s/nstlim=(\d+)/nstlim=$newnsteps/;
                    $set_nstlim = 1;
                }
                #if (/ig=([-]?\d+)/) {                    
                #    print "Info: Setting the random seed to $ig\n";
                #    $line =~ s/ig=([-]?\d+)/ig=$ig/;
                #    $set_ig = 1;
                #}
                
                $config_tmp .= $line . "\n";                
            }            
            close (FI);
            
            if (!$set_nstlim) {
                print "Error: Failed to set nstlim in the configuration file\n";
                exit 0;
            }
            
            #if (!$set_ig) {
            #    print "Error: Failed to set ig in the configuration file\n";
            #    exit 0;
            #}
            
            open(FO, ">$INPUT.$step.$config_ext") or die "Error: Failed to open configuration file $INPUT.$step.$config_ext for writing\n";
            print FO $config_tmp;
            close(FO);                        
            
            #If the rst restart and output files are available for the current step and are OK - continue the crashed run
            print "Info: The simulation for the current step $steps[$stepid] will be continued from $INPUT.$steps[$stepid].$rst_ext\_bkp$last_rst_bkp\n";
            $command .= " -c $INPUT.$steps[$stepid].$rst_ext\_bkp$last_rst_bkp"; # <<< COMMAND                        
        }
        #For non [free-equil, freeMD and AMD] runs or in case the restart rst file does not exist - use the output from the previous step as the input (the default)
        else {
            if (-f "$INPUT.$steps[$stepid-1].$rst_ext" && -s "$INPUT.$steps[$stepid-1].$rst_ext") {
                $command .= " -c $INPUT.$steps[$stepid-1].$rst_ext"; # <<< COMMAND
            }
            else {
                print "Error: File $INPUT.$steps[$stepid-1].$rst_ext does not exist or has zero size\n";
                exit 1;
            }
        }
    }
        
    #-------------------------------------------------------------------------------
    
    #Select output files
    $command .= " -o $INPUT.$step.$out_ext -r $INPUT.$step.$rst_ext -inf $INPUT.$step.$mdinfo_ext";
    
    #-------------------------------------------------------------------------------
    
    #Now check if we need the position restraints
    open(FI, "$INPUT.$step.$config_ext") or die "Error: Failed to open file $INPUT.$step.$config_ext for reading\n";
    my $addrestraints = 0;
    while (<FI>) {if (/restraintmask/) {$addrestraints = 1; last;}}
    close(FI);
    if ($addrestraints && $stepid == 0) {$command .= " -ref $input_inpcrd";}
    elsif ($addrestraints && $stepid != 0) {$command .= " -ref $INPUT.".$steps[$stepid-1].".$rst_ext";}    
    
    #-------------------------------------------------------------------------------
    
    #Now check for md trajectory flag
    if ($step ne $em1 && $step ne $em2) {
        $command .= " -x $INPUT.$step.$mdcrd_ext";
    }
    
    #-------------------------------------------------------------------------------
    #Remove "logfile if necessery"
    if (-f "logfile") {
        print "Info: Removing logfile from the previous iteration\n";
        unlink("logfile");
    }    
    
    #-------------------------------------------------------------------------------
    
    #Execute
    
    if ($DUMMY) {
        print "Info: Dummy mode is ON\n";
        print "Info: The input/output/config files have been updated\n";
        print "Info: To run the simulation execute the command manually:\n";
        print "COMMAND $command\n";
        print "Info: Exiting\n";
        exit 0;
    }
    else {
        print "Info: Forking [$command]\n";
        
        my $pid = fork();
        if (not defined $pid) {
            print "Error: Failed to fork a thread\n";
            exit 1;
        }
        #Execute in the forked thread
        if (not $pid) {
            system($command);
            exit;
        }
        
        #-------------------------------------------------------------------------------
        
        #Monitor activity in the main process    
        while (1) {            
            sleep(1);    
        
            my $res = waitpid($pid, WNOHANG); #Request status code from the child thread    
                
            my $print_info="";            
            
            if (($step eq $em1 || $step eq $em2) && -f "$INPUT.$step.$out_ext") {
                my $read = 0;
                my @data = ();
                open(FI, "$INPUT.$step.$out_ext") or die "Error: Failed to open file $INPUT.$step.$mdinfo_ext for reading\n";
                while (<FI>) {
                    chomp;                                            	
                    if (/NSTEP\s+ENERGY\s+RMS\s+GMAX\s+NAME\s+NUMBER/) {$read = 1; next;}
                    if ($read) {
                        $read = 0;
                        @data = split(/\s+/);                                        
                    }
                }
                close (FI);
                if (@data == 0) {
                    $print_info = "Waiting for the information update ...";            
                }
                else {
                    $print_info = "Current progress: Nstep=$data[1] | Energy=$data[2] | RMS=$data[3]";
                }
                print $print_info . "       \r";
            }
            elsif (($step ne $em1 && $step ne $em2) && -f "$INPUT.$step.$mdinfo_ext") {           
                open(FI, "$INPUT.$step.$mdinfo_ext") or die "Error: Failed to open file $INPUT.$step.$mdinfo_ext for reading\n";
                while (<FI>) {
                    chomp;                                
    
                    if (/Total steps :\s+(.*)\s+\| Completed :\s+(.*)\s+\| Remaining/) {
                        $print_info = "Current progress: Total steps=$1 | Completed=$2";
                    }
                    if (/Estimated time remaining:\s+(.*)/) {
                        $print_info .= " | Estimated time remaining=$1";
                    }
                }
                close(FI);
                if ($print_info eq "") {
                    $print_info = "Waiting for the information update ...";
                }
                
                print $print_info . "       \r";
            }
            else {
                print "Waiting for the output files to appear ... \r";            
            }            
            
            #The child has failed
            if ($res == -1) {
                print "Error: Some error occurred while forking leap". $? ."\n";
                exit();
            }
            #The child has completed
            if ($res) {
                if ($LOMONOSOV) {
                    print "\nInfo: Task has been submitted to Lomonosov-1/2\n";
                    print "Info: Exiting\n";
                    exit 0;
                }
                print "\nInfo: Calculations have been completed\n";
                last;
            }
            #The child is still running
            if (!$res) {next;} 
        }
        #-------------------------------------------------------------------------------
        #Done
    }    
}

#Read the average energies from the output file of the previous step and update the config file of the current step
sub updateAMD {
    my $out_prev = shift;
    my $conf_amd = shift;
    
    print "Info: Updating the configuration for the AMD step\n";
    
    my $etot = 0;
    my $edih = 0;
    
    print "Info: Reading the average total potential energy and average dihedral energy from the file $out_prev\n";
    my $read = 0;
    open(FI, $out_prev) or die "Error: Failed to open output file $out_prev for reading\n";
    while (<FI>) {
        chomp;
        if (/A V E R A G E S\s+O V E R\s+\d+\s+S T E P S/) {            
            $read = 1;
            next;
        }
        if (!$read) {
            next;
        }
        
        if (/Etot\s+=\s+(.*)\s+EKtot/) {
            if ($etot) {print "Error: Duplicate average value for Etot in file $out_prev\n"; exit 0;}
            $etot = $1;
            $etot =~ s/\s+//g;
        }
        if (/DIHED\s+=\s+(.*)$/) {
            if ($edih) {print "Error: Duplicate average value for Edihed in file $out_prev\n"; exit 0;}
            
            $edih = $1;
            $edih =~ s/\s+//g;
        }
        if (/------------------------------------------------------------------------------/){
            $read = 0;
            next;
        }                
    }
    
    print "Info: Reading the number of amino acids in the protein and the number of atoms in the molecular system from $INPUT.pdb\n";
    my $pdb = readpdb("$INPUT.pdb");
    my $atoms = scalar keys %{$pdb};
    
    my %residues = ();
    my $noWAT_noNa_noCl = 1;
    my $pdb2 = pdbReformat($pdb, $noWAT_noNa_noCl);
    foreach my $chain (keys %{$pdb2}) {
        foreach my $resname (keys %{$pdb2->{$chain}}) {
            foreach my $resid (keys %{$pdb2->{$chain}->{$resname}}) {                
                $residues{$chain."__".$resname."__".$resid} = 1;
            }
        }
    }
    my $resnum = scalar keys %residues;        
    
    print "Info: The collected values are: Etot=$etot, EDihed=$edih, Residues(non water and ions)=$resnum, Atoms(total)=$atoms\n";
    
    if (!$etot || !$edih || !$atoms || !$resnum) {
        print "Error: Failed to collect one or more of the required parameters\n";
        exit 1;
    }    
    
    my $EthreshP = $etot + (0.16 * $atoms);
    my $alphaP = 0.16 * $atoms;
    my $EthreshD = $edih + (4 * $resnum);
    my $alphaD = (4 * $resnum)/5;        
    
    print "Info: The estimated AMD parameters are: EthreshP=$EthreshP, alphaP=$alphaP, EthreshD=$EthreshD, alphaD=$alphaD\n";
    
    if (!$EthreshP || !$alphaP || !$EthreshD || !$alphaD) {
        print "Error: Failed to estimate one or more of the required parameters\n";
        exit 1;
    }
    #if ($EthreshP > 0 || $alphaP < 0 || $EthreshD  < 0 || $alphaD < 0) {
    if ($EthreshP > 0 || $alphaP < 0 || $alphaD < 0) {
        print "Error: One or more of the estimated parameters are on the wrong scale\n";
        exit 1;
    }
    
    print "Info: Modifying the configuration for the AMD step $conf_amd\n";
    my $newconf = "";
    
    my $mod_ethD = 0; my $mod_ethP = 0; my $mod_alphaD = 0; my $mod_alphaP = 0;
    
    open(FI, $conf_amd) or die "Error: Failed to open file $conf_amd for reading\n";    
    while (<FI>) {
        chomp;
        my $line = $_;
        if ($line =~ /ethreshd=XXXX/) {
            $line =~ s/ethreshd=XXXX/ethreshd=$EthreshD/; $mod_ethD = 1;
        }
        if ($line =~ /ethreshp=XXXX/) {
            $line =~ s/ethreshp=XXXX/ethreshp=$EthreshP/; $mod_ethP = 1;
        }
        if ($line =~ /alphad=XXXX/) {
            $line =~ s/alphad=XXXX/alphad=$alphaD/; $mod_alphaD = 1;
        }
        if ($line =~ /alphap=XXXX/) { 
            $line =~ s/alphap=XXXX/alphap=$alphaP/; $mod_alphaP = 1;
        }        
        $newconf .= $line . "\n";  
    }    
    close(FI);
    
    if (!$mod_ethD || !$mod_ethP || !$mod_alphaD || !$mod_alphaP) {
        print "Error: Failed to modify at least one of the AMD parameters in the configuration file (ethD $mod_ethD; ethP $mod_ethP; alphaD $mod_alphaD; alphaP $mod_alphaP;)\n";
        print "Error: Maybe the configuration file $conf_amd was modified before. Re-run with amdonly=true or rewrite the configuration using print=$amd\n";
        exit 1;
    }
    
    open(FO, ">$conf_amd") or die "Error: Failed to open file $conf_amd for writing\n";
    print FO $newconf;
    close(FO);
    
    if (!-s $conf_amd) {
        print "Error: Updating the configuration file $conf_amd has failed\n";
        exit 1;
    }    
}

#Check if output files for particular step are available
sub checkOutput {
    my $step = shift;
    
    print "Info: Checking for output files\n";
    
    my $allpresent = 1;
    
    #First check if files do exist and have non-zero size
    foreach my $ext (@extentions) {
        
        #Skip mdcrd for energy minimizations
        if ($step eq $em1 || $step eq $em2) {
            if ($ext eq $mdcrd_ext) {
                next;
            }            
        }        
        
        my $file = "$INPUT.$step.$ext";
        if (!-f $file || !-s $file) {$allpresent = 0; print "Info: File $file ... N/A\n";}
        else {
            #Check the content of the output file
            if ($ext eq $out_ext) {
                open (FI, $file) || die "Error: Failed to open file $file for reading\n";
                my $total_wall_time: = 0;
                while (<FI>) {
                    if (/Total wall time/ || /FINAL RESULTS/ || /Final Performance Info/) {
                        $total_wall_time=1;
                        last;
                    }                    
                }                
                close(FI);
                
                if (!$total_wall_time) {
                    $allpresent = 0;
                    print "Info: File $file ... INCOMPLETE\n";
                    next;
                }                
            }                        
            print "Info: File $file ... OK\n";
        }                
    }
    
    if ($step eq $amd) {
        if (-f "amd.log") {
            print "Info: File amd.log ... OK\n";
        }
        else {
            $allpresent = 0;
            print "Info: File amd.log ... N/A\n";
        }
    }
    
    return $allpresent;
}

#Attempt to create PDB from the results of the current step (will not terminate of failure to create the file, only on the existence of previous file)
sub preparePDB {
    my $step = shift;
    my $mode = shift;
    
    if (-f $input_prmtop && -s $input_prmtop && -f "$INPUT.$step.$rst_ext" && -s "$INPUT.$step.$rst_ext") {        
    
        print "Info: Checking for the final PDB file after step $step\n";    
        
        my $file = "$INPUT.$step.$pdb_ext";
        print "Info: File $file ... ";
    
        if (!-f $file || !-s $file) {
            system("$AMBPDB -p $input_prmtop -c $INPUT.$step.$rst_ext > $file 2> /dev/null");
            if (!-f $file || !-s $file) {
                print "FAILED creating the file\n";
            }
            else {
                print "created from $input_prmtop and $INPUT.$step.$rst_ext\n";
            }
        }
        else {
            print "file exists\n";
            
            if ($mode == 2) {
                print "Error: File $file can be from the previous (unfinished) run\n";
                print "Error: Remove manually and restart to update\n";
                exit 1;
            }                    
        }
        
    }
}


#Clean output files for particular step are available
sub bkpOutput {
    my $step = shift;
    
    print "Info: Creating backups for old files\n";
    
    #First check if files do exist and have non-zero size
    
    #Pass one - select a common i for all files to be backuped
    my $last_bkp_i = 0;
    opendir (DIR, "./") or die "Error: Failed to list files in the current folder\n";
    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);
        if ($file =~ /$INPUT\.$step.*\_bkp(\d+)/) {
            if ($last_bkp_i < $1) {
                $last_bkp_i = $1;
            }                
        }
        if ($step eq $amd && $file =~ /"amd.log\_bkp(\d+)"/) {
            if ($last_bkp_i < $1) {
                $last_bkp_i = $1;
            }                
        }        
    }
    closedir(DIR);
    
    #Pass two - copy if needed    
    my $i = $last_bkp_i+1;
    foreach my $ext (@extentions) {                           
        my $file_original = "$INPUT.$step.$ext";
        if (-f $file_original) {
            my $file_copy = "$file_original\_bkp$i";            
            move($file_original, $file_copy);
            print "Info: Creating a backup copy of the file $file_original as $file_copy\n";
        }                    
    }
    if ($step eq $amd && -f "amd.log") {
        move("amd.log", "amd.log\_bkp$i");
        print "Info: Creating a backup copy of the file amd.log as amd.log\_bkp$i\n";
    }
}

sub printToFile {
    my $file = shift;
    my $content = shift;
    open(FO, ">$file") or die "Error: Failed to open configuration file $file for writing\n";
    print FO $content;
    close(FO);
    
    if (! -f $file || ! -s $file) {
        print "Error: Failed to print configuration file $file\n";
        exit 1;
    }
    
}

#Read in a pdb file
# version 2. Support added for PDB with more than 100 000 atoms
sub readpdb {
    my $input = shift;
    
    my %pdb = ();
    open(FI, $input) or die "ERROR: Can not open file $input for reading\n\n";
    
    my $last_num = 0;
    
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
                if ($num < $last_num) {
                    $num = $last_num;
                }                
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
            
            $last_num = $num;                    
        }
    }
    close(FI);
    
    return \%pdb;
}

sub calcWater {
    my $pdb = shift;    
    my $water = 0;    
    my $resid_prev = -1;
    
    foreach my $num (sort {$a <=> $b} keys %{$pdb}) {        
        my $chain = $pdb->{$num}->{"chainid"};
        my $resname = $pdb->{$num}->{"res"};
        my $resid = $pdb->{$num}->{"resn"};
        my $atom = $pdb->{$num}->{"atom"};
        
        my $x = $pdb->{$num}->{"xcor"};
        my $y = $pdb->{$num}->{"ycor"};
        my $z = $pdb->{$num}->{"zcor"};
                
        my @data = ($x, $y, $z);
        
        if ($resid_prev != $resid && ($resname eq "WAT" || $resname eq "HOH" || $resname =~ /^TIP/)) {
            $water++;
        }
        $resid_prev = $resid;
    }
    
    return $water;
}

sub pdbReformat {
    my $pdb = shift;
    my $noWAT_noNa_noCl = shift;
    
    my %hash = ();
    foreach my $num (sort {$a <=> $b} keys %{$pdb}) {
        #print "$num\n";
        my $chain = $pdb->{$num}->{"chainid"};
        my $resname = $pdb->{$num}->{"res"};
        my $resid = $pdb->{$num}->{"resn"};
        my $atom = $pdb->{$num}->{"atom"};
        
        my $x = $pdb->{$num}->{"xcor"};
        my $y = $pdb->{$num}->{"ycor"};
        my $z = $pdb->{$num}->{"zcor"};
                
        my @data = ($x, $y, $z);

        #Filter out Water (WAT) and ions (Na+ and Cl-)        
        if ($noWAT_noNa_noCl) {                    
            if ($resname eq "Na+" || $resname eq "Cl-" || $resname eq "WAT" ) {
                next;
            }
        }
        
        if ($hash{$chain}{$resname}{$resid}{$atom}) {
            print "Warning: Overwriting duplicate coordinates for atom $chain:$resname:$resid:$atom}\n";
        }
        
        $hash{$chain}{$resname}{$resid}{$atom} = \@data;        
    }
    
    return \%hash;
}

sub printConfig {
    
    my $step = shift;
    
    my $EQTIME_PER_SUBSTEP = $EQTIME / $EQSTEPS;
                
    my $configuration = "";
    
    if ($step eq $em1) {            
            
        $configuration = "Step-1: Minimize only the water and ions, restraining the heavy atoms of protein and ligands at constant $RESTRAINTS_CONSTANT kcal/mol-A^2
&cntrl  
imin=1,       ! Flag to run minimization - Single point energy calculation
irest=0,      ! Flag to restart a simulation -  Do not restart the simulation; instead, run as a new simulation  
ntx=1,        ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates, but no velocities, will be read    
ntmin=1,      ! For ncyc cycles the steepest descent method is used then conjugate gradient is switched on
maxcyc=5000,  ! Maximum number of minimization cycles to allow
ncyc=2500,    ! The method of minimization will be switched from SD to CG after ncyc cycles
cut=$CUTOFF,  ! The nonbonded cutoff, in Angstroms
nsnb=20,      ! Frequency at which the non-bonded list is updated
ntb=1,        ! Impose periodic boundaries at constant volume
ntp=0,        ! Flag for constant pressure dynamics - No pressure scaling  
ntc=1,        ! Flag for SHAKE to perform bond length constraints - SHAKE is not performed
ntf=1,        ! Force evaluation - complete interaction is calculated
igb=0,        ! Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models. Default is 0.
iwrap=1,      ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,       ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,     ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory  
ntpr=10,      ! Print the progress of the minimization to output file every ntpr steps
nmropt = 0,   ! No nmr-type analysis will be done.  
ntr=1,        ! Restrain atoms using a harmonic potential
restraintmask='!(:WAT,Na+,Cl-) & !\@H=', ! Apply restraints to all non-hydrogen atoms of all residues except water and ions
restraint_wt=$RESTRAINTS_CONSTANT,        ! Force constant (kcal/(mol Angstroms^2))
&end
/
";
    } elsif ($step eq $em1_vacuum) {            
            
        $configuration = "Step-1: Minimize in vacuum only the water and ions, restraining the heavy atoms of protein and ligands at constant $RESTRAINTS_CONSTANT kcal/mol-A^2
&cntrl  
imin=1,       ! Flag to run minimization - Single point energy calculation
irest=0,      ! Flag to restart a simulation -  Do not restart the simulation; instead, run as a new simulation  
ntx=1,        ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates, but no velocities, will be read    
ntmin=1,      ! For ncyc cycles the steepest descent method is used then conjugate gradient is switched on
maxcyc=5000,  ! Maximum number of minimization cycles to allow
ncyc=2500,    ! The method of minimization will be switched from SD to CG after ncyc cycles
cut=$CUTOFF,  ! The nonbonded cutoff, in Angstroms
nsnb=20,      ! Frequency at which the non-bonded list is updated
ntb=0,        ! Do not impose periodic boundaries at constant volume
ntp=0,        ! Flag for constant pressure dynamics - No pressure scaling  
ntc=1,        ! Flag for SHAKE to perform bond length constraints - SHAKE is not performed
ntf=1,        ! Force evaluation - complete interaction is calculated
igb=6,        ! No continuum solvent model is used  
ntxo=2,       ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,     ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory  
ntpr=10,      ! Print the progress of the minimization to output file every ntpr steps
nmropt = 0,   ! No nmr-type analysis will be done.  
ntr=1,        ! Restrain atoms using a harmonic potential
restraintmask='!(:WAT,Na+,Cl-) & !\@H=', ! Apply restraints to all non-hydrogen atoms of all residues except water and ions
restraint_wt=$RESTRAINTS_CONSTANT,        ! Force constant (kcal/(mol Angstroms^2))
&end
/
";
    }elsif ($step eq $water) {
    
    my $nstwater = 20/$DT_STEP; #i.e., 20ps
    
        $configuration = "Step-2: Let water and ions move (NTP, $TEMPERATURE K), restraining all atoms of protein and ligand at 10 kcal/mol-A^2
&cntrl  
imin=0,         ! Flag to run minimization - No minimization  
irest=0,        ! Flag to restart a simulation -  Do not restart the simulation; instead, run as a new simulation  
ntx=1,          ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates, but no velocities, will be read  
nstlim=$nstwater,   ! Number of MD-steps to be performed
t=0.0,          ! The time at the start (psec) this is for your own reference and is not critical.
dt=$DT_STEP,       ! The time step (psec).
ntc=2,          ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained  
ntf=2,          ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,    ! The nonbonded cutoff, in Angstroms
nsnb=20,        ! Frequency at which the non-bonded list is updated
ntb=2,          ! Impose periodic boundaries at constant pressure 
ntt=1,          ! Switch for temperature scaling - Constant temperature, using the weak-coupling algorithm
temp0=300.0,    ! Reference temperature at which the system is to be kept, if ntt > 0.
tempi=200.0,    ! Initial temperature
tautp=0.5,      ! Time constant, in ps, for heat bath coupling for the system, if ntt = 1 (smaller constant gives tighter coupling) 
ntp=1 ,         ! Flag for constant pressure dynamics - md with isotropic position scaling
taup=1.0,       ! Pressure relaxation time (in ps), when NTP > 0. The recommended is [1 - 5]. (use larger values to disturb Newton's equations as little as possible)
nscm=2500,      ! Flag for the removal of translational and rotational center-of-mass
iwrap=1,        ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ioutfm=1,       ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory  
ntxo=2,         ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ntpr=1000,      ! Print the progress of the minimization to output file every ntpr steps
ntwx=10000,     ! Every ntwx steps, the coordinates will be written to the mdcrd file
ntwr=10000,     ! Every ntwr steps during dynamics, the “restrt” file will be written
nmropt = 0,     ! No nmr-type analysis will be done.
ntr=1,          ! Flag for restraining specified atoms in Cartesian space using a harmonic potential
restraintmask='!(:WAT,Na+,Cl-)' , ! String that specifies the restrained atoms when ntr=1.
restraint_wt=10.0,                ! The weight (in kcal/mol−Å2) for the positional restraints.
&end
/
";
    } elsif ($step eq $em2) {
            
        $configuration = "Step-3: Minimize water and protein
&cntrl  
imin=1,       ! Flag to run minimization - Single point energy calculation
irest=0,      ! Flag to restart a simulation -  Do not restart the simulation; instead, run as a new simulation  
ntx=1,        ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates, but no velocities, will be read    
ntmin=1,      ! For ncyc cycles the steepest descent method is used then conjugate gradient is switched on
maxcyc=10000, ! Maximum number of minimization cycles to allow
ncyc=5000,    ! The method of minimization will be switched from SD to CG after ncyc cycles
cut=$CUTOFF,  ! The nonbonded cutoff, in Angstroms  
nsnb=20,      ! Frequency at which the non-bonded list is updated  
ntb=1,        ! Impose periodic boundaries at constant volume
ntp=0,        ! Flag for constant pressure dynamics - No pressure scaling  
ntc=1,        ! Flag for SHAKE to perform bond length constraints - SHAKE is not performed
ntf=1,        ! Force evaluation - complete interaction is calculated
igb=0,        ! Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models. Default is 0.  
iwrap=1,      ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,       ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,     ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=10,      ! Print the progress of the minimization to output file every ntpr steps
nmropt = 0,   ! No nmr-type analysis will be done.
&end
/
";
    } elsif ($step eq $heat) {
        my $heat_steps=100/$DT_STEP; #i.e., 100ps
        
        $configuration = "Step-4: Heating from 0 to $TEMPERATURE K (NVT) restraining heavy atoms of the protein and ligands at constant $RESTRAINTS_CONSTANT kcal/mol-A^2
&cntrl     
imin=0,             ! Flag to run minimization - No minimization
irest=0,            ! Flag to restart a simulation -  Do not restart the simulation; instead, run as a new simulation
ntx=1,              ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates, but no velocities, will be read
nstlim=$heat_steps, ! Number of MD-steps to be performed
dt=$DT_STEP,        ! The time step (psec). For temperatures above 300K, the step size should be reduced.
ntc=2,              ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained
ntf=2,              ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,        ! The nonbonded cutoff, in Angstroms
nsnb=20,            ! Frequency at which the non-bonded list is updated    
ntb=1,              ! Impose periodic boundaries at constant volume
ntp=0,              ! Flag for constant pressure dynamics - No pressure scaling
!ntt=1,             ! Switch for temperature scaling - Constant temperature, using the weak-coupling algorithm
!tautp=0.5,         ! Time constant, in ps, for heat bath coupling for the system, if ntt = 1 (smaller constant gives tighter coupling)
ntt=3,              ! Langevin thermostat
gamma_ln=2.0,       ! Collision frequency in ps−1 for Langevin thermostat
ig=-1,              ! Random seed for Langevin thermostat will be based on the current date and time
tempi=0.0,          ! Initial temperature
temp0=$TEMPERATURE, ! Reference temperature at which the system is to be kept, if ntt > 0.
igb=0,              ! Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models. Default is 0.    
nscm=500,           ! Flag for the removal of translational and rotational center-of-mass  
iwrap=1,            ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,             ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,           ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=10000,         ! Print the progress of the run to output file every ntpr steps
ntwr=50000,         ! Every ntwr steps during dynamics, the “restrt” file will be written
ntwx=50000,         ! Every ntwx steps, the coordinates will be written to the mdcrd file
ntr=1,              ! Restrain atoms using a harmonic potential
restraintmask='!(:WAT,Na+,Cl-) & !\@H=', ! Apply restraints to all non-hydrogen atoms of all residues except water and ions
restraint_wt=$RESTRAINTS_CONSTANT,       ! Force constant (kcal/(mol Angstroms^2))
nmropt=1            ! NMR restraints will be read (See TEMP0 control below)
&end

&wt
TYPE='TEMP0',          ! Varies the target temperature TEMP0
istep1=0,              ! Initial step
istep2=$heat_steps,    ! Final step
value1=0,              ! Initial temp0 (K)
value2=$TEMPERATURE, / ! Final temp0 (K)
&end  
&wt TYPE='END' &end     ! End of varying conditions
/  
";
    } elsif ($step =~ /^$equil\__(.*)/) {    
        
        my $restraint_constant = $1;
        $restraint_constant =~ s/\s+//g;
        
        if ($restraint_constant==0) {
            
            my $nsteps = sprintf("%.0f", $EQTIME_PER_SUBSTEP/$DT_STEP);
            
            $configuration = "Step-5: Equilibration at the target temperature ($TEMPERATURE K) in the NPT ensemble with no restraints
total_steps=$nsteps
&cntrl    
imin=0,            ! Flag to run minimization - No minimization
irest=1,           ! Flag to restart a simulation - Restart the simulation, reading coordinates and velocities from a previously saved restart file
ntx=5,             ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates and velocities will be read
nstlim=$nsteps,    ! Number of MD-steps to be performed
dt=$DT_STEP,       ! The time step (psec). For temperatures above 300K, the step size should be reduced.
ntc=2,             ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained
ntf=2,             ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,       ! The nonbonded cutoff, in Angstroms
ntb=2,             ! Impose periodic boundaries at constant pressure 
ntp=1,             ! Flag for constant pressure dynamics - md with isotropic position scaling
taup=2.0,          ! Pressure relaxation time (in ps), when NTP > 0. The recommended is [1 - 5]. (use larger values to disturb Newton's equations as little as possible)
ntt=3,             ! Langevin thermostat
gamma_ln=2.0,      ! Collision frequency in ps−1 for Langevin thermostat
ig=-1,             ! Random seed for Langevin thermostat will be based on the current date and time
temp0=$TEMPERATURE,! Reference temperature at which the system is to be kept, if ntt > 0.  
iwrap=1,           ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,            ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,          ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=100000,       ! Print the progress of the run to output file every ntpr steps
ntwr=100000,       ! Every ntwr steps during dynamics, the “restrt” file will be written
ntwx=100000,       ! Every ntwx steps, the coordinates will be written to the mdcrd file  
&end
/
";    
        }        
        elsif ($restraint_constant>0 && $restraint_constant <1) {

        my $nsteps = sprintf("%.0f", $EQTIME_PER_SUBSTEP/$DT_STEP);
        
        $configuration = "Step-5: Equilibration at the target temperature ($TEMPERATURE K) in the NPT ensemble with restraints at constant $restraint_constant kcal/mol-A^2
&cntrl    
imin=0,            ! Flag to run minimization - No minimization
irest=1,           ! Flag to restart a simulation - Restart the simulation, reading coordinates and velocities from a previously saved restart file
ntx=5,             ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates and velocities will be read
nstlim=$nsteps,    ! Number of MD-steps to be performed
dt=$DT_STEP,       ! The time step (psec). For temperatures above 300K, the step size should be reduced.
ntc=2,             ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained
ntf=2,             ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,       ! The nonbonded cutoff, in Angstroms
ntb=2,             ! Impose periodic boundaries at constant pressure
ntp=1,             ! Flag for constant pressure dynamics - md with isotropic position scaling
taup=1.0,          ! Pressure relaxation time (in ps), when NTP > 0. The recommended is [1 - 5]. (use larger values to disturb Newton's equations as little as possible)
ntt=3,             ! Langevin thermostat
gamma_ln=2.0,      ! Collision frequency in ps−1 for Langevin thermostat
ig=-1,             ! Random seed for Langevin thermostat will be based on the current date and time
temp0=$TEMPERATURE,! Reference temperature at which the system is to be kept, if ntt > 0.  
iwrap=1,           ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,            ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,          ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=10000,        ! Print the progress of the run to output file every ntpr steps
ntwr=100000,       ! Every ntwr steps during dynamics, the “restrt” file will be written
ntwx=100000,       ! Every ntwx steps, the coordinates will be written to the mdcrd file  
ntr=1,             ! Restrain atoms using a harmonic potential
restraintmask='!(:WAT,Na+,Cl-) & !\@H=', ! Apply restraints to all non-hydrogen atoms of all residues except water and ions
restraint_wt=$restraint_constant,        ! Force constant (kcal/(mol Angstroms^2))
&end
/
";                    
	}
        elsif ($restraint_constant>=1) {

        my $nsteps = sprintf("%.0f", $EQTIME_PER_SUBSTEP/$DT_STEP);
        
        $configuration = "Step-5: Equilibration at the target temperature ($TEMPERATURE K) in the NPT ensemble with restraints at constant $restraint_constant kcal/mol-A^2
&cntrl    
imin=0,            ! Flag to run minimization - No minimization
irest=1,           ! Flag to restart a simulation - Restart the simulation, reading coordinates and velocities from a previously saved restart file
ntx=5,             ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates and velocities will be read
nstlim=$nsteps,    ! Number of MD-steps to be performed
dt=$DT_STEP,       ! The time step (psec). For temperatures above 300K, the step size should be reduced.
ntc=2,             ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained
ntf=2,             ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,       ! The nonbonded cutoff, in Angstroms
ntb=2,             ! Impose periodic boundaries at constant pressure
ntp=1,             ! Flag for constant pressure dynamics - md with isotropic position scaling
taup=1.0,          ! Pressure relaxation time (in ps), when NTP > 0. The recommended is [1 - 5]. (use larger values to disturb Newton's equations as little as possible)
ntt=3,             ! Langevin thermostat
gamma_ln=2.0,      ! Collision frequency in ps−1 for Langevin thermostat
ig=-1,             ! Random seed for Langevin thermostat will be based on the current date and time
temp0=$TEMPERATURE,! Reference temperature at which the system is to be kept, if ntt > 0.  
iwrap=1,           ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,            ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,          ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=10000,        ! Print the progress of the run to output file every ntpr steps
ntwr=100000,       ! Every ntwr steps during dynamics, the “restrt” file will be written
ntwx=100000,       ! Every ntwx steps, the coordinates will be written to the mdcrd file  
ntr=1,             ! Restrain atoms using a harmonic potential
restraintmask='!(:WAT,Na+,Cl-) & !\@H=', ! Apply restraints to all non-hydrogen atoms of all residues except water and ions
restraint_wt=$restraint_constant,        ! Force constant (kcal/(mol Angstroms^2))
&end
/
";                    
        }
        elsif ($restraint_constant<0) {print "Error: The restraints constant must be positive or zero ($restraint_constant)\n"; exit 1;}        
    } elsif ($step eq $free) {
        
        my $write = 100000;
        #my $nsteps = 500000000 * 0.002/$DT_STEP; #i.e., a microsecond MD trajectory by default
        my $nsteps = sprintf("%.0f", $RUNTIME/$DT_STEP);
        
        if ($nsteps > 999999999) {                                    
            $nsteps = 999999999;
            print "Warning: The length of the $step has changed and was set to the top-limit number of steps $nsteps (~".sprintf("%.2f", $nsteps*$DT_STEP/1000)." ns)\n\n";
        } 
                                                        #The nsteps must be 9-digit only. You can run several consequtive MDs and then
                                                        #merge them together to get a longer trajectory
        
        $configuration = "Step-6: Free MD simulation in the NVT ensemble at $TEMPERATURE K
total_steps=$nsteps
&cntrl    
imin=0,            ! Flag to run minimization - No minimization
irest=1,           ! Flag to restart a simulation - Restart the simulation, reading coordinates and velocities from a previously saved restart file
ntx=5,             ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates and velocities will be read
nstlim=$nsteps,    ! Number of MD-steps to be performed
dt=$DT_STEP,       ! The time step (psec). For temperatures above 300K, the step size should be reduced.
ntc=2,             ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained
ntf=2,             ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,       ! The nonbonded cutoff, in Angstroms
ntb=1,             ! Impose periodic boundaries at constant volume
ntp=0,             ! Flag for constant pressure dynamics - No pressure scaling
ntt=3,             ! Langevin thermostat
gamma_ln=2.0,      ! Collision frequency in ps−1 for Langevin thermostat
ig=-1,             ! Random seed for Langevin thermostat will be based on the current date and time
temp0=$TEMPERATURE ! Reference temperature at which the system is to be kept, if ntt > 0.  
iwrap=1,           ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,            ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,          ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=$write,       ! Print the progress of the run to output file every ntpr steps
ntwr=$write,       ! Every ntwr steps during dynamics, the “restrt” file will be written
ntwx=$write,       ! Every ntwx steps, the coordinates will be written to the mdcrd file
&end
/
";
    } elsif ($step eq $amd) {
        my $write = 100000;
        my $nsteps = sprintf("%.0f", $RUNTIME/$DT_STEP);
        
        if ($nsteps > 999999999) {                                    
            $nsteps = 999999999;
            print "Warning: The length of the $step has changed and was set to the top-limit number of steps $nsteps (~".sprintf("%.2f", $nsteps*$DT_STEP/1000)." ns)\n\n";
        }
        
        $configuration = "Step-7: Accelerated MD in the NVT ensemble at $TEMPERATURE K
total_steps=$nsteps
&cntrl    
imin=0,            ! Flag to run minimization - No minimization
irest=1,           ! Flag to restart a simulation - Restart the simulation, reading coordinates and velocities from a previously saved restart file
ntx=5,             ! Option to read the initial coordinates, velocities and box size from the inpcrd file - Coordinates and velocities will be read
nstlim=$nsteps,    ! Number of MD-steps to be performed
dt=$DT_STEP,       ! The time step (psec). For temperatures above 300K, the step size should be reduced.
ntc=2,             ! Flag for SHAKE to perform bond length constraints - bonds involving hydrogen are constrained
ntf=2,             ! Force evaluation - bond interactions involving H-atoms omitted (use with ntc=2)
cut=$CUTOFF,       ! The nonbonded cutoff, in Angstroms
ntb=1,             ! Impose periodic boundaries at constant volume
ntp=0,             ! Flag for constant pressure dynamics - No pressure scaling
ntt=3,             ! Langevin thermostat
gamma_ln=2.0,      ! Collision frequency in ps−1 for Langevin thermostat
ig=-1,             ! Random seed for Langevin thermostat will be based on the current date and time
temp0=$TEMPERATURE,! Reference temperature at which the system is to be kept, if ntt > 0.
iwrap=1,           ! The coordinates written to the restart and trajectory files will be wrapped into a primary box
ntxo=2,            ! Format of the final coordinates, velocities, and box size written to the restart file - NetCDF file
ioutfm=1,          ! The format of coordinate and velocity trajectory files - Binary NetCDF trajectory
ntpr=$write,       ! Print the progress of the run to output file every ntpr steps
ntwr=$write,       ! Every ntwr steps during dynamics, the “restrt” file will be written
ntwx=$write,       ! Every ntwx steps, the coordinates will be written to the mdcrd file
iamd=3,            ! Accelerated MD - boost the whole potential with an extra boost to the torsions
ethreshd=XXXX,
alphad=XXXX,
ethreshp=XXXX,
alphap=XXXX,
&end
/
";
    }
    else {
        print "Error: Unknown step $step\n";
        exit 1;
    }
    printToFile("$INPUT.$step.$config_ext", $configuration);
}
