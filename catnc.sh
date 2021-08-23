#!/bin/bash

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

echo ""
echo "Concatenate NetCDF/Out/aMDlog Amber trajectory files from crushed/re-started runs"
echo "v. 2.0 2019-12-15"
echo ""

if [ ! -n "$1" ]; then
    echo
    echo "Usage: $0 <step> "
    echo "Usage: $0 <step> force #force concatenation of the incomplete run"
    echo
    echo "Example: $0 step6_free force"
    echo "Example: $0 step7_amd force"
    echo
    exit
fi

FORCE=0
if [ -n "$2" ] && [ "$2" == "force" ]; then
    FORCE=1
fi

#METHODS

#Compare aMD log with the out file
compare()
{
    VAR1=`cat $AMDLOG | grep -v '#' | awk '{print $2}'`
    VAR2=`cat $OUT | grep NSTEP | awk '{var=substr($0, 9, 9); gsub(/ /, "", var); print var;}'`
    if [ "$VAR1" != "$VAR2" ];then
        echo "Error: The sequence of steps in files $OUT and $AMDLOG is different"
        exit
    fi
}


#MAIN

CPPTRAJ=`which cpptraj`
if [ ! -x "$CPPTRAJ" ];then
    echo "Error: Failed to execute cpptraj"
    exit
fi

STEP=$1

BACKUP="files_bkp"
if [ -d "$BACKUP" ];then
    echo "Error: Backup folder $BACKUP already exists"
    exit
fi

echo "Input: Using stepid $STEP"
echo

PRMTOP=`ls | grep \.prmtop$`
if [ ! -f "$PRMTOP" ];then
    echo "Error: Failed to select prmtop file"
    exit
fi
echo "Info: Using prmtop file $PRMTOP"

MASK=`echo $PRMTOP | sed s/\.prmtop$//`
echo "Info: Using filemask for other files $MASK"

LASTOUT="$MASK.$STEP.out"
if [ ! -f "$LASTOUT" ];then
    echo "Error: Failed to select the last OUT file with name $LASTOUT"
    exit
fi
echo "Info: Last OUT file is $LASTOUT"

if [ "`cat $LASTOUT | grep 'Total wall time'`" == "" ] ;then
    echo "Error: The last OUT file $LASTOUT is incomplete"
    
    if [ "$FORCE" == "0" ];then
        exit
    fi  
    
    echo "Info: User has forced override of this error"
fi

LASTNC="$MASK.$STEP.nc"
if [ ! -f "$LASTNC" ];then
    echo "Error: Failed to select the last NC file with name $LASTNC"
    exit
fi
echo "Info: Last NC file is $LASTNC"

#Now for aMD check the amd.log file
if [ "$STEP" == "step7_amd" ];then
    LASTAMDLOG="amd.log"
    if [ ! -f "$LASTAMDLOG" ];then
        echo "Error: aMD log file $LASTAMDLOG is not present"
        exit
    fi
    
    OUT=$LASTOUT
    AMDLOG=$LASTAMDLOG
    compare
    
    echo "Info: Last aMD log file is $LASTAMDLOG"
fi

echo "Info: Searching for part-files in NC and OUT format"
COMMAND_NC=""
COMMAND_OUT=""
COMMAND_AMD=""

#RST files are indicators of a non-null run (i.e., if rst file was created than at least some
#steps were calculated 
lastrst=`ls | grep "$MASK.$STEP.rst_bkp" | sed s/.*rst.bkp// | sort -V | tail -n 1`
if [ "$lastrst" != "" ];then
    echo "Info: Last available RST file is $MASK.$STEP.rst_bkp$lastrst"
else
    echo "Info: No RST backup files were found"
    echo "Info: Seems like there is no work to do in this folder"
    echo "Done!"
    echo
    exit
fi

i=0
while [ "$i" -lt "$lastrst" ]
do
    let i=i+1
    
    #First check the rst file
    PARTRST="$MASK.$STEP.rst_bkp$i"
    if [ ! -f "$PARTRST" ];then
        echo "Info: $PARTRST ... missing"
        continue
    fi
    echo "Info: $PARTRST ... available"
    
    #Next check the out file
    PARTOUT="$MASK.$STEP.out_bkp$i"
    if [ ! -f "$PARTOUT" ];then
        echo "Info: $PARTOUT ... missing"
        continue
    fi
    echo "Info: $PARTOUT ... available"

    #If we are here than the out file must not be empty
    if [ "`cat $PARTOUT | grep 'NSTEP =  '`" == "" ];then
        echo "Error: OUT file for part #$i $PARTOUT does not contain step information"
        exit
    fi

    #Finally, if we are here than the nc file must be present
    PARTNC="$MASK.$STEP.nc_bkp$i"
    if [ ! -f "$PARTNC" ];then
        echo "Error: NC file for part #$i is not present"
        exit 
    fi
    
    #Now for aMD check the amd.log file
    if [ "$STEP" == "step7_amd" ];then
        PARTAMDLOG="amd.log_bkp$i"
        if [ ! -f "$PARTAMDLOG" ];then
            echo "Error: aMD log file $PARTAMDLOG is not present"
            exit
        fi
        
        OUT=$PARTOUT
        AMDLOG=$PARTAMDLOG
        compare
        
        echo "Info: $PARTAMDLOG ... available"
    fi
    
    if [ "$STEP" == "step7_amd" ];then
        echo "Info: Part #$i - concatenating files $PARTNC, $PARTOUT and $PARTAMDLOG"
        COMMAND_AMD="$COMMAND_AMD $PARTAMDLOG"
    else
        echo "Info: Part #$i - concatenating files $PARTNC and $PARTOUT"
    fi
    
    COMMAND_NC="$COMMAND_NC -y $PARTNC"
    COMMAND_OUT="$COMMAND_OUT $PARTOUT"
done

TEMPNC="TEMP.nc"
if [ -f "$TEMPNC" ];then
    echo "Error: Temporary NC file $TEMPNC already exists"
    exit
fi

TEMPOUT="TEMP.out"
if [ -f "$TEMPOUT" ];then
    echo "Error: Temporary OUT file $TEMPOUT already exists"
    exit
fi

TEMPAMD="TEMP.amd"
if [ -f "$TEMPAMD" ];then
    echo "Error: Temporary AMDLOG file $TEMPAMD already exists"
    exit
fi

COMMAND_NC="$CPPTRAJ -p $PRMTOP $COMMAND_NC -y $LASTNC -x $TEMPNC"
COMMAND_OUT="cat $COMMAND_OUT $LASTOUT"
COMMAND_AMD="cat $COMMAND_AMD $LASTAMDLOG"

echo "Info: Concatenating trajectories"
echo "Info: Executing command [$COMMAND_NC]"
echo
echo "-------------------cpptraj log starts-------------------"
$COMMAND_NC
echo "-------------------cpptraj  log  ends-------------------"
echo

echo "Info: Concatenating OUT files"
echo "Info: Executing command [$COMMAND_OUT]"
$COMMAND_OUT > $TEMPOUT
echo 

if [ "$STEP" == "step7_amd" ];then
    echo "Info: Concatenating aMD log files"
    echo "Info: Executing command [$COMMAND_AMD] and postprocessing"
    $COMMAND_AMD | awk 'BEGIN {H=0; S=0}
    /^  #/ {if (H==0 || H==1) {H=1; print;}; next;}
    // {if (H==1 || H==0) {H=2}; S=S+$1; $2=S; printf("%12s%10s%22s%22s%22s%22s%22s%22s\n", $1, $2, $3, $4, $5, $6, $7, $8)}
    ' > $TEMPAMD
    echo     
fi

#Making a backup of old part-files
echo "Info: Backing up obsolete data to $BACKUP"
mkdir $BACKUP
j=1
while [ "$j" -le "$i" ]
do
    file="$MASK.$STEP.last.pdb_bkp$j"; if [ -f $file ]; then mv $file $BACKUP/; fi
    file="$MASK.$STEP.mdinfo_bkp$j"; if [ -f $file ]; then mv $file $BACKUP/; fi
    file="$MASK.$STEP.rst_bkp$j"; if [ -f $file ]; then mv $file $BACKUP/; fi
    mv "$MASK.$STEP.out_bkp$j" $BACKUP/
    mv "$MASK.$STEP.nc_bkp$j" $BACKUP/
    if [ "$STEP" == "step7_amd" ];then
	mv amd.log_bkp$j $BACKUP/
    fi
    let j=j+1
done
mv $LASTNC $BACKUP
mv $LASTOUT $BACKUP

if [ "$STEP" == "step7_amd" ];then
    mv $LASTAMDLOG $BACKUP
fi

echo "Info: Storing final NC and OUT files"
mv $TEMPNC $LASTNC
mv $TEMPOUT $LASTOUT

if [ "$STEP" == "step7_amd" ];then
    mv $TEMPAMD $LASTAMDLOG
fi

echo ""
echo "Info: Checking the final trajectory file $LASTNC"
echo -n "Info: "`$CPPTRAJ -p $PRMTOP -y $LASTNC -tl`
echo ""
echo "Done!"
echo ""
