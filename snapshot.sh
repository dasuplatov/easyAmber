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

if [ ! -n "$4" ] && [ ! -n "$6" ];then
    echo ""
    echo "This script will print snapshots from Amber NetCDF trajectory"
    echo "v. 2.0 2019-12-15"
    echo
    echo "NB!: Frame numbering usually should be +1 to what VMD tells you"
    echo 
    echo "Write PDB and RST for a selected frame"
    echo "Usage: $0 <step> <input: solvent/nosolvent> <output: solvent/nosolvent> <frame num>"
    echo 
    echo "Write PDB batch for within a selected range"
    echo "Usage: $0 <step> <input: solvent/nosolvent> <output: solvent/nosolvent> <start> <stop> <increment>"
    echo ""
    echo "                 <input: solvent/nosolvent> - activates one of the two filename/filetype scenarios (parm7+netcdf or pdb+netcdf)"
    echo "                 <output: solvent/nosolvent> - when input=solvent you have an option of stripping water from the output"
    echo 
    exit
fi

CPPTRAJ=`which cpptraj`
if [ ! -x "$CPPTRAJ" ];then
    echo "Error: Failed to execute cpptraj"
    exit
fi

STEP=$1

INPUT_SOLVENT=""
if [ "$2" == "solvent" ];then
    INPUT_SOLVENT="1"
elif [ "$2" == "nosolvent" ];then
    INPUT_SOLVENT="0"
else 
    echo "Error: Input Solvent option \"$2\" was not recognized"
    exit
fi

OUTPUT_SOLVENT=""
if [ "$3" == "solvent" ];then
    OUTPUT_SOLVENT="1"
elif [ "$3" == "nosolvent" ];then
    OUTPUT_SOLVENT="0"
else 
    echo "Error: Output Solvent option \"$3\" was not recognized"
    exit
fi

if [ "$INPUT_SOLVENT" == "0" ] && [ "$OUTPUT_SOLVENT" == "1" ]; then
    echo "Error: You have indicated input_nosolvent and output_solvent which is obviously not possible"
    exit 1
fi

START=0
STOP=0
INCREMENT=0

if [ -n "$6" ];then
    START=$4
    STOP=$5
    INCREMENT=$6
else
    START=$4
    STOP=$4
    INCREMENT=1
fi

STRUCTURE=""
TRJ=""
MASK=""
MODEL=""

if [ "$INPUT_SOLVENT" == "1" ];then
    STRUCTURE=`ls | grep ".*.prmtop$"`
    MASK=`echo $STRUCTURE | sed s/.prmtop//` 
    TRJ="$MASK.$STEP.nc"
    MODEL="$MASK.pdb"
elif [ "$INPUT_SOLVENT" == "0" ];then
    STRUCTURE=`ls | grep ".*.nowater.pdb$"`
    MASK=`echo $STRUCTURE | sed s/.nowater.pdb//` 
    TRJ="$MASK.$STEP.nowater.nc"
    MODEL="$STRUCTURE"
else 
    echo "Error: Unknown value of the input_solvent flag ($INPUT_SOLVENT)"
    exit
fi

if [ -f "$STRUCTURE" ];then
    echo "Input: Amber structure parameters file: $STRUCTURE"
else
    echo "Error: Amber structure parameters file not found"
    exit
fi

if [ -f "$MODEL" ];then
    echo "Input: Amber PDB model file: $MODEL"
else
    echo "Error: Amber PDB model file $MODEL not found"
    exit
fi

if [ -f "$TRJ" ];then
    echo "Input: NetCDF trajectory file $TRJ"
else
    echo "Error: Trajectory file $TRJ not found"
    exit
fi

OUTPUT="$MASK.$STEP.snapshots_"$START"_"$STOP"_"$INCREMENT".pdb"
OUTPUT_RST="$MASK.$STEP.snapshots_"$START".rst"

if [ -f "$OUTPUT" ];then
    echo "Error: The pdb output file $OUTPUT already exists"
    exit
fi

echo "Input: pdb/prmtop  $STRUCTURE"
echo "Input: nc          $TRJ"
echo "Input: step        $STEP"
echo "Input: start       $START"
echo "Input: stop        $STOP"
echo "Input: increment   $INCREMENT"
echo "Input: out solvent $OUTPUT_SOLVENT"
echo "Input: output      $OUTPUT"
if [ "$START" == "$STOP" ];then
    echo "Input: output rst  $OUTPUT_RST"
fi
echo

cpptrajfile="cpptraj_$STEP"_"snapshots_"$START"_"$STOP"_"$INCREMENT".txt"
echo "trajin $TRJ $START $STOP $INCREMENT
reference $MODEL parm $STRUCTURE [ref]
rms align ref [ref] @CA @CA
" > $cpptrajfile

if [ "$OUTPUT_SOLVENT" == "0" ];then
    echo "strip :WAT,Na+,Cl-" >> $cpptrajfile
fi

echo "trajout $OUTPUT pdb" >> $cpptrajfile

if [ "$START" == "$STOP" ];then
    echo "trajout $OUTPUT_RST restart onlyframes" >> $cpptrajfile
fi

$CPPTRAJ -i $cpptrajfile -p $STRUCTURE

echo "Info: The output pdb has been written to $OUTPUT"
if [ "$START" == "$STOP" ];then
    echo "Info: The output rst has been written to $OUTPUT_RST"
fi
echo "Done!"
echo ""
