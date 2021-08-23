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

if [ ! -n "$1" ];then
    echo ""
    echo "This script will strip water molecules from the input trajectory"
    echo "and write a new trajectory to a separate file"
    echo "v. 2.0 2019-12-15"
    echo
    echo "Usage: $0 <step>"
    echo ""
    echo "Example: $0 step6_free"
    echo "Example: $0 step7_amd"
    echo 
    exit
fi

CPPTRAJ=`which cpptraj`
if [ ! -x "$CPPTRAJ" ];then
    echo "Error: Failed to execute cpptraj"
    exit
fi

STEP=$1

PRMTOP=`ls | grep '\.prmtop$'`

if [ ! -f "$PRMTOP" ];then
    echo "Error: The prmtop file $PRMTOP does not exist"
    exit
fi

MASK=`echo $PRMTOP | sed s/.prmtop$//`

PDB="$MASK.pdb"

if [ ! -f "$PDB" ];then
    echo "Error: The pdb file $PDB does not exist"
    exit
fi

PDBOUT="$MASK.nowater.pdb"

if [ -f "$PDBOUT" ];then
    echo "Error: The output pdb file $PDBOUT already exists"
    exit
fi

NC="$MASK.$STEP.nc"

if [ ! -f "$NC" ];then
    echo "Error: The nc trajectory file $NC does not exist"
    exit
fi

NCOUT="$MASK.$STEP.nowater.nc"

if [ -f "$NCOUT" ];then
    echo "Error: The output nc trajectory file $NCOUT already exists"
    exit
fi

echo "Input: prmtop     $PRMTOP"
echo "Input: nc         $NC"
echo "Input: step       $STEP"
echo "Input: nc output  $NCOUT"
echo "Input: pdb output $PDBOUT"
echo

#strip :WAT,Na+,Cl-

#To autoimage a complex molecular system select anchor manually
#autoimage [<mask> | anchor <mask> [fixed <mask>] [mobile <mask>]]
#autoimage :1-1000
#More at https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml
echo "trajin $NC
trajout $NCOUT
strip :WAT
center
autoimage
" > "cpptraj_"$STEP"_nc_nowater.txt"

$CPPTRAJ -i  "cpptraj_"$STEP"_nc_nowater.txt" -p $PRMTOP

#strip :WAT,Na+,Cl-
echo "trajin $PDB
trajout $PDBOUT
strip :WAT
center
autoimage
" > "cpptraj_"$STEP"_pdb_nowater.txt"

$CPPTRAJ -i "cpptraj_"$STEP"_pdb_nowater.txt" -p $PRMTOP

echo "Info: The output nc  has been written to $NCOUT"
echo "Info: The output pdb has been written to $PDBOUT"
echo "Done!"
echo ""
