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

echo
echo "Convert PDB water to AMBER water"
echo "Version 2.0 2019-12-15"
echo 
    
if [ ! -n "$2" ];then
    echo "Usage: $0 <input.pdb> <output.pdb>"    
    echo ""
    exit
fi

INPUT=$1
OUTPUT=$2

echo "Info: Reading input PDB data from $INPUT"

awk 'BEGIN {RESI=""; RESN=""; RESI_PREV=""; H=0; COUNT_RES=0; COUNT_ATOMS=0; LINES_UNTOUCHED=0;}
/ATOM|HETATM/ {
    RESI=substr($0,23,4);
    RESN=substr($0,18,3);
    if (RESN=="WAT" || RESN=="HOH") {
        BEFORE_ATOM=substr($0,1,12)
        ATOM=substr($0,13,4)
        AFTER_ATOM=substr($0,17)
        
        gsub(" ","", ATOM)
        
        if (RESI_PREV != RESI) {
           H=0
           COUNT_RES=COUNT_RES+1
        }
        
        if (ATOM=="HW" && RESI_PREV==RESI && H==0) {
            ATOM=" H1 "
            H=H+1
            COUNT_ATOMS=COUNT_ATOMS+1
        }
        else if (ATOM=="HW" && RESI_PREV==RESI && H==1) {
            ATOM=" H2 "
            H=H+1
            COUNT_ATOMS=COUNT_ATOMS+1
        }
        else if (ATOM=="OW") {
            ATOM=" O  "
            COUNT_ATOMS=COUNT_ATOMS+1
        }
        else {
            print "Warning: RESN="RESN" RESI="RESI" ATOM="ATOM" H="H" (<- hydrogens count in the residue) RESI_PREV="RESI_PREV""  | "cat >&2"
        }
        
        RESI_PREV = RESI
        print BEFORE_ATOM""ATOM""AFTER_ATOM
    }
    else {
        print;
        ATOM_LINES_UNTOUCHED=ATOM_LINES_UNTOUCHED+1
    }
    next;
}
// {
     print;
     NONATOM_LINES_UNTOUCHED=NONATOM_LINES_UNTOUCHED+1;
}
END {
    print "Info: Renamed "COUNT_ATOMS" atoms in "COUNT_RES" solvent residues"  | "cat >&2"
    print "Info: Printed "ATOM_LINES_UNTOUCHED" original lines with atom coordinates of non-water residues"  | "cat >&2"
    print "Info: Printed "NONATOM_LINES_UNTOUCHED" original lines with other data"  | "cat >&2"

}' $INPUT > $OUTPUT

echo "Info: Printed the output PDB data to $OUTPUT"
#echo 
#echo "NB!"
#echo "Don't forget to add 'HOH=TP4' or 'WAT=TP4' line to the leap configuration file"
#echo
echo "Done!"
echo ""
