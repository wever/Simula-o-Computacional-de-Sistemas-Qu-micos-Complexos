#!/bin/bash 

source /scratch/covidufscar/weverson.gomes3/gromacs_19/bin/GMXRC


find . -iname  pullx*.xvg > pullx.dat
find . -iwholename  '*npt/topol.tpr' > topol.dat

gmx_19 wham -it topol.dat  -ix pullx.dat  -o -hist ; qtgrace profile.xvg 
