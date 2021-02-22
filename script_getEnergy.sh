#!/bin/bash 

source /scratch/covidufscar/weverson.gomes3/gromacs_19/bin/GMXRC


fullPathFiles=$(find . -type f -name "*npt.*edr" | sort)
rootPath=`pwd`
GMX=gmx_19

echo "ALL FILES: " ${fullPathFiles}

folder=data    

mkdir ${folder}
property="Enthalpy"

n=0
for d in ${fullPathFiles}; 
do
    dirPath=$(dirname "${d}")
    fileName=$(basename -- "${d}")
    outFilename="${fileName%.*}"

    echo "$dirPath" "${outFilename}"; 
    
    cd ${dirPath}
    echo "enthalpy" | $GMX energy -f ${fileName} -o ${rootPath}/${folder}/${n}.${property} -xvg none

    $GMX distance -s topol.tpr -f ${outFilename}.xtc -select 'com of group "r_1" plus com of group "r_10"' -xvg none -n ${rootPath}/residues1_10.ndx -oav ${rootPath}/$folder/${n}.dist
    cd ${rootPath}
    
    n=$((n+1))
done

