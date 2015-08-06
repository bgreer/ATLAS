#!/bin/bash
# make a grid of points and print list to file
# $1 = output file
# $2 = clon
# $3 = clat
# $4 = lonrn
# $5 = latrn
# $6 = spacing
# $7 = max tiles per sub-file (0 = unlimited)

# examples:
# standard dense-pack
# ./grid.sh grid_file 180 0 90 90 7.5 0
# ultra-dense-pack
# ./grid.sh grid_file 180 0 90 90 0.25 7500

file=$1
clon=$2
clat=$3
lonrn=$4
latrn=$5
spacing=$6

numlon=`echo "scale=0; ${lonrn}/${spacing}+1" | bc -l 2>/dev/null`
numlat=`echo "scale=0; ${latrn}/${spacing}+1" | bc -l 2>/dev/null`

numpos=`echo "scale=0; ${numlon}*${numlat}" | bc -l 2>/dev/null`

tempfile="grid_tempfile"
rm -f ${tempfile}
touch ${tempfile}

echo "numlon,numlat = ${numlon},${numlat}"
echo "numpos = ${numpos}"

# set up division of tiles
if [ $7 -eq "0" ]
then
  echo "Outputting to single file."
  numfiles=1
  tilesperfile=${numpos}
else
  echo "Splitting output into multiple files."
  numfiles=`echo "scale=0; (${numpos} + $7 - 1)/$7" | bc -l 2>/dev/null`
  echo "Number of files: ${numfiles}"
  tilesperfile=`echo "scale=0; ${numpos}/${numfiles}" | bc -l 2>/dev/null`
  echo "Number of tiles per file (approx): ${tilesperfile}"
fi

currfile=0
ind=0

percent="0"
# header
#echo ${numpos} ${numlon} ${numlat} >> ${file}
for (( i=1; i<=${numlon}; i++ ))
do
	thislon=`echo "scale=5; ${i}*${spacing} + ${clon} - (${lonrn}/2)" | bc -l 2>/dev/null`
#	check=`echo "scale=0; (${i})%(${numlon}/10)" | bc -l 2>/dev/null`
#	if [ ${check} -eq "0" ]
#	then
#		percent=`echo "scale=0; ${percent}+10" | bc -l 2>/dev/null`
#		echo "${percent}%"
#	fi
	for (( j=1; j<=${numlat}; j++ ))
	do
		thislat=`echo "scale=5; ${j}*${spacing} + ${clat} - (${latrn}/2)" | bc -l 2>/dev/null`
		echo ${thislon} ${thislat} >> ${tempfile}
    ind=$((${ind}+1))
    # check for full file
    if [ ${ind} -ge ${tilesperfile} ]
    then
      if [ $7 -ne "0" ]
      then
        echo "Writing to file ${currfile}"
        # make file
        rm -f ${file}_${currfile}
        touch ${file}_${currfile}
        echo "${ind} 1 1" >> ${file}_${currfile}
        cat ${tempfile} >> ${file}_${currfile}
        # clear temp file
        rm -f ${tempfile}
        touch ${tempfile}
        ind=0
        currfile=$((${currfile}+1))
      fi
    fi
	done
done

# if no splitting files, wrap up grid file
if [ $7 -eq "0" ]
then
  echo ${numpos} ${numlon} ${numlat} >> ${file}
  cat ${tempfile} >> ${file}
fi

rm -f ${tempfile}
