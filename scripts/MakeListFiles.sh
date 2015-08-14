#!/bin/bash

# Input:
#   Directory of files
#   Number of dops needed in the end
#	Filename (doplist)
#	Cadence in seconds
# Output:
#	list of dopplergram filenames

# example:
# ./MakeListFiles.sh DOPS/ 2048 doplist 45

source float.sh # to do floating point math

echo "MakeList Start (Files)"
echo

maindir=$1
numneeded=$2
doplist=$3
cadence=$4

num=`ls ${maindir}/*.fits | wc | awk '{print $1}'`

echo "num in dir = ${num}"
echo "num needed = ${numneeded}"

if [ ${num} -eq ${numneeded} ]
then
	echo "all dops present, making list"
	rm -f ${doplist}
	touch ${doplist}
	for (( i=1; i<=${num}; i++ ))
	do
		item=$(ls ${maindir}/*.fits | head -n ${i} | tail -n 1)
#		echo '"'${maindir}/${item}'"' >> ${doplist}
		echo '"'${item}'"' >> ${doplist}
	done
fi
prevtime=-1
if [ ${num} -lt ${numneeded} ]
then
	echo "not all dops present, checking times"
	rm -f ${doplist}
	touch ${doplist}
	for (( i=1; i<=${num}; i++ ))
	do
		item=$(ls ${maindir}/*.fits | head -n ${i} | tail -n 1)
#		echo "./DopTime ${item}"
		time0=`./DopTime ${item} | head -n 2 | tail -n 1 | sed -e 's/[.]/-/g' | sed -e 's/_/ /g' | sed -e 's/TAI/ /g'`
#		echo $time0
		currtime=`date -d "${time0}" +"%s"` # unix time in secs
		if [ ${prevtime} -gt 0 ]
		then
			while [ $((${currtime}-${prevtime})) -gt ${cadence} ]
			do
				echo "diff = "$((${currtime}-${prevtime}))
				prevtime=$((${prevtime}+${cadence}))
				echo '"0"' >> ${doplist}
			done
		fi
		echo '"'${item}'"' >> ${doplist}
		prevtime=${currtime}
	done
fi

# do final check
numfinal=`wc ${doplist} | awk '{print $1}'`
if [ ${numfinal} -ne ${numneeded} ]
then
	echo "error, needed ${numneeded}, only made ${numfinal}"
fi
