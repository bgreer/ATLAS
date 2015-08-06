#!/bin/bash
# use this script to add quotes around every line in a file

# $1 is the input file
# $2 is the output file

numlines=`wc $1 | awk '{print $1}'`

echo "${numlines} lines found, putting quotes around them."

rm -f $2
touch $2

for (( i=1; i<=${numlines}; i++ ))
do
	line=`head -n ${i} $1 | tail -n 1`
	echo "\"${line}\"" >> $2
done
