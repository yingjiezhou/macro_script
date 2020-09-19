#!/bin/bash

workdir="/star/u/pengliu/gpfs01/He4/R14AuAu200GeV/low"
rootfiledir="${workdir}/output1"
cshfiledir="${workdir}/input/part1/runinfo/csh"
sessfiledir="${workdir}/input/part1"

id="9BC05E069298FC662492BE0AA53E0773"

n=0
i=0
failedID=""

while [ $i -lt 8845 ]; do
	
	if [ ! -e ${rootfiledir}/${id}_$i.deuteron.histo.root ]; then
		echo "${id}_$i.deuteron.histo.root doesn't exist.............."

		failedID="${failedID}${i},"

		let "n+=1"
	fi
    
	let "i+=1"
done

echo "failed job ID:"
echo "${failedID%?}"

if [ $n -gt 0 ]; then
	
	star-submit -r  ${failedID%?}  ${sessfiledir}/${id}.session.xml

fi

echo "Total number of files missing = $n"
