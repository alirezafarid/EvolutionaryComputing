#!/bin/bash

javac -Xlint:unchecked  -Xlint:deprecation  -cp contest.jar player99.java individual.java CustomComparator.java
jar cmf MainClass.txt submission.jar player99.class individual.class CustomComparator.class 


if [ ! -f $1 ]; then 
	touch $1 
fi

if [ ! -d "results" ]; then
	mkdir results	
fi
 
 for i in `seq 1 $2`;
        do
       	java  -jar testrun.jar   -submission=player99  -evaluation=$3   -seed=$i >> $1
done    
 
mv $1 results/

python scripts/get_avg_EA.py $1 $3

rm results/$1
