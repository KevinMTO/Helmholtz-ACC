#!/bin/bash


for d in $(find . -maxdepth 2 -type d)
do
    for file in $d/*
    do  
        if [ ! -d "$file" ]; then
                    touch temp.txt
                    grep MFlops $file> temp.txt

                    echo
                    echo "processed"
                    echo "$file"
                    echo

                    count=0
                    total=0
                    totalSQ=0 

                    for i in $( awk '{ print $3; }' temp.txt )
                    do 
                        total=$(echo "$total+$i" | bc )
                        ((count++))
                    done

                    AVG=$(echo "scale=3; ($total / $count)" | bc)

                    sumofdiff=0
                    for i in $( awk '{ print $3; }' temp.txt )
                    do 
                        diff=$(echo "$i-$AVG" | bc )
                        diff2=$(echo "$diff*$diff" | bc )
                        sumofdiff=$(echo "$sumofdiff+$diff2" | bc )
                    done

                    VARIANCE=$(echo "scale=3; ($sumofdiff / $count)" | bc)

                    echo "${file##*/} VARIANCE ${VARIANCE} AVG ${AVG} " >> $d/sumup.txt

                    rm -f temp.txt
        fi
    done 
done 

