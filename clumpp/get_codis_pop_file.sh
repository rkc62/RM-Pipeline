#!/bin/bash

out="codis/input/codis_popfile.txt"
> ${out}	


for i in {1..100}
do
    struc_result="../structure/str/output/str_${i}.str_f"
    awk 'NR==31,NR==56' ${struc_result} >> ${out}
done


sed 's/ 0:/ 26:/g' ${out} > tmp
mv tmp ${out}