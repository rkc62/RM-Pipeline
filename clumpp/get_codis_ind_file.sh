#!/bin/bash

out="codis/input/str_indfile.txt"
> ${out}

for i in {1..100}
do
    struc_result="../structure/codis/output/codis_${i}.str_f"
    awk 'NR==82,NR==2585' ${struc_result} >> ${out}
done

awk '{$2=gensub(/[^0-9]+/,"","g",$2); print}' "${out}" > out1
awk -v c=4 '{if($c == 0){$c = 26} print $0}' out1 > out2

rm out1

mv out2 ${out}
