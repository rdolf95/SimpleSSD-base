#!/bin/bash
#echo 'alibaba_30 alibaba_130 alibaba_150 alibaba_230 alibaba_401 3000'

#trace_names=("alibaba_30" "alibaba_130" "alibaba_150" "alibaba_230" "alibaba_401")
#trace_names=("alibaba_12" "alibaba_100" "alibaba_130" "alibaba_507" "alibaba_538" "alibaba_730" "alibaba_731" "alibaba_743" "alibaba_792")
#trace_names=("alibaba_507")
trace_names=("alibaba_124" "alibaba_806" "alibaba_727" "alibaba_804" "alibaba_810")


# erase exisitng text of logs
for (( i = 0 ; i < ${#trace_names[@]} ; i++ )) ; do
    echo  >$(pwd)/log/stdout/no_refresh/"${trace_names[$i]}_3000.txt"
    echo  >$(pwd)/log/stderr/no_refresh/"${trace_names[$i]}_3000.txt"
done


for (( i = 0 ; i < ${#trace_names[@]} ; i++ )) ; do
    python3 script_ali_no_refresh.py ${trace_names[$i]} 3000 >>$(pwd)/log/stdout/no_refresh/"${trace_names[$i]}_3000.txt" 2>>$(pwd)/log/stderr/no_refresh/"${trace_names[$i]}_3000.txt" &
done


#alibaba_430