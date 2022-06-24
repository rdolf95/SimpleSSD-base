#!/bin/bash
echo  >$(pwd)/log/stdout/pe3000.txt
echo  >$(pwd)/log/stdout/pe5000.txt
echo 'alibaba_150 alibaba_230 3000 5000'
python3 script_ali.py alibaba_150 alibaba_230 3000 >>$(pwd)/log/stdout/pe3000.txt 2>>$(pwd)/log/stderr/pe3000.txt &
python3 script_ali.py alibaba_150 alibaba_230 5000 >>$(pwd)/log/stdout/pe5000.txt 2>>$(pwd)/log/stderr/pe5000.txt &
#alibaba_430