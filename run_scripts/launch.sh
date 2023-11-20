#!/bin/bash

python EXECPROC/Preproc.py
python EXECPROC/Preproc_pdaf.py

cp HGS/Grokfiles/wells/wells_flow_1.inc HGS/Grokfiles/wells_flow.inc
python EXECPROC/Spinup.py


for i in $(seq 1 1 95) 
do 
cp HGS/Grokfiles/wells/wells_flow_$i.inc HGS/Grokfiles/wells_flow.inc

python EXECPROC/Execproc.py
python EXECPROC/obshgs2obspest_Sequential.py

preplot HGS/Flowo.pm

cp HGS/Flowo.head_pm.0001 HGS/IC/Ini_pm_SS
cp HGS/Flowo.head_olf.0001 HGS/IC/Ini_olf_SS

done
#SLEEP 3
#)

