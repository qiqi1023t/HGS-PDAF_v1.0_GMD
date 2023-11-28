#!/bin/bash

python EXECPROC/Preproc.py
python EXECPROC/Preproc_pdaf.py

cp HGS/Grokfiles/wells/wells_flow_1.inc HGS/Grokfiles/wells_flow.inc
python EXECPROC/Spinup.py

