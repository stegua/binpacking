#!/bin/bash

for MYFILE in `ls /data/gualandi/MyDATA/gap/*$1`; do 
   export MYFILE; 
   qsub my_pbs -v MYFILE; 
done;
