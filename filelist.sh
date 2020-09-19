#!/bin/bash
rm pico_total.list runnumber.list
touch pico_nfs.list pico_local.list runnumber.list
get_file_list.pl -keys path,filename -cond production=P16ij,trgsetupname=AuAu_200_production_2016,filetype=daq_reco_picoDst,filename~st_physics,storage=nfs -limit 0 -delim / -distinct >> pico_nfs.list
get_file_list.pl -keys path,filename -cond production=P16ij,trgsetupname=AuAu_200_production_2016,filetype=daq_reco_picoDst,filename~st_physics,storage=local -limit 0 -delim / -distinct >> pico_local.list
sort -t '/' -k 10 pico_nfs.list > pico_nfs_sorted.list
sort -t '/' -k 12 pico_local.list > pico_local_sorted.list
cat pico_nfs.list pico_local.list > pico_total.list
cut -d '/' -f 10 pico_nfs_sorted.list | sort | uniq >> runnumber.list
cut -d '/' -f 12 pico_local_sorted.list | sort | uniq >> runnumber.list
rm pico_nfs.list pico_local.list pico_nfs_sorted.list pico_local_sorted.list

get_file_list.pl -keys 'runnumber' -cond production=P16id,trgsetupname=AuAu_200_production_high_2014, filetype=daq_reco_picoDst,filename~st_physics,storage=nfs -limit 0 >& runNumber_ht_high

#=========================================================
#Matt run14 mid, P16id, SL18f
get_file_list.pl -keys path,filename -cond production=P16id,filetype=daq_reco_picoDst,trgsetupname=AuAu_200_production_mid_2014,library=SL18f,storage!=HPSS -delim / -distinct | tee pico_nfs.list

#low, -limit 0 mean unlimited file number
get_file_list.pl -keys runnumber -cond production=P16id,filetype=daq_reco_picoDst,trgsetupname=AuAu_200_production_low_2014,library=SL18f,storage!=HPSS -limit 0 -delim / -distinct | tee pico_nfs_low_number_18h.list

#Run14 only have mid and low
get_file_list.pl -keys path,filename -cond production=P16id,filetype=daq_reco_picoDst,trgsetupname=AuAu_200_production_low_2014||AuAu_200_production_mid_2014||AuAu_200_production_high_2014||AuAu_200_production_2014,library=SL18f,storage!=HPSS -delim / -distinct | tee pico_nfs.list

# will not use 18ih for my Run14 analysis
get_file_list.pl -keys path,filename -cond production=P18ih,filetype=daq_reco_picoDst,trgsetupname=AuAu_200_production_low_2014,library=SL18h,storage!=HPSS -limit 0 -delim / -distinct | tee pico_nfs_low.list

get_file_list.pl -keys path,filename -cond 'trgsetupname=AuAu_200_production_mid\low_2014,filename~st_physics, production=p16id,library=SL18f, filetype=daq_reco_picoDst,storage!=HPSS' -limit 0

get_file_list.pl -keys path,filename -cond production=P12id,collision=pp200,trgsetupname=pp200_production_2012,library=SL18f,filename~st_physics,filetype=daq_reco_picoDst,storage!=HPSS -limit 0 -delim / -distinct  > tmpppp.list
