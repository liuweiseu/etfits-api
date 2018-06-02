#! /bin/bash

export data_dir=/home/jeffc/build/serendip6/src

for fitsfile in `ls ${data_dir}/*.working`; do
    #echo $fitsfile
	if [ ! -e $fitsfile.plotted ]
	then
		python ~/build/etfits_api/plot/s6_plot_waterfall.py $fitsfile
		python ~/build/etfits_api/plot/s6_plot_power_histogram.py $fitsfile
		python ~/build/etfits_api/plot/s6_plot_baseband.py $fitsfile
		touch $fitsfile.plotted
	fi
done
