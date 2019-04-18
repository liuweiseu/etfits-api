#! /bin/bash

for fitsfile in `ls *.fits`; do
    echo $fitsfile
	if [ ! -e $fitsfile.plotted ]
	then
		echo ...waterfall
		python /usr/local/bin/s6_plot_waterfall.py $fitsfile
		echo ...power histogram
		python /usr/local/bin/s6_plot_power_histogram.py $fitsfile
		echo ...baseband
		python /usr/local/bin/s6_plot_baseband.py $fitsfile
		echo ...meanpower
		python /usr/local/bin/s6_plot_meanpower.py $fitsfile
		touch $fitsfile.plotted
	fi
done
