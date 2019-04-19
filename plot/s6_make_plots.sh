#! /bin/bash
fitsfile=$1
echo ${fitsfile}
if [ ! -e $1.plotted ]
then
	echo ...waterfall
	python /usr/local/bin/s6_plot_waterfall.py ${fitsfile}
	echo ...power histogram
	python /usr/local/bin/s6_plot_power_histogram.py ${fitsfile}
	echo ...baseband
	python /usr/local/bin/s6_plot_baseband.py ${fitsfile}
	echo ...meanpower
	python /usr/local/bin/s6_plot_meanpower.py ${fitsfile}
	touch $1.plotted
fi


