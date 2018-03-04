
if [ $# -lt 2 ]; then
	echo "usage: test-psfrec.sh serun exposurename [ ccd ]"
	echo "\te.g. test-psfrec.sh wlse0001 decam--25--12-i-7 35"
	exit 45
fi

config=$WL_DIR/etc/wl.config

serun=$1
exposurename=$2

if [ $# -gt 2 ]; then
	ccds=$3
else
	ccds=$(seq 1 62)
fi

dir=$DESDATA/wlbnl/$serun/$exposurename


for ccd in $ccds; do

	ccd=$( printf "%02d" $ccd)

	root=$dir/${serun}_${exposurename}_$ccd

	psf_file=${root}_psf.fits
	fitpsf_file=${root}_fitpsf.fits

	out_file=${root}_checkpsf.rec

	echo "output file: $out_file"
	test-psfrec $config $psf_file $fitpsf_file > $out_file
done
