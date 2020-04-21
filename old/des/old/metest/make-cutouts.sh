if [[ $# -lt 2 ]]; then
    echo "usage: make-cutouts.sh run s2nmin"
    exit 1
fi

run=$1
s2nmin=$2

d=$DESDATA/wlbnl/$run/collated/cutouts
infile=$d/s2n-gt-$s2nmin.dat
outd=$d/s2n-gt-$s2nmin

if [[ ! -e $outd ]]; then
    mkdir $outd
fi

front='(NR > 1) {print $7,1,$5,$4,51,"'
back='/" $1 "-" $2 ".png"}'
awkprog="${front}${outd}${back}"
echo "$awkprog"
head -1000 $infile | awk "$awkprog" | imcutout -r -8,50000

phpshow=$outd/phpshow.php
if [[ ! -e $phpshow ]]; then
    cp -v ~/www/tmp/phpshow/phpshow.php $outd/
fi
