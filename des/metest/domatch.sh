merun=$1
serun=$2

rad=2
sef=$DESDATA/wlbnl/$serun/collated/$serun-radec.dat
mef=$DESDATA/wlbnl/$merun/collated/$merun-radec.dat
out=$DESDATA/wlbnl/$merun/collated/match-$serun-$merun.dat

echo "will write to $out"
pv $sef | smatch -v -r $rad $mef > $out
