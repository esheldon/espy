function usage {
    echo "
    usage: domatch.sh [-s] merun serun

        -s use stars for the SE data
        -h for this help message
"
}


stars="N"
while getopts "s" Option; do
    case $Option in
        s) stars="Y";;
        h) usage
           exit 1 ;;
        [?]) usage
             exit 1 ;;
    esac
done
shift $(($OPTIND - 1))


if [[ $# -lt 2 ]]; then
    usage
    exit 1
fi

merun=$1
serun=$2

rad=1
mef=$DESDATA/wlbnl/$merun/collated/$merun-radec.dat
if [[ $stars == "Y" ]]; then
    sef=$DESDATA/wlbnl/$serun/collated/$serun-stars-radec.dat
    out=$DESDATA/wlbnl/$merun/collated/match-$serun-stars-$merun.dat
else
    sef=$DESDATA/wlbnl/$serun/collated/$serun-radec.dat
    out=$DESDATA/wlbnl/$merun/collated/match-$serun-$merun.dat
fi

echo "will write to $out"
pv $sef | smatch -v -r $rad $mef > $out
