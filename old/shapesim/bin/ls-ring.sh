# list the yaml files sorted by second index
function usage {
    echo "
    usage: ls-ring.sh [-h -r] run

        -r reverse the order by s/n
        -t list bytrial
        -h for this help message
"
}


reverse="N"
bytrial="n"
while getopts "hrt" Option; do
    case $Option in
        r) reverse="Y" ;;
        t) bytrial="Y" ;;
        h) usage
           exit 1 ;;
        [?]) usage
             exit 1 ;;
    esac
done
shift $(($OPTIND - 1))

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

run=$1

if [[ $reverse == "Y" ]]; then
    for i in 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00; do
        if [[ $bytrial == "Y" ]]; then
            ls $run-*-0$i-*.yaml 2> /dev/null | sort
        else
            ls $run-*-0$i.yaml 2> /dev/null | sort
        fi
    done
else
    for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19; do
        if [[ $bytrial == "Y" ]]; then
            ls $run-*-0$i-*.yaml 2> /dev/null | sort
        else
            ls $run-*-0$i.yaml 2> /dev/null | sort
        fi
    done
fi
