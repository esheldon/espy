# list the yaml files in order fastest to slowest
run=$1
for i in 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00; do
    ls $run-*-0$i.yaml | sort
done
