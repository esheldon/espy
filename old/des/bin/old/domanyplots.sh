if [ $# -lt 1 ]; then
    echo domanyplots.sh serun
    exit 45
fi
serun=$1

# combined plots from all exposures
echo "Doing size mag"
echo "----------------------------------------------------------"
python plot_size_mag.py $serun
echo "Doing compare sizes"
echo "----------------------------------------------------------"
python compare-sizes.py $serun
echo "Doing whisker plots combined over all exposures"
echo "----------------------------------------------------------"
python plot_checkpsf.py --serun=$serun --types=allstarbin,allpsfbin,allbothbin

