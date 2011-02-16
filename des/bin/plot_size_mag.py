import des
import sys

if len(sys.argv) < 2:
    sys.stdout.write('plot_size_mag.py run\n')
    sys.exit(45)

run = sys.argv[1]
s = des.select.SizeMagSelector(run)
s.plot()
