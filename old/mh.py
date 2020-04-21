import math
import sys

def sizes(rib_len=10.0, spine_len=10.0):
    # there are two end-to-end
    total_rib_len=rib_len*2
    rib_width=1.0             # inches

    # end-to-end
    total_spine_len = spine_len
    spine_width=1.25        # inches

    # arclen = total_rib_len = pi*r
    # so projected onto the ground is 2*r
    height = total_rib_len/math.pi
    rib_ground_footlen = 2*height
    spine_ground_footlen = total_spine_len
    sys.stdout.write('height is: %s\n' % (height,))
    sys.stdout.write('Footprint is: %sX%s = %s\n' % (rib_ground_footlen,
                                                     spine_ground_footlen,
                                                     rib_ground_footlen*spine_ground_footlen))

    my_area = 7.5*12.7
    sly_area = 11.8*10.0
    sys.stdout.write('\nOur tent areas:\n')
    sys.stdout.write('  Mine: 7.5x12.7 %s\n' % my_area)
    sys.stdout.write('  slys: 11.8x10: %s\n' % sly_area)
    sys.stdout.write('  Total area: %s\n' % (my_area+sly_area))

    # treating as squares which is more right since that is how
    # they stack.  0.1 for thickness of the plastic
    rib_area_max = ((rib_width+0.1)/12.0)**2
    spine_area_max = ((spine_width+0.1)/12.0)**2

    cubic_feet_max = (rib_len*6*rib_area_max) + (spine_len*2*spine_area_max)
    sys.stdout.write('\nVolume in storage max: %s cubit feet\n' %
                     (cubic_feet_max,))

sizes()

