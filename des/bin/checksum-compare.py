"""
 Compare our checksum list on files to the master list
"""
import sys
def load_files(fname, is_master=False):
    if is_master:
        namestart=2
    else:
        namestart=27

    flist={}
    with open(fname) as fobj:
        for line in fobj:
            md5sum,name=line.split()
            name=name[namestart:]

            flist[name] = md5sum
    return flist        

def compare_md5sums(master_flist, test_flist):
    for name,md5sum in test_flist.iteritems():
        try:
            master_md5sum=master_flist[name]
            if md5sum != master_md5sum:
                print name,'has md5sum',md5sum,'instead of',master_md5sum
        except:
            print '%s not found in master list' % name
def main():
    if len(sys.argv) != 3:
        print "usage: checksum-compare.py master_list test_list"
        sys.exit(45)

    master_file=sys.argv[1]
    test_file=sys.argv[2]

    master_flist = load_files(master_file, is_master=True)
    test_flist = load_files(test_file, is_master=False)

    compare_md5sums(master_flist, test_flist)

main()
