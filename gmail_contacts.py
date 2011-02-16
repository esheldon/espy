import csv
import sys

if __name__=="__main__":
    """
    Convert a gmail.csv contacts file into a mutt alias file.  Attempt to
    create unique a unique alias for each entry. This is created from the
    name entry, filling spaces with underscores, or from the email address
    when that is empty.
    """
    if (len(sys.argv) < 2):
        sys.stdout.write('-syntax: gmail_contacts gmail.csv\n')
    else:
        reader = csv.reader(open(sys.argv[1],"r"))
        first=True
        for row in reader:
            if first:
                sys.stdout.write('#\n')
                sys.stdout.write('# ' + row[0]+'   '+row[1]+'\n')
                sys.stdout.write('# --------------------------------------------\n')
                first=False
            else:
                # Name and email address
                email = row[1]
                name = row[0]
                if name == '':
                    name=email

                # Create an alias
                alias = name.replace(' ','_')
                email = '<'+email+'>'
                
                sys.stdout.write('alias  %-40s  %-40s  %-40s\n' %
                                 (alias,name,email))


