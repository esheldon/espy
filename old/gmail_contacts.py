"""
    %prog style filename

The contacts go on stdout
"""
import csv
import sys

from optparse import OptionParser
parser=OptionParser(__doc__)

def extract_mutt(fname):
    """
    Convert a gmail.csv contacts file into a mutt alias file.  Attempt to
    create unique a unique alias for each entry. This is created from the
    name entry, filling spaces with underscores, or from the email address
    when that is empty.
    """
    reader = csv.reader(open(fname))
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


def extract_sup(fname):
    """
    Convert a gmail.csv contacts file into a sup contacts.
    First Name,Middle Name,Last Name,Title,Suffix,Initials,Web Page,Gender,Birthday,Anniversary,Location,Language,Internet Free Busy,Notes,E-mail Address
    """
    reader = csv.reader(open(fname))
    first=True
    for row in reader:
        if first:
            first=False
            continue

        # Name and email address
        email = row[1]
        first_name = row[0].strip()
        middle_name = row[1].strip()
        last_name = row[2].strip()
        email = row[14].strip()

        email = email.replace('"','').strip()
        if email == '':
            continue

        
        full_name=[]
        if first_name != '':
            full_name.append(first_name)
        if middle_name != '':
            full_name.append(middle_name)
        if last_name != '':
            full_name.append(last_name)

        full_name=' '.join(full_name)

        if email.lower() == full_name.lower():
            full_name=''

        if full_name != '':
            email = '<%s>' % email
            print full_name,email
        else:
            print email


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    style=args[0]
    fname=args[1]

    if style=='mutt':
        extract_mutt(fname)
    elif style=='sup':
        extract_sup(fname)
    else:
        raise ValueError("bad style: %s" % style)

main()
