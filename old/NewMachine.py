import sys

NameDict = {'greeting     ciao':0,
            'greeting       yo':0,
            'geek          sol':0,
            'geek         halo':0,
            'geek         lens':0,
            'geek        metal':0,
            'drinks        rum':0,
            'drinks       beer':0,
            'drinks        gin':0,
            'drinks        rye':0,
            'misc       falcon':0,
            'misc        wotan':0,
            'misc         doom':0,
            'misc        sloth':0,
            'misc         silk':0}

Andreas = {'greeting     ciao':1,
           'greeting       yo':1,
           'geek          sol':0,
           'geek         halo':0,
           'geek         lens':-1,
           'geek        metal':0,
           'drinks        rum':0,
           'drinks       beer':0,
           'drinks        gin':0,
           'drinks        rye':-1,
           'misc       falcon':0,
           'misc        wotan':1,
           'misc         doom':1,
           'misc        sloth':1,
           'misc         silk':-1}

Erin = {'greeting     ciao':1,
        'greeting       yo':0,
        'geek          sol':1,
        'geek         halo':1,
        'geek         lens':0,
        'geek        metal':1,
        'drinks        rum':-1,
        'drinks       beer':0,
        'drinks        gin':1,
        'drinks        rye':0,
        'misc       falcon':1,
        'misc        wotan':0,
        'misc         doom':1,
        'misc        sloth':0,
        'misc         silk':1}


Michael = {'greeting     ciao':0,
           'greeting       yo':0,
           'geek          sol':-1,
           'geek         halo':1,
           'geek         lens':-1,
           'geek        metal':0,
           'drinks        rum':-1,
           'drinks       beer':-1,
           'drinks        gin':-1,
           'drinks        rye':-1,
           'misc       falcon':-1,
           'misc        wotan':-1,
           'misc         doom':-1,
           'misc        sloth':0,
           'misc         silk':0}

Beth = {'greeting     ciao':1,
        'greeting       yo':1,
        'geek          sol':0,
        'geek         halo':1,
        'geek         lens':0,
        'geek        metal':1,
        'drinks        rum':-1,
        'drinks       beer':-1,
        'drinks        gin':-1,
        'drinks        rye':-1,
        'misc       falcon':-1,
        'misc        wotan':-1,
        'misc         doom':-1,
        'misc        sloth':+1,
        'misc         silk':+1}

Roman = {'greeting     ciao':0,
         'greeting       yo':0,
         'geek          sol':0,
         'geek         halo':1,
         'geek         lens':0,
         'geek        metal':0,
         'drinks        rum':0,
         'drinks       beer':0,
         'drinks        gin':0,
         'drinks        rye':0,
         'misc       falcon':0,
         'misc        wotan':0,
         'misc         doom':0,
         'misc        sloth':0,
         'misc         silk':0}


if __name__ == "__main__":
    sys.stdout.write("%17s %10s\n" % ("name","score"))
    sys.stdout.write("-"*28+'\n')
    for name in sorted(NameDict):
        NameDict[name] += Andreas[name]
        NameDict[name] += Erin[name]
        NameDict[name] += Michael[name]
        NameDict[name] += Beth[name]
        NameDict[name] += Roman[name]
        sys.stdout.write("%17s %10s\n" % (name,NameDict[name]))
