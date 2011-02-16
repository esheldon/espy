# This demonstrates that scalar (including strings) and
# typle parameters cannot be modified, but lists and
# dictionaries can. 

# Note: You cannot pass an undefined parameter either.
# Thus you cannot do the IDL style "return all optional
# outputs through keywords" thing unless the user first
# defines them.  A workaround coulde be passing a dictionary
# or returning a tuple with the dictionary of optional things
# as a member.

def test_scalar(scalar=3):
    scalar = 5

def test_str(str=""):
    str = "goodbye"

def test_list(list=[]):
    list.append(35)

def test_tuple(typle=()):
    tuple = (5,)

def test_dict(dict={}):
    if not 'a' in dict:
        dict['a'] = 3



    
if __name__ == "__main__":

    scalar = 27
    print "scalar = ",scalar
    test_scalar(scalar)
    print "Now scalar = ",scalar

    str = "hello world"
    print "\nstr = ",str
    test_str(str)
    print "Now str = ",str

    list = [83]
    print "\nlist = ",list
    test_list(list)
    print "Now list = ",list

    tuple = (66,)
    print "\ntuple = ",tuple
    test_tuple(tuple)
    print "Now tuple = ",tuple



    dict = {'c':35}
    print "\ndict = ",dict
    test_dict(dict)
    print "Now dict = ", dict
