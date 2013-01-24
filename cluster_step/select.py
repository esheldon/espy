from numpy import where
class Selector(object):
    def __init__(self):
        pass

    def get_logic(self, data, setname):
        from esutil.numpy_util import strmatch

        logic = (data['flags']==0) | (data['flags'] == 2**16)

        logic=logic & (data['s2n_w'] > 10) & (data['s2n_w'] < 1.e6)

        # has minimal effect for use1
        #logic=logic & (data['Ts2n'] > 2)   & (data['Ts2n'] < 1.e6)
        # has minimal effect for use1
        #logic=logic & (data['Tmean'] > 2)  & (data['Ts2n'] < 1.e6)

        if setname=='use0':
            # default set
            pass
        elif setname=='use1':

            isdev=strmatch(data['model'],'gdev')
            isexp=strmatch(data['model'],'gexp')

            Ts2n20=(data['Ts2n'] > 20)

            logic=logic & (isexp | (isdev & Ts2n20) )

        elif setname=='use2':

            isexp=strmatch(data['model'],'gexp')
            logic = logic & isexp

        elif setname=='use3':
            Ts2n20=(data['Ts2n'] > 20)
            logic = logic & Ts2n20

        elif setname=='use4':
            isexp=strmatch(data['model'],'gexp')
            Ts2n20=(data['Ts2n'] > 20)
            logic = logic & isexp & Ts2n20

        elif setname=='use5':
            # same as use1 but with Ts2n > 3.5
            isdev=strmatch(data['model'],'gdev')
            isexp=strmatch(data['model'],'gexp')

            Ts2n20=(data['Ts2n'] > 20)
            Ts2n3_5=(data['Ts2n'] > 3.5)

            logic=logic & ((isexp & Ts2n3_5) | (isdev & Ts2n20) )

        else:
            raise ValueError("bad use name: '%s'" % setname)

        return logic

    def select(self, data, setname):
        logic=self.get_logic(data, setname)
        w,=where(logic)
        return w

#parser.add_option('--s2n',default='10,200', help="s/n range, %default")
#parser.add_option('--Ts2n',default='2,1.e6', help="Ts2n range, %default")
#parser.add_option('--sratio',default='1.0,1.e6',
#                  help='sratio range, %default')
#parser.add_option('--Tmean',default='2,1.e6',
#                  help='Tmean range, %default')
#parser.add_option('--mag',default='0,100',
#                  help='mag range, %default')
#parser.add_option('--prob',default='0,1',
#                  help='fit_prob range, %default')
#parser.add_option('-t','--type',default=None,
#                  help="limit to objects best fit by this model")


