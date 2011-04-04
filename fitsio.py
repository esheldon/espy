import sys

class FITS(dict):
    def __init__(self, filename):
        self['filename'] = filename
        self.fobj = open(filename, 'r')
        self.hdus = []

    def read_header(self, ext):
        if ext != 0:
            raise ValueError("just doing primary hdu now")
        self.fobj.seek(0)
        mult=2880

        cards = []
        h = FITSHeader()
        while True:
            tmp = self.fobj.read(mult)
            h.load_cards(tmp)
            if h.end_card != None:
                return h


class FITSHeader(dict):
    def __init__(self):
        self.cards=[]
        self.header_text=''
        self.end_card=None

    def load_cards(self, text):
        """
        Load cards from the input text.  If a keyword is END then
        stop and return, other wise process until the end of the chunk.

        the text chunk should be 2880 characters in size.
        """

        self.header_text += text

        tsize = len(text)
        if tsize != 2880:
            raise ValueError("Each header chunk must be 2880 chars in size")
        cardsize=80

        # this is 36
        ncards_max = 2880/cardsize

        for i in xrange(ncards_max):
            i1 = i*cardsize
            i2 = (i+1)*cardsize
            tmp = text[i1:i2]

            card = FITSCard(tmp)

            if card.type == 'end':
                self.end_card=i

            self.cards.append(card)

    def __repr__(self):
        rep=""
        for c in self.cards:
            rep+=c.text + '\n'
            if c.type == 'end':
                return rep


class FITSCard:
    def __init__(self, text):
        self.text=text

        self.type=None
        self.key=None
        self.val=None
        self.comment=None

        self._process_text()

    def _process_text(self):

        card_text = self.text.strip()
        if card_text == '':
            return

        # deal with special fields that don't have a = sign
        keyfield=card_text[0:8]
        if keyfield[0:4] == 'END':
            self.type='end'
            return
        if keyfield[0:7] == 'COMMENT':
            # note we are going back to unstripped version
            card.type='comment'
            card.val = self.text[7:]
            return
        elif keyfield[0:7] == 'HISTORY':
            # note we are going back to unstripped version
            card.type='history'
            card.val = self.text[7:]
            return
        elif keyfield.strip() == '':
            # standard allows key to be blank but values to be
            # there, but I'm going to skip it
            return

        if card_text[8] != '=':
            # this is not a keyword-value field, probably
            self.type = 'unknown'
            return

        self.key = card_text[0:8].strip()
        if self.key[0:5] == 'TTYPE':
            self.type = 'ttype'
        elif self.key[0:5] == 'TFORM':
            self.type = 'tform'
        elif self.key[0:4] == 'TDIM':
            self.type = 'tdim'
        else:
            self.type = 'keyword'

        val = card_text[9:].strip()

        comment_id = val.find('/')
        if comment_id != -1:
            # this will be '' if the / is at the last char
            self.comment = val[comment_id+1:]
            val = val[0:comment_id]
        self.valstring = val.strip()
        try:
            self.val = eval(self.valstring)
        except:
            self.val = self.valstring


    def __repr__(self):
        rep="""type:    {type}
key:     {key}
val:     {val}
comment: 
{comment}""".format(type=self.type,
                   key=self.key,
                   val=self.val,
                   comment=self.comment)
        return rep

if __name__=='__main__':
    from optparse import OptionParser
    parser=OptionParser(__doc__)
    
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) > 0:
        filename=args[0]
        if len(args) > 1:
            ext=int(args[1])
        else:
            ext=0

        f = FITS(filename)
        f.read_header_text(ext)

