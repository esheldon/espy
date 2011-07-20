def tselect(t):
    import sdsspy
    from esutil.numpy_util import where1

    plogic = (t['resolve_status'] & sdsspy.flagval('resolve_status','survey_primary')) != 0
    blended_logic = (t['objc_flags'] & sdsspy.flagval('object1','blended')) == 0
    nodeblend_logic = (t['objc_flags'] & sdsspy.flagval('object1','nodeblend')) != 0

    pw = where1( plogic )
    aw = where1( plogic & (blended_logic | nodeblend_logic) )

    print 'primary:',pw.size
    print 'primary + blend cuts:',aw.size
def qsort(lst):
    n = len(lst)
    if n <= 1:
        return lst

    pivot = lst[n/2]
    del lst[n/2]

    low = []
    high =[]

    for el in lst:
        if el < pivot:
            low.append(el)
        else:
            high.append(el)
    
    return qsort(low) + [pivot] + qsort(high)
