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
