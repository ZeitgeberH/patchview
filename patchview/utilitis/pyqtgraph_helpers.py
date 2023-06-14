def getAllChildIndexFromItem(item):
    '''Get all child from item'''
    children = []
    for i in range(item.childCount()):
        children.append(item.child(i))
    return children

def getAllSiblings(sel, includeSel=False):
    '''get all siblings of sel (not including sel)'''
    p = sel.parent()
    children = getAllChildIndexFromItem(p)
    if not includeSel:
        children.remove(sel)
    return children