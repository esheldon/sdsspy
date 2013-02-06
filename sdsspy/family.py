"""
classes
-------
Family:
    A class to find the family of deblended objects

functions
---------
get_family
    Given the input ids or catalog, find the family
get_parent
    Given the id info, find the parent and return its
    photoObj structure
"""

from numpy import array, where, concatenate
import sdsspy

class Family(object):
    """
    Reconstruct the family list for an object

    examples
    --------
    fam=Family(struct, index)

    Get the indices into struct for each type of family member

        children=fam.get_children()
        parent=fam.get_parent()
        grandparent=fam.get_grandparent()

        ra=struct['ra'][children]

    """
    def __init__(self, data, index, verbose=False):
        self.data=data
        self.index=index
        self.verbose=verbose
        self.obj=self.data[index]

        wfield,=where(data['field'] == self.obj['field'])
        self.wfield=wfield

        self.find_family()

    def __repr__(self):
        lines=[]
        lines.append("  children: %s" % self.wchild)
        lines.append("  parent:   %s" % self.wparent)
        lines.append("  grandparent:   %s" % self.wgrandparent)
        return "\n".join(lines)

    def get_all(self):
        return concatenate( (self.wgrandparent, self.wparent, self.wchild) )

    def get_children(self):
        return self.wchild
    def get_parent(self):
        return self.wparent
    def get_grandparent(self):
        return self.wgrandparent

    def empty_array(self):
        return array([], dtype='i8')

    def find_family(self):
        """
        work our way up to a parent, possibly a grandparent parent
        """

        obj=self.data[self.index]

        fdata=self.data[self.wfield]

        findex,=where(fdata['id']==obj['id'])

        # indices are in the field at this point
        wgrandparent=self.empty_array()
        if obj['nchild'] > 0:
            # this is a parent

            if obj['nchild']==1:
                # it is a grandparent parent
                wgrandparent=findex.copy()
                # the real parent is the one who has me listed as parent
                wparent,=where(fdata['parent']==obj['id'])
                if wparent.size==0:
                    # the faint version simply isn't here
                    # we can go no further
                    wchild=self.empty_array()
                else:
                    # now the children have this parent listed
                    wchild,=where(fdata['parent']==fdata['id'][wparent])

                mess='is grandparent parent with %s children' % obj['nchild']
            else:
                # normal parent
                wparent=findex.copy()

                # I may have a grandparent version
                wgrandparent,=where(fdata['id'] == obj['parent'])

                # children list me as parent
                wchild,=where(fdata['parent']==obj['id'])
                mess='is parent with %s children' % obj['nchild']
        else:
            # this is a child
            if obj['nchild'] != 0:
                # this also has children: thus it is is the faint version of a
                # grandparent object.  We call it the parent
                wparent=findex.copy()

                # the grandparent must be the parent of this object
                wgrandparent,=where(fdata['id']==obj['parent'])

                # and children are now its children
                wchild,=where(fdata['parent']==obj['id'])
                mess='is child of %s with %s children' % (obj['parent'],obj['nchild'])
            else:
                # normal child
                wparent,=where(fdata['id']==obj['parent'])
                wchild,=where(fdata['parent']==obj['parent'])

                # we can try to find a grandparent parent if the parent
                # is found
                if wparent.size > 0:
                    if fdata['parent'][wparent] != -1:
                        # a grandparent parent does exist
                        wgrandparent,=where(fdata['id']==fdata['parent'][wparent])
                mess="is child with %s siblings" % (wchild.size-1)
                       

        if wgrandparent.size > 0:
            wgrandparent=self.wfield[wgrandparent]
        if wparent.size > 0:
            wparent=self.wfield[wparent]
        if wchild.size > 0:
            wchild=self.wfield[wchild]

        self.wgrandparent=wgrandparent
        self.wparent=wparent
        self.wchild=wchild

        if self.verbose:
            print mess

        if self.verbose > 1:
            print 'grandparent:  ',self.wgrandparent
            print 'parent:       ',self.wparent
            print 'children:     ',self.wchild


def get_parent(**keys):
    """
    Read in a photoObj structure for the parent of the input object

    You must enter as keywords objid or photoid or run,rerun,camcol,field,id or
    a structure with those.
    """
    fam=get_family(**keys)
    parent=fam.get_parent()
    if parent.size == 0:
        raise ValueError("Could not read parent")
    return fam.data[parent[0]]

def get_family(**keys):
    """
    Given the input ids, find the family.
    
    Returns a Family object.  Reads the photoobj files to get
    the data and find the family.

    You must enter as keywords objid or photoid or run,rerun,camcol,field,id

    """
    from .util import get_id_info, get_photoid

    data=sdsspy.read('photoObj',lower=True, **keys)
    if data is None:
        raise ValueError("failed to read any photoObj files for this id")

    ids=get_id_info(**keys)
    pids=get_photoid(data)

    index,=where(pids == ids['photoid'])
    if index.size == 0:
        raise ValueError("result doesn't contain object?")

    fam=Family(data, index[0])

    return fam

