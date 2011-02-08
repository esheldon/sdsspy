"""
    Module:
        sdsspy.files

    Description:
        A set of functions for dealing with SDSS files, including reading files,
        finding file lists, getting run metadata, converting SDSS id numbers 
        to strings, etc.

    Some Useful Functions. 

        read(ftype, run=None, camcol=None, field=None, id=None, **keys)
            Read an SDSS file given the input file type and id info

        filename(ftype, run=None, camcol=None, field=None, **keys)
            Generate an SDSS file name from input file type and id info

        filedir(ftype, run=None, camcol=None, **keys)
            Generate an SDSS dir path from input file type and id info

        filespec(ftype)
            Return the path specificition for the in input file type.  This information
            read from
                $SDSSPY_DIR/share/sdssFileTypes.par

        file_list(ftype, run=None, camcol=None, field=None, **keys)
            Get lists of SDSS files based on id info.

        runlist()
            Return the sdss runlist.  This is read from $PHOTO_REDUX/runList.par
            This data is cached and only read once, so is preferable to using
            sdsspy.files.read('runList')

        expand_sdssvars(string, **keys):
            Convert environment variables and sdss vars, e.g. $RUNNUM, in file
            specifications to values.  Used to generate file locations.

        Routines for converting SDSS id numbers to strings:
            These functions return scalars for scalar input, otherwise lists. 

            run2string(runs)
                Convert input run to a string padded to 6 zeros.
            field2string(fields)
                Convert input field to a string padded to 4 zeros.
            id2string(ids)
                Convert input id to string
            camcol2string(camcols, band=None)
                Convert input camcol to a string.
            filter2string(bands)
                Return alphabetic representation of filter
                  fstr = filter2string(2)      # returns 'r'
                  fstr = filter2string('2')    # returns 'r'
                  fstr = filter2string('r')    # returns 'r'
                  fstr = filter2string('R')    # returns 'r'

            rerun2string(reruns)
                Convert input rerun to string
            stripe2string(stripes)
                Convert input stripe to string

    Dependencies:
        numpy
        esutil for reading files and struct processing.


    Modification History:
        Documented: 2007-06-01, Erin Sheldon
"""
import os
import sys
from sys import stdout

import re
import glob

import numpy

import esutil
from esutil.ostools import path_join

import sdsspy


def read(ftype, run=None, camcol=None, field=None, id=None, **keys):
    """
    Module:
        sdsspy.files
    Name:
        read
    Purpose:
        Read an SDSS file given the input file type and id info
    Calling Sequence:
        data=read(ftype, run=None, camcol=None, field=None, id=None, **keys)

    Example:
        data=read('psField', 756, 3, 125)
        imdict=read('fpAtlas', 756, 3, 125, 125, trim=True)


    Inputs:
        ftype: 
            A file type, such as 'psField', 'fpAtlas', 'calibObj.gal'.  See
            $SDSSPY_DIR/share/sdssFileTypes.par for a list of types.

            The ftype is case-insensitive.

    Keyword Inputs:

        NOTE: For run,camcol,field, *one and only one* of these can be a
        sequence.  You can also send globs such as '*' for these.

        run: 
            SDSS run id, needed for most files.
        camcol:
            SDSS camera column id
        field:
            SDSS field identifier
        id:
            SDSS id within a field.  E.g. the atlas reader requires an id.
        filter:
            SDSS filter id, such as 'u','g','r','i','z'.  You can also send
            an index [0,5] representing those.
        rerun:
            You normally don't have to specify the rerun, it can typically
            be determined from the run list.

        rows:
            A subset of the rows to return for tables.  Ignored if reading
            multiple files.
        columns: 
            A subset of columns to read.

        combine:
            When reading multiple files, the result is usually a list of
            results, one per file.  If combine=True, these are combined
            into a single structure.
        ensure_native:
            Ensure the data is in native byte order.
        lower, upper:
            If lower, make sure all field names are in lower case. Upper
            is the opposite.  Useful for reading outputs from IDL mwrfits
            which is all caps.
        trim:
            Trim atlas images.  See sdsspy.atlas.read_atlas for more info.
        ext:
            Which extension to read.  Some files, like psField files, are
            multi-extension.

        verbose:  
            Print the names of files being read.
    Other Keywords:
        These will probably be specific to a given file type. As an example:
        id=: SDSS id within a field.  Some user-defined file type may depend
            on the id

    """

    flist = file_list(ftype, run, camcol, field, id=id, **keys)

    lftype = ftype.lower()
    if lftype == 'fpatlas':
        return _read_atlas(flist, id=id, **keys)
    else:
        if is_yanny(flist[0]):
            return _read_yanny(flist, **keys)

        import esutil as eu

        if len(flist) == 1:
            flist=flist[0]

        ext=keys.get('ext', None)
        if ext is None:
            fs = filespec(ftype)
            if fs['ext'] > -1:
                ext=fs['ext']
                keys['ext'] = ext

        if lftype == 'psField' and ext != 6:
            return _read_psfield(flist, ext, **keys)

        return eu.io.read(flist, **keys)

def filename(ftype, run=None, camcol=None, field=None, **keys):
    """
    Module:
        sdsspy.files
    Name:
        filename
    Purpose:
        Generate an SDSS file name from input file type and id info
    Calling Sequence:
        fname=filename(ftype, run=None, camcol=None, field=None, **keys)

    Example:
        filename('psField', 756, 3, 125)


    Inputs:
        ftype: 
            A file type, such as 'psField', 'fpAtlas', 'calibObj.gal'.  See
            $SDSSPY_DIR/share/sdssFileTypes.par for a list of types.  
            
            The ftype is case-insensitive.

    Keyword Inputs:

        NOTE: run,camcol,field can also be globs such as '*'

        run: 
            SDSS run id, needed for most files.
        camcol:
            SDSS camera column id

        field:
            SDSS field identifier
        filter:
            SDSS filter id, such as 'u','g','r','i','z'.  You can also send
            an index [0,5] representing those.

        rerun:
            You normally don't have to specify the rerun, it can typically
            be determined from the run list.

    Other Keywords:
        These will probably be specific to a given file type. As an example:
        id=: SDSS id within a field.  Some user-defined file type may depend
            on the id

    """
    fs=FileSpec()
    return fs.filename(ftype, run, camcol, field, **keys)

def filedir(ftype, run=None, camcol=None, **keys):
    """
    Module:
        sdsspy.files
    Name:
        filedir
    Purpose:
        Generate an SDSS dir path from input file type and id info
    Calling Sequence:
        fname=filedir(ftype, run=None, camcol=None, **keys)

    Example:
        filedir('psField', 756, 3)


    Inputs:
        ftype: 
            A file type, such as 'psField', 'fpAtlas', 'calibObj.gal'.  See
            $SDSSPY_DIR/share/sdssFileTypes.par for a list of types.  
            
            The ftype is case-insensitive.

    Keyword Inputs:
        run: 
            SDSS run id, needed for most files.
        camcol:
            SDSS camera column id

        rerun:
            You normally don't have to specify the rerun, it can typically
            be determined from the run list.

    Other Keywords:
        These will probably be specific to a given file type

    """

    fs=FileSpec()
    return fs.dir(ftype, run, camcol, **keys)
    
def filespec(ftype):
    """
    Module:
        sdsspy.files
    Name:
        filespec
    Purpose:
        Return the path specificition for the in input file type.  This information
        read from
            $SDSSPY_DIR/share/sdssFileTypes.par
    Inputs:
        ftype:  The file type. The ftype is case-insensitive.

    Output:
        A dictionary with the file specification:
            'dir': A directory pattern, e.g. '$PHOTO_REDUX/$RERUN/$RUNNUM/objcs/$COL'
            'name': A name pattern, e.g. 'fpAtlas-$RUNSTR-$COL-$FIELDSTR.fit'
            'ftype': The file type, e.g. fpAtlas
            'ext': The extension.  Supported for fits files.
    """
    fs = FileSpec()
    return fs.filespec(ftype)


def is_yanny(fname):
    return fname[-4:] == '.par'

def _read_yanny(fname, **keys):
    if is_sequence(fname):
        if len(fname) > 1:
            raise ValueError("Only read one .par file at a time")
        fname=fname[0]
    verbose = keys.get('verbose',False)
    if verbose:
        stdout.write("Reading file: '%s'\n" % fname)
        
    if 'runList' in fname or 'sdssFileTypes' in fname:
        return sdsspy.yanny.readone(fname)
    elif 'sdssMaskbits' in fname:
        return sdsspy.yanny.read(fname)
    else:
        raise ValueError("Unexpected yanny .par file: '%s'" % fname)


def _read_psfield(fname, ext, **keys):
    """

    This reader is designed to read extensions from a single file.  If you want
    to read multiple files extension 6 then just use the main read() function.

    """
    
    if is_sequence(ext):
        n_ext=len(ext)
        if n_ext == 1:
            ext=ext[0]
    else:
        n_ext=1

    if n_ext == 1:
        data = eu.io.read(fname, ext=ext, **keys)
    else:
        data = []
        for text in ext:
            tdata = eu.io.read(fname, ext=text, **keys)
            data.append(tdata)

    return data


def _read_atlas(flist, **keys):
    import atlas
    trim = keys.get('trim',False)
    id = keys.get('id', None)
    verbose=keys.get('verbose', False)

    if id is None:
        raise ValueError("You must enter id(s) to read fpAtlas files")

    if len(flist) > 1:
        raise ValueError("You can only read from one atlas image at a "
                         "time, %s requested" % len(flist))
    fname=flist[0]

    if is_sequence(id):
        if len(id) > 1:
            raise ValueError("You can only read one object at a time "
                             "from atlas images, %s requested" % len(id))
        id=id[0]

    if verbose:
        print "Reading id %s from '%s'" % (id,fname)
    imdict=atlas.read_atlas(fname, id, trim=trim)
    return imdict

class FileSpec:
    def __init__(self, reload=False):
        self.load(reload=reload)

    def load(self, reload=False):
        if not hasattr(FileSpec, '_filetypes') or reload:
            if 'SDSSPY_DIR' not in os.environ:
                raise ValueError('SDSSPY_DIR environment var must be set')
            d=os.environ['SDSSPY_DIR']
            d = os.path.join(d, 'share')
            f = os.path.join(d, 'sdssFileTypes.par')

            self._filetypes = sdsspy.yanny.readone(f)

            self._ftypes_lower = self._filetypes['ftype'].copy()
            for i in xrange(self._ftypes_lower.size):
                self._ftypes_lower[i] = self._ftypes_lower[i].lower()
    def reload(self):
        self.load(reload=True)

    def filespec(self, ftype):
        w,=numpy.where(self._ftypes_lower == ftype.lower())
        if w.size == 0:
            raise ValueError("File type '%s' is unknown" % ftype)

        fs = {'ftype': str(self._filetypes['ftype'][w][0]),
              'dir':   str(self._filetypes['dir'][w][0]),
              'name':  str(self._filetypes['name'][w][0]),
              'ext':   int(self._filetypes['ext'][w][0])}
        return fs

    def dir_pattern(self, ftype):
        fs = self.filespec(ftype)
        return fs['dir']
    def name_pattern(self, ftype):
        fs = self.filespec(ftype)
        return fs['name']
    def file_pattern(self, ftype):
        fs = self.filespec(ftype)
        d = fs['dir']
        name=fs['name']
        path = os.path.join(d,name)
        return path

    def dir(self, ftype, run=None, camcol=None, **keys):
        p = self.dir_pattern(ftype)
        d = expand_sdssvars(p, run=run, camcol=camcol, **keys)
        return d
    def name(self, ftype, run=None, camcol=None, field=None, **keys):
        n = self.name_pattern(ftype)
        n = expand_sdssvars(n, run=run, camcol=camcol, field=field, **keys)
        return n
    def filename(self, ftype, run=None, camcol=None, field=None, **keys):
        f = self.file_pattern(ftype)
        f = expand_sdssvars(f, run=run, camcol=camcol, field=field, **keys)
        return f

_file_list_cache={}
def file_list(ftype, run=None, camcol=None, field=None, **keys):
    """
    Module:
        sdsspy.files
    Name:
        file_list
    Purpose:
        Get lists of SDSS files based on id info.
    Calling Sequence:
        fl = file_list(ftype, run=None, camcol=None, field=None, **keys)

    Examples:
        # get atlas file for a given run, camcol, and field
        f=file_list('fpAtlas', 756, 3, 125)
        # get for a range of field
        f=file_list('fpAtlas', 756, 3, range(125,135))
        # get for all fields in the column
        f=file_list('fpAtlas', 756, 3, '*')
        # get for a given field and all columns
        f=file_list('fpAtlas', 756, '*', 125)

    Inputs:
        ftype: 
            A file type, such as 'psField', 'fpAtlas', 'calibObj.gal'.  See
            $SDSSPY_DIR/share/sdssFileTypes.par for a list of types.

            The ftype is case-insensitive.

    Keyword Inputs:

        NOTE: For run,camcol,field, *one and only one* of these can be a
        sequence.  You can also send globs such as '*' for these.

        run: 
            SDSS run id, needed for most files.
        camcol:
            SDSS camera column id
        field:
            SDSS field identifier
        id:
            SDSS id within a field.  E.g. the atlas reader requires an id.
        filter:
            SDSS filter id, such as 'u','g','r','i','z'.  You can also send
            an index [0,5] representing those.
        rerun:
            You normally don't have to specify the rerun, it can typically
            be determined from the run list.

    Output:
        A list of files.  The output is always a list even if only one file
        is returned.

    """
    glob_pattern = keys.get('glob', None)
    if glob_pattern is None:
        if ftype is None:
            raise ValueError('send filetype and some id info or the full pattern on glob= keyword')

        
        run_is_sequence, camcol_is_sequence, field_is_sequence =_check_id_sequences(run,camcol,field)

        if run_is_sequence:
            flist=[]
            for trun in run:
                flist += file_list(ftype, trun, camcol, field, **keys)
            return flist

        if camcol_is_sequence:
            flist=[]
            for tcamcol in camcol:
                flist += file_list(ftype, run, tcamcol, field, **keys)
            return flist

        if field_is_sequence:
            flist=[]
            for tfield in field:
                flist += file_list(ftype, run, camcol, tfield, **keys)
            return flist

        fs=FileSpec()
        glob_pattern = fs.filename(ftype, run=run, camcol=camcol, field=field, **keys)

    # add * at the end to catch .gz files
    glob_pattern += '*'

    if glob_pattern not in _file_list_cache:
        flist = glob.glob(glob_pattern)
        _file_list_cache[glob_pattern] = flist
    else:
        flist = _file_list_cache[glob_pattern]
    return flist

def runlist():
    """
    Module:
        sdsspy.files
    Name:
        runlist
    Purpose:
        Return the sdss runlist.  This is read from $PHOTO_REDUX/runList.par
        This data is cached and only read once, so is preferable to using
        sdsspy.files.read('runList')

    Calling Sequence:
        rl = sdsspy.files.runlist()

    Outputs:
        A numpy array with fields. This is the typedef from the yanny .par file:
            typedef struct {
                int       run;        # Run number
                char      rerun[];    # Relevant rerun directory
                int       exist;      # Does a directory exist (has it ever been queued) (0/1)?
                int       done;       # Is this field done (0-Not done/1-Done)?
                int       calib;      # Is this field calibrated (0-No/1-Yes)?
                                      # All the fields below are only updated if the run is done
                                      # Otherwise, they are set to zero.
                int       startfield; # Starting field - determined from fpObjc present
                int       endfield;   # Ending field
                char      machine[];  # The machine the data is on
                char      disk[];     # The disk this run is on....
            } RUNDATA;

    """
    rl = RunList()
    return rl.data


class SimpleFileCache:
    """
    This is a base class for file caches.  Inherit from this and implement a
    .load() method.
    """
    def __init__(self, reload=False):
        self.load(reload=reload)

    def load(self, reload=False):
        """
        Over-ride this method, and replace SimpleCache below with your
        class name.
        """
        if not hassattr(SimpleFileCache, 'data') or reload:
            SimpleFileCache.data=None

    def reload(self):
        self.load(reload=True)

class RunList(SimpleFileCache):
    """
    This is a cache for the runList
    """
    def load(self, reload=False):
        if not hasattr(RunList, 'data') or reload:
            RunList.data = read('runList')


def _check_id_sequences(run,camcol,field):
    run_is_sequence=False
    camcol_is_sequence=False
    field_is_sequence=False

    nseq=0
    if is_sequence(run):
        nseq += 1
        run_is_sequence=True
    if is_sequence(camcol):
        nseq += 1
        camcol_is_sequence=True
    if is_sequence(field):
        nseq += 1
        field_is_sequence=True

    # only one of run/camcol/field can be a sequence.
    if nseq > 1:
        raise ValueError("only one of run/camcol/field can be a sequence")

    return run_is_sequence, camcol_is_sequence, field_is_sequence 

def is_sequence(var):
    if isinstance(var, (str, unicode)):
        return False

    try:
        l = len(var)
        return True
    except:
        return False


def find_rerun(run):
    rl = runlist()
    w,=numpy.where(rl['run'] == run)
    if w.size == 0:
        raise ValueError("Run %s not found in runList.par" % run)
    return rl['rerun'][w[0]]

# convert sdss numbers to strings in file names and such
def stripe2string(stripes):
    """
    ss = stripe2String(stripes)
    Return the string version of the stripe.  9->'09'
    Range checking is applied.
    """
    return tostring(stripes, 0, 99)

def run2string(runs):
    """
    rs = run2string(runs)
    Return the string version of the run.  756->'000756'
    Range checking is applied.
    """
    return tostring(runs,0,999999)


def rerun2string(reruns):
    """
    rrs = rerun2string(reruns)
    Return the string version of the rerun.  No zfill is used.
    """
    return tostring(reruns)

def camcol2string(camcols):
    """
    cs = camcol2string(camcols)
    Return the string version of the camcol.  1 -> '1'
    Range checking is applied.
    """
    return tostring(camcols,1,6)
    
def field2string(fields):
    """
    fs = field2string(field)
    Return the string version of the field.  25->'0025'
    Range checking is applied.
    """
    return tostring(fields,0,9999)

def id2string(ids):
    """
    istr = id2string(ids)
    Return the string version of the id.  25->'00025'
    Range checking is applied.
    """
    return tostring(ids,0,99999)



filter_dict = {0: 'u',
             1: 'g',
             2: 'r',
             3: 'i',
             4: 'z',
             'u':'u',
             'g':'g',
             'r':'r',
             'i':'i',
             'z':'z'}

def filter2string(filter):
    """
    fstr = filter2string(filter)
    Return alphabetic representation of filter
      fstr = filter2string(2)      # returns 'r'
      fstr = filter2string('2')    # returns 'r'
      fstr = filter2string('r')    # returns 'r'
      fstr = filter2string('R')    # returns 'r'
    """

    if not numpy.isscalar(filter):
        # Scalar pars cannot be modified
        return [filter2string(f) for f in filter]

    if filter == '*':
        return '*'

    if filter not in filter_dict:
        raise ValueError("bad filter indicator: %s" % filter)

    return filter_dict[filter]

def tostring(val, nmin=None, nmax=None):
    
    if not numpy.isscalar(val):
        return [tostring(v,nmin,nmax) for v in val]

    if isinstance(val, (str,unicode)):
        return val

    if nmin is not None:
        if val < nmin:
            raise ValueError("Number ranges below min value of %s\n" % nmin)
    if nmax is not None:
        if val > nmax:
            raise ValueError("Number ranges higher than max value of %s\n" % nmax)


    if nmax is not None:
        nlen = len(str(nmax))
        vstr = str(val).zfill(nlen)
    else:
        vstr = str(val)

    return vstr






def expand_sdssvars(string_in, **keys):

    string = string_in

    # this will expand all environment variables, e.g. $PHOTO_SWEEP
    # if they don't exist, the result will be incomplete

    string = os.path.expandvars(string)

    if string.find('$RUNNUM') != -1:
        run=keys.get('run', None)
        if run is None:
            raise ValueError("run keyword must be sent: '%s'" % string)
        string = string.replace('$RUNNUM', tostring(run))
    if string.find('$RUNSTR') != -1:
        run=keys.get('run', None)
        if run is None:
            raise ValueError("run keyword must be sent: '%s'" % string)
        string = string.replace('$RUNSTR', run2string(run))

    if string.find('$RERUN') != -1:
        rerun=keys.get('rerun', None)
        if rerun is None:
            # try to determine from run number
            run = keys.get('run',None)
            if run is not None:
                try:
                    rerun = find_rerun(run)
                except:
                    rerun = None
        if rerun is None:
            # try to determine rerun
            raise ValueError("rerun keyword must be sent: '%s'" % string)
        string = string.replace('$RERUN', tostring(rerun))


    if string.find('$COL') != -1:
        camcol=keys.get('camcol', None)
        if camcol is None:
            raise ValueError("camcol keyword must be sent: '%s'" % string)
        string = string.replace('$COL', camcol2string(camcol))

    if string.find('$FIELDSTR') != -1:
        field=keys.get('field', None)
        if field is None:
            raise ValueError("field keyword must be sent: '%s'" % string)
        string = string.replace('$FIELDSTR', field2string(field))


    if string.find('$IDSTR') != -1:
        id=keys.get('id',None)
        if id is None:
            raise ValueError("id keyword must be sent: '%s'" % string)
        string = string.replace('$IDSTR', id2string(id))

    if string.find('$ID') != -1:
        id=keys.get('id',None)
        if id is None:
            raise ValueError("id keyword must be sent: '%s'" % string)
        string = string.replace('$ID', tostring(id))

    if string.find('$FILTER') != -1:
        fiter=keys.get('filter',None)
        if filter is None:
            raise ValueError("filter keyword must be sent: '%s'" % string)
        string = string.replace('$FILTER', filter2string(band))

    if string.find('$TYPE') != -1:
        type=keys.get('type',None)
        if type is None:
            raise ValueError("type keyword must be sent: '%s'" % string)
        string = string.replace('$TYPE', tostring(type))

    # see if there are any leftover un-expanded variables.  If so
    # raise an exception
    if string.find('$') != -1:
        raise ValueError("There were unexpanded variables: '%s'" % string)

    return string

