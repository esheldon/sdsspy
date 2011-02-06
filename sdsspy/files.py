"""
    Module:
        sdsspy.files

    Description:
        A set of functions for dealing with SDSS files, including reading files,
        finding file lists, getting run metadata, converting SDSS id numbers 
        to strings, etc.

    Classes:
        FileList:
            A class for genering file lists.  This reads actual files from
            disk.  If you want to do more complex selections, such as using the
            resolve and calibration information, see the sdsspy.window module.


    Functions. 
            read(type, 
                 run=, 
                 rerun=301, 
                 camcol='*', 
                 field='*', 
                 band='*', 
                 type=,
                 nocombine=)

                Read from SDSS files.

            filename(type, 
                     run=, 
                     rerun=301, 
                     camcol='*', 
                     field='*', 
                     band='*', 
                     add_dir=)

                Generate SDSS file names from id info.

            filedir(subdir, 
                    run=, 
                    rerun=301, 
                    camcol='*', 
                    field='*',
                    band='*')

                    
                Generate an SDSS directory name from id info.

        expand_sdssvars:
            Convert environment variables and sdss vars, e.g. $RUNNUM, in file
            specifications to values.  Used to generate file locations.

        Routines for converting SDSS id numbers to strings:
            These functions return scalars for scalar input, otherwise lists. 

            int2string(numbers, min, max)
                Convert numbers to strings with zfill between min and max
            stripe2string(stripes)
                Convert input stripe(s) to string(s)
            run2string(runs)
                Convert input run(s) to string(s)
            rerun2string(reruns)
                Convert input rerun(s) to string(s)
            camcol2string(camcols, band=None)
                Convert input camcol(s) to string(s). If band is sent it
                is incorporated into the string as in SDSS files.
            field2string(fields)
                Convert input field(s) to string(s)
            id2string(ids)
                Convert input id(s) to string(s)
            filter2string(bands)
                Return alphabetic representation of band; works on both string
                and numerical input.

    Dependencies:
        numpy
        esutil


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

_DEFAULT_RERUN=301

debug=False


def filespec(ftype):
    fs = FileSpec()
    return fs.filespec(ftype)

def filename(ftype, run=None, camcol=None, field=None, **keys):
    fs=FileSpec()
    return fs.filename(ftype, run, camcol, field, **keys)

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
def file_list(ftype=None, run=None, camcol=None, field=None, **keys):

    glob_pattern = keys.get('glob', None)
    if glob_pattern is None:
        if ftype is None:
            raise ValueError('send filetype and some id info or the full pattern on glob= keyword')

        
        run_is_sequence, camcol_is_sequence, field_is_sequence =_check_sequences(run,camcol,field)

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

    if glob_pattern not in _file_list_cache or reload:
        flist = glob.glob(glob_pattern)
        _file_list_cache[glob_pattern] = flist
    else:
        flist = _file_list_cache[glob_pattern]
    return flist

def runlist():
    rl = RunList()
    return rl.runlist()

class RunList:
    """
    This is basically a cache for the runList
    """
    def __init__(self, reload=False):
        self.load(reload=reload)
    def load(self, reload=False):
        if not hasattr(RunList, '_runlist') or reload:
            f=file_list('runList')
            if len(f) == 0:
                fname = filename('runList')
                raise ValueError("file not found: '%s'" % fname)
            self._runlist = sdsspy.yanny.readone(f[0], defchar=10)

    def reload(self):
        self.load(reload=True)

    def runlist(self):
        return self._runlist

def _check_sequences(run,camcol,field):
    run_is_sequence=False
    camcol_is_sequence=False
    field_is_sequence=False

    nseq=0
    if _is_sequence(run):
        nseq += 1
        run_is_sequence=True
    if _is_sequence(camcol):
        nseq += 1
        camcol_is_sequence=True
    if _is_sequence(field):
        nseq += 1
        field_is_sequence=True

    # only one of run/camcol/field can be a sequence.
    if nseq > 1:
        raise ValueError("only one of run/camcol/field can be a sequence")

    return run_is_sequence, camcol_is_sequence, field_is_sequence 

def _is_sequence(var):
    try:
        l = len(var)
        return True
    except:
        return False


def read(filetype, 
         run=None, 
         rerun=_DEFAULT_RERUN, 
         camcol=None,
         field='*', 
         **keywords):
    """
    res = sdsspy.files.read(filetype, 
                            run=None, 
                            rerun=301,
                            camcol=None, 
                            field='*',
                            columns=None,
                            combine=True, 
                            verbose=False, 
                            header=False, 
                            ensure_native=False,
                            **file_keywords)

    Inputs:
        filetype: File type, e.g. 'calibobj'
    Keywords:
        run: 
            The SDSS run.  Most files require this.
        rerun: 
            The SDSS rerun.  Default is 301
        camcol: 
            Camera column.  Required for many SDSS files.  Can be '*'
        field: 
            The field to read.  Can be '*' or a sequence or array.  Default is
            '*' meaning "read all fields"
        band: 
            bandpass. Required for some files.  Can also be '*'
        type: 
            The type, e.g. 'gal' for calibObj.  Can also be '*'
        combine: 
            Results are initially kept in a list of arrays, with each list
            element corresponding to a single file read.  By default these are
            combined into one big array.  If this keyword is False this
            combination step is not performed, saving some memory but making
            use difficult.

        rows: 
            A subset of rows to read from a file.  Ignored if reading multiple
            files.
        columns: 
            A subset of columns to read.
        ensure_native:
            For some file formats, e.g. .fits and .rec, you can force the 
            byte ordering to native by setting this keyword.
        verbose:  
            Print the names of files being read.
    Revision History:
        Documented: 2007-06-01, Erin Sheldon, NYU
    """

    # get some keywords
    rows=keywords.get('rows',None)
    columns=keywords.get('columns',None)
    verbose=keywords.get('verbose',False)
    combine=keywords.get('combine',True)
    ensure_native=keywords.get('ensure_native',False)
    
    # Note, we always get the full field, if this filetype is split by fields.
    # then we can extract a subset later

    ftype=filetype.lower()
    fl = FileList(ftype, run, rerun=rerun, camcol=camcol, **keywords)
    flist = fl.flist(subfields=field)
    
    data = esutil.io.read(flist, 
                          ext=filespec[ftype]['hdu'],
                          rows=rows,
                          columns=columns, 
                          view=numpy.ndarray, 
                          verbose=verbose, 
                          lower=True,
                          ensure_native=ensure_native,
                          combine=combine)

    # for some types we can split by field after the fact,e.g. calibobj
    # only try this if the data on disk are not split by field already.
    post_fieldsplit_types = ['calibobj']
    if field != '*' \
            and not filespec[ftype]['byfield'] \
            and ftype in post_fieldsplit_types:

        data = extract_subfields(data, field, verbose=verbose)

    return data


def extract_subfields(data, fields, verbose=False):
    if 'field' in data.dtype.names:
        if verbose:
            stdout.write("Extracting a subset of fields\n")
        h,rev=esutil.stat.histogram(data['field'],min=0,rev=True)

        f2get = ensure_sequence(fields)
        keep = numpy.zeros(data.size)
        for f in f2get:
            if rev[f] != rev[f+1]:
                w=rev[ rev[f]:rev[f+1] ]
                keep[w] = 1
        wkeep, = numpy.where(keep != 0)
        if wkeep.size == 0:
            raise ValueError("None of the requested fields matched")

        data_keep = data[wkeep]
    else:
        data_keep = data

    return data_keep



def filedir(filetype, run=None, rerun=_DEFAULT_RERUN, camcol='*', **keywords):
    """
    dir = sdsspy.files.filedir(filetype, run=None, **keywords)

    E.g. For calibobj, all that is required is the filetype
    and rerun.

       filedir = filedir('calibobj', rerun=301)

    Note also, the rerun can be omitted if using the default
    which is 301

       filedir = filedir('calibobj')

    Keywords:
        run
        rerun
        camcol
        field
        id
        band
        type

    """

    ftype=filetype.lower()
    if ftype not in filespec:
        raise ValueError("Unsupported file type: '%s'" % ftype)

    dir = expand_sdssvars(filespec[ftype]['dir'], 
                          run=run, 
                          rerun=rerun, 
                          camcol=camcol,
                          **keywords)
    return dir




def filename_old(filetype_input, 
             run=None,
             rerun=_DEFAULT_RERUN,
             camcol='*', 
             field='*',
             add_dir=True, 
             **keywords):
    """
    name = sdsspy.files.filename(filetype, run, rerun=301, add_dir=True, **keywords)

    Keywords:
        run (can either be argument or keyword)
        rerun (default 301)
        camcol
        field
        id
        band
        type

    Get an SDSS file name, such as calibobj, etc, given the relevant info.
    run,rerun,camcol,band must be scalars.  
    
    Certain file names requre more info than others. asTrans only requires the
    run, but calibObj needs run,camcol (assuming default rerun). The directory
    is also added unless add_dir = False


    # example: calibobj for run 756, camcol 1 type 'gal' and default
    #     rerun of 301
    >>> filename('calibobj', 756, camcol=1, type='gal')
    '/global/data/boss/sweeps/2010-01-11/301/calibObj-000756-1-gal.fits.gz'
    """

    ftype=filetype_input.lower()

    if ftype not in filespec:
        raise ValueError("Unsupported file type: '%s'" % ftype)

    name = filespec[ftype]['name']
    name = expand_sdssvars(name, 
                           run=run, 
                           rerun=rerun, 
                           camcol=camcol,
                           field=field,
                           **keywords)

    if add_dir:
        dir = filedir(ftype, 
                      run=run, 
                      rerun=rerun, 
                      camcol=camcol,
                      field=field, 
                      **keywords)
        name = os.path.join(dir, name)
    
    return name



                 

__files2fields_pattern_compiled = re.compile('-[0-9][0-9][0-9][0-9]$')
__files2fields_errmess = \
      'Incorrectly formatted name %s.  Must end in -iiii.extension\n'
def files2fields(filenames):
    """
    field = sdsspy.files.files2fields(filenames)
    Extract the field number from field-by-field type of files such
    as a tsObj file that end with front-iiii.extension
    """
    if numpy.isscalar(filenames):
        filenames = [filenames]

    fields = []
    search = __files2fields_pattern_compiled.search

    for name in filenames:
        if name.find('.') != -1:
            sp = name.split('.')[0]
            mo = search(sp)
            if mo:
                field = int( sp[mo.start()+1:] )
                fields.append(field)
            else:
                sys.stdout.write(__files2fields_errmess % s)
                return []
        else:
            sys.stdout.write(__files2fields_errmess % s)
            return []

    return fields


class FileListOld:
    """
    Class:
        FileList
    Purpose:
        A class to find sdss file lists on disk.  Lists for previous calls
        are cached internally for speed.

    Construction and loading of lists:
        By default only the file type is needed, in which case a glob pattern
        with '*' for all file elements is created.  Subsets are determined
        through the keyword.

    Required Arguments for Construction or load():
        filetype:  An SDSS filetype.  Supported types:
            calibObj

    Keywords for Construction and load():
        run
        rerun (default 301)
        camcol
        field
        id
        band
        type (e.g. 'gal' for calibObj)
        
    Examples:
        import sdsspy
        ft=sdsspy.files.FileList('calibobj', run=756, type='gal')
        flist = ft.flist()
        for f in flist:
            print f

        ft.load('fpAtlas', run=94, camcol=3)
        for f in ft.flist():
            print f
    """
    def __init__(self, 
                 filetype, 
                 run=None, 
                 rerun=_DEFAULT_RERUN, 
                 camcol='*', 
                 field='*', 
                 id='*', 
                 band='*',
                 type='*',
                 **keywords):

        if not hasattr(FileList, '_flist'):
            # this class variable will hold all previously requested file 
            # lists
            FileList._flist={} 

        self.load(filetype, run, rerun, camcol, field, id, band, type)

    def load(self, 
             filetype, 
             run=None, 
             rerun=301, 
             camcol='*', 
             field='*', 
             id='*', 
             band='*', 
             type='*'):

        self._filetype=filetype.lower()
        self._run=run
        self._rerun=rerun
        self._camcol=camcol
        self._field=field
        self._id=id
        self._band=band
        self._type=type
        self._key = self.makekey(filetype,run,rerun,camcol,field,id,band,type)

        if self._key not in FileList._flist:
            
            # now load the file list
            pattern = filename(filetype,
                               run=run,rerun=rerun,camcol=camcol,field=field,id=id,
                               band=band,type=type)
            if filespec[self._filetype]['check_compress']:
                pattern += '*'

            flist = glob.glob(pattern)
            if len(flist) == 0:
                raise RuntimeError("no files matched pattern: %s" % pattern)
            FileList._flist[self._key] = flist

            self._used_cache = False
        else:
            self._used_cache = True

    def used_cache(self):
        return self._used_cache

    def reload(self):
        if self._key in FileList._flist:
            del FileList._flist[self._key]
        self.load(self._filetype,
                  self._run,
                  self._rerun,
                  self._camcol,
                  self._field,
                  self._id,
                  self._band,
                  self._type)

    def flist(self, subfields='*'):
        if self._key not in FileList._flist:
            raise ValueError("File has not been loaded")

        if subfields == '*' or not filespec[self._filetype]['byfield']:
            # just return everything if '*' or if the data is not split
            # by field anyway
            return FileList._flist[self._key]
        else:
            # if this file type is split by fields, extract the desired subset
            if isinstance(subfields,(str,unicode)):
                if subfields != '*':
                    raise ValueError("String not supported for field keyword "
                                     "except for '*'")
            else:
                if isinstance(subfields,(list,tuple,numpy.ndarray)):
                    f2get = subfields
                else:
                    f2get = [subfields]

                keep = []
                this_flist=FileList._flist[self._key]
                for f in f2get:
                    tname = filename(self._filetype,self._run,
                                     rerun=self._rerun,
                                     camcol=self._camcol,
                                     field=f,
                                     id=self._id,
                                     band=self._band,
                                     type=self._type)

                    if tname in this_flist:
                        keep.append(tname)

                if len(keep) == 0:
                    raise ValueError("None of the requested fields were found")

                return keep 



    def clear(self):
        FileList._flist={}

    def makekey(self,filetype,run,rerun,camcol,field,id,band,type):
        key='%s-%s-%s-%s-%s-%s-%s-%s' % (filetype,run,rerun,camcol,field,id,band,type)
        return key


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

def ensure_sequence(var):
    if isinstance(var, (list,tuple,numpy.ndarray)):
        return var
    else:
        return [var]

