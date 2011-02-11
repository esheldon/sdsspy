import distutils.sysconfig
from distutils.core import setup, Extension, Command
import numpy
import os
import sys
import glob


data_files=[]

# data files go under PREFIX/share
data_files = [('share',['sdsspy/share/sdssMaskbits.par','sdsspy/share/sdssFileTypes.par'])]

# if you run python setup.py with_ups sdsspy then the PREFIX/ups dir and table
# file will be written

class AddUPS(Command):
    _data_files = data_files
    user_options=[]
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        # create the ups table

        main_libdir=distutils.sysconfig.get_python_lib()
        pylib_install_subdir = main_libdir.replace(distutils.sysconfig.PREFIX+os.sep,'')

        pylib_install_subdir = pylib_install_subdir.replace('dist-packages','site-packages')

        upstext="""
setupOptional("numpy")
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
        """ % pylib_install_subdir 

        upsdir='ups'
        if not os.path.exists(upsdir):
            sys.stdout.write("Creating ups dir: %s\n" % upsdir)
            os.makedirs(upsdir)

        upsname=os.path.join(upsdir,'sdsspy.table')
        sys.stdout.write('Writing ups table file: %s\n' % upsname)
        upsfile=open(upsname, 'w')
        upsfile.write(upstext)
        upsfile.close()


        AddUPS._data_files.append(('ups',['ups/sdsspy.table']))


# extensions
#
# the atlas reader is in c

#c_files = glob.glob('sdsspy/atlas/*.c')
atlas_sources = glob.glob('sdsspy/atlas/*.c')
#atlas_sources = [f for f in c_files if 'py_psf' not in f]
#psf_sources   = [f for f in c_files if 'py_atlas' not in f]

"""
print 'atlas sources:',atlas_sources
print
print 'psf sources:',psf_sources
"""

atlas_module = Extension('sdsspy.atlas._py_atlas', sources=atlas_sources)
#psf_module = Extension('sdsspy.atlas._py_psf', sources=psf_sources)

# required to compile the atlas reader
os.environ['CFLAGS'] = '-DLINKAGE -DCHECK_LEAKS -DSTAND_ALONE -DSDSS_LITTLE_ENDIAN'

packages = ['sdsspy', 'sdsspy.atlas', 'sdsspy.sandbox']
setup(name='sdsspy',
      version='0.1',
      cmdclass={"with_ups": AddUPS},
      description='A set of python routines for working with SDSS data', 
      packages=packages,
      #ext_modules=[atlas_module, psf_module], 
      ext_modules=[atlas_module], 
      data_files=data_files,
      include_dirs=numpy.get_include())


