import distutils.sysconfig
from distutils.core import setup, Extension, Command
import numpy
import os
import sys
import glob


data_files=[]

data_files = [('share',['share/sdssMaskbits.par','share/sdss_filetypes.par'])]
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


# the atlas reader is in c
atlas_sources=glob.glob('sdsspy/atlas/*.c')
atlas_module = Extension('sdsspy.atlas._py_atlas', sources=atlas_sources)

os.environ['CFLAGS'] = '-DLINKAGE -DCHECK_LEAKS -DSTAND_ALONE -DSDSS_LITTLE_ENDIAN'

packages = ['sdsspy', 'sdsspy.atlas']
setup(name='sdsspy',
      version='0.1',
      cmdclass={"with_ups": AddUPS},
      description='A set of python routines for working with SDSS data', 
      packages=packages,
      ext_modules=[atlas_module], 
      data_files=data_files,
      include_dirs=numpy.get_include())


