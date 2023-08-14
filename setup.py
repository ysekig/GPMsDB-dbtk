# setup.py
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'GPMsDB_dbtk', 'VERSION'), 'r') as f:
        return f.readline().strip()

setup(
      cmdclass    = {'build_ext':build_ext},
      )

setup(
    name='GPMsDB_dbtk',
    python_requires='>=3.6',
    version=version(),
    author='Yuji Sekiguchi',
    author_email='y.sekiguchi@aist.go.jp',
    packages=['GPMsDB_dbtk', 'GPMsDB_dbtk.util'],
    scripts=['bin/GPMsDB_dbtk'],
    package_data={'GPMsDB_dbtk': ['VERSION']},
    url='https://github.com/ysekig/GPMsDB-dbtk',
    mdclass={'build_ext': build_ext},
    description='Toolkit for bacterial and archaeal identification based on MALDI-TOF MS peak lists.',
    install_requires=['biolib>=0.1.0', 'numpy']
)
