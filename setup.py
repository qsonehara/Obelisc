from Cython.Distutils import build_ext
from Cython.Build import cythonize
from setuptools import setup, Extension, find_packages
import numpy

try:
    from Cython.Distutils import build_ext
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

if USE_CYTHON:
    ext = '.pyx'
    cmdclass = {'build_ext': build_ext}
    ext_modules = cythonize([Extension("LOCH_MappingTools", ["obelisc/LOCH_MappingTools.pyx"], include_dirs=[numpy.get_include()])])
else:
    ext = '.cpp'
    cmdclass = {}
    ext_modules = [Extension("LOCH_MappingTools", ["obelisc/LOCH_MappingTools.cpp"], include_dirs=[numpy.get_include()])]

setup(
    name='Obelisc',
    version='0.0.1',
    description='SNP streak-based nonparametric linkage analysis tool',
    long_description="README.md",
    author='Kyuto Sonehara',
    author_email='qsonehara@sg.med.osaka-u.ac.jp',
    url='https://github.com/qsonehara/Obelisc',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=["numpy","pandas","matplotlib","pysnptools","argparse"],
    entry_points={
        "console_scripts": [
            "obelisc = obelisc.main:main"
        ]
    },
    cmdclass=cmdclass,
    ext_modules=ext_modules
)