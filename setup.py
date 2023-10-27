import platform
from setuptools import setup, Extension
from Cython.Distutils import build_ext

NAME = "evenlyspacedstreamlines"
VERSION = "0.1.0"
DESCR = "Generate evenly-spaced streamlines from an orientation field on a triangulated 3D surface"
KEYWORDS = "vector,field,visualization,surface,streamline"
URL = "http://github.com/jacquemv/evenlyspacedstreamlines"
REQUIRES = ['numpy', 'cython']
AUTHOR = "Vincent Jacquemet"
EMAIL = "vincent.jacquemet@umontreal.ca"
LICENSE = "MIT"
SRC_DIR = "evenlyspacedstreamlines"
PACKAGES = [SRC_DIR]

if platform.system() == 'Windows':
    compiler_args = ['/openmp', '/O2']
    linker_args = []
else:
    compiler_args = ['-fopenmp', '-Ofast']
    linker_args = ['-fopenmp']
    
ext = Extension(SRC_DIR + ".runengine",
                sources=[SRC_DIR + "/engine.cpp", SRC_DIR + "/runengine.pyx"],
                libraries=[],
                extra_compile_args=compiler_args,
                extra_link_args=linker_args,
                language="c++",
                include_dirs=[SRC_DIR])
ext.cython_directives = {'language_level': "3"}

EXTENSIONS = [ext]

setup(install_requires=REQUIRES,
      packages=PACKAGES,
      zip_safe=False,
      name=NAME,
      version=VERSION,
      description=DESCR,
      keywords=KEYWORDS,
      long_description=open('README.md', 'r').read(),
      long_description_content_type='text/markdown',
      author=AUTHOR,
      author_email=EMAIL,
      url=URL,
      license=LICENSE,
      cmdclass={"build_ext": build_ext},
      ext_modules=EXTENSIONS
)
