from setuptools import setup, Extension
from Cython.Distutils import build_ext

NAME = "evenlyspacedstreamlines"
VERSION = "0.1.0"
DESCR = "Generate evenly spaced streamlines from an orientation field on a triangulated 3D surface"
URL = "http://"
REQUIRES = ['numpy', 'cython']

AUTHOR = "Vincent Jacquemet"
EMAIL = "vincent.jacquemet@umontreal.ca"

LICENSE = "MIT"

SRC_DIR = "evenlyspacedstreamlines"
PACKAGES = [SRC_DIR]

ext = Extension("runengine",
                sources=[SRC_DIR + "/engine.cpp", SRC_DIR + "/runengine.pyx"],
                libraries=[],
                extra_compile_args=['-Ofast', '-fopenmp'],
                extra_link_args=['-fopenmp'],
                language="c++",
                include_dirs=[SRC_DIR])

EXTENSIONS = [ext]

setup(install_requires=REQUIRES,
      packages=PACKAGES,
      zip_safe=False,
      name=NAME,
      version=VERSION,
      description=DESCR,
      author=AUTHOR,
      author_email=EMAIL,
      url=URL,
      license=LICENSE,
      cmdclass={"build_ext": build_ext},
      ext_modules=EXTENSIONS
)
