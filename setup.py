from setuptools import setup, find_packages
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
setup(
    name='ggchem',
    version='1.0',
    description='Gas-Grain CHEMical (GGCHEM) code for interstelar clouds',
    long_description=str(open(path.join(here, "src/longdescription.txt")).read()),
    # The project's main homepage.
    url='will-be-done',
    # Author details
    author='Jixing Ge',
    author_email='gejixing666@gmail.com',
    # Choose your license
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    py_modules=["src/__init__","src/ggchemlib","src/ggchemGUI","src/ggfuncs"],
    install_requires=['matplotlib','numpy','scipy','progressbar','numba']
)
