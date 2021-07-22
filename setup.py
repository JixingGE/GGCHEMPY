from setuptools import setup, find_packages
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
setup(
    name='ggchempy',
    version='1.0',
    description='Gas-Grain CHEMical (GGCHEM) code for interstelar clouds',
    long_description=str(open(path.join(here, "longdescription.txt")).read()),
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
    py_modules=["__init__","ggchempylib","ggfuncs"],
    #data_files=[('in', ['in/network2.txt', 'in/ed.txt', 'in/iabun.txt'])],
    install_requires=['matplotlib','numpy','scipy','progressbar','numba']
)
