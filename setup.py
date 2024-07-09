from setuptools import setup, find_packages
import os
 
from ex-fem import __version__
 
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
 
def get_version():
    with open("ex-fem/__version__.py", "r") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip("\"'")
    raise RuntimeError("Unable to find version string.")
 
setup(
    name='ex-fem',
    version=get_version(),  
 
    url='https://https://github.com/kinfungchan/ex-fem-user',
    author='Kin Fung Chan',
    author_email='kin.chan@eng.ox.ac.uk',
 
    packages=find_packages(),
 
    setup_requires=[
        'numpy',
    ],
 
)