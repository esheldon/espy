import os
import glob
from setuptools import setup, find_packages

scripts = glob.glob('bin/*')
scripts = [s for s in scripts if '~' not in s]

setup(
    name="espy",
    version="0.9.1",
    packages=find_packages(),
    scripts=scripts,
    author='Erin Sheldon',
    author_email='erin.sheldon@gmail.com',
    url='https://github.com/esheldon/espy',
)
