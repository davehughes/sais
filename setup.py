from setuptools import setup, Extension
import subprocess

# This is a huge hack, but just makes sure the shared library is built during install
status = subprocess.call(' '.split('cd sais-lite && make dll'), shell=True)
if status:
    raise Exception("Failed to build sais-lite DLL.")

setup(name='sais',
      version='0.1.0',
      description='Python wrapper around Yuta Mori\'s implementation of '
                  'SA-IS suffix array construction.'
    )
