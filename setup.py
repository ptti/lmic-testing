#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='rdt',
      version='0.1',
      description='Modelling of rapid diagnostic testing',
      author=['William Waites'],
      author_email='ww@okapi.cc',
      keywords=['COVID-19'],
      classifiers=[
          # How mature is this project? Common values are
          #   3 - Alpha
          #   4 - Beta
          #   5 - Production/Stable
          'Development Status :: 4 - Beta',

          # Intended audience
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',

          # License
          'License :: OSI Approved :: GNU General Public License (GPL)',
          # Specify the Python versions you support here. In particular,
          # ensure that you indicate whether you support Python 2, Python 3
          # or both.
          'Programming Language :: Python :: 3',
      ],
      license='GPLv3',
      packages=find_packages(),
      install_requires=[
          'matplotlib',
          'pyabc'
      ],
      python_requires='>=3.1.*',
      entry_points={
          'console_scripts': [
              'rdt_sim = rdt.sim:command',
          ],
      },
      package_data={
          'kappa' : [ '*.ka' ],
      }
)
