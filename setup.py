#! /usr/bin/env python3

import glob
import os
import setuptools
import sys

if sys.version_info < (3, 2):
    sys.exit('Minimum supported Python version is 3.2')

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
#with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

# Get the current version
exec(open("src/rase/version.py").read())

rase_pys = list(map(lambda x: x.split('/')[-1].replace(".py", ""), glob.glob('src/rase/rase_*.py')))
rase_pys_strs = list(map(lambda x: '{z}.py = rase.{z}:main'.format(z=x), rase_pys))

#print(['rase = rase.rase:main'] + rase_pys_strs, file=sys.stderr)

setuptools.setup(
    name='rase',
    version=VERSION,
    description='description',
    #long_description=long_description,
    url='https://github.com/karel-brinda/rase',
    author='Karel Brinda',
    author_email='kbrinda@hsph.harvard.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: Unix',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        'wheel',
    ],
    package_data={
        'src': [
            '*.py',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    scripts=glob.glob("scripts/*"),
    entry_points={
        'console_scripts': ['rase = rase.rase:main'] + rase_pys_strs,
    },
)
