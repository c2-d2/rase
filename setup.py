# see https://github.com/pypa/sampleproject

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
exec(open("rase/version.py").read())

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
    packages=["rase"],
    install_requires=[
        'wheel',
    ],
    package_data={
        'rase': [
            '*.py',
        ],
    },
    scripts=glob.glob("scripts/*"),
    entry_points={
        'console_scripts': [
            'rase = rase.rase:main',
        ],
    },
)
