""""
MaterialsCoord: NearNeighbor benchmarking.
"""
from setuptools import setup, find_packages


with open('README.md', 'r') as file:
    long_description = file.read()

setup(
    name='materialscoord',
    version='0.1.0',
    url='https://github.com/hillarypan/MaterialsCoord',
    license='A modified BSD license',
    author='Hillary Pan',
    author_email='hp393@cornell.edu',
    description='Materials coordination benchmarking',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Information Technology',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering',
        'Topic :: Other/Nonlisted Topic',
        'Operating System :: OS Independent',
        ],
    keywords='crystal-structure crystallography benchmark',
    test_suite='nose.collector',
    packages=find_packages(),
    install_requires=['spglib', 'numpy', 'pymatgen>=2017.12.30',
                      'monty', 'pandas'],
    extras_require={'docs': ['sphinx', 'sphinx-argparse',
                             'sphinx-autodoc-typehints', 'm2r'],
                    'tests': ['nose', 'coverage', 'coveralls']},
    package_data={'materialscoord': ["structures/*/*"]},
    data_files=['LICENSE'],
    )
