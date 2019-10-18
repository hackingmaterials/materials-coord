""""
MaterialsCoord: NearNeighbor benchmarking.
"""
from pathlib import Path
from setuptools import setup, find_packages
from materialscoord import __version__ as version

long_description = Path("README.md").read_text()
reqs_raw = Path("requirements.txt").read_text()
reqs_list = [r for r in reqs_raw.split("\n")]

setup(
    name='materialscoord',
    version=version,
    url='https://github.com/hackingmaterials/MaterialsCoord',
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
    packages=find_packages(),
    install_requires=reqs_list,
    extras_require={'tests': ['tests', 'coverage', 'coveralls']},
    package_data={'materialscoord': ["structures/*/*"]},
    data_files=['LICENSE'],
)
