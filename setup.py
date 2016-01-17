from __future__ import print_function
import os

from setuptools import setup, find_packages

readme_dir = os.path.dirname(__file__)
readme_filename = os.path.join(readme_dir, 'README.md')

try:
    with open(readme_filename, 'r') as f:
        readme = f.read()
except:
    readme = ""

try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except:
    print(
        "Conversion of long_description from MD to reStructuredText failed...")


if __name__ == '__main__':
    setup(
        name='pytcga',
        version="0.0.4",
        description="Store and query public TCGA data",
        author="Arun Ahuja",
        author_email="aahuja11 {at} gmail {dot}",
        url="https://github.com/arahuja/pytcga",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'pandas >=0.13.1',
            'nose >=1.3.6',
            'beautifulsoup4',
            'requests',
        ],
        long_description=readme,
        packages=find_packages(exclude=["test", "tests"]),
    )
