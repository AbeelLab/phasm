import os
from setuptools import setup, find_packages


# Nice hack to retrieve docstring and version without importing the module
__version__ = None
__doc__ = ''
read_doc = 0  # Not started, in progress, done
init_file = os.path.join(os.path.dirname(__file__), 'phasm', '__init__.py')
for line in open(init_file).readlines():
    if line.startswith('version_info') or line.startswith('__version__'):
        exec(line.strip())
    elif line.startswith('"""'):
        if read_doc == 0:
            read_doc = 1
            line = line.lstrip('"')
        elif read_doc == 1:
            read_doc = 2
    if read_doc == 1:
        __doc__ += line


setup(
    name='phasm',
    version=__version__,
    packages=find_packages(),

    author='Lucas van Dijk',
    author_email='info@lucasvandijk.nl',
    url='https://bitbucket.org/tudelft-bioinformatics/phasm',

    description='Haplotype-aware de novo genome assembler for aneuploid '
                'organisms from long read data.',
    long_description=__doc__.strip(),
    license='MIT',

    # Dependencies
    install_requires=[
        'networkx>=1.9',
    ],
    setup_requires=['pytest-runner'],
    test_requires=['pytest', 'pytest-cov']
)
