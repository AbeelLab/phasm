import os
import sys
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.errors import CompileError


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


__version__ = '0.0.1'


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        'phasm.overlapper',
        [
            'src/overlapper.cpp',
            'src/phasm.cpp'
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            os.path.join(os.path.dirname(__file__), 'vendor')
        ],
        language='c++',
        extra_compile_args=[
            '-fopenmp',
            '-O3',
            '-DNDEBUG',
            '-DSEQAN_HAS_OPENMP=1'
        ],
        extra_link_args=[
            '-fopenmp',
            '-lrt',
            '-lpthread'
        ]
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++14 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' %
                        self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' %
                        self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


setup(
    name='phasm',
    version=__version__,
    packages=find_packages(),

    author='Lucas van Dijk',
    author_email='info@lucasvandijk.nl',
    url='https://github.com/lrvdijk/phasm',

    description='Haplotype-aware de novo genome assembler for aneuploid '
                'and polyploid organisms from long read data.',
    long_description=__doc__.strip(),
    license='MIT',

    # Dependencies
    install_requires=[
        'pybind11>=1.7',
        'networkx>=1.9',
        'numpy>=1.11',
        'scipy>=0.16'
    ],
    setup_requires=['pytest-runner'],
    test_requires=['pytest', 'pytest-cov'],

    # C++ extension
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,

    # Entry points
    entry_points={
        'console_scripts': [
            'phasm = phasm.cli.assembler:main',
            'phasm-convert = phasm.cli.convert:main'
        ]
    }
)
