#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
from os.path import dirname, join

from setuptools import find_packages, setup


def read(*names, **kwargs):
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()


# activate once added, do not remove
long_description = '{}\n{}'.format(
    read('README.md'),
    read('CHANGELOG.md'),
    )


setup(
    name='haddock3',
    version='3.0.0',
    description='Haddock 3.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='Apache License 2.0',
    author='HADDOCK',
    author_email='A.M.J.J.Bonvin@uu.nl',
    url='https://github.com/haddocking/haddock3',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    # py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3.9',
        ],
    project_urls={
        'webpage': 'https://github.com/haddocking/haddock3',
        'Documentation': 'https://github.com/haddocking/haddock3#readme',
        'Changelog': '',
        'Issue Tracker': 'https://github.com/haddocking/haddock3/issues',
        'Discussion Forum': 'https://github.com/haddocking/haddock3/issues',
        },
    keywords=[
        'Structural Biology',
        'Biochemistry',
        'Docking',
        'Protein docking',
        'Proteins',
        ],
    python_requires='>=3.9, <3.10',
    install_requires=[
        # not added on purpose
        ],
    extras_require={
        },
    setup_requires=[
        ],
    entry_points={
        'console_scripts': [
            'haddock3 = haddock.clis.cli:maincli',
            "haddock3-mpitask = haddock.clis.cli_mpi:maincli",
            'haddock3-bm = haddock.clis.cli_bm:maincli',
            'haddock3-cfg = haddock.clis.cli_cfg:maincli',
            'haddock3-clean = haddock.clis.cli_clean:maincli',
            'haddock3-copy = haddock.clis.cli_cp:maincli',
            'haddock3-dmn = haddock.clis.cli_dmn:maincli',
            'haddock3-pp = haddock.clis.cli_pp:maincli',
            'haddock3-score = haddock.clis.cli_score:maincli',
            'haddock3-unpack = haddock.clis.cli_unpack:maincli',
            ]
        },
    # cmdclass={'build_ext': optional_build_ext},
    # ext_modules=[
    #    Extension(
    #        splitext(relpath(path, 'src').replace(os.sep, '.'))[0],
    #        sources=[path],
    #        include_dirs=[dirname(path)]
    #    )
    #    for root, _, _ in os.walk('src')
    #    for path in glob(join(root, '*.c'))
    # ],
    )
