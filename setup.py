#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup


def read(*names, **kwargs):
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()


# activate once added, do not remove
#long_description = '{}\n{}'.format(
#    read('README.rst'),
#    read(join('docs', 'CHANGELOG.rst')),
#    )


setup(
    name='haddock3',
    version='0.0.3-alpha',
    description='Haddock 3.',
    long_description='',#long_description,
    long_description_content_type='text/x-rst',
    license='Apache License 2.0',
    author='HADDOCK',
    author_email='A.M.J.J.Bonvin@uu.nl',
    url='https://github.com/haddocking/haddock3',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    #py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        #'Development Status :: 4 - Beta',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Microsoft',
        'Programming Language :: Python :: 3.8',
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
    python_requires='>=3.8, <4',
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
