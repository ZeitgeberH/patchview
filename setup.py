#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import find_packages, setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = []
with open('requirements.txt', 'r') as f:
   requirements += [l for l in f.read().split('\n') if l.strip()]

## forked version of NeuroM   
requirements +=[
    'NeuroM @ git+https://github.com/ZeitgeberH/NeuroM@patchview#egg=NeuroM'
]
setup_requirements = ['pytest-runner']

test_requirements = ['pytest']

#with open('requirements_dev.txt', 'r') as f:
#    dev_requirements = [l for l in f.read().split('\n') if l.strip()]
dev_requirements = [] 

setup(
    author="Ming Hu",
    author_email='',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: patch-clamp users',
        'License :: OSI Approved :: BSD 3-Clause License',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "Operating System :: OS Independent",
    ],
    description="Patchview perform data analysis and visualization on whole-cell recording data, including firing pattern analysis, event analysis, synatpic connection detection, morphorlocial analysis and more.",
    install_requires=requirements,
    extras_require={'dev': dev_requirements},
    license="BSD 3-Clause license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='patchview',
    name='patchview',
    packages=find_packages(include=['patchview','patchview/Data',\
        'patchview/HekaIO','patchview/neurom','patchview/utilitis']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/zeitgeberH/patchview',
    version='0.2.5',
    zip_safe=False,
    entry_points={
        'gui_scripts': [
            'patchview=patchview:main'
        ],
    },
)
