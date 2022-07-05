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

setup_requirements = ['pytest-runner']

test_requirements = ['pytest']

#with open('requirements_dev.txt', 'r') as f:
#    dev_requirements = [l for l in f.read().split('\n') if l.strip()]
dev_requirements = [] 

setup(
    author="Ming Hu",
    author_email='ming.hu@bcm.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
    ],
    description="Patchview perform data analysis and visualization on whole-cell recording data, including firing pattern analysis, event analysis, synatpic connection detection, morphorlocial analysis and more.",
    install_requires=requirements,
    extras_require={'dev': dev_requirements},
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='patchview',
    name='patchview',
    packages=find_packages(include=['src','src.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/zeitgeberH/patchview',
    version='1.0',
    zip_safe=False,
)
