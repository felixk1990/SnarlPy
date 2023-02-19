# @Author:  Felix Kramer
# @Date:   2021-01-14T22:42:23+01:00
# @Email:  kramer@mpi-cbg.de
# @Project: go-with-the-flow
# @Last modified by:    Felix Kramer
# @Last modified time: 2021-11-02T11:26:43+01:00
# @License: MIT


import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "entanglement",
    version = "0.0.1",
    author = "felixk1990",
    author_email = "felixuwekramer@protonmail.com",
    description = "A repository holding the structure for the 'entanglement' package as well example, experiments and galleries.",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/felixk1990/entanglement-analysis",
    packages=setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires = '>=3.6',
)

from setuptools import setup, find_packages
