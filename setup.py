from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="ctdna-lengths",
    version="0.0.1",
    author="Per Høgfeldt",
    description="A script for extracting the length distribution of ctDNA WGS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #TODO: fix url
    url="https://github.com/hogfeldt/",
    #packages=find_packages("src"),
    #package_dir={"": "src"},
    #test_suite="test",
    install_requires=[
        "pysam",
        "numpy",
    ],
    #TODO fix endpoints
    entry_points={
        "console_scripts": [
            "ctdna-lengths=scripts.ctdna_lengths:main",
        ]
    },
    include_package_data=True,
    python_requires=">=3.6",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
