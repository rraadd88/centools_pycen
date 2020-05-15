import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="centools_pycen",
    version="0.0.1",
    author="",
    author_email="",
    description="CEN-tools - pyCEN package",
    long_description="PyCEN is the Python implementation of the CEN-tools web interface (http://cen-tools.com/) to analyse contexts in pre-defined or user-defined data in detail.",
    long_description_content_type="text/markdown",
    url="https://gitlab.ebi.ac.uk/petsalakilab/centools_pycen/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
