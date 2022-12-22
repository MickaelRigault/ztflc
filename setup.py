#! /usr/bin/env python
#

DESCRIPTION = "ztflc: Force photometry lc fitter"
LONG_DESCRIPTION = """ Force photometry lc fitter"""

DISTNAME = "ztflc"
AUTHOR = "Mickael Rigault"
MAINTAINER = "Mickael Rigault"
MAINTAINER_EMAIL = "m.rigault@ipnl.in2p3.fr"
URL = "https://github.com/MickaelRigault/ztflc/"
LICENSE = "BSD (3-clause)"
DOWNLOAD_URL = "https://github.com/MickaelRigault/ztflc/tarball/0.2"
VERSION = "0.2.8"

try:
    from setuptools import setup, find_packages

    _has_setuptools = True
except ImportError:
    from distutils.core import setup

    _has_setuptools = False


def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import ztfquery
    except ImportError:
        install_requires.append("ztfquery")

    try:
        import pandas
    except ImportError:
        install_requires.append("pandas")

    try:
        import iminuit
    except ImportError:
        install_requires.append("iminuit")

    try:
        import numpy
    except ImportError:
        install_requires.append("numpy")

    return install_requires


if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        # This should be updated if new submodules are added
        packages = ["ztflc"]

    setup(
        name=DISTNAME,
        author=AUTHOR,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        install_requires=install_requires,
        scripts=["bin/forcephoto.py"],
        packages=packages,
        include_package_data=True,
        # package_data={'pysedm': ['data/*.*']},
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.5",
            "License :: OSI Approved :: BSD License",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
    )
