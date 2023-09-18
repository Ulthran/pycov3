"""Python setup.py for pycov3 package"""
import io
import os
from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("pycov3", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="pycov3",
    version=read("pycov3", "VERSION"),
    description="Generate cov3 files used in DEMIC",
    url="https://github.com/Ulthran/pycov3/",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="Charlie Bushman",
    packages=find_packages(exclude=["tests", ".github"]),
    install_requires=read_requirements("requirements.txt"),
    entry_points={"console_scripts": ["pycov3 = pycov3.__main__:main"]},
    extras_require={"test": read_requirements("requirements-test.txt")},
)
