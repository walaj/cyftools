from setuptools import setup, find_packages

setup(
    name="pysift",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "pysift = pysift.core:main",
        ],
    },
    # additional metadata
)
