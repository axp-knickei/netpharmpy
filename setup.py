from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="netpharmpy",
    version="1.0.0",
    description="Network pharmacology pipeline for drug discovery and target prediction",
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.8",
)
