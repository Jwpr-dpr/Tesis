from setuptools import setup, find_packages

setup(
    name="Photonics",
    version="0.1",
    packages=find_packages(),
    description="Computational package to perform calculations on photonic crystals in python",
    author="Jhon Wilmer Pino Román, Hernán Alejandro Gomez",
    author_email="Jwprdpr@gmail.com , hagomez@udemedellin.edu.co",
    long_description=open('README.md').read(),
    install_requires = ['numpy','scipy','matplotlib'],
    url="https://github.com/Jwpr-dpr/Tesis.git"

)
    