import setuptools

with open("README.md", "r") as fh:
    description = fh.read()
from setuptools import setup, find_packages
setuptools.setup(
    name="IRonPy",
    version="0.3.88",
    author="berre",
    packages=find_packages(),
    package_data={'': ['*.dll'], '': ['src/*.so'], '': ['data/*.csv']},
    zip_safe=False,
    include_package_data=True,
    author_email="berreergunn@gmail.com",
    description="IRon Package",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/berreergun/IRonPy",
    license='MIT',
    python_requires='>=2',
    install_requires=["pywin32 >= 1.0;platform_system=='Windows'","matplotlib",'numpy','pandas','sklearn']
)