import setuptools
from gas_dynamics.__about__ import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name="gas_dynamics",
      version=__version__,
      author="Fernando de la Fuente",
      author_email="FernandoAdelaFuente@gmail.com",
      description="Package containing functions for working with compressible flow",
      long_description=long_description,
      long_description_content_type="text/markdown",
      project_urls = {"ReadtheDocs" : "http://gas-dynamics.readthedocs.io", "Github" : "http://github.com/fernancode/gas_dynamics"},
      packages=setuptools.find_packages(),
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent"
      ],
      license="MIT",
      python_requires=">=3.6",
      tests_require=['pytest'],
      setup_requires=["numpy==1.19.3", "pytest-runner"],
      install_requires=["numpy==1.19.3", "scipy", "matplotlib==3.3.2"]
)