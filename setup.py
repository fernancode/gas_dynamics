import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name="gas-dynamics",
      version="0.2.1",
      author="Fernando de la Fuente",
      author_email="FernandoAdelaFuente@gmail.com",
      description="Package containing functions for working with compressible flow",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="http://github.com/fernancode/gas-dynamics",
      packages=setuptools.find_packages(),
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent"
      ],
      license="MIT",
      python_requires=">=3.6",
)