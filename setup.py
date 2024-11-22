import setuptools
with open('requirements.txt') as f:
    required = f.read().splitlines()
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="miniQuant",                     # This is the name of the package
    version="1.0",                        # The initial release version
    author="Haoran Li",                     # Full name of the author
    description="M͟i͟xed Bayesian n̲etwork for i̲soform quantification (miniQuant) provides a highly-accurate bioinformatics tool for transcript abundance estimation.",
    long_description=long_description,      # Long description read from the the readme file
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),    # List of all python modules to be installed
    license_files = ('LICENSE',),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],                                      # Information to filter the project on PyPi website
    python_requires='>=3.8',                # Minimum version requirement of the package
    py_modules=["miniQuant"],             # Name of the python package
    package_dir={'':'isoform_quantification'},     # Directory of the source code of the package
    install_requires=required                    # Install other dependencies if any
)
