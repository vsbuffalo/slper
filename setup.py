with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="slper",
    version="0.0.1",
    author="Vince Buffalo",
    author_email="vsbuffaloAAAAA@gmail.com",
    description="SLiM Helper",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vsbuffalo/slper",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
