import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='imgMS',
    version='0.2.3',
    author='nikadilli',
    author_email='nikadilli@gmail.com',
    url="https://github.com/nikadilli/imgMS",
    packages=setuptools.find_packages(),
    scripts=[],
    license='LICENSE.txt',
    description='Package for data reduction of LA-ICP-MS data.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "matplotlib",
        "pandas",
        "numpy",
        "patsy",
        "Pillow",
        "scikit-learn",
        "scipy",
        "sklearn",
        "statsmodels",
        "xlrd",
        "XlsxWriter"
    ],
)
