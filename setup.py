from setuptools import setup

setup(
    name='imgMS',
    version='0.0.1',
    author='nikadilli',
    author_email='nikadilli@gmail.com',
    packages=['imgMS'],
    scripts=[],
    license='LICENSE.txt',
    description='Package for data reduction of LA-ICP-MS data.',
    long_description=open('README.md').read(),
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
