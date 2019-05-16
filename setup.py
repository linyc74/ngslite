from setuptools import setup
import ngslite


setup(
    name='ngslite',
    version=ngslite.__version__,
    description='Light-weight functions for next-generation sequencing (NGS) data analysis',
    url='https://github.com/linyc74/ngslite',
    author='Yu-Cheng Lin',
    author_email='yclin.python@gmail.com',
    license='MIT',
    packages=['ngslite'],
    python_requires='>3.6',
    install_requires=['numpy', 'pandas', 'biopython', 'scipy'],
    zip_safe=False
)


# In command line type in: "python setup.py sdist"

