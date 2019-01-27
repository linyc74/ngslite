from setuptools import setup

setup(name='ngslite',
      version='0.5',
      description='Light-weight functions for next-generation sequencing (NGS) data analysis',
      url='https://github.com/linyc74/ngslite',
      author='Yu-Cheng Lin',
      author_email='yclin.python@gmail.com',
      license='MIT',
      packages=['ngslite'],
      python_requires='>3.6',
      install_requires=['numpy'],
      zip_safe=False)

# In command line type in: "python setup.py sdist"

