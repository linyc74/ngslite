from setuptools import setup


def main():
    setup(
        name='ngslite',
        version='1.1.4',
        description='Light-weight functions for next-generation sequencing (NGS) data analysis',
        url='https://github.com/linyc74/ngslite',
        author='Yu-Cheng Lin',
        author_email='yclin.python@gmail.com',
        license='MIT',
        packages=['ngslite'],
        python_requires='>3.6',
        install_requires=[
            'numpy>=1.18',
            'matplotlib>=3.2',
            'pandas>=1.0',
            'biopython>=1.77',
            'scipy>=1.5',
            'dna_features_viewer>=3.0',
        ],
        zip_safe=False
    )


if __name__ == '__main__':
    main()

    # In command line type in: "python setup.py sdist"
