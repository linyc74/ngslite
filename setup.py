from setuptools import setup


def main():
    setup(
        name='ngslite',
        version='1.2.0',
        description='Light-weight functions for next-generation sequencing (NGS) data analysis',
        url='https://github.com/linyc74/ngslite',
        author='Yu-Cheng Lin',
        author_email='yclin.python@gmail.com',
        license='MIT',
        packages=['ngslite'],
        python_requires='>3.6',
        install_requires=['numpy>=1.18'],
        zip_safe=False
    )


if __name__ == '__main__':
    main()

    # In command line type in: "python setup.py sdist"
