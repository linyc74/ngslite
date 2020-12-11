from setuptools import setup, find_packages


def get_version() -> str:
    with open('ngslite/__init__.py') as fh:
        for line in fh:
            if line.startswith('__version__'):
                return line.split('=')[1].strip()[1:-1]


def main():

    setup(
        name='ngslite',
        version=get_version(),
        description='Light-weight functions for next-generation sequencing (NGS) data analysis',
        url='https://github.com/linyc74/ngslite',
        author='Yu-Cheng Lin',
        author_email='yclin.python@gmail.com',
        license='MIT',
        packages=find_packages(),
        python_requires='>3.6',
        install_requires=['numpy', 'pandas'],
        zip_safe=False
    )


if __name__ == '__main__':
    main()

