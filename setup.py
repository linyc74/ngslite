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
        install_requires=['numpy', 'pandas'],
        zip_safe=False
    )


if __name__ == '__main__':
    main()

    # python setup.py sdist
    # cd sdist
    # twine upload --repository pypi ngslite-[VERSION].tar.gz

