from setuptools import find_packages, setup

setup(
    name='qualiloop',
    packages=find_packages(),
    version='0.1.0',
    description='My first Python library',
    author='Lilian Denzler',
    license='UCL',
    install_requires=[],#here put bioptools etc
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
