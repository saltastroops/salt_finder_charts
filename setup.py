"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['numpy', 'aplpy', 'astropy', 'ephem', 'pillow']

setup_requirements = ['pytest-runner']

test_requirements = ['pytest']

setup(
    author="SAAO/SALT",
    author_email='salthelp@salt.ac.za',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Generate finder charts for the Southern African Large Telescope (SALT)",
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='salt_finder_charts',
    name='salt_finder_charts',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/saltastroops/salt_finder_charts',
    version='0.1.0',
    zip_safe=False,
)
