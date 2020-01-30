
from os import path
from setuptools import setup
import referenceseeker


# Get the long description from the README file
setup_dir = path.abspath(path.dirname(__file__))
with open(path.join(setup_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='referenceseeker',
    version=referenceseeker.__version__,
    description='ReferenceSeeker: rapid determination of appropriate reference genomes.',
    keywords=['bioinformatics', 'ngs', 'wgs', 'microbial genomics'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPLv3',
    author='Oliver Schwengers',
    author_email='oliver.schwengers@computational.bio.uni-giessen.de',
    url='https://github.com/oschwengers/referenceseeker',
    packages=['referenceseeker'],
    python_requires='>=3.5',
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'biopython >= 1.71'
    ],
    entry_points={
        'console_scripts': [
            'referenceseeker=referenceseeker.referenceseeker:main',
            'referenceseeker_db=referenceseeker.database:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Bug Reports': 'https://github.com/oschwengers/referenceseeker/issues',
        'Source': 'https://github.com/oschwengers/referenceseeker'
    },
)
