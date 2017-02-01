from setuptools import setup

versionfile = 'elisa/version.py'
exec(compile(open(versionfile, 'rb').read(), versionfile, 'exec'))

setup(
    name='liulab_elisa',
    version=__version__,
    url='https://github.com/WendyLiuLab/elisascripts',
    license='MIT',
    author='Tim D. Smith',
    author_email='elisa@tim-smith.us',
    description='Analyses single-cell ELISA image stacks.',
    packages=['elisa', 'elisa.test'],
    package_data={'elisa.test': ['*.tif']},
    platforms='any',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
    ],
    install_requires=[
        # 'cairo', but it isn't on pypi
        'ijroi',
        'numpy',
        'matplotlib',
        'pandas',
        'pillow',
        'scipy',
        'shapely',
        'tifffile',
        'typing',
    ],
    entry_points={
        'console_scripts': [
            'elisa_id_singlets = elisa.id_singlets:main',
            'elisa_normalize_bg = elisa.normalize_bg:main',
            'elisa_register = elisa.register:main',
            'elisa_stitch = elisa.stitch:main',
        ]
    },
)
