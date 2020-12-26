from setuptools import setup, find_packages, Extension
setup(
    name='mykmeanspp',
    version='0.1.0',
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            'mykmeanspp',
            # the files to compile into our module relative to ``setup.py``
            ['kmeans.c'],
        ),
    ]
)
