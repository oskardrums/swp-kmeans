from setuptools import setup, Extension

setup(
    name='mykmeanspp',
    version='0.1.0',
    ext_modules=[
        Extension('mykmeanspp', ['kmeans.c']),
    ]
)
