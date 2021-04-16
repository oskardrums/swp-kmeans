from setuptools import setup, Extension

setup(
    name='clustering',
    version='0.1.0',
    ext_modules=[
        Extension('clustering', ['kmeans.c', 'mat.c', 'clustering.c', 'nsc.c', 'log.c', 'jaccard.c']),
    ]
)
