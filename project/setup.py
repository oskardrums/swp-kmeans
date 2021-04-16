from setuptools import setup, Extension

setup(
    name='clustering_project',
    version='0.1.0',
    ext_modules=[
        Extension(
            'clustering',
            sources=['kmeans.c', 'mat.c', 'clustering.c', 'nsc.c', 'log.c', 'jaccard.c']
        ),
    ]
)
