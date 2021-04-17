from setuptools import setup, Extension

setup(
    name='clustering_project',
    version='0.1.0',
    ext_modules=[
        Extension(
            'clustering',
            sources=['kmeans.c', 'mat.c', 'clustering.c', 'nsc.c', 'jaccard.c'],
            extra_link_args=["-flto"],
            extra_compile_args=["-std=c99", "-Wall", "-Wextra", "-Wpedantic", "-Werror", "-flto"]
        ),
    ]
)
