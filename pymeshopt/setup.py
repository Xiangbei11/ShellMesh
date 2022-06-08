from setuptools import setup
from setuptools.extension import Extension
import numpy as np
ext_modules = [Extension("pymeshopt",sources=["pymeshopt.pyx","cmeshopt.cpp"],language="c++")]
setup(
    name = 'pymeshopt',
    ext_modules = ext_modules,
    include_dirs=[np.get_include()],
    setup_requires  = ["numpy"],
    install_requires = ["numpy"],
)