from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "eigen_backend",                # Name of the Python module
        ["eigen_backend.cpp"],          # Source file
        include_dirs=["./Eigen"],       # Relative path to the Eigen headers
    ),
]

setup(
    name="eigen_backend",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
