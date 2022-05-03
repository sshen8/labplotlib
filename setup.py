import setuptools

setuptools.setup(
    name="labplotlib",
    version="1.0.0",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    install_requires=["numpy", "Bio", "matplotlib", "pandas", "shapely"],
)
