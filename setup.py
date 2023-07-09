import setuptools

setuptools.setup(
    name="labplotlib",
    version="1.2.0",
    packages=setuptools.find_packages(),
    package_data={"labplotlib": ["blood/data/blank_cells.eps.xml"]},
    python_requires='>=3.6',
    install_requires=["numpy", "Bio", "matplotlib", "pandas", "shapely"],
)
