#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()


install_requires = [
    "geopandas",
    'vtools3>=3.8',         # moved unit_conversions to vtools, so dependent on the version where we did that
    "gdal",
    "libgdal",  # libgdal is for calling gdal_translate in clip_dems
    "xarray",
    "netcdf4",
    "scipy",
    "matplotlib",
    "statsmodels",
    "palettable",
    "pyyaml",
    "boto3",  # only used for download_hrrr
    "statsmodels",
    "beautifulsoup4",
    "pyproj",
    "nodepy",
    "sympy",  # sympy is dependency of nodepy
    "shapely>=2.0",
    "rasterstats",
    "param",
]

requirements = install_requires

setup_requirements = [
    "pytest-runner",
]

test_requirements = [
    "pytest>=3",
]

setup(
    author="California Department of Water Resources",
    author_email="Eli.Ateljevich@water.ca.gov, Kijin.Nam@water.ca.gov",
    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="Python package for SCHISM",
    entry_points={
        "console_scripts": [
            "batch_metrics=schimpy.batch_metrics:main",
            "clip_dems=schimpy.clip_dems:main",
            "contour_smooth=schimpy.contour_smooth:main",
            "convert_mesh=schimpy.convert_mesh:main",
            "convert_polygons=schimpy.convert_polygons:main",
            "convert_linestrings=schimpy.convert_linestrings:main",
            "combine_consume=schimpy.combine_consume:main",
            "prepare_schism=schimpy.prepare_schism:main",
            "create_vgrid_lsc2=schimpy.create_vgrid_lsc2:main",
            "schism_hotstart=schimpy.schism_hotstart:main",
            "split_quad=schimpy.split_quad:main",
            "model_time=schimpy.model_time:main",
            "gen_elev2d=schimpy.gen_elev2d:main",
            "small_areas=schimpy.small_areas:main",
            "station=schimpy.station:main",
            "create_hotstart=schimpy.schism_hotstart:main",
            "create_nudging=schimpy.nudging:main",
            "interpolate_structure=schimpy.interpolate_structure:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="schimpy",
    name="schimpy",
    packages=find_packages(include=["schimpy", "schimpy.*"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/CADWRDeltaModeling/schimpy",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
)
