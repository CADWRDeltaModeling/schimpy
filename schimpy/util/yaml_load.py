from schimpy.schism_yaml import load, dump
from pathlib import Path


def yaml_from_file(filename, envvar=None):
    """
    Load a YAML file and return its contents.

    Parameters
    ----------
    filename : str|Path
        The path to the YAML file.

    Returns
    -------
    dict
        The contents of the YAML file as a dictionary.
    """
    with open(filename, "r") as file:
        return load(file, envvar=envvar)


def yaml_to_yaml(infile, outfile, envvar=None):
    """
    Load a YAML file and write its contents to another YAML file.

    Parameters
    ----------
    infile : str|Path
        The path to the input YAML file.
    outfile : str|Path
        The path to the output YAML file.
    envvar : dict, optional
        Environment variables to substitute in the YAML file.
    """
    data = yaml_from_file(infile, envvar=envvar)
    with open(outfile, "w") as file:
        file.write(data)
