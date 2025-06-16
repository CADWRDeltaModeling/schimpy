from schimpy.schism_yaml import load, dump
from pathlib import Path
import pandas as pd
import string


def csv_from_file(filename, envvar=None, **kwargs):
    """
    Load a CSV file and substitute environment variables in string fields.

    Parameters
    ----------
    filename : str
        Path to the CSV file.
    envvar : dict, optional
        Dictionary of variables to substitute (e.g., {'calsim_dss': 'blah.dss'}).
    kwargs : passed to pd.read_csv

    Returns
    -------
    pd.DataFrame
        DataFrame with substitutions applied.
    """
    df = pd.read_csv(filename, **kwargs)
    if envvar is None:
        return df

    # Substitute in column names
    df.columns = [
        string.Template(str(col)).safe_substitute(**envvar) for col in df.columns
    ]

    # Substitute in index (if it's string/object)
    if df.index.dtype == "object":
        df.index = [
            string.Template(str(idx)).safe_substitute(**envvar) for idx in df.index
        ]

    # Substitute in all string/object cells
    for col in df.select_dtypes(include=["object"]).columns:
        df[col] = (
            df[col]
            .astype(str)
            .apply(lambda x: string.Template(x).safe_substitute(**envvar))
        )

    return df


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
