from schimpy.schism_yaml import load, load_raw
import pandas as pd
import string
import yaml
import io


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


class NamedStringIO(io.StringIO):
    def __init__(self, value, name="in_memory.yaml"):
        super().__init__(value)
        self.name = name


def yaml_from_dict(input_dict, envvar=None):
    """
    Convert a dictionary to a YAML string with environment variable substitution.

    Parameters
    ----------
    input_dict : dict
        The dictionary to convert to YAML.
    envvar : dict, optional
        Environment variables to substitute in the YAML output.

    Returns
    -------
    str
        The YAML representation of the dictionary.
    """
    yaml_str = yaml.safe_dump(input_dict)
    stream = NamedStringIO(yaml_str, name="in_memory.yaml")

    return load(stream)


def yaml_from_file(filename, envvar=None, raw=False):
    """
    Load a YAML file and return its contents.

    Parameters
    ----------
    filename : str|Path
        The path to the YAML file.
    envvar : dict, optional
        Variables to substitute into the YAML.
    raw : bool, optional
        If True, leave every scalar as a string instead of resolving YAML
        types. Needed where the literal text of a value matters.

    Returns
    -------
    dict
        The contents of the YAML file as a dictionary.
    """
    with open(filename, "r") as file:
        if raw:
            return load_raw(file, envvar=envvar)
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
        file.write(yaml.safe_dump(data))

def write_yaml(data, outfile):
    """
    Write a dictionary to a YAML file.

    Parameters
    ----------
    data : dict
        The data to write to the YAML file.
    outfile : str|Path
        The path to the output YAML file.
    """
    with open(outfile, "w") as file:
        file.write(yaml.safe_dump(data))
