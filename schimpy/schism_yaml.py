"""A customized version of YAML parser for SCHISM
It stores document in an ordered dict and supports variable substitution.
"""

import yaml
from yaml.reader import *
from yaml.scanner import *
from yaml.parser import *
from yaml.composer import *
from yaml.constructor import *
from yaml.resolver import *
from yaml.nodes import *
from yaml.dumper import *
from yaml.loader import *
from yaml.tokens import *
import collections
import string
import re
import os
import warnings
import sys

__all__ = ["load", "dump"]
include_keywords = [
    "include",
]
substitute_keywords = [
    "config",
    "env",
]
yaml_extensions = [
    ".yaml",
]


def check_env(env):
    """Check the integrity of environment variables"""
    msg = ""
    n_errors = 0
    for k, v in env.items():
        if "$" in k:
            msg += "A Name of environment variable cannot be a variable: " "%s\n" % k
            n_errors += 1
        if not isinstance(v, (str, float, int)):
            msg += "Environment variables do not accept non-scalar variables: %s\n" % k
            n_errors += 1
    if n_errors > 0:
        raise ValueError(msg)


def substitute_env(env):
    """Substitute environmental variables. Automatically modifies env in place."""
    if env is None:
        return
    check_env(env)
    while True:
        count_substitutions = 0
        for k, v in env.items():
            template = string.Template(v)
            substituted = template.safe_substitute(**env)
            if v != substituted:
                env[k] = substituted
                count_substitutions += 1
        if count_substitutions == 0:
            break


class SubstituteComposer(Composer):
    """Composer with substitution"""

    def __init__(self, env=None):
        self.env = env
        super(SubstituteComposer, self).__init__()

    def compose_scalar_node(self, anchor):
        if self.env is None:
            return super(SubstituteComposer, self).compose_scalar_node(anchor)
        event = self.get_event()
        tag = event.tag
        if tag is None or tag == "!":
            tag = self.resolve(ScalarNode, event.value, event.implicit)
        template = string.Template(event.value)
        value = template.safe_substitute(**self.env)
        if "$" in value:
            raise ComposerError(
                "Expected a substitution",
                event.start_mark,
                "No corresponding config variable",
            )
        node = ScalarNode(
            tag, value, event.start_mark, event.end_mark, style=event.style
        )
        if anchor is not None:
            self.anchors[anchor] = node
        return node


class SubstituteConstructor(SafeConstructor):
    """Customized constructor for including.
    If a key value is one of include_keywords, a following YAML file
    will be included.
    """

    def __init__(self, env=None):
        self.env = env
        super(SubstituteConstructor, self).__init__()
        # super(SubstituteConstructor, self).add_constructor('!include',
        # SubstituteConstructor.include)

    def construct_pairs(self, node, deep=False):
        """Overridden construct_pairs function with substitution"""
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None,
                None,
                "expected a mapping node, but found %s" % node.id,
                node.start_mark,
            )
        pairs = []
        for key_node, value_node in node.value:
            if key_node.value in include_keywords:
                included = self.include(value_node)
                if isinstance(included, dict):
                    pairs.extend(iter(included.items()))
                else:
                    raise ValueError("Included YAML must be a dictionary")
            else:
                key = self.construct_object(key_node, deep=deep)
                value = self.construct_object(value_node, deep=deep)
                pairs.append((key, value))
        return pairs

    def include(self, node):
        """Process a node with !include tag.
        This code is copied from Stackoverflow.
        """
        if isinstance(node, yaml.ScalarNode):
            return self.extractFile(self.construct_scalar(node))
        elif isinstance(node, yaml.SequenceNode):
            fns = self.construct_sequence(node)
            if ".yaml" in fns[0]:
                result = {}  # create dictionary to store all values in .yaml list
                for filename in self.construct_sequence(node):
                    tmp_rslt = self.extractFile(filename)
                    for key in tmp_rslt.keys():
                        if key in result.keys():
                            result[key] = result[key] + tmp_rslt[key]
                        else:
                            result[key] = tmp_rslt[key]
            else:
                result = []
                for filename in self.construct_sequence(node):
                    result += self.extractFile(filename)
            return result
        elif isinstance(node, yaml.MappingNode):
            result = {}
            for k, v in self.construct_mapping(node).items():
                result[k] = self.extractFile(v)
            return result
        else:
            print("Error: unrecognized node type in the include constructor")
            raise ConstructorError

    def extractFile(self, filename):
        """Load an yaml file to include"""
        filepath = os.path.join(self.root, filename)
        if os.path.exists(filepath):
            with open(filepath, "r") as f:
                loader = SubstituteLoader(f, self.env)
                return loader.get_single_data()
        else:
            print("Error: cannot find the file:" + filename)
            raise ValueError("Cannot find included file: {}".format(filename))


SubstituteConstructor.add_constructor("!include", SubstituteConstructor.include)


class SubstituteLoader(
    Reader, Scanner, Parser, SubstituteComposer, SubstituteConstructor, Resolver
):
    """Raw Loader with substitute."""

    def __init__(self, stream, env):
        """Constructor"""
        self.root = os.path.split(stream.name)[0]
        Reader.__init__(self, stream)
        Scanner.__init__(self)
        Parser.__init__(self)
        SubstituteComposer.__init__(self, env)
        SubstituteConstructor.__init__(self, env)
        Resolver.__init__(self)


class SubstituteRawLoader(
    Reader, Scanner, Parser, SubstituteComposer, SubstituteConstructor, BaseResolver
):
    """Raw Loader with substitute."""

    def __init__(self, stream, env):
        """Constructor"""
        self.root = os.path.split(stream.name)[0]
        Reader.__init__(self, stream)
        Scanner.__init__(self)
        Parser.__init__(self)
        SubstituteComposer.__init__(self, env)
        SubstituteConstructor.__init__(self, env)
        BaseResolver.__init__(self)


class RawLoader(Reader, Scanner, Parser, Composer, SubstituteConstructor, BaseResolver):
    """Raw Loader, No automatic resolving."""

    def __init__(self, stream):
        self.root = os.path.split(stream.name)[0]
        Reader.__init__(self, stream)
        Scanner.__init__(self)
        Parser.__init__(self)
        Composer.__init__(self)
        SubstituteConstructor.__init__(self)
        BaseResolver.__init__(self)


def dict_constructor(loader, node):
    """Constructor with OrderedDict"""
    return collections.OrderedDict(loader.construct_pairs(node))


def dict_representer(dumper, data):
    """Dumper with OrderedDict"""
    return dumper.represent_mapping(
        BaseResolver.DEFAULT_MAPPING_TAG, iter(data.items())
    )


# Add constructor and representer
yaml.add_representer(collections.OrderedDict, dict_representer, Dumper=SafeDumper)
yaml.add_constructor(
    BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor, Loader=SubstituteLoader
)
yaml.add_constructor(
    BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor, Loader=RawLoader
)


def check_unsubstituted_vars(data, path=None, found=None):
    """Recursively check for unsubstituted {var} in any string in the data structure."""
    import re

    if found is None:
        found = []
    if path is None:
        path = []
    if isinstance(data, dict):
        for k, v in data.items():
            check_unsubstituted_vars(v, path + [str(k)], found)
    elif isinstance(data, (list, tuple)):
        for idx, v in enumerate(data):
            check_unsubstituted_vars(v, path + [str(idx)], found)
    elif isinstance(data, str):
        matches = re.findall(r"\{[^{}]+\}", data)
        if matches:
            found.append((".".join(path), matches))
    return found


def collapsed_nested_matches(node):
    for key, item in node.items():
        if isinstance(item, dict):
            if key in item.keys():
                node[key] = item[key]
            collapsed_nested_matches(item)


def load(stream, envvar=None, collapse_nested=True):
    """Load a schism YAML"""
    # First round to get environmental variables
    loader = RawLoader(stream)
    try:
        data = loader.get_single_data()
    finally:
        loader.dispose()
    stream.seek(0)

    # Second round with substitution
    env = {}
    for key in substitute_keywords:
        if key in data:
            env.update(data[key])
    # If envvar is provided, update the environment variables
    # with the provided variables
    if envvar is not None:
        env.update(envvar)
    substitute_env(env)
    loader = SubstituteLoader(stream, env)
    try:
        data = loader.get_single_data()
        unsub = check_unsubstituted_vars(data)
        if unsub:
            msg = [f"Unsubstituted variables found in YAML:"]
            for path, matches in unsub:
                msg.append(f"  {path}: {matches}")
            msg.append(r"This could be an issue of using {var} instead of ${var}")
            warnings.warn("\n".join(msg))
            
        collapsed_nested_matches(data)
        return data
    finally:
        loader.dispose()


def load_raw(stream):
    """Load a schism YAML"""
    # First round to get environmental variables
    loader = RawLoader(stream)
    try:
        data = loader.get_single_data()
    finally:
        loader.dispose()
    stream.seek(0)

    # Second round with substitution
    env = {}
    for key in substitute_keywords:
        if key in data:
            env.update(data[key])
    substitute_env(env)
    loader = SubstituteRawLoader(stream, env)
    try:
        data = loader.get_single_data()
        return data
    finally:
        loader.dispose()


def dump(data, stream=None, Dumper=Dumper, **kwds):
    """
    Serialize a Python object into a YAML stream.
    If stream is None, return the produced string instead.
    """
    return yaml.dump(data, stream, Dumper=Dumper, **kwds)


def safe_dump(data, stream=None, **kwds):
    """
    Serialize a Python object into a YAML stream.
    If stream is None, return the produced string instead.
    """
    return yaml.safe_dump(data, stream, **kwds)


import click


@click.command()
@click.argument("filename", type=click.Path(exists=True))
def schism_yaml_cli(filename):
    """Load and display a YAML file that follows SCHISM YAML format."""
    with open(filename, "r") as f:
        data_ = load(f)
        print(yaml.safe_dump(data_))


if __name__ == "__main__":
    schism_yaml_cli()
