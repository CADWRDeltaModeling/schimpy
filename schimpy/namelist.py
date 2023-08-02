def parse(file_content):
    """
    Here's a simple implementation of a parser for the Fortran namelist format 
    If there's a line starting with !, the parser will store it as the full_line_comment and attach it to the next key-value pair it encounters. 
    If there's an inline comment starting with !, it will be stored as inline_comment.

    For the most part this parser is comment preserving; you might notice some whitespace and blank line(s) being eliminated in the round trip through write method below
    """
    namelists = {}
    current_namelist = None
    full_line_comment = ""
    for line in file_content.splitlines():
        line = line.strip()
        if line.startswith("!"):
            if len(full_line_comment) == 0:
                full_line_comment = line[1:]
            else:
                full_line_comment = full_line_comment + '\n' + line[1:]
        elif line.startswith("&"):
            current_namelist = line[1:].strip()
            namelists[current_namelist] = {"comment": full_line_comment}
            full_line_comment = ""
        elif line.startswith("/"):
            current_namelist = None
        elif current_namelist is not None:
            line_parts = line.split("!", maxsplit=1)
            key_value_pair = line_parts[0].strip()
            inline_comment = line_parts[1] if len(line_parts) == 2 else ""
            if not key_value_pair:
                continue
            key, value = key_value_pair.split("=")
            key = key.strip()
            value = value.strip()
            try:
                value = int(value)
            except:
                try:
                    value=float(value)
                except:
                    pass
            namelists[current_namelist][key] = {
                "value": value, "full_line_comment": full_line_comment, "inline_comment": inline_comment}
            full_line_comment = ""
    return namelists


def write(namelist_data):
    """
    writes out the namelist dictionary from the parse method to a string. The comments are preserved but not the whitespace or indentations
    """
    lines = []
    for namelist_name, namelist_values in namelist_data.items():
        comment = namelist_values.get("comment")
        if comment and len(comment) != 0:
            lines.append("!{}".format(comment.replace("\n", "\n!")))
        lines.append("&{}".format(namelist_name))
        for key, value_data in namelist_values.items():
            if key == "comment":  # already taken care of above
                continue
            if "full_line_comment" in value_data:  # do this first
                full_line_comment = value_data.get("full_line_comment")
                if len(full_line_comment) != 0:
                    lines.append("!{}".format(
                        full_line_comment.replace("\n", "\n!")))
            if key == "full_line_comment":
                continue  # took care of it above
            lines.append("  {} = {}".format(key, value_data["value"]))
            inline_comment = value_data.get("inline_comment")
            if inline_comment and len(inline_comment) != 0:
                lines[-1] = lines[-1] + \
                    (" !{}".format(inline_comment.replace('\n', '')))
        lines.append("/")
    return "\n".join(lines)


class TopLevelNamelist:
    """
    TopLevelNamelist is a class for the top-level elements of the dictionary, and Namelist is a class for all other nested elements. 
    The TopLevelNamelist class creates attributes for each key-value pair in the dictionary and, if the value is itself a dictionary, creates a new Namelist object with the nested dictionary. 
    """

    def __init__(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                setattr(self, k, Namelist(v))
            else:
                setattr(self, k, v)


class Namelist:
    """
    The Namelist class simply creates an attribute for each key-value pair in the dictionary.
    """

    def __init__(self, d):
        for k, v in d.items():
            setattr(self, k, v)

    def __eq__(self, other):
        return self.value == other


def map_to_object(d):
    """
    In this example, TopLevelNamelist is a class for the top-level elements of the dictionary, and Namelist is a class for all other nested elements. 
    The TopLevelNamelist class creates attributes for each key-value pair in the dictionary and, if the value is itself a dictionary, 
    creates a new Namelist object with the nested dictionary. The Namelist class simply creates an attribute for each key-value pair in the dictionary.

    The map_to_object function takes a dictionary d as input and returns a new TopLevelNamelist object with the values from d.

    This way, you can easily map the dictionary returned by parse_namelist to objects, with different classes for the top-level and nested elements, 
    and access the values in the dictionary as attributes of the objects.
    """
    return TopLevelNamelist(d)


def map_to_dict(obj):
    """
    This function takes an object of type TopLevelNamelist as input and returns a dictionary representation of the object. 
    The function uses the vars built-in function to get a dictionary of attributes and values for the object, and then iterates over the key-value pairs in the dictionary.
    If the value is an instance of the Namelist class, it recursively maps the Namelist object to a dictionary using the map_to_dict function. If the value is not an instance of the Namelist class, it simply adds the key-value pair to the dictionary.
    This way, you can easily map an object of type TopLevelNamelist back to a dictionary representation, preserving the structure of the original dictionary.
    """
    d = {}
    for k, v in vars(obj).items():
        if isinstance(v, Namelist):
            d[k] = map_to_dict(v)
        else:
            d[k] = v
    return d
