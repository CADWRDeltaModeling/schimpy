def parse(file_content):
    """
    Here's a simple implementation of a parser for the Fortran namelist format
    If there's a line starting with !, the parser will store it as the full_line_comment and attach it to the next key-value pair it encounters.
    If there's an inline comment starting with !, it will be stored as inline_comment.

    For the most part this parser is comment preserving; you might notice some whitespace and blank line(s) being eliminated in the round trip through write method below.
    NOTE: If a comment block is preceded by a blank line, this fact is now preserved and used when writing back out.
    """
    namelists = {}
    current_namelist = None

    full_line_comment = ""
    full_line_comment_leading_blank = False  # NEW: did a blank line precede this comment block?
    prev_blank = False  # NEW: was the previous raw line blank?

    for raw_line in file_content.splitlines():
        stripped = raw_line.strip()

        # Track blanks explicitly
        if stripped == "":
            prev_blank = True
            continue

        if stripped.startswith("!"):
            # Full-line comment, possibly spanning multiple lines
            if len(full_line_comment) == 0:
                full_line_comment = stripped[1:]
                full_line_comment_leading_blank = prev_blank
            else:
                full_line_comment = full_line_comment + "\n" + stripped[1:]
                # keep whatever leading_blank state we already had
            prev_blank = False

        elif stripped.startswith("&"):
            current_namelist = stripped[1:].strip()
            namelists[current_namelist] = {
                "comment": full_line_comment,
            }
            # NEW: remember if the namelist comment was preceded by a blank line
            if len(full_line_comment) != 0:
                namelists[current_namelist]["comment_leading_blank"] = full_line_comment_leading_blank
            full_line_comment = ""
            full_line_comment_leading_blank = False
            prev_blank = False

        elif stripped.startswith("/"):
            current_namelist = None
            full_line_comment = ""
            full_line_comment_leading_blank = False
            prev_blank = False

        elif current_namelist is not None:
            # Key/value line (possibly with inline comment)
            line_parts = stripped.split("!", maxsplit=1)
            key_value_pair = line_parts[0].strip()
            inline_comment = line_parts[1] if len(line_parts) == 2 else ""
            if not key_value_pair:
                prev_blank = False
                continue

            key, value = key_value_pair.split("=")
            key = key.strip()
            value = value.strip()

            try:
                value = int(value)
            except Exception:
                try:
                    value = float(value)
                except Exception:
                    pass

            entry = {
                "value": value,
                "full_line_comment": full_line_comment,
                "inline_comment": inline_comment,
            }
            # NEW: remember if the key's full-line comment was preceded by a blank line
            if len(full_line_comment) != 0:
                entry["full_line_comment_leading_blank"] = full_line_comment_leading_blank

            namelists[current_namelist][key] = entry

            full_line_comment = ""
            full_line_comment_leading_blank = False
            prev_blank = False

        else:
            # Outside any namelist and not a comment; just reset flags
            full_line_comment = ""
            full_line_comment_leading_blank = False
            prev_blank = False

    return namelists


def write(namelist_data):
    """
    writes out the namelist dictionary from the parse method to a string. The comments are preserved but not the whitespace or indentations.

    If a comment block (either namelist-level or per-key full-line comment) was preceded by a blank line in the original file,
    a single blank line is written before that comment block.
    """
    lines = []
    for namelist_name, namelist_values in namelist_data.items():
        comment = namelist_values.get("comment")
        comment_leading_blank = namelist_values.get("comment_leading_blank", False)

        # Namelist-level comment (possibly multi-line)
        if comment and len(comment) != 0:
            if comment_leading_blank:
                lines.append("")  # restore the blank line before the comment
            lines.append("!{}".format(comment.replace("\n", "\n!")))

        lines.append("&{}".format(namelist_name))

        for key, value_data in namelist_values.items():
            if key in ("comment", "comment_leading_blank"):
                continue  # handled above

            # Handle per-key full-line comment
            if isinstance(value_data, dict) and "full_line_comment" in value_data:
                full_line_comment = value_data.get("full_line_comment", "")
                if len(full_line_comment) != 0:
                    leading_blank = value_data.get("full_line_comment_leading_blank", False)
                    if leading_blank:
                        lines.append("")  # restore blank before the comment
                    lines.append("!{}".format(full_line_comment.replace("\n", "\n!")))

            if key in ("full_line_comment", "full_line_comment_leading_blank"):
                continue  # internal metadata, not an actual namelist key

            if not isinstance(value_data, dict) or "value" not in value_data:
                continue  # skip any non-standard entries

            # Key/value line
            lines.append("  {} = {}".format(key, value_data["value"]))

            inline_comment = value_data.get("inline_comment")
            if inline_comment and len(inline_comment) != 0:
                lines[-1] = lines[-1] + (
                    " !{}".format(inline_comment.replace("\n", ""))
                )

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
