def parse(file_content):
    namelists = {}
    current_namelist = None
    full_line_comment = None
    for line in file_content.splitlines():
        line = line.strip()
        if line.startswith("!"):
            full_line_comment = line[1:].strip()
        elif line.startswith("&"):
            current_namelist = line[1:].strip()
            namelists[current_namelist] = {"comment": full_line_comment}
            full_line_comment = None
        elif line.startswith("/"):
            current_namelist = None
        elif current_namelist is not None:
            line_parts = line.split("!", maxsplit=1)
            key_value_pair = line_parts[0].strip()
            inline_comment = line_parts[1].strip() if len(line_parts) == 2 else None
            if not key_value_pair:
                continue
            key, value = key_value_pair.split("=")
            key = key.strip()
            value = value.strip()
            namelists[current_namelist][key] = {"value": value, "full_line_comment": full_line_comment, "inline_comment": inline_comment}
            full_line_comment = None
    return namelists

def write(namelist_data):
    lines = []
    for namelist_name, namelist_values in namelist_data.items():
        lines.append("&{}".format(namelist_name))
        comment = namelist_values.get("comment")
        if comment:
            lines.append("! {}".format(comment))
        for key, value_data in namelist_values.items():
            if key == "comment":
                continue
            lines.append("  {} = {}".format(key, value_data["value"]))
            inline_comment = value_data.get("inline_comment")
            if inline_comment:
                lines.append("  ! {}".format(inline_comment))
        lines.append("/")
    return "\n".join(lines)