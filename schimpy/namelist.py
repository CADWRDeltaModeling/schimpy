def parse(file_content):
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
            namelists[current_namelist][key] = {"value": value, "full_line_comment": full_line_comment, "inline_comment": inline_comment}
            full_line_comment = ""
    return namelists

def write(namelist_data):
    lines = []
    for namelist_name, namelist_values in namelist_data.items():
        comment = namelist_values.get("comment")
        if comment and len(comment) != 0:
            lines.append("!{}".format(comment.replace("\n","\n!")))
        lines.append("&{}".format(namelist_name))
        for key, value_data in namelist_values.items():
            if key == "comment": # already taken care of above
                continue 
            if "full_line_comment" in value_data: # do this first
                full_line_comment = value_data.get("full_line_comment")
                if len(full_line_comment) != 0:
                    lines.append("!{}".format(full_line_comment.replace("\n","\n!")))
            if key == "full_line_comment":
                continue # took care of it above
            lines.append("  {} = {}".format(key, value_data["value"]))
            inline_comment = value_data.get("inline_comment")
            if inline_comment and len(inline_comment) != 0 :
                lines[-1] = lines[-1] + (" !{}".format(inline_comment.replace('\n','')))
        lines.append("/")
    return "\n".join(lines)