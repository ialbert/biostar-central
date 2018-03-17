import re


def count_diffs(diff_list, str_format=False):
    nadd, nsub = 0, 0
    for diff in diff_list:
        if diff.startswith("-") and not diff.startswith("---"):
            nsub += 1
        elif diff.startswith("+") and not diff.startswith("+++"):
            nadd += 1
    if str_format:
        nadd = f"{nadd} line" if nadd == 1 else f"{nadd} lines"
        nsub = f"{nsub} line" if nsub == 1 else f"{nsub} lines"

    return nadd, nsub


def add_color(diff, back, strong):
    return f"<span style='background-color:{back};'><strong style='color:{strong};'>{diff}</strong></span>"


def color_diffs(diff_list):
    colored_diffs = []
    line = 1
    nadd, nsub = count_diffs(diff_list, str_format=True)

    for diff in diff_list:

        if diff.startswith("-"):
            diff = f"{diff.strip()} {nsub}\n" if diff.startswith("---") else diff
            diff = add_color(diff=diff, back="#ff6666", strong="#661400")

        elif diff.startswith("+"):
            diff = f"{diff.strip()} {nadd}\n" if diff.startswith("+++") else diff
            diff = add_color(diff=diff, back="#99ff99", strong="#008000")

        elif diff.startswith("@@"):
            # Extract line number from a string like: @@ -17,7 +17,7 @@

            line = re.search(r"\d+,", diff)
            line = int(line.group(0).split(",")[0]) + 1
            colored_diffs.append("\n" + diff + "..........\n")
            continue

        # Add the line number
        color_line = f"<em style='color: gray;'>{line}</em>\t"
        colored_diffs.append(color_line + diff)
        line += 1

    return colored_diffs
