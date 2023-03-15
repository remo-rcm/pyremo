output_file_patterns = {
    "e": "e{expid:06d}{type}{middle}{date}{end}.{suffix}",
    "t": "e{expid:06d}{type}{date}{end}.{suffix}",
}


def get_file_pattern(
    type, expid=None, date=None, code=None, middle=None, end=None, suffix="nc"
):
    wild = "*"
    if middle is None:
        middle = ""
    if expid is None:
        expid = wild
    else:
        expid = f"{int(expid):06d}"
    if type in ["e", "n"] and code is None:
        middle = "_*_"
    elif type in ["e", "n"]:
        middle = f"_c{int(code):03d}_"
    else:
        middle = ""
    if end is None:
        end = wild
    if date is None:
        date = wild
    return "e{expid}{type}{middle}{date}{end}.{suffix}".format(
        expid=expid, type=type, middle=middle, date=date, end=end, suffix=suffix
    )
