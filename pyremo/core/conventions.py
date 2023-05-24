def output_pattern(
    type, expid=None, date=None, code=None, middle=None, end=None, suffix=".nc"
):
    """Create a REMO output file pattern

    Creates a REMO output file pattern from different arguments. If arguments are
    missing, they will be replaced with wildcards.

    Parameters
    ----------
    type :
        Type of output file, should be in ``["e", "t", "f", "g", "n", "p"]``.
    expid : str or int
        Experiment ID.
    date : str
        Date in the filename.
    code : str or int
        Code in e or n files. Will be ignored for other filetypes.
    middle : str
        Some string in the middle of filename. For e and n files this would be the code
        automatically.
    end : str
        Some string in the end of the filename. Could be used for a wildcard.
    suffix : str
        File suffix.

    Returns
    -------
    Filename patter : str

    """
    wild = "*"
    if middle is None:
        middle = ""
    if expid is None:
        expid = wild
    else:
        expid = f"{int(expid):06d}"
    if type in ["e", "n", "p"]:
        if code is None:
            middle = "_c*_"
        else:
            middle = f"_c{int(code):03d}_"
    else:
        middle = ""
    if date is None:
        date = wild
    if type in ["t", "f", "g"] and len(date) < 10 and end is None:
        end = wild
    if type in ["e", "n", "p"] and len(date) < 6 and end is None:
        end = wild
    if end is None:
        end = ""
    return "e{expid}{type}{middle}{date}{end}{suffix}".format(
        expid=expid, type=type, middle=middle, date=date, end=end, suffix=suffix
    )
