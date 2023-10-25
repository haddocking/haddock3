"""haddock3-restraints validate_tbl subcommand.

The validate_tbl subcommand validates an input TBL file.
Usage:
    haddock3-restraints validate_tbl <tbl_file> [--silent] [--quick]
"""
from pathlib import Path

import re


def add_validate_tbl_arguments(validate_tbl_subcommand):
    """Add arguments to the score subcommand."""
    validate_tbl_subcommand.add_argument(
        "tbl_file",
        type=str,
        help="TBL file to be validated",
        )

    validate_tbl_subcommand.add_argument(
        "--pcs",
        help="PCS mode",
        action='store_true',
        )

    validate_tbl_subcommand.add_argument(
        "--quick",
        help="Check global formatting before going line by line "
        "(opening/closing parenthesis and quotation marks",
        action='store_true',
        )
    
    validate_tbl_subcommand.add_argument(
        "--silent",
        help="Only output errors, do not output TBL file at the end",
        action='store_true',
        )

    return validate_tbl_subcommand


def check_parenthesis(file):
    open_parenthesis = re.compile('[(]')
    close_parenthesis = re.compile('[)]')
    quotation_marks = re.compile('[\"]')
    opened = 0
    closed = 0
    quote = 0
    for _match in open_parenthesis.finditer(file):
        opened += 1
    for _match in close_parenthesis.finditer(file):
        closed += 1
    for _match in quotation_marks.finditer(file):
        quote += 1
    if opened != closed:
        raise Exception("Problem with TBL file parentheses ({:d} opening for {:d} "
                        "closing parentheses)".format(opened, closed))
    if quote % 2 != 0:
        raise Exception("Problem with TBL file, odd number of quotation marks "
                        "({:d} quotation marks)".format(quote))


def validate_tbldata(restraints, pcs=False):
    # List of selection keywords
    selectors = ('name', 'resn', 'atom', 'resi', 'attr', 'segi', 'chem', 'id',
                 'byres', 'not')
    connectors = ('and', 'or', 'not', 'byres')
    output = ""
    parentmatch = re.compile('[()]')
    # Global mode is activated outside any assign statement
    mode = "global"
    # Remove any carriage return/new line caracters
    lines = restraints.replace('\r', '').split("\n")
    # Line number
    lnr = 0
    # Temporary line storage for future output
    tmp_output = None
    tmp = None
    level = None
    s = None
    selections = None
    postselections = None
    numbers = None
    types = None
    lastassign = None
    for line in lines:
        lnr += 1
        # Take everything which is before putative "!" (comment) caracter
        if line.find("!") > -1:
            line = line[:line.find("!")]
        # Remove whitespaces and merge with previous line if in OR statement
        # for AIR restraints
        if mode != "format":
            line = line.strip()
        else:
            line = tmp + line.strip()
            mode = "postassign"
        # End of line
        if not len(line):
            continue
        # Check if all "" are closed
        if line.count('"') % 2 > 0:
            raise Exception('Unclosed " at line {:d} for line: {}'.format(lnr, line))
        if mode in ("global", "postglobal"):
            # Start of an assign statement
            if line.lower().startswith("assi"):
                if mode == "postglobal":
                    output += "\n!\n"
                    output += tmp_output
                mode = "assign"
                selections = []
                # assign is the only character of the line,check following lines
                if line.find("(") == -1:
                    continue
                # Reset temporary buffer
                tmp_output = None
            # Check for "OR" restraint
            elif mode == "postglobal" and line.lower().startswith("or"):
                mode = "postassign"
                selections = []
                line = line[len("or"):]
                # Case where "OR" statement is the only one on the line
                # (rest of selection at the next line)
                if line == "":
                    continue
            # We are not treating an assign params (postglobal) or at the
            # end (global) and no "OR" restraint is present -> ERROR
            else:
                raise Exception("Invalid TBL file: Unknown statement (line {:d}): {}".format(lnr, line))
        # We check if the selection is made over two lines
        # (thanks to "segid" keyword)
        if mode == "postassign":
            if line.count("segid") == 1:
                mode = "format"
                tmp = line
                continue

        matched = True
        while matched:
            matched = False
            # We are looking for parenthesis as start of the selections
            if mode in ("assign", "postassign"):
                # Ambiguous restraint for two different pairs of atoms
                # (ex: THR 1 B <-> ALA 2 A OR GLY 2 B <-> ASP 10 A )
                pos = line.find("(")
                if pos != -1:
                    matched = True
                    line = line[pos + 1:]
                    # We look for an "OR" selection
                    if mode == "postassign":
                        mode = "postsel"
                    else:
                        mode = "sel"
                    lastassign = lnr
                    s = ""
                    level = 1
            # Get the structural selections
            if mode in ("sel", "postsel"):
                # Detect opening and closing parenthesis
                for match in parentmatch.finditer(line):
                    if match.group() == "(":
                        level += 1
                    else:
                        level -= 1
                    # End of the parenthesis content
                    if level == 0:
                        if mode == "postsel":
                            mode = "postassign"
                        else:
                            mode = "assign"
                        matched = True
                        # print l[:match.start()]
                        # Get back the parentheses content and check
                        # for selection keywords
                        idx_connectors = []
                        sel = line[:match.start()]
                        # We can directly check for the 1st keyword presence
                        if len(sel) > 0:
                            syntax_ok = False
                            for s in selectors:
                                if sel.strip().startswith(s):
                                    syntax_ok = True
                            if not syntax_ok and not \
                                    sel.strip().startswith('byres') and not \
                                    sel.strip().startswith('not'):
                                raise Exception("1) Missing or wrong keyword in-term {} (stopped at line {:d})".
                                                format(sel, lnr))
                        # Get indexes of connectors (AND, OR, etc.)
                        for c in connectors:
                            idx_connectors.append([m.start()+len(c) for m in
                                                   re.finditer(c, sel)])
                        # Flatten the list
                        tmp_list = [item for sublist in idx_connectors for
                                    item in sublist]
                        idx_connectors = tmp_list
                        # Check that each connector is followed by a keyword
                        for i in idx_connectors:
                            new_sel = sel[i:].strip().strip("(").strip()
                            if len(new_sel) > 0:
                                syntax_ok = False
                                for s in selectors:
                                    if new_sel.startswith(s):
                                        syntax_ok = True
                                        break
                                if not syntax_ok:
                                    raise Exception("Missing or wrong keyword in-term {} (stopped at line {:d})".
                                                    format(sel, lnr))
                        s += sel
                        selections.append(s)
                        s = None
                        # Get the rest of the line
                        line = line[match.end():]
                        # Go to the process of the selection
                        break
                # No parenthesis, we get the whole line
                if level > 0:
                    if line != "":
                        s += line + "\n\t"
                    else:
                        # To avoid multiple blank lines
                        if repr(s) == "''":
                            # print repr(l), repr(s), selections
                            s += "\n\t"
        # TBD
        if mode in ("sel", "postsel"):
            continue
        # Selection parsed for "OR" line
        if mode == "postassign":
            if len(selections) != postselections:
                raise Exception("Invalid TBL file: wrong number of selections: in-term {:d}, cross-term {:d} "
                                "(stopped at line {:d})".format(postselections, len(selections), lnr))
            tmp_output += " or"
            for s in selections:
                tmp_output += "\t({})\n".format(s)
            # We let the possibility for other "OR"
            mode = "postglobal"
        if len(line) == 0:
            continue
        # Define the distance restraints type to adapt the parsing
        if mode == "assign":
            mode = "numbers"
            if pcs:
                if len(selections) == 5:
                    types = (" {:.3f}", " {:.3f}")
                else:
                    raise Exception("Invalid TBL file: wrong number of selections (must be 5 in PCS mode)")
            else:
                if len(selections) == 2:
                    types = (" {:.3f}", " {:.3f}", " {:.3f}")
                elif len(selections) == 4:
                    types = (" {:.3f}", " {:.3f}", " {:.3f}", " {:d}")
                elif len(selections) == 5:
                    raise Exception("Invalid TBL file: wrong number of selections (can be 5 only in PCS mode)")
                elif len(selections) == 6:
                    types = (" {:.3f}", " {:.3f}")
                else:
                    check_parenthesis(restraints)
                    raise Exception("Invalid TBL file: wrong number of selections (must be 2, 4 or 6)")
            postselections = len(selections)
            numbers = []
        # Distance restraints parsing
        if mode == "numbers":
            ll = line.split()
            for num in ll:
                if len(numbers) == len(types):
                    break
                numbers.append(float(num))
            if len(numbers) == len(types):
                tmp_output = "assign "
                for s in selections:
                    tmp_output += "\t({})\n".format(s)
                tmp_output = tmp_output[:-len("\n")]
                for n, t in zip(numbers, types):
                    tmp_output += t.format(n)
                tmp_output += "\n"
                mode = "postglobal"
    # "OR" lines have been parsed, we store the selections
    if mode == "postglobal":
        output += "!\n"
        output += tmp_output
        mode = "global"
    # If mode is not back to global, something has not been processed properly
    if mode != "global":
        raise Exception("Invalid TBL file: Malformed ASSIGN statement (line {:d}), use --quick to check for "
                        "putative syntax issues".format(lastassign))
    if not len(output.strip()):
        raise Exception("Invalid or empty TBL file")

    # Remove extra lines before each assign statement
    output = output.replace("\n\n", "\n")
    # Remove extra line at the beginning of the file
    if output.startswith("\n"):
        output = output.replace("\n", "", 1)
    return output


def validate_tbl(tbl_file, pcs, quick=False, silent=False):
    """Get the passive residues."""
    if quick:
        tbldata = open(tbl_file).read()
        # Check the parenthesis and quotation marks opening/closure
        check_parenthesis(tbldata)
    
    if Path(tbl_file).exists():
        tbldata = open(tbl_file).read()
        # Parse and process the restraints
        if silent:
            validate_tbldata(tbldata, pcs)
        else:
            print(validate_tbldata(tbldata, pcs))
    else:
        raise Exception(f"TBL file {tbl_file} does not exist, check the path")
    
    return
