"""PDB Input/Output"""


_to_remove = ['REMAR', 'CTERB', 'CTERA', 'NTERA', 'NTERB', 'CONECT']
_to_rename = {'WAT ':'TIP3', 'HSD':'HIS', 'HSE':'HIS', 'HID':'HIS', 'HIE':'HIS',
              ' 0.00969':' 0.00   '}

def sanitize(input_file_name, output_file_name):
    """Remove problematic portions of a PDB file."""
    with open(output_file_name, 'w') as output_handler:
        with open(input_file_name) as input_handler:
            for line in input_handler:
                # Ignoring lines containing any tag from _to_remove
                if not any([tag in line for tag in _to_remove]):
                    for tag, new_tag in _to_rename.items():
                        line = line.replace(tag, new_tag)
                    output_handler.write(line)

