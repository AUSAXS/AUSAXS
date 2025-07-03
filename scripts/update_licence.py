import os
import sys

directory = "include"
if len(sys.argv) == 2:
    directory = sys.argv[1]

default_license = \
"""\
// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje\
"""

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith(".cpp") or file.endswith(".h"): # only cpp files
            file_path = os.path.join(root, file)
            with open(file_path, 'r+') as f:
                content = f.read()

                if "License" in content:
                    # delete the old license - first 4 lines
                    content = content.split('\n')
                    content = content[4:]
                    content = '\n'.join(content)

                    # check for empty lines at the beginning
                    i = 0
                    while i < len(content) and len(content[i].split()) == 0:
                        i += 1

                    content = content[i:]

                # insert the new license
                content = default_license + '\n\n' + content

                # write at the beginning of the file
                f.seek(0)
                f.write(content)
                f.truncate()