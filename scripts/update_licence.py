import os
import sys

directory = "source"
if len(sys.argv) == 2:
    directory = sys.argv[1]

default_license = \
"""/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/"""

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith(".cpp"):
            file_path = os.path.join(root, file)
            with open(file_path, 'r+') as f:
                content = f.read()

                if "License" in content:
                    continue

                content = default_license + '\n\n' + content

                f.seek(0)
                f.write(content)