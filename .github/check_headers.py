import os
import sys

directories = ["source", "include", "executable", "scripts"]
suffixes = (".cpp", ".h", ".py")

flag_fail = False
for directory in directories:
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(suffixes):
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as f:
                    line = f.readline().strip()
                    if not "SPDX-License-Identifier: LGPL-3.0-or-later" in line:
                        print(f"Missing or imporperly formatted license header in {file_path}")
                        flag_fail = True
                    
                    line = f.readline()
                    tokens = line.split()
                    if len(tokens) < 3 or tokens[0] not in ["//", "#"] or tokens[1] != "Author:":
                        print(f"Missing or improperly formatted author line in {file_path}")
                        flag_fail = True

if flag_fail:
    print("Some files are missing the license header or author line.")
    exit(1)

print("All files have a license header and author line.")
exit(0)