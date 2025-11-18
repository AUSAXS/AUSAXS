import os

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
                    this_fail = False
                    if not "SPDX-License-Identifier: LGPL-3.0-or-later" in line:
                        flag_fail = True
                        this_fail = True

                    line = f.readline()
                    tokens = line.split()
                    if len(tokens) < 3 or tokens[0] not in ["//", "#"] or tokens[1] != "Author:":
                        flag_fail = True
                        this_fail = True

                    if this_fail:
                        print(f"Check license header and author line in {file_path}")

if flag_fail:
    print("Some files are missing the license header or author line.")
    exit(1)

print("All files have a license header and author line.")
exit(0)