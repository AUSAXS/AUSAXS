import os

directories = ["source", "include", "executable", "scripts"]
suffixes = (".cpp", ".h")

# the only files allowed to touch <stdexcept>: they define the ausaxs::except:: hierarchy itself
allowlist = [
    os.path.join("include", "core", "utility", "Exceptions.h"),
    os.path.join("source", "core", "utility", "Exceptions.cpp"),
    os.path.join("source", "gpu", "SYCL", "SYCLBackend.cpp"),
]

flag_fail = False
for directory in directories:
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(suffixes):
                file_path = os.path.join(root, file)
                if file_path in allowlist:
                    continue
                with open(file_path, 'r') as f:
                    for line_no, line in enumerate(f, start=1):
                        if "#include <stdexcept>" in line:
                            flag_fail = True
                            print(f"{file_path}:{line_no}: includes <stdexcept>; throw ausaxs::except::* instead of a raw std:: exception")
                        if "throw std::" in line:
                            flag_fail = True
                            print(f"{file_path}:{line_no}: throws a standard exception; throw ausaxs::except::* instead")

if flag_fail:
    print("Some files bypass the ausaxs::except:: hierarchy.")
    exit(1)

print("No files bypass the ausaxs::except:: hierarchy.")
exit(0)
