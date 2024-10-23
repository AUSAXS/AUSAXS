import sys
import os

def generate_header(input_files, output_header):
    header_content = "#pragma once\n\n"
    header_content += "#include <plotting.cpp>\n"
    header_content += "#include <array>\n\n"

    for input_file in input_files:
        variable_name = os.path.basename(input_file).replace(".", "_")
        with open(input_file, 'rb') as file:
            content = file.read()
            array_content = ", ".join(f"{byte}" for byte in content)
            
            header_content += (
                f"inline constexpr std::array<char, {len(content)}> resources::{variable_name} = {{\n"
                f"    {array_content}\n"
                f"    }};\n\n"
            )

    with open(output_header, 'w') as header_file:
        header_file.write(header_content)
    print(f"Header file '{output_header}' generated successfully.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python generate_header.py <output_header> <input_file1> <input_file2> ...")
        sys.exit(1)

    output_header = sys.argv[1]
    input_files = sys.argv[2:]
    generate_header(input_files, output_header)