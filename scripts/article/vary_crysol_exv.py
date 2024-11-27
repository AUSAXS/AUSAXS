import math
import subprocess
import re
import os

def parse_default_data(filename):
    """Parse the default data file and return a dictionary of element data."""
    element_data = {}
    with open(filename, 'r') as f:
        lines = f.readlines()

    current_element = None
    data_block = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue  # Skip empty lines
        if line.startswith('data_'):
            if current_element and data_block:
                element_data[current_element] = data_block
            current_element = line[5:]  # Get the element id after 'data_'
            data_block = {}
        elif line.startswith('#'):
            if current_element and data_block:
                element_data[current_element] = data_block
            current_element = None
            data_block = {}
        elif current_element:
            # Split line into key and value
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
            key = key.strip()
            value = value.strip()
            # Remove '_atsas_atomic_group.' prefix from key
            if key.startswith('_atsas_atomic_group.'):
                key = key[len('_atsas_atomic_group.'):]
            # Remove surrounding quotes from value if any
            if value.startswith("'") and value.endswith("'"):
                value = value[1:-1]
            # Convert numerical values to appropriate types
            if key in ['displaced_volume', 'radius_vdw', 'number_of_electrons', 'weight', 'scattering_length_bound_coherent']:
                try:
                    if '.' in value or 'e' in value.lower():
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass  # Keep as string if conversion fails
            data_block[key] = value
    # After loop, make sure to add the last element
    if current_element and data_block:
        element_data[current_element] = data_block
    return element_data

def calculate_radius(volume):
    """Calculate the spherical radius from volume."""
    if volume <= 0:
        return 0.0
    radius = ((3 * volume) / (4 * math.pi)) ** (1/3)
    return radius

def parse_input_file(filename):
    """Parse the input file and return a dictionary of methods and their volume sets."""
    methods = {}
    with open(filename, 'r') as file:
        lines = file.readlines()

    current_method = None
    headers = []
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('###'):
            # Start of a new method
            current_method = line.strip('# ').strip()
            methods[current_method] = []
            headers = []
        elif current_method and line:
            if line.startswith('chi2'):
                # Header line
                headers = line.split()
            else:
                # Data line
                values = line.split()
                if len(values) != len(headers):
                    continue  # Skip malformed lines
                volume_set = {}
                chi2_value = values[0]
                volume_set['chi2'] = float(chi2_value)
                for header, value in zip(headers[1:], values[1:]):
                    volume_set[header] = float(value)
                methods[current_method].append(volume_set)
    return methods

def generate_parameters_file(volume_set, element_data, method_name):
    """Generate the parameters file for CRYSOL based on the volume set."""
    output_lines = []
    for element, volume in volume_set.items():
        if element == 'chi2':
            continue  # Skip chi2 in volume set
        element_info = element_data.get(element, {})
        description = element_info.get('description', 'Unknown')
        number_of_electrons = element_info.get('number_of_electrons', 'Unknown')
        weight = element_info.get('weight', 'Unknown')
        weight_source = element_info.get('weight.source', 'Unknown')
        scattering_length = element_info.get('scattering_length_bound_coherent', 'Unknown')
        scattering_length_units = element_info.get('scattering_length_bound_coherent.units', 'Unknown')

        # Use the volume from the volume set, and calculate radius
        displaced_volume = volume
        radius_vdw = calculate_radius(volume)

        # Skip elements with zero volume or radius
        if displaced_volume <= 0 or radius_vdw <= 0:
            continue

        # Generate the data section
        data_section = f"data_{element}\n"
        data_section += f"_atsas_atomic_group.id                                       {element}\n"
        data_section += f"_atsas_atomic_group.description                              {description}\n"
        data_section += f"_atsas_atomic_group.displaced_volume                         {displaced_volume:.3f}\n"
        data_section += f"_atsas_atomic_group.displaced_volume.units                   angstrom^3\n"
        data_section += f"_atsas_atomic_group.displaced_volume.source                  '{method_name} method'\n"
        data_section += f"_atsas_atomic_group.radius_vdw                               {radius_vdw:.3f}\n"
        data_section += f"_atsas_atomic_group.radius_vdw.units                         angstrom\n"
        data_section += f"_atsas_atomic_group.radius_vdw.source                        '{method_name} method'\n"
        data_section += f"_atsas_atomic_group.number_of_electrons                      {number_of_electrons}\n"
        data_section += f"_atsas_atomic_group.weight                                   {weight}\n"
        data_section += f"_atsas_atomic_group.weight.source                            '{weight_source}'\n"
        data_section += f"_atsas_atomic_group.scattering_length_bound_coherent         {scattering_length}\n"
        data_section += f"_atsas_atomic_group.scattering_length_bound_coherent.units   {scattering_length_units}\n"
        data_section += "#\n"
        output_lines.append(data_section)
    return output_lines

def run_crysol(pdb_file, data_file):
    """Run CRYSOL with the specified parameters and capture the output."""
    # Create the temp/crysol directory if it doesn't exist
    temp_dir = os.path.join('temp', 'crysol')
    os.makedirs(temp_dir, exist_ok=True)

    # Get the base filenames
    pdb_filename = os.path.abspath(pdb_file)
    data_filename = os.path.abspath(data_file)

    command = ['crysol', pdb_filename, data_filename, '-cst']
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=temp_dir)
    if process.returncode != 0:
        print(f"CRYSOL failed with error code {process.returncode}")
        print(process.stderr)
        return None
    
    output = process.stdout
    return output

def extract_chi2(crysol_output):
    """Extract the chi2 value from CRYSOL output."""
    # Updated regex pattern
    pattern = r"Chi-square of fit[\s\.]*:\s*([0-9\.Ee+\-]+)"
    match = re.search(pattern, crysol_output)
    if match:
        try:
            chi2_value = float(match.group(1))
            return chi2_value
        except ValueError:
            print(f"Unable to convert chi2 value '{match.group(1)}' to float.")
            return None
    else:
        print("Chi2 value not found in CRYSOL output.")
        return None

start_index = {}
def main():
    home = os.path.expanduser("~")
    input_filename = 'output/vary_exv_tables/SASDJY3.txt'
    default_data_filename = home+'/tools/ATSAS/share/atsas/data/atomic_groups.cif.original'
    pdb_file = 'data/SASDJY3/SASDJY3.pdb'
    data_file = 'data/SASDJY3/SASDJY3.dat'

    element_data = parse_default_data(default_data_filename)
    methods = parse_input_file(input_filename)

    chi2_values = []

    for method_name, volume_sets in methods.items():
        print(f"Processing method: {method_name}")
        start_index[method_name] = len(chi2_values)
        for idx, volume_set in enumerate(volume_sets):
            if (idx > 100): break
            # Generate the parameters file for this volume set
            parameters_content = generate_parameters_file(volume_set, element_data, method_name)
            parameters_filename = home+'/tools/ATSAS/share/atsas/data/atomic_groups.cif'
            with open(parameters_filename, 'w') as outfile:
                outfile.write('\n'.join(parameters_content))
            # Run CRYSOL with the parameters file
            crysol_output = run_crysol(pdb_file, data_file)
            if crysol_output is None:
                continue  # Skip if CRYSOL failed
            # Extract the chi2 value
            chi2_value = extract_chi2(crysol_output)
            print(f"  Processed volume set {idx+1}/{len(volume_sets)} with chi2 = {chi2_value}")
            if chi2_value is not None:
                chi2_values.append(chi2_value)
            else:
                print(f"Chi2 not found for volume set {idx+1}")

    # Save the chi2 values to a file
    with open('chi2_values.txt', 'w') as f:
        for method_name, idx in start_index.items():
            f.write(f"\n### {method_name}\n")
            for chi2 in chi2_values[idx:]:
                f.write(f"{chi2}\n")
    print("Chi2 values saved to chi2_values.txt")

if __name__ == '__main__':
    main()