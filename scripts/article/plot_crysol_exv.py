import matplotlib.pyplot as plt

# Read Data
with open('chi2_values.txt') as f:
    data = f.read().splitlines()

# Initialize variables
methods = []
current_method = None
current_chi_values = []

# Parse the data
for line in data:
    line = line.strip()
    if line.startswith('###'):
        # Save the previous method's data
        if current_method is not None:
            methods.append((current_method, current_chi_values))
        # Start a new method
        current_method = line.strip('#').strip()
        current_chi_values = []
    elif line:
        try:
            chi_value = float(line)
            current_chi_values.append(chi_value)
        except ValueError:
            pass  # Ignore non-numeric lines
# Add the last method's data
if current_method is not None:
    methods.append((current_method, current_chi_values))

# Plot the density histogram
fig, ax = plt.subplots()
for method, chi_values in methods:
    ax.hist(chi_values, bins=20, alpha=0.5, label=method, density=True)
ax.set_xlabel('$\chi^2_r$')
ax.set_ylabel('Density')
ax.legend()
plt.show()
