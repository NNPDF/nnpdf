import yaml

# Input and output filenames
input_file = "kinematics_XDQ.yaml"
output_file = "kinematics_XUB.yaml"

# New x mid values
new_x_mids = [
    0.01, 0.0129155, 0.01668101, 0.02154435, 0.02782559, 0.03593814,
    0.04641589, 0.05994843, 0.07742637, 0.1, 0.18, 0.26, 0.34, 0.42,
    0.5, 0.58, 0.66, 0.74, 0.82, 0.9
]

# Load YAML
with open(input_file, "r") as f:
    data = yaml.safe_load(f)

# Modify Q2 and x mid values
for i, bin_entry in enumerate(data["bins"]):
    # 1) Replace Q2 mid with 10000
    bin_entry["Q2"]["mid"] = 10000.0
    # 2) Replace x mid with the corresponding new value (if available)
    if i < len(new_x_mids):
        bin_entry["x"]["mid"] = new_x_mids[i]

# Save new YAML
with open(output_file, "w") as f:
    yaml.dump(data, f, sort_keys=False)

print(f"✅ Modified YAML saved to {output_file}")