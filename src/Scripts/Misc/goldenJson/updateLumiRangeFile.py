import json

# File paths
input_file = "lumi_range_EIGHTEEN.json"
output_file = "data_range_EIGHTEEN.json"

# Load the input JSON file
with open(input_file, "r") as infile:
    lumi_data = json.load(infile)

# Transform the data
converted_data = {}
for run, ranges in lumi_data.items():
    converted_ranges = []
    for lumi_range in ranges:
        if "-" in lumi_range:
            start, end = map(int, lumi_range.split("-"))
            converted_ranges.append([start, end])
        else:
            # Single lumi section (e.g., "1830")
            lumi = int(lumi_range)
            converted_ranges.append([lumi, lumi])
    converted_data[run] = converted_ranges

# Save the transformed data to the output JSON file
with open(output_file, "w") as outfile:
    json.dump(converted_data, outfile, indent=2)

print(f"Converted lumi range saved to {output_file}")