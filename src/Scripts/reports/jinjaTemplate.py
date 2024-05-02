from jinja2 import Environment, FileSystemLoader
from coffea.util import load
import argparse
import importlib


# Start the parser
parser = argparse.ArgumentParser(description="Fill up the unfilled tex files")


# Add arguments to the parser
parser.add_argument('-i', '--input', type=str, required=True, help="Input the input tex file")
parser.add_argument('-o', '--output', type=str, required=True, help="Input the output tex file")
parser.add_argument('-c', '--config', type=str, required=True, help="config python file containing data dictionary")



args = parser.parse_args()

inputFile = args.input
outputFile = args.output
configFile = args.config

module = importlib.import_module(configFile)
data = module.Data

# Configure Jinja with custom block start and end strings
env = Environment(loader=FileSystemLoader('.'), variable_start_string='((*', variable_end_string='*))', block_start_string='((%', block_end_string='%))', comment_start_string='<!--', comment_end_string='-->')

# print(inputFile, outputFile, configFile)

template = env.get_template(inputFile)


# Render the template with data
output = template.render(data)

# Save the rendered LaTeX code
with open(outputFile, 'w') as f:
    f.write(output)
