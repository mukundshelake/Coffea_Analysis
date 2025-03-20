import argparse
import json
import os


def create_output_path(input_file, output_directory):
    # Extract the last 4 subfolders from the input file path
    subfolders = os.path.normpath(input_file).split(os.sep)[-4:-1]
    # Create the new output file path
    output_file_path = os.path.join(output_directory, *subfolders, os.path.splitext(os.path.basename(input_file))[0] + '_output.h5')
    return output_file_path

parser = argparse.ArgumentParser(description="Process some eras.")

# Define the allowed choices
allowed_eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']

# Add the --eras argument with choices and default value
parser.add_argument(
    '-e', '--eras', 
    choices=allowed_eras, 
    nargs='*', 
    default=allowed_eras,
    help="Specify one or more eras. Allowed values are: 'UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018'. If not provided, all eras will be used by default."
)

parser.add_argument(
    '-c', '--channels',
    nargs='+',  # Accept one or more values
    type=str,
    default=[],  # Default to an empty list if no arguments are provided
    help="Specify one or more channels to process. By default, all channels"
)


# Add the --sample flag argument
parser.add_argument(
    '-s', '--sample',
    action='store_true',
    help="If provided, the sample mode will be enabled."
)

parser.add_argument(
    '--skimmed',
    action='store_true',
    help="If provided, skimmed datasets will be used for analysis."
)

parser.add_argument(
    '--onlyData',
    action='store_true',
    help="If provided, only the data will be processed."
)

parser.add_argument(
    '--onlyMC',
    action='store_true',
    help="If provided, only the MC will be processed."
)

# Add the --output argument
parser.add_argument(
    '-t','--timestamp',
    type=str,
    default='timestamp',
    help="Specify the timestamp'."
)

args = parser.parse_args()

# Display the parsed arguments
print(f"Selected eras: {args.eras}")
print(f"Sample mode: {args.sample}")
print(f"Timestamp for book keeping: {args.timestamp}")

if len(args.channels) > 0:
    print(f"Channels: {args.channels}")
else:
    print("Channels: All")


outputDir = f"outputs/{args.timestamp}"
datasetFlag = 'data'

if args.sample:
    datasetFlag = 'sample'

if args.skimmed:
    datasetFlag = 'skimmed_' + datasetFlag

fileset = {}
for era in args.eras:
    with open(f'../../Datasets/{datasetFlag}Files_{era}.json', 'r') as json_file:
        dicti = json.load(json_file)
        for pr in dicti['Data_mu']:
            if args.onlyMC:
                print("Working on only the MC, ignoring data")
                continue
            if len(args.channels) > 0:
                if pr not in args.channels:
                    # print(f'skipping {era}_{pr} as not in list')
                    continue
            datasetName = f'{era}_{pr}'
            fileset[datasetName] = {"files": dicti['Data_mu'][pr]}
        for pr in dicti['MC_mu']:
            if args.onlyData:
                print("Working on only the data, ignoring MC")
                continue
            if len(args.channels) > 0:
                if pr not in args.channels:
                    # print(f'skipping {era}_{pr} as not in list')
                    continue
            datasetName = f'{era}_{pr}'
            fileset[datasetName] = {"files": dicti['MC_mu'][pr]}

output_directory = os.path.join('outputs', args.timestamp)
os.makedirs(output_directory, exist_ok=True)

input_files = []

for dataset in fileset:
    input_files.extend(fileset[dataset]['files'])

for input_file in input_files:
    output_file = create_output_path(input_file, output_directory)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    process_file(input_file, output_directory)