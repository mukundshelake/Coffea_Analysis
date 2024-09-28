import uproot
import awkward as ak
import json
import pandas as pd
for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
    with open(f'../../Datasets/sampleFiles_{era}.json', 'r') as json_file:
        dicti = json.load(json_file)
        for dataMC in dicti:
            if 'Data' in dataMC:
                continue
            for pr in dicti[dataMC]:
                datasetName = f'{era}_{pr}'
                for file in dicti[dataMC][pr]:
                    # print(f"Transferring {file}")
                    # Using pexpect to run the scp command and handle password prompt
                    dataFile = uproot.open(file)
                    tree = dataFile["Events"]
                    branch = tree["GenPart_pdgId"].array()
                    third_element = branch[:, 2]
                    fourth_element = branch[:, 3]

                    # Combine 3rd and 4th elements into pairs
                    df = pd.DataFrame({"third": ak.to_numpy(third_element), "fourth": ak.to_numpy(fourth_element)})

                    unique_pairs_df = df.drop_duplicates()
                    # Find unique pairs
                    unique_pairs = ak.Array(unique_pairs_df.to_records(index=False))

                    # Print the unique pairs
                    print(f"for {era}/{pr} the unique pairs are:")
                    for pair in unique_pairs:
                        print(pair)
                    # print(unique_pairs)

