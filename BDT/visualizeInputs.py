import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import joblib
import numpy as np
from sklearn.impute import SimpleImputer
import os

timestamp = 'Nov17'
outputDir = f'outputs/{timestamp}'
era = 'UL2016preVFP'

parquet_files = [f for f in os.listdir(outputDir) if f.endswith('AL.parquet') and f.startswith('events_')]
df = pd.concat([pd.read_parquet(os.path.join(outputDir, f)) for f in parquet_files], ignore_index=True)


# parquetFile = f"events_{era}.parquet"
print(len(df))

# # Load the Parquet file
# df = pd.read_parquet(f"{outputDir}/{parquetFile}")

# Sample a fraction of the DataFrame
# df = df.sample(frac=0.1, random_state=42)

# plot the ttbarpz distribution when y value is 0
if 'ttbarpz' in df.columns:
    plt.hist(df[df['y'] == -1]['ttbarpz'], bins=100, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['ttbarpz'], bins=100, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['ttbarpz'], bins=100, color = 'red', label='qq', density=True, histtype='step')
    plt.xlabel('ttbarpz')
    plt.xlim(-3000, 3000)
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/ttbarpz.png")
    # reset the plt
    plt.clf()

# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
if 'FW1' in df.columns:
    plt.hist(df[df['y'] == -1]['FW1'], bins=100, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['FW1'], bins=100, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['FW1'], bins=100, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('FW1')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/FW1.png")
    plt.clf()



if 'FW2' in df.columns:
    plt.hist(df[df['y'] == -1]['FW2'], bins=100, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['FW2'], bins=100, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['FW2'], bins=100, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('FW2')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/FW2.png")
    plt.clf()



if 'FW3' in df.columns:
    plt.hist(df[df['y'] == -1]['FW3'], bins=100, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['FW3'], bins=100, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['FW3'], bins=100, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('FW3')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/FW3.png")
    plt.clf()


if 'FW4' in df.columns:
    plt.hist(df[df['y'] == -1]['FW4'], bins=100, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['FW4'], bins=100, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['FW4'], bins=100, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('FW4')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/FW4.png")
    plt.clf()






if 'nJet' in df.columns:
    plt.hist(df[df['y'] == -1]['nJet'], bins=10, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['nJet'], bins=10, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['nJet'], bins=10, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    plt.xlim(3, 12)
    plt.xlabel('nJet')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/nJet.png")
    plt.clf()



if 'Jet_HT' in df.columns:
    plt.hist(df[df['y'] == -1]['Jet_HT'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Jet_HT'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Jet_HT'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlim(0, 2000)
    plt.xlabel('Jet_HT')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Jet_HT.png")
    plt.clf()




if 'pT_sum' in df.columns:
    plt.hist(df[df['y'] == -1]['pT_sum'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['pT_sum'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['pT_sum'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('pT_sum')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/pT_sum.png")
    plt.clf()



if 'Sxx' in df.columns:
    plt.hist(df[df['y'] == -1]['Sxx'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Sxx'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Sxx'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('Sxx')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Sxx.png")
    plt.clf()




if 'Syy' in df.columns:
    plt.hist(df[df['y'] == -1]['Syy'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Syy'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Syy'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('Syy')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Syy.png")
    plt.clf()


if 'Szz' in df.columns:
    plt.hist(df[df['y'] == -1]['Szz'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Szz'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Szz'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('Szz')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Szz.png")
    plt.clf()



if 'Sxy' in df.columns:
    plt.hist(df[df['y'] == -1]['Sxy'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Sxy'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Sxy'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('Sxy')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Sxy.png")
    plt.clf()




if 'Sxz' in df.columns:
    plt.hist(df[df['y'] == -1]['Sxz'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Sxz'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Sxz'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('Sxz')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Sxz.png")
    plt.clf()



if 'Syz' in df.columns:
    plt.hist(df[df['y'] == -1]['Syz'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['Syz'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['Syz'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('Syz')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/Syz.png")
    plt.clf()

if 'S' in df.columns:
    plt.hist(df[df['y'] == -1]['S'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['S'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['S'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('S')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/S.png")
    plt.clf()

if 'A' in df.columns:
    plt.hist(df[df['y'] == -1]['A'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['A'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['A'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('A')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/A.png")
    plt.clf()

if 'AL' in df.columns:
    plt.hist(df[df['y'] == -1]['AL'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['AL'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['AL'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('AL')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/AL.png")
    plt.clf()


if 'P' in df.columns:
    plt.hist(df[df['y'] == -1]['P'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['P'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['P'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('P')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/P.png")
    plt.clf()

if 'p2in' in df.columns:
    plt.hist(df[df['y'] == -1]['p2in'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['p2in'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['p2in'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('p2in')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/p2in.png")
    plt.clf()

if 'p2out' in df.columns:
    plt.hist(df[df['y'] == -1]['p2out'], bins=20, color='blue', label='qg', density=True, histtype='step')
    plt.hist(df[df['y'] == 1]['p2out'], bins=20, color = 'black', label='gg', density=True, histtype='step')
    plt.hist(df[df['y'] == 0]['p2out'], bins=20, color = 'red', label='qq', density=True, histtype='step')
    # set x-limit
    # plt.xlim(-200000, 200000)
    plt.xlabel('p2out')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()
    # save the figure
    plt.savefig(f"{outputDir}/p2out.png")
    plt.clf()
