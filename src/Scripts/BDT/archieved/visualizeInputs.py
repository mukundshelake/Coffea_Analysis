import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import joblib
import numpy as np
from sklearn.impute import SimpleImputer
import os

timestamp = 'Nov11'
outputDir = f'outputs/{timestamp}'
era = 'UL2016preVFP'
parquetFile = f"events_{era}.parquet"

# Load the Parquet file
df = pd.read_parquet(f"{outputDir}/{parquetFile}")

# Sample a fraction of the DataFrame
# df = df.sample(frac=0.1, random_state=42)

# plot the ttbarpz distribution when y value is 0
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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


# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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






# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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




# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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




# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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




# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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



# plot similar for FW1 binned between -1 * 10^6 to 1 * 10^6 with 100 bins
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