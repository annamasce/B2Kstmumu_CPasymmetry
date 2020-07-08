import ROOT
import matplotlib
ROOT.gROOT.SetBatch(True)
import numpy as np
import progressbar
import pandas as pd
import matplotlib.pyplot as plt

# Jpsi_data = '/vols/lhcb/amascell/rootFiles/data/B2KstJpsi_run2_sWeight_newBDT.root'
Jpsi_data = '/home/anna/master_thesis/data/withBDT/B2KstJpsi_2016_sWeight_newBDT.root'


fileName = Jpsi_data
treeName = 'DecayTree'
d = ROOT.RDataFrame(treeName, fileName)
print('Applying Jpsi weights')
d_weights = d.Define('weight', 'sWeight')

dict_K = d_weights.AsNumpy(['K_P', 'K_ETA', 'weight'])

print(dict_K, type(dict_K))

num_bin = 30

a = dict(
    K_P=dict_K['K_P'],
    K_ETA=dict_K['K_ETA'],
    weight=dict_K['weight'])

event_df = pd.DataFrame(a)

# sort dataframe
index_column = 'K_P'
value_column = 'weight'
event_df.hist(column=[index_column], bins=num_bin, weights=event_df[value_column])
plt.show()
mean = event_df[value_column].sum() / num_bin
sorted_events = event_df.sort_values(by=[index_column], ignore_index=True)


bins = [pd.DataFrame() for _ in range(num_bin)]
bar = progressbar.ProgressBar(max_value=len(event_df))

i = 0
for b in range(num_bin):
    while True:
        try:
            bins[b] = bins[b].append(sorted_events.loc[i], ignore_index=True)
            i += 1
            bar.update(i)
            if bins[b][value_column].sum() > mean:
                break
        except KeyError:
            break

sums = [el[value_column].sum() for el in bins]
low_limits = [el[index_column][0] for el in bins]
low_limits.append(float(bins[-1].tail(1)[index_column]))
print('mean: {} sums:'.format(round(mean, 3)), sums)
print(low_limits)
event_df.hist(column=[index_column], bins=low_limits, weights=event_df[value_column])
plt.show()
