###
# File created by Leonardo Cencetti on 7/7/20
###

import random
import numpy as np
import progressbar
import pandas as pd
import matplotlib.pyplot as plt

list_size = 10000
mu = 0
sigma = 10
num_bin = 20

a = dict(
    K_P=np.random.normal(mu, sigma, 10000),
    K_ETA=np.random.normal(mu, sigma, 10000),
    weight=np.random.uniform(0, 1, 10000))

event_df = pd.DataFrame(a)

lower_bounds = []
# sort dataframe
index_column = 'K_P'
value_column = 'weight'
event_df.hist(column=[index_column], bins=num_bin, weights=event_df[value_column])
plt.show()
mean = event_df[value_column].sum() / num_bin
sorted_events = event_df.sort_values(by=[index_column], ignore_index=True)


bins = [pd.DataFrame() for _ in range(num_bin)]

i = 0
for b in progressbar.progressbar(range(num_bin)):
    while True:
        try:
            bins[b] = bins[b].append(sorted_events.loc[i], ignore_index=True)
            i += 1
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
