#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')


# create test data
np.random.seed(19680801)
data = [sorted(np.random.normal(0, std, 100)) for std in range(1, 5)]
print(data)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), sharey=True)
quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)

inds = np.arange(1, len(medians) + 1)

ax1.set_title('Default violin plot')
ax1.set_ylabel('Observed values')
ax1.violinplot(data, showmeans=True, showmedians=False)
ax1.scatter(inds, medians, color='white', s=5, zorder=3)
ax1.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)

quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)

inds = np.arange(1, len(medians) + 1)

ax2.set_title('Default violin plot')
ax2.set_ylabel('Observed values')
ax2.violinplot(data, showmeans=True, showmedians=False)
ax2.scatter(inds, medians, color='white', s=5, zorder=3)
ax2.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)


labels = ['A', 'B', 'C', 'D']
for ax in [ax1, ax2]:
    set_axis_style(ax, labels)

plt.subplots_adjust(bottom=0.15, wspace=0.05)
plt.show()