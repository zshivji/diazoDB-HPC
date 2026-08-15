import matplotlib.pyplot as plt
import numpy as np

category_names = ['I', 'II', 'III', 'III-anf', 'III-vnf', 'IVa']

results = {'DiazoDB\n7911': np.array([2064, 4178, 338, 271, 216, 168]), 
            'NFix\n4488': np.array([4085, 0, 0, 261, 56, 0]),
            'NSDB\n1536': np.array([139, 181, 31, 17, 17, 0]), 
            'Cyano\n1342': np.array([1310, 9, 0, 23, 0, 0])}


labels = list(results.keys())
data = np.array(list(results.values()))
data_cum = data.cumsum(axis=1)
category_colors = plt.colormaps['RdYlGn'](
    np.linspace(0.15, 0.85, data.shape[1]))[::-1]

fig, ax = plt.subplots(figsize=(9.2, 3))
ax.invert_yaxis()
ax.xaxis.set_visible(False)
ax.set_xlim(0, np.sum(data, axis=1).max())

for i, (colname, color) in enumerate(zip(category_names, category_colors)):
    widths = data[:, i]
    starts = data_cum[:, i] - widths
    # reduce spacing between bars
    rects = ax.barh(labels, widths, left=starts, height=0.7,
                    label=colname, color=color)

    r, g, b, _ = color
    text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
    # ax.bar_label(rects, label_type='center', color=text_color)
ax.legend(ncols=len(category_names), bbox_to_anchor=(0, 1),
            loc='lower left', fontsize='small')

fig.savefig('../diazoDB-comparison/diazoDB-comp-counts.svg', format='svg',
            dpi=6000, transparent=True)
plt.show()
