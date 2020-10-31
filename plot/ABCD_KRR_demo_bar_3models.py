import pandas as pd
import matplotlib.pyplot as plt

sail_blue = [0, 0.1255, 0.2471, 1]
mint = [0.6784, 0.9373, 0.8196, 1]
grey = [0.9, 0.9, 0.9, 1]
Nbehavior = 36

df = pd.DataFrame(dict(
	allAArandWA_WAbetter = [20, 12, 21],
	allAArandWA_sigdiff = [20+5, 12+20, 21+6],
	allAArandWA_all = [Nbehavior] * 3,
	randWA_WAbetter = [26, 13, 25],
	randWA_sigdiff = [26+2, 13+14, 25+2],
	randWA_all = [Nbehavior] * 3,
	allAA_WAbetter = [17, 11, 18],
	allAA_sigdiff = [17+6, 11+18, 18+8],
	allAA_all = [Nbehavior] * 3))

fig = plt.figure(figsize=(15, 8))
# transparent backgroung
ax = plt.gca()
ax.patch.set_alpha(0)
# invisible right and upper axes
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

# tick labels
plt.xticks([0.85, 1.85, 2.85], ['Predictive COD', 'Pearson\'s correlation', 'MSE'])
ax.tick_params(axis='both', which='major', labelsize=18, direction='out', width=2)

allAArandWA_list = [plt.bar([1, 2, 3], df.allAArandWA_all, align='edge', width=0.2, linewidth=0, color=grey),
	plt.bar([1, 2, 3], df.allAArandWA_sigdiff, align='edge', width=0.2, linewidth=0, color=mint),
	plt.bar([1, 2, 3], df.allAArandWA_WAbetter, align='edge', width=0.2, linewidth=0, color=sail_blue)]

randWA_list = [plt.bar([1-0.25, 2-0.25, 3-0.25], df.randWA_all, align='edge', width=0.2, linewidth=0, color=grey),
	plt.bar([1-0.25, 2-0.25, 3-0.25], df.randWA_sigdiff, align='edge', width=0.2, linewidth=0, color=mint),
	plt.bar([1-0.25, 2-0.25, 3-0.25], df.randWA_WAbetter, align='edge', width=0.2, linewidth=0, color=sail_blue)]

#allAA_list = [plt.bar([1-0.5, 2-0.5, 3-0.5], df.allAA_all, align='edge', width=0.2, linewidth=0, color=grey),
#	plt.bar([1-0.5, 2-0.5, 3-0.5], df.allAA_sigdiff, align='edge', width=0.2, linewidth=0, color=mint),
#	plt.bar([1-0.5, 2-0.5, 3-0.5], df.allAA_WAbetter, align='edge', width=0.2, linewidth=0, color=sail_blue)]
h1 = plt.bar([1-0.5, 2-0.5, 3-0.5], df.allAA_all, align='edge', width=0.2, linewidth=0, color=grey)
h2 = plt.bar([1-0.5, 2-0.5, 3-0.5], df.allAA_sigdiff, align='edge', width=0.2, linewidth=0, color=mint)
h3 = plt.bar([1-0.5, 2-0.5, 3-0.5], df.allAA_WAbetter, align='edge', width=0.2, linewidth=0, color=sail_blue)

plt.text(0.6, 38, 'All-AA\nmodel', fontsize=14, ha='center')
plt.text(0.85, 40, 'Random-WA\nmodel', fontsize=14, ha='center')
plt.text(1.1, 38, 'Combined\nmodel', fontsize=14, ha='center')

plt.text(1.6, 38, 'All-AA\nmodel', fontsize=14, ha='center')
plt.text(1.85, 40, 'Random-WA\nmodel', fontsize=14, ha='center')
plt.text(2.1, 38, 'Combined\nmodel', fontsize=14, ha='center')

plt.text(2.6, 38, 'All-AA\nmodel', fontsize=14, ha='center')
plt.text(2.85, 40, 'Random-WA\nmodel', fontsize=14, ha='center')
plt.text(3.1, 38, 'Combined\nmodel', fontsize=14, ha='center')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])

plt.ylabel('Number of behaviors', fontsize=20)
plt.xlabel('Accuracy metric', fontsize=20)
plt.legend((h3, h2, h1), ('WA better than AA', 'AA better than WA', 'No significant difference'), fontsize=16, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')

fname = '/Users/jli/Documents/Research/my_projects/fairAI/ABCD_race/compare_3models'
plt.savefig(fname+'.png', format='png', transparent=True)
plt.savefig(fname+'.eps', format='eps', transparent=True)