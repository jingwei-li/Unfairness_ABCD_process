import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--metric', required=True, 
	help='The accuracy metric to be plotted. Choose from pCOD, corr and MSE')
args = parser.parse_args()
metric = args.metric

sail_blue = [0, 0.1255, 0.2471, 1]
mint = [0.6784, 0.9373, 0.8196, 1]
grey = [0.9, 0.9, 0.9, 1]
Nbehavior = 36

df = pd.DataFrame(dict(
	allAArandWA_WAbetter = [24, 11],
	allAArandWA_sigdiff = [24+5, 11+20],
	allAArandWA_all = [Nbehavior] * 2,
	randWA_WAbetter = [26, 13],
	randWA_sigdiff = [26+2, 13+16],
	randWA_all = [Nbehavior] * 2,
	allAA_WAbetter = [19, 11],
	allAA_sigdiff = [19+8, 11+22],
	allAA_all = [Nbehavior] * 2))

fig = plt.figure(figsize=(6.5, 8))
# transparent backgroung
ax = plt.gca()
ax.patch.set_alpha(0)
# invisible right and upper axes
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

# tick labels
plt.xticks([1, 2, 3], ['All-AA\nmodel', 'Random-WA\nmodel', 'Combined\nmodel'])
ax.tick_params(axis='both', which='major', labelsize=14, direction='out', width=2)

# metric conditions
if metric == 'pCOD':
	idx = 0
	xstr = 'predictive COD'
elif metric == 'corr':
	idx = 1
	xstr = 'Pearson\'s r'
else:
	print('Unknown metric')

# create bars
bw = 0.7
allAArandWA_list = [plt.bar(3-bw/2, df.allAArandWA_all[idx], align='edge', width=bw, linewidth=0, color=grey),
	plt.bar(3-bw/2, df.allAArandWA_sigdiff[idx], align='edge', width=bw, linewidth=0, color=mint),
	plt.bar(3-bw/2, df.allAArandWA_WAbetter[idx], align='edge', width=bw, linewidth=0, color=sail_blue)]

randWA_list = [plt.bar(2-bw/2, df.randWA_all[idx], align='edge', width=bw, linewidth=0, color=grey),
	plt.bar(2-bw/2, df.randWA_sigdiff[idx], align='edge', width=bw, linewidth=0, color=mint),
	plt.bar(2-bw/2, df.randWA_WAbetter[idx], align='edge', width=bw, linewidth=0, color=sail_blue)]

h1 = plt.bar(1-bw/2, df.allAA_all[idx], align='edge', width=bw, linewidth=0, color=grey)
h2 = plt.bar(1-bw/2, df.allAA_sigdiff[idx], align='edge', width=bw, linewidth=0, color=mint)
h3 = plt.bar(1-bw/2, df.allAA_WAbetter[idx], align='edge', width=bw, linewidth=0, color=sail_blue)

box = ax.get_position()
ax.set_position([box.x0+0.1, box.y0, box.width * 0.7, box.height])

plt.ylabel('Number of behaviors', fontsize=20)
plt.xlabel('Accuracy metric: ' + xstr, fontsize=19)
plt.legend((h3, h2, h1), ('WA better than AA', 'AA better than WA', 'No significant difference'), fontsize=16, frameon=False, bbox_to_anchor=(-0.2, 1.16), loc='upper left')

fname = '/Users/jli/Documents/Research/my_projects/fairAI/ABCD_race/compare_3models_' + metric + '_reg_y_FC'
plt.savefig(fname+'.png', format='png', transparent=True)
plt.savefig(fname+'.eps', format='eps', transparent=True)