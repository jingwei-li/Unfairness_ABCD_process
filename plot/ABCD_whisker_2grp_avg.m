function ABCD_whisker_2grp_avg(avg_data, colormat, y_label, x_labels, outdir, outstem)

	% ABCD_whisker_2grp_avg(avg_data, colormat, y_label, x_labels, outdir, outstem)
	%
	% Inputs:
	%   - avg_data: 3 x #behaviors matrix. Average accuracy of group 1, group 2, 
	%               and difference, respectively.
	%   - colormat: 3 x 3 color mat, each row is the RGB of one group.
	%   - y_label: Y-axis label.
	%   - x_labels: cell with 3 entries. Name for each group and the difference.
	%   - outdir: output directory, full path.
	%   - outstem: output filename stem. The full output filename will be
	%              <outdir>/Mean_<outstem>.png(.eps)

	addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

	f = figure('visible', 'off');
	aboxplot(avg_data, 'colormap', colormat, 'colorrev', 1);
	hold on
	xlimit = get(gca, 'xlim');
	plot(xlimit, [0 0], ':k');
	hold off
	
	pf = get(gcf, 'position');
	set(gcf, 'position', [0 0 300 800]);
	set(gca, 'position', [0.3 0.3 0.6 0.6])
	ylm = get(gca, 'ylim');
	set(gca, 'ytick', ylm(1):0.2:ylm(2))
	yl = ylabel(y_label, 'fontsize', 16, 'linewidth', 2);
	
	set(gca, 'XTickLabel', x_labels, 'fontsize', 16, 'linewidth', 2)
	rotateXLabels( gca(), 45 );
	set(gca, 'tickdir', 'out', 'box', 'off');
	
	outname = fullfile(outdir, ['Mean_' outstem]);
	export_fig(outname, '-png', '-nofontswap', '-a1');
	set(gcf, 'color', 'w');
	hgexport(f, outname)
	close

	rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))
end