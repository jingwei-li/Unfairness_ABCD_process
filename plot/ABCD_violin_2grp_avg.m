function ABCD_violin_2grp_avg(avg_data, colormat, y_label, x_labels, outdir, outstem)

	% ABCD_violin_2grp_avg(avg_data, colormat, y_label, x_labels, outdir, outstem)
	%
	% Inputs:
	%   - avg_data: 3 x #splits matrix. Average accuracy of group 1, group 2, 
	%               and difference, respectively.
	%   - colormat: 3 x 3 color mat, each row is the RGB of one group.
	%   - y_label: Y-axis label.
	%   - x_labels: cell with 3 entries. Name for each group and the difference.
	%   - outdir: output directory, full path.
	%   - outstem: output filename stem. The full output filename will be
	%              <outdir>/Mean_<outstem>.png(.eps)


	f = figure('visible', 'off');
	hold on
    vio = violinplot(avg_data, [], [], 'ShowMean', true);
    for i = 1:length(vio)
        vio(i).ViolinPlot.LineWidth = 2;
        vio(i).ScatterPlot.Marker = '.';
        vio(i).MedianPlot.SizeData = 15;
        vio(i).ViolinPlot.FaceColor = colormat(i,:);
    end
	xlimit = get(gca, 'xlim');
	plot(xlimit, [0 0], ':k');
	hold off
	
	pf = get(gcf, 'position');
	set(gcf, 'position', [0 0 300 800]);
	set(gca, 'position', [0.3 0.3 0.6 0.6])
	ylm = get(gca, 'ylim');
	%set(gca, 'ytick', ylm(1):0.2:ylm(2))
	yl = ylabel(y_label, 'fontsize', 16, 'linewidth', 2);
	
	set(gca, 'XTickLabel', x_labels, 'fontsize', 16, 'linewidth', 2)
	rotateXLabels( gca(), 45 );
	set(gca, 'tickdir', 'out', 'box', 'off');
	
	outname = fullfile(outdir, ['Mean_' outstem]);
	export_fig(outname, '-png', '-nofontswap', '-a1');
	set(gcf, 'color', 'w');
	hgexport(f, outname)
	close

end