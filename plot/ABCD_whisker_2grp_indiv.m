function ABCD_whisker_2grp_indiv(data, colormat, y_label, legends, ...
	tit, colloq_nm, sigdiff_idx, outdir, outstem, metric)

	% ABCD_whisker_2grp_indiv(data, colormat, y_label, legends, ...
	%	  tit, colloq_nm, sigdiff_idx, outdir, outstem)
	%
	% Inputs:
	%   - data: #splits x 3 x #behaviors matrix. Accuracy of group 1, 
	%                group 2, and the difference.
	%   - colormat: 3 x 3 color mat, each row is the RGB of one group.
	%   - y_label: Y-axis label.
	%   - legends: cell with 3 entries. Legend for each group and the difference.
	%   - colloq_nm: colloquial behavioral names.
	%   - sigdiff_idx: indices of behaviors with significant group difference.
	%   - outdir: output directory, full path.
	%   - outstem: output name relative to outdir, without extension.
	%

	addpath(genpath( '/home/jingweil/storage/from_HOME/code/plotting_functions/'))
	nbhvr = size(data, 3);

	if(~exist('metric', 'var'))
		metric = 'predictive_COD';
	end

	f = figure('visible', 'off');
	aboxplot(data, 'colormap', colormat)
	hold on
	xlimit = get(gca, 'xlim');
	plot(xlimit, [0 0], ':k');
	hold off
	
	pf = get(gcf, 'position');
	set(gcf, 'position', [0 0 100+90*nbhvr 900])
	set(gca, 'position', [0.35 0.4 0.6 0.5])
	
	if(strfind(metric, 'COD'))
		ylm = get(gca, 'ylim');
		if(ylm(1)<-1)
			ylm(1) = -1;
			warning('There are behaviors with accuracy lower than -1.')
		end
		if(ylm(2)>2)
			ylm(2) = 1;
			warning('There are values higher than 1.')
		end
		set(gca, 'ylim', ylm, 'ytick', [ylm(1):0.2:ylm(2)]);
	end
	yl = ylabel(y_label);
	set(yl, 'fontsize', 16, 'linewidth', 2)
	
	l = legend(legends);
	set(l, 'fontsize', 12, 'linewidth', 2, 'location', 'best', 'box', 'off')
	
	set(gca, 'xticklabel', colloq_nm, 'fontsize', 16, 'linewidth', 2);
	rotateXLabels( gca(), 45 );
	set(gca, 'tickdir', 'out', 'box', 'off')
	
	% plot * on the significant behaviors
	if(~isempty(sigdiff_idx))
		ylimvals = get(gca, 'ylim');
		ylimvals_new = ylimvals;
		ylimvals_new(2) = (ylimvals(2)-ylimvals(1))*0.08 + ylimvals(2);
		set(gca, 'ylim', ylimvals_new)
		lp = get(l, 'position');
		if(lp(2)>0.5)
			lp(2) = lp(2) - 0.03;
		end
		set(l, 'position', lp)
		text(sigdiff_idx, repmat(ylimvals(2), size(sigdiff_idx)), '*');
	end
	
	if(~isempty(tit))
		title(tit, 'fontsize', 16)
	end

	outname = fullfile(outdir, [outstem ]);
	export_fig(outname, '-png', '-nofontswap', '-a1');
	set(gcf, 'color', 'w');
	hgexport(f, outname)
	close

	rmpath(genpath( '/home/jingweil/storage/from_HOME/code/plotting_functions/'))
end