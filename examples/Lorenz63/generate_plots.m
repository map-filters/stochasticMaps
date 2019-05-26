% declare labels for lines in plot
algorithms = {'EnKF','Linear','Linear + 1 RBF','Linear + 2 RBF'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT DATA                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filter_all = cell(length(M_vect),4);
filter_out = cell(length(M_vect),4);

for j = 1:length(M_vect)
	M = M_vect(j);

	% load enkf data
	load(['data/enkf_pert_M' num2str(M)],'filter');
	[~, opt_idx] = min([filter.RMSE_mean]);
	filter_all{j,1} = filter;
	filter_out{j,1} = filter(opt_idx);

	% load order 1 data
	load(['data/SM_order1_M' num2str(M)],'filter');
	[~, opt_idx] = min([filter.RMSE_mean]);
 	filter_all{j,2} = filter;
	filter_out{j,2} = filter(opt_idx);

	% load order 2 data
	load(['data/SM_order2_M' num2str(M)],'filter');
	[~, opt_idx] = min([filter.RMSE_mean]);
 	filter_all{j,3} = filter;
	filter_out{j,3} = filter(opt_idx);

	% load order 3 data
	load(['data/SM_order3_M' num2str(M)],'filter');
	[~, opt_idx] = min([filter.RMSE_mean]);
	filter_all{j,4} = filter;
	filter_out{j,4} = filter(opt_idx);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CONVERGENCE OF STATISTICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
hold on
for k = 1:size(filter_out,2)
    filt_array = cell2mat(filter_out(:,k));
	plot(M_vect, [filt_array.RMSE_mean], '-o', 'markerFaceColor', 'auto', ...
		'MarkerSize', 7, 'DisplayName', algorithms{k},'LineWidth',2.5)
end
plot(M_vect, ones(length(M_vect),1)*model.sigma_y, '--k', 'DisplayName', '$\sigma_{y}$')
set(gca,'FontSize',17)
xlabel('Ensemble size $M$','FontSize',20)
ylabel('Average RMSE','FontSize',20)
set(gca,'XTick',[10,60,100,200,400,600]);
ylim([0.25,0.9])
xlim([min(M_vect), max(M_vect)])
h = legend('show');
set(h,'FontSize',16)
set(gca,'LineWidth',2)
hold off
print('-depsc','figures/L63_mean_RMSE')

figure()
hold on
for k=1:size(filter_out,2)
    filt_array = cell2mat(filter_out(:,k));
    plot(M_vect, [filt_array.spread_mean], '-o', 'markerFaceColor', 'auto', ...
    	'MarkerSize', 7, 'DisplayName', algorithms{k}, 'LineWidth', 2.5)
end
set(gca,'FontSize',17)
xlabel('Ensemble size $M$', 'FontSize', 20)
ylabel('Average spread', 'FontSize', 20)
set(gca,'XTick',[10,60,100,200,400,600]);
xlim([min(M_vect), max(M_vect)])
h = legend('show');
set(h,'FontSize',16);
set(gca,'LineWidth',2)
hold off
print('-depsc','figures/L63_mean_spread')

figure()
hold on
for k=1:size(filter_out,2)
    filt_array = cell2mat(filter_out(:,k));
    plot(M_vect, [filt_array.covprob_mean], '-o', 'markerFaceColor', 'auto', ...
    	'MarkerSize', 7, 'DisplayName', algorithms{k}, 'LineWidth', 2.5)
end
set(gca,'FontSize',17)
xlabel('Ensemble size $M$', 'FontSize', 20)
ylabel('Average coverage probability', 'FontSize', 20)
set(gca,'XTick',[10,60,100,200,400,600]);
ylim([0.8,1.0])
xlim([min(M_vect), max(M_vect)])
h = legend('show','location','southeast');
set(h,'FontSize',16)
set(gca,'LineWidth',2)
hold off
print('-depsc','figures/L63_mean_covprob')

% -- END OF FILE --
