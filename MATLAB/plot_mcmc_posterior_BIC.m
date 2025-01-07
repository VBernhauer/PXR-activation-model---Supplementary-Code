function [] = plot_mcmc_posterior_BIC()

    %%% command:
    %%% plot_mcmc_posterior_BIC()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    
    n = 227;
    k = 11;
    k_no_res_dep = 10; 
    
    loglik_load = load('mlvalues.mat');
    loglik = loglik_load.mlvalue_array(:,:);
    BIC = k * log(n) - 2 * loglik;
    
    loglik_no_res_dep = load('../MATLAB_no_resouce_depletion/mlvalues.mat');
    loglik_no_res_dep = loglik_no_res_dep.mlvalue_array(:,:);   
    BIC_no_res_dep = k_no_res_dep * log(n) - 2 * loglik_no_res_dep;

    % length(BIC)
    % length(BIC_no_res_dep)

    figure(1);
    ax = gca;
    set(gcf,'Position', [300, 100, 600, 300]);
    set(ax,...
        'box','on',...
        'FontSize',10);
    set(gca,'TickLength',[0.025, 0.01])
    hold on;
    histogram(BIC,...
          'FaceColor', 'none',...
          'EdgeColor', [0 0 0],...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    histogram(BIC_no_res_dep,...
          'FaceColor', 'none',...
          'EdgeColor', [0 1 1],...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    l = legend('$k_\mathrm{r} > 0$','$k_\mathrm{r} = 0$','Interpreter','latex');
    set(l,'Location','northeast');
    set(l,'FontSize',10);   
    xlabel('Bayesian information criterion','FontSize',14);
    ylabel('Normalized frequency','FontSize',14);


    %%% save figure %%%
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    savefig(strcat('figures/loglikelihood.fig'));
    exportgraphics(gcf,'figures/loglikelihood.png');
    % exportgraphics(gcf,'../../LaTeX/figures/loglikelihood.eps','ContentType','vector');

    end