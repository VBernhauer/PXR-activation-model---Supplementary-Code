function [] = plot_mcmc_posteriors_activity()

    %%% command:
    %%% plot_mcmc_posteriors_activity()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    chains_1OHMid3 = [];
    for jj = 1:5
        jjchains = load(strcat('./chains_1OHMid3/chains_1OHMid3_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = size(jjchains,1)/2:size(jjchains,1)
            if mod(kk,2)==0
                chains_1OHMid3 = [chains_1OHMid3;jjchains(kk,:)];
            end
        end
    end
    
    chains_activity_cyp3a4 = [];
    for jj = 1:5
        jjchains = load(strcat('./chains_activity_cyp3a4/chains_activity_cyp3a4_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = size(jjchains,1)/2:size(jjchains,1)
            if mod(kk,2)==0
                chains_activity_cyp3a4 = [chains_activity_cyp3a4;jjchains(kk,:)];
            end
        end
    end
    
    chains_activity_cyp2c9 = [];
    for jj = 1:5
        jjchains = load(strcat('./chains_activity_cyp2c9/chains_activity_cyp2c9_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = size(jjchains,1)/2:size(jjchains,1)
            if mod(kk,2)==0
                chains_activity_cyp2c9 = [chains_activity_cyp2c9;jjchains(kk,:)];
            end
        end
    end
    
    chains_activity_cyp2b6 = [];
    for jj = 1:5
        jjchains = load(strcat('./chains_activity_cyp2b6/chains_activity_cyp2b6_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = size(jjchains,1)/2:size(jjchains,1)
            if mod(kk,2)==0
                chains_activity_cyp2b6 = [chains_activity_cyp2b6;jjchains(kk,:)];
            end
        end
    end

    a{1} = '$k_{\mathrm{met}_\mathrm{cyp3a4,ss}}$';
    a{2} = '$k_{\mathrm{met}_\mathrm{cyp3a4,sm}}$';
    a{3} = '$k_{\mathrm{met}_\mathrm{cyp2c9,sm}}$';
    a{4} = '$k_{\mathrm{met}_\mathrm{cyp2b6,sm}}$';
    
    %%% Aggregated posterior plots %%%
    fcolor = {[0 0 0],[1 0 0],[0 0 1],[0 1 1],[1 0 1]};
    figlabs = {'(a)','(b)'};
    xlims = {[-5.5 -2.0],[-0.6 1.2]};
    figure(2);
    tiledplot = tiledlayout(1,1,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 600, 400]);
    mm = 1;
    ax(mm) = nexttile(mm);
    set(ax(mm),...
        'box','on',...
        'FontSize',12,...
        'XLim',xlims{mm}); 
    set(gca,'TickLength',[0.025, 0.01])
%     text(0.9,0.9,figlabs{mm},...
%                 'Units','Normalized',...
%                 'HorizontalAlignment','center',...
%                 'FontSize',14,...
%                 'FontWeight','bold');
    hold on;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = histogram(log10(chains_1OHMid3(:,1)),...
          'FaceColor', 'none',...
          'EdgeColor', fcolor{1},...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    h2 = histogram(log10(chains_activity_cyp3a4(:,1)),...
          'FaceColor', 'none',...
          'EdgeColor', fcolor{2},...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    h3 = histogram(log10(chains_activity_cyp2c9(:,1)),...
          'FaceColor', 'none',...
          'EdgeColor', fcolor{3},...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    h4 = histogram(log10(chains_activity_cyp2b6(:,1)),...
          'FaceColor', 'none',...
          'EdgeColor', fcolor{4},...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    l=legend(a{1},a{2},a{3},a{4});
    set(l,'Location','northwest');
    set(l,'FontSize',14);
    set(l,'Interpreter','latex');
    
    xlabel(tiledplot,'Log_{10} value','FontSize',14);
    ylabel(tiledplot,'Normalized frequency','FontSize',14);
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';
    

    savefig(strcat('figures/posteriors_activity.fig'));
    exportgraphics(gcf,'figures/posteriors_activity.png');
%     exportgraphics(gcf,'figures/posteriors_activity.pdf','ContentType','vector');
%     exportgraphics(gcf,'../LaTeX/figures/posteriors_activity.eps','ContentType','vector');
                     

    end