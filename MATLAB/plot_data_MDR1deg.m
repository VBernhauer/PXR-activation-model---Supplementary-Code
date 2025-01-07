function [] = plot_data_MDR1deg()

    %%% command:
    %%% plot_data_MDR1deg()
   
    clc;
    close all;
    clear all;
    set(0,'DefaultFigureVisible','on');
        
    % load data and parameters
    metadata = load('datamat.mat');
    
    %%% CYPs %%%
    data = metadata.dataMDR1deg;
    time = metadata.time_MDR1deg;
    timeXaxis = [0,8,12,24,48,72,120,168];
    YLims = [0.6 1.4];    
    YTick = 0.7:0.2:1.3;    
    YTickLabel = {'0.7','0.9','1.1','1.3'};
    ylabels = ' mRNA degradation';
    Text  = 'MDR1';
    

    %%%% Each donor in different color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    facecol = {'red','blue','green'};
    names = {'Donor 1','Donor 2','Donor 3'};
    markers = {'o','diamond','square'};
    markersize = 50;
    markerfa = {0.5,0.5,0.5};
    mcolors = [0 0 0];
    
    tiledplot = tiledlayout(1,1,'TileSpacing','Compact');
%     set(gcf, 'Position',  [500, 100, 500, 350]);
    ax(1) = nexttile(1);
    set(ax(1),...
        'box','on',...
        'XLim',[-0.05*max(timeXaxis) 1.05*max(timeXaxis)],...
        'XTick',timeXaxis,...
        'XTickLabel',timeXaxis,...
        'XTickLabelRotation',45,...
        'YLim',YLims,...
        'YTick',YTick,...
        'YTickLabel',YTickLabel,...
        'FontSize',10);
        set(gca,'TickLength',[0.025, 0.01])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kk = 1:size(data,1)
            scatter(time,data(kk,:),markersize,... 
                'Marker',markers{kk}, ...
                'MarkerEdgeColor',mcolors,...
                'MarkerFaceColor',facecol{kk},...
                'MarkerFaceAlpha',markerfa{kk});
        end
        hold on;
        plot(0:1:168,ones(1,length(0:1:168)),'--','Color','black');
%         yline(1,'--');
    title(strcat('\rm',Text,ylabels),'FontSize',10);                  
    leg = legend(ax(1),names,'Location','SouthOutside','FontSize',10,'Orientation','Horizontal');
    leg.Layout.Tile = 'North';
    xlabel(tiledplot,'Time (hours)','FontSize',14);
%     ylabel(tiledplot,'Fold change to corresponding control','FontSize',14);
    ylabel(tiledplot,['Fold change',newline,'to corresponding control'],'FontSize',14);  
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';


    %%% save figure %%%
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    savefig(strcat('figures/data_MDR1deg.fig'));
    exportgraphics(gcf,'figures/data_MDR1deg.png');
    % exportgraphics(gcf,'../../LaTeX/figures/data_MDR1deg.eps','ContentType','vector');
            
    end