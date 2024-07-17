function [] = plot_mcmc_posteriors()

    %%% command:
    %%% plot_mcmc_posteriors()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    chains = [];
    for ii = 1:5
        for jj = 1:5
            jjchains = load(strcat('./chains_',num2str(ii),'/chains_',num2str(jj),'.mat'));
            jjchains = jjchains.chains(:,:);
            for kk = size(jjchains,1)/2:size(jjchains,1)
                if mod(kk,50)==0
                    chains = [chains;jjchains(kk,:)];
                end
            end
        end
    end
    
    a{1} = 'k_\mathrm{pxr,max}';
    a{2} = 'k_\mathrm{r}';
    a{3} = 'k_\mathrm{pxr,deg}';
    a{4} = 'k_{mRNA_\mathrm{cyp3a4}^\mathrm{fold}}';
    a{5} = 'k_{mRNA_\mathrm{cyp3a4,deg}}';
    a{6} = 'k_{mRNA_\mathrm{cyp2c9}^\mathrm{fold}}';
    a{7} = 'k_{mRNA_\mathrm{cyp2c9,deg}}';
    a{8} = 'k_{mRNA_\mathrm{cyp2b6}^\mathrm{fold}}';
    a{9} = 'k_{mRNA_\mathrm{cyp2b6,deg}}';
    a{10} = 'k_{mRNA_\mathrm{mdr1}^\mathrm{fold}}';
    a{11} = 'k_{mRNA_\mathrm{mdr1,deg}}';
    labs = {a{1},...
            a{2},...
            a{3},...
            a{4},...
            a{5},...
            a{6},...
            a{7},...
            a{8},...
            a{9},...
            a{10},...
            a{11}};

%     ax=figure(1);
%     ecornerplot(log10(chains),'ks',true,'names',labs);
        
    ax=figure(1);
    clf;
    nn = 0;
    for aa = 1:length(labs)
        for bb = 1:length(labs)
            nn = nn+1;
            if aa<bb
                ax(nn)=subaxis(length(labs),length(labs),nn,'sh',0.01,'sv',0.01);    
                set(ax(nn),'Visible','off');
            end
            if aa == bb
                %%% subaxis source %%%
                % Aslak Grinsted (2021). Subaxis - Subplot (https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot), 
                % MATLAB Central File Exchange. Retrieved July 26, 2021.
                ax(nn)=subaxis(length(labs),length(labs),nn,'sh',0.01,'sv',0.01);
                hold on
                histogram(log10(chains(:,bb)),...
                          'FaceColor', [0 0 0],...
                          'EdgeColor', [0 0 0],...
                          'FaceAlpha', 0.25,...
                          'EdgeAlpha', 1,...
                          'Normalization', 'probability')
                if aa == length(labs)
                    xlabel(strcat('$',labs{bb},'$'),'Interpreter','latex')
                    a = gca;
                    a.FontSize = 8; 
                    set(get(a,'XLabel'),'Rotation',0);
                    xtickangle(45);
                    ax(nn).XAxis.TickLabelFormat = '%.2f';
                end
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                if aa < length(labs)
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                end            
            end
            if aa>bb
                %%% subaxis source %%%
                % Aslak Grinsted (2021). Subaxis - Subplot (https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot), 
                % MATLAB Central File Exchange. Retrieved July 26, 2021.
                ax(nn)=subaxis(length(labs),length(labs),nn,'sh',0.01,'sv',0.01);
                hold on
                scatter(log10(chains(:,bb)), log10(chains(:,aa)),3,'filled','black',...
                        'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05)
                if aa == length(labs)
                    xlabel(strcat('$',labs{bb},'$'),'Interpreter','latex')
                    a = gca;
                    a.FontSize = 8; 
                    xtickangle(45);
                    ax(nn).XAxis.TickLabelFormat = '%.2f';
                end
                if bb == 1
                    ylabel(strcat('$',labs{aa},'$'),'Interpreter','latex')
                    a = gca;
                    a.FontSize = 8; 
                    set(get(gca,'ylabel'),'rotation',90);
                    hs = get(gca, 'YLabel');
                    pos = get(hs, 'Position');
                    if aa == 7
                        pos(1) = pos(1)+0.4;
                    else
                        pos(1) = pos(1)+0.3;
                    end
                    set(hs, 'Position', pos);
                    ax(nn).YAxis.TickLabelFormat = '%.2f';
                    ytickangle(45)
                end
                if aa < length(labs)
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                end
                if bb > 1
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                end
            end
        end
    end
    set(gcf,'PaperOrientation','landscape');
    
    savefig(strcat('figures/corner.fig'));
%     exportgraphics(gcf,'figures/corner.png');
%     exportgraphics(gcf,'../LaTeX/figures/corner.eps','ContentType','vector'); 


    loga{1} = '$\log_{10}k_\mathrm{pxr,max}$';
    loga{2} = '$\log_{10}k_\mathrm{r}$';
    loga{3} = '$\log_{10}k_\mathrm{pxr,deg}$';
    loga{4} = '$\log_{10}k_{mRNA_\mathrm{cyp3a4}^\mathrm{fold}}$';
    loga{5} = '$\log_{10}k_{mRNA_\mathrm{cyp3a4,deg}}$';
    loga{6} = '$\log_{10}k_{mRNA_\mathrm{cyp2c9}^\mathrm{fold}}$';
    loga{7} = '$\log_{10}k_{mRNA_\mathrm{cyp2c9,deg}}$';
    loga{8} = '$\log_{10}k_{mRNA_\mathrm{cyp2b6}^\mathrm{fold}}$';
    loga{9} = '$\log_{10}k_{mRNA_\mathrm{cyp2b6,deg}}$';
    loga{10} = '$\log_{10}k_{mRNA_\mathrm{mdr1}^\mathrm{fold}}$';
    loga{11} = '$\log_{10}k_{mRNA_\mathrm{mdr1,deg}}$';
    
    %%% Aggregated posterior plots %%%
    fcolor = {{[0 0 0],[1 1 0],[0 1 0]};...
              {[1 0 0],[0 0 1],[0 1 1],[1 0 1]};...
              {[1 0 0],[0 0 1],[0 1 1],[1 0 1]}};
    figlabs = {'(a)','(b)','(c)'};
    xlims = {[-3.0 -0.4],[-2.3 -0.5],[-2.1 -1.2]};
    figure(2);
    tiledplot = tiledlayout(1,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1700, 350]);
    mm = 0;
    for aa = 1:3
        mm = mm+1;
        ax(mm) = nexttile(mm);
        set(ax(mm),...
            'box','on',...
            'FontSize',12,...
            'XLim',xlims{mm});
        set(gca,'TickLength',[0.025, 0.01])
        text(0.9,0.9,figlabs{mm},...
                    'Units','Normalized',...
                    'HorizontalAlignment','center',...
                    'FontSize',14,...
                    'FontWeight','bold');
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch mm
            case 1
                i = 0;
                for nn = 1:3
                    i = i + 1;
                    histogram(log10(chains(:,nn)),...
                          'FaceColor', 'none',...
                          'EdgeColor', fcolor{mm}{i},...
                          'FaceAlpha', 1,...
                          'EdgeAlpha', 1,...
                          'Normalization', 'probability',...
                          'DisplayStyle', 'stairs',...
                          'LineWidth',1.5);
                hold on;
                end
                l=legend(loga{1},loga{2},loga{3});
%                     l=legend({strcat('$\Log_{10}$',labs{1}),strcat('Log_{10}',labs{2}),strcat('Log_{10}',labs{3})});
                set(l,'Location','northwest');
                set(l,'FontSize',12);
                set(l,'interpreter','latex');
%                 set(l,'NumColumns',2);
            case 2
                i = 0;
                for nn = 2:5
                    i = i + 1;
                    histogram(log10(chains(:,2*nn)),...
                          'FaceColor', 'none',...
                          'EdgeColor', fcolor{mm}{i},...
                          'FaceAlpha', 1,...
                          'EdgeAlpha', 1,...
                          'Normalization', 'probability',...
                          'DisplayStyle', 'stairs',...
                          'LineWidth',1.5);
                hold on;
                end
                l=legend(loga{4},loga{6},loga{8},loga{10});
%                     l=legend(strcat('Log_{10}',labs{4}),strcat('Log_{10}',labs{6}),strcat('Log_{10}',labs{8}),strcat('Log_{10}',labs{10}));
                set(l,'Location','northwest');
                set(l,'FontSize',12);
                set(l,'interpreter','latex');
%                 set(l,'NumColumns',2);
            case 3
                i = 0;
                for nn = 2:5
                    i = i + 1;
                    histogram(log10(chains(:,2*nn+1)),...
                          'FaceColor', 'none',...
                          'EdgeColor', fcolor{mm}{i},...
                          'FaceAlpha', 1,...
                          'EdgeAlpha', 1,...
                          'Normalization', 'probability',...
                          'DisplayStyle', 'stairs',...
                          'LineWidth',1.5);
                hold on;
                end
                l=legend(loga{5},loga{7},loga{9},loga{11});
%                     l=legend(strcat('Log_{10}',labs{5}),strcat('Log_{10}',labs{7}),strcat('Log_{10}',labs{9}),strcat('Log_{10}',labs{11}));
                set(l,'Location','northwest');
                set(l,'FontSize',12);
                set(l,'interpreter','latex');
%                 set(l,'NumColumns',2);
        end
    end
    xlabel(tiledplot,'Log_{10} value','FontSize',14);
    ylabel(tiledplot,'Normalized frequency','FontSize',14);
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';


    savefig(strcat('figures/posteriors.fig'));
    exportgraphics(gcf,'figures/posteriors.png');
%     exportgraphics(gcf,'figures/posteriors.pdf','ContentType','vector');
%     exportgraphics(gcf,'../LaTeX/figures/posteriors.eps','ContentType','vector');
                     

    end