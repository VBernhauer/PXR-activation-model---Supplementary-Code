function [] = plot_CYP3A4_validation()

    %%% command:
    %%% plot_CYP3A4_validation()  
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
        
    metadata = load('datamat.mat');
    
    pars = readtable('maxLikValues.txt');
    pars = pars.Var2(2:end);
       
    Dose    = [0.1,1,10,30];
    
    tmax_valid    = metadata.tmax_h;
    T_valid       = metadata.t_h;


    %%% CYP3A4 original %%%
    Rif = {1,10};
    data = metadata.data_l;
    tmax = metadata.tmax_l;      
    T = metadata.t_l;   
    ylabels = {'CYP3A4 fold mRNA expression, 1\muM';...
               'CYP3A4 fold mRNA expression, 10\muM'};
    mcolors = {[0 0 1];...
               [1 0 0]};

    data = vertcat(data{:});


    %%% CYP3A4 validation %%%
    time_valid          = metadata.time_valid;
    dataCYP3A4_valid    = metadata.dataCYP3A4_valid;
    
    timeXaxis   = [0,24,48,72,120,168];
    YLims       = [0.4 4.6];    
    YTick       = 1:1:4;     
    yaxislabels = "";
              
    
    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    colors = [[0 1 1];[0 0 1];[1 0 0];[1 0 1]];
    figlabs = {'(a)','(b)','(c)'};

    for ii = 1:length(Dose)
        Legend{ii} = strcat(num2str(Dose(ii)),' \muM');
    end

    markersize = 40;
    ftsize = 16;

    tiledplot = tiledlayout(1,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1500, 400]);
    for mm = 1:3
        ax(mm) = nexttile(mm);
        set(ax(mm),...
            'box','on',...
            'XLim',[-0.05*max(timeXaxis) 1.05*max(timeXaxis)],...
            'XTick',timeXaxis,...
            'XTickLabel',timeXaxis,...
            'XTickLabelRotation',45,...
            'YLim',YLims,...
            'YTick',YTick,...
            'FontSize',ftsize);
        set(gca,'TickLength',[0.015, 0.01])

        xlabel('Time (hours)','FontSize',14);
        ylabel(textwrap(yaxislabels,35),'FontSize',ftsize);    

        text(0.1,0.9,figlabs{mm},...
        'Units','Normalized',...
        'HorizontalAlignment','center',...
        'FontSize',ftsize,...
        'FontWeight','bold');

        hold on;
        if mm == 1
            for i = 1:length(Dose)
                output = model(pars,Dose(i),tmax_valid,T_valid);
                if i == length(Dose) 
                    output_24h_30uM = model(pars,Dose(4),tmax_valid,24);
                    disp(['CYP3A4 mRNA fold predicted: ' num2str(output_24h_30uM(1))])
                    disp(['CYP3A4 mRNA fold measured: ' num2str(dataCYP3A4_valid(4,1))])
                end
                plot(ax(1),T_valid,output,'Color',colors(i,:),'LineWidth',1,'LineStyle','-')
                hold on;                
            end           
            for kk = 1:size(dataCYP3A4_valid,1)
                scatter(time_valid,dataCYP3A4_valid(kk,:),markersize,... 
                    'Marker','v', ...
                    'MarkerEdgeColor',[0 0 0],...
                    'MarkerFaceColor',colors(kk,:),...
                    'MarkerFaceAlpha',0.5);
            end
            title(strcat('\rm','CYP3A4 fold mRNA expression'),'FontSize',ftsize);
        else
            rr = mm - 1;
            %%% shaded area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            outputvec = [];
            for jj = 1:5
                for oo = 1:5
                    chains = load(strcat('./chains_',num2str(jj),'/chains_',num2str(oo),'.mat'));
                    chains = chains.chains(:,:);
                    for nn = size(chains,1)/2:size(chains,1)
                        if mod(nn,100)==0
                            output = model(chains(nn,:),Rif{rr},tmax,T);
                            outputvec = [outputvec;output];
                        end
                    end
                end
            end
            outputvecmin = min(outputvec,[],1);
            outputvecmax = max(outputvec,[],1);

            outputvec_025prc = prctile(outputvec,2.5,1);
            outputvec_975prc = prctile(outputvec,97.5,1);

            patch(ax(mm),[T,fliplr(T)],[outputvecmin,fliplr(outputvecmax)],mcolors{rr},'FaceAlpha',0.1,'EdgeColor',mcolors{rr},'EdgeAlpha',0.1);
            hold on;
            patch(ax(mm),[T,fliplr(T)],[outputvec_025prc,fliplr(outputvec_975prc)],mcolors{rr},'FaceAlpha',0.25,'EdgeColor',mcolors{rr},'EdgeAlpha',0.25);
            hold on;
            % plot(ax(mm),T,outputvecmin,'Color',[mcolors{rr} 0.1],'LineWidth',0.1,'LineStyle','-')
            % hold on;
            % plot(ax(mm),T,outputvecmax,'Color',[mcolors{rr} 0.1],'LineWidth',0.1,'LineStyle','-')
            % hold on;
            % plot(ax(mm),T,outputvec_025prc,'Color',[mcolors{rr} 0.2],'LineWidth',0.1,'LineStyle','-')
            % hold on;
            % plot(ax(mm),T,outputvec_975prc,'Color',[mcolors{rr} 0.2],'LineWidth',0.1,'LineStyle','-')
            % hold on;
            %%% maximum likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            MLpars = readmatrix('maxLikValues.txt');
            outputmaxlik = model(MLpars(2:end,2),Rif{rr},tmax,T);
            plot(ax(mm),T,outputmaxlik,'Color',mcolors{rr},'LineWidth',1,'LineStyle','-')

            for kk = 1:size(dataCYP3A4_valid,1)
                if kk == mm
                    scatter(time_valid,dataCYP3A4_valid(kk,:),markersize,... 
                        'Marker','v', ...
                        'MarkerEdgeColor',[0 0 0],...
                        'MarkerFaceColor',colors(kk,:),...
                        'MarkerFaceAlpha',0.5);
                end
            end
            title(strcat('\rm',ylabels{rr}),'FontSize',ftsize);
        end
    end

    leg = legend(ax(1),Legend{:,:},'Location','northeast','FontSize',12,'Orientation','Vertical');
    title(leg,'RIF concentration, $L_\mathrm{pxr}$','Interpreter','latex')
    % xlabel(tiledplot,'Time (hours)','FontSize',14);    
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    %%% save figure %%%
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    savefig(strcat('figures/CYP3A4_validation.fig'));
    exportgraphics(gcf,'figures/CYP3A4_validation.png');
    % exportgraphics(gcf,'../../LaTeX/figures/CYP3A4_validation.eps','ContentType','vector');
             
    
    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% output of the model %%%
    function [output] = model(pars,Dose,tmax,T) 
        solution = ode23s(@ode,[0 tmax],...
                               [0 1],...
                               [],...
                               pars,...
                               Dose);
        output = deval(solution,T,2);

    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x,pars,Xint)        
        dxdt = zeros(2,1);

        dxdt(1) = pars(1)*(1-x(1))*Xint*exp(-pars(2)*t) - pars(3)*x(1);   % activated PXR
        dxdt(2) = pars(4)*x(1) + pars(5)*(1-x(2));                          % CYP3A4
    end

end