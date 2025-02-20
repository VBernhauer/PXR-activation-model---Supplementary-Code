function [] = plot_mcmc_kinetics_activity(cyp)

    %%% command: plot_mcmc_kinetics_activity(cyp), where cyp = 1 or 2 or 3
   
    close all;
    set(0,'DefaultFigureVisible','on');

    metadata = load('datamat_activity.mat');
    MLpars = readmatrix('maxLikValues.txt');
    mlpars = MLpars(2:end,2);
    
    switch cyp
        case 1
            data    = metadata.data_act{1};
            stdev   = metadata.stdev_act{1};
            k_cyp           = mlpars(4);
            k_cypdeg        = mlpars(5);
            substrate_total = 10*80; 
            label = '3a4';
            ylabelname = '1-OH-midazolam';
        case 2
            data    = metadata.data_act{2};
            stdev   = metadata.stdev_act{2};
            k_cyp           = mlpars(6);
            k_cypdeg        = mlpars(7);
            substrate_total = 10*80; 
            label = '2c9';
            ylabelname = '4-OH-diclofenac';
        case 3
            data    = metadata.data_act{3};
            stdev   = metadata.stdev_act{3};
            k_cyp           = mlpars(8);
            k_cypdeg        = mlpars(9); 
            substrate_total = 40*80; 
            label = '2b6';
            ylabelname = 'OH-bupropion';
    end
    
    k_pxr           = mlpars(1);
    k_r             = mlpars(2);
    k_pxrdeg        = mlpars(3);

    time    = metadata.time_act;
    tmax    = 240;%metadata.tmax_act;
    Rif = 10;
    Sub = substrate_total;
    
    T = 0:0.1:4;
    T_act = 0:0.1:240; %metadata.t_act;
    
    timeXaxis = [0,4];
    
    timeXaxisAct = [0,24,48,72,120,168,240];
    
    YLims = [-5 45];    
    YTick = 0:10:40; 
    
    YLimsPred = [-50 750];    
    YTickPred = 0:100:700; 
    
    Dose = [0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50];
    
    
    titles  = {'24 h rifampicin pre-treatment';...
             '48 h rifampicin pre-treatment';...
             '72 h rifampicin pre-treatment';...
             '120 h rifampicin pre-treatment';...
             '24 h rifampicin pre-treatment';...
             '48 h rifampicin pre-treatment';...
             '72 h rifampicin pre-treatment';...
             '120 h rifampicin pre-treatment'};

    %%%% Each donor in different color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colors = jet(length(Dose));
    for ii = 1:length(Dose)
        Legend{ii} = strcat(num2str(Dose(ii)),' \muM');
    end
    figlabs = {'(A)','(B)','(C)','(D)',...
               '(E)','(F)','(G)','(H)',...
               '(I)','(J)'};
           
    facecol = {'red','green'};
    names = {'Donor 1','Donor 3'};
    markers = {'o','square'};
    markersize = 50;
    markerfa = {0.5,0.5,0.5};
    mcolors = [0 0 0];
    
    figure(1);
    clc;
    tiledplot = tiledlayout(3,4);
    set(gcf, 'Position',  [500, 100, 1250, 800]);
    mm = 0;
    for bb = 1:2
        for aa = 1:size(time,2)
            if bb == 1 || bb == 2
                mm = mm+1;
                ax(mm) = nexttile(mm);
                set(ax(mm),...
                    'box','on',...
                    'XLim',[-0.15*max(timeXaxis) 1.15*max(timeXaxis)],...
                    'XTick',timeXaxis,...
                    'XTickLabel',timeXaxis,...
                    'XTickLabelRotation',45,...
                    'YLim',YLims,...
                    'YTick',YTick,...
                    'FontSize',10);
                set(gca,'TickLength',[0.025, 0.01])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hold on;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if bb == 1
                    for kk = 1:size(data,1)
                        scatter(4,data(kk,aa),markersize,...
                                'Marker',markers{kk},...
                                'MarkerEdgeColor',mcolors,...
                                'MarkerFaceColor',facecol{kk},...
                                'MarkerFaceAlpha',markerfa{kk});
                        hold on;
                    end 
                end
                
                if bb == 2
                    CYP = model();
                    %%% shaded area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    outputvec = [];
                    for jj = 1:5
                        chains = load(strcat('./chains_activity_cyp',label,'/chains_activity_cyp',label,'_',num2str(jj),'.mat'));
                        chains = chains.chains(:,:);
                        for nn = size(chains,1)/2:size(chains,1)
                            if mod(nn,5)==0
                                output = substrate_total*(1 - exp(-chains(nn,:)*T*CYP(aa)));
                                outputvec = [outputvec;output];
                            end
                        end
                    end
                    outputvecmin = min(outputvec,[],1);
                    outputvecmax = max(outputvec,[],1);

                    outputvec_025prc = prctile(outputvec,2.5,1);
                    outputvec_975prc = prctile(outputvec,97.5,1);

                    patch(ax(mm),[T,fliplr(T)],[outputvecmin,fliplr(outputvecmax)],mcolors,'FaceAlpha',0.1,'EdgeColor',mcolors,'EdgeAlpha',0.1);
                    patch(ax(mm),[T,fliplr(T)],[outputvec_025prc,fliplr(outputvec_975prc)],mcolors,'FaceAlpha',0.25,'EdgeColor',mcolors,'EdgeAlpha',0.25);

                    % plot(ax(mm),T,outputvecmin,'Color',[mcolors 0.1],'LineWidth',0.1,'LineStyle','-')
                    % plot(ax(mm),T,outputvecmax,'Color',[mcolors 0.1],'LineWidth',0.1,'LineStyle','-')
                    hold on;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% maximum likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    mlpars_act = readmatrix(strcat('maxLikValues_activity_cyp',label,'.txt'));
                    outputmaxlik = substrate_total*(1 - exp(-mlpars_act(2:end,2)*T*CYP(aa))); 
                    plot(ax(mm),T,outputmaxlik,'Color',mcolors,'LineWidth',1,'LineStyle','-')
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    errorbar(4,nanmean(data(:,aa)),stdev(1,aa),... 
                    'Marker','diamond', ...
                    'MarkerSize',7,...
                    'MarkerEdgeColor',[0 0 0],...
                    'MarkerFaceColor','white',...
                    'LineStyle','none',...
                    'LineWidth',0.5,...
                    'Color',[0 0 0],...
                    'CapSize',7);
                end                
            end
            xlabel('Time (hours)','FontSize',10);
            title(strcat('\rm',titles{mm}),'FontSize',10,'Position',[2.3,37.5,0]);
            text(0.1,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',10,...
                'FontWeight','Normal');
            hold on;
        end
    end
            
    ax(9) = nexttile(9,[1 2]);
    set(ax(9),...
        'box','on',...
        'XLim',[-0.05*max(timeXaxisAct) 1.05*max(timeXaxisAct)],...
        'XTick',timeXaxisAct,...
        'XTickLabel',timeXaxisAct,...
        'YLim',YLimsPred,...
        'YTick',YTickPred,...
        'XTickLabelRotation',45,...
        'FontSize',10);
    set(gca,'TickLength',[0.5*0.025, 0.5*0.01])
    hold on;
    %%% shaded area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outputvec = [];
    for jj = 1:5
        chains = load(strcat('./chains_activity_cyp',label,'/chains_activity_cyp',label,'_',num2str(jj),'.mat'));
        chains = chains.chains(:,:);
        for nn = size(chains,1)/2:size(chains,1)
            if mod(nn,5)==0
                output = model_activity(chains(nn,:));
                outputvec = [outputvec;output];
            end
        end
    end
    outputvecmin = min(outputvec,[],1);
    outputvecmax = max(outputvec,[],1);

    outputvec_025prc = prctile(outputvec,2.5,1);
    outputvec_975prc = prctile(outputvec,97.5,1);

    patch(ax(9),[T_act,fliplr(T_act)],[outputvecmin,fliplr(outputvecmax)],mcolors,'FaceAlpha',0.1,'EdgeColor',mcolors,'EdgeAlpha',0.1);
    patch(ax(9),[T_act,fliplr(T_act)],[outputvec_025prc,fliplr(outputvec_975prc)],mcolors,'FaceAlpha',0.25,'EdgeColor',mcolors,'EdgeAlpha',0.25);

    % plot(ax(9),T_act,outputvecmin,'Color',[mcolors 0.05],'LineWidth',0.1,'LineStyle','-')
    % plot(ax(9),T_act,outputvecmax,'Color',[mcolors 0.05],'LineWidth',0.1,'LineStyle','-')
    hold on;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% maximum likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mlpars_act = readmatrix(strcat('maxLikValues_activity_cyp',label,'.txt'));
    outputmaxlik = model_activity(mlpars_act(2:end,2));
    plot(ax(9),T_act,outputmaxlik,'Color',mcolors,'LineWidth',1,'LineStyle','-')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    xlabel('Time post rifampicin treatment (hours)','FontSize',10);
    text(0.1,0.9,figlabs{9},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',10,...
                'FontWeight','Normal');
    hold on;    

    ax(10) = nexttile(11,[1 2]);
    set(ax(10),...
        'box','on',...
        'XLim',[-0.05*max(timeXaxisAct) 1.05*max(timeXaxisAct)],...
        'XTick',timeXaxisAct,...
        'XTickLabel',timeXaxisAct,...
        'YLim',YLimsPred,...
        'YTick',YTickPred,...
        'XTickLabelRotation',45,...
        'FontSize',10);
    set(gca,'TickLength',[0.5*0.025, 0.5*0.01])
    hold on;
    mlpars_act = readmatrix(strcat('maxLikValues_activity_cyp',label,'.txt'));
    for i = 1:length(Dose)
        output = model_dose(mlpars_act(2:end,2),Dose(i));
        plot(ax(10),T_act,output,'Color',colors(i,:),'LineWidth',1,'LineStyle','-')
        hold on;
    end      
    xlabel('Time post rifampicin treatment (hours)','FontSize',10);
    text(0.1,0.9,figlabs{10},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',10,...
                'FontWeight','Normal');
    hold on;

        
    leg1 = legend(ax(4),names,'Location','eastoutside','FontSize',10,'Orientation','Vertical');
    leg2 = legend(ax(10),Legend{:,:},'Location','eastoutside','FontSize',10,'Orientation','Vertical');
    title(leg2,'RIF concentration, $L_\mathrm{pxr}$','Interpreter','latex')
    ylabel(tiledplot,strcat(ylabelname,' (pmol  per well)'),'FontSize',14);
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    %%% save figure %%%
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    savefig(strcat('figures/activity_cyp',label,'.fig'));
    exportgraphics(gcf,strcat('figures/activity_cyp',label,'.png'));
    % exportgraphics(gcf,strcat('../../LaTeX/figures/activity_cyp',label,'.eps'),'ContentType','vector');

    %%% output of the model %%%
    function [output] = model()         
        solution = ode23s(@ode,[0 tmax],...
                       [0 1],...
                       []);          
        output = deval(solution,time,2);          
    end

    function [output] = model_activity(pars)         
        solution = ode23s(@ode_activity,[0 tmax],...
                       [0 1 Sub 0],...
                       [],...
                       pars);          
        output = deval(solution,T_act,4);          
    end

    function [output] = model_dose(pars,D)         
        solution = ode23s(@ode_dose,[0 tmax],...
                       [0 1 Sub 0],...
                       [],...
                       pars,...
                       D);          
        output = deval(solution,T_act,4);          
    end

    %%% ODE system dose %%%
    function [dxdt] = ode_dose(t,x,pars,D)        
        dxdt = zeros(4,1);
        
        dxdt(1) = k_pxr*(1-x(1))*D*exp(-k_r*t) - k_pxrdeg*x(1);          % activated PXR
        dxdt(2) = k_cyp*x(1) + k_cypdeg*(1-x(2));                          % CYP
        dxdt(3) = - pars*x(2)*x(3);
        dxdt(4) = pars*x(2)*x(3);                                              % activity
    end

    %%% ODE system activity %%%
    function [dxdt] = ode_activity(t,x,pars)        
        dxdt = zeros(4,1);
        
        dxdt(1) = k_pxr*(1-x(1))*Rif*exp(-k_r*t) - k_pxrdeg*x(1);          % activated PXR
        dxdt(2) = k_cyp*x(1) + k_cypdeg*(1-x(2));                          % CYP
        dxdt(3) = - pars*x(2)*x(3);
        dxdt(4) = pars*x(2)*x(3);                                              % activity
    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x)        
        dxdt = zeros(2,1);

        dxdt(1) = k_pxr*(1-x(1))*Rif*exp(-k_r*t) - k_pxrdeg*x(1);          % activated PXR
        dxdt(2) = k_cyp*x(1) + k_cypdeg*(1-x(2));                          % CYP
    end
            
end