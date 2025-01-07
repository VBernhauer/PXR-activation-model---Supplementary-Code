function [] = posterior_values_ci()

    %%% command:
    %%% posterior_values_ci()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    metadata = load('datamat.mat');

    data    = {metadata.data_deg;metadata.data_l;metadata.data_h};
    stdev   = {metadata.stdev_deg;metadata.stdev_l;metadata.stdev_h};
    time    = {metadata.time_deg;metadata.time_l;metadata.time_h};
    tmax    = {metadata.tmax_l;metadata.tmax_h};
    Rif     = {1;10};

    data = vertcat(data{:});
    stdev = vertcat(stdev{:});
 
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    chains = [];
    for kk = 1:5
        for jj = 1:5
            jjchains = load(strcat('./chains_',num2str(kk),'/chains_',num2str(jj),'.mat'));
            jjchains = jjchains.chains(:,:);
            burnin = size(jjchains,1)/2;
            chains = [chains; jjchains(burnin+1:end,:)];
        end
    end
    
    mlvalue = -99999.99;
    mlpars = zeros(1,size(chains,2));
    mlvalue_array = [];
    for ll = 1:size(chains,1)
        if mod(ll,100)==0
            disp(ll);
        end
        llvalue = loglike(chains(ll,:));
        mlvalue_array = [mlvalue_array; llvalue];
        if llvalue > mlvalue
            mlvalue = llvalue;
            mlpars = chains(ll,:);
        end
    end
    %%% save chains %%%
    save(strcat('./mlvalues.mat'),'mlvalue_array');
            
    parsmean = mean(chains);
    parsmedian = median(chains);
    parsquantile = [];
    for kk = 1:size(chains,2)
        parsquantile = [parsquantile; quantile(chains(:,kk),[0.05 0.95])];
    end
    
    
    a{1} = 'k_{pxr}';
    a{2} = 'k_{res}';
    a{3} = 'k_{pxr,deg}';
    a{4} = 'k_{mRNA_{cyp3a4}^{fold}}';
    a{5} = 'k_{mRNA_{cyp3a4,deg}}';
    a{6} = 'k_{mRNA_{cyp2c9}^{fold}}';
    a{7} = 'k_{mRNA_{cyp2c9,deg}}';
    a{8} = 'k_{mRNA_{cyp2b6}^{fold}}';
    a{9} = 'k_{mRNA_{cyp2b6,deg}}';
    a{10} = 'k_{mRNA_{mdr1}^{fold}}';
    a{11} = 'k_{mRNA_{mdr1,deg}}';
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
        
    names = {'parameter','mean','median','95CI lower', '95CI upper'};
    
    fid = fopen('posteriorValues.txt','w');

    fprintf(fid, '%2s %2s %2s %2s %2s\n', names{:});
    for kk = 1:length(labs)
        fprintf(fid,'%0s %.8f %.8f %.8f %.8f\n',labs{kk},parsmean(kk),parsmedian(kk),parsquantile(kk,1),parsquantile(kk,2));
    end
    fclose(fid);
    
    fidml = fopen('maxLikValues.txt','w');

    fprintf(fidml, '%2s %.8f \n', 'MLvalue', mlvalue);
    for kk = 1:length(labs)
        fprintf(fidml,'%0s %.8f\n',labs{kk},mlpars(kk));
    end
    fclose(fidml);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars)        
        [solution] = model(pars);
        value = 0;
        for ii=1:length(data)
            for mm=1:size(data{ii},1)
                idx=find(~isnan(data{ii}(mm,:)));
                value = value+sum(lognormpdf(data{ii}(mm,idx),solution{ii}(1,idx),stdev{ii}(1,idx)));
            end        
        end
    end

    %%% output of the model %%%
    function [output] = model(pars)            
        %%% degradation of protein mRNA %%%
        output{1} = exp(-pars(5) .* time{1});
        output{2} = exp(-pars(7) .* time{1});
        output{3} = exp(-pars(9) .* time{1});

        solution_l = ode23s(@ode,[0 tmax{1}],...
                       [0 1 1 1 1],...
                       [],...
                       pars,...
                       Rif{1});
        solution_h = ode23s(@ode,[0 tmax{2}],...
                       [0 1 1 1 1],...
                       [],...
                       pars,...
                       Rif{2});           
        %%% mRNA kinetics for measured proteins %%%
        % low
        mm = 4;
        for ii = 2:5
            output{mm} = deval(solution_l,time{2},ii);
            mm = mm + 1;
        end
        % high
        for ii = 2:5
            output{mm} = deval(solution_h,time{3},ii);
            mm = mm + 1;
        end   
    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x,pars,Xint)        
        dxdt = zeros(5,1);

        dxdt(1) = pars(1)*(1-x(1))*Xint*exp(-pars(2)*t) - pars(3)*x(1);    % activated PXR
        dxdt(2) = pars(4)*x(1) + pars(5)*(1-x(2));                         % CYP3A4
        dxdt(3) = pars(6)*x(1) + pars(7)*(1-x(3));                         % CYP2C9
        dxdt(4) = pars(8)*x(1) + pars(9)*(1-x(4));                         % CYP2B6
        dxdt(5) = pars(10)*x(1) + pars(11)*(1-x(5));                       % MDR1
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end