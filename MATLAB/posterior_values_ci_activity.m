function [] = posterior_values_ci_activity(cyp)

    %%% command:
    %%% posterior_values_ci_activity(cyp), where cyp = 1 or 2 or 3
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    metadata = load('datamat_activity.mat');
    
    switch cyp
        case 1
            data    = metadata.data_act{1};
            stdev   = mean(metadata.stdev_act{1});
            substrate_total = 10*80; 
            label = '3a4';
            a{1} = 'k_{met,cyp3a4}';  
            labs = {a{1}}; 
        case 2
            data    = metadata.data_act{2};
            stdev   = mean(metadata.stdev_act{2});
            substrate_total = 10*80; 
            label = '2c9';
            a{1} = 'k_{met,cyp2c9}';  
            labs = {a{1}}; 
        case 3
            data    = metadata.data_act{3};
            stdev   = mean(metadata.stdev_act{3});
            substrate_total = 40*80; 
            label = '2b6';
            a{1} = 'k_{met,cyp2b6}';  
            labs = {a{1}}; 
    end
    
    time    = metadata.time_act;
    tmax    = metadata.tmax_act;
    tt = 4;
    
    MLpars = readmatrix('maxLikValues.txt');
    mlpars = MLpars(2:end,2);
    Rif = 10;
    
    chains = [];
    for jj = 1:5
        jjchains = load(strcat('./chains_activity_cyp',label,'/chains_activity_cyp',label,'_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        burnin = size(jjchains,1)/2;
        chains = [chains; jjchains(burnin+1:end,:)];
    end
    
    mlvalue = -99999.99;
    mlparams = zeros(1,size(chains,2));
    mlvalue_array = [];
    for ll = 1:size(chains,1)
        if mod(ll,100)==0
            disp(ll);
        end
        llvalue = loglike(chains(ll,:));
        mlvalue_array = [mlvalue_array; llvalue];
        if llvalue > mlvalue
            mlvalue = llvalue;
            mlparams = chains(ll,:);
        end
    end
    %%% save chains %%%
    save(strcat('./mlvalues_activity_cyp',label,'.mat'),'mlvalue_array');
            
    parsmean = mean(chains);
    parsmedian = median(chains);
    parsquantile = [];
    for kk = 1:size(chains,2)
        parsquantile = [parsquantile; quantile(chains(:,kk),[0.05 0.95])];
    end
            
    names = {'parameter','mean','median','95CI lower', '95CI upper'};
    
    fid = fopen(strcat('posteriorValues_activity_cyp',label,'.txt'),'w');

    fprintf(fid, '%2s %2s %2s %2s %2s\n', names{:});
    for kk = 1:length(labs)
        fprintf(fid,'%0s %.8f %.8f %.8f %.8f\n',labs{kk},parsmean(kk),parsmedian(kk),parsquantile(kk,1),parsquantile(kk,2));
    end
    fclose(fid);
    
    fidml = fopen(strcat('maxLikValues_activity_cyp',label,'.txt'),'w');

    fprintf(fidml, '%2s %.8f \n', 'MLvalue', mlvalue);
    for kk = 1:length(labs)
        fprintf(fidml,'%0s %.8f\n',labs{kk},mlparams(kk));
    end
    fclose(fidml);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars)        
        CYP = model();   
        value = 0;
        OH = substrate_total*(1 - exp(-pars*tt*CYP));
        for ii=1:size(data,2)  
            idx=~isnan(data(:,ii));
            value = value+sum(lognormpdf(data(idx,ii),OH(1,ii),stdev));     
        end
    end

    %%% output of the model %%%
    function [output] = model()         
        solution = ode23s(@ode,[0 tmax],...
                       [0 1 ],...
                       []);          
        output = deval(solution,time,2);          
    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x)        
        dxdt = zeros(2,1);

        dxdt(1) = mlpars(1)*(1-x(1))*Rif*exp(-mlpars(2)*t) - mlpars(3)*x(1);   % activated PXR
        dxdt(2) = mlpars(4)*x(1) + mlpars(5)*(1-x(2));                         % CYP3A4  
    end

end