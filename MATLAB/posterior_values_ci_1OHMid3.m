function [] = posterior_values_ci_1OHMid3()

    %%% command:
    %%% posterior_values_ci_1OHMid3()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    metadata = load('datamat_activity.mat');
    
    data    = metadata.data_3{2};
    stdev   = mean(metadata.stdev_3{2});
    substrate_total = 800; 
    time    = metadata.time_act;
    tmax    = metadata.tmax_act;
    tt = 4;
    
    MLpars = readmatrix('maxLikValues.txt');
    mlpars = MLpars(2:end,2);
    Rif = 10;
    
    chains = [];
    for jj = 1:5
        jjchains = load(strcat('./chains_1OHMid3','/chains_1OHMid3_',num2str(jj),'.mat'));
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
    save(strcat('./mlvalues_1OHMid3XXX.mat'),'mlvalue_array');
            
    parsmean = mean(chains);
    parsmedian = median(chains);
    parsquantile = [];
    for kk = 1:size(chains,2)
        parsquantile = [parsquantile; quantile(chains(:,kk),[0.05 0.95])];
    end
    
    
    a{1} = 'k_{met,cyp3a4}';  
    labs = {a{1}}; 
        
    names = {'parameter','mean','median','95CI lower', '95CI upper'};
    
    fid = fopen('posteriorValues_1OHMid3XXX.txt','w');

    fprintf(fid, '%2s %2s %2s %2s %2s\n', names{:});
    for kk = 1:length(labs)
        fprintf(fid,'%0s %.8f %.8f %.8f %.8f\n',labs{kk},parsmean(kk),parsmedian(kk),parsquantile(kk,1),parsquantile(kk,2));
    end
    fclose(fid);
    
    fidml = fopen('maxLikValues_1OHMid3.txt','w');

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
        OH = substrate_total.*(1 - exp(-pars*tt*CYP));
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