function [results] = mcmc_activity(nstart,nend,pars,cyp)

    metadata = load('datamat_activity.mat');
    MLpars = readmatrix('maxLikValues.txt');
    mlpars = MLpars(2:end,2);
    
    switch cyp
        case 1
            data            = metadata.data_act{1};
            stdev           = mean(metadata.stdev_act{1});
            k_cyp           = mlpars(4);
            k_cypdeg        = mlpars(5);
            substrate_total = 10*80; 
            label = '3a4';
        case 2
            data            = metadata.data_act{2};
            stdev           = mean(metadata.stdev_act{2});
            k_cyp           = mlpars(6);
            k_cypdeg        = mlpars(7);
            substrate_total = 10*80; 
            label = '2c9';
        case 3
            data            = metadata.data_act{3};
            stdev           = mean(metadata.stdev_act{3});
            k_cyp           = mlpars(8);
            k_cypdeg        = mlpars(9); 
            substrate_total = 40*80; 
            label = '2b6';
    end

    k_pxr           = mlpars(1);
    k_r             = mlpars(2);
    k_pxrdeg        = mlpars(3);
    time            = metadata.time_act;
    tmax            = metadata.tmax_act;
    tt              = 4;
    
    Rif             = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iround=nstart:nend
        nSamples = 5*10^3;
        chains = slicesample(pars(iround,:),nSamples,...
        'burnin',0,...
        'thin',1,...
        'logpdf',@logprob);
        results{iround} = chains;
        
        %%% save chains %%%
        save(strcat('./chains_activity_cyp',label,'/chains_activity_cyp',label,'_',num2str(iround),'.mat'),'chains');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% helper functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
    end

    %%% log-probability %%%
    function [value] = logprob(pars)  
        lp = prior(pars);
        if ~isfinite(lp)
            value = -inf;
        else
            value = lp+loglike(pars);
        end   
    end

    %%% priors %%%
    function [flag] = prior(pars)
        if all(pars>0)
            flag = 0;
        else
            flag = -inf;
        end    
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars)        
        CYP = model();
        value = 0;
        OH = substrate_total*(1 - exp(-pars*tt*CYP));
        for jj=1:size(data,2)  
            idx=~isnan(data(:,jj));
            value = value+sum(lognormpdf(data(idx,jj),OH(1,jj),stdev));    
        end
    end

    %%% output of the model %%%
    function [output] = model()         
        solution = ode23s(@ode,[0 tmax],...
                       [0 1],...
                       []);          
        output = deval(solution,time,2);        
    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x)        
        dxdt = zeros(2,1);

        dxdt(1) = k_pxr*(1-x(1))*Rif*exp(-k_r*t) - k_pxrdeg*x(1);          % activated PXR
        dxdt(2) = k_cyp*x(1) + k_cypdeg*(1-x(2));                          % CYP
    end

end