function [results] = mcmc_1OHMid3(nstart,nend,pars)

    metadata = load('datamat_activity.mat');
    
    data    = metadata.data_3{1};
    stdev   = mean(metadata.stdev_3{1});
    substrate_total = 800;
    time    = metadata.time_act;
    tmax    = metadata.tmax_act;
    tt = 4;
    
    MLpars = readmatrix('maxLikValues.txt');
    mlpars = MLpars(2:end,2);
    Rif = 10;

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
        save(strcat('./chains_1OHMid3/chains_1OHMid3_',num2str(iround),'.mat'),'chains');
        
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
        OH = substrate_total.*(1 - exp(-pars*tt*CYP));
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

        dxdt(1) = mlpars(1)*(1-x(1))*Rif*exp(-mlpars(2)*t) - mlpars(3)*x(1);   % activated PXR
        dxdt(2) = mlpars(4)*x(1) + mlpars(5)*(1-x(2));                         % CYP3A4  
    end

end