function [sim, rt, rt_X] = KLS2_simulate_trans(m,s,r,rt,rt_X,trans,Yguess,Pguess,sim_0)

%% Compute stationary measure across states

    et = trans.et;
    et_X = trans.et_X;
    N = trans.N;
    sim = sim_0;
    
    %Measure of firms type 1
    M = KLS2_measure_trans(r,N,trans.apt_ind,trans.measure);   
    sim.M = M;
    
    %Measure of firms type 2
    M_X = KLS2_measure_trans(r,N,trans.apt_ind_X,trans.measure_X);   
    sim.M_X = M_X;
    
    
    sim.measure_end = squeeze(sim.M(:,:,N));
    sim.measure_X_end =  squeeze(sim.M_X(:,:,N)); 
    sim.measure_all_end =  sim.measure_end*(1-m.Xshare)+m.Xshare*sim.measure_X_end ;
     
for n=2:N
    
    %Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pfv(n)/Pguess(n); 
    m.xi_real_lag=m.Pfv(n-1)/Pguess(n-1); 
    sim.xi_real(n,1) = m.xi_real;
    
	
	if s.high_z_ForDebt == 0
		if n==2
			FXdebt =m.lambda_period_2 + (1-m.lambda_period_2)*(m.xi_real/m.xi_real_lag); 
			FXdebt_X=m.lambda_period_2_X+(1-m.lambda_period_2_X)*(m.xi_real/m.xi_real_lag);        
		else
		
			FXdebt =m.lambda + (1-m.lambda)*(m.xi_real/m.xi_real_lag); 
			FXdebt_X=m.lambda_X+(1-m.lambda_X)*(m.xi_real/m.xi_real_lag);
		
			if s.UIPdeviations == 1;
			
				FXdebt =m.lambda_v(n) + (1-m.lambda_v(n))*(m.xi_real/m.xi_real_lag); 
				FXdebt_X=m.lambda_X_v(n)+(1-m.lambda_X_v(n))*(m.xi_real/m.xi_real_lag);
				
			end
		
		end
	
	elseif s.high_z_ForDebt == 1
	
	    if n==2
			FXdebt =m.lambda_period_2_mat + (1-m.lambda_period_2_mat).*(m.xi_real/m.xi_real_lag); 
			FXdebt_X=m.lambda_X_period_2_mat+(1-m.lambda_X_period_2_mat).*(m.xi_real/m.xi_real_lag);        
		else
			FXdebt =m.lambda_mat + (1-m.lambda_mat).*(m.xi_real/m.xi_real_lag); 
			FXdebt_X=m.lambda_X_mat+(1-m.lambda_X_mat).*(m.xi_real/m.xi_real_lag);
		end
	
	end
   
    measure = squeeze(sim.M(:,:,n));
    measure_X = squeeze(sim.M_X(:,:,n));
    
    measure_lag = squeeze(sim.M(:,:,n-1));
    measure_X_lag = squeeze(sim.M_X(:,:,n-1));    
    
%% Constrained

    if n>2
        
    % Extensive margin:
    
        % unconstrained export policies:
        e_unc = rt{n}.pi_u_x-m.F_base>= rt{n}.pi_u_nx;
        if s.model == 1  || s.model == 2  || s.model == 4
            e_X_unc = rt_X{n}.pi_u_x-m.F_X_base>= rt_X{n}.pi_u_nx;
        elseif s.model == 3
            e_X_unc = rt_X{n}.pi_u_x>= rt_X{n}.pi_u_nx;
        end
        
        % Measure of firms constrained on the extensive margin:
        sim.ext_const(n,1) = sum(sum(measure.*(e_unc-rt{n}.e)));
        sim.ext_const_X(n,1) = sum(sum(measure_X.*(e_X_unc-rt_X{n}.e)));
        sim.ext_const_all(n,1) = sim.ext_const(n,1)*(1-m.Xshare) + sim.ext_const_X(n,1)*m.Xshare;        

    % Intensive margin:        

        % error margin (the same expressions can be numericallu slghtly different; matters for indicators; 'ind')
        err = 1e-10;
        
        % in terms of capital (indicator)
        sim.int_const_k_ind(n,1) = sum(sum(measure.*( rt{n}.e.*((rt{n}.k_u_x-err)>rt{n}.k) + (1-rt{n}.e).*((rt{n}.k_u_nx-err)>rt{n}.k) )));
        sim.int_const_k_ind_X(n,1) = sum(sum(measure_X.*( rt_X{n}.e.*((rt_X{n}.k_u_x-err)>rt_X{n}.k)+(1-rt_X{n}.e).*((rt_X{n}.k_u_nx-err)>rt_X{n}.k) )));
        sim.int_const_k_ind_all(n,1) = sim.int_const_k_ind(n,1)*(1-m.Xshare) + sim.int_const_k_ind_X(n,1)*m.Xshare;
        
        % in terms of capital (percentages)
        sim.int_const_k_prc(n,1) = sum(sum(measure.*( rt{n}.e.*(rt{n}.k./rt{n}.k_u_x)+(1-rt{n}.e).*(rt{n}.k./rt{n}.k_u_nx) )));
        sim.int_const_k_prc_X(n,1) = sum(sum(measure_X.*( rt_X{n}.e.*(rt_X{n}.k./rt_X{n}.k_u_x)+(1-rt_X{n}.e).*(rt_X{n}.k./rt_X{n}.k_u_nx) )));
        sim.int_const_k_prc_all(n,1) = sim.int_const_k_prc(n,1)*(1-m.Xshare) +  sim.int_const_k_prc_X(n,1)*m.Xshare; 
        
        % in terms of profits (indicator)
        sim.int_const_pi_ind(n,1) = sum(sum(measure.*(rt{n}.e.*((rt{n}.pi_u_x-err)>rt{n}.pi)+(1-rt{n}.e).*((rt{n}.pi_u_nx-err)>rt{n}.pi))));
        sim.int_const_pi_ind_X(n,1) = sum(sum(measure_X.*( rt_X{n}.e.*((rt_X{n}.pi_u_x-err)>rt_X{n}.pi)+(1-rt_X{n}.e).*((rt_X{n}.pi_u_nx-err)>rt_X{n}.pi) )));
        sim.int_const_pi_ind_all(n,1) = sim.int_const_pi_ind(n,1)*(1-m.Xshare) + sim.int_const_pi_ind_X(n,1)*m.Xshare;
        
        % in terms of capital (percentages)
        sim.int_const_pi_prc(n,1) = sum(sum(measure.*(rt{n}.e.*(rt{n}.pi./rt{n}.pi_u_x)+(1-rt{n}.e).*(rt{n}.pi./rt{n}.pi_u_nx))));
        sim.int_const_pi_prc_X(n,1) = sum(sum(measure_X.*( rt_X{n}.e.*(rt_X{n}.pi./rt_X{n}.pi_u_x)+(1-rt_X{n}.e).*(rt_X{n}.pi./rt_X{n}.pi_u_nx) )));
        sim.int_const_pi_prc_all(n,1) = sim.int_const_pi_prc(n,1)*(1-m.Xshare) + sim.int_const_pi_prc_X(n,1)*m.Xshare;        
        
        
    end

%% Investment

    if n>1 && n<N        
        
        sim.inv_agg(n,1) = (sum(sum(squeeze(sim.M(:,:,n+1)).*rt{n+1}.k)) - (1-m.delta_v(n))*sum(sum(squeeze(sim.M(:,:,n)).*rt{n}.k)))*(1-m.Xshare);
        sim.inv_agg_X(n,1) = (sum(sum(squeeze(sim.M_X(:,:,n+1)).*rt_X{n+1}.k)) - (1-m.delta_v(n))*sum(sum(squeeze(sim.M_X(:,:,n)).*rt_X{n}.k)))*m.Xshare;
        sim.inv_agg_all(n,1)=sim.inv_agg(n,1)+sim.inv_agg_X(n,1);
        
    elseif n==N %For last period (n=N)   
        % We assume that we reached steady state when n=N
        sim.inv_agg(n,1) = (sum(sum(squeeze(sim.M(:,:,n)).* (m.delta_v(n)*rt{n}.k) )) )*(1-m.Xshare);
        sim.inv_agg_X(n,1) = (sum(sum(squeeze(sim.M_X(:,:,n)).* (m.delta_v(n)*rt_X{n}.k) )) )*m.Xshare;
        sim.inv_agg_all(n,1)=sim.inv_agg(n,1)+sim.inv_agg_X(n,1);
        
    end

%% Exporters and non-exporters
        
    sim.share_x(n,1) = sum(sum(measure.*et(:,:,n)))*(1-m.Xshare) + m.Xshare*sum(sum(measure_X.*et_X(:,:,n))); 
    sim.share_x_composition(n,1) =m.Xshare*sum(sum(measure_X.*et_X(:,:,n)))/sim.share_x(n,1); % Share of exporters type 2
    sim.share_nx(n,1) = sum(sum(measure.*(1-et(:,:,n))))*(1-m.Xshare) + sum(sum(measure_X.*(1-et_X(:,:,n))))*m.Xshare;    
    sim.share_nx_composition(n,1) = sum(sum(measure_X.*(1-et_X(:,:,n))))*m.Xshare /sim.share_nx(n,1); % Share of non-exporters type 2
    
%% Compute price and quantity indexes

%Domestic sales
    if s.model==1  || s.model == 2  || s.model == 4
        Pd_sigma = sum( measure(rt{n}.pd>0).*(rt{n}.pd(rt{n}.pd>0).^(1-m.sigma)) )*(1-m.Xshare)...
                   + sum( measure_X(rt_X{n}.pd>0).*(rt_X{n}.pd(rt_X{n}.pd>0).^(1-m.sigma)) )*(m.Xshare);
        Yd_sigma = sum( measure(rt{n}.yd>0).*(rt{n}.yd(rt{n}.yd>0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare)...
                   +sum( measure_X(rt_X{n}.yd>0).*(rt_X{n}.yd(rt_X{n}.yd>0).^((m.sigma-1)/m.sigma)) )*(m.Xshare);

        sim.Pd(n,1) = Pd_sigma^(1/(1-m.sigma));
        sim.Yd(n,1) = Yd_sigma^(m.sigma/(m.sigma-1));     
        sim.PdYd(n,1) = sum( measure(rt{n}.pd>0).*(rt{n}.pd(rt{n}.pd>0).*rt{n}.yd(rt{n}.pd>0)))*(1-m.Xshare)...
                        + sum( measure_X(rt_X{n}.pd>0).*(rt_X{n}.pd(rt_X{n}.pd>0).*rt_X{n}.yd(rt_X{n}.pd>0)))*(m.Xshare);
    elseif s.model==3
        Pd_sigma = sum( measure(rt{n}.pd>0).*(rt{n}.pd(rt{n}.pd>0).^(1-m.sigma)) )*(1-m.Xshare);
        Yd_sigma = sum( measure(rt{n}.yd>0).*(rt{n}.yd(rt{n}.yd>0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare);

        sim.Pd(n,1) = Pd_sigma^(1/(1-m.sigma));
        sim.Yd(n,1) = Yd_sigma^(m.sigma/(m.sigma-1));     
        sim.PdYd(n,1) = sum(measure(rt{n}.pd>0).*(rt{n}.pd(rt{n}.pd>0).*rt{n}.yd(rt{n}.pd>0)))*(1-m.Xshare);
    end

    if s.model==3 
    %Exporters vs non-exporters among firms that can sell in both markets
        Yd_sigma_NonX = sum( measure(rt{n}.yd>0 & rt{n}.e==0).*(rt{n}.yd(rt{n}.yd>0 & rt{n}.e==0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare);
        sim.Yd_NonX(n,1) = Yd_sigma_NonX^(m.sigma/(m.sigma-1));    
        sim.PdYd_NonX(n,1) = sum( measure(rt{n}.pd>0 & rt{n}.e==0).*(rt{n}.pd(rt{n}.pd>0 & rt{n}.e==0).*rt{n}.yd(rt{n}.pd>0 & rt{n}.e==0)))*(1-m.Xshare);

        Yd_sigma_X = sum( measure(rt{n}.yd>0 & rt{n}.e==1).*(rt{n}.yd(rt{n}.yd>0 & rt{n}.e==1).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare);

        sim.Yd_X(n,1) = Yd_sigma_X^(m.sigma/(m.sigma-1));    
        sim.PdYd_X(n,1) = sum( measure(rt{n}.pd>0 & rt{n}.e==1).*(rt{n}.pd(rt{n}.pd>0 & rt{n}.e==1).*rt{n}.yd(rt{n}.pd>0 & rt{n}.e==1)))*(1-m.Xshare);
    end
        
%Exports    
    Px_sigma = sum( measure(rt{n}.pf>0).*(rt{n}.pf(rt{n}.pf>0).^(1-m.sigma)) )*(1-m.Xshare)...
        + sum( measure_X(rt_X{n}.pf>0).*(rt_X{n}.pf(rt_X{n}.pf>0).^(1-m.sigma)) )*m.Xshare;
    
    Yx_sigma = sum( measure(rt{n}.yf>0).*(rt{n}.yf(rt{n}.yf>0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare);
    Yx_sigma_X = sum( measure_X(rt_X{n}.yf>0).*(rt_X{n}.yf(rt_X{n}.yf>0).^((m.sigma-1)/m.sigma)) )*m.Xshare;
    Yx_sigma_all = Yx_sigma_X + Yx_sigma;
    
    Px_sigma_zero = Px_sigma<=0;
    Yx_sigma_zero_all = Yx_sigma_all<=0;
    sim.Yx_sigma(n,1) = (1-Yx_sigma_zero_all)*(((1-Yx_sigma_zero_all)*Yx_sigma_all + Yx_sigma_zero_all)^(m.sigma/(m.sigma-1)));                    
   
    Yx_sigma_zero = Yx_sigma<=0;
    Yx_sigma_zero_X = Yx_sigma_X<=0;
    sim.Yx(n,1) = (1-Yx_sigma_zero)*(((1-Yx_sigma_zero)*Yx_sigma + Yx_sigma_zero)^(m.sigma/(m.sigma-1)));       
    sim.Yx_X(n,1) = (1-Yx_sigma_zero_X)*(((1-Yx_sigma_zero_X)*Yx_sigma_X + Yx_sigma_zero_X)^(m.sigma/(m.sigma-1))); 
    
    sim.Yx_sumall(n,1) = (1-m.Xshare)*sum(measure(rt{n}.yf>0).*rt{n}.yf(rt{n}.yf>0)) + (m.Xshare)*sum(measure_X(rt_X{n}.yf>0).*rt_X{n}.yf(rt_X{n}.yf>0));
    
    sim.Px(n,1) = m.xi_real*(1-Px_sigma_zero)*(((1-Px_sigma_zero)*Px_sigma + Px_sigma_zero)^(1/(1-m.sigma)));
    sim.PxYx(n,1) = m.xi_real*sum( measure(rt{n}.pf>0).*( rt{n}.pf(rt{n}.pf>0).*rt{n}.yf(rt{n}.pf>0)))*(1-m.Xshare);    
    sim.PxYx_X(n,1) = m.xi_real*sum( measure_X(rt_X{n}.pf>0).*( rt_X{n}.pf(rt_X{n}.pf>0).*rt_X{n}.yf(rt_X{n}.pf>0)))*m.Xshare;

    sim.PxYx_all(n,1) = sim.PxYx_X(n,1) + sim.PxYx(n,1);

    sim.PxYx_all_USD(n,1) = sim.PxYx_all(n,1)/m.xi_real;
    sim.PxYx_USD(n,1) = sim.PxYx(n,1)/m.xi_real;
    sim.PxYx_X_USD(n,1) = sim.PxYx_X(n,1)/m.xi_real;           

    
%Imports for a small open economy    
    sim.Pm(n,1) =  m.xi_real*m.pm_v(n); %in units of domestic final goods.
    sim.Ym(n,1) = m.omega_m^m.sigma*(m.A_v(n)^(m.sigma-1))*Yguess(n)*(sim.Pm(n)^(-m.sigma));
    sim.Pm_sigma = (m.omega_m^m.sigma)*(sim.Pm(n)^(1-m.sigma));
    sim.Ym_sigma = m.omega_m*sim.Ym(n,1)^((m.sigma-1)/m.sigma);
    sim.PmYm(n,1) = sim.Pm(n)*sim.Ym(n,1);   
      
    
%% Compute aggregate variables needed to evaluate market clearing conditions

%Labor and capital    
    sim.K(n,1) = sum(sum(measure.*rt{n}.k))*(1-m.Xshare);
    sim.K_X(n,1) = m.Xshare*sum(sum(measure_X.*rt_X{n}.k));
    sim.K_all(n,1) = sim.K(n,1) + sim.K_X(n,1);
    
    sim.N(n,1) = sum(sum(measure.*rt{n}.n))*(1-m.Xshare);
    sim.N_X(n,1) = m.Xshare*sum(sum(measure_X.*rt_X{n}.n));
    
%Fixed costs   
    if s.model==1  || s.model == 2  || s.model == 4
        sim.FC(n,1) = sum(sum(measure.*rt{n}.e)).*m.F_base*(1-m.Xshare) + sum(sum(measure_X.*rt_X{n}.e)).*m.F_X_base*(m.Xshare);
    elseif s.model==3
        sim.FC(n,1) = sum(sum(measure.*rt{n}.e)).*m.F_base*(1-m.Xshare);
    end

    % Total labor
    sim.N_all(n,1) = sim.N(n,1) + sim.N_X(n,1)+sim.FC(n,1)*(1-s.fcost_fgoods);    
    
    
%Final good output and price
    %Tradable good 
    sim.Pcpi(n,1) = (sim.Pd(n,1)^(1-m.sigma) + max(sim.Pm_sigma,0))^(1/(1-m.sigma));
    sim.Ycpi(n,1) = m.A_v(n)*(sim.Yd(n,1)^((m.sigma-1)/m.sigma) + max(sim.Ym_sigma,0))^(m.sigma/(m.sigma-1));        
    
    %Final good
    sim.Y(n,1) = sim.Ycpi(n,1); 
    sim.P(n,1) = Pguess(n);     
        
%Consumption    
    sim.C(n,1) = sum(sum(measure.*rt{n}.c))*(1-m.Xshare);
    sim.C_X(n,1) = m.Xshare*sum(sum(measure_X.*rt_X{n}.c)); 
    sim.C_all(n,1) = sim.C(n,1) + sim.C_X(n,1);
    
%% Market clearing conditions

%Labor
    sim.n_supply(n,1) = 1;
    sim.n_demand(n,1) = sim.N_all(n,1); %Total labor demand = Demand from entrepreneurs  + fixed export costs + non-tradable goods     
    sim.mc_n(n,1) = 10*log((1+sim.n_demand(n,1))/(1+sim.n_supply(n,1)));        
                                     
%Assets 
%This market clearing condition is the result of reformulating the debt market clearing condition
%In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below.   
    sim.a_supply(n,1) = sum(sum(measure.*r.a_grid_mat))*(1-m.Xshare)+m.Xshare*sum(sum(measure_X.*r.a_grid_mat)); 
    sim.a_demand(n,1) = sim.K_all(n,1);
    sim.mc_a(n,1) = 10*log((1+sim.a_demand(n,1))/(1+sim.a_supply(n,1)));          

%Final goods   
%Total production of final goods is used for: (i) consumption, (ii) investment, (iii) fixed costs
%Consumption and fixed costs expenditures are straightforward to compute
%Aggregate investment is derived as follows: nu*k_new + integral[k' - (1-delta)k] 
%Now, integral[(1-delta)k]=(1-delta)kagg and [nu*k_new + integral(k')]=kagg
%Then, aggregate investment equals delta*kagg
    sim.y_supply(n,1) = sim.Y(n,1);
    sim.y_demand(n,1) = sim.C_all(n,1)+ sim.inv_agg_all(n,1)+ sim.FC(n,1)*s.fcost_fgoods ;  
    
    sim.mc_y(n,1) = 10*log((1+sim.y_demand(n,1))/(1+sim.y_supply(n,1)));
    
%Beliefs about tradable good
%To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
%In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
    sim.mc_p_belief(n,1) = 10*log((1+sim.Pcpi(n,1))/(1+m.A_v(n)));        
    sim.mc_y_belief(n,1) = 10*log((1+sim.Y(n,1))/(1+Yguess(n)));                                                       
    
%% Display statistics    
 
%Setup
    rt{n}.d = (1+m.rv(n))*(rt{n}.k - r.a_grid_mat);  % This is the 'expected debt' at period n+1
    rt_X{n}.d = (1+m.rv(n))*(rt_X{n}.k - r.a_grid_mat);  % This is the 'expected debt' at period n+1
    
    rt{n}.d_repayment =FXdebt.*(1+m.rv(n)).*(rt{n}.k - r.a_grid_mat);  % This is the amount you have to repay at n+1
    rt_X{n}.d_repayment = FXdebt_X.*(1+m.rv(n)).*(rt_X{n}.k - r.a_grid_mat);  % This is the amount you have to repay at n+1
    
	
%     rt{n}.d_cost =(FXdebt-1)*(1+m.rv(n))*(rt{n}.k - rt{n}.a);
%     rt_X{n}.d_cost = (FXdebt_X-1)*(1+m.rv(n))*(rt_X{n}.k - rt_X{n}.a);     
 
    sim.credit_all(n,1) = (sum(sum(measure.*max(rt{n}.d,0)))*(1-m.Xshare)+m.Xshare*sum(sum(measure_X.*max(rt_X{n}.d,0))) )./ (1+m.rv(n)) ;   %Amount borrowed in current period (sum of positive d / 1+r), domestic units
    sim.credit_repayment_all(n,1) = sum(sum(measure.*max(rt{n}.d_repayment,0)))*(1-m.Xshare)+m.Xshare*sum(sum(measure_X.*max(rt_X{n}.d_repayment,0))); %Credit repayment in current period (sum of positive d * FXdebt), domestic units  
    
    sim.d_agg(n,1) = sum(sum(measure.*rt{n}.d))*(1-m.Xshare)+m.Xshare*sum(sum(measure_X.*rt_X{n}.d));
    sim.d_agg_repayment(n,1) = sum(sum(measure.*rt{n}.d_repayment))*(1-m.Xshare)+m.Xshare*sum(sum(measure_X.*rt_X{n}.d_repayment));
    
    %Low export cost firms (type 2)
    sim.credit_X(n,1) = m.Xshare*sum(sum(measure_X.*max(rt_X{n}.d,0))); 
    sim.credit_repayment_X(n,1) = m.Xshare*sum(sum(measure_X.*max(rt_X{n}.d_repayment,0))); 
    sim.d_agg_X(n,1)=m.Xshare*sum(sum(measure_X.*rt_X{n}.d));
    sim.d_agg_X_repayment(n,1)=m.Xshare*sum(sum(measure_X.*rt_X{n}.d_repayment));


	if s.high_z_ForDebt == 0
		sim.r_tilde(n,1) = (1+m.rv(n))*FXdebt - 1;
		sim.r_tilde_X(n,1) = (1+m.rv(n))*FXdebt_X - 1;
	elseif s.high_z_ForDebt == 1
			sim.r_tilde(n,1) = sum(sum(measure.*(1+m.rv(n)).*FXdebt)) - 1;
		sim.r_tilde_X(n,1) = sum(sum(measure_X.*(1+m.rv(n)).*FXdebt_X)) - 1;
   	end	
		
	if n==2
		sim.r_effective_avg(n,1) = (1-m.Xshare)*(1+sim.r_tilde(1))*(1 + sum(sum(sum(measure.*rt{1}.lm))))...
								+ m.Xshare*(1+sim.r_tilde_X(1))*(1 + sum(sum(sum(measure_X.*rt_X{1}.lm))));     
	else
		sim.r_effective_avg(n,1) = (1-m.Xshare)*(1+sim.r_tilde(n))*(1 + sum(sum(sum(measure.*rt{n}.lm))))...
								+ m.Xshare*(1+sim.r_tilde_X(n))*(1 + sum(sum(sum(measure_X.*rt_X{n}.lm))));     
	end

	
    
    sim.GDP(n,1) = sim.PdYd(n,1) + sim.PxYx_all(n,1);    
    sim.TFP(n,1) = sim.GDP(n,1)/((sim.K_all(n,1)^m.alpha)*(sim.N_all(n,1)^(1-m.alpha)));
    
    sim.NX_GDP(n,1) = (sim.PxYx_all(n,1)-sim.PmYm(n,1))/sim.GDP(n,1);    
    
    % Is this correct?
    sim.NFA_GDP(n,1) = -sim.d_agg_repayment(n,1)/sim.GDP(n,1); %Is this the right way to measure it? Or should we take into account changes in the RER? (as in d_agg_repayment)
    
    sales_x = m.xi_real*rt{n}.pf.*rt{n}.yf;
    sales_x_X = m.xi_real*rt_X{n}.pf.*rt_X{n}.yf;    
        
    sales_d = rt{n}.pd.*rt{n}.yd;
    sales = sales_d+sales_x; %r.pf is denominated in foreign currency, so we adjust it    
    
    x_share = sales_x./sales;    
      
    sales_d_X =rt_X{n}.pd.*rt_X{n}.yd;
    sales_X = sales_d_X+sales_x_X; 
    x_share_X = sales_x_X./sales_X;     

    if s.model==1  || s.model == 2  || s.model == 4
        sim.sales_avg_nx(n,1) = ((1-m.Xshare)*sum(sum(measure.*sales.*(1-rt{n}.e))))+(m.Xshare)*sum(sum(measure_X.*sales_X.*(1-rt_X{n}.e)) )/sim.share_nx(n,1); 
    elseif s.model==3
        sim.sales_avg_nx(n,1) = (1-m.Xshare)*sum(sum(measure.*sales.*(1-rt{n}.e)))/sim.share_nx(n,1);  
    end
    sim.sales_avg_x(n,1) = ((1-m.Xshare)*sum(sum(measure.*sales.*rt{n}.e))+m.Xshare* sum(sum(measure_X.*sales_X.*rt_X{n}.e))   )/sim.share_x(n,1);  
    sim.exp_avg_x(n,1) = ((1-m.Xshare)*sum(sum(measure.*sales_x.*rt{n}.e))+m.Xshare* sum(sum(measure_X.*sales_x_X.*rt_X{n}.e))   )/sim.share_x(n,1);  

%    
    
%Statistics to be displayed    
     sim.X_GDP(n,1) = sim.PxYx_all(n,1)/sim.GDP(n,1);     
     sim.credit_gdp(n,1) = sim.credit_repayment_all(n,1)/sim.GDP(n,1); % is this correct?
     sim.xpremium_sales(n,1) = sim.sales_avg_x(n,1)/sim.sales_avg_nx(n,1);

     %Average export intensity among exporters
     sim.x_share_avg = ( (1-m.Xshare)*sum(sum(measure(x_share>0).*x_share(x_share>0))) +m.Xshare* sum(sum(measure_X(x_share_X>0).*x_share_X(x_share_X>0)))  )/sim.share_x(n,1);
     
     sim.X_D(n,1) = sim.PxYx_all(n,1)/sim.PdYd(n,1);       
  
%% Real variables
%Laspeyres with prices fixed at pre-devaluation values

%Real GDP
    xi_real_ss = m.Pfv(1)/m.Pt(1);
    sales_d = rt{1}.pd.*rt{n}.yd;
    sales_x = xi_real_ss*rt{1}.pf.*rt{n}.yf;
    sales = sales_d + sales_x;
    
    sales_d_X = rt_X{1}.pd.*rt_X{n}.yd;
    sales_x_X = xi_real_ss*rt_X{1}.pf.*rt_X{n}.yf;
    sales_X = sales_d_X + sales_x_X;    
    
    sim.GDP_Laspeyres(n,1) = (1-m.Xshare)*sum(sum(measure.*sales)) + m.Xshare*sum(sum(measure_X.*sales_X));
    sim.X_Laspeyres(n,1) = (1-m.Xshare)*sum(sum(measure.*sales_x)) + m.Xshare*sum(sum(measure_X.*sales_x_X));
    sim.D_Laspeyres(n,1) = (1-m.Xshare)*sum(sum(measure.*sales_d)) + m.Xshare*sum(sum(measure_X.*sales_d_X));
    
    sim.X_Laspeyres_Type1(n,1) = (1-m.Xshare)*sum(sum(measure.*sales_x));
    sim.X_Laspeyres_Type2(n,1) = m.Xshare*sum(sum(measure_X.*sales_x_X));
    
%% Solution-related statistics

    sim.a_min_share(n,1) = sum(measure(1,:));
    sim.a_max_share(n,1) = sum(measure(end,:));
    
end % Transition Loop

sim.GDP_Laspeyres(1,1) = sim.GDP(1);
sim.X_Laspeyres(1,1) = sim.PxYx_all(1);
sim.D_Laspeyres(1,1) = sim.PdYd(1);

sim.X_Laspeyres_Type1(1,1) = sim.PxYx(1);
sim.X_Laspeyres_Type2(1,1) = sim.PxYx_X(1);

end

