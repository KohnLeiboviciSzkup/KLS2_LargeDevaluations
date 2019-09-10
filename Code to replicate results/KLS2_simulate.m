function [sim,r,r_X,s] = KLS2_simulate(m,s,r,r_X)

%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag = m.Pf_lag/m.P_lag; 
	
	if s.high_z_ForDebt == 0
		FXdebt = m.lambda + (1-m.lambda)*(m.xi_real/m.xi_real_lag);
		FXdebt_X = m.lambda_X + (1-m.lambda_X)*(m.xi_real/m.xi_real_lag);
    elseif s.high_z_ForDebt == 1
		FXdebt = m.lambda_mat + (1-m.lambda_mat).*(m.xi_real/m.xi_real_lag);
		FXdebt_X = m.lambda_X_mat + (1-m.lambda_X_mat).*(m.xi_real/m.xi_real_lag);	
	
	end
	
	
    %% Compute stationary measure across states

    [sim.measure,sim.diff_M] = KLS2_measure(s,r);  % Firms type 1
    [sim.measure_X,sim.diff_M_X] = KLS2_measure(s,r_X); % Firms type 2
   
    %Measure across all firms
    sim.measure_all = sim.measure*(1-m.Xshare) + m.Xshare *sim.measure_X ;

%% Exporters and non-exporters
        
    sim.share_x = sum(sum(sim.measure.*r.e))*(1-m.Xshare) + m.Xshare*sum(sum(sim.measure_X.*r_X.e)); 
    sim.share_x_composition =m.Xshare*sum(sum(sim.measure_X.*r_X.e))/sim.share_x; % Share of exporters type 2
    sim.share_nx = sum(sum(sim.measure.*(1-r.e)))*(1-m.Xshare) + sum(sum(sim.measure_X.*(1-r_X.e)))*m.Xshare;    
    sim.share_nx_composition = m.Xshare*sum(sum(sim.measure_X.*(1-r_X.e)))/sim.share_nx; % Share of non-exporters type 2
                
%% Productivity

    sim.z_avg = sum(sum(sim.measure,1).*r.z_grid'); %This should equal the average productivity computed using r.z_pi
    sim.z_avg_X = sum(sum(sim.measure_X,1).*r.z_grid');

%% Constrained Measures

    % Extensive margin:

    % unconstrained export policies:
    FC = m.F_base*(s.fcost_fgoods+(1-s.fcost_fgoods)*m.w);
    FC_X = m.F_X_base*(s.fcost_fgoods+(1-s.fcost_fgoods)*m.w);
    
    r.e_unc = r.pi_u_x-FC>= r.pi_u_nx;
    if s.model == 1 || s.model == 2 || s.model == 4
        r_X.e_unc = r_X.pi_u_x-FC_X>= r_X.pi_u_nx;
    elseif s.model == 3
        r_X.e_unc = r_X.pi_u_x>= r_X.pi_u_nx;
    end

    % Measure of firms constrained on the extensive margin:
    sim.ext_const = sum(sum(sim.measure.*(r.e_unc-r.e)));
    sim.ext_const_X = sum(sum(sim.measure_X.*(r_X.e_unc-r_X.e)));
    sim.ext_const_all = sim.ext_const*(1-m.Xshare) + sim.ext_const_X*m.Xshare;

% Intensive margin:

    % error margin (the same expressions can be numericallu slghtly different; matters for indicators; 'ind')
    err = 1e-10;

    % in terms of capital (indicator)
    sim.int_const_k_ind = sum(sum(sim.measure.*( r.e.*((r.k_u_x-err)>r.k) + (1-r.e).*((r.k_u_nx-err)>r.k) )));
    sim.int_const_k_ind_X = sum(sum(sim.measure_X.*( r_X.e.*((r_X.k_u_x-err)>r_X.k)+(1-r_X.e).*((r_X.k_u_nx-err)>r_X.k) )));
    sim.int_const_k_ind_all = sim.int_const_k_ind*(1-m.Xshare) + sim.int_const_k_ind_X*m.Xshare;

    % in terms of capital (percentages)
    sim.int_const_k_prc = sum(sum(sim.measure.*( r.e.*(r.k./r.k_u_x)+(1-r.e).*(r.k./r.k_u_nx) )));
    sim.int_const_k_prc_X = sum(sum(sim.measure_X.*( r_X.e.*(r_X.k./r_X.k_u_x)+(1-r_X.e).*(r_X.k./r_X.k_u_nx) )));
    sim.int_const_k_prc_all = sim.int_const_k_prc*(1-m.Xshare) +  sim.int_const_k_prc_X*m.Xshare;

    % in terms of profits (indicator)
    sim.int_const_pi_ind = sum(sum(sim.measure.*(r.e.*((r.pi_u_x-err)>r.pi)+(1-r.e).*((r.pi_u_nx-err)>r.pi))));
    sim.int_const_pi_ind_X = sum(sum(sim.measure_X.*( r_X.e.*((r_X.pi_u_x-err)>r_X.pi)+(1-r_X.e).*((r_X.pi_u_nx-err)>r_X.pi) )));
    sim.int_const_pi_ind_all = sim.int_const_pi_ind*(1-m.Xshare) + sim.int_const_pi_ind_X*m.Xshare;

    % in terms of capital (percentages)
    sim.int_const_pi_prc = sum(sum(sim.measure.*(r.e.*(r.pi./r.pi_u_x)+(1-r.e).*(r.pi./r.pi_u_nx))));
    sim.int_const_pi_prc_X = sum(sum(sim.measure_X.*( r_X.e.*(r_X.pi./r_X.pi_u_x)+(1-r_X.e).*(r_X.pi./r_X.pi_u_nx) )));
    sim.int_const_pi_prc_all = sim.int_const_pi_prc*(1-m.Xshare) + sim.int_const_pi_prc_X*m.Xshare;

% Constrained exporters:

    % in terms of capital (indicator)
    sim.int_const_k_ind_x = sum(sum( sim.measure.*( r.e.*((r.k_u_x-err)>r.k) ))) /  sum(sum( sim.measure.*r.e))  ;
    sim.int_const_k_ind_X_x = sum(sum( sim.measure_X.*( r_X.e.*((r_X.k_u_x-err)>r_X.k) ) )) /  sum(sum( sim.measure_X.*r_X.e));
    sim.int_const_k_ind_all_x = sim.int_const_k_ind_x*(1-m.Xshare) + sim.int_const_k_ind_X_x*m.Xshare;

    % in terms of capital (percentages)
    sim.int_const_k_prc_x = sum(sum(sim.measure.*( r.e.*(r.k./r.k_u_x) ))) /  sum(sum( sim.measure.*r.e)) ;
    sim.int_const_k_prc_X_x = sum(sum(sim.measure_X.*( r_X.e.*(r_X.k./r_X.k_u_x ) ))) /  sum(sum( sim.measure_X.*r_X.e)) ;
    sim.int_const_k_prc_all_x = sim.int_const_k_prc_x*(1-m.Xshare) +  sim.int_const_k_prc_X_x*m.Xshare;

     % in terms of profits (indicator)
    sim.int_const_pi_ind_x = sum(sum(sim.measure.*(r.e.*((r.pi_u_x-err)>r.pi) ))) / sum(sum(sim.measure.*(r.e)));
    sim.int_const_pi_ind_X_x = sum(sum(sim.measure_X.*( r_X.e.*((r_X.pi_u_x-err)>r_X.pi) ))) / sum(sum(sim.measure_X.*( r_X.e)));
    sim.int_const_pi_ind_all_x = sim.int_const_pi_ind_x*(1-m.Xshare) + sim.int_const_pi_ind_X_x*m.Xshare;

    % in terms of capital (percentages)
    sim.int_const_pi_prc_x = sum(sum(sim.measure.*(r.e.*(r.pi./r.pi_u_x)))) / sum(sum(sim.measure.*(r.e))) ;
    sim.int_const_pi_prc_X_x = sum(sum(sim.measure_X.*( r_X.e.*(r_X.pi./r_X.pi_u_x) )))  / sum(sum(sim.measure_X.*( r_X.e))) ;
    sim.int_const_pi_prc_all_x = sim.int_const_pi_prc_x*(1-m.Xshare) + sim.int_const_pi_prc_X_x*m.Xshare;   
        
% Constrained non-exporters:

    % in terms of capital (indicator)
    sim.int_const_k_ind_nx = sum(sum( sim.measure.*( (1-r.e).*((r.k_u_nx-err)>r.k) ))) /  sum(sum( sim.measure.*(1-r.e)))  ;
    sim.int_const_k_ind_X_nx = sum(sum( sim.measure_X.*( (1-r_X.e).*((r_X.k_u_nx-err)>r_X.k) ) )) /  sum(sum( sim.measure_X.*(1-r_X.e)));
    sim.int_const_k_ind_all_nx = sim.int_const_k_ind_nx*(1-m.Xshare) + sim.int_const_k_ind_X_nx*m.Xshare;

    % in terms of capital (percentages)
    sim.int_const_k_prc_nx = sum(sum(sim.measure.*( (1-r.e).*(r.k./r.k_u_nx) ))) /  sum(sum( sim.measure.*(1-r.e))) ;
    sim.int_const_k_prc_X_nx = sum(sum(sim.measure_X.*( (1-r_X.e).*(r_X.k./r_X.k_u_nx ) ))) /  sum(sum( sim.measure_X.*(1-r_X.e))) ;
    sim.int_const_k_prc_all_nx = sim.int_const_k_prc_nx*(1-m.Xshare) +  sim.int_const_k_prc_X_nx*m.Xshare;

     % in terms of profits (indicator)
    sim.int_const_pi_ind_nx = sum(sum(sim.measure.*((1-r.e).*((r.pi_u_nx-err)>r.pi) ))) / sum(sum(sim.measure.*(1-r.e)));
    sim.int_const_pi_ind_X_nx = sum(sum(sim.measure_X.*( (1-r_X.e).*((r_X.pi_u_nx-err)>r_X.pi) ))) / sum(sum(sim.measure_X.*(1- r_X.e)));
    sim.int_const_pi_ind_all_nx = sim.int_const_pi_ind_nx*(1-m.Xshare) + sim.int_const_pi_ind_X_nx*m.Xshare;

    % in terms of capital (percentages)
    sim.int_const_pi_prc_nx = sum(sum(sim.measure.*((1-r.e).*(r.pi./r.pi_u_nx)))) / sum(sum(sim.measure.*(1-r.e))) ;
    sim.int_const_pi_prc_X_nx = sum(sum(sim.measure_X.*( (1-r_X.e).*(r_X.pi./r_X.pi_u_nx) )))  / sum(sum(sim.measure_X.*(1- r_X.e))) ;
    sim.int_const_pi_prc_all_nx = sim.int_const_pi_prc_nx*(1-m.Xshare) + sim.int_const_pi_prc_X_nx*m.Xshare;   
           
    
    
%% Compute price and quantity indexes

%Domestic sales

    if s.model==1  || s.model == 2 || s.model == 4
        sim.Pd_sigma = sum( sim.measure(r.pd>0).*(r.pd(r.pd>0).^(1-m.sigma)) )*(1-m.Xshare)...
                       + sum( sim.measure_X(r_X.pd>0).*(r_X.pd(r_X.pd>0).^(1-m.sigma)) )*(m.Xshare);
        
        sim.Yd_sigma = sum( sim.measure(r.yd>0).*(r.yd(r.yd>0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare)...
                       +sum( sim.measure_X(r_X.yd>0).*(r_X.yd(r_X.yd>0).^((m.sigma-1)/m.sigma)) )*(m.Xshare);
                           
        sim.Pd = sim.Pd_sigma^(1/(1-m.sigma));
        sim.Yd = sim.Yd_sigma^(m.sigma/(m.sigma-1));     
        
        sim.PdYd =sum( sim.measure(r.pd>0).*(r.pd(r.pd>0).*r.yd(r.pd>0)))*(1-m.Xshare)...
                  + sum( sim.measure_X(r_X.pd>0).*(r_X.pd(r_X.pd>0).*r_X.yd(r_X.pd>0)))*(m.Xshare);
              
    elseif s.model==3
        
        sim.Pd_sigma = sum( sim.measure(r.pd>0).*(r.pd(r.pd>0).^(1-m.sigma)) )*(1-m.Xshare);
        sim.Yd_sigma = sum( sim.measure(r.yd>0).*(r.yd(r.yd>0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare);
        sim.Pd = ((sim.Pd_sigma )^(1/(1-m.sigma)));
        sim.Yd = ((sim.Yd_sigma )^(m.sigma/(m.sigma-1)));     
        sim.PdYd =sum( sim.measure(r.pd>0).*(r.pd(r.pd>0).*r.yd(r.pd>0)))*(1-m.Xshare);
    
    end

%Exports    
    sim.Px_sigma = sum( sim.measure(r.pf>0).*(r.pf(r.pf>0).^(1-m.sigma)) )*(1-m.Xshare)...
                 + sum( sim.measure_X(r_X.pf>0).*(r_X.pf(r_X.pf>0).^(1-m.sigma)) )*m.Xshare;
    
    Yx_sigma = sum( sim.measure(r.yf>0).*(r.yf(r.yf>0).^((m.sigma-1)/m.sigma)) )*(1-m.Xshare);
    Yx_sigma_X = sum( sim.measure_X(r_X.yf>0).*(r_X.yf(r_X.yf>0).^((m.sigma-1)/m.sigma)) )*m.Xshare;             
    sim.Yx_sigma = Yx_sigma_X + Yx_sigma;             
    Px_sigma_zero = sim.Px_sigma<=0;
    Yx_sigma_zero = sim.Yx_sigma<=0;
    sim.Px = m.xi_real*(1-Px_sigma_zero)*(((1-Px_sigma_zero)*sim.Px_sigma + Px_sigma_zero)^(1/(1-m.sigma)));
    sim.Yx_all = (1-Yx_sigma_zero)*(((1-Yx_sigma_zero)*sim.Yx_sigma + Yx_sigma_zero)^(m.sigma/(m.sigma-1)));                    
      
    Yx_sigma_zero = Yx_sigma<=0;
    Yx_sigma_zero_X = Yx_sigma_X<=0;
    sim.Yx = (1-Yx_sigma_zero)*(((1-Yx_sigma_zero)*Yx_sigma + Yx_sigma_zero)^(m.sigma/(m.sigma-1)));       
    sim.Yx_X = (1-Yx_sigma_zero_X)*(((1-Yx_sigma_zero_X)*Yx_sigma_X + Yx_sigma_zero_X)^(m.sigma/(m.sigma-1)));     
    
    sim.PxYx = m.xi_real*sum( sim.measure(r.pf>0).*(r.pf(r.pf>0).*r.yf(r.pf>0)))*(1-m.Xshare);    
    sim.PxYx_X = m.xi_real*sum( sim.measure_X(r_X.pf>0).*(r_X.pf(r_X.pf>0).*r_X.yf(r_X.pf>0)))*m.Xshare;
    sim.PxYx_all = sim.PxYx + sim.PxYx_X;
      
    sim.PxYx_all_USD = sim.PxYx_all/m.xi_real;   
    sim.PxYx_USD = sim.PxYx/m.xi_real;
    sim.PxYx_X_USD = sim.PxYx_X/m.xi_real;    

%Imports
    sim.Pm =  m.xi_real*m.Pm; %in units of domestic final goods.
    sim.Ym = m.omega_m^m.sigma*(m.A^(m.sigma-1))*m.Y*(sim.Pm^(-m.sigma));
    sim.Pm_sigma = (m.omega_m^m.sigma)*(sim.Pm^(1-m.sigma));
    sim.Ym_sigma = m.omega_m*sim.Ym^((m.sigma-1)/m.sigma);
    sim.PmYm=sim.Pm*sim.Ym;
       
%% Compute aggregate variables needed to evaluate market clearing conditions

%Labor and capital
    sim.K_all = sum(sum(sim.measure.*r.k))*(1-m.Xshare)+sum(sum(sim.measure_X.*r_X.k))*m.Xshare;
    sim.K = sum(sum(sim.measure.*r.k))*(1-m.Xshare);
    sim.K_X = sum(sum(sim.measure_X.*r_X.k))*m.Xshare;
    
    sim.inv_agg = m.delta*sim.K;
    sim.inv_agg_X = m.delta*sim.K_X;
    sim.inv_agg_all = sim.inv_agg + sim.inv_agg_X;
    
    sim.N = sum(sum(sim.measure.*r.n))*(1-m.Xshare);
    sim.N_X = sum(sum(sim.measure_X.*r_X.n))*m.Xshare;        
    
    if s.model == 1  || s.model == 2 || s.model == 4
        sim.FC =  m.F_base * sum(sum(sim.measure.*r.e))*(1-m.Xshare) + m.F_X_base * sum(sum(sim.measure_X.*r_X.e))*m.Xshare; % sim.FC in units according to s.fcost_fgoods
    elseif s.model == 3
        sim.FC = m.F_base * sum(sum(sim.measure.*r.e))*(1-m.Xshare) ;   
    end
        
    % Total labor
    sim.N_all = sim.N+sim.N_X+ sim.FC * (1-s.fcost_fgoods);
        
%Final good output and price
    %Tradable good 
    sim.Pcpi = (sim.Pd^(1-m.sigma) + max(sim.Pm_sigma,0) )^(1/(1-m.sigma)); % This should be equal to 1, Price of the domestic final good in units of the domestic final good
    sim.Ycpi = m.A*(sim.Yd^((m.sigma-1)/m.sigma) + max(sim.Ym_sigma,0) )^(m.sigma/(m.sigma-1));        

    %Final good
    sim.Y = sim.Ycpi; 
    sim.P = m.P;
    sim.xi_real = m.xi_real;
    
%Consumption    
    sim.C_all = sum(sum(sim.measure.*r.c))*(1-m.Xshare)+sum(sum(sim.measure_X.*r_X.c))*m.Xshare;                      
    sim.C = sum(sum(sim.measure.*r.c))*(1-m.Xshare);
    sim.C_X = sum(sum(sim.measure_X.*r_X.c))*m.Xshare;   
    
%% Market clearing conditions

%Labor
    sim.n_supply = 1;
    sim.n_demand = sim.N_all;
    sim.mc_n = 10*log((1+sim.n_demand)/(1+sim.n_supply));        
                                  
%Assets 
%This market clearing condition is the result of reformulating the debt market clearing condition
%In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below. 
    sim.a_supply = sum(sum(sim.measure.*r.a_grid_mat))*(1-m.Xshare) + sum(sum(sim.measure_X.*r.a_grid_mat))*m.Xshare; %sim.measure contains all types of agents: type 1 and 2 (only exporters)  
    sim.a_demand = sim.K_all;
    sim.mc_a = 10*log((1+sim.a_demand)/(1+sim.a_supply));          

%Final goods   
%Total production of final goods is used for: (i) consumption, (ii) investment, (iii) fixed costs
%Consumption and fixed costs expenditures are straightforward to compute
%Aggregate investment is derived as follows: nu*k_new + integral[k' - (1-delta)k] 
%Now, integral[(1-delta)k]=(1-delta)kagg and [nu*k_new + integral(k')]=kagg
%Then, aggregate investment equals delta*kagg
    sim.y_supply = sim.Y;    
    sim.y_demand = sim.C_all + m.delta*sim.K_all + sim.FC*s.fcost_fgoods;
    sim.mc_y = 10*log((1+sim.y_demand)/(1+sim.y_supply)); 
   
%Beliefs about tradable good
%To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
%In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
    sim.mc_p_belief = 10*log((1+sim.Pcpi)/(1+m.A)); %Now what we called sim.Pcpi should be equal to 1, see the notes. Once we divide by P, P enters only inside Xi=P*/P x Xi_nom.
                                                    %sim.Pcpi is not the price index but an object (function of prices) that has to equal 1 in equilibrium                                                       
    sim.mc_y_belief = 10*log((1+sim.Ycpi)/(1+m.Y));                                                       
                        
%% Display statistics    

%Setup
    %a = k - d/(1+r)
    r.d = (1+m.r)*(r.k - r.a_grid_mat); % This is the 'expected debt' at period n+1
    r_X.d = (1+m.r)*(r_X.k - r_X.a_grid_mat);  % This is the 'expected debt' at period n+1
    
        
    sim.credit_all = (sum(sum(sim.measure.*max(r.d,0)))*(1-m.Xshare) + m.Xshare*sum(sum(sim.measure_X.*max(r_X.d,0))))/(1+m.r); %Amount borrowed in current period (sum of positive d / 1+r), domestic units
    sim.credit_repayment_all = sum(sum(sim.measure.*max(r.d.*FXdebt,0)))*(1-m.Xshare) + m.Xshare*sum(sum(sim.measure_X.*max(r_X.d.*FXdebt_X,0))); %Credit repayment in current period (sum of positive d * FXdebt), domestic units  
               	   
			   
    sim.d_agg = (sum(sum(sim.measure.*r.d))*(1-m.Xshare) + m.Xshare*sum(sum(sim.measure_X.*r_X.d)))/(1+m.r); 
    sim.d_agg_X = m.Xshare*sum(sum(sim.measure_X.*r_X.d))/(1+m.r);
    sim.d_agg_repayment = (sum(sum(sim.measure.*r.d.*FXdebt))*(1-m.Xshare) + m.Xshare*sum(sum(sim.measure_X.*r_X.d.*FXdebt_X)));     
    
    
    sim.r_tilde = m.r;
    sim.r_tilde_X = m.r;
 
    sim.r_effective_avg = (1-m.Xshare)*(1+sim.r_tilde)*(1 + sum(sum(sum(sim.measure.*r.lm))))...
                            + m.Xshare*(1+sim.r_tilde_X)*(1 + sum(sum(sum(sim.measure_X.*r_X.lm))));     
     
    sim.GDP = sim.PdYd + sim.PxYx_all; %GDP denominated in units of final goods, sim.Px is already denominated in domestic currency
    sim.TFP = sim.GDP/((sim.K_all^m.alpha)*(sim.N_all^(1-m.alpha)));
    
    sim.NX_GDP = (sim.PxYx_all-sim.PmYm)/sim.GDP;
    sim.NFA_GDP = -sim.d_agg_repayment/sim.GDP; %Is this the right way to measure it? Or should we take into account changes in the RER? (as in d_agg_repayment)
    
    sales_x =  m.xi_real*r.pf.*r.yf; %r.pf is denominated in units of the foreign good, so we adjust it
    sales_d =r.pd.*r.yd;
    sales = sales_d+sales_x; 
    x_share = sales_x./sales;   
    
    sales_x_X =  m.xi_real*r_X.pf.*r_X.yf; %r.pf is denominated in foreign currency, so we adjust it
    sales_d_X =r_X.pd.*r_X.yd;
    sales_X = sales_d_X+sales_x_X; 
    x_share_X = sales_x_X./sales_X;   
    
    ln_sales = log(sales);
    ln_sales_X = log(sales_X);
    ln_sales_mean = sum(sum(sim.measure.*log(sales)))*(1-m.Xshare) + sum(sum(sim.measure_X.*log(sales_X)))*m.Xshare;
    sim.ln_sales_sd = sqrt(sum(sum(sim.measure .*(ln_sales - ln_sales_mean).^2))*(1-m.Xshare)+sum(sum(sim.measure_X.*(ln_sales_X - ln_sales_mean).^2))*m.Xshare);
    
    
    sim.sales_avg_nx = ((1-m.Xshare)*sum(sum(sim.measure.*sales.*(1-r.e))) + (m.Xshare)*sum(sum(sim.measure_X.*sales_X.*(1-r_X.e))))/sim.share_nx;  
    sim.labor_avg_nx = ((1-m.Xshare)*sum(sum(sim.measure.*r.n.*(1-r.e)))+ (m.Xshare)*sum(sum(sim.measure_X.*r_X.n.*(1-r_X.e))))/sim.share_nx;    
        
    sim.sales_avg_x = ((1-m.Xshare)*sum(sum(sim.measure.*sales.*r.e))+m.Xshare*sum(sum(sim.measure_X.*sales_X)))/ sim.share_x;                  
    sim.labor_avg_x = (sim.FC*(1-s.fcost_fgoods)+(1-m.Xshare)*sum(sum(sim.measure.*r.n.*r.e))+m.Xshare*sum(sum(sim.measure_X.*r_X.n)))/sim.share_x;      
    
    sim.sales_avg_type1 = sum(sum(sim.measure.*sales*(1-m.Xshare) + m.Xshare*sim.measure_X.*sales_X.*(1-r_X.e) )) / sum(sum((1-m.Xshare)*sim.measure + m.Xshare*sim.measure_X.*(1-r_X.e) ));  
    sim.sales_avg_type2 = sum(sum(sim.measure_X.*sales_X.*r_X.e*m.Xshare )) / sum(sum(sim.measure_X.*r_X.e*m.Xshare ));  
    
    sim.sales_avg_type1_2 = sum(sum(sim.measure.*sales));  
    sim.sales_avg_type2_2 = sum(sum(sim.measure_X.*sales_X));  
    
    sim.salesx_avg_type1 = sum(sum(sim.measure.*r.e.*sales_x)) / sum(sum(sim.measure.*r.e));  
    sim.salesx_avg_type2 = sum(sum(sim.measure_X.*sales_x_X.*r_X.e )) / sum(sum(sim.measure_X.*r_X.e ));  
    
    sim.salesd_avg_type1 = sum(sum(sim.measure.*sales_d*(1-m.Xshare) + m.Xshare*sim.measure_X.*sales_d_X.*(1-r_X.e) )) / sum(sum((1-m.Xshare)*sim.measure + m.Xshare*sim.measure_X.*(1-r_X.e) ));  
    sim.salesd_avg_type2 = sum(sum(sim.measure_X.*sales_d_X.*r_X.e*m.Xshare )) / sum(sum(sim.measure_X.*r_X.e*m.Xshare ));  
    
    sim.salesd_avg_type1_2 = sum(sum(sim.measure.*sales_d));  
    sim.salesd_avg_type2_2 = sum(sum(sim.measure_X.*sales_d_X));
    
%Size distribution
    sales_bothgroups = [sales(:);sales_X(:)];
    measure_bothgroups = [sim.measure(:)*(1-m.Xshare);sim.measure_X(:)*m.Xshare];
    
    [sales_bothgroups_sorted,sales_bothgroups_sorted_ind] = sort(sales_bothgroups);
    measure_bothgroups_sorted = measure_bothgroups(sales_bothgroups_sorted_ind);
    measure_bothgroups_sorted_cumsum = cumsum(measure_bothgroups_sorted);
    
    sales_p25_ind = find(measure_bothgroups_sorted_cumsum>=0.25,1,'first');
    sim.sales_p25 = sales_bothgroups_sorted(sales_p25_ind);    
    
    sales_p50_ind = find(measure_bothgroups_sorted_cumsum>=0.5,1,'first');
    sim.sales_p50 = sales_bothgroups_sorted(sales_p50_ind);    
        
    sales_p75_ind = find(measure_bothgroups_sorted_cumsum>=0.75,1,'first');
    sim.sales_p75 = sales_bothgroups_sorted(sales_p75_ind);
    
    sim.sales_share_top25 = sum(measure_bothgroups_sorted(sales_p75_ind:end).*sales_bothgroups_sorted(sales_p75_ind:end))/...
                        sum(measure_bothgroups_sorted.*sales_bothgroups_sorted);

%Size distribution for Exporters
    %sales_x_bothgroups = [sales_x(:);sales_x_X(:)];
    measure_x_bothgroups = [sim.measure(:).*r.e(:)*(1-m.Xshare);sim.measure_X(:).*r_X.e(:)*m.Xshare];
    
    %[sales_x_bothgroups_sorted,sales_x_bothgroups_sorted_ind] = sort(sales_x_bothgroups);
    measure_x_bothgroups_sorted = measure_x_bothgroups(sales_bothgroups_sorted_ind);
    measure_x_bothgroups_sorted_cumsum = cumsum(measure_x_bothgroups_sorted);
    
    x_sales_p25_ind = find(measure_x_bothgroups_sorted_cumsum>=0.25*sim.share_x,1,'first');
    sim.x_sales_p25 = sales_bothgroups_sorted(x_sales_p25_ind);    
    
    x_sales_p50_ind = find(measure_x_bothgroups_sorted_cumsum>=0.50*sim.share_x,1,'first');
    sim.x_sales_p50 = sales_bothgroups_sorted(x_sales_p50_ind);    
      
    x_sales_p75_ind = find(measure_x_bothgroups_sorted_cumsum>=0.75*sim.share_x,1,'first');
    sim.x_sales_p75 = sales_bothgroups_sorted(x_sales_p75_ind);   
                    
%Size distribution for Non-Exporters
    %sales_d_bothgroups = [sales_d(:);sales_d_X(:)];
    measure_d_bothgroups = [sim.measure(:).*(1-r.e(:))*(1-m.Xshare);sim.measure_X(:).*(1-r_X.e(:))*m.Xshare];
    
    %[sales_d_bothgroups_sorted,sales_d_bothgroups_sorted_ind] = sort(sales_d_bothgroups);
    measure_d_bothgroups_sorted = measure_d_bothgroups(sales_bothgroups_sorted_ind);
    measure_d_bothgroups_sorted_cumsum = cumsum(measure_d_bothgroups_sorted);
    
    nx_sales_p25_ind = find(measure_d_bothgroups_sorted_cumsum>=0.25*(1-sim.share_x),1,'first');
    sim.nx_sales_p25 = sales_bothgroups_sorted(nx_sales_p25_ind);    
    
    nx_sales_p50_ind = find(measure_d_bothgroups_sorted_cumsum>=0.50*(1-sim.share_x),1,'first');
    sim.nx_sales_p50 = sales_bothgroups_sorted(nx_sales_p50_ind);    
      
    nx_sales_p75_ind = find(measure_d_bothgroups_sorted_cumsum>=0.75*(1-sim.share_x),1,'first');
    sim.nx_sales_p75 = sales_bothgroups_sorted(nx_sales_p75_ind);    
    
                    
%Statistics to be displayed    
    sim.X_GDP = sim.PxYx_all/sim.GDP;  
    sim.X_X_GDP = sim.PxYx_X/sim.GDP; 
    
    sim.credit_gdp = sim.credit_repayment_all/sim.GDP;
    sim.xpremium_labor = sim.labor_avg_x/sim.labor_avg_nx;
    sim.xpremium_sales = sim.sales_avg_x/sim.sales_avg_nx;
    sim.xpremium_sales_med = sim.x_sales_p50/sim.nx_sales_p50;
    
    sim.type2premium_sales = sim.sales_avg_type2/sim.sales_avg_type1;
    sim.type2premium_sales_2 = sim.sales_avg_type2_2/sim.sales_avg_type1_2;
    sim.type2premium_salesx = sim.salesx_avg_type2/sim.salesx_avg_type1;
    sim.type2premium_salesd = sim.salesd_avg_type2/sim.salesd_avg_type1;
    sim.type2premium_salesd_2 = sim.salesd_avg_type2_2/sim.salesd_avg_type1_2;
    
    
    sim.x_share_avg =( (1-m.Xshare)* sum(sim.measure(x_share>0).*x_share(x_share>0))+m.Xshare* sum(sim.measure_X(x_share_X>0).*x_share_X(x_share_X>0))   ) / ((1-m.Xshare)*sum(sim.measure(x_share>0)) +m.Xshare*sum(sim.measure_X(x_share_X>0))  );
    sim.x_share_avg_type1 =  sum(sim.measure(x_share>0).*x_share(x_share>0))/sum(sim.measure(x_share>0));
    sim.x_share_avg_type2 =sum(sim.measure_X(x_share_X>0).*x_share_X(x_share_X>0))/sum(sim.measure_X(x_share_X>0));
    
    sim.X_D = sim.PxYx_all/sim.PdYd;  

    sim.k_wagebill = sim.a_demand/(m.w*sim.n_supply); %Includes fixed costs in case when they are denominated in labor units
    
%% Solution-related statistics

    sim.a_min_share = sum(sim.measure(1,:));
    sim.a_max_share = sum(sim.measure(end,:)); 
    
end

