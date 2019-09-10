
%% This code saves output from KLS2_main

    sim_fun.w=Prices(2,:); % Wages
    
    % Share of exporters and non-exporters
    share_nx_all = sim_fun.share_nx;
    share_nx_all(1) = sim_0.share_nx;
    share_x_all = sim_fun.share_x;
    share_x_all(1) = sim_0.share_x;
    
    % Export elasticity
    XElasticity_log_value = (log(sim_fun.PxYx_all_USD) - log(sim_fun.PxYx_all_USD(1)))./(log(sim_fun.xi_real) - log(sim_fun.xi_real(1)));        
    XElasticity_log_value(1) = 0;
    
    XElasticity_log = (log(sim_fun.X_Laspeyres) - log(sim_fun.X_Laspeyres(1)))./(log(sim_fun.xi_real) - log(sim_fun.xi_real(1)));        
    XElasticity_log(1) = 0;
    
    XElasticity = ((sim_fun.X_Laspeyres - sim_fun.X_Laspeyres(1))/sim_fun.X_Laspeyres(1))./((sim_fun.xi_real - sim_fun.xi_real(1))/sim_fun.xi_real(1));        
    XElasticity(1) = 0;

% Moments:

%     Share of exporters with X/Y>0.6
%     Avg. export intensity, X/Y<0.6
%     Avg. export intensity, X/Y>0.6
%     Share of exporters
%     Std. dev. Log(sales)
%     Share of sales by largest 25% firms
%     NX/GDP
%     Credit/GDP

    moments_0 = [sim_0.share_x_composition sim_0.x_share_avg_type1 sim_0.x_share_avg_type2 sim_0.share_x sim_0.ln_sales_sd sim_0.sales_share_top25 sim_0.NX_GDP sim_0.credit_gdp];
    moments_end = [sim_end.share_x_composition sim_end.x_share_avg_type1 sim_end.x_share_avg_type2 sim_end.share_x sim_end.ln_sales_sd sim_end.sales_share_top25 sim_end.NX_GDP sim_end.credit_gdp];
   
    

    sim_fun.inv_agg_all(1)=m.delta*sim_0.K_all;
        
    Results= [sim_fun.share_x, sim_fun.share_x_composition, sim_fun.share_nx, sim_fun.share_nx_composition, sim_fun.ext_const, sim_fun.ext_const_X, sim_fun.ext_const_all, sim_fun.int_const_k_ind, sim_fun.int_const_k_ind_X, sim_fun.int_const_k_ind_all, sim_fun.int_const_k_prc, sim_fun.int_const_k_prc_X, sim_fun.int_const_k_prc_all, sim_fun.int_const_pi_ind, sim_fun.int_const_pi_ind_X, sim_fun.int_const_pi_ind_all, sim_fun.int_const_pi_prc, sim_fun.int_const_pi_prc_X, sim_fun.int_const_pi_prc_all, sim_fun.Pd, sim_fun.Yd, sim_fun.PdYd, sim_fun.Yx_sigma, sim_fun.Px, sim_fun.Yx, sim_fun.Yx_X, sim_fun.PxYx, sim_fun.PxYx_X, sim_fun.PxYx_all, sim_fun.PxYx_all_USD, sim_fun.PxYx_USD, sim_fun.PxYx_X_USD, sim_fun.Pm, sim_fun.Ym, sim_fun.PmYm, sim_fun.K_all, sim_fun.K, sim_fun.K_X, sim_fun.inv_agg, sim_fun.inv_agg_X, sim_fun.inv_agg_all, sim_fun.N_X, sim_fun.N, sim_fun.FC, sim_fun.N_all, sim_fun.Pcpi, sim_fun.Ycpi, sim_fun.Y, sim_fun.P, sim_fun.xi_real, sim_fun.C_all, sim_fun.C, sim_fun.C_X, sim_fun.n_supply, sim_fun.n_demand, sim_fun.mc_n, sim_fun.a_supply, sim_fun.a_demand, sim_fun.mc_a, sim_fun.y_supply, sim_fun.y_demand, sim_fun.mc_y, sim_fun.mc_p_belief, sim_fun.mc_y_belief, sim_fun.credit_all, sim_fun.credit_repayment_all, sim_fun.d_agg, sim_fun.d_agg_X, sim_fun.d_agg_repayment, sim_fun.r_tilde, sim_fun.r_tilde_X, sim_fun.r_effective_avg, sim_fun.GDP, sim_fun.TFP, sim_fun.NX_GDP, sim_fun.NFA_GDP, sim_fun.sales_avg_nx, sim_fun.sales_avg_x, sim_fun.X_GDP, sim_fun.credit_gdp, sim_fun.xpremium_sales, sim_fun.X_D, sim_fun.a_min_share, sim_fun.a_max_share, sim_fun.Yx_sumall, sim_fun.credit_X, sim_fun.credit_repayment_X, sim_fun.d_agg_X_repayment, sim_fun.exp_avg_x, sim_fun.GDP_Laspeyres, sim_fun.X_Laspeyres, sim_fun.D_Laspeyres, sim_fun.X_Laspeyres_Type1, sim_fun.X_Laspeyres_Type2, sim_fun.w' m.pm_v m.zagg_v m.rv]';

	% Series in Results:
%	1    sim_fun.share_x
%	2	 sim_fun.share_x_composition
%	3	 sim_fun.share_nx
%	4	 sim_fun.share_nx_composition
%	5	 sim_fun.ext_const
%	6	 sim_fun.ext_const_X
%	7	 sim_fun.ext_const_all
%	8	 sim_fun.int_const_k_ind
%	9	 sim_fun.int_const_k_ind_X
%	10	 sim_fun.int_const_k_ind_all
%	11	 sim_fun.int_const_k_prc
%	12	 sim_fun.int_const_k_prc_X
%	13	 sim_fun.int_const_k_prc_all
%	14	 sim_fun.int_const_pi_ind
%	15	 sim_fun.int_const_pi_ind_X
%	16	 sim_fun.int_const_pi_ind_all
%	17	 sim_fun.int_const_pi_prc
%	18	 sim_fun.int_const_pi_prc_X
%	19	 sim_fun.int_const_pi_prc_all
%	20	 sim_fun.Pd
%	21	 sim_fun.Yd
%	22	 sim_fun.PdYd
%	23	 sim_fun.Yx_sigma
%	24	 sim_fun.Px
%	25	 sim_fun.Yx
%	26	 sim_fun.Yx_X
%	27	 sim_fun.PxYx
%	28	 sim_fun.PxYx_X
%	29	 sim_fun.PxYx_all
%	30	 sim_fun.PxYx_all_USD
%	31	 sim_fun.PxYx_USD
%	32	 sim_fun.PxYx_X_USD
%	33	 sim_fun.Pm
%	34	 sim_fun.Ym
%	35	 sim_fun.PmYm
%	36	 sim_fun.K_all
%	37	 sim_fun.K
%	38	 sim_fun.K_X
%	39	 sim_fun.inv_agg
%	40	 sim_fun.inv_agg_X
%	41	 sim_fun.inv_agg_all
%	42	 sim_fun.N_X
%	43	 sim_fun.N
%	44	 sim_fun.FC
%	45	 sim_fun.N_all
%	46	 sim_fun.Pcpi
%	47	 sim_fun.Ycpi
%	48	 sim_fun.Y
%	49	 sim_fun.P
%	50	 sim_fun.xi_real
%	51	 sim_fun.C_all
%	52	 sim_fun.C
%	53	 sim_fun.C_X
%	54	 sim_fun.n_supply
%	55	 sim_fun.n_demand
%	56	 sim_fun.mc_n
%	57	 sim_fun.a_supply
%	58	 sim_fun.a_demand
%	59	 sim_fun.mc_a
%	60	 sim_fun.y_supply
%	61	 sim_fun.y_demand
%	62	 sim_fun.mc_y
%	63	 sim_fun.mc_p_belief
%	64	 sim_fun.mc_y_belief
%	65	 sim_fun.credit_all
%	66	 sim_fun.credit_repayment_all
%	67	 sim_fun.d_agg
%	68	 sim_fun.d_agg_X
%	69	 sim_fun.d_agg_repayment
%	70	 sim_fun.r_tilde
%	71	 sim_fun.r_tilde_X
%	72	 sim_fun.r_effective_avg
%	73	 sim_fun.GDP
%	74	 sim_fun.TFP
%	75	 sim_fun.NX_GDP
%	76	 sim_fun.NFA_GDP
%	77	 sim_fun.sales_avg_nx
%	78	 sim_fun.sales_avg_x
%	79	 sim_fun.X_GDP
%	80	 sim_fun.credit_gdp
%	81	 sim_fun.xpremium_sales
%	82	 sim_fun.X_D
%	83	 sim_fun.a_min_share
%	84	 sim_fun.a_max_share
%	85	 sim_fun.Yx_sumall
%	86	 sim_fun.credit_X
%	87	 sim_fun.credit_repayment_X
%	88	 sim_fun.d_agg_X_repayment
%	89	 sim_fun.exp_avg_x
%	90	 sim_fun.GDP_Laspeyres
%	91	 sim_fun.X_Laspeyres
%	92	 sim_fun.D_Laspeyres
%	93	 sim_fun.X_Laspeyres_Type1
%	94	 sim_fun.X_Laspeyres_Type2
%	95	 sim_fun.w
%	96	 m.pm_v
%	97	 m.zagg_v
%	98	 m.rv

   
    periods = 10;

    series{1} = log(sim_fun.xi_real(1:periods))-log(sim_fun.xi_real(1));
    series{2} = XElasticity_log(1:periods);
    series{3} = log(sim_fun.w(1:periods)')-log(sim_fun.w(1));
    series{4} = log(sim_fun.GDP_Laspeyres(1:periods))-log(sim_fun.GDP_Laspeyres(1));
    series{5} = log(sim_fun.PmYm(1:periods))-log(sim_fun.PmYm(1));
    series{6} = log(sim_fun.X_Laspeyres(1:periods))-log(sim_fun.X_Laspeyres(1));
    series{7} = sim_fun.NX_GDP(1:periods); 
    series{8} = log(sim_fun.C_all(1:periods))-log(sim_fun.C_all(1));
    series{9} = [sim_0.K_all*m.delta;sim_fun.inv_agg_all(2:periods)]-sim_0.K_all*m.delta;
    series{10} = sim_fun.share_x(1:periods);
    series{11} = log([sim_0.credit_repayment_all;sim_fun.credit_repayment_all(2:periods)])-log(sim_0.credit_repayment_all);
    series{12} = log([sim_0.credit_repayment_all/sim_fun.xi_real(1);sim_fun.credit_repayment_all(2:periods)./sim_fun.xi_real(2:periods)])-log(sim_0.credit_repayment_all/sim_fun.xi_real(1));
    series{13} = sim_fun.r_effective_avg(1:periods);
    series{14} = XElasticity_log_value(1:periods);
        
    output = [];
    for i=1:numel(series)
        output = [output;series{i}];
    end  
    

    
    