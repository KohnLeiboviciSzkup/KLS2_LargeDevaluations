function [sim, r, r_X, s] = KLS2_simulate_shocks(m,s,r,r_X)

tic;

%% Generate random shocks and simulate policy functions

%Initalize objects
    ind_z = zeros(s.N,s.T*(1+s.burn));
    
%Initialize productivity realizations at stationary distribution
%Input -> z_pi, stationary distribution of productivity shocks, vector of size [s.z_grid_size x 1], computed by approx_shocks.m when run by model_solve.m
%Output -> ind_z(:,1), vector of indexes across productivity states for each firm, size [s.N x 1]
    
    k=1;
    for i=1:length(r.log_z_grid)
        v(k:k+floor(r.z_pi(i)*s.N)-1)=i;  
        k=length(v)+1;
    end
    v(length(v):s.N)=round(mean(v));
    ind_z(:,1)=v';
    clear v;

%Compute productivity realizations after initial state
%Input -> z.P, transition matrix computed by rowenhorst.m
%Output -> ind_z(:,2:end)
    C = cumsum(r.z_P,2);
    R = rand(s.N,s.T*(1+s.burn));
    for j=2:s.T*(1+s.burn);
        f = repmat(R(:,j),1,s.z_grid_size);
        ind_z(:,j) = 1+sum( f > C(ind_z(:,j-1),:) ,2);
        clear f;
    end
    clear R C;
    
    ind_z = ind_z';
    
    
%Alternative way of simulating shocks (based on compecon function ddpsimul.m)
    %ind_z = ddpsimul(r.z_P,v,s.T*(1+s.burn)-1)';

    
%%  Generating Policies
    
    e   =  r.e;
    e_X =  r_X.e;
    a   = r.a_grid_mat;
    a_X =  r.a_grid_mat; 
    ap_ind =r.ap_ind;
    ap_ind_X =r_X.ap_ind;
    yd   = r.yd;
    yd_X =  r_X.yd;
    yf   = r.yf;
    yf_X = r_X.yf;
    pd   = r.pd;
    pd_X =  r_X.pd;
    pf   =  r.pf;
    pf_X = r_X.pf;
    k   =   r.k;
    k_X =   r_X.k;
    n   =   r.n;
    n_X =   r_X.n;
    c   =   r.c;
    c_X =   r_X.c;
    lm  =   r.lm;
    lm_X  =   r_X.lm;

%%    
    
    %Initialize objects
    ind_a = zeros(s.T*(1+s.burn)+1,s.N);
    exporter = zeros(s.T*(1+s.burn)+1,s.N);   

    %First period initialization
    ind_a(1,1:s.N) = 1; %We initialize every firm with lowest assets
        
    N1=round(s.N*(1-m.Xshare)); % First N1 firms are type 1 (high trade costs), the remaining are type 2.
    N2=s.N-N1;
    
    %Compute asset and export states for t>1    
    for t=1:s.T*(1+s.burn)  

        
        %we map ind_a and ind_z to the policy functions using ind        
        ind = s.a_grid_size*(ind_z(t,:)-1)+ind_a(t,:);  

        
        %Update state variables    
        exporter(t+1,1:N1) = e(ind(1:N1)); % r.exporter is a decision rule given states today, sim.exporter updates simulation   
        exporter(t+1,N1+1:end) = e_X(ind(N1+1:end)); % r.exporter is a decision rule given states today, sim.exporter updates simulation   
       
        ind_a(t+1,1:N1)= ap_ind(ind(1:N1));
        ind_a(t+1,N1+1:end)= ap_ind_X(ind(N1+1:end));

    end
    
    ind_z_sim=ind_z;
    a_sim=r.a_grid(ind_a); %Asset values
    exporter_decision_sim = exporter(2:end,:); %Export decisions (as opposed to status)
    clear ind  ind_z
    % First N1 firms are type 1, firms N1+1 to N are type 2 (X)
    
    % Burn first s.burn*s.T time periods
    ind_a = ind_a(s.T*s.burn+1:end-1,:);
    ind_a_sim = ind_a;
    clear ind_a
    
    a_sim = a_sim(s.T*s.burn+1:end-1,:);
    a_sim_end = a_sim(end,:);
    ind_z_sim=ind_z_sim(s.T*s.burn+1:end,:);
    exporter_sim = exporter(s.T*s.burn+1:end-1,:);
    clear exporter
    
    sim.exporter_end = exporter_sim(end,:);
    exporter_decision_sim = exporter_decision_sim(s.T*s.burn+1:end,:);
    
    % saving for transition
    sim.exporter_decision_T = exporter_decision_sim(end-1,:)+1;

%% Compute Statistics    
    
     %Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
     m.xi_real = m.Pf/m.P;
     sim.xi_real = m.xi_real;
     m.xi_real_lag = m.Pf_lag/m.P_lag; 
     FXdebt = m.lambda + (1-m.lambda)*(m.xi_real/m.xi_real_lag);
     FXdebt_X = m.lambda_X + (1-m.lambda_X)*(m.xi_real/m.xi_real_lag);
    
    %we map ind_a and ind_z to the policy functions using ind        
    ind=s.a_grid_size*(ind_z_sim-1)+ind_a_sim;
    % sim.ind=ind;
     
    % Saving indices for transition:
    sim.ind_z_end = ind_z_sim(end,:); 
    sim.ind_a_end = ind_a_sim(end,:);     
    ind_end = ind(end,:); 
    clear ind_z_sim ind_a_sim ind
    
    
    
    
    %Productivity     
% HERE %     sim.z = r.z_grid(ind_z);
     
    %Quantities        
%     yf_sim = zeros(s.T,s.N);
%     yf_sim(:,1:N1) = yf(ind(:,1:N1)).*exporter_decision_sim(:,1:N1);
%     yf_sim(:,N1+1:end) = yf_X(ind(:,N1+1:end)).*exporter_decision_sim(:,N1+1:end);
%     yf_sim_end = yf_sim(end,:);
%     clear yf_sim
    yf_sim_end = zeros(1,s.N);
    yf_sim_end(1,1:N1) = yf(ind_end(1:N1)).*exporter_decision_sim(end,1:N1);
    yf_sim_end(1,N1+1:end) = yf_X(ind_end(N1+1:end)).*exporter_decision_sim(end,N1+1:end);


    
%     yd_sim = zeros(s.T,s.N);
%     yd_sim(:,1:N1) = yd(ind(:,1:N1));
%     yd_sim(:,N1+1:end) = yd_X(ind(:,N1+1:end));
%     yd_end=yd_sim(end,:); 
%     clear yd_sim
    yd_end = zeros(1,s.N);
    yd_end(:,1:N1) = yd(ind_end(1:N1));
    yd_end(:,N1+1:end) = yd_X(ind_end(N1+1:end));

    
%     pf_sim = zeros(s.T,s.N);
%     pf_sim(:,1:N1) = pf(ind(:,1:N1)).*exporter_decision_sim(:,1:N1);
%     pf_sim(:,N1+1:end) = pf_X(ind(:,N1+1:end)).*exporter_decision_sim(:,N1+1:end);
%     pf_sim_end = pf_sim(end,:);
%     clear pf_sim
    pf_sim_end = zeros(1,s.N);
    pf_sim_end(1,1:N1) = pf(ind_end(1:N1)).*exporter_decision_sim(end,1:N1);
    pf_sim_end(1,N1+1:end) = pf_X(ind_end(N1+1:end)).*exporter_decision_sim(end,N1+1:end);

    
%     pd_sim = zeros(s.T,s.N);
%     pd_sim(:,1:N1) = pd(ind(:,1:N1));
%     pd_sim(:,N1+1:end) = pd_X(ind(:,N1+1:end));
%     pd_end=pd_sim(end,:);
%     clear pd_sim
    pd_end = zeros(1,s.N);
    pd_end(1,1:N1) = pd(ind_end(1:N1));
    pd_end(1,N1+1:end) = pd_X(ind_end(N1+1:end));
 

%     k_sim = zeros(s.T,s.N);
%     k_sim(:,1:N1) = k(ind(:,1:N1));
%     k_sim(:,N1+1:end) = k_X(ind(:,N1+1:end));
%     k_sim_end = k_sim(end,:);
%     clear k_sim
    k_sim_end = zeros(1,s.N);
    k_sim_end(1,1:N1) = k(ind_end(1:N1));
    k_sim_end(1,N1+1:end) = k_X(ind_end(N1+1:end));
    
    
%     n_sim = zeros(s.T,s.N);
%     n_sim(:,1:N1) = n(ind(:,1:N1));
%     n_sim(:,N1+1:end) = n_X(ind(:,N1+1:end));
%     n_sim_end = n_sim(end,:);
%     clear n_sim
    n_sim_end = zeros(1,s.N);
    n_sim_end(1,1:N1) = n(ind_end(1:N1));
    n_sim_end(1,N1+1:end) = n_X(ind_end(N1+1:end));

    
%     c_sim = zeros(s.T,s.N);
%     c_sim(:,1:N1) = c(ind(:,1:N1));
%     c_sim(:,N1+1:end) = c_X(ind(:,N1+1:end));
%     c_sim_end = c_sim(end,:);
%     clear c_sim
    c_sim_end = zeros(1,s.N);
    c_sim_end(1,1:N1) = c(ind_end(1:N1));
    c_sim_end(1,N1+1:end) = c_X(ind_end(N1+1:end));
 
    
%     d_sim = zeros(s.T,s.N);
%     d = (1+m.r)*(k - a);        %a = k - d/(1+r) This is the 'expected debt' at period n+1
%     d_X = (1+m.r)*(k_X - a_X);  %a = k - d/(1+r) This is the 'expected debt' at period n+1
%     d_sim(:,1:N1) = d(ind(:,1:N1));
%     d_sim(:,N1+1:end) = d_X(ind(:,N1+1:end));
%     d_end = d_sim(end,:);
%     clear d_sim
    d_end = zeros(1,s.N);
    d = (1+m.r)*(k - a);        %a = k - d/(1+r) This is the 'expected debt' at period n+1
    d_X = (1+m.r)*(k_X - a_X);  %a = k - d/(1+r) This is the 'expected debt' at period n+1
    d_end(1,1:N1) = d(ind_end(1:N1));
    d_end(1,N1+1:end) = d_X(ind_end(N1+1:end));


%     lm_sim = zeros(s.T,s.N);
%     lm_sim(1,1:N1) = lm(ind(:,1:N1));
%     lm_sim(1,N1+1:end) = lm_X(ind(:,N1+1:end));
%     sim.lm_end=lm_sim(end,:);
%     clear lm_sim    
    sim.lm_end = zeros(1,s.N);
    sim.lm_end(1,1:N1) = lm(ind_end(1:N1));
    sim.lm_end(1,N1+1:end) = lm_X(ind_end(N1+1:end));


 %% Market Clearing conditions for the LAST simulated period


 % Domestic sales
     sim.Pd_sigma = sum(pd_end(pd_end>0).^(1-m.sigma))/s.N; 
     sim.Yd_sigma =  sum(yd_end(yd_end>0).^((m.sigma-1)/m.sigma))/s.N; 
     sim.Pd = sim.Pd_sigma^(1/(1-m.sigma));
     sim.Yd = sim.Yd_sigma^(m.sigma/(m.sigma-1));
     sim.PdYd = sum(  pd_end.*yd_end)./s.N;
         
 
 % Exports    
 
     sim.Px_sigma = sum(pf_sim_end(pf_sim_end>0).^(1-m.sigma))/s.N;     
     sim.Yx_sigma = sum(yf_sim_end(yf_sim_end>0).^((m.sigma-1)/m.sigma))/s.N;
     sim.Px = m.xi_real* sim.Px_sigma^(1/(1-m.sigma));
     sim.Yx_all = m.xi_real *sim.Yx_sigma^(m.sigma/(m.sigma-1));
     sim.PxYx_all = m.xi_real*sum(pf_sim_end.*yf_sim_end)./s.N;
     sim.PxYx_all_USD = sim.PxYx_all/m.xi_real;      
     
     sim.PxYx = m.xi_real*sum( pf_sim_end(1:N1).*yf_sim_end(1:N1) )./s.N;  % Divided by the TOTAL number firms
     sim.PxYx_X = m.xi_real*sum(  pf_sim_end(N1+1:end).*yf_sim_end(N1+1:end) )./s.N; % Divided by the TOTAL number firms
     
 % Imports
    sim.Pm =  m.xi_real*m.Pm; %in units of domestic final goods.
    sim.Ym = m.omega_m^m.sigma*(m.A^(m.sigma-1))*m.Y*(sim.Pm^(-m.sigma));
    sim.Pm_sigma = (m.omega_m^m.sigma)*(sim.Pm^(1-m.sigma));
    sim.Ym_sigma = m.omega_m*sim.Ym^((m.sigma-1)/m.sigma);
    sim.PmYm=sim.Pm*sim.Ym;
     
 % Compute aggregate variables needed to evaluate market clearing conditions
 
 %Capital
     sim.K_all = sum(k_sim_end)/s.N;
    
 %Fixed costs   
     sim.FC = (1/s.N)*( sum( exporter_decision_sim(end,1:N1).* m.F_base ) +  sum(exporter_decision_sim(end,N1+1:end).* m.F_X_base ) );
    
 
 % Labor
     sim.N_all=sum(n_sim_end)/s.N+(1-s.fcost_fgoods)*sim.FC;
          
 % Final good output and price
     %Tradable good 
     sim.Pcpi = (sim.Pd^(1-m.sigma) + sim.Pm_sigma )^(1/(1-m.sigma));
     sim.Ycpi = m.A*(sim.Yd^((m.sigma-1)/m.sigma) + sim.Ym_sigma)^(m.sigma/(m.sigma-1));        
 
     %Final good
     sim.Y = sim.Ycpi; 
     sim.P = m.P;
 
 % Consumption
     
     sim.C = sum(c_sim_end)/s.N;
      
 %% Market clearing conditions
 
 %Labor
     sim.n_supply = 1;
     sim.n_demand = sim.N_all; %Total labor demand = Demand from entrepreneurs      + fixed export costs
     sim.mc_n = 10*log((1+sim.n_demand)/(1+sim.n_supply));        
                                   
 %Assets 
 %This market clearing condition is the result of reformulating the debt market clearing condition
 %In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below. 
     sim.a_supply = sum(a_sim_end)/s.N;
     sim.a_demand = sim.K_all;   % we already divided K by N
     sim.mc_a = 10*log((1+sim.a_demand)/(1+sim.a_supply));          
 
 %Final goods   
 %Total production of final goods is used for: (i) consumption, (ii) investment, (iii) fixed costs
 %Consumption and fixed costs expenditures are straightforward to compute
 %Aggregate investment is derived as follows: nu*k_new + integral[k' - (1-delta)k] 
 %Now, integral[(1-delta)k]=(1-delta)kagg and [nu*k_new + integral(k')]=kagg
 %Then, aggregate investment equals delta*kagg
     sim.y_supply = sim.Y;
     sim.y_demand = sim.C + m.delta*sim.K_all + sim.FC*s.fcost_fgoods ; 
     sim.mc_y = 10*log((1+sim.y_demand)/(1+sim.y_supply)); 
     
 %Beliefs about tradable good
 %To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
 %In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
     sim.mc_p_belief = 10*log((1+sim.Pcpi)/(1+m.A)); %Now what we called sim.Pcpi should be equal to 1, see the notes. Once we divide by P, P enters only inside Xi=P*/P x Xi_nom.
                                                     %sim.Pcpi is not the price index but an object (function of prices) that has to equal 1 in equilibrium                                                       
     sim.mc_y_belief = 10*log((1+sim.Ycpi)/(1+m.Y));                                                       
           

     toc;

%% Display statistics      

    sim.credit_all = (sum(max(d_end,0))/(1+m.r))/s.N; %Amount borrowed in current period (sum of positive d / 1+r), domestic units
    sim.credit_repayment_all = (sum( max(d_end(1:N1),0))*FXdebt +sum( max(d_end(N1+1:end),0))*FXdebt_X   )/s.N;  %Credit repayment in current period (sum of positive d * FXdebt), domestic units  
    
    sim.d_agg =  (sum(d_end)/(1+m.r))/s.N;   %d_agg_all
    sim.d_agg_X =   (sum(d_end(N1+1:end))/(1+m.r))/N2;  
    sim.d_agg_repayment = (sum( d_end(1:N1))*FXdebt +sum( d_end(N1+1:end))*FXdebt_X )/s.N; %d_agg_repayment_all

    sim.r_tilde = m.r;
    sim.r_tilde_X = m.r;
        
    sim.r_effective_avg = ( (1+sim.r_tilde)*(1 + sum(sim.lm_end(1:N1) )) + (1+sim.r_tilde_X)*(1 + sum(sim.lm_end(N1+1:end))) )/s.N ; 

    
    %Sales     
    sim.sales_x_end = m.xi_real*pf_sim_end.*yf_sim_end; %pf_sim is denominated in foreign currency, so we adjust it
    sim.sales_d_end = pd_end.*yd_end;
    sim.sales_end = sim.sales_d_end+sim.sales_x_end;
    x_share_end = sim.sales_x_end./sim.sales_end; %Export intensity  
    
    
    ln_sales = log(sim.sales_end);
    ln_sales_mean = sum(log(sim.sales_end))/s.N;
    sim.ln_sales_sd = sqrt((1/s.N)*sum((ln_sales(end,:) - ln_sales_mean).^2));
    
    
    
    sim.GDP =  sum(sim.sales_end)./s.N;  
    sim.TFP = sim.GDP/((sim.K_all^m.alpha)*(sim.N_all^(1-m.alpha)));
    
    sim.NX_GDP = (sim.PxYx-sim.PmYm)/sim.GDP;
    sim.NFA_GDP = -sim.d_agg_repayment/sim.GDP;
    
        
    sim.sales_avg_nx = sum(sim.sales_end.*(1-exporter_decision_sim(end,:)))/sum(1-exporter_decision_sim(end,:));
    sim.labor_avg_nx = sum(n_sim_end.*(1-exporter_decision_sim(end,:)))/sum(1-exporter_decision_sim(end,:));
           
    sim.sales_avg_x = sum(sim.sales_end.*exporter_decision_sim(end,:))/sum(exporter_decision_sim(end,:));
    sim.labor_avg_x = (s.N*sim.FC*(1-s.fcost_fgoods) +sum(n_sim_end.*exporter_decision_sim(end,:))  ) /sum(exporter_decision_sim(end,:));
    
    sim.sales_avg_type1  = (sum(sim.sales_end(1:N1))+sum(sim.sales_end(N1+1:end).*(exporter_decision_sim(end,N1+1:end)==0) )  ) /  (N1+sum(exporter_decision_sim(end,N1+1:end)==0));
    sim.sales_avg_type2  = sum(sim.sales_end(N1+1:end).*(exporter_decision_sim(end,N1+1:end)==1) )   /  sum(exporter_decision_sim(end,N1+1:end)==1);      
     
    sim.sales_avg_type1_2 =sum(sim.sales_end(1:N1))/N1; 
    sim.sales_avg_type2_2 =sum(sim.sales_end(N1+1:end))/N2; 
    
    sim.salesx_avg_type1 = sum(sim.sales_x_end(1:N1))/ sum(exporter_decision_sim(end,1:N1)==1) ; 
    sim.salesx_avg_type2 = sum(sim.sales_x_end(N1+1:end))/ sum(exporter_decision_sim(end,N1+1:end)==1) ; 
    
    sim.salesd_avg_type1 = (sum(sim.sales_d_end(1:N1)) + sum(sim.sales_d_end(N1+1:end).*(exporter_decision_sim(end,N1+1:end)==0) ) ) /  (N1+sum(exporter_decision_sim(end,N1+1:end)==0)); 
    sim.salesd_avg_type2 = sum(sim.sales_d_end(N1+1:end).*(exporter_decision_sim(end,N1+1:end)==1) )   /  sum(exporter_decision_sim(end,N1+1:end)==1);
    
    sim.salesd_avg_type1_2 =sum(sim.sales_d_end(1:N1))/N1; 
    sim.salesd_avg_type2_2 =sum(sim.sales_d_end(N1+1:end))/N2; 
    
    
    % Firm-level statistics
    sim.share_exporters = sum(exporter_decision_sim(end,:))/s.N; % share of exporters
    sim.share_x=sim.share_exporters;
    sim.share_nonexporters = 1-sim.share_exporters;
    sim.share_x_composition = sum(exporter_decision_sim(end,N1+1:end))/sum(exporter_decision_sim(end,:)); % Share of exporters type 2
    sim.share_nx_composition =sum((1-exporter_decision_sim(end,N1+1:end)))/sim.share_nonexporters; % Share of non-exporters type 2
        
    sim.share_exporters_vec = sum(exporter_decision_sim,2)/s.N;
    
    %Export entry and exit rates
    sim.share_starters =   sum(exporter_decision_sim.*(1-exporter_sim),2) ./ sum(1-exporter_sim,2);   % Vector of share of starters by period
    sim.share_stoppers =   sum((1-exporter_decision_sim).*exporter_sim,2) ./ sum(exporter_sim,2) ;   % Vector of share of stoppers by period
    clear exporter_sim
    
    %Export intensity
    temp = x_share_end(exporter_decision_sim(end,:)==1);
    sim.export_sales_ratio_med =median(temp(:));        
    clear temp     
    
    % saving for transition:
     sim.exporter_decision_T = exporter_decision_sim(end-1,:);
     y1 = prctile(sim.sales_end(exporter_decision_sim(end,:)==1),[25 50 75]);
     y0 = prctile(sim.sales_end(exporter_decision_sim(end,:)==0),[25 50 75]); 
     exporter_decision_sim_end = exporter_decision_sim(end,:);
     clear exporter_decision_sim
     
    %Size distribution
    y_all = prctile(sim.sales_end,[25 50 75]); 
    sim.sales_p25=y_all(1);
    sim.sales_p50=y_all(2);
    sim.sales_p75=y_all(3);
    sim.sales_share_top25 = sum(sim.sales_end(sim.sales_end>sim.sales_p75))/sum(sim.sales_end);

    
    %Size distribution for Exporters
    sim.x_sales_p25=y1(1);
    sim.x_sales_p50=y1(2);
    sim.x_sales_p75=y1(3);
    
    %Size distribution for Non-Exporters
    sim.nx_sales_p25=y0(1);
    sim.nx_sales_p50=y0(2);
    sim.nx_sales_p75=y0(3);
    
      
    %Statistics to be displayed    
    sim.X_GDP = sim.PxYx/sim.GDP;  
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
     
    sim.x_share_avg = sum(x_share_end)/  sum(x_share_end>0);  %conditional on exporting
    sim.x_share_avg_type1 =  sum(x_share_end(1:N1))/  sum(x_share_end(1:N1)>0); 
    sim.x_share_avg_type2 =  sum(x_share_end(N1+1:end))/ sum(x_share_end(N1+1:end)>0);

    sim.X_D = sim.PxYx/sim.PdYd ;  

    sim.k_wagebill = sim.a_demand/(m.w*sim.n_supply);  %without fixed costs
    
    
    
%% Constrained (by productivity quartile):

    pi = r.pi;
    pi_X = r_X.pi;
    pi_sim_end = zeros(1,s.N);
    pi_sim_end(1,1:N1) = pi(ind_end(1:N1));
    pi_sim_end(1,N1+1:end) = pi_X(ind_end(N1+1:end));  

    % unconstrained export policies:
    e_unc = r.pi_u_x-m.F_base>= r.pi_u_nx;
    e_unc_X = r_X.pi_u_x-m.F_X_base>= r_X.pi_u_nx;

    % optimal choice of exporting
    exporter_decision_unc = zeros(1,s.N);
    exporter_decision_unc(1,1:N1) = e_unc(ind_end(1:N1));
    exporter_decision_unc(1,N1+1:end) = e_unc_X(ind_end(N1+1:end));
    
    % percentiles
%     indz_p99 = prctile(sim.ind_z_end(1,:),99);
    indz_p75 = prctile(sim.ind_z_end(1,:),75);
    indz_p50 = prctile(sim.ind_z_end(1,:),50);
    indz_p25 = prctile(sim.ind_z_end(1,:),25);
    
    % categories of productivity
    z_cat{1,1} = sim.ind_z_end(1,:)<=indz_p25;
    z_cat{2,1} = sim.ind_z_end(1,:)>indz_p25 & sim.ind_z_end(1,:)<=indz_p50;
    z_cat{3,1} = sim.ind_z_end(1,:)>indz_p50 & sim.ind_z_end(1,:)<=indz_p75;
    z_cat{4,1} = sim.ind_z_end(1,:)>indz_p75;
    N_z_cat{1} = sum(z_cat{1,1}==1);
    N_z_cat{2} = sum(z_cat{2,1}==1);
    N_z_cat{3} = sum(z_cat{3,1}==1);
    N_z_cat{4} = sum(z_cat{4,1}==1);
    
    % unconstrained capital:
    k_unc   = r.k_u_nx.*(1-e) + r.k_u_x.*e;
    k_unc_X = r_X.k_u_nx.*(1-e_X) + r_X.k_u_x.*e_X;
    k_unc_sim_end = zeros(1,s.N);
    k_unc_sim_end(1,1:N1) = k_unc(ind_end(1:N1));
    k_unc_sim_end(1,N1+1:end) = k_unc_X(ind_end(N1+1:end));
 
    % unconstrained profits:
    pi_unc   = r.pi_u_nx.*(1-e)+r.pi_u_x.*e;
    pi_unc_X = r_X.pi_u_nx.*(1-e_X)+r_X.pi_u_x.*e_X;
    pi_unc_sim_end = zeros(1,s.N);
    pi_unc_sim_end(1,1:N1) = pi_unc(ind_end(1:N1));
    pi_unc_sim_end(1,N1+1:end) = pi_unc_X(ind_end(N1+1:end));
    
    % aggregate:
    sim.ext_k_const = sum(k_sim_end./k_unc_sim_end)/s.N;
    sim.ext_pi_const = sum(pi_sim_end./pi_unc_sim_end)/s.N;
    sim.ext_k_const_ind = sum(k_sim_end<k_unc_sim_end)/s.N;
    sim.ext_pi_const_ind = sum(pi_sim_end<pi_unc_sim_end)/s.N;
    sim.ext_const_total = sum(exporter_decision_unc - exporter_decision_sim_end)/s.N;
    
    % by productivity:
    for i=1:numel(z_cat)
        
        % constrained on the extensive margin
        sim.ext_const(i,1) = sum(exporter_decision_unc(z_cat{i}==1) - exporter_decision_sim_end(z_cat{i}==1))/N_z_cat{i};
        
        % constrained on the intensive margin (profits)
        sim.int_const_pi_prc(i,1) = sum(pi_sim_end(z_cat{i}==1)./pi_unc_sim_end(z_cat{i}==1))/N_z_cat{i};
        sim.int_const_pi_ind(i,1) = sum(pi_sim_end(z_cat{i}==1)<pi_unc_sim_end(z_cat{i}==1))/N_z_cat{i};
        
        % constrained on the inetnsive margin
        sim.int_const_k_prc(i,1) = sum(k_sim_end(z_cat{i}==1)./k_unc_sim_end(z_cat{i}==1))/N_z_cat{i}; 
        sim.int_const_k_ind(i,1) = sum(k_sim_end(z_cat{i}==1)<k_unc_sim_end(z_cat{i}==1))/N_z_cat{i}; 
        
    end
    
end

