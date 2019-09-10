function [sim, rt, rt_X] = KLS2_simulate_shocks_trans(m,s,r,rt,rt_X,trans,Yguess,Pguess,sim_0)

%% Compute stationary measure across states

    T_trans = trans.N;
    
    sim_0.yf=[];
    sim_0.yd=[];
    sim_0.pf=[];
    sim_0.pd=[];
    sim_0.k=[];
    sim_0.n=[];
    sim_0.c=[];
    sim_0.d=[];
    sim_0.sales_d=[];
    sim_0.sales_f=[];
    sim_0.sales=[];
    
    sim_end.yf=[];
    sim_end.yd=[];
    sim_end.pf=[];
    sim_end.pd=[];
    sim_end.k=[];
    sim_end.n=[];
    sim_end.c=[];
    sim_end.d=[];
    sim_end.sales_d=[];
    sim_end.sales_f=[];
    sim_end.sales=[];
    

    %sim = sim_0;
    
    % There are two type of firms: firms that can export and produce
    % domestically (type 1) and firms that only export (type 2)
    N1=round(s.N*(1-m.Xshare)); % First N1 firms are type 1, the remaining are type 2.
    N2=s.N-N1;

    
%% Generate random shocks and simulate policy functions

%Initalize objects
    ind_z = zeros(s.N,T_trans);
    ind_a = zeros(s.N,T_trans);

    ind_z(:,1) = sim_0.ind_z_end;    % productivity state in the last period of SS
    ind_a(:,1) = sim_0.ind_a_end;     % median(sim_0.ind_a(end-10:end,:),1)';    % assets state in the last period of SS
    %end-10:

    ind_k(:,1) = sim_0.exporter_decision_T+1; %ones(s.N,1); %Status does not matter
      
 
%Compute productivity realizations after initial state
%Input -> z.P, transition matrix computed by rowenhorst.m
%Output -> ind_z(:,2:end)
    C = cumsum(r.z_P,2);
    R = rand(s.N,T_trans);
    for j=2:T_trans;
        f = repmat(R(:,j),1,s.z_grid_size);
        ind_z(:,j) = 1+sum( f > C(ind_z(:,j-1),:) ,2);
        clear f;
    end
    clear R C;
    
   ind_z = ind_z';
   ind_a = ind_a';    % assets state in the last period of SS
   ind_k = ind_k';  
  
    ind_trans=zeros(T_trans,s.N);
    k_state=zeros(T_trans,s.N);
    
%%  Generating Policies
    
   
    % states of trans.et: assets,productivity, time (a,z,t)
     et = trans.et; 
     apt = trans.apt;
     apt_ind = trans.apt_ind;
     kt = trans.kt;
     et_X = trans.et_X;
     apt_X = trans.apt_X;
     apt_ind_X  = trans.apt_ind_X;
     kt_X = trans.kt_X;

     % initialization:
    exporter_decision = zeros(T_trans,s.N);

   
     %Compute asset and export states for t>1    
    for t=1:T_trans  
        
        ap_ind = squeeze(apt_ind(:,:,t));  
        e  = squeeze(et(:,:,t));    
        k  = squeeze(kt(:,:,t));  
        ap_ind_X = squeeze(apt_ind_X(:,:,t));  
        e_X  = squeeze(et_X(:,:,t)); 
        k_X  = squeeze(kt_X(:,:,t));  
        
        %Construct 3-dimensional indexes from unidimensional indexes ind_a, ind_z, and ind_k     
        %The policy functions are multidimensional matrices -> we map ind_a, ind_z, and ind_k to the policy functions using ind        
        ind = s.a_grid_size*(ind_z(t,:)-1)+ind_a(t,:);  
        ind=ind';
        
        ind_trans(t,:)=ind;
     
        %Update state variables             
        ind_a(t+1,1:N1)= ap_ind(ind(1:N1));
        ind_a(t+1,N1+1:end)= ap_ind_X(ind(N1+1:end));  
        
        k_state(t,1:N1) = k(ind(1:N1));
        k_state(t,N1+1:end) = k_X(ind(N1+1:end));
        
        exporter_decision(t,1:N1) = e(ind(1:N1));
        exporter_decision(t,N1+1:end) = e_X(ind(N1+1:end));         
        
        z(t,:) = rt{1}.z_grid_mat(ind);
        a(t,:) = rt{1}.a_grid_mat(ind);
                
    end    
    
    %Initialize objects
    yf_trans  = zeros(T_trans,s.N);
    yd_trans  = zeros(T_trans,s.N);
    pf_trans  = zeros(T_trans,s.N);
    pd_trans  = zeros(T_trans,s.N);
    pf1_trans = zeros(T_trans,s.N);
    pd1_trans = zeros(T_trans,s.N);    
    k_trans   = zeros(T_trans,s.N);
    n_trans   = zeros(T_trans,s.N);
    c_trans   = zeros(T_trans,s.N);
    debt_trans = zeros(T_trans,s.N);
    lm_trans = zeros(T_trans,s.N);

    % prices in the initial steady state:
    pd1 = rt{1}.pd;
    pf1 = rt{1}.pf;
    pd1_X = rt_X{1}.pd;
    pf1_X = rt_X{1}.pf;
        
    
%% Here starts the loop to compute market clearing conditions and other
%  statistics
    

% initializations:
sales_x  = zeros(T_trans,s.N);
sales =  zeros(T_trans,s.N);
sales_x_Laspeyres  = zeros(T_trans,s.N);
sales_d_Laspeyres  = zeros(T_trans,s.N);
sales_Laspeyres =  zeros(T_trans,s.N);

for t=1:T_trans   %t=2:T_trans
    
%Real exchange rate:
    sim.xi_real(t) = m.Pfv(t)/Pguess(t); %Pf=1 for all n
    
    if t==1
        FXdebt = 1;
        FXdebt_X = 1;
    elseif t==2
        FXdebt = m.lambda_period_2 + (1-m.lambda_period_2)*(sim.xi_real(t)/sim.xi_real(t-1)); 
        FXdebt_X = m.lambda_period_2_X+(1-m.lambda_period_2_X)*(sim.xi_real(t)/sim.xi_real(t-1));        
    else
        FXdebt = m.lambda + (1-m.lambda)*(sim.xi_real(t)/sim.xi_real(t-1)); 
        FXdebt_X = m.lambda_X+(1-m.lambda_X)*(sim.xi_real(t)/sim.xi_real(t-1));
    end

%% Investment

    if t>1 && t<T_trans  
        inv(t,:) = k_state(t+1,:) - (1-m.delta)*k_state(t,:);
        sim.inv_agg(t,1) = sum(k_state(t+1,:) - (1-m.delta)*k_state(t,:) )/s.N;         
    else       
       inv(t,:) = m.delta*k_state(t,:);
       sim.inv_agg(t,1) = sum(m.delta*k_state(t,:))/s.N;       
    end
    

%% Exporters and non-exporters

%%  Generating Policies

    e = rt{t}.e;
    ap_ind = rt{t}.ap_ind;
    yd = rt{t}.yd;
    yf= rt{t}.yf;
    pd= rt{t}.pd;
    pf = rt{t}.pf;
    k  = rt{t}.k;
    n  = rt{t}.n;
    c = rt{t}.c; 
    debt = rt{t}.debt;
    lm = rt{t}.lm;
    
    e_X = rt_X{t}.e;
    ap_ind_X = rt_X{t}.ap_ind;
    yd_X = rt_X{t}.yd;
    yf_X = rt_X{t}.yf;
    pd_X = rt_X{t}.pd;
    pf_X = rt_X{t}.pf;
    k_X  = rt_X{t}.k;
    n_X = rt_X{t}.n;
    c_X  = rt_X{t}.c;
    debt_X = rt_X{t}.debt; 
    lm_X = rt_X{t}.lm;
    
    ind=ind_trans(t,:);
            
    %Quantities            
    yf_trans(t,1:N1) = yf(ind(1:N1)).*exporter_decision(t,1:N1);
    yf_trans(t,N1+1:end) = yf_X(ind(N1+1:end)).*exporter_decision(t,N1+1:end);
    yd_trans(t,1:N1) = yd(ind(1:N1));
    yd_trans(t,N1+1:end) = yd_X(ind(N1+1:end));
       
    % Prices:
    pf_trans(t,1:N1) = pf(ind(1:N1)).*exporter_decision(t,1:N1);
    pf_trans(t,N1+1:end) = pf_X(ind(N1+1:end)).*exporter_decision(t,N1+1:end);
    pd_trans(t,1:N1) = pd(ind(1:N1));
    pd_trans(t,N1+1:end) = pd_X(ind(N1+1:end));    

    % Keeping prices constant:
    pf1_trans(t,1:N1) = pf1(ind(1:N1)).*exporter_decision(t,1:N1);
    pf1_trans(t,N1+1:end) = pf1_X(ind(N1+1:end)).*exporter_decision(t,N1+1:end);
    pd1_trans(t,1:N1) = pd1(ind(1:N1));
    pd1_trans(t,N1+1:end) = pd1_X(ind(N1+1:end));      
    
    % Capital:
    k_trans(t,1:N1) = k(ind(1:N1));
    k_trans(t,N1+1:end) = k_X(ind(N1+1:end));
    
    % Labor:
    n_trans(t,1:N1) = n(ind(1:N1));
    n_trans(t,N1+1:end) = n_X(ind(N1+1:end));    
    
    % Consumption:
    c_trans(t,1:N1) = c(ind(1:N1));
    c_trans(t,N1+1:end) = c_X(ind(N1+1:end));       
    
    % Foreign Sales 
    sales_x(t,:) = sim.xi_real(t).*pf_trans(t,:).*yf_trans(t,:);
    sim.sales_x_total(t,1) = sum(sales_x(t,:))/s.N;
    
    % Foreign Sales Laspeyers:
    sales_x_Laspeyres(t,:) = sim.xi_real(1).*pf1_trans(t,:).*yf_trans(t,:);
    sim.sales_x_total_Laspeyres(t,1) = sum(sales_x_Laspeyres(t,:))/s.N;    
    
    % Total Sales:
    sales(t,:) = pd_trans(t,:).*yd_trans(t,:) + sim.xi_real(t).*pf_trans(t,:).*yf_trans(t,:);
    sim.sales_total(t,1) = sum(sales(t,:))/s.N; 
        
    % Total Sales Laspeyers:
    sales_Laspeyres(t,:) = pd1_trans(t,:).*yd_trans(t,:) + sim.xi_real(1).*pf1_trans(t,:).*yf_trans(t,:);
    sim.sales_total_Laspeyres(t,1) = sum(sales_Laspeyres(t,:))/s.N;     
    
    %Domestic Sales Laspeyres
    sales_d_Laspeyres(t,:) = pd1_trans(t,:).*yd_trans(t,:);
    sim.sales_d_total_Laspeyres(t,1) = sum(sales_d_Laspeyres(t,:))/s.N;    
    
    % Debt
    debt_trans(t,1:N1) = debt(ind(1:N1));
    debt_trans(t,N1+1:end) = debt_X(ind(N1+1:end));
    
    % Lagrange Multipliers:
    lm_trans(t,1:N1) = lm(ind(1:N1));
    lm_trans(t,N1+1:end) = lm_X(ind(N1+1:end));
    
    sim.share_x(t,1) = sum(exporter_decision(t,:))/s.N;
    sim.share_d(t,1) = sum(1-exporter_decision(t,:))/s.N;   
    
    if t>1
        sim.entry(t,1) = sum(exporter_decision(t,:).*(1-exporter_decision(t-1,:)))/sum((1-exporter_decision(t-1,:)));
        sim.exit(t,1)  = sum((1-exporter_decision(t,:)).*exporter_decision(t-1,:))/sum(exporter_decision(t-1,:));   
    end

%% Compute aggregate variables needed to evaluate market clearing conditions

%Domestic sales
    Pd_sigma = sum(pd_trans(t,:).^(1-m.sigma))/s.N;
    Yd_sigma = sum(yd_trans(t,:).^((m.sigma-1)/m.sigma))/s.N;  

    sim.Pd(t,1) = Pd_sigma^(1/(1-m.sigma));
    sim.Yd(t,1) = Yd_sigma^(m.sigma/(m.sigma-1));     

%Imports for a small open economy    
    sim.Pm(t,1) =  sim.xi_real(t)*m.pm_v(t); %in units of domestic final goods.
    sim.Ym(t,1) = m.omega_m^m.sigma*(m.A_v(t)^(m.sigma-1))*Yguess(t)*(sim.Pm(t)^(-m.sigma));
    sim.Pm_sigma = (m.omega_m^m.sigma)*(sim.Pm(t)^(1-m.sigma));
    sim.Ym_sigma = m.omega_m*sim.Ym(t,1)^((m.sigma-1)/m.sigma);

%Labor and capital    
    sim.K_all(t,1) = sum(k_trans(t,:))/s.N;
    
%Fixed costs   
    if s.model==1  || s.model == 2  || s.model == 4
        sim.FC(t,1) = (sum(exporter_decision(t,1:N1))/s.N).*m.F_base + (sum(exporter_decision(t,N1+1:end))/s.N).*m.F_X_base;
    elseif s.model==3
        sim.FC(t,1) = (sum(exporter_decision(t,1:N1))/s.N).*m.F_base;
    end

    % Total labor
    sim.N_all(t,1) = sum(n_trans(t,:))/s.N + sim.FC(t,1)*(1-s.fcost_fgoods); 
        
%Consumption    
    sim.C_all(t,1) = sum(c_trans(t,:))/s.N;    
  
%Final good
    sim.Pcpi(t,1) = (sim.Pd(t,1)^(1-m.sigma) + max(sim.Pm_sigma,0))^(1/(1-m.sigma));
    sim.Ycpi(t,1) = m.A_v(t)*(sim.Yd(t,1)^((m.sigma-1)/m.sigma) + max(sim.Ym_sigma,0))^(m.sigma/(m.sigma-1));            
    
    sim.Y(t,1) = sim.Ycpi(t,1); 
    sim.P(t,1) = Pguess(t);       
    
%% Market clearing conditions    
    
%Labor
    sim.n_supply(t,1) = 1;
    sim.n_demand(t,1) = sim.N_all(t,1); %Total labor demand = Demand from entrepreneurs  + fixed export costs + non-tradable goods     
    sim.mc_n(t,1) = 10*log((1+sim.n_demand(t,1))/(1+sim.n_supply(t,1)));        
    
%Final goods   
%Total production of final goods is used for: (i) consumption, (ii) investment, (iii) fixed costs
%Consumption and fixed costs expenditures are straightforward to compute
%Aggregate investment is derived as follows: nu*k_new + integral[k' - (1-delta)k] 
%Now, integral[(1-delta)k]=(1-delta)kagg and [nu*k_new + integral(k')]=kagg
%Then, aggregate investment equals delta*kagg
    sim.y_supply(t,1) = sim.Y(t,1);
    sim.y_demand(t,1) = sim.C_all(t,1)+ sim.inv_agg(t,1)+ sim.FC(t,1)*s.fcost_fgoods;      
    sim.mc_y(t,1) = 10*log((1+sim.y_demand(t,1))/(1+sim.y_supply(t,1)));    
    
%Beliefs about tradable good
%To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
%In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
    sim.mc_p_belief(t,1) = 10*log((1+sim.Pcpi(t,1))/(1+m.A_v(t)));        
    sim.mc_y_belief(t,1) = 10*log((1+sim.Y(t,1))/(1+Yguess(t)));           
  
end    
    
%% Decomposition

for t=1:T_trans

    sim.Xtotal(t,1) = sum(sales_x_Laspeyres(t,:))/s.N;
    sim.Dtotal(t,1) = sum(sales_d_Laspeyres(t,:))/s.N;
    
    if t>1
        new_exporters = exporter_decision(t,:)==1 & exporter_decision(t-1,:)==0;
        sim.Xnew(t,1) = sum(sales_x_Laspeyres(t,new_exporters))/s.N;
        
        cont_back = exporter_decision(t,:)==1 & exporter_decision(t-1,:)==1;       
        sim.Xcont_back(t,1) = sum(sales_x_Laspeyres(t,cont_back))/s.N;
        sim.Dcont_back(t,1) = sum(sales_d_Laspeyres(t,cont_back))/s.N;
        sim.Ycont_back(t,1) = sum(sales_Laspeyres(t,cont_back))/s.N;
    end
    
    if t<T_trans
       exit_exporters = exporter_decision(t,:)==1 & exporter_decision(t+1,:)==0; 
       sim.Xexit(t,1) = sum(sales_x_Laspeyres(t,exit_exporters))/s.N;
       
       cont_fwd = exporter_decision(t,:)==1 & exporter_decision(t+1,:)==1;
       sim.Xcont_fwd(t,1) = sum(sales_x_Laspeyres(t,cont_fwd))/s.N;
       sim.Dcont_fwd(t,1) = sum(sales_d_Laspeyres(t,cont_fwd))/s.N;
       sim.Ycont_fwd(t,1) = sum(sales_Laspeyres(t,cont_fwd))/s.N;
    end

end

for t=2:T_trans
    sim.growth_Xtotal(t,1) = (sim.Xtotal(t,1)-sim.Xtotal(t-1,1))/sim.Xtotal(t-1,1);
    sim.growth_Xextensive(t,1) = (sim.Xnew(t,1)-sim.Xexit(t-1,1))/sim.Xtotal(t-1,1);
    sim.growth_Xintensive(t,1) = (sim.Xcont_back(t,1)-sim.Xcont_fwd(t-1,1))/sim.Xtotal(t-1,1);
    sim.growth_Xintensive_Y(t,1) = (sim.Ycont_back(t,1)-sim.Ycont_fwd(t-1,1))/sim.Xtotal(t-1,1);
    sim.growth_Xintensive_D(t,1) = -(sim.Dcont_back(t,1)-sim.Dcont_fwd(t-1,1))/sim.Xtotal(t-1,1);
end


%% Decomposition - Relative to initial SS

for t=1:T_trans

    sim.Ytotal2(t,1) = sum(sales_Laspeyres(t,:))/s.N;
    sim.Xtotal2(t,1) = sum(sales_x_Laspeyres(t,:))/s.N;
    sim.Dtotal2(t,1) = sum(sales_d_Laspeyres(t,:))/s.N;
    
    if t>1
        new_exporters = exporter_decision(t,:)==1 & exporter_decision(1,:)==0;
        sim.Xnew2(t,1) = sum(sales_x_Laspeyres(t,new_exporters))/s.N;
        
        cont_back = exporter_decision(t,:)==1 & exporter_decision(1,:)==1;       
        sim.Xcont_back2(t,1) = sum(sales_x_Laspeyres(t,cont_back))/s.N;
        sim.Dcont_back2(t,1) = sum(sales_d_Laspeyres(t,cont_back))/s.N;
        sim.Ycont_back2(t,1) = sum(sales_Laspeyres(t,cont_back))/s.N;
    end
    
   exit_exporters = exporter_decision(1,:)==1 & exporter_decision(t,:)==0; 
   sim.Xexit2(t,1) = sum(sales_x_Laspeyres(1,exit_exporters))/s.N;

   cont_fwd = exporter_decision(1,:)==1 & exporter_decision(t,:)==1;
   sim.Xcont_fwd2(t,1) = sum(sales_x_Laspeyres(1,cont_fwd))/s.N;
   sim.Dcont_fwd2(t,1) = sum(sales_d_Laspeyres(1,cont_fwd))/s.N;
   sim.Ycont_fwd2(t,1) = sum(sales_Laspeyres(1,cont_fwd))/s.N;

end

for t=2:T_trans
    sim.growth_Xtotal2(t,1) = (sim.Xtotal2(t,1)-sim.Xtotal2(1,1))/sim.Xtotal2(1,1);
    sim.growth_Xextensive2(t,1) = (sim.Xnew2(t,1)-sim.Xexit2(t,1))/sim.Xtotal2(1,1);
    sim.growth_Xintensive2(t,1) = (sim.Xcont_back2(t,1)-sim.Xcont_fwd2(t,1))/sim.Xtotal2(1,1);
    sim.growth_Xintensive_Y2(t,1) = (sim.Ycont_back2(t,1)-sim.Ycont_fwd2(t,1))/sim.Xtotal2(1,1);
    sim.growth_Xintensive_D2(t,1) = -(sim.Dcont_back2(t,1)-sim.Dcont_fwd2(t,1))/sim.Xtotal2(1,1);
    
    sim.growth_Xintensive_Y3(t,1) = (sim.Ycont_back2(t,1)-sim.Ycont_fwd2(t,1))/sim.Ytotal2(1,1);
    sim.growth_Xintensive_D3(t,1) = -(sim.Dcont_back2(t,1)-sim.Dcont_fwd2(t,1))/sim.Dtotal2(1,1);
    sim.growth_Xintensive_Y3_wt(t,1) = sim.Ytotal2(1,1)/sim.Xtotal2(1,1);
    sim.growth_Xintensive_D3_wt(t,1) = sim.Dtotal2(1,1)/sim.Xtotal2(1,1);
end


%% Regressions


dln_sales_x_rel1994(1,:) = log(sales_x_Laspeyres(2,:)) - log(sales_x_Laspeyres(1,:));
dln_sales_x_rel1994(2,:) = log(sales_x_Laspeyres(3,:)) - log(sales_x_Laspeyres(1,:));
dln_sales_x_rel1994(3,:) = log(sales_x_Laspeyres(4,:)) - log(sales_x_Laspeyres(1,:));
dln_sales_x_rel1994(4,:) = log(sales_x_Laspeyres(5,:)) - log(sales_x_Laspeyres(1,:));
dln_sales_x_rel1994(5,:) = log(sales_x_Laspeyres(6,:)) - log(sales_x_Laspeyres(1,:));

%Low X/Y
    dummy_1995_lowXY = zeros(5,s.N);
    dummy_1995_lowXY(1,1:N1) = 1;

    dummy_1996_lowXY = zeros(5,s.N);
    dummy_1996_lowXY(2,1:N1) = 1;

    dummy_1997_lowXY = zeros(5,s.N);
    dummy_1997_lowXY(3,1:N1) = 1;
    
    dummy_1998_lowXY = zeros(5,s.N);
    dummy_1998_lowXY(4,1:N1) = 1;
    
    dummy_1999_lowXY = zeros(5,s.N);
    dummy_1999_lowXY(5,1:N1) = 1;    

%High X/Y
    dummy_1995_highXY = zeros(5,s.N);
    dummy_1995_highXY(1,N1+1:end) = 1;

    dummy_1996_highXY = zeros(5,s.N);
    dummy_1996_highXY(2,N1+1:end) = 1;

    dummy_1997_highXY = zeros(5,s.N);
    dummy_1997_highXY(3,N1+1:end) = 1;
    
    dummy_1998_highXY = zeros(5,s.N);
    dummy_1998_highXY(4,N1+1:end) = 1;
    
    dummy_1999_highXY = zeros(5,s.N);
    dummy_1999_highXY(5,N1+1:end) = 1;       
    
%Constant
    constant = ones(5,s.N);
    
%Valid observations
    valid = isfinite(dln_sales_x_rel1994);
    
    sim.dln_sales_x_rel1994 = dln_sales_x_rel1994(valid);
    sim.dummy_1995_lowXY = dummy_1995_lowXY(valid);
    sim.dummy_1996_lowXY = dummy_1996_lowXY(valid);
    sim.dummy_1997_lowXY = dummy_1997_lowXY(valid);
    sim.dummy_1998_lowXY = dummy_1998_lowXY(valid);
    sim.dummy_1999_lowXY = dummy_1999_lowXY(valid);
    sim.dummy_1995_highXY = dummy_1995_highXY(valid);
    sim.dummy_1996_highXY = dummy_1996_highXY(valid);
    sim.dummy_1997_highXY = dummy_1997_highXY(valid);
    sim.dummy_1998_highXY = dummy_1998_highXY(valid);
    sim.dummy_1999_highXY = dummy_1999_highXY(valid);
    sim.constant = constant(valid);
    
    sim.regcoeffs = regress(sim.dln_sales_x_rel1994,[sim.constant...
        sim.dummy_1996_lowXY sim.dummy_1997_lowXY sim.dummy_1998_lowXY sim.dummy_1999_lowXY...
        sim.dummy_1995_highXY sim.dummy_1996_highXY sim.dummy_1997_highXY sim.dummy_1998_highXY sim.dummy_1999_highXY]);
    
    sim.dln_sales_x_rel1994_lowXY = sim.regcoeffs(1) + [0 sim.regcoeffs(2) sim.regcoeffs(3) sim.regcoeffs(4) sim.regcoeffs(5)];
    sim.dln_sales_x_rel1994_highXY = sim.regcoeffs(1) + [sim.regcoeffs(6) sim.regcoeffs(7) sim.regcoeffs(8) sim.regcoeffs(9) sim.regcoeffs(10)];
    
%Aggregate exports by export intensity category
for t=1:6
    sim.agg_salesx_lowXY(t,1) = sum(sales_x_Laspeyres(t,1:N1))/s.N;
    sim.agg_salesx_highXY(t,1) = sum(sales_x_Laspeyres(t,N1+1:end))/s.N;
end
sim.agg_salesx_lowXY = sim.agg_salesx_lowXY/sim.agg_salesx_lowXY(1,1);
sim.agg_salesx_highXY = sim.agg_salesx_highXY/sim.agg_salesx_highXY(1,1);

%% Exports growth by debt

%Debt vs. savings
    firms_debt1994 = debt_trans(1,:)>0;
    firms_savings1994 = debt_trans(1,:)<=0;
    sim.N_debt = sum(debt_trans(1,:)>0);
    sim.N_savings = sum(debt_trans(1,:)<=0);
    
%Debt: binding vs. not binding
    firms_debt1994binding = debt_trans(1,:)>0 & lm_trans(1,:)>0;
    firms_debt1994notbinding = debt_trans(1,:)>0 & lm_trans(1,:)==0;

%Debt vs. savings, by export status  
    firms_debt1994_x = debt_trans(1,:)>0 & exporter_decision(1,:)==1;
    firms_savings1994_x = debt_trans(1,:)<=0 & exporter_decision(1,:)==1;
    firms_debt1994_nx = debt_trans(1,:)>0 & exporter_decision(1,:)==0;
    firms_savings1994_nx = debt_trans(1,:)<=0 & exporter_decision(1,:)==0;
    sim.N_debt_x = sum(debt_trans(1,:)>0 & exporter_decision(1,:)==1);
    sim.N_savings_x = sum(debt_trans(1,:)<=0 & exporter_decision(1,:)==1);
    sim.N_debt_nx = sum(debt_trans(1,:)>0 & exporter_decision(1,:)==01);
    sim.N_savings_nx = sum(debt_trans(1,:)<=0 & exporter_decision(1,:)==0);    
    
%Debt vs. savings, by export intensity categories     
    firms_debt1994_lowXY = firms_debt1994_x;
    firms_debt1994_lowXY(1,N1+1:end) = 0;
    firms_savings1994_lowXY = firms_savings1994_x;
    firms_savings1994_lowXY(1,N1+1:end) = 0;    
    sim.N_debt_low_XY = sum(firms_debt1994_lowXY);
    sim.N_savings_low_XY = sum(firms_savings1994_lowXY);
    
    firms_debt1994_highXY = firms_debt1994_x;
    firms_debt1994_highXY(1,1:N1) = 0;
    firms_savings1994_highXY = firms_savings1994_x;
    firms_savings1994_highXY(1,1:N1) = 0;   
    sim.N_debt_high_XY = sum(firms_debt1994_highXY);
    sim.N_savings_high_XY = sum(firms_savings1994_highXY);  
    
%Regression
    dln_sales_rel1994(1,:) = log(sales_Laspeyres(2,:)) - log(sales_Laspeyres(1,:));
    dln_sales_rel1994(2,:) = log(sales_Laspeyres(3,:)) - log(sales_Laspeyres(1,:));
    dln_sales_rel1994(3,:) = log(sales_Laspeyres(4,:)) - log(sales_Laspeyres(1,:));
    dln_sales_rel1994(4,:) = log(sales_Laspeyres(5,:)) - log(sales_Laspeyres(1,:));
    dln_sales_rel1994(5,:) = log(sales_Laspeyres(6,:)) - log(sales_Laspeyres(1,:));
    
    dln_sales_x_rel1994(1,:) = log(sales_x_Laspeyres(2,:)) - log(sales_x_Laspeyres(1,:));
    dln_sales_x_rel1994(2,:) = log(sales_x_Laspeyres(3,:)) - log(sales_x_Laspeyres(1,:));
    dln_sales_x_rel1994(3,:) = log(sales_x_Laspeyres(4,:)) - log(sales_x_Laspeyres(1,:));
    dln_sales_x_rel1994(4,:) = log(sales_x_Laspeyres(5,:)) - log(sales_x_Laspeyres(1,:));
    dln_sales_x_rel1994(5,:) = log(sales_x_Laspeyres(6,:)) - log(sales_x_Laspeyres(1,:));

    iy_rel1994(1,:) = inv(2,:)./sales(2,:) - inv(1,:)./sales(1,:);
    iy_rel1994(2,:) = inv(3,:)./sales(3,:) - inv(1,:)./sales(1,:);
    iy_rel1994(3,:) = inv(4,:)./sales(4,:) - inv(1,:)./sales(1,:);
    iy_rel1994(4,:) = inv(5,:)./sales(5,:) - inv(1,:)./sales(1,:);
    iy_rel1994(5,:) = inv(6,:)./sales(6,:) - inv(1,:)./sales(1,:);
    
    z_1994(1,:) = z(1,:);
    z_1994(2,:) = z(1,:);
    z_1994(3,:) = z(1,:);
    z_1994(4,:) = z(1,:);
    z_1994(5,:) = z(1,:);
    
    a_1994(1,:) = a(1,:);
    a_1994(2,:) = a(1,:);
    a_1994(3,:) = a(1,:);
    a_1994(4,:) = a(1,:);
    a_1994(5,:) = a(1,:);    
    
    %Debt
        dummy_1995_debt = zeros(5,s.N);
        dummy_1995_debt(1,firms_debt1994) = 1;

        dummy_1996_debt = zeros(5,s.N);
        dummy_1996_debt(2,firms_debt1994) = 1;

        dummy_1997_debt = zeros(5,s.N);
        dummy_1997_debt(3,firms_debt1994) = 1;

        dummy_1998_debt = zeros(5,s.N);
        dummy_1998_debt(4,firms_debt1994) = 1;

        dummy_1999_debt = zeros(5,s.N);
        dummy_1999_debt(5,firms_debt1994) = 1;    

    %Savings
        dummy_1995_savings = zeros(5,s.N);
        dummy_1995_savings(1,firms_savings1994) = 1;

        dummy_1996_savings = zeros(5,s.N);
        dummy_1996_savings(2,firms_savings1994) = 1;

        dummy_1997_savings = zeros(5,s.N);
        dummy_1997_savings(3,firms_savings1994) = 1;

        dummy_1998_savings = zeros(5,s.N);
        dummy_1998_savings(4,firms_savings1994) = 1;

        dummy_1999_savings = zeros(5,s.N);
        dummy_1999_savings(5,firms_savings1994) = 1;       

    %High export intensity
        dummy_highXY = zeros(5,s.N);
        dummy_highXY(:,N1+1:end) = 1;
        
    %Constant
        constant = ones(5,s.N);

        z = z(2:6,:);
        
    %Valid observations - Exports
        valid = isfinite(dln_sales_x_rel1994);
                
        sim.dln_sales_x_rel1994 = dln_sales_x_rel1994(valid);
        sim.dummy_1995_debt = dummy_1995_debt(valid);
        sim.dummy_1996_debt = dummy_1996_debt(valid);
        sim.dummy_1997_debt = dummy_1997_debt(valid);
        sim.dummy_1998_debt = dummy_1998_debt(valid);
        sim.dummy_1999_debt = dummy_1999_debt(valid);
        sim.dummy_1995_savings = dummy_1995_savings(valid);
        sim.dummy_1996_savings = dummy_1996_savings(valid);
        sim.dummy_1997_savings = dummy_1997_savings(valid);
        sim.dummy_1998_savings = dummy_1998_savings(valid);
        sim.dummy_1999_savings = dummy_1999_savings(valid);
        sim.constant = constant(valid);
        sim.ln_z = log(z(valid));
        sim.ln_a = log(a(valid));
        sim.ln_z_1994 = log(z_1994(valid));
        sim.ln_a_1994 = log(a_1994(valid));        
        sim.dummy_highXY = dummy_highXY(valid);
        
        sim.regcoeffs2 = regress(sim.dln_sales_x_rel1994,[sim.constant ...
            sim.dummy_1996_savings sim.dummy_1997_savings sim.dummy_1998_savings sim.dummy_1999_savings...
            sim.dummy_1995_debt sim.dummy_1996_debt sim.dummy_1997_debt sim.dummy_1998_debt sim.dummy_1999_debt...
            sim.ln_z_1994]);

        sim.dln_sales_x_rel1994_savings = sim.regcoeffs2(1) + [0 sim.regcoeffs2(2) sim.regcoeffs2(3) sim.regcoeffs2(4) sim.regcoeffs2(5)];
        sim.dln_sales_x_rel1994_debt = sim.regcoeffs2(1) + [sim.regcoeffs2(6) sim.regcoeffs2(7) sim.regcoeffs2(8) sim.regcoeffs2(9) sim.regcoeffs2(10)];
    
    %Valid observations - Total sales
        %valid = isfinite(dln_sales_rel1994);
        valid = isfinite(dln_sales_rel1994) & isfinite(dln_sales_x_rel1994);
                
        sim.dln_sales_rel1994 = dln_sales_rel1994(valid);
        sim.dummy_1995_debt = dummy_1995_debt(valid);
        sim.dummy_1996_debt = dummy_1996_debt(valid);
        sim.dummy_1997_debt = dummy_1997_debt(valid);
        sim.dummy_1998_debt = dummy_1998_debt(valid);
        sim.dummy_1999_debt = dummy_1999_debt(valid);
        sim.dummy_1995_savings = dummy_1995_savings(valid);
        sim.dummy_1996_savings = dummy_1996_savings(valid);
        sim.dummy_1997_savings = dummy_1997_savings(valid);
        sim.dummy_1998_savings = dummy_1998_savings(valid);
        sim.dummy_1999_savings = dummy_1999_savings(valid);
        sim.constant = constant(valid);
        sim.ln_z = log(z(valid));
        sim.ln_a = log(a(valid));
        sim.ln_z_1994 = log(z_1994(valid));
        sim.ln_a_1994 = log(a_1994(valid));   
        sim.dummy_highXY = dummy_highXY(valid);
        
        sim.regcoeffs3 = regress(sim.dln_sales_rel1994,[sim.constant ...
            sim.dummy_1996_savings sim.dummy_1997_savings sim.dummy_1998_savings sim.dummy_1999_savings...
            sim.dummy_1995_debt sim.dummy_1996_debt sim.dummy_1997_debt sim.dummy_1998_debt sim.dummy_1999_debt...
            sim.ln_z_1994]);

        sim.dln_sales_rel1994_savings = sim.regcoeffs3(1) + [0 sim.regcoeffs3(2) sim.regcoeffs3(3) sim.regcoeffs3(4) sim.regcoeffs3(5)];
        sim.dln_sales_rel1994_debt = sim.regcoeffs3(1) + [sim.regcoeffs3(6) sim.regcoeffs3(7) sim.regcoeffs3(8) sim.regcoeffs3(9) sim.regcoeffs3(10)];
       
    %Valid observations - Investment
        %valid = isfinite(iy_rel1994);
         valid = isfinite(iy_rel1994) & isfinite(dln_sales_x_rel1994);
         
        sim.iy_rel1994 = iy_rel1994(valid);
        sim.dummy_1995_debt = dummy_1995_debt(valid);
        sim.dummy_1996_debt = dummy_1996_debt(valid);
        sim.dummy_1997_debt = dummy_1997_debt(valid);
        sim.dummy_1998_debt = dummy_1998_debt(valid);
        sim.dummy_1999_debt = dummy_1999_debt(valid);
        sim.dummy_1995_savings = dummy_1995_savings(valid);
        sim.dummy_1996_savings = dummy_1996_savings(valid);
        sim.dummy_1997_savings = dummy_1997_savings(valid);
        sim.dummy_1998_savings = dummy_1998_savings(valid);
        sim.dummy_1999_savings = dummy_1999_savings(valid);
        sim.constant = constant(valid);
        sim.ln_z = log(z(valid));
        sim.ln_a = log(a(valid));
        sim.ln_z_1994 = log(z_1994(valid));
        sim.ln_a_1994 = log(a_1994(valid));   
        sim.dummy_highXY = dummy_highXY(valid);
        
        sim.regcoeffs4 = regress(sim.iy_rel1994,[sim.constant ...
            sim.dummy_1996_savings sim.dummy_1997_savings sim.dummy_1998_savings sim.dummy_1999_savings...
            sim.dummy_1995_debt sim.dummy_1996_debt sim.dummy_1997_debt sim.dummy_1998_debt sim.dummy_1999_debt...
            sim.ln_z_1994]);

        sim.iy_rel1994_savings = sim.regcoeffs4(1) + [0 sim.regcoeffs4(2) sim.regcoeffs4(3) sim.regcoeffs4(4) sim.regcoeffs4(5)];
        sim.iy_rel1994_debt = sim.regcoeffs4(1) + [sim.regcoeffs4(6) sim.regcoeffs4(7) sim.regcoeffs4(8) sim.regcoeffs4(9) sim.regcoeffs4(10)];
                    
    %Transpose
        sim.table_reg = [sim.iy_rel1994_debt' sim.iy_rel1994_savings' sim.dln_sales_rel1994_debt' sim.dln_sales_rel1994_savings' sim.dln_sales_x_rel1994_debt' sim.dln_sales_x_rel1994_savings'];
        
        
    
for t=1:T_trans   
    
    %Debt vs. savings:
    
        debt_share(t,1) = sum(debt_trans(t,:)>0)/s.N;
        debt_savings(t,1) = -sum(debt_trans(t,debt_trans(t,:)>0))/sum(debt_trans(t,debt_trans(t,:)<=0));
          
        X_debt1994(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994))/s.N; 
        X_savings1994(t,1) = sum(sales_x_Laspeyres(t,firms_savings1994))/s.N;  
        
        Y_debt1994(t,1) = sum(sales_Laspeyres(t,firms_debt1994))/s.N; 
        Y_savings1994(t,1) = sum(sales_Laspeyres(t,firms_savings1994))/s.N;     
  
        C_debt1994(t,1) = sum(c_trans(t,firms_debt1994))/s.N; 
        C_savings1994(t,1) = sum(c_trans(t,firms_savings1994))/s.N; 
        
        sharex_debt1994(t,1) = mean(exporter_decision(t,firms_debt1994));
        sharex_savings1994(t,1) = mean(exporter_decision(t,firms_savings1994));
        
        if t>1 && t<T_trans        
            inv_agg_debt1994(t,1) = sum(k_state(t+1,firms_debt1994) - (1-m.delta)*k_state(t,firms_debt1994) )/s.N; 
            inv_agg_savings1994(t,1) = sum(k_state(t+1,firms_savings1994) - (1-m.delta)*k_state(t,firms_savings1994) )/s.N;   
        else       
           inv_agg_debt1994(t,1) = sum(m.delta*k_state(t,firms_debt1994))/s.N;       
           inv_agg_savings1994(t,1) = sum(m.delta*k_state(t,firms_savings1994))/s.N; 
        end            
        
        sim.debt_share = debt_share;
        sim.debt_savings = debt_savings;
        sim.X_debt1994 = X_debt1994;
        sim.X_savings1994 = X_savings1994;
        sim.Y_debt1994 = Y_debt1994;
        sim.Y_savings1994 = Y_savings1994;
        sim.C_debt1994 = C_debt1994;
        sim.C_savings1994 = C_savings1994;
        sim.sharex_debt1994 = sharex_debt1994;
        sim.sharex_savings1994 = sharex_savings1994;
        sim.inv_agg_debt1994 = inv_agg_debt1994;
        sim.inv_agg_savings1994 = inv_agg_savings1994;        
        

    %Debt: binding vs. not binding        
        X_debt1994binding(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994binding))/s.N; 
        X_debt1994notbinding(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994notbinding))/s.N; 
         
        Y_debt1994binding(t,1) = sum(sales_Laspeyres(t,firms_debt1994binding))/s.N; 
        Y_debt1994notbinding(t,1) = sum(sales_Laspeyres(t,firms_debt1994notbinding))/s.N;          
  
        C_debt1994binding(t,1) = sum(c_trans(t,firms_debt1994binding))/s.N; 
        C_debt1994notbinding(t,1) = sum(c_trans(t,firms_debt1994notbinding))/s.N; 
        
        sharex_debt1994binding(t,1) = mean(exporter_decision(t,firms_debt1994binding));
        sharex_debt1994notbinding(t,1) = mean(exporter_decision(t,firms_debt1994notbinding));
      
        if t>1 && t<T_trans        
            inv_agg_debt1994binding(t,1) = sum(k_state(t+1,firms_debt1994binding) - (1-m.delta)*k_state(t,firms_debt1994binding) )/s.N;     
            inv_agg_debt1994notbinding(t,1) = sum(k_state(t+1,firms_debt1994notbinding) - (1-m.delta)*k_state(t,firms_debt1994notbinding) )/s.N;     
        else        
           inv_agg_debt1994binding(t,1) = sum(m.delta*k_state(t,firms_debt1994binding))/s.N; 
           inv_agg_debt1994notbinding(t,1) = sum(m.delta*k_state(t,firms_debt1994notbinding))/s.N;  
        end   

        sim.table_debtbinding_all = [X_debt1994binding, X_debt1994notbinding, Y_debt1994binding, Y_debt1994notbinding, C_debt1994binding, C_debt1994notbinding, sharex_debt1994binding, sharex_debt1994notbinding, inv_agg_debt1994binding, inv_agg_debt1994notbinding];
        
    %Debt vs. savings, exporters  
        sim.debt_share_x(t,1) = sum(debt_trans(t,firms_debt1994_x)>0)/(sum(firms_debt1994_x)+sum(firms_savings1994_x));
        sim.debt_savings_x(t,1) = -sum(debt_trans(t,firms_debt1994_x))/sum(debt_trans(t,firms_savings1994_x));
          
        sim.X_debt1994_x(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994_x))/s.N; 
        sim.X_savings1994_x(t,1) = sum(sales_x_Laspeyres(t,firms_savings1994_x))/s.N;  
        
        sim.Y_debt1994_x(t,1) = sum(sales_Laspeyres(t,firms_debt1994_x))/s.N; 
        sim.Y_savings1994_x(t,1) = sum(sales_Laspeyres(t,firms_savings1994_x))/s.N;     

        sim.C_debt1994_x(t,1) = sum(c_trans(t,firms_debt1994_x))/s.N; 
        sim.C_savings1994_x(t,1) = sum(c_trans(t,firms_savings1994_x))/s.N;          
        
        sim.sharex_debt1994_x(t,1) = mean(exporter_decision(t,firms_debt1994_x));
        sim.sharex_savings1994_x(t,1) = mean(exporter_decision(t,firms_savings1994_x));
        
        if t>1 && t<T_trans        
            sim.inv_agg_debt1994_x(t,1) = sum(k_state(t+1,firms_debt1994_x) - (1-m.delta)*k_state(t,firms_debt1994_x) )/s.N; 
            sim.inv_agg_savings1994_x(t,1) = sum(k_state(t+1,firms_savings1994_x) - (1-m.delta)*k_state(t,firms_savings1994_x) )/s.N;   
        else       
           sim.inv_agg_debt1994_x(t,1) = sum(m.delta*k_state(t,firms_debt1994_x))/s.N;       
           sim.inv_agg_savings1994_x(t,1) = sum(m.delta*k_state(t,firms_savings1994_x))/s.N; 
        end         
        

    %Debt vs. savings, non-exporters  
        sim.debt_share_nx(t,1) = sum(debt_trans(t,firms_debt1994_nx)>0)/(sum(firms_debt1994_nx)+sum(firms_savings1994_nx));
        sim.debt_savings_nx(t,1) = -sum(debt_trans(t,firms_debt1994_nx))/sum(debt_trans(t,firms_savings1994_nx));
          
        sim.X_debt1994_nx(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994_nx))/s.N; 
        sim.X_savings1994_nx(t,1) = sum(sales_x_Laspeyres(t,firms_savings1994_nx))/s.N;  
        
        sim.Y_debt1994_nx(t,1) = sum(sales_Laspeyres(t,firms_debt1994_nx))/s.N; 
        sim.Y_savings1994_nx(t,1) = sum(sales_Laspeyres(t,firms_savings1994_nx))/s.N;     
 
        sim.C_debt1994_nx(t,1) = sum(c_trans(t,firms_debt1994_nx))/s.N; 
        sim.C_savings1994_nx(t,1) = sum(c_trans(t,firms_savings1994_nx))/s.N;         
        
        sim.sharex_debt1994_nx(t,1) = mean(exporter_decision(t,firms_debt1994_nx));
        sim.sharex_savings1994_nx(t,1) = mean(exporter_decision(t,firms_savings1994_nx));
        
        if t>1 && t<T_trans        
            sim.inv_agg_debt1994_nx(t,1) = sum(k_state(t+1,firms_debt1994_nx) - (1-m.delta)*k_state(t,firms_debt1994_nx) )/s.N; 
            sim.inv_agg_savings1994_nx(t,1) = sum(k_state(t+1,firms_savings1994_nx) - (1-m.delta)*k_state(t,firms_savings1994_nx) )/s.N;   
        else       
           sim.inv_agg_debt1994_nx(t,1) = sum(m.delta*k_state(t,firms_debt1994_nx))/s.N;       
           sim.inv_agg_savings1994_nx(t,1) = sum(m.delta*k_state(t,firms_savings1994_nx))/s.N; 
        end           

    %Debt vs. savings, EXPORTERS by export intensity categories  

        %Debt vs. savings, high export intensity
        sim.debt_share_highXY(t,1) = sum(debt_trans(t,firms_debt1994_highXY)>0)/(sum(firms_debt1994_highXY)+sum(firms_savings1994_highXY));
        sim.debt_savings_highXY(t,1) = -sum(debt_trans(t,firms_debt1994_highXY))/sum(debt_trans(t,logical(1-firms_debt1994_highXY)));
          
        sim.X_debt1994_highXY(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994_highXY))/s.N; 
        sim.X_savings1994_highXY(t,1) = sum(sales_x_Laspeyres(t,firms_savings1994_highXY))/s.N;  
        
        sim.Y_debt1994_highXY(t,1) = sum(sales_Laspeyres(t,firms_debt1994_highXY))/s.N; 
        sim.Y_savings1994_highXY(t,1) = sum(sales_Laspeyres(t,firms_savings1994_highXY))/s.N;     
  
        sim.C_debt1994_highXY(t,1) = sum(c_trans(t,firms_debt1994_highXY))/s.N; 
        sim.C_savings1994_highXY(t,1) = sum(c_trans(t,firms_savings1994_highXY))/s.N;             
        
        sim.sharex_debt1994_highXY(t,1) = mean(exporter_decision(t,firms_debt1994_highXY));
        sim.sharex_savings1994_highXY(t,1) = mean(exporter_decision(t,firms_savings1994_highXY));
        
        if t>1 && t<T_trans        
            sim.inv_agg_debt1994_highXY(t,1) = sum(k_state(t+1,firms_debt1994_highXY) - (1-m.delta)*k_state(t,firms_debt1994_highXY) )/s.N; 
            sim.inv_agg_savings1994_highXY(t,1) = sum(k_state(t+1,firms_savings1994_highXY) - (1-m.delta)*k_state(t,firms_savings1994_highXY) )/s.N;   
        else       
           sim.inv_agg_debt1994_highXY(t,1) = sum(m.delta*k_state(t,firms_debt1994_highXY))/s.N;       
           sim.inv_agg_savings1994_highXY(t,1) = sum(m.delta*k_state(t,firms_savings1994_highXY))/s.N; 
        end         
        

    %Debt vs. savings, low export intensity  
        sim.debt_share_lowXY(t,1) = sum(debt_trans(t,firms_debt1994_lowXY)>0)/(sum(firms_debt1994_lowXY)+sum(firms_savings1994_lowXY));
        sim.debt_savings_lowXY(t,1) = -sum(debt_trans(t,firms_debt1994_lowXY))/sum(debt_trans(t,logical(1-firms_debt1994_lowXY)));
          
        sim.X_debt1994_lowXY(t,1) = sum(sales_x_Laspeyres(t,firms_debt1994_lowXY))/s.N; 
        sim.X_savings1994_lowXY(t,1) = sum(sales_x_Laspeyres(t,firms_savings1994_lowXY))/s.N;  
        
        sim.Y_debt1994_lowXY(t,1) = sum(sales_Laspeyres(t,firms_debt1994_lowXY))/s.N; 
        sim.Y_savings1994_lowXY(t,1) = sum(sales_Laspeyres(t,firms_savings1994_lowXY))/s.N;     

        sim.C_debt1994_lowXY(t,1) = sum(c_trans(t,firms_debt1994_lowXY))/s.N; 
        sim.C_savings1994_lowXY(t,1) = sum(c_trans(t,firms_savings1994_lowXY))/s.N;     
        
        sim.sharex_debt1994_lowXY(t,1) = mean(exporter_decision(t,firms_debt1994_lowXY));
        sim.sharex_savings1994_lowXY(t,1) = mean(exporter_decision(t,firms_savings1994_lowXY));
        
        if t>1 && t<T_trans        
            sim.inv_agg_debt1994_lowXY(t,1) = sum(k_state(t+1,firms_debt1994_lowXY) - (1-m.delta)*k_state(t,firms_debt1994_lowXY) )/s.N; 
            sim.inv_agg_savings1994_lowXY(t,1) = sum(k_state(t+1,firms_savings1994_lowXY) - (1-m.delta)*k_state(t,firms_savings1994_lowXY) )/s.N;   
        else       
           sim.inv_agg_debt1994_lowXY(t,1) = sum(m.delta*k_state(t,firms_debt1994_lowXY))/s.N;       
           sim.inv_agg_savings1994_lowXY(t,1) = sum(m.delta*k_state(t,firms_savings1994_lowXY))/s.N; 
        end       
        
end



%% More Micro-Level Statistics by for median productivity exporter (Debt vs Savings)

%  Productivity groups (Debt vs Savings, Exporters)
    debt_share_zcat = zeros(T_trans,1);
    debt_savings_zcat = zeros(T_trans,1);
    D_debt1994_zcat = zeros(T_trans,1);
    D_savings1994_zcat = zeros(T_trans,1);    
    X_debt1994_zcat = zeros(T_trans,1);
    X_savings1994_zcat = zeros(T_trans,1);
    Y_debt1994_zcat = zeros(T_trans,1);
    Y_savings1994_zcat = zeros(T_trans,1);
    C_debt1994_zcat = zeros(T_trans,1);
    C_savings1994_zcat = zeros(T_trans,1);
    sharex_debt1994_zcat = zeros(T_trans,1);
    sharex_savings1994_zcat = zeros(T_trans,1);
    inv_agg_debt1994_zcat = zeros(T_trans,1);
    inv_agg_savings1994_zcat = zeros(T_trans,1);
    
    indz_p50_exporter = prctile(ind_z(1,exporter_decision(1,:)==1),50);
    z_cat = ind_z(1,:)==indz_p50_exporter;

    for t=1:T_trans     
        
    %Debt vs. savings
        debt_share_zcat(t) = sum(debt_trans(t,firms_debt1994_x==1 & z_cat==1)>0)/(sum(firms_debt1994_x.*z_cat)+sum(firms_savings1994_x.*z_cat));
        debt_savings_zcat(t) = -sum(debt_trans(t,firms_debt1994_x==1 & z_cat==1))/sum(debt_trans(t,firms_savings1994_x==1 & z_cat==1));
          
        D_debt1994_zcat(t) = sum(sales_d_Laspeyres(t,firms_debt1994_x==1 & z_cat==1))/s.N; 
        D_savings1994_zcat(t) = sum(sales_d_Laspeyres(t,firms_savings1994_x==1 & z_cat==1))/s.N;          
        
        X_debt1994_zcat(t) = sum(sales_x_Laspeyres(t,firms_debt1994_x==1 & z_cat==1))/s.N; 
        X_savings1994_zcat(t) = sum(sales_x_Laspeyres(t,firms_savings1994_x==1 & z_cat==1))/s.N;  
        
        Y_debt1994_zcat(t) = sum(sales_Laspeyres(t,firms_debt1994_x==1 & z_cat==1))/s.N; 
        Y_savings1994_zcat(t) = sum(sales_Laspeyres(t,firms_savings1994_x==1 & z_cat==1))/s.N;     

        C_debt1994_zcat(t) = sum(c_trans(t,firms_debt1994_x==1 & z_cat==1))/s.N; 
        C_savings1994_zcat(t) = sum(c_trans(t,firms_savings1994_x==1 & z_cat==1))/s.N;          
        
        sharex_debt1994_zcat(t) = mean(exporter_decision(t,firms_debt1994_x==1 & z_cat==1));
        sharex_savings1994_zcat(t) = mean(exporter_decision(t,firms_savings1994_x==1 & z_cat==1));
        
        if t>1 && t<T_trans        
            inv_agg_debt1994_zcat(t) = sum(k_state(t+1,firms_debt1994_x==1 & z_cat==1) - (1-m.delta)*k_state(t,firms_debt1994_x==1 & z_cat==1) )/s.N; 
            inv_agg_savings1994_zcat(t) = sum(k_state(t+1,firms_savings1994_x==1 & z_cat==1) - (1-m.delta)*k_state(t,firms_savings1994_x==1 & z_cat==1) )/s.N;   
        else       
           inv_agg_debt1994_zcat(t) = sum(m.delta*k_state(t,firms_debt1994_x==1 & z_cat==1))/s.N;       
           inv_agg_savings1994_zcat(t) = sum(m.delta*k_state(t,firms_savings1994_x==1 & z_cat==1))/s.N; 
        end         

    end   
    
    % number of firms with debt and savings in each category
    sim.N_zcat_debt_x(1) = sum(debt_trans(t,firms_debt1994_x==1 & z_cat==1)>0);
    sim.N_zcat_sav_x(1) = sum(debt_trans(t,firms_savings1994_x==1 & z_cat==1)>0);

     
    sim.table_zcat_x = [X_debt1994_zcat(:), X_savings1994_zcat(:), Y_debt1994_zcat(:), Y_savings1994_zcat(:), inv_agg_debt1994_zcat(:), inv_agg_savings1994_zcat(:), debt_share_zcat(:), D_debt1994_zcat(:), D_savings1994_zcat(:)];           
    


%% Output:

        % saving results for all firms:
        sim.table_debtsavings1994 = [X_debt1994(1:10), X_savings1994(1:10),Y_debt1994(1:10), Y_savings1994(1:10), inv_agg_debt1994(1:10), inv_agg_savings1994(1:10), inv_agg_debt1994(1:10)./Y_debt1994(1:10), inv_agg_savings1994(1:10)./Y_savings1994(1:10), debt_share(1:10)];        
        sim.table_debtsavings1994_avg = [X_debt1994(1:10)*(s.N/sim.N_debt), X_savings1994(1:10)*(s.N/sim.N_savings),Y_debt1994(1:10)*(s.N/sim.N_debt), Y_savings1994(1:10)*(s.N/sim.N_savings), inv_agg_debt1994(1:10)*(s.N/sim.N_debt), inv_agg_savings1994(1:10)*(s.N/sim.N_savings), debt_share(1:10)]; 
        sim.table_debtsavings1994_all = [debt_share, debt_savings, X_debt1994, X_savings1994, Y_debt1994, Y_savings1994, C_debt1994, C_savings1994, sharex_debt1994, sharex_savings1994, inv_agg_debt1994, inv_agg_savings1994];
        
        % saving results for exporters:
        avg_debt = (s.N/sim.N_debt_x);
        avg_sav = (s.N/sim.N_savings_x);
        sim.table_exporters1994 = [sim.X_debt1994_x(1:10), sim.X_savings1994_x(1:10), sim.Y_debt1994_x(1:10), sim.Y_savings1994_x(1:10), sim.inv_agg_debt1994_x(1:10), sim.inv_agg_savings1994_x(1:10), sim.inv_agg_debt1994_x(1:10)./sim.Y_debt1994_x(1:10), sim.inv_agg_savings1994_x(1:10)./sim.Y_savings1994_x(1:10), sim.debt_share_x(1:10)];  
        sim.table_exporters1994_avg = [sim.X_debt1994_x(1:10)*avg_debt, sim.X_savings1994_x(1:10)*avg_sav, sim.Y_debt1994_x(1:10)*avg_debt, sim.Y_savings1994_x(1:10)*avg_sav, sim.inv_agg_debt1994_x(1:10)*avg_debt, sim.inv_agg_savings1994_x(1:10)*avg_sav, sim.debt_share_x(1:10)];  
        sim.table_exporters1994_all = [sim.debt_share_x, sim.debt_savings_x, sim.X_debt1994_x, sim.X_savings1994_x, sim.Y_debt1994_x, sim.Y_savings1994_x, sim.C_debt1994_x, sim.C_savings1994_x, sim.sharex_debt1994_x, sim.sharex_savings1994_x, sim.inv_agg_debt1994_x, sim.inv_agg_savings1994_x];
        
        % saving results for non-exporters:
        avg_debt = (s.N/sim.N_debt_nx);
        avg_sav = (s.N/sim.N_savings_nx);        
        sim.table_nonexporters1994 = [sim.X_debt1994_nx(1:10), sim.X_savings1994_nx(1:10), sim.Y_debt1994_nx(1:10), sim.Y_savings1994_nx(1:10), sim.inv_agg_debt1994_nx(1:10), sim.inv_agg_savings1994_nx(1:10), sim.inv_agg_debt1994_nx(1:10)./sim.Y_debt1994_nx(1:10), sim.inv_agg_savings1994_nx(1:10)./sim.Y_savings1994_nx(1:10) ,sim.debt_share_nx(1:10)];
        sim.table_nonexporters1994_avg = [sim.X_debt1994_nx(1:10)*avg_debt, sim.X_savings1994_nx(1:10)*avg_sav, sim.Y_debt1994_nx(1:10)*avg_debt, sim.Y_savings1994_nx(1:10)*avg_sav, sim.inv_agg_debt1994_nx(1:10)*avg_debt, sim.inv_agg_savings1994_nx(1:10)*avg_sav, sim.debt_share_nx(1:10)];                
        sim.table_nonexporters1994_all = [sim.debt_share_nx, sim.debt_savings_nx, sim.X_debt1994_nx, sim.X_savings1994_nx, sim.Y_debt1994_nx, sim.Y_savings1994_nx, sim.C_debt1994_nx, sim.C_savings1994_nx, sim.sharex_debt1994_nx, sim.sharex_savings1994_nx, sim.inv_agg_debt1994_nx, sim.inv_agg_savings1994_nx];
        
        % saving results for highXY:
        avg_debt = (s.N/sim.N_debt_high_XY);
        avg_sav = (s.N/sim.N_savings_high_XY); 
        sim.table_highXY1994 = [sim.X_debt1994_highXY(1:10), sim.X_savings1994_highXY(1:10), sim.Y_debt1994_highXY(1:10), sim.Y_savings1994_highXY(1:10), sim.inv_agg_debt1994_highXY(1:10), sim.inv_agg_savings1994_highXY(1:10),  sim.inv_agg_debt1994_highXY(1:10)./sim.Y_debt1994_highXY(1:10), sim.inv_agg_savings1994_highXY(1:10)./sim.Y_savings1994_highXY(1:10), sim.debt_share_highXY(1:10)];
        sim.table_highXY1994_avg = [sim.X_debt1994_highXY(1:10)*avg_debt, sim.X_savings1994_highXY(1:10)*avg_sav, sim.Y_debt1994_highXY(1:10)*avg_debt, sim.Y_savings1994_highXY(1:10)*avg_sav, sim.inv_agg_debt1994_highXY(1:10)*avg_debt, sim.inv_agg_savings1994_highXY(1:10)*avg_sav, sim.debt_share_highXY(1:10)];
        sim.table_highXY1994_all = [sim.debt_share_highXY, sim.debt_savings_highXY, sim.X_debt1994_highXY, sim.X_savings1994_highXY, sim.Y_debt1994_highXY, sim.Y_savings1994_highXY, sim.C_debt1994_highXY, sim.C_savings1994_highXY, sim.sharex_debt1994_highXY, sim.sharex_savings1994_highXY, sim.inv_agg_debt1994_highXY, sim.inv_agg_savings1994_highXY];
                
        % saving results for lowXY:
        avg_debt = (s.N/sim.N_debt_low_XY);
        avg_sav = (s.N/sim.N_savings_low_XY);        
        sim.table_lowXY1994 = [sim.X_debt1994_lowXY(1:10), sim.X_savings1994_lowXY(1:10), sim.Y_debt1994_lowXY(1:10), sim.Y_savings1994_lowXY(1:10), sim.inv_agg_debt1994_lowXY(1:10), sim.inv_agg_savings1994_lowXY(1:10), sim.inv_agg_debt1994_lowXY(1:10)./sim.Y_debt1994_lowXY(1:10), sim.inv_agg_savings1994_lowXY(1:10)./sim.Y_savings1994_lowXY(1:10), sim.X_debt1994_lowXY(1:10)];
        sim.table_lowXY1994_avg = [sim.X_debt1994_lowXY(1:10)*avg_debt, sim.X_savings1994_lowXY(1:10)*avg_sav, sim.Y_debt1994_lowXY(1:10)*avg_debt, sim.Y_savings1994_lowXY(1:10)*avg_sav, sim.inv_agg_debt1994_lowXY(1:10)*avg_debt, sim.inv_agg_savings1994_lowXY(1:10)*avg_sav, sim.X_debt1994_lowXY(1:10)];
        sim.table_lowXY1994_all = [sim.debt_share_lowXY, sim.debt_savings_lowXY, sim.X_debt1994_lowXY, sim.X_savings1994_lowXY, sim.Y_debt1994_lowXY, sim.Y_savings1994_lowXY, sim.C_debt1994_lowXY, sim.C_savings1994_lowXY, sim.sharex_debt1994_lowXY, sim.sharex_savings1994_lowXY, sim.inv_agg_debt1994_lowXY, sim.inv_agg_savings1994_lowXY];
        
        
%% Aguiar Figures

indsales_p90 = prctile(sales_Laspeyres(1,:),90);
indsales_p60 = prctile(sales_Laspeyres(1,:),40);

sales_cat{1} = sales_Laspeyres(1,:)>=0; %All firms
sales_cat{2} = sales_Laspeyres(1,:)>=indsales_p60;
sales_cat{3} = sales_Laspeyres(1,:)>=indsales_p90;
sample = sales_cat{2};


sales_growth_x = zeros(T_trans ,1);
sales_growth_nx = zeros(T_trans,1);

for t=2:T_trans  
    sales_growth_x(t,1) = sum(log(sales_Laspeyres(t,exporter_decision(1,sample)==1)) - log(sales_Laspeyres(1,exporter_decision(1,sample)==1)))/s.N;
    sales_growth_nx(t,1) = sum(log(sales_Laspeyres(t,exporter_decision(1,sample)==0)) - log(sales_Laspeyres(1,exporter_decision(1,sample)==0)))/s.N;
end

inv_x = zeros(T_trans,1);
inv_nx = zeros(T_trans,1);

for t=2:T_trans  

    if t>1 && t<T_trans        
        inv_x(t,1) = sum(k_state(t+1,exporter_decision(1,sample)==1) - (1-m.delta)*k_state(t,exporter_decision(1,sample)==1) )/s.N; 
        inv_nx(t,1) = sum(k_state(t+1,exporter_decision(1,sample)==0) - (1-m.delta)*k_state(t,exporter_decision(1,sample)==0) )/s.N;   
    else     
       inv_x(t,1) = sum(m.delta*k_state(t,exporter_decision(1,sample)==1))/s.N;       
       inv_nx(t,1) = sum(m.delta*k_state(t,exporter_decision(1,sample)==0))/s.N;    
    end

end  

for t=1:T_trans  
    sales_total_Laspeyres_x(t,:) = sum(sales_Laspeyres(t,exporter_decision(1,sample)==1))/s.N;   
    sales_total_Laspeyres_nx(t,:) = sum(sales_Laspeyres(t,exporter_decision(1,sample)==0))/s.N;   
end

sim.Aguiar = [sales_growth_x sales_growth_nx inv_x inv_nx sales_total_Laspeyres_x sales_total_Laspeyres_nx];

%% Saving debt matrix 
sim.debt_trans = debt_trans;
sim.lm_trans = lm_trans;


%% Clear

clear   firms_debt1994binding firms_debt1994notbinding firms_debt1994_x firms_savings1994_x...
        firms_savings1994_x firms_debt1994_nx firms_savings1994_nx firms_debt1994_lowXY...
        firms_savings1994_lowXY firms_debt1994_highXY firms_savings1994_highXY...
        firms_debt1994 firms_savings1994
    
    
    


    