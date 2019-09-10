function r = KLS2_dynamicproblem_typeX_trans(m,s,r,vp,period2)

                  
   FC = 0;
               

%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag=m.Pf_lag/m.P_lag;
    if period2==2
        FXdebt =m.lambda_period_2_X + (1-m.lambda_period_2_X)*((m.xi_real/m.xi_real_lag)); 
    elseif period2==1
        FXdebt = m.lambda_X + (1-m.lambda_X)*(m.xi_real/m.xi_real_lag);    
    end

%% Dynamic problem: value functions and policy functions.    


    v_new = ((r.pi_x.^(1-m.gamma))./(1-m.gamma))./(1-m.beta);
    ap_ind = zeros(size(v_new));    
   
   
%Value function iteration algorithm           
    vp_t=vp';            
    %c_aux =  m.w*(1-m.zeta_w) +  (1-m.zeta_pi)*(r.pi_x- FC)  + r.a_grid'*ones(1,s.z_grid_size)*(1+m.r)*FXdebt;   
    
    c_aux =  m.w*(1-m.zeta_w) +  (1-m.zeta_pi)*( r.pd_x.*r.yd_x + m.xi_real*r.yf_x.*r.pf_x - m.w*r.n_x - FC)...
                              - (1+m.r)*FXdebt*(r.k_x - r.a_grid'*ones(1,s.z_grid_size) ) +  (1-m.delta)*r.k_x ;

    
%     if s.timing==2 && period2==2
%         
%         c_aux =  (m.P_lag/m.P) *(m.w*(1-m.zeta_w) +  (1-m.zeta_pi)*( r.pd_x.*r.yd_x + m.xi_real*r.yf_x.*r.pf_x - m.w*r.n_x - FC))...
%                               - (1+m.r)*FXdebt*(r.k_x - r.a_grid'*ones(1,s.z_grid_size) ) +  (1-m.delta)*r.k_x+ s.flag_T*m.Y*(1-(m.P_lag/m.P));
%         
%     end
    
    
    for j = 1:s.z_grid_size
        for i = 1:s.a_grid_size
             c = c_aux(i,j) - r.a_grid; 
             neg_c_indexes = c<=0; %Indices of a' for which consumption<0
             u = ((c).^(1-m.gamma))./(1-m.gamma); 
             u(neg_c_indexes) = -1e50; %Only allow for positive consumption
             [v,index] = max(u + m.beta*r.z_P(j,:)*vp_t);

             v_new(i,j) = v;                  
             ap_ind(i,j)=index;                           
        end
    end
    
%Store output from value function iteration            
    r.ap = r.a_grid(ap_ind);
    r.ap_ind = ap_ind;
    r.c = c_aux - r.ap;
    r.v = v_new;
    r.e=ones(s.a_grid_size,s.z_grid_size);
        
%% Store output to be used in simulation
     
    r.pd = zeros(size(r.e));
    r.yd = zeros(size(r.e));

    r.pf = r.pf_x;
    r.yf = r.yf_x;

    r.k = r.k_x;
    r.n = r.n_x;
    r.pi = r.pi_x;

    if period2==1
        r.lm = r.lm_x;
        r.const = r.lm>0;  
    end
    
end