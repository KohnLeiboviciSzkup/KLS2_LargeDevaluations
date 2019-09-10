function r = KLS2_dynamicproblem_trans(m,s,r,vp,period2,varargin)

    if isempty(varargin) || varargin{1}==1 %High export cost firms                
        FC = m.F_base*(s.fcost_fgoods+(1-s.fcost_fgoods)*m.w);
        if period2==2
		
		
            if s.high_z_ForDebt == 0  
				lambda=m.lambda_period_2;
			elseif s.high_z_ForDebt == 1
				lambda=m.lambda_period_2_mat;
			end
			
			
        elseif period2==1
		    
			if s.high_z_ForDebt == 0  
				lambda=m.lambda;
			elseif s.high_z_ForDebt == 1
				lambda=m.lambda_mat;
			end
		
        end
    elseif varargin{1}==2 %Low export cost firms                
        FC = m.F_X_base*(s.fcost_fgoods+(1-s.fcost_fgoods)*m.w);
        if period2==2
      
            if s.high_z_ForDebt == 0  
				lambda=m.lambda_period_2_X;
			elseif s.high_z_ForDebt == 1
				lambda=m.lambda_X_period_2_mat;
			end			
			
			
        elseif period2==1
            
			
			if s.high_z_ForDebt == 0  
				lambda=m.lambda_X;
			elseif s.high_z_ForDebt == 1
				lambda=m.lambda_X_mat;
			end
			
        end
    end
       
%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)

	FXdebt = lambda + (1-lambda).*(m.xi_real/m.xi_real_lag);    
   
%% Dynamic problem: value functions and policy functions.    

    v_new = ((r.pi_nx.^(1-m.gamma))./(1-m.gamma))./(1-m.beta);
    ap_ind = zeros(s.a_grid_size,s.z_grid_size);

%Export policy function       
    
    e = zeros(s.a_grid_size,s.z_grid_size);
    for j = 1:s.z_grid_size
        for i = 1:s.a_grid_size               
                e(i,j) = ( r.pi_x(i,j)- FC >= r.pi_nx(i,j) );
        end
    
    end
    
%Value function iteration algorithm    
    vp_t=vp';           
     
    c_aux =  m.w*(1-m.zeta_w) +  (1-m.zeta_pi)*(e.* (r.pd_x.*r.yd_x + m.xi_real*r.yf_x.*r.pf_x - m.w*r.n_x - FC)...
                              +   (1- e).*(  r.pd_nx.*r.yd_nx   - m.w*r.n_nx ) ) ...
                              - (1+m.r).*FXdebt.*((e.*r.k_x + (1- e).*r.k_nx) - r.a_grid'*ones(1,s.z_grid_size) )...
                              +  (1-m.delta)*( e.*r.k_x + (1- e).*r.k_nx);
    

    
    
    for j = 1:s.z_grid_size
        for i = 1:s.a_grid_size
             c = c_aux(i,j) - r.a_grid; 
             neg_c_indexes = c<=0; %indices of a' for which consumption<0
             u = ((c).^(1-m.gamma))./(1-m.gamma); 
             u(neg_c_indexes) = -1e50; % Only allow for positive consumption
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
    r.e=e;
        
%% Store output to be used in simulation
     
    r.pd = (1-r.e).*r.pd_nx + r.e.*r.pd_x;
    r.yd = (1-r.e).*r.yd_nx + r.e.*r.yd_x;

    r.pf = (1-r.e).*r.pf_nx + r.e.*r.pf_x;
    r.yf = (1-r.e).*r.yf_nx + r.e.*r.yf_x;

    r.k = (1-r.e).*r.k_nx + r.e.*r.k_x;
    r.n = (1-r.e).*r.n_nx + r.e.*r.n_x;
    r.pi = (1-r.e).*r.pi_nx + r.e.*r.pi_x;

%     if period2==1 %Not period 2
%         r.lm = (1-r.e).*r.lm_nx + r.e.*r.lm_x;
%         r.const = r.lm>0;  
%     end
    
    r.lm = (1-r.e).*r.lm_nx + r.e.*r.lm_x;
    r.const = r.lm>0;  

    r.debt = (r.k-r.a_grid_mat)*(1+m.r);    
    
end