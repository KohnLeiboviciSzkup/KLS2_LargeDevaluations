function r = KLS2_dynamicproblem_typeX(m,s,r)

%Export costs (we assume that these firms pay no fixed cost)
    FC = 0;   
    
%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag = m.xi_real; %m.Pf_lag/m.P_lag;
    FXdebt = m.lambda + (1-m.lambda)*(m.xi_real/m.xi_real_lag); 
    
%% Dynamic problem: value functions and policy functions.    

%Initialize solution objects 
    %Value function
        v_new = ((r.pi_x.^(1-m.gamma))./(1-m.gamma))./(1-m.beta);    
        v_old = v_new;  
        
    %Policy functions
        e = ones(s.a_grid_size,s.z_grid_size);
        ap = zeros(s.a_grid_size,s.z_grid_size);
        ap_ind = zeros(s.a_grid_size,s.z_grid_size); 

    %Value function iteration algorithm    
        iter = 0;    
        v_diff = 1;    
 
    %Asset grid
        r.agrid_T = r.a_grid';          
        
         
    while v_diff>s.eps

        iter = iter + 1;
        v_old_T = v_old';
        
        for j = 1:s.z_grid_size                
            for i = 1:s.a_grid_size
                                              
                c = m.w*(1-m.zeta_w) +  (1-m.zeta_pi)*(r.pi_x(i,j)- FC)  + r.a_grid(i)*(1+m.r)*FXdebt - r.a_grid ; 
                u = ((c).^(1-m.gamma))./(1-m.gamma); 
                neg_c_indexes = c<=0; %indices of a' for which consumption<0
                u(neg_c_indexes) = -1e50; % Only allow for positive consumption
                [v,index] = max(u + m.beta*r.z_P(j,:)*v_old_T);
                
                v_new(i,j) = v;                 
                ap(i,j) = r.a_grid(index); 
                ap_ind(i,j)=index;                                       
              
            end
        end          
        
        %Consumption and utility 
        c = m.w*(1-m.zeta_w) +  (1-m.zeta_pi)* (r.pi_x-FC)  + r.agrid_T*ones(1,s.z_grid_size).*(1+m.r).*FXdebt - ap;
        u = ((c).^(1-m.gamma))./(1-m.gamma);  

        %Accelerator
        for g=1:s.ac_iter            
            for j = 1:s.z_grid_size
                   v_new(:,j) =  u(:,j) + (m.beta*r.z_P(j,:)*v_new(ap_ind(:,j),:)')';                 
            end   
        end

        %Convergence                
        v_diff = max(max(abs(v_new-v_old)));
        v_old  = v_new;         

        if mod(iter,100)==0
            disp(['Delta V: ' num2str(v_diff)]);
        end        
    end
    
%Store output from value function iteration       
    r.v = v_new;
    r.ap = ap;
    r.ap_ind = ap_ind;
    r.c =c;
    r.e =e;
    
%% Store output to be used in simulation

    r.pd = zeros(size(r.e));
    r.yd = zeros(size(r.e));
    
    r.pf = r.pf_x;    
    r.yf = r.yf_x;  
    
    r.k = r.k_x;        
    r.n = r.n_x;        
    r.pi = r.pi_x;  
    
    r.lm = r.lm_x;        
    r.const = r.lm>0;
    
      
%Fixed costs
 % We assume that these firms pay no fixed cost TC=0
    r.F_mat = zeros(s.a_grid_size,s.z_grid_size); 
    
    %Adjust fixed costs if denominated in units of labor    
    if s.fcost_fgoods==0
        r.F_mat = m.w*r.F_mat;
    end 
        
end