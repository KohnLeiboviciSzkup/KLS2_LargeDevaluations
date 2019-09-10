function r = KLS2_dynamicproblem(m,s,r,varargin)

    if isempty(varargin) || varargin{1}==1 %High export cost firms                
        FC = m.F_base*(s.fcost_fgoods+(1-s.fcost_fgoods)*m.w);
    elseif varargin{1}==2 %Low export cost firms                
        FC = m.F_X_base*(s.fcost_fgoods+(1-s.fcost_fgoods)*m.w);
    end

%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag = m.xi_real; % For now, in SS
    FXdebt = m.lambda + (1-m.lambda)*(m.xi_real/m.xi_real_lag);    
	
	if s.high_z_ForDebt == 1
		FXdebt = m.lambda_mat + (1-m.lambda_mat).*(m.xi_real/m.xi_real_lag);    
   	end	
	
   
%% Dynamic problem: value functions and policy functions    

%Initialize solution objects 
%Value function
    v_new = ((r.pi_nx.^(1-m.gamma))./(1-m.gamma))./(1-m.beta);   
    v_old = v_new;         

%Policy functions
    e = zeros(size(v_new));
    ap = zeros(size(v_new));
    ap_ind = zeros(size(v_new));    

%Exports
    for i = 1:s.a_grid_size
        for j = 1:s.z_grid_size                
            e(i,j) = r.pi_x(i,j)-FC >= r.pi_nx(i,j);
        end
    end
    
%Productivity process
    r.zP_T = r.z_P';

%Asset grid
    r.agrid_T = r.a_grid';        

%Value function iteration algorithm    
iter = 0;    
v_diff = 1;
    
    while v_diff>s.eps

        iter = iter + 1;
        v_old_T = v_old';
        
        for i = 1:s.a_grid_size
            for j = 1:s.z_grid_size  

				
                if s.high_z_ForDebt == 0                              
					c = m.w*(1-m.zeta_w) + (1-m.zeta_pi)*(e(i,j)*(r.pi_x(i,j)-FC) + (1-e(i,j))*r.pi_nx(i,j)) + r.a_grid(i)*(1+m.r)*FXdebt - r.a_grid;
                elseif s.high_z_ForDebt == 1 
					c = m.w*(1-m.zeta_w) + (1-m.zeta_pi)*(e(i,j)*(r.pi_x(i,j)-FC) + (1-e(i,j))*r.pi_nx(i,j)) + r.a_grid(i)*(1+m.r)*FXdebt(i,j) - r.a_grid;
				end
				
				u = (c.^(1-m.gamma))./(1-m.gamma); 
                neg_c_indexes = c<=0; %Indices of a' for which consumption<0                
                u(neg_c_indexes) = -1e50; %Only allow for positive consumption
                [v,index] = max(u + m.beta*r.z_P(j,:)*v_old_T);
                
                v_new(i,j) = v;                 
                ap(i,j) = r.a_grid(index); 
                ap_ind(i,j) = index;                                               
               
            end
        end          
        
        %Consumption and utility 
        c = m.w*(1-m.zeta_w) + (1-m.zeta_pi)*(e.*(r.pi_x-FC) + (1-e).*r.pi_nx) + r.agrid_T*ones(1,s.z_grid_size).*(1+m.r).*FXdebt - ap;
        u = (c.^(1-m.gamma))./(1-m.gamma);  

        %Accelerator
        for g = 1:s.ac_iter
            for j = 1:s.z_grid_size           
                v_new(:,j) = u(:,j) + m.beta*v_new(ap_ind(:,j),:)*r.zP_T(:,j);               
            end   
        end

        %Convergence
        v_diff = max(max(abs(v_new-v_old)));
        v_old = v_new;         

        if mod(iter,100)==0
            disp(['Delta V: ' num2str(v_diff)]);
        end        
    end
    
%Store output from value function iteration    
    r.v = v_new;
    r.ap = ap;    
    r.ap_ind = ap_ind;    
    r.c = c;    
    r.e = e;
    
%% Store output to be used in simulation
     
    r.pd = (1-r.e).*r.pd_nx + r.e.*r.pd_x;
    r.yd = (1-r.e).*r.yd_nx + r.e.*r.yd_x;

    r.pf = (1-r.e).*r.pf_nx + r.e.*r.pf_x;
    r.yf = (1-r.e).*r.yf_nx + r.e.*r.yf_x;

    r.k = (1-r.e).*r.k_nx + r.e.*r.k_x;
    r.n = (1-r.e).*r.n_nx + r.e.*r.n_x;
    r.pi = (1-r.e).*r.pi_nx + r.e.*r.pi_x;

    r.lm = (1-r.e).*r.lm_nx + r.e.*r.lm_x;
    r.const = r.lm>0;  
    
    r.debt = (r.k-r.a_grid_mat)*(1+m.r);
    
end