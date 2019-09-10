function r = KLS2_staticproblem(m,s,r,varargin)
    
    if isempty(varargin) || varargin{1}==1 %High export cost firms          
        tau = m.tau;     
		if s.high_z_ForDebt == 0
			lambda=m.lambda;
		elseif s.high_z_ForDebt == 1
			lambda=m.lambda_mat;
		end
		
    elseif varargin{1}==2 %Low export cost firms                
        tau = m.tau_X;  
		if s.high_z_ForDebt == 0
			lambda=m.lambda_X;
		elseif s.high_z_ForDebt == 1
			lambda=m.lambda_X_mat;
		end		
		
    end

%% Useful objects

%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag= m.Pf_lag/m.P_lag; 
    FXdebt =lambda + (1-lambda).*((m.xi_real/m.xi_real_lag));     
    
    %% Exporters
    %Solution for every state (a,z), conditional on e'=1

    %Aggregate constants
    Phi0_x = 1 + (m.xi_real^m.sigma)*(tau^(1-m.sigma))*(m.Yf/m.Y);    
    Phi1_x = (m.sigma/(m.sigma-1))*m.w*(1/(1-m.alpha))*(m.Y^(-1/m.sigma))*(Phi0_x^(m.alpha/(1-m.alpha)));        
    Phi2_x = ((m.w*m.alpha/(1-m.alpha))^(m.alpha-1))*(1/Phi0_x);  

    %Capital     

	if s.high_z_ForDebt == 0	
		%Unconstrained choice of capital      
		if (1+m.r)*FXdebt-(1-m.delta)>0
			r.k_u_x = ( (r.z_grid_mat.^((1-m.sigma)/m.sigma)) .* Phi1_x .* (Phi2_x.^((1/m.sigma)+(m.alpha/(1-m.alpha)))) .* ...
				((1+m.r)*FXdebt-(1-m.delta))^((1-m.alpha)*((1/m.sigma)+m.alpha/(1-m.alpha))) ).^(-m.sigma);
		elseif (1+m.r)*FXdebt-(1-m.delta)<=0
			r.k_u_x = 1e4*ones(s.a_grid_size,s.z_grid_size); %i.e. we set it at constrained level (see NotesApr13.txt)
			r.infinitecapital = 1;
		end

		%Constrained choice of capital
		r.k_c_x = r.a_grid_mat*(1+(m.theta/((1+m.r)*FXdebt-m.theta)));

		%Capital choice
		r.k_x = zeros(s.a_grid_size,s.z_grid_size);
		if m.theta<(1+m.r)*FXdebt
			r.k_x = min(r.k_u_x,r.k_c_x);
		else
			r.k_x(r.a_grid_mat_p) = r.k_u_x(r.a_grid_mat_p);
			r.k_x(r.a_grid_mat_np) = max(r.k_u_x(r.a_grid_mat_np),r.k_c_x(r.a_grid_mat_np));
		end
		
	elseif s.high_z_ForDebt == 1
	
		%Unconstrained choice of capital        
		ind_temp = (1+m.r)*FXdebt-(1-m.delta)>0;
		r.k_u_x = zeros(s.a_grid_size,s.z_grid_size);
		r.k_u_x(ind_temp==1) = ( (r.z_grid_mat(ind_temp==1).^((1-m.sigma)/m.sigma)) .* Phi1_x .* (Phi2_x.^((1./m.sigma)+(m.alpha./(1-m.alpha)))) .* ...
											 ((1+m.r).*FXdebt(ind_temp==1)-(1-m.delta)).^((1-m.alpha)*((1./m.sigma)+m.alpha./(1-m.alpha))) ).^(-m.sigma);
		ind_mat = ones(s.a_grid_size,s.z_grid_size);    
		r.k_u_x(ind_temp==0) =  1e4*ind_mat(ind_temp==0);                              
		r.infinitecapital   = zeros (s.a_grid_size,s.z_grid_size);    
		r.infinitecapital(ind_temp==0) = 1;


		%Constrained choice of capital
		r.k_c_x = r.a_grid_mat.*(1+(m.theta./((1+m.r).*FXdebt-m.theta)));

		%Capital choice
		r.k_x = zeros(s.a_grid_size,s.z_grid_size);
		ind_temp = m.theta<(1+m.r)*FXdebt;
		r.k_x(ind_temp==1) = min( r.k_u_x(ind_temp==1),r.k_c_x(ind_temp==1)) ;
		r.k_x(ind_temp==0) = r.k_u_x(ind_temp==0);

	end	

    %Indicator matrix for binding constraint
    r.const_x = r.k_x ~= r.k_u_x;

    %Output values, given capital stock
    r.yd_x = ( Phi1_x.*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_x.^(m.alpha/(m.alpha-1))) ).^((m.sigma*(m.alpha-1))/((1-m.alpha)+m.alpha*m.sigma));
    r.yf_x = (m.xi_real^m.sigma)*((tau)^(-m.sigma))*(m.Yf/m.Y).*r.yd_x;
    r.y_x = (r.yd_x+tau*r.yf_x);        

    %Labor, given output values and capital stock (without F)
    r.n_x = (r.y_x.^(1/(1-m.alpha))) .* (r.z_grid_mat.^(1/(m.alpha-1))) .* (r.k_x.^(m.alpha/(m.alpha-1)));       

    %Prices
    r.pd_x = ((r.yd_x/m.Y).^(-1/m.sigma));
    r.pf_x = ((r.yf_x/m.Yf).^(-1/m.sigma));                

    %Profits (in units of the final good)
    r.pi_x = r.pd_x.*r.yd_x + m.xi_real*r.yf_x.*r.pf_x - ((1+m.r).*FXdebt-(1-m.delta)).*r.k_x - m.w*r.n_x;
	
    %Lagrange multiplier and interest rate
    %Test: Check whether expression equals zero when firms are unconstrained    
    lm_x_raw = (((m.w*m.alpha)/(1-m.alpha))*(r.y_x.^(1/(1-m.alpha))).*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_x.^(1/(m.alpha-1))) - (1+m.r).*FXdebt+(1-m.delta))./((1+m.r).*FXdebt-m.theta); 
    r.lm_x = r.const_x.*lm_x_raw;

    %% Non-exporters  
    %Solution for every state (a,z), conditional on e'=0

    %Aggregate constants
    Phi0_nx = 1;    
    Phi1_nx = (m.sigma/(m.sigma-1))*m.w*(1/(1-m.alpha))*(m.Y^(-1/m.sigma))*(Phi0_nx^(m.alpha/(1-m.alpha)));        
    Phi2_nx = ((m.w*m.alpha/(1-m.alpha))^(m.alpha-1))*(1/Phi0_nx);  

    %Capital   

	if s.high_z_ForDebt == 0		
		%Unconstrained choice of capital
		if (1+m.r)*FXdebt-(1-m.delta)>0
			r.k_u_nx = ( (r.z_grid_mat.^((1-m.sigma)/m.sigma)) .* Phi1_nx .* (Phi2_nx.^((1/m.sigma)+(m.alpha/(1-m.alpha)))) .* ... 
				((1+m.r)*FXdebt-(1-m.delta))^((1-m.alpha)*((1/m.sigma)+m.alpha/(1-m.alpha))) ).^(-m.sigma);
		elseif (1+m.r)*FXdebt-(1-m.delta)<=0                        
			r.k_u_nx = 1e4*ones(s.a_grid_size,s.z_grid_size); %i.e. we set it at constrained level (see NotesApr13.txt)
			r.infinitecapital = 1;        
		end

		%Constrained choice of capital
		r.k_c_nx = r.a_grid_mat*(1+(m.theta/((1+m.r)*FXdebt-m.theta)));

		%Capital choice
		r.k_nx = zeros(s.a_grid_size,s.z_grid_size);
		if m.theta<(1+m.r)*FXdebt
			r.k_nx = min(r.k_u_nx,r.k_c_nx);
		else
			r.k_nx(r.a_grid_mat_p) = r.k_u_nx(r.a_grid_mat_p);
			r.k_nx(r.a_grid_mat_np) = max(r.k_u_nx(r.a_grid_mat_np),r.k_c_nx(r.a_grid_mat_np));
		end
		
	elseif s.high_z_ForDebt == 1
	
		%Unconstrained choice of capital
		ind_temp = (1+m.r)*FXdebt-(1-m.delta)>0;
		r.k_u_nx = zeros(s.a_grid_size,s.z_grid_size);
		r.k_u_nx(ind_temp==1) = ( (r.z_grid_mat(ind_temp==1).^((1-m.sigma)./m.sigma)) .* Phi1_nx .* (Phi2_nx.^((1./m.sigma)+(m.alpha./(1-m.alpha)))) .* ... 
												  ((1+m.r).*FXdebt(ind_temp==1)-(1-m.delta)).^((1-m.alpha)*((1/m.sigma)+m.alpha./(1-m.alpha))) ).^(-m.sigma);    
		ind_mat = ones(s.a_grid_size,s.z_grid_size);                                      
		r.k_u_nx(ind_temp==0) = 1e4*ind_mat(ind_temp==0);

		%Constrained choice of capital
		r.k_c_nx = r.a_grid_mat.*(1+(m.theta./((1+m.r).*FXdebt-m.theta)));

		%Capital choice
		r.k_nx = zeros(s.a_grid_size,s.z_grid_size);
		ind_temp = m.theta<(1+m.r).*FXdebt;
		r.k_nx(ind_temp==1) = min(r.k_u_nx(ind_temp==1),r.k_c_nx(ind_temp==1));
		r.k_nx(ind_temp==0) = r.k_u_nx(ind_temp==0);

	end	


    %Indicator matrix for binding constraint            
    r.const_nx = r.k_nx ~= r.k_u_nx;     

    %Output values, given capital stock
    r.yd_nx = ( Phi1_nx.*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_nx.^(m.alpha/(m.alpha-1))) ).^((m.sigma*(m.alpha-1))/((1-m.alpha)+m.alpha*m.sigma));
    r.yf_nx = zeros(s.a_grid_size,s.z_grid_size);
    r.y_nx = (r.yd_nx+tau*r.yf_nx);        

    %Labor, given output values and capital stock
    r.n_nx = (r.y_nx.^(1/(1-m.alpha))) .* (r.z_grid_mat.^(1/(m.alpha-1))) .* (r.k_nx.^(m.alpha/(m.alpha-1)));

    %Prices
    r.pd_nx = ((r.yd_nx/m.Y).^(-1/m.sigma));
    r.pf_nx = zeros(s.a_grid_size,s.z_grid_size);                     

    %Profits
    r.pi_nx = r.pd_nx.*r.yd_nx - ((1+m.r).*FXdebt-(1-m.delta)).*r.k_nx - m.w*r.n_nx;
	
    %Lagrange multiplier and interest rate
    %Test: Check whether expression equals zero when firms are unconstrained         
    lm_nx_raw = (((m.w*m.alpha)/(1-m.alpha))*(r.y_nx.^(1/(1-m.alpha))).*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_nx.^(1/(m.alpha-1))) - (1+m.r).*FXdebt+(1-m.delta))./((1+m.r).*FXdebt-m.theta); 
    r.lm_nx = r.const_nx.*lm_nx_raw;                   

%% Unconstrained profits

%EXPORTERS

    %Output values, given capital stock
    yd_u_x = (Phi1_x.*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_u_x.^(m.alpha/(m.alpha-1))) ).^((m.sigma*(m.alpha-1))/((1-m.alpha)+m.alpha*m.sigma));
    yf_u_x = (m.xi_real^m.sigma)*((tau)^(-m.sigma))*(m.Yf/m.Y).*yd_u_x;
    r.y_u_x = (yd_u_x+tau*yf_u_x);        

    %Labor, given output values and capital stock
    r.n_u_x = (r.y_u_x.^(1/(1-m.alpha))) .* (r.z_grid_mat.^(1/(m.alpha-1))) .* (r.k_u_x.^(m.alpha/(m.alpha-1)));       

    %Prices
    pd_u_x = ((yd_u_x/m.Y).^(-1/m.sigma));
    pf_u_x = ((yf_u_x/m.Yf).^(-1/m.sigma));                

    %Profits
    r.pi_u_x = pd_u_x.*yd_u_x + m.xi_real*pf_u_x.*yf_u_x - ((1+m.r).*FXdebt-(1-m.delta)).*r.k_u_x - m.w*r.n_u_x;
	

%NON-EXPORTERS
    %Output values, given capital stock
    yd_u_nx = ( Phi1_nx.*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_u_nx.^(m.alpha/(m.alpha-1))) ).^((m.sigma*(m.alpha-1))/((1-m.alpha)+m.alpha*m.sigma));
    yf_u_nx = zeros(s.a_grid_size,s.z_grid_size);
    r.y_u_nx = (yd_u_nx+tau*yf_u_nx);        

    %Labor, given output values and capital stock
    r.n_u_nx = (r.y_u_nx.^(1/(1-m.alpha))) .* (r.z_grid_mat.^(1/(m.alpha-1))) .* (r.k_u_nx.^(m.alpha/(m.alpha-1)));

    %Prices
    pd_u_nx = ((yd_u_nx/m.Y).^(-1/m.sigma));
    pf_u_nx = zeros(s.a_grid_size,s.z_grid_size);                     

    %Profits
    r.pi_u_nx = pd_u_nx.*yd_u_nx - ((1+m.r).*FXdebt-(1-m.delta)).*r.k_u_nx - m.w*r.n_u_nx;  
	    
    
end