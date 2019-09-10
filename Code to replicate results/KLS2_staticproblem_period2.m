%This code solves the static problem for period 2 of the transition
%Here the capital stock is taken as given, as it was chosen in period 2

function r = KLS2_staticproblem_period2(m,s,r,rt,varargin)
    

    if isempty(varargin) || varargin{1}==1
          
        tau=m.tau;
        
		if s.high_z_ForDebt == 0
			lambda_period_2=m.lambda_period_2;
		elseif s.high_z_ForDebt == 1
			lambda_period_2=m.lambda_period_2_mat;
		end
        
    elseif varargin{1}==2  %Type 2 firms
                
        tau=m.tau_X;
        
		if s.high_z_ForDebt == 0
			lambda_period_2=m.lambda_period_2_X;
		elseif s.high_z_ForDebt == 1
			lambda_period_2=m.lambda_X_period_2_mat;
		end
        
    end

%% Useful objects
           
    %Impact of changes in exchange rate on debt returns (Real exchange rate only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag=m.Pf_lag/m.P_lag; 
    FXdebt =lambda_period_2 + (1-lambda_period_2).*((m.xi_real/m.xi_real_lag));        
    	
	
%% Capital as a state variable:
r.k = rt{1}.k;
%r.k2 = (1-rt{1}.e).*rt{1}.k_nx + rt{1}.e.*rt{1}.k_x;
  
%% Exporters
%Solution for every state (a,e,z), conditional on e'=1

%Solution objects
    Phi = (1+ (1/tau^(m.sigma-1))*(m.xi_real^m.sigma)*(m.Yf/m.Y));

%Phi
    Phi_x = ((1-m.alpha)*((m.sigma-1)/m.sigma)*(1/m.w)*(m.Y^(1/m.sigma)))^(m.sigma/(m.alpha*m.sigma+1-m.alpha));

%Capital stock 
    r.k_x = rt{1}.k;       
    
%Labor    
    r.n_x = (Phi^(1/(m.alpha*m.sigma+1-m.alpha)))*Phi_x*(r.z_grid_mat.*(r.k.^(m.alpha))).^((m.sigma-1)/(m.alpha*m.sigma+1-m.alpha));  

%Output values are fixed
    r.yd_x = ((1/Phi)^((m.alpha*m.sigma)/(m.alpha*m.sigma+1-m.alpha))) .* (Phi_x.^(1-m.alpha)) .* (r.z_grid_mat.*(r.k.^(m.alpha))).^(m.sigma/(m.alpha*m.sigma+1-m.alpha));   
    r.yf_x = ((1/tau^(m.sigma))*(m.xi_real^m.sigma)*(m.Yf/m.Y)) ./ (Phi^((m.alpha*m.sigma)/(m.alpha*m.sigma+1-m.alpha))) ...
             .* (Phi_x.^(1-m.alpha)) .* (r.z_grid_mat.*(r.k.^(m.alpha))).^(m.sigma/(m.alpha*m.sigma+1-m.alpha));   
    r.y_x  = (r.yd_x+tau*r.yf_x);        

%Prices adjust (change in demand)
    r.pd_x = (r.yd_x/m.Y).^(-1/m.sigma);
    r.pf_x = (r.yf_x/m.Yf).^(-1/m.sigma);                
    r.pf_x(r.yf_x==0) = 0;

%Profits
    r.pi_x = r.pd_x.*r.yd_x + m.xi_real*r.pf_x.*r.yf_x - ((1+m.r).*FXdebt-(1-m.delta)).*r.k_x - m.w*r.n_x;
	r.pi_x_wig = r.pd_x.*r.yd_x + m.xi_real*r.pf_x.*r.yf_x - m.w*r.n_x;

%Lagrange multiplier and interest rate
    r.lm_x = rt{1}.lm_x;  
        
%% Non-exporters  
%Solution for every state (a,e,z), conditional on e'=0

%Phi
    Phi_nx = ((1-m.alpha)*((m.sigma-1)/m.sigma)*(1/m.w)*(m.Y^(1/m.sigma)))^(m.sigma/(m.alpha*m.sigma+1-m.alpha));

%Capital stock 
    r.k_nx = rt{1}.k;     
    
%Labor    
    r.n_nx = Phi_nx*(r.z_grid_mat.*(r.k.^(m.alpha))).^((m.sigma-1)/(m.alpha*m.sigma+1-m.alpha)); 

%Output values, given capital stock
    r.yd_nx = (Phi_nx.^(1-m.alpha))*(r.z_grid_mat.*(r.k.^(m.alpha))).^(m.sigma/(m.alpha*m.sigma+1-m.alpha));        
    r.yf_nx = zeros(s.a_grid_size,s.z_grid_size);
    r.y_nx = (r.yd_nx+tau*r.yf_nx);        

%Prices
    r.pd_nx = ((r.yd_nx/m.Y).^(-1/m.sigma));
    r.pf_nx = zeros(s.a_grid_size,s.z_grid_size);                     

%Profits
    r.pi_nx = r.pd_nx.*r.yd_nx - ((1+m.r).*FXdebt-(1-m.delta)).*r.k_nx - m.w*r.n_nx;
	r.pi_nx_wig = r.pd_nx.*r.yd_nx - m.w*r.n_nx;

%Lagrange multiplier and interest rate
    r.lm_nx = rt{1}.lm_nx;   
        
%% Unconstrained profits: 

%Productivity
   % r.z_grid_mat = ones(s.a_grid_size,1)*r.z_grid'; 


% The issue is that these exporters should have chosen an unconstrained
% amount of capital for today yesterday when they did not know the
% change in m.r.
       
%EXPORTERS:   
%     %Constants:
%     Phi0_x = 1 + (m.xi_real_lag^m.sigma)*(tau^(1-m.sigma))*(m.Yf/m.Y_lag);    
%     Phi1_x = (m.sigma/(m.sigma-1))*m.w_lag*(1/(1-m.alpha))*(m.Y_lag^(-1/m.sigma))*(Phi0_x^(m.alpha/(1-m.alpha)));        
%     Phi2_x = ((m.w_lag*m.alpha/(1-m.alpha))^(m.alpha-1))*(1/Phi0_x);  
%         
%     %Unconstrained Capital (chosen for old interest rate):
%     r.k_u_x = ( (r.z_grid_mat.^((1-m.sigma)/m.sigma)) .* Phi1_x .* (Phi2_x.^((1/m.sigma)+(m.alpha/(1-m.alpha)))) .* ((1+m.r_lag)*1-(1-m.delta))^((1-m.alpha)*((1/m.sigma)+m.alpha/(1-m.alpha))) ).^(-m.sigma);
%       
%     %Choice of Labor   
%     n_u_x = (Phi2.^(1/(m.alpha*m.sigma+(1-m.alpha)))).*Phi.*(r.z_grid_mat.*(r.k_u_x.^(m.alpha))).^((m.sigma-1)/(m.alpha*m.sigma+1-m.alpha));          
% 
%     %Output values
%    yd_u_x = (1/(Phi2^((m.alpha*m.sigma)/(m.alpha*m.sigma+(1-m.alpha))))) * (Phi^(1-m.alpha)) .* (r.z_grid_mat.*(r.k_u_x.^(m.alpha))).^(m.sigma/(m.alpha*m.sigma+1-m.alpha));      
%    yf_u_x = 1/tau*(Phi2-1)*yd_u_x;      
% 
%     %Prices
%     pd_u_x = (yd_u_x/m.Y).^(-1/m.sigma);
%     pf_u_x = (yf_u_x/m.Yf).^(-1/m.sigma);                
%         
%     %Profits
%     r.pi_u_x = pd_u_x.*yd_u_x + m.xi_real*pf_u_x.*yf_u_x - ((1+m.r)*FXdebt-(1-m.delta))*r.k_u_x - m.w*n_u_x;        
%         
% %NON-EXPORTERS:    
%     %Constants:
%     Phi0_nx = 1;    
%     Phi1_nx = (m.sigma/(m.sigma-1))*m.w_lag*(1/(1-m.alpha))*(m.Y_lag^(-1/m.sigma))*(Phi0_nx^(m.alpha/(1-m.alpha)));        
%     Phi2_nx = ((m.w_lag*m.alpha/(1-m.alpha))^(m.alpha-1))*(1/Phi0_nx);  
%     
%     %Unconstrained Capital (chosen for old interest rate):   
%     r.k_u_nx = ( (r.z_grid_mat.^((1-m.sigma)/m.sigma)) .* Phi1_nx .* (Phi2_nx.^((1/m.sigma)+(m.alpha/(1-m.alpha)))) .* ((1+m.r_lag)*1-(1-m.delta))^((1-m.alpha)*((1/m.sigma)+m.alpha/(1-m.alpha))) ).^(-m.sigma);
%           
%     %Choice of Labor     
%     n_u_nx = Phi*(r.z_grid_mat.*(r.k_u_nx.^(m.alpha))).^((m.sigma-1)/(m.alpha*m.sigma+1-m.alpha));
%                 
%     %Output values, given capital stock
%     yd_u_nx = r.z_grid_mat.*(r.k_u_nx.^m.alpha).*(n_u_nx.^(1-m.alpha));   
%     yf_u_nx = zeros(2,s.a_grid_size,s.z_grid_size);    
% 
%     %Prices
%     pd_u_nx = ((yd_u_nx/m.Y).^(-1/m.sigma));
%     pf_u_nx = zeros(2,s.a_grid_size,s.z_grid_size);                     
%         
%     %Profits
%     r.pi_u_nx = pd_u_nx.*yd_u_nx - ((1+m.r)*FXdebt-(1-m.delta))*r.k_u_nx - m.w*n_u_nx;  
%                     
end