function [r] = KLS2_staticproblem_typeX(m,s,r)

%% Useful objects

%Assets
    r.a_grid_mat = r.a_grid'*ones(1,s.z_grid_size);         
    
%Productivity
    r.z_grid_mat = ones(s.a_grid_size,1)*r.z_grid'; 

%Indicator matrices for asset grid
    r.a_grid_mat_np = r.a_grid_mat<=0; %Non-positive assets
    r.a_grid_mat_p = r.a_grid_mat>0; %Positive assets 
         
%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)   
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag = m.Pf_lag/m.P_lag;
    FXdebt = m.lambda_X + (1-m.lambda_X)*((m.xi_real/m.xi_real_lag));  
      
%% Only Exporters
%Solution for every state (a,z)

%Aggregate constants
    Phi1_x = (m.sigma/(m.sigma-1))*m.w*(1/(1-m.alpha))*(m.Yf^(-1/m.sigma))*(m.tau_X^(1/(1-m.alpha)))/m.xi_real;        
    Phi2_x = ( (m.w*m.alpha/(1-m.alpha))^(m.alpha-1))/m.tau_X;  

%Capital        
    %Unconstrained choice of capital
        if (1+m.r)*FXdebt-(1-m.delta)>0
            r.k_u_x = ( (r.z_grid_mat.^((1-m.sigma)/m.sigma)) .* Phi1_x .* (Phi2_x.^((1/m.sigma)+(m.alpha/(1-m.alpha)))) .* ((1+m.r)*FXdebt-(1-m.delta))^((1-m.alpha)*((1/m.sigma)+m.alpha/(1-m.alpha))) ).^(-m.sigma);
        elseif (1+m.r)*FXdebt-(1-m.delta)<=0  
            r.k_u_x = 1e4*ones(size(r.a_grid_mat));
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

    %Indicator matrix for binding constraint
        r.const_x = r.k_x ~= r.k_u_x;

%Output values, given capital stock        
    r.yd_x = zeros(s.a_grid_size,s.z_grid_size);
    r.yf_x = ( Phi1_x.*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_x.^(m.alpha/(m.alpha-1))) ).^((m.sigma*(m.alpha-1))/((1-m.alpha)+m.alpha*m.sigma));
    r.y_x = (r.yd_x+m.tau_X*r.yf_x);        

%Labor, given output values and capital stock
    r.n_x = (r.y_x.^(1/(1-m.alpha))) .* (r.z_grid_mat.^(1/(m.alpha-1))) .* (r.k_x.^(m.alpha/(m.alpha-1)));       

%Prices
     r.pd_x =zeros(s.a_grid_size,s.z_grid_size);
     r.pf_x = ((r.yf_x/m.Yf).^(-1/m.sigma));                

%Profits (in units of the final good)
    r.pi_x = m.xi_real*r.yf_x.*r.pf_x - ((1+m.r)*FXdebt-(1-m.delta))*r.k_x - m.w*r.n_x;
    r.pi_x_wig = m.xi_real*r.yf_x.*r.pf_x - m.w*r.n_x;

%Lagrange multiplier and interest rate
%Test: Check whether expression equals zero when firms are unconstrained    
    lm_x_raw = (((m.w*m.alpha)/(1-m.alpha))*(r.y_x.^(1/(1-m.alpha))).*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_x.^(1/(m.alpha-1))) - (1+m.r)*FXdebt+(1-m.delta))/((1+m.r)*FXdebt-m.theta); 
    r.lm_x = r.const_x.*lm_x_raw;
          
%% Unconstrained profits: 
   
%Output:
    yd_u_x = zeros(s.a_grid_size,s.z_grid_size);
    yf_u_x = ( Phi1_x.*(r.z_grid_mat.^(1/(m.alpha-1))).*(r.k_u_x.^(m.alpha/(m.alpha-1))) ).^((m.sigma*(m.alpha-1))/((1-m.alpha)+m.alpha*m.sigma));
    r.y_u_x = (yd_u_x + m.tau_X*yf_u_x);        

%Labor, given output values and capital stock
    r.n_u_x = (r.y_u_x.^(1/(1-m.alpha))) .* (r.z_grid_mat.^(1/(m.alpha-1))) .* (r.k_u_x.^(m.alpha/(m.alpha-1)));       

% Prices
    pd_u_x =zeros(s.a_grid_size,s.z_grid_size);
    pf_u_x = (yf_u_x/m.Yf).^(-1/m.sigma);                

%Profits (in units of the final good)
    r.pi_u_x = m.xi_real*yf_u_x.*pf_u_x - ((1+m.r)*FXdebt-(1-m.delta))*r.k_u_x - m.w*r.n_u_x; 
    
% Additional variables needed for measures of constrained:
    r.k_u_nx = zeros(s.a_grid_size,s.z_grid_size);
    r.pi_u_nx = zeros(s.a_grid_size,s.z_grid_size);
    
end