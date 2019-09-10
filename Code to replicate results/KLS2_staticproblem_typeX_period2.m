function [r] = KLS2_staticproblem_typeX_period2(m,s,r,rt)

%% Useful objects

         
%Real exchange rate (only appears in FXdebt as m.xi_real/m.xi_lag_real)
    m.xi_real = m.Pf/m.P;
    m.xi_real_lag=m.Pf_lag/m.P_lag; 
    FXdebt =m.lambda_period_2_X + (1-m.lambda_period_2_X)*((m.xi_real/m.xi_real_lag));        
    
%% Only Exporters
%Solution for every state (a,e,z)

%Capital stock and labor
    r.k_x = rt{1}.k_x;   


%Output values
    r.yf_x = (((m.sigma/(m.sigma-1))*m.w/(1-m.alpha)).^((m.sigma*(m.alpha-1))/(1-m.alpha+m.alpha*m.sigma))) .* ...
             (m.tau_X^(-m.sigma/(1-m.alpha+m.alpha*m.sigma))) .* ...
             (r.z_grid_mat.^(m.sigma/(1-m.alpha+m.alpha*m.sigma))) .* ...
             (r.k_x.^((m.sigma*m.alpha)/(1-m.alpha+m.alpha*m.sigma))) .* ...
             m.Yf.^((1-m.alpha)/(1-m.alpha+m.alpha*m.sigma)) .* ...
             (m.xi_real.^((m.sigma*(1-m.alpha))/(1-m.alpha+m.alpha*m.sigma)));         
    r.yd_x = zeros(s.a_grid_size,s.z_grid_size); 
    r.y_x = (r.yd_x+m.tau_X*r.yf_x);        

    r.n_x = ((m.tau_X*r.yf_x)./(r.z_grid_mat.*(r.k_x.^m.alpha))).^(1/(1-m.alpha));  
    
%Prices
    r.pd_x = zeros(size(r.yf_x));
    r.pf_x = (r.yf_x/m.Yf).^(-1/m.sigma);                
    r.pf_x(r.yf_x==0) = 0;

%Profits
    r.pi_x = m.xi_real*r.yf_x.*r.pf_x - ((1+m.r)*FXdebt-(1-m.delta))*r.k_x - m.w*r.n_x;
    r.pi_x_wig = m.xi_real*r.pf_x.*r.yf_x - m.w*r.n_x;    

%Lagrange multiplier and interest rate
%This is not working for period 2:
    r.lm_x = zeros(s.a_grid_size,s.z_grid_size); 
            

end