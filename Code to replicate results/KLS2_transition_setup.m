%Period 1 = Initial steady state
%Period 2 = Shock realized, full path of exchange rate becomes known, transition begins
%Period N = Final steady state

%Length of the transition
trans.N = 25; 

%Shock to collateral constraint
theta_old = m.theta;
theta_new = m.theta; 
m.theta_v = ones(trans.N,1); 
m.theta_v(1:2) = theta_old;
m.theta_v(3:end) = theta_new;  

if s.financial_crisis==1
   m.theta_v(3:6) = [0.27558498	0.22725906	0.30416216	0.42243861];
   m.theta_v(7:end) = m.theta_v(6);
end


%Real interest rate shock
r_old = m.r;
r_new = m.r; 
m.rv = zeros(trans.N,1); 
m.rv(1:2) = r_old;

if s.model==1
    m.rv(3:6) = [0.15163496     0.18177558     0.11478128    0.097790579];
elseif s.model==2
    m.rv(3:6) = [0.099033836    0.091638678    0.094480253    0.092716351];
   
   if s.FFshocks == 1 %In the workspace for baseline, this is r_tilde
        m.rv(3:6) = [0.098331621	0.098664491	0.101663141	0.099531417];
   end

elseif s.model==3
    m.rv(3:6) = [0.15446371     0.17665159     0.10696911    0.085073114];
elseif s.model==4
    m.rv(3:6) = [0.14681557     0.17725749     0.11304378      0.1000285];
end

if s.irp==1 && s.model==1
    m.rv(3:6) = [0.08097009    0.081485339    0.087172864    0.087531566];
end

m.rv(7:end) = m.rv(6);
    
if s.irf_pm==1 || s.financial_crisis==1
   m.rv(3:end) = m.r;
end






% Beta shock
beta_old = m.beta;
beta_new = m.beta; 
m.beta_v = zeros(trans.N,1); 
m.beta_v(1) = beta_old;
m.beta_v(2:end) = beta_new;

% delta shock
delta_old = m.delta;
delta_new = m.delta; 
m.delta_v = zeros(trans.N,1); 
m.delta_v(1) = delta_old;
m.delta_v(2:end) = delta_new;

%Foreign CPI shock
Pf_old = m.Pf;
Pf_new = m.Pf; 
m.Pfv = zeros(trans.N,1); 
m.Pfv(1) = Pf_old; 
m.Pfv(2) = Pf_old + (Pf_new-Pf_old)/2;
m.Pfv(3:end) = Pf_new;

%Shock to pm 
if s.model==1
    pm_path = [1   0.55493726     0.61656867     0.75728231     0.78231349];
elseif s.model==2
    pm_path = [1   0.51496959     0.60601235      0.7631675     0.78336051];
    
    if s.FFshocks == 1
        pm_path = [1   0.55493726     0.61656867     0.75728231     0.78231349];
    end
elseif s.model==3
    pm_path = [1   0.58388002     0.61894929     0.74326028      0.7763353];
elseif s.model==4
    pm_path = [1   0.54728636     0.61840267      0.7608604     0.78255951];

end

if s.irp==1 && s.model==1
    pm_path  = [1   0.5545584     0.64208869     0.77585961     0.78848017];
end

if s.financial_crisis==1
    
     pm_path = [1	0.64336679	0.68235635	0.73848649	0.76185425];	
end

if s.irf_pm==1
   pm_path = [1 1/1.8 1/1.8 1/1.8 1/1.8]; 
end

pm_old = m.Pm;
pm_new = m.Pm;
m.pm_v = zeros(trans.N,1); 
m.pm_v(1) = pm_old;    
m.pm_v(2) = pm_path(2); 
m.pm_v(3) = pm_path(3); 
m.pm_v(4) = pm_path(4); 
m.pm_v(5) = pm_path(5); 
m.pm_v(6:end) = pm_path(5);
    
% Foreign output shock
% Negative shock -> devaluation
% Positive shock -> apreciation
Yf_old = m.Yf;
Yf_new = m.Yf;
m.Yfv = zeros(trans.N,1); 
m.Yfv(1) = Yf_old;
m.Yfv(2) = Yf_old + (Yf_new-Yf_old)/2;
m.Yfv(3:end) = Yf_new;    

%Share of debt denominated in domestic final goods
m.lambda_period_2= m.lambda;
m.lambda_period_2_X = m.lambda_X;  

if s.UIPdeviations == 1

	m.lambda_old = m.lambda;
	m.lambda_new = 0;  
	m.lambda_v = zeros(trans.N,1); 
	m.lambda_v(1) = m.lambda_old;    
	m.lambda_v(2) = m.lambda_old; 
	m.lambda_period_2 = m.lambda_old; 
	m.lambda_v(3:end) = m.lambda_new;

	m.lambda_X_old = m.lambda;
	m.lambda_X_new = 0;  
	m.lambda_X_v = zeros(trans.N,1); 
	m.lambda_X_v(1) = m.lambda_X_old;    
	m.lambda_X_v(2) = m.lambda_X_old; 
    m.lambda_period_2_X = m.lambda_X_old; 
	m.lambda_X_v(3:end) = m.lambda_X_new;

end
    

%Shock to labor income 
%Labor income = (1-zeta_w)*w
zeta_w_old = m.zeta_w;
zeta_w_new = m.zeta_w;
m.zeta_w_v = zeros(trans.N,1); 
m.zeta_w_v(1) = zeta_w_old;
m.zeta_w_v(2:end) = zeta_w_new; 

%Shock to dividend income
%Dividend income = (1-zeta_pi)*pi
zeta_pi_old = m.zeta_pi;
zeta_pi_new = m.zeta_pi;
m.zeta_pi_v = zeros(trans.N,1); 
m.zeta_pi_v(1) = zeta_pi_old;
m.zeta_pi_v(2:end) = zeta_pi_new; 

%Shock to productivity of final good producer
A_old = m.A;
A_new = m.A;
m.A_v = ones(trans.N,1); 
m.A_v(1) = A_old;
m.A_v(2:end) = A_new; 


if s.model==1
    zagg_path = [1          0.97730792     0.99095589      1.0119136      1.0325626];
elseif s.model==2
    zagg_path = [1          0.98132361     0.99772036      1.0152815      1.0437912];
 
    if s.FFshocks == 1
       
        zagg_path = [1          0.97730792     0.99095589      1.0119136      1.0325626];
    end
   
elseif s.model==3
    zagg_path = [1          0.9508588     0.97499091     0.99634435      1.0254037];
elseif s.model==4
    zagg_path = [1          0.9808992     0.99527242      1.0151195      1.0353159];

end

if s.irp==1 && s.model==1
    zagg_path  = [1          0.88479301     0.96040402     0.99564284      1.0319007];
end


if s.financial_crisis==1
    
    zagg_path = [ 1     0.9270488	0.98180721	1.038341	1.0465126];
    
end

if s.irf_pm==1
    zagg_path = [1           1     1      1      1];
end




zagg_old = m.zagg;
zagg_new = m.zagg;
m.zagg_v = ones(trans.N,1); 
m.zagg_v(1) = zagg_old;
m.zagg_v(2:5) = zagg_path(2:5);     
m.zagg_v(6:end) = zagg_path(5);    

