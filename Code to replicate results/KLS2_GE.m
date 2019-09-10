function [mcc,m,r,r_X,s,sim] = KLS2_GE(x,m,s,r)

%% Aggregate prices and quantities

% Guessed prices           
    m.w = exp(x(1));        
    m.Y = exp(x(2));    
    m.P = exp(x(3));     

% Lagged prices
% Steady state: Prices are constant
    m.Pf_lag = m.Pf;
    m.P_lag  = m.P ;
      
%% Solution        
        
% Display
    disp('------------------------------------------------');
    disp(['Guesses: Y=' num2str(m.Y,'%10.8g') ', w=' num2str(m.w,'%10.8g') ', r=' num2str(m.r,'%10.8g') ', P=' num2str(m.P,'%10.8g')]);                

% Static problem    
    % Adjust aggregate output
    m.Y = m.Y*m.A^(m.sigma-1);
    
    % Initialization
    r_X = r;
    
    % High export cost firms
    r = KLS2_staticproblem(m,s,r);  
        
    % Low export cost firms
    if s.model==1 || s.model == 2 || s.model == 4 %Model with two types of firms
	
		r_X = KLS2_staticproblem(m,s,r_X,2); 
    
	elseif s.model==3 %Baseline model, with one type of firm that only exports
        r_X = KLS2_staticproblem_typeX(m,s,r_X);
    end

    % Restore aggregate output
    m.Y = m.Y/(m.A^(m.sigma-1));

% Dynamic problem 
    %High export cost firms
    r = KLS2_dynamicproblem(m,s,r);        
        
    %Low export cost firms
	if s.model==1  || s.model == 2 || s.model == 4 %Model with two types of firms
		r_X = KLS2_dynamicproblem(m,s,r_X,2);
	elseif s.model==3 %Baseline model, with one type of firm that only exports  
		r_X = KLS2_dynamicproblem_typeX(m,s,r_X);
	end
	
% Simulation    
    if s.flag_simulate == 0
        [sim,r,r_X,s] = KLS2_simulate(m,s,r,r_X);                  
    elseif s.flag_simulate == 1                
        [sim,r,r_X,s] = KLS2_simulate_shocks(m,s,r,r_X);                 
    end

% Market clearing conditions    
    mcc = [sim.mc_n sim.mc_y sim.mc_y_belief]; 
    sim.mcc = mcc;     

% Display
if s.display==1
    disp(['GE: Y_MCC=' num2str(sim.mc_y) ', N_MCC=' num2str(sim.mc_n) ', A_MCC=' num2str(sim.mc_a) ', P_CES=' num2str(sim.mc_p_belief) ', Y_CES=' num2str(sim.mc_y_belief)]);
    disp(' ');
    disp(['X/GDP: ' num2str(sim.X_GDP)]);  
    disp(['X_X/GDP: ' num2str(sim.X_X_GDP)]);
    disp(['Avg exp_int type 1: ' num2str(sim.x_share_avg_type1)]);
    disp(['Avg exp_int type 2: ' num2str(sim.x_share_avg_type2)]);
    disp(['Share of exporters type 2: ' num2str(sim.share_x_composition)]);                         
    disp(['Share of exporters: ' num2str(sim.share_x)]);                     
    disp(['Credit/GDP: ' num2str(sim.credit_gdp)]);
    disp(['Sd of log sales: ' num2str(sim.ln_sales_sd)]);
    disp(['Share of sales top 25: ' num2str(sim.sales_share_top25)]);
    disp(['Exporter sales premium: ' num2str(sim.xpremium_sales)]); 
    disp(['Exporter labor premium: ' num2str(sim.xpremium_labor)]); 
    disp(['Exporter med sales premium: ' num2str(sim.xpremium_sales_med)]); 
    disp(['Type2 vs Type1 sales premium: ' num2str(sim.type2premium_sales)]); 
    disp(['Type2 vs Type1 sales premium 2: ' num2str(sim.type2premium_sales_2)]);
    disp(['Type2 vs Type1 sales_x premium: ' num2str(sim.type2premium_salesx)]); 
    disp(['Type2 vs Type1 sales_d premium: ' num2str(sim.type2premium_salesd)]); 
    disp(['Type2 vs Type1 sales_d premium 2: ' num2str(sim.type2premium_salesd_2)]);
end          
    
end