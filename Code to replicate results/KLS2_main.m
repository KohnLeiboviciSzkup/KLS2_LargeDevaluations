%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to replicate results in 
% "Financial Frictions and Export Dynamics in Large Devaluations"
% Journal of International Economics, 2019
% David Kohn, Fernando Leibovici and Michal Szkup 
% August 2019
%
% ---- Main file  ------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Start

% Clear memory and workspace, close all windows
    clear;
    close all;
    clc;

% Timer
    tic;

% Directory settings
    results.dir = pwd;
    results.dir_results = strcat(results.dir,'\results\');

%% Control panel

% Model choice
    s.model = 1; %=1 if baseline model with two types of firms
                 %=2 if frictionless model with two types of firms
                 %=3 if baseline model with one type of firms that only exports
                 %=4 if baseline model with one type of firms

% pm shocks
    s.irf_pm = 0; %=0 if baseline experiment 
	              %=1 if pm shock only, impulse response (use s.model=1)

    s.lambda_flag = 0; %=0 if baseline calibration
                       %=1 if lambda=0,lambdaX=0
                       %=2 if lambda=1,lambdaX=1
                       %=3 if lambda=0.5,lambdaX=0

% Alternative experiments 				  
    s.UIPdeviations = 0;%=0 if baseline experiment				   
	                    %=1 if UIP deviations experiment (use s.model=1)

    s.high_z_ForDebt = 0;%=0 if baseline experiment				   
	                     %=1 if high z firms hold high foreign debt experiment (use s.model=1)						
  
    s.financial_crisis = 0;%=0 if baseline experiment				   
    	                   %=1 if financial crisis shock (theta shock)
  
    s.irp = 0; %=0 if baseline experiment				   
    	       %=1 if IRP experiment
    
    s.FFshocks = 0; %=0 if baseline experiment				   
                    %=1 if Shocks of baseline experiment (with s.model==2)
               
% Saves workspace after running code
    s.save_workspace = 0;

% Save prices
    s.save_prices = 0;
    
                              
% Steady-state in general equilibrium?
    s.GE = 1; %=0 if partial equilibrium
              %=1 if general equilibrium        

% Transition
    s.transition = 1;
    s.transition_GE = 1;    

% Load prices    
    s.load_prices = 1;
    
% Denomination of fixed costs (includes sunk costs)
    s.fcost_fgoods = 0; 
        %=0 if fixed costs denominated in units of labor
        %=1 if fixed costs denominated in units of the final good

    
% Simulation
    s.flag_simulate = 0; % s.flag_simulate = 0, simulates by iterating on measure of firms,
                         % s.flag_simulate = 1, simulates by generating sequence of random shocks
 
% Display
    s.display=1; %Display output on every GE iteration
              

    
%% Calibration
    
% Calibrated parameters
    if s.model==1 %Baseline model with two types of firms

        m.Xshare = 0.044243; %Fraction of firms that only exports                                     
        m.tau = 5.7089; %Iceberg trade cost 
        m.tau_X  = 1.7578; %Trade cost of firms that only export
        m.theta = 0.49375; %Collateral constraint
        F=0.043183; %Fixed export cost
        m.F_base = F;
        m.F_X_base = F; 
        m.log_z_sigma = 0.25693; %Standard deviation of log-normal productivity distribution                             
        m.log_z_rho = 0.87598; %Persistence productivity 
        m.lambda = 0.45;
        m.lambda_X = 0.45;
        m.beta = 0.85404 ; %Discount factor 
                
        
        
    elseif s.model==2 % Fricitonless model with two types of firms    
        
        %Frictionless
        m.Xshare = 0.043508; %Fraction of firms that only exports                                     
        m.tau = 5.5844; %Iceberg trade cost 
        m.tau_X  = 1.7123; %Trade cost of firms that only export
        m.theta = 1e8; %Collateral constraint
        F=0.041353; %Fixed export cost
        m.F_base = F;
        m.F_X_base = F; 
        m.log_z_sigma = 0.25198; %Standard deviation of log-normal productivity distribution                          
        m.log_z_rho = 0.85726; %Persistence productivity   
        m.lambda = 1; %0.55;
        m.lambda_X = 1; %0.55;
        m.beta = 0.89392; %Discount factor        
    
    elseif s.model==3 % Baseline model with one type of firms that only exports
      
        %Financial frictions
        m.Xshare = 0.32; %Fraction of firms that only exports
        m.tau = 100000; %Iceberg trade cost 
        m.tau_X  = 4.3096; %Trade cost of firms that only export        
        m.theta =  0.46827; %0.363478667137117; %Collateral constraint
        m.F_base = 100; %Fixed export cost for type 1 firms
        m.F_X_base = 0; %Fixed export cost for type 2 firms
        m.log_z_sigma = 0.25546; %Standard deviation of log-normal productivity distribution   
        m.log_z_rho = 0.88509; %Persistence productivity 
        m.lambda = 0.45;
        m.lambda_X = 0.45;
        m.beta = 0.84823  ; %Discount factor       
        
    elseif s.model==4 % Baseline model with one type of firms 
        
        m.Xshare = 0; %Fraction of firms that only exports                                     
        m.tau = 4.4371; %Iceberg trade cost 
        m.tau_X  = 4.4371; %Trade cost of firms that only export
        m.theta = 0.50228; %Collateral constraint
        F=0.039384; %Fixed export cost
        m.F_base = F;
        m.F_X_base = F; 
        m.log_z_sigma = 0.26299; %Standard deviation of log-normal productivity distribution 
        m.log_z_rho = 0.9159; %Persistence productivity 
        m.lambda = 0.45;
        m.lambda_X = 0.45;
        m.beta = 0.83507  ; %Discount factor  
        
    end
    
    % Overwrite values of lambda
    if s.lambda_flag==1
       m.lambda = 0;
       m.lambda_X = 0;
    elseif s.lambda_flag==2
       m.lambda = 1;
       m.lambda_X = 1;            
    elseif s.lambda_flag==3    
       m.lambda = 0.5;
       m.lambda_X = 0;              
    end
	
	if s.high_z_ForDebt == 1
	   m.lambda = 1;
       m.lambda_X = 1;            
	end
    
    
% Pre-assigned
    
    m.r = 0.08; %Interest rate
    m.sigma = 4; %Elasticity of substitution across varieties
    m.gamma = 2; %Risk aversion
    m.alpha = 0.33; %Capital share in production 
    m.delta = 0.06; %Capital depreciation rate               
    m.z_mean = 1; %Average productivity level of log-normal distribution
    
% Additional parameters    
    m.zeta_w = 0;
    m.zeta_pi = 0;
    m.A = 1;
    m.zagg = 1;
    
% Rest of the world
    m.Yf = 10; %Output in the rest of the world                                
    m.Pf = 1; %Price index in the rest of the world 
    m.Pm = 1; %Price index of goods imported from the rest of the world
    m.omega_m = 1; %Measure of varieties produced by the rest of the world 
    
% Productivity mean  
    m.log_z_mu = log(m.z_mean)-(m.log_z_sigma^2)*(1/((1-m.log_z_rho^2)))*(1/2); %Normalize average productivity to 1, log_z_mu is the mean of log_z

%% Setup shocks to parameters
    KLS2_transition_setup;   
    
%% Solution options

% Productivity 
    s.z_grid_size = 100; %Productivity grid size
    s.z_grid_power = 1; %Curvature parameter to control distance across grid points

% Assets 
    s.a_grid_size = 250; %Asset grid size
    s.a_grid_size_negative = 0; %Number of negative asset values, if assets can be negative

    s.a_grid_power_pos = 2; %Curvature parameter to control distance across grid points -- for a>=0 (for a<0, grid spacing is linear)
    s.a_grid_power_neg = 1; 

    s.a_grid_ub = 200; %Upper bound on asset grid, CHECK WHETHER IT BINDS!            
    s.a_grid_lb_pos = 1e-3;         
    s.a_grid_ub_neg = -1e-3;        
    s.a_grid_lb_neg = -25;           

% Convergence parameters    
    s.eps = 1e-8; 
    s.eps_sm = 1e-8; %Stationary measure

% Value function accelerator
    s.ac_iter = 8; %Number of iterations                      

% Optimizer options        
    s.options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'Algorithm',{'levenberg-marquardt',1});                                                                                                        
 
%% Simulation by generating sequence of shocks (if s.flag_simulate==1)    
    
    if s.flag_simulate == 1

        %Simulations    
        s.T = 40; % Number of periods
        s.N = 1000000; % Number of firms         
        s.burn = 1; % Number of time periods burnt = s.burn*s.T
        s.seed = 88;
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',s.seed)); % Sets seed for randomizations       

    end

%% Setup asset and productivity grids
    KLS2_setup_grids;                                              

	if s.high_z_ForDebt == 1 % Needs to be AFTER KLS2_setup_grids
	
		% Foreign debt holding
		z_bar_debt = 0.74; % threshold for productivity

		% matrix of debt holding for low export intensity firms
		m.lambda_mat = m.lambda.*ones(s.a_grid_size,s.z_grid_size);    
		ind_debt = find(cumsum(r.z_pi)>z_bar_debt,1,'first');
		m.lambda_mat(:,ind_debt+1:end)= 0;    
		m.lambda_period_2_mat = m.lambda_mat;
		clear ind_debt 

		% matrix of debt holding for high export intensity firms
		m.lambda_X_mat = m.lambda_X.*ones(s.a_grid_size,s.z_grid_size);    
		ind_debt = find(cumsum(r.z_pi)>z_bar_debt,1,'first');
		m.lambda_X_mat(:,ind_debt+1:end)= 0;  
		m.lambda_X_period_2_mat =  m.lambda_X_mat;
		clear ind_debt 
		
	end
          
%% STEP 1: Solve initial steady state

% Initial prices and aggregate quantities     

    if s.load_prices==0

        m.Y = 1.9139;
        m.P = 0.65814;
        m.w = 0.95052;

        if s.model == 2 % Frictionless model

           m.Y = 2.243605626578187;
           m.P = 0.645988153380563;
           m.w=1.111350641402580;

        end

    elseif s.load_prices==1

        if s.model == 1
            load KLS2_prices.mat;
            if s.irf_pm==1 && s.lambda_flag==0
                load KLS2_prices_pm.mat;
                              
            elseif s.irf_pm==1 && s.lambda_flag==1
                load KLS2_prices_pm_lambda0.mat;
            elseif s.irf_pm==1 && s.lambda_flag==2
                load KLS2_prices_pm_lambda1.mat;
            elseif s.irf_pm==1 && s.lambda_flag==3
                load KLS2_prices_pm_lambda05lambdaX0.mat;
            end
			
			if s.UIPdeviations == 1
				load KLS2_prices_UIPdeviations.mat;
			elseif s.high_z_ForDebt == 1 
				load KLS2_prices_high_z_ForDebt.mat;			
            elseif s.financial_crisis==1
                load KLS2_prices_fincrisis.mat;
            elseif s.irp==1
                load KLS2_prices_IRP.mat;    
   
            end
			
            
        elseif s.model == 2
           load KLS2_prices_noFF.mat;
           
           if s.irf_pm==1 && s.lambda_flag==0
                load KLS2_prices_noFF_pm.mat;
                
           end
           
           
        elseif s.model == 3
           load KLS2_prices_noReallocation.mat;
        elseif s.model == 4
           load KLS2_prices_OneType.mat;
        end
        
       

        
        

        m.Y = Prices(1,1);
        m.P = Prices(3,1);
        m.w = Prices(2,1);
    end


% Solve GE    
        if s.GE == 1       
            solver.x0 = [m.w m.Y m.P];
            [solver.z,solver.fval,solver.exitflag] = fsolve(@(x)KLS2_GE(x,m,s,r),log(solver.x0),s.options);                
            [solver.mcc_0,m_0,r_0,r_X_0,s_0,sim_0] = KLS2_GE(solver.z,m,s,r);      
        elseif s.GE ==0       
            solver.x0 = [m.w m.Y m.P];   
            [solver.mcc_0,m_0,r_0,r_X_0,s_0,sim_0] = KLS2_GE(log(solver.x0),m,s,r);  
        end

%% STEP 2: Solve final steady state       
     
    if s.irp==1
        
        m.lambda = 1;
        m.lambda_X = 1;
        
    end

% Initial prices and aggregate quantities    
    
    if s.load_prices==0

        m.Y = 2.131319415265014;
        m.P = 0.540185587060054;
        m.w = 1.065995774679367;

        if s.model == 2 % Frictionless model

            m.Y = 2.710544000000001;
            m.P = 0.579754040000000;
            m.w = 1.179494200000000;

        end

    elseif s.load_prices==1

        m.Y = Prices(1,end);
        m.P = Prices(3,end);
        m.w = Prices(2,end);

    end

% Shocks
    m.r = m.rv(end);
    m.Pf = m.Pfv(end);
    m.zeta_w = m.zeta_w_v(end);
    m.zeta_pi = m.zeta_pi_v(end);
    m.A = m.A_v(end);
    m.Yf = m.Yfv(end);
    m.beta=m.beta_v(end);
    m.Pm = m.pm_v(end);
    m.theta = m.theta_v(end);
    r.z_grid = m.zagg_v(end)*r.z_grid_original;
	
	if s.UIPdeviations == 1
	
		m.lambda = m.lambda_v(end);
		m.lambda_X = m.lambda_X_v(end);
	
	end
        
% Solve GE    
    if s.GE == 1     
        solver.x0 = [m.w m.Y m.P];    
        [solver.z,solver.fval,solver.exitflag] = fsolve(@(x)KLS2_GE(x,m,s,r),log(solver.x0),s.options);                
        [solver.mcc_end,m_end,r_end,r_X_end,s_end,sim_end] = KLS2_GE(solver.z,m,s,r);      
    elseif s.GE ==0             
        solver.x0 = [m.w m.Y m.P];    
        [solver.mcc_end,m_end,r_end,r_X_end,s_end,sim_end] = KLS2_GE(log(solver.x0),m,s,r); 
    end
    
%% Initialize transition objects

    trans.vt = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.apt = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.apt_ind = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.kt = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.et = zeros(s.a_grid_size,s.z_grid_size,trans.N);

    trans.vt_X = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.apt_X = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.apt_ind_X = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.kt_X = zeros(s.a_grid_size,s.z_grid_size,trans.N);
    trans.et_X = zeros(s.a_grid_size,s.z_grid_size,trans.N);

    m.Yt = zeros(trans.N,1); 
    m.wt = zeros(trans.N,1); 
    m.Pt = zeros(trans.N,1);     
        
% Store results into transition objects    
    trans.vt(:,:,1)=r_0.v;
    trans.apt(:,:,1)=r_0.ap;
    trans.apt_ind(:,:,1)=r_0.ap_ind;
    trans.et(:,:,1)=r_0.e;
    trans.kt(:,:,1)=r_0.k;
    trans.vt_X(:,:,1)=r_X_0.v;
    trans.apt_X(:,:,1)=r_X_0.ap;
    trans.apt_ind_X(:,:,1)=r_X_0.ap_ind;
    trans.et_X(:,:,1)=r_X_0.e;
    trans.kt_X(:,:,1)=r_X_0.k;    
    
    trans.vt(:,:,trans.N)=r_end.v;
    trans.apt(:,:,trans.N)=r_end.ap;
    trans.apt_ind(:,:,trans.N)=r_end.ap_ind;
    trans.et(:,:,trans.N)=r_end.e;
    trans.kt(:,:,trans.N)=r_end.k;
    trans.vt_X(:,:,trans.N)=r_X_end.v;
    trans.apt_X(:,:,trans.N)=r_X_end.ap;
    trans.apt_ind_X(:,:,trans.N)=r_X_end.ap_ind;
    trans.et_X(:,:,trans.N)=r_X_end.e;
    trans.kt_X(:,:,trans.N)=r_X_end.k;
     
% Store results from initial and final steady states
    rt{1,1} = r_0;
    rt{trans.N,1} = r_end;

    rt_X{1,1} = r_X_0;
    rt_X{trans.N,1} = r_X_end;

    trans.w0 = m_0.w;
    
% Initialize measure
    if s.flag_simulate==0
        
        trans.measure = sim_0.measure; 
        trans.measure_X = sim_0.measure_X;
    end

    
%% STEP 3: guess sequence of aggregate prices for N period, with N sufficiently large

    m.Yt(1)=m_0.Y;
    m.wt(1)=m_0.w;
    m.Pt(1)=m_0.P;

    m.Yt(trans.N)=m_end.Y;
    m.wt(trans.N)=m_end.w;
    m.Pt(trans.N)=m_end.P;     

    if s.transition==1 
        
        
        
        if s.load_prices == 0

            % Initial guess using linear interpolation
            m.Yt(2)=m_0.Y;
            m.wt(2)=m_0.w;
            m.Pt(2)=m_0.P;    

            stepY=(m.Yt(trans.N)-m.Yt(1))/(trans.N-2);
            stepw=(m.wt(trans.N)-m.wt(1))/(trans.N-2);
            stepP=(m.Pt(trans.N)-m.Pt(1))/(trans.N-2);

            for i=3:trans.N-1        
                m.Yt(i)=m.Yt(i-1)+stepY;
                m.wt(i)=m.wt(i-1)+stepw;
                m.Pt(i)=m.Pt(i-1)+stepP;       
            end     

        elseif s.load_prices == 1

            m.Yt(2:end-1) = Prices(1,2:end-1);
            m.wt(2:end-1) = Prices(2,2:end-1);
            m.Pt(2:end-1) = Prices(3,2:end-1);            

        end
    
        % Guessed prices and quantities, input for solver
        Guess = [m.Yt(2:trans.N-1)' m.wt(2:trans.N-1)' m.Pt(2:trans.N-1)'];
       
%% STEP 4: Solving for the GE prices and quantities   

        if s.transition_GE==0
            Prices_sol = Guess;    
            [mc, m, r, sim_fun, rt, rt_X, trans] = KLS2_transition(Prices_sol,m,r,s,rt,rt_X,trans,sim_0);    
        else
            solve_prices = @(Prices)KLS2_transition(Prices,m,r,s,rt,rt_X,trans,sim_0); %Function that takes prices as inputs and market clearing values as output
            [Prices_sol, mcc_sol, exit_sol] = fsolve(solve_prices,Guess,s.options);
            [mc, m, r, sim_fun, rt, rt_X, trans] = KLS2_transition(Prices_sol,m,r,s,rt,rt_X,trans,sim_0);
        end


        Prices(1,:) = [m.Yt(1) Prices_sol(1:trans.N-2) m.Yt(trans.N)];
        Prices(2,:) = [m.wt(1) Prices_sol(trans.N-2+1:trans.N-2+trans.N-2) m.wt(trans.N)];
        Prices(3,:) = [m.Pt(1) Prices_sol(trans.N-2+trans.N-2+1:end) m.Pt(trans.N)];
            
    % Display convergence statistics      
         disp(' ');
         disp(' ');
    
        if s.flag_simulate==0
           diff_measures = sum((sim_fun.measure_end(:)-sim_end.measure(:)).^2);
           disp(['Convergence of measures: ' num2str(diff_measures)]);

        end

        diff_output_initial = sim_0.Y(1)/m_end.Y;
        diff_output = sim_fun.Y(trans.N-1)/m_end.Y;

        disp(['Convergence in output: ' num2str(diff_output)]);
        disp(['Initial diff in output: ' num2str(diff_output_initial)]);       
        
    %% STEP 5: Output and Figures   

        % Output:
        if s.flag_simulate==0
            KLS2_output     
        end
        

        %End timer    
        elapsed_time = toc;
    
    end  % if s.transition==1     
     
    if s.save_workspace==1
        rt = [];
        rt_X = [];
        sim_fun.M = [];
        sim_fun.M_X = [];
        r_0 = [];
        r_X_0 = [];
        r_end = [];
        r_X_end = [];
        trans = [];

        time_tag = clock; %year, month, day, hour, minute
        filename = ['workspace_' num2str(time_tag(1)) '_' num2str(time_tag(2)) '_' num2str(time_tag(3)) '_' num2str(time_tag(4))...
           '_' num2str(time_tag(5)) '.mat'];
        save(filename);        
    end

    if s.save_prices==1
        if s.model==1
            
            if s.irf_pm == 0 && s.UIPdeviations == 0 && s.high_z_ForDebt == 0 && s.financial_crisis==0
                save('KLS2_prices.mat','Prices');
            elseif s.irf_pm == 1 %pm_shock only (impulse response)
                 if s.lambda_flag == 0  %=0 if baseline calibration
                    save('KLS2_prices_pm.mat','Prices');
                elseif s.lambda_flag == 1  %=1 if lambda=0,lambdaX=0
                    save('KLS2_prices_pm_lambda0.mat','Prices');
                elseif s.lambda_flag == 2 %=2 if lambda=1,lambdaX=1
                    save('KLS2_prices_pm_lambda1.mat','Prices');
                elseif s.lambda_flag == 3   %=3 if lambda=0.5,lambdaX=0 
                    save('KLS2_prices_pm_lambda05lambdaX0.mat','Prices');
                end
            elseif s.UIPdeviations == 1
				save('KLS2_prices_UIPdeviations.mat','Prices');
			elseif s.high_z_ForDebt == 1 
                save('KLS2_prices_high_z_ForDebt.mat','Prices');
            elseif s.financial_crisis==1
                save('KLS2_prices_fincrisis.mat','Prices'); 
			elseif s.irp==1
                save('KLS2_prices_IRP.mat','Prices'); 
            
            end
			
        elseif s.model==2
            
            if s.irf_pm==0
                save('KLS2_prices_noFF.mat','Prices');
            elseif s.irf_pm==1 && s.lambda_flag==0
                save('KLS2_prices_noFF_pm.mat','Prices');
            end
            
        elseif s.model==3
            save('KLS2_prices_noReallocation.mat','Prices');
        elseif s.model==4
            save('KLS2_prices_OneType.mat','Prices');
        
        end
    end
    