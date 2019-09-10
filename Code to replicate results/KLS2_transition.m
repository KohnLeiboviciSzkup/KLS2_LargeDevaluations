function [mc m r sim_fun rt rt_X trans] = KLS2_transition(Guess,m,r,s,rt,rt_X,trans,sim_0)

  
    Yguess = [m.Yt(1) Guess(1:trans.N-2) m.Yt(trans.N)];
    wguess = [m.wt(1) Guess(trans.N-2+1:trans.N-2+trans.N-2) m.wt(trans.N)];
    Pguess = [m.Pt(1) Guess(trans.N-2+trans.N-2+1:end) m.Pt(trans.N)];    

    r_initial = r; 

    for t=trans.N-1:-1:2    

        %% Shocks

        % shock to interest rate
        m.r = m.rv(t);
        m.r_lag = m.rv(t-1);

        % shock to foreign price
        m.Pf_lag = m.Pfv(t-1);
        m.Pf = m.Pfv(t);

        % shock to foreign demand
        m.Yf = m.Yfv(t);

        % Shock to pm
        m.Pm = m.pm_v(t);
        
        % Shock to zeta_w and zeta_pi
        m.zeta_w = m.zeta_w_v(t);
        m.zeta_pi = m.zeta_pi_v(t);

        % shock to the production of domestic final goods
        m.A = m.A_v(t);

        % shock to collateral constraint
        m.theta = m.theta_v(t);

        % shock to discount factor
        m.beta = m.beta_v(t);        
        
        % shock to depreciation rate
        m.delta = m.delta_v(t);
		
		% fraction of foreign debt
		if s.UIPdeviations == 1;
			m.lambda = m.lambda_v(t);
			m.lambda_X = m.lambda_X_v(t);
		end
        
        %aggregate productivity shock
        r.z_grid = m.zagg_v(t)*r.z_grid_original;
        r.z_grid_mat = ones(s.a_grid_size,1)*r.z_grid'; 

        r_initial.z_grid = m.zagg_v(t)*r_initial.z_grid_original;
        r_initial.z_grid_mat = ones(s.a_grid_size,1)*r_initial.z_grid';         
        
    %% GE prices
    
        m.Y = Yguess(t);
        m.Y_lag = Yguess(t-1);             

        m.w = wguess(t);
        m.w_lag = wguess(t-1);             

        m.P = Pguess(t);
        m.P_lag = Pguess(t-1);               

        m.xi_real = m.Pf/m.P;
        m.xi_real_lag=m.Pf_lag/m.P_lag;            
        
    %Value function
        vp=squeeze(trans.vt(:,:,t+1));  
        vp_X=squeeze(trans.vt_X(:,:,t+1));

    %% Solve static and dynamic problems   
    
        % Period 2
        if t==2 
                                    
                           
            m.Y = m.Y*m.A^(m.sigma-1);

            %Static problem                
            r_temp = KLS2_staticproblem_period2(m,s,r_initial,rt);

            if s.model==1  || s.model == 2  || s.model == 4

               r_X_temp = KLS2_staticproblem_period2(m,s,r_initial,rt_X,2); 

            elseif s.model==3

               r_X_temp = KLS2_staticproblem_typeX_period2(m,s,r_initial,rt_X); 

            end

            m.Y = m.Y/(m.A^(m.sigma-1));             

            %Dynamic problem
            r_temp2 = KLS2_dynamicproblem_trans(m,s,r_temp,vp,2);

            if s.model==1  || s.model == 2  || s.model == 4

                r_X_temp2 = KLS2_dynamicproblem_trans(m,s,r_X_temp,vp_X,2,2);

            elseif s.model==3

                r_X_temp2 =  KLS2_dynamicproblem_typeX_trans(m,s,r_X_temp,vp_X,2);

            end                  
                                                
                
              

        else

        %Static problem
            m.Y = m.Y*m.A^(m.sigma-1);
            
            r_temp = KLS2_staticproblem(m,s,r_initial);

            if s.model==1  || s.model == 2  || s.model == 4
                
               r_X_temp = KLS2_staticproblem(m,s,r_initial,2);
                
            elseif s.model==3
                
               r_X_temp = KLS2_staticproblem_typeX(m,s,r_initial);
                
            end            
           
            m.Y = m.Y/(m.A^(m.sigma-1));

        %Dynamic problem
            r_temp2 = KLS2_dynamicproblem_trans(m,s,r_temp,vp,1);

            if s.model==1  || s.model == 2  || s.model == 4
                
               r_X_temp2 = KLS2_dynamicproblem_trans(m,s,r_X_temp,vp_X,1,2);
                               
            elseif s.model==3
                
               r_X_temp2 = KLS2_dynamicproblem_typeX_trans(m,s,r_X_temp,vp_X,1);
                
            end            
            
        end        

        
        % Storing output
        trans.vt(:,:,t) = r_temp2.v;
        trans.apt(:,:,t) = r_temp2.ap;
        trans.et(:,:,t) = r_temp2.e; 
        trans.apt_ind(:,:,t) = r_temp2.ap_ind;
        trans.kt(:,:,t) = r_temp2.k; 
        rt{t,1} = r_temp2;

        trans.vt_X(:,:,t) = r_X_temp2.v;
        trans.apt_X(:,:,t) = r_X_temp2.ap;
        trans.et_X(:,:,t) = r_X_temp2.e;   
        trans.apt_ind_X(:,:,t) = r_X_temp2.ap_ind;
        trans.kt_X(:,:,t) = r_X_temp2.k; 
        rt_X{t,1} = r_X_temp2;

    end


    if s.flag_simulate == 0

        [sim_fun, rt, rt_X] = KLS2_simulate_trans(m,s,r,rt,rt_X,trans,Yguess,Pguess,sim_0);

    elseif s.flag_simulate == 1

       [sim_fun, rt, rt_X] = KLS2_simulate_shocks_trans(m,s,r,rt,rt_X,trans,Yguess,Pguess,sim_0);        

    end

    % Market clearing conditions
     
        mc = [sim_fun.mc_n(2:trans.N-1)' sim_fun.mc_y(2:trans.N-1)' sim_fun.mc_y_belief(2:trans.N-1)']';
     

    end


