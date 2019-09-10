%Productivity
    %Rouwenhorst discretization method
    if s.z_grid_size~=1
        [r.log_z_grid r.z_P r.z_pi] = KLS2_rouwenhorst(m.log_z_rho,m.log_z_sigma,m.log_z_mu,s.z_grid_size);
        r.z_grid = exp(r.log_z_grid');  
        r.z_grid_original = exp(r.log_z_grid');  
    else
        r.z_P = 1;
        r.z_pi = 1;
        r.z_grid = 1;
    end

    %Productivity process statistics
    r.z_min = min(r.z_grid);
    r.z_max = max(r.z_grid);
    r.z_mean = sum(r.z_grid.*r.z_pi);                    
                
%Asset grid                    
    %Check whether assets can be negative (ie. if theta>1+r)
    %If they can, construct grid with negative asset values -- otherwise, construct grid with positive asset values only                       
    FXdebt = 1; %m.lambda + (1-m.lambda)*(m.xi/m.xi_lag), FXdebt is =1 in steady state since xi=xi_lag     
    if m.theta<=(1+r_old)*FXdebt || s.a_grid_size_negative==0 %If assets cannot be negative           
        s.a_grid_size_neg = 0; %Number of negative gridpoints
        s.a_grid_lb = s.a_grid_lb_pos; %Asset grid's lower bound
    else %If assets can be negative
        s.a_grid_size_neg = s.a_grid_size_negative; %Number of negative gridpoints
        s.a_grid_lb = s.a_grid_lb_neg; %Asset grid's lower bound
    end
    s.a_grid_size_pos = s.a_grid_size - s.a_grid_size_neg; %Number of positive gridpoints       

    %Compute it here since asset grid depends on NBL, which is a function of prices
    %If s.a_grid_power>1, allocate more points towards lower values 
    %If s.a_grid_power=1, equally spaced points
    %If s.a_grid_power<1, allocate more points towards higher values
    %Typically, will want to have more points toward lower values of the asset state since policy functions have more curvature at
    %those points, and become more linear for higher values    
    r.a_grid_pos = linspace(0,1,s.a_grid_size_pos);
    r.a_grid_pos = r.a_grid_pos.^s.a_grid_power_pos;
    r.a_grid_pos = r.a_grid_pos.*(s.a_grid_ub-s.a_grid_lb_pos) + s.a_grid_lb_pos; 
    if s.a_grid_size_neg>0
        r.a_grid_neg = linspace(0,1,s.a_grid_size_neg);
        r.a_grid_neg = r.a_grid_neg.^s.a_grid_power_neg;
        r.a_grid_neg = r.a_grid_neg.*(s.a_grid_ub_neg-s.a_grid_lb_neg) + s.a_grid_lb_neg; 
        r.a_grid = [r.a_grid_neg r.a_grid_pos];            
    elseif s.a_grid_size_neg==0
        r.a_grid = r.a_grid_pos;
    end  
    
%Grids in matrix form
    %Assets
        r.a_grid_mat = r.a_grid'*ones(1,s.z_grid_size);         
        r.a_grid_mat_np = r.a_grid_mat<=0; %Non-positive assets
        r.a_grid_mat_p = r.a_grid_mat>0; %Positive assets 
        
    %Productivity
        r.z_grid_mat = ones(s.a_grid_size,1)*r.z_grid'; 
        r.z_grid = m.zagg*r.z_grid_original;
    