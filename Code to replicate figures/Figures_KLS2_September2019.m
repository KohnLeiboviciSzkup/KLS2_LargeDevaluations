%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to replicate the figures in 
% "Financial Frictions and Export Dynamics in Large Devaluations"
% Journal of International Economics, 2019
% David Kohn, Fernando Leibovici and Michal Szkup 
% August 2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Setup

    close all;
    clc;
    clear;

    sav = 0; % if sav = 1 then figures are saved 
             % if sav = 0 then figures are not saved     

    f=14;   % Fontsize for every plot
    lw=3;   % Linewidth for every plot
         

%% Figure 1 in Paper and Online Appendix
      
    % Load data of rer and real exports and use it to date devaluations,
    % and to produce figures of export dynamics around devaluation episodes.
    % A devaluation is an INCREASE in rer
    
    load data_yearly         
    
    detrending = 3;   % if detrending = 0 then no detrending
                      % if detrending = 1 then log-linear detrending
                      % if detrending = 2 then subtract pre-dev log-growth  
                      % if detrending = 3 then subtract average log-growth

    w = 3;   % look at exports w_before periods before to w periods after devaluation date
    w_before = 1;
    k = 1;   % initial year of the data (row number), k=1 corresponds to 1980
    n = w+w_before+1;

    rexpo  = log(r_expo(k:end,:));
    rer    = log(rer(k:end,:));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find large devaluation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define threshold increase in rer
    d_rer_x = 20;  % A large devaluation is an increase in rer of more than d_rer_x % 

    % Calculate log change in rer (rer in logs already!)
    d_rer=rer(2:end,:)-rer(1:end-1,:);

    % Find large devaluations -only increases in rer- (put a 1 when there is a large devaluation)
    dev=(100*d_rer>=d_rer_x);

    %Restrict to devaluation with w_before periods before
    if w_before>=1
        dev(1:w_before,:)=0;
    end

    % Control ONLY for consecutive periods (date devaluation at the first of two consecutive drops)
    for j=1:size(dev,2)
        for i=1:size(dev,1)-w
             if dev(i,j)==1 
                 dev(i+1,j)=0;
             end
        end
    end

    % Number of large devaluations in the sample
    ndev = sum(sum(dev));

    % Time
    t=1980:2013;
    t=t(k:end);
    t=t';
    dev=[zeros(1,size(dev,2)); dev]; % For 1980 observation that we lost when we calculated changes
    num_t = numel(t);
    
     
    % Countries
    num_countries = numel(countries);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct exports series around a window of devaluation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Indicator for devaluation period
    dev_w=dev;
    for j=1:num_countries
        for i=w_before+1:num_t-w   
             if dev(i,j)==1
                 dev_w(i-w_before:i+w,j) = 1;             
             end    
        end          
    end
    dev_period = dev;
    dev = dev_w;
    clear dev_w

       
    % Applying the indicator for devaluation dev:
    rexpo_dev     = rexpo.*dev;
    rer_dev       = rer.*dev;
    
    % Differencing data (data already in logs!):
    d_rer     = [ zeros(1,size(d_rer,2)) ; d_rer].*dev;
    d_rexpo   = rexpo(2:end,:)-rexpo(1:end-1,:);
    d_rexpo   = [zeros(1,size(d_rer,2)) ; d_rexpo].*dev;
      
    % (Log) Linear detrending    (variables in logs)
    rexpo_cl = zeros(num_t,num_countries);
    for i=1:num_countries
        rexpo_cl(isnan(rexpo(:,i))==0,i) = detrend(rexpo(isnan(rexpo(:,i))==0,i));
    end
 
    % Applying the indicator for devaluation dev:
    rexpo_cl_dev   = rexpo_cl.*dev;

    % Get all devaluation data together into simpler matrices    
    w_predev = 4;

    % Initialization:
    rexpo_dev_w     =  zeros(w+w_before+1,1);
    rer_dev_w       =  zeros(w+w_before+1,1);
    d_rer_dev     =  zeros(w+w_before+1,1);
    d_rexpo_dev    =  zeros(w+w_before+1,1);
    rexpo_cl_dev_w    =  zeros(w+w_before+1,1);
    d_rexpo_predev = NaN(w_predev,1);
    countries_dev = {};
    m=1; %initialization
    for j=1:num_countries
        ii=1;
        for i=1:num_t
            if dev_period(ii,j)==1 %dev(ii,j)>0 %what about those in first few periods?, %NAN does not count
               countries_dev{1,m} =  countries{1,j}; 
               rexpo_dev_w(:,m)     =  rexpo_dev(ii-w_before:ii+w,j);
               rer_dev_w(:,m)       =  rer_dev(ii-w_before:ii+w,j);
               d_rer_dev(:,m)     =  d_rer(ii-w_before:ii+w,j);
               d_rexpo_dev(:,m)   =  d_rexpo(ii-w_before:ii+w,j);
               rexpo_cl_dev_w(:,m)    =  rexpo_cl_dev(ii-w_before:ii+w,j);
               predev_ind = min(max((ii-w_before)-w_predev,1),w_predev);
               if ii>w_before+1
                   temp_ind = 1;
                   for h=1:predev_ind 
                        temp_ind2 = ii-w_before-temp_ind+1;
                        d_rexpo_predev(h,m) = rexpo(temp_ind2,j) - rexpo(temp_ind2-1,j);
                        temp_ind = temp_ind + 1;
                   end
               end
               m = m+1; % update column index
            end
            d_rexpo_predev(d_rexpo_predev==0) = NaN;
            d_rexpo_predev_mean = nanmean(d_rexpo_predev);
            %ii = min(ii+w+w_before+w_before+1,size(dev,1)); %at least w_before apart from last episode...
            ii=ii+1; %update row index

            if ii>=num_t-w+1 
                   break
            end

        end
        d_rexpo_dev2(:,j) = rexpo(2:end,j) - rexpo(1:end-1,j);
    end
    d_rexpo_dev3=nanmean(d_rexpo_dev2);
    dev_country=max(dev);
    temp=(d_rexpo_dev3.*dev_country);
    temp=temp(temp>0);
    d_rexpo_mean_allperiods=[temp(1:4) temp(6:7) temp(7) temp(7) temp(9:11) temp(11)]; 
    %without temp(5) Japan, without temp(8) Russia
    
    rexpo_dev = rexpo_dev_w;
    rexpo_cl_dev=rexpo_cl_dev_w;
    rer_dev = rer_dev_w;

    clear rexpo_dev_w rer_dev_w rexpo_cl_dev_w 

    % Current devaluations in order from 1:m:
    % Argentina	2002
    % Brazil	1999
    % Iceland 	2008/2009(solo 2008)
    % Indonesia	1998
    % Korea	1998
    % Malaysia	1998
    % Mexico	1982
    % Mexico	1986
    % Mexico	1995 (9th devaluation)
    % Russia	1999
    % Turkey	2001
    % Venezuela	2002
    % Venezuela	2010

    % Note: we drop Russia (some data is missing) and Japan (year of
    % devaluation close to end of sample, so not enough years afterwards
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute variables to plot in figures %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find countries with missing exports observations
    temp = sum(isnan(rexpo_dev))==0; % gets rid of Russian devaluation (before we got rid of Japan with windows)
    rexpo_dev     =  rexpo_dev(:,temp);
    rer_dev       =  rer_dev(:,temp);
    d_rer_dev       =  d_rer_dev(:,temp);
    d_rexpo_dev     =  d_rexpo_dev(:,temp);
    d_rexpo_predev_mean     = d_rexpo_predev_mean(:,temp);
    d_rexpo_predev_mean(isnan(d_rexpo_predev_mean)==1) = 0;
    rexpo_cl_dev   =  rexpo_cl_dev(:,temp);
 
    % Real exchange rate 
    d_rer_mean = nanmean(d_rer_dev'); 
    d_rer_med = nanmedian(d_rer_dev');
    rel_rer = rer_dev - ones(n,1)*rer_dev(1,:);
    rel_rer_mean = nanmean(rel_rer');
    rel_rer_med = nanmedian(rel_rer');   
        
    % Real exports   
    
    % Average growth by country
    d_rexpo_raw = rexpo(:,2:end)-rexpo(:,1:end-1);
    d_rexpo_mean_country = nanmean(d_rexpo_raw');

    % Annual growth
    d_rexpo_mean = nanmean(d_rexpo_dev');
    d_rexpo_med = nanmedian(d_rexpo_dev');

    % Growth relative to pre-devaluation period
    rel_rexpo = rexpo_dev - ones(n,1)*rexpo_dev(1,:);
    rel_rexpo_mean = nanmean(rel_rexpo');
    rel_rexpo_med = nanmedian(rel_rexpo');    

    % Growth relative to pre-devaluation period subtracting pre-dev avg growth
    rel_rexpo_detrend = rexpo_dev - ones(n,1)*rexpo_dev(1,:) - (0:n-1)'*d_rexpo_predev_mean;
    rel_rexpo_detrend_mean = nanmean(rel_rexpo_detrend');
    rel_rexpo_detrend_med = nanmedian(rel_rexpo_detrend');   
    rexpo_elasticity_detrend = rel_rexpo_detrend./rel_rer;
    rexpo_elasticity_detrend_med = nanmedian(rexpo_elasticity_detrend,2) ;
    rexpo_elasticity_detrend_med(1)=0;     

    % Growth relative to pre-devaluation period subtracting avg growth
    rel_rexpo_detrend_avg = rexpo_dev - ones(n,1)*rexpo_dev(1,:) - (0:n-1)'*d_rexpo_mean_allperiods;
    rel_rexpo_detrend_avg_mean = nanmean(rel_rexpo_detrend_avg');
    rel_rexpo_detrend_avg_med = nanmedian(rel_rexpo_detrend_avg');   
    rexpo_elasticity_detrend_avg = rel_rexpo_detrend_avg./rel_rer;
    rexpo_elasticity_detrend_avg_med = nanmedian(rexpo_elasticity_detrend_avg,2) ;
    rexpo_elasticity_detrend_avg_med(1)=0;    
    
    % Elasticity to RER relative to pre-devaluation period
    rexpo_elasticity = rel_rexpo./rel_rer;
    rexpo_elasticity(1,:) = 0;
    rexpo_elasticity_mean = nanmean(rexpo_elasticity');
    rexpo_elasticity_med = nanmedian(rexpo_elasticity');
    rexpo_elasticity_l = (rexpo_cl_dev - ones(n,1)*rexpo_cl_dev(1,:))./rel_rer;
    rexpo_elasticity_l(1,:) = 0;
    rexpo_elasticity_l_mean = nanmean(rexpo_elasticity_l');
    rexpo_elasticity_l_med = nanmedian(rexpo_elasticity_l');          
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figures Data for Paper %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    periods = [-w_before:1:w];

    % RER
    figure(101)
    plot(periods,rel_rer_med,'LineWidth',lw); grid on; 
    axis([periods(1) periods(end) (min(rel_rer_med)-0.1) (max(rel_rer_med)+0.1) ])
    line([0,0],[-25,45],'linestyle','--','Color','r','LineWidth',lw)
    set(gca,'XTick',-1:1:3,'YTick',0:0.1:0.5);
    box off
    xlabel('Years after devaluation','Fontsize',f)
    title({'Real Exchange Rate';'(log change from pre-devaluation year, median)'},'Fontsize',f)
    set(gca,'Fontsize',f)

    if sav==1

        % saving RER figure:
        Filename = 'KLS2_data_RER.eps'; 
        saveas(gcf,Filename,'epsc')  

       % fixing the dotted line:
        fix_dottedline(Filename)

    end
    
    if detrending == 0
    
        % Exports:
        figure(102)
        plot(periods,rexpo_elasticity_med,'LineWidth',lw); grid on; 
        axis([periods(1) periods(end) (min(rexpo_elasticity_med)-0.1) (max(rexpo_elasticity_med)+0.1) ])
        line([0,0],[-25,45],'linestyle','--','Color','r','LineWidth',lw)
        set(gca,'XTick',-1:1:3,'YTick',0:0.4:1.6);
        box off
        xlabel('Years after Devaluation','Fontsize',f)
        title('Elasticity of Exports to RER (median)','Fontsize',f)
        set(gca,'Fontsize',f)

        if sav==1

            % saving real Exports figure:
            Filename = 'KLS2_data_XElasticity.eps'; 
            saveas(gcf,Filename,'epsc')  
           % fixing the dotted line:
            fix_dottedline(Filename)
      
        end

       
    elseif detrending == 1    
       
        % Exports: 
        figure(102)
        plot(periods,rexpo_elasticity_l_med,'LineWidth',lw); grid on; 
        axis([periods(1) periods(end) (min(rexpo_elasticity_l_med)-0.1) (max(rexpo_elasticity_l_med)+0.1) ])
        line([0,0],[-25,45],'linestyle','--','Color','r','LineWidth',lw)
        set(gca,'XTick',-1:1:3,'YTick',-0.2:0.2:0.8);
        box off
        xlabel('Years after Devaluation','Fontsize',f)
        title('Elasticity of Exports to RER (median)','Fontsize',f)
        set(gca,'Fontsize',f)

        if sav==1

            % saving real Exports figure:
            Filename = 'KLS2_data_XElasticity.eps'; 
            saveas(gcf,Filename,'epsc')  
           % fixing the dotted line:
            fix_dottedline(Filename)
       
        end

      
    
    elseif detrending == 2        % Subtracting Pre-Dev log-growth
       
        % Exports: 
        figure(102)
        
        plot(periods,rexpo_elasticity_detrend_med','LineWidth',lw); grid on; 
        line([0,0],[-1,1],'linestyle','--','Color','r','LineWidth',lw)
        axis([periods(1) periods(end) -0.1 0.2 ])
        set(gca,'XTick',-1:1:3,'YTick',-0.2:0.1:0.2);
        box off
        xlabel('Years after Devaluation','Fontsize',f)
        title('Elasticity of Exports to RER (median)','Fontsize',f)
        set(gca,'Fontsize',f)

        if sav==1

            Filename = 'KLS2_data_XElasticity.eps'; 
            saveas(gcf,Filename,'epsc')  

           % fixing the dotted line:
            fix_dottedline(Filename)       
        end


    elseif detrending == 3        % Subtracting Average log-growth
       
        % Exports: 
        figure(102)
        plot(periods,rexpo_elasticity_detrend_avg_med,'LineWidth',lw); grid on; 
        line([0,0],[-25,45],'linestyle','--','Color','r','LineWidth',lw)
        axis([periods(1) periods(end) -0.1 0.8 ])
        set(gca,'XTick',-1:1:3,'YTick',-0.2:0.2:0.6);
        box off
        xlabel('Years after Devaluation','Fontsize',f)
        title('Elasticity of Exports to RER (median)','Fontsize',f)
        set(gca,'Fontsize',f)

        if sav==1


            Filename = 'KLS2_data_XElasticity.eps'; 
            saveas(gcf,Filename,'epsc')  

           % fixing the dotted line:
            fix_dottedline(Filename) 
        end

    end    
    

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figures in the Paper
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clearvars -except sav f lw

         
    N=20; % Number of transition periods

    nT=5;
    t_talk = 1:nT;
    t_talk2 = 1:12;
    t_talk_label = t_talk - 2;
    t_talk2_label = t_talk2 - 2;


    load KLS2_Workspaces_September2019
    
% Workspaces
% Mex_3shocks_FF
% Mex_3shocks_NoFF
% Mex_Pm_FF
% Mex_Pm_NoFF
% Mex_Pm_lambda0
% Mex_Pm_lambda1
% Mex_Pm_lambda05lambdax0
% Mex_high_z_ForDebt
% Mex_UIPdeviations
% Mex_3shocks_FF_FinancialCrisis
% Mex_3shocks_FF_IRP
% Mex_3shocks_FF_PmRShocks
% Mex_3shocks_FF_PmZShocks
% Mex_3shocks_FF_OnlyPmShock
% Mex_3shocks_NoFF_BaselineShocks
% Mex_3shocks_NoReallocation
% Mex_3shocks_OneType
% MicroEvidence_FF


% Series in each workspace:
%	1    sim_fun.share_x
%	2	 sim_fun.share_x_composition
%	3	 sim_fun.share_nx
%	4	 sim_fun.share_nx_composition
%	5	 sim_fun.ext_const
%	6	 sim_fun.ext_const_X
%	7	 sim_fun.ext_const_all
%	8	 sim_fun.int_const_k_ind
%	9	 sim_fun.int_const_k_ind_X
%	10	 sim_fun.int_const_k_ind_all
%	11	 sim_fun.int_const_k_prc
%	12	 sim_fun.int_const_k_prc_X
%	13	 sim_fun.int_const_k_prc_all
%	14	 sim_fun.int_const_pi_ind
%	15	 sim_fun.int_const_pi_ind_X
%	16	 sim_fun.int_const_pi_ind_all
%	17	 sim_fun.int_const_pi_prc
%	18	 sim_fun.int_const_pi_prc_X
%	19	 sim_fun.int_const_pi_prc_all
%	20	 sim_fun.Pd
%	21	 sim_fun.Yd
%	22	 sim_fun.PdYd
%	23	 sim_fun.Yx_sigma
%	24	 sim_fun.Px
%	25	 sim_fun.Yx
%	26	 sim_fun.Yx_X
%	27	 sim_fun.PxYx
%	28	 sim_fun.PxYx_X
%	29	 sim_fun.PxYx_all
%	30	 sim_fun.PxYx_all_USD
%	31	 sim_fun.PxYx_USD
%	32	 sim_fun.PxYx_X_USD
%	33	 sim_fun.Pm
%	34	 sim_fun.Ym
%	35	 sim_fun.PmYm
%	36	 sim_fun.K_all
%	37	 sim_fun.K
%	38	 sim_fun.K_X
%	39	 sim_fun.inv_agg
%	40	 sim_fun.inv_agg_X
%	41	 sim_fun.inv_agg_all
%	42	 sim_fun.N_X
%	43	 sim_fun.N
%	44	 sim_fun.FC
%	45	 sim_fun.N_all
%	46	 sim_fun.Pcpi
%	47	 sim_fun.Ycpi
%	48	 sim_fun.Y
%	49	 sim_fun.P
%	50	 sim_fun.xi_real
%	51	 sim_fun.C_all
%	52	 sim_fun.C
%	53	 sim_fun.C_X
%	54	 sim_fun.n_supply
%	55	 sim_fun.n_demand
%	56	 sim_fun.mc_n
%	57	 sim_fun.a_supply
%	58	 sim_fun.a_demand
%	59	 sim_fun.mc_a
%	60	 sim_fun.y_supply
%	61	 sim_fun.y_demand
%	62	 sim_fun.mc_y
%	63	 sim_fun.mc_p_belief
%	64	 sim_fun.mc_y_belief
%	65	 sim_fun.credit_all
%	66	 sim_fun.credit_repayment_all
%	67	 sim_fun.d_agg
%	68	 sim_fun.d_agg_X
%	69	 sim_fun.d_agg_repayment
%	70	 sim_fun.r_tilde
%	71	 sim_fun.r_tilde_X
%	72	 sim_fun.r_effective_avg
%	73	 sim_fun.GDP
%	74	 sim_fun.TFP
%	75	 sim_fun.NX_GDP
%	76	 sim_fun.NFA_GDP
%	77	 sim_fun.sales_avg_nx
%	78	 sim_fun.sales_avg_x
%	79	 sim_fun.X_GDP
%	80	 sim_fun.credit_gdp
%	81	 sim_fun.xpremium_sales
%	82	 sim_fun.X_D
%	83	 sim_fun.a_min_share
%	84	 sim_fun.a_max_share
%	85	 sim_fun.Yx_sumall
%	86	 sim_fun.credit_X
%	87	 sim_fun.credit_repayment_X
%	88	 sim_fun.d_agg_X_repayment
%	89	 sim_fun.exp_avg_x
%	90	 sim_fun.GDP_Laspeyres
%	91	 sim_fun.X_Laspeyres
%	92	 sim_fun.D_Laspeyres
%	93	 sim_fun.X_Laspeyres_Type1
%	94	 sim_fun.X_Laspeyres_Type2
%	95	 sim_fun.w
%	96	 m.pm_v
%	97	 m.zagg_v
%	98	 m.rv


%  Series in the data (deviations from SS)

    data_RER = [0.0000
                0.4054
                0.299
                0.1539
                0.1376 ]';

    data_GDP = [0.0000
                -0.0835
                -0.0505
                -0.0074
                0.0144  ]';        

    data_Inv = [ 0.0000
                -0.0333
                -0.0188
                -0.0066
                 0.014     ]';                


    % Data from Mexican Devaluation, Real Exports log growth subtracting mean
    % log growth

    data_Xreal     = [0
                      0.130
                      0.211
                      0.253
                      0.270     ]';

        
%% Figure 2 in the Paper
  
    % RER
    figure(201);
    plot(t_talk_label,log(Mex_3shocks_FF(50,t_talk))-log(Mex_3shocks_FF(50,1)),'r-',t_talk_label,log(Mex_3shocks_NoFF(50,t_talk))-log(Mex_3shocks_NoFF(50,1)),'k-',t_talk_label,data_RER,'b','LineWidth',lw);    
    h=title({'Real Exchange Rate';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:0.5);
    axis([t_talk_label(1) t_talk_label(end) 0 0.45]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    
    
    if sav==1

        Filename = 'KLS2_RER_noFFnew'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end

    % Real GDP    
    figure(202);
    plot(t_talk_label,log(Mex_3shocks_FF(90,t_talk))-log(Mex_3shocks_FF(90,1)),'r-',t_talk_label,log(Mex_3shocks_NoFF(90,t_talk))-log(Mex_3shocks_NoFF(90,1)),'k-',t_talk_label,data_GDP,'b','LineWidth',lw);    
    h=title({'Real GDP';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.1 0.1]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_GDP_noFFnew'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
    
    
    % Investment_Output
    
    Inv_GDP_FF = (Mex_3shocks_FF(41,:)./Mex_3shocks_FF(73,:));
    Inv_GDP_NoFF = (Mex_3shocks_NoFF(41,:)./Mex_3shocks_NoFF(73,:));
 
    
    figure(203);
    plot(t_talk_label,Inv_GDP_FF(1,t_talk)-Inv_GDP_FF(1,1),'r-',t_talk_label,Inv_GDP_NoFF(1,t_talk)-Inv_GDP_NoFF(1,1),'k-',t_talk_label,data_Inv,'b','LineWidth',lw);    
    h=title({'Investment over GDP';'(change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.06 0.06]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_Inv_noFFnew'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
    
    % Exports elasticity

    %	91	 sim_fun.X_Laspeyres
    %	50	 sim_fun.xi_real

    XElast_FF = log(Mex_3shocks_FF(91,:)./Mex_3shocks_FF(91,1))./log(Mex_3shocks_FF(50,:)./Mex_3shocks_FF(50,1));
    XElast_NoFF = log(Mex_3shocks_NoFF(91,:)./Mex_3shocks_NoFF(91,1))./log(Mex_3shocks_NoFF(50,:)./Mex_3shocks_NoFF(50,1));
    XElast_Data =data_Xreal./data_RER;
    XElast_Data(1) = 0;
    XElast_FF(1)=0;
    XElast_NoFF(1)=0;


    %Data vs SMM-Trans, 5 periods   
    figure(204);
    plot(t_talk_label,XElast_FF(1,t_talk),'r-',t_talk_label,XElast_NoFF(1,t_talk),'k-',t_talk_label,XElast_Data,'b','LineWidth',lw);
    h=title('Elasticity of Exports to RER');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:1:3.0);
    axis([t_talk_label(1)-0.01 t_talk_label(end) 0 3]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    
    
    if sav==1

        Filename = 'KLS2_XElasticity3_data_noFFnew'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
    
%% Figure 3 in the Paper: Micro Level Evidence, Median productivity exporters

MicroEvidence_FF(:,1:4)=log(MicroEvidence_FF(:,1:4))-log(MicroEvidence_FF(1,1:4));
MicroEvidence_FF(:,5:8)=(MicroEvidence_FF(:,5:8))-(MicroEvidence_FF(1,5:8));
MicroEvidence_FF(:,9:11)=log(MicroEvidence_FF(:,9:11))-log(MicroEvidence_FF(1,9:11));

% Data: Changes relative to pre-devaluation SS
% Columns:
%(1) X (debt>0), log change
%(2) X (debt<0), log change
%(3) Y (debt>0), log change
%(4) Y (debt<0), log change
%(5) Inv (debt>0), change
%(6) Inv (debt<0), change
%(7) Inv/Y (debt>0), change
%(8) Inv/Y (debt<0), change
%(9) # firms (debt>0), log change
%(10) Domestic sales (debt>0), log change
%(11) Domestic (debt<0), log change



    
%     MicroEvidence_FF =[
%     0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
%     0.673248	0.701009	0.129046	0.251123	-0.007140	0.004642	-0.071584	0.092430	-0.135261	-0.184193	-0.127853
%     0.534876	0.740107	-0.037745	0.191759	0.016025	0.011175	0.167607	0.277805	-0.227799	-0.376576	-0.328321
%     0.643018	0.818877	0.009049	0.246753	0.012555	-0.000308	0.122342	-0.030931	-0.241403	-0.390095	-0.312988
%     0.707155	0.798439	0.046970	0.223019	0.011050	-0.000190	0.100666	-0.026098	-0.293481	-0.380318	-0.342426
%     0.753759	0.771184	0.077870	0.193523	0.009844	-0.001134	0.084172	-0.048669	-0.334757	-0.367034	-0.375825
%     0.778883	0.771938	0.087685	0.175522	0.005917	-0.001392	0.047005	-0.054588	-0.370753	-0.374971	-0.427485
%     0.778027	0.741040	0.085988	0.148473	0.005268	-0.002171	0.041332	-0.074751	-0.399996	-0.377660	-0.447482
%     0.774146	0.655916	0.080837	0.079447	0.002632	-0.002641	0.017902	-0.086442	-0.423444	-0.384312	-0.487823
%     0.768775	0.571035	0.071682	0.017125	0.000442	-0.002824	-0.001371	-0.090956	-0.434162	-0.397968	-0.512011
%     ];
% 


    figure(301);
    plot(t_talk2_label(1:5),MicroEvidence_FF(t_talk2(1:5),7),'k-',t_talk2_label(1:5),MicroEvidence_FF(t_talk2(1:5),8),'r--','LineWidth',lw);
    h=ylabel('Change from s.s.');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );        
    h=title('Investment over Output');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );     
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:1:t_talk2_label(5)],'YTick',-0.1:0.1:0.6);
    %axis([t_talk2_label(1)-0.01 t_talk2_label(5) -0.1 0.3]);
    axis([t_talk2_label(1)-0.01 t_talk2_label(5) -0.1 0.6]);
    
    grid on;
    hleg = legend('Firms with Debt','Firms with Savings','Location','NorthEast');            
    thh = findobj(hleg,'type','text');
    set(thh,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_MicroEvidence_FF_1'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end

    figure(302);
    plot(t_talk2_label(1:6),MicroEvidence_FF(t_talk2(1:6),3),'k-',t_talk2_label(1:6),MicroEvidence_FF(t_talk2(1:6),4),'r--','LineWidth',lw);
    h=title('Output');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );   
    h=ylabel('Log change from s.s.');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:1:t_talk2_label(5)],'YTick',0:0.1:0.25);
    %axis([t_talk2_label(1)-0.01 t_talk2_label(5) -0.05 0.30]);
    axis([t_talk2_label(1)-0.01 t_talk2_label(5) -0.06 0.30]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize',f);  
    
    if sav==1

        Filename = 'KLS2_MicroEvidence_FF_2'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end
   
    figure(303);
    plot(t_talk2_label(1:6),MicroEvidence_FF(t_talk2(1:6),1),'k-',t_talk2_label(1:6),MicroEvidence_FF(t_talk2(1:6),2),'r--','LineWidth',lw);
    h=ylabel('Log change from s.s.');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );        
    h=title('Exports');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );     
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:1:t_talk2_label(5)],'YTick',0:0.2:0.9);
    axis([t_talk2_label(1)-0.01 t_talk2_label(5) -0.01 0.9]);
    grid on;
    
    if sav==1

        Filename = 'KLS2_MicroEvidence_FF_3'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end
    
%% Figure 4 in the Paper
          
      % RER, Inv/output, RGDP, XElast
     NoReallocation_3shocks(1,1:10)=Mex_3shocks_NoReallocation(50,1:10);
     NoReallocation_3shocks(2,1:10)=Mex_3shocks_NoReallocation(41,1:10)./Mex_3shocks_NoReallocation(73,1:10);
     NoReallocation_3shocks(3,1:10)=Mex_3shocks_NoReallocation(90,1:10);
     NoReallocation_3shocks(4,1:10)=(log(Mex_3shocks_NoReallocation(91,1:10))-log(Mex_3shocks_NoReallocation(91,1)))./(log(Mex_3shocks_NoReallocation(50,1:10))-log(Mex_3shocks_NoReallocation(50,1)));
     NoReallocation_3shocks(4,1)=0;
     NoReallocation_3shocks=NoReallocation_3shocks';
     
     
     OneType_3shocks(1,1:10)=Mex_3shocks_OneType(50,1:10);
     OneType_3shocks(2,1:10)=Mex_3shocks_OneType(41,1:10)./Mex_3shocks_OneType(73,1:10);
     OneType_3shocks(3,1:10)=Mex_3shocks_OneType(90,1:10);
     OneType_3shocks(4,1:10)=(log(Mex_3shocks_OneType(91,1:10))-log(Mex_3shocks_OneType(91,1)))./(log(Mex_3shocks_OneType(50,1:10))-log(Mex_3shocks_OneType(50,1)));
     OneType_3shocks(4,1)=0;
     OneType_3shocks=OneType_3shocks';
     
     
     
     
     figure(401);
     plot(t_talk_label,XElast_FF(1,t_talk)./XElast_FF(1,5),'r-',t_talk_label,OneType_3shocks(t_talk,4)./OneType_3shocks(5,4),'k-',t_talk_label,NoReallocation_3shocks(t_talk,4)./NoReallocation_3shocks(5,4),'c-',t_talk_label,XElast_Data./XElast_Data(end),'b','LineWidth',lw);
     h=title('Elasticity of Exports to RER');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.5:3.0);
     axis([t_talk_label(1)-0.01 t_talk_label(end) 0 1.05]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     h=ylabel('Export elasticity as % of final s.s.');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f )
     hleg = legend('Baseline','One type','No reallocation', 'Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);    
     
     if sav==1
 
         Filename = 'KLS2_XElasticity3_data_OneTypeNoReallocation2'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
     
     
     
%% Figure 5 in the Paper:  Export elasticity by export-intensity groups
 
    % Results from the model: 
    %Column 1 is Low X/Y FF model, Column 2 is High X/Y FF model 
    %Column 3 is Low X/Y NoFF model, Column 4 is High X/Y NoFF model
    
%     model_validation=[
% 0.000	0.000	0.000	0.000
% 0.994	0.332	0.975	0.362
% 0.720	0.519	0.564	0.567
% 0.389	0.328	0.354	0.356
% 0.337	0.291	0.277	0.281
% 0.348	0.297	0.286	0.291
% ];

  model_validation(:,1)=log(Mex_3shocks_FF(93,1:6))-log(Mex_3shocks_FF(93,1));
  model_validation(:,2)=log(Mex_3shocks_FF(94,1:6))-log(Mex_3shocks_FF(94,1));
  model_validation(:,3)=log(Mex_3shocks_NoFF(93,1:6))-log(Mex_3shocks_NoFF(93,1));
  model_validation(:,4)=log(Mex_3shocks_NoFF(94,1:6))-log(Mex_3shocks_NoFF(94,1));
  

  % Results from the data (PPI-Adjusted)
  % Column 1 is Low X/Y, Column 2 is High X/Y for a regression without controls
  % Column 3 is Low X/Y, Column 4 is High X/Y for a regression with
  % controls and industry FE
    

 data_validation=[
0.000	0.000	0.000	0.000
0.621	0.313	0.672	0.346
0.788	0.221	0.834	0.252
0.796	0.178	0.854	0.212
0.873	0.175	0.917	0.218
0.761	0.015	0.816	0.076
];
 
    
    figure(501);
    plot(t_talk2_label(1:6),model_validation(:,1),'k-',t_talk2_label(1:6),model_validation(:,2),'r--','LineWidth',lw);
    h=ylabel('Log change from s.s.');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );        
    h=title('Model');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );     
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:1:t_talk2_label(5)],'YTick',0:0.25:1.2);
    axis([t_talk2_label(1)-0.01 t_talk2_label(5) 0 1.05]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize',f);  
    hleg = legend('Low X/Y firms','High X/Y firms','Location','NorthEast');            
    thh = findobj(hleg,'type','text');
    set(thh,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_Validation_Model'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end   
    
    figure(502);
    plot(t_talk2_label(1:6),data_validation(:,3),'k-',t_talk2_label(1:6),data_validation(:,4),'r--','LineWidth',lw);
    h=ylabel('Log change from 1994');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );        
    h=title('Data');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );     
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:1:t_talk2_label(5)],'YTick',0:0.25:1);
    axis([t_talk2_label(1)-0.01 t_talk2_label(5) 0 1.05]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize',f);  
    hleg = legend('Low X/Y firms','High X/Y firms','Location','Northwest');            
    thh = findobj(hleg,'type','text');
    set(thh,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_Validation_Data'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end  
    

    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONLINE APPENDIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%% Figure 3 in the Online Appendix: Shocks

    shock_baseline_A = Mex_3shocks_FF(97,1:5);
    shock_baseline_pm = Mex_3shocks_FF(96,1:5);
    shock_baseline_r = Mex_3shocks_FF(98,2:6);  % lagged one period
    
    shock_noff_A = Mex_3shocks_NoFF(97,1:5);
    shock_noff_pm = Mex_3shocks_NoFF(96,1:5);
    shock_noff_r = Mex_3shocks_NoFF(98,2:6);  % lagged one period
    
    
    figure(1301);
    subplot(3,1,1); plot(t_talk_label,shock_baseline_pm-shock_baseline_pm(1),'r-',t_talk_label,shock_noff_pm-shock_noff_pm(1),'k-','LineWidth',lw);    
    h=title({'Price of imports, p_m'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.5:0.25:2);
    axis([t_talk_label(1) t_talk_label(end) -0.5 0.05]);
    grid on;
    hleg = legend('Baseline','Frictionless','Location','SouthEast','Orientation','horizontal');
    set(legend, 'Position', [0.515 0.03 0.0001 0.0001])
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    
        
    subplot(3,1,2); plot(t_talk_label,shock_baseline_A-shock_baseline_A(1),'r-',t_talk_label,shock_noff_A-shock_noff_A(1),'k-','LineWidth',lw);    
    h=title({'Aggregate productivity, A'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.05:0.05:2);
    axis([t_talk_label(1) t_talk_label(end) -0.05 0.05]);
    grid on;
        
    subplot(3,1,3); plot(t_talk_label,shock_baseline_r-shock_baseline_r(1),'r-',t_talk_label,shock_noff_r-shock_noff_r(1),'k-','LineWidth',lw);    
    h=title({'Real interest rate, r'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.05:2);
    axis([t_talk_label(1) t_talk_label(end) 0 0.12]);
    grid on;

    if sav==1      
        saveas(gcf,'KLS2_Shocks2', 'epsc');
        saveas(gcf,'KLS2_Shocks2', 'fig');
    end
     
%% Figure 4 in the Online Appendix: Shocks decomposition
    % RER, RGDP, inv/GDP log changes wrt SS (for inv/GDP is changes)
        
    Shocks3(:,1)= log(Mex_3shocks_FF(50,1:5))-log(Mex_3shocks_FF(50,1)); 
    Shocks3(:,2)= log(Mex_3shocks_FF(90,1:5))-log(Mex_3shocks_FF(90,1)); 
    Shocks3(:,3)= Mex_3shocks_FF(41,1:5)./Mex_3shocks_FF(73,1:5)-Mex_3shocks_FF(41,1)./Mex_3shocks_FF(73,1); 
        
    PmShock(:,1)= log(Mex_3shocks_FF_OnlyPmShock(50,1:5))-log(Mex_3shocks_FF_OnlyPmShock(50,1)); 
    PmShock(:,2)= log(Mex_3shocks_FF_OnlyPmShock(90,1:5))-log(Mex_3shocks_FF_OnlyPmShock(90,1)); 
    PmShock(:,3)= Mex_3shocks_FF_OnlyPmShock(41,1:5)./Mex_3shocks_FF_OnlyPmShock(73,1:5)-Mex_3shocks_FF_OnlyPmShock(41,1)./Mex_3shocks_FF_OnlyPmShock(73,1); 
          
    PmZShock(:,1)= log(Mex_3shocks_FF_PmZShocks(50,1:5))-log(Mex_3shocks_FF_PmZShocks(50,1)); 
    PmZShock(:,2)= log(Mex_3shocks_FF_PmZShocks(90,1:5))-log(Mex_3shocks_FF_PmZShocks(90,1)); 
    PmZShock(:,3)= Mex_3shocks_FF_PmZShocks(41,1:5)./Mex_3shocks_FF_PmZShocks(73,1:5)-Mex_3shocks_FF_PmZShocks(41,1)./Mex_3shocks_FF_PmZShocks(73,1); 

    PmRShock(:,1)= log(Mex_3shocks_FF_PmRShocks(50,1:5))-log(Mex_3shocks_FF_PmRShocks(50,1)); 
    PmRShock(:,2)= log(Mex_3shocks_FF_PmRShocks(90,1:5))-log(Mex_3shocks_FF_PmRShocks(90,1)); 
    PmRShock(:,3)= Mex_3shocks_FF_PmRShocks(41,1:5)./Mex_3shocks_FF_PmRShocks(73,1:5)-Mex_3shocks_FF_PmRShocks(41,1)./Mex_3shocks_FF_PmRShocks(73,1); 
    
        
    figure(1401);
    subplot(3,1,1); plot(t_talk_label,Shocks3(:,1),'r-',t_talk_label,PmShock(:,1),'k-',t_talk_label,PmRShock(:,1),'b-','LineWidth',lw);    
    h=title({'Real exchange rate'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:2);
    axis([t_talk_label(1) t_talk_label(end) 0 0.45]);
    grid on;
    hleg = legend('Baseline','Only p_m shocks','Only p_m and r shocks','Location','SouthEast','Orientation','horizontal');
    set(legend, 'Position', [0.515 0.04 0.00001 0.00001])
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    subplot(3,1,2); plot(t_talk_label,Shocks3(:,2),'r-',t_talk_label,PmShock(:,2),'k-',t_talk_label,PmRShock(:,2),'b-','LineWidth',lw);    
    h=title({'Real GDP'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:2);
    axis([t_talk_label(1) t_talk_label(end) -0.1 0.1]);
    grid on;

    subplot(3,1,3); plot(t_talk_label,Shocks3(:,3),'r-',t_talk_label,PmShock(:,3),'k-',t_talk_label,PmRShock(:,3),'b-','LineWidth',lw);    
    h=title({'Investment / GDP'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.3:0.15:2);
    axis([t_talk_label(1) t_talk_label(end) -0.25 0.4]);
    grid on;

    if sav==1    
        saveas(gcf,'KLS2_Shocks4', 'epsc');
        saveas(gcf,'KLS2_Shocks4', 'fig');      
    end
%% Figure 5 in the Online Appendix: No FF model with same shocks as in FF model 
    %m.pm_v, m.zagg_v, m.rv N=10	
    ShocksFF_NoFF1=Mex_3shocks_NoFF_BaselineShocks(96:98,:)';

    % RER, Inv/output, RGDP, XElast
    ShocksFF_NoFF1_m(:,1)=Mex_3shocks_NoFF_BaselineShocks(50,1:10);
    ShocksFF_NoFF1_m(:,2)=Mex_3shocks_NoFF_BaselineShocks(41,1:10)./Mex_3shocks_NoFF_BaselineShocks(73,1:10);
    ShocksFF_NoFF1_m(:,3)=Mex_3shocks_NoFF_BaselineShocks(90,1:10);
    ShocksFF_NoFF1_m(:,4)=(log(Mex_3shocks_NoFF_BaselineShocks(91,1:10))-log(Mex_3shocks_NoFF_BaselineShocks(91,1)))./(log(Mex_3shocks_NoFF_BaselineShocks(50,1:10))-log(Mex_3shocks_NoFF_BaselineShocks(50,1)));
    ShocksFF_NoFF1_m(1,4)=0;
    
    
    
    
  
    %RER
    figure(1501);
    plot(t_talk_label,log(Mex_3shocks_FF(50,t_talk))-log(Mex_3shocks_FF(50,1)),'r-',t_talk_label,log(ShocksFF_NoFF1_m(t_talk,1))-log(ShocksFF_NoFF1_m(1,1)),'k-',t_talk_label,data_RER,'b','LineWidth',lw);    
    h=title({'Real Exchange Rate';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:0.5);
    axis([t_talk_label(1) t_talk_label(end) 0 0.45]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_RER_FFshocks'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end

    %Real GDP    
    figure(1502);
    plot(t_talk_label,log(Mex_3shocks_FF(90,t_talk))-log(Mex_3shocks_FF(90,1)),'r-',t_talk_label,log(ShocksFF_NoFF1_m(t_talk,3))-log(ShocksFF_NoFF1_m(1,3)),'k-',t_talk_label,data_GDP,'b','LineWidth',lw);    
    h=title({'Real GDP';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.1 0.1]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_GDP_FFshocks'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
    
    % Investment_Output
    Inv_GDP_FF = (Mex_3shocks_FF(41,:)./Mex_3shocks_FF(73,:));

    figure(1503);
    plot(t_talk_label,Inv_GDP_FF(1,t_talk)-Inv_GDP_FF(1,1),'r-',t_talk_label,ShocksFF_NoFF1_m(t_talk,2)-ShocksFF_NoFF1_m(1,2),'k-',t_talk_label,data_Inv,'b','LineWidth',lw);    
    h=title({'Investment over GDP';'(change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.1 0.06]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);

    if sav==1

        Filename = 'KLS2_Inv_FFshocks'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end

     
     % Exports elasticity

     %	91	 sim_fun.X_Laspeyres
     %	50	 sim_fun.xi_real

     XElast_FF = log(Mex_3shocks_FF(91,:)./Mex_3shocks_FF(91,1))./log(Mex_3shocks_FF(50,:)./Mex_3shocks_FF(50,1));
     XElast_Data =data_Xreal./data_RER;
     XElast_Data(1) = 0;
     XElast_FF(1)=0;
 
     
     % Data vs SMM-Trans, 5 periods   
     figure(1504);
     plot(t_talk_label,XElast_FF(1,t_talk),'r-',t_talk_label,ShocksFF_NoFF1_m(t_talk,4),'k-',t_talk_label,XElast_Data,'b','LineWidth',lw);
     h=title('Elasticity of Exports to RER');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:1:3.0);
     axis([t_talk_label(1)-0.01 t_talk_label(end) 0 3]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Baseline','Frictionless','Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);    
     
     if sav==1
 
         Filename = 'KLS2_XElasticity3_data_FFshocks'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
  
%% Figure 6 in the Online Appendix:  Current account  (trade balance)in terms of GDP for both data and model
   
    data_tb = [-0.026240473;0.042354178;0.027330479;0.011417082;-0.002782437]';
    % 1994, 1995, and so on...
    %(Net exports of goods and services/GDP, current USD, World Bank data)

    baseline_tb = Mex_3shocks_FF(75,1:5)';
    noff_tb = Mex_3shocks_NoFF(75,1:5)';
    
    


    figure(1601);
    plot(t_talk_label,baseline_tb,'r-',t_talk_label,data_tb,'k-','LineWidth',lw);    
    h=title({'Trade balance in terms of GDP'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.05:0.05:2);
    axis([t_talk_label(1) t_talk_label(end) -0.08 0.08]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Data','Location','SouthEast','Orientation','horizontal');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    


    if sav==1
        saveas(gcf,'KLS2_tb', 'epsc');
        saveas(gcf,'KLS2_tb', 'fig');
    end
          
%% Figure 7 in the Online Appendix: No Reallocation
      
    % RER, Inv/output, RGDP, XElast
    NoReallocation_3shocks(:,1)=Mex_3shocks_NoReallocation(50,1:10);
    NoReallocation_3shocks(:,2)=Mex_3shocks_NoReallocation(41,1:10)./Mex_3shocks_NoReallocation(73,1:10);
    NoReallocation_3shocks(:,3)=Mex_3shocks_NoReallocation(90,1:10);
    NoReallocation_3shocks(:,4)=(log(Mex_3shocks_NoReallocation(91,1:10))-log(Mex_3shocks_NoReallocation(91,1)))./(log(Mex_3shocks_NoReallocation(50,1:10))-log(Mex_3shocks_NoReallocation(50,1)));
    NoReallocation_3shocks(1,4)=0;
    
        
    %RER
    figure(1701);
    plot(t_talk_label,log(Mex_3shocks_FF(50,t_talk))-log(Mex_3shocks_FF(50,1)),'r-',t_talk_label,log(NoReallocation_3shocks(t_talk,1))-log(NoReallocation_3shocks(1,1)),'k-',t_talk_label,data_RER,'b','LineWidth',lw);    
    h=title({'Real Exchange Rate';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:0.5);
    axis([t_talk_label(1) t_talk_label(end) 0 0.45]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','No Reallocation','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_RER_NoReallocation'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end


    %Real GDP    
    figure(1702);
    plot(t_talk_label,log(Mex_3shocks_FF(90,t_talk))-log(Mex_3shocks_FF(90,1)),'r-',t_talk_label,log(NoReallocation_3shocks(t_talk,3))-log(NoReallocation_3shocks(1,3)),'k-',t_talk_label,data_GDP,'b','LineWidth',lw);    
    h=title({'Real GDP';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.1 0.1]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','No Reallocation','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_GDP_NoReallocation'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
     
     
     % Investment_Output

     
     figure(1703);
     plot(t_talk_label,Inv_GDP_FF(1,t_talk)-Inv_GDP_FF(1,1),'r-',t_talk_label,NoReallocation_3shocks(t_talk,2)-NoReallocation_3shocks(1,2),'k-',t_talk_label,data_Inv,'b','LineWidth',lw);    
     h=title({'Investment over GDP';'(change from steady state)'});
     to = findobj( h, 'type', 'text' );
     set(to,'fontsize', f );    
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
     axis([t_talk_label(1) t_talk_label(end) -0.1 0.06]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Baseline','No Reallocation','Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);
     
     if sav==1
 
         Filename = 'KLS2_Inv_NoReallocation'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
     
   
     
    % Exports elasticity
 
        
     figure(1704);
     plot(t_talk_label,XElast_FF(1,t_talk)./XElast_FF(1,5),'r-',t_talk_label,NoReallocation_3shocks(t_talk,4)./NoReallocation_3shocks(5,4),'k-',t_talk_label,XElast_Data./XElast_Data(end),'b','LineWidth',lw);
     h=title('Elasticity of Exports to RER');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.25:3.0);
     axis([t_talk_label(1)-0.01 t_talk_label(end) 0 1.1]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Baseline','No Reallocation','Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);    
     
     if sav==1
 
         Filename = 'KLS2_XElasticity3_data_NoReallocation2'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
          
     
%% Figure 8 in the Online Appendix: One Type
     
    % RER, Inv/output, RGDP, XElast
    OneType_3shocks(:,1)=Mex_3shocks_OneType(50,1:10);
    OneType_3shocks(:,2)=Mex_3shocks_OneType(41,1:10)./Mex_3shocks_OneType(73,1:10);
    OneType_3shocks(:,3)=Mex_3shocks_OneType(90,1:10);
    OneType_3shocks(:,4)=(log(Mex_3shocks_OneType(91,1:10))-log(Mex_3shocks_OneType(91,1)))./(log(Mex_3shocks_OneType(50,1:10))-log(Mex_3shocks_OneType(50,1)));
    OneType_3shocks(1,4)=0;
%     OneType_3shocks=[1.585947751	2.263110053	2.084623696	1.835771635	1.805946763	1.809778704	1.811452456	1.812173319	1.812631647	1.812710078
%                     0.075139093	0.031795152	0.055646837	0.067583348	0.086273283	0.07927159	0.075566371	0.073751253	0.073010892	0.072099623
%                     2.183508655	2.085541497	2.10915021	2.169678528	2.213981544	2.222532068	2.227177468	2.230050086	2.233362808	2.234916232
%                     0	2.275312163	2.353749237	2.454046395	2.456785664	2.50514961	2.521982603	2.526264714	2.527174445	2.524656404]';

    %RER
    figure(1801);
    plot(t_talk_label,log(Mex_3shocks_FF(50,t_talk))-log(Mex_3shocks_FF(50,1)),'r-',t_talk_label,log(OneType_3shocks(t_talk,1))-log(OneType_3shocks(1,1)),'k-',t_talk_label,data_RER,'b','LineWidth',lw);    
    h=title({'Real Exchange Rate';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:0.5);
    axis([t_talk_label(1) t_talk_label(end) 0 0.45]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','One Type','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_RER_OneType'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end


    %Real GDP    
    figure(1802);
    plot(t_talk_label,log(Mex_3shocks_FF(90,t_talk))-log(Mex_3shocks_FF(90,1)),'r-',t_talk_label,log(OneType_3shocks(t_talk,3))-log(OneType_3shocks(1,3)),'k-',t_talk_label,data_GDP,'b','LineWidth',lw);    
    h=title({'Real GDP';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.1 0.1]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','One Type','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_GDP_OneType'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
     
     
     %Investment_Output

     
     figure(1803);
     plot(t_talk_label,Inv_GDP_FF(1,t_talk)-Inv_GDP_FF(1,1),'r-',t_talk_label,OneType_3shocks(t_talk,2)-OneType_3shocks(1,2),'k-',t_talk_label,data_Inv,'b','LineWidth',lw);    
     h=title({'Investment over GDP';'(change from steady state)'});
     to = findobj( h, 'type', 'text' );
     set(to,'fontsize', f );    
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
     axis([t_talk_label(1) t_talk_label(end) -0.1 0.06]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Baseline','One Type','Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);
     
     if sav==1
 
         Filename = 'KLS2_Inv_OneType'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
     
   
     
    % Exports elasticity
 
       
     figure(1804);
     plot(t_talk_label,XElast_FF(1,t_talk)./XElast_FF(1,5),'r-',t_talk_label,OneType_3shocks(t_talk,4)./OneType_3shocks(5,4),'k-',t_talk_label,XElast_Data./XElast_Data(end),'b','LineWidth',lw);
     h=title('Elasticity of Exports to RER');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.25:3.0);
     axis([t_talk_label(1)-0.01 t_talk_label(end) 0 1.1]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Baseline','One Type','Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);    
     
     if sav==1
 
         Filename = 'KLS2_XElasticity3_data_OneType2'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
     
%% Figure 9 in the Online Appendix: Firm exports growth by X/Sales and 
    % imported intermediates/intermediates
    
    % Data
    % LowX,LowM	HighX,LowM	LowX,HighM	HighX,HighM
	Xelast_Mint=[ 0.00000000	0.00000000	0.00000000	0.00000000
        0.63679564	0.31719444	0.63679564	0.31719444
        0.90274144	0.49658614	0.92277384	0.51661854
        0.96034254	0.47691664	0.90750874	0.42408284
        0.94386684	0.53764624	1.04610784	0.63988724
        ];

     figure(1901);
     plot(1994:1:1998,Xelast_Mint(t_talk,1),'b-',1994:1:1998,Xelast_Mint(t_talk,2),'c-',1994:1:1998,Xelast_Mint(t_talk,3),'r-',1994:1:1998,Xelast_Mint(t_talk,4),'k-','LineWidth',lw);
     h=title('Elasticity of Exports to RER');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     set(gca,'FontSize',f,'XTick',1994:1:1998,'YTick',0:0.2:1.2);
     axis([1994-0.01 1998 0 1.2]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Low X, Low M','High X, Low M','Low X, High M','High X, High M','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);    
     
     if sav==1
 
         Filename = 'KLS2_XElasticity_ImportIntensity1'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end	

    
    

%% Figure 10 in the Online Appendix: Firm exports growth by X/Sales and 
    % imported intermediates/sales
    
     % Data
     % LowX,LowM	HighX,LowM	LowX,HighM	HighX,HighM
     Xelast_Mint2=[ 0.00000000	0.00000000	0.00000000	0.00000000
        0.67026435	0.34668065	0.67026435	0.34668065
        0.84695945	0.50349425	0.94183605	0.59837085
        0.90543725	0.48431295	0.92283605	0.50171175
        0.92120045	0.54406105	1.03313705	0.65599765];

     figure(11001);
     plot(1994:1:1998,Xelast_Mint2(t_talk,1),'b-',1994:1:1998,Xelast_Mint2(t_talk,2),'c-',1994:1:1998,Xelast_Mint2(t_talk,3),'r-',1994:1:1998,Xelast_Mint2(t_talk,4),'k-','LineWidth',lw);
     h=title('Elasticity of Exports to RER');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     set(gca,'FontSize',f,'XTick',1994:1:1998,'YTick',0:0.2:1.2);
     axis([1994-0.01 1998 0 1.2]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Low X, Low M','High X, Low M','Low X, High M','High X, High M','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);    
     
     if sav==1
 
         Filename = 'KLS2_XElasticity_ImportIntensity2'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end	
    
%% Figure 11 in the Online Appendix: External Finance Dependence and Exports Growth (real, in domestic
    % goods)
    
    % Data
    data_EFD=[0.14	nan	0.08	-0.45	0.40	0.03	-0.14	-0.08	0.28	0.24	0.18	0.20	0.25	0.22	0.33	0.23	1.14	-0.15	0.53	0.06	0.09	0.01	0.24	0.45	0.77	0.31	0.96	0.47
            0.56	0.48	0.57	0.67	0.79	1.00	0.64	0.47	0.66	0.48	1.39	0.46	0.71	0.48	0.49	0.76	0.27	0.38	0.38	1.08	1.18	0.99	0.93	0.56	0.48	0.60	0.12	0.57
            0.50	0.29	0.73	0.62	0.97	1.64	0.75	0.56	0.97	1.23	0.73	0.50	-0.09	0.89	0.57	0.98	0.41	0.51	0.20	0.81	0.40	0.56	1.15	1.22	0.62	0.68	1.46	0.73
            559	175	146	7	376	201	92	97	82	95	188	196	180	313	37	85	361	53	51	237	74	38	381	252	186	225	39	53
            ]';

    [b1,b1int,r1,r1int,stats1] = regress(data_EFD(:,2),[ones(28,1) data_EFD(:,1)]);
    [b2,b2int,r2,r2int,stats2] = regress(data_EFD(:,3),[ones(28,1) data_EFD(:,1)]);
    b1=round(b1*100)/100;
    b2=round(b2*100)/100;
    stats1=round(stats1*100)/100;
    stats2=round(stats2*100)/100;
    
    figure(11101);
    scatter(data_EFD(:,1),data_EFD(:,2),data_EFD(:,4)./5,'k','filled');
    h=title('Real Exports 1-Yr Growth');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    h=xlabel('External Finance Dependence');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hold on;
    line([min(data_EFD(:,1)),max(data_EFD(:,1))],[min(data_EFD(:,1))*b1(2)+b1(1),max(data_EFD(:,1))*b1(2)+b1(1)],'color','k','linewidth',lw)
    h=text(0.68,0.85,['y=',num2str(b1(1)),num2str(b1(2)),'x']);
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    
    
    h=text(0.68,0.7,['R^2=',num2str(stats1(1))]);
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',-0.5:0.5:1.5,'YTick',0:0.4:2);
    axis([-0.5 1.5 -0.2 1.8]);
        
    if sav==1

        Filename = 'KLS2_EFD_1'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end

    figure(11102);
    scatter(data_EFD(:,1),data_EFD(:,3),data_EFD(:,4)./5,'k','filled');
     h=title('Real Exports 5-Yr Growth');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    h=xlabel('External Finance Dependence');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hold on;
    line([min(data_EFD(:,1)),max(data_EFD(:,1))],[min(data_EFD(:,1))*b2(2)+b2(1),max(data_EFD(:,1))*b2(2)+b2(1)],'color','k','linewidth',lw)
    h=text(0.7,1.1,['y=',num2str(b2(1)),'+',num2str(b2(2)),'x']);
     to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    
    h=text(0.7,0.95,['R^2=',num2str(stats2(1))]);
     to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',-0.5:0.5:1.5,'YTick',0:0.4:2);
    axis([-0.5 1.5 -0.2 1.8]);
        
    
    if sav==1

        Filename = 'KLS2_EFD_2'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end
   
   
%% Figure 12 in the Online Appendix: XElasticity in IRP model
        
    % 3 Shocks, IRP lambda 1, calibrated shocks
    Xelast_IRP=(log(Mex_3shocks_FF_IRP(91,1:10))-log(Mex_3shocks_FF_IRP(91,1)))./(log(Mex_3shocks_FF_IRP(50,1:10))-log(Mex_3shocks_FF_IRP(50,1)));
    Xelast_IRP(1)=0;
    Xelast_IRP=Xelast_IRP';

    figure(11201);
    plot(t_talk_label,XElast_FF(1,t_talk),'r-',t_talk_label,Xelast_IRP(t_talk,1),'k-',t_talk_label,XElast_Data,'b','LineWidth',lw);
    h=title('Elasticity of Exports to RER');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:1:3.0);
    axis([t_talk_label(1)-0.01 t_talk_label(end) 0 3]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','IRP','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_XElasticity3_data_IRP2'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end
        
     
%% Figure 13 in the Online Appendix: Distribution of foreign debt

    % Exports elasticity:

    %	91	 sim_fun.X_Laspeyres
    %	50	 sim_fun.xi_real

    XElast_Pm_FF = log(Mex_Pm_FF(91,:)./Mex_Pm_FF(91,1))./log(Mex_Pm_FF(50,:)./Mex_Pm_FF(50,1));
    XElast_Pm_FF(1)=0;
    XElast_Pm_lambda0 = log(Mex_Pm_lambda0(91,:)./Mex_Pm_lambda0(91,1))./log(Mex_Pm_lambda0(50,:)./Mex_Pm_lambda0(50,1));
    XElast_Pm_lambda0(1)=0;
    XElast_Pm_lambda1 = log(Mex_Pm_lambda1(91,:)./Mex_Pm_lambda1(91,1))./log(Mex_Pm_lambda1(50,:)./Mex_Pm_lambda1(50,1));
    XElast_Pm_lambda1(1)=0;
    XElast_Pm_lambda05lambdax0 = log(Mex_Pm_lambda05lambdax0(91,:)./Mex_Pm_lambda05lambdax0(91,1))./log(Mex_Pm_lambda05lambdax0(50,:)./Mex_Pm_lambda05lambdax0(50,1));
    XElast_Pm_lambda05lambdax0(1)=0;


    figure(11301);
    plot( t_talk2_label,XElast_Pm_lambda0(1,t_talk2)/XElast_Pm_lambda0(1,end),'k-',t_talk2_label,XElast_Pm_FF(1,t_talk2)/XElast_Pm_FF(1,end),'r-',t_talk2_label,XElast_Pm_lambda05lambdax0(1,t_talk2)/XElast_Pm_lambda05lambdax0(1,end),'c-',t_talk2_label,XElast_Pm_lambda1(1,t_talk2)/XElast_Pm_lambda1(1,end),'b-','LineWidth',lw);
    h=title('Elasticity of Exports to RER');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:2:t_talk2_label(end)],'YTick',0:0.2:2.6);
    axis([t_talk2_label(1)-0.01 t_talk2_label(end) 0 1.1]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );    
    h=ylabel('Export elasticity as % of final s.s.');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );      
    hleg = legend('All foreign debt','Baseline','Foreign-debt for high export intensity','No foreign debt ','Location','SouthEast');            
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_XElasticity_foreigndebt_3c'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
              
    end
    
    

    figure(11302);
    plot(t_talk2_label,100*(log(Mex_Pm_lambda0(95,t_talk2))-log(Mex_Pm_lambda0(95,1))),'k-',t_talk2_label,100*(log(Mex_Pm_FF(95,t_talk2))-log(Mex_Pm_FF(95,1))),'r-',t_talk2_label,100*(log(Mex_Pm_lambda05lambdax0(95,t_talk2))-log(Mex_Pm_lambda05lambdax0(95,1))),'c-', t_talk2_label,100*(log(Mex_Pm_lambda1(95,t_talk2))-log(Mex_Pm_lambda1(95,1))),'b-','LineWidth',lw);
    h=ylabel('% deviation from s.s.');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );        
    h=title('Real wage');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );     
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:2:t_talk2_label(end)],'YTick',0:5:25);
    axis([t_talk2_label(1)-0.01 t_talk2_label(end) 0 25]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize',f);    
    hleg = legend('All foreign debt','Baseline','Foreign-debt for high export intensity','No foreign debt ','Location','SouthEast');            
    thh = findobj(hleg,'type','text');
    set(thh,'fontsize',f);
    if sav==1

        Filename = 'KLS2_XElasticity_foreigndebt_explanation_1'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end    
    
    figure(11303);
    plot(t_talk2_label,100*(log(Mex_Pm_lambda0(42,t_talk2))-log(Mex_Pm_lambda0(42,1))),'k-',t_talk2_label,100*(log(Mex_Pm_FF(42,t_talk2))-log(Mex_Pm_FF(42,1))),'r-',t_talk2_label,100*(log(Mex_Pm_lambda05lambdax0(42,t_talk2))-log(Mex_Pm_lambda05lambdax0(42,1))),'c-',t_talk2_label,100*(log(Mex_Pm_lambda1(42,t_talk2))-log(Mex_Pm_lambda1(42,1))),'b-','LineWidth',lw);   
    h=title('Labor of firms with high export intensity');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );   
    h=ylabel('% deviation from s.s.');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',[-1 t_talk2_label(1)+1:2:t_talk2_label(end)],'YTick',0:15:60);
    axis([t_talk2_label(1)-0.01 t_talk2_label(end) 0 60]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize',f);    
     hleg = legend('All foreign debt','Baseline','Foreign-debt for high export intensity','No foreign debt ','Location','SouthEast');            
    thh = findobj(hleg,'type','text');
    set(thh,'fontsize',f);
    
    
    if sav==1

        Filename = 'KLS2_XElasticity_foreigndebt_explanation_2'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end
         
     
%% Figure 14 in the Online Appendix: Financial Crisis Shock 
      
    %m.pm_v, m.zagg_v, m.theta N=10	
    Theta_3shocks=Mex_3shocks_FF_FinancialCrisis(96:97,1:10)';
    Theta_3shocks(1:10,3)=[  0.49375	0.49375	0.27558498	0.22725906	0.30416216	0.42243861	0.42243861	0.42243861	0.42243861	0.42243861]';


      % RER, Inv/output, RGDP, XElast
     Theta_3shocks_m(1,1:10)=Mex_3shocks_FF_FinancialCrisis(50,1:10);
     Theta_3shocks_m(2,1:10)=Mex_3shocks_FF_FinancialCrisis(41,1:10)./Mex_3shocks_FF_FinancialCrisis(73,1:10);
     Theta_3shocks_m(3,1:10)=Mex_3shocks_FF_FinancialCrisis(90,1:10);
     Theta_3shocks_m(4,1:10)=(log(Mex_3shocks_FF_FinancialCrisis(91,1:10))-log(Mex_3shocks_FF_FinancialCrisis(91,1)))./(log(Mex_3shocks_FF_FinancialCrisis(50,1:10))-log(Mex_3shocks_FF_FinancialCrisis(50,1)));
     Theta_3shocks_m(4,1)=0;
     Theta_3shocks_m=Theta_3shocks_m';                
                    

    figure(11401);
    subplot(3,1,1); plot(t_talk_label,shock_baseline_pm,'r-',t_talk_label,Theta_3shocks(1:5,1),'k-','LineWidth',lw);    
    h=title({'Price of imports, p_m'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.25:2);
    axis([t_talk_label(1) t_talk_label(end) 0.5 1.05]);
    grid on;
    hleg = legend('Baseline','Financial Crisis','Location','SouthEast','Orientation','horizontal');
    set(legend, 'Position', [0.515 0.03 0.0001 0.0001])
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    subplot(3,1,2); plot(t_talk_label,shock_baseline_A,'r-',t_talk_label,Theta_3shocks(1:5,2),'k-','LineWidth',lw);    
    h=title({'Aggregate productivity, A'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0.8:0.1:2);
    axis([t_talk_label(1) t_talk_label(end) 0.8 1.1]);
    grid on;
    set(th,'fontsize',f);      

    subplot(3,1,3); plot(t_talk_label,0.4938*ones(5,1),'r-',t_talk_label,Theta_3shocks(1:5,3),'k-','LineWidth',lw);    
    h=title({'Collateral, \theta'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:2);
    axis([t_talk_label(1) t_talk_label(end) 0.201 0.6]);
    grid on;

    if sav==1
        saveas(gcf,'KLS2_Shocks_theta', 'epsc');
        saveas(gcf,'KLS2_Shocks_theta', 'fig');
    end
    
%% Figure 15 in the Online Appendix: Financial Crisis Shock     

    % RER
    figure(11501);
    plot(t_talk_label,log(Mex_3shocks_FF(50,t_talk))-log(Mex_3shocks_FF(50,1)),'r-',t_talk_label,log( Theta_3shocks_m(t_talk,1))-log( Theta_3shocks_m(1,1)),'k-',t_talk_label,data_RER,'b','LineWidth',lw);    
    h=title({'Real Exchange Rate';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:0.1:0.5);
    axis([t_talk_label(1) t_talk_label(end) 0 0.45]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Financial Crisis','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_RER_Theta'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end


    % Real GDP    
    figure(11502);
    plot(t_talk_label,log(Mex_3shocks_FF(90,t_talk))-log(Mex_3shocks_FF(90,1)),'r-',t_talk_label,log( Theta_3shocks_m(t_talk,3))-log( Theta_3shocks_m(1,3)),'k-',t_talk_label,data_GDP,'b','LineWidth',lw);    
    h=title({'Real GDP';'(log change from steady state)'});
    to = findobj( h, 'type', 'text' );
    set(to,'fontsize', f );    
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.2:0.1:0.1);
    axis([t_talk_label(1) t_talk_label(end) -0.25 0.1]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Financial Crisis','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);
    
    if sav==1

        Filename = 'KLS2_GDP_Theta'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
    
     % Investment_Output
     figure(11503);
     plot(t_talk_label,Inv_GDP_FF(1,t_talk)-Inv_GDP_FF(1,1),'r-',t_talk_label, Theta_3shocks_m(t_talk,2)- Theta_3shocks_m(1,2),'k-',t_talk_label,data_Inv,'b','LineWidth',lw);    
     h=title({'Investment over GDP';'(change from steady state)'});
     to = findobj( h, 'type', 'text' );
     set(to,'fontsize', f );    
     set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',-0.1:0.05:0.1);
     axis([t_talk_label(1) t_talk_label(end) -0.1 0.06]);
     grid on;
     h=xlabel('Years');
     to = findobj( h, 'type', 'text' );
     set( to, 'fontsize', f );
     hleg = legend('Baseline','Financial Crisis','Data','Location','SouthEast');
     th = findobj(hleg,'type','text');
     set(th,'fontsize',f);
     
     if sav==1
 
         Filename = 'KLS2_Inv_Theta'; 
         saveas(gcf,Filename, 'epsc');
         saveas(gcf,Filename, 'fig');
                 
     end
     
  
    % Exports elasticity
 
    figure(11504);
    plot(t_talk_label,XElast_FF(1,t_talk),'r-',t_talk_label, Theta_3shocks_m(t_talk,4),'k-',t_talk_label,XElast_Data,'b','LineWidth',lw);
    h=title('Elasticity of Exports to RER');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:1:3.0);
    axis([t_talk_label(1)-0.01 t_talk_label(end) 0 3]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Baseline','Financial Crisis','Data','Location','SouthEast');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_XElasticity3_data_Theta'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');

    end
    
    
    
%% Figure 16 in the Online Appendix : Debt held by more productive firms

    XElast_UIPdev = log(Mex_UIPdeviations(91,:)./Mex_UIPdeviations(91,1))./log(Mex_UIPdeviations(50,:)./Mex_UIPdeviations(50,1));
    XElast_UIPdev(1)=0;
    
    figure(11601);
    p = plot(t_talk_label',[ XElast_Data' XElast_FF(1:nT)' XElast_UIPdev(1:nT)'],'LineWidth',lw); grid on; 
    set(p(3),'Color','k')
    h=title('Elasticity of Exports to RER');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:1:3.0);
    axis([t_talk_label(1)-0.01 t_talk_label(end) 0 3]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Data','Baseline','Endogenous Debt Choice','Location','NorthWest');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_UIPdeviation_Elasticity'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end

%% Figure 17 in the Online Appendix: Optimal choice of debt
     
    XElast_zFD = log(Mex_high_z_ForDebt(91,:)./Mex_high_z_ForDebt(91,1))./log(Mex_high_z_ForDebt(50,:)./Mex_high_z_ForDebt(50,1));
    XElast_zFD(1)=0;
    
    
    figure(11701);
    p = plot(t_talk_label',[ XElast_Data' XElast_FF(1:nT)' XElast_zFD(1:nT)'],'LineWidth',lw); grid on; 
    set(p(3),'Color','k')
    h=title('Elasticity of Exports to RER');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    set(gca,'FontSize',f,'XTick',t_talk_label(1):1:t_talk_label(end),'YTick',0:1:3.0);
    axis([t_talk_label(1)-0.01 t_talk_label(end) 0 3]);
    grid on;
    h=xlabel('Years');
    to = findobj( h, 'type', 'text' );
    set( to, 'fontsize', f );
    hleg = legend('Data','Baseline','\lambda=0 for high productivity','Location','NorthWest');
    th = findobj(hleg,'type','text');
    set(th,'fontsize',f);    

    if sav==1

        Filename = 'KLS2_FD_lambda_045_0_Elasticity'; 
        saveas(gcf,Filename, 'epsc');
        saveas(gcf,Filename, 'fig');
                
    end
     


         
       


     
     