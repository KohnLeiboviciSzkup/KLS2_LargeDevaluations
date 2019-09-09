//========================================================================//
//initialize 

clear
pause on
set more off

//========================================================================//
//control panel 

//INEGI
	local INEGI = 1
	
//options
	local balanced 		 = 1
	local ppi_adjustment = 1
	
//years
	local year_start 	= 1994
	local year_end 		= 1999
	local finalSS_start = 1999
	local finalSS_end 	= 1999
	
//log file
	local logfile 		= 1
	local logfile_name 	= "KLS2_2019_02_22"

//========================================================================//
//directories and log file

if `INEGI'==0 {
	local directory_input 	= "D:\Dropbox\Shared\Large devaluations\JIErevision\Data\INEGI_EIA19942003\stata\19_02_21"	
	local directory_output 	= "D:\Dropbox\Shared\Large devaluations\JIErevision\Data\INEGI_EIA19942003\stata\19_02_21"	

	cd "`directory_output'"
	
	//use "`directory_input'\ejem_eia_94_03.dta", clear
	use "`directory_input'\ejem_eia_94_03_ext20170120.dta", clear
}
else if `INEGI'==1 {
	local directory_input = "C:\Local_D\USUARIOS_FF\Fernando_Leibovici_LM_356"
	local directory_output ="C:\Local_D\USUARIOS_FF\Fernando_Leibovici_LM_356"

	use "`directory_input'\lad_eia_1994_2003_f.dta", clear
}

if `logfile'==1 {
	log using `logfile_name'.log, replace
}

//========================================================================//
//rename variables

ren per year
ren clase industry_id

if `INEGI'==0 {
	ren control plant_id //called "folio" in the data description spreadsheet
}
else if `INEGI'==1 {
	ren IDE_EIM plant_id
}

//destring
	destring year		, replace
	destring industry_id, replace
	destring plant_id	, replace

//adjust year
if `INEGI'==0 {
	gen year_temp		= 1900+year if year>20
	replace year_temp	= 2000+year if year<20
	drop year
	ren year_temp year
}	

//SWITCH THIS ON TO GET THE ERROR MESSAGE WE GOT BACK FROM INEGI
/* 	gen year_temp=1900+year if year>20
	replace year_temp=2000+year if year<20
	drop year
	ren year_temp year	 */
	
//sales, output, labor, capital, inventories
	if `INEGI'==0 {
		local pre = "v"
	}
	if `INEGI'==1 {
		local pre = "V"
	}

	//sales and output
		ren `pre'25 output_value
		ren `pre'26 sales_h
		ren `pre'27 sales_x
		ren `pre'28 sales

	//labor
		ren `pre'1 workers
		ren `pre'3 hours
		ren `pre'4 wagebill
		
	//intermediate inputs
		ren `pre'6 inputs_domestic
		ren `pre'7 inputs_imported
		
	//inventories
		ren `pre'35 inventory_goods_start 	//inventories of finished goods at beginning of the year
		ren `pre'36 inventory_goods_end		//inventories of finished goods at end of the year
		ren `pre'39 inventory_inputs_start 	//inventories of inputs at beginning of the year
		ren `pre'40 inventory_inputs_end 	//inventories of inputs at end of the year
		
	//fixed assets
		//stock of capital at the start of the year
			ren `pre'43 k_machines_start 
			ren `pre'44 k_structures_start
			ren `pre'45 k_land_start
			ren `pre'46 k_transport_start
			ren `pre'47 k_other_start
	
		//increases of capital stock
			//purchases, domestic new
				ren `pre'49 k_machines_buy_h_new
				ren `pre'50 k_structures_buy_h_new
				ren `pre'51 k_transport_buy_h_new
				ren `pre'52 k_other_buy_h_new
				
			//purchases, domestic used
				ren `pre'54 k_machines_buy_h_used
				ren `pre'55 k_structures_buy_h_used
				ren `pre'56 k_land_buy_h_used
				ren `pre'57 k_transport_buy_h_used
				ren `pre'58 k_other_buy_h_used
				
			//purchases, imported
				ren `pre'60 k_machines_buy_m
				ren `pre'61 k_structures_buy_m
				ren `pre'62 k_transport_buy_m
				ren `pre'63 k_other_buy_m		
				
			//produced internally
				ren `pre'65 k_machines_produced
				ren `pre'66 k_structures_produced
				ren `pre'67 k_transport_produced
				ren `pre'68 k_other_produced	
				
		//reductions of capital stock		
			//depreciation 1
			//"bajas del ejercicio a valor historico"
			//approximate translation: depreciation over this period, at historic value
				ren `pre'70 k_machines_dep1
				ren `pre'71 k_structures_dep1
				ren `pre'72 k_land_dep1
				ren `pre'73 k_transport_dep1
				ren `pre'74 k_other_dep1
				
			//depreciation 2
			//"depreciacion historica asignada del ejercicio"
				ren `pre'82 k_machines_dep2
				ren `pre'83 k_structures_dep2
				ren `pre'84 k_transport_dep2
				ren `pre'85 k_other_dep2
			
			//sales
				ren `pre'92 k_machines_sales
				ren `pre'93 k_structures_sales
				ren `pre'94 k_land_sales
				ren `pre'95 k_transport_sales
				ren `pre'96 k_other_sales	

//all
	gen all = 1

//========================================================================//
//producer price index adjustment
//data from IFS (IMF)

if `ppi_adjustment'==1 {

	replace output_value = output_value/1.39 if year==1995
	replace output_value = output_value/1.86 if year==1996
	replace output_value = output_value/2.18 if year==1997
	replace output_value = output_value/2.53 if year==1998
	replace output_value = output_value/2.89 if year==1999

	replace sales_h = sales_h/1.39 if year==1995
	replace sales_h = sales_h/1.86 if year==1996
	replace sales_h = sales_h/2.18 if year==1997
	replace sales_h = sales_h/2.53 if year==1998
	replace sales_h = sales_h/2.89 if year==1999
	
	replace sales_x = sales_x/1.39 if year==1995
	replace sales_x = sales_x/1.86 if year==1996
	replace sales_x = sales_x/2.18 if year==1997
	replace sales_x = sales_x/2.53 if year==1998
	replace sales_x = sales_x/2.89 if year==1999

	replace sales = sales/1.39 if year==1995
	replace sales = sales/1.86 if year==1996
	replace sales = sales/2.18 if year==1997
	replace sales = sales/2.53 if year==1998
	replace sales = sales/2.89 if year==1999	
	
	replace wagebill = wagebill/1.39 if year==1995
	replace wagebill = wagebill/1.86 if year==1996
	replace wagebill = wagebill/2.18 if year==1997
	replace wagebill = wagebill/2.53 if year==1998
	replace wagebill = wagebill/2.89 if year==1999		
	
	replace inputs_domestic = inputs_domestic/1.39 if year==1995
	replace inputs_domestic = inputs_domestic/1.86 if year==1996
	replace inputs_domestic = inputs_domestic/2.18 if year==1997
	replace inputs_domestic = inputs_domestic/2.53 if year==1998
	replace inputs_domestic = inputs_domestic/2.89 if year==1999	

	replace inputs_imported = inputs_imported/1.39 if year==1995
	replace inputs_imported = inputs_imported/1.86 if year==1996
	replace inputs_imported = inputs_imported/2.18 if year==1997
	replace inputs_imported = inputs_imported/2.53 if year==1998
	replace inputs_imported = inputs_imported/2.89 if year==1999		
	
	replace inventory_goods_start = inventory_goods_start/1.39 if year==1995
	replace inventory_goods_start = inventory_goods_start/1.86 if year==1996
	replace inventory_goods_start = inventory_goods_start/2.18 if year==1997
	replace inventory_goods_start = inventory_goods_start/2.53 if year==1998
	replace inventory_goods_start = inventory_goods_start/2.89 if year==1999	
	
	replace inventory_goods_end = inventory_goods_end/1.39 if year==1995
	replace inventory_goods_end = inventory_goods_end/1.86 if year==1996
	replace inventory_goods_end = inventory_goods_end/2.18 if year==1997
	replace inventory_goods_end = inventory_goods_end/2.53 if year==1998
	replace inventory_goods_end = inventory_goods_end/2.89 if year==1999	
	
	replace inventory_inputs_start = inventory_inputs_start/1.39 if year==1995
	replace inventory_inputs_start = inventory_inputs_start/1.86 if year==1996
	replace inventory_inputs_start = inventory_inputs_start/2.18 if year==1997
	replace inventory_inputs_start = inventory_inputs_start/2.53 if year==1998
	replace inventory_inputs_start = inventory_inputs_start/2.89 if year==1999	
	
	replace inventory_inputs_end = inventory_inputs_end/1.39 if year==1995
	replace inventory_inputs_end = inventory_inputs_end/1.86 if year==1996
	replace inventory_inputs_end = inventory_inputs_end/2.18 if year==1997
	replace inventory_inputs_end = inventory_inputs_end/2.53 if year==1998
	replace inventory_inputs_end = inventory_inputs_end/2.89 if year==1999	
	
}
	
//========================================================================//
//setup main variables	

//panel data
	xtset plant_id year, yearly					
	
//new variables	
	//exporter dummy	
		gen exporter = 0
		replace exporter = 1 if sales_x > 0					
		
	//export intensity
		gen x_share = 0
		replace x_share = sales_x/sales if exporter==1			

		gen x_share2 = .
		replace x_share2 = sales_x/sales if exporter==1			
		
	//industry 2-digits
		tostring industry_id, g(industry_id_str)
		gen industry2d_id = substr(industry_id_str,1,2)
		destring industry2d_id, replace	
		drop industry_id_str
			
	//industry 3-digits
		tostring industry_id, g(industry_id_str)
		gen industry3d_id = substr(industry_id_str,1,3)
		destring industry3d_id, replace		
		drop industry_id_str		
		
	//final steady-state
		gen finalSS = 0
		replace finalSS = 1 if year>=`finalSS_start' & year<=`finalSS_end'
		
	//capital stock
		gen k = k_machines_start + k_structures_start + k_land_start + k_transport_start + k_other_start
		
	//investment
		//approach #1: back out investment given observed capital stock and depreciation rate
			local delta = 0.06
			gen i_gross = F.k - (1-`delta')*k if F.year==year+1 & F.k~=0 & F.k~=.
			gen i_net 	= F.k - k if F.year==year+1 & F.k~=0 & F.k~=.
		
		//approach #2: directly add up variables
			gen i_gross_increase 	= (k_machines_buy_h_new + k_structures_buy_h_new + k_transport_buy_h_new + k_other_buy_h_new) + (k_machines_buy_h_used + k_structures_buy_h_used + k_land_buy_h_used + k_transport_buy_h_used + k_other_buy_h_used) + (k_machines_buy_m + k_structures_buy_m + k_transport_buy_m + k_other_buy_m) + (k_machines_produced + k_structures_produced + k_transport_produced + k_other_produced)
			gen i_gross_decrease 	= (k_machines_dep1 + k_structures_dep1 + k_land_dep1 + k_transport_dep1 + k_other_dep1) + (k_machines_sales + k_structures_sales + k_land_sales + k_transport_sales + k_other_sales)
			gen i_net_decrease 		= (k_machines_sales + k_structures_sales + k_land_sales + k_transport_sales + k_other_sales)			
			
			gen i_gross_direct 		= i_gross_increase - i_gross_decrease
			gen i_net_direct 		= i_gross_increase - i_net_decrease
		
		//ppi-adjustment
			if `ppi_adjustment'==1 {
				replace k = k/1.39 if year==1995
				replace k = k/1.86 if year==1996
				replace k = k/2.18 if year==1997
				replace k = k/2.53 if year==1998
				replace k = k/2.89 if year==1999			
			
				replace i_gross = i_gross/1.39 if year==1995
				replace i_gross = i_gross/1.86 if year==1996
				replace i_gross = i_gross/2.18 if year==1997
				replace i_gross = i_gross/2.53 if year==1998
				replace i_gross = i_gross/2.89 if year==1999
				
				replace i_net = i_net/1.39 if year==1995
				replace i_net = i_net/1.86 if year==1996
				replace i_net = i_net/2.18 if year==1997
				replace i_net = i_net/2.53 if year==1998
				replace i_net = i_net/2.89 if year==1999				
				
				replace i_gross_increase = i_gross_increase/1.39 if year==1995
				replace i_gross_increase = i_gross_increase/1.86 if year==1996
				replace i_gross_increase = i_gross_increase/2.18 if year==1997
				replace i_gross_increase = i_gross_increase/2.53 if year==1998
				replace i_gross_increase = i_gross_increase/2.89 if year==1999

				replace i_gross_decrease = i_gross_decrease/1.39 if year==1995
				replace i_gross_decrease = i_gross_decrease/1.86 if year==1996
				replace i_gross_decrease = i_gross_decrease/2.18 if year==1997
				replace i_gross_decrease = i_gross_decrease/2.53 if year==1998
				replace i_gross_decrease = i_gross_decrease/2.89 if year==1999

				replace i_net_decrease = i_net_decrease/1.39 if year==1995
				replace i_net_decrease = i_net_decrease/1.86 if year==1996
				replace i_net_decrease = i_net_decrease/2.18 if year==1997
				replace i_net_decrease = i_net_decrease/2.53 if year==1998
				replace i_net_decrease = i_net_decrease/2.89 if year==1999

				replace i_gross_direct = i_gross_direct/1.39 if year==1995
				replace i_gross_direct = i_gross_direct/1.86 if year==1996
				replace i_gross_direct = i_gross_direct/2.18 if year==1997
				replace i_gross_direct = i_gross_direct/2.53 if year==1998
				replace i_gross_direct = i_gross_direct/2.89 if year==1999				
				
				replace i_net_direct = i_net_direct/1.39 if year==1995
				replace i_net_direct = i_net_direct/1.86 if year==1996
				replace i_net_direct = i_net_direct/2.18 if year==1997
				replace i_net_direct = i_net_direct/2.53 if year==1998
				replace i_net_direct = i_net_direct/2.89 if year==1999					
			}		

	//inventories (at the start of the period)
		gen inventories 		= inventory_goods_start
		gen inventories_sales 	= inventories/sales
		
		gen inventories_inputs 				= inventory_inputs_start
		gen inventoriesinputs_intermediates = inventory_inputs_start/(inputs_domestic+inputs_imported)
		
	//share of imported intermediates
		gen intermediates_mshare = inputs_imported/(inputs_domestic+inputs_imported)
		
		gen importedintermediates_wagebill = inputs_imported/wagebill
		
		gen importedintermediates_sales = inputs_imported/sales
		
	//tag plants
		egen tag_plant = tag(plant_id)
		
//create missing values
	replace output_value	= . if output_value	<=0
	replace sales			= . if sales		<=0	
	
	replace sales_h			= . if sales_h	<0
	replace sales_x			= . if sales_x	<0
	replace sales_x			= . if sales_x==0 & sales_h==0 
	replace sales_h			= . if sales_x==0 & sales_h==0 
			
	replace workers			= . if workers	<=0	
	replace hours			= . if hours	<=0	
	replace wagebill		= . if wagebill	<=0
	
	replace k				= . if k			<=0
	replace i_net			= . if i_net		<-k
	replace i_net_direct	= . if i_net_direct	<-k
	
	replace inventories						= . if inventories						<0
	replace inventories_sales				= . if inventories_sales				<0		
	replace inventories_inputs 				= . if inventories_inputs				<0
	replace inventoriesinputs_intermediates	= . if inventoriesinputs_intermediates	<0
	replace intermediates_mshare			= . if intermediates_mshare				<0
	replace importedintermediates_wagebill	= . if importedintermediates_wagebill	<0
	
//log
	gen ln_wagebill = ln(wagebill)
	
//clean data
	if `year_end'~=2003 {
		drop if year>`year_end'
	}
		
	//drop firms with missing values of key variables
		drop if sales_h	== . 
		drop if sales_x	== .
		drop if sales	== . 

	//drop firms with missing observations within observed period
	//missing observations = a (firm,year) pair such that (firm,year-k) and (firm,year+t) exist for some t,k
		bysort plant_id: egen year_max = max(year)
		bysort plant_id: egen year_min = min(year)
		bysort plant_id: egen obs_num  = count(year)
		gen year_num = year_max-year_min+1
		
		gen missing = 0
		replace missing = 1 if year_num~=obs_num
		drop if missing == 1

	//balanced panel
		if `balanced'==1 {
			keep if year_min==`year_start' & year_max==`year_end'
		}

//========================================================================//
//construct additional variables
		
//export intensity categories 
//classify firms based on export intensity in year 1994			
	gen xy_cat_temp 	= 0 if x_share==0                & x_share~=. & year==1994 	//Non-exporters
	replace xy_cat_temp = 1 if x_share<=0.6 & x_share~=0 & x_share~=. & year==1994	//X/Y<=0.6
	replace xy_cat_temp = 2 if x_share>0.6               & x_share~=. & year==1994	//X/Y>0.6
	bysort plant_id: egen xy_cat = mean(xy_cat_temp)	
	
	gen x_share_1994_temp = x_share if year==1994
	bysort plant_id: egen x_share_1994 = mean(x_share_1994_temp)
	
//size categories
//classify firms based on sales in year 1994	
	egen sales_p10temp = pctile(sales) if year==1994, p(10)
	egen sales_p25temp = pctile(sales) if year==1994, p(25)
	egen sales_p50temp = pctile(sales) if year==1994, p(50)	
	egen sales_p75temp = pctile(sales) if year==1994, p(75)	
	egen sales_p90temp = pctile(sales) if year==1994, p(90)
	
	bysort plant_id: egen sales_p10 = mean(sales_p10temp)
	bysort plant_id: egen sales_p25 = mean(sales_p25temp)
	bysort plant_id: egen sales_p50 = mean(sales_p50temp)
	bysort plant_id: egen sales_p75 = mean(sales_p75temp)
	bysort plant_id: egen sales_p90 = mean(sales_p90temp)
	
	gen size_cat_temp 	  = .
	replace size_cat_temp = 1 if year==1994 & sales<=sales_p25
	replace size_cat_temp = 2 if year==1994 & sales>sales_p25 & sales<=sales_p50
	replace size_cat_temp = 3 if year==1994 & sales>sales_p50 & sales<=sales_p75
	replace size_cat_temp = 4 if year==1994 & sales>sales_p75 & sales~=.
	bysort plant_id: egen size_cat = mean(size_cat_temp)		
		
//growth rates relative to 1994	
	//exports
		gen ln_sales_x_1994_temp 				= ln(sales_x) if year==1994
		bysort plant_id: egen ln_sales_x_1994 	= mean(ln_sales_x_1994_temp)
				
		gen ln_sales_x 			= ln(sales_x)
		gen dln_sales_x_rel1994 = ln_sales_x - ln_sales_x_1994
		
	//domestic sales
		gen ln_sales_h_1994_temp 				= ln(sales_h) if year==1994
		bysort plant_id: egen ln_sales_h_1994 	= mean(ln_sales_h_1994_temp)
		
		gen ln_sales_h 			= ln(sales_h)
		gen dln_sales_h_rel1994 = ln_sales_h - ln_sales_h_1994	
		
	//sales
		gen ln_sales_1994_temp 				= ln(sales) if year==1994
		bysort plant_id: egen ln_sales_1994 = mean(ln_sales_1994_temp)
		
		gen ln_sales 		  = ln(sales)
		gen dln_sales_rel1994 = ln_sales - ln_sales_1994	

	//capital
		gen ln_k_1994_temp 				= ln(k) if year==1994
		bysort plant_id: egen ln_k_1994 = mean(ln_k_1994_temp)
		
		gen ln_k 		  = ln(k)
		gen dln_k_rel1994 = ln_k - ln_k_1994	
	
//years after devaluation
	gen years_after_devaluation = year - 1994	
	
//========================================================================//
//summary of dataset and variables

//dataset
	describe, f
	xtdescribe, p(9)
	tab obs_num if tag_plant==1 //number of plants by length of observed spells

//industries
	tab industry_id
	tab industry3d_id
	tab industry2d_id
	
//summary of raw variables
	sum output_value 					if year==1994, detail
	sum sales 							if year==1994, detail
	sum sales_x 						if year==1994, detail
	sum sales_h 						if year==1994, detail
	sum workers 						if year==1994, detail
	sum hours 							if year==1994, detail
	sum wagebill 						if year==1994, detail
	sum k 								if year==1994, detail
	sum i_net 							if year==1994, detail
	sum i_gross 						if year==1994, detail
	sum i_net_direct 					if year==1994, detail
	sum i_gross_direct 					if year==1994, detail
	sum inventories 					if year==1994, detail
	sum inventories_sales 				if year==1994, detail
	sum inventories_inputs 				if year==1994, detail
	sum inventoriesinputs_intermediates if year==1994, detail
	sum intermediates_mshare 			if year==1994, detail
	sum importedintermediates_wagebill 	if year==1994, detail
	
//summary of constructed variables
	sum exporter 	if year==1994, detail
	sum x_share 	if year==1994, detail
	sum x_share2 	if year==1994, detail
	
//========================================================================//
//PAPER
//========================================================================//

//calibration targets, period 1994
//used to compute table 1 of the paper

//share of exporters
	sum exporter if year==1994, detail
		
//share of firms with high export intensity		
	tab xy_cat if year==1994
	
//average export intensity by export intensity category
	sum x_share2 if xy_cat==1 & year==1994, detail
	sum x_share2 if xy_cat==2 & year==1994, detail
	
//share of total sales by largest 25% of firms
	table size_cat if year==1994, c(sum sales count sales)
	
//sd log sales	
	sum ln_sales if year==1994, detail	
	
//========================================================================//
//validation evidence: average firm-level adjustment by export intensity
//used to compute figure 5

//construct dummies
	gen year1995_xylow = 0
	replace year1995_xylow = 1 if year==1995 & xy_cat==1

	gen year1996_xylow = 0
	replace year1996_xylow = 1 if year==1996 & xy_cat==1

	gen year1997_xylow = 0
	replace year1997_xylow = 1 if year==1997 & xy_cat==1

	gen year1998_xylow = 0
	replace year1998_xylow = 1 if year==1998 & xy_cat==1

	gen year1999_xylow = 0
	replace year1999_xylow = 1 if year==1999 & xy_cat==1

	gen year1995_xyhigh = 0
	replace year1995_xyhigh = 1 if year==1995 & xy_cat==2

	gen year1996_xyhigh = 0
	replace year1996_xyhigh = 1 if year==1996 & xy_cat==2

	gen year1997_xyhigh = 0
	replace year1997_xyhigh = 1 if year==1997 & xy_cat==2

	gen year1998_xyhigh = 0
	replace year1998_xyhigh = 1 if year==1998 & xy_cat==2

	gen year1999_xyhigh = 0
	replace year1999_xyhigh = 1 if year==1999 & xy_cat==2

//regressions
	reg dln_sales_x_rel1994 year1996_xylow year1997_xylow year1998_xylow year1999_xylow year1995_xyhigh year1996_xyhigh year1997_xyhigh year1998_xyhigh year1999_xyhigh if year>1994, robust

	reg dln_sales_x_rel1994 year1996_xylow year1997_xylow year1998_xylow year1999_xylow year1995_xyhigh year1996_xyhigh year1997_xyhigh year1998_xyhigh year1999_xyhigh inventories_sales inventoriesinputs_intermediates importedintermediates_wagebill if year>1994, robust

	reg dln_sales_x_rel1994 year1996_xylow year1997_xylow year1998_xylow year1999_xylow year1995_xyhigh year1996_xyhigh year1997_xyhigh year1998_xyhigh year1999_xyhigh i.industry2d_id inventories_sales inventoriesinputs_intermediates importedintermediates_wagebill if year>1994, robust

	reg dln_sales_x_rel1994 year1996_xylow year1997_xylow year1998_xylow year1999_xylow year1995_xyhigh year1996_xyhigh year1997_xyhigh year1998_xyhigh year1999_xyhigh i.years_after_devaluation#i.industry2d_id inventories_sales inventoriesinputs_intermediates importedintermediates_wagebill if year>1994, robust		
	
//========================================================================//
//validation evidence: decomposition relative to pre-devaluation
//used to compute table 3 

//Exporter in period t, but non-exporter in 1994
	gen decomp_newexporter2 	= 0 //Discard sums corresponding to these values
	replace decomp_newexporter2 = 1 if exporter==1 & L.exporter==0 & year==1995
	replace decomp_newexporter2 = 1 if exporter==1 & L2.exporter==0 & year==1996
	replace decomp_newexporter2 = 1 if exporter==1 & L3.exporter==0 & year==1997
	replace decomp_newexporter2 = 1 if exporter==1 & L4.exporter==0 & year==1998
	replace decomp_newexporter2 = 1 if exporter==1 & L5.exporter==0 & year==1999

//Exporter in 1994, but non-exporter in period t
	gen decomp_futurenonx2_1994_1995 		= 0 //Discard sums corresponding to these values
	replace decomp_futurenonx2_1994_1995 	= 1 if exporter==1 & F.exporter==0 & year==1994

	gen decomp_futurenonx2_1994_1996 		= 0 //Discard sums corresponding to these values
	replace decomp_futurenonx2_1994_1996 	= 1 if exporter==1 & F2.exporter==0 & year==1994

	gen decomp_futurenonx2_1994_1997 		= 0 //Discard sums corresponding to these values
	replace decomp_futurenonx2_1994_1997 	= 1 if exporter==1 & F3.exporter==0 & year==1994

	gen decomp_futurenonx2_1994_1998 		= 0 //Discard sums corresponding to these values
	replace decomp_futurenonx2_1994_1998 	= 1 if exporter==1 & F4.exporter==0 & year==1994

	gen decomp_futurenonx2_1994_1999 		= 0 //Discard sums corresponding to these values
	replace decomp_futurenonx2_1994_1999 	= 1 if exporter==1 & F5.exporter==0 & year==1994

//Exporter in 1994 and in period t, values after the devaluation
	gen decomp_contexporter_back2 		= 0 //Discard sums corresponding to these values
	replace decomp_contexporter_back2 	= 1 if exporter==1 & L.exporter==1 & year==1995
	replace decomp_contexporter_back2 	= 1 if exporter==1 & L2.exporter==1 & year==1996
	replace decomp_contexporter_back2 	= 1 if exporter==1 & L3.exporter==1 & year==1997	
	replace decomp_contexporter_back2 	= 1 if exporter==1 & L4.exporter==1 & year==1998	
	replace decomp_contexporter_back2 	= 1 if exporter==1 & L5.exporter==1 & year==1999
	
//Exporter in 1994 and in period t, values in 1994
	gen decomp_contx_fwd2_1994_1995 	= 0 //Discard sums corresponding to these values
	replace decomp_contx_fwd2_1994_1995 = 1 if exporter==1 & F.exporter==1 & year==1994

	gen decomp_contx_fwd2_1994_1996 	= 0 //Discard sums corresponding to these values
	replace decomp_contx_fwd2_1994_1996 = 1 if exporter==1 & F2.exporter==1 & year==1994
	
	gen decomp_contx_fwd2_1994_1997 	= 0 //Discard sums corresponding to these values
	replace decomp_contx_fwd2_1994_1997 = 1 if exporter==1 & F3.exporter==1 & year==1994

	gen decomp_contx_fwd2_1994_1998 	= 0 //Discard sums corresponding to these values
	replace decomp_contx_fwd2_1994_1998 = 1 if exporter==1 & F4.exporter==1 & year==1994

	gen decomp_contx_fwd2_1994_1999 	= 0 //Discard sums corresponding to these values
	replace decomp_contx_fwd2_1994_1999 = 1 if exporter==1 & F5.exporter==1 & year==1994	

//All exports
	table year, c(sum sales_x) //Total exports by year

	table year decomp_newexporter2, c(sum sales_x) f(%20.8g) //Total exports by firms that export in current period but did not export in the previous period 
	
	table year decomp_futurenonx2_1994_1995, c(sum sales_x) f(%20.8g) //Total exports by firms that export in current period but do not export in the following period 
	table year decomp_futurenonx2_1994_1996, c(sum sales_x) f(%20.8g) //Total exports by firms that export in current period but do not export in the following period 
	table year decomp_futurenonx2_1994_1997, c(sum sales_x) f(%20.8g) //Total exports by firms that export in current period but do not export in the following period 
	table year decomp_futurenonx2_1994_1998, c(sum sales_x) f(%20.8g) //Total exports by firms that export in current period but do not export in the following period 
	table year decomp_futurenonx2_1994_1999, c(sum sales_x) f(%20.8g) //Total exports by firms that export in current period but do not export in the following period 
	
	table year decomp_contexporter_back2	, c(sum sales_x) f(%20.8g) //Total sales by firms that export in current period and in the previous period
	table year decomp_contx_fwd2_1994_1995	, c(sum sales_x) f(%20.8g) //Total sales by firms that export in current period and in the following period	
	table year decomp_contx_fwd2_1994_1996	, c(sum sales_x) f(%20.8g) //Total sales by firms that export in current period and in the following period	
	table year decomp_contx_fwd2_1994_1997	, c(sum sales_x) f(%20.8g) //Total sales by firms that export in current period and in the following period	
	table year decomp_contx_fwd2_1994_1998	, c(sum sales_x) f(%20.8g) //Total sales by firms that export in current period and in the following period	
	table year decomp_contx_fwd2_1994_1999	, c(sum sales_x) f(%20.8g) //Total sales by firms that export in current period and in the following period	
	
//========================================================================//
//APPENDIX
//========================================================================//

//export intensity distribution
//used to compute fig 2 and table 3 

//bins
	gen xy_category 	= 0
	replace xy_category = 1 	if x_share2>0 & x_share2<=0.1
	replace xy_category = 2 	if x_share2>0.1 & x_share2<=0.2
	replace xy_category = 3 	if x_share2>0.2 & x_share2<=0.3
	replace xy_category = 4 	if x_share2>0.3 & x_share2<=0.4
	replace xy_category = 5 	if x_share2>0.4 & x_share2<=0.5
	replace xy_category = 6 	if x_share2>0.5 & x_share2<=0.6
	replace xy_category = 7 	if x_share2>0.6 & x_share2<=0.7
	replace xy_category = 8 	if x_share2>0.7 & x_share2<=0.8
	replace xy_category = 9 	if x_share2>0.8 & x_share2<=0.9	
	replace xy_category = 10 	if x_share2>0.9 & x_share2<=1

//histograms	
//used to compute fig 2 of the paper
	tab xy_category if year==1994
	tab xy_category if year==1999
	tab xy_category	
		
//statistics by export intensity bin
//used to compute table 3 
	table xy_category 	if year==1994, c(sum sales sum sales_x sum sales_h mean x_share count exporter)	
	table xy_cat 		if year==1994, c(sum sales sum sales_x sum sales_h mean x_share count exporter)	
	
//========================================================================//
//validation evidence: exports growth and external finance dependence
	
table industry3d_id year, c(count plant_id sum sales sum sales_x sum sales_h)
table industry3d_id year, c(mean exporter mean x_share sum exporter)
table industry3d_id year, c(sum workers sum hours sum wagebill sum k)
table industry3d_id year, c(sum i_net sum i_gross sum i_net_direct sum i_gross_direct)		
		
//========================================================================//
//Imported intermediates

//Summary statistics on use of intermediate inputs
// 	=> All firms
// 	=> Exporters vs. non-exporters
// 	=> High X/Y vs. low X/Y

	sum intermediates_mshare if year==1994, detail
	sum intermediates_mshare if year==1994 & exporter==0, detail
	sum intermediates_mshare if year==1994 & exporter==1, detail
	sum intermediates_mshare if year==1994 & xy_cat==0, detail
	sum intermediates_mshare if year==1994 & xy_cat==1, detail
	sum intermediates_mshare if year==1994 & xy_cat==2, detail	
	
	sum importedintermediates_wagebill if year==1994, detail
	sum importedintermediates_wagebill if year==1994 & exporter==0, detail
	sum importedintermediates_wagebill if year==1994 & exporter==1, detail
	sum importedintermediates_wagebill if year==1994 & xy_cat==0, detail
	sum importedintermediates_wagebill if year==1994 & xy_cat==1, detail
	sum importedintermediates_wagebill if year==1994 & xy_cat==2, detail		
	
	sum importedintermediates_sales if year==1994, detail
	sum importedintermediates_sales if year==1994 & exporter==0, detail
	sum importedintermediates_sales if year==1994 & exporter==1, detail
	sum importedintermediates_sales if year==1994 & xy_cat==0, detail
	sum importedintermediates_sales if year==1994 & xy_cat==1, detail
	sum importedintermediates_sales if year==1994 & xy_cat==2, detail		

//Extend specification from the paper by controlling also for import intensity category x year dummies
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.intermediates_mshare_high#i.years_after_devaluation if year>1994, robust
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.industry2d_id inventories_sales inventoriesinputs_intermediates i.intermediates_mshare_high#i.years_after_devaluation if year>1994, robust
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.years_after_devaluation#i.industry2d_id inventories_sales inventoriesinputs_intermediates i.intermediates_mshare_high#i.years_after_devaluation if year>1994, robust

	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.importedint_wagebill_high#i.years_after_devaluation if year>1994, robust
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.industry2d_id inventories_sales inventoriesinputs_intermediates i.importedint_wagebill_high#i.years_after_devaluation if year>1994, robust
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.years_after_devaluation#i.industry2d_id inventories_sales inventoriesinputs_intermediates i.importedint_wagebill_high#i.years_after_devaluation if year>1994, robust	

	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.importedint_sales_high#i.years_after_devaluation if year>1994, robust
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.industry2d_id inventories_sales inventoriesinputs_intermediates i.importedint_sales_high#i.years_after_devaluation if year>1994, robust
	reg dln_sales_x_rel1994 i.xy_cat#i.years_after_devaluation i.years_after_devaluation#i.industry2d_id inventories_sales inventoriesinputs_intermediates i.importedint_sales_high#i.years_after_devaluation if year>1994, robust		

//========================================================================//
//end

if `logfile'==1 {
	log close			
}

//========================================================================//