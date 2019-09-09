//========================================================================================================//

clear
set more off
pause on

//local directory_input = "C:\Users\David\Dropbox\Compartida\Compartida con Fer y Mike\Large devaluations\SubmissionApril2017\Data\March2017\KLS2 Section 2 Data"
local directory_input = "C:\Dropbox\Shared\Large devaluations\SubmissionApril2017\Data\March2017\KLS2 Section 2 Data"
//local directory_input = "C:\Users\Michal\Dropbox\Large devaluations\Data\WoldBank Enterprise Survey"
use "`directory_input'\comprehensive.dta", clear

//========================================================================================================//
//Setup variables

/*
The weights in the more recent Enterprise Surveys data are probability weights. Using these weights allows inference on the population of non-agricultural private firms (that meet the Enterprise Surveys eligibility criteria) in a country. In Stata, a survey design should be declared before performing any analysis. Specifically, this command should be used: svyset idstd [pweight=wt], strata(strata). The survey command svy should also be used in calculating any statistics to be interpreted for the population of non-agricultural private firms. For statistics related to specific types of firms, analysts should use the subpopulation option in Stata.

http://www.ats.ucla.edu/stat/stata/faq/svy_stata_subpop.htm
*/

//Devaluation countries
	gen devaluation = 0
	replace devaluation = 1 if country=="Argentina2006" || country=="Brazil2003" || country=="Indonesia2003" || country=="Malaysia2002"  || country=="Mexico2006" || country=="Turkey2002"|| country=="Turkey2004"|| country=="Turkey2005"	
		
	gen devaluation_country = ""
	replace devaluation_country = country if devaluation==1
	
//IDs		  
    ren idstd id_firm 
	encode country, g(survey_id)
	
//Weight
	gen weight = wt
	replace weight = 1 if wt==.
		
//Sales
	gen ln_sales = ln(c274a1y) //Total sales one year ago (000s LCU)
	
//Export intensity
	gen x_share = c211a2/100 //% of sales exporter directly

//Exporter
	gen exporter_old = exporter
	drop exporter
	
	gen exporter_wb = .
	replace exporter_wb = 0 if exporter_old==2
	replace exporter_wb = 1 if exporter_old==1
	
	gen exporter = 0 if x_share==0
	replace exporter = 1 if x_share>0 & x_share<=1
		
//Foreign debt
    gen fdebt_share = c230/100 //Share of total borrowing denominated in foreign currency
	
	gen fdebt_dummy = .
	replace fdebt_dummy = 0 if fdebt_share==0
	replace fdebt_dummy = 1 if fdebt_share>0 & fdebt_share<=1
	
	gen fdebt_share_nozeros = .
	replace fdebt_share_nozeros = fdebt_share if fdebt_dummy==1
	
//Sample selection	
	keep if sector==1 //Restrict attention to manufactures
	
//========================================================================================================//
//Foreign-denominated debt & financial constraints

//drop if income==1 //Excluding low income

//Foreign-denominated debt
	table country if devaluation==1 [pw=weight], c(mean fdebt_dummy)
	table country exporter if devaluation==1 [pw=weight], c(mean fdebt_dummy)

	table country if devaluation==1 [pw=weight], c(mean fdebt_share_nozeros)
	table country exporter if devaluation==1 [pw=weight], c(mean fdebt_share_nozeros)

	//Counts
		table country if devaluation==1, c(count fdebt_dummy)
		table country exporter if devaluation==1, c(count fdebt_dummy)

		table country if devaluation==1, c(count fdebt_share_nozeros)
		table country exporter if devaluation==1, c(count fdebt_share_nozeros)
	
//Financial constraints
	gen accesstofinance_mobstacle = .
	replace accesstofinance_mobstacle = 0 if c218k<=1 & c218k~=.
	replace accesstofinance_mobstacle = 1 if c218k>1 & c218k~=.

	gen costfinance_mobstacle = .
	replace costfinance_mobstacle = 0 if c218l<=1 & c218l~=.
	replace costfinance_mobstacle = 1 if c218l>1 & c218l~=.	
	
	table country if devaluation==1 [pw=weight], c(mean accesstofinance_mobstacle)
	table country exporter if devaluation==1 [pw=weight], c(mean accesstofinance_mobstacle)

	table country if devaluation==1 [pw=weight], c(mean costfinance_mobstacle)
	table country exporter if devaluation==1 [pw=weight], c(mean costfinance_mobstacle)

	//Counts
		table country if devaluation==1, c(count accesstofinance_mobstacle)
		table country exporter if devaluation==1, c(count accesstofinance_mobstacle)

		table country if devaluation==1, c(count costfinance_mobstacle)
		table country exporter if devaluation==1, c(count costfinance_mobstacle)	
	
//Interaction between foreign-denominate debt and credit constraints
	table country fdebt_dummy if devaluation==1 [pw=weight], c(count accesstofinance_mobstacle)
	table country fdebt_dummy if devaluation==1 [pw=weight], c(count costfinance_mobstacle)
	
	//Counts
		table country fdebt_dummy if devaluation==1, c(count accesstofinance_mobstacle)
		table country fdebt_dummy if devaluation==1, c(count costfinance_mobstacle)	
		
//========================================================================================================//
//Foreign-denominated debt & financial constraints
//Statistics by size

gen workers = c262a1y

gen workers_cat = .
replace workers_cat = 1 if workers<=25 & workers~=.
replace workers_cat = 2 if workers>25 & workers<=100 & workers~=.
replace workers_cat = 3 if workers>100 & workers<=250 & workers~=.
replace workers_cat = 4 if workers>250 & workers~=.

table country workers_cat if devaluation==1 [pw=weight], c(mean fdebt_dummy)
table country workers_cat if devaluation==1 [pw=weight], c(mean fdebt_share_nozeros)

table country workers_cat if devaluation==1 [pw=weight], c(mean accesstofinance_mobstacle)
table country workers_cat if devaluation==1 [pw=weight], c(mean costfinance_mobstacle)
	
//Counts	
	table country workers_cat if devaluation==1, c(count fdebt_dummy)
	table country workers_cat if devaluation==1, c(count fdebt_share_nozeros)

	table country workers_cat if devaluation==1, c(count accesstofinance_mobstacle)
	table country workers_cat if devaluation==1, c(count costfinance_mobstacle)

//========================================================================================================//
//Foreign-denominated debt & financial constraints
//Statistics by export intensity

gen x_share_cat = .
replace x_share_cat = 0 if x_share==0
replace x_share_cat = 1 if x_share>0 & x_share<=0.60
replace x_share_cat = 2 if x_share>0.60 & x_share<=1

table country x_share_cat if devaluation==1 [pw=weight], c(mean fdebt_dummy)
table country x_share_cat if devaluation==1 [pw=weight], c(mean fdebt_share_nozeros)

table country x_share_cat if devaluation==1 [pw=weight], c(mean accesstofinance_mobstacle)
table country x_share_cat if devaluation==1 [pw=weight], c(mean costfinance_mobstacle)

//Counts
	table country x_share_cat if devaluation==1, c(count fdebt_dummy)
	table country x_share_cat if devaluation==1, c(count fdebt_share_nozeros)

	table country x_share_cat if devaluation==1, c(count accesstofinance_mobstacle)
	table country x_share_cat if devaluation==1, c(count costfinance_mobstacle)

//========================================================================================================//

