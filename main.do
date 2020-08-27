clear 
set more off 

*Packages 
*ssc install sxpose   
*ssc install dmariano


*ADD main directory
global main_dir "C:\Users\Prerana Hiriyur\Desktop\EC4304 Project STATA"
 

global data "$main_dir\data"
global graphs "$main_dir\graphs"

cd "$data"


/********
USA DATA
*********/

*USA CPI Data
preserve 
	clear
	freduse CPIAUCNS
	rename CPIAUCNS usa_cpi_nsa
	label var usa_cpi_nsa "USA Monthly CPI, NSA"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "usa_cpi.dta", replace 
restore

*US Fed funds rate 
preserve 
	clear 
	freduse FEDFUNDS
	rename FEDFUNDS usa_ir
	label var usa_ir "Effective Federal Funds Rate"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "us_fed_funds.dta", replace
restore 

*merge usa datasets 
use usa_cpi.dta, clear
merge 1:1 date using us_fed_funds.dta, gen(merge_fed_funds)
save usa.dta, replace 

/**********
JAPAN DATA
**********/
*JAPAN/USA Ex rate
preserve 
	clear
	freduse EXJPUS
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	rename EXJPUS jpy_exrate_nsa 
	label var jpy_exrate_nsa "Japan / U.S. Foreign Exchange Rate Monthly, NSA"
	save "jpy_exrate.dta", replace
restore

*Japan CPI Data
preserve
	clear 
	freduse JPNCPIALLMINMEI
	rename JPNCPIALLMINMEI jpy_cpi_nsa
	label var jpy_cpi_nsa "Japan Monthly CPI, NSA"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "jpy_cpi.dta", replace 
restore

*Japan Interbank rates - Immediate
preserve
	clear 
	freduse IRSTCI01JPM156N
	rename IRSTCI01JPM156N jpy_ir
	label var jpy_ir "Immediate Rates: Less than 24 Hours: Call Money/Interbank Rate for Japan"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "jpy_interbank_imm.dta", replace
restore 



*merge japan datasets
use jpy_exrate.dta, clear 
merge 1:1 date using jpy_cpi.dta, gen(merge_jpy_cpi)
merge 1:1 date using jpy_interbank_imm.dta, gen(merge_jpy_intb_imm)
sort date 
order date jpy_exrate_nsa jpy_cpi_nsa jpy_ir 
save jpy.dta, replace 

 

/*************************************************************
Dataset for EUR USD - Germany as Proxy for Economic variables 
**************************************************************/
*EUR USD Ex rate data
preserve 
	clear
	freduse EXUSEU
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	rename EXUSEU eur_exrate_nsa 
	replace eur_exrate_nsa = 1/eur_exrate_nsa
	label var eur_exrate_nsa "Euro / U.S. Foreign Exchange Rate Monthly, NSA"
	save eur_exrate.dta, replace
restore 

/************
GERMANY DATA
************/

*Germany CPI Data
preserve
	clear 
	freduse DEUCPIALLMINMEI
	rename DEUCPIALLMINMEI ger_cpi_nsa
	label var ger_cpi_nsa "Germany Monthly CPI, NSA"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "ger_cpi.dta", replace 
restore

*Germany Interbank rates - Immediate 
preserve
	clear 
	freduse IRSTCI01DEM156N
	rename IRSTCI01DEM156N ger_ir
	label var ger_ir "Immediate Rates: Less than 24 Hours: Call Money/Interbank Rate for Germany"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "ger_interbank_imm.dta", replace
restore 


/**************
EURO AREA DATA
***************/

*Harmonized Euro Area CPI Index 
preserve 
	clear 
	freduse CP0000EZ19M086NEST
	rename CP0000EZ19M086NEST eur_cpi_nsa
	label var eur_cpi_nsa "Harmonized Euro Area Monthly CPI (19 Countries), NSA"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "eur_cpi.dta", replace 
restore

*Euro Area Interbank Rates - Immediate 
preserve
	clear 
	freduse IRSTCI01EZM156N 
	rename IRSTCI01EZM156N eur_ir
	label var eur_ir "Immediate Rates: Less than 24 Hours: Call Money/Interbank Rate for Euro Area"
	drop date 
	gen date = mofd(daten)
	format date %tm 
	drop daten
	save "eur_interbank_imm.dta", replace
restore 

*merge eur datasets 
use eur_exrate.dta, clear 
foreach v in eur ger{
preserve 
	merge 1:1 date using "`v'_cpi.dta", gen(merge_`v'_cpi)
	merge 1:1 date using "`v'_interbank_imm.dta", gen(merge_`v'_intb_imm)
	sort date 
	order date eur_exrate_nsa `v'_cpi_nsa `v'_ir 
	save eur_`v'.dta, replace 
restore
}

/******************************
Importing Production Indicator
*******************************/ 
preserve 
	clear 
	import excel ppl_nsa, clear
	drop if _n <= 6
	drop A
	rename B country
	keep if _n == 1 | country == "Japan" | country == "United States" | country == "Germany"
	replace country = "UnitedStates" if country == "United States" /*just for convenience cos stata vars cant have spaces*/
	sxpose, clear firstnames /*transposes the dataset*/
	drop if _n <= 2 /*strings*/
	rename Country date_str
	rename Japan jpy_prod 
	label var jpy_prod "Japan Production Index, NSA"
	rename Germany ger_prod 
	label var ger_prod "Germany Production Index, NSA"
	rename UnitedStates usa_prod
	label var usa_prod "USA Production Index, NSA"
	*dealing with year and quarter date values in the middle of months
	*logic: if the 5th place of date string M it is monthly
	gen m = substr(date_str, 5,1)
	*dropping if not a monthly value 
	drop if m != "M"
	drop m
	gen date = monthly(date_str, "YM")
	format date %tm
	drop date_str
	destring jpy_prod usa_prod ger_prod, replace force
	save prod_indx.dta, replace 
restore 


/***************************
VIX monthly from Yahoo.com
****************************/
preserve
	clear
	import delimited vix
	gen daten = daily(date, "YMD")
	format daten %d
	drop date
	keep daten close 
	rename close vix
	label var vix "VIX Closing Values Monthly"
	gen date = mofd(daten)
	format date %tm 
	drop daten 
	*yahoo has 2 values of most recent data, drop duplicate
	duplicates drop date, force
	save vix.dta, replace 
restore


/****************************
merging all datasets together
*****************************/ 
use jpy.dta, clear
merge 1:1 date using eur_eur.dta, gen(merge_eur)
merge 1:1 date eur_exrate_nsa using eur_ger.dta, gen(merge_ger) 
merge 1:1 date using usa.dta, gen(merge_usa)
merge 1:1 date using prod_indx.dta, gen(merge_prod)
merge 1:1 date using vix.dta, gen(merge_vix)

save final_dataset_full.dta, replace 

*Finally choose our sample 
keep if date >= tm(2001m1) /*data available for all variables*/
*keep if date < tm(2018m8) /*data not available for jpy_prod after this*/

drop merge*
save final_dataset.dta, replace 


/************************
Forecasting and analysis
*************************/ 
use final_dataset.dta, clear

tsset date

/***********
random walk
***********/ 

*Using a rolling window scheme the optimal forecast for random walk is the last period value
foreach v in eur jpy {
gen f_rw_`v' = L.`v'_exrate_nsa
gen f_rw_`v'_err = `v'_exrate_nsa - f_rw_`v'
egen rmsfe_rw_`v' = sd(f_rw_`v'_err)
gen ub_rw_`v' = f_rw_`v' + 0.674*rmsfe_rw_`v'
gen lb_rw_`v' = f_rw_`v' + 0.674*rmsfe_rw_`v'
egen rmspe_rw_`v' = sd(f_rw_`v'_err/`v'_exrate_nsa)
gen mspe_rw_`v' = 100*rmspe_rw_`v'^2 
egen mfe_rw_`v' = mean(f_rw_`v'_err*100/`v'_exrate_nsa)
}


/*********************************************
Taylor's rule & Incorporating VIX *For JPYUSD* 
**********************************************/

*creating the delta e to use taylor's rule from Ince (2013)
gen ljpy_exrate_nsa = log(jpy_exrate_nsa)
gen delta_s_jpy = D.ljpy_exrate_nsa
label var delta_s_jpy "Log differential of JPY USD Monthly Exchange Rate"

*components of taylor's rule
*creating inflation variables from CPI 
gen jpyinfl = D.jpy_cpi_nsa/ L.jpy_cpi_nsa
gen usainfl = D.usa_cpi_nsa / L.usa_cpi_nsa

*production gap variables
*on original industrial production
tsfilter hp jpy_prod_hp = jpy_prod, smooth(14400)
tsfilter hp usa_prod_hp = usa_prod, smooth(14400)

*to check whether cyclical component effectively captured
pergram jpy_prod_hp, xline(.03125 .16667)

*output gap in actual terms
gen jpy_prod_gap = jpy_prod - jpy_prod_hp
gen usa_prod_gap = usa_prod - usa_prod_hp

*create difference in vix 
gen diff_vix = D.vix

*our 3 models 

/**********************************************
Taylor's rule without coeff constraints + vix
**********************************************/
global prd = tm(2009m10) /*cut dataset to construct model selection + POOS*/

*construct forecasts of exchange rate from forecasts of delta s 

*Taylor's Rule 
reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap if date < $prd
predict f_s_tr, xb
predict err_s_tr, resid  
gen X = exp(f_s_tr)
gen f_jpyexrate_tr = X*L.jpy_exrate_nsa
label var f_jpyexrate_tr "Exchange rate forecast using TR, JPYUSD"

*Taylor's Rule + 1 lag of VIX
reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L.diff_vix if date < $prd
predict f_s_tr_vix1, xb 
predict err_s_vix1, resid  
gen Xvix1 = exp(f_s_tr_vix1)
gen f_jpyexrate_tr_vix1 = Xvix1*L.jpy_exrate_nsa
label var f_jpyexrate_tr_vix1 "Exchange rate forecast using TR + 1 lag of VIX, JPYUSD" 

*Taylor's Rule + 2 lag of VIX
reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L(1/2).diff_vix if date < $prd
predict f_s_tr_vix2, xb 
predict err_s_vix2, resid  
gen Xvix2 = exp(f_s_tr_vix2)
gen f_jpyexrate_tr_vix2 = Xvix2*L.jpy_exrate_nsa
label var f_jpyexrate_tr_vix2 "Exchange rate forecast using TR + 2 lags of VIX, JPYUSD" 

*Taylor's Rule + 3 lag of VIX
reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L(1/3).diff_vix if date < $prd
predict f_s_tr_vix3, xb 
predict err_s_vix3, resid  
gen Xvix3 = exp(f_s_tr_vix3)
gen f_jpyexrate_tr_vix3 = Xvix3*L.jpy_exrate_nsa
label var f_jpyexrate_tr_vix3 "Exchange rate forecast using TR + 3 lags of VIX, JPYUSD" 


*examine residuals
foreach v in tr vix1 vix2 vix3 {
tsline err_s_`v'
graph export "$graphs\\err_s_`v'_jpy.png", width(1450) height(1054) replace
corrgram err_s_`v'
ac err_s_`v'
graph export "$graphs\\ac_`v'_jpy.png", width(1450) height(1054) replace
}
*errors are serially autocorrelated ---- biased forecasts 

*Using newey to ensure HAC errors while testing for joint significance of coefficients
newey delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L(1/3).diff_vix if date < $prd, lag(4)

*granger causality test for all variables 
foreach var in L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap{
testparm `var'
}

testparm L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap

forval i = 1/3{
testparm L(1/`i').diff_vix
}



*PLS
gen p = _n
gen poos_pls = 1 if p >= 107 & p <161
gen trhat =.
gen vix1hat =.
gen vix2hat =.

gen trhat_exrate =.
gen vix1hat_exrate =.
gen vix2hat_exrate =.

local i = 0
forvalues t =107/160{
local i = `i' + 1
*creating pseudo OOS for PLS evaluation: TR only
qui reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap if p >= 3 + `i' & p <= `t'
predict fit_tr, xb
gen X_tr = exp(fit_tr)
gen pls_oos_jpy_tr = X_tr*L.jpy_exrate_nsa


*creating pseudo OOS for PLS evaluation: TR + 1 lag of VIX
qui reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L.diff_vix if p >= 3 + `i' & p <= `t'
predict fit_vix1, xb
gen X_vix1 = exp(fit_vix1)
gen pls_oos_jpy_vix1 = X_vix1*L.jpy_exrate_nsa

*creating pseudo OOS for PLS evaluation: TR + 2 lags of VIX
qui reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L(1/2).diff_vix if p >= 3 + `i' & p <= `t'
predict fit_vix2, xb
gen X_vix2 = exp(fit_vix2)
gen pls_oos_jpy_vix2 = X_vix2*L.jpy_exrate_nsa

*log differential forecast
replace trhat = fit_tr if p==(`t'+1)
replace vix1hat = fit_vix1 if p==(`t'+1)
replace vix2hat = fit_vix2 if p==(`t'+1)

*creating ex rate forecast from delta 
replace trhat_exrate = pls_oos_jpy_tr if p==(`t'+1)
replace vix1hat_exrate = pls_oos_jpy_vix1 if p==(`t'+1)
replace vix2hat_exrate = pls_oos_jpy_vix2 if p==(`t'+1)


drop fit* pls* X_* 
}

*using PLS instead of AIC because dependent variable (RW vs others) are different
*creating poos errors for PLS
gen tr_oos_err = trhat_exrate - jpy_exrate_nsa
gen vix1_oos_err = vix1hat_exrate - jpy_exrate_nsa
gen vix2_oos_err = vix2hat_exrate - jpy_exrate_nsa

*calculation of PLS for 3 models
foreach v in tr vix1 vix2{
egen `v'_oos_pls = sd(`v'_oos_err), by(poos_pls)
replace `v'_oos_pls =. if poos_pls ==.
}

*calculation of PLS for Random Walk 
egen rw_oos_pls = sd(f_rw_jpy_err), by(poos_pls)
replace rw_oos_pls =. if poos_pls==.

drop trhat vix1hat vix2hat tr_oos_err vix1_oos_err vix2_oos_err 
*using PLS against the models, ranked: RW < TR+2VIX < TR+1VIX < TR

*now, since PLS are quite similar, we want to do a forecast combination of the models
*using Granger-Ramanathan because possibility of bias with correlated errors
*using all available data until hold out period, 2014m4


global m = tm(2014m4)

*in-sample forecasts for TR
qui reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap if date <= $m
predict fit_tr, xb
gen X_tr = exp(fit_tr)
gen f_jpy_tr = X_tr*L.jpy_exrate_nsa

*in-sample forecasts for TR+VIX1
qui reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L.diff_vix if date<=$m
predict fit_vix1, xb
gen X_vix1 = exp(fit_vix1)
gen f_jpy_vix1 = X_vix1*L.jpy_exrate_nsa

*in-sample forecasts for TR+VIX2
qui reg delta_s_jpy  L.jpyinfl L.usainfl L.jpy_prod_gap L.usa_prod_gap L(1/2).diff_vix if date<=$m
predict fit_vix2, xb
gen X_vix2 = exp(fit_vix2)
gen f_jpy_vix2 = X_vix2*L.jpy_exrate_nsa

*Granger Ramanathan Forecast combination regression (because correlated forecast, biased/ unbiased, constant ?)
constraint 1 f_rw_jpy+trhat_exrate+vix1hat_exrate+vix2hat_exrate =1
*with constant, accounting for biased forecasts
cnsreg jpy_exrate_nsa f_rw_jpy trhat_exrate vix1hat_exrate vix2hat_exrate if date > $prd & date <= $m, constraints(1)

*without constant
cnsreg jpy_exrate_nsa f_rw_jpy trhat_exrate vix1hat_exrate vix2hat_exrate if date > $prd & date <= $m, constraints(1) nocons

*repeat after dropping tr, vix1
constraint 1 f_rw_jpy+vix2hat_exrate = 1
cnsreg jpy_exrate_nsa f_rw_jpy vix2hat_exrate if date > $prd & date <= $m, constraints(1) 
   
**********************************************************************************************
*OOS Forecasts (Assuming availabiity of forecasted values of taylor's rule components and Vix)
**********************************************************************************************


*1 month ahead direct forecast for TR + 2 lags of VIX
gen vix2_1mf =.
gen vix2_1mfL =.
gen vix2_1mfU =.

local k= 2 /* start of rolling window*/ 
local m = 160 /* end of rolling window*/

forval i=1/51 {
local j = `i' + 1
local k = `k' + 1
local m = `m' + 1

qui reg delta_s_jpy  L(`i').jpyinfl L(`i').usainfl L(`i').jpy_prod_gap L(`i').usa_prod_gap L(`i'/`j').diff_vix if p >= `k' & p < `m'
predict vix2_1mfit`i', xb
predict vix2_1m_se`i', stdf


replace vix2_1mf = vix2_1mfit`i' if p == `m' 
replace vix2_1mfL = vix2_1mfit`i' - vix2_1m_se`i'*0.674 if p == `m' 
replace vix2_1mfU = vix2_1mfit`i' + vix2_1m_se`i'*0.674 if p == `m' 

drop vix2_1mfit`i' vix2_1m_se`i'
}

label variable vix2_1mf "1-month ahead direct forecast"
label variable vix2_1mfL "1-month ahead direct forecast (50% lower bound)"
label variable vix2_1mfU "1-month ahead direct forecast (50% upper bound)"

tsline vix2_1mf delta_s_jpy vix2_1mfU vix2_1mfL if p > 160 & p < 212, lp(solid solid dash dash) legend(rows(4)) title("1 Month Ahead Direct Forecasts") subtitle("Using rolling window estimates") xtitle("") ytitle("JPY USD Log Difference")
graph export "$graphs\\1m_logdf_jpyusd.png", width(1450) height(1054) replace

*6 month ahead direct forecast for TR + 2 lags of VIX
gen vix2_6mf =.
gen vix2_6mfL =.
gen vix2_6mfU =.

local k= 7 /* start of rolling window*/ 
local m = 166 /* end of rolling window*/

forval i=1/45 {
local j = `i' + 1
local k = `k' + 1
local m = `m' + 1

qui reg delta_s_jpy  L(`i').jpyinfl L(`i').usainfl L(`i').jpy_prod_gap L(`i').usa_prod_gap L(`i'/`j').diff_vix if p >= `k' & p < `m'
predict vix2_6mfit`i', xb
predict vix2_6m_se`i', stdf


replace vix2_6mf = vix2_6mfit`i' if p == `m'
replace vix2_6mfL = vix2_6mfit`i' - vix2_6m_se`i'*0.674 if p == `m'
replace vix2_6mfU = vix2_6mfit`i' + vix2_6m_se`i'*0.674 if p == `m'

drop vix2_6mfit`i' vix2_6m_se`i'
}
label variable vix2_6mf "6-month ahead direct forecast"
label variable vix2_6mfL "6-month ahead direct forecast (50% lower bound)"
label variable vix2_6mfU "6-month ahead direct forecast (50% upper bound)"

tsline vix2_6mf delta_s_jpy vix2_6mfU vix2_6mfL if p > 166 & p < 212, lp(solid solid dash dash) legend(rows(4)) title("6 Month Ahead Direct Forecasts") subtitle("Using rolling window estimates") xtitle("") ytitle("JPY USD Log Difference")
graph export "$graphs\\6m_logdf_jpyusd.png", width(1450) height(1054) replace

*12 month ahead direct forecast for TR + 2 lags of VIX
gen vix2_12mf =.
gen vix2_12mfL =.
gen vix2_12mfU =.

local k= 14 /* start of rolling window*/ 
local m = 172 /* end of rolling window*/

forval i=1/39 {
local j = `i' + 1
local k = `k' + 1
local m = `m' + 1

qui reg delta_s_jpy  L(`i').jpyinfl L(`i').usainfl L(`i').jpy_prod_gap L(`i').usa_prod_gap L(`i'/`j').diff_vix if p >= `k' & p < `m'
predict vix2_12mfit`i', xb
predict vix2_12m_se`i', stdf


replace vix2_12mf = vix2_12mfit`i' if p == `m'
replace vix2_12mfL = vix2_12mfit`i' - vix2_12m_se`i'*0.674 if p == `m'
replace vix2_12mfU = vix2_12mfit`i' + vix2_12m_se`i'*0.674 if p == `m' 

drop vix2_12mfit`i' vix2_12m_se`i'
}

label variable vix2_12mf "12-month ahead direct forecast"
label variable vix2_12mfL "12-month ahead direct forecast (50% lower bound)"
label variable vix2_12mfU "12-month ahead direct forecast (50% upper bound)"

tsline vix2_12mf delta_s_jpy vix2_12mfU vix2_12mfL if p > 172 & p < 212, lp(solid solid dash dash) legend(rows(4)) title("12 Month Ahead Direct Forecasts") subtitle("Using rolling window estimates") xtitle("") ytitle("JPY USD Log Difference")
graph export "$graphs\\12m_logdf_jpyusd.png", width(1450) height(1054) replace

*transforming back the dependent variable to actual exchange rate
foreach i in 1 6 12{
gen X_vix2_`i' = exp(vix2_`i'mf)
gen X_vix2_`i'U = exp(vix2_`i'mfU)
gen X_vix2_`i'L = exp(vix2_`i'mfL)
gen exr_vix2_`i'mf = X_vix2_`i'*L.jpy_exrate_nsa /*using observed lagged value*/
label var exr_vix2_`i'mf "`i' Month Ahead Point Forecast, Transformed"
gen exr_vix2_`i'mfU = X_vix2_`i'U*L.jpy_exrate_nsa 
label var exr_vix2_`i'mfU "`i' Month Ahead Upper Bound, Transformed"
gen exr_vix2_`i'mfL = X_vix2_`i'L*L.jpy_exrate_nsa 
label var exr_vix2_`i'mfL "`i' Month Ahead Lower Bound, Transformed"
tsline jpy_exrate_nsa exr_vix2_`i'mf exr_vix2_`i'mfU exr_vix2_`i'mfL if p > 160 & p < 213, lp(solid solid dash dash) legend(rows(4)) title("`i' Month Ahead Direct Forecasts - Transformed") subtitle("Using rolling window estimates") xtitle("") ytitle("JPY USD")
graph export "$graphs\\`i'm_jpyusd.png", width(1450) height(1054) replace
}


/*Forecast evaluation*/
*chosen model based on above tests: TR + 2 lags of VIX, compare against random walk
*Using Diebold Mariano, with rolling window forecasts we can use the test

foreach i in 1m 6m 12m{
gen f_vix2_jpy_`i'err = jpy_exrate_nsa - exr_vix2_`i'f
gen d_jpy_`i'=f_rw_jpy_err^2 - f_vix2_jpy_`i'err^2  /*generate differential under squared loss*/

tsline d_jpy_`i' if p>160
ac d_jpy_`i'
corrgram d_jpy_`i'
}
*cannot reject the null at Q-statistic in sqrt(51/45), errors are robust

dmariano jpy_exrate_nsa f_rw_jpy exr_vix2_1mf, crit(MSE) kernel(bartlett) maxlag(0)

foreach i in 1m 6m{
reg d_jpy_`i', r

}
*for 12m, errors are not robust so we use newey with 3 lags
newey d_jpy_12m, lag(3)

**DM: if the coeff is negative, RW is better
**DM_1m: RW is better
**DM_6m: RW is better
**DM_12m: RW is better

*do the correction for small sample size
disp -2.61*sqrt(1+(1/51)*(1-2*1)+(1/51^2)*(1*(1-1)))
disp t(50,-2.5842851)*2
disp -1.79*sqrt(1+(1/45)*(1-2*6)+(1/45^2)*(6*(6-1)))
disp t(44,-1.5710963)*2
disp -1.50*sqrt(1+(1/39)*(1-2*12)+(1/39^2)*(12*(12-1)))
disp t(38,-1.0575175)*2

*can only reject null at RW & 1-month forecast comparison, the rest are inconclusive
*RW beats 1-month ahead forecasts for JPYUSD
drop f_s_tr_vix3 err_s_vix3 Xvix3
drop ljpy_exrate_nsa delta_s_jpy jpyinfl usainfl jpy_prod_hp usa_prod_hp jpy_prod_gap usa_prod_gap diff_vix f_s_tr err_s_tr X f_jpyexrate_tr f_s_tr_vix1 err_s_vix1 Xvix1 f_jpyexrate_tr_vix1 f_s_tr_vix2 err_s_vix2 Xvix2 f_jpyexrate_tr_vix2 p poos_pls trhat_exrate vix1hat_exrate vix2hat_exrate tr_oos_pls vix1_oos_pls vix2_oos_pls rw_oos_pls fit_tr X_tr f_jpy_tr fit_vix1 X_vix1 f_jpy_vix1 fit_vix2 X_vix2 f_jpy_vix2 vix2_1mf vix2_1mfL vix2_1mfU vix2_6mf vix2_6mfL vix2_6mfU vix2_12mf vix2_12mfL vix2_12mfU X_vix2_1 X_vix2_1U X_vix2_1L exr_vix2_1mf exr_vix2_1mfU exr_vix2_1mfL X_vix2_6 X_vix2_6U X_vix2_6L exr_vix2_6mf exr_vix2_6mfU exr_vix2_6mfL X_vix2_12 X_vix2_12U X_vix2_12L exr_vix2_12mf exr_vix2_12mfU exr_vix2_12mfL
drop f_vix2_jpy_1merr d_jpy_1m f_vix2_jpy_6merr d_jpy_6m f_vix2_jpy_12merr d_jpy_12m

/*********************************************
Taylor's rule + incorporating VIX *For EURUSD* 
**********************************************/

*creating the delta e to use taylor's rule from Ince (2013)
*data from FRED gives USDEUR, using EURUSD for ease of analysis
*using Germany as the proxy for Euro currency 
*complete dataset ends at 2018m12 
gen leur_exrate_nsa = log(eur_exrate_nsa) 
gen delta_s_eur = D.leur_exrate_nsa
label var delta_s_eur "Log differential of EURUSD Monthly Exchange Rate"

*components of taylor's rule
*creating inflation variables from CPI 
gen gerinfl = D.ger_cpi_nsa/ L.ger_cpi_nsa
gen usainfl = D.usa_cpi_nsa / L.usa_cpi_nsa

*production gap variables
*on original industrial production, Germany as proxy for Euro
tsfilter hp ger_prod_hp = ger_prod, smooth(14400)
tsfilter hp usa_prod_hp = usa_prod, smooth(14400)

*to check whether cyclical component effectively captured
pergram ger_prod_hp, xline(.03125 .16667)

*output gap in actual terms
gen ger_prod_gap = ger_prod - ger_prod_hp
gen usa_prod_gap = usa_prod - usa_prod_hp

*create difference in vix 
gen diff_vix = D.vix

*our 3 models 

/**********************************************
Taylor's rule without coeff constraints + vix
**********************************************/
global prd2 = tm(2009m12) /*cut dataset to construct model selection + POOS*/

*construct forecasts of exchange rate from forecasts of delta s 
*taylor's Rule
reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap if date < $prd2
predict f_s_tr, xb
predict err_s_tr, resid  
gen X = exp(f_s_tr)
gen f_eurexrate_tr = X*L.eur_exrate_nsa
label var f_eurexrate_tr "Exchange rate forecast using TR, EURUSD"

*Taylor's Rule + 1 lag of VIX
reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L.diff_vix if date < $prd2
predict f_s_tr_vix1, xb 
predict err_s_vix1, resid  
gen Xvix1 = exp(f_s_tr_vix1)
gen f_eurexrate_tr_vix1 = Xvix1*L.eur_exrate_nsa
label var f_eurexrate_tr_vix1 "Exchange rate forecast using TR + 1 lag of VIX, EURUSD" 

*Taylor's Rule + 2 lags of VIX
reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L(1/2).diff_vix if date < $prd2
predict f_s_tr_vix2, xb 
predict err_s_vix2, resid  
gen Xvix2 = exp(f_s_tr_vix2)
gen f_eurexrate_tr_vix2 = Xvix2*L.eur_exrate_nsa
label var f_eurexrate_tr_vix2 "Exchange rate forecast using TR + 2 lags of VIX, EURUSD" 

*Taylor's rule + 3 lags of VIX
reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L(1/3).diff_vix if date < $prd2
predict f_s_tr_vix3, xb 
predict err_s_vix3, resid  
gen Xvix3 = exp(f_s_tr_vix3)
gen f_eurexrate_tr_vix3 = Xvix3*L.eur_exrate_nsa
label var f_eurexrate_tr_vix3 "Exchange rate forecast using TR + 3 lags of VIX, EURUSD" 

*examine residuals
foreach v in tr vix1 vix2 vix3{
tsline err_s_`v'
graph export "$graphs\\err_s_`v'_eur.png", width(1450) height(1054) replace
corrgram err_s_`v'
ac err_s_`v'
graph export "$graphs\\ac_`v'_eur.png", width(1450) height(1054) replace
}
*sqrt(108) \approx 10, thus from Ljung Box Q test, errors have serial autocorrelation  ---- biased forecasts 

*Using newey to ensure HAC errors while testing for joint significance of coefficients
newey delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L(1/3).diff_vix if date < $prd2, lag(5)

*granger causality test for all variables 
foreach var in L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap{
testparm `var'
}

testparm L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap

forval i = 1/3{
testparm L(1/`i').diff_vix
}



*PLS
gen p = _n
gen poos_pls = 1 if p >= 109 & p <163
gen trhat =.
gen vix1hat =.
gen vix2hat =.

gen trhat_exrate =.
gen vix1hat_exrate =.
gen vix2hat_exrate =.

local i = 0
forvalues t =109/163{
local `i' = `i' + 1
*creating pseudo OOS for PLS evaluation: TR only
qui reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap if p >= 3 + `i' & p <= `t'
predict fit_tr, xb
gen X_tr = exp(fit_tr)
gen pls_oos_eur_tr = X_tr*L.eur_exrate_nsa


*creating pseudo OOS for PLS evaluation: TR + 1 lag of VIX
qui reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L.diff_vix if p >= 3 + `i' & p <= `t'
predict fit_vix1, xb
gen X_vix1 = exp(fit_vix1)
gen pls_oos_eur_vix1 = X_vix1*L.eur_exrate_nsa

*creating pseudo OOS for PLS evaluation: TR + 2 lags of VIX
qui reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L(1/2).diff_vix if p >= 3 + `i' & p <= `t'
predict fit_vix2, xb
gen X_vix2 = exp(fit_vix2)
gen pls_oos_eur_vix2 = X_vix2*L.eur_exrate_nsa

replace trhat = fit_tr if p==(`t'+1)
replace vix1hat = fit_vix1 if p==(`t'+1)
replace vix2hat = fit_vix2 if p==(`t'+1)

*creating ex rate forecast from delta 
replace trhat_exrate = pls_oos_eur_tr if p==(`t'+1)
replace vix1hat_exrate = pls_oos_eur_vix1 if p==(`t'+1)
replace vix2hat_exrate = pls_oos_eur_vix2 if p==(`t'+1)


drop fit* pls* X_* 
}

*using PLS instead of AIC because dependent variable (RW vs others) are different
*creating oos errors for PLS
gen tr_oos_err = trhat_exrate - eur_exrate_nsa
gen vix1_oos_err = vix1hat_exrate - eur_exrate_nsa
gen vix2_oos_err = vix2hat_exrate - eur_exrate_nsa

*calculation of PLS for 3 models
foreach v in tr vix1 vix2{
egen `v'_oos_pls = sd(`v'_oos_err), by(poos_pls)
replace `v'_oos_pls =. if poos_pls ==.
}

*calculation of PLS for Random Walk 
egen rw_oos_pls = sd(f_rw_eur_err), by(poos_pls)
replace rw_oos_pls =. if poos_pls==.

drop trhat vix1hat vix2hat tr_oos_err vix1_oos_err vix2_oos_err 
*using PLS against the models, ranked: RW < TR+1VIX < TR+2VIX < TR

*now, since PLS are quite similar, we want to do a forecast combination of the models
*using Granger-Ramanathan because possibility of bias with correlated errors
*using all available data until hold out period, 2014m6

global m = tm(2014m6)

*in-sample forecasts for TR
qui reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap if date <= $m
predict fit_tr, xb
gen X_tr = exp(fit_tr)
gen f_eur_tr = X_tr*L.eur_exrate_nsa

*in-sample forecasts for TR+VIX1
qui reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L.diff_vix if date<=$m
predict fit_vix1, xb
gen X_vix1 = exp(fit_vix1)
gen f_eur_vix1 = X_vix1*L.eur_exrate_nsa

*in-sample forecasts for TR+VIX2
qui reg delta_s_eur  L.gerinfl L.usainfl L.ger_prod_gap L.usa_prod_gap L(1/2).diff_vix if date<=$m
predict fit_vix2, xb
gen X_vix2 = exp(fit_vix2)
gen f_eur_vix2 = X_vix2*L.eur_exrate_nsa

*Granger Ramanathan Forecast combination regression (because correlated forecast, biased/ unbiased, constant ?)
constraint 1 f_rw_eur+trhat_exrate+vix1hat_exrate+vix2hat_exrate =1
*with constant, accounting for biased forecasts
cnsreg eur_exrate_nsa f_rw_eur trhat_exrate vix1hat_exrate vix2hat_exrate if date > $prd2 & date <= $m, constraints(1)

*repeat after dropping tr and vix2 model
constraint 1 f_rw_eur+vix1hat_exrate = 1
cnsreg eur_exrate_nsa f_rw_eur vix1hat_exrate if date > $prd2 & date <= $m, constraints(1)

*picks RW + VIX1
   
**********************************************************************************************
*OOS Forecasts (Assuming availabiity of forecasted values of taylor's rule components)
**********************************************************************************************

*1 month ahead direct forecast for TR+1 lag of VIX only
gen vix1_1mf =.
gen vix1_1mfL =.
gen vix1_1mfU =.

local k= 2 /* start of rolling window*/ 
local m = 162 /* end of rolling window*/

forval i=1/55 {

local k = `k' + 1
local m = `m' + 1

qui reg delta_s_eur  L(`i').gerinfl L(`i').usainfl L(`i').ger_prod_gap L(`i').usa_prod_gap L(`i').diff_vix if p >= `k' & p < `m'
predict vix1_1mfit`i', xb
predict vix1_1m_se`i', stdf


replace vix1_1mf = vix1_1mfit`i' if p == `m' 
replace vix1_1mfL = vix1_1mfit`i' - vix1_1m_se`i'*0.674 if p == `m' 
replace vix1_1mfU = vix1_1mfit`i' + vix1_1m_se`i'*0.674 if p == `m' 

drop vix1_1mfit`i' vix1_1m_se`i'
}

label variable vix1_1mf "1-month ahead direct forecast"
label variable vix1_1mfL "1-month ahead direct forecast (50% lower bound)"
label variable vix1_1mfU "1-month ahead direct forecast (50% upper bound)"

tsline vix1_1mf delta_s_eur vix1_1mfU vix1_1mfL if p > 162 & p < 218, lp(solid solid dash dash) legend(rows(4)) title("1 Month Ahead Direct Forecasts") subtitle("Using rolling window estimates") xtitle("") ytitle("EUR USD Log Difference")
graph export "$graphs\\1m_logdf_eurusd.png", width(1450) height(1054) replace

*6 month ahead direct forecast for TR + 2 lags of VIX
gen vix1_6mf =.
gen vix1_6mfL =.
gen vix1_6mfU =.

local k= 8 /* start of rolling window*/ 
local m = 168 /* end of rolling window*/

forval i=1/49 {

local k = `k' + 1
local m = `m' + 1

qui reg delta_s_eur  L(`i').gerinfl L(`i').usainfl L(`i').ger_prod_gap L(`i').usa_prod_gap L(`i').diff_vix if p >= `k' & p < `m'
predict vix1_6mfit`i', xb
predict vix1_6m_se`i', stdf


replace vix1_6mf = vix1_6mfit`i' if p == `m'
replace vix1_6mfL = vix1_6mfit`i' - vix1_6m_se`i'*0.674 if p == `m'
replace vix1_6mfU = vix1_6mfit`i' + vix1_6m_se`i'*0.674 if p == `m'

drop vix1_6mfit`i' vix1_6m_se`i'
}
label variable vix1_6mf "6-month ahead direct forecast"
label variable vix1_6mfL "6-month ahead direct forecast (50% lower bound)"
label variable vix1_6mfU "6-month ahead direct forecast (50% upper bound)"

tsline vix1_6mf delta_s_eur vix1_6mfU vix1_6mfL if p > 168 & p < 218, lp(solid solid dash dash) legend(rows(4)) title("6 Month Ahead Direct Forecasts") subtitle("Using rolling window estimates") xtitle("") ytitle("EUR USD Log Difference")
graph export "$graphs\\6m_logdf_eurusd.png", width(1450) height(1054) replace

*12 month ahead direct forecast for TR
gen vix1_12mf =.
gen vix1_12mfL =.
gen vix1_12mfU =.

local k= 14 /* start of rolling window*/ 
local m = 174 /* end of rolling window*/

forval i=1/43 {

local k = `k' + 1
local m = `m' + 1

qui reg delta_s_eur  L(`i').gerinfl L(`i').usainfl L(`i').ger_prod_gap L(`i').usa_prod_gap L(`i').diff_vix if p >= `k' & p < `m'
predict vix1_12mfit`i', xb
predict vix1_12m_se`i', stdf


replace vix1_12mf = vix1_12mfit`i' if p == `m'
replace vix1_12mfL = vix1_12mfit`i' - vix1_12m_se`i'*0.674 if p == `m'
replace vix1_12mfU = vix1_12mfit`i' + vix1_12m_se`i'*0.674 if p == `m' 

drop vix1_12mfit`i' vix1_12m_se`i'
}

label variable vix1_12mf "12-month ahead direct forecast"
label variable vix1_12mfL "12-month ahead direct forecast (50% lower bound)"
label variable vix1_12mfU "12-month ahead direct forecast (50% upper bound)"

tsline vix1_12mf delta_s_eur vix1_12mfU vix1_12mfL if p > 174 & p < 218, lp(solid solid dash dash) legend(rows(4)) title("12 Month Ahead Direct Forecasts") subtitle("Using rolling window estimates") xtitle("") ytitle("EUR USD Log Difference")
graph export "$graphs\\12m_logdf_eurusd.png", width(1450) height(1054) replace

*transforming back the dependent variable to actual exchange rate
foreach i in 1 6 12{
gen X_vix1_`i' = exp(vix1_`i'mf)
gen X_vix1_`i'U = exp(vix1_`i'mfU)
gen X_vix1_`i'L = exp(vix1_`i'mfL)
gen exr_vix1_`i'mf = X_vix1_`i'*L.eur_exrate_nsa /*using observed lag value*/
label var exr_vix1_`i'mf "`i' Month Ahead Point Forecast, Transformed"
gen exr_vix1_`i'mfU = X_vix1_`i'U*L.eur_exrate_nsa
label var exr_vix1_`i'mfU "`i' Month Ahead Upper Bound, Transformed"
gen exr_vix1_`i'mfL = X_vix1_`i'L*L.eur_exrate_nsa
label var exr_vix1_`i'mfL "`i' Month Ahead Lower Bound, Transformed"
tsline eur_exrate_nsa exr_vix1_`i'mf exr_vix1_`i'mfU exr_vix1_`i'mfL if p > 162 & p < 219, lp(solid solid dash dash) legend(rows(4)) title("`i' Month Ahead Direct Forecasts - Transformed") subtitle("Using rolling window estimates") xtitle("") ytitle("EUR USD")
graph export "$graphs\\`i'm_eurusd.png", width(1450) height(1054) replace
}

/*Forecast evaluation*/
*chosen model based on above tests: TR, compare against random walk
*Using Diebold Mariano, with rolling window forecasts we can use the test
*as we have already determined, the models all have errors that are serially
*autocorrelated, so we will use newey 

foreach i in 1m 6m 12m{
gen f_vix1_eur_`i'err = eur_exrate_nsa - exr_vix1_`i'f
gen d_eur_`i'=f_rw_eur_err^2 - f_vix1_eur_`i'err^2  /*generate differential under squared loss*/

tsline d_eur_`i' if p>160
ac d_eur_`i'
corrgram d_eur_`i'
}
*cannot reject the null at Q-statistic in sqrt(55/49/43), errors are robust

dmariano eur_exrate_nsa f_rw_eur exr_vix1_1mf, crit(MSE) kernel(bartlett) maxlag(0)

foreach i in 1m 6m 12m{
reg d_eur_`i', r
}

**DM: if the coeff is negative, RW is better
**DM_1m: RW is better
**DM_6m: RW is better
**DM_12m: RW is better
*but t/n: all coeffs are VERY close to zero (esp the 1-month ahead forecasts)

*do the correction for small sample size
disp -0.92*sqrt(1+(1/55)*(1-2*1)+(1/55^2)*(1*(1-1)))
disp t(54,-.55488574)*2
disp -0.15*sqrt(1+(1/49)*(1-2*6)+(1/49^2)*(6*(6-1)))
disp t(48,-.79892681)*2
disp -1.13*sqrt(1+(1/43)*(1-2*12)+(1/43^2)*(12*(12-1)))
disp t(42,-.50540144)*2

*cannot reject the null for ALL forecasts
*RW & VIX1 model forecast performance are not significantly different from each other

