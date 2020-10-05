clear all
. use https://dl.dropboxusercontent.com/u/359550/booij2016/work_base
/*
multiple line comments like this */

. summarize male age ist_raw fes_raw cito_raw treat matched

// Create new variables:


gen ist_normsq = ist_norm^2
label var ist_normsq "normalized IST score squared"

gen ist_normcub = ist_norm^3
label var ist_normcub "normalized IST score cubed"

gen istnorm_z1 = ist_norm*z

gen istnorm_z0 = ist_norm*(z-1)

gen istnormsq_z1 = ist_normsq*z

gen istnormsq_z0 = ist_normsq*(z-1)

egen fes = std(fes_raw)
label var fes "standardized FES score"

egen cito = std(cito_raw) , std(.87)
label var cito "standardized CITO score"

egen math = std(reswis_raw) 
label var math "standardized average math score grades 2 - 6"

egen language = std(reslang_raw) 
label var language "standardized average language score grades 2 - 6"

egen other = std(resother_raw) 
label var other "standardized average other subjects score grades 2 - 6"

generate lnw = ln(startwage_raw) 


******************************************************************************

// TABLE 1 

label var male "Male"
label var age "Age"
label var treat "GT Program"
label var matched "Matched record"
label var drop "Retention"
label var reswis_raw "GPA math"
label var reslang_raw "GPA language"
label var resother_raw "GPA other"
label var nsubj56 "Avg no. subjects"
label var sctrack "Avg no. science subjects"
label var iwisB "Advanced math"
label var matchedhe "Matched record HE"
label var chosesc "Chose science field"
label var switch1yr "Switched fields"
label var startwage_raw "Predicted earnings"

estpost summarize male age ist_raw fes_raw cito_raw treat matched drop ///
                  reswis_raw reslang_raw resother_raw nsubj56 sctrack iwisB ///
				  matchedhe chosesc switch1yr startwage_raw , quietly

esttab , varwidth(40) cell((min(fmt(%9.2f)) max(fmt(%9.2f)) ///
				  mean(fmt(%9.2f)) sd(fmt(%9.2f)))) label nonumber ///
				  title("Table 1a")

estpost tabulate cohort , quietly

esttab , cell(b (label (N))) nonumber label title("Table 1b")

****************************************************************
// TABLE 2

// - Baseline Quadratic
qui reg treat z ist_norm ist_normsq male age cito fes i.cohort , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_quad

// - No controls
qui reg treat z ist_norm ist_normsq i.cohort , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_quad_nox

// - Linear
qui reg treat z ist_norm male age cito fes i.cohort , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_linear

// - Cubic
qui reg treat z ist_norm ist_normsq ist_normcub male age cito fes i.cohort , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_cub

// - Quadratic split
qui reg treat z istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 male age cito fes i.cohort , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_quadsplt

// - Zoom
qui reg treat z istnorm_z1 istnorm_z0 male age cito fes i.cohort if ist_norm>=-18 & ist_norm<=18 , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_zoom

// - Donut
qui reg treat z ist_norm ist_normsq male age cito fes i.cohort if ist_norm>4 | ist_norm<-4 , cluster(group)
test z
estadd scalar Pvalue = r(p)
estadd scalar Fstat = r(F)
estimates store fs_donut

// - Balancing: Male
qui reg male z ist_norm ist_normsq i.cohort , robust cluster(group)
estimates store bal_male

// - Balancing: Age
qui reg age z ist_norm ist_normsq i.cohort , robust cluster(group)
estimates store bal_age

// - Balancing: FES
qui reg fes z ist_norm ist_normsq i.cohort , robust cluster(group)
estimates store bal_fes

// - Balancing: CITO
qui reg cito z ist_norm ist_normsq i.cohort , robust cluster(group)
estimates store bal_cito

estadd ysumm: * 

esttab fs_quad fs_quad_nox fs_linear fs_cub fs_quadsplt fs_zoom fs_donut bal_male bal_age bal_fes bal_cito ///
	   , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   scalars(ymean ysd Pvalue Fstat r2) obslast sfmt(%9.2f) ///
       drop(ist_norm* istnorm* ist_normcub male age cito fes *cohort _cons) ///
       title("Table 2") ///
	   mtitles("Quadratic" "No controls" "Linear" "Cubic" "Quad. Split" "Zoom" "Donut" "Male" "Age" "FES" "CITO") 
	   
estimates clear
   
**************************************

// TABLE 3 

// - Matched
qui reg treat z matched male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls matched (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 260.4 
estimates store rf_matched

// - Retention
qui reg treat z drop male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls drop (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 261.9 
estimates store rf_drop

// - GPA Math
qui reg treat z math male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls math (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 250.7 
estimates store rf_GPAmath

// - GPA Language
qui reg treat z language male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls language (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 256.0 
estimates store rf_GPAlang

// - GPA Other
qui reg treat z other male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls other (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 252.5 
estimates store rf_GPAother

// - No. Subjects
qui reg treat z nsubj56 male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 133.1 
estimates store rf_subj

// - No. Science Subjects
qui reg treat z sctrack male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 132.1  
estimates store rf_scsubj

// - Advanced Math
qui reg treat z iwisB male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
test z

qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estadd ysumm
test treat
estadd scalar Pvalue = r(p)
estadd scalar FS_Fstat = 131.7  
estimates store rf_advmath

esttab rf_matched rf_drop rf_GPAmath rf_GPAlang rf_GPAother rf_subj rf_scsubj rf_advmath ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux scalars(ymean ysd Pvalue FS_Fstat r2_a) ///
       sfmt(%9.2f) drop(*cohort ist_norm* _cons) obslast ///
       title("Table 3") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

*******************************************************************************
	 
// TABLE 4  

// - Females only:

// 		- Sample Selection: Matched
qui ivregress 2sls matched (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_matched_fem

// 		- Sample Selection: Retention
qui ivregress 2sls drop (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_drop_fem

// 		- GPA Math
qui ivregress 2sls math (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_GPAmath_fem

// 		- GPA Language
qui ivregress 2sls language (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_GPAlang_fem

// 		- GPA Other
qui ivregress 2sls other (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_GPAother_fem

// 		- No. Subjects
qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_subj_fem

// 		- No. Science Subjects
qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_scsubj_fem

// 		- Advanced Math
qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_advmath_fem


// - Males only:

// 		- Sample Selection: Matched
qui ivregress 2sls matched (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_matched_mal

// 		- Sample Selection: Retention
qui ivregress 2sls drop (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_drop_mal

// 		- GPA Math
qui ivregress 2sls math (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_GPAmath_mal

// 		- GPA Language
qui ivregress 2sls language (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_GPAlang_mal

// 		- GPA Other
qui ivregress 2sls other (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_GPAother_mal

// 		- No. Subjects
qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_subj_mal

// 		- No. Science Subjects
qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_scsubj_mal

// 		- Advanced Math
qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_advmath_mal

// Table 4 (in three parts):

esttab rf_matched rf_drop rf_GPAmath rf_GPAlang rf_GPAother rf_subj rf_scsubj rf_advmath ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
       drop(male age fes cito *cohort ist_norm* _cons) noobs ///
       title("Table 4a") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

esttab rf_matched_fem rf_drop_fem rf_GPAmath_fem rf_GPAlang_fem rf_GPAother_fem rf_subj_fem rf_scsubj_fem rf_advmath_fem ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
       drop(male age fes cito *cohort ist_norm* _cons) noobs ///
       title("Table 4b") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

esttab rf_matched_mal rf_drop_mal rf_GPAmath_mal rf_GPAlang_mal rf_GPAother_mal rf_subj_mal rf_scsubj_mal rf_advmath_mal ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
       drop(male age fes cito *cohort ist_norm* _cons) noobs ///
       title("Table 4c") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

estimates clear

*******************************************************************************
	   
// TABLE 5

// Baseline (N=2438)

// 		- Matched HE
qui ivregress 2sls matchedhe (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_matchedhe

// 		- Sc. field
qui ivregress 2sls chosesc (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_chosesc

//		- FY Switch
qui ivregress 2sls switch1yr (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_switch1yr

//		- log(w)
qui ivregress 2sls lnw (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_lnw

// Females (N=1135)

// 		- Matched HE
qui ivregress 2sls matchedhe (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_matchedhe_fem

// 		- Sc. field
qui ivregress 2sls chosesc (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_chosesc_fem

//		- FY Switch
qui ivregress 2sls switch1yr (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_switch1yr_fem

//		- log(w)
qui ivregress 2sls lnw (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==0 , cluster(group)
estimates store rf_lnw_fem

// Males only (N=1303)

// 		- Matched HE
qui ivregress 2sls matchedhe (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_matchedhe_mal

// 		- Sc. field
qui ivregress 2sls chosesc (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_chosesc_mal

//		- FY Switch
qui ivregress 2sls switch1yr (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_switch1yr_mal

//		- log(w)
qui ivregress 2sls lnw (treat = z) male age fes cito i.cohort ist_norm ist_normsq if male==1 , cluster(group)
estimates store rf_lnw_mal

// Table 5 (in three parts):

esttab rf_matchedhe rf_chosesc rf_switch1yr rf_lnw ///
	   , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
       drop(male age fes cito *cohort ist_norm* _cons) ///
       title("Table 5a") ///
	   mtitles("Matched HE" "Sc. Field" "FY Switch" "log(w)")
	   
   
esttab rf_matchedhe_fem rf_chosesc_fem rf_switch1yr_fem rf_lnw_fem ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   drop(male age fes cito *cohort ist_norm* _cons) ///
       title("Table 5b") ///
	   mtitles("Matched HE" "Sc. Field" "FY Switch" "log(w)")

esttab rf_matchedhe_mal rf_chosesc_mal rf_switch1yr_mal rf_lnw_mal ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   drop(male age fes cito *cohort ist_norm* _cons) ///
       title("Table 5c") ///
	   mtitles("Matched HE" "Sc. Field" "FY Switch" "log(w)")

estimates clear	   

*******************************************************************************
	   
// TABLE A1

// Baseline

qui ivregress 2sls matched (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_matched

qui ivregress 2sls drop (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_drop

qui ivregress 2sls math (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_GPAmath

qui ivregress 2sls language (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_GPAlang

qui ivregress 2sls other (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_GPAother

qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_subj

qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_scsubj

qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort ist_norm ist_normsq , cluster(group)
estimates store rf_advmath

// Quadratic Split

qui ivregress 2sls matched (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_matched_qs

qui ivregress 2sls drop (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_drop_qs

qui ivregress 2sls math (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_GPAmath_qs

qui ivregress 2sls language (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_GPAlang_qs

qui ivregress 2sls other (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_GPAother_qs

qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_subj_qs

qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_scsubj_qs

qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 istnormsq_z1 istnormsq_z0 , cluster(group)
estimates store rf_advmath_qs

// Zoom

qui ivregress 2sls matched (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_matched_zm

qui ivregress 2sls drop (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_drop_zm

qui ivregress 2sls math (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_GPAmath_zm

qui ivregress 2sls language (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_GPAlang_zm

qui ivregress 2sls other (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_GPAother_zm

qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_subj_zm

qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_scsubj_zm

qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 , cluster(group)
estimates store rf_advmath_zm

// Donut

qui ivregress 2sls matched (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_matched_dn

qui ivregress 2sls drop (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_drop_dn

qui ivregress 2sls math (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_GPAmath_dn

qui ivregress 2sls language (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_GPAlang_dn

qui ivregress 2sls other (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_GPAother_dn

qui ivregress 2sls nsubj56 (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_subj_dn

qui ivregress 2sls sctrack (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_scsubj_dn

qui ivregress 2sls iwisB (treat = z) male age fes cito i.cohort ist_norm ist_normsq if ist_norm>4 | ist_norm<-4 , cluster(group)
estimates store rf_advmath_dn

// Table A1 (in four parts):

esttab rf_matched rf_drop rf_GPAmath rf_GPAlang rf_GPAother rf_subj rf_scsubj rf_advmath ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   drop(male age fes cito *cohort ist_norm* _cons) noobs ///
       title("Table A1a") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")
	   
esttab rf_matched_qs rf_drop_qs rf_GPAmath_qs rf_GPAlang_qs rf_GPAother_qs rf_subj_qs rf_scsubj_qs rf_advmath_qs ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   drop(male age fes cito *cohort istnorm* _cons) noobs ///
       title("Table A1b") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

esttab rf_matched_zm rf_drop_zm rf_GPAmath_zm rf_GPAlang_zm rf_GPAother_zm rf_subj_zm rf_scsubj_zm rf_advmath_zm ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   drop(male age fes cito *cohort istnorm* _cons) noobs ///
       title("Table A1c") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

esttab rf_matched_dn rf_drop_dn rf_GPAmath_dn rf_GPAlang_dn rf_GPAother_dn rf_subj_dn rf_scsubj_dn rf_advmath_dn ///
       , b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) staraux ///
	   drop(male age fes cito *cohort ist_norm* _cons) noobs ///
       title("Table A1d") ///
	   mtitles("Matched" "Retention" "Math" "Language" "Other" "No. subjects" "No. science subjects" "Adv. math")

estimates clear

*******************************************************************************

// FIGURE 1

qui reg treat z ist_norm 
predict fs_lin

qui reg treat z ist_norm ist_normsq  
predict fs_quad
	   
qui reg treat z ist_norm ist_normsq ist_normcub 
predict fs_cub

qui reg treat z istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 
predict fs_zoom

qui reg treat z ist_norm ist_normsq if ist_norm<-4 | ist_norm>4  
predict fs_donut
 
sum ist_norm
di (70+43)/4 // yields number of bins, based on range size divided by bin size
xtile x = ist_norm , n(28)
egen m1 = mean(treat), by(x) // take mean of each bin
egen m2 = tag(m1) // yields one observation per group

twoway (scatter m1 ist_norm if m2==1 & (inrange(ist_norm, -70,-0.1)), mcolor(white) mlcolor(black) ysc(noextend))  /// Local average 
       (scatter m1 ist_norm if m2==1 & (inrange(ist_norm, 0,43)), mcolor(white) mlcolor(black) ysc(noextend))  /// 
       (lfit fs_lin ist_norm if fs_lin>0 & (inrange(ist_norm, -70,-0.1)), lpattern(dash) lcolor(black)) /// Linear 
       (lfit fs_lin ist_norm if fs_lin<1 & (inrange(ist_norm, 0,43)), lpattern(dash) lcolor(black)) /// 
	   (qfit fs_quad ist_norm if fs_quad>0 & (inrange(ist_norm, -70,-0.1)), lcolor(black)) /// Quadratic 
	   (qfit fs_quad ist_norm if fs_quad<1 & (inrange(ist_norm, 0,43)), lcolor(black)) /// 
	   (fpfit fs_cub ist_norm if fs_cub>0 & (inrange(ist_norm, -70,-0.1)), est(degree(3)) lpattern(dot) lcolor(black)) /// Cubic 
	   (fpfit fs_cub ist_norm if fs_cub<1 & (inrange(ist_norm, 0,43)), est(degree(3)) lpattern(dot) lcolor(black))   /// 
	   (lfit fs_zoom ist_norm if fs_zoom>0 & (inrange(ist_norm, -18,-0.1)), lcolor(blue) lwidth(medthick)) /// Zoom 
	   (lfit fs_zoom ist_norm if fs_zoom<1 & (inrange(ist_norm, 0,18)), lcolor(blue) lwidth(medthick)) /// 
	   (qfit fs_donut ist_norm if fs_donut>0 & (inrange(ist_norm, -70,-4)), lcolor(gray) lwidth(medthick)) /// Donut 
	   (qfit fs_donut ist_norm if fs_donut<1 & (inrange(ist_norm, 4,43)), lcolor(gray) lwidth(medthick)) /// 
	   , xline(0 , lpattern(dash) lcolor(red)) ///
	   title("Figure 1") ytitle("Fraction Assigned") xtitle("") ///
	   legend(ring(0) pos(11) row(6) order(1 3 5 7 9 11) label(1 "Local average") ///
	   label(3 "Linear fit") label(5 "Quadratic fit") label(7 "Cubic fit") ///
	   label(9 "Zoom") label(11 "Donut")) name(panel1,replace) nodraw
	   
histogram ist_norm , width(4) start(-72) fraction color(black) fcolor(gray) ///
          xtitle("IST distance to threshold") ///
		  xline(0 , lpattern(dash)) name(panel2,replace) nodraw

graph combine panel1 panel2 , xcommon row(2) imargin(0 0 0 0) iscale(0.7) 


*******************************************************************************

// FIGURE 2 

estimates clear
drop m1 m2

//  - FES:

qui reg fes z ist_norm  
predict fes_lin

qui reg fes z ist_norm ist_normsq 
predict fes_quad
	   
qui reg fes z ist_norm ist_normsq ist_normcub 
predict fes_cub

qui reg fes z istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 
predict fes_zoom

qui reg fes z ist_norm ist_normsq if ist_norm<-4 | ist_norm>4 
predict fes_donut
 
egen m1 = mean(fes), by(x) 
egen m2 = tag(m1) 


twoway (scatter m1 ist_norm if m2==1 & (inrange(ist_norm, -65,-0.1)), mcolor(white) mlcolor(black))  /// Local average 
       (scatter m1 ist_norm if m2==1 & (inrange(ist_norm, 0,50)), mcolor(white) mlcolor(black))  /// 
       (lfit fes_lin ist_norm if (inrange(ist_norm, -65,-0.1)), lpattern(dash) lcolor(black)) /// Linear 
       (lfit fes_lin ist_norm if (inrange(ist_norm, 0,50)), lpattern(dash) lcolor(black)) /// 
	   (qfit fes_quad ist_norm if (inrange(ist_norm, -65,-0.1)), lcolor(black)) /// Quadratic 
	   (qfit fes_quad ist_norm if (inrange(ist_norm, 0,50)), lcolor(black)) /// 
	   (fpfit fes_cub ist_norm if (inrange(ist_norm, -65,-0.1)), est(degree(3)) lpattern(dot) lcolor(black)) /// Cubic 
	   (fpfit fes_cub ist_norm if (inrange(ist_norm, 0,50)), est(degree(3)) lpattern(dot) lcolor(black))   /// 
	   (lfit fes_zoom ist_norm if (inrange(ist_norm, -18,-0.1)), lcolor(blue) lwidth(medthick)) /// Zoom 
	   (lfit fes_zoom ist_norm if (inrange(ist_norm, 0,18)), lcolor(blue) lwidth(medthick)) /// 
	   (qfit fes_donut ist_norm if (inrange(ist_norm, -65,-4)), lcolor(gray) lwidth(medthick)) /// Donut
	   (qfit fes_donut ist_norm if (inrange(ist_norm, 4,50)), lcolor(gray) lwidth(medthick)) /// 
	   , xline(0 , lpattern(dash) lcolor(red)) yscale(r(-2 2))  ///
	   title("FES") ytitle("Standardized pre-test score") xtitle("IST distance to threshold") ///
	   legend(ring(0) pos(11) row(6) order(1 3 5 7 9 11) label(1 "Local average") ///
	   label(3 "Linear fit") label(5 "Quadratic fit") label(7 "Cubic fit") ///
	   label(9 "Zoom") label(11 "Donut")) name(panelfes,replace) nodraw

//  CITO:

qui reg cito z ist_norm 
predict cito_lin

qui reg cito z ist_norm ist_normsq 
predict cito_quad
	   
qui reg cito z ist_norm ist_normsq ist_normcub 
predict cito_cub

qui reg cito z istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 
predict cito_zoom

qui reg cito z ist_norm ist_normsq if ist_norm<-4 | ist_norm>4  
predict cito_donut
 
egen m3 = mean(cito), by(x) 
egen m4 = tag(m1) 

twoway (scatter m3 ist_norm if m4==1 & (inrange(ist_norm, -65,-0.1)), mcolor(white) mlcolor(black))  /// Local average 
       (scatter m3 ist_norm if m4==1 & (inrange(ist_norm, 0,50)), mcolor(white) mlcolor(black))  /// 
       (lfit cito_lin ist_norm if (inrange(ist_norm, -65,-0.1)), lpattern(dash) lcolor(black)) /// Linear 
       (lfit cito_lin ist_norm if (inrange(ist_norm, 0,50)), lpattern(dash) lcolor(black)) /// 
	   (qfit cito_quad ist_norm if (inrange(ist_norm, -65,-0.1)), lcolor(black)) /// Quadratic 
	   (qfit cito_quad ist_norm if (inrange(ist_norm, 0,50)), lcolor(black)) /// 
	   (fpfit cito_cub ist_norm if (inrange(ist_norm, -65,-0.1)), est(degree(3)) lpattern(dot) lcolor(black)) /// Cubic 
	   (fpfit cito_cub ist_norm if (inrange(ist_norm, 0,50)), est(degree(3)) lpattern(dot) lcolor(black))   /// 
	   (lfit cito_zoom ist_norm if (inrange(ist_norm, -18,-0.1)), lcolor(blue) lwidth(medthick)) /// Zoom 
	   (lfit cito_zoom ist_norm if (inrange(ist_norm, 0,18)), lcolor(blue) lwidth(medthick)) /// 
	   (qfit cito_donut ist_norm if (inrange(ist_norm, -65,-4)), lcolor(gray) lwidth(medthick)) /// Donut 
	   (qfit cito_donut ist_norm if (inrange(ist_norm, 4,50)), lcolor(gray) lwidth(medthick)) /// 
	   , xline(0 , lpattern(dash) lcolor(red)) yscale(r(-2 2))  ///
	   title("CITO") ytitle("Standardized pre-test score") xtitle("IST distance to threshold") ///
	   legend(off) name("panelcito") nodraw

graph combine panelfes panelcito , ycommon imargin(0 0 0 0)


*******************************************************************************

// FIGURE 3 

estimates clear
drop m1 m2 m3 m4

foreach i in math language other {
qui reg `i' z ist_norm 
predict lin

qui reg `i' z ist_norm ist_normsq 
predict quad

qui reg `i' z ist_norm ist_normsq ist_normcub 
predict cub

qui reg `i' z istnorm_z1 istnorm_z0 if ist_norm>=-18 & ist_norm<=18 
predict zoom

qui reg `i' z ist_norm ist_normsq if ist_norm<-4 | ist_norm>4  
predict donut

egen m1= mean(`i'), by(x)
egen m2=tag(m1)

twoway (scatter m1 ist_norm if m2==1 & (inrange(ist_norm, -70,0)), mcolor(white) mlcolor(black) ysc(noextend))  /// Local average 
       (scatter m1 ist_norm if m2==1 & (inrange(ist_norm, 0,50)), mcolor(white) mlcolor(black) ysc(noextend))  /// 
       (lfit lin ist_norm if (inrange(ist_norm, -70,0)), lpattern(dash) lcolor(black)) /// Linear 
       (lfit lin ist_norm if (inrange(ist_norm, 0,50)), lpattern(dash) lcolor(black)) /// 
	   (qfit quad ist_norm if (inrange(ist_norm, -70,0)), lcolor(black)) /// Quadratic 
	   (qfit quad ist_norm if (inrange(ist_norm, 0,50)), lcolor(black)) /// 
	   (fpfit cub ist_norm if (inrange(ist_norm, -70,0)), lpattern(dot) lcolor(black)) /// Cubic 
	   (fpfit cub ist_norm if (inrange(ist_norm, 0,50)), est(degree(3)) lpattern(dot) lcolor(black))   /// 
	   (lfit zoom ist_norm if (inrange(ist_norm, -18,0)), lcolor(blue) lwidth(medthick)) /// Zoom 
	   (lfit zoom ist_norm if (inrange(ist_norm, 0,18)), lcolor(blue) lwidth(medthick)) /// 
	   (qfit donut ist_norm if (inrange(ist_norm, -70,-4)), lcolor(gray) lwidth(medthick)) /// Donut 
	   (qfit donut ist_norm if (inrange(ist_norm, 4,50)), lcolor(gray) lwidth(medthick)) /// 
	   , xline(0 , lpattern(dash) lcolor(red)) ///
	   ytitle("Standardized test score") xtitle("IST distance to threshold") ///
	   legend(ring(0) pos(11) row(6) order(1 3 5 7 9 11) label(1 "Local average") ///
	   label(3 "Linear fit") label(5 "Quadratic fit") label(7 "Cubic fit") ///
	   label(9 "Zoom") label(11 "Donut")) name("`i'",replace) nodraw
	   drop lin quad cub zoom donut m1 m2
	   }   


graph combine math language other , ycommon iscale(0.5) row(1)

