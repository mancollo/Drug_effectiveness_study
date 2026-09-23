*** ------------------------------------------------------------------------- ***
*** STEP 1: CREATE RETENTION INDICATORS                                       ***
*** ------------------------------------------------------------------------- ***
bysort unique_id: egen present_wave2 = max(cond(Surveycode == 2, 1, 0))
bysort unique_id: egen present_wave3 = max(cond(Surveycode == 3, 1, 0))
bysort unique_id: egen present_wave4 = max(cond(Surveycode == 4, 1, 0))

*** ------------------------------------------------------------------------- ***
*** STEP 2: CONVERT STRINGS & CATEGORIZE AGE                                 ***
*** ------------------------------------------------------------------------- ***
capture encode County, gen(county_num)
capture encode sex, gen(sex_num)
capture encode School_code, gen(school_num)

gen age_cat = .
replace age_cat = 1 if age >= 6  & age <= 10
replace age_cat = 2 if age >= 11 & age <= 14
replace age_cat = 3 if age >= 15 & age <= 18
label define age_lbl 1 "6-10 years" 2 "11-14 years" 3 "15-18 years", replace
label values age_cat age_lbl

*** ------------------------------------------------------------------------- ***
*** STEP 3: ESTIMATE RETENTION PROBABILITIES & IPW                           ***
*** ------------------------------------------------------------------------- ***
logit present_wave2 i.county_num i.age_cat i.sex_num i.Hkinfect i.Asinfect i.Ttinfect if Surveycode == 1
predict p_ret2 if Surveycode == 1

logit present_wave3 i.county_num i.age_cat i.sex_num i.Hkinfect i.Asinfect i.Ttinfect if Surveycode == 1
predict p_ret3 if Surveycode == 1

logit present_wave4 i.county_num i.age_cat i.sex_num i.Hkinfect i.Asinfect i.Ttinfect if Surveycode == 1
predict p_ret4 if Surveycode == 1

bysort unique_id: egen p_wave2 = max(p_ret2)
bysort unique_id: egen p_wave3 = max(p_ret3)
bysort unique_id: egen p_wave4 = max(p_ret4)

gen ipw = 1 if Surveycode == 1
replace ipw = 1 / p_wave2 if Surveycode == 2
replace ipw = 1 / p_wave3 if Surveycode == 3
replace ipw = 1 / p_wave4 if Surveycode == 4

*** ========================================================================= ***
*** MAIN STATISTICAL ANALYSIS (VCE(ROBUST) & ERROR-HANDLED)        ***
*** ========================================================================= ***

*** ========================================================================= ***
*** 1.0 PREVALENCE AT EACH WAVE BY SUBGROUPS                                  ***
*** ========================================================================= ***

*** ========================================================================= ***
*** 1.1 BASELINE (Surveycode == 1)- Table 3                                   ***
*** ========================================================================= ***

*** --- HOOKWORM (Hkinfect) --- ***
bysort county_num: tab Hkinfect if Surveycode == 1
bysort school_num: tab Hkinfect if Surveycode == 1
bysort age_cat:    tab Hkinfect if Surveycode == 1
bysort sex_num:    tab Hkinfect if Surveycode == 1
bysort county_num: binreg Hkinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort school_num: binreg Hkinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort age_cat:    binreg Hkinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort sex_num:    binreg Hkinfect [pw = ipw] if Surveycode == 1, rr vce(robust)

*** --- ASCARIS LUMBRICOIDES (Asinfect) --- ***
bysort county_num: tab Asinfect if Surveycode == 1
bysort school_num: tab Asinfect if Surveycode == 1
bysort age_cat:    tab Asinfect if Surveycode == 1
bysort sex_num:    tab Asinfect if Surveycode == 1
bysort county_num: poisson Asinfect [pw = ipw] if Surveycode == 1, irr vce(robust) // Use poisson to handle high prevalence and weights
bysort school_num: poisson Asinfect [pw = ipw] if Surveycode == 1, irr vce(robust) // Use poisson to handle high prevalence and weights
bysort age_cat:    binreg Asinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort sex_num:    binreg Asinfect [pw = ipw] if Surveycode == 1, rr vce(robust)

*** --- TRICHURIS TRICHIURA (Ttinfect) --- ***
bysort county_num: tab Ttinfect if Surveycode == 1
bysort school_num: tab Ttinfect if Surveycode == 1
bysort age_cat:    tab Ttinfect if Surveycode == 1
bysort sex_num:    tab Ttinfect if Surveycode == 1
bysort county_num: binreg Ttinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort school_num: binreg Ttinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort age_cat:    binreg Ttinfect [pw = ipw] if Surveycode == 1, rr vce(robust)
bysort sex_num:    binreg Ttinfect [pw = ipw] if Surveycode == 1, rr vce(robust)


*** ========================================================================= ***
*** 1.2 FOLLOW-UP 1 (Surveycode == 2) - Table 4                               ***
*** ========================================================================= ***

*** --- HOOKWORM (Hkinfect) --- ***
bysort county_num: tab Hkinfect if Surveycode == 2
bysort school_num: tab Hkinfect if Surveycode == 2
bysort age_cat:    tab Hkinfect if Surveycode == 2
bysort sex_num:    tab Hkinfect if Surveycode == 2
bysort county_num: binreg Hkinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort school_num: binreg Hkinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort age_cat:    binreg Hkinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort sex_num:    binreg Hkinfect [pw = ipw] if Surveycode == 2, rr vce(robust)

*** --- ASCARIS LUMBRICOIDES (Asinfect) --- ***
bysort county_num: tab Asinfect if Surveycode == 2
bysort school_num: tab Asinfect if Surveycode == 2
bysort age_cat:    tab Asinfect if Surveycode == 2
bysort sex_num:    tab Asinfect if Surveycode == 2
bysort county_num: binreg Asinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort school_num: binreg Asinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort age_cat:    binreg Asinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort sex_num:    binreg Asinfect [pw = ipw] if Surveycode == 2, rr vce(robust)

*** --- TRICHURIS TRICHIURA (Ttinfect) --- ***
bysort county_num: tab Ttinfect if Surveycode == 2
bysort school_num: tab Ttinfect if Surveycode == 2
bysort age_cat:    tab Ttinfect if Surveycode == 2
bysort sex_num:    tab Ttinfect if Surveycode == 2
bysort county_num: binreg Ttinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort school_num: binreg Ttinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort age_cat:    binreg Ttinfect [pw = ipw] if Surveycode == 2, rr vce(robust)
bysort sex_num:    binreg Ttinfect [pw = ipw] if Surveycode == 2, rr vce(robust)


*** ========================================================================= ***
*** 1.3 FOLLOW-UP 2 (Surveycode == 3)  - Table 5                              ***
*** ========================================================================= ***

*** --- HOOKWORM (Hkinfect) --- ***
bysort county_num: tab Hkinfect if Surveycode == 3
bysort school_num: tab Hkinfect if Surveycode == 3
bysort age_cat:    tab Hkinfect if Surveycode == 3
bysort sex_num:    tab Hkinfect if Surveycode == 3
bysort county_num: binreg Hkinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort school_num: binreg Hkinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort age_cat:    binreg Hkinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort sex_num:    binreg Hkinfect [pw = ipw] if Surveycode == 3, rr vce(robust)

*** --- ASCARIS LUMBRICOIDES (Asinfect) --- ***
bysort county_num: tab Asinfect if Surveycode == 3
bysort school_num: tab Asinfect if Surveycode == 3
bysort age_cat:    tab Asinfect if Surveycode == 3
bysort sex_num:    tab Asinfect if Surveycode == 3
bysort county_num: binreg Asinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort school_num: binreg Asinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort age_cat:    binreg Asinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort sex_num:    binreg Asinfect [pw = ipw] if Surveycode == 3, rr vce(robust)

*** --- TRICHURIS TRICHIURA (Ttinfect) --- ***
bysort county_num: tab Ttinfect if Surveycode == 3
bysort school_num: tab Ttinfect if Surveycode == 3
bysort age_cat:    tab Ttinfect if Surveycode == 3
bysort sex_num:    tab Ttinfect if Surveycode == 3
bysort county_num: binreg Ttinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort school_num: binreg Ttinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort age_cat:    binreg Ttinfect [pw = ipw] if Surveycode == 3, rr vce(robust)
bysort sex_num:    binreg Ttinfect [pw = ipw] if Surveycode == 3, rr vce(robust)


*** ========================================================================= ***
*** 1.4 FOLLOW-UP 3 (Surveycode == 4) - Table 6                               ***
*** ========================================================================= ***

*** --- HOOKWORM (Hkinfect) --- ***
bysort county_num: tab Hkinfect if Surveycode == 4
bysort school_num: tab Hkinfect if Surveycode == 4
bysort age_cat:    tab Hkinfect if Surveycode == 4
bysort sex_num:    tab Hkinfect if Surveycode == 4
bysort county_num: binreg Hkinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort school_num: binreg Hkinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort age_cat:    binreg Hkinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort sex_num:    binreg Hkinfect [pw = ipw] if Surveycode == 4, rr vce(robust)

*** --- ASCARIS LUMBRICOIDES (Asinfect) --- ***
bysort county_num: tab Asinfect if Surveycode == 4
bysort school_num: tab Asinfect if Surveycode == 4
bysort age_cat:    tab Asinfect if Surveycode == 4
bysort sex_num:    tab Asinfect if Surveycode == 4
bysort county_num: binreg Asinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort school_num: binreg Asinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort age_cat:    binreg Asinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort sex_num:    binreg Asinfect [pw = ipw] if Surveycode == 4, rr vce(robust)

*** --- TRICHURIS TRICHIURA (Ttinfect) --- ***
bysort county_num: tab Ttinfect if Surveycode == 4
bysort school_num: tab Ttinfect if Surveycode == 4
bysort age_cat:    tab Ttinfect if Surveycode == 4
bysort sex_num:    tab Ttinfect if Surveycode == 4
bysort county_num: binreg Ttinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort school_num: binreg Ttinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort age_cat:    binreg Ttinfect [pw = ipw] if Surveycode == 4, rr vce(robust)
bysort sex_num:    binreg Ttinfect [pw = ipw] if Surveycode == 4, rr vce(robust)

*** ========================================================================= ***
*** 2. OVERALL PREVALENCE PER WAVE                                            ***
*** ========================================================================= ***
bysort Surveycode: binreg Sthinfect [pw = ipw], rr vce(robust)
bysort Surveycode: binreg Hkinfect [pw = ipw], rr vce(robust)
bysort Surveycode: binreg Asinfect [pw = ipw], rr vce(robust)
bysort Surveycode: binreg Ttinfect [pw = ipw], rr vce(robust)

*** ========================================================================= ***
*** 3. RELATIVE REDUCTION IN PREVALENCE (VS BASELINE)                         ***
*** ========================================================================= ***

* --- Any STH (Sthinfect) --- *
poisson Sthinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (RelRed_FU1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Sthinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (RelRed_FU2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Sthinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (RelRed_FU3: (1 - exp(_b[4.Surveycode])) * 100)


* --- Hookworm (Hkinfect) --- *
poisson Hkinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (RelRed_FU1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (RelRed_FU2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (RelRed_FU3: (1 - exp(_b[4.Surveycode])) * 100)


* --- Ascaris Lumbricoides (Asinfect) --- *
poisson Asinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (RelRed_FU1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (RelRed_FU2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (RelRed_FU3: (1 - exp(_b[4.Surveycode])) * 100)


* --- Trichuris Trichiura (Ttinfect) --- *
poisson Ttinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (RelRed_FU1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (RelRed_FU2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (RelRed_FU3: (1 - exp(_b[4.Surveycode])) * 100)

*** ========================================================================= ***
*** 4. MEAN INTENSITY OF STH INFECTIONS                                       ***
*** ========================================================================= ***
bysort Surveycode: nbreg hkepg [pw = ipw], irr vce(robust)
bysort Surveycode: nbreg asepg [pw = ipw], irr vce(robust)
bysort Surveycode: nbreg ttepg [pw = ipw], irr vce(robust)


*** ========================================================================= ***
*** 5. RELATIVE REDUCTION IN MEAN INTENSITY (OVERALL)                         ***
*** ========================================================================= ***

* Hookworm (hkepg)
nbreg hkepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_FU1: (1 - exp(_b[2.Surveycode])) * 100)

nbreg hkepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_FU2: (1 - exp(_b[3.Surveycode])) * 100)

nbreg hkepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_FU3: (1 - exp(_b[4.Surveycode])) * 100)

* Ascaris Lumbricoides (asepg)
nbreg asepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_FU1: (1 - exp(_b[2.Surveycode])) * 100)

nbreg asepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_FU2: (1 - exp(_b[3.Surveycode])) * 100)

nbreg asepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_FU3: (1 - exp(_b[4.Surveycode])) * 100)

* Trichuris Trichiura (ttepg)
nbreg ttepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_FU1: (1 - exp(_b[2.Surveycode])) * 100)

nbreg ttepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_FU2: (1 - exp(_b[3.Surveycode])) * 100)

nbreg ttepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_FU3: (1 - exp(_b[4.Surveycode])) * 100)

*-------------------------------------------------------------------------------
* EGG REDUCTION RATES (ERR) WITH 95% CIs
* Subgroups: Overall, County, School, Age Category, Sex
*-------------------------------------------------------------------------------

*===============================================================================
* 1. HOOKWORM
*===============================================================================

* --- HOOKWORM: FOLLOW-UP 1 ---
* Overall
poisson hkepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (ERR_Overall: (1 - exp(_b[2.Surveycode])) * 100)

* County
quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[2.Surveycode])) * 100)

* School
quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

* Age Category
quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust) iter(100)
nlcom (ERR_Age2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[2.Surveycode])) * 100)


* Sex
quietly nbreg hkepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[2.Surveycode])) * 100)


* --- HOOKWORM: FOLLOW-UP 2 ---
* Overall
quietly nbreg hkepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[3.Surveycode])) * 100)

* County
quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[3.Surveycode])) * 100)

* School
quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[3.Surveycode])) * 100)

* Age Category
quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust) iter(100)
nlcom (ERR_Age2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[3.Surveycode])) * 100)

* Sex
quietly nbreg hkepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[3.Surveycode])) * 100)


* --- HOOKWORM: FOLLOW-UP 3 ---
* Overall
quietly nbreg hkepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[4.Surveycode])) * 100)

* County
quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[4.Surveycode])) * 100)

* School
quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[4.Surveycode])) * 100)

* Age Category
quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust) iter(100)
nlcom (ERR_Age2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[4.Surveycode])) * 100)


* Sex
quietly nbreg hkepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg hkepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[4.Surveycode])) * 100)


*===============================================================================
* 2. ASCARIS LUMBRICOIDES
*===============================================================================

* --- ASCARIS LUMBRICOIDES: FOLLOW-UP 1 ---
* Overall
quietly nbreg asepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[2.Surveycode])) * 100)

* County
quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[2.Surveycode])) * 100)

* School
quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[2.Surveycode])) * 100)

* Age Category
quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 4 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age4: (1 - exp(_b[2.Surveycode])) * 100)

* Sex
quietly nbreg asepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[2.Surveycode])) * 100)


* --- ASCARIS LUMBRICOIDES: FOLLOW-UP 2 ---
* Overall
quietly nbreg asepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[3.Surveycode])) * 100)

* County
quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[3.Surveycode])) * 100)

* School
quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[3.Surveycode])) * 100)

* Age Category
quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust) iter(100)
nlcom (ERR_Age3: (1 - exp(_b[3.Surveycode])) * 100)


* Sex
quietly nbreg asepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[3.Surveycode])) * 100)


* --- ASCARIS LUMBRICOIDES: FOLLOW-UP 3 ---
* Overall
quietly nbreg asepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[4.Surveycode])) * 100)

* County
quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[4.Surveycode])) * 100)

* School
quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[4.Surveycode])) * 100)

* Age Category
quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[4.Surveycode])) * 100)


* Sex
quietly nbreg asepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg asepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[4.Surveycode])) * 100)


*===============================================================================
* 3. TRICHURIS TRICHIURA
*===============================================================================

* --- TRICHURIS TRICHIURA: FOLLOW-UP 1 ---
* Overall
quietly nbreg ttepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[2.Surveycode])) * 100)

* County
quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[2.Surveycode])) * 100)

* School
quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[2.Surveycode])) * 100)

* Age Category
quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age2: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 4 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Age4: (1 - exp(_b[2.Surveycode])) * 100)

* Sex
quietly nbreg ttepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[2.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 2), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[2.Surveycode])) * 100)


* --- TRICHURIS TRICHIURA: FOLLOW-UP 2 ---
* Overall
quietly nbreg ttepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[3.Surveycode])) * 100)

* County
quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[3.Surveycode])) * 100)

* School
quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[3.Surveycode])) * 100)

* Age Category
quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age2: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 3), irr vce(robust) iter(100)
nlcom (ERR_Age3: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 4 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Age4: (1 - exp(_b[3.Surveycode])) * 100)

* Sex
quietly nbreg ttepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[3.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 3), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[3.Surveycode])) * 100)


* --- TRICHURIS TRICHIURA: FOLLOW-UP 3 ---
* Overall
quietly nbreg ttepg i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Overall: (1 - exp(_b[4.Surveycode])) * 100)

* County
quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_C3: (1 - exp(_b[4.Surveycode])) * 100)

* School
quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch3: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch4: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch5: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch6: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch7: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Sch8: (1 - exp(_b[4.Surveycode])) * 100)

* Age Category
quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age1: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age2: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Age3: (1 - exp(_b[4.Surveycode])) * 100)

* Sex
quietly nbreg ttepg i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Male: (1 - exp(_b[4.Surveycode])) * 100)

quietly nbreg ttepg i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 4), irr vce(robust)
nlcom (ERR_Female: (1 - exp(_b[4.Surveycode])) * 100)

*** ========================================================================= ***
*** CURE RATES (CR) WITH 95% CONFIDENCE INTERVALS    ***
*** ========================================================================= ***

*** ========================================================================= ***
*** 1. HOOKWORM (Hkinfect)                                                    ***
*** ========================================================================= ***

* --- HOOKWORM: FOLLOW-UP 1 ---
* Overall
poisson Hkinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[2.Surveycode])) * 100)

* County
poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[2.Surveycode])) * 100)

* School
poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[2.Surveycode])) * 100)

* Age Category
poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[2.Surveycode])) * 100)

* Sex
poisson Hkinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[2.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[2.Surveycode])) * 100)


* --- HOOKWORM: FOLLOW-UP 2 ---
* Overall
poisson Hkinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[3.Surveycode])) * 100)

* County
poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[3.Surveycode])) * 100)

* School
poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[3.Surveycode])) * 100)

* Age Category
poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[3.Surveycode])) * 100)

* Sex
poisson Hkinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[3.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[3.Surveycode])) * 100)


* --- HOOKWORM: FOLLOW-UP 3 ---
* Overall
poisson Hkinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[4.Surveycode])) * 100)

* County
poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[4.Surveycode])) * 100)

* School
poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[4.Surveycode])) * 100)

* Age Category
poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[4.Surveycode])) * 100)

* Sex
poisson Hkinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[4.Surveycode])) * 100)

poisson Hkinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[4.Surveycode])) * 100)

*** ========================================================================= ***
*** 2. ASCARIS LUMBRICOIDES (Asinfect)                                        ***
*** ========================================================================= ***

* --- ASCARIS: FOLLOW-UP 1 ---
* Overall
poisson Asinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[2.Surveycode])) * 100)

* County
poisson Asinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[2.Surveycode])) * 100)

* School
poisson Asinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[2.Surveycode])) * 100)

* Age Category
poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[2.Surveycode])) * 100)

* Sex
poisson Asinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[2.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[2.Surveycode])) * 100)


* --- ASCARIS: FOLLOW-UP 2 ---
* Overall
poisson Asinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[3.Surveycode])) * 100)

* County
poisson Asinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[3.Surveycode])) * 100)

* School
poisson Asinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[3.Surveycode])) * 100)

* Age Category
poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[3.Surveycode])) * 100)

* Sex
poisson Asinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[3.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[3.Surveycode])) * 100)

* --- ASCARIS: FOLLOW-UP 3 ---
* Overall
poisson Asinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[4.Surveycode])) * 100)

* County
poisson Asinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[4.Surveycode])) * 100)


* School
poisson Asinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch3: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch4: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch5: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch6: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 4), vce(robust)irr
nlcom (CR_Sch7: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[4.Surveycode])) * 100)


* Age Category
poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[4.Surveycode])) * 100)


* Sex
poisson Asinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[4.Surveycode])) * 100)

poisson Asinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[4.Surveycode])) * 100)

*** ========================================================================= ***
*** 3. TRICHURIS TRICHIURA (Ttinfect)                                       ***
*** ========================================================================= ***

* --- TRICHURIS TRICHIURA: FOLLOW-UP 1 ---
* Overall
poisson Ttinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[2.Surveycode])) * 100)

* County
poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[2.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[2.Surveycode])) * 100)

* School
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[2.Surveycode])) * 100)

* Age Category
poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[2.Surveycode])) * 100)
poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[2.Surveycode])) * 100)

* Sex
poisson Ttinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[2.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 2), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[2.Surveycode])) * 100)

* --- TRICHURIS TRICHIURA: FOLLOW-UP 2 ---
* Overall
poisson Ttinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[3.Surveycode])) * 100)

* County
poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[3.Surveycode])) * 100)

* School
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[3.Surveycode])) * 100)

* Age Category
poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[3.Surveycode])) * 100)

* Sex
poisson Ttinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[3.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 3), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[3.Surveycode])) * 100)


* --- TRICHURIS TRICHIURA: FOLLOW-UP 3 ---
* Overall
poisson Ttinfect i.Surveycode [pw = ipw] if (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Overall: (1 - exp(_b[4.Surveycode])) * 100)

* County
poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if county_num == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_C3: (1 - exp(_b[4.Surveycode])) * 100)

* School
poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch3: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 4 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch4: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 5 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch5: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 6 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch6: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 7 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch7: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if school_num == 8 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Sch8: (1 - exp(_b[4.Surveycode])) * 100)

* Age Category
poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age1: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age2: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if age_cat == 3 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Age3: (1 - exp(_b[4.Surveycode])) * 100)

* Sex
poisson Ttinfect i.Surveycode [pw = ipw] if sex_num == 1 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Female: (1 - exp(_b[4.Surveycode])) * 100)

poisson Ttinfect i.Surveycode [pw = ipw] if sex_num == 2 & (Surveycode == 1 | Surveycode == 4), vce(robust) irr
nlcom (CR_Male: (1 - exp(_b[4.Surveycode])) * 100)
