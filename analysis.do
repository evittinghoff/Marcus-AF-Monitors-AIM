cap program drop analysis
program define analysis
syntax [, table1 srlag bp_stats bp_km bp_cco scram_cco peth plotdata reviews ]
	
	if "`table1'"~="" {
		use ./dta/baseline, clear
		keep id ptid srsens age male raceth ///
			smokstat htn dm cad hf bmi ///
			amd drd prp dsp flc qnd stl dft prc 
		qui save ./dta/table1, replace

		qui use ./dta/af_episodes, clear
		qui collapse days (count) af_n=onset, by(id device)
		qui collapse (sum) days af_n, by(id)
		qui merge 1:1 id using ./dta/table1, nogen
		qui replace af_n = 0 if missing(af_n)
		qui recode af_n 2/max=1, gen(af_any)
		label var days "total days of device use"
		label var af_n "# AF episodes detected"
		label var af_any "any AF episodes"
		label val af_any yesno
		qui save ./dta/table1, replace
		
		* number of times button pushed
		qui use ./dta/triggers, clear
		egen ntriggers = rownonmiss(dttm*)
		label var ntriggers "# of times button pushed"
		keep id ntriggers
		qui merge 1:1 id using ./dta/table1, nogen
		qui save ./dta/table1, replace
		
		* WRISTAS data
		qui use ./dta/wristas
		collapse (sum) n_events = event, by(id)
		qui merge 1:1 id using ./dta/table1, nogen
		qui gen in_wristas = cond(~missing(n_events),1,0)
		qui gen any_wrev = cond(n_events>0,1,0) if in_wristas
		label var in_wristas "any WrisTAS data"
		label var any_wrev "any WrisTAS-detected drinking event"
		label val in_wristas any_wrev yesno
		qui save ./dta/table1, replace
		
		* SCRAM data
		qui use ./dta/scram, clear
		collapse (count) scev=sttm, by(id)
		qui merge 1:1 id using ./dta/table1, nogen
		qui gen in_scram = cond(~missing(scev),1,0)
		qui gen any_scev = cond(scev>0,1,0) if in_scram
		label var in_scram "any SCRAM data"
		label var any_scev "any SCRAM-detected drinking event"
		label val in_scram any_scev yesno
		qui save ./dta/table1, replace
		
		* PETH data
		qui use ./dta/peth, clear
		collapse peth_amt (max) peth_det, by(id)
		qui merge 1:1 id using ./dta/table1, nogen
		qui gen in_peth = cond(~missing(peth_amt),1,0)
		label var in_peth "any PeTH data"
		label var peth_amt "average of 2- and 4-week PeTH values"
		label var peth_det "any positive PeTH"
		label val in_peth peth_det yesno
		qui save ./dta/table1, replace
		
		di _n as result "coincidentally, column N's match"
		tab srsens af_any
		foreach bv in af_any srsens {
			local tl: var label `bv'
			maketable, by(`bv') total(before) ///
				vars(age conts\age contn\male cat\raceth cat\ ///
					smokstat cat\htn cat\dm cat\cad cat\hf cat\ ///
					amd cat\drd cat\prp cat\dsp cat\flc cat\qnd cat\stl cat\ ///
					dft cat\prc cat\ntriggers conts\days conts\af_any cat\af_n conts\ ///
					in_wristas cat\any_wrev cat\in_scram cat\any_scev cat\ ///
					in_peth cat\peth_amt conts\peth_det cat) ///
				title(Table 1: sample characteristics by `tl')
		}
		
	}

	if "`srlag'"~="" {
		use ./dta/aflags, clear
		qui keep if ~missing(etohaf04) | ~missing(etohaf403)
		qui replace etohaf04 = "0"+etohaf04 if substr(etohaf04,2,1)==":"
		qui replace etohaf403 = "0"+etohaf403 if substr(etohaf403,2,1)==":"
		qui gen lag = real(substr(etohaf04,1,2))+real(substr(etohaf04,4,2))/60 
		qui replace lag = real(substr(etohaf403,1,2))+real(substr(etohaf403,4,2))/60 ///
			if missing(lag) & ~missing(etohaf403)
		label var lag "Q1: self-reported hours from ETOH consumption to AF onset"
		printlabel lag, space
		tabstat lag, by(week) s(n mean p5 p25 p50 p75 p95) format(%8.3g)
	}
	
	if "`bp_stats'"~="" {
		* duration of monitoring
		di _n as result "Q2 statistics"
		qui use ./dta/af_episodes, clear
		qui collapse days, by(id device)
		qui collapse (sum) days, by(id)
		label var days "total days of device use, from AF patch ID data"
		printlabel days, space
		tabstat days, s(n mean sd p50 p25 p75 min max) format(%8.3g)
		
		* number of times button pushed
		qui use ./dta/triggers, clear
		egen ntriggers = rownonmiss(dttm*)
		label var ntriggers "# of times button pushed"
		printlabel ntriggers, space
		tabstat ntriggers, s(n mean sd p50 p25 p75 min max) format(%8.3g)
		keep id dttm*
		qui reshape long dttm, i(id) j(pushno)
		qui gen dt = dofc(dttm)
		keep id dt
		qui duplicates drop
		qui collapse (count) ntdays=dt, by(id)
		label var ntdays"# of days button pushed"
		printlabel ntdays, space
		tabstat ntdays, s(n mean sd p50 p25 p75 min max) format(%8.3g)
		
		* events and triggers
		use ./dta/af_episodes, clear
		qui keep if ~missing(onset)
		qui duplicates drop id onset, force
		preserve
		collapse (count) n_af=onset, by(id)
		label var n_af "# AF episodes"
		printlabel n_af, space
		tabstat n_af, s(n mean sd p50 p25 p75 min max) format(%8.3g)
		restore
		preserve
		qui gen afdt = dofc(onset)
		keep id afdt
		qui duplicates drop
		collapse (count) n_afdays = afdt, by(id)
		label var n_afdays "# days with AF episodes"
		printlabel n_afdays, space
		tabstat n_afdays, s(n mean sd p50 p25 p75 min max) format(%8.3g)
		restore
		printlabel duration, space
		tabstat duration, s(n mean sd p50 p25 p75 min max) format(%8.3g)
	
		* linked triggers and episodes
		qui merge m:1 id using ./dta/triggers, nogen 
		qui drop if missing(onset)
		qui gen minlag = 5000
		qui gen lag = .
		foreach trigger of varlist dttm* {
			qui replace lag = (onset-`trigger')/3600000 ///
				if inrange((onset-`trigger')/3600000,0,minlag)
			qui replace minlag = lag if lag<minlag
			}
		qui recode lag (min/1=1 "<=1 hour") (1/3=2 ">1-3 hours") ///
			(3/6=3 ">3-6 hours") (6/12=4 ">6-12 hours") (12/24=5 ">12-24 hours") ///
			(24/48=6 ">24-48 hours") (.=8 "no earlier push") (48/max=7 ">48"), ///
			gen(lagcat)
		label var lagcat "Hours from most recent trigger push to AF onset"
		tab lagcat
	}
	
	* KM plots for time from button pushes to AF events
	if "`bp_km'"~="" {
		use ./dta/triggers, clear
		keep id ptid dttm*
		qui reshape long dttm, i(id) j(trno)
		qui drop if missing(dttm)
		drop trno
		qui duplicates drop
		qui bysort ptid (dttm): gen trno = _n
		qui tsset ptid trno
		qui gen censlag = (f.dttm-dttm)/3600000
		
		preserve
		use ./dta/af_episodes, clear
		keep id onset
		qui duplicates drop
		bysort id (onset): gen evno = _n
		qui sum evno
		local maxevno = r(max)
		qui reshape wide onset, i(id) j(evno)
		qui save events, replace
		restore
		
		qui merge m:1 id using events, nogen
		erase events.dta
		qui gen evlag = .
		qui gen onset = .
		forvalues i = 1/`maxevno' {
			qui replace onset = onset`i' ///
				if missing(evlag) & dttm<onset`i' & ~missing(onset`i')
			qui replace evlag = (onset`i'-dttm)/3600000 ///
				if missing(evlag) & dttm<onset`i' & ~missing(onset`i')
		}
		qui gen event = cond(evlag<censlag,1,0)
		qui gen lag = min(evlag, censlag)
		qui replace event = 0 if lag>48
		qui recode lag 48/max=48
		qui stset lag,f(event)
		sts graph, failure per(100) ///
			xtitle("{bf:Time from trigger push to AF onset (hours)}") ///
			ytitle("{bf:Cumulative incidence (%)}") ///
			title("") ///
			ylabel(0(5)30, angle(horizontal) format(%8.0f)) ///
			xlabel(0(6)48) ///
			graphregion(color(white))
		qui graph export ./figures/lwzio_km.pdf, replace
		
		qui replace event = cond(evlag<48,1,0)
		qui replace lag = min(evlag, 48)
		qui stset lag,f(event)
		sts graph, failure per(100) ///
			xtitle("{bf:Time from trigger push to AF onset (hours)}") ///
			ytitle("{bf:Cumulative incidence (%)}") ///
			title("") ///
			ylabel(0(5)30, angle(horizontal) format(%8.0f)) ///
			xlabel(0(6)48) ///
			graphregion(color(white))
		qui graph export ./figures/lwzio_km2.pdf, replace
	}
	
	* case crossover using button pushes
	if "`bp_cco'"~="" {
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		format onsetdt %td
		* keep first episode on each date
		qui bysort id onsetdt (onset): keep if _n==1
		qui gen af = cond(~missing(onset),1,0)
		qui merge m:1 id using ./dta/baseline, nogen keepusing(id srsens)
		qui merge m:1 id using ./dta/triggers, nogen
		* rationalize start and stop days
		qui bysort id source: egen fonset = min(onset)
		qui bysort id source: egen lonset = max(onset)
		qui gen fdttm = .
		qui gen ldttm = .
		forvalues i = 1/4 {
			qui replace fdttm = fdttm`i' if source==`i'
			qui replace ldttm = ldttm`i' if source==`i'
		}
		qui gen start = min(dateon,fonset,fdttm)-1
		qui gen stop = max(dateoff,lonset,ldttm)+1
		format start stop fonset lonset fdttm ldttm %tc
		qui gen day = onsetdt-dofc(start)+1
		qui replace days = dofc(stop)-dofc(start)+1
		qui bysort id source: egen fday = min(day)
		qui expand days+1 if day==fday, gen(cntrl)
		qui bysort id source cntrl: replace day = _n if cntrl==1
		qui bysort id source day (cntrl): keep if _n==1
		qui replace af = 0 if cntrl==1
		* impute average onset time on AF days for control days
		qui gen onset_hr = hh(onset)+mm(onset)/60+ss(onset)/3600 if af==1
		qui bysort id: egen m_onset_hr = mean(onset_hr)
		qui replace onset = 3600000*(24*(dofc(start)+day-1)+m_onset_hr) if af==0
		qui gen dow = dow(dofc(onset))
		qui gen week = ceil(day/7)
		qui egen group = group(id source)
		qui gen nt4 = 0
		qui gen nt12 = 0
		qui gen nt24 = 0
		qui gen nt48 = 0
		foreach hh in 4 12 24 48 {
			foreach trigger of varlist dttm* {
				qui replace nt`hh' = nt`hh'+1 if inrange(onset-`trigger',0,`hh'*3600000)
			}
			qui recode nt`hh' (2/max=1 ">=1"), gen(nt`hh'any)
			qui recode nt`hh' (3/max=2 ">=2"), gen(nt`hh'cat)
			qui recode nt`hh' (6/max=5 ">=5"), gen(nt`hh'cat2)
			forvalues i = 0/1 {
				qui gen nt`hh'_`i' = nt`hh'*cond(srsens==`i',1,0)
				qui gen nt`hh'any_`i' = nt`hh'any*cond(srsens==`i',1,0)
				qui gen nt`hh'cat_`i' = nt`hh'cat*cond(srsens==`i',1,0)
				qui gen nt`hh'cat2_`i' = nt`hh'cat2*cond(srsens==`i',1,0)
			}
			label var nt`hh' "# button pushes, last `hh' hours"
			label var nt`hh'any "any button pushes, last `hh' hours"
			label var nt`hh'cat "# button pushes, last `hh' hours"
			label var nt`hh'cat2 "# button pushes, last `hh' hours"
			label var nt`hh'_1 "# button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'any_1 "any button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'cat_1 "# button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'cat2_1 "# button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'_0 "# button pushes, last `hh' hrs, ETOH insensitive"
			label var nt`hh'any_0 "any button pushes, last `hh' hrs, ETOH insensitive"
			label var nt`hh'cat_0 "# button pushes, last `hh' hrs, ETOH insensitive"
			label var nt`hh'cat2_0 "# button pushes, last `hh' hrs, ETOH insensitive"
		}
		label var af "AF episode"
		label var dow "day of week"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		label val *any* yesno
		label var week "follow-up week"
		bvmv af i.nt4cat i.dow, ///
			model(clogit) modelopts(group(group)) id(id)
		bvmv af nt4any_1 nt4any_0 i.dow, ///
			model(clogit) modelopts(group(group)) id(id)
		qui testparm nt4any_*, equal
		pvalue r(p)
		di as result ///
			"  equality of button push effects, by self-reported sensitivity: `r(pvalue)'"
		bvmv af i.nt4cat2 i.dow, ///
			model(clogit) modelopts(group(group)) id(id)
		bvmv af i.nt12cat i.dow, ///
			model(clogit) modelopts(group(group)) id(id)
		bvmv af i.nt24cat i.dow, ///
			model(clogit) modelopts(group(group)) id(id)
		bvmv af i.nt48cat i.dow, ///
			model(clogit) modelopts(group(group)) id(id)
	}
	
	* case crossover using SCRAM data
	if "`scram_cco'"~="" {
		qui use ./dta/scram, clear
		qui bysort id: egen fsdt = min(dt)
		qui bysort id: egen lsdt = max(dt)
		keep id epnum sttm_ebrac pk_ebrac fsdt lsdt
		qui drop if missing(sttm_ebrac)
		qui bysort id: gen nr = _N
		qui sum nr
		local mnr = r(max)
		qui reshape wide sttm_ebrac pk_ebrac, i(id) j(epnum)
		qui save tmp, replace
		
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		* keep first episode on each date
		qui bysort id onsetdt (onset): keep if _n==1
		qui gen af = cond(~missing(onset),1,0)
		qui merge m:1 id using tmp, nogen
		qui drop if ~inrange(onsetdt,fsdt,lsdt)
		* rationalize start and stop days
		qui bysort id source: egen fonset = min(onset)
		qui bysort id source: egen lonset = max(onset)
		qui gen start = min(dateon,fonset)-1
		qui gen stop = max(dateoff,lonset)+1
		qui gen day = onsetdt-dofc(start)+1
		qui replace days = dofc(stop)-dofc(start)+1
		qui bysort id source: egen fday = min(day)
		qui expand days+1 if day==fday, gen(cntrl)
		qui bysort id source cntrl: replace day = _n if cntrl==1
		qui bysort id source day (cntrl): keep if _n==1
		qui replace af = 0 if cntrl==1
		* impute average onset time on AF days for control days
		qui gen onset_hr = hh(onset)+mm(onset)/60+ss(onset)/3600 if af==1
		qui bysort id: egen m_onset_hr = mean(onset_hr)
		qui replace onset = 3600000*(24*(dofc(start)+day-1)+m_onset_hr) if af==0
		qui gen dow = dow(dofc(onset))
		qui gen week = ceil(day/7)
		qui egen group = group(id source)
		qui merge m:1 id using tmp, nogen
		erase tmp.dta
		foreach hh in 12  {
			qui gen mpk_ebrac`hh' = 0
			forvalues i = 1/`mnr' {
				qui replace mpk_ebrac`hh' = max(mpk_ebrac`hh',pk_ebrac`i') ///
					if inrange(onset-sttm_ebrac`i',0,`hh'*3600000)
			}
			qui recode mpk_ebrac`hh' (0=0 "0") (0/.05=1 ">0-0.05") ///
						(.05/.15=2 ">0.05-0.15") (.15/max=3 ">0.15"), gen(mpk_ebrac`hh'cat)
			qui replace mpk_ebrac`hh' = mpk_ebrac`hh'*10
			label var mpk_ebrac`hh' "SCRAM peak eBrAC (per 0.1), last `hh' hours"
			label var mpk_ebrac`hh'cat "SCRAM peak eBrAC, last `hh' hours"
		}
		label var af "AF episode"
		label var dow "day of week"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		label var week "follow-up week"
		foreach hh in 12 {
			foreach tt in _ebrac {
				bvmv af mpk`tt'`hh' i.dow, ///
					model(clogit) modelopts(group(group)) id(id)
				bvmv af i.mpk`tt'`hh'cat i.dow, ///
					model(clogit) modelopts(group(group)) id(id)
			}
		}
	}
	
	if "`peth'"~="" {
		* trigger pushes by date
		use ./dta/triggers, clear
		keep id dttm*
		qui reshape long dttm, i(id) j(epno)
		qui gen onsetdt = dofc(dttm)
		qui collapse (count) bp_n=dttm, by(id onsetdt)
		qui save tmp, replace
		
		* af episodes by week
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		qui replace dateon = dofc(dateon)
		qui bysort id: egen fdateon = min(dateon)
		qui collapse (count) af_n = onset, by(id fdateon onsetdt)
		qui merge 1:1 id onsetdt using tmp, nogen
		qui save tmp, replace
			
		* SCRAM
		qui use ./dta/scram, clear
		qui gen onsetdt = dofc(sttm_ebrac)
		keep id onsetdt auc_ebrac
		collapse (sum) auc_ebrac, by(id onsetdt)
		qui merge 1:1 id onsetdt using tmp, nogen
		
		qui bysort id: egen fonsetdt = min(onsetdt)
		qui gen start = min(fonsetdt, fdateon)
		qui gen day = onsetdt-start+1
		qui recode day 1/14=2 15/28=4 29/max=5, gen(week)
		
		* count # of episode and trigger pushes on each date
		qui collapse (sum) af_n bp_n auc_ebrac, by(id week)
		qui save tmp, replace
	
		qui merge 1:1 id week using ./dta/peth, nogen //keep(1 3)
		erase tmp.dta
		qui keep if inlist(week,2,4) 
		qui replace peth_amt = 0 if peth_det==0 & missing(peth_amt)
		foreach x in af bp {
			qui replace `x'_n = 0 if missing(`x'_n)
		}
		qui recode af_n 2/max=1, gen(af_any)
		qui recode auc_ebrac (0=0 "0") (0/10=1 ">0/10") (10/25=2 ">10-25") ///
			(25/max=3 ">25"), gen(auc_ebraccat)
		label var week "study week"
		label var af_n "# AF episodes in two-week period"
		label var af_any "any AF episodes in two-week period"
		label var bp_n "# button pushes in two-week period"
		label var auc_ebrac "total SCRAM eBrAC AUC in two-week period"
		label var auc_ebraccat "total SCRAM eBrAC AUC in two-week period"
		label define week 2 "1-2" 4 "3-4", replace
		label val week week
		
		maketable, by(peth_det) npv ///
			vars(af_n conts\af_n contn\bp_n conts\bp_n contn\ ///
				auc_ebrac conts\auc_ebrac contn\auc_ebraccat cat) ///
			title(# AF episodes, button pushes, SCRAM AUC by PETH detectability) ///
			subtitle(unit of analysis is two-week period)
			
		bvmv peth_det bp_n auc_ebrac i.auc_ebraccat, ///
			model(logistic) modelopts(vce(cluster ptid)) id(ptid) ///
			bvonly bvn bvadj(i.week) bvtag(Minimally** adjusted) ///
			subtitle(**for study week, using robust SEs)
			
		bvmv peth_det bp_n auc_ebrac i.week, ///
			model(logistic) modelopts(vce(cluster ptid)) id(ptid) ///
			subtitle(using robust SEs)
			
		bvmv peth_det bp_n i.auc_ebraccat i.week, ///
			model(logistic) modelopts(vce(cluster ptid)) id(ptid) ///
			subtitle(using robust SEs)
			
		bvmv af_n peth_det peth_amt auc_ebrac i.auc_ebraccat, ///
			model(nbreg) modelopts(vce(cluster ptid)) id(ptid) ///
			bvonly bvn bvadj(i.week) bvtag(Minimally** adjusted) ///
			subtitle(**for study week, using robust SEs)
			
		* Figure 1 boxplot
		label var bp_n "{bf:Button Presses}"
		label var auc_ebrac "{bf:Blood Alcohol Concentration AUC}"
		qui recode peth_det (1=0 "{bf:PEth Positive}") (0=1 "{bf:PEth Negative}"), ///
			gen(peth_pos)
		graph box bp_n, over(peth_pos, label(labsize(vlarge))) ///
			ytitle(,size(vlarge)) graphregion(color(white)) ///
			saving(box_bpn,replace)
		graph box auc_ebrac, over(peth_pos, label(labsize(vlarge))) ///
			ytitle(,size(vlarge)) graphregion(color(white)) ///
			saving(box_auc,replace)
		graph combine box_bpn.gph box_auc.gph, ///
			rows(1) ysize(4) xsize(7) graphregion(color(white))
		graph export ./figures/figure1.pdf, replace

	}
			
	* plot data for Sean
	if "`plotdata'"~="" {
		use ./dta/triggers, clear
		keep id dttm*
		qui reshape long dttm, i(id) j(trno)
		qui drop if missing(dttm)
		qui duplicates drop id dttm, force
		
		preserve
		use ./dta/af_episodes, clear
		keep id onset
		qui duplicates drop
		bysort id (onset): gen evno = _n
		qui sum evno
		local maxevno = r(max)
		qui reshape wide onset, i(id) j(evno)
		qui save events, replace
		restore
		
		qui merge m:1 id using events, nogen
		qui gen hours = .
		forvalues i = 1/`maxevno' {
			qui replace hours = (onset`i'-dttm)/3600000 ///
				if missing(hours) & ~missing(onset`i') & ///
					inrange((onset`i'-dttm)/3600000,0,12)
		}
		qui keep if ~missing(hours)
		keep id dttm hours
		qui export delimited id dttm hours using ./csv/plotdata_new_1, replace
					
		use ./dta/scram, clear
		keep id sttm pk_ebrac
		
		qui merge m:1 id using events, nogen
		erase events.dta
		qui gen hours = .
		qui gen pkebrac = .
		forvalues i = 1/`maxevno' {
			qui replace hours = (onset`i'-sttm)/3600000 ///
				if missing(hours) & ~missing(onset`i') ///
					& inrange((onset`i'-sttm)/3600000,0,12)
			qui replace pkebrac = pk_ebrac if ~missing(hours)
		}
		qui keep if ~missing(hours)
		keep id sttm hours pkebrac
		qui export delimited id sttm hours pkebrac using ./csv/plotdata_new_2, ///
			replace
	}
	
	if "`reviews'"~="" {
		* trigger pushes by date
		use ./dta/triggers, clear
		keep id dttm*
		qui reshape long dttm, i(id) j(epno)
		qui gen onsetdt = dofc(dttm)
		qui collapse (count) bp_n=dttm, by(id onsetdt)
		qui save tmp, replace
			
		* SCRAM
		qui use ./dta/scram, clear
		qui gen onsetdt = dofc(sttm_ebrac)
		keep id onsetdt auc_ebrac
		collapse (sum) auc_ebrac, by(id onsetdt)
		qui merge 1:1 id onsetdt using tmp, nogen
		
		qui merge m:1 id using ./dta/insample, nogen keep(3)
		
		di _n as result "Correlation of button pushes and SCRAM, by date"
		pwcorr bp_n auc_ebrac
		spearman bp_n auc_ebrac, stats(rho p)
		
		qui bysort id: egen start = min(onsetdt)
		qui gen day = onsetdt-start+1
		qui recode day 1/14=2 15/28=4 29/max=5, gen(week)
		
		* count # of episode and trigger pushes overall
		preserve
		qui collapse (sum) bp_n auc_ebrac, by(id)
		qui save tmp, replace
				
		di _n as result "Correlation of button pushes and SCRAM, overall"
		pwcorr bp_n auc_ebrac, sig
		spearman bp_n auc_ebrac
		restore
		preserve
		qui collapse (sum) bp_n auc_ebrac, by(id week)
		qui save tmp, replace
		
		qui merge 1:1 id week using ./dta/peth, nogen //keep(1 3)
		erase tmp.dta
		qui keep if inlist(week,2,4) 
		qui replace peth_amt = 0 if peth_det==0 & missing(peth_amt)
		foreach x in bp {
			qui replace `x'_n = 0 if missing(`x'_n)
		}
		qui recode auc_ebrac (0=0 "0") (0/10=1 ">0/10") (10/25=2 ">10-25") ///
			(25/max=3 ">25"), gen(auc_ebraccat)
		qui recode peth_amt (0=0 "0") (0/19=1 "Q1 (>0-19)") (19/37=2 "Q2 (>19-37)") ///
			(37/81=3 "Q3 (>37-81)") (81/max=4 "Q4 (>81)"), gen(peth_amtcat)
		qui recode bp_n (0=0 "0") (0/6=1 "Q1 (1-6)") (7/12=2 "Q2 (7-12)") ///
			(8/20=3 "Q3 (8-21)") (21/max=4 "Q4 (>21)"), gen(bp_ncat)
			
		label var week "study week"
		label var bp_n "# button pushes in two-week period"
		label var bp_ncat "# button pushes in two-week period"
		label var auc_ebrac "total SCRAM eBrAC AUC in two-week period"
		label var auc_ebraccat "total SCRAM eBrAC AUC in two-week period"
		label var peth_amtcat "PeTH amount in two-week period"
		label define week 2 "1-2" 4 "3-4", replace
		label val week week
		
		di _n as result "distribution of PETH readings"
		tab peth_det
		di _n as result "PETH amounts for readings marked as detectable"
		tabstat peth_amt if peth_det==1, s(n mean sd min p25 p50 p75 max) ///
			format(%8.2f)

		di _n as result "Correlations of PeTH amounts with SCRAM AUC"
		spearman peth_amt auc_ebrac
		tab peth_amtcat auc_ebraccat, row col nokey 
		di _n as result "Correlations of PeTH amounts with button pushes"
		spearman peth_amt bp_n
		tab peth_amtcat bp_ncat, row col nokey
		restore
		
		preserve
		qui replace bp_n = round(bp_n)
		qui recode bp_n (0=0 "0") (1 2=1 "1-2") (3 4=2 "3-4") (5/max=3 ">=5"), ///
			gen(bp_ncat)
		label var bp_ncat "# button pushes in 24-hour period"
		tab bp_ncat 
		restore
		
		* duration of wearing monitor
		preserve
		use ./dta/coverage, clear
		qui merge 1:1 id using ./dta/insample, nogen keep(3)
		qui gen days1 = max(0,(5-lw201))*14/4
		qui gen days2 = max(0,(5-lw401))*14/4
		qui egen days = rowtotal(days1 days2 zio*days)
		qui recode days (0/7=1 "0-7") (7/14=2 ">7-14") (14/21=3 ">14-21") ///
			(21/max=4 ">21-28"), gen(dayscat)
		label var dayscat "days of monitor use"
		tab dayscat
		restore
		
		* any use of alcohol among 56 with AF
		use ./dta/triggers, clear
		keep id dttm*
		qui merge 1:1 id using ./dta/insample, nogen keep(3)
		qui reshape long dttm, i(id) j(epno)
		qui gen onsetdt = dofc(dttm)
		qui collapse (count) bp_n=dttm, by(id onsetdt)
		qui save tmp, replace
		
		* af episodes by week
		use ./dta/af_episodes, clear
		qui merge m:1 id using ./dta/insample, nogen keep(3)
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		qui replace dateon = dofc(dateon)
		qui bysort id: egen fdateon = min(dateon)
		qui collapse (count) af_n = onset, by(id fdateon onsetdt)
		qui merge 1:1 id onsetdt using tmp, nogen
		erase tmp.dta
		qui collapse (sum) af_n bp_n, by(id)
		di _n as result "# ppts who had any AF and consumed any alcohol"
		count if af_n>0 & bp_n>0 & ~missing(af_n) & ~missing(bp_n)
		foreach x in af bp {
			qui recode `x'_n 2/max=1, gen(`x'_any)
		}
		label define yesno 0 "no" 1 "yes"
		label val *_any yesno
		label var bp_any "any button pushes"
		label var af_any "any AF episodes"
		tab bp_any af_any, exact
		
		* heterogeneity of AF response across ppts
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		format onsetdt %td
		* keep first episode on each date
		qui bysort id onsetdt (onset): keep if _n==1
		qui gen af = cond(~missing(onset),1,0)
		qui merge m:1 id using ./dta/baseline, nogen keepusing(id srsens)
		qui merge m:1 id using ./dta/triggers, nogen
		* rationalize start and stop days
		qui bysort id source: egen fonset = min(onset)
		qui bysort id source: egen lonset = max(onset)
		qui gen fdttm = .
		qui gen ldttm = .
		forvalues i = 1/4 {
			qui replace fdttm = fdttm`i' if source==`i'
			qui replace ldttm = ldttm`i' if source==`i'
		}
		qui gen start = min(dateon,fonset,fdttm)-1
		qui gen stop = max(dateoff,lonset,ldttm)+1
		format start stop fonset lonset fdttm ldttm %tc
		qui gen day = onsetdt-dofc(start)+1
		qui replace days = dofc(stop)-dofc(start)+1
		qui bysort id source: egen fday = min(day)
		qui expand days+1 if day==fday, gen(cntrl)
		qui bysort id source cntrl: replace day = _n if cntrl==1
		qui bysort id source day (cntrl): keep if _n==1
		qui replace af = 0 if cntrl==1
		qui replace onsetdt = dofc(start)+day-1
		qui drop if missing(day)
		* impute average onset time on AF days for control days
		qui gen onset_hr = hh(onset)+mm(onset)/60+ss(onset)/3600 if af==1
		qui bysort id: egen m_onset_hr = mean(onset_hr)
		qui replace onset = 3600000*(24*(dofc(start)+day-1)+m_onset_hr) if af==0
		qui gen dow = dow(onsetdt)
		qui gen week = ceil(day/7)
		qui egen group = group(id source)
		qui bysort group: egen in_cco = max(af)
		qui gen nt4 = 0
		foreach hh in 4 {
			foreach trigger of varlist dttm* {
				qui replace nt`hh' = nt`hh'+1 if inrange(onset-`trigger',0,`hh'*3600000)
			}
			qui recode nt`hh' (2/max=1 ">=1"), gen(nt`hh'any)
			qui recode nt`hh' (3/max=2 ">=2"), gen(nt`hh'cat)
			qui recode nt`hh' (6/max=5 ">=5"), gen(nt`hh'cat2)
			forvalues i = 0/1 {
				qui gen nt`hh'_`i' = nt`hh'*cond(srsens==`i',1,0)
				qui gen nt`hh'any_`i' = nt`hh'any*cond(srsens==`i',1,0)
				qui gen nt`hh'cat_`i' = nt`hh'cat*cond(srsens==`i',1,0)
				qui gen nt`hh'cat2_`i' = nt`hh'cat2*cond(srsens==`i',1,0)
			}
			label var nt`hh' "# button pushes, last `hh' hours"
			label var nt`hh'any "any button pushes, last `hh' hours"
			label var nt`hh'cat "# button pushes, last `hh' hours"
			label var nt`hh'cat2 "# button pushes, last `hh' hours"
			label var nt`hh'_1 "# button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'any_1 "any button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'cat_1 "# button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'cat2_1 "# button pushes, last `hh' hrs, ETOH sensitive"
			label var nt`hh'_0 "# button pushes, last `hh' hrs, ETOH insensitive"
			label var nt`hh'any_0 "any button pushes, last `hh' hrs, ETOH insensitive"
			label var nt`hh'cat_0 "# button pushes, last `hh' hrs, ETOH insensitive"
			label var nt`hh'cat2_0 "# button pushes, last `hh' hrs, ETOH insensitive"
		}
		label var af "AF episode"
		label var dow "day of week"
		label var in_cco "in case-crossover sample"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		label val *any* yesno
		label var week "follow-up week"
		bvmv af nt4any_1 nt4any_0 i.dow, ///
			model(clogit) modelopts(group(group)) id(id) prntonly(nt4any_*)
		qui testparm nt4any_*, equal
		pvalue r(p)
		di as result ///
			"  equality of button push effects: `r(pvalue)'" _n(2) ///
			"ICC of for participant effects, melogit model, all ppts"
		bvmv af i.nt4any i.dow, ///
			model(melogit) modelopts(|| id:, or) id(id) prntonly(i.nt4any)
		estat icc
		di _n as result "ICC of for participant effects, melogit model, CCO ppts"
		bvmv af i.nt4any i.dow if in_cco==1, ///
			model(melogit) modelopts(|| id:, or) id(id) prntonly(i.nt4any)
		estat icc
				
		* trigger pushes by date
		use ./dta/triggers, clear
		keep id dttm*
		qui reshape long dttm, i(id) j(epno)
		qui gen onsetdt = dofc(dttm)
		qui collapse (count) bp_n=dttm, by(id onsetdt)
		di _n as result "drinking episodes per day, all days"
		tabstat bp_n, s(n mean sd p25 p50 p75 max) format(%8.3g)
		di _n as result "drinking episodes per day, days with any drinking"
		tabstat bp_n if bp_n>0, s(n mean sd p25 p50 p75 max) format(%8.3g)
		
		* interaction of lag and number of button pushes
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		format onsetdt %td
		* keep first episode on each date
		qui bysort id onsetdt (onset): keep if _n==1
		qui gen af = cond(~missing(onset),1,0)
		qui merge m:1 id using ./dta/baseline, nogen keepusing(id srsens)
		qui merge m:1 id using ./dta/triggers, nogen
		* rationalize start and stop days
		qui bysort id source: egen fonset = min(onset)
		qui bysort id source: egen lonset = max(onset)
		qui gen fdttm = .
		qui gen ldttm = .
		forvalues i = 1/4 {
			qui replace fdttm = fdttm`i' if source==`i'
			qui replace ldttm = ldttm`i' if source==`i'
		}
		qui gen start = min(dateon,fonset,fdttm)-1
		qui gen stop = max(dateoff,lonset,ldttm)+1
		format start stop fonset lonset fdttm ldttm %tc
		qui gen day = onsetdt-dofc(start)+1
		qui replace days = dofc(stop)-dofc(start)+1
		qui bysort id source: egen fday = min(day)
		qui expand days+1 if day==fday, gen(cntrl)
		qui bysort id source cntrl: replace day = _n if cntrl==1
		qui bysort id source day (cntrl): keep if _n==1
		qui replace af = 0 if cntrl==1
		qui replace onsetdt = dofc(start)+day-1
		qui drop if missing(day)
		* impute average onset time on AF days for control days
		qui gen onset_hr = hh(onset)+mm(onset)/60+ss(onset)/3600 if af==1
		qui bysort id: egen m_onset_hr = mean(onset_hr)
		qui replace onset = 3600000*(24*(dofc(start)+day-1)+m_onset_hr) if af==0
		qui gen dow = dow(onsetdt)
		qui gen week = ceil(day/7)
		qui egen group = group(id source)
		qui bysort group: egen in_cco = max(af)
		qui gen nt04 = 0
		qui gen nt48 = 0
		qui gen nt812 = 0
		foreach s in 0 4 8 {
			local e = `s'+4
			local hh = "`s'`e'"
			foreach trigger of varlist dttm* {
				qui replace nt`hh' = nt`hh'+1 ///
					if inrange(onset-`trigger',`s'*3600000,`e'*3600000)
			}
			qui recode nt`hh' (2/max=1 ">=1"), gen(nt`hh'any)
			qui recode nt`hh' (3/max=2 ">=2"), gen(nt`hh'cat)
			label var nt`hh' "# button pushes, last `s'-`e' hours"
			label var nt`hh'any "any button pushes, last `s'-`e' hours"
			label var nt`hh'cat "# button pushes, last `s'-`e' hours"
		}
		label var in_cco "in case-crossover sample"
		label var af "AF episode"
		label var dow "day of week"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		label val *any* yesno
		label var week "follow-up week"
		bvmv af nt04 nt48 nt812, ///
			model(clogit) modelopts(group(group)) id(id) prntonly(nt*)
		qui testparm nt*, equal
		pvalue r(p)
		di as result ///
			"  equality of button push effects: `r(pvalue)'" 
			
		* button pushes by dow and time of day
		use ./dta/triggers, clear
		keep id dttm* source* days*
		qui reshape long dttm source days, i(id) j(epno)
		rename dttm onset
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		format onsetdt %td
		qui gen dow = dow(onsetdt)
		label var dow "day of week"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		qui gen time = ceil((hh(onset)+1)/4)
		label define time 1 "12am-4am" 2 "4am-8am" 3 "8am-12pm" 4 "12pm-4pm" ///
			5 "4pm-8pm" 6 "8pm-12am"
		label val time time
		label var time "time of day"
		di _n as result "overall distribution of button push times"
		tab dow time, row col nokey
		
		* individual distributions
		qui collapse (count) bp_n=onset, by(id onsetdt dow time source days)
		format bp_n %8.0f
		label var bp_n "# button presses"
		qui bysort id source (onsetdt): gen recno = _n
		qui expand 6*(days+1) if recno==1, gen(extra)
		qui bysort id source onsetdt extra: ///
			replace onsetdt = onsetdt+int(_n/6) if extra==1
		qui bysort id source onsetdt extra: replace time = _n if extra==1
		qui bysort id source onsetdt time (extra): keep if _n==1
		qui replace bp_n = 0 if extra==1
		qui replace dow = dow(onsetdt) if extra==1
		di _n as result "participant #s of button presses by time of day"
		tabstat bp_n, by(time) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
		di _n as result "participant #s of button presses in 4-hour windows, by day of week"
		tabstat bp_n, by(dow) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
		qui collapse (sum) bp_n, by(id onsetdt dow)
		di _n as result "participant #s of button presses in 24-hour windows, by day of week"
		tabstat bp_n, by(dow) s(mean sd p25 p50 p75 p95 max) format(%8.3g)

		* scram events by dow and time of day	
		use ./dta/scram, clear
		keep id sttm
		rename sttm onset
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		format onsetdt %td
		qui gen dow = dow(onsetdt)
		label var dow "day of week"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		qui gen time = ceil((hh(onset)+1)/4)
		label define time 1 "12am-4am" 2 "4am-8am" 3 "8am-12pm" 4 "12pm-4pm" ///
			5 "4pm-8pm" 6 "8pm-12am"
		label val time time
		label var time "time of day"
		di _n as result "overall distribution of SCRAM events"
		tab dow time, row col nokey
		
		* individual distributions
		qui gen days = 28
		qui collapse (count) sc_n=onset, by(id onsetdt dow time days)
		format sc_n %8.0f
		label var sc_n "# SCRAM episodes"
		qui bysort id (onsetdt): gen recno = _n
		qui expand 6*(days+1) if recno==1, gen(extra)
		qui bysort id onsetdt extra: ///
			replace onsetdt = onsetdt+int(_n/6) if extra==1
		qui bysort id onsetdt extra: replace time = _n if extra==1
		drop if time>6
		qui bysort id onsetdt time (extra): keep if _n==1
		qui replace sc_n = 0 if extra==1
		qui replace dow = dow(onsetdt) if extra==1
		di _n as result "participant #s of SCRAM episodes by time of day"
		tabstat sc_n, by(time) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
		di _n as result "participant #s of SCRAM episodes in 4-hour windows, by day of week"
		tabstat sc_n, by(dow) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
		qui collapse (sum) sc_n, by(id onsetdt dow)
		di _n as result "participant #s of SCRAM episodes in 24-hour windows, by day of week"
		tabstat sc_n, by(dow) s(mean sd p25 p50 p75 p95 max) format(%8.3g)

		* AF episodes by dow and time of day	
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui gen onsetdt = dofc(onset)
		format onsetdt %td
		qui gen dow = dow(onsetdt)
		label var dow "day of week"
		label define dow 0 "Sunday" 1 "Monday" 2 "Tuesday" 3 "Wednesday" ///
			4 "Thursday" 5 "Friday" 6 "Saturday"
		label val dow dow
		qui gen time = ceil((hh(onset)+1)/4)
		label define time 1 "12am-4am" 2 "4am-8am" 3 "8am-12pm" 4 "12pm-4pm" ///
			5 "4pm-8pm" 6 "8pm-12am"
		label val time time
		label var time "time of day"
		di _n as result "overall distribution of AF episode times"
		qui keep if ~missing(dow) & ~missing(time)
		tab dow time, row col nokey
		
		* individual distributions
		qui collapse (count) af_n=onset, by(id onsetdt dow time source days)
		format af_n %8.0f
		label var af_n "# AF episodes"
		qui bysort id source (onsetdt): gen recno = _n
		qui expand 6*(days+1) if recno==1, gen(extra)
		qui bysort id source onsetdt extra: ///
			replace onsetdt = onsetdt+int(_n/6) if extra==1
		qui bysort id source onsetdt extra: replace time = _n if extra==1
		qui bysort id source onsetdt time (extra): keep if _n==1
		qui replace af_n = 0 if extra==1
		qui replace dow = dow(onsetdt) if extra==1
		di _n as result "participant #s of AF episodes by time of day"
		tabstat af_n, by(time) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
		di _n as result "participant #s of AF episodes in 4-hour windows, by day of week"
		tabstat af_n, by(dow) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
		qui collapse (sum) af_n, by(id onsetdt dow)
		di _n as result "participant #s of AF episodes in 24-hour windows, by day of week"
		tabstat af_n, by(dow) s(mean sd p25 p50 p75 p95 max) format(%8.3g)
	}
	
end

cd "~/work/parnassus/marcus/projects/monitors/programs"
cd ./dta/
!gunzip -f *.dta.gz
cd ../

set tracedepth 1
set trace off
cap log close
log using analysis, replace
analysis, reviews //table1 scram_cco peth
log close

cd ./dta
!gzip -f *.dta
cd ../
