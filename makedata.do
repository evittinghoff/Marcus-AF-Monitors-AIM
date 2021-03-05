
cap program drop date_period
program define date_period
	syntax varlist 
	tempvar hour min sec
	gen `varlist'_dt = dofc(`varlist')
	format `varlist'_dt %td
	gen `hour' = hh(`varlist')
	gen `min' = mm(`varlist')
	gen `sec' = ss(`varlist')
	gen `varlist'_period = ceil((`hour'*3600+`min'*60+`sec')/1800) `if' `in'
	end
	
cap program drop makedata
program define makedata
syntax [, main af_episodes scram wristas peth plotdata ]

	if "`main'"~="" {
		qui import delimited ///
			using "./csv/EtOHInAFEventMonitor_DATA_2020-08-13_1109.csv", ///
			clear

		qui encode redcap, gen(tmp)
		qui recode tmp (1=0 "enrollment") (2=2 "week 2") (3=4 "week 4") ///
			(4=5 "early termination"), gen(week)
		drop redcap tmp
		label var week "study week"
		rename dataanalysisyn insample
		
		drop *notes* *comp* *trigger*rhythm *trigger*other *onset *af*_*
		rename *___* *_*
		rename *trigger* *dttm*
				
		* baseline and birth date
		set trace off
		qui gen dob = date(demo01,"MDY",2020)
		qui gen bldt = date(lifewatch_start,"MDY",2020)
		qui replace bldt = dofc(clock(zio1_start,"MDYhm",2020)) if missing(bldt)
		format dob bldt %td
			
		qui gen male = cond(demo02==1,1,0) if ~missing(demo02)
		label var male "male sex"
		
		* race ethnicity
		qui gen raceth = 1 if demo11_5==1
		qui replace raceth = 2 if demo11_2==1
		qui replace raceth = 3 if demo11_3==1
		qui replace raceth = 4 if demo11_6==1
		qui egen nr = rownonmiss(demo11_?)
		qui replace raceth = 5 if missing(raceth) & nr>0
		drop nr
		label define raceth 1 "white" 2 "Asian" 3 "Black" 4 "Latino" 5 "other"
		label val raceth raceth
		label var raceth "race/ethnicity"
		
		rename physex08 wtkg
		rename physex09 htcm
		rename pmh01 htn
		rename pmh04 cad
		rename pmh06 hf
		rename pmh09 dm
		rename pmh30 bmi
		rename pmh33 smokstat
		rename pmh35 packyrs
		label var htn "known HTN"
		label var cad "known CAD"
		label var hf "known history of CHF"
		label var dm "known diabetes"
		label var bmi "body mass index"
		label var smokstat "smoking history"
		label var packyrs "pack-years of smoking"
		
		rename aad01 bb
		rename aad05 amd
		rename aad08 drd
		rename aad11 prp
		rename aad14 dsp
		rename aad17 flc
		rename aad20 qnd
		rename aad23 stl
		rename aad26 dft
		rename aad29 prc
		label var bb "beta-blockers"
		label var amd "amiodarone"
		label var drd "dronedarone"
		label var prp "propafenone"
		label var dsp "disopiramide"
		label var flc "flecanide"
		label var qnd "quinidine"
		label var stl "sotalol"
		label var dft "dofetilide"
		label var prc "procainimide"
		
		qui gen srsens = cond(etohaf01==1,1,0)
		label var srsens "AF sensitive to drinking by self-report"
		label val srsens srsens
		
		label define yesno 1 "yes" 0 "no"
		label val ///
			male htn cad hf dm bb amd drd prp dsp flc qnd stl dft prc srsens ///
			yesno
		label define smokstat 1 "never" 2 "former" 3 "current"
		label val smokstat smokstat
		
		* put id, baseline date, birthdate on every record
		rename record_id ptid
		foreach x in id bldt dob raceth male smokstat packyrs ///
			htn cad hf dm bmi amd drd prp dsp flc qnd stl dft prc srsens {
			forvalues i = 1/3 {
				qui replace `x' = `x'[_n-`i'] if missing(`x') & ~missing(`x'[_n-`i']) & ///
					ptid==ptid[_n-`i']
				qui replace `x' = `x'[_n+`i'] if missing(`x') & ~missing(`x'[_n+`i']) & ///
					ptid==ptid[_n+`i']
				qui replace insample = insample[_n-`i'] ////
					if missing(insample) & ~missing(insample[_n-`i']) & ///
					ptid==ptid[_n-`i']
			}
		}
		
		qui keep if insample==1
		
		* age 
		qui gen age = int((bldt-dob)/365.25)
		label var age "age (years)"
		
		preserve
		keep id
		qui duplicates drop
		qui save ./dta/insample, replace
		restore
		
		* baseline variables
		preserve
		qui keep if week==0
	
		keep ptid id-tricussurgreplaced heh_yn-scram_serial wristas_device_number ///
			wristas_v1_activation ilr0* srsens age male raceth smokstat packyrs ///
			htn cad hf dm bmi amd drd prp dsp flc qnd stl dft prc
		rename v1_* *
		*des
		*sum
		qui save ./dta/baseline, replace
		restore
		
		* duration of coverage
		preserve
		keep ptid id lw201 lw401 zio1days zio2days zio3days
		qui egen nnum = rownonmiss(lw*01 *days)
		qui keep if nnum>0
		collapse (lastnm) lw*01 *days, by(id)
		qui egen ziodays = rowtotal(zio*days)
		label define howoften 1 "nearly always" 2 "most of the time" ///
			3 "~half the time" 4 "infrequently"
		label val lw* howoften
		label var lw201 "LW use, week 2"
		label var lw401 "LW use, week 4"
		label var ziodays "total days of Zio patch use"
		label var zio1days "days of first Zio patch use"
		label var zio2days "days of second Zio patch use"
		label var zio3days "days of third Zio patch use"
		qui save ./dta/coverage, replace
		restore
		
		* AF lags by self report
		preserve
		keep id ptid week etohaf*
		qui save ./dta/aflags, replace
		restore
		
		* triggers
		foreach i in 2 4 5 {
			preserve
			qui keep if week==`i'
			keep id ptid week lwdttm? lwdttm?? zio?dttm? zio?dttm??
			if `i'~=2 qui append using ./dta/triggers, force
			qui save ./dta/triggers, replace
			restore
		}
		
		* date-times
		use ./dta/triggers, clear
		foreach x of varlist lwdttm? lwdttm?? zio?dttm? zio?dttm?? {
			rename `x' tmp
			local type: type tmp
			if "`type'"=="byte" qui tostring tmp, force replace
			qui gen `x' = clock(tmp,"MD20Yhm") if ~missing(tmp)
			format `x' %tc
			drop tmp
		}	
		* change LW time zone from Central to Pacific
		foreach x of varlist lwdttm? lwdttm?? {
			qui replace `x' = `x'-2*3600000
		}
		
		* make a single sequence of triggers
		qui reshape long lwdttm zio1dttm zio2dttm zio3dttm, i(id week) j(trno)
		egen nnm = rownonmiss(lwdttm zio1dttm zio2dttm zio3dttm)
		qui drop if nnm==0
		drop nnm
		rename lwdttm dttm1
		rename zio1dttm dttm2
		rename zio2dttm dttm3
		rename zio3dttm dttm4		
		forvalues i = 1/4 {
			qui bysort id: egen fdttm`i' = min(dttm`i')
			qui bysort id: egen ldttm`i' = max(dttm`i')
			label var fdttm`i' "first trigger from device `i'"
			label var ldttm`i' "last trigger from device `i'"
		}
		qui replace ldttm3 = ldttm3-365*24*3600000 if id=="094_RC"
		format fdttm* ldttm* %tc
		qui reshape long dttm, i(id week trno) j(source)
		qui drop if missing(dttm)
		drop week
		qui duplicates drop id dttm, force
		
		* add days of coverage with each device
		qui merge m:1 id using ./dta/coverage, nogen
		qui gen days = max(0,(5-lw201))*14/4 if source==1
		qui replace days = max(0,(5-lw401))*14/4 if source==1 & missing(days)
		forvalues i = 1/3 {
			qui replace days = zio`i'days if source==`i'+1
		} 
		qui bysort id (dttm): replace trno = _n
		qui reshape wide dttm source days, i(id) j(trno)
		label define source 1 "LW" 2 "Zio patch 1" 3 "Zio patch 2" 4 "Zio patch 3"
		label val source* source
		format dttm* %tC
		forvalues i = 1/89 {
			label var dttm`i' "date/time of triggering event `i'"
			label var source`i' "source device triggering event `i'"
		}
		qui replace dttm41 = dttm40+3600000 if id=="094_RC"
		qui replace dttm1 = dttm1+365*24*3600000 if inlist(id,"050_MDS","061_WDK")
		qui compress
		*des
		*sum
		qui save ./dta/triggers, replace
	}
	
	if "`af_episodes'"~="" {
		qui import excel using "./excel/AF Episodes v4 (2).xlsx", ///
			clear firstrow cellrange(a1:h5781)
		rename *, lower
		rename af* *
		rename patchid device
		rename onset tmp
		qui gen onset = clock(tmp,"DMYhms")
		qui replace onset = clock(tmp,"MDYhm") if missing(onset)
		qui replace onset = clock(tmp,"MD20Yhms") if missing(onset)
		format onset %tc
		duplicates drop device onset, force
		format onset %tC

		drop tmp type
		rename duration tmp
		qui replace tmp = "0"+tmp if substr(tmp,3,1)==":"
		qui replace tmp = "00"+tmp if substr(tmp,2,1)==":"
		qui gen duration = real(substr(tmp,1,3))+real(substr(tmp,5,2))/60+ ///
			real(substr(tmp,8,2))/3600
		drop tmp
		qui save ./dta/af_episodes, replace

		* append patient IDs
		qui import excel using "./excel/AF Burden and Patch ID.xlsx", ///
			clear firstrow cellrange(a1:m108)
		rename studyid id
		rename patch? device?
		rename patch?dateon dateon?
		rename patch?dateoff dateoff?
		rename patch?afburden burden?
		keep id device? dateon? dateoff? burden?
		foreach x in dateon2 dateoff2 {
			rename `x' tmp
			qui gen `x' = clock(tmp,"MD20Yhm")
			qui replace `x' = cofd(date(tmp,"MD20Y")) if missing(`x')
			format `x' %tc
			drop tmp
		}
		qui destring burden3, replace
		qui reshape long device dateon dateoff burden, i(id) j(devno)
		qui replace dateoff = cofd(mdy(8,1,2017)) ///
			if dofc(dateon)==mdy(7,27,2017) & dofc(dateoff)==mdy(8,1,2016)
		qui replace dateoff = cofd(mdy(9,26,2017)) ///
			if dofc(dateon)==mdy(9,12,2017) & dofc(dateoff)==mdy(9,26,2016)
		qui replace dateon = cofd(mdy(10,26,2017)) ///
			if dofc(dateon)==mdy(10,26,2027) & dofc(dateoff)==mdy(11,9,2017)
		qui replace dateon = cofd(mdy(11,9,2017)) ///
			if dofc(dateon)==mdy(11,9,2022) &dofc(dateoff)==mdy(11,22,2017)
		qui replace dateon = dateon-365*24*3600000 ///
			if id=="036_CAS" & year(dofc(dateon))==2017
		qui replace dateoff = dateoff-365*24*3600000 ///
			if id=="036_CAS" & year(dofc(dateoff))==2017
		qui gen days = dofc(dateoff)-dofc(dateon)+1
		qui drop if missing(device)
		qui gen source = 1 if substr(device,1,9)=="LifeWatch"
		qui bysort id (devno): replace source = _n+1 if substr(device,1,1)=="N"
		label define source 1 "LW" 2 "Zio patch 1" 3 "Zio patch 2" 4 "Zio patch 3"
		label val source* source
		format dateon dateoff %tc
		qui merge 1:m device using ./dta/af_episodes, nogen keep(1 3)
		qui merge m:1 id using ./dta/insample, keep(3) nogen
		label var device "LW/Zio device ID"
		label var beats "AF beats"
		label var minhr "min AF beats/hour"
		label var maxhr "max AF beats/hour"
		label var avghr "average AF beats/hour"
		label var onset "AF episode onset time"
		label var duration "AF episode duration (hours)"
		label var days "days of device use"
		label var devno "device #"
		label var source "device"
		qui compress
		qui save ./dta/af_episodes, replace
	}
	
	if "`scram'"~="" {
		use ./dta/scram_epdatat_ebrac, clear
		rename *time* *tm*
		qui gen dt = dofc(sttm)
		preserve
		use ./dta/scram_daydatat, clear
		rename *, lower
		rename date dt
		keep id dt
		qui save tmp, replace
		restore
		qui merge m:1 id dt using tmp, nogen keepusing(id dt)
		erase tmp.dta
		qui merge m:1 id using ./dta/insample, nogen keep(3)
		label var id "id"
		label var epnum "episode #"
		label var stday "day of episode start"
		label var sttm "date/time of episode start"
		label var dt "date of episode start"
		label var endday "day of episode end"
		label var endtm "date/time of episode end"
		label var pkday "day of episode TAC peak"
		label var pktm "date/time of episode TAC peak"
		label var pktac "episode TAC peak"
		label var avtac "average episode TAC"
		label var absrate "absorption rate"
		label var elimrate "elimination rate"
		label var auc "AUC for episode"
		label var tmhrs "episode duration (hours)"
		label var amsconf "episode AMS confirmable"
		label var tamper "epiosde during AMS suggested tamper"
		label var sttm_ebrac "date/time eBrAC first detectable"
		label var stdf_ebrac "eBrAC to TAC latency"
		label var pk_ebrac "peak eBrAC value"
		label var pktm_ebrac "date/time of eBrAC peak"
		label var pklat_ebrac "eBrAC start to peak (hours)"
		label var pktdf_ebrac "eBrAC to TAC peak latency"
		label var auc_ebrac "eBrAC AUC"
		label var pktac "peak TAC"
		format dt %td
		qui compress
		des
		sum
		qui save ./dta/scram, replace
	}
	
	if "`wristas'"~="" {
		use ./dta/wristas_combined_05122020, clear
		rename *, lower
		rename temp* temp
		rename date_time dttm
		rename gp event_n
		rename drink_event event
		qui replace id = upper(substr(id,2,.))
		qui gen bac = alcohol/1000
		label var dttm "date/time of event"
		label var event_n "event #"
		label var event "alcohol event"
		label var bac "blood alcohol concentration"
		label define yesno 1 "yes" 0 "no"
		label val event yesno
		*bysort id (evnet_n dttm): list evnt_n dttm evnt alc temp, noobs clean
		keep id dttm event_n event bac
		order id dttm event_n event bac
		qui merge m:1 id using ./dta/insample, nogen keep(3)
		qui compress
		des
		sum
		qui save ./dta/wristas, replace
	}
	
	if "`peth'"~="" {
		use "./dta/PETH_PATCH 10202020.dta", clear
		rename record_id ptid
		drop id
		rename studyid id
		rename peth_result peth_det
		rename peth_quant peth_amt
		qui encode redcap, gen(tmp)
		qui recode tmp (1=2 "week 2") (2=4 "week 4") (3=5 "early termination"), ///
			gen(week)
		keep id ptid week peth_amt peth_det
		label var id "study ID"
		label var ptid "numeric study ID"
		label var week "study week"
		label var peth_det "detectable PeTH"
		label var peth_amt "PeTH quantitation"
		
		label define yesno 1 "yes" 0 "no"
		label val peth_det yesno
		
		aorder
		order ptid id week peth_det peth_amt
		duplicates drop 
		sort ptid week
		des
		sum
		qui compress
		qui save ./dta/peth, replace
	}

	if "`plotdata'"~="" {
		* AF episodes from patch
		use ./dta/af_episodes, clear
		qui duplicates drop id onset, force
		qui bysort id (onset): gen afno = _n
		qui sum afno
		local mafno = r(max)
		qui bysort id: egen fdateon = min(dateon)
		qui bysort id: egen ldateoff = max(dateoff)
		keep id fdateon ldateoff onset afno
		qui reshape wide onset, i(id) j(afno)
		
		* add triggers
		qui merge 1:1 id using ./dta/triggers, nogen keep(1 3) keepusing(id dttm*)
		
		* add SCRAM data
		preserve
		qui use ./dta/scram, clear
		keep id pktm_ebrac pk_ebrac auc_ebrac
		qui duplicates drop
		qui bysort id (pktm): gen scno = _n
		qui sum scno
		local mscno = r(max)
		qui reshape wide pktm_ebrac pk_ebrac auc_ebrac, i(id) j(scno)
		qui save tmp, replace
		restore
		qui merge 1:1 id using tmp, nogen keep(1 3)
		erase tmp.dta
		
		* rationalize start and stop days
		qui egen fonset = rowmin(onset*)
		qui egen lonset = rowmax(onset*)
		qui egen fdttm = rowmin(dttm*)
		qui egen ldttm = rowmax(dttm*)
		qui egen fpktm = rowmin(pktm*)
		qui egen lpktm = rowmax(pktm*)
		qui gen start = min(fdateon,fonset,fdttm,fpktm)
		qui gen stop = max(ldateoff,lonset,ldttm,lpktm)
		
		* numeric id
		qui encode id, gen(ptid)
		
		* map date-times to date and period
		foreach x of varlist start stop onset* dttm* pktm* {
			qui date_period `x'
		}

		qui gen days = stop_dt-start_dt+1
		qui expand days
		qui bysort ptid: gen day = _n
		qui bysort ptid: gen dt = start_dt+day-1
		qui drop if day>28
		qui gen week = ceil(day/7)
		format dt start_dt stop_dt %td
		format start stop %tc
		qui expand 48
		qui bysort ptid dt: gen period = _n
		
		qui drop if ~inrange(dt,start_dt,stop_dt) ///
			| (dt==start_dt & start_period>period) ///
			| (dt==stop_dt & stop_period<period)
			
		* any AF
		qui gen anyaf = 0 
		forvalues i = 1/`mafno' {
			qui replace anyaf = 1 if onset`i'_dt==dt & onset`i'_period==period
		}
		* # button pushes
		qui gen nbp = 0
		forvalues i = 1/89 {
			qui replace nbp = nbp+1 if dttm`i'_dt==dt & dttm`i'_period==period
		}
		* SCRAM eBrAC and AUC
		qui gen pk_ebrac = .
		qui gen auc_ebrac = .
		forvalues i = 1/`mscno' {
			qui replace pk_ebrac = max(pk_ebrac,pk_ebrac`i') ///
				if pktm_ebrac`i'_dt==dt & pktm_ebrac`i'_period==period
			qui replace auc_ebrac = auc_ebrac+auc_ebrac`i' ///
				if pktm_ebrac`i'_dt==dt & pktm_ebrac`i'_period==period ///
					& ~missing(auc_ebrac) & ~missing(auc_ebrac`i')
			qui replace auc_ebrac = auc_ebrac`i' ///
				if pktm_ebrac`i'_dt==dt & pktm_ebrac`i'_period==period ///
					& missing(auc_ebrac) & ~missing(auc_ebrac`i')
		}
		* PeTH
		preserve
		use ./dta/peth, clear
		keep id week peth*
		qui drop if week==5
		collapse peth_amt (max) peth_det, by(id week)
		qui reshape wide peth_det peth_amt, i(id) j(week)
		qui save tmp, replace
		restore
		qui merge m:1 id using tmp, nogen keep(1 3)
		erase tmp.dta
		qui gen pdet = .
		qui gen pamt = .
		foreach x in det amt {
			foreach j in 2 4 {
			qui replace p`x' = peth_`x'`j' if week<=`j' & ~missing(peth_`x'`j')
			}
		}
		keep id day period anyaf nbp pk_ebrac auc_ebrac pdet pamt
		order id day period anyaf nbp pk_ebrac auc_ebrac pdet pamt
		sort id day period
		label var day "study day"
		label var period "30-minute period within day"
		label var anyaf "any AF in 30-minute period"
		label var nbp "# button pushes in 30-minute period"
		label var pk_ebrac "peak SCRAM eBrAC in period"
		label var auc_ebrac "total SCRAM eBrAC AUC in period"
		label var pdet "detectable PeTH at end of 2-week period"
		label var pamt "PeTH amount at end of 2-week period"
		des
		sum
		qui export delimited id day period anyaf nbp pk_ebrac auc_ebrac pdet pamt ///
			using ./csv/plotdata3.csv, replace
		}
	
end

cd "~/work/parnassus/marcus/projects/monitors/programs"
cd ./csv
!gunzip -f *.csv.gz
cd ../excel/
!gunzip -f *.xlsx.gz
cd ../dta/
!gunzip -f *.dta.gz
cd ../

set tracedepth 1
set trace off
cap log close
log using makedata, replace
makedata, main //scram wristas peth plotdata
log close

cd ./dta
!gzip -f *.dta
cd ../Excel
!gzip -f *.xlsx
cd ../csv
!gzip -f *.csv
cd ../
