// generate_truth.do
// Generates "true" cause-specific and all-cause survival quantities
// for each scenario, by simulating a very large population (size $Ntruth)
// and recording the resulting survival probabilities to file.

// NOTE: other-cause (background) mortality rates are taken from
// popmort_males for all simulated individuals, i.e. male population
// mortality rates are used regardless of simulated sex to get reference
// adjusted all-cause survival using male other-cause mortality rates. 
// To use sex-specific background mortality instead, the merge step below would 
// need to be split/merged by sex against the corresponding popmort file(s).

// Requires globals.do to have been run first (defines $Nscenarios, $Ntruth, $maxt).
// Requires the popmort_males file to be available in the working directory.

set seed 215789

forvalues scen = 1/$Nscenarios {
	foreach val in a b c {
		// load scenario-specific parameters (lambda, gamma, hazard ratios, etc.)
		do "${ROOT}/Do/Scenarios/Scenario`scen'/`val'.do"
		di "Scenario `scen'`val'"
		quietly {
		clear
		set obs $Ntruth
		gen id = _n
		gen sex = cond(runiform() < 0.5, 1, 2) //approx. equal sex distr.
		gen sexc = sex-1
		gen agediag = floor(rnormal(${meanage}, ${sdage}))
		replace agediag = min(agediag, 99)
		gen agediagc = agediag-${meanage}
		gen datediag = mdy(1,1,2020)
		gen yydx = year(datediag) //everyone is diagnosed in 2020
		
		// generate cause-specific survival times
		survsim t_cancer d_cancer, distribution(weibull) ///
		lambdas($lambda) gammas($gamma) maxt($maxt) ///
		covariates(agediagc ${ageloghr} sexc ${sexloghr})
		
// Generate time to death from other causes using population mortality
// rates. NOTE: this currently uses popmort_males, i.e. all simulated individuals
// (regardless of their simulated `sex`, to use males as the reference population 
// for reference adjusted measures) are assigned  background mortality rates for males.
// This can be changed to any population mortality file to use other rates as the reference.
		gen _age =.
		gen _year=.
		forvalues i=0/`=$maxt-1' {
			capture drop _merge prob //drop _merge and prob, but if an error occurs keep running the program
			replace _age = min(floor(agediag + `i'), 99) //increase age by one
			replace _year = min(yydx + `i', 2022) //increase year by one, same as pmmaxyear option
			quietly merge m:1 _age _year using popmort_males, keep(matched master) //merge with popmort file
			gen t_other`i' = (-log(runiform()))/(-ln(prob)) //generate uniform number and transform to exponential(1) r.v., -ln(prob) is the hazard rate, so this is the exponential survival time with hazard -ln(prob). Generates a time to death from expected mortality
			replace t_other`i'=1 if t_other`i'>=1 //if survival time is longer than a year then set it to 1, i.e. the person survived the whole year.
		}
		
		gen t_other= . //set variable timeback to missing for everyone
		forvalues i = 0/`=$maxt-1' {
			quietly replace t_other=t_other`i'+`i' if t_other`i'< 1 & t_other == . 
			} //total survival time is nr. of full years survived plus a fraction of one year survived. Only replace if there is no previously recorded survival time.
		replace t_other = $maxt if t_other==.
			
		drop t_other? t_other?? _age _year _merge rate 
		
		// All-cause survival time is the minimum of cause-specific and
		// other-cause survival times.
		gen time = min(t_other, t_cancer) //all-cause survival time
		gen died = time < ${maxt} //death indicator, if the survival time is <=maxt
		//replace time = ${maxt} if died == 0 //censor at time maxt
		replace time = ceil(time*365.241)/365.241
		
		// Compute true net survival (based on cause-specific survival) and
		// true reference adjusted all-cause survival, by sex, at
		// 1, 5, 10 and 15 years. Write results as global macros to a .do
		// file so they can be loaded later for comparison against estimates.
		file open truth_file using "${ROOT}/Truths/Scenario`scen'/`val'.do", write replace
		
		foreach t in 1 5 10 15{
			forvalues s=1/2 {
				count if sex==`s'
				local N = r(N)
				count if (t_cancer >= `t') & sex==`s'
				local true_net`t's`s' = r(N)/`N'
				count if (time >= `t') & sex==`s'
				local true_ra`t's`s' = r(N)/`N'
				file write truth_file "global true_net_`t's`s' =  `true_net`t's`s''" _n
				file write truth_file "global true_ra_`t's`s' =  `true_ra`t's`s''" _n
			}
		}	
		
		file close truth_file
	}
}
}
