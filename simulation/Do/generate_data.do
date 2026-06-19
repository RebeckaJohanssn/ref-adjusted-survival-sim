// generate_data.do
// Generates one set of simulated datasets for a single scenario/variant,
// writing $Nsim simulated samples (each of size $Nobs) to disk.
//
// Designed to be called as a separate Stata process per scenario/variant,
// e.g. from generate_data_parallel.do, with three positional arguments:
//   1: ROOT directory of the project
//   2: scenario number
//   3: scenario variant (a/b/c)
// Can also be run directly for a single scenario/variant by setting
// ROOT/scen/val manually below instead of passing them as arguments.
//
// NOTE on background mortality: unlike generate_truth.do (which uses
// popmort_males for all simulated individuals), this script merges on sex
// as well as age/year, using popmort2022, i.e. sex-specific mortality
// rates ARE used here.
//
// Requires globals.do to have been run (defines $Nsim, $Nobs, $maxt, etc.)
// and the relevant Scenario`scen'/`val'.do file to exist.
// Requires popmort2022 (population mortality rates by sex/age/year) to be
// available in the working directory or on the adopath.
// Output directory ${ROOT}/Simulated_Data/Scenario`scen'/`val'/ must exist
// before running, or the save command will fail.

global ROOT `1'
local scen `2'
local val `3'

//set seed 136677
// NOTE: no seed is set here. When run in parallel via generate_data_parallel.do,
// each spawned process is a separate Stata session with its own random state.

cd $ROOT
do "${ROOT}/globals.do"

do "${ROOT}/Do/Scenarios/Scenario`scen'/`val'.do"
forvalues j = 1/$Nsim {
	display "Dataset `j'"
	quietly {
	clear
	set obs $Nobs
	gen id = _n
	gen sex = cond(runiform() < 0.5, 1, 2) //approx. equal sex distr.
	gen sexc = sex-1
	gen agediag = floor(rnormal(${meanage}, ${sdage}))
	gen agediagc = agediag-$meanage
	replace agediag = min(agediag, 99)
	gen datediag = mdy(1,1,2020) //everyone is diagnosed in 2020
	gen yydx = year(datediag)
	egen agegroup = cut(agediag), at(0 45 55 65 75 130) icodes 
	tab agegroup, gen(agegroup) //create agegroups

	// generate cause-specific survival times
	survsim t_cancer d_cancer, distribution(weibull) ///weibull distribution
	lambdas($lambda) gammas($gamma) maxt($maxt) ///
	covariates(agediagc ${ageloghr} sexc ${sexloghr})
	
	// generate time to death from other causes using sex-specific
	// population mortality rates (popmort2022), updating age/year each
	// follow-up year as in generate_truth.do
	gen _age =.
	gen _year=.
	forvalues i=0/`=$maxt-1' {
		capture drop _merge prob //drop _merge and prob, but if an error occurs keep running the program
		replace _age = min(floor(agediag + `i'), 99)  //increase age by one
		replace _year = min(yydx + `i', 2022) //increase year by one, same as pmaxyear option
		quietly merge m:1 _age _year sex using popmort2022, keep(matched master) //merge with popmort file
		gen t_other`i' = (-log(runiform()))/(-ln(prob)) 
		replace t_other`i'=1 if t_other`i'>=1 //if survival time is longer than a year then set it to 1, i.e. the person survived the whole year.
	}

	gen t_other=. //set variable t_other to missing for everyone
	forvalues i = 0/`=$maxt-1' {
		quietly replace t_other=t_other`i'+`i' if t_other`i'<1 & t_other == . 
		} //total survival time is nr. of full years survived plus a fraction of one year survived. Only replace if there is no previously recorded survival time.

	replace t_other = $maxt if t_other==.

	drop t_other? t_other?? _age _year _merge rate 

	// all-cause survival time and outcome
	gen time = min(t_other, t_cancer) //all-cause survival time
	gen died = time < ${maxt} //death indicator, if the survival time is <=maxt
	replace time = ceil(time*365.241)/365.241

	// save simulated dataset
	save "${ROOT}/Simulated_Data/Scenario`scen'/`val'/simdata`j'", replace
	}
}
