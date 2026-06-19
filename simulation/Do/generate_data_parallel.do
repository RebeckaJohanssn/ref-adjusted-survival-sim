// generate_data_parallel.do
// Launches generate_data.do as separate parallel Stata processes — one per
// scenario/variant combination — using winexec to spawn each as a background
// job. This is much faster than running scenarios sequentially, but:
//   - it is Windows-only (winexec is a Windows-specific command)
//   - it relies on $statapath pointing to a valid Stata executable
//     (set in globals.do)
//   - each process does NOT wait for the previous one to finish, so make
//     sure your machine has enough cores/memory to run $Nscenarios * 3
//     instances of Stata simultaneously
//
// Requires globals.do to have been run first (defines $Nscenarios, $statapath).

set seed 9731
forvalues s = 1/$Nscenarios {
	foreach val in a b c {
		winexec "${statapath}" /e run "${ROOT}/Do/generate_data.do" "${ROOT}" `s' `val'
	}
}
