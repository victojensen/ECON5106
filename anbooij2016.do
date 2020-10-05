clear all
use work_base

// 
global lopts lc(black black black blue gs11) lp(dash solid dot solid solid)

// define variable lists
local xvars male age cito fes
local outcomes matched drop reswis_raw reslang_raw resother_raw nsubj56 sctrack iwisB 
local heoutcomes matchedhe chosesc switch1yr startwage_raw 

// define specifications
global average i.xpos
global linear c.ist_norm
global quadratic c.ist_norm##c.ist_norm
global quad_split c.ist_norm##c.ist_norm c.z#(c.ist_norm##c.ist_norm)
global cubic c.ist_norm##c.ist_norm##c.ist_norm
global zoom $linear if inrange(ist_norm,-18,18)
global donut $quadratic if !inrange(ist_norm,-4,4)

// discretize the support of the running variable for the graphs
g x    = autocode(ist_norm + 1, 40, -80, 80)
g xpos = x + 100 // stata cannot make dummies (i.x) of negative values

// Table 1
sum male `xvars' ist_raw treat `outcomes' `heoutcomes'
tab cohort 

// in the rest of the analysis use log(wage)
replace startwage_raw = log(startwage_raw) 
// and normalized test scores
qui foreach v of var cito fes res* {
	forv c=1998/2010 {
		sum `v' if cohort==`c'
		replace `v' = (`v' - r(mean)) / r(sd) if cohort==`c'
	}
	replace `v' = 0 if missing(`v')
}

// little program that makes simple RD graph
program rdgraph
	args y
	macro shift 1 // chop off y from the arguments so that we can pass the rest to twoway below
	quietly {
		// residualize the outcome
		sum `y'
		local m = r(mean)
		reg `y' i.cohort `xvars'
		predict p`y', res
		replace p`y' = p`y' + `m'
		
		foreach spec in average linear quadratic cubic zoom donut {
			reg p`y' z $`spec'
			predict y_`spec' if e(sample)
			label var y_`spec' `spec'
		}
	}
	twoway (sc y_a x, ms(Oh)) ///
		(line y_l-y_d ist_norm if ist_norm>=0, $lopts sort) ///
		(line y_l-y_d ist_norm if ist_norm< 0, $lopts sort), ///
		legend(order(1 2 3 4 5 6)) legend(pos(11) ring(0)) ///
		ytitle(`y') xline(0) name(rd_`y', replace) `*'
	capture drop y_* p`y'
end

// Figure 1
rdgraph treat xsca(alt) 
twoway hist ist_norm, start(-80) width(4) xline(0) ysca(reverse) ylabel(0) fysize(30) name(hist_ist, replace)
gr combine rd_treat hist_ist, col(1) xcommon imargin(b=0 t=0)

// Table 2 - first stage
reg treat z i.cohort $quadratic, cluster(group)
est sto nox
qui foreach spec in quadratic linear cubic quad_split zoom donut {
	reg treat z i.cohort `xvars' $`spec', cluster(group)
	est sto `spec'
}
esttab *, keep(z) b(2) se stat(N r2) mtitles

// Table 2 - balancing
est clear
qui foreach y in `xvars' {
	reg `y' z i.cohort $quadratic, cluster(group)
	est sto `y'
}
esttab *, keep(z) b(2) se stat(N r2) mtitles

// balancing F test
// using suest
est clear
qui foreach y in `xvars' {
	reg `y' z i.cohort $quadratic
	est sto `y'
}
qui suest `xvars', cluster(group)
test z
// or alternatively using reg
qui reg z `xvars' i.cohort $quadratic, cluster(group)
testparm `xvars'

// Figure 2 - balancing
rdgraph fes
rdgraph cito

// Figure 3 - reduced forms
rdgraph reswis_raw 
rdgraph reslang_raw 
rdgraph resother_raw

// Table 4 - high school
qui foreach sample in 1 male==0 male==1 {
	est clear
	foreach y of local outcomes {
		ivregress 2sls `y' (treat = z) i.cohort `xvars' $quadratic if `sample', cluster(group)
		est sto `y'
	}
	noi esttab *, keep(treat) b(2) se mtitles title(`sample')
}

// Table 5 - higher education
qui foreach sample in 1 male==0 male==1 {
	est clear
	foreach y of local heoutcomes {
		ivregress 2sls `y' (treat = z) i.cohort `xvars' $quadratic if `sample', cluster(group)
		est sto `y'
	}
	noi esttab *, keep(treat) b(2) se mtitles title(`sample')
}

// Table A1 - specification checks
qui foreach spec in quadratic quad_split zoom donut {
	est clear
	foreach y of local outcomes {
		ivregress 2sls `y' (treat = z) i.cohort `xvars' $`spec', cluster(group)
		est sto `y'
	}
	noi esttab *, keep(treat) b(2) se mtitles title(`spec')
}


// Extra: estimate x0
g xc = .
g r2 = 0
qui forv c=1998/2010 {
	forv x=90/120 {
		g p = ist_raw>=`x'
		reg treat p $quadratic if cohort==`c'
		cap replace xc = `x' if cohort==`c' & e(r2) > r2
		cap replace r2 = max(e(r2), r2) if cohort==`c'
		drop p
	}
}
// check that they are the same
tab x0 xc
