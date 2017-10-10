* Demand Example.
/*
- Uses data generated in OtherComps.jl
- which saves a file "demand_example.csv"
- contains a two firm market
- sequentially does demand computation while one firm goes from l1 to l2 to l3
*/


import delimited /Users/austinbean/Desktop/dynhosp/demand_example.csv, clear


* 1 is the firm changing levels, 2 is the firm not doing so.  
collapse (mean) mn_1 = main  mn_2 = alt (count) ocount_1 = main ocount_2 = alt (sd) sd_1 = main sd_2 = alt  (p5) p5_1 = main  p5_2 = alt  (p95) p95_1 = main p95_2 = alt, by(lev)

reshape long mn_ ocount_ sd_ p5_ p95_, i(lev) j(nn)


* Replace the x-values and replot...
replace lev = 1.5 if lev == 1 & nn == 1
replace lev = 2.5 if lev == 2 & nn == 1
replace lev = 3.5 if lev == 3 & nn == 1

graph twoway (bar mn_ lev if nn == 1, barw(0.4)) (bar mn_ lev if nn == 2, barw(0.4)) (rcap p95_ p5_ lev if nn == 1, lcolor(orange)) (rcap p95_ p5_ lev if nn == 2, lcolor(orange)), xlabel(1 "1" 1.5 "1" 2 "1" 2.5 "2" 3 "1" 3.5 "3") xline(1.75) xline(2.75) yscale(range(0 1200)) ytitle("Patients Admitted") xtitle("Facility Levels") graphregion(color(white)) legend(lab(1 "Investing") lab( 2 "Not Investing") lab(3 "5-95 %-tile") lab( 4 "5-95 %-tile" ) ) title("Simulated Demand Response to Firm Investment") subtitle("Grayson County, TX") note("Texoma Med. Center - Denison, TX"  "Wilson N Jones Med. Center - Sherman, TX")




graph save Graph "/Users/austinbean/Desktop/dynhosp/DemandSim.gph"
