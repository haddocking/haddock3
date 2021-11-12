remark parameter file for dipolar coupling calculations
set message off echo off end

bonds     OO  XX               1000.0      3.00
bonds     OO  YY               1000.0      3.00
bonds     OO  ZZ               1000.0      3.00

angle     XX OO YY             1000.0      90.0
angle     XX OO ZZ             1000.0      90.0
angle     YY OO ZZ             1000.0      90.0

NONBonded  XX		0.01	0.01	0.01 	0.01
NONBonded  YY		0.01	0.01	0.01 	0.01
NONBonded  ZZ		0.01	0.01	0.01 	0.01
NONBonded  OO		0.01	0.01	0.01 	0.01

set message on echo on end
