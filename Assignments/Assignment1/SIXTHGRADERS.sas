data sixgraders;
input x1 x2 @@;
cards;
28.87  12.38  12.20  11.33
20.10  10.34  22.55  10.51
69.05  14.10  14.30  11.72
65.40  14.20  31.79  12.33
29.59  12.00  11.60  11.38
44.82  12.82  68.47  14.02
77.37  13.70  42.64  12.95
24.67  11.53  16.70  12.05
65.01  13.11  86.27  15.37
9.99   11.23  76.73  14.74
;

/*Question 2*/
/*Part (a)*/
data new_sixgraders;
set sixgraders;
arcsinx1 = arsin(sqrt(x1/100));
sqrtx1 = sqrt(x1);
	
proc univariate data=new_sixgraders normal;
	var x1 arcsinx1 sqrtx1;
	symbol1 value=dot color=blue;
	qqplot x1 arcsinx1 sqrtx1 /normal(mu=est sigma=est color=red w=2) square;
run;

/*Part (b): Sqrt is the best*/
proc iml;
	use new_sixgraders;
	read all var{sqrtx1 x2} into x;
	close new_sixgraders;
	
	n = nrow(x);
	p = ncol(x);
	xbar = x[:,]`;
	S = (x`*x - n*xbar*xbar`)/(n-1);
	mu0 = {0.5, 12};
	mu0[1] = sqrt(0.5);
	
	tsq = n*(xbar-mu0)`*inv(S)*(xbar-mu0);
	alpha = 0.05;
	cv = (n-1)*p/(n-p)*finv(1-alpha,p,n-p);
	pval = 1-probf((n-p)/(p*(n-1))*tsq,p,n-p);
	print tsq[format=10.4] cv[format=10.4] pval;
run;