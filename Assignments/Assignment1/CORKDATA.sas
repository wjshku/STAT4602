data cork;
input N E S W @@;
cards;
72 66  76  77  91  79  100 75
60 53  66  63  56  68  47  50 
56 57  64  58  79  65  70  61 
41 29  36  38  81  80  68  58 
32 32  35  36  78  55  67  60 
30 35  34  26  46  38  37  38 
39 39  31  27  39  35  34  37 
42 43  31  25  32  30  30  32 
37 40  31  25  60  50  67  54 
33 29  27  36  35  37  48  39 
32 30  34  28  39  36  39  31 
63 45  74  63  50  34  37  40 
54 46  60  52  43  37  39  50 
47 51  52  45  48  54  57  43 
;


/* Question 1*/
/*Part (a)*/
/*Q-Q Plot for all 4 variables*/
proc univariate data=cork normal;
var N E S W;
symbol1 value=dot color=blue;
qqplot N E S W /normal(mu=est sigma=est color=red w=2) square;
run;
/*ChiSQ Plot for multivariate normality testing */
%include ’STAT4602A1/chisq_plot.sas’;
%chisq_plot(cork,N E S W);

/*Part (b)*/
proc iml;
	use cork;
	read all var{N E S W} into x1;
	close cork;
	
	n1 = nrow(x1);
	mean1 = x1[:,]`;
	W1 = (x1`*x1 - n1*mean1*mean1`);
	S1 = W1/(n1-1);
	print n1[format=5.0] mean1[format=10.4] S1[format=10.4] W1[format=10.4];
	alpha = 0.05;
	c = {1 -1 0 0,
		 0  1 -1 0,
		 0 0 1 -1};
	q = nrow(c);
	tsq = n1*(c*mean1)`*inv(c*s1*c`)*(c*mean1);
	cv = (n1-1)*q/(n1-q)*finv(1-alpha,q,n1-q);
	pval = 1-probf(tsq*(n1-q)/(q*(n1-1)),q,n1-q);
	print tsq[format=10.4] cv[format=10.4] pval[format=10.4];
	
/*Part (c) (i)*/
proc iml;
	use cork;
	read all var{N E S W} into x1;
	close cork;
	
	n = nrow(x1);
	mean1 = x1[:,]`;
	W1 = (x1`*x1 - n*mean1*mean1`);
	S = W1/(n-1);
	print n[format=5.0] mean1[format=10.4] S[format=10.4] W1[format=10.4];
	alpha = 0.05;
	c = {1 0 -1 0,
		 0  1 0 -1};
	q = nrow(c);
	tsq = n*(c*mean1)`*inv(c*S*c`)*(c*mean1);
	cv = (n-1)*q/(n-q)*finv(1-alpha,q,n-q);
	pval = 1-probf(tsq*(n-q)/(q*(n-1)),q,n-q);
	print tsq[format=10.4] cv[format=10.4] pval[format=10.4];
	
	/*Part (c) (ii)*/
	alpha = 0.05;
	a = {1,-1,1,-1};
	merror = abs(tinv(alpha/(2*1),n-1))*sqrt(a`*S*a/n);
	ci1 = (a`*mean1-merror)||(a`*mean1+merror);
	
	ci = ci1;
	r = {"Contrast"};
	c = {"Lower Bound" "Upper Bound"};
	print ci[label="95% C.I" rowname=r colname=c format=10.4];

	
	