data plant;
   infile 'd:\plant.txt' delimiter='09'x;
   /* delimiter='09'x means that the observations are seperated by tab */
   input number $ x1 x2;
   run;

/* (a) */
proc univariate data=plant normal;
   var x1 x2;
   symbol1 value=dot color=blue;
   qqplot x1 x2 /normal(mu=est sigma=est color=red w=2) square;
   run;

/* (b) */
%include 'd:\chisq_plot.sas'; /* include the macro */
%chisq_plot(plant,x1 x2); /* run the macro */

/* (c) */
proc iml;
   use plant;
   read all var {x1 x2} into x;
   close plant;

   /* (i) */
   mu = {50,1000};
   n = nrow(x);
   xbar = x[:,]`;
   s = (x`*x-n*xbar*xbar`)/(n-1);
   tsq = n*(xbar-mu)`*inv(s)*(xbar-mu);
   print 'Hotelling T-square Statistic',, tsq;

   /* (ii) */
   x2 = 200;
   ex1x2 = xbar[1]+s[1,2]*inv(s[2,2])*(x2-xbar[2]);
   print 'MLE of E(x1|x2=200)',, ex1x2;

quit;
