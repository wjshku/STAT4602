/* ChiSQ Plot for multivariate normality testing */
/* PROC PLOT, INSTEAD OF PROC GPLOT, IS USED IN THIS MACRO */

%macro chisq_plot(dname,group);
proc iml;
   use &dname;
   read all var {&group} into x;
   close &dname;
   
   n = nrow(x);
   p = ncol(x);
   xbar = x[:,]`;
   s = (x`*x-n*xbar*xbar`)/(n-1);
   free di2 xi2;
   do i = 1 to n;
      di2 = di2//((x[i,]`-xbar)`*inv(s)*(x[i,]`-xbar));
      xi2 = xi2//cinv((i-.5)/n,p);
     /* // vertical concatenation */
     /*CINV(P, df, nc) P-quantile of chi-square distribution (inverse of CDF) */
   end;
   maxx = int(max(xi2))+1;
   call symput('maxp',left(char(maxx)));  
   di2 = di2//maxx;
   call sort(di2,1);
   create result var {xi2 di2};
   append;
   close result;
quit;

data result;
   set result;
run;

proc plot data=result ;
   plot di2*xi2='*' ;
  title 'Plot of Chi-square Distribution vs Ordered Squared Mahalanobis Distance';
  title2 'x axis: Ordered Squared Mahalanobis Distance';
  title3 'y axis: Chi square';
   run;
quit;
proc datasets lib=work;
   delete result;
   run;
quit;
%mend;
