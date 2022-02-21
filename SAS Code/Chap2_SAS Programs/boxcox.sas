     /* Two macros, BOXCOX1 uses REG and BOXCOX2 uses GLM */

*----- Box-Cox Transforms for Dependent Variable in Regression -----*;
* BoxCox macro with positional parameters separated by commas.
* Invoke the macro by the following macro call :
*       %BOXCOX1(dname, yname, xname, lambda1, lambda2, width)
* where
* dname=name of SAS dataset,
* yname=name of dependent variable,
* xname=names of independent variables,
* lambda1=lower lambda value,
* lambda2=upper lambda values,
* width=interval width between lambdas.
* Example:     %boxcox1(c.class, weight, age height, -3, 3, 0.5)
* To get translated statements of the macro, add the following
* statement before the %BOXCOX1(  ,  ,  ,  ,  ,  )  call above :
*       OPTIONS MPRINT;
* The macro will print out the lambda values together with the
* corresponding root mean square error and the estimates. It
* also provides a hi-resolution smooth-plot of RMSE vs lambda.
*------------------------------------------------------------------*;

%macro BOXCOX1(dname,yname,xname,lambda1,lambda2,width);

 %let FLAG=0; 
 data _tmp_;
 set &dname;
 if (&yname <=0) then do;
   put "ERROR: Non-positive value of dependent variable &yname"
   "(obs=" _n_ ") " &yname;
   call symput('FLAG','1');
   end;
 run;

 %if (&FLAG) %then %do;
   %put BOXCOX: Cannot compute Box-Cox transformation;
   %goto done;
   %end;

  /* Find the number of intervals between lambda1 & lambda2 */
 data _null_;
 steps=round((&lambda2-&lambda1)/&width,1);
 call symput('steps', left(steps));


   /* Find the geometric mean and transformed Y for lambdas */
 data LOG;
 set &dname (keep= &yname &xname);
 logy=log(&yname);

 proc means noprint; var logy;
 output out=MEANLOG mean=meanlogy;

 data BOXCOX  (drop=meanlogy gmean );
 set LOG;
 if _n_=1 then do; set MEANLOG; gmean=exp(meanlogy); end;
 retain gmean ;
 array bc bc0-bc&steps;
 do over bc; lambda=&lambda1+(_i_-1)*&width;
    if lambda=0 then bc=gmean*logy;
    else  bc=(&yname**lambda-1)/lambda/gmean**(lambda-1);
 end;
        /* Regressions using all transformed Y */
 proc reg outest=RMSE (drop=_model_ _type_  _depvar_ bc0-bc&steps);
  model bc0-bc&steps = &xname /noprint;

        /* Collate RMSE with lambdas */
 data LAMBDA;
  lambda=&lambda1+(_n_-1)*&width; set RMSE;
 proc print; run;

      /* Next step reqiures hi-resolution graphics  */
 goptions gunit=cm ftitle=swissb htitle=0.6 
                               ftext=swissb htext=0.42;
 proc gplot ;
  symbol1 i=spline v=dot c=black h=0.2;
  axis1 label=(a=90);
  axis2 label=('Box-Cox Power (' f=cgreek 'l' f=swiss ')' );
  plot _rmse_*lambda=1 / vaxis=axis1 haxis=axis2 hminor=1 frame
                                          des="Box-Cox plot of RMSE * Lambda for &yname";
  label _rmse_ = "Root Mean Squared Error";
  title " Box-Cox Transform for Dependent Variable: &yname";
 run;
 
     /* Clean up datasets */
 %done:
 proc datasets mt=data nolist; delete LOG MEANLOG BOXCOX RMSE LAMBDA _tmp_;
  run;
 title1;  title2;
%mend;

*----- Box-Cox Transforms for Dependent Variable in ANOVA --------*;
* BoxCox macro with positional parameters separated by commas.
* Invoke the macro by the following macro call :
*    %BOXCOX2(dname, class, yname, effects, lambda1, lambda2, width )
* where
* dname= name of SAS dataset,
* class= variables to be used as class variables to form effects,
* yname= name of dependent variable,
* effects= design effects of the model,
* lambda1= lower lambda value,
* lambda2= upper lambda value,
* width= interval width between lambdas.
* Example:   %boxcox2(A.CROPS, A B, YIELD, X A|B X*A, -3, 3, 0.5 )
* To get translated statements of the macro, add the following
* statement before the %BOXCOX2(  ,  ,  ,  ,  ,  ,  )  call above :
*       OPTIONS MPRINT;
* The macro will print out the lambda values together with the
* corresponding root mean square error and the estimates.  It
* also provides a hi-resolution smooth-plot of RMSE vs lambda.
*------------------------------------------------------------------*;

%macro BOXCOX2(dname,class,yname,effects,lambda1,lambda2,width);

     /* Find the number of intervals between lambda1 & lambda2 */
 data _null_;
 steps=round((&lambda2-&lambda1)/&width,1);
 call symput('steps', left(steps));

     /* Find the geometric mean and transformed Y for lambdas */
 data LOG;
 set &dname;
 logy=log(&yname);

 proc means noprint; var logy;
 output out=MEANLOG mean=meanlogy;

 data BOXCOX  (drop=meanlogy gmean );
 set LOG;
 if _n_=1 then do; set MEANLOG; gmean=exp(meanlogy); end;
 retain gmean ;
 array bc bc0-bc&steps;
 do over bc; lambda=&lambda1+(_i_-1)*&width;
    if lambda=0 then bc=gmean*logy;
    else  bc=(&yname**lambda-1)/lambda/gmean**(lambda-1);
 end;
         /* Linear models using all transformed Y */
 proc glm noprint outstat=RMSE (keep=_TYPE_ DF SS);
  class &class;
  model bc0-bc&steps = &effects ;

         /* Collate RMSE with lambdas */
 data LAMBDA (drop= _TYPE_ DF SS);
  lambda=&lambda1+(_n_-1)*&width; set RMSE (where=(_TYPE_='ERROR'));
  RMSE=sqrt(SS/DF);
 proc print; run;

         /* Next step reqiures hi-resolution graphics  */
 proc gplot ;
  symbol1 i=spline v=dot c=black h=0.2;
  axis1 label=(a=90);
  axis2 label=('Box-Cox Power (' f=cgreek 'l' f=swiss ')' );
  plot rmse*lambda=1 / vaxis=axis1 haxis=axis2 hminor=1 frame
                                          des="Box-Cox plot of RMSE * Lambda for &yname";
  label rmse = "Root Mean Squared Error";
  title " Box-Cox Transform for Dependent Variable: &yname";
 run;
          /* Clean up datasets */
 proc datasets mt=data nolist; delete LOG MEANLOG BOXCOX RMSE LAMBDA;
  run;
 title1;  title2;

symbol1 i=none;
%mend;
