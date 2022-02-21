data fev;
infile "STAT4602A1/FEV.DAT" delimiter=" ";
input Fev Pat Time;
run;

Data sex;
input Sex $;
datalines;
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
M
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
F
;
run;

Data fevmerge;
MERGE fev sex;
run;

Data fevm;
set fevmerge;
If Sex = "M" Then output;
run;

Data fevf;
set fevmerge;
If Sex = "F" Then output;
run;

/*Part (a)*/
title "Mean Response(Male)";
proc sgplot data=fevm;
   vline Time / response=Fev group=Time stat=mean limitstat=stderr;
   yaxis label='Mean +/- SEM';
run;


title "Mean Response(Female)";
proc sgplot data=fevf;
   vline Time / response=Fev group=Time stat=mean limitstat=stderr;
   yaxis label='Mean +/- SEM';
run;

/*Part (b)*/
proc sql noprint;
    create table have_aggregated as 
        select Fev
             , Sex
             , Time
             , Pat
        from fevm
        group by Pat
        order by Pat
    ;
quit;

proc transpose data=fevm
               out=fevmtrans(drop=_NAME_);
    id Time;
    by Pat;
    var Fev;
run;

proc transpose data=fevf
               out=fevftrans(drop=_NAME_);
    id Time;
    by Pat;
    var Fev;
run;

/*(i)*/
proc iml;
	use fevmtrans;
	read all var{"5" "10" "30" "60" "120"} into m;
	close fevmtrans;
	
	use fevftrans;
	read all var{"5" "10" "30" "60" "120"} into f;
	close fevftrans;
	
	mn = nrow(m);
	p = ncol(m);
	mbar = m[:,]`;
	Sm = (m`*m - mn*mbar*mbar`)/(mn-1);

	
	fn = nrow(f);
	p = ncol(f);
	fbar = f[:,]`;
	Sf = (f`*f - fn*fbar*fbar`)/(fn-1);
	
	S = ((mn-1)*Sm + (fn-1)*Sf)/(mn+fn-2);
	tsq = mn*fn/(mn+fn)*(mbar-fbar)`*inv(S)*(mbar-fbar);
	
	alpha = 0.05;
	cv = (mn+fn-2)*p/(mn+fn-p-1)*finv(1-alpha,p,mn+fn-p-1);
	pval = 1-probf((mn+fn-1-p)/(p*(mn+fn-2))*tsq,p,mn+fn-p-1);
	print tsq[format=10.4] cv[format=10.4] pval;
run;

/*(i)*/
proc iml;
	use fevmtrans;
	read all var{"5" "10" "30" "60" "120"} into m;
	close fevmtrans;
	
	use fevftrans;
	read all var{"5" "10" "30" "60" "120"} into f;
	close fevftrans;
	
	mn = nrow(m);
	p = ncol(m);
	mbar = m[:,]`;
	Sm = (m`*m - mn*mbar*mbar`)/(mn-1);

	
	fn = nrow(f);
	p = ncol(f);
	fbar = f[:,]`;
	Sf = (f`*f - fn*fbar*fbar`)/(fn-1);
	
	C = {1 -1 0 0 0,
	     0 1 -1 0 0,
	     0 0 1 -1 0,
	     0 0 0 1 -1};
	S = ((mn-1)*Sm + (fn-1)*Sf)/(mn+fn-2);
	tsq = mn*fn/(mn+fn)*(C*(mbar-fbar))`*inv(C*S*C`)*(C*(mbar-fbar));
	
	q = nrow(C);
	alpha = 0.05;
	cv = (mn+fn-2)*q/(mn+fn-q-1)*finv(1-alpha,q,mn+fn-q-1);
	pval = 1-probf((mn+fn-q-1)/(q*(mn+fn-2))*tsq,q,mn+fn-q-1);
	print tsq[format=10.4] cv[format=10.4] pval;
run;