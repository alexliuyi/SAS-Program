
/************************************************/
/******** Marginal Structure Model **************/
/************************************************/

%let tdvar= f_invisit f_ervisit; /* Time-Varing Variables */

%let fxvar= indexage male female 
	    b_Obesity b_Diabetes b_DrugAbuse; /* Fixed Variables */
	
%let outvar= f_invisit

data all_bipolar; 
set n.all_bipolar; 
	vis_p=vis+1; 
run;

proc sort data=all; by enrolid vis; run;

data MSM;
merge all(in=a keep=enrolid vis   cohort &outvar &fxvar where=(vis>1))
      all(in=b keep=enrolid vis_p cohort &tdvar  rename=(vis_p=vis cohort=cohort_p f_invisit=f_invisit_p f_ervisit=f_ervisit_p));
by enrolid vis;
if a and b;
run;



/* MSM Model */

/* Step 1. Computing Treatment Selction Weights */
/* treatment selection weights:  numerator calculation (probability of treatment with baseline covariates) */
proc logistic data=MSM;
	class cohort_p (ref='Drug01') ;
	model cohort = cohort_p &fxvar /link=glogit;
	output out=predtrt0(where=(cohort=_level_)) pred=predtrt0;
run;
     
/* treatment selection weights:  denominator calculation (probability of treatment with baseline covariates and time-dependent covariates) */
proc logistic data=MSM;
	class cohort_p (ref='Drug01') ;
	model cohort = cohort_p &fxvar &tdvar f_invisit_p f_ervisit_p /link=glogit;
	output out=predtrt1(where=(cohort=_level_)) pred=predtrt1;
run;
     




/* Step 2. Computing Stablized IPTW */ 
/* Observations where the time-varying factors are strong predictors of the current treatment selection are down-weighted in the analyes, 
	b/c such observations are over-represented in the observed data*/
proc sql;
/*ratio of probabilities for treatment*/
create table WEIGHT as
select distinct *, predtrt0/predtrt1 as predtrt
from predtrt1(keep=enrolid vis predtrt1)
	 natural full join
	 predtrt0(keep=enrolid vis predtrt0)
order enrolid, vis;
quit;

data MSM_weight;
merge MSM(in=a) Weight(in=b keep=enrolid vis predtrt);
by enrolid vis;
if a and b;
if first.enrolid then stabwt=predtrt;
					else stabwt=predtrt*dum;
retain dum;
dum=stabwt;
keep enrolid vis cohort cohort_p &fxvar &outvar stabwt f_medical_cost f_mh_medical_cost f_incost_new f_MH_incost_new;
run;


/* Step 3. Final Model: Weighted repeated measures analysis */

%macro Datacombine (outdata);

	proc sort data=LSMeans out=MS ; by cohort_p; run;
	proc sort data=Estimate out=ES ( where=(level1^=' ')); by level1; run;

	data OUTPUT;
	retain cohort_p MU STDERRMU ProbZ;
	merge ES(in=a keep=level1 probz rename=(level1=cohort_p)) MS(in=b keep=cohort_p MU STDERRMU);
	by cohort_p;
	if a and b;
	run;

	
data final;
retain variable;
merge 	OUTPUT(where=(cohort_p='Drug01') 		rename=(MU=MU0 STDERRMU=STDERRMU0 ProbZ=ProbZ0)) 
		OUTPUT(where=(cohort_p='Drug02') 	rename=(MU=MU1 STDERRMU=STDERRMU1 ProbZ=ProbZ1)) 
		OUTPUT(where=(cohort_p='Drug03') 	rename=(MU=MU2 STDERRMU=STDERRMU2 ProbZ=ProbZ2)) 
		OUTPUT(where=(cohort_p='Drug04') 	rename=(MU=MU3 STDERRMU=STDERRMU3 ProbZ=ProbZ3)) 
		OUTPUT(where=(cohort_p='Drug05') 	rename=(MU=MU4 STDERRMU=STDERRMU4 ProbZ=ProbZ4)) 
		OUTPUT(where=(cohort_p='Drug06') 	rename=(MU=MU5 STDERRMU=STDERRMU5 ProbZ=ProbZ5)) 
		OUTPUT(where=(cohort_p='Drug07') 	rename=(MU=MU6 STDERRMU=STDERRMU6 ProbZ=ProbZ6)) 
		;
drop cohort_p;
length variable $50;
variable="&outcome";
DROP ProbZ0;
run;

proc append data=final base=&outdata force;
run;

proc export data=&outdata outfile="E:\..\final_method.csv" dbms=csv replace; run;

%mend;


/* Health Care Utlization: Any visit */
ods select none;

proc datasets; delete final_binary;  quit;
%macro logitmodel (outcome);
	
	ods output  GEEEmpPEst=Estimate;
	ods output 	LSMeans=LSMeans;
	proc genmod data=MSM_weight descending;
	class  enrolid cohort_p(ref='0Lurasidone');
	weight stabwt;
	model &outcome = cohort_p &fxvar 	  					 
		  /  dist=bin link=logit type3;
	repeated subject = enrolid/type=ar(1) ;
	lsmeans cohort_p/ilink ;
	run; 
	ods output close;

	%datacombine(outdata=final_binary);

%mend;
%logitmodel(outcome=f_invisit);


/* Health Care Utlization: # of visit */
proc datasets; delete final_count;  quit;
%macro ngmodel (outcome=);

	ods output  GEEEmpPEst=Estimate;
	ods output 	LSMeans=LSMeans;
	proc genmod data=MSM_weight;
	class  enrolid cohort_p(ref='Drug01');
	weight stabwt;
	model &outcome = cohort_p &fxvar 					 
		  /  dist=negbin type3;
	repeated subject = enrolid/type=ar(1) ;
	lsmeans cohort_p/ilink ;
	run; 
	ods output close;

	%datacombine(outdata=final_count);

%mend;
%ngmodel(outcome=f_los  );




/* Health Care Cost */
proc datasets; delete final_cost;  quit;
%macro gamodel (outcome=);

	ods output  GEEEmpPEst=Estimate;
	ods output 	LSMeans=LSMeans;
	proc genmod data=MSM_weight;
		a=_mean_;
		b=_resp_;
		d=b/a+log(a);
		variance var=a**2;
		deviance dev=d;
	class  enrolid cohort_p(ref='Drug01');
	weight stabwt;
	model &outcome = cohort_p &fxvar			  					 
		  /  link=log type3;
	repeated subject = enrolid/type=ar(1) ;
	lsmeans cohort_p/ilink ;
	run; 
	ods output close;

	%datacombine(outdata=final_cost);

%mend;

%gamodel(outcome=f_incost );

