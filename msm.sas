
/************************************************/
/******** Marginal Structure Model **************/
/************************************************/

%let tdvar=	 f_invisit f_ervisit; /* Time-Varing Variables */

%let fxvar= indexage male female 
			b_Pervasive b_Disruptive b_Attention b_Depression b_Anxiety b_Adjustment b_Obesity b_Diabetes b_DrugAbuse
			comed_Antidepressants comed_Anxiolytics comed_Other; /* Fixed Variables */
	
%let outvar= f_los f_invisit f_ervisit f_officevisit f_outotvisit f_pharmvisit f_out_anyvisit
	 		 f_incost--f_pharmcost f_n_invisit f_n_ervisit--f_n_pharmvisit f_n_out_anyvisit f_out_total_cost f_total_cost f_incost_new
			 f_MH_los f_MH_invisit f_MH_ervisit f_MH_officevisit f_MH_outotvisit f_MH_pharmvisit f_MH_out_anyvisit
	 		 f_MH_incost--f_MH_pharmcost f_MH_n_invisit f_MH_n_ervisit--f_MH_n_pharmvisit f_MH_n_out_anyvisit f_MH_out_total_cost f_MH_total_cost f_MH_incost_new;  /* Outcomes*/


data all_bipolar; 
set n.all_bipolar; 
	vis_p=vis+1; 
if f_invisit=1 then f_incost_new = f_incost;
if f_MH_invisit=1 then f_MH_incost_new = f_MH_incost;
run;

proc sort data=all_bipolar; by enrolid vis; run;

data MSM;
merge all_bipolar(in=a keep=enrolid vis   cohort &outvar &fxvar where=(vis>1))
      all_bipolar(in=b keep=enrolid vis_p cohort &tdvar  rename=(vis_p=vis cohort=cohort_p f_invisit=f_invisit_p f_ervisit=f_ervisit_p));
by enrolid vis;
if a and b;

f_medical_cost = f_out_total_cost + f_incost;
f_MH_medical_cost = f_MH_out_total_cost + f_MH_incost;
run;



/* MSM Model */

/* Step 1. Computing Treatment Selction Weights */
/* treatment selection weights:  numerator calculation (probability of treatment with baseline covariates) */
proc logistic data=MSM;
	class cohort_p (ref='0Lurasidone') ;
	model cohort = cohort_p &fxvar /link=glogit;
	output out=predtrt0(where=(cohort=_level_)) pred=predtrt0;
run;
     
/* treatment selection weights:  denominator calculation (probability of treatment with baseline covariates and time-dependent covariates) */
proc logistic data=MSM;
	class cohort_p (ref='0Lurasidone') ;
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
merge 	OUTPUT(where=(cohort_p='0Lurasidone') 		rename=(MU=MU0 STDERRMU=STDERRMU0 ProbZ=ProbZ0)) 
		OUTPUT(where=(cohort_p='1Quetiapine') 	rename=(MU=MU1 STDERRMU=STDERRMU1 ProbZ=ProbZ1)) 
		OUTPUT(where=(cohort_p='2Risperidone') 	rename=(MU=MU2 STDERRMU=STDERRMU2 ProbZ=ProbZ2)) 
		OUTPUT(where=(cohort_p='3Aripiprazole') rename=(MU=MU3 STDERRMU=STDERRMU3 ProbZ=ProbZ3)) 
		OUTPUT(where=(cohort_p='4Olanzapine') 	rename=(MU=MU4 STDERRMU=STDERRMU4 ProbZ=ProbZ4)) 
		OUTPUT(where=(cohort_p='5No/Mini') 	rename=(MU=MU5 STDERRMU=STDERRMU5 ProbZ=ProbZ5)) 
		OUTPUT(where=(cohort_p='6Other_Treat') 	rename=(MU=MU6 STDERRMU=STDERRMU6 ProbZ=ProbZ6)) 
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
	class  enrolid cohort_p(ref='0Lurasidone');
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
	class  enrolid cohort_p(ref='0Lurasidone');
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






/* -------------------------Main Part End Here------------------------------------------------------- */


/* Health Care Cost without No/Min group*/
proc datasets; delete final_cost;  quit;
%macro gamodel (outcome=);

	ods output  GEEEmpPEst=Estimate;
	ods output 	LSMeans=LSMeans;
	proc genmod data=MSM_weight;
	where cohort_p not in ('5No/Mini'); 
		a=_mean_;
		b=_resp_;
		d=b/a+log(a);
		variance var=a**2;
		deviance dev=d;
	class  enrolid cohort_p(ref='0Lurasidone');
	weight stabwt;
	model &outcome = cohort_p &fxvar			  					 
		  /  link=log type3;
	repeated subject = enrolid/type=ar(1) ;
	lsmeans cohort_p/ilink ;
	run; 
	ods output close;

	%datacombine(outdata=final_cost);

%mend;

%gamodel(outcome=f_outotcost);				 /* Error: Error in computing the variance function.*/
%gamodel(outcome=f_medical_cost);			 /* Error: Error in computing the variance function.*/
%gamodel(outcome=f_total_cost);				 /* Error: Error in computing the variance function.*/
%gamodel(outcome=f_MH_outotcost);			 /* Error: Error in computing the variance function.*/
%gamodel(outcome=f_MH_medical_cost);		 /* ERROR: Negative variance from user-defined variance function.
												ERROR: Error in computing the variance function.*/
%gamodel(outcome=f_MH_total_cost);			 /* ERROR: Termination due to Floating Point Exception*/



/* GLM for INCOST and ERCOST */
%let fxvar= indexage male female
			b_Pervasive b_Disruptive b_Attention b_Depression b_Anxiety b_Adjustment b_Obesity b_Diabetes b_DrugAbuse
			comed_Antidepressants comed_Anxiolytics comed_Other;

proc summary data=All_bipolar nway;
class enrolid &fxvar;
var f_incost f_out_total_cost f_ercost f_officecost f_outotcost  f_pharmcost f_total_cost
	f_MH_incost f_MH_out_total_cost f_MH_ercost f_MH_officecost f_MH_outotcost  f_MH_pharmcost f_MH_total_cost;
output out=alex(drop=_:) sum=;
run;

proc means data=alex noprint;
var &fxvar;
output out=meanvar(drop=_type_ _freq_) mean=;
run;

data alex_sum;
merge alex
	  all_bipolar(keep=enrolid cohort vis where=(vis=1));
by enrolid;
f_medical_cost = f_incost + f_out_total_cost;
f_MH_medical_cost = f_MH_incost + f_MH_out_total_cost;
run;

data alex_final;
set alex_sum meanvar end=k;
if cohort='0Lurasidone'   then gp = 0;
if cohort='1Quetiapine'   then gp = 1;
if cohort='2Risperidone'  then gp = 2;
if cohort='3Aripiprazole' then gp = 3;
if cohort='4Olanzapine'   then gp = 4;
if cohort='5No/Mini'      then gp = 5;
if cohort='6Other_Treat'  then gp = 6;
output;
if k then do; gp=0; cohort='0Lurasidone';   output; end;
if k then do; gp=1; cohort='1Quetiapine';   output; end;
if k then do; gp=2; cohort='2Risperidone';  output; end;
if k then do; gp=3; cohort='3Aripiprazole'; output; end;
if k then do; gp=4; cohort='4Olanzapine';   output; end;
if k then do; gp=5; cohort='5No/Mini'; 		output; end;
if k then do; gp=6; cohort='6Other_Treat';  output; end;
run;


%macro glm(dep=);
proc genmod data=alex_final;
a=_mean_;
b=_resp_;
d=b/a+log(a);
variance var=a**2;
deviance dev=d;
class  enrolid cohort(ref='0Lurasidone');
model &dep = &fxvar cohort/link=log;
output out=est(keep=cohort &dep.2 enrolid where=(enrolid=. and cohort^=' ')) p=&dep.2;
ods output parameterestimates = odds(keep=level1 ProbChiSq rename=(level1=cohort ProbChiSq=&dep.) where=(cohort^=''));
run; 

proc sort data=est(rename=(&dep.2=&dep));
by cohort;
run;
proc sort data=odds;
by cohort;
run;

proc transpose data=est out=est1;
id cohort;
var &dep;
run;

proc transpose data=odds out=odds1;
id cohort;
var &dep;
run;

data result;
retain _name_ _0Lurasidone _1Quetiapine _1p _2Risperidone _2p _3Aripiprazole _3p _4Olanzapine _4p _5No_Mini _5p _6Other_Treat _6p;
merge est1(in=a drop=_label_)
	  odds1(drop=_label_ _0Lurasidone rename=(_1Quetiapine=_1p _2Risperidone=_2p _3Aripiprazole=_3p _4Olanzapine=_4p _5No_Mini=_5p _6Other_Treat=_6p));
by _name_;
retain _name_ _0Lurasidone _1Quetiapine _1p _2Risperidone _2p _3Aripiprazole _3p _4Olanzapine _4p _5No_Mini _5p _6Other_Treat _6p;
run;

proc append base=table_all data=result force;
run;
%mend;

ods select none;

%glm(dep=f_incost); 
%glm(dep=f_out_total_cost); 
%glm(dep=f_ercost); 
%glm(dep=f_officecost); 
%glm(dep=f_outotcost); 
%glm(dep=f_medical_cost); 
%glm(dep=f_pharmcost); 
%glm(dep=f_total_cost); 

%glm(dep=f_MH_incost); 
%glm(dep=f_MH_out_total_cost); 
%glm(dep=f_MH_ercost); 
%glm(dep=f_MH_officecost); 
%glm(dep=f_MH_outotcost); 
%glm(dep=f_MH_medical_cost); 
%glm(dep=f_MH_pharmcost); 
%glm(dep=f_MH_total_cost); 

proc export data=Table_all outfile="E:\Yigong Zhou\APs\RESULTS\MSM output\glm results.csv" dbms=csv replace; run;



%let fxvar= indexage male
			b_Pervasive b_Disruptive b_Attention b_Depression b_Anxiety b_Adjustment b_Obesity b_Diabetes b_DrugAbuse
			comed_Antidepressants comed_Anxiolytics comed_Other;

/* Odds Ratio for invisit */
ods select all;
ods output  GEEEmpPEst=Estimate;
ods output 	LSMeans=LSMeans;
proc genmod data=MSM_weight descending;
class  enrolid cohort_p(ref='0Lurasidone');
weight stabwt;
model f_invisit = cohort_p &fxvar			  					 
	  /  dist=bin link=logit type3;
repeated subject = enrolid/type=ar(1) ;

lsmeans cohort_p/ilink;
estimate "Odds for age"    indexage 1/exp;
estimate "Odds for male"   male 1/exp;
estimate "Odds for L vs Que" cohort_p 1 0 0 0 0 0 -1/exp;
estimate "Odds for L vs Ris" cohort_p 0 1 0 0 0 0 -1/exp;
estimate "Odds for L vs Ari" cohort_p 0 0 1 0 0 0 -1/exp;
estimate "Odds for L vs Ola" cohort_p 0 0 0 1 0 0 -1/exp;
estimate "Odds for L vs Min" cohort_p 0 0 0 0 1 0 -1/exp;
estimate "Odds for L vs Oth" cohort_p 0 0 0 0 0 1 -1/exp;
estimate "Odds for b_Pervasive"   	b_Pervasive 1/exp;
estimate "Odds for b_Disruptive"    b_Disruptive 1/exp;
estimate "Odds for b_Attention"   	b_Attention 1/exp;
estimate "Odds for b_Depression"  	b_Depression 1/exp;
estimate "Odds for b_Anxiety"   	b_Anxiety 1/exp;
estimate "Odds for b_Adjustment"    b_Adjustment 1/exp;
estimate "Odds for b_Obesity"   	b_Obesity 1/exp;
estimate "Odds for b_Diabetes"   	b_Diabetes 1/exp;
estimate "Odds for b_DrugAbuse"   	b_DrugAbuse 1/exp;
estimate "Odds for comed_Antidepressants"   comed_Antidepressants 1/exp;
estimate "Odds for comed_Anxiolytics" 	  	comed_Anxiolytics 1/exp;
estimate "Odds for comed_Other"   			comed_Other 1/exp;
run; 
