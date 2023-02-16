/* Logistic models of individual survival */

proc logistic data = RT_pop plots(only)=effectplot(clband);
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =   exp pgland pgland*exp;
run;

/* no block or exp*pgland effects */
proc logistic data = RT_pop plots(only)=effectplot(clband);
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =   exp pgland ;
run;


proc logistic data = RT_pop plots(only)=effectplot(clband);
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =  oid_moy  pgland  pgland*oid_moy;
run;

/* pas d'effet bloc */
proc logistic data = RT_pop plots(only)=effectplot(clband);
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =  oid_moy  pgland ;
run;


proc logistic data = RT_pop plots(only)=effectplot(clband);
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =  oid_moy  H09v ;
run;

proc logistic data = RT_pop plots(only)=effectplot(clband); /* exp for graphical output */
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =   H09v exp  ;
run;

proc logistic data = RT_pop plots(only)=effectplot(clband); /* exp for graphical output */
class  bloc exp (ref='exp');
model survie2017 (Event = '1') =   H09v exp  H09v*exp;
run;

proc logistic data = RT_pop plots(only)=effectplot(clband); 
class  bloc exp (ref='exp') gel_13 (ref='0');
model survie2017 (Event = '1') =  oid_moy  gel_13 H09v   ;
run;

/* with family effects */
proc logistic data = RT_pop ; 
class  bloc exp (ref='exp') gel_13(ref='0') fameff;
model survie2017 (Event = '1') =   bloc  exp fameff fameff*exp   pgland gel_13;
run;

/* with heterozygosity */

proc logistic data = RT_pop_tot plots(only)=effectplot(clband);
class  bloc exp (ref='exp')fameff;
model survie2017 (Event = '1') =   exp fameff pgland  gel_13 PHt PHt*exp;
units PHt = 0.1;
run;

proc logistic data = RT_pop_tot plots(only)=effectplot(clband);
class  bloc exp (ref='exp') fameff;
model survie2017 (Event = '1') =   exp fameff pgland  gel_13 PHt ;
units PHt = 0.1;
run;



/* Path analysis */

proc calis data=RT_pop; 
  path                      
   survie2017 <- H12v ,
   survie2017 <- oid_moy,
   survie2017 <- gel_13,
   H12v <- oid_moy,
   H12v <- pgland,
   gel_13 <- oid_moy,
   gel_13 <- H12v ;
effpart                   /* for direct and indirect effects */
  survie2017 <-  H12v oid_moy ,
  H12v <- pgland oid_moy ;
run;






/* Family effects */

proc sort data = RT_pop;
by exp;
run;

proc glm data = RT_pop; by exp;
class bloc exp fameff;
model oid_moy = bloc fameff/SS3;
run;

proc glm data = RT_pop; by exp;
class bloc exp fameff;
model H17v = bloc  fameff/SS3;
run;


proc glm data = RT_pop;
class bloc exp fameff;
model pgland = bloc exp fameff/SS3;
run;
proc glm data = RT_pop;
class bloc exp fameff;
model H12v = bloc exp fameff/SS3;
run;




