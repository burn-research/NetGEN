#include "udf.h"
#include "prop.h"
#include "pdf_transport.h"
#include "flamelet.h"
#include "math.h"

/* This UDF will reproduce the steady PASR model:

* To use from Fluent:
 * 1) set up Fluent just like EDC model (import Chemkin mechanism)
 * 2) build and load the UDF library
 * 3) Set the UDF in the menu models -> species models -> user defined EDC scales
 * 4) Allocate 10 UDM

 Follow further instructions in the end-user input section


/* ------------------------- start user input ------------------------- */
#define C1          2.1377
#define C2          0.4083
#define EDC_ACC     0    /* acceleration factor: 0=slow+stable -> 1=fast+unstable */
#define VERBOSITY   0

/*Definition of Coefficient from Case 1 Z. Li et al. for Dynamic PaSR*/
#define Cp1JM 1.7     
#define Cp2JM 1.45
#define Cd1JM 1.0
#define Cd2JM 0.9    


#define SMALL2 1e-8
#define SMALL3 1e-3

/* -------------------------  end user input  ------------------------- */

/* Define counter to determine C_mix, Ferrarotti et al. 2018 PROCI*/
 int counter=1;    /* 1 for static Cmix, 2 for fractal, 3 for dynamic JM*/

/*Define parameters*/
real Cmix=0.5;  /*To be modified if counter=1, so if the static Cmix approach is chosen*/
real D=4;       /*To be modified if counter=2, so if the static fractal approach is chosen*/

/*Define species to calculate tau chem as in EDC*/
#define num_species 4					/* Define number of species to be used for tau c calculation */
char *strs[num_species]={"ch4", "o2", "co2","h2o"}; 	/* definition of species that will be used for tau chemistry calculation*/
 

/* ---------------- Dynamic PaSR usage instructions -------------------- */

/* Enable 3 UDS: UDS-0 is mean mixture fraction, UDS-1 is mixture fraction variance, UDS-2 is scalar dissipation*/
/* Activate the source term for UDS1 (mixture craction variance) and UDS2 (chi) */
/* Set the UDS diffusivity to user-defined

/*Boundary conditions for UDS are:
 For UDS-0 (which is mean mixture fraction): 1 at fuel inlet and 0 at oxidizer inlet, 0 flux at walls.
 For UDS-1 (which is variance of mixture fraction) :  0 value at inlets and 0 flux at walls
 For UDS-2 (which is scalar dissipation of mixture fraction) :  0 value at inlets and 0 flux at walls*/

/* Procedure to follow in Fluent :
 /* 1) Deactivate all the equations except from UDS-0, UDS-1, UDS-2
 /* 2) Patch UDS-2 as 0.1 in the fluid domain
 /* 3) Run the simulations until residuals for all the UDS are below 1e-4
 /* 4) At this point it is possible to activate all the equations again */


 /* Useful parameters are stored in:
 /* UDM-1 = Minimum timescale
    UDM-2 = Chemical timescale
    UDM-3 = Mixing timescale
    UDM-4 = Mass fraction of the cell reactive zone
    UDM-5 = Mixture fraction variance source term
    UDM-6 = Mixture fraction dissipation rate source term
    UDM-7 = Cmix in dynamic formulation
    UDM-8 = Integral timescale
    UDM-9 = Kolmogorov timescale */

/* --------------------- end of the instructions -------------------- */

/* UDS Diffusivity */
DEFINE_DIFFUSIVITY(uds_diffusivity, c, t, i)
{
  return (C_K_L(c,t)/C_CP(c,t)+C_MU_T(c,t));
}

/*source term for mixture fraction variance*/
DEFINE_SOURCE(uds_variance_source,c,t,dS,eqn)/*attach to the UDS-1 source term in fluid cell zone*/
{
  real prod, diss, DV[ND_ND], source;
  NV_V(DV, =, C_UDSI_G(c,t,0));
  prod = 2.0 * C_MU_T(c,t) / M_species_sct * NV_DOT(DV, DV);
  diss = C_R(c,t)*C_UDSI(c,t,2);
  C_UDMI(c,t,5) = prod-diss;
  dS[eqn] = 0;
  return (prod - diss);
}

/*source term for scalar dissipation*/
DEFINE_SOURCE(uds_chi_source,c,t,dS,eqn)/*attach to the UDS-2 source term in fluid cell zone*/
{
  real prod, diss, DV[ND_ND];
  real Diss1=0.0;
  real Diss2=0.0;
  real Prod1=0.0;
  real Prod2=0.0;
  real Pf=0.0;
  real Pk=0.0;
  real Rt=0.0;
  real iRt=0.0;
  real var_diss=0.0;
  real alpha_s = 1;
  real source_chi=0.0;
  real source_limit =1e12;

  
  NV_V(DV, =, C_UDSI_G(c,t,0));
  
 
        var_diss=C_UDSI(c,t,1)/MAX(C_UDSI(c,t,2),SMALL3);
        iRt=var_diss/MAX(TRB_TIM_SCAL(c,t),SMALL3);
        Pk=C_MU_T(c,t)*C_STRAIN_RATE_MAG(c,t)*C_STRAIN_RATE_MAG(c,t); /*mu_t*S^2 */
        Pf=2*C_MU_T(c,t)/M_species_sct*NV_DOT(DV, DV); /*2*rho*D_t*gradZ^2 */
       
        Prod1 = Cp1JM/MAX(C_UDSI(c,t,1),SMALL3)*Pf*iRt; /*Cp1*Pf*chi/var*Rt^-1 */
        Prod2 = Cp2JM*Pk/MAX(C_K(c,t),SMALL3);  /* Cp2*Pk*Chi/k */
        Diss1 = Cd1JM*C_R(c,t)/MAX(C_UDSI(c,t,1),SMALL3);/* Cd1*rho*Chi^2/Zvar */
        Diss2 = Cd2JM*C_R(c,t)/MAX(TRB_TIM_SCAL(c,t),SMALL3); /*Cd2*rho*epsilon/k*Chi */ 

        prod = Prod1*C_UDSI(c,t,2) + Prod2*C_UDSI(c,t,2);
        diss = Diss1*C_UDSI(c,t,2)*C_UDSI(c,t,2)+ Diss2*C_UDSI(c,t,2);
        dS[eqn] =(Prod1 + Prod2) - (2*Diss1*C_UDSI(c,t,2) + Diss2);
        source_chi = alpha_s * (prod - diss);

	/* Limit the source term if necessary for numerical stability reasons */
        if (source_chi > source_limit)
        {
        source_chi = source_limit;
        } 
        if (source_chi < -source_limit)
        {
        source_chi = -source_limit;
        }
        C_UDMI(c,t,6) = source_chi;
   

  return (source_chi);
}




/* calculate tau_c*/
#define max_time_scale 100

real Cell_Chemical_Time_Scale_2016(cell_t c, Thread *t, char **specie_names, int number_spe)
{

Material *m, *sp;

real temp, pressure, result, interm;

real dydt[MAX_PDF_SPECIES], y[MAX_PDF_SPECIES];


int ns,  ind;
  
m = THREAD_MATERIAL(t);
  
mixture_species_loop(m,sp,ns)
{
	y[ns] = C_YI(c,t,ns);
}

temp = C_T(c,t);

pressure = op_pres + C_P(c,t);

/*Message("specie index %i name %s\n",ind,MIXTURE_SPECIE_NAME(m,ind));*/


 compute_net_react_rate_helper(y, ABS_P(C_P(c,t),op_pres), C_T(c,t), dydt, NULL, NULL, -1);

/*Message("after isat call\n");*/


result = 0.0;
          
for (ns=0; ns<number_spe; ns++)
{
 
    ind = mixture_specie_index(m, specie_names[ns]);
	
	if (fabs(dydt[ind]) < 1e-08) /* avoid division by 0 */
        interm = 0.0; /* indicator that specie k has an infinite time scale, otherwise it will false the MAX results */
    else
        interm = (y[ind]/fabs(dydt[ind]));

    result = MAX(result, interm);
                
               
}
   
if (result == 0.0) /* in this case, interm has been equal to 0 for all the species, meaning that all the species have an infinite time scale */
{
    result = max_time_scale;
}
   
/* clip the maximum value of result to max_time_scale, this is not mandatory... */
result = MIN(result,max_time_scale);



/*Message("result %g\n \n",result);*/
	
return result;

}


DEFINE_EDC_MDOT(edc_mdot, c, t, mdot, calc_tau, tau)
 {
  real ted = 0.0;
  real nu, tau_i, tau_k;
  real tau_m;
  real tau_c;
  real gamma;
  
  if (calc_tau)     
  { 

    if( M_turb_model==K_OMEGA || 
	M_turb_model==K_OMEGA_EASM ||
	M_turb_model==WJ_BSL_EARSM ||
	M_turb_model==SST || 
	M_turb_model==TRANS_SST || 
	M_turb_model==DES_SST || 
	M_turb_model==SAS_SST || 
	M_turb_model==RSM_K_OMEGA)
    	{
      	   real Cmu = (rp_kw_easm || rp_kw_wj_bsl_earsm) ? M_kw_easm_Cmu : M_kw_beta_star_inf;
      	   ted = MAX( Cmu*C_K(c,t)*C_O(c,t), 1.e-3);
	}

    else
    ted = MAX( C_D(c,t), 1.e-3);
    
    nu=C_MU_L(c,t)/C_R(c,t);
    tau_k = pow(nu/ted, 0.5); 	/*Kolmogorov time-scale*/
    tau_i = C_K(c,t)/ted;       /*Integral time-scale*/
    C_UDMI(c,t,8)=Cmix*tau_i;
    C_UDMI(c,t,9)=tau_k;
    
    if(counter==1)
      {
     /*static cmix*/ 
      tau_m=Cmix*tau_i;
    
      }
     else if (counter==2)
      {
      /*fractal*/
      real Ret, k2, alpha, alpha_exp, Cmix_eq_fractal;
      Ret=C_MU_T(c,t)/C_MU_L(c,t);
      k2=pow(C_K(c,t),2);
      alpha=(3*(D-3))/(1+D);
      alpha_exp=(1-alpha)/2;
      Cmix_eq_fractal=pow((ted*C_MU_L(c,t)/k2),alpha_exp);
      tau_m=tau_i*pow((ted*C_MU_L(c,t)/k2),alpha_exp);
      }
     else 
      	{
      /*dynamic*/
      real Cmix_eq_dynamic;
      tau_m=MAX(C_UDSI(c,t,1),SMALL2)/MAX(C_UDSI(c,t,2),SMALL2);
      	if(tau_m>tau_i)
      		{tau_m=tau_i;}
     	 if(tau_m>0.1)
		{tau_m=0.1;}
      	if(tau_m<tau_k)
		{
     		tau_m=tau_k;
      		}
	Cmix_eq_dynamic=tau_m/MAX(tau_i,SMALL2);
      	C_UDMI(c,t,7)=Cmix_eq_dynamic;
	}
  
      tau_c=Cell_Chemical_Time_Scale_2016(c, t, strs, num_species);
      C_UDMI(c,t,2)=tau_c;
      C_UDMI(c,t,1)=MIN(tau_m,tau_c);
      C_UDMI(c,t,3)=tau_m;
      
      *tau = tau_m;
  }
      


 gamma = C_UDMI(c,t,2)/(C_UDMI(c,t,2)+C_UDMI(c,t,3));
 C_UDMI(c,t,4)=gamma;
 
 *mdot =gamma;
 
 }
