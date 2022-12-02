/***********************************************************************
   UDF for the:
Extraction of mass flowrate across neighbor cells
************************************************************************/

/*include file needed to attach the udf to fluent */
#include "udf.h"

/*include files for writing output text files */
#include <stdlib.h>
#include <stdio.h>
#include "udf.h"
#include "mem.h"
#include "stdio.h"
#include "sg.h"
#include "sg_pdf.h"

DEFINE_ON_DEMAND(calc_mdot)
{
  Domain *domain;             	/* Initialize domain 		*/
  domain = Get_Domain(1);     	/* Get domain 			*/

  face_t f;			/* Initialize face index 	*/		
  Thread *ft;			/* Initialize face thread 	*/
  real xc_0[ND_ND];		/* Initialize centroid vector and velocity vector */
  real Af[ND_ND];		/* Face area vector 		*/
  real A;			/* Face area value */

  double mi, qi;		/* Initialize mass flowrate and heat flux */

  real u, v, h, rho;		/* Initialize quantities for heat flux */
  real vel[ND_ND];		/* Initialize velocity vector */

  /* Initialize file to print mass flowrates */	
  FILE *fp;
  fp = fopen("Neighbours_cell_flows","wb");
  fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\n", "x", "y", "z" ,"c_id0", "c_id1", "mass_flowrate");

  /* Initialize for exporting single cell indexes */
  cell_t c;
  Thread *ct;
  int id, zone_id;

  /* Initialize for exporting mass flowrates loop */
  int c0_id, c1_id;  	 	/* Initialize cell indexes */
  cell_t c0, c1;		/* Initialize cell indexes */
  Thread *c0_t, *c1_t;		/* Initialize cell thread  */
  
  /* Initialize file to print cell indexes and centroids */
  FILE *fid;
  fid = fopen("Cell_indexes", "wb");
  fprintf(fid, "%s\t%s\t%s\t%s\n", "ID", "x", "y", "z");

  /* Initialize a file to print the boundary cells and their associated boundary with face area */
  FILE *fib;
  fib = fopen("Boundary_cells", "wb");
  fprintf(fid, "%s\t%s\t%s\n", "ID", "BC", "Flowrate");
  
  /* Loop all over the cell in the domain */
  thread_loop_c(ct, domain)
  {
   	begin_c_loop(c, ct)
	{
		id = C_ID(c, ct);
		C_CENTROID(xc_0, c, ct);
		fprintf(fid, "%d\t%f\t%f\t%f\n", id, xc_0[0], xc_0[1], xc_0[2]);
	}
	end_c_loop(c, ct)
   }
   fclose(fid);	
  

  /* Begin a loop over the faces to get neighbors and mass flowrates */
  thread_loop_f(ft, domain)  /* loop over the faces in the domain */
  {
	begin_f_loop(f, ft)
	{
		c0 = F_C0(f, ft);		/* Get the internal cell */  	
		c0_t = F_C0_THREAD(f, ft);	/* Get the internal cell thread */
		c0_id = C_ID(c0, c0_t);		/* Get the internal cell index */

		c1_t = F_C1_THREAD(f, ft);	/* get the external cell */

		C_CENTROID(xc_0, c0, c0_t);	/* Get the cell centroid coordinates */

		if(c1_t == NULL)
		{
			c1_id = -1;	/* C1 does not exist since C0 is on a boundary face */
			mi = 0.0;	
		}
		else
		{
			c1 = F_C1(f, ft);		/* Get neighbor cell index */
			c1_id = C_ID(c1, c1_t);		/* Get neighbor cell index */

			/* Identify the geometry of the problem */
			if(ND_ND == 2)
			{
				mi = F_FLUX(f, ft)*2*M_PI; 	/* get the flux across the face */
			}
			else if(ND_ND == 3)
			{
				mi = F_FLUX(f, ft);
			}
		}	
		fprintf(fp, "%f\t%f\t%f\t%d\t%d\t%E\n", xc_0[0], xc_0[1], xc_0[2], c0_id, c1_id, mi);
	}
	end_f_loop(f,ft)
  }
fclose(fp);

  /* Begin a loop over the faces to print boundary cells and their associated BC's index */
  thread_loop_f(ft, domain)
  {
  	begin_f_loop(f, ft)
	{
		if(BOUNDARY_FACE_THREAD_P(ft))
		{ 
		zone_id = THREAD_ID(ft);				/* If it's a boundary cell print the id */
		c0 = F_C0(f, ft);					/* Get cell index */
		c0_t = F_C0_THREAD(f, ft);				/* Get cell thread */
		c0_id = C_ID(c0, c0_t);					/* Get cell index */
		F_AREA(Af, f, ft);					/* Get face area vector */
		A = NV_DOT(Af, Af);					/* Get face area value */

		u = C_U(c0, c0_t);
		v = C_V(c0, c0_t);
		h = C_H(c0, c0_t);
		rho = C_R(c0, c0_t);

		vel[0] = u;
		vel[1] = v;

		/* Identify the geometry of the problem */
		if(ND_ND == 2)
		{
			mi = F_FLUX(f, ft)*2*M_PI;				/* Get the flux across the cell */
			qi = rho*h*NV_DOT(vel, Af);				/* Get heat flux across the cell */
		}
		else if(ND_ND == 3)
		{
			mi = F_FLUX(f, ft);
			qi = 0;
		}
		fprintf(fib, "%d\t%d\t%E\t%E\t%E\n", c0_id, zone_id, mi, A, qi);	/* Print the information */
		}
	}
	end_f_loop(f, ft)
  }
  fclose(fib);
}





