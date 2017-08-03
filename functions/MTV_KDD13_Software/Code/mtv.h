#include <math.h>
#include <time.h>			
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include "matrix.h"
#include <omp.h>
#define sqrt2 1.414213562373095
#define sqrt3 1.732050807568877
bool stopPri3D(double *X1, double *X2, double *X3, double *Z, int m, int n, int d, double absTol, double relTol, int display)
{
	double rk = 0, epri = 0, ex = 0, ez = 0;
	int i,j,l,k, l2, ind1, ind2, ind3;

	l = m*n;
	l2 = n*d;
	#pragma omp parallel for shared(X1,X2,X3,Z) private(i,j,k,ind1,ind2,ind3) reduction(+: rk,ex,ez)
    for(i = 0; i<n; i++)
		for(j=0;j<m;j++)
			for(k=0;k<d;k++){
					ind1 = k*l + i*m + j;
					ind2 = k*l+ j*n + i;
					ind3 = l2*j+d*i + k;
		            rk += (X1[ind1]-Z[ind1])*(X1[ind1]-Z[ind1]);
					rk += (X2[ind2]-Z[ind1])*(X2[ind2]-Z[ind1]);
					rk += (X3[ind3]-Z[ind1])*(X3[ind3]-Z[ind1]);

		            ex += X1[ind1]*X1[ind1] + X2[ind2]*X2[ind2] + X3[ind3]*X3[ind3];
		            ez += Z[ind1]*Z[ind1];
	        }
	rk = sqrt(rk);
	ez = sqrt3*sqrt(ez);
	ex = sqrt(ex);
	epri = ex>ez?ex:ez;
	i = m>n?m:n;
	epri *= relTol;
	epri += sqrt3*sqrt(1.0*m*n*d)*absTol;

	if(display >0)
	   printf("PriError and Resideual is %f %f\n", (float) epri, (float) rk);
	return (rk <= epri);
}

bool stopDual3D(double *Z, double *Zo, double *U1, double *U2, double *U3, int m, int n, int d, 
				double rho, double absTol, double relTol, int display)
{
	double sk = 0, ed = 0;
	int i;
	#pragma omp parallel for shared(U1,U2,U3,Z,Zo) private(i) reduction(+: sk,ed)
    for (i=0;i<m*n*d;i++)
	{
		sk += (Z[i]-Zo[i])*(Z[i]-Zo[i]);
		ed += U1[i]*U1[i] + U2[i]*U2[i] + U3[i]*U3[i];
	}
	sk =rho*sqrt3*sqrt(sk);
	ed = sqrt(ed);
    ed *= relTol;
	i = m>n?m:n;
	ed += sqrt3*sqrt(1.0*m*n*d)*absTol;
	if (display > 0)
	   printf("DualError and Resideual is %f %f\n", (float) ed, (float) sk);
	return (sk <= ed);
}

bool stopPri(double *X, double *X1, double *Z, int m, int n, double absTol, double relTol,int display)
{
	double rk = 0, epri = 0, ex = 0, ez = 0;
	int i,j,l,k;

	#pragma omp parallel for shared(X,X1,Z) private(i,j,l,k) reduction(+: rk,ex,ez)
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
	    {
			l = i*n+j;
			k = i+j*m;
		    rk += (X[l]-Z[l])*(X[l]-Z[l]) + (X1[k]-Z[l])*(X1[k]-Z[l]);
		    ex += X[l]*X[l] + X1[k]*X1[k];
		    ez += Z[l]*Z[l];
	    }
	rk = sqrt(rk);
	ez = sqrt2*sqrt(ez);
	ex = sqrt(ex);
	epri = ex>ez?ex:ez;
	i = m>n?m:n;
	epri *= relTol;
	epri += sqrt2*i*absTol;
	
	if (display > 0)
	   printf("PriError and Resideual is %f %f\n", (float) epri, (float) rk);
	return (rk <= epri);
}

bool stopDual(double *Z, double *Zo, double *U, double *U1, int m, int n, double rho, double absTol, double relTol, int display)
{
	double sk = 0, ed = 0;
	int i;
	#pragma omp parallel for shared(U,U1,Z,Zo) private(i) reduction(+: sk,ed)
    for (i=0;i<m*n;i++)
	{
		sk += (Z[i]-Zo[i])*(Z[i]-Zo[i]);
		ed += U[i]*U[i] + U1[i]*U1[i];
	}
	sk =rho*sqrt2*sqrt(sk);
	ed = sqrt(ed);
    ed *= relTol;
	i = m>n?m:n;
	ed += sqrt2*i*absTol;
	if (display > 0)
	   printf("DualError and Resideual is %f %f\n", (float) ed, (float) sk);
	return (sk < ed);
}
// The function TV1D_denoise is written by Dr.Laurent Condat, PhD, CNRS research fellow in France.
// For more details of this function, please refer to http://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/research.html
void TV1D_denoise(double* input, double* output, const int width, const double lambda) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		double umin=lambda, umax=-lambda;	/*u is the dual variable*/
		double vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const double twolambda=2.0*lambda;	/*auxiliary variable*/
		const double minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1); 
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}

			if(k>width-1)
				break;

			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
}


/*
 Function:  mtv
 Dec. 18, completely decompable
 Author:    Sen Yang, Arizona State University
 */
size_t mtv(double* X, const double* Y, const int m, const int n, double lam, double rho,
                     const int maxIter, double absTol, double relTol, const int display)
{
	double *X1, *U, *vtemp, *Z, *Zo, *U1,flam2,td, rhoinv,obj;
	int iter,i,j;

	X1 = (double *) malloc(m*n*sizeof(double));
	Z = (double *) malloc(m*n*sizeof(double));
	Zo = (double *) malloc(m*n*sizeof(double));
	U  = (double *) calloc(m*n, sizeof(double));
	U1  = (double *) calloc(m*n, sizeof(double));

	if (m>n)
	   vtemp = (double *) malloc(m*sizeof(double));
	else
	   vtemp = (double *) malloc(n*sizeof(double));

    rhoinv = 1/rho;

	for (i = 0; i < n; i++)
		TV1D_denoise(&X[m*i], &Z[m*i], m, lam);
	
	for (i = 0; i<m; i++)
	{
		for (j = 0; j<n; j++)
	        vtemp[j] = Z[i+j*m] ; 
		TV1D_denoise(vtemp, &X1[n*i], n, lam);
	}

	td = 1.0/(1 + 2* rho);
	for (i = 0; i<m; i++)
		for (j = 0; j<n; j++)
			Z[i+j*m] = (Y[i+j*m] + 2*rho*X1[i*n+j])*td;

	memcpy(Zo, Z, m*n*sizeof(double));

	for (iter = 0; iter < maxIter; iter++)
	{
		flam2 = lam * rhoinv;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j<m; j++)
			   vtemp[j] = Z[m*i+j] - U[m*i+j] * rhoinv;
			TV1D_denoise(vtemp, &X[m*i], m, flam2);
		}
		
		for (i = 0; i<m; i++)
		{
			for (j = 0; j<n; j++)
		        vtemp[j] = Z[i+j*m] - U1[i+j*m] * rhoinv; 
			TV1D_denoise(vtemp, &X1[n*i], n, flam2);
		}
		
		td = 1.0/(1 + 2* rho);
		for(i = 0; i<n; i++)
            for(j = 0; j<m;j++)
			Z[i*m+j] = (Y[i*m+j] + U[i*m+j] + U1[i*m+j] 
                        + rho*(X1[i+j*n] + X[i*m+j]))*td;

		for (i = 0; i<n; i++)
		   for (j=0; j<m;j++)
		   {
			   U[i*m + j] += rho * (X[i*m + j] -Z[i*m + j]);
			   U1[i*m + j] += rho * (X1[i + j*n] -Z[i*m + j]);
		   }

		if(!stopDual(Z,Zo,U,U1,(int)m,(int)n,rho,absTol,relTol,display))
		{
		    memcpy(Zo, Z, m*n*sizeof(double));
			continue;
        }		
		if(stopPri(X,X1,Z,(int)m,(int)n,absTol,relTol,display))
			break;
        memcpy(Zo, Z, m*n*sizeof(double));
	}
	
	for (i = 0; i<n; i++)
		for (j=0; j<m;j++)
        {
          X[i*m+j] += X1[i + j*n];
          X[i*m+j]/=2;
        }
	free(X1);
	free(vtemp);
	free(U);
	free(Z);
	free(U1);
	free(Zo);
	return (size_t) iter;

}

int mtvp(double* X, const double* Y, const int m, const int n, double lam2, double rho,
                      const int maxIter, const double absTol, const double relTol, const int display)
{
	double *X1, *U, *T1,*Z, *U1, *Zo,flam2, td, rhoinv;
	int iter;
    int i,j;

	X1 = (double *) malloc(m*n*sizeof(double));
	U  = (double *) calloc(m*n, sizeof(double));
    T1 = (double *) malloc(m*n*sizeof(double));
	Z = (double *) malloc(m*n*sizeof(double));
    Zo = (double *) malloc(m*n*sizeof(double));
	U1  = (double *) calloc(m*n, sizeof(double));

    rhoinv = 1/rho;

	//omp_set_num_threads(num_threads);
	//tempTime = omp_get_wtime();
	#pragma omp parallel shared(Z,X,T1,X1) private(i,j)
	{
	   //printf("\nTreads %d on CPU %d  ",omp_get_thread_num(),sched_getcpu());
	   #pragma omp for
	   for (i = 0; i < n; i++)
		   TV1D_denoise(&X[m*i], &Z[m*i], m, lam2);

	   #pragma omp for	
	   for (i = 0; i<m; i++)
		   for (j = 0; j<n; j++)
		    T1[n*i+j] = Z[i+j*m];

       #pragma omp for
       for (i = 0; i<m; i++)
	       TV1D_denoise(&T1[n*i], &X1[n*i], n, lam2);
	}
	//timeP += omp_get_wtime() - tempTime;


	td = 1.0/(1 + 2* rho);
	//tempTime = omp_get_wtime();
	#pragma omp parallel for shared(Z,Y,X1) private(i,j)
	for (i = 0; i<m; i++)
		for (j = 0; j<n; j++)
			Z[i+j*m] = (Y[i+j*m] + 2*rho*X1[i*n+j])*td;
	//timeP += omp_get_wtime() - tempTime;
	memcpy(Zo, Z, m*n*sizeof(double));
	flam2 = lam2 * rhoinv;
	
   
   for (iter = 0; iter < maxIter; iter++)
   {
	   //tempTime = omp_get_wtime();
	   //omp_set_num_threads(num_threads);
	   #pragma omp parallel shared(T1,X1,Y,U,U1,Z,X) private(i,j)
	   {
        #pragma omp for 
        for(i = 0; i<n; i++)
        	for (j = 0; j<m; j++)
			   T1[m*i+j] = Z[m*i+j] - U[m*i+j]*rhoinv;
        
        #pragma omp for
		for (i = 0; i < n; i++)
			TV1D_denoise(&T1[m*i], &X[m*i], m, flam2);
		
        #pragma omp for
		for (i = 0; i<m; i++)
			for (j = 0; j<n; j++)
		        T1[n*i+j] = Z[i+j*m] - U1[i+j*m]*rhoinv;
        
        #pragma omp for
        for (i = 0; i<m; i++)
			TV1D_denoise(&T1[n*i], &X1[n*i], n, flam2);
        #pragma omp for
		for(i = 0; i<n; i++)
            for(j = 0; j<m;j++)
			Z[i*m+j] = (Y[i*m+j] + U[i*m+j] + U1[i*m+j] 
                        + rho*(X1[i+j*n] + X[i*m+j]))*td;
		
        #pragma omp for
		for (i = 0; i<n; i++)
		   for (j=0; j<m;j++)
			   U[i*m + j] += rho * (X[i*m + j] - Z[i*m + j]);

        #pragma omp for
		for (i = 0; i<n; i++)
		   for (j=0; j<m;j++)
			   U1[i*m + j] += rho * (X1[i + j*n] - Z[i*m + j]);
        
	    }
	  // printf("iteration idex is %d\n",iter);
		if(!stopDual(Z,Zo,U,U1,(int)m,(int)n,rho,absTol,relTol,display))
		{
		    memcpy(Zo, Z, m*n*sizeof(double));
			continue;
        }		
		if(stopPri(X,X1,Z,(int)m,(int)n,absTol,relTol,display))
			break;
	   //timeP += omp_get_wtime() - tempTime;
        memcpy(Zo, Z, m*n*sizeof(double));
	}
    //tempTime = omp_get_wtime();
    #pragma omp parallel for shared(X,X1) private(i,j)
	for (i = 0; i<n; i++)
		for (j=0; j<m;j++)
	   {
          X[i*m+j] += X1[i + j*n];
          X[i*m+j]/=2;
        }
	//timeP += omp_get_wtime() - tempTime;
	//printf(" timeP is %f ",(float)timeP);
	free(X1);
	free(T1);
	free(U);
	free(Z);
	free(U1);
	free(Zo);
	return iter;
}

size_t mtvp3(double *X1, double *Y, int m, int n, int d,
                    double lam2, double rho, int maxIter, double absTol, double relTol, int display)
{
	double *X2, *X3, *U1, *U2, *U3, *Z, *Zo, *T, td, rhoinv,flam1,flam2;
	int iter, i,j, k, ind1,ind2,ind3;
    const int l = (int)(m*n);
	const int l2 = (int)(d*n);

	rhoinv = 1/rho;
	X2 = (double *) malloc(m*n*d*sizeof(double));
	X3 = (double *) malloc(m*n*d*sizeof(double));
	U1 = (double *) malloc(m*n*d*sizeof(double));
	U2 = (double *) malloc(m*n*d*sizeof(double));
	U3 = (double *) malloc(m*n*d*sizeof(double));
	Z  = (double *) malloc(m*n*d*sizeof(double));
	Zo  = (double *) malloc(m*n*d*sizeof(double));
	T  = (double *) malloc(m*n*d*sizeof(double));

	// compute intinial solutions
	memcpy(U1, Y, m*n*d*sizeof(double));
	td = 1.0/(1 + 3* rho);
	
	//omp_set_num_threads(num_threads);
	#pragma omp parallel shared(U2,Y,U1,U3,X1,X2,X3,Z) private(i,j,k)
	{
		//printf("\nTreads %d on CPU %d  ",omp_get_thread_num(),sched_getcpu());
		#pragma omp for
	    for(k=0; k<d; k++)
		for(j=0; j<m;j++)
			for(i=0;i<n;i++)
			U2[k*l+ j*n + i] = Y[k*l + i*m + j];		
	    #pragma omp for
	    for(k=0; k<d; k++)
		for(j=0; j<m;j++)
			for(i=0;i<n;i++)
			U3[l2*j+d*i + k] = Y[k*l+i*m+j];

	    #pragma omp for
	    for (j = 0; j < n*d; j++)
		    TV1D_denoise(&U1[m*j], &X1[m*j], m, 3*lam2);

	    #pragma omp for
	    for (i = 0; i < m*d; i++)
		    TV1D_denoise(&U2[n*i], &X2[n*i], n, 3*lam2);
			
	    #pragma omp for
	    for (k = 0; k < m*n; k++)
		    TV1D_denoise(&U3[d*k], &X3[d*k], d, 3*lam2);

	    #pragma omp for
	    for(i = 0; i<n; i++)
		for(j=0;j<m;j++)
			for(k=0;k<d;k++){
				X1[k*l + i*m + j] += X2[k*l+ j*n + i] + X3[l2*j+d*i + k];
				X1[k*l + i*m + j] /=3.;
				}

		#pragma omp for
		for(i = 0; i<n*m*d; i++)
			Z[i] = (Y[i] + 3*rho*X1[i])*td;
	}

	memset(U1,0,m*n*d*sizeof(double));
	memset(U2,0,m*n*d*sizeof(double));
	memset(U3,0,m*n*d*sizeof(double));
	memcpy(Zo,Z,m*n*d*sizeof(double));

	flam2 = lam2 * rhoinv;
	for(iter = 0; iter<maxIter; iter++)
	{
		// compute X1
		//omp_set_num_threads(num_threads);
    	#pragma omp parallel shared(Z,Y,U1,U2,U3,X1,X2,X3,T) private(i,j,k,ind1,ind2,ind3)
		{
			#pragma omp for
		    for(i=0;i<m*n*d;i++)
			    T[i] = Z[i] - U1[i] * rhoinv;

		    #pragma omp for
		    for(i=0;i<d*n;i++)
				TV1D_denoise(&T[m*i], &X1[m*i], (int)m, flam2);
        
		// compute X2
		    #pragma omp for
		    for(k=0; k<d; k++)
		       for(j=0; j<m;j++)
			   for(i=0;i<n;i++){
				ind1 = k*l+ j*n + i;
			    T[ind1] = Z[k*l + i*m + j] - U2[ind1]*rhoinv;
			}

            #pragma omp for
		    for(i=0;i<d*m;i++)
			   TV1D_denoise(&T[n*i], &X2[n*i], (int)n, flam2);
        
		// compute X3
		    #pragma omp for
		    for(i = 0; i<n; i++)
			for(j=0;j<m;j++)
				for(k=0;k<d;k++){
					ind1 = l2*j+d*i + k;
					T[ind1] = Z[k*l+i*m+j] - U3[ind1]*rhoinv;
				}
		    #pragma omp for
  		    for(i=0;i<n*m;i++)
			    TV1D_denoise(&T[d*i], &X3[d*i], d, flam2);

		    #pragma omp for
            for(i = 0; i<n; i++)
			for(j=0;j<m;j++)
				for(k=0;k<d;k++){
					ind1 = k*l + i*m + j;
					ind2 = k*l+ j*n + i;
					ind3 = l2*j+d*i + k;

		           Z[ind1] = (Y[ind1] + U1[ind1] + U2[ind2] + U3[ind3] 
                     + rho*(X1[ind1] + X2[ind2] + X3[ind3]))*td;

					U1[ind1] += rho * (X1[ind1] - Z[ind1]);
					U2[ind2] += rho * (X2[ind2] - Z[ind1]);
					U3[ind3] += rho * (X3[ind3] - Z[ind1]);
				}
		}
		if(!stopDual3D(Z,Zo,U1,U2,U3, (int)m,(int)n, (int)d, rho,absTol,relTol,display))
		{
		    memcpy(Zo, Z, m*n*d*sizeof(double));
			continue;
        }		
		if(stopPri3D(X1,X2,X3,Z,(int)m,(int)n,(int)d,absTol,relTol,display))
			break;
        memcpy(Zo, Z, m*n*d*sizeof(double));
	}
	#pragma omp parallel for shared(X1,X2,X3) private(i,j,k,ind1,ind2,ind3)
    for(i = 0; i<n; i++)
		for(j=0;j<m;j++)
			for(k=0;k<d;k++){
				ind1 = k*l + i*m + j;
				ind2 = k*l+ j*n + i;
				ind3 = l2*j+d*i + k;
				X1[ind1] +=  X2[ind2] + X3[ind3];
				X1[ind1] /=3.0;
			}

	free(T);
	free(X2);
	free(X3);
	free(Z);
	free(U1);
	free(U2);
	free(U3);
    free(Zo);
	return iter;
}
