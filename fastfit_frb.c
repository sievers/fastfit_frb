#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include <omp.h>
#include <assert.h>
#include <string.h>


#define DM0 4150
//#define SCATPOW -4.4
#define SCATPOW 4.0
//#define SCATPOW 0.0

//gcc-4.8 -O3 -fopenmp -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so
//gcc-5 -O3 -march=native -fPIC --shared fastfit_frb.c -o libfastfit_frb.so

float erfcx_jls(float x)
{

  
  static float pp[12]={1.24191724424408,
		       0.069834734623543,
		       -0.0358068266053169,
		       0.00756104756501004,
		       -0.000646903838206737,
		       -0.000188296978740529,
		       0.000100748126939714,
		       -2.31471757832635e-05,
		       1.33483371978706e-06,
		       1.29465312953991e-06,
		       -6.44459659333348e-07,
		       1.60424433876047e-07};
  static float pp2[9]={0.663675454804955,
		       -0.278469065276596,
		       0.0490341793971719,
		       -0.00759975102038169,
		       0.00106541074275422,
		       -0.000137503478564939,
		       1.65405454409028e-05,
		       -1.87148312096791e-06,
		       2.00570753007308e-07};
  
  
  if (x>1) {
    float xinv=1/x;
    //float xinv=1e-5*x; //in here for speed testing purposes. Doesn't seem to make much difference
    float xx=2*xinv-1;

    
    float tot=pp[0]+pp[1]*xx;
    float tm=1.0;
    float tc=xx;
    float tn;
    float xxx=2*xx;
    for (int j=2;j<12;j++) {
      //for (int j=2;j<7;j++) {
      tn=xxx*tc-tm;
      tm=tc;
      tc=tn;
      tot+=pp[j]*tn;
    }
    
    return tot/(1+2*x);
    //return tot;
  } else {
    if (x<0)
      return erfc(x)*exp(x*x);
    
    float xx=2*x-1;      
    float tot=pp2[0]+pp2[1]*xx;
    float tm=1.0;
    float tc=xx;
    float tn;
    for (int j=2;j<9;j++) {
      tn=2*xx*tc-tm;
      tm=tc;
      tc=tn;
      tot+=pp2[j]*tn;
    }
    return tot;
  }
  
  
  
  
}

/*--------------------------------------------------------------------------------*/
void get_derivs_1d_real_safe(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0)
{
#if 0
  for (int i=0;i<ny;i++) {
    f[i]=2.0*y[i];
    dfda[i]=3.0*y[i];
    dfdb[i]=4.3*y[i];
    dfdx0[i]=12.353*y[i];
  }
  return;
#endif

  float rt_b=sqrt(b);
  float b_inv=1/b;
  float rt_b_inv=1/rt_b;
  float rt_pi_inv=1/sqrt(M_PI);
  for (int i=0;i<ny;i++) {
    float u0=(x0-y[i]+a*b/2)*rt_b_inv;
    if (u0<0) {
      float to_exp=a*(x0-y[i]+0.25*a*b);
      float expval=exp(to_exp);
      float g=0.5*a*expval;
      float my_erfc=erfc(u0);
      float u0_exp=exp(-u0*u0);
      f[i]=g*my_erfc;
      
      float my_erfc_deriv=-u0_exp*rt_b*rt_pi_inv;
      //float g_deriv=g*(1+a*(x0-y[i]+0.5*a*b));
      float g_deriv=(0.5+0.5*a*(0.5*a*b+x0-y[i]))*expval;
      dfda[i]=g_deriv*my_erfc+g*my_erfc_deriv;   

      float g_deriv_b=0.25*a*a*g;
      float my_erfc_deriv_b=-(0.25*a*rt_b_inv-0.5*(x0-y[i])*rt_b_inv*b_inv)*u0_exp*2*rt_pi_inv;
      dfdb[i]=g_deriv_b*my_erfc+g*my_erfc_deriv_b;
      float g_deriv_x0=g*a;
      float my_erfc_deriv_x0=-u0_exp*rt_b_inv*2*rt_pi_inv;
      dfdx0[i]=g_deriv_x0*my_erfc+g*my_erfc_deriv_x0;
      
    } else {
      float to_exp=a*(x0-y[i]+a*b*0.25);
      float my_erfc_part=erfcx_jls(u0);
      float g_part=0.5*a;
      float exp_vec=exp(-u0*u0+to_exp);
      f[i]=my_erfc_part*g_part*exp_vec;
      
      float my_erfc_deriv_part=-rt_b*rt_pi_inv;
      //float g_deriv_part=(1+a*(x0-y[i]+0.5*a*b))*g_part;
      float g_deriv_part=0.5*(1+a*(a*b/2+x0-y[i]));
      dfda[i]=(g_deriv_part*my_erfc_part+g_part*my_erfc_deriv_part)*exp_vec;
      
      float g_deriv_b_part=0.25*a*a*g_part;
      float my_erfc_deriv_b_part=-(0.25*a*rt_b_inv-0.5*(x0-y[i])*b_inv*rt_b_inv)*2*rt_pi_inv;
      dfdb[i]=(g_deriv_b_part*my_erfc_part+g_part*my_erfc_deriv_b_part)*exp_vec;
      float g_deriv_x0_part=g_part*a;
      float my_erfc_deriv_x0_part=-2*rt_b_inv*rt_pi_inv;
      dfdx0[i]=(g_deriv_x0_part*my_erfc_part+g_part*my_erfc_deriv_x0_part)*exp_vec;
    }
  }
}

/*--------------------------------------------------------------------------------*/

void get_derivs_1d_real_neg(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0)
{
  float rt_b=sqrt(b);
  float b_inv=1/b;
  float rt_b_inv=1/rt_b;
  float rt_pi_inv=1/sqrt(M_PI);
  for (int i=0;i<ny;i++) {
    float u0=(x0-y[i]+a*b/2)*rt_b_inv;
    float to_exp=a*(x0-y[i]+a*b*0.25);
    float my_erfc_part=erfcx_jls(u0);
    float g_part=0.5*a;
    float exp_vec=exp(-u0*u0+to_exp);
    f[i]=my_erfc_part*g_part*exp_vec;

    float my_erfc_deriv_part=-rt_b*rt_pi_inv;
    float g_deriv_part=(1+a*(x0-y[i]+0.5*a*b))*g_part;
    dfda[i]=(g_deriv_part*my_erfc_part+g_part*my_erfc_deriv_part)*exp_vec;
    float g_deriv_b_part=0.25*a*a*g_part;
    float my_erfc_deriv_b_part=-(0.25*a*rt_b_inv-0.5*(x0-y[i])*b_inv*rt_b_inv)*2*rt_pi_inv;
    dfdb[i]=(g_deriv_b_part*my_erfc_part+g_part*my_erfc_deriv_b_part)*exp_vec;
    float g_deriv_x0_part=g_part*a;
    float my_erfc_deriv_x0_part=-2*rt_b_inv*rt_pi_inv;
    dfdx0[i]=(g_deriv_x0_part*my_erfc_part+g_part*my_erfc_deriv_x0_part)*exp_vec;
  }
}
/*--------------------------------------------------------------------------------*/

void get_derivs_1d_real_pos(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0)
{
  float rt_b=sqrt(b);
  float b_inv=1/b;
  float rt_b_inv=1/rt_b;
  float rt_pi_inv=1/sqrt(M_PI);
  for (int i=0;i<ny;i++) {
    float u0=(x0-y[i]+a*b/2)*rt_b_inv;
    //if (u0<0)
    //continue;
    float to_exp=a*(x0-y[i]+0.25*a*b);
    float g=0.5*a*exp(to_exp);
    float my_erfc=erfc(u0);
    float u0_exp=exp(-u0*u0);
    float my_erfc_deriv=-u0_exp*rt_b*rt_pi_inv;
    float g_deriv=g*(1+a*(x0-y[i]+0.5*a*b));
    f[i]=g*my_erfc;
    dfda[i]=g_deriv*my_erfc+g*my_erfc_deriv;   
    float g_deriv_b=0.25*a*a*g;
    float my_erfc_deriv_b=-(0.25*a*rt_b_inv-0.5*(x0-y[i])*rt_b_inv*b_inv)*u0_exp*2*rt_pi_inv;
    dfdb[i]=g_deriv_b*my_erfc+g*my_erfc_deriv_b;
    float g_deriv_x0=g*a;
    float my_erfc_deriv_x0=-u0_exp*rt_b_inv*2*rt_pi_inv;
    dfdx0[i]=g_deriv_x0*my_erfc+g*my_erfc_deriv_x0;
  }
}
/*--------------------------------------------------------------------------------*/
float *fvector(int n)
{
  float *v=(float *)malloc(sizeof(float)*n);
  if (v) {
    memset(v,0,sizeof(float)*n);
    return v;
  } else {
    printf("Malloc failure in fvector when requesting %d elements\n",n);
    return NULL;
  }

}

/*--------------------------------------------------------------------------------*/

void get_frb_derivs_compact(float *t_offset, int nt, float dt, float *freq, int nfreq, float a, float b, float t0, float dm, float amp, float alpha, float ref_freq,  float *f, float *dfda, float *dfdb, float *dfdt, float *dfdm, float *dfdamp, float *dfdind)
// scattering kernel is exp(-a*t), so usual scattering time is 1/a
// gaussian profile is exp(-b*(t-t0)^2), so usual gaussian sigma is sqrt(2*b)
{
  
  //printf("dt is %12.5f\n",dt);
  float *t_tmp=fvector(nt);
  float *dfda_tmp=fvector(nt);
  float *dfdb_tmp=fvector(nt);
  float *dfdx_tmp=fvector(nt);
  float *f_tmp=fvector(nt);
  float *lags=fvector(nfreq);
  double ref_lag=DM0*(1/(double)ref_freq/(double)ref_freq);
  for (int chan=0;chan<nfreq;chan++)
    lags[chan]=DM0/((double)freq[chan]*(double)freq[chan])-ref_lag;
  
  for (int chan=0;chan<nfreq;chan++) {
    float t_arr=t0+dm*lags[chan];
    
    for (int i=0;i<nt;i++) {
      // t_tmp[i]=(t0+lags[chan]-t_offset[chan])+i*dt;
      t_tmp[i]=t_offset[chan]+i*dt;
    }
    
    float tmp=freq[chan]/ref_freq;
    float logfac=log(tmp); //need this for spectral index fitting
    float scatfac=pow(tmp,SCATPOW);
    float a_eff=a*scatfac;  //this is the actual scattering at this frequency
    //get the derivatives for the channel as they depend on time, scattering, and intrinsic FWHM
    get_derivs_1d_real_safe(t_tmp,nt,a_eff,b,t0+dm*lags[chan],f_tmp,dfda_tmp,dfdb_tmp,dfdx_tmp);
    float nu_fac=pow(freq[chan]/ref_freq,alpha);
    float amp_fac=nu_fac*amp;

    for (int i=0;i<nt;i++) {
      int ii=chan*nt+i;
      
      f[ii]=f_tmp[i]*amp_fac;
      dfda[ii]=dfda_tmp[i]*amp_fac*scatfac;
      dfdb[ii]=dfdb_tmp[i]*amp_fac;
      dfdt[ii]=dfdx_tmp[i]*amp_fac;
      dfdm[ii]=dfdx_tmp[i]*amp_fac*lags[chan];
      dfdamp[ii]=f_tmp[i]*nu_fac;
      dfdind[ii]=f_tmp[i]*amp*logfac;
      
    }
    //void get_derivs_1d_real_safe(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0)            
  }
  

  free(lags);
  free(t_tmp);
  free(dfda_tmp);
  free(dfdb_tmp);
  free(dfdx_tmp);
  free(f_tmp);
}

/*================================================================================*/
#if 0
int main(int argc, char *argv[])
{
  int npt=10000000;
  float dt=2e-6;
  float *y=fvector(npt);
  for (int i=0;i<npt;i++)
    y[i]=dt*i;
  float t0=5;
  float a=0.5;
  float b=0.25;
  float *f=fvector(npt);
  float *dfda=fvector(npt);
  float *dfdb=fvector(npt);
  float *dfdx0=fvector(npt);
  //double t1=omp_get_wtime();
  get_derivs_1d_real_pos(y,npt,a,b,t0,f,dfda,dfdb,dfdx0);
  //double t2=omp_get_wtime();
  //printf("Took %12.5f seconds to run, time per point was %14.4e\n",t2-t1,(t2-t1)/npt);
  if (npt<1001) {
    FILE *outfile=fopen("fastfit_dump.txt","w");
    for (int i=0;i<npt;i++)
      fprintf(outfile,"%12.5f %14.7e %14.7e %14.7e %14.7e\n",y[i],f[i],dfda[i],dfdb[i],dfdx0[i]);
    fclose(outfile);
  }
  //t1=omp_get_wtime();
  get_derivs_1d_real_safe(y,npt,a,b,t0,f,dfda,dfdb,dfdx0);
  //t2=omp_get_wtime();
  //printf("Took %12.5f seconds to run, time per point was %14.4e\n",t2-t1,(t2-t1)/npt);
  if (npt<1001) {
    FILE *outfile=fopen("fastfit_dump_neg.txt","w");
    for (int i=0;i<npt;i++)
      fprintf(outfile,"%12.5f %14.7e %14.7e %14.7e %14.7e\n",y[i],f[i],dfda[i],dfdb[i],dfdx0[i]);
    fclose(outfile);
  }
  

}

#endif
