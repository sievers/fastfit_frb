float erfcx_jls(float x);
void get_derivs_1d_real_safe(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0);
void get_derivs_1d_real_neg(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0);
void get_derivs_1d_real_pos(float *y,int ny, float a,float b,float x0,float *f, float *dfda, float *dfdb, float *dfdx0);
void get_frb_derivs_compact(float *t_offset, int nt, float dt, float *freq, int nfreq, float a, float b, float t0, float dm, float amp, float alpha, float ref_freq,  float *f, float *dfda, float *dfdb, float *dfdt, float *dfdm, float *dfdamp, float *dfdind);

