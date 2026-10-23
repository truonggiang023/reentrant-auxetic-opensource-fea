/* Finite-strain multiplicative Hencky J2 material for CalculiX 2.20.
 * Stress = second Piola-Kirchhoff; tangent = dS/dGreen strain, tensor shear.
 * Native decks use nine descriptors and an immutable REENTRANT_MATERIAL_DATA
 * sidecar; full inline tables below are only for the C material-point API and
 * <=21-property coupons. Stock CCX 2.20 must never receive larger tables.
 * State: inverse plastic metric (six), alpha, omega, post-initiation u, D.
 * No mesh/topology mutation: residual stiffness retains failed elements.
 * Properties: E,nu,uf,residual,L,nh,nf,eta_cutoff, hardening(sigma,alpha),
 *             fracture locus(eta,epsilon_f), viscosity_time. Units: mm, N, MPa.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static double *external_tables=NULL;
static int external_nh=0,external_nf=0;
__attribute__((constructor)) static void identify(void){
 puts("REENTRANT_MATERIAL_ABI=3; multiplicative Hencky J2; residual elements");
 const char *path=getenv("REENTRANT_MATERIAL_DATA");
 if(!path)return; /* Inline <=21-property coupon tests need no sidecar. */
 FILE *f=fopen(path,"r");char magic[64];
 if(!f||fscanf(f,"%63s %d %d",magic,&external_nh,&external_nf)!=3||
    strcmp(magic,"REENTRANT_TABLES_V1")||external_nh<2||external_nf<2||
    external_nh+external_nf>10000){fprintf(stderr,"Invalid REENTRANT material sidecar\n");exit(201);}
 external_tables=malloc(2*(external_nh+external_nf)*sizeof(double));
 if(!external_tables){fprintf(stderr,"Cannot allocate REENTRANT tables\n");exit(201);}
 for(int i=0;i<2*(external_nh+external_nf);i++)
  if(fscanf(f,"%lf",external_tables+i)!=1||!isfinite(external_tables[i])){
   fprintf(stderr,"Invalid REENTRANT table value\n");exit(201);}
 fclose(f);
 for(int k=0;k<external_nh;k++){
  double sigma=external_tables[2*k],alpha=external_tables[2*k+1];
  if(sigma<=0||(k==0&&alpha!=0)||(k>0&&(alpha<=external_tables[2*k-1]||sigma<external_tables[2*k-2]))){
   fprintf(stderr,"Invalid hardening table in REENTRANT sidecar\n");exit(201);}
 }
 for(int k=0;k<external_nf;k++){
  const double*t=external_tables+2*external_nh;
  if(t[2*k+1]<=0||(k>0&&t[2*k]<=t[2*k-2])){
   fprintf(stderr,"Invalid fracture locus in REENTRANT sidecar\n");exit(201);}
 }
}
extern void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
static const int vi[6]={0,1,2,0,0,1}, vj[6]={0,1,2,1,2,2};
static void unpack(const double*v,double*a){
 memset(a,0,9*sizeof(double));
 for(int k=0;k<6;k++)a[3*vi[k]+vj[k]]=a[3*vj[k]+vi[k]]=v[k];
}
static void pack(const double*a,double*v){for(int k=0;k<6;k++)v[k]=a[3*vi[k]+vj[k]];}
static void mul(const double*a,const double*b,double*c){
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){c[3*i+j]=0;for(int k=0;k<3;k++)c[3*i+j]+=a[3*i+k]*b[3*k+j];}
}
static void spectral(const double*z,const double*w,double*a){
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){a[3*i+j]=0;for(int k=0;k<3;k++)a[3*i+j]+=z[i+3*k]*w[k]*z[j+3*k];}
}
static int eigen(const double*a,double*z,double*w){
 memcpy(z,a,9*sizeof(double));int n=3,ld=3,lw=32,info;double work[32];char v='V',u='U';
 dsyev_(&v,&u,&n,z,&ld,w,work,&lw,&info);return info;
}
static double interp(double x,const double*t,int n,int xcol){
 if(x<=t[xcol])return t[1-xcol];
 for(int k=1;k<n;k++)if(x<=t[2*k+xcol]){
 double x0=t[2*k-2+xcol],x1=t[2*k+xcol];
 return t[2*k-2+1-xcol]+(x-x0)/(x1-x0)*(t[2*k+1-xcol]-t[2*k-2+1-xcol]);}
 return t[2*(n-1)+1-xcol];
}
int reentrant_response_dt(const double*C,const double*old,const double*p,int elastic,double dt,double*S,double*state){
 double E=p[0],nu=p[1],G=E/(2*(1+nu)),K=E/(3*(1-2*nu));int nh=(int)p[5],nf=(int)p[6];
 const double*hard=p+8,*locus=hard+2*nh;
 double z[9],w[3],U[9],invU[9],cp[9],temp[9],B[9];
 if(eigen(C,z,w)||w[0]<=0)return 1;
 double lnJ=0;for(int i=0;i<3;i++){lnJ+=.5*log(w[i]);w[i]=sqrt(w[i]);}
 spectral(z,w,U);for(int i=0;i<3;i++)w[i]=1/w[i];spectral(z,w,invU);
 unpack(old,cp);if(cp[0]+cp[4]+cp[8]<1e-12){memset(cp,0,sizeof(cp));cp[0]=cp[4]=cp[8]=1;}
 mul(U,cp,temp);mul(temp,U,B);
 if(eigen(B,z,w)||w[0]<=0)return 2;
 double eps[3],dev[3],tau[3],q=0;
 for(int i=0;i<3;i++){eps[i]=.5*log(w[i]);dev[i]=eps[i]-lnJ/3; q+=dev[i]*dev[i];}
 q=2*G*sqrt(1.5*q);
 double J=exp(lnJ);
 double dp=0,y=J*interp(old[6],hard,nh,1);
 if(!elastic&&q>y){double lo=0,hi=q/(3*G);
 for(int k=0;k<60;k++){double m=.5*(lo+hi);if(q-3*G*m>J*interp(old[6]+m,hard,nh,1))lo=m;else hi=m;} dp=.5*(lo+hi);}
 double scale=q>1e-14?fmax(0,1-3*G*dp/q):1;
 double eta=q*scale>1e-12?K*lnJ/(q*scale):0;
 double domega=eta>p[7]?dp/interp(eta,locus,nf,0):0;
 double omega=old[7]+domega;
 double fraction=domega>0?fmin(1,fmax(0,(omega-1)/domega)):(old[7]>=1?1:0);
 double post=old[8]+p[4]*dp*fraction;
 double target=fmin(1-p[3],post/p[2]),viscosity=p[8+2*nh+2*nf];
 double D=elastic?old[9]:fmax(old[9],(viscosity*old[9]+dt*target)/(viscosity+dt));
 for(int i=0;i<3;i++){tau[i]=(1-D)*(K*lnJ+2*G*scale*dev[i]);w[i]=exp(2*(lnJ/3+scale*dev[i]));}
 double T[9],Smat[9],Bn[9],cpn[9];spectral(z,tau,T);
 mul(invU,T,temp);mul(temp,invU,Smat);pack(Smat,S);
 spectral(z,w,Bn);mul(invU,Bn,temp);mul(temp,invU,cpn);pack(cpn,state);
 state[6]=old[6]+dp;state[7]=omega;state[8]=post;state[9]=D;
 return 0;
}
int reentrant_response(const double*C,const double*old,const double*p,int elastic,double*S,double*state){
 return reentrant_response_dt(C,old,p,elastic,1.,S,state);
}
/* Native CalculiX interface, 32-bit integer ABI. Arrays are Fortran column major. */
void umat_(const char*name,const int*iel,const int*ip,const int*kode,const double*p,
 const double*e,const double*e0,const double*beta,const double*F0,const double*J0,
 const double*F1,const double*J1,const int*thermal,const double*temp,const double*dt,
 const double*time,const double*ttime,const int*icmd,const int*elastic,const int*mi,
 const int*ns,const double*s0,double*s1,double*stress,double*tangent,
 const int*ori,const double*gauss,const double*orab,double*newdt,const int*ipkon,int len){
 (void)name;(void)e0;(void)beta;(void)F0;(void)J0;(void)F1;(void)J1;(void)thermal;(void)temp;
 (void)icmd;(void)dt;(void)time;(void)ttime;(void)ori;(void)gauss;(void)orab;(void)ipkon;(void)len;
 int count=-(*kode)-100;
 if(p[0]<=0||p[1]<=-1||p[1]>=.5||p[2]<=0||p[3]<=0||p[3]>=1||p[4]<=0){
  fprintf(stderr,"Invalid REENTRANT physical parameters\n");exit(201);}

 if(*ns<10||count<9||count>21){fprintf(stderr,"REENTRANT requires 10 states and <=21 deck constants\n");exit(201);}
 double *expanded=NULL;
 if(p[5]<0){
  if(count!=9||!external_tables||-(int)p[5]!=external_nh||(int)p[6]!=external_nf){
   fprintf(stderr,"REENTRANT sidecar missing or inconsistent; set REENTRANT_MATERIAL_DATA\n");exit(201);}
  int n=9+2*(external_nh+external_nf);
  expanded=malloc(n*sizeof(double));
  if(!expanded)exit(201);
  memcpy(expanded,p,8*sizeof(double));expanded[5]=external_nh;
  memcpy(expanded+8,external_tables,2*(external_nh+external_nf)*sizeof(double));
  expanded[n-1]=p[8];p=expanded;
 }else if(count!=9+2*(int)p[5]+2*(int)p[6]){
  fprintf(stderr,"Invalid inline REENTRANT table counts\n");exit(201);
 }
 int off=(*iel-1)*(*ns)*mi[0]+(*ip-1)*(*ns);
 double C[9],old[10],state[10];unpack(e,C);for(int i=0;i<9;i++)C[i]*=2;C[0]++;C[4]++;C[8]++;
 memcpy(old,s0+off,sizeof(old));
 if(reentrant_response_dt(C,old,p,*elastic,*dt,stress,state)){*newdt=.25;free(expanded);return;}
 memcpy(s1+off,state,sizeof(state));
 /* Do not request cutbacks merely because a Newton trial has large damage.
  * Global line search must be allowed to correct trial states; xstateini
  * remains committed. Request pnewdt only for invalid constitutive kinematics. */
 /* Return a valid tangent even for stress-only requests. CCX 2.20 also
  * consults this array in forced-displacement predictors and hourglass paths. */
 double A[6][6];const double h=1e-7;
 for(int j=0;j<6;j++){
  double Cp[9],Cm[9],sp[6],sm[6],unused[10];memcpy(Cp,C,sizeof(C));memcpy(Cm,C,sizeof(C));
  Cp[3*vi[j]+vj[j]]+=2*h;Cm[3*vi[j]+vj[j]]-=2*h;
  if(j>=3){Cp[3*vj[j]+vi[j]]+=2*h;Cm[3*vj[j]+vi[j]]-=2*h;}
  if(reentrant_response_dt(Cp,old,p,*elastic,*dt,sp,unused)||reentrant_response_dt(Cm,old,p,*elastic,*dt,sm,unused)){*newdt=.25;free(expanded);return;}
  for(int i=0;i<6;i++)A[i][j]=(sp[i]-sm[i])/(2*h*(j>=3?2:1));
 }
 /* Native CCX interface accepts only symmetric material tangents. Use a
  * symmetric numerical algorithmic approximation; stress/state remain exact. */
 int k=0;for(int j=0;j<6;j++)for(int i=0;i<=j;i++)tangent[k++]=.5*(A[i][j]+A[j][i]);
 free(expanded);
}
