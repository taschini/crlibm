#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "trigpi.h"

void sinpi_accurate(double *polyTblh, double *polyTblm, double *polyTbll,int index, int quadrant, int sgn, double y)
{
  double y2h, y2m, y2l;
  Mul123(&y2h,&y2m,&y2l,y,y,0.0);
  
  double sah, sam, sal, cah, cal, cam;
  double t1h, t1l, t1m, t2h, t2m, t2l, t3h, t3m, t3l, t4h, t4m, t4l, t5h, t5m, t5l, t6h, t6m, t6l, t7h, t7m, t7l, t8h, t8m, t8l, t9h, t9m, t9l, t10h, t10m, t10l, t11h, t11m, t11l, t12h, t12m, t12l, t13h, t13m, t13l, t14h, t14m, t14l, t15h, t15m, t15l, t16h, t16m, t16l, t17h, t17m, t17l, t18h, t18m, t18l, t19h, t19m, t19l, t20h, t20m, t20l;
  double temph, tempm, templ;
  double polysinh, polysinm, polysinl, polycosh, polycosm, polycosl;
  double cats_h, cats_m, cats_l, sacs_h, sacs_l, sacs_m;
  double resh, resm, resl;
  sah=sincosTable[index<<1].hi;
  sam=sincosTable[index<<1].mi;
  sal=sincosTable[index<<1].lo;
  cah=sincosTable[(index<<1)+1].hi;
  cam=sincosTable[(index<<1)+1].mi;
  cal=sincosTable[(index<<1)+1].lo;
  Mul12(&t1h,&t1m,c10h.d+y*(c11h.d+y*(c12h.d+y*c13h.d)),y);
  Add22(&t2h,&t2m,c9h.d,c9l.d,t1h,t1m);
  MulAdd212(&t3h,&t3m,c8h.d,c8l.d,y,t2h,t2m);
  MulAdd212(&t4h,&t4m,c7h.d,c7l.d,y,t3h,t3m);
  MulAdd212(&t5h,&t5m,c6h.d,c6l.d,y,t4h,t4m);
  MulAdd212(&t6h,&t6m,c5h.d,c5l.d,y,t5h,t5m);
  MulAdd212(&t7h,&t7m,c4h.d,c4l.d,y,t6h,t6m);
  MulAdd212(&t8h,&t8m,c3h.d,c3l.d,y,t7h,t7m);
  MulAdd212(&t9h,&t9m,c2h.d,c2l.d,y,t8h,t8m);
  MulAdd212(&t10h,&t10m,c1h.d,c1l.d,y,t9h,t9m);
  Mul123(&t11h,&t11m,&t11l,y,t10h,t10m);
  Add233Cond(&polycosh,&polycosm,&polycosl,c0h.d,c0l.d,t11h,t11m,t11l);



  Mul133(&t1h,&t1m,&t1l,s20h.d,y2h,y2m,y2l);
  Add133Cond(&t2h,&t2m,&t2l,s18h.d,t1h,t1m,t1l);
  Mul33(&t3h,&t3m,&t3l,t2h,t2m,t2l,y2h,y2m,y2l);
  Add133Cond(&t4h,&t4m,&t4l,s16h.d,t3h,t3m,t3l);
  Mul33(&t5h,&t5m,&t5l,y2h,y2m,y2l,t4h,t4m,t4l);
  Add133Cond(&t6h,&t6m,&t6l,s14h.d,t5h,t5m,t5l);
  Mul33(&t7h,&t7m,&t7l,t6h,t6m,t6l,y2h,y2m,y2l);
  Add133Cond(&t8h,&t8m,&t8l,s12h.d,t7h,t7m,t7l);
  Mul33(&t9h,&t9m,&t9l,y2h,y2m,y2l,t8h,t8m,t8l);
  Renormalize3(&temph,&tempm,&templ,t9h,t9m,t9l);
  Add133Cond(&t10h,&t10m,&t10l,s10h.d,temph,tempm,templ);
  Mul33(&t11h,&t11m,&t11l,y2h,y2m,y2l,t10h,t10m,t10l);
  Add233Cond(&t12h,&t12m,&t12l,s8h.d,s8l.d,t11h,t11m,t11l);
  Mul33(&t13h,&t13m,&t13l,y2h,y2m,y2l,t12h,t12m,t12l);
  Renormalize3(&temph,&tempm,&templ,t13h,t13m,t13l);
  Add233Cond(&t14h,&t14m,&t14l,s6h.d,s6l.d,temph,tempm,templ);
  Renormalize3(&temph,&tempm,&templ,t14h,t14m,t14l);
  Mul33(&t15h,&t15m,&t15l,y2h,y2m,y2l,temph,tempm,templ);
  Renormalize3(&temph,&tempm,&templ,t15h,t15m,t15l);
  Add233Cond(&t16h,&t16m,&t16l,s4h.d,s4l.d,temph,tempm,templ);
  Mul33(&t17h,&t17m,&t17l,t16h,t16m,t16l,y2h,y2m,y2l);
  Add33(&t18h,&t18m,&t18l,s2h.d,s2m.d,s2l.d,t17h,t17m,t17l);
  Mul33(&t19h,&t19m,&t19l,t18h,t18m,t18l,y2h,y2m,y2l);
  Add33(&t20h,&t20m,&t20l,s0h.d,s0m.d,s0l.d,t19h,t19m,0.0);
  Renormalize3(&temph,&tempm,&templ,t20h,t20m,t20l);
  Mul133(&polysinh,&polysinm,&polysinl,y,temph,tempm,0.0);
  if(quadrant==0)
  {
    Mul33(&sacs_h,&sacs_m,&sacs_l,sah,sam,sal,polycosh,polycosm,polycosl);
    Mul33(&cats_h,&cats_m,&cats_l,cah,cam,cal,polysinh,polysinm,polysinl);
    Add33(&resh,&resm,&resl,sacs_h,sacs_m,sacs_l,cats_h,cats_m,cats_l);
  }
  if(quadrant==1)
  {
    Mul33(&sacs_h,&sacs_m,&sacs_l,cah,cam,cal,polycosh,polycosm,polycosl);
    Mul33(&cats_h,&cats_m,&cats_l,sah,sam,sal,polysinh,polysinm, polysinl);
    Add33(&resh,&resm,&resl,sacs_h,sacs_m,sacs_l,cats_h,cats_m,cats_l);
  }
  if(sgn==-1)
  {
    resh=-resh;
    resl=-resl;
    resm=-resm;
  }
  Renormalize3(polyTblh,polyTblm,polyTbll,resh,resm,resl);
}

double sinpi_rn(double x){
  double tsl,tsh,tsxh,tsxl,ts2,csh,csl,csxh,csxl,y,y2,rem,kover256,sah,sal,cah,cal,res,resh,resm,resl,sacs_h,sacs_l,cats_h,cats_l;
  long int absxhi;
  int sgn,index,quadrant,infzero;
  db_number k,k2,x_split,xd;
  infzero=0;
  if(x<0) { x=-x; infzero=1; }
  x_split.d=x; xd.d=x;
  absxhi = x_split.i[HI] & 0x7fffffff;
  /* SPECIAL CASES: x=(Nan, Inf) sin(x)=Nan */
  if (absxhi>=0x7ff00000) {
    x_split.l=0xfff8000000000000LL;
    return x_split.d - x_split.d;
  }
  
  /* We multiply x by 2^(-1074) in two stpes (x*2^(-950)*2^(-124)); 
     the result of this operation is a denormalized number */
  k.d=x;
  k.i[HI]+=(-950)<<20;
  k.d=k.d*TWOPOWERM124;
  k2.i[HI]=k.i[HI];
  k2.i[LO]=k.i[LO];
  /* The sign of sin(Pi*x) depends of the sign and of the parity of x. */
  if(infzero==0)
  {
    if (k2.i[LO]%2==0) sgn=1; else sgn=-1;
  }
   else
  {
    if (k2.i[LO]%2==0) sgn=-1; else sgn=1;
  }
  /* If x<XMAX_RETURN_PIX_FOR_SIN we can return sin(Pi*x)=Pi*x */
  if (absxhi <XMAX_RETURN_PIX_FOR_SIN)
  {
    Mul133(&resh,&resm,&resl,x,pi_hi.d,pi_mi.d,pi_lo.d);    
    RoundToNearest3(&res,resh,resm,resl);
    if (infzero==1) return -res; else return res;
  }
  /* x is an integer if it is greater than this bound, so we return 0 */
  if(xd.i[HI]>=0x43400000)
  {
    if (infzero==1) return -0.0; else return 0.0;
  }
  /* We multiply by 2^1074 in two steps; k2 is the nearest integer of x */
  k2.d=k2.d*TWOPOWER124;
  k2.d=k2.d*TWOPOWER950;

  rem=x-k2.d;
  if (rem>0.5) rem=x-(k2.d+1.0);
  if(x-k2.d<0) sgn=-sgn;
  if (rem<0.0) rem=-rem;
  if (rem==0.5) 
  {
    if(sgn==1) return 1.0; else return -1.0;
  }

  /* Now, we have 0<=x<0.5 */
  xd.d=rem+1.0;
  index=(xd.i[HI]>>12)&0x3f;
  if (rem<0.25) 
  {
    quadrant=0; 
    kover256=index*ONEOVER256;
  } 
   else 
  { 
    quadrant=1; 
    index=63-index;
    kover256=0.5-index*ONEOVER256;
  }
  /* We have 0<=y<4/512 */
  y=rem-kover256;

  sah=sincosTable[index<<1].hi;
  sal=sincosTable[index<<1].mi;
  cah=sincosTable[(index<<1)+1].hi;
  cal=sincosTable[(index<<1)+1].mi;

  if(y!=0.0)
  {
    /* We begin the polynomial evaluation */
    y2=y*y;
    ts2=y2*((s2h.d + y2*(s4h.d + y2*(s6h.d+y2*(s8h.d+y2*(s10h.d+y2*(s12h.d+y2*(s14h.d+y2*s16h.d))))))));
    Add22(&tsxh,&tsxl,s0h.d,s0m.d,ts2,0.0);
    Mul122(&tsh,&tsl,y,tsxh,tsxl);
    Mul12(&csxh,&csxl,(c1h.d+y*(c2h.d+y*(c3h.d+y*(c4h.d+y*(c5h.d+y*(c6h.d+y*(c7h.d+y*c8h.d))))))),y);
    Add12(csh,csl,c0h.d,csxh);
    if(quadrant==0)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,sah,sal);
      Mul22(&cats_h,&cats_l,cah,cal,tsh,tsl);
      Add22Cond(&resh,&resl,cats_h,cats_l,sacs_h,sacs_l);
    }
    if(quadrant==1)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,cah,cal);
      Mul22(&cats_h,&cats_l,sah,sal,tsh,tsl);
      Add22Cond(&resh,&resl,sacs_h,sacs_l,cats_h,cats_l);
    }
  }
   else
  {
    if(quadrant==0)
    {
      resh=sah;
      resl=sal;
    }
    if(quadrant==1)
    {
      resh=cah;
      resl=cal;
    }
  }
  if(resh == (resh + (resl * RN_CST_SINCOS)))
  {
    if(sgn==-1) return -resh; else return resh;
  }
   else 
  {
    sinpi_accurate(&resh,&resm,&resl,index,quadrant,sgn,y);
    RoundToNearest3(&res,resh,resm,resl);
    return res;
  }
}

double sinpi_rd(double x){
  double tsl,tsh,tsxh,tsxl,ts2,csh,csl,csxh,csxl,y,y2,rem,kover256,sah,sal,cah,cal,res,resh,resm,resl,sacs_h,sacs_l,cats_h,cats_l;
  long int absxhi;
  int sgn,index,quadrant,infzero;
  db_number k,k2,x_split,xd;
  infzero=0;
  if(x<0) { x=-x; infzero=1; }
  x_split.d=x; xd.d=x;
  absxhi = x_split.i[HI] & 0x7fffffff;
  /* SPECIAL CASES: x=(Nan, Inf) sin(x)=Nan */
  if (absxhi>=0x7ff00000) {
    x_split.l=0xfff8000000000000LL;
    return x_split.d - x_split.d;
  }
  
  /* We multiply x by 2^(-1074) in two stpes (x*2^(-950)*2^(-124)); 
     the result of this operation is a denormalized number */
  k.d=x;
  k.i[HI]+=(-950)<<20;
  k.d=k.d*TWOPOWERM124;
  k2.i[HI]=k.i[HI];
  k2.i[LO]=k.i[LO];
  /* The sign of sin(Pi*x) depends of the sign and of the parity of x. */
  if(infzero==0)
  {
    if (k2.i[LO]%2==0) sgn=1; else sgn=-1;
  }
   else
  {
    if (k2.i[LO]%2==0) sgn=-1; else sgn=1;
  }
  /* If x<XMAX_RETURN_PIX_FOR_SIN we can return sin(Pi*x)=Pi*x */
  if (absxhi <XMAX_RETURN_PIX_FOR_SIN)
  {
    Mul133(&resh,&resm,&resl,x,pi_hi.d,pi_mi.d,pi_lo.d);    
    RoundDownwards3(&res,resh,resm,resl);
    if (infzero==1) return -res; else return res;
  }
  /* x is an integer if it is greater than this bound, so we return 0 */
  if(xd.i[HI]>=0x43400000)
  {
    if (infzero==1) return -0.0; else return 0.0;
  }
  /* We multiply by 2^1074 in two steps; k2 is the nearest integer of x */
  k2.d=k2.d*TWOPOWER124;
  k2.d=k2.d*TWOPOWER950;

  rem=x-k2.d;
  if (rem>0.5) rem=x-(k2.d+1.0);
  if(x-k2.d<0) sgn=-sgn;
  if (rem<0.0) rem=-rem;
  if (rem==0.5) 
  {
    if(sgn==1) return 1.0; else return -1.0;
  }

  /* Now, we have 0<=x<0.5 */
  xd.d=rem+1.0;
  index=(xd.i[HI]>>12)&0x3f;
  if (rem<0.25) 
  {
    quadrant=0; 
    kover256=index*ONEOVER256;
  } 
   else 
  { 
    quadrant=1; 
    index=63-index;
    kover256=0.5-index*ONEOVER256;
  }
  /* We have 0<=y<4/512 */
  y=rem-kover256;

  sah=sincosTable[index<<1].hi;
  sal=sincosTable[index<<1].mi;
  cah=sincosTable[(index<<1)+1].hi;
  cal=sincosTable[(index<<1)+1].mi;

  if(y!=0.0)
  {
    /* We begin the polynomial evaluation */
    y2=y*y;
    ts2=y2*((s2h.d + y2*(s4h.d + y2*(s6h.d+y2*(s8h.d+y2*(s10h.d+y2*(s12h.d+y2*(s14h.d+y2*s16h.d))))))));
    Add22(&tsxh,&tsxl,s0h.d,s0m.d,ts2,0.0);
    Mul122(&tsh,&tsl,y,tsxh,tsxl);
    Mul12(&csxh,&csxl,(c1h.d+y*(c2h.d+y*(c3h.d+y*(c4h.d+y*(c5h.d+y*(c6h.d+y*(c7h.d+y*c8h.d))))))),y);
    Add12(csh,csl,c0h.d,csxh);
    if(quadrant==0)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,sah,sal);
      Mul22(&cats_h,&cats_l,cah,cal,tsh,tsl);
      Add22Cond(&resh,&resl,cats_h,cats_l,sacs_h,sacs_l);
    }
    if(quadrant==1)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,cah,cal);
      Mul22(&cats_h,&cats_l,sah,sal,tsh,tsl);
      Add22Cond(&resh,&resl,sacs_h,sacs_l,cats_h,cats_l);
    }
  }
   else
  {
    if(quadrant==0)
    {
      resh=sah;
      resl=sal;
    }
    if(quadrant==1)
    {
      resh=cah;
      resl=cal;
    }
  }
  if(sgn==-1) { 
    resh=-resh;
    resl=-resl;
  }
  TEST_AND_RETURN_RD(resh,resl,RDROUNDCST);
  sinpi_accurate(&resh,&resm,&resl,index,quadrant,sgn,y);
  RoundDownwards3(&res,resh,resm,resl);
  return res;
}

double sinpi_ru(double x){
  double tsl,tsh,tsxh,tsxl,ts2,csh,csl,csxh,csxl,y,y2,rem,kover256,sah,sal,cah,cal,res,resh,resm,resl,sacs_h,sacs_l,cats_h,cats_l;
  long int absxhi;
  int sgn,index,quadrant,infzero;
  db_number k,k2,x_split,xd;
  infzero=0;
  if(x<0) { x=-x; infzero=1; }
  x_split.d=x; xd.d=x;
  absxhi = x_split.i[HI] & 0x7fffffff;
  /* SPECIAL CASES: x=(Nan, Inf) sin(x)=Nan */
  if (absxhi>=0x7ff00000) {
    x_split.l=0xfff8000000000000LL;
    return x_split.d - x_split.d;
  }
  
  /* We multiply x by 2^(-1074) in two stpes (x*2^(-950)*2^(-124)); 
     the result of this operation is a denormalized number */
  k.d=x;
  k.i[HI]+=(-950)<<20;
  k.d=k.d*TWOPOWERM124;
  k2.i[HI]=k.i[HI];
  k2.i[LO]=k.i[LO];
  /* The sign of sin(Pi*x) depends of the sign and of the parity of x. */
  if(infzero==0)
  {
    if (k2.i[LO]%2==0) sgn=1; else sgn=-1;
  }
   else
  {
    if (k2.i[LO]%2==0) sgn=-1; else sgn=1;
  }
  /* If x<XMAX_RETURN_PIX_FOR_SIN we can return sin(Pi*x)=Pi*x */
  if (absxhi <XMAX_RETURN_PIX_FOR_SIN)
  {
    Mul133(&resh,&resm,&resl,x,pi_hi.d,pi_mi.d,pi_lo.d);    
    RoundUpwards3(&res,resh,resm,resl);
    if (infzero==1) return -res; else return res;
  }
  /* x is an integer if it is greater than this bound, so we return 0 */
  if(xd.i[HI]>=0x43400000)
  {
    if (infzero==1) return -0.0; else return 0.0;
  }
  /* We multiply by 2^1074 in two steps; k2 is the nearest integer of x */
  k2.d=k2.d*TWOPOWER124;
  k2.d=k2.d*TWOPOWER950;

  rem=x-k2.d;
  if (rem>0.5) rem=x-(k2.d+1.0);
  if(x-k2.d<0) sgn=-sgn;
  if (rem<0.0) rem=-rem;
  if (rem==0.5) 
  {
    if(sgn==1) return 1.0; else return -1.0;
  }

  /* Now, we have 0<=x<0.5 */
  xd.d=rem+1.0;
  index=(xd.i[HI]>>12)&0x3f;
  if (rem<0.25) 
  {
    quadrant=0; 
    kover256=index*ONEOVER256;
  } 
   else 
  { 
    quadrant=1; 
    index=63-index;
    kover256=0.5-index*ONEOVER256;
  }
  /* We have 0<=y<4/512 */
  y=rem-kover256;

  sah=sincosTable[index<<1].hi;
  sal=sincosTable[index<<1].mi;
  cah=sincosTable[(index<<1)+1].hi;
  cal=sincosTable[(index<<1)+1].mi;

  if(y!=0.0)
  {
    /* We begin the polynomial evaluation */
    y2=y*y;
    ts2=y2*((s2h.d + y2*(s4h.d + y2*(s6h.d+y2*(s8h.d+y2*(s10h.d+y2*(s12h.d+y2*(s14h.d+y2*s16h.d))))))));
    Add22(&tsxh,&tsxl,s0h.d,s0m.d,ts2,0.0);
    Mul122(&tsh,&tsl,y,tsxh,tsxl);
    Mul12(&csxh,&csxl,(c1h.d+y*(c2h.d+y*(c3h.d+y*(c4h.d+y*(c5h.d+y*(c6h.d+y*(c7h.d+y*c8h.d))))))),y);
    Add12(csh,csl,c0h.d,csxh);
    if(quadrant==0)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,sah,sal);
      Mul22(&cats_h,&cats_l,cah,cal,tsh,tsl);
      Add22Cond(&resh,&resl,cats_h,cats_l,sacs_h,sacs_l);
    }
    if(quadrant==1)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,cah,cal);
      Mul22(&cats_h,&cats_l,sah,sal,tsh,tsl);
      Add22Cond(&resh,&resl,sacs_h,sacs_l,cats_h,cats_l);
    }
  }
   else
  {
    if(quadrant==0)
    {
      resh=sah;
      resl=sal;
    }
    if(quadrant==1)
    {
      resh=cah;
      resl=cal;
    }
  }
  if(sgn==-1) { 
    resh=-resh;
    resl=-resl;
  }
  TEST_AND_RETURN_RU(resh,resl,RDROUNDCST);
  sinpi_accurate(&resh,&resm,&resl,index,quadrant,sgn,y);
  RoundUpwards3(&res,resh,resm,resl);
  return res;
}

double sinpi_rz(double x){
  double tsl,tsh,tsxh,tsxl,ts2,csh,csl,csxh,csxl,y,y2,rem,kover256,sah,sal,cah,cal,res,resh,resm,resl,sacs_h,sacs_l,cats_h,cats_l;
  long int absxhi;
  int sgn,index,quadrant,infzero;
  db_number k,k2,x_split,xd;
  infzero=0;
  if(x<0) { x=-x; infzero=1; }
  x_split.d=x; xd.d=x;
  absxhi = x_split.i[HI] & 0x7fffffff;
  /* SPECIAL CASES: x=(Nan, Inf) sin(x)=Nan */
  if (absxhi>=0x7ff00000) {
    x_split.l=0xfff8000000000000LL;
    return x_split.d - x_split.d;
  }
  
  /* We multiply x by 2^(-1074) in two stpes (x*2^(-950)*2^(-124)); 
     the result of this operation is a denormalized number */
  k.d=x;
  k.i[HI]+=(-950)<<20;
  k.d=k.d*TWOPOWERM124;
  k2.i[HI]=k.i[HI];
  k2.i[LO]=k.i[LO];
  /* The sign of sin(Pi*x) depends of the sign and of the parity of x. */
  if(infzero==0)
  {
    if (k2.i[LO]%2==0) sgn=1; else sgn=-1;
  }
   else
  {
    if (k2.i[LO]%2==0) sgn=-1; else sgn=1;
  }
  /* If x<XMAX_RETURN_PIX_FOR_SIN we can return sin(Pi*x)=Pi*x */
  if (absxhi <XMAX_RETURN_PIX_FOR_SIN)
  {
    Mul133(&resh,&resm,&resl,x,pi_hi.d,pi_mi.d,pi_lo.d);    
    RoundTowardsZero3(&res,resh,resm,resl);
    if (infzero==1) return -res; else return res;
  }
  /* x is an integer if it is greater than this bound, so we return 0 */
  if(xd.i[HI]>=0x43400000)
  {
    if (infzero==1) return -0.0; else return 0.0;
  }
  /* We multiply by 2^1074 in two steps; k2 is the nearest integer of x */
  k2.d=k2.d*TWOPOWER124;
  k2.d=k2.d*TWOPOWER950;

  rem=x-k2.d;
  if (rem>0.5) rem=x-(k2.d+1.0);
  if(x-k2.d<0) sgn=-sgn;
  if (rem<0.0) rem=-rem;
  if (rem==0.5) 
  {
    if(sgn==1) return 1.0; else return -1.0;
  }

  /* Now, we have 0<=x<0.5 */
  xd.d=rem+1.0;
  index=(xd.i[HI]>>12)&0x3f;
  if (rem<0.25) 
  {
    quadrant=0; 
    kover256=index*ONEOVER256;
  } 
   else 
  { 
    quadrant=1; 
    index=63-index;
    kover256=0.5-index*ONEOVER256;
  }
  /* We have 0<=y<4/512 */
  y=rem-kover256;

  sah=sincosTable[index<<1].hi;
  sal=sincosTable[index<<1].mi;
  cah=sincosTable[(index<<1)+1].hi;
  cal=sincosTable[(index<<1)+1].mi;

  if(y!=0.0)
  {
    /* We begin the polynomial evaluation */
    y2=y*y;
    ts2=y2*((s2h.d + y2*(s4h.d + y2*(s6h.d+y2*(s8h.d+y2*(s10h.d+y2*(s12h.d+y2*(s14h.d+y2*s16h.d))))))));
    Add22(&tsxh,&tsxl,s0h.d,s0m.d,ts2,0.0);
    Mul122(&tsh,&tsl,y,tsxh,tsxl);
    Mul12(&csxh,&csxl,(c1h.d+y*(c2h.d+y*(c3h.d+y*(c4h.d+y*(c5h.d+y*(c6h.d+y*(c7h.d+y*c8h.d))))))),y);
    Add12(csh,csl,c0h.d,csxh);
    if(quadrant==0)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,sah,sal);
      Mul22(&cats_h,&cats_l,cah,cal,tsh,tsl);
      Add22Cond(&resh,&resl,cats_h,cats_l,sacs_h,sacs_l);
    }
    if(quadrant==1)
    {
      Mul22(&sacs_h,&sacs_l,csh,csl,cah,cal);
      Mul22(&cats_h,&cats_l,sah,sal,tsh,tsl);
      Add22Cond(&resh,&resl,sacs_h,sacs_l,cats_h,cats_l);
    }
  }
   else
  {
    if(quadrant==0)
    {
      resh=sah;
      resl=sal;
    }
    if(quadrant==1)
    {
      resh=cah;
      resl=cal;
    }
  }
  if(sgn==-1) { 
    resh=-resh;
    resl=-resl;
  }
  TEST_AND_RETURN_RZ(resh,resl,RDROUNDCST);
  sinpi_accurate(&resh,&resm,&resl,index,quadrant,sgn,y);
  RoundTowardsZero3(&res,resh,resm,resl);
  return res;
}


