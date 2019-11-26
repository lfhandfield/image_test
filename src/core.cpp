/*
 * core.cpp
 *
 * Copyright (C) 2013 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "core.h"

namespace LFHPrimitive{


LFH_GOLD double lngamma(double xx){ double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j; y=x=xx; tmp=x+5.5f;tmp-=(x+0.5)*log(tmp);ser= 1.000000000190015;
    for(j=0;j<=5;j++)ser += cof[j]/++y;
    return(-tmp+log(2.5066282746310005*ser/x));
}
LFH_GOLD double lngammaexp(double xx){ double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j; y=x=exp(xx); tmp=x+5.5f;tmp-=(x+0.5)*log(tmp);ser= 1.000000000190015;
    for(j=0;j<=5;j++)ser += cof[j]/++y;
    return(-tmp+log(2.5066282746310005*ser/x));
}
LFH_GOLD void gamma(Magscaled<double> &fout, double xx){double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j; y=x=exp(xx); tmp=x+5.5f;tmp-=(x+0.5)*log(tmp);ser= 1.000000000190015;
    for(j=0;j<=5;j++)ser += cof[j]/++y;
    fout.value = 2.5066282746310005*ser/x;
    fout.exponent = -tmp;
}
LFH_GOLD double lngammadiff(double a, double b){
    double tmp,y,ser,ser2;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=a;
    tmp = -log_x_p1(b/a) + (a+0.5) * log_x_p1(b/ (a+5.5f)) + b * (log(a+b+5.5f) -1.0);
    ser2= 1.000000000190015;
    ser= 0;
    for(j=0;j<=5;j++) {
        ser2 += cof[j]/++y;
        ser -= cof[j] /(y*(y+b));
    }
    return(tmp+log_x_p1(ser * b/ser2)); // (sum C_i / (i+x)) / (sum C_i / (i+x+b))
}


LFH_GOLD double lngammadiffexp(double a, double b){
    double tmp,y,x,ser,ser2;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=exp(a);
    x = exp(b);
    //tmp = -log_x_p1(x/y) + (y+0.5) * log_x_p1(x/ (y+5.5f)) + x * (log(y+x+5.5f) -1.0);
    tmp = -log_e_x_p1(b-a) + (y+0.5) * log_e_x_p1(b - log(y+5.5f)) + x * (log(y+5.5f+x) -1.0);
    ser2= 1.000000000190015;
    ser= 0;

    for(j=0;j<=5;j++) {
        ser2 += cof[j]/++y;
        ser += cof[j] /(y*(y+x));
    }
    //printf("%e %e %e\n", ser,ser2,ser/ser2);

    return(tmp+log_1_me_x(b+ log(ser/ser2))); // (sum C_i / (i+x)) / (sum C_i / (i+x+b))
}
LFH_GOLD double lngammasum(double xx, double pivot){double x,y,tmp,ser;double x2,y2,tmp2,ser2;int j;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    y=x=pivot - xx;
    tmp=x+5.5f;
    tmp-=(x+0.5)*log(tmp);
    y2=x2=pivot + xx;
    tmp2=x2+5.5f;
    tmp2-=(x2+0.5)*log(tmp2);
    ser= 1.000000000190015; ser2= 1.000000000190015;
    for(j=0;j<=5;j++){
        ser += cof[j]/++y;
        ser2 += cof[j]/++y2;
    }
    return(-tmp+log(2.5066282746310005*ser/x)-tmp2+log(2.5066282746310005*ser2/x2));
}

double lngammachoose(double x, double y){  // log(gamma(y+1)) - log(gamma(x+1)) - log(gamma(y-x+1))
    y -= x;
    if (x < y) {return lngammadiff(x +1.0,y) - lngamma(y+1.0);}
    else {return lngammadiff(y+1.0,x) - lngamma(x+1.0);}
}


LFH_GOLD double lnbeta(const double a, const double b) {
return (b > a) ? lngamma(a) - lngammadiff(b,a):lngamma(b) - lngammadiff(a,b);}
LFH_GOLD double lnbetaexp(const double a, const double b) {
return (b > a) ? lngammaexp(a) - lngammadiffexp(b,a):lngammaexp(b) - lngammadiffexp(a,b);}
inline double gamma(double xx){return exp(lngamma(xx));}
LFH_GOLD double logdensity_chisquarre(double x, double free){ double fout = log(0.25) - x + (free - 2.0) *log(x * 0.5); fout *= 0.5; fout -= lngamma(free * 0.5); return fout;}
LFH_GOLD double Prob_NBdistrib(uint32_t x, double r, double p){return exp(lngamma(r + x) - lngamma((double)(x+1)) - lngamma(r)) * pow(p, (double)x) *pow(1.0 - p, (double)r);}




double log_exp_m1(double xx){ // log(e^x - 1), where x>0
    if (xx < 0.002) return (xx <= 0.0) ? nan("") : log(xx * (1.0 + 0.5*xx*(1.0 + xx/3.0*(1.0 + 0.25*xx*(1.0 + 0.2*xx)))));
    else return (xx < 2.0) ? log(exp(xx) - 1.0) : xx + log(1.0 - exp(-xx));
}

double d_log_exp_m1_dx(double xx){ // 1 / (1 - e^-x) , where x>0
    if (xx < 0.002f) return (xx <= 0.0) ? nan("") : 1.0d / (xx * (1.0 - 0.5*xx*(1.0 - xx/3.0*(1.0 - 0.25*xx*(1.0 - 0.2*xx)))));
    else return 1.0 / (1.0 - exp(-xx));
}

double d_tanh(double xx){
    double tmp;
    if (fabs(xx) < 8.0){
        tmp = tanh(xx);
        return 1.0 - tmp * tmp;
    }
    tmp = exp(-2.0 * fabs(xx));
    return 8.0 * tmp * (0.5 - tmp);
}

double d2_tanh(double xx){
    double tmp;
    if (fabs(xx) < 8.0){
        tmp = tanh(xx);
        return -2.0 * tmp * (1.0 - tmp * tmp);
    }
    tmp = exp(-2.0 * fabs(xx));
    return tanh(xx) * -16.0 * tmp * (0.5 - tmp);
}



/*
double lnBesselScaled_Kn(double n, double x){ // Modified bessel function of second , multiplied by  * |x|^{n}
  n = abs(n);
  x = abs(x);

  if(n == 0) {
    return log(gsl_sf_bessel_K0_scaled(x));
  }
  else if(n == 1) {
    return log(gsl_sf_bessel_K1_scaled(x));
  }
  else if(x <= 5.0) {
    return bessel_Kn_scaled_small_x(n, x, result);
  }
  else if(GSL_ROOT3_DBL_EPSILON * x > 0.25 * (n*n + 1)) {
    return gsl_sf_bessel_Knu_scaled_asympx_e((double)n, x, result);
  }
  else if(GSL_MIN(0.29/(n*n), 0.5/(n*n + x*x)) < GSL_ROOT3_DBL_EPSILON) {
    return gsl_sf_bessel_Knu_scaled_asymp_unif_e((double)n, x, result);
  }
  else {

    double two_over_x = 2.0/x;
    gsl_sf_result r_b_jm1;
    gsl_sf_result r_b_j;
    int stat_0 = gsl_sf_bessel_K0_scaled_e(x, &r_b_jm1);
    int stat_1 = gsl_sf_bessel_K1_scaled_e(x, &r_b_j);
    double b_jm1 = r_b_jm1.val;
    double b_j   = r_b_j.val;
    double b_jp1;
    int j;

    for(j=1; j<n; j++) {
      b_jp1 = b_jm1 + j * two_over_x * b_j;
      b_jm1 = b_j;
      b_j   = b_jp1;
    }

    result->val  = b_j;
    result->err  = n * (fabs(b_j) * (fabs(r_b_jm1.err/r_b_jm1.val) + fabs(r_b_j.err/r_b_j.val)));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
}*/

double log_x_p1(const double x){
    if (fabs(x) < exp(-8.0)) return x*(1.0 - x *( 0.5 - x * ((1.0/3) - x * (0.25 - 0.2 *x))));
    else return log(1.0 + x);
}
double log_x_1plus_overx(const double x){
    if (fabs(x) < exp(-8.0)) return 1.0 - x *( 0.5 - x * ((1.0/3) - x * (0.25 - 0.2 *x)));
    else return log(1.0 + x) / x;
}
double log_e_x_1plus_overe_x(const double x){
    double tmp;
    if (x < -8.0) {
        tmp = exp(x);
        return 1.0 - tmp *( 0.5 - tmp * ((1.0/3) - tmp * (0.25 - 0.2 *tmp)));
    }else{
        tmp = exp(-x);
        return (log(1.0 + tmp) + x) * tmp;
    }
}

double logistic(const double x){
    if (x > 18.75f) return 1.0 - exp(-x) + exp(-2.0*x);
    else if (x < -18.375f) return exp(x) - exp(2.0*x);
    else return 1.0 / (1.0 + exp(-x));
}
double d_logistic_dx(double x){
    x = fabs(x);
    if (x > 18.75f) return exp(-x) - 2.0*exp(-2.0*x) + 3.0*exp(-3.0*x);
    x = exp(-x);
    return x / ((1.0 + x) *(1.0 + x));
}

double powerllogistic(const double x, const double y){ // log((1.0 + e^(-x))) ^ -y  = -y * log( 1 + exp(x))
    if (x < -36.75f) return y * x;
    double t = exp(-x);
    if  (x > 7.25f){
		const double c2 =  1.0/3.0;
		const double c4 =  1.0/5.0;
		const double c5 =  1.0/6.0;
		const double c6 =  1.0/7.0;
		const double c8 =  1.0/9.0;
		const double c9 =  1.0/10.0;
		return y * t*(-1.0 + t * (0.5 - t*(c2 - t*(0.25 - t*(c4 - t*(c5 - t*(c6 - t*(0.125 - t*(c8 - t*c9)))))))));
    }else return log(t + 1.0) * -y;
}
double powerllogit(const double x, const double y){// computes -log(x^(1/-y) - 1.0) = -log( e^(log(x)/-y) - 1.0)
    double lx = log(x) / y;  log(x) > 36.75 * y;
    if (lx > 36.75f) return -lx;
    else if (lx > 9.1188196555451620800313608440928e-4) return -log( exp(lx) -1.0);
    const double c2 =  1.0/3.0;
    const double c4 =  1.0/5.0;
    const double c5 =  1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c8 =  1.0/9.0;
    const double c9 =  1.0/10.0;
    return -log(lx*(1.0 + lx * (0.5 + lx*(c2 + lx*(0.25 + lx*(c4 + lx*(c5 + lx*(c6 + lx*(0.125 + lx*(c8 + lx*c9))))))))));
}

double log_e_x_p1(const double x){
    if (x > 36.75f) return x;
    double t = exp(x);
    if (x < -18.375f) return t;
    else if  (x < -7.25f){
		const double c2 =  1.0/3.0;
		const double c4 =  1.0/5.0;
		const double c5 =  1.0/6.0;
		const double c6 =  1.0/7.0;
		const double c8 =  1.0/9.0;
		const double c9 =  1.0/10.0;
		return t*(1.0 - t * (0.5 - t*(c2 - t*(0.25 - t*(c4 - t*(c5 - t*(c6 - t*(0.125 - t*(c8 - t*c9)))))))));
    }else return log(t + 1.0);
}
double log_e_x_m1(const double x){
    if (x > 36.75f) return x;
    else if (x > 9.1188196555451620800313608440928e-4) return log( exp(x) -1.0);
    const double c2 =  1.0/3.0;
    const double c4 =  1.0/5.0;
    const double c5 =  1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c8 =  1.0/9.0;
    const double c9 =  1.0/10.0;
    return log(x*(1.0 + x * (0.5 + x*(c2 + x*(0.25 + x*(c4 + x*(c5 + x*(c6 + x*(0.125 + x*(c8 + x*c9))))))))));
}
double log_1_me_x(double x){ // log( 1 - exp(x))
    if (x < -18.50f){
        x = exp(x);
        return x * (-1.0 + 0.5 * x);
    }else if (x < -9.1188196555451620800313608440928e-4) return log(1.0 - exp(x));
    const double c2 =  1.0/3.0;
    const double c4 =  1.0/5.0;
    const double c5 =  1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c8 =  1.0/9.0;
    const double c9 =  1.0/10.0;
    return log(-x*(1.0 + x * (0.5 + x*(c2 + x*(0.25 + x*(c4 + x*(c5 + x*(c6 + x*(0.125 + x*(c8 + x*c9))))))))));
    //return log(-x) + 0.5 * x * (1.0 + (1.0/12.0) * x * (1.0 - x * x * (1.0 / 120.0) * (1.0 - x * x * (1 / 63.0))));
}
double log_1_me_x_v2(double x){ // log( 1 - exp(x))
    if (x < -18.50f){
        x = exp(x);
        return x * (1.0 - 0.5 * x);
    }else if (x < -9.1188196555451620800313608440928e-4) return log(1.0 - exp(x));
    return log(-x) + 0.5 * x * (1.0 + (1.0/12.0) * x * (1.0 - x * x * (1.0 / 120.0) * (1.0 - x * x * (1 / 63.0))));
}
double log_1_me_nx(const double x){
    if (x > 36.75f) return exp(-x);
    else if (x > 9.1188196555451620800313608440928e-4) return log(1.0 - exp(-x));
    const double c2 =  1.0/3.0;
    const double c4 =  1.0/5.0;
    const double c5 =  1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c8 =  1.0/9.0;
    const double c9 =  1.0/10.0;
    return log(x*(1.0 - x * (0.5 - x*(c2 - x*(0.25 - x*(c4 - x*(c5 - x*(c6 - x*(0.125 - x*(c8 - x*c9))))))))));
}
double d_log_1_me_nx_dx(const double x){
    if (x > 36.75f) return -exp(-x);
    else if (x > 9.1188196555451620800313608440928e-4) return -exp(-x) / (1.0 - exp(-x));
    const double c2 =  1.0/12.0;
    const double c3 =  1.0/60.0;
    const double c4 =  1.0/42.0;
    const double c5 =  1.0/40.0;
    double sx = x*x;
    return (1.0 / x) - 0.5 - c2 * x * (1.0 + c3 * sx * (1.0 + c4* sx * (1.0 + c5* sx)));
}

double e_x_1minus(const double x){
    if (fabs(x) > 0.03125f) return exp(x) - 1.0f;
    return x*(1.0 + x * (0.5 + x*(TEMPLATE_CONSTANTS<3>::ONE_OVER_X + x*(0.25 + x*(TEMPLATE_CONSTANTS<5>::ONE_OVER_X + x*(TEMPLATE_CONSTANTS<6>::ONE_OVER_X + x*(TEMPLATE_CONSTANTS<7>::ONE_OVER_X + x*(0.125 + x*(TEMPLATE_CONSTANTS<9>::ONE_OVER_X + x*TEMPLATE_CONSTANTS<10>::ONE_OVER_X)))))))));
}
double e_x_1minus_overx(const double x){
    if (fabs(x) > 0.03125f) return (exp(x) - 1.0) /x;
    return 1.0 + x * (0.5 + x*(TEMPLATE_CONSTANTS<3>::ONE_OVER_X + x*(0.25 + x*(TEMPLATE_CONSTANTS<5>::ONE_OVER_X + x*(TEMPLATE_CONSTANTS<6>::ONE_OVER_X + x*(TEMPLATE_CONSTANTS<7>::ONE_OVER_X + x*(0.125 + x*(TEMPLATE_CONSTANTS<9>::ONE_OVER_X + x*TEMPLATE_CONSTANTS<10>::ONE_OVER_X))))))));
}

double one_over_e_x_1minus(const double x){ // (domain x > 0), 1 / (exp(x) - 1)
    double tmp;
    if (x > 18.375f){
        tmp = exp(-x);
        return tmp * (1.0+tmp) ;
    }else if (fabs(x) > 0.03125f) return 1.0 / (exp(x) - 1.0);
    tmp =  x*(1.0 + x * (0.5 + x*(TEMPLATE_CONSTANTS<3>::ONE_OVER_X + x*(0.25 + x*(TEMPLATE_CONSTANTS<5>::ONE_OVER_X + x*(TEMPLATE_CONSTANTS<6>::ONE_OVER_X + x*(TEMPLATE_CONSTANTS<7>::ONE_OVER_X + x*(0.125 + x*(TEMPLATE_CONSTANTS<9>::ONE_OVER_X + x*TEMPLATE_CONSTANTS<10>::ONE_OVER_X)))))))));
    return 1.0 / tmp;
}

double logit(const double x){return log(x) - log(1.0 - x);}
double logit_signif(const double x, const double e){return-(log_e_x_m1(-log(x)-e));}
double logit(const Magscaled<double> x){return-(log_e_x_m1(-log(x.value)-x.exponent));}

double log_average(const double log_x, const double log_y){return ((log_x > log_y) ? log_x + log_e_x_p1(log_y - log_x) : log_y + log_e_x_p1(log_x - log_y)) - M_LN2;}

double logtransformed_sum(const double x,const double y){ // log(exp(x) + exp(y))
    if (x > y + 0.03125) return x + log(1 - exp(y-x));
    else if (y > x + 0.03125) return y + log(1 - exp(x-y));
    else {
        double diff = (x-y); diff = diff * diff * 0.125;
        return 0.5 * (x + y) + log(2) + diff * (1.0 - diff * (3.0 - 5.625 * diff));
    }
}

double sampleGaussian(){return(sqrt(-2 * log((1 + rand()) / ((float)RAND_MAX))) * cos(M_PI *2 * rand() / ((float)RAND_MAX)));}
double sampleGaussian(double mu, double sigma){return sampleGaussian() * sigma + mu;}
LFH_GOLD double randExponential(){return log(1.0 + RAND_MAX) - log(0.5 + rand());}
LFH_GOLD double randUniform(){return (0.5 + rand()) / (1.0 + RAND_MAX);}

static double erfc8_sum(double x){
  static double P[] = {2.97886562639399288862, 7.409740605964741794425, 6.1602098531096305440906, 5.019049726784267463450058, 1.275366644729965952479585264, 0.5641895835477550741253201704};
  static double Q[] = {3.3690752069827527677, 9.608965327192787870698, 17.08144074746600431571095, 12.0489519278551290360340491, 9.396034016235054150430579648, 2.260528520767326969591866945};
  double num=0.0, den=0.0;
  int i;
  num = P[5];
  den = x + Q[5];
  for (i=4; i>=0; --i) {den = x*den + Q[i];num = x*num + P[i];
  }
return num/den;}

double cumNormal(double x){return 0.5 * erfc(-sqrt(0.5) * x);}
double lncumNormal(double x){
    double x2;
    if (x >= 6.0) { x2 = x*x;
            x2 = exp(-0.5 * (x*x + log(M_PI * 2.0))) * (-1.0 + (1.0 - (3.0 - (15.0 - (105.0 - (945.0 - (10395.0 - 135135.0/ x2)/ x2)/ x2)/ x2) / x2) / x2) / x2) / x;
            return x2 * (1.0 - x2 * (0.5 - x2 * ((1.0/3.0) - x2*0.25)));
    }else if (x >= -32.0) return log(0.5 * erfc(-sqrt(0.5) * x));
    x2 = x*x;
    return -0.5 * (x2 + log(2.0*M_PI)) - log(-x - (1.0 - (2.0 - (10.0 - (74.0 - (706.0 - (8162.0 - (110410.0 - 1708394.0 / x2) / x2) / x2) / x2)/ x2) / x2) / x2 ) / x);
}
void cumNormal(Magscaled<double> &fout, double x){
    double x2;
    if (x >= 6.0) { x2 = x*x;
        x2 = exp(-0.5 * (x*x + log(M_PI * 2.0))) * (-1.0 + (1.0 - (3.0 - (15.0 - (105.0 - (945.0 - (10395.0 - 135135.0/ x2)/ x2)/ x2)/ x2) / x2) / x2) / x2) / x;
        fout.exponent = x2 * (1.0 - x2 * (0.5 - x2 * ((1.0/3.0) - x2*0.25)));
        fout.value = 1.0;
    }else if (x >= -32.0){
        fout.exponent = 0.0f;
        fout.value = 0.5 * erfc(-sqrt(0.5) * x);
    }else{
        x2 = x*x;
        fout.value = 1.0 / (-x - (1.0 - (2.0 - (10.0 - (74.0 - (706.0 - (8162.0 - (110410.0 - 1708394.0 / x2) / x2) / x2) / x2)/ x2) / x2) / x2 ) / x);
        fout.exponent = -0.5 * (x2 + log(2.0*M_PI));
    }
}
double logitcumNormal(double x){
    double x2;
    if (fabs(x) > 32.0){
        double x2 = x*x;
        double y = fabs(x);
        return 0.5 * (x * y + log(2.0*M_PI)) - log(y) - log(1.0 - (1.0 - (2.0 - (10.0 - (74.0 - (706.0 - (8162.0 - (110410.0 - 1708394.0 / x2) / x2) / x2) / x2)/ x2) / x2) / x2 ) / x2);
    }else{
        double x2 = erfc(-sqrt(0.5) * x);
        return log(0.5 * x2) - log(1.0 - 0.5 * x2);
    }
}



//= (-0.5 * x2 + log(...))
//=  - log(erfc(sqrt(0.5)*x))
//log(erfc(-sqrt(0.5) * x)) - log(2.0);

// 0.5 / erfc( sqrt(0.5) * x) \sim e^0.5*x^2 ( sqrt(0.125 * pi) * x + sqrt(0.125 * pi) / x - sqrt(pi) / x^3
// log(0.5 / erfc( sqrt(0.5) * x)) \sim log(e^0.5*x^2 ( sqrt(0.125 * pi) * x + sqrt(0.125 * pi) / x - sqrt(pi) / x^3)
// log(0.5 / erfc( sqrt(0.5) * x)) \sim 0.5*x^2  + log(( sqrt(0.125 * pi) * x + sqrt(0.125 * pi) / x - sqrt(pi) / x^3)
// log(2.0 * erfc( sqrt(0.5) * x)) \sim -0.5*x^2  - log(( sqrt(0.125 * pi) * x + sqrt(0.125 * pi) / x - sqrt(pi) / x^3)

// 1 2 10 74 706 8162 110410 1708394

double d_lncumNormal_dx(double x){
    if (x >= -32.0) return exp(-0.5* x*x) / (erfc(-sqrt(0.5) * x) * sqrt(0.5 * M_PI));
    else{
        double x2 = x*x;
        return -( x + (1.0 - (2.0 - (10.0 - (74.0 - (706.0 - (8162.0 - 110410.0 / x2) / x2) / x2)/ x2) / x2) / x2 ) / x);
    }
}
double d2_lncumNormal_dx2(double x){
    if (x >= -32.0) {
        double e = exp(-0.5* x*x) * sqrt(2.0 / M_PI) / erfc(-sqrt(0.5) * x);
        return -(x + e) * e;
    }else{
        double x2 = x*x;
        return -1.0 + (1.0 - (6.0 - (50.0 - (518.0 - (6354.0 - (89782.0 - 1435330.0 / x2) / x2) / x2)/ x2) / x2) / x2 ) / x2;
   }
}

double logP_to_Normal( const double mu, const double sigmasquare, const double logP ){

}
double logP_to_Exponential(const double scale, const double logP){return(-scale * logP);}
double logP_to_Gamma(const double scale, const double shape, const double logP){

}
double HypergeomericLogitPval(double k, uint32_t n, uint32_t K, uint32_t N){ // gets the probability masses that accounts for >= probfraction (or all if no threshold provided)
    uint32_t base = (uint32_t) floor(k +0.5);
    return HypergeomericLogitPval(base, n,K,N, k + 0.5 - base);
}
double HypergeomericLogitPval(uint32_t k, uint32_t n, uint32_t K, uint32_t N, double topfr){ // gets the probability masses that accounts for >= probfraction (or all if no threshold provided)
    uint32_t mode_p = ((n + 1) *(K +1)) / (N +2);
    bool flip;
    if ((flip = (k > mode_p))){k = n - k; K = N - K; topfr = 1.0 - topfr;}
    // now k <= mode_p;
    double logprob = lngammachoose(k, K) + lngammachoose(n - k, N - K) -  lngammachoose(n, N);
    double dblk = (double)k;
    double tmp;
    double frac = 1.0;
    double probsum = topfr;
    do{
        tmp = probsum;
        frac *= dblk;
        frac *= dblk + N - K - n;
        dblk-= 1.0;
        frac /= -dblk + K;
        frac /= -dblk + n;
        probsum += frac;
     //  printf("frac: %e\n", frac);
    }while(tmp != probsum);
    probsum = -log(probsum) - logprob;
return (flip) ? log_e_x_m1(probsum): -log_e_x_m1(probsum);}
// sqrt(2) * InverseErfc(2 / (1+exp(x)))^2
double logitPval_to_Zscore(double x){
    if (fabs(x) > 2.78136714){ // sqrt((-1/2)*(Log[2Pi] - 2*Log[1 + E^x]  + Log[2*Log[1 + E^x] - Log[2Pi]]))
        double term;
        if (x < 0.0){
            term = 2.0 * (-x + log_e_x_p1(x)) - log(M_2PI);
            return -sqrt( term - log(term));
        }else{
            term = 2.0 * (x + log_e_x_p1(-x)) - log(M_2PI);
            return sqrt( term - log(term));
        }
        /*		const double c2 =  1.0/3.0;
		const double c4 =  1.0/5.0;
		const double c5 =  1.0/6.0;
		const double c6 =  1.0/7.0;
		const double c8 =  1.0/9.0;
		const double c9 =  1.0/10.0;
		return t*(1.0 - t * (0.5 - t*(c2 - t*(0.25 - t*(c4 - t*(c5 - t*(c6 - t*(0.125 - t*(c8 - t*c9)))))))));
		*/
    }else{
        const double c1f = sqrt(M_PI * 2.0) * 0.25;
        const double c2 = (-4.0 + M_PI) / 48.0;
        const double c3 = (64.0 - M_PI * (40.0 - 7.0 * M_PI)) / 7680.0;
        const double c4 = (-2176.0 + M_PI * (2464.0 - M_PI * (980.0 - 127.0 * M_PI))) / 2580480.0;
        const double c5 = (126976.0 - M_PI * (225280.0 - M_PI * (150528.0 - M_PI * (42672.0 - 4369.0 * M_PI)))) / 1486356480.0;
//        const double c4 = (-2176.0 + M_PI * (2464.0 - M_PI * (980.0 - 127.0 * M_PI))) / 2580480.0;
//        const double c5 = (126976.0 - M_PI * (225280.0 - M_PI * (150528.0 - M_PI * (42672.0 - 4369.0 * M_PI)))) / 1486356480.0;
        double sqr = x * x;
        return x * c1f * (1.0 + sqr * (c2 + sqr * (c3 + sqr * (c4 + sqr * c5))));
    }
}



//(-1/2)*(Log[2Pi]+ Log[ (2*Log[1 + E^x] - Log[2Pi]) / (1 + E^x)^2])



} // end of namespace

 //   #include "Advanced.cpp"
