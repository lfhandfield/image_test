
// contains modified code from: (optional, needs #define GNU_SCIENTIFIC_LIBRARY in primitive.h)

/* gsl_machine.h
 * specfunc/exp.c
 * specfunc/gamma.c
 * specfunc/gamma_inc.c
 * specfunc/beta.c
 * specfunc/beta_inc.c
 * specfunc/log.c
 * specfunc/log.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
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



#define OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
//#define UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)
#define UNDERFLOW_ERROR(result) /* ignored */

#define INTERNAL_OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; return GSL_EOVRFLW; } while(0)

// #define INTERNAL_UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; return GSL_EUNDRFLW; } while(0)
#define INTERNAL_UNDERFLOW_ERROR(result) /* ignored */

#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)
#define DOMAIN_ERROR_MSG(msg, result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ((msg), GSL_EDOM); } while(0)
#define DOMAIN_ERROR_E10(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; (result)->e10 = 0 ; GSL_ERROR ("domain error", GSL_EDOM); } while(0)
#define OVERFLOW_ERROR_E10(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; (result)->e10 = 0; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
// #define UNDERFLOW_ERROR_E10(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; (result)->e10 = 0; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)
#define UNDERFLOW_ERROR_E10(result) /* ignored */

#define OVERFLOW_ERROR_2(r1,r2) do { (r1)->val = GSL_POSINF; (r1)->err = GSL_POSINF; (r2)->val = GSL_POSINF ; (r2)->err=GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
// #define UNDERFLOW_ERROR_2(r1,r2) do { (r1)->val = 0.0; (r1)->err = GSL_DBL_MIN; (r2)->val = 0.0 ; (r2)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)
#define UNDERFLOW_ERROR_2(r1,r2) /* ignored */
#define DOMAIN_ERROR_2(r1,r2) do { (r1)->val = GSL_NAN; (r1)->err = GSL_NAN;  (r2)->val = GSL_NAN; (r2)->err = GSL_NAN;  GSL_ERROR ("domain error", GSL_EDOM); } while(0)
#define GSL_SF_FACT_NMAX 170
#define GSL_SF_GAMMA_XMAX  171.0
#define GSL_SF_DOUBLEFACT_NMAX 297
// #define CHECK_UNDERFLOW(r) if (fabs((r)->val) < GSL_DBL_MIN) GSL_ERROR("underflow", GSL_EUNDRFLW);
#define CHECK_UNDERFLOW(r) /* ignored */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)
#define GSL_ERROR(reason, gsl_errno) \
do { \
gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
return gsl_errno ; \
} while (0)
/* GSL_ERROR_VAL: call the error handler, and return the given value */
#define GSL_ERROR_VAL(reason, gsl_errno, value) \
do { \
gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
return value ; \
} while (0)
/* GSL_ERROR_VOID: call the error handler, and then return
 (for void functions which still need to generate an error) */
#define GSL_ERROR_VOID(reason, gsl_errno) \
do { \
gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
return ; \
} while (0)
/* GSL_ERROR_NULL suitable for out-of-memory conditions */
#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)
#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02
#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)
#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)
#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01
#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)
/* !MACHINE CONSTANTS! */
/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON
typedef struct
	{
		long double dat[2];
	}
	gsl_complex_long_double;
typedef struct
	{
		double dat[2];
	}
	gsl_complex;
typedef struct
	{
		float dat[2];
	}
	gsl_complex_float;
enum {
	GSL_SUCCESS  = 0,
	GSL_FAILURE  = -1,
	GSL_CONTINUE = -2,  /* iteration has not converged */
	GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
	GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
	GSL_EFAULT   = 3,   /* invalid pointer */
	GSL_EINVAL   = 4,   /* invalid argument supplied by user */
	GSL_EFAILED  = 5,   /* generic failure */
	GSL_EFACTOR  = 6,   /* factorization failed */
	GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
	GSL_ENOMEM   = 8,   /* malloc failed */
	GSL_EBADFUNC = 9,   /* problem with user-supplied function */
	GSL_ERUNAWAY = 10,  /* iterative process is out of control */
	GSL_EMAXITER = 11,  /* exceeded max number of iterations */
	GSL_EZERODIV = 12,  /* tried to divide by zero */
	GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
	GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
	GSL_EUNDRFLW = 15,  /* underflow */
	GSL_EOVRFLW  = 16,  /* overflow  */
	GSL_ELOSS    = 17,  /* loss of accuracy */
	GSL_EROUND   = 18,  /* failed because of roundoff error */
	GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
	GSL_ENOTSQR  = 20,  /* matrix not square */
	GSL_ESING    = 21,  /* apparent singularity detected */
	GSL_EDIVERGE = 22,  /* integral or series is divergent */
	GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
	GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
	GSL_ECACHE   = 25,  /* cache limit exceeded */
	GSL_ETABLE   = 26,  /* table limit exceeded */
	GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
	GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
	GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
	GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
	GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
	GSL_EOF      = 32   /* end of file */
} ;
#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif
#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif
#define GSL_POSZERO (+0)
#define GSL_NEGZERO (-0)
/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */
/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))
void no_error_handler (const char *reason, const char *file, int line, int gsl_errno);
void gsl_error (const char * reason, const char * file, int line, int gsl_errno){

	fprintf (stderr, "%s.\n",reason);
	fflush (stderr);

	abort ();
}
void no_error_handler(const char *reason, const char *file, int line, int gsl_errno){
	/* do nothing */
	reason = 0;
	file = 0;
	line = 0;
	gsl_errno = 0;
	return;
}
struct cheb_series_struct {
	double * c;   /* coefficients                */
	int order;    /* order of expansion          */
	double a;     /* lower interval point        */
	double b;     /* upper interval point        */
	int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;
static inline int	cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result){
	int j;
	double d  = 0.0;
	double dd = 0.0;

	double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
	double y2 = 2.0 * y;

	double e = 0.0;

	for(j = cs->order; j>=1; j--) {
		double temp = d;
		d = y2*d - dd + cs->c[j];
		e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
		dd = temp;
	}

	{
		double temp = d;
		d = y*d - dd + 0.5 * cs->c[0];
		e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
	}

	result->val = d;
	result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

	return GSL_SUCCESS;
}
int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result){
	/* synthetic extended precision constants */
	const double P1 = 4 * 7.8539812564849853515625e-01;
	const double P2 = 4 * 3.7748947079307981766760e-08;
	const double P3 = 4 * 2.6951514290790594840552e-15;
	const double TwoPi = 2*(P1 + P2 + P3);

	const double y = GSL_SIGN(theta) * 2 * floor(fabs(theta)/TwoPi);
	double r = ((theta - y*P1) - y*P2) - y*P3;

	if(r >  M_PI) { r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
	else if (r < -M_PI) r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */

	result->val = r;

	if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
		result->val = GSL_NAN;
		result->err = GSL_NAN;
		GSL_ERROR ("error", GSL_ELOSS);
	}
	else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val - theta);
		return GSL_SUCCESS;
	}
	else {
		double delta = fabs(result->val - theta);
		result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
		return GSL_SUCCESS;
	}
}

double log_1plusx_mx_e(double x){
	if(x <= -1.0) myexit("invalid Domain for log_1plusx_mx_e\n");
	else if(fabs(x) < GSL_ROOT5_DBL_EPSILON) {
		const double c1 = -0.5;
		const double c2 =  1.0/3.0;
		const double c3 = -1.0/4.0;
		const double c4 =  1.0/5.0;
		const double c5 = -1.0/6.0;
		const double c6 =  1.0/7.0;
		const double c7 = -1.0/8.0;
		const double c8 =  1.0/9.0;
		const double c9 = -1.0/10.0;
		const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
		return x*x * (c1 + x*(c2 + x*(c3 + x*(c4 + x*t))));
	}
	else if(fabs(x) < 0.5) {
		const double logdata[] = {
			-1.12100231323744103373737274541,
			0.19553462773379386241549597019,
			-0.01467470453808083971825344956,
			0.00166678250474365477643629067,
			-0.00018543356147700369785746902,
			0.00002280154021771635036301071,
			-2.8031253116633521699214134172e-06,
			3.5936568872522162983669541401e-07,
			-4.6241857041062060284381167925e-08,
			6.0822637459403991012451054971e-09,
			-8.0339824424815790302621320732e-10,
			1.0751718277499375044851551587e-10,
			-1.4445310914224613448759230882e-11,
			1.9573912180610336168921438426e-12,
			-2.6614436796793061741564104510e-13,
			3.6402634315269586532158344584e-14,
			-4.9937495922755006545809120531e-15,
			6.8802890218846809524646902703e-16,
			-9.5034129794804273611403251480e-17,
			1.3170135013050997157326965813e-17
		};
		const Tuple<double, 20u> log_cheb = Tuple<double, 20>(logdata);

		return x*x * cheb_eval(log_cheb, 0.5*(8.0*x + 1.0)/(x+2.0));
	}
	else return log(1.0 + x) - x;
}

static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.04955262679620434040357683080,
  0.00449293488768382749558001242,
 -0.00129194104658496953494224761,
 -0.00001836389292149396270416979,
  0.00002211114704099526291538556,
 -5.23337485234257134673693179020e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.85627169231704599442806370690e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.29599561220523396007359328540e-19
};
static cheb_series erfc_xlt1_cs = {erfc_xlt1_data, 19, -1, 1, 12};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.013343124200271211203618353102,
  0.003824682739750469767692372556,
 -0.001058699227195126547306482530,
  0.000283859419210073742736310108,
 -0.000073906170662206760483959432,
  0.000018725312521489179015872934,
 -4.62530981164919445131297264430e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.07462122724551777372119408710e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.65174789720310713757307724790e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.22171792645348091472914001250e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.06258237507038698392265499770e-14,
  9.89705409478327321641264227110e-15,
 -1.90685978789192181051961024995e-15,
  3.50826648032737849245113757340e-16
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.003736240359381998520654927536,
 -0.000916623948045470238763619870,
  0.000199094325044940833965078819,
 -0.000040276384918650072591781859,
  7.76515264697061049477127605790e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.61833026634844152345304095560e-8,
  8.00253111512943601598732144340e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.78022521563251805044056974560e-11,
  6.17253683874528285729910462130e-12,
 -9.96019290955316888445830597430e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.92607828989125810013581287560e-15,
 -6.07970619384160374392535453420e-16,
  9.12600607264794717315507477670e-17
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

static double erfc8_sum(double x)
{
	/* estimates erfc(x) valid for 8 < x < 100 */
	/* This is based on index 5725 in Hart et al */

	static double P[] = {
		2.97886562639399288862,
		7.409740605964741794425,
		6.1602098531096305440906,
		5.019049726784267463450058,
		1.275366644729965952479585264,
		0.5641895835477550741253201704
	};
	static double Q[] = {
		3.3690752069827527677,
		9.608965327192787870698,
		17.08144074746600431571095,
		12.0489519278551290360340491,
		9.396034016235054150430579648,
		2.260528520767326969591866945,
		1.0
	};
	double num=0.0, den=0.0;
	int i;

	num = P[5];
	for (i=4; i>=0; --i) {
		num = x*num + P[i];
	}
	den = Q[6];
	for (i=5; i>=0; --i) {
		den = x*den + Q[i];
	}

	return num/den;
}

inline static double erfc8(double x){double e;e = erfc8_sum(x);e *= exp(-x*x);return e;}
inline static double log_erfc8(double x){double e; e = erfc8_sum(x); e = log(e) - x*x; return e;}
double erfc(double x)
{
	const double ax = fabs(x);
	double e_val;

	/* CHECK_POINTER(result) */

	if(ax <= 1.0) {
		double t = 2.0*ax - 1.0;
		static double erfc_xlt1_data[20] = {
			1.06073416421769980345174155056,
			-0.42582445804381043569204735291,
			0.04955262679620434040357683080,
			0.00449293488768382749558001242,
			-0.00129194104658496953494224761,
			-0.00001836389292149396270416979,
			0.00002211114704099526291538556,
			-5.23337485234257134673693179020e-7,
			-2.78184788833537885382530989578e-7,
			1.41158092748813114560316684249e-8,
			2.72571296330561699984539141865e-9,
			-2.06343904872070629406401492476e-10,
			-2.14273991996785367924201401812e-11,
			2.22990255539358204580285098119e-12,
			1.36250074650698280575807934155e-13,
			-1.95144010922293091898995913038e-14,
			-6.85627169231704599442806370690e-16,
			1.44506492869699938239521607493e-16,
			2.45935306460536488037576200030e-18,
			-9.29599561220523396007359328540e-19
		};
		e_val = cheb_eval(Tuple<double,20>(erfc_xlt1_data) , t);
	}
	else if(ax <= 5.0) {
		double ex2 = exp(-x*x);
		double t = 0.5*(ax-3.0);
		static double erfc_x15_data[25] = {
			0.44045832024338111077637466616,
			-0.143958836762168335790826895326,
			0.044786499817939267247056666937,
			-0.013343124200271211203618353102,
			0.003824682739750469767692372556,
			-0.001058699227195126547306482530,
			0.000283859419210073742736310108,
			-0.000073906170662206760483959432,
			0.000018725312521489179015872934,
			-4.62530981164919445131297264430e-6,
			1.11558657244432857487884006422e-6,
			-2.63098662650834130067808832725e-7,
			6.07462122724551777372119408710e-8,
			-1.37460865539865444777251011793e-8,
			3.05157051905475145520096717210e-9,
			-6.65174789720310713757307724790e-10,
			1.42483346273207784489792999706e-10,
			-3.00141127395323902092018744545e-11,
			6.22171792645348091472914001250e-12,
			-1.26994639225668496876152836555e-12,
			2.55385883033257575402681845385e-13,
			-5.06258237507038698392265499770e-14,
			9.89705409478327321641264227110e-15,
			-1.90685978789192181051961024995e-15,
			3.50826648032737849245113757340e-16
		};
		e_val = ex2 * cheb_eval(Tuple<double,25>(erfc_x15_data), t);
	}
	else if(ax < 10.0) {
		double exterm = exp(-x*x) / ax;
		double t = (2.0*ax - 15.0)/5.0;

		static double erfc_x510_data[20] = {
			1.11684990123545698684297865808,
			0.003736240359381998520654927536,
			-0.000916623948045470238763619870,
			0.000199094325044940833965078819,
			-0.000040276384918650072591781859,
			7.76515264697061049477127605790e-6,
			-1.44464794206689070402099225301e-6,
			2.61311930343463958393485241947e-7,
			-4.61833026634844152345304095560e-8,
			8.00253111512943601598732144340e-9,
			-1.36291114862793031395712122089e-9,
			2.28570483090160869607683087722e-10,
			-3.78022521563251805044056974560e-11,
			6.17253683874528285729910462130e-12,
			-9.96019290955316888445830597430e-13,
			1.58953143706980770269506726000e-13,
			-2.51045971047162509999527428316e-14,
			3.92607828989125810013581287560e-15,
			-6.07970619384160374392535453420e-16,
			9.12600607264794717315507477670e-17
		};

		e_val = exterm * cheb_eval(Tuple<double,20>(erfc_x510_data) , t);
	}
	else {
		e_val = erfc8(ax);
	}

	return (x < 0.0) ? 2.0 - e_val : e_val;
}

int gsl_sf_erfc_e(double x, gsl_sf_result * result)
{
  const double ax = fabs(x);
  double e_val, e_err;

  /* CHECK_POINTER(result) */

  if(ax <= 1.0) {
    double t = 2.0*ax - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&erfc_xlt1_cs, t, &c);
    e_val = c.val;
    e_err = c.err;
  }
  else if(ax <= 5.0) {
    double ex2 = exp(-x*x);
    double t = 0.5*(ax-3.0);
    gsl_sf_result c;
    cheb_eval_e(&erfc_x15_cs, t, &c);
    e_val = ex2 * c.val;
    e_err = ex2 * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON);
  }
  else if(ax < 10.0) {
    double exterm = exp(-x*x) / ax;
    double t = (2.0*ax - 15.0)/5.0;
    gsl_sf_result c;
    cheb_eval_e(&erfc_x510_cs, t, &c);
    e_val = exterm * c.val;
    e_err = exterm * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON + GSL_DBL_EPSILON);
  }
  else {
    e_val = erfc8(ax);
    e_err = (x*x + 1.0) * GSL_DBL_EPSILON * fabs(e_val);
  }

  if(x < 0.0) {
    result->val  = 2.0 - e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    result->val  = e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}

int gsl_sf_log_erfc_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x*x < 10.0*GSL_ROOT6_DBL_EPSILON) {
    const double y = x / M_SQRTPI;
    /* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
    const double c3 = (4.0 - M_PI)/3.0;
    const double c4 = 2.0*(1.0 - M_PI/3.0);
    const double c5 = -0.001829764677455021;  /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
    const double c6 =  0.02629651521057465;   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
    const double c7 = -0.01621575378835404;
    const double c8 =  0.00125993961762116;
    const double c9 =  0.00556964649138;
    const double c10 = -0.0045563339802;
    const double c11 =  0.0009461589032;
    const double c12 =  0.0013200243174;
    const double c13 = -0.00142906;
    const double c14 =  0.00048204;
    double series = c8 + y*(c9 + y*(c10 + y*(c11 + y*(c12 + y*(c13 + c14*y)))));
    series = y*(1.0 + y*(1.0 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*series)))))));
    result->val = -2.0 * series;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }else if(x > 8.0) {
    result->val = log_erfc8(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }else {
    gsl_sf_result result_erfc;
    gsl_sf_erfc_e(x, &result_erfc);
    result->val  = log(result_erfc.val);
    result->err  = fabs(result_erfc.err / result_erfc.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

double gsl_sf_log_erfc(double x){gsl_sf_result result; gsl_sf_log_erfc_e(x, &result); return result.val;}



double polygamma0(double xx){
	return( 32768.0f * (lngamma(xx * (1.0f + (1.0f /65536.0f))) - lngamma(xx * (1.0f - (1.0f /65536.0f)))) / xx );
}

double gammastar(const double x)
{
	if(x <= 0.0) myexit("invalid Domain for gammastar\n");

	else if(x < 0.5) {
		double lg = lngamma(x);
		const double lx = log(x);
		const double c  = 0.5*(log(2.0f*M_PI));
		return lg - (x-0.5)*lx + x - c;
	}
	else if(x < 2.0) {
		const double t = 4.0/3.0*(x-0.5) - 1.0;

		const double gstar_a_data[30] = {
			2.16786447866463034423060819465,
			-0.05533249018745584258035832802,
			0.01800392431460719960888319748,
			-0.00580919269468937714480019814,
			0.00186523689488400339978881560,
			-0.00059746524113955531852595159,
			0.00019125169907783353925426722,
			-0.00006124996546944685735909697,
			0.00001963889633130842586440945,
			-6.3067741254637180272515795142e-06,
			2.0288698405861392526872789863e-06,
			-6.5384896660838465981983750582e-07,
			2.1108698058908865476480734911e-07,
			-6.8260714912274941677892994580e-08,
			2.2108560875880560555583978510e-08,
			-7.1710331930255456643627187187e-09,
			2.3290892983985406754602564745e-09,
			-7.5740371598505586754890405359e-10,
			2.4658267222594334398525312084e-10,
			-8.0362243171659883803428749516e-11,
			2.6215616826341594653521346229e-11,
			-8.5596155025948750540420068109e-12,
			2.7970831499487963614315315444e-12,
			-9.1471771211886202805502562414e-13,
			2.9934720198063397094916415927e-13,
			-9.8026575909753445931073620469e-14,
			3.2116773667767153777571410671e-14,
			-1.0518035333878147029650507254e-14,
			3.4144405720185253938994854173e-15,
			-1.0115153943081187052322643819e-15
		};
		const Tuple<double, 30> log_cheb_a = Tuple<double, 30>(gstar_a_data);
		return cheb_eval(log_cheb_a, t);
	}
	else if(x < 10.0) {
		const double t = 0.25*(x-2.0) - 1.0;

		static double gstar_b_data[] = {
			0.0057502277273114339831606096782,
			0.0004496689534965685038254147807,
			-0.0001672763153188717308905047405,
			0.0000615137014913154794776670946,
			-0.0000223726551711525016380862195,
			8.0507405356647954540694800545e-06,
			-2.8671077107583395569766746448e-06,
			1.0106727053742747568362254106e-06,
			-3.5265558477595061262310873482e-07,
			1.2179216046419401193247254591e-07,
			-4.1619640180795366971160162267e-08,
			1.4066283500795206892487241294e-08,
			-4.6982570380537099016106141654e-09,
			1.5491248664620612686423108936e-09,
			-5.0340936319394885789686867772e-10,
			1.6084448673736032249959475006e-10,
			-5.0349733196835456497619787559e-11,
			1.5357154939762136997591808461e-11,
			-4.5233809655775649997667176224e-12,
			1.2664429179254447281068538964e-12,
			-3.2648287937449326771785041692e-13,
			7.1528272726086133795579071407e-14,
			-9.4831735252566034505739531258e-15,
			-2.3124001991413207293120906691e-15,
			2.8406613277170391482590129474e-15,
			-1.7245370321618816421281770927e-15,
			8.6507923128671112154695006592e-16,
			-3.9506563665427555895391869919e-16,
			1.6779342132074761078792361165e-16,
			-6.0483153034414765129837716260e-17
		};
		const Tuple<double, 30> log_cheb_b = Tuple<double, 30>(gstar_b_data);
		double c = cheb_eval(log_cheb_b, t);
		return c/(x*x) + 1.0 + 1.0/(12.0*x);

	}
	else if(x < 1.0/GSL_ROOT4_DBL_EPSILON) {
		const double y = 1.0/(x*x);
		const double c0 =  1.0/12.0;
		const double c1 = -1.0/360.0;
		const double c2 =  1.0/1260.0;
		const double c3 = -1.0/1680.0;
		const double c4 =  1.0/1188.0;
		const double c5 = -691.0/360360.0;
		const double c6 =  1.0/156.0;
		const double c7 = -3617.0/122400.0;
		const double ser = c0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
		return exp(ser/x);
	}
	else if(x < 1.0/ExCo<double>::epsilon()) {
		/* Use Stirling formula for Gamma(x).
		 */
		const double xi = 1.0/x;
		return(1.0 + xi/12.0*(1.0 + xi/24.0*(1.0 - xi*(139.0/180.0 + 571.0/8640.0*xi))));
	}
	else return(1.0);
}

double incgamma_dom(const double a, const double x) {
	if(a < 10.0) {
		double lnr;
		double lg;
		lg = lngamma(a+1.0);
		lnr = a * log(x) - x - lg;
		return exp(lnr);
	}else{
		double gstar;
		double ln_term;
		double term1;
		if (x < 0.5*a) {
			double u = x/a;
			double ln_u = log(u);
			ln_term = ln_u - u + 1.0;
		} else {
			double mu = (x-a)/a;
			ln_term = log_1plusx_mx_e(mu);  /* log(1+mu) - mu */
		};
		gstar = gammastar(a);
		term1 = exp(a*ln_term)/sqrt(2.0*M_PI*a);
		return  term1/gstar;
	}
}

double
gamma_inc_P_series(const double a, const double x)
{
	const int nmax = 5000;
	double dd = incgamma_dom(a, x);
	double sum  = 1.0;
	double term = 1.0;
	int n;
	for(n=1; n<nmax; n++) {
		term *= x/(a+n);
		sum  += term;
		if(fabs(term/sum) < ExCo<double>::epsilon()) break;
	}
	return dd * sum;
}

double
gamma_inc_F_CF(const double a, const double x)
{
	const int    nmax  =  5000;
	const double smallr =  pow(ExCo<double>::epsilon(), 3.0f);

	double hn = 1.0;           /* convergent */
	double Cn = 1.0 / smallr;
	double Dn = 1.0;
	int n;

	/* n == 1 has a_1, b_1, b_0 independent of a,x,
	 so that has been done by hand                */
	for ( n = 2 ; n < nmax ; n++ ) {
		double an;
		double delta;

		if (n & 1)
			an = 0.5*(n-1)/x;
		else
			an = (0.5*n-a)/x;

		Dn = 1.0 + an * Dn;
		if ( fabs(Dn) < smallr )
			Dn = smallr;
		Cn = 1.0 + an/Cn;
		if ( fabs(Cn) < smallr )
			Cn = smallr;
		Dn = 1.0 / Dn;
		delta = Cn * Dn;
		hn *= delta;
		if(fabs(delta-1.0) < ExCo<double>::epsilon()) break;
	}
	return hn;
}


double gamma_inc_Q_CF(const double a, const double x){
	double D = incgamma_dom(a,x);
	double F = gamma_inc_F_CF(a, x);
	return D * (a/x) * F;
}

double
gamma_inc_Q_series(const double a, const double x)
{
	double term1;  /* 1 - x^a/Gamma(a+1) */
	double sum;    /* 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ... */
	double term2;  /* a temporary variable used at the end */

	{
		/* Evaluate series for 1 - x^a/Gamma(a+1), small a */
		const double pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
		const double lnx  = log(x);
		const double el   = M_EULER+lnx;
		const double c1 = -el;
		const double c2 = M_PI*M_PI/12.0 - 0.5*el*el;
		const double c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
		const double c4 = -0.04166666666666666667
		* (-1.758243446661483480 + lnx)
		* (-0.764428657272716373 + lnx)
		* ( 0.723980571623507657 + lnx)
		* ( 4.107554191916823640 + lnx);
		const double c5 = -0.0083333333333333333
		* (-2.06563396085715900 + lnx)
		* (-1.28459889470864700 + lnx)
		* (-0.27583535756454143 + lnx)
		* ( 1.33677371336239618 + lnx)
		* ( 5.17537282427561550 + lnx);
		const double c6 = -0.0013888888888888889
		* (-2.30814336454783200 + lnx)
		* (-1.65846557706987300 + lnx)
		* (-0.88768082560020400 + lnx)
		* ( 0.17043847751371778 + lnx)
		* ( 1.92135970115863890 + lnx)
		* ( 6.22578557795474900 + lnx);
		const double c7 = -0.00019841269841269841
		* (-2.5078657901291800 + lnx)
		* (-1.9478900888958200 + lnx)
		* (-1.3194837322612730 + lnx)
		* (-0.5281322700249279 + lnx)
		* ( 0.5913834939078759 + lnx)
		* ( 2.4876819633378140 + lnx)
		* ( 7.2648160783762400 + lnx);
		const double c8 = -0.00002480158730158730
		* (-2.677341544966400 + lnx)
		* (-2.182810448271700 + lnx)
		* (-1.649350342277400 + lnx)
		* (-1.014099048290790 + lnx)
		* (-0.191366955370652 + lnx)
		* ( 0.995403817918724 + lnx)
		* ( 3.041323283529310 + lnx)
		* ( 8.295966556941250 + lnx);
		const double c9 = -2.75573192239859e-6
		* (-2.8243487670469080 + lnx)
		* (-2.3798494322701120 + lnx)
		* (-1.9143674728689960 + lnx)
		* (-1.3814529102920370 + lnx)
		* (-0.7294312810261694 + lnx)
		* ( 0.1299079285269565 + lnx)
		* ( 1.3873333251885240 + lnx)
		* ( 3.5857258865210760 + lnx)
		* ( 9.3214237073814600 + lnx);
		const double c10 = -2.75573192239859e-7
		* (-2.9540329644556910 + lnx)
		* (-2.5491366926991850 + lnx)
		* (-2.1348279229279880 + lnx)
		* (-1.6741881076349450 + lnx)
		* (-1.1325949616098420 + lnx)
		* (-0.4590034650618494 + lnx)
		* ( 0.4399352987435699 + lnx)
		* ( 1.7702236517651670 + lnx)
		* ( 4.1231539047474080 + lnx)
		* ( 10.342627908148680 + lnx);

		term1 = a*(c1+a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))));
	}

	{
		// Evaluate the sum.
		const int nmax = 5000;
		double t = 1.0;
		int n;
		sum = 1.0;

		for(n=1; n<nmax; n++) {
			t *= -x/(n+1.0);
			sum += (a+1.0)/(a+n+1.0)*t;
			if(fabs(t/sum) < ExCo<double>::epsilon()) break;
		}
	}

	term2 = (1.0 - term1) * a/(a+1.0) * x * sum;
	return term1 + term2;
}

double gamma_inc_Q_asymp_unif(const double a, const double x) {
	const double rta = sqrt(a);
	const double eps = (x-a)/a;
	double ln_term = log_1plusx_mx_e(eps);  /* log(1+eps) - eps */
	const double eta  = ((eps) >= 0.0 ? 1 : -1) * sqrt(-2.0*ln_term);
	double R;
	double c0, c1;
	if(fabs(eps) < GSL_ROOT5_DBL_EPSILON) {
		c0 = -1.0/3.0 + eps*(1.0/12.0 - eps*(23.0/540.0 - eps*(353.0/12960.0 - eps*589.0/30240.0)));
		c1 = -1.0/540.0 - eps/288.0;
	}
	else {
		const double rt_term = sqrt(-2.0 * ln_term/(eps*eps));
		const double lam = x/a;
		c0 = (1.0 - 1.0/rt_term)/eps;
		c1 = -(eta*eta*eta * (lam*lam + 10.0*lam + 1.0) - 12.0 * eps*eps*eps) / (12.0 * eta*eta*eta*eps*eps*eps);
	}
	R = exp(-0.5*a*eta*eta)* M_SQRT_1_2 /(sqrt(M_PI)*rta) * (c0 + c1/a);
	return 0.5 * erfc(eta*rta * M_SQRT_1_2) + R;
}

/* Q large x asymptotic */
double
gamma_inc_Q_large_x(const double a, const double x)
{
	const int nmax = 5000;
	double D = incgamma_dom(a,x);
	double sum  = 1.0;
	double term = 1.0;
	double last = 1.0;
	int n;
	for(n=1; n<nmax; n++) {
		term *= (a-n)/x;
		if(fabs(term/last) > 1.0) break;
		if(fabs(term/sum)  < ExCo<double>::epsilon()) break;
		sum  += term;
		last  = term;
	}
	return D * (a/x) * sum;
}


/* Evaluate the continued fraction for exprel.
 * [Abramowitz+Stegun, 4.2.41]
 */
static
int	exprel_n_CF(const int N, const double x, gsl_sf_result * result){
	const double RECUR_BIG = GSL_SQRT_DBL_MAX;
	const int maxiter = 5000;
	int n = 1;
	double Anm2 = 1.0;
	double Bnm2 = 0.0;
	double Anm1 = 0.0;
	double Bnm1 = 1.0;
	double a1 = 1.0;
	double b1 = 1.0;
	double a2 = -x;
	double b2 = N+1;
	double an, bn;

	double fn;

	double An = b1*Anm1 + a1*Anm2;   /* A1 */
	double Bn = b1*Bnm1 + a1*Bnm2;   /* B1 */

	/* One explicit step, before we get to the main pattern. */
	n++;
	Anm2 = Anm1;
	Bnm2 = Bnm1;
	Anm1 = An;
	Bnm1 = Bn;
	An = b2*Anm1 + a2*Anm2;   /* A2 */
	Bn = b2*Bnm1 + a2*Bnm2;   /* B2 */

	fn = An/Bn;

	while(n < maxiter) {
		double old_fn;
		double del;
		n++;
		Anm2 = Anm1;
		Bnm2 = Bnm1;
		Anm1 = An;
		Bnm1 = Bn;
		an = ( (n & 1) ? ((n-1)/2)*x : -(N+(n/2)-1)*x );
		bn = N + n - 1;
		An = bn*Anm1 + an*Anm2;
		Bn = bn*Bnm1 + an*Bnm2;

		if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
			An /= RECUR_BIG;
			Bn /= RECUR_BIG;
			Anm1 /= RECUR_BIG;
			Bnm1 /= RECUR_BIG;
			Anm2 /= RECUR_BIG;
			Bnm2 /= RECUR_BIG;
		}

		old_fn = fn;
		fn = An/Bn;
		del = old_fn/fn;

		if(fabs(del - 1.0) < 2.0*GSL_DBL_EPSILON) break;
	}

	result->val = fn;
	result->err = 2.0*(n+1.0)*GSL_DBL_EPSILON*fabs(fn);

	if(n == maxiter)
		GSL_ERROR ("error", GSL_EMAXITER);
	else
		return GSL_SUCCESS;
}
/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
int	gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result){
	const double adx = fabs(dx);


	if(x + adx > GSL_LOG_DBL_MAX) {
		OVERFLOW_ERROR(result);
	}
	else if(x - adx < GSL_LOG_DBL_MIN) {
		UNDERFLOW_ERROR(result);
	}
	else {
		const double ex  = exp(x);
		const double edx = exp(adx);
		result->val  = ex;
		result->err  = ex * ((GSL_DBL_EPSILON > edx - 1.0/edx) ? GSL_DBL_EPSILON  : edx - 1.0/edx);
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
return GSL_FAILURE;
}
int gsl_sf_exp_e(const double x, gsl_sf_result * result){
	if(x > GSL_LOG_DBL_MAX) {
		OVERFLOW_ERROR(result);
	}
	else if(x < GSL_LOG_DBL_MIN) {
		UNDERFLOW_ERROR(result);
	}
	else {
		result->val = exp(x);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
return GSL_FAILURE;
}
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result){
	if(x > ExCo<int>::mkMaximum()-1) {
		OVERFLOW_ERROR_E10(result);
	}
	else if(x < ExCo<int>::mkMinimum()+1) {
		UNDERFLOW_ERROR_E10(result);
	}
	else {
		const int N = (int) floor(x/M_LN10);
		result->val = exp(x-N*M_LN10);
		result->err = 2.0 * (fabs(x)+1.0) * GSL_DBL_EPSILON * fabs(result->val);
		result->e10 = N;
		return GSL_SUCCESS;
	}
return GSL_FAILURE;
}
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result){
	const double ay  = fabs(y);

	if(y == 0.0) {
		result->val = 0.0;
		result->err = 0.0;
		return GSL_SUCCESS;
	}
	else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
			&& (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
			) {
		const double ex = exp(x);
		result->val = y * ex;
		result->err = (2.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		const double ly  = log(ay);
		const double lnr = x + ly;

		if(lnr > GSL_LOG_DBL_MAX - 0.01) {
			OVERFLOW_ERROR(result);
		}
		else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
			UNDERFLOW_ERROR(result);
		}
		else {
			const double sy   = GSL_SIGN(y);
			const double M    = floor(x);
			const double N    = floor(ly);
			const double a    = x  - M;
			const double b    = ly - N;
			const double berr = 2.0 * GSL_DBL_EPSILON * (fabs(ly) + fabs(N));
			result->val  = sy * exp(M+N) * exp(a+b);
			result->err  = berr * fabs(result->val);
			result->err += 2.0 * GSL_DBL_EPSILON * (M + N + 1.0) * fabs(result->val);
			return GSL_SUCCESS;
		}
	}
return GSL_FAILURE;
}
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result){
	const double ay  = fabs(y);

	if(y == 0.0) {
		result->val = 0.0;
		result->err = 0.0;
		result->e10 = 0;
		return GSL_SUCCESS;
	}
	else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
			&& (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
			) {
		const double ex = exp(x);
		result->val = y * ex;
		result->err = (2.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
		result->e10 = 0;
		return GSL_SUCCESS;
	}
	else {
		const double ly  = log(ay);
		const double l10_val = (x + ly)/M_LN10;

		if(l10_val > ExCo<int>::mkMaximum()-1) {
			OVERFLOW_ERROR_E10(result);
		}
		else if(l10_val < ExCo<int>::mkMinimum()+1) {
			UNDERFLOW_ERROR_E10(result);
		}
		else {
			const double sy  = GSL_SIGN(y);
			const int    N   = (int) floor(l10_val);
			const double arg_val = (l10_val - N) * M_LN10;
			const double arg_err = 2.0 * GSL_DBL_EPSILON * fabs(ly);

			result->val  = sy * exp(arg_val);
			result->err  = arg_err * fabs(result->val);
			result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
			result->e10 = N;

			return GSL_SUCCESS;
		}
	}
return GSL_FAILURE;
}
int gsl_sf_exp_mult_err_e(const double x, const double dx,
						  const double y, const double dy,
						  gsl_sf_result * result){
	const double ay  = fabs(y);

	if(y == 0.0) {
		result->val = 0.0;
		result->err = fabs(dy * exp(x));
		return GSL_SUCCESS;
	}
	else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
			&& (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
			) {
		double ex = exp(x);
		result->val  = y * ex;
		result->err  = ex * (fabs(dy) + fabs(y*dx));
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		const double ly  = log(ay);
		const double lnr = x + ly;

		if(lnr > GSL_LOG_DBL_MAX - 0.01) {
			OVERFLOW_ERROR(result);
		}
		else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
			UNDERFLOW_ERROR(result);
		}
		else {
			const double sy  = GSL_SIGN(y);
			const double M   = floor(x);
			const double N   = floor(ly);
			const double a   = x  - M;
			const double b   = ly - N;
			const double eMN = exp(M+N);
			const double eab = exp(a+b);
			result->val  = sy * eMN * eab;
			result->err  = eMN * eab * 2.0*GSL_DBL_EPSILON;
			result->err += eMN * eab * fabs(dy/y);
			result->err += eMN * eab * fabs(dx);
			return GSL_SUCCESS;
		}
	}
return GSL_FAILURE;
}
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx,
							  const double y, const double dy,
							  gsl_sf_result_e10 * result){
	const double ay  = fabs(y);

	if(y == 0.0) {
		result->val = 0.0;
		result->err = fabs(dy * exp(x));
		result->e10 = 0;
		return GSL_SUCCESS;
	}
	else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
			&& (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
			) {
		const double ex = exp(x);
		result->val  = y * ex;
		result->err  = ex * (fabs(dy) + fabs(y*dx));
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		result->e10 = 0;
		return GSL_SUCCESS;
	}
	else {
		const double ly  = log(ay);
		const double l10_val = (x + ly)/M_LN10;

		if(l10_val > ExCo<int>::mkMaximum()-1) {
			OVERFLOW_ERROR_E10(result);
		}
		else if(l10_val < ExCo<int>::mkMinimum()+1) {
			UNDERFLOW_ERROR_E10(result);
		}
		else {
			const double sy  = GSL_SIGN(y);
			const int    N   = (int) floor(l10_val);
			const double arg_val = (l10_val - N) * M_LN10;
			const double arg_err = dy/fabs(y) + dx + 2.0*GSL_DBL_EPSILON*fabs(arg_val);

			result->val  = sy * exp(arg_val);
			result->err  = arg_err * fabs(result->val);
			result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
			result->e10 = N;

			return GSL_SUCCESS;
		}
	}
return GSL_FAILURE;
}

int gsl_sf_expm1_e(const double x, gsl_sf_result * result){
	const double cut = 0.002;

	if(x < GSL_LOG_DBL_MIN) {
		result->val = -1.0;
		result->err = GSL_DBL_EPSILON;
		return GSL_SUCCESS;
	}
	else if(x < -cut) {
		result->val = exp(x) - 1.0;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < cut) {
		result->val = x * (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < GSL_LOG_DBL_MAX) {
		result->val = exp(x) - 1.0;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		OVERFLOW_ERROR(result);
	}
}
double gsl_sf_expm1 (double x){gsl_sf_result_struct res; gsl_sf_expm1_e(x,&res); return res.val;}
int gsl_sf_exprel_e(const double x, gsl_sf_result * result){
	const double cut = 0.002;

	if(x < GSL_LOG_DBL_MIN) {
		result->val = -1.0/x;
		result->err = GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < -cut) {
		result->val = (exp(x) - 1.0)/x;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < cut) {
		result->val = (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < GSL_LOG_DBL_MAX) {
		result->val = (exp(x) - 1.0)/x;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		OVERFLOW_ERROR(result);
	}
}
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result){
	const double cut = 0.002;

	if(x < GSL_LOG_DBL_MIN) {
		result->val = -2.0/x*(1.0 + 1.0/x);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < -cut) {
		result->val = 2.0*(exp(x) - 1.0 - x)/(x*x);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < cut) {
		result->val = (1.0 + 1.0/3.0*x*(1.0 + 0.25*x*(1.0 + 0.2*x*(1.0 + 1.0/6.0*x))));
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < GSL_LOG_DBL_MAX) {
		result->val = 2.0*(exp(x) - 1.0 - x)/(x*x);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		OVERFLOW_ERROR(result);
	}
}
/* Chebyshev coefficients for
 * x^2(Gamma*(x) - 1 - 1/(12x)), x = 4(t+1)+2, -1 < t < 1
 */
static double gstar_b_data[] = {
0.0057502277273114339831606096782,
0.0004496689534965685038254147807,
-0.0001672763153188717308905047405,
0.0000615137014913154794776670946,
-0.0000223726551711525016380862195,
8.0507405356647954540694800545e-06,
-2.8671077107583395569766746448e-06,
1.0106727053742747568362254106e-06,
-3.5265558477595061262310873482e-07,
1.2179216046419401193247254591e-07,
-4.1619640180795366971160162267e-08,
1.4066283500795206892487241294e-08,
-4.6982570380537099016106141654e-09,
1.5491248664620612686423108936e-09,
-5.0340936319394885789686867772e-10,
1.6084448673736032249959475006e-10,
-5.0349733196835456497619787559e-11,
1.5357154939762136997591808461e-11,
-4.5233809655775649997667176224e-12,
1.2664429179254447281068538964e-12,
-3.2648287937449326771785041692e-13,
7.1528272726086133795579071407e-14,
-9.4831735252566034505739531258e-15,
-2.3124001991413207293120906691e-15,
2.8406613277170391482590129474e-15,
-1.7245370321618816421281770927e-15,
8.6507923128671112154695006592e-16,
-3.9506563665427555895391869919e-16,
1.6779342132074761078792361165e-16,
-6.0483153034414765129837716260e-17
};
static cheb_series gstar_b_cs = {
gstar_b_data,
29,
-1, 1,
18
};
/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
0.99999999999980993227684700473478,
676.520368121885098567009190444019,
-1259.13921672240287047156078755283,
771.3234287776530788486528258894,
-176.61502916214059906584551354,
12.507343278686904814458936853,
-0.13857109526572011689554707,
9.984369578019570859563e-6,
1.50563273514931155834e-7
};
int	gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta){

	if(zr != 0.0 || zi != 0.0) {
		const double ax = fabs(zr);
		const double ay = fabs(zi);
		const double min = (ax< ay) ? ax : ay;
		const double max = (ax> ay) ? ax : ay;
		lnr->val = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
		lnr->err = 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
		theta->val = atan2(zi, zr);
		theta->err = GSL_DBL_EPSILON * fabs(lnr->val);
		return GSL_SUCCESS;
	}
	else {
		DOMAIN_ERROR_2(lnr, theta);
	}
return GSL_FAILURE;
}
/* complex version of Lanczos method; this is not safe for export
 * since it becomes bad in the left half-plane
 */
int	lngamma_lanczos_complex(double zr, double zi, gsl_sf_result * yr, gsl_sf_result * yi){
	int k;
	gsl_sf_result log1_r,    log1_i;
	gsl_sf_result logAg_r,   logAg_i;
	double Ag_r, Ag_i;
	double yi_tmp_val, yi_tmp_err;

	zr -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

	Ag_r = lanczos_7_c[0];
	Ag_i = 0.0;
	for(k=1; k<=8; k++) {
		double R = zr + k;
		double I = zi;
		double a = lanczos_7_c[k] / (R*R + I*I);
		Ag_r +=  a * R;
		Ag_i -=  a * I;
	}

	gsl_sf_complex_log_e(zr + 7.5, zi, &log1_r,  &log1_i);
	gsl_sf_complex_log_e(Ag_r, Ag_i,   &logAg_r, &logAg_i);

	/* (z+0.5)*log(z+7.5) - (z+7.5) + 0.5f * log(2.0 * M_PI) + log(Ag(z)) */
	yr->val = (zr+0.5)*log1_r.val - zi*log1_i.val - (zr+7.5) + 0.5f * log(2.0 * M_PI) + logAg_r.val;
	yi->val = zi*log1_r.val + (zr+0.5)*log1_i.val - zi + logAg_i.val;
	yr->err = 4.0 * GSL_DBL_EPSILON * fabs(yr->val);
	yi->err = 4.0 * GSL_DBL_EPSILON * fabs(yi->val);
	yi_tmp_val = yi->val;
	yi_tmp_err = yi->err;
	gsl_sf_angle_restrict_symm_err_e(yi_tmp_val, yi);
	yi->err += yi_tmp_err;
	return GSL_SUCCESS;
}
/* Lanczos method for real x > 0;
 * gamma=7, truncated at 1/(z+8)
 * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
 */
static
int	lngamma_lanczos(double x, gsl_sf_result * result){
	int k;
	double Ag;
	double term1, term2;

	x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

	Ag = lanczos_7_c[0];
	for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

	/* (x+0.5)*log(x+7.5) - (x+7.5) + 0.5f * log(2.0 * M_PI) + log(Ag(x)) */
	term1 = (x+0.5)*log((x+7.5)/M_E);
	term2 = 0.5f * log(2.0 * M_PI) + log(Ag);
	result->val  = term1 + (term2 - 7.0);
	result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2) + 7.0);
	result->err += GSL_DBL_EPSILON * fabs(result->val);

	return GSL_SUCCESS;
}
/* x = eps near zero
 * gives double-precision for |eps| < 0.02
 */
static
int	lngamma_sgn_0(double eps, gsl_sf_result * lng, double * sgn){
	/* calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps/2 */
	const double c1  = -0.07721566490153286061;
	const double c2  = -0.01094400467202744461;
	const double c3  =  0.09252092391911371098;
	const double c4  = -0.01827191316559981266;
	const double c5  =  0.01800493109685479790;
	const double c6  = -0.00685088537872380685;
	const double c7  =  0.00399823955756846603;
	const double c8  = -0.00189430621687107802;
	const double c9  =  0.00097473237804513221;
	const double c10 = -0.00048434392722255893;
	const double g6  = c6+eps*(c7+eps*(c8 + eps*(c9 + eps*c10)));
	const double g   = eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*g6)))));

	/* calculate Gamma(eps) eps, a positive quantity */
	const double gee = g + 1.0/(1.0+eps) + 0.5*eps;

	lng->val = log(gee/fabs(eps));
	lng->err = 4.0 * GSL_DBL_EPSILON * fabs(lng->val);
	*sgn = GSL_SIGN(eps);

	return GSL_SUCCESS;
}
/* This gets bad near the negative half axis. However, this
 * region can be avoided by use of the reflection formula, as usual.
 * Only the first two terms of the series are kept.
 */
#if 0
static
int	lngamma_complex_stirling(const double zr, const double zi, double * lg_r, double * arg){
	double re_zinv,  im_zinv;
	double re_zinv2, im_zinv2;
	double re_zinv3, im_zinv3;
	double re_zhlnz, im_zhlnz;
	double r, lnr, theta;
	gsl_sf_complex_log_e(zr, zi, &lnr, &theta);  /* z = r e^{i theta} */
	r = exp(lnr);
	re_zinv =  (zr/r)/r;
	im_zinv = -(zi/r)/r;
	re_zinv2 = re_zinv*re_zinv - im_zinv*im_zinv;
	re_zinv2 = 2.0*re_zinv*im_zinv;
	re_zinv3 = re_zinv2*re_zinv - im_zinv2*im_zinv;
	re_zinv3 = re_zinv2*im_zinv + im_zinv2*re_zinv;
	re_zhlnz = (zr - 0.5)*lnr - zi*theta;
	im_zhlnz = zi*lnr + zr*theta;
	*lg_r = re_zhlnz - zr + 0.5*log(2.0 * M_PI) + re_zinv/12.0 - re_zinv3/360.0;
	*arg  = im_zhlnz - zi + 1.0/12.0*im_zinv - im_zinv3/360.0;
	return GSL_SUCCESS;
}
#endif /* 0 */
inline static int lngamma_2_pade(const double eps, gsl_sf_result * result) {
	/* Use (2,2) Pade for Log[Gamma[2+eps]]/eps
	 * plus a correction series.
	 */
	const double n1 = 1.000895834786669227164446568;
	const double n2 = 4.209376735287755081642901277;
	const double d1 = 2.618851904903217274682578255;
	const double d2 = 10.85766559900983515322922936;
	const double num = (eps + n1) * (eps + n2);
	const double den = (eps + d1) * (eps + d2);
	const double pade = 2.85337998765781918463568869 * num/den;
	const double c0 =  0.0001139406357036744;
	const double c1 = -0.0001365435269792533;
	const double c2 =  0.0001067287169183665;
	const double c3 = -0.0000693271800931282;
	const double c4 =  0.0000407220927867950;
	const double eps5 = eps*eps*eps*eps*eps;
	const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
	result->val = eps * (pade + corr);
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
}
/* Chebyshev fit for f(y) = Re(Psi(1+Iy)) + M_EULER - y^2/(1+y^2) - y^2/(2(4+y^2))
 * 1 < y < 10
 *   ==>
 * y(x) = (9x + 11)/2,  -1 < x < 1
 * x(y) = (2y - 11)/9
 *
 * g(x) := f(y(x))
 */
static double r1py_data[] = {
1.59888328244976954803168395603,
0.67905625353213463845115658455,
-0.068485802980122530009506482524,
-0.005788184183095866792008831182,
0.008511258167108615980419855648,
-0.004042656134699693434334556409,
0.001352328406159402601778462956,
-0.000311646563930660566674525382,
0.000018507563785249135437219139,
0.000028348705427529850296492146,
-0.000019487536014574535567541960,
8.0709788710834469408621587335e-06,
-2.2983564321340518037060346561e-06,
3.0506629599604749843855962658e-07,
1.3042238632418364610774284846e-07,
-1.2308657181048950589464690208e-07,
5.7710855710682427240667414345e-08,
-1.8275559342450963966092636354e-08,
3.1020471300626589420759518930e-09,
6.8989327480593812470039430640e-10,
-8.7182290258923059852334818997e-10,
4.4069147710243611798213548777e-10,
-1.4727311099198535963467200277e-10,
2.7589682523262644748825844248e-11,
4.1871826756975856411554363568e-12,
-6.5673460487260087541400767340e-12,
3.4487900886723214020103638000e-12,
-1.1807251417448690607973794078e-12,
2.3798314343969589258709315574e-13,
2.1663630410818831824259465821e-15
};
static cheb_series r1py_cs = {
r1py_data,
29,
-1,1,
18
};
/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.         to  1.00000D+00
 with weighted error   2.03E-17
 log weighted error  16.69
 significant figures required  16.39
 decimal places required  17.37

 Series for APSI       on the interval  0.         to  2.50000D-01
 with weighted error   5.54E-17
 log weighted error  16.26
 significant figures required  14.42
 decimal places required  16.86

 */
static double psics_data[23] = {
-.038057080835217922,
.491415393029387130,
-.056815747821244730,
.008357821225914313,
-.001333232857994342,
.000220313287069308,
-.000037040238178456,
.000006283793654854,
-.000001071263908506,
.000000183128394654,
-.000000031353509361,
.000000005372808776,
-.000000000921168141,
.000000000157981265,
-.000000000027098646,
.000000000004648722,
-.000000000000797527,
.000000000000136827,
-.000000000000023475,
.000000000000004027,
-.000000000000000691,
.000000000000000118,
-.000000000000000020
};
static double apsics_data[16] = {
-.0204749044678185,
-.0101801271534859,
.0000559718725387,
-.0000012917176570,
.0000000572858606,
-.0000000038213539,
.0000000003397434,
-.0000000000374838,
.0000000000048990,
-.0000000000007344,
.0000000000001233,
-.0000000000000228,
.0000000000000045,
-.0000000000000009,
.0000000000000002,
-.0000000000000000
};
static cheb_series psi_cs = {
psics_data,
22,
-1, 1,
17
};
static cheb_series apsi_cs = {
apsics_data,
15,
-1, 1,
9
};
#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
0.0,  /* Infinity */              /* psi(0) */
-M_EULER,                          /* psi(1) */
0.42278433509846713939348790992,  /* ...    */
0.92278433509846713939348790992,
1.25611766843180047272682124325,
1.50611766843180047272682124325,
1.70611766843180047272682124325,
1.87278433509846713939348790992,
2.01564147795560999653634505277,
2.14064147795560999653634505277,
2.25175258906672110764745616389,
2.35175258906672110764745616389,
2.44266167997581201673836525479,
2.52599501330914535007169858813,
2.60291809023222227314862166505,
2.67434666166079370172005023648,
2.74101332832746036838671690315,
2.80351332832746036838671690315,
2.86233685773922507426906984432,
2.91789241329478062982462539988,
2.97052399224214905087725697883,
3.02052399224214905087725697883,
3.06814303986119666992487602645,
3.11359758531574212447033057190,
3.15707584618530734186163491973,
3.1987425128519740085283015864,
3.2387425128519740085283015864,
3.2772040513135124700667631249,
3.3142410883505495071038001619,
3.3499553740648352213895144476,
3.3844381326855248765619282407,
3.4177714660188582098952615740,
3.4500295305349872421533260902,
3.4812795305349872421533260902,
3.5115825608380175451836291205,
3.5409943255438998981248055911,
3.5695657541153284695533770196,
3.5973435318931062473311547974,
3.6243705589201332743581818244,
3.6506863483938174848844976139,
3.6763273740348431259101386396,
3.7013273740348431259101386396,
3.7257176179372821503003825420,
3.7495271417468059598241920658,
3.7727829557002943319172153216,
3.7955102284275670591899425943,
3.8177324506497892814121648166,
3.8394715810845718901078169905,
3.8607481768292527411716467777,
3.8815815101625860745049801110,
3.9019896734278921969539597029,
3.9219896734278921969539597029,
3.9415975165651470989147440166,
3.9608282857959163296839747858,
3.9796962103242182164764276160,
3.9982147288427367349949461345,
4.0163965470245549168131279527,
4.0342536898816977739559850956,
4.0517975495308205809735289552,
4.0690389288411654085597358518,
4.0859880813835382899156680552,
4.1026547480502049565823347218,
4.1190481906731557762544658694,
4.1351772229312202923834981274,
4.1510502388042361653993711433,
4.1666752388042361653993711433,
4.1820598541888515500147557587,
4.1972113693403667015299072739,
4.2121367424746950597388624977,
4.2268426248276362362094507330,
4.2413353784508246420065521823,
4.2556210927365389277208378966,
4.2697055997787924488475984600,
4.2835944886676813377364873489,
4.2972931188046676391063503626,
4.3108066323181811526198638761,
4.3241399656515144859531972094,
4.3372978603883565912163551041,
4.3502848733753695782293421171,
4.3631053861958823987421626300,
4.3757636140439836645649474401,
4.3882636140439836645649474401,
4.4006092930563293435772931191,
4.4128044150075488557724150703,
4.4248526077786331931218126607,
4.4367573696833950978837174226,
4.4485220755657480390601880108,
4.4601499825424922251066996387,
4.4716442354160554434975042364,
4.4830078717796918071338678728,
4.4942438268358715824147667492,
4.5053549379469826935258778603,
4.5163439489359936825368668713,
4.5272135141533849868846929582,
4.5379662023254279976373811303,
4.5486045001977684231692960239,
4.5591308159872421073798223397,
4.5695474826539087740464890064,
4.5798567610044242379640147796,
4.5900608426370772991885045755,
4.6001618527380874001986055856
};
#define PSI_1_TABLE_NMAX 100
static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
0.0,  /* Infinity */              /* psi(1,0) */
M_PI*M_PI/6.0,                    /* psi(1,1) */
0.644934066848226436472415,       /* ...      */
0.394934066848226436472415,
0.2838229557371153253613041,
0.2213229557371153253613041,
0.1813229557371153253613041,
0.1535451779593375475835263,
0.1331370146940314251345467,
0.1175120146940314251345467,
0.1051663356816857461222010,
0.0951663356816857461222010,
0.0869018728717683907503002,
0.0799574284273239463058557,
0.0740402686640103368384001,
0.0689382278476838062261552,
0.0644937834032393617817108,
0.0605875334032393617817108,
0.0571273257907826143768665,
0.0540409060376961946237801,
0.0512708229352031198315363,
0.0487708229352031198315363,
0.0465032492390579951149830,
0.0444371335365786562720078,
0.0425467743683366902984728,
0.0408106632572255791873617,
0.0392106632572255791873617,
0.0377313733163971768204978,
0.0363596312039143235969038,
0.0350841209998326909438426,
0.0338950603577399442137594,
0.0327839492466288331026483,
0.0317433665203020901265817,
0.03076680402030209012658168,
0.02984853037475571730748159,
0.02898347847164153045627052,
0.02816715194102928555831133,
0.02739554700275768062003973,
0.02666508681283803124093089,
0.02597256603721476254286995,
0.02531510384129102815759710,
0.02469010384129102815759710,
0.02409521984367056414807896,
0.02352832641963428296894063,
0.02298749353699501850166102,
0.02247096461137518379091722,
0.02197713745088135663042339,
0.02150454765882086513703965,
0.02105185413233829383780923,
0.02061782635456051606003145,
0.02020133322669712580597065,
0.01980133322669712580597065,
0.01941686571420193164987683,
0.01904704322899483105816086,
0.01869104465298913508094477,
0.01834810912486842177504628,
0.01801753061247172756017024,
0.01769865306145131939690494,
0.01739086605006319997554452,
0.01709360088954001329302371,
0.01680632711763538818529605,
0.01652854933985761040751827,
0.01625980437882562975715546,
0.01599965869724394401313881,
0.01574770606433893015574400,
0.01550356543933893015574400,
0.01526687904880638577704578,
0.01503731063741979257227076,
0.01481454387422086185273411,
0.01459828089844231513993134,
0.01438824099085987447620523,
0.01418415935820681325171544,
0.01398578601958352422176106,
0.01379288478501562298719316,
0.01360523231738567365335942,
0.01342261726990576130858221,
0.01324483949212798353080444,
0.01307170929822216635628920,
0.01290304679189732236910755,
0.01273868124291638877278934,
0.01257845051066194236996928,
0.01242220051066194236996928,
0.01226978472038606978956995,
0.01212106372098095378719041,
0.01197590477193174490346273,
0.01183418141592267460867815,
0.01169577311142440471248438,
0.01156056489076458859566448,
0.01142844704164317229232189,
0.01129931481023821361463594,
0.01117306812421372175754719,
0.01104961133409026496742374,
0.01092885297157366069257770,
0.01081070552355853781923177,
0.01069508522063334415522437,
0.01058191183901270133041676,
0.01047110851491297833872701,
0.01036260157046853389428257,
0.01025632035036012704977199,  /* ...        */
0.01015219706839427948625679,  /* psi(1,99)  */
0.01005016666333357139524567   /* psi(1,100) */
};
/* digamma for x both positive and negative; we do both
 * cases here because of the way we use even/odd parts
 * of the function
 */
static int psi_x(const double x, gsl_sf_result * result) {
	const double y = fabs(x);

	if(x == 0.0 || x == -1.0 || x == -2.0) {
		DOMAIN_ERROR(result);
	}
	else if(y >= 2.0) {
		const double t = 8.0/(y*y)-1.0;
		gsl_sf_result result_c;
		cheb_eval_e(&apsi_cs, t, &result_c);
		if(x < 0.0) {
			const double s = sin(M_PI*x);
			const double c = cos(M_PI*x);
			if(fabs(s) < 2.0*GSL_SQRT_DBL_MIN) {
				DOMAIN_ERROR(result);
			}
			else {
				result->val  = log(y) - 0.5/x + result_c.val - M_PI * c/s;
				result->err  = M_PI*fabs(x)*GSL_DBL_EPSILON/(s*s);
				result->err += result_c.err;
				result->err += GSL_DBL_EPSILON * fabs(result->val);
				return GSL_SUCCESS;
			}
		} else {
			result->val  = log(y) - 0.5/x + result_c.val;
			result->err  = result_c.err;
			result->err += GSL_DBL_EPSILON * fabs(result->val);
			return GSL_SUCCESS;
		}
	}
	else { /* -2 < x < 2 */
		gsl_sf_result result_c;
		if(x < -1.0) { /* x = -2 + v */
			const double v  = x + 2.0;
			const double t1 = 1.0/x;
			const double t2 = 1.0/(x+1.0);
			const double t3 = 1.0/v;
			cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);

			result->val  = -(t1 + t2 + t3) + result_c.val;
			result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)) + fabs(x/(t3*t3)));
			result->err += result_c.err;
			result->err += GSL_DBL_EPSILON * fabs(result->val);
			return GSL_SUCCESS;
		}
		else if(x < 0.0) { /* x = -1 + v */
			const double v  = x + 1.0;
			const double t1 = 1.0/x;
			const double t2 = 1.0/v;
			cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);

			result->val  = -(t1 + t2) + result_c.val;
			result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)));
			result->err += result_c.err;
			result->err += GSL_DBL_EPSILON * fabs(result->val);
			return GSL_SUCCESS;
		}
		else if(x < 1.0) { /* x = v */
			const double t1 = 1.0/x;
			cheb_eval_e(&psi_cs, 2.0*x-1.0, &result_c);

			result->val  = -t1 + result_c.val;
			result->err  = GSL_DBL_EPSILON * t1;
			result->err += result_c.err;
			result->err += GSL_DBL_EPSILON * fabs(result->val);
			return GSL_SUCCESS;
		}
		else { /* x = 1 + v */
			const double v = x - 1.0;
			return cheb_eval_e(&psi_cs, 2.0*v-1.0, result);
		}
	}
}
/* psi(z) for large |z| in the right half-plane; [Abramowitz + Stegun, 6.3.18]
 static
 gsl_complex
 psi_complex_asymp(gsl_complex z)
 {
 // coefficients in the asymptotic expansion for large z;   * let w = z^(-2) and write the expression in the form   *   ln(z) - 1/(2z) - 1/12 w (1 + c1 w + c2 w + c3 w + ... )
 static const double c1 = -0.1;
 static const double c2 =  1.0/21.0;
 static const double c3 = -0.05;

 gsl_complex zi = gsl_complex_inverse(z);
 gsl_complex w  = gsl_complex_mul(zi, zi);
 gsl_complex cs;

 // Horner method evaluation of term in parentheses
 gsl_complex sum;
 sum = gsl_complex_mul_real(w, c3/c2);
 sum = gsl_complex_add_real(sum, 1.0);
 sum = gsl_complex_mul_real(sum, c2/c1);
 sum = gsl_complex_mul(sum, w);
 sum = gsl_complex_add_real(sum, 1.0);
 sum = gsl_complex_mul_real(sum, c1);
 sum = gsl_complex_mul(sum, w);
 sum = gsl_complex_add_real(sum, 1.0);

 // correction added to log(z)
 cs = gsl_complex_mul(sum, w);
 cs = gsl_complex_mul_real(cs, -1.0/12.0);
 cs = gsl_complex_add(cs, gsl_complex_mul_real(zi, -0.5));

 return gsl_complex_add(gsl_complex_log(z), cs);
 }*/
/* psi(z) for complex z in the right half-plane
 static int	 psi_complex_rhp(
 gsl_complex z,
 gsl_sf_result * result_re,
 gsl_sf_result * result_im
 )
 {
 int n_recurse = 0;
 int i;
 gsl_complex a;

 if(GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0)
 {
 result_re->val = 0.0;
 result_im->val = 0.0;
 result_re->err = 0.0;
 result_im->err = 0.0;
 return GSL_EDOM;
 }

 // compute the number of recurrences to apply
 if(GSL_REAL(z) < 20.0 && fabs(GSL_IMAG(z)) < 20.0)
 {
 const double sp = sqrt(20.0 + GSL_IMAG(z));
 const double sn = sqrt(20.0 - GSL_IMAG(z));
 const double rhs = sp*sn - GSL_REAL(z);
 if(rhs > 0.0) n_recurse = ceil(rhs);
 }

 // compute asymptotic at the large value z + n_recurse
 a = psi_complex_asymp(gsl_complex_add_real(z, n_recurse));

 result_re->err = 2.0 * GSL_DBL_EPSILON * fabs(GSL_REAL(a));
 result_im->err = 2.0 * GSL_DBL_EPSILON * fabs(GSL_IMAG(a));

 // descend recursively, if necessary
 for(i = n_recurse; i >= 1; --i)
 {
 gsl_complex zn = gsl_complex_add_real(z, i - 1.0);
 gsl_complex zn_inverse = gsl_complex_inverse(zn);
 a = gsl_complex_sub(a, zn_inverse);

 // accumulate the error, to catch cancellations
 result_re->err += 2.0 * GSL_DBL_EPSILON * fabs(GSL_REAL(zn_inverse));
 result_im->err += 2.0 * GSL_DBL_EPSILON * fabs(GSL_IMAG(zn_inverse));
 }

 result_re->val = GSL_REAL(a);
 result_im->val = GSL_IMAG(a);

 result_re->err += 2.0 * GSL_DBL_EPSILON * fabs(result_re->val);
 result_im->err += 2.0 * GSL_DBL_EPSILON * fabs(result_im->val);

 return GSL_SUCCESS;
 }*/
/* generic polygamma; assumes n >= 0 and x > 0
 */
int gsl_sf_psi_e(const double x, gsl_sf_result * result) {
	return psi_x(x, result);
}
/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
1.00000000000000000000000000000,
0.083333333333333333333333333333,
-0.00138888888888888888888888888889,
0.000033068783068783068783068783069,
-8.2671957671957671957671957672e-07,
2.0876756987868098979210090321e-08,
-5.2841901386874931848476822022e-10,
1.3382536530684678832826980975e-11,
-3.3896802963225828668301953912e-13,
8.5860620562778445641359054504e-15,
-2.1748686985580618730415164239e-16,
5.5090028283602295152026526089e-18,
-1.3954464685812523340707686264e-19,
3.5347070396294674716932299778e-21,
-8.9535174270375468504026113181e-23
};
int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result){

	if(s <= 1.0 || q <= 0.0) {
		DOMAIN_ERROR(result);
	}
	else {
		const double max_bits = 54.0;
		const double ln_term0 = -s * log(q);

		if(ln_term0 < GSL_LOG_DBL_MIN + 1.0) {
			UNDERFLOW_ERROR(result);
		}
		else if(ln_term0 > GSL_LOG_DBL_MAX - 1.0) {
			OVERFLOW_ERROR (result);
		}
		else if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25)) {
			result->val = pow(q, -s);
			result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
			return GSL_SUCCESS;
		}
		else if(s > 0.5*max_bits && q < 1.0) {
			const double p1 = pow(q, -s);
			const double p2 = pow(q/(1.0+q), s);
			const double p3 = pow(q/(2.0+q), s);
			result->val = p1 * (1.0 + p2 + p3);
			result->err = GSL_DBL_EPSILON * (0.5*s + 2.0) * fabs(result->val);
			return GSL_SUCCESS;
		}
		else {
			/* Euler-Maclaurin summation formula
			 * [Moshier, p. 400, with several typo corrections]
			 */
			const int jmax = 12;
			const int kmax = 10;
			int j, k;
			const double pmax  = pow(kmax + q, -s);
			double scp = s;
			double pcp = pmax / (kmax + q);
			double ans = pmax*((kmax+q)/(s-1.0) + 0.5);

			for(k=0; k<kmax; k++) {
				ans += pow(k + q, -s);
			}
			for(j=0; j<=jmax; j++) {
				double delta = hzeta_c[j+1] * scp * pcp;
				ans += delta;
				if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
				scp *= (s+2*j+1)*(s+2*j+2);
				pcp /= (kmax + q)*(kmax + q);
			}

			result->val = ans;
			result->err = 2.0 * (jmax + 1.0) * GSL_DBL_EPSILON * fabs(ans);
			return GSL_SUCCESS;
		}
	}
    return GSL_FAILURE;
}
double gsl_sf_hzeta_e_s_is_2(const double q){
    /* Euler-Maclaurin summation formula
     * [Moshier, p. 400, with several typo corrections]
     */
    const int jmax = 12;
    const int kmax = 10;
    int j, k;
    const double pmax  = pow(kmax + q, -2.0);
    double scp = 2.0;
    double pcp = pmax / (kmax + q);
    double ans = pmax*(kmax+q + 0.5);

    for(k=0; k<kmax; k++) {
        ans += pow(k + q, -2.0);
    }
    for(j=0; j<=jmax; j++) {
        double delta = hzeta_c[j+1] * scp * pcp;
        ans += delta;
        if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
        scp *= (2*j+3)*(2*j+4);
        pcp /= (kmax + q)*(kmax + q);
    }
    return ans;
}
static int psi_n_xg0(const int n, const double x, gsl_sf_result * result){
	if(n == 0) {
		return gsl_sf_psi_e(x, result);
	}
	else {
		/* Abramowitz + Stegun 6.4.10 */
		gsl_sf_result ln_nf;
		gsl_sf_result hzeta;
		int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
		int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
		int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
											hzeta.val, hzeta.err,
											result);
		if((n & 1) ==0) result->val = -result->val;
		return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
	}
}
/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
int gsl_sf_psi_int_e(const int n, gsl_sf_result * result){

	if(n <= 0) {
		DOMAIN_ERROR(result);
	}
	else if(n <= PSI_TABLE_NMAX) {
		result->val = psi_table[n];
		result->err = GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		/* Abramowitz+Stegun 6.3.18 */
		const double c2 = -1.0/12.0;
		const double c3 =  1.0/120.0;
		const double c4 = -1.0/252.0;
		const double c5 =  1.0/240.0;
		const double ni2 = (1.0/n)*(1.0/n);
		const double ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2*c5)));
		result->val  = log(n) - 0.5/n + ser;
		result->err  = GSL_DBL_EPSILON * (fabs(log(n)) + fabs(0.5/n) + fabs(ser));
		result->err += GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
}
int gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result) {
	const double ay = fabs(y);


	if(ay > 1000.0) {
		/* [Abramowitz+Stegun, 6.3.19] */
		const double yi2 = 1.0/(ay*ay);
		const double lny = log(ay);
		const double sum = yi2 * (1.0/12.0 + 1.0/120.0 * yi2 + 1.0/252.0 * yi2*yi2);
		result->val = lny + sum;
		result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lny) + fabs(sum));
		return GSL_SUCCESS;
	}
	else if(ay > 10.0) {
		/* [Abramowitz+Stegun, 6.3.19] */
		const double yi2 = 1.0/(ay*ay);
		const double lny = log(ay);
		const double sum = yi2 * (1.0/12.0 +
								  yi2 * (1.0/120.0 +
										 yi2 * (1.0/252.0 +
												yi2 * (1.0/240.0 +
													   yi2 * (1.0/132.0 + 691.0/32760.0 * yi2)))));
		result->val = lny + sum;
		result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lny) + fabs(sum));
		return GSL_SUCCESS;
	}
	else if(ay > 1.0){
		const double y2 = ay*ay;
		const double x  = (2.0*ay - 11.0)/9.0;
		const double v  = y2*(1.0/(1.0+y2) + 0.5/(4.0+y2));
		gsl_sf_result result_c;
		cheb_eval_e(&r1py_cs, x, &result_c);
		result->val  = result_c.val - M_EULER + v;
		result->err  = result_c.err;
		result->err += 2.0 * GSL_DBL_EPSILON * (fabs(v) + M_EULER + fabs(result_c.val));
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		result->err *= 5.0; /* FIXME: losing a digit somewhere... maybe at x=... ? */
		return GSL_SUCCESS;
	}
	else {
		/* [Abramowitz+Stegun, 6.3.17]
		 *
		 * -M_EULER + y^2 Sum[1/n 1/(n^2 + y^2), {n,1,M}]
		 *   +     Sum[1/n^3, {n,M+1,Infinity}]
		 *   - y^2 Sum[1/n^5, {n,M+1,Infinity}]
		 *   + y^4 Sum[1/n^7, {n,M+1,Infinity}]
		 *   - y^6 Sum[1/n^9, {n,M+1,Infinity}]
		 *   + O(y^8)
		 *
		 * We take M=50 for at least 15 digit precision.
		 */
		const int M = 50;
		const double y2 = y*y;
		const double c0 = 0.00019603999466879846570;
		const double c2 = 3.8426659205114376860e-08;
		const double c4 = 1.0041592839497643554e-11;
		const double c6 = 2.9516743763500191289e-15;
		const double p  = c0 + y2 *(-c2 + y2*(c4 - y2*c6));
		double sum = 0.0;
		double v;

		int n;
		for(n=1; n<=M; n++) {
			sum += 1.0/(n * (n*n + y*y));
		}

		v = y2 * (sum + p);
		result->val  = -M_EULER + v;
		result->err  = GSL_DBL_EPSILON * (M_EULER + fabs(v));
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
}
int gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result) {
	if(n <= 0) {
		DOMAIN_ERROR(result);
	}
	else if(n <= PSI_1_TABLE_NMAX) {
		result->val = psi_1_table[n];
		result->err = GSL_DBL_EPSILON * result->val;
		return GSL_SUCCESS;
	}
	else {
		/* Abramowitz+Stegun 6.4.12
		 * double-precision for n > 100
		 */
		const double c0 = -1.0/30.0;
		const double c1 =  1.0/42.0;
		const double c2 = -1.0/30.0;
		const double ni2 = (1.0/n)*(1.0/n);
		const double ser =  ni2*ni2 * (c0 + ni2*(c1 + c2*ni2));
		result->val = (1.0 + 0.5/n + 1.0/(6.0*n*n) + ser) / n;
		result->err = GSL_DBL_EPSILON * result->val;
		return GSL_SUCCESS;
	}
}
int gsl_sf_psi_1_e(const double x, gsl_sf_result * result) {

	if(x == 0.0 || x == -1.0 || x == -2.0) {
		DOMAIN_ERROR(result);
	}
	else if(x > 0.0)
	{
		return psi_n_xg0(1, x, result);
	}
	else if(x > -5.0)
	{
		/* Abramowitz + Stegun 6.4.6 */
		int M = -(int)floor(x);
		double fx = x + M;
		double sum = 0.0;
		int m;

		if(fx == 0.0)
			DOMAIN_ERROR(result);

		for(m = 0; m < M; ++m)
			sum += 1.0/((x+m)*(x+m));

		{
			int stat_psi = psi_n_xg0(1, fx, result);
			result->val += sum;
			result->err += M * GSL_DBL_EPSILON * sum;
			return stat_psi;
		}
	}
	else
	{
		/* Abramowitz + Stegun 6.4.7 */
		const double sin_px = sin(M_PI * x);
		const double d = M_PI*M_PI/(sin_px*sin_px);
		gsl_sf_result r;
		int stat_psi = psi_n_xg0(1, 1.0-x, &r);
		result->val = d - r.val;
		result->err = r.err + 2.0*GSL_DBL_EPSILON*d;
		return stat_psi;
	}
}
int gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result) {

	if(n == 0)
	{
		return gsl_sf_psi_e(x, result);
	}
	else if(n == 1)
	{
		return gsl_sf_psi_1_e(x, result);
	}
	else if(n < 0 || x <= 0.0) {
		DOMAIN_ERROR(result);
	}
	else {
		gsl_sf_result ln_nf;
		gsl_sf_result hzeta;
		int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
		int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
		int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
											hzeta.val, hzeta.err,
											result);
		if((n & 1) ==0) result->val = -result->val;
		return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
	}
}
/* x near a negative integer
 * Calculates sign as well as log(|gamma(x)|).
 * x = -N + eps
 * assumes N >= 1
 */
static
int	lngamma_sgn_sing(int N, double eps, gsl_sf_result * lng, double * sgn){
	if(eps == 0.0) {
		lng->val = 0.0;
		lng->err = 0.0;
		*sgn = 0.0;
		GSL_ERROR ("error", GSL_EDOM);
	}
	else if(N == 1) {
		/* calculate series for
		 * g = eps gamma(-1+eps) + 1 + eps/2 (1+3eps)/(1-eps^2)
		 * double-precision for |eps| < 0.02
		 */
		const double c0 =  0.07721566490153286061;
		const double c1 =  0.08815966957356030521;
		const double c2 = -0.00436125434555340577;
		const double c3 =  0.01391065882004640689;
		const double c4 = -0.00409427227680839100;
		const double c5 =  0.00275661310191541584;
		const double c6 = -0.00124162645565305019;
		const double c7 =  0.00065267976121802783;
		const double c8 = -0.00032205261682710437;
		const double c9 =  0.00016229131039545456;
		const double g5 = c5 + eps*(c6 + eps*(c7 + eps*(c8 + eps*c9)));
		const double g  = eps*(c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*g5)))));

		/* calculate eps gamma(-1+eps), a negative quantity */
		const double gam_e = g - 1.0 - 0.5*eps*(1.0+3.0*eps)/(1.0 - eps*eps);

		lng->val = log(fabs(gam_e)/fabs(eps));
		lng->err = 2.0 * GSL_DBL_EPSILON * fabs(lng->val);
		*sgn = ( eps > 0.0 ? -1.0 : 1.0 );
		return GSL_SUCCESS;
	}
	else {
		double g;

		/* series for sin(Pi(N+1-eps))/(Pi eps) modulo the sign
		 * double-precision for |eps| < 0.02
		 */
		const double cs1 = -1.6449340668482264365;
		const double cs2 =  0.8117424252833536436;
		const double cs3 = -0.1907518241220842137;
		const double cs4 =  0.0261478478176548005;
		const double cs5 = -0.0023460810354558236;
		const double e2  = eps*eps;
		const double sin_ser = 1.0 + e2*(cs1+e2*(cs2+e2*(cs3+e2*(cs4+e2*cs5))));

		/* calculate series for ln(gamma(1+N-eps))
		 * double-precision for |eps| < 0.02
		 */
		double aeps = fabs(eps);
		double c1, c2, c3, c4, c5, c6, c7;
		double lng_ser;
		gsl_sf_result c0;
		gsl_sf_result psi_0;
		gsl_sf_result psi_1;
		gsl_sf_result psi_2;
		gsl_sf_result psi_3;
		gsl_sf_result psi_4;
		gsl_sf_result psi_5;
		gsl_sf_result psi_6;
		psi_2.val = 0.0;
		psi_3.val = 0.0;
		psi_4.val = 0.0;
		psi_5.val = 0.0;
		psi_6.val = 0.0;
		gsl_sf_lnfact_e(N, &c0);
		gsl_sf_psi_int_e(N+1, &psi_0);
		gsl_sf_psi_1_int_e(N+1, &psi_1);
		if(aeps > 0.00001) gsl_sf_psi_n_e(2, N+1.0, &psi_2);
		if(aeps > 0.0002)  gsl_sf_psi_n_e(3, N+1.0, &psi_3);
		if(aeps > 0.001)   gsl_sf_psi_n_e(4, N+1.0, &psi_4);
		if(aeps > 0.005)   gsl_sf_psi_n_e(5, N+1.0, &psi_5);
		if(aeps > 0.01)    gsl_sf_psi_n_e(6, N+1.0, &psi_6);
		c1 = psi_0.val;
		c2 = psi_1.val/2.0;
		c3 = psi_2.val/6.0;
		c4 = psi_3.val/24.0;
		c5 = psi_4.val/120.0;
		c6 = psi_5.val/720.0;
		c7 = psi_6.val/5040.0;
		lng_ser = c0.val-eps*(c1-eps*(c2-eps*(c3-eps*(c4-eps*(c5-eps*(c6-eps*c7))))));

		/* calculate
		 * g = ln(|eps gamma(-N+eps)|)
		 *   = -ln(gamma(1+N-eps)) + ln(|eps Pi/sin(Pi(N+1+eps))|)
		 */
		g = -lng_ser - log(sin_ser);

		lng->val = g - log(fabs(eps));
		lng->err = c0.err + 2.0 * GSL_DBL_EPSILON * (fabs(g) + fabs(lng->val));

		*sgn = ( (N & 1) ? -1.0 : 1.0 ) * ( eps > 0.0 ? 1.0 : -1.0 );

		return GSL_SUCCESS;
	}
}
inline
static
int	lngamma_1_pade(const double eps, gsl_sf_result * result){
	/* Use (2,2) Pade for Log[Gamma[1+eps]]/eps
	 * plus a correction series.
	 */
	const double n1 = -1.0017419282349508699871138440;
	const double n2 =  1.7364839209922879823280541733;
	const double d1 =  1.2433006018858751556055436011;
	const double d2 =  5.0456274100274010152489597514;
	const double num = (eps + n1) * (eps + n2);
	const double den = (eps + d1) * (eps + d2);
	const double pade = 2.0816265188662692474880210318 * num / den;
	const double c0 =  0.004785324257581753;
	const double c1 = -0.01192457083645441;
	const double c2 =  0.01931961413960498;
	const double c3 = -0.02594027398725020;
	const double c4 =  0.03141928755021455;
	const double eps5 = eps*eps*eps*eps*eps;
	const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
	result->val = eps * (pade + corr);
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
}
int gsl_sf_lngamma_e(double x, gsl_sf_result * result) {

	if(fabs(x - 1.0) < 0.01) {
		/* Note that we must amplify the errors
		 * from the Pade evaluations because of
		 * the way we must pass the argument, i.e.
		 * writing (1-x) is a loss of precision
		 * when x is near 1.
		 */
		int stat = lngamma_1_pade(x - 1.0, result);
		result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
		return stat;
	}
	else if(fabs(x - 2.0) < 0.01) {
		int stat = lngamma_2_pade(x - 2.0, result);
		result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
		return stat;
	}
	else if(x >= 0.5) {
		return lngamma_lanczos(x, result);
	}
	else if(x == 0.0) {
		DOMAIN_ERROR(result);
	}
	else if(fabs(x) < 0.02) {
		double sgn;
		return lngamma_sgn_0(x, result, &sgn);
	}
	else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
		/* Try to extract a fractional
		 * part from x.
		 */
		double z  = 1.0 - x;
		double s  = sin(M_PI*z);
		double as = fabs(s);
		if(s == 0.0) {
			DOMAIN_ERROR(result);
		}
		else if(as < M_PI*0.015) {
			/* x is near a negative integer, -N */
			if(x < ExCo<int>::mkMinimum() + 2.0) {
				result->val = 0.0;
				result->err = 0.0;
				GSL_ERROR ("error", GSL_EROUND);
			}
			else {
				int N = -(int)(x - 0.5);
				double eps = x + N;
				double sgn;
				return lngamma_sgn_sing(N, eps, result, &sgn);
			}
		}
		else {
			gsl_sf_result lg_z;
			lngamma_lanczos(z, &lg_z);
			result->val = log(M_PI) - (log(as) + lg_z.val);
			result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg_z.err;
			return GSL_SUCCESS;
		}
	}
	else {
		/* |x| was too large to extract any fractional part */
		result->val = 0.0;
		result->err = 0.0;
		GSL_ERROR ("error", GSL_EROUND);
	}
}
int	gsl_sf_exprel_n_e(const int N, const double x, gsl_sf_result * result){
	if(N < 0) {
		DOMAIN_ERROR(result);
	}
	else if(x == 0.0) {
		result->val = 1.0;
		result->err = 0.0;
		return GSL_SUCCESS;
	}
	else if(fabs(x) < GSL_ROOT3_DBL_EPSILON * N) {
		result->val = 1.0 + x/(N+1) * (1.0 + x/(N+2));
		result->err = 2.0 * GSL_DBL_EPSILON;
		return GSL_SUCCESS;
	}
	else if(N == 0) {
		return gsl_sf_exp_e(x, result);
	}
	else if(N == 1) {
		return gsl_sf_exprel_e(x, result);
	}
	else if(N == 2) {
		return gsl_sf_exprel_2_e(x, result);
	}
	else {
		if(x > N && (-x + N*(1.0 + log(x/N)) < GSL_LOG_DBL_EPSILON)) {
			/* x is much larger than n.
			 * Ignore polynomial part, so
			 * exprel_N(x) ~= e^x N!/x^N
			 */
			gsl_sf_result lnf_N;
			double lnr_val;
			double lnr_err;
			double lnterm;
			gsl_sf_lnfact_e(N, &lnf_N);
			lnterm = N*log(x);
			lnr_val  = x + lnf_N.val - lnterm;
			lnr_err  = GSL_DBL_EPSILON * (fabs(x) + fabs(lnf_N.val) + fabs(lnterm));
			lnr_err += lnf_N.err;
			return gsl_sf_exp_err_e(lnr_val, lnr_err, result);
		}
		else if(x > N) {
			/* Write the identity
			 *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
			 * then use the asymptotic expansion
			 * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
			 */
			double ln_x = log(x);
			gsl_sf_result lnf_N;
			double lg_N;
			double lnpre_val;
			double lnpre_err;
			gsl_sf_lnfact_e(N, &lnf_N);    /* log(N!)       */
			lg_N  = lnf_N.val - log(N);       /* log(Gamma(N)) */
			lnpre_val  = x + lnf_N.val - N*ln_x;
			lnpre_err  = GSL_DBL_EPSILON * (fabs(x) + fabs(lnf_N.val) + fabs(N*ln_x));
			lnpre_err += lnf_N.err;
			if(lnpre_val < GSL_LOG_DBL_MAX - 5.0) {
				int stat_eG;
				gsl_sf_result bigG_ratio;
				gsl_sf_result pre;
				int stat_ex = gsl_sf_exp_err_e(lnpre_val, lnpre_err, &pre);
				double ln_bigG_ratio_pre = -x + (N-1)*ln_x - lg_N;
				double bigGsum = 1.0;
				double term = 1.0;
				int k;
				for(k=1; k<N; k++) {
					term *= (N-k)/x;
					bigGsum += term;
				}
				stat_eG = gsl_sf_exp_mult_e(ln_bigG_ratio_pre, bigGsum, &bigG_ratio);
				if(stat_eG == GSL_SUCCESS) {
					result->val  = pre.val * (1.0 - bigG_ratio.val);
					result->err  = pre.val * (2.0*GSL_DBL_EPSILON + bigG_ratio.err);
					result->err += pre.err * fabs(1.0 - bigG_ratio.val);
					result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
					return stat_ex;
				}
				else {
					result->val = 0.0;
					result->err = 0.0;
					return stat_eG;
				}
			}
			else {
				OVERFLOW_ERROR(result);
			}
		}
		else if(x > -10.0*N) {
			return exprel_n_CF(N, x, result);
		}
		else {
			/* x -> -Inf asymptotic:
			 * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
			 *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
			 */
			double sum  = 1.0;
			double term = 1.0;
			int k;
			for(k=1; k<N; k++) {
				term *= (N-k)/x;
				sum  += term;
			}
			result->val = -N/x * sum;
			result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
			return GSL_SUCCESS;
		}
	}
}
/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/
/*
 int	gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  ){
  if(x >= 0.0)
  {
    gsl_complex z = gsl_complex_rect(x, y);
    return psi_complex_rhp(z, result_re, result_im);
  }
  else
  {
    // reflection formula [Abramowitz+Stegun, 6.3.7]
    gsl_complex z = gsl_complex_rect(x, y);
    gsl_complex omz = gsl_complex_rect(1.0 - x, -y);
    gsl_complex zpi = gsl_complex_mul_real(z, M_PI);
    gsl_complex cotzpi = gsl_complex_cot(zpi);
    int ret_val = psi_complex_rhp(omz, result_re, result_im);
    if(GSL_IS_REAL(GSL_REAL(cotzpi)) && GSL_IS_REAL(GSL_IMAG(cotzpi)))
    {
      result_re->val -= M_PI * GSL_REAL(cotzpi);
      result_im->val -= M_PI * GSL_IMAG(cotzpi);
      return ret_val;
    }
    else
    {
      GSL_ERROR("singularity", GSL_EDOM);
    }
  }
}*/
/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
/*
double gsl_sf_psi_int(const int n){
  EVAL_RESULT(gsl_sf_psi_int_e(n, &result));
}
double gsl_sf_psi(const double x){
  EVAL_RESULT(gsl_sf_psi_e(x, &result));
}
double gsl_sf_psi_1piy(const double x){
  EVAL_RESULT(gsl_sf_psi_1piy_e(x, &result));
}
double gsl_sf_psi_1_int(const int n){
  EVAL_RESULT(gsl_sf_psi_1_int_e(n, &result));
}
double gsl_sf_psi_1(const double x){
  EVAL_RESULT(gsl_sf_psi_1_e(x, &result));
}
double gsl_sf_psi_n(const int n, const double x){
  EVAL_RESULT(gsl_sf_psi_n_e(n, x, &result));
}*/
int	gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result){
	const double adx = fabs(dx);


	if(x + adx > ExCo<int>::mkMaximum() - 1) {
		OVERFLOW_ERROR_E10(result);
	}
	else if(x - adx < ExCo<int>::mkMinimum() + 1) {
		UNDERFLOW_ERROR_E10(result);
	}
	else {
		const int    N  = (int)floor(x/M_LN10);
		const double ex = exp(x-N*M_LN10);
		result->val = ex;
		result->err = ex * (2.0 * GSL_DBL_EPSILON * (fabs(x) + 1.0) + adx);
		result->e10 = N;
		return GSL_SUCCESS;
	}
return GSL_FAILURE;
}
static double lopx_data[21] = {
2.16647910664395270521272590407,
-0.28565398551049742084877469679,
0.01517767255690553732382488171,
-0.00200215904941415466274422081,
0.00019211375164056698287947962,
-0.00002553258886105542567601400,
2.9004512660400621301999384544e-06,
-3.8873813517057343800270917900e-07,
4.7743678729400456026672697926e-08,
-6.4501969776090319441714445454e-09,
8.2751976628812389601561347296e-10,
-1.1260499376492049411710290413e-10,
1.4844576692270934446023686322e-11,
-2.0328515972462118942821556033e-12,
2.7291231220549214896095654769e-13,
-3.7581977830387938294437434651e-14,
5.1107345870861673561462339876e-15,
-7.0722150011433276578323272272e-16,
9.7089758328248469219003866867e-17,
-1.3492637457521938883731579510e-17,
1.8657327910677296608121390705e-18
};
static cheb_series lopx_cs = {
lopx_data,
20,
-1, 1,
10
};
int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result){

	if(x <= -1.0) {
		DOMAIN_ERROR(result);
	}
	else if(fabs(x) < GSL_ROOT6_DBL_EPSILON) {
		const double c1 = -0.5;
		const double c2 =  1.0/3.0;
		const double c3 = -1.0/4.0;
		const double c4 =  1.0/5.0;
		const double c5 = -1.0/6.0;
		const double c6 =  1.0/7.0;
		const double c7 = -1.0/8.0;
		const double c8 =  1.0/9.0;
		const double c9 = -1.0/10.0;
		const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
		result->val = x * (1.0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*t)))));
		result->err = GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(fabs(x) < 0.5) {
		double t = 0.5*(8.0*x + 1.0)/(x+2.0);
		gsl_sf_result c;
		cheb_eval_e(&lopx_cs, t, &c);
		result->val = x * c.val;
		result->err = fabs(x * c.err);
		return GSL_SUCCESS;
	}
	else {
		result->val = log(1.0 + x);
		result->err = GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
}
static struct {int n; double f; long i; } fact_table[GSL_SF_FACT_NMAX + 1] = {
{ 0,  1.0,     1L     },
{ 1,  1.0,     1L     },
{ 2,  2.0,     2L     },
{ 3,  6.0,     6L     },
{ 4,  24.0,    24L    },
{ 5,  120.0,   120L   },
{ 6,  720.0,   720L   },
{ 7,  5040.0,  5040L  },
{ 8,  40320.0, 40320L },
{ 9,  362880.0,     362880L    },
{ 10, 3628800.0,    3628800L   },
{ 11, 39916800.0,   39916800L  },
{ 12, 479001600.0,  479001600L },
{ 13, 6227020800.0,                               0 },
{ 14, 87178291200.0,                              0 },
{ 15, 1307674368000.0,                            0 },
{ 16, 20922789888000.0,                           0 },
{ 17, 355687428096000.0,                          0 },
{ 18, 6402373705728000.0,                         0 },
{ 19, 121645100408832000.0,                       0 },
{ 20, 2432902008176640000.0,                      0 },
{ 21, 51090942171709440000.0,                     0 },
{ 22, 1124000727777607680000.0,                   0 },
{ 23, 25852016738884976640000.0,                  0 },
{ 24, 620448401733239439360000.0,                 0 },
{ 25, 15511210043330985984000000.0,               0 },
{ 26, 403291461126605635584000000.0,              0 },
{ 27, 10888869450418352160768000000.0,            0 },
{ 28, 304888344611713860501504000000.0,           0 },
{ 29, 8841761993739701954543616000000.0,          0 },
{ 30, 265252859812191058636308480000000.0,        0 },
{ 31, 8222838654177922817725562880000000.0,       0 },
{ 32, 263130836933693530167218012160000000.0,     0 },
{ 33, 8683317618811886495518194401280000000.0,    0 },
{ 34, 2.95232799039604140847618609644e38,  0 },
{ 35, 1.03331479663861449296666513375e40,  0 },
{ 36, 3.71993326789901217467999448151e41,  0 },
{ 37, 1.37637530912263450463159795816e43,  0 },
{ 38, 5.23022617466601111760007224100e44,  0 },
{ 39, 2.03978820811974433586402817399e46,  0 },
{ 40, 8.15915283247897734345611269600e47,  0 },
{ 41, 3.34525266131638071081700620534e49,  0 },
{ 42, 1.40500611775287989854314260624e51,  0 },
{ 43, 6.04152630633738356373551320685e52,  0 },
{ 44, 2.65827157478844876804362581101e54,  0 },
{ 45, 1.19622220865480194561963161496e56,  0 },
{ 46, 5.50262215981208894985030542880e57,  0 },
{ 47, 2.58623241511168180642964355154e59,  0 },
{ 48, 1.24139155925360726708622890474e61,  0 },
{ 49, 6.08281864034267560872252163321e62,  0 },
{ 50, 3.04140932017133780436126081661e64,  0 },
{ 51, 1.55111875328738228022424301647e66,  0 },
{ 52, 8.06581751709438785716606368564e67,  0 },
{ 53, 4.27488328406002556429801375339e69,  0 },
{ 54, 2.30843697339241380472092742683e71,  0 },
{ 55, 1.26964033536582759259651008476e73,  0 },
{ 56, 7.10998587804863451854045647464e74,  0 },
{ 57, 4.05269195048772167556806019054e76,  0 },
{ 58, 2.35056133128287857182947491052e78,  0 },
{ 59, 1.38683118545689835737939019720e80,  0 },
{ 60, 8.32098711274139014427634118320e81,  0 },
{ 61, 5.07580213877224798800856812177e83,  0 },
{ 62, 3.14699732603879375256531223550e85,  0 },
{ 63, 1.982608315404440064116146708360e87,  0 },
{ 64, 1.268869321858841641034333893350e89,  0 },
{ 65, 8.247650592082470666723170306800e90,  0 },
{ 66, 5.443449390774430640037292402480e92,  0 },
{ 67, 3.647111091818868528824985909660e94,  0 },
{ 68, 2.480035542436830599600990418570e96,  0 },
{ 69, 1.711224524281413113724683388810e98,  0 },
{ 70, 1.197857166996989179607278372170e100,  0 },
{ 71, 8.504785885678623175211676442400e101,  0 },
{ 72, 6.123445837688608686152407038530e103,  0 },
{ 73, 4.470115461512684340891257138130e105,  0 },
{ 74, 3.307885441519386412259530282210e107,  0 },
{ 75, 2.480914081139539809194647711660e109,  0 },
{ 76, 1.885494701666050254987932260860e111,  0 },
{ 77, 1.451830920282858696340707840860e113,  0 },
{ 78, 1.132428117820629783145752115870e115,  0 },
{ 79, 8.946182130782975286851441715400e116,  0 },
{ 80, 7.156945704626380229481153372320e118,  0 },
{ 81, 5.797126020747367985879734231580e120,  0 },
{ 82, 4.753643337012841748421382069890e122,  0 },
{ 83, 3.945523969720658651189747118010e124,  0 },
{ 84, 3.314240134565353266999387579130e126,  0 },
{ 85, 2.817104114380550276949479442260e128,  0 },
{ 86, 2.422709538367273238176552320340e130,  0 },
{ 87, 2.107757298379527717213600518700e132,  0 },
{ 88, 1.854826422573984391147968456460e134,  0 },
{ 89, 1.650795516090846108121691926250e136,  0 },
{ 90, 1.485715964481761497309522733620e138,  0 },
{ 91, 1.352001527678402962551665687590e140,  0 },
{ 92, 1.243841405464130725547532432590e142,  0 },
{ 93, 1.156772507081641574759205162310e144,  0 },
{ 94, 1.087366156656743080273652852570e146,  0 },
{ 95, 1.032997848823905926259970209940e148,  0 },
{ 96, 9.916779348709496892095714015400e149,  0 },
{ 97, 9.619275968248211985332842594960e151,  0 },
{ 98, 9.426890448883247745626185743100e153,  0 },
{ 99, 9.332621544394415268169923885600e155,  0 },
{ 100, 9.33262154439441526816992388563e157,  0 },
{ 101, 9.42594775983835942085162312450e159,  0 },
{ 102, 9.61446671503512660926865558700e161,  0 },
{ 103, 9.90290071648618040754671525458e163,  0 },
{ 104, 1.02990167451456276238485838648e166,  0 },
{ 105, 1.08139675824029090050410130580e168,  0 },
{ 106, 1.146280563734708354534347384148e170,  0 },
{ 107, 1.226520203196137939351751701040e172,  0 },
{ 108, 1.324641819451828974499891837120e174,  0 },
{ 109, 1.443859583202493582204882102460e176,  0 },
{ 110, 1.588245541522742940425370312710e178,  0 },
{ 111, 1.762952551090244663872161047110e180,  0 },
{ 112, 1.974506857221074023536820372760e182,  0 },
{ 113, 2.231192748659813646596607021220e184,  0 },
{ 114, 2.543559733472187557120132004190e186,  0 },
{ 115, 2.925093693493015690688151804820e188,  0 },
{ 116, 3.393108684451898201198256093590e190,  0 },
{ 117, 3.96993716080872089540195962950e192,  0 },
{ 118, 4.68452584975429065657431236281e194,  0 },
{ 119, 5.57458576120760588132343171174e196,  0 },
{ 120, 6.68950291344912705758811805409e198,  0 },
{ 121, 8.09429852527344373968162284545e200,  0 },
{ 122, 9.87504420083360136241157987140e202,  0 },
{ 123, 1.21463043670253296757662432419e205,  0 },
{ 124, 1.50614174151114087979501416199e207,  0 },
{ 125, 1.88267717688892609974376770249e209,  0 },
{ 126, 2.37217324288004688567714730514e211,  0 },
{ 127, 3.01266001845765954480997707753e213,  0 },
{ 128, 3.85620482362580421735677065923e215,  0 },
{ 129, 4.97450422247728744039023415041e217,  0 },
{ 130, 6.46685548922047367250730439554e219,  0 },
{ 131, 8.47158069087882051098456875820e221,  0 },
{ 132, 1.11824865119600430744996307608e224,  0 },
{ 133, 1.48727070609068572890845089118e226,  0 },
{ 134, 1.99294274616151887673732419418e228,  0 },
{ 135, 2.69047270731805048359538766215e230,  0 },
{ 136, 3.65904288195254865768972722052e232,  0 },
{ 137, 5.01288874827499166103492629211e234,  0 },
{ 138, 6.91778647261948849222819828311e236,  0 },
{ 139, 9.61572319694108900419719561353e238,  0 },
{ 140, 1.34620124757175246058760738589e241,  0 },
{ 141, 1.89814375907617096942852641411e243,  0 },
{ 142, 2.69536413788816277658850750804e245,  0 },
{ 143, 3.85437071718007277052156573649e247,  0 },
{ 144, 5.55029383273930478955105466055e249,  0 },
{ 145, 8.04792605747199194484902925780e251,  0 },
{ 146, 1.17499720439091082394795827164e254,  0 },
{ 147, 1.72724589045463891120349865931e256,  0 },
{ 148, 2.55632391787286558858117801578e258,  0 },
{ 149, 3.80892263763056972698595524351e260,  0 },
{ 150, 5.71338395644585459047893286526e262,  0 },
{ 151, 8.62720977423324043162318862650e264,  0 },
{ 152, 1.31133588568345254560672467123e267,  0 },
{ 153, 2.00634390509568239477828874699e269,  0 },
{ 154, 3.08976961384735088795856467036e271,  0 },
{ 155, 4.78914290146339387633577523906e273,  0 },
{ 156, 7.47106292628289444708380937294e275,  0 },
{ 157, 1.17295687942641442819215807155e278,  0 },
{ 158, 1.85327186949373479654360975305e280,  0 },
{ 159, 2.94670227249503832650433950735e282,  0 },
{ 160, 4.71472363599206132240694321176e284,  0 },
{ 161, 7.59070505394721872907517857094e286,  0 },
{ 162, 1.22969421873944943411017892849e289,  0 },
{ 163, 2.00440157654530257759959165344e291,  0 },
{ 164, 3.28721858553429622726333031164e293,  0 },
{ 165, 5.42391066613158877498449501421e295,  0 },
{ 166, 9.00369170577843736647426172359e297,  0 },
{ 167, 1.50361651486499904020120170784e300,  0 },
{ 168, 2.52607574497319838753801886917e302,  0 },
{ 169, 4.26906800900470527493925188890e304,  0 },
{ 170, 7.25741561530799896739672821113e306,  0 },
/*
 { 171, 1.24101807021766782342484052410e309,  0 },
 { 172, 2.13455108077438865629072570146e311,  0 },
 { 173, 3.69277336973969237538295546352e313,  0 },
 { 174, 6.42542566334706473316634250653e315,  0 },
 { 175, 1.12444949108573632830410993864e318,  0 },
 { 176, 1.97903110431089593781523349201e320,  0 },
 { 177, 3.50288505463028580993296328086e322,  0 },
 { 178, 6.23513539724190874168067463993e324,  0 },
 { 179, 1.11608923610630166476084076055e327,  0 },
 { 180, 2.00896062499134299656951336898e329,  0 },
 { 181, 3.63621873123433082379081919786e331,  0 },
 { 182, 6.61791809084648209929929094011e333,  0 },
 { 183, 1.21107901062490622417177024204e336,  0 },
 { 184, 2.22838537954982745247605724535e338,  0 },
 { 185, 4.12251295216718078708070590390e340,  0 },
 { 186, 7.66787409103095626397011298130e342,  0 },
 { 187, 1.43389245502278882136241112750e345,  0 },
 { 188, 2.69571781544284298416133291969e347,  0 },
 { 189, 5.09490667118697324006491921822e349,  0 },
 { 190, 9.68032267525524915612334651460e351,  0 },
 { 191, 1.84894163097375258881955918429e354,  0 },
 { 192, 3.54996793146960497053355363384e356,  0 },
 { 193, 6.85143810773633759312975851330e358,  0 },
 { 194, 1.32917899290084949306717315158e361,  0 },
 { 195, 2.59189903615665651148098764559e363,  0 },
 { 196, 5.08012211086704676250273578535e365,  0 },
 { 197, 1.00078405584080821221303894971e368,  0 },
 { 198, 1.98155243056480026018181712043e370,  0 },
 { 199, 3.94328933682395251776181606966e372,  0 },
 { 200, 7.88657867364790503552363213932e374,  0 }
 */
};
static struct {int n; double f; long i; } doub_fact_table[GSL_SF_DOUBLEFACT_NMAX + 1] = {
{ 0,  1.000000000000000000000000000,    1L    },
{ 1,  1.000000000000000000000000000,    1L    },
{ 2,  2.000000000000000000000000000,    2L    },
{ 3,  3.000000000000000000000000000,    3L    },
{ 4,  8.000000000000000000000000000,    8L    },
{ 5,  15.00000000000000000000000000,    15L   },
{ 6,  48.00000000000000000000000000,    48L   },
{ 7,  105.0000000000000000000000000,    105L  },
{ 8,  384.0000000000000000000000000,    384L  },
{ 9,  945.0000000000000000000000000,    945L  },
{ 10, 3840.000000000000000000000000,    3840L   },
{ 11, 10395.00000000000000000000000,    10395L  },
{ 12, 46080.00000000000000000000000,       46080L       },
{ 13, 135135.0000000000000000000000,       135135L      },
{ 14, 645120.00000000000000000000000,      645120L      },
{ 15, 2.02702500000000000000000000000e6,   2027025L     },
{ 16, 1.03219200000000000000000000000e7,   10321920L    },
{ 17, 3.4459425000000000000000000000e7,    34459425L    },
{ 18, 1.85794560000000000000000000000e8,   185794560L   },
{ 19, 6.5472907500000000000000000000e8,            0 },
{ 20, 3.7158912000000000000000000000e9,            0 },
{ 21, 1.37493105750000000000000000000e10,          0 },
{ 22, 8.1749606400000000000000000000e10,           0 },
{ 23, 3.1623414322500000000000000000e11,           0 },
{ 24, 1.96199055360000000000000000000e12,          0 },
{ 25, 7.9058535806250000000000000000e12,           0 },
{ 26, 5.1011754393600000000000000000e13,           0 },
{ 27, 2.13458046676875000000000000000e14,          0 },
{ 28, 1.42832912302080000000000000000e15,          0 },
{ 29, 6.1902833536293750000000000000e15,           0 },
{ 30, 4.2849873690624000000000000000e16,           0 },
{ 31, 1.91898783962510625000000000000e17,          0 },
{ 32, 1.37119595809996800000000000000e18,          0 },
{ 33, 6.3326598707628506250000000000e18,           0 },
{ 34, 4.6620662575398912000000000000e19,           0 },
{ 35, 2.21643095476699771875000000000e20,          0 },
{ 36, 1.67834385271436083200000000000e21,          0 },
{ 37, 8.2007945326378915593750000000e21,           0 },
{ 38, 6.3777066403145711616000000000e22,           0 },
{ 39, 3.1983098677287777081562500000e23,           0 },
{ 40, 2.55108265612582846464000000000e24,          0 },
{ 41, 1.31130704576879886034406250000e25,          0 },
{ 42, 1.07145471557284795514880000000e26,          0 },
{ 43, 5.6386202968058350994794687500e26,           0 },
{ 44, 4.7144007485205310026547200000e27,           0 },
{ 45, 2.53737913356262579476576093750e28,          0 },
{ 46, 2.16862434431944426122117120000e29,          0 },
{ 47, 1.19256819277443412353990764062e30,          0 },
{ 48, 1.04093968527333324538616217600e31,          0 },
{ 49, 5.8435841445947272053455474391e31,           0 },
{ 50, 5.2046984263666662269308108800e32,           0 },
{ 51, 2.98022791374331087472622919392e33,          0 },
{ 52, 2.70644318171066643800402165760e34,          0 },
{ 53, 1.57952079428395476360490147278e35,          0 },
{ 54, 1.46147931812375987652217169510e36,          0 },
{ 55, 8.6873643685617511998269581003e36,           0 },
{ 56, 8.1842841814930553085241614926e37,           0 },
{ 57, 4.9517976900801981839013661172e38,           0 },
{ 58, 4.7468848252659720789440136657e39,           0 },
{ 59, 2.92156063714731692850180600912e40,       0 },
{ 60, 2.84813089515958324736640819942e41,       0 },
{ 61, 1.78215198865986332638610166557e42,       0 },
{ 62, 1.76584115499894161336717308364e43,       0 },
{ 63, 1.12275575285571389562324404931e44,       0 },
{ 64, 1.13013833919932263255499077353e45,       0 },
{ 65, 7.2979123935621403215510863205e45,        0 },
{ 66, 7.4589130387155293748629391053e46,        0 },
{ 67, 4.8896013036866340154392278347e47,        0 },
{ 68, 5.0720608663265599749067985916e48,        0 },
{ 69, 3.3738248995437774706530672060e49,        0 },
{ 70, 3.5504426064285919824347590141e50,        0 },
{ 71, 2.39541567867608200416367771623e51,       0 },
{ 72, 2.55631867662858622735302649017e52,       0 },
{ 73, 1.74865344543353986303948473285e53,       0 },
{ 74, 1.89167582070515380824123960272e54,       0 },
{ 75, 1.31149008407515489727961354964e55,       0 },
{ 76, 1.43767362373591689426334209807e56,       0 },
{ 77, 1.00984736473786927090530243322e57,       0 },
{ 78, 1.12138542651401517752540683649e58,       0 },
{ 79, 7.9777941814291672401518892225e58,        0 },
{ 80, 8.9710834121121214202032546920e59,        0 },
{ 81, 6.4620132869576254645230302702e60,        0 },
{ 82, 7.3562883979319395645666688474e61,        0 },
{ 83, 5.3634710281748291355541151243e62,        0 },
{ 84, 6.1792822542628292342360018318e63,        0 },
{ 85, 4.5589503739486047652209978556e64,        0 },
{ 86, 5.3141827386660331414429615754e65,        0 },
{ 87, 3.9662868253352861457422681344e66,        0 },
{ 88, 4.6764808100261091644698061863e67,        0 },
{ 89, 3.5299952745484046697106186396e68,        0 },
{ 90, 4.2088327290234982480228255677e69,        0 },
{ 91, 3.2122956998390482494366629620e70,        0 },
{ 92, 3.8721261107016183881809995223e71,        0 },
{ 93, 2.98743500085031487197609655470e72,       0 },
{ 94, 3.6397985440595212848901395509e73,        0 },
{ 95, 2.83806325080779912837729172696e74,       0 },
{ 96, 3.4942066022971404334945339689e75,        0 },
{ 97, 2.75292135328356515452597297515e76,       0 },
{ 98, 3.4243224702511976248246432895e77,        0 },
{ 99, 2.72539213975072950298071324540e78,       0 },
{ 100, 3.4243224702511976248246432895e79,       0 },
{ 101, 2.75264606114823679801052037785e80,      0 },
{ 102, 3.4928089196562215773211361553e81,       0 },
{ 103, 2.83522544298268390195083598919e82,      0 },
{ 104, 3.6325212764424704404139816015e83,       0 },
{ 105, 2.97698671513181809704837778865e84,      0 },
{ 106, 3.8504725530290186668388204976e85,       0 },
{ 107, 3.1853757851910453638417642339e86,       0 },
{ 108, 4.1585103572713401601859261374e87,       0 },
{ 109, 3.4720596058582394465875230149e88,       0 },
{ 110, 4.5743613929984741762045187512e89,       0 },
{ 111, 3.8539861625026457857121505465e90,       0 },
{ 112, 5.1232847601582910773490610013e91,       0 },
{ 113, 4.3550043636279897378547301176e92,       0 },
{ 114, 5.8405446265804518281779295415e93,       0 },
{ 115, 5.0082550181721881985329396352e94,       0 },
{ 116, 6.7750317668333241206863982681e95,       0 },
{ 117, 5.8596583712614601922835393732e96,       0 },
{ 118, 7.9945374848633224624099499564e97,       0 },
{ 119, 6.9729934618011376288174118541e98,       0 },
{ 120, 9.5934449818359869548919399477e99,       0 },
{ 121, 8.4373220887793765308690683435e100,      0 },
{ 122, 1.17040028778399040849681667362e102,       0 },
{ 123, 1.03779061691986331329689540625e103,       0 },
{ 124, 1.45129635685214810653605267528e104,       0 },
{ 125, 1.29723827114982914162111925781e105,       0 },
{ 126, 1.82863340963370661423542637086e106,       0 },
{ 127, 1.64749260436028300985882145742e107,       0 },
{ 128, 2.34065076433114446622134575470e108,       0 },
{ 129, 2.12526545962476508271787968008e109,       0 },
{ 130, 3.04284599363048780608774948111e110,       0 },
{ 131, 2.78409775210844225836042238090e111,       0 },
{ 132, 4.0165567115922439040358293151e112,        0 },
{ 133, 3.7028500103042282036193617666e113,        0 },
{ 134, 5.3821859935336068314080112822e114,        0 },
{ 135, 4.9988475139107080748861383849e115,        0 },
{ 136, 7.3197729512057052907148953438e116,        0 },
{ 137, 6.8484210940576700625940095873e117,        0 },
{ 138, 1.01012866726638733011865555744e119,       0 },
{ 139, 9.5193053207401613870056733264e119,        0 },
{ 140, 1.41418013417294226216611778042e121,       0 },
{ 141, 1.34222205022436275556779993902e122,       0 },
{ 142, 2.00813579052557801227588724819e123,       0 },
{ 143, 1.91937753182083874046195391280e124,       0 },
{ 144, 2.89171553835683233767727763739e125,       0 },
{ 145, 2.78309742114021617366983317355e126,       0 },
{ 146, 4.2219046860009752130088253506e127,        0 },
{ 147, 4.0911532090761177752946547651e128,        0 },
{ 148, 6.2484189352814433152530615189e129,        0 },
{ 149, 6.0958182815234154851890356000e130,        0 },
{ 150, 9.3726284029221649728795922783e131,        0 },
{ 151, 9.2046856051003573826354437561e132,        0 },
{ 152, 1.42463951724416907587769802630e134,       0 },
{ 153, 1.40831689758035467954322289468e135,       0 },
{ 154, 2.19394485655602037685165496051e136,       0 },
{ 155, 2.18289119124954975329199548675e137,       0 },
{ 156, 3.4225539762273917878885817384e138,        0 },
{ 157, 3.4271391702617931126684329142e139,        0 },
{ 158, 5.4076352824392790248639591467e140,        0 },
{ 159, 5.4491512807162510491428083336e141,        0 },
{ 160, 8.6522164519028464397823346347e142,        0 },
{ 161, 8.7731335619531641891199214170e143,        0 },
{ 162, 1.40165906520826112324473821082e145,       0 },
{ 163, 1.43002077059836576282654719098e146,       0 },
{ 164, 2.29872086694154824212137066574e147,       0 },
{ 165, 2.35953427148730350866380286512e148,       0 },
{ 166, 3.8158766391229700819214753051e149,        0 },
{ 167, 3.9404222333837968594685507847e150,        0 },
{ 168, 6.4106727537265897376280785126e151,        0 },
{ 169, 6.6593135744186166925018508262e152,        0 },
{ 170, 1.08981436813352025539677334714e154,       0 },
{ 171, 1.13874262122558345441781649128e155,       0 },
{ 172, 1.87448071318965483928245015709e156,       0 },
{ 173, 1.97002473472025937614282252992e157,       0 },
{ 174, 3.2615964409499994203514632733e158,        0 },
{ 175, 3.4475432857604539082499394274e159,        0 },
{ 176, 5.7404097360719989798185753611e160,        0 },
{ 177, 6.1021516157960034176023927864e161,        0 },
{ 178, 1.02179293302081581840770641427e163,       0 },
{ 179, 1.09228513922748461175082830877e164,       0 },
{ 180, 1.83922727943746847313387154568e165,       0 },
{ 181, 1.97703610200174714726899923887e166,       0 },
{ 182, 3.3473936485761926211036462131e167,        0 },
{ 183, 3.6179760666631972795022686071e168,        0 },
{ 184, 6.1592043133801944228307090322e169,        0 },
{ 185, 6.6932557233269149670791969232e170,        0 },
{ 186, 1.14561200228871616264651187999e172,       0 },
{ 187, 1.25163882026213309884380982464e173,       0 },
{ 188, 2.15375056430278638577544233437e174,       0 },
{ 189, 2.36559737029543155681480056857e175,       0 },
{ 190, 4.0921260721752941329733404353e176,        0 },
{ 191, 4.5182909772642742735162690860e177,        0 },
{ 192, 7.8568820585765647353088136358e178,        0 },
{ 193, 8.7203015861200493478863993359e179,        0 },
{ 194, 1.52423511936385355864990984535e181,       0 },
{ 195, 1.70045880929340962283784787050e182,       0 },
{ 196, 2.98750083395315297495382329688e183,       0 },
{ 197, 3.3499038543080169569905603049e184,        0 },
{ 198, 5.9152516512272428904085701278e185,        0 },
{ 199, 6.6663086700729537444112150067e186,        0 },
{ 200, 1.18305033024544857808171402556e188,       0 },
{ 201, 1.33992804268466370262665421635e189,       0 },
{ 202, 2.38976166709580612772506233164e190,       0 },
{ 203, 2.72005392664986731633210805920e191,       0 },
{ 204, 4.8751138008754445005591271565e192,        0 },
{ 205, 5.5761105496322279984808215214e193,        0 },
{ 206, 1.00427344298034156711518019425e195,       0 },
{ 207, 1.15425488377387119568553005492e196,       0 },
{ 208, 2.08888876139911045959957480403e197,       0 },
{ 209, 2.41239270708739079898275781478e198,       0 },
{ 210, 4.3866663989381319651591070885e199,        0 },
{ 211, 5.0901486119543945858536189892e200,        0 },
{ 212, 9.2997327657488397661373070276e201,        0 },
{ 213, 1.08420165434628604678682084470e203,       0 },
{ 214, 1.99014281187025170995338370390e204,       0 },
{ 215, 2.33103355684451500059166481610e205,       0 },
{ 216, 4.2987084736397436934993088004e206,        0 },
{ 217, 5.0583428183525975512839126509e207,        0 },
{ 218, 9.3711844725346412518284931849e208,        0 },
{ 219, 1.10777707721921886373117687056e210,       0 },
{ 220, 2.06166058395762107540226850068e211,       0 },
{ 221, 2.44818734065447368884590088393e212,       0 },
{ 222, 4.5768864963859187873930360715e213,        0 },
{ 223, 5.4594577696594763261263589712e214,        0 },
{ 224, 1.02522257519044580837604008002e216,       0 },
{ 225, 1.22837799817338217337843076851e217,       0 },
{ 226, 2.31700301993040752692985058084e218,       0 },
{ 227, 2.78841805585357753356903784452e219,       0 },
{ 228, 5.2827668854413291614000593243e220,        0 },
{ 229, 6.3854773479046925518730966640e221,        0 },
{ 230, 1.21503638365150570712201364459e223,       0 },
{ 231, 1.47504526736598397948268532937e224,       0 },
{ 232, 2.81888441007149324052307165546e225,       0 },
{ 233, 3.4368554729627426721946568174e226,        0 },
{ 234, 6.5961895195672941828239876738e227,        0 },
{ 235, 8.0766103614624452796574435210e228,        0 },
{ 236, 1.55670072661788142714646109101e230,       0 },
{ 237, 1.91415665566659953127881411447e231,       0 },
{ 238, 3.7049477293505577966085773966e232,        0 },
{ 239, 4.5748344070431728797563657336e233,        0 },
{ 240, 8.8918745504413387118605857518e234,        0 },
{ 241, 1.10253509209740466402128414180e236,       0 },
{ 242, 2.15183364120680396827026175195e237,       0 },
{ 243, 2.67916027379669333357172046456e238,       0 },
{ 244, 5.2504740845446016825794386748e239,        0 },
{ 245, 6.5639426708018986672507151382e240,        0 },
{ 246, 1.29161662479797201391454191399e242,       0 },
{ 247, 1.62129383968806897081092663913e243,       0 },
{ 248, 3.2032092294989705945080639467e244,        0 },
{ 249, 4.0370216608232917373192073314e245,        0 },
{ 250, 8.0080230737474264862701598667e246,        0 },
{ 251, 1.01329243686664622606712104019e248,       0 },
{ 252, 2.01802181458435147454008028642e249,       0 },
{ 253, 2.56362986527261495194981623168e250,       0 },
{ 254, 5.1257754090442527453318039275e251,        0 },
{ 255, 6.5372561564451681274720313908e252,        0 },
{ 256, 1.31219850471532870280494180544e254,       0 },
{ 257, 1.68007483220640820876031206743e255,       0 },
{ 258, 3.3854721421655480532367498580e256,        0 },
{ 259, 4.3513938154145972606892082546e257,        0 },
{ 260, 8.8022275696304249384155496309e258,        0 },
{ 261, 1.13571378582320988503988335446e260,       0 },
{ 262, 2.30618362324317133386487400329e261,       0 },
{ 263, 2.98692725671504199765489322224e262,       0 },
{ 264, 6.0883247653619723214032673687e263,        0 },
{ 265, 7.9153572302948612937854670389e264,        0 },
{ 266, 1.61949438758628463749326912007e266,       0 },
{ 267, 2.11340038048872796544071969939e267,       0 },
{ 268, 4.3402449587312428284819612418e268,        0 },
{ 269, 5.6850470235146782270355359914e269,        0 },
{ 270, 1.17186613885743556369012953528e271,       0 },
{ 271, 1.54064774337247779952663025366e272,       0 },
{ 272, 3.1874758976922247332371523360e273,        0 },
{ 273, 4.2059683394068643927077005925e274,        0 },
{ 274, 8.7336839596766957690697974006e275,        0 },
{ 275, 1.15664129333688770799461766294e277,       0 },
{ 276, 2.41049677287076803226326408256e278,       0 },
{ 277, 3.2038963825431789511450909263e279,        0 },
{ 278, 6.7011810285807351296918741495e280,        0 },
{ 279, 8.9388709072954692736948036845e281,        0 },
{ 280, 1.87633068800260583631372476186e283,       0 },
{ 281, 2.51182272495002686590823983534e284,       0 },
{ 282, 5.2912525401673484584047038284e285,        0 },
{ 283, 7.1084583116085760305203187340e286,        0 },
{ 284, 1.50271572140752696218693588728e288,       0 },
{ 285, 2.02591061880844416869829083919e289,       0 },
{ 286, 4.2977669632255271118546366376e290,        0 },
{ 287, 5.8143634759802347641640947085e291,        0 },
{ 288, 1.23775688540895180821413535163e293,       0 },
{ 289, 1.68035104455828784684342337075e294,       0 },
{ 290, 3.5894949676859602438209925197e295,        0 },
{ 291, 4.8898215396646176343143620089e296,        0 },
{ 292, 1.04813253056430039119572981576e298,       0 },
{ 293, 1.43271771112173296685410806860e299,       0 },
{ 294, 3.08150963985904315011544565835e300,       0 },
{ 295, 4.2265172478091122522196188024e301,        0 },
{ 296, 9.1212685339827677243417191487e302,        0 },
{ 297, 1.25527562259930633890922678431e304,       0 },
/*
 { 298, 2.71813802312686478185383230631e305,       0 },
 { 299, 3.7532741115719259533385880851e306,        0 },
 { 300, 8.1544140693805943455614969189e307,  }
 */
};
/* Chebyshev coefficients for Gamma*(3/4(t+1)+1/2), -1<t<1
 */
static double gstar_a_data[30] = {
2.16786447866463034423060819465,
-0.05533249018745584258035832802,
0.01800392431460719960888319748,
-0.00580919269468937714480019814,
0.00186523689488400339978881560,
-0.00059746524113955531852595159,
0.00019125169907783353925426722,
-0.00006124996546944685735909697,
0.00001963889633130842586440945,
-6.3067741254637180272515795142e-06,
2.0288698405861392526872789863e-06,
-6.5384896660838465981983750582e-07,
2.1108698058908865476480734911e-07,
-6.8260714912274941677892994580e-08,
2.2108560875880560555583978510e-08,
-7.1710331930255456643627187187e-09,
2.3290892983985406754602564745e-09,
-7.5740371598505586754890405359e-10,
2.4658267222594334398525312084e-10,
-8.0362243171659883803428749516e-11,
2.6215616826341594653521346229e-11,
-8.5596155025948750540420068109e-12,
2.7970831499487963614315315444e-12,
-9.1471771211886202805502562414e-13,
2.9934720198063397094916415927e-13,
-9.8026575909753445931073620469e-14,
3.2116773667767153777571410671e-14,
-1.0518035333878147029650507254e-14,
3.4144405720185253938994854173e-15,
-1.0115153943081187052322643819e-15
};
static cheb_series gstar_a_cs = {
gstar_a_data,
29,
-1, 1,
17
};
/* series for gammastar(x)
 * double-precision for x > 10.0
 */
static
int	gammastar_ser(const double x, gsl_sf_result * result){
	/* Use the Stirling series for the correction to Log(Gamma(x)),
	 * which is better behaved and easier to compute than the
	 * regular Stirling series for Gamma(x).
	 */
	const double y = 1.0/(x*x);
	const double c0 =  1.0/12.0;
	const double c1 = -1.0/360.0;
	const double c2 =  1.0/1260.0;
	const double c3 = -1.0/1680.0;
	const double c4 =  1.0/1188.0;
	const double c5 = -691.0/360360.0;
	const double c6 =  1.0/156.0;
	const double c7 = -3617.0/122400.0;
	const double ser = c0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
	result->val = exp(ser/x);
	result->err = 2.0 * GSL_DBL_EPSILON * result->val * ((1.0 > ser/x) ? 1.0 : ser/x);
	return GSL_SUCCESS;
}
/* Chebyshev expansion for log(gamma(x)/gamma(8))
 * 5 < x < 10
 * -1 < t < 1
 */
static double gamma_5_10_data[24] = {
-1.5285594096661578881275075214,
4.8259152300595906319768555035,
0.2277712320977614992970601978,
-0.0138867665685617873604917300,
0.0012704876495201082588139723,
-0.0001393841240254993658962470,
0.0000169709242992322702260663,
-2.2108528820210580075775889168e-06,
3.0196602854202309805163918716e-07,
-4.2705675000079118380587357358e-08,
6.2026423818051402794663551945e-09,
-9.1993973208880910416311405656e-10,
1.3875551258028145778301211638e-10,
-2.1218861491906788718519522978e-11,
3.2821736040381439555133562600e-12,
-5.1260001009953791220611135264e-13,
8.0713532554874636696982146610e-14,
-1.2798522376569209083811628061e-14,
2.0417711600852502310258808643e-15,
-3.2745239502992355776882614137e-16,
5.2759418422036579482120897453e-17,
-8.5354147151695233960425725513e-18,
1.3858639703888078291599886143e-18,
-2.2574398807738626571560124396e-19
};
static const cheb_series gamma_5_10_cs = {
gamma_5_10_data,
23,
-1, 1,
11
};
/* gamma(x) for x >= 1/2
 * assumes x >= 1/2
 */
static
int	gamma_xgthalf(const double x, gsl_sf_result * result){

	if(x == 0.5) {
		result->val = 1.77245385090551602729817;
		result->err = GSL_DBL_EPSILON * result->val;
		return GSL_SUCCESS;
	} else if (x <= (GSL_SF_FACT_NMAX + 1.0) && x == floor(x)) {
		int n = (int) floor (x);
		result->val = fact_table[n - 1].f;
		result->err = GSL_DBL_EPSILON * result->val;
		return GSL_SUCCESS;
	}
	else if(fabs(x - 1.0) < 0.01) {
		/* Use series for Gamma[1+eps] - 1/(1+eps).
		 */
		const double eps = x - 1.0;
		const double c1 =  0.4227843350984671394;
		const double c2 = -0.01094400467202744461;
		const double c3 =  0.09252092391911371098;
		const double c4 = -0.018271913165599812664;
		const double c5 =  0.018004931096854797895;
		const double c6 = -0.006850885378723806846;
		const double c7 =  0.003998239557568466030;
		result->val = 1.0/x + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*c7))))));
		result->err = GSL_DBL_EPSILON;
		return GSL_SUCCESS;
	}
	else if(fabs(x - 2.0) < 0.01) {
		/* Use series for Gamma[1 + eps].
		 */
		const double eps = x - 2.0;
		const double c1 =  0.4227843350984671394;
		const double c2 =  0.4118403304264396948;
		const double c3 =  0.08157691924708626638;
		const double c4 =  0.07424901075351389832;
		const double c5 = -0.00026698206874501476832;
		const double c6 =  0.011154045718130991049;
		const double c7 = -0.002852645821155340816;
		const double c8 =  0.0021039333406973880085;
		result->val = 1.0 + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
		result->err = GSL_DBL_EPSILON;
		return GSL_SUCCESS;
	}
	else if(x < 5.0) {
		/* Exponentiating the logarithm is fine, as
		 * long as the exponential is not so large
		 * that it greatly amplifies the error.
		 */
		gsl_sf_result lg;
		lngamma_lanczos(x, &lg);
		result->val = exp(lg.val);
		result->err = result->val * (lg.err + 2.0 * GSL_DBL_EPSILON);
		return GSL_SUCCESS;
	}
	else if(x < 10.0) {
		/* This is a sticky area. The logarithm
		 * is too large and the gammastar series
		 * is not good.
		 */
		const double gamma_8 = 5040.0;
		const double t = (2.0*x - 15.0)/5.0;
		gsl_sf_result c;
		cheb_eval_e(&gamma_5_10_cs, t, &c);
		result->val  = exp(c.val) * gamma_8;
		result->err  = result->val * c.err;
		result->err += 2.0 * GSL_DBL_EPSILON * result->val;
		return GSL_SUCCESS;
	}
	else if(x < GSL_SF_GAMMA_XMAX) {
		/* We do not want to exponentiate the logarithm
		 * if x is large because of the inevitable
		 * inflation of the error. So we carefully
		 * use pow() and exp() with exact quantities.
		 */
		double p = pow(x, 0.5*x);
		double e = exp(-x);
		double q = (p * e) * p;
		double pre = M_SQRT2 * M_SQRTPI * q/sqrt(x);
		gsl_sf_result gstar;
		int stat_gs = gammastar_ser(x, &gstar);
		result->val = pre * gstar.val;
		result->err = (x + 2.5) * GSL_DBL_EPSILON * result->val;
		return stat_gs;
	}
	else {
		OVERFLOW_ERROR(result);
	}
}
/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double * sgn){
	if(fabs(x - 1.0) < 0.01) {
		int stat = lngamma_1_pade(x - 1.0, result_lg);
		result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
		*sgn = 1.0;
		return stat;
	}
	else if(fabs(x - 2.0) < 0.01) {
		int stat = lngamma_2_pade(x - 2.0, result_lg);
		result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
		*sgn = 1.0;
		return stat;
	}
	else if(x >= 0.5) {
		*sgn = 1.0;
		return lngamma_lanczos(x, result_lg);
	}
	else if(x == 0.0) {
		*sgn = 0.0;
		DOMAIN_ERROR(result_lg);
	}
	else if(fabs(x) < 0.02) {
		return lngamma_sgn_0(x, result_lg, sgn);
	}
	else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
		/* Try to extract a fractional
		 * part from x.
		 */
		double z = 1.0 - x;
		double s = sin(M_PI*x);
		double as = fabs(s);
		if(s == 0.0) {
			*sgn = 0.0;
			DOMAIN_ERROR(result_lg);
		}
		else if(as < M_PI*0.015) {
			/* x is near a negative integer, -N */
			if(x < ExCo<int>::mkMinimum() + 2.0) {
				result_lg->val = 0.0;
				result_lg->err = 0.0;
				*sgn = 0.0;
				GSL_ERROR ("error", GSL_EROUND);
			}
			else {
				int N = -(int)(x - 0.5);
				double eps = x + N;
				return lngamma_sgn_sing(N, eps, result_lg, sgn);
			}
		}
		else {
			gsl_sf_result lg_z;
			lngamma_lanczos(z, &lg_z);
			*sgn = (s > 0.0 ? 1.0 : -1.0);
			result_lg->val = log(M_PI) - (log(as) + lg_z.val);
			result_lg->err = 2.0 * GSL_DBL_EPSILON * fabs(result_lg->val) + lg_z.err;
			return GSL_SUCCESS;
		}
	}
	else {
		/* |x| was too large to extract any fractional part */
		result_lg->val = 0.0;
		result_lg->err = 0.0;
		*sgn = 0.0;
		GSL_ERROR ("error", GSL_EROUND);
	}
}
int	gsl_sf_gamma_e(const double x, gsl_sf_result * result){
	if(x < 0.5) {
		int rint_x = (int)floor(x+0.5);
		double f_x = x - rint_x;
		double sgn_gamma = ( ((rint_x & 1) ==0) ? 1.0 : -1.0 );
		double sin_term = sgn_gamma * sin(M_PI * f_x) / M_PI;

		if(sin_term == 0.0) {
			DOMAIN_ERROR(result);
		}
		else if(x > -169.0) {
			gsl_sf_result g;
			gamma_xgthalf(1.0-x, &g);
			if(fabs(sin_term) * g.val * GSL_DBL_MIN < 1.0) {
				result->val  = 1.0/(sin_term * g.val);
				result->err  = fabs(g.err/g.val) * fabs(result->val);
				result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
				return GSL_SUCCESS;
			}
			else {
				UNDERFLOW_ERROR(result);
			}
		}
		else {
			/* It is hard to control it here.
			 * We can only exponentiate the
			 * logarithm and eat the loss of
			 * precision.
			 */
			gsl_sf_result lng;
			double sgn;
			int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
			int stat_e   = gsl_sf_exp_mult_err_e(lng.val, lng.err, sgn, 0.0, result);
			return GSL_ERROR_SELECT_2(stat_e, stat_lng);
		}
	}
	else {
		return gamma_xgthalf(x, result);
	}
return GSL_FAILURE;
}
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result) {
	if(x <= 0.0) {
		DOMAIN_ERROR(result);
	}
	else if(x < 0.5) {
		gsl_sf_result lg;
		const int stat_lg = gsl_sf_lngamma_e(x, &lg);
		const double lx = log(x);
		const double c  = 0.5*log(2.0f * M_PI);
		const double lnr_val = lg.val - (x-0.5)*lx + x - c;
		const double lnr_err = lg.err + 2.0 * GSL_DBL_EPSILON *((x+0.5)*fabs(lx) + c);
		const int stat_e  = gsl_sf_exp_err_e(lnr_val, lnr_err, result);
		return GSL_ERROR_SELECT_2(stat_lg, stat_e);
	}
	else if(x < 2.0) {
		const double t = 4.0/3.0*(x-0.5) - 1.0;
		return cheb_eval_e(&gstar_a_cs, t, result);
	}
	else if(x < 10.0) {
		const double t = 0.25*(x-2.0) - 1.0;
		gsl_sf_result c;
		cheb_eval_e(&gstar_b_cs, t, &c);
		result->val  = c.val/(x*x) + 1.0 + 1.0/(12.0*x);
		result->err  = c.err/(x*x);
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if(x < 1.0/GSL_ROOT4_DBL_EPSILON) {
		return gammastar_ser(x, result);
	}
	else if(x < 1.0/GSL_DBL_EPSILON) {
		/* Use Stirling formula for Gamma(x).
		 */
		const double xi = 1.0/x;
		result->val = 1.0 + xi/12.0*(1.0 + xi/24.0*(1.0 - xi*(139.0/180.0 + 571.0/8640.0*xi)));
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		result->val = 1.0;
		result->err = 1.0/x;
		return GSL_SUCCESS;
	}
}
int	gsl_sf_gammainv_e(const double x, gsl_sf_result * result){

	if (x <= 0.0 && x == floor(x)) { /* negative integer */
		result->val = 0.0;
		result->err = 0.0;
		return GSL_SUCCESS;
	} else if(x < 0.5) {
		gsl_sf_result lng;
		double sgn;
		int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
		if(stat_lng == GSL_EDOM) {
			result->val = 0.0;
			result->err = 0.0;
			return GSL_SUCCESS;
		}
		else if(stat_lng != GSL_SUCCESS) {
			result->val = 0.0;
			result->err = 0.0;
			return stat_lng;
		}
		else {
			return gsl_sf_exp_mult_err_e(-lng.val, lng.err, sgn, 0.0, result);
		}
	}
	else {
		gsl_sf_result g;
		int stat_g = gamma_xgthalf(x, &g);
		if(stat_g == GSL_EOVRFLW) {
			UNDERFLOW_ERROR(result);
		}
		else {
			result->val  = 1.0/g.val;
			result->err  = fabs(g.err/g.val) * fabs(result->val);
			result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
			CHECK_UNDERFLOW(result);
			return GSL_SUCCESS;
		}
	}
return GSL_FAILURE;
}
/*
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg) {
	if(zr <= 0.5) {
		// Transform to right half plane using reflection; in fact we do a little better by stopping at 1/2.
		double x = 1.0-zr;
		double y = -zi;
		gsl_sf_result a, b;
		gsl_sf_result lnsin_r, lnsin_i;

		int stat_l = lngamma_lanczos_complex(x, y, &a, &b);
		int stat_s = gsl_sf_complex_logsin_e(M_PI*zr, M_PI*zi, &lnsin_r, &lnsin_i);

		if(stat_s == GSL_SUCCESS) {
			int stat_r;
			lnr->val = log(M_PI) - lnsin_r.val - a.val;
			lnr->err = lnsin_r.err + a.err + 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
			arg->val = -lnsin_i.val - b.val;
			arg->err = lnsin_i.err + b.err + 2.0 * GSL_DBL_EPSILON * fabs(arg->val);
			stat_r = gsl_sf_angle_restrict_symm_e(&(arg->val));
			return GSL_ERROR_SELECT_2(stat_r, stat_l);
		}
		else {
			DOMAIN_ERROR_2(lnr,arg);
		}
	}
	else {
		// otherwise plain vanilla Lanczos
		return lngamma_lanczos_complex(zr, zi, lnr, arg);
	}
}*/
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result){

	if(x < 0.0 || n < 0) DOMAIN_ERROR(result);
	else if(n == 0) {
		result->val = 1.0;
		result->err = 0.0;
	}else if(n == 1) {
		result->val = x;
		result->err = 0.0;
	}else if(x == 0.0) {
		result->val = 0.0;
		result->err = 0.0;
	}else{
		const double log2pi = log(2.0 * M_PI);
		const double ln_test = n*(log(x)+1.0) + 1.0 - (n+0.5)*log(n+1.0) + 0.5*log2pi;

		if (ln_test < GSL_LOG_DBL_MIN+1.0) {UNDERFLOW_ERROR(result);
		}else if(ln_test > GSL_LOG_DBL_MAX-1.0) {
			OVERFLOW_ERROR(result);
		}else{
			double product = 1.0;
			int k;
			for(k=1; k<=n; k++) {
				product *= (x/k);
			}
			result->val = product;
			result->err = n * GSL_DBL_EPSILON * product;
			CHECK_UNDERFLOW(result);
			return GSL_SUCCESS;
		}
	}
    return GSL_SUCCESS;
}
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result){

	if(n < 18) {
		result->val = fact_table[n].f;
		result->err = 0.0;
		return GSL_SUCCESS;
	}
	else if(n <= GSL_SF_FACT_NMAX){
		result->val = fact_table[n].f;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		OVERFLOW_ERROR(result);
	}
}
/*
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result){
	if(n < 26) {
		result->val = doub_fact_table[n].f;
		result->err = 0.0;
		return GSL_SUCCESS;
	}
	else if(n <= GSL_SF_DOUBLEFACT_NMAX){
		result->val = doub_fact_table[n].f;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		OVERFLOW_ERROR(result);
	}
}*/
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result){

	if(n <= GSL_SF_FACT_NMAX){
		result->val = log(fact_table[n].f);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else {
		gsl_sf_lngamma_e(n+1.0, result);
		return GSL_SUCCESS;
	}
}
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result){

	if(n <= GSL_SF_DOUBLEFACT_NMAX){
		result->val = log(doub_fact_table[n].f);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
	else if((n & 1)) {
		gsl_sf_result lg;
		gsl_sf_lngamma_e(0.5*(n+2.0), &lg);
		result->val = 0.5*(n+1.0) * M_LN2 - 0.5*log(M_PI) + lg.val;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
		return GSL_SUCCESS;
	}
	else {
		gsl_sf_result lg;
		gsl_sf_lngamma_e(0.5*n+1.0, &lg);
		result->val = 0.5*n*M_LN2 + lg.val;
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
		return GSL_SUCCESS;
	}
}
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result){

	if(m > n) {
		DOMAIN_ERROR(result);
	}
	else if(m == n || m == 0) {
		result->val = 0.0;
		result->err = 0.0;
		return GSL_SUCCESS;
	}
	else {
		gsl_sf_result nf;
		gsl_sf_result mf;
		gsl_sf_result nmmf;
		if(m*2 > n) m = n-m;
		gsl_sf_lnfact_e(n, &nf);
		gsl_sf_lnfact_e(m, &mf);
		gsl_sf_lnfact_e(n-m, &nmmf);
		result->val  = nf.val - mf.val - nmmf.val;
		result->err  = nf.err + mf.err + nmmf.err;
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
}
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result){
	if(m > n) {
		DOMAIN_ERROR(result);
	}
	else if(m == n || m == 0) {
		result->val = 1.0;
		result->err = 0.0;
		return GSL_SUCCESS;
	}
	else if (n <= GSL_SF_FACT_NMAX) {
		result->val = (fact_table[n].f / fact_table[m].f) / fact_table[n-m].f;
		result->err = 6.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	} else {
		if(m*2 < n) m = n-m;

		if (n - m < 64)  /* compute product for a manageable number of terms */
		{
			double prod = 1.0;
			unsigned int k;

			for(k=n; k>=m+1; k--) {
				double tk = (double)k / (double)(k-m);
				if(tk > GSL_DBL_MAX/prod) {
					OVERFLOW_ERROR(result);
				}
				prod *= tk;
			}
			result->val = prod;
			result->err = 2.0 * GSL_DBL_EPSILON * prod * fabs(n-m);
			return GSL_SUCCESS;
		}
		else
		{
			gsl_sf_result lc;
			const int stat_lc = gsl_sf_lnchoose_e (n, m, &lc);
			const int stat_e  = gsl_sf_exp_err_e(lc.val, lc.err, result);
			return GSL_ERROR_SELECT_2(stat_lc, stat_e);
		}
	}
}

static int	beta_cont_frac(const double a, const double b, const double x, gsl_sf_result * result){
	const unsigned int max_iter = 512;        /* control iterations      */
	const double cutoff = 2.0 * ExCo<double>::delta();  /* control the zero cutoff */
	unsigned int iter_count = 0;
	double cf;

	/* standard initialization for continued fraction */
	double num_term = 1.0;
	double den_term = 1.0 - (a+b)*x/(a+1.0);
	if (fabs(den_term) < cutoff) den_term = cutoff;
	den_term = 1.0/den_term;
	cf = den_term;

	while(iter_count < max_iter) {
		const int k  = iter_count + 1;
		double coeff = k*(b-k)*x/(((a-1.0)+2*k)*(a+2*k));
		double delta_frac;

		/* first step */
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		if(fabs(den_term) < cutoff) den_term = cutoff;
		if(fabs(num_term) < cutoff) num_term = cutoff;
		den_term  = 1.0/den_term;

		delta_frac = den_term * num_term;
		cf *= delta_frac;

		coeff = -(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1.0));

		/* second step */
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		if(fabs(den_term) < cutoff) den_term = cutoff;
		if(fabs(num_term) < cutoff) num_term = cutoff;
		den_term = 1.0/den_term;

		delta_frac = den_term*num_term;
		cf *= delta_frac;

		if(fabs(delta_frac-1.0) < 2.0*GSL_DBL_EPSILON) break;

		++iter_count;
	}

	result->val = cf;
	result->err = iter_count * 4.0 * GSL_DBL_EPSILON * fabs(cf);

	if(iter_count >= max_iter) {
        printf("beta_cont_frac could not converge in time, %e %e %e as input\n", a, b, x);
		GSL_ERROR ("error, too many iteration", GSL_EMAXITER);
	}else return 0;
}

int
gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result)
{
  if(a < 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == 0.0)
  {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x <= 0.5*a) {
    /* If the series is quick, do that. It is
     * robust and simple.
     */
    result->val  = 1.0 - gamma_inc_P_series(a, x);
    result->err  = 0.0;
    return GSL_SUCCESS;
  }
  else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
    /* Then try the difficult asymptotic regime.
     * This is the only way to do this region.
     */
    result->val = gamma_inc_Q_asymp_unif(a, x);
    return GSL_SUCCESS;
  }
  else if(a < 0.2 && x < 5.0) {
    /* Cancellations at small a must be handled
     * analytically; x should not be too big
     * either since the series terms grow
     * with x and log(x).
     */
    result->val = gamma_inc_Q_series(a, x);
    return GSL_SUCCESS;
  }
  else if(a <= x) {
    if(x <= 1.0e+06) {
      /* Continued fraction is excellent for x >~ a.
       * We do not let x be too large when x > a since
       * it is somewhat pointless to try this there;
       * the function is rapidly decreasing for
       * x large and x > a, and it will just
       * underflow in that region anyway. We
       * catch that case in the standard
       * large-x method.
       */
      result->val = gamma_inc_Q_CF(a, x);
    }
    else {
      result->val = gamma_inc_Q_large_x(a, x);
    }
    return GSL_SUCCESS;
  }
  else {
    if(x > a - sqrt(a)) {
      /* Continued fraction again. The convergence
       * is a little slower here, but that is fine.
       * We have to trade that off against the slow
       * convergence of the series, which is the
       * only other option.
       */
      result->val =  gamma_inc_Q_CF(a, x);
      return GSL_SUCCESS;
    }
    else {
      result->val  = gamma_inc_P_series(a, x);
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result)
{
  if(a <= 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 20.0 || x < 0.5*a) {
    /* Do the easy series cases. Robust and quick.
     */
    result->val = gamma_inc_P_series(a, x); return GSL_SUCCESS;
  }
  else if(a > 1.0e+06 && (x-a)*(x-a) < a) {
    /* Crossover region. Note that Q and P are
     * roughly the same order of magnitude here,
     * so the subtraction is stable.
     */
    result->val = gamma_inc_Q_asymp_unif(a, x);
    return GSL_SUCCESS;
  }
  else if(a <= x) {
    /* Q <~ P in this area, so the
     * subtractions are stable.
     */
    gsl_sf_result Q;
    int stat_Q;
    if(a > 0.2*x) {
      result->val = 1.0 - gamma_inc_Q_CF(a, x);
    }   else {
      result->val =  1.0 - gamma_inc_Q_large_x(a, x);
    }
    return GSL_SUCCESS;
  }
  else {
    if((x-a)*(x-a) < a) {
      /* This condition is meant to insure
       * that Q is not very close to 1,
       * so the subtraction is stable.
       */
      result->val = 1.0 -  gamma_inc_Q_CF(a, x);
    }
    else {
      result->val  = gamma_inc_P_series(a, x);
    }
    return GSL_SUCCESS;
  }
}

double gsl_sf_gamma_inc_P(const double a, const double x){
  gsl_sf_result result;
  gsl_sf_gamma_inc_P_e(a, x, &result);
  return result.val;
}

double gsl_sf_gamma_inc_Q(const double a, const double x)
{ gsl_sf_result result;
  gsl_sf_gamma_inc_Q_e(a, x, &result);
  return result.val;
}

static double isnegint (const double x) {return (x < 0) && (x == floor(x));}
int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn) {
	if(x == 0.0 || y == 0.0) {
		*sgn = 0.0;
		DOMAIN_ERROR(result);
	} else if (isnegint(x) || isnegint(y)) {
		*sgn = 0.0;
		DOMAIN_ERROR(result); // not defined for negative integers
	}

	// See if we can handle the postive case with min/max < 0.2

	if (x > 0 && y > 0) {
		const double max = (x > y) ? x : y;
		const double min = (x < y) ? x : y;
		const double rat = min/max;

		if(rat < 0.2) {
			// min << max, so be careful with the subtraction

			double lnpre_val;
			double lnpre_err;
			double lnpow_val;
			double lnpow_err;
			double t1, t2, t3;
			gsl_sf_result lnopr;
			gsl_sf_result gsx, gsy, gsxy;
			gsl_sf_gammastar_e(x, &gsx);
			gsl_sf_gammastar_e(y, &gsy);
			gsl_sf_gammastar_e(x+y, &gsxy);
			gsl_sf_log_1plusx_e(rat, &lnopr);
			lnpre_val = log(gsx.val*gsy.val/gsxy.val * M_SQRT2*M_SQRTPI);
			lnpre_err = gsx.err/gsx.val + gsy.err/gsy.val + gsxy.err/gsxy.val;
			t1 = min*log(rat);
			t2 = 0.5*log(min);
			t3 = (x+y-0.5)*lnopr.val;
			lnpow_val  = t1 - t2 - t3;
			lnpow_err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
			lnpow_err += fabs(x+y-0.5) * lnopr.err;
			result->val  = lnpre_val + lnpow_val;
			result->err  = lnpre_err + lnpow_err;
			result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
			*sgn = 1.0;
			return 0;
		}
	}

	// General case - Fallback
	{
		gsl_sf_result lgx, lgy, lgxy;
		double sgx, sgy, sgxy, xy = x+y;
		int stat_gx  = gsl_sf_lngamma_sgn_e(x, &lgx, &sgx);
		int stat_gy  = gsl_sf_lngamma_sgn_e(y, &lgy, &sgy);
		int stat_gxy = gsl_sf_lngamma_sgn_e(xy, &lgxy, &sgxy);
		*sgn = sgx * sgy * sgxy;
		result->val  = lgx.val + lgy.val - lgxy.val;
		result->err  = lgx.err + lgy.err + lgxy.err;
		result->err += GSL_DBL_EPSILON * (fabs(lgx.val) + fabs(lgy.val) + fabs(lgxy.val));
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_ERROR_SELECT_3(stat_gx, stat_gy, stat_gxy);
	}
}
int gsl_sf_lnbeta_e(const double x, const double y, gsl_sf_result * result) {
	double sgn;
	int status = gsl_sf_lnbeta_sgn_e(x,y,result,&sgn);
	if (sgn == -1) {
		DOMAIN_ERROR(result);
	}
	return status;
}
int gsl_sf_log_e(const double x, gsl_sf_result * result) {
	if(x <= 0.0) {
		DOMAIN_ERROR(result);
	}
	else {
		result->val = log(x);
		result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return GSL_SUCCESS;
	}
}
int	gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result){
	if(a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0) {
		DOMAIN_ERROR(result);
	}
	else if(x == 0.0) {
		result->val = 0.0;
		result->err = 0.0;
		return 0;
	}
	else if(x == 1.0) {
		result->val = 1.0;
		result->err = 0.0;
		return 0;
	}
	else {
		gsl_sf_result ln_beta;
		gsl_sf_result ln_x;
		gsl_sf_result ln_1mx;
		gsl_sf_result prefactor;
		const int stat_ln_beta = gsl_sf_lnbeta_e(a, b, &ln_beta);
		const int stat_ln_1mx = gsl_sf_log_1plusx_e(-x, &ln_1mx);
		const int stat_ln_x = gsl_sf_log_e(x, &ln_x);
		const int stat_ln = GSL_ERROR_SELECT_3(stat_ln_beta, stat_ln_1mx, stat_ln_x);

		const double ln_pre_val = -ln_beta.val + a * ln_x.val + b * ln_1mx.val;
		const double ln_pre_err =  ln_beta.err + fabs(a*ln_x.err) + fabs(b*ln_1mx.err);
		const int stat_exp = gsl_sf_exp_err_e(ln_pre_val, ln_pre_err, &prefactor);

		if(stat_ln != 0) {
			result->val = 0.0;
			result->err = 0.0;
			GSL_ERROR ("error", GSL_ESANITY);
		}

		if(x < (a + 1.0)/(a+b+2.0)) {
			/* Apply continued fraction directly. */
			gsl_sf_result cf;
			const int stat_cf = beta_cont_frac(a, b, x, &cf);
			int stat;
			result->val = prefactor.val * cf.val / a;
			result->err = (fabs(prefactor.err * cf.val) + fabs(prefactor.val * cf.err))/a;

			stat = GSL_ERROR_SELECT_2(stat_exp, stat_cf);
			if(stat == 0) {
				CHECK_UNDERFLOW(result);
			}
			return stat;
		}
		else {
			/* Apply continued fraction after hypergeometric transformation. */
			gsl_sf_result cf;
			const int stat_cf = beta_cont_frac(b, a, 1.0-x, &cf);
			int stat;
			const double term = prefactor.val * cf.val / b;
			result->val  = 1.0 - term;
			result->err  = fabs(prefactor.err * cf.val)/b;
			result->err += fabs(prefactor.val * cf.err)/b;
			result->err += 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(term));
			stat = GSL_ERROR_SELECT_2(stat_exp, stat_cf);
			if(stat == 0) {
				CHECK_UNDERFLOW(result);
			}
			return stat;
		}
	}
}
// gsl_sf_log_beta_inc_e added by LFH, computed "error" is off



double gsl_sf_lnbeta(const double x, const double y)
{ gsl_sf_result result;
  gsl_sf_lnbeta_e(x, y, &result);
  return result.val;
}
static double
beta_cont_frac_wr (const double a, const double b, const double x,
                const double epsabs)
{
  const unsigned int max_iter = 512;    /* control iterations      */
  const double cutoff = 2.0 * GSL_DBL_MIN;      /* control the zero cutoff */
  unsigned int iter_count = 0;
  double cf;

  /* standard initialization for continued fraction */
  double num_term = 1.0;
  double den_term = 1.0 - (a + b) * x / (a + 1.0);

  if (fabs (den_term) < cutoff)
    den_term = GSL_NAN;

  den_term = 1.0 / den_term;
  cf = den_term;

  while (iter_count < max_iter)
    {
      const int k = iter_count + 1;
      double coeff = k * (b - k) * x / (((a - 1.0) + 2 * k) * (a + 2 * k));
      double delta_frac;

      /* first step */
      den_term = 1.0 + coeff * den_term;
      num_term = 1.0 + coeff / num_term;

      if (fabs (den_term) < cutoff)
        den_term = GSL_NAN;

      if (fabs (num_term) < cutoff)
        num_term = GSL_NAN;

      den_term = 1.0 / den_term;

      delta_frac = den_term * num_term;
      cf *= delta_frac;

      coeff = -(a + k) * (a + b + k) * x / ((a + 2 * k) * (a + 2 * k + 1.0));

      /* second step */
      den_term = 1.0 + coeff * den_term;
      num_term = 1.0 + coeff / num_term;

      if (fabs (den_term) < cutoff)
        den_term = GSL_NAN;

      if (fabs (num_term) < cutoff)
        num_term = GSL_NAN;

      den_term = 1.0 / den_term;

      delta_frac = den_term * num_term;
      cf *= delta_frac;

      if (fabs (delta_frac - 1.0) < 2.0 * GSL_DBL_EPSILON)
        break;

      if (cf * fabs (delta_frac - 1.0) < epsabs)
        break;

      ++iter_count;
    }

  if (iter_count >= max_iter)
    return GSL_NAN;

  return cf;
}


static double
beta_inc_AXPY (const double A, const double Y,
               const double a, const double b, const double x)
{
  if (x == 0.0) return Y;
  else if (x == 1.0) return A + Y;
  else if (a > 1e5 && b < 10 && x > a / (a + b))
    {
      /* Handle asymptotic regime, large a, small b, x > peak [AS 26.5.17] */
      double N = a + (b - 1.0) / 2.0;
      return A * gsl_sf_gamma_inc_Q (b, -N * log (x)) + Y;
    }
  else if (b > 1e5 && a < 10 && x < b / (a + b))
    {
      /* Handle asymptotic regime, small a, large b, x < peak [AS 26.5.17] */
      double N = b + (a - 1.0) / 2.0;
      return A * gsl_sf_gamma_inc_P (a, -N * log1p (-x)) + Y;
    }
  else
    {
      double ln_beta = gsl_sf_lnbeta (a, b);
      double ln_pre = -ln_beta + a * log (x) + b * log1p (-x);

      double prefactor = exp (ln_pre);

      if (x < (a + 1.0) / (a + b + 2.0)){
          /* Apply continued fraction directly. */
          double epsabs = fabs (Y / (A * prefactor / a)) * GSL_DBL_EPSILON;
          double cf = beta_cont_frac_wr (a, b, x, epsabs);
          return A * (prefactor * cf / a) + Y;
      }else{
          /* Apply continued fraction after hypergeometric transformation. */
          double epsabs = fabs ((A + Y) / (A * prefactor / b)) * GSL_DBL_EPSILON;
          double cf = beta_cont_frac_wr (b, a, 1.0 - x, epsabs);
          double term = prefactor * cf / b;
          return (A == -Y) ? -A * term : A * (1 - term) + Y;
        }
    }
}


int	gsl_sf_log_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result){
	if(a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0) {
		DOMAIN_ERROR(result);
	}
	else if(x == 0.0) {
		result->val = GSL_NEGINF;
		result->err = 0.0;
		return 0;
	}
	else if(x == 1.0) {
		result->val = 0.0;
		result->err = 0.0;
		return 0;
	}else{
		gsl_sf_result ln_beta;
		gsl_sf_result ln_x;
		gsl_sf_result ln_1mx;

		const int stat_ln_beta = gsl_sf_lnbeta_e(a, b, &ln_beta);
		const int stat_ln_1mx = gsl_sf_log_1plusx_e(-x, &ln_1mx);
		const int stat_ln_x = gsl_sf_log_e(x, &ln_x);
		const int stat_ln = GSL_ERROR_SELECT_3(stat_ln_beta, stat_ln_1mx, stat_ln_x);

		const double ln_pre_val = -ln_beta.val + a * ln_x.val + b * ln_1mx.val;

		if(stat_ln != 0) {result->val = 0.0;result->err = 0.0;GSL_ERROR ("error", GSL_ESANITY);}

		if(x < (a + 1.0)/(a+b+2.0)) {
			/* Apply continued fraction directly. */
			gsl_sf_result cf;
			const int stat_cf = beta_cont_frac(a, b, x, &cf);
			int stat;
			const int stat_ln_xcf = gsl_sf_log_e(cf.val, result);
			result->val += ln_pre_val - log(a);
			result->err += ln_beta.err + fabs(a*ln_x.err) + fabs(b*ln_1mx.err);
			stat = GSL_ERROR_SELECT_2(stat_cf,stat_ln_xcf);
		//	if(stat == 0) {CHECK_UNDERFLOW(result);}
			return stat;
		}else{
			/* Apply continued fraction after hypergeometric transformation. */
			gsl_sf_result cf;
			gsl_sf_result prefactor;

			const double ln_pre_err =  ln_beta.err + fabs(a*ln_x.err) + fabs(b*ln_1mx.err);
			const int stat_exp = gsl_sf_exp_err_e(ln_pre_val, ln_pre_err, &prefactor);

			const int stat_cf = beta_cont_frac(b, a, 1.0-x, &cf);
			int stat;
			const double term = prefactor.val * cf.val / b;
			const int stat_ln_1mx_new = gsl_sf_log_1plusx_e(-term, result);
			result->err += fabs(prefactor.err * cf.val)/b;
			result->err += fabs(prefactor.val * cf.err)/b;

			stat = GSL_ERROR_SELECT_3(stat_exp, stat_cf, stat_ln_1mx_new);
		//	if(stat == 0) {CHECK_UNDERFLOW(result);}
			return stat;
		}
	}
}



double incgamma_frac(const double a, const double x) {
	if(a < 0.0 || x < 0.0) myexit("invalid Domain for incgamma_frac\n");
	else if(x == 0.0) return(1.0);
	else if(a == 0.0) return(0.0);
	else if(x <= 0.5*a) {
		// If the series is quick, do that. It is robust and simple.
		return 1.0f - gamma_inc_P_series(a, x);
	}
	else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
		/* Then try the difficult asymptotic regime.
		 * This is the only way to do this region.
		 */
		return gamma_inc_Q_asymp_unif(a, x);
	}
	else if(a < 0.2 && x < 5.0) {
		/* Cancellations at small a must be handled
		 * analytically; x should not be too big
		 * either since the series terms grow
		 * with x and log(x).
		 */
		return gamma_inc_Q_series(a, x);
	}
	else if(a <= x) {
		if(x <= 1.0e+06) {
			/* Continued fraction is excellent for x >~ a.
			 * We do not let x be too large when x > a since
			 * it is somewhat pointless to try this there;
			 * the function is rapidly decreasing for
			 * x large and x > a, and it will just
			 * underflow in that region anyway. We
			 * catch that case in the standard
			 * large-x method.
			 */
			return gamma_inc_Q_CF(a, x);
		}
		else {
			return gamma_inc_Q_large_x(a, x);
		}
	}
	else {
		if(x > a - sqrt(a)) {
			/* Continued fraction again. The convergence
			 * is a little slower here, but that is fine.
			 * We have to trade that off against the slow
			 * convergence of the series, which is the
			 * only other option.
			 */
			return gamma_inc_Q_CF(a, x);
		}
		else {
			return 1.0f - gamma_inc_P_series(a, x);
		}
	}
}

/* LFH log transformed */
double log_incgamma_dom(const double a, const double x) {
	if(a < 10.0) {
		double lnr;
		double lg;
		lg = lngamma(a+1.0);
		lnr = a * log(x) - x - lg;
		return lnr;
	}else{
		double gstar;
		double ln_term;

		if (x < 0.5*a) {
			double u = x/a;
			double ln_u = log(u);
			ln_term = ln_u - u + 1.0;
		} else {
			double mu = (x-a)/a;
			ln_term = log_1plusx_mx_e(mu);  /* log(1+mu) - mu */
		};
		gstar = gammastar(a);
		return  a*ln_term - log(gstar) - log(sqrt(2.0*M_PI*a));
	}
}

/* LFH log transformed */

double
log_gamma_inc_P_series(const double a, const double x)
{
	const int nmax = 5000;
	double dd = log_incgamma_dom(a, x);
	double sum  = 1.0;
	double term = 1.0;
	int n;
	for(n=1; n<nmax; n++) {
		term *= x/(a+n);
		sum  += term;
		if(fabs(term/sum) < ExCo<double>::epsilon()) break;
	}
	return dd + log(sum);
}

double log_gamma_inc_F_CF(const double a, const double x) {
	const int    nmax  =  5000;
	const double smallr =  pow(ExCo<double>::epsilon(), 3.0f);

	double hn = 0.0;           /* convergent */
	double Cn = 1.0 / smallr;
	double Dn = 1.0;
	int n;

	/* n == 1 has a_1, b_1, b_0 independent of a,x,
	 so that has been done by hand                */
	for ( n = 2 ; n < nmax ; n++ ) {
		double an;
		double delta;

		if (n & 1)
			an = 0.5*(n-1)/x;
		else
			an = (0.5*n-a)/x;

		Dn = 1.0 + an * Dn;
		if ( fabs(Dn) < smallr )
			Dn = smallr;
		Cn = 1.0 + an/Cn;
		if ( fabs(Cn) < smallr )
			Cn = smallr;
		Dn = 1.0 / Dn;
		delta = Cn * Dn;
		hn += log(fabs(delta));
		if(fabs(delta-1.0) < ExCo<double>::epsilon()) break;
	}
	return hn;
}




/* Q large x asymptotic, LFH log transformed */
double
log_gamma_inc_Q_large_x(const double a, const double x)
{
	const int nmax = 5000;
	double logD = log_incgamma_dom(a,x);
	double sum  = 1.0;
	double term = 1.0;
	double last = 1.0;
	int n;
	for(n=1; n<nmax; n++) {
		term *= (a-n)/x;
		if(fabs(term/last) > 1.0) break;
		if(fabs(term/sum)  < ExCo<double>::epsilon()) break;
		sum  += term;
		last  = term;
	}
	return logD + log(a) - log(x) + log(sum);
}

/* LFH log transformed */
double log_gamma_inc_Q_CF(const double a, const double x){
	double logD = log_incgamma_dom(a,x);
	double logF = log_gamma_inc_F_CF(a, x);
	return logD + log(a) - log(x) + logF;
}
/* LFH log transformed */

// computes log(1+x)
double gslfh_sf_log_1plusx(const double x){
	if (fabs(x) >= 0.5f) return log(1.0 + x);
	else if(fabs(x) >= GSL_ROOT6_DBL_EPSILON) {
		double t = 0.5*(8.0*x + 1.0)/(x+2.0);
		gsl_sf_result c;
		cheb_eval_e(&lopx_cs, t, &c);
		return  x * c.val;
	}else if(fabs(x) > -1.0f) {
		const double c1 = -0.5;
		const double c2 =  1.0/3.0;
		const double c3 = -1.0/4.0;
		const double c4 =  1.0/5.0;
		const double c5 = -1.0/6.0;
		const double c6 =  1.0/7.0;
		const double c7 = -1.0/8.0;
		const double c8 =  1.0/9.0;
		const double c9 = -1.0/10.0;
		const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
		return x * (1.0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*t)))));
	}else return GSL_NEGINF;
}

double incgamma_log_frac(const double a, const double x) {

	if(a < 0.0 || x < 0.0) myexit("invalid Domain for incgamma_log_frac\n");
	else if(x == 0.0) return(0.0);
	else if(a == 0.0) return(GSL_NEGINF);
	else if(x <= 0.5*a) {
		// If the series is quick, do that. It is robust and simple.
		return gslfh_sf_log_1plusx(-gamma_inc_P_series(a, x));
	}
	else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
		/* Then try the difficult asymptotic regime.
		 * This is the only way to do this region.
		 */
		return log(gamma_inc_Q_asymp_unif(a, x));
	}
	else if(a < 0.2 && x < 5.0) {
		/* Cancellations at small a must be handled
		 * analytically; x should not be too big
		 * either since the series terms grow
		 * with x and log(x).
		 */
		return log(gamma_inc_Q_series(a, x));
	}
	else if(a <= x) {
		if(x <= 1.0e+06) {
			/* Continued fraction is excellent for x >~ a.
			 * We do not let x be too large when x > a since
			 * it is somewhat pointless to try this there;
			 * the function is rapidly decreasing for
			 * x large and x > a, and it will just
			 * underflow in that region anyway. We
			 * catch that case in the standard
			 * large-x method.
			 */
			return log_gamma_inc_Q_CF(a, x);
		}else{
			return log_gamma_inc_Q_large_x(a, x);
		}
	}
	else {
		if(x > a - sqrt(a)) {
			/* Continued fraction again. The convergence
			 * is a little slower here, but that is fine.
			 * We have to trade that off against the slow
			 * convergence of the series, which is the
			 * only other option.
			 */
			return log_gamma_inc_Q_CF(a, x);
		}
		else {
			return gslfh_sf_log_1plusx(-gamma_inc_P_series(a, x));
		}
	}
}

/* LFH log transformed */
double incgamma_1minus_frac(const double a, const double x) {
	if(a < 0.0 || x < 0.0) myexit("invalid Domain for incgamma_log_1minus_frac\n");
	else if(x == 0.0) return(0.0);
	else if(a == 0.0) return(1.0);
	else if(x <= 0.5*a) {
		// If the series is quick, do that. It is robust and simple.
		return gamma_inc_P_series(a, x);
	}
	else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
		/* Then try the difficult asymptotic regime.
		 * This is the only way to do this region.
		 */
		return 1.0f - gamma_inc_Q_asymp_unif(a, x);
	}
	else if(a < 0.2 && x < 5.0) {
		/* Cancellations at small a must be handled
		 * analytically; x should not be too big
		 * either since the series terms grow
		 * with x and log(x).
		 */
		return 1.0f - gamma_inc_Q_series(a, x);
	}
	else if(a <= x) {
		if(x <= 1.0e+06) {
			/* Continued fraction is excellent for x >~ a.
			 * We do not let x be too large when x > a since
			 * it is somewhat pointless to try this there;
			 * the function is rapidly decreasing for
			 * x large and x > a, and it will just
			 * underflow in that region anyway. We
			 * catch that case in the standard
			 * large-x method.
			 */
			return 1.0f - gamma_inc_Q_CF(a, x);
		}
		else {
			return 1.0f - gamma_inc_Q_large_x(a, x);
		}
	}
	else {
		if(x > a - sqrt(a)) {
			/* Continued fraction again. The convergence
			 * is a little slower here, but that is fine.
			 * We have to trade that off against the slow
			 * convergence of the series, which is the
			 * only other option.
			 */
			return 1.0f - gamma_inc_Q_CF(a, x);
		}
		else {
			return gamma_inc_P_series(a, x);
		}
	}
}

double incgamma_log_1minus_frac(const double a, const double x) {
	if(a < 0.0 || x < 0.0) myexit("invalid Domain for incgamma_log_1minus_frac\n");
	else if(x == 0.0) return(GSL_NEGINF);
	else if(a == 0.0) return(0.0);
	else if(x <= 0.5*a) {
		// If the series is quick, do that. It is robust and simple.
		return log_gamma_inc_P_series(a, x);
	}
	else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
		/* Then try the difficult asymptotic regime.
		 * This is the only way to do this region.
		 */
		return gslfh_sf_log_1plusx(-gamma_inc_Q_asymp_unif(a, x));
	}
	else if(a < 0.2 && x < 5.0) {
		/* Cancellations at small a must be handled
		 * analytically; x should not be too big
		 * either since the series terms grow
		 * with x and log(x).
		 */
		return gslfh_sf_log_1plusx(-gamma_inc_Q_series(a, x));
	}
	else if(a <= x) {
		if(x <= 1.0e+06) {
			/* Continued fraction is excellent for x >~ a.
			 * We do not let x be too large when x > a since
			 * it is somewhat pointless to try this there;
			 * the function is rapidly decreasing for
			 * x large and x > a, and it will just
			 * underflow in that region anyway. We
			 * catch that case in the standard
			 * large-x method.
			 */
			return gslfh_sf_log_1plusx(-gamma_inc_Q_CF(a, x));
		}
		else {
			return gslfh_sf_log_1plusx(-gamma_inc_Q_large_x(a, x));
		}
	}
	else {
		if(x > a - sqrt(a)) {
			/* Continued fraction again. The convergence
			 * is a little slower here, but that is fine.
			 * We have to trade that off against the slow
			 * convergence of the series, which is the
			 * only other option.
			 */
			return gslfh_sf_log_1plusx(-gamma_inc_Q_CF(a, x));
		}
		else {
			return log_gamma_inc_P_series(a, x);
		}
	}
}


 int gsl_sf_beta_e(const double x, const double y, gsl_sf_result * result) {
	if((x > 0 && y > 0) && x < 50.0 && y < 50.0) {
		// Handle the easy case
		gsl_sf_result gx, gy, gxy;
		gsl_sf_gamma_e(x, &gx);
		gsl_sf_gamma_e(y, &gy);
		gsl_sf_gamma_e(x+y, &gxy);
		result->val  = (gx.val*gy.val)/gxy.val;
		result->err  = gx.err * fabs(gy.val/gxy.val);
		result->err += gy.err * fabs(gx.val/gxy.val);
		result->err += fabs((gx.val*gy.val)/(gxy.val*gxy.val)) * gxy.err;
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		return 0;
	}
	else if (isnegint(x) || isnegint(y)) {
		DOMAIN_ERROR(result);
	} else if (isnegint(x+y)) {  // infinity in the denominator
		result->val = 0.0;
		result->err = 0.0;
		return 0;
	} else {
		gsl_sf_result lb;
		double sgn;
		int stat_lb = gsl_sf_lnbeta_sgn_e(x, y, &lb, &sgn);
		if(stat_lb == 0) {
			int status = gsl_sf_exp_err_e(lb.val, lb.err, result);
			result->val *= sgn;
			return status;
		}
		else {
			result->val = 0.0;
			result->err = 0.0;
			return stat_lb;
		}
	}
}

double gsl_sf_beta(const double x, const double y)
{ gsl_sf_result result;
  gsl_sf_beta_e(x, y, &result);
  return result.val;
}
/*

 double Bessel0(double x)  {
 double m1, m2, n1, n2;
 int    n;

 //Labels: e100, e200
 double a,a1,a2,b,c,e,e1,e2,e3,e4,x1;
 int    m;
 //Calculate P and Q polynomials; P0(x)=m1, P1(x)=m2,
 //Q0(x)=n1, Q1(x)=n2
 a = 1; a1 = 1; a2 = 1; b = 1; c = 1;
 e1 = 1e6;
 m = -1;
 x1 = 1.0 / (8.0 * x);
 x1 = x1 * x1;
 m1 = 1.0; m2 = 1.0;
 n1 = -1.0 / (8.0 * x);
 n2 = -3.0 * n1;
 n = 0;
 while(true){
 m = m + 2; a = a * m * m;
 m = m + 2; a = a * m * m;
 c = c * x1;
 a1 = a1 * a2; a2 = a2 + 1.0; a1 = a1 * a2;
 a2 = a2 + 1.0; e2 = a * c / a1;
 e4 = 1.0 + (m + 2) / m + (m + 2) * (m + 2) / (a2 * 8 * x) + (m + 2) * (m + 4) / (a2 * 8 * x);
 e4 = e4 * e2;
 //Test for divergence
 if (fabs(e4) > e1) break;
 e1 = fabs(e2);
 m1 = m1 - e2;
 m2 = m2 + e2 * (m + 2) / m;
 n1 = n1 + e2 * (m + 2) * (m + 2) / (a2 * 8 * x);
 n2 = n2 - e2 * (m + 2) * (m + 4) / (a2 * 8 * x);
 n++;
 //Test for convergence criterion
 if (e1 < pow(2.0f,-22.0f)) break;
 }
 a = M_PI;
 e = e2;
 b = sqrt(2.0 / (a * x));
 return(b * (m1 * cos(x - a / 4) - n1 * sin(x - a / 4.0)));
 }

 double Bessel1(double x)  {
 double m1, m2, n1, n2;
 int    n;

 //Labels: e100, e200
 double a,a1,a2,b,c,e,e1,e2,e3,e4,x1;
 int    m;
 //Calculate P and Q polynomials; P0(x)=m1, P1(x)=m2,
 //Q0(x)=n1, Q1(x)=n2
 a = 1; a1 = 1; a2 = 1; b = 1; c = 1;
 e1 = 1e6;
 m = -1;
 x1 = 1.0 / (8.0 * x);
 x1 = x1 * x1;
 m1 = 1.0; m2 = 1.0;
 n1 = -1.0 / (8.0 * x);
 n2 = -3.0 * n1;
 n = 0;
 while(true){
 m = m + 2; a = a * m * m;
 m = m + 2; a = a * m * m;
 c = c * x1;
 a1 = a1 * a2; a2 = a2 + 1.0; a1 = a1 * a2;
 a2 = a2 + 1.0; e2 = a * c / a1;
 e4 = 1.0 + (m + 2) / m + (m + 2) * (m + 2) / (a2 * 8 * x) + (m + 2) * (m + 4) / (a2 * 8 * x);
 e4 = e4 * e2;
 //Test for divergence
 if (fabs(e4) > e1) break;
 e1 = fabs(e2);
 m1 = m1 - e2;
 m2 = m2 + e2 * (m + 2) / m;
 n1 = n1 + e2 * (m + 2) * (m + 2) / (a2 * 8 * x);
 n2 = n2 - e2 * (m + 2) * (m + 4) / (a2 * 8 * x);
 n++;
 //Test for convergence criterion
 if (e1 < e3) break;
 }
 e200: a = M_PI;
 e = e2;
 b = sqrt(2.0 / (a * x));
 return(b * (m2 * cos(x - 3 * a / 4.0) - n2 * sin(x - 3 * a / 4.0)));
 }

 */

 int
gsl_sf_bessel_K_scaled_steed_temme_CF2(const double nu, const double x,
                                       double * K_nu, double * K_nup1,
                                       double * Kp_nu)
{
  const int maxiter = 10000;

  int i = 1;
  double bi = 2.0*(1.0 + x);
  double di = 1.0/bi;
  double delhi = di;
  double hi    = di;

  double qi   = 0.0;
  double qip1 = 1.0;

  double ai = -(0.25 - nu*nu);
  double a1 = ai;
  double ci = -ai;
  double Qi = -ai;

  double s = 1.0 + Qi*delhi;

  for(i=2; i<=maxiter; i++) {
    double dels;
    double tmp;
    ai -= 2.0*(i-1);
    ci  = -ai*ci/i;
    tmp  = (qi - bi*qip1)/ai;
    qi   = qip1;
    qip1 = tmp;
    Qi += ci*qip1;
    bi += 2.0;
    di  = 1.0/(bi + ai*di);
    delhi = (bi*di - 1.0) * delhi;
    hi += delhi;
    dels = Qi*delhi;
    s += dels;
    if(fabs(dels/s) < GSL_DBL_EPSILON) break;
  }

  hi *= -a1;

  *K_nu   = sqrt(M_PI/(2.0*x)) / s;
  *K_nup1 = *K_nu * (nu + x + 0.5 - hi)/x;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;
  if(i == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}

static double g1_dat[14] = {
  -1.14516408366268311786898152867,
   0.00636085311347084238122955495,
   0.00186245193007206848934643657,
   0.000152833085873453507081227824,
   0.000017017464011802038795324732,
  -6.4597502923347254354668326451e-07,
  -5.1819848432519380894104312968e-08,
   4.5189092894858183051123180797e-10,
   3.2433227371020873043666259180e-11,
   6.8309434024947522875432400828e-13,
   2.8353502755172101513119628130e-14,
  -7.9883905769323592875638087541e-16,
  -3.3726677300771949833341213457e-17,
  -3.6586334809210520744054437104e-20
};
static cheb_series g1_cs = {
  g1_dat,
  13,
  -1, 1,
  7
};

/* nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu]) */
static double g2_dat[15] =
{
  1.882645524949671835019616975350,
 -0.077490658396167518329547945212,
 -0.018256714847324929419579340950,
  0.0006338030209074895795923971731,
  0.0000762290543508729021194461175,
 -9.5501647561720443519853993526e-07,
 -8.8927268107886351912431512955e-08,
 -1.9521334772319613740511880132e-09,
 -9.4003052735885162111769579771e-11,
  4.6875133849532393179290879101e-12,
  2.2658535746925759582447545145e-13,
 -1.1725509698488015111878735251e-15,
 -7.0441338200245222530843155877e-17,
 -2.4377878310107693650659740228e-18,
 -7.5225243218253901727164675011e-20
};
static cheb_series g2_cs = {
  g2_dat,
  14,
  -1, 1,
  8
};
static
int
gsl_sf_temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2)
{
  const double anu = fabs(nu);    /* functions are even */
  const double x = 4.0*anu - 1.0;
  gsl_sf_result r_g1;
  gsl_sf_result r_g2;
  cheb_eval_e(&g1_cs, x, &r_g1);
  cheb_eval_e(&g2_cs, x, &r_g2);
  *g1 = r_g1.val;
  *g2 = r_g2.val;
  *g_1mnu = 1.0/(r_g2.val + nu * r_g1.val);
  *g_1pnu = 1.0/(r_g2.val - nu * r_g1.val);
  return GSL_SUCCESS;
}


int
gsl_sf_bessel_K_scaled_temme(const double nu, const double x,
                             double * K_nu, double * K_nup1, double * Kp_nu)
{
  const int max_iter = 15000;

  const double half_x    = 0.5 * x;
  const double ln_half_x = log(half_x);
  const double half_x_nu = exp(nu*ln_half_x);
  const double pi_nu   = M_PI * nu;
  const double sigma   = -nu * ln_half_x;
  const double sinrat  = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));
  const double sinhrat = (fabs(sigma) < GSL_DBL_EPSILON ? 1.0 : sinh(sigma)/sigma);
  const double ex = exp(x);

  double sum0, sum1;
  double fk, pk, qk, hk, ck;
  int k = 0;
  int stat_iter;

  double g_1pnu, g_1mnu, g1, g2;
  int stat_g = gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);

  fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = 0.5/half_x_nu * g_1pnu;
  qk = 0.5*half_x_nu * g_1mnu;
  hk = pk;
  ck = 1.0;
  sum0 = fk;
  sum1 = hk;
  while(k < max_iter) {
    double del0;
    double del1;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    hk  = -k*fk + pk;
    del0 = ck * fk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if(fabs(del0) < 0.5*fabs(sum0)*GSL_DBL_EPSILON) break;
  }

  *K_nu   = sum0 * ex;
  *K_nup1 = sum1 * 2.0/x * ex;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;

  stat_iter = ( k == max_iter ? GSL_EMAXITER : GSL_SUCCESS );
  return GSL_ERROR_SELECT_2(stat_iter, stat_g);
}

double gsl_sf_bessel_Knu_scaled_e10_e(double nu, double x){
  /* CHECK_POINTER(result) */
  x =fabs(x);
  nu =fabs(nu);
  // if (nu == 0) TODO!

    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */
    double K_mu, K_mup1, Kp_mu;
    double K_nu, K_nup1, K_num1;
    int n, e10 = 0;

    if(x < 2.0) {
      gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    } else {
      gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }

    /* recurse forward to obtain K_num1, K_nu */
    K_nu   = K_mu;
    K_nup1 = K_mup1;

    for(n=0; n<N; n++) {
      K_num1 = K_nu;
      K_nu   = K_nup1;
      /* rescale the recurrence to avoid overflow */
      if (fabs(K_nu) > GSL_SQRT_DBL_MAX) {
        double p = floor(log(fabs(K_nu))/M_LN10);
        double factor = pow(10.0, p);
        K_num1 /= factor;
        K_nu /= factor;
        e10 += p;
      };
      K_nup1 = 2.0*(mu+n+1)/x * K_nu + K_num1;
    }
    return K_nu;
}

 // assumes n >= 1.0f!
 /*
static
double bessel_Kn_small_x(const double n, const double x) {
  int k;
  double y = 0.25 * x * x;
  double ln_x_2 = log(0.5*fabs(x));
  double k_term;
  double term1, sum1, ln_pre1;
  double term2, sum2, pre2;

  double npk_fact = lngamma(n);
  ln_pre1 = -n*ln_x_2 + npk_fact;
  npk_fact = exp(npk_fact);
  sum1 = 1.0;
  k_term = 1.0;
  for(k=1; k<=n-1; k++) {
    k_term *= -y/(k * (n-k));
    sum1 += k_term;
  }
  term1 = 0.5 * exp(ln_pre1) * sum1;

  pre2 = 0.5 * exp(n*ln_x_2);
  if(pre2 > 0.0) {
    const int KMAX = 20;

    double yk = 1.0;
    double k_fact  = 1.0;
    double psi_kp1 = -M_EULER;
    double psi_npkp1;
    psi_npkp1 = d_lngamma_dx(n) + 1.0/n;
    sum2 = (psi_kp1 + psi_npkp1 - 2.0*ln_x_2)/npk_fact;
    for(k=1; k<KMAX; k++) {
      psi_kp1   += 1.0/k;
      psi_npkp1 += 1.0/(n+k);
      k_fact    *= k;
      npk_fact *= n+k;
      yk *= y;
      k_term = yk*(psi_kp1 + psi_npkp1 - 2.0*ln_x_2)/(k_fact*npk_fact);
      sum2 += k_term;
    }
    term2 = ( GSL_IS_ODD(n) ? -1.0 : 1.0 ) * pre2 * sum2;
  }
  else {
    term2 = 0.0;
  }

  return term1 + term2;
}*/


double BesselJ0 (double X) {
	/***********************************************************************
	 This subroutine calculates the First Kind Bessel Function of
	 order 0, for any real number X. The polynomial approximation by
	 series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	 REFERENCES:
	 M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	 C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	 VOL.5, 1962.
	 ************************************************************************/
	const double
	P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
	P4=-0.2073370639E-5, P5= 0.2093887211E-6,
	Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
	Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
	R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
	R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
	S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
	S4= 59272.64853, S5=267.8532712, S6=1.0;
	double
	AX,FR,FS,Z,FP,FQ,XX,Y, TMP;

	if (X==0.0) return 1.0;
	AX = fabs(X);
	if (AX < 8.0) {
		Y = X*X;
		FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
		FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
		TMP = FR/FS;
	}
	else {
		Z = 8./AX;
		Y = Z*Z;
		XX = AX-0.785398164;
		FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
		FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
		TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
	}
	return TMP;
}

double Sign(double X, double Y) {
	if (Y<0.0) return (-fabs(X));
	else return (fabs(X));
}

double BesselJ1 (double X) {
	/**********************************************************************
	 This subroutine calculates the First Kind Bessel Function of
	 order 1, for any real number X. The polynomial approximation by
	 series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	 REFERENCES:
	 M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	 C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	 VOL.5, 1962.
	 ***********************************************************************/
	const double
	P1=1.0, P2=0.183105E-2, P3=-0.3516396496E-4, P4=0.2457520174E-5,
	P5=-0.240337019E-6,  P6=0.636619772,
	Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
	Q4=-0.88228987E-6, Q5= 0.105787412E-6,
	R1= 72362614232.0, R2=-7895059235.0, R3=242396853.1,
	R4=-2972611.439,   R5=15704.48260,  R6=-30.16036606,
	S1=144725228442.0, S2=2300535178.0, S3=18583304.74,
	S4=99447.43394,    S5=376.9991397,  S6=1.0;

	double AX,FR,FS,Y,Z,FP,FQ,XX, TMP;

	AX = fabs(X);
	if (AX < 8.0) {
		Y = X*X;
		FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
		FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
		TMP = X*(FR/FS);
	}
	else {
		Z = 8.0/AX;
		Y = Z*Z;
		XX = AX-2.35619491;
		FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
		FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
		TMP = sqrt(P6/AX)*(cos(XX)*FP-Z*sin(XX)*FQ)*Sign(S6,X);
	}
	return TMP;
}

double BesselJN (int N, double X) {
	/************************************************************************
	 This subroutine calculates the first kind modified Bessel function
	 of integer order N, for any REAL X. We use here the classical
	 recursion formula, when X > N. For X < N, the Miller's algorithm
	 is used to avoid overflows.
	 -----------------------------
	 REFERENCE:
	 C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
	 MATHEMATICAL TABLES, VOL.5, 1962.
	 *************************************************************************/
	const int IACC = 40;
	const double BIGNO = 1e10,  BIGNI = 1e-10;

	double TOX,BJM,BJ,BJP,SUM,TMP;
	int J, JSUM, M;

	if (N == 0) return BesselJ0(X);
	if (N == 1) return BesselJ1(X);
	if (X == 0.0) return 0.0;

	TOX = 2.0/X;
	if (X > 1.0*N) {
		BJM = BesselJ0(X);
		BJ  = BesselJ1(X);
		for (J=1; J<N; J++) {
			BJP = J*TOX*BJ-BJM;
			BJM = BJ;
			BJ  = BJP;
		}
		return BJ;
	}
	else {
		M = (int) (2*((N+floor(sqrt(1.0*(IACC*N))))/2));
		TMP = 0.0;
		JSUM = 0;
		SUM = 0.0;
		BJP = 0.0;
		BJ  = 1.0;
		for (J=M; J>0; J--) {
			BJM = J*TOX*BJ-BJP;
			BJP = BJ;
			BJ  = BJM;
			if (fabs(BJ) > BIGNO) {
				BJ  = BJ*BIGNI;
				BJP = BJP*BIGNI;
				TMP = TMP*BIGNI;
				SUM = SUM*BIGNI;
			}
			if (JSUM != 0)  SUM += BJ;
			JSUM = 1-JSUM;
			if (J == N)  TMP = BJP;
		}
		SUM = 2.0*SUM-BJ;
		return (TMP/SUM);
	}
}

double BesselI0(double X) {
	double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
	P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067429;
	P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
	Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
	Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
	Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
	if (fabs(X) < 3.75) {
		Y=(X/3.75)*(X/3.75);
		return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
	}
	else {
		AX=fabs(X);
		Y=3.75/AX;
		BX=exp(AX)/sqrt(AX);
		AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
		return (AX*BX);
	}
}

double BesselI1(double X) {
	double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
	P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
	P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
	Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
	Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
	Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
	if (fabs(X) < 3.75) {
		Y=(X/3.75)*(X/3.75);
		return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
	}
	else {
		AX=fabs(X);
		Y=3.75/AX;
		BX=exp(AX)/sqrt(AX);
		AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
		return (AX*BX);
	}
}

double BesselIN(int N, double X) {
	/*
	 !     This subroutine calculates the first kind modified Bessel function
	 !     of integer order N, for any REAL X. We use here the classical
	 !     recursion formula, when X > N. For X < N, the Miller's algorithm
	 !     is used to avoid overflows.
	 !     REFERENCE:
	 !     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
	 !     MATHEMATICAL TABLES, VOL.5, 1962.
	 */

	int IACC = 40;
	double BIGNO = 1e10, BIGNI = 1e-10;
	double TOX, BIM, BI, BIP, BSI;
	int J, M;

	if (N==0)  return (BesselI0(X));
	if (N==1)  return (BesselI1(X));
	if (X==0.0) return 0.0;

	TOX = 2.0/X;
	BIP = 0.0;
	BI  = 1.0;
	BSI = 0.0;
	M = (int) (2*((N+floor(sqrt(IACC*N)))));
	for (J = M; J>0; J--) {
		BIM = BIP+J*TOX*BI;
		BIP = BI;
		BI  = BIM;
		if (fabs(BI) > BIGNO) {
			BI  = BI*BIGNI;
			BIP = BIP*BIGNI;
			BSI = BSI*BIGNI;
		}
		if (J==N)  BSI = BIP;
	}
	return (BSI*BesselI0(X)/BI);
}

/* cdf/inverse_normal.c
 * cdf/gauss.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

 /* cdf/inverse_normal.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

static double sliwa_rat_eval (const double a[], const size_t na, const double b[], const size_t nb, const double x){
  size_t i, j;
  double u, v, r;
  u = a[na - 1];
  for (i = na - 1; i > 0; i--) u = x * u + a[i - 1];
  v = b[nb - 1];
  for (j = nb - 1; j > 0; j--) v = x * v + b[j - 1];
  r = u / v;
  return r;
}

static double sliwa_small (double q){
  const double a[8] = { 3.387132872796366608, 133.14166789178437745,
    1971.5909503065514427, 13731.693765509461125,
    45921.953931549871457, 67265.770927008700853,
    33430.575583588128105, 2509.0809287301226727
  };
  const double b[8] = { 1.0, 42.313330701600911252,
    687.1870074920579083, 5394.1960214247511077,
    21213.794301586595867, 39307.89580009271061,
    28729.085735721942674, 5226.495278852854561
  };
  double r = 0.180625 - q * q;
  return q * sliwa_rat_eval (a, 8, b, 8, r);
}

static double sliwa_intermediate (double r){
  const double a[] = { 1.42343711074968357734, 4.6303378461565452959,
    5.7694972214606914055, 3.64784832476320460504,
    1.27045825245236838258, 0.24178072517745061177,
    0.0227238449892691845833, 7.7454501427834140764e-4
  };
  const double b[] = { 1.0, 2.05319162663775882187,
    1.6763848301838038494, 0.68976733498510000455,
    0.14810397642748007459, 0.0151986665636164571966,
    5.475938084995344946e-4, 1.05075007164441684324e-9
  };
  return sliwa_rat_eval (a, 8, b, 8, (r - 1.6));
}

static double sliwa_tail (double r) {
  const double a[] = { 6.6579046435011037772, 5.4637849111641143699,
    1.7848265399172913358, 0.29656057182850489123,
    0.026532189526576123093, 0.0012426609473880784386,
    2.71155556874348757815e-5, 2.01033439929228813265e-7
  };
  const double b[] = { 1.0, 0.59983220655588793769,
    0.13692988092273580531, 0.0148753612908506148525,
    7.868691311456132591e-4, 1.8463183175100546818e-5,
    1.4215117583164458887e-7, 2.04426310338993978564e-15
  };
  return sliwa_rat_eval (a, 8, b, 8, (r - 5.0));
}

double sliwa_gsl_cdf_ugaussian_Pinv (const double P){
  double r, x, pp;
  double dP = P - 0.5;
  if (P == 1.0) return GSL_POSINF;
  else if (P == 0.0) return GSL_NEGINF;
  if (fabs (dP) <= 0.425) return sliwa_small (dP);
  pp = (P < 0.5) ? P : 1.0 - P;
  r = sqrt (-log (pp));
  x = (r <= 5.0) ? sliwa_intermediate (r) : sliwa_tail (r);
  return (P < 0.5) ? -x : x;
}

double handfield_gsl_cdf_ugaussian_logPinv (const double P){
  double r, x;
  if (P >= 0.0) return GSL_POSINF;
  double dP = exp(P) - 0.5;
  if (fabs (dP) <= 0.425) return sliwa_small (dP);
  if (P < -0.7){
      r = sqrt(-P);
      return -((r <= 5.0) ? sliwa_intermediate (r) : sliwa_tail (r));
  }else{
      r = sqrt( -log((P < -0.002) ? (1 - exp(P)) :  -P * (1.0 + 0.5*P*(1.0 + P/3.0*(1.0 + 0.25*P*(1.0 + 0.2*P)))) ) );
      return   (r <= 5.0) ? sliwa_intermediate (r) : sliwa_tail (r);
  }
}



#define GAUSS_SCALE (16.0)
#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif
#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI (M_2_SQRTPI * M_SQRT1_2 / 2.0)
#endif
#define SQRT32 (4.0 * M_SQRT2)

static double
get_del (double x, double rational)
{
  double xsq = 0.0;
  double del = 0.0;
  double result = 0.0;

  xsq = floor (x * GAUSS_SCALE) / GAUSS_SCALE;
  del = (x - xsq) * (x + xsq);
  del *= 0.5;

  result = exp (-0.5 * xsq * xsq) * exp (-1.0 * del) * rational;

  return result;
}

static double
gauss_small (const double x)
{
  unsigned int i;
  double result = 0.0;
  double xsq;
  double xnum;
  double xden;

  const double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };

  xsq = x * x;
  xnum = a[4] * xsq;
  xden = xsq;

  for (i = 0; i < 3; i++)
    {
      xnum = (xnum + a[i]) * xsq;
      xden = (xden + b[i]) * xsq;
    }

  result = x * (xnum + a[3]) / (xden + b[3]);

  return result;
}

/*
 * Normal cdf for 0.66291 < fabs(x) < sqrt(32).
 */
static double
gauss_medium (const double x)
{
  unsigned int i;
  double temp = 0.0;
  double result = 0.0;
  double xnum;
  double xden;
  double absx;

  const double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };

  absx = fabs (x);

  xnum = c[8] * absx;
  xden = absx;

  for (i = 0; i < 7; i++)
    {
      xnum = (xnum + c[i]) * absx;
      xden = (xden + d[i]) * absx;
    }

  temp = (xnum + c[7]) / (xden + d[7]);

  result = get_del (x, temp);

  return result;
}
/*
 * Normal cdf for
 * {sqrt(32) < x < GAUSS_XUPPER} union { GAUSS_XLOWER < x < -sqrt(32) }.
 */
static double
gauss_large (const double x)
{
  int i;
  double result;
  double xsq;
  double temp;
  double xnum;
  double xden;
  double absx;

  const double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  absx = fabs (x);
  xsq = 1.0 / (x * x);
  xnum = p[5] * xsq;
  xden = xsq;

  for (i = 0; i < 4; i++)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }

  temp = xsq * (xnum + p[4]) / (xden + q[4]);
  temp = (M_1_SQRT2PI - temp) / absx;

  result = get_del (x, temp);

  return result;
}
double gsl_cdf_ugaussian_P (const double x){
  double result;
  double absx = fabs (x);
  if (absx < GSL_DBL_EPSILON / 2) return 0.5;
  else if (absx < 0.66291) return 0.5 + gauss_small (x);
  else if (absx < SQRT32){
      result = gauss_medium (x);
      return (x > 0.0) ? 1.0 - result : result;
  }
  else if (x > 8.572) return 1.0;
  else if (x < -37.519) return 0.0;
  else{
      result = gauss_large (x);

      if (x > 0.0)
        {
          result = 1.0 - result;
        }
    }
  return result;
}


/* cdf/tdist.c
 *
 * Copyright (C) 2002 Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the Student's t cumulative distribution function using
 * the method detailed in
 *
 * W.J. Kennedy and J.E. Gentle, "Statistical Computing." 1980.
 * Marcel Dekker. ISBN 0-8247-6898-1.
 *
 * G.W. Hill and A.W. Davis. "Generalized asymptotic expansions
 * of Cornish-Fisher type." Annals of Mathematical Statistics,
 * vol. 39, 1264-1273. 1968.
 *
 * G.W. Hill. "Algorithm 395: Student's t-distribution," Communications
 * of the ACM, volume 13, number 10, page 617. October 1970.
 *
 * G.W. Hill, "Remark on algorithm 395: Student's t-distribution,"
 * Transactions on Mathematical Software, volume 7, number 2, page
 * 247. June 1982.
 */

static double
poly_eval (const double c[], unsigned int n, double x)
{
  unsigned int i;
  double y = c[0] * x;

  for (i = 1; i < n; i++)
    {
      y = x * (y + c[i]);
    }

  y += c[n];

  return y;
}
static double
cornish_fisher (double t, double n)
{
  const double coeffs6[10] = {
    0.265974025974025974026,
    5.449696969696969696970,
    122.20294372294372294372,
    2354.7298701298701298701,
    37625.00902597402597403,
    486996.1392857142857143,
    4960870.65,
    37978595.55,
    201505390.875,
    622437908.625
  };
  const double coeffs5[8] = {
    0.2742857142857142857142,
    4.499047619047619047619,
    78.45142857142857142857,
    1118.710714285714285714,
    12387.6,
    101024.55,
    559494.0,
    1764959.625
  };
  const double coeffs4[6] = {
    0.3047619047619047619048,
    3.752380952380952380952,
    46.67142857142857142857,
    427.5,
    2587.5,
    8518.5
  };
  const double coeffs3[4] = {
    0.4,
    3.3,
    24.0,
    85.5
  };

  double a = n - 0.5;
  double b = 48.0 * a * a;

  double z2 = a * log1p (t * t / n);
  double z = sqrt (z2);

  double p5 = z * poly_eval (coeffs6, 9, z2);
  double p4 = -z * poly_eval (coeffs5, 7, z2);
  double p3 = z * poly_eval (coeffs4, 5, z2);
  double p2 = -z * poly_eval (coeffs3, 3, z2);
  double p1 = z * (z2 + 3.0);
  double p0 = z;

  double y = p5;
  y = (y / b) + p4;
  y = (y / b) + p3;
  y = (y / b) + p2;
  y = (y / b) + p1;
  y = (y / b) + p0;

  if (t < 0)
    y *= -1;

  return y;
}
double gsl_cdf_tdist_P (const double x, const double nu){
  double P;
  double x2 = x * x;
  if (nu > 30 && x2 < 10 * nu){
      double u = cornish_fisher (x, nu);
      P = gsl_cdf_ugaussian_P (u);
      return P;
    }
  if (x2 < nu){
      double u = x2 / nu;
      double eps = u / (1 + u);
      return beta_inc_AXPY ((x >= 0) ? 0.5 : -0.5, 0.5, 0.5, nu / 2.0, eps);
  }else{
      double v = nu / (x * x);
      double eps = v / (1 + v);
      if (x >= 0) return  beta_inc_AXPY (-0.5, 1.0, nu / 2.0, 0.5, eps);
      else return beta_inc_AXPY (0.5, 0.0, nu / 2.0, 0.5, eps);
  }
}

double gsl_cdf_ugaussian_Q (const double x){
  double result;
  double absx = fabs (x);
  if (absx < GSL_DBL_EPSILON / 2) return 0.5;
  else if (absx < 0.66291){
      result = gauss_small (x);
      return (x < 0.0) ? fabs(result) + 0.5 : 0.5 - result;
  }else if (absx < SQRT32){
      result = gauss_medium (x);
      return  (x < 0.0) ? 1.0 - result : result;
  }else if (x > 37.519) return 0.0;
  else if (x < -8.572) return 1.0;
  else{
      result = gauss_large (x);
      return (x < 0.0)?result = 1.0 - result : result;
  }
}

double gsl_cdf_tdist_Q (const double x, const double nu) {
  double Q;
  double x2 = x * x;
  if (nu > 30 && x2 < 10 * nu){
      double u = cornish_fisher (x, nu);
      return gsl_cdf_ugaussian_Q (u);
  }
  if (x2 < nu){
      double u = x2 / nu;
      double eps = u / (1 + u);
      return beta_inc_AXPY ( (x >= 0) ? -0.5 : 0.5 , 0.5, 0.5, nu / 2.0, eps);
    }else{
      double v = nu / (x * x);
      double eps = v / (1 + v);
        if (x >= 0)return beta_inc_AXPY (0.5, 0.0, nu / 2.0, 0.5, eps);
        else return beta_inc_AXPY (-0.5, 1.0, nu / 2.0, 0.5, eps);
    }
}


/* cdf/tdistinv.c
 *
 * Copyright (C) 2007, 2010 Brian Gough
 * Copyright (C) 2002 Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


static double
inv_cornish_fisher (double z, double nu)
{
  double a = 1 / (nu - 0.5);
  double b = 48.0 / (a * a);

  double cf1 = z * (3 + z * z);
  double cf2 = z * (945 + z * z * (360 + z * z * (63 + z * z * 4)));

  double y = z - cf1 / b + cf2 / (10 * b * b);

  double t = GSL_SIGN (z) * sqrt (nu * expm1 (a * y * y));

  return t;
}

double
gsl_ran_tdist_pdf (const double x, const double nu)
{
  double p;

  double lg1 = lngamma (nu / 2);
  double lg2 = lngamma ((nu + 1) / 2);

  p = ((exp (lg2 - lg1) / sqrt (M_PI * nu))
       * pow ((1 + x * x / nu), -(nu + 1) / 2));
  return p;
}



double
gsl_cdf_tdist_Pinv (const double P, const double nu)
{
  double x, ptail;

  if (P == 1.0)
    {
      return GSL_POSINF;
    }
  else if (P == 0.0)
    {
      return GSL_NEGINF;
    }

  if (nu == 1.0)
    {
      x = tan (M_PI * (P - 0.5));
      return x;
    }
  else if (nu == 2.0)
    {
      x = (2 * P - 1) / sqrt (2 * P * (1 - P));
      return x;
    }

  ptail = (P < 0.5) ? P : 1 - P;

  if (sqrt (M_PI * nu / 2) * ptail > pow (0.05, nu / 2))
    {
      double xg = sliwa_gsl_cdf_ugaussian_Pinv (P);
      x = inv_cornish_fisher (xg, nu);
    }
  else
    {
      /* Use an asymptotic expansion of the tail of integral */

      double beta = gsl_sf_beta (0.5, nu / 2);

      if (P < 0.5)
        {
          x = -sqrt (nu) * pow (beta * nu * P, -1.0 / nu);
        }
      else
        {
          x = sqrt (nu) * pow (beta * nu * (1 - P), -1.0 / nu);
        }

      /* Correct nu -> nu/(1+nu/x^2) in the leading term to account
         for higher order terms. This avoids overestimating x, which
         makes the iteration unstable due to the rapidly decreasing
         tails of the distribution. */

      x /= sqrt (1 + nu / (x * x));
    }

  {
    double dP, phi;
    unsigned int n = 0;

  start:
    dP = P - gsl_cdf_tdist_P (x, nu);
    phi = gsl_ran_tdist_pdf (x, nu);

    if (dP == 0.0 || n++ > 32)
      goto end;

    {
      double lambda = dP / phi;
      double step0 = lambda;
      double step1 = ((nu + 1) * x / (x * x + nu)) * (lambda * lambda / 4.0);

      double step = step0;

      if (fabs (step1) < fabs (step0))
        {
          step += step1;
        }

      if (P > 0.5 && x + step < 0)
        x /= 2;
      else if (P < 0.5 && x + step > 0)
        x /= 2;
      else
        x += step;

      if (fabs (step) > 1e-10 * fabs (x))
        goto start;
    }

  end:
    if (fabs(dP) > GSL_SQRT_DBL_EPSILON * P)
      {
        GSL_ERROR_VAL("inverse failed to converge", GSL_EFAILED, GSL_NAN);
      }

    return x;
  }
}


double handfield_gsl_cdf_tdist_logPinv(const double PP, const double nu){
  double x, ptail;
  if (PP == 0.0) return GSL_POSINF;
  else if (PP == GSL_NEGINF) return GSL_NEGINF;
  if (nu == 1.0) return tan (M_PI * (exp(PP) - 0.5));
  else if (nu == 2.0) return (2 * exp(PP) - 1) / sqrt (2 * exp(PP) * (1 - exp(PP)));
  double P = exp(PP);
  ptail = (P < 0.5) ? P : 1 - P;

  if (sqrt (M_PI * nu / 2) * ptail > pow (0.05, nu / 2)){
      double xg = sliwa_gsl_cdf_ugaussian_Pinv (P);
      x = inv_cornish_fisher (xg, nu);
  }else{
      /* Use an asymptotic expansion of the tail of integral */

      double beta = gsl_sf_beta (0.5, nu / 2);

      if (P < 0.5) x = -sqrt (nu) * pow (beta * nu * P, -1.0 / nu);
      else x = sqrt (nu) * pow (beta * nu * (1 - P), -1.0 / nu);

      /* Correct nu -> nu/(1+nu/x^2) in the leading term to account
         for higher order terms. This avoids overestimating x, which
         makes the iteration unstable due to the rapidly decreasing
         tails of the distribution. */

      x /= sqrt (1 + nu / (x * x));
    }

  {
    double dP, phi;
    unsigned int n = 0;

  start:
    dP = P - gsl_cdf_tdist_P (x, nu);
    phi = gsl_ran_tdist_pdf (x, nu);

    if (dP == 0.0 || n++ > 32)
      goto end;

    {
      double lambda = dP / phi;
      double step0 = lambda;
      double step1 = ((nu + 1) * x / (x * x + nu)) * (lambda * lambda / 4.0);
      double step = step0;

      if (fabs (step1) < fabs (step0)) step += step1;
      if ((P > 0.5 && x + step < 0)||(P < 0.5 && x + step > 0)) x /= 2;
      else x += step;

      if (fabs (step) > 1e-10 * fabs (x))
        goto start;
    }

  end:
    if (fabs(dP) > GSL_SQRT_DBL_EPSILON * P)
      {
        GSL_ERROR_VAL("inverse failed to converge", GSL_EFAILED, GSL_NAN);
      }

    return x;
  }
}

double
gsl_cdf_tdist_Qinv (const double Q, const double nu)
{
  double x, qtail;

  if (Q == 0.0)return GSL_POSINF;
  else if (Q == 1.0) return GSL_NEGINF;
  else if (nu == 1.0) return  tan (M_PI * (0.5 - Q));
  else if (nu == 2.0) return (1 - 2 * Q) / sqrt (2 * Q * (1 - Q));

  qtail = (Q < 0.5) ? Q : 1 - Q;

  if (sqrt (M_PI * nu / 2) * qtail > pow (0.05, nu / 2))
    {
      double xg = sliwa_gsl_cdf_ugaussian_Pinv (-Q);
      x = inv_cornish_fisher (xg, nu);
    }
  else
    {
      /* Use an asymptotic expansion of the tail of integral */

      double beta = gsl_sf_beta (0.5, nu / 2);

      if (Q < 0.5) x = sqrt (nu) * pow (beta * nu * Q, -1.0 / nu);
      else x = -sqrt (nu) * pow (beta * nu * (1 - Q), -1.0 / nu);

      /* Correct nu -> nu/(1+nu/x^2) in the leading term to account
         for higher order terms. This avoids overestimating x, which
         makes the iteration unstable due to the rapidly decreasing
         tails of the distribution. */

      x /= sqrt (1 + nu / (x * x));
    }

  {
    double dQ, phi;
    unsigned int n = 0;

  start:
    dQ = Q - gsl_cdf_tdist_Q (x, nu);
    phi = gsl_ran_tdist_pdf (x, nu);

    if (dQ == 0.0 || n++ > 32) goto end;

    {
      double lambda = - dQ / phi;
      double step0 = lambda;
      double step1 = ((nu + 1) * x / (x * x + nu)) * (lambda * lambda / 4.0);

      double step = step0;

      if (fabs (step1) < fabs (step0)) step += step1;

      if (Q < 0.5 && x + step < 0) x /= 2;
      else if (Q > 0.5 && x + step > 0) x /= 2;
      else x += step;

      if (fabs (step) > 1e-10 * fabs (x))
        goto start;
    }
  }

end:

  return x;
}
