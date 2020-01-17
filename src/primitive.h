/*
 * primitive.h
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

// optionnal external dependency, which uses a modified file "GSLfunc.hpp"
// MAF,



#ifndef _defined_Primitive
#define _defined_Primitive

#define GNU_SCIENTIFIC_LIBRARY
#ifndef LFH_HAS_RUNNINGTIME_STATISTICS
#define LFH_HAS_RUNNINGTIME_STATISTICS true
#endif

#include "bastructs.h"

namespace LFHPrimitive{
	void LFHBreakpoint();
template <class A,class B> class AbstractKeyIterator{
    public:
    A curkey;
    const B* target;
    AbstractKeyIterator(){}
    AbstractKeyIterator(const B& _target) : target(&_target){}
    operator const A& () const{return(curkey);}
    const A& operator()()const{return(curkey);}

    virtual bool first()=0;
    virtual bool last()=0;
    virtual bool next()=0;
    virtual bool prev()=0;
};
template <class A,class B> class AbstractIterator{
public:
    AbstractIterator(){}
    A* first(B & ob){init(ob); return(next(ob));}
    virtual void init( B &)=0;
    virtual A* next( B &)=0;
};
enum functionname{
    FUNCNAME_NULL =0,
    FUNCNAME_SQUARRE,
    FUNCNAME_INVERSE,
    FUNCNAME_SIN,
    FUNCNAME_COS,
    FUNCNAME_LNGAMMA,
    FUNCNAME_POLYGAMMA_0,
    FUNCNAME_SINH,
    FUNCNAME_COSH,
    FUNCNAME_BESSEL_J0,
    FUNCNAME_BESSEL_J1,
    FUNCNAME_BESSEL_I0,
    FUNCNAME_BESSEL_I1
};



// for ordering, partial ordering is used, where either min and max needs to be strictly larger to for "<=" "=>" operator to be false, "a>b" is now !(a<=b)
// hence, equality "==" implies that both min and max are equal, but SETCMP_IS_NOT_SUBSET, SETCMP_IS_NOT_SUPERSET can be set

// x-y-Y-X  16 128    ==> ==
// x-y-X-Y  32 128    ==> >
// x-X-y-Y  32 128 1  ==> >
// y-x-X-Y  32 64     ==> ==
// y-x-Y-X  16 64     ==> <
// y-Y-x-X  16 64 1   ==> <

// y-x-X-Y  =>  "x==y"
// z-x-Z-X  ==> "z<x"

// is "z<y" guarantied? nope

// y-z-Z-Y  ==> "z == y"
// z-y-Z-Y  ==> "z<y"

template<class A, unsigned int value = ((sizeof(A) <= 4) ? 0 : 0) > // enum are 32 bit maximum... lame way to check
class ExCo_Enumflag{    public:    static const int flag = 0;};
template<class A> class ExCo_Enumflag<A, 1>{public: static const int flag = 1;};
// list of classes that are definitly not enums:
template<class A> class ExCo_Enumflag<vector<A>, 1>{public: static const int flag = 0;};
template<class A> class ExCo_Enumflag<A*, 1>{public: static const int flag = 0;};
template<class A, unsigned int SIZE> class ExCo_Enumflag<A[SIZE], 1>{public: static const int flag = 0;};
template<class A, class B> class ExCo_Enumflag<pair<A,B>, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<string, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<double, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<float, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<bool, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<char, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<unsigned char, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<short, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<unsigned short, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<int, 1>{public: static const int flag = 0;};
template< > class ExCo_Enumflag<unsigned int, 1>{public: static const int flag = 0;};
template<class A, unsigned int SIZE> class ExCo_Enumflag<Tuple<A,SIZE>, 1>{public: static const int flag = 0;};
template<class A, class B> class ExCo_Enumflag<KeyElem<A,B>, 1>{public: static const int flag = 0;};
template<class A, unsigned int SIZE> class ExCo_Enumflag<WeightElem<A,SIZE>, 1>{public: static const int flag = 0;};
template<class A, unsigned int SIZE> class ExCo_Enumflag<Trianglix<A,SIZE>, 1>{public: static const int flag = 0;};
template<class A> class ExCo_Enumflag<Vector<A>, 1>{public: static const int flag = 0;};
template<class A> class ExCo_Enumflag<AppendixPtr<A>, 1>{public: static const int flag = 0;};
template<class A> class ExCo_Enumflag<RessourcePtr<A>, 1>{public: static const int flag = 0;};


template<class NODE, unsigned int DIMS = 3u, class INTCLASS = uint32_t> class GaussCloud;
template<class C, class V>
class GaussScope{ // the weight of vector under a weight
public:
	double w,w2;
    C mean;
	V var;
    typedef GaussScope<C,V> SAFETYPE;

    GaussScope(): w(0),w2(0){}
    GaussScope(const C &pos, const V &_var, double _w = 1.0f, double _w2 = 0.0f):w(_w),w2(_w2){mean = pos * _w;var = var * _w; if (w2 == 0.0f) w2 = _w*_w; }

	GaussScope<C,V> operator+(const GaussScope<C,V>& other)const {return GaussScope<C,V>(mean+ other.mean, var+ other.var, w+other.w,w2+other.w2);}
	GaussScope<C,V> operator-(const GaussScope<C,V>& other)const {return GaussScope<C,V>(mean- other.mean, var+ other.var, w+other.w,w2+other.w2);}
	GaussScope<C,V>& operator+=(const GaussScope<C,V>& other){
		w += other.w;
        w2 += other.w2;
		mean += other.mean;
        var += other.var;
		return(*this);
	}

	GaussScope<C,V>& operator-=(const GaussScope<C,V>& other){ // TODO
		w += other.w;
        w2 += other.w2;
		mean += other.mean;
        var += other.var;
		return(*this);
	}

	C getMean() const{return( mean / w);}
	GaussScope<C,V>& toZero();
};
template <class A, class B> struct IteratorScope : public AbstractIterator<A,B> {
	enum {valid = false};
};
	template <class Node, int base_resol, int incr_resol, int nb_resol, class intclass, int nbdim, int flag> class MultiResolHashTable;
//	static void loadBMP(const char* const path,DataGrid<Tuple<unsigned char, 3>,2> &im);
//	static void saveBMP(const char* const path,const DataGrid<Tuple<unsigned char, 3>,2> &im);

	template<class prec, unsigned int nbchan, Tuple_flag flag> static void loadWAV(const char* path,DataGrid<Tuple<prec, nbchan>,1> &im);
	template<class prec, unsigned int nbchan, Tuple_flag flag> static void saveWAV(const char* path,DataGrid<Tuple<prec, nbchan>,1> &im);

 //   template<class TYPE> static void FFT_pow2(TYPE* array, unsigned char mag, bool inverse = false);

	template<class TYPE> void pow2_BLUE_FFT_routine(TYPE* data, Complex<double>* blue, unsigned char mag, unsigned int subsize); //  1 <= mult <= 7
	template<class TYPE> void pow2_BLUE_IFFT_routine(TYPE* data, Complex<double>* blue, unsigned char mag, unsigned int subsize); //  1 <= mult <= 7

	string to_string(int i);
	char* cloneString(const char* const what);
	template<class TYPE> void pow2_FFT_routine(TYPE* data, unsigned char mag); //output may need swaping
	template<class TYPE> void pow2_IFFT_routine(TYPE* data, unsigned char mag); //output may need swaping
	template<class TYPE> void pow2_IFFT_routine_swapped(TYPE* data, unsigned char mag); // assumes input is swaped, output does not need swaping
	template<class TYPE> void pow2_IFFT_routine_swapped_shift(TYPE* data, unsigned char mag, unsigned int size); // assumes input is swaped, output does not need swaping, only writes the size value, shifted to fit at the end of the array
	template<class TYPE> void comp_FFT_routine(TYPE* data, unsigned char mag, unsigned int mult);
    template<class TYPE> void comp_IFFT_routine(TYPE* data, unsigned char mag, unsigned int mult);
    template<class TYPE> static void pow2_bitswap_permutation(TYPE* array, unsigned char mag);

	// special fouriers
	template<class TYPE> void pow2_FFT_routine_zerohalf(TYPE* data, unsigned char mag); // ignores the 2nd half, assumes they are zeroes
	template<class prec, int nbchan, Tuple_flag flag> static void rescale(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);
	template<class prec, int nbchan, Tuple_flag flag> static void lowpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);
	template<class prec, int nbchan, Tuple_flag flag> static void highpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);
	template<class prec, int nbchan, Tuple_flag flag> static void lowandhighpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double low, double high);

	template<class prec, int nbchan, Tuple_flag flag> static void fourrierfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);



	template<class prec, int nbchan, Tuple_flag flag> static void extractFreqfromWAV(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);


	class ArgumentParser;
	template<int TEMPVAL>
	class templateInt{
		static const int value = TEMPVAL;
	};

	extern WarH static_warning_handdle;
	extern map< void* , pair< unsigned int, unsigned int > > man_in_the_middle;
  //  extern AliasBank static_alias_bank;
enum ODEfunctions{
    LFHP_ODE_NONE,
    LFHP_ODE_2ND_COS, // k*f(x) * cos = f"(x);
    LFHP_ODE_2ND_SIN, // k*f(x) * sin = f"(x);
    LFHP_ODE_2ND_CYLINDER //  f"(x) = k*( sin(f(x)) * u * cos(f(x)));
};
enum LFHwarnings{
    LFH_WARNING_NONE =0,
    LFH_WARNING_NAN =1,
    LFH_WARNING_MAXITE =2,
    LFH_WARNING_UNEXPECTED_INPUT=3,
    LFH_WARNING_MATRIX_SIZES_MISMATCH=4,
    LFH_ERROR_CANT_OPEN_FILE=5,
    LFH_ERROR_CANT_ALLOCATE_MEMORY=6
};
class WarH{
	int lastw;
	int wcount;
	map<int, char*> warns;
	int level;
	void reportlast(){
		if (lastw == 0) return;
		if (wcount > 1) fprintf(stderr,"%ix:", wcount);

		if (lastw < 10){
			switch((LFHwarnings)lastw){
				case LFH_WARNING_NAN: fprintf(stderr,"NaN obtained!\n"); break;
				case LFH_WARNING_MAXITE: fprintf(stderr,"Max Iteration Reached!\n"); break;
				case LFH_WARNING_UNEXPECTED_INPUT:  fprintf(stderr,"Function received unexpected input!\n"); break;
				case LFH_ERROR_CANT_OPEN_FILE:  fprintf(stderr,"Could not open file! (critical error)\n"); break;
				case LFH_ERROR_CANT_ALLOCATE_MEMORY:  myexit("Could allocate memory! (critical error)\n");
				default:  fprintf(stderr,"Unknown Warning used!\n");
			}
		}else fprintf(stderr,"some warning\n");
		fflush(stderr);
	}
public:
	WarH(): lastw(0), level(0){

	}
	~WarH(){ reportlast();
	}
	int makeWarning(char*)const{
		return(0);
	}
	void report(char* arg)const{if (level) fprintf(stderr,"some warning\n");}
	void operator<<(LFHwarnings what){if (level) this->Warning((int) what);}

	void operator()(bool test, LFHwarnings what){if (level) {if (test) this->Warning((int) what);}}
	void Warning(int w){
	    if (level) {
            if (lastw == w) wcount++;
            else{
                reportlast();
                lastw = w;
                wcount++;
            }
		}
	}
};

template<bool a> class LinkAssert{public: char junk[ a ? 1 : -1]; };

	// Polyvalent Reader Writer, which call convertions from what is read into the desired output
	class WiseIOUnit{
	public:
		template<class T> static T getFrom(const char* path);
	};

	// event type, entity array
/*
	template<class A>
	class Oper1{
	public:
		static const bool NBARG = 1;
		static const bool NBCONSTARG = 0;
		typedef A& ARG1;
		typedef void RETT;
		virtual void operator()(A &) const =0;

		template <class CA> void operator()(CA & ca) const{
			if (IteratorScope< A , CA >::valid) {
			IteratorScope< A , CA > ite;
			A* pt;
			for( pt = ite.first(ca);pt != NULL;pt = ite.next(ca)) (*this)(*pt);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}

		template <class CA> void apply(CA & ca) const{
			if (IteratorScope< A , CA >::valid) {
				IteratorScope< A , CA > ite;
				A* pt;
				for( pt = ite.first(ca);pt != NULL;pt = ite.next(ca)) (*this)(*pt);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}
	};

	template<class A, class B>
	class Oper2{
	public:
		static const bool IsPOD = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

		typedef A ARG1;
		typedef B ARG2;



		virtual void operator()(A& , B& ) const =0;
		template <class CA,class CB> void operator()(CA & ca, CB & cb) const{ ca(*this,cb);}
		template <class CA,class CB> void operator()(CA & ca, const CB & cb) const{ ca(*this,cb);}


		template <class CA> void apply(CA & ca , const B & b) const{
			if (IteratorScope< A , CA >::valid) {
				IteratorScope< A , CA > ite;
				A* pt;
				for( pt = ite.first(ca);pt != NULL;pt = ite.next(ca)) (*this)(*pt, b);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}

		template <class CA,class CB> void apply(CA & ca , CB & cb) const{ // assumes containners are of the same size!
			if ((IteratorScope< A , CA >::valid)&&(IteratorScope< B , CB >::valid)) {
				IteratorScope< A , CA > itea;
				IteratorScope< B , CB > iteb;
				A* pa; B* pb;
				for( pa = itea.first(ca) , pb = iteb.first(cb) ;pa != NULL;pa = itea.next(ca),pb = iteb.next(cb)) (*this)(*pa, *pb);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}
		template <class CA,class CB> void apply(CA & ca , const CB & cb) const{ // assumes containners are of the same size!
			if ((IteratorScope< A , CA >::valid)&&(IteratorScope< B , CB >::valid)) {
				IteratorScope< A , CA > itea;
				IteratorScope< B , CB > iteb;
				A* pa; B* pb;
				for( pa = itea.first(ca) , pb = iteb.first(cb) ;pa != NULL;pa = itea.next(ca),pb = iteb.next(cb)) (*this)(*pa, *pb);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}
	};

	template<class A, class B, class C>
	class Oper3{
	public:
		virtual void operator()(A &, B &, C &) const =0;
		template <class CA,class CB,class CC> void operator()(CA & ca, CB & cb, CC & cc) const{ ca(*this,cb,cc);}
		template <class CA,class CB,class CC> void operator()(CA & ca, CB & cb, const CC & cc) const{ ca(*this,cb,cc);}
		template <class CA,class CB,class CC> void operator()(CA & ca, const CB & cb, const CC & cc) const{ ca(*this,cb,cc);}
	};
*/
	template<int order> double cheb_eval(const Tuple<double, order> cs, const double x);

template<class O, class I, class S>
class ScaleConvert{
public:
    static const int NBARG = 1;
    typedef I ARG1;
    typedef O OUTPUT;
    S scale;
    ScaleConvert(const S& _scale): scale(_scale){}
    O operator()(const I& i) const{return (O)(i * scale);}
};




// functions!
double deflatedNBterm_approx(const double r, const double m); // testest
double deflatedNBterm(const double r, const double m); // r - log(1 + e^m) -log(exp(log(1+e^-m)e^r) -1)
double deflatedNBterm_dr(const double r, const double m); // r - log(1 + e^m) -log(exp(log(1+e^-m)e^r) -1)
double deflatedNBterm_dm(const double r, const double m); // r - log(1 + e^m) -log(exp(log(1+e^-m)e^r) -1)
double deflatedNBterm_dr2(const double r, const double m); // r - log(1 + e^m) -log(exp(log(1+e^-m)e^r) -1)
double deflatedNBterm_dm2(const double r, const double m); // r - log(1 + e^m) -log(exp(log(1+e^-m)e^r) -1)
LFH_GOLD    void deflatedNBterm_wr_f_dr_dm(double* fout3, const uint32_t x, const double r, const double m); // derivative +-0.01% accurate, does *not* contain log(!x)

void deflatedNB_getMeanAndVar(double* fout2, const double r, const double m);
void deflatedNB_getMeanAndVarSilly(double* fout2, const double r, const double m); // testing function
double deflatedNB_pvalue_silly(const uint32_t x, const double r, const double m); // is slow as x increases
double deflatedNB_pvalue(const uint32_t x, const double r, const double m);
double deflatedNB_logpvalue(const uint32_t x, const double r, const double m);
double NB_logitpvalue(const uint32_t x, const double r, const double m, const double p0);
double NB_pvalue(const uint32_t x, const double r, const double m);
double NB_logpvalue(const uint32_t x, const double r, const double m);
uint32_t deflatedNB_sample(const double r, const double m, uint32_t maxval); // is slow as x increases
double deflatedNB_pvalueSilly(const uint32_t x, const double r, const double m); // testing function
void deflatedNB_Zscore_d2(double* fout6, const uint32_t x, const double r, const double m); // is slow as x increases
void deflatedNBterm_wr_f_d2(double* fout6, const uint32_t x, const double r, const double m); // f, df/dr, df/dm, d2f/dr2, d2f/dm2

LFH_GOLD    double quasilimit_func(const double r); // x - log(log(1.0+exp(x))) >= 0
#ifdef GNU_SCIENTIFIC_LIBRARY
double log_1plusx_mx_e(const double x); // Log(1 + x) - x

Vector<Tuple<double, 5u> > derivativeAnalyser(double (&fnc)(const double),double (&der)(const double),bool posonly =false, bool show=true); // x, F[x], F'[x], relerror -1 (2 values, adapted and fixed steps)

// find intersection with 0 for a monotonic function, min > max allowed, but range spawned must contain no Nan, if no zeroes, return closest
template<class F, class I> I monotonicSolver(const F& f, const I &minimum, const I & maximum);


LFH_GOLD    double d_lngamma_dx(double x); // aka psi
LFH_GOLD    double d2_lngamma_dx2(double x);
LFH_GOLD    double d_lngamma_dx_intdiff(double x, int k); // aka psi(x+k) - psi(x)
LFH_GOLD	double Pvalue_to_stdnorm(double x);
LFH_GOLD	double LogPvalue_to_stdnorm(double x);
LFH_GOLD	double logdensity_chisquarre(double x, double free);
LFH_GOLD	double Pvalue_chisquarre_Ptail(double x, double free);
LFH_GOLD	double Pvalue_chisquarre_Ntail(double x, double free);
LFH_GOLD	double Pvalue_Gamma_Ptail(double x, double k, double theta);
LFH_GOLD	double Pvalue_Gamma_Ntail(double x, double k, double theta);
LFH_GOLD	double Pvalue_GammaRate_Ptail(double x, double alpha, double beta);
LFH_GOLD	double Pvalue_GammaRate_Ntail(double x, double alpha, double beta);
LFH_GOLD	double Pvalue_SumLogPvalues_Ptail(double sum, int nbpvals);
LFH_GOLD	double Pvalue_Beta_Ntail(double x, double a, double b);
LFH_GOLD	double Pvalue_Beta_Ptail(double x, double a, double b);
LFH_GOLD	double Pvalue_Fdistrib_Ntail(double x, double a, double b);
LFH_GOLD	double Pvalue_Fdistrib_Ptail(double x, double a, double b);
LFH_GOLD    double Pvalue_Tdistrib_Ptail(double x, double v);
LFH_GOLD	double Pvalue_Tdistrib_Ntail(double x, double v);
LFH_GOLD	double Pvalue_to_Tdistrib(double x, double v);
LFH_GOLD	double logPvalue_to_Tdistrib(double x, double v);

LFH_GOLD	double Pvalue_NBdistrib_LT(uint32_t x, double r, double p);
LFH_GOLD	double Pvalue_NBdistrib_GE(uint32_t x, double r, double p);
LFH_GOLD	double Pvalue_NBdistrib_LE(uint32_t x, double r, double p);
LFH_GOLD	double Pvalue_NBdistrib_GT(uint32_t x, double r, double p);
LFH_GOLD	double Pvalue_NBdistrib_LH(uint32_t x, double r, double p);
LFH_GOLD	double Pvalue_NBdistrib_GH(uint32_t x, double r, double p);

LFH_GOLD	double LogPvalue_NBdistrib_LE(uint32_t x, double r, double p);
LFH_GOLD	double LogPvalue_NBdistrib_GT(uint32_t x, double r, double p);
LFH_GOLD	double LogPvalue_NBdistrib_LT(uint32_t x, double r, double p);
LFH_GOLD	double LogPvalue_NBdistrib_GE(uint32_t x, double r, double p);
LFH_GOLD	double LogPvalue_NBdistrib_LH(uint32_t x, double r, double p);
LFH_GOLD	double LogPvalue_NBdistrib_GH(uint32_t x, double r, double p);
LFH_GOLD	double NBdistrib_exppara(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogProb_NBdistrib_exppara(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_GE(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LE(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_GT(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LT(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_GH(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LH(uint32_t x, double log_r, double logit_p);
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LH_exact(uint32_t x, double log_r, double logit_p); // mean = exp(log_r + logit_p)
double lfh_gsl_sf_log_beta_logit_inc_e(const double a, const double b, const double logit_x);
double log_beta_cont_frac(const double a, const double b,const double x);
LFH_GOLD	double LogPvalue_chisquarre_Ptail(double x, double free);
LFH_GOLD	double LogPvalue_chisquarre_Ntail(double x, double free);
LFH_GOLD	double LogPvalue_Gamma_Ptail(double x, double k, double theta);
LFH_GOLD	double LogPvalue_Gamma_Ntail(double x, double k, double theta);
LFH_GOLD	double LogPvalue_GammaRate_Ptail(double x, double alpha, double beta);
LFH_GOLD	double LogPvalue_GammaRate_Ntail(double x, double alpha, double beta);
LFH_GOLD	double LogPvalue_SumLogPvalues_Ptail(double sum, int nbpvals);
LFH_GOLD	double LogPvalue_Beta_Ntail(double x, double a, double b);
LFH_GOLD	double LogPvalue_Beta_Ptail(double x, double a, double b);
LFH_GOLD	double LogPvalue_Fdistrib_Ntail(double x, double a, double b);
LFH_GOLD	double LogPvalue_Fdistrib_Ptail(double x, double a, double b);

double BesselJ0(double);
double BesselJ1(double);
double BesselI0(double);
double BesselI1(double);

#endif

	double gammastar(double);
	double incgamma_dom(double, double); // D(a,x) := x^a e^(-x) / Gamma(a+1), domminant part
	double incgamma_frac(double, double); // CDP looking incomplete gamma
	double polygamma0(double);
	double sinc(double);
	double sinc_d(double);

	double mean_of_clamped_gamma(double k, double theta, double max);
	double mean_of_clamped_gamma_rate(double k, double beta, double max);
//    double mean_of_clamped_gamma2(double mean, double var, double max);

	void addRuler(DataGrid<unsigned char, 3> &RGBimage, double bounds[4], int nbsep[2]);


	double sampleGaussian();

	double distanceToEllipse(double squared_dist1, double squared_dist2, double squared_width, double squarred_focal_dist);
	double CubicRealRoot(double*, bool want_middle); // smallest root returned otherwise
    double QuarticRealRoot(double*); // todo
    double QuarticRealRootPos(double*); // find minimum positive root to quartic polynomial, return -1 if none is found

	double CubicInterpolationRoot(double* , bool is_monotonic);

    void RGB2Dvec(unsigned char *, double *); // reads 2 double, write 3 chars

enum struct_property{
	STRUCTPROP_IS_VALID=1, // is special values exists, then this checks for them
	STRUCTPROP_HAS_VALID_INVERSE=2, // check if 1.0f / X is a special value!
	STRUCTPROP_CAN_COMMUTE=4, // check is for all Y, YX = XY
};



void GPParamSearch(const Vector<KeyElem<double, double> >&, double &mean, double &var_noise, double &var_signal, double &scale);

// abstract class

	// treats A as an enum
template<typename A>
class EnumBox{
    public:
    static const bool IsPOD = true;
    static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

    A data;
    EnumBox(){}
    EnumBox(const A& f_init):data((A)f_init){}

    operator const A&() const{ return data;}
    operator A(){ return data;}

    A& operator()(){return data;}
    A operator()()const{return data;}

    bool operator>(const EnumBox<A> &o)const { return data > o.data;}
    bool operator>=(const EnumBox<A> &o)const{ return data >= o.data;}
    bool operator<(const EnumBox<A> &o)const{ return data < o.data;}
    bool operator<=(const EnumBox<A> &o)const{ return data <= o.data;}
    bool operator==(const EnumBox<A> &o)const{ return data == o.data;}
    bool operator!=(const EnumBox<A> &o)const{ return data != o.data;}

    bool operator>(const A &o)const { return data > o;}
    bool operator>=(const A &o)const{ return data >= o;}
    bool operator<(const A &o)const{ return data < o;}
    bool operator<=(const A &o)const{ return data <= o;}
    bool operator==(const A &o)const{ return data == o;}
    bool operator!=(const A &o)const{ return data != o;}

    EnumBox<A>& operator=(const EnumBox<A> &o) {data = o.data; return *this;}
    EnumBox<A>& operator=(const A &o) {data = o; return *this;}

    EnumBox<A>& toMemmove(EnumBox<A> &o) {data = o.data; return *this;}
    EnumBox<A>& toMemmove(A &o) {data = o; return *this;}

    EnumBox<A>& toZero() {data = (A)0; return *this;}
    EnumBox<A>& toOne() {data = (A)1; return *this;}

    void show(FILE* f = stdout, int level =0)const{fprintf(f,"enum%i", (int)data); if (level == 0) fprintf(f,"\n");}
};
template <class C, unsigned int DIM>
class ConstGrid{
	public:
		virtual C operator()(const typename Tuple<unsigned int, DIM>::TUPLE_TYPE &coor) const =0;
		virtual void getDims(typename Tuple<unsigned int, DIM>::TUPLE_TYPE &o_dims) const =0;
        virtual ~ConstGrid(){};
	};

/*
template<class I>
class Permutation{
public:
	I* map;
	Permutation();
	int operator[](const int &offset) const;
	Permutation<I> operator-() const;
	Permutation<I> operator+(const Permutation<size> &) const;
	template<class C, Tuple_flag Cflag, unsigned int esize> Tuple<C, size+esize, Cflag> operator+(const Tuple<C, size+esize, Cflag>&) const;
};*/

template<int buffersize>
class SuperString{
public:
	char buffer[buffersize];
	vector<int> chunks;
	SuperString();
	char* operator[](int chunkId);
	void setChunk(int chunkId, char* string, int stringlength =0);
	char* operator()();
};
enum tiffflag{
    TIFFFLAG_Compression = 259,
    TIFFFLAG_PHOTOMETRICINTERPRETATION = 262, // NEEDED , short
    TIFFFLAG_DOCUMENTNAME= 269, // ASCII
    TIFFFLAG_PAGENAME = 285, // ASCII
    TIFFFLAG_XPOSITION = 286, // rationnal
    TIFFFLAG_YPOSITION = 287, // rationnal
    TIFFFLAG_PAGENUMBER = 297 // short
    };

template<class C,unsigned int size, Tuple_flag Cflag> class IsLFHPrimitive< Tuple<C, size,Cflag> > {public: enum{ans = true};};

template<class A, class B, unsigned int size, Tuple_flag Cflag>
class STDRETTYPE2<Tuple<A,size,Cflag> , Tuple<B,size,Cflag > >{
public:
	typedef Tuple<typename STDRETTYPE2<A,B>::DEF_TYPE, size, Cflag> DEF_TYPE; // default
	typedef Tuple<typename STDRETTYPE2<A,B>::PLUS_TYPE, size, Cflag> PLUS_TYPE; // operator+
	typedef Tuple<typename STDRETTYPE2<A,B>::MINU_TYPE, size, Cflag> MINU_TYPE; // operator-
	typedef Tuple<typename STDRETTYPE2<A,B>::PROD_TYPE, size, Cflag> PROD_TYPE; // operator*
	typedef Tuple<typename STDRETTYPE2<A,B>::DIVI_TYPE, size, Cflag> DIVI_TYPE; // operator/
	typedef Tuple<typename STDRETTYPE2<A,B>::OP2_TYPE, size, Cflag> OP2_TYPE; // operator()
};

template<class A, class B, unsigned int size, Tuple_flag Cflag>
class STDRETTYPE2<Tuple<A,size,Cflag> , B >{
public:
	typedef Tuple<typename STDRETTYPE2<A,B>::DEF_TYPE, size, Cflag> DEF_TYPE; // default
	typedef Tuple<typename STDRETTYPE2<A,B>::PLUS_TYPE, size, Cflag> PLUS_TYPE; // operator+
	typedef Tuple<typename STDRETTYPE2<A,B>::MINU_TYPE, size, Cflag> MINU_TYPE; // operator-
	typedef Tuple<typename STDRETTYPE2<A,B>::PROD_TYPE, size, Cflag> PROD_TYPE; // operator*
	typedef Tuple<typename STDRETTYPE2<A,B>::DIVI_TYPE, size, Cflag> DIVI_TYPE; // operator/
	typedef Tuple<typename STDRETTYPE2<A,B>::OP2_TYPE, size, Cflag> OP2_TYPE; // operator()
};



class ArgumentParser{
	void help_routine(char c);
public:
    virtual ~ArgumentParser(){}
	void readList(char* const, Vector<unsigned int> & _out);
	void readList(char* const, vector<unsigned int> & _out);
	int operator()(int argv, char * const * args, bool has_prog_name_as_argument = true);
	int operator()(const char * argslist);

	virtual void nbaddtoken(char const * const token, int& min, int& max)=0;
	virtual void store(char* const * token, int nbtoken)=0; // return nb token used
	virtual int defstore(char* const * token, int nbtoken)=0; // return nb token used
	virtual void help() {printf("help unavailable!\n");}
};

template<int nbstate>
class HMM{
public:
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;

	//typedef LFHPrimitive::YESNO<true> IsPOD;
	TMatrix<double,nbstate,nbstate> transition;
	Tuple<double,nbstate> boundary;
	WeightElem< TMatrix<double,nbstate,nbstate>, 1 > learned_transition; // EM step buffer;


	template <Tuple_flag TF> void init(double transitrate, Tuple<double, nbstate, TF> const &);

	void EMinit();
	template<int lenght> void runHMM(Tuple<Tuple<double, nbstate>, lenght> &likelyhoods, double* weights = NULL);
	void runHMM(Tuple<double, nbstate>* likelyhoods, int lenght, double* weights = NULL);
	template< unsigned int nbdimss> DataGrid< Tuple<double, nbstate>, nbdimss> runHMM(DataGrid< Tuple<double, nbstate>, nbdimss> const &likelyhoods, int direction, bool learn = false);
	template< unsigned int nbdimss> DataGrid< Tuple<double, nbstate>, nbdimss> runHMM_meta(DataGrid< Tuple<double, nbstate>, nbdimss> const &likelyhoods, int direction, double* weights = NULL);

	void swapstates(int a,int b);

	double EMfinit();

	void show(FILE* f= stdout)const;
	ERRCODE save(FILE* f) const;
	ERRCODE load(FILE* f, unsigned int size = 0);
};

template<class A,class B>
class Convert{
	public:
	void operator()(A  &, B &) const;
	void operator()(A  &, const B &) const;
};

template<class A>
class getNorm{
public:
	void operator()(double &, A &) const;
};

// if any negative, its zero
template<class A>
class getPositiveNorm{
public:
	void operator()(double &, A &) const;
};

// template<class C> class Square : LFHDECL_OPER2(C,C);

/*
template <class C, int SIZ>
class VectorArray{
public:
	uint32_t asize;
	C* darray;
    Tuple<uint32_t, SIZ> indexes;
    Tuple<uint32_t, SIZ> sizes;

    VectorArray();
    uint32_t getSize(int index) const{return sizes[index];}
    C& pushIn(int index);

    void show(FILE* f = stdout, int level = 0) const;
};*/

// class which stores structured data
// data may be owned or not
// data may have a addresslink
template <class C, unsigned int nbdim>
class DataGrid{ //,: public ConstGrid<C,nbdim>  public TypedFunctor<C,const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE& >{
    static const bool IsPOD = false;
    static const unsigned int type_dimentions = nbdim;
	public:


	typedef C INNER_TYPE;
	typedef typename MT_IFTYPE<nbdim-1, Tuple<unsigned int, nbdim> , unsigned int >::TYPE ITERATOR_TYPE;
    typedef DataGrid<C, nbdim> SAFETYPE;


	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
    template <class S> class SUBS_INNER{public: typedef DataGrid<S, nbdim> type;};

	class KeyIterator : public AbstractKeyIterator< Tuple<unsigned int, nbdim> , DataGrid<C, nbdim> >{
	public:
		KeyIterator(const DataGrid<C, nbdim>& tar) :AbstractKeyIterator< Tuple<unsigned int, nbdim> , DataGrid<C, nbdim> >(tar){}
		bool first(){ExOp::toZero((*this).curkey); return((*this).target->data != NULL);}
		bool last(){for(unsigned int i=0;i<nbdim;i++) (*this).curkey[i] = (*this).target->dims[i]-1;return((*(Vector<C>*)this).getSize() != 0);}
		bool next(){
			unsigned int dir=0;
			for(dir=0;dir< nbdim;dir++) if ((*this).curkey[dir] == (*this).target->dims[dir]-1) (*this).curkey[dir] =0; else {(*this).curkey[dir]++; break;}
			return(dir < nbdim);
		}
		bool prev(){
			unsigned int dir=0;
			for(dir=0;dir< nbdim;dir++) if ((*this).curkey[dir] == 0) (*this).curkey[dir] =(*this).target->dims[dir]-1; else {(*this).curkey[dir]--; break;}
			return(dir < nbdim);
		}
		Tuple<unsigned int, nbdim+1> padwith(const unsigned int pad) const{
			Tuple<unsigned int, nbdim+1> fout;
			memcpy(&(fout[1]), (*this).curkey.data, sizeof(unsigned int) * nbdim);
			fout[0] = pad;
			return(fout);
		}
		void write_to(unsigned int* target) const {memcpy(target, (*this).curkey.data, sizeof(unsigned int) * nbdim);}
	};

	class Iterator{
    public:
        Tuple<uint32_t, nbdim> coor;
        C* cur;
        C* maxval;
        DataGrid<C,nbdim>& target;
        Iterator(DataGrid<C,nbdim>& trg): target(trg){}
        operator bool (){cur = target.data; maxval = cur + target.totsize(); coor.toZero(); return (target.totsize()!=0u);}
        bool operator++(int){uint32_t i; cur++; for(i=0;i<coor.getSize();i++) {if ((++coor[i]) < target.dims[i]) break; coor[i]=0;} return (i < coor.getSize());}
        Tuple<uint32_t, nbdim> operator()()const{return coor;}
        C* operator->(){return cur;}
        C& operator*(){return *cur;}
	};
    DataGrid<C,nbdim>::Iterator getIterator(){return Iterator(*this);}


    ExCoMeMdEcLaRe( LFHCONCAT2(DataGrid<C,nbdim>) );
	C* data;
	uint32_t dims[nbdim];
	DataGrid(void* owner = NULL);
	DataGrid(char* path);
	DataGrid(unsigned int* dims, void* owner = NULL);
	DataGrid(const DataGrid<C,nbdim>&);
	template<unsigned int fsize> DataGrid(const DataGrid<Tuple<C,fsize> ,nbdim-1>&);
	template<unsigned int fsize> DataGrid(const DataGrid< C[fsize]  ,nbdim-1>&);
	template<unsigned int fsize> DataGrid(const DataGrid< DataGrid<C,fsize> , nbdim-fsize >&);
	~DataGrid();

    DataGrid<C, nbdim>& toZero();
    DataGrid<C, nbdim>& toOne();
    DataGrid<C, nbdim>& toRand();
    DataGrid<C,nbdim>& toMemmove(DataGrid<C,nbdim>& source);
	DataGrid<C,nbdim-1> makeSlice(uint32_t which_dim, int coor)const;
	DataGrid<C,nbdim+1> fromSlice(uint32_t which_dim) const;
    DataGrid<C,nbdim> mkconcatenate(const DataGrid<C,nbdim> &, unsigned int direction)const;
	DataGrid<C,nbdim> selectedSlices(unsigned int direction, Vector<unsigned int> coor); // TODO
    template<unsigned int fsize> DataGrid<Tuple<C,fsize>,nbdim-1> mkTupleGrid(const Tuple<unsigned int,fsize> &which) const;

	KeyIterator getKeyIterator() const{return(KeyIterator(*this));}

	//operators

	DataGrid<C,nbdim>& operator=(const DataGrid<C,nbdim>& other);

	template<class O> DataGrid<C,nbdim>& operator=(const DataGrid<O,nbdim>& other);
	template<unsigned int S, Tuple_flag flag> DataGrid<C,nbdim>& operator=(const DataGrid<Tuple<C, S, flag>,nbdim-1>& other);
	template<class O, unsigned int S, Tuple_flag flag> DataGrid<C,nbdim>& operator=(const DataGrid<Tuple<O, S, flag>,nbdim-1>& other);
	template<Tuple_flag flag> DataGrid<C,nbdim>& operator=(const DataGrid<Tuple<C, 0u, flag>,nbdim-1>& other){LFH_FUNCTION_MISSING;}
	template<class O, Tuple_flag flag> DataGrid<C,nbdim>& operator=(const DataGrid<Tuple<O, 0u, flag>,nbdim-1>& other){LFH_FUNCTION_MISSING;}

	DataGrid<C,nbdim> mkOrthoFlip(unsigned int which_dim_source, unsigned int which_dim_target) const;
	const DataGrid<C,nbdim>& toOrthoFlip(unsigned int which_dim_source, unsigned int which_dim_target){DataGrid<C,nbdim> tmp = this->mkOrthoFlip(which_dim_source, which_dim_target);(*this) = tmp; return (*this);}
	DataGrid<C,nbdim> mkDimFlip(unsigned int which_dim) const;
	const DataGrid<C,nbdim>& toDimFlip(unsigned int which_dim);

	unsigned int totsize() const; // total size of array

	//C& operator[](unsigned int);
	//C& operator[](int);
	//C operator[](unsigned int) const;
	//C operator[](int) const;

	C& operator()(unsigned int*);
	C& operator()(int*);
	C operator()(unsigned int*) const;
	C operator()(int*) const;

	template<class O> DataGrid<O,nbdim> operator()(O (*fnct)(const C&)) const;
	DataGrid<C,nbdim>& operator()(C& (*fnct)(C&));


	typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1>,C>::TYPE operator()(unsigned int) const;
	typename MT_IFTYPE<nbdim-2, DataGrid<C,nbdim-2>,C>::TYPE operator()(unsigned int,unsigned int) const;
	typename MT_IFTYPE<nbdim-3, DataGrid<C,nbdim-3>,C>::TYPE operator()(unsigned int,unsigned int,unsigned int) const;

	C linearInterpolation(const Tuple<double,nbdim>&) const;

	C& operator()(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE & );
	C& operator()(const typename Tuple<int,nbdim>::TUPLE_TYPE&);
	C operator()(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE&) const; // constgrid func
	C operator()(const typename Tuple<int,nbdim>::TUPLE_TYPE&) const; // constgrid func
	typename Tuple<unsigned int, nbdim>::TUPLE_TYPE getDims() const; // constgrid func

	void drawLine(const Tuple<unsigned int,nbdim> &a, const Tuple<unsigned int,nbdim> &b, const C&);
	template<unsigned int size> void drawLine(const Tuple<unsigned int,nbdim-1> &a, const Tuple<unsigned int,nbdim-1> &b, const Tuple<C,size> &);

	void drawChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
	void drawCharFlip(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
	void drawUpChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
	void drawDownChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
	template<unsigned int size> void drawChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);
	template<unsigned int size> void drawCharFlip(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);
	template<unsigned int size> void drawUpChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);
	template<unsigned int size> void drawDownChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);

	Tuple<unsigned int, nbdim> addressof(C const * const item) const;

	template <int ssize> C* operator()(const Tuple<unsigned int,ssize>&);
	template <int ssize> C* operator()(const Tuple<int,ssize>&);
	template <int ssize> const C * const operator()(const Tuple<unsigned int,ssize>&) const;
	template <int ssize> const C * const operator()(const Tuple<int,ssize>&) const;

	void setRow(const C* const,const Tuple<unsigned int, nbdim>& coor, const int direction);

	template<class O> DataGrid<C,nbdim>& operator+=(KeyElem< Tuple<double,nbdim>, O > const & other);

	template<class O> DataGrid<C,nbdim>& operator+=(DataGrid<O,nbdim> const & other);
	template<class O> DataGrid<C,nbdim>& operator-=(DataGrid<O,nbdim>  const & other);
	template<class O> DataGrid<C,nbdim>& operator*=(DataGrid<O,nbdim>  const & other);
	template<class O> DataGrid<C,nbdim>& operator/=(DataGrid<O,nbdim>  const & other);
	template<class O> DataGrid<C,nbdim>& operator+=(O const & other);
	template<class O> DataGrid<C,nbdim>& operator-=(O const & other);
	template<class O> DataGrid<C,nbdim>& operator*=(O const & other);
	template<class O> DataGrid<C,nbdim>& operator/=(O const & other);



	template<class A, class B> DataGrid<C,nbdim>& toAddMult(A const & a, B const &b) {return (*this += (a*b));}
	template<class A, class B> DataGrid<C,nbdim>& toAddMult(DataGrid<A,nbdim> const &, B const &);
	//////////////////////////////
	// custom constructors
	DataGrid<C,nbdim> Crop(const Tuple<unsigned int, nbdim> &min, const Tuple<unsigned int, nbdim> &max) const;

	DataGrid<C,2> pseudoInverse(DataGrid<C,nbdim>*L = NULL, DataGrid<C,nbdim>*R = NULL) const;// return pseudo inverse, if factors R and/or L are provided, returns L * ((this)^-1) * R
	DataGrid<C,nbdim> Inverse() const; // assumes the Matrix is inversible!!!

	void initeye();

	void update_pseudoInverse(const DataGrid<C,nbdim>& target);
	void update_Inverse(const DataGrid<C,nbdim>& target);// assumes the Matrix is still inversible!!!

	void makeFreqConvolutionMap(int size, double low, double high);

	template <class D> DataGrid(DataGrid<D,nbdim> const & other);

	void setSizes(unsigned int const * const dims);
	DataGrid<C,nbdim>& toSizes(unsigned int const * const dims){this->setSizes(dims); return *this;}
	unsigned int const * const getSizes() const;
	DataGrid<C,nbdim>& setSizes(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE &dims);
//	Tuple<unsigned int,nbdim> getSizes() const;

    DataGrid< Tuple<C, nbdim> ,nbdim> makeGradient() const;
    DataGrid<double,nbdim> makeGradientNorm() const;


    DataGrid<unsigned int, nbdim> makeWatershedSortIndex(unsigned int &nbpix, bool min_gradient_bassin = true, bool tree_skip_save = true) const;
	DataGrid<unsigned int, nbdim> makeSlicesIndexes(const unsigned int nbslices) const{Vector<double> daslices; for(unsigned int i=0;i < nbslices-1;i++) daslices.push_back(((double)i+1)/nbslices); return makeSlicesIndexes(daslices);}
	DataGrid<unsigned int, nbdim> makeSlicesIndexes(const Vector<double> &fractions)const;

	template<class D> void matrixleftmultiply(const ConstGrid<D,2>& R);
	void settomatrixMultiply(const DataGrid<C,nbdim>& L, const DataGrid<C,nbdim>& R);
	void settomatrixMultiplyofTransposed(const DataGrid<C,nbdim>& L, const DataGrid<C,nbdim>& Transposed_R);

	void leftHouseHolderMultiply_lame(const double * const vec, int lenght,bool hint);
	void rightHouseHolderMultiply_lame(const double * const vec, int lenght,bool hint);
	void leftHouseHolderMultiply_lame(const double * const vec, int lenght, const double& sqrt_den,bool hint);
	void rightHouseHolderMultiply_lame(const double * const vec, int lenght, const double& sqrt_den,bool hint);
	void leftHouseHolderMultiply(const double * const vec, int lenght, const double& sqrt_den,bool hint);
	void rightHouseHolderMultiply(const double * const vec, int lenght, const double& sqrt_den,bool hint);
	void dualHouseHolderMultiply(const double * const vec, int lenght, const double& sqrt_den,bool hint); // multiply on both sides, assume squarre symetric Matrix!

	DataGrid<C,nbdim>& convolvecircle_zeroborder(double radius);

	const DataGrid<C,nbdim>& blur(const GaussianDistribution<nbdim> &);
	const DataGrid<C,nbdim>& blur_zeroborder(const GaussianDistribution<nbdim> &);
	DataGrid<C,nbdim>& blur_crude(const GaussianDistribution<nbdim> &);

	DataGrid<C,nbdim>& blur_crude_dir(unsigned int dir, double wind);


	template<class D> void leftHouseHolderMultiply(const D * const vec, int lenght, const D& sqrt_den,bool hint);
	template<class D> void rightHouseHolderMultiply(const D * const vec, int lenght, const D& sqrt_den,bool hint);

	// row columns operations

	void leftoffdiagelimination(double factor, int from, int to);// add columns (top) to // add columns
	void rightoffdiagelimination(double factor, int from, int to);

	void positivedefinite_inverse_and_involution(DataGrid<C, nbdim> &f_out, const DataGrid<C, nbdim> &invol);
	void solveLinearSystem(const C* const target, C * _out);
	void LinearTransform(const C* const target, C * _out);

	// this computes (R*M*(R^-1))^-1, where M is a large Matrix, and R a reduction
	void downsampled_Matrix(const DataGrid<double, nbdim> &f_inner,const DataGrid<double, nbdim> &f_reduction);

	void Matrix_factor_out_from_symetric_Matrix(DataGrid<C, nbdim> &io_symmetric_Matrix);

	void makeDiagonalizer(const ConstGrid<double,2> &symmetric_Matrix, Vector<double> &eigenval);
	void makeDiagonalizer_ofinverse(const ConstGrid<double,2> &symmetric_Matrix, Vector<double> &eigenval);

	template<unsigned int size> void makeDiagonalizer(const TMatrix<double,size,size> &symmetric_Matrix, Vector<double> &eigenval);
	template<unsigned int size> void makeDiagonalizer_ofinverse(const TMatrix<double,size,size> &symmetric_Matrix, Vector<double> &eigenval);

	void solve_sym_and_project(DataGrid<double, 2> &project, const ConstGrid<double,2> &symmetric_Matrix); // TODO // computes : project dot symm_Matrix^-1 dot (this)

	void makeInverse(const ConstGrid<double,2> &symmetric_Matrix);

	void LinearSystem_fixpoint_solve(const DataGrid<C, nbdim>&target, double (*weights)(int a, int b));
	void LinearSystem_residual(const DataGrid<C, nbdim>&guess, const DataGrid<C, nbdim>&target, const ConstGrid<double,2>&);

	DataGrid<int, nbdim> directNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&), bool equalitycount = false);
	DataGrid<int, nbdim> indirectNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&), bool equalitycount = false);
	DataGrid<int, nbdim> knightNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&), bool equalitycount = false);

	vector< Tuple<unsigned int,nbdim > > getCoor_Mathching(const C&)const;

	unsigned int get_directNeightbor(const Tuple<unsigned int , nbdim> & query, Tuple< Tuple<unsigned int , nbdim>,  nbdim * 2> &_out )const;
	unsigned int get_indirectNeightbor(const Tuple<unsigned int , nbdim> & query, Tuple< Tuple<unsigned int , nbdim>,  (unsigned int)(TEMPLATE_INT_POWER<3,nbdim>::ans) -1u > &_out )const;
	unsigned int get_knightNeightbor(const Tuple<unsigned int , nbdim> & query, Tuple< Tuple<unsigned int , nbdim>, nbdim * (nbdim -1) * 4  > &_out )const;

	unsigned int getDimSize(unsigned int dim);
	void getDims(unsigned int* out_dim);
	void exponentialBlur(double std_dev, bool circular);
	template<unsigned int size, Tuple_flag Cflag> DataGrid<Tuple<C, size,Cflag>, nbdim> exponentialBlur(Tuple<double,size,Cflag> apertures, bool circular);
	template<unsigned int size, Tuple_flag Cflag> DataGrid<Tuple<C, size,Cflag>, nbdim> crudeGaussianBlur(Tuple<double,size,Cflag> apertures, bool circular);
	DataGrid<C, nbdim> crudeGaussianBlur(double aperture, bool circular = false);

	DataGrid<C,nbdim>& toresize(unsigned int* newdims);
	DataGrid<C,nbdim> resize(unsigned int* newdims) const;
	DataGrid<C,nbdim> resize_dir(unsigned int dims, unsigned int size) const;

	DataGrid<C,nbdim>& toblur(double aperture){for(unsigned int i=0;i<nbdim;i++) this->toblur_dir(aperture,i); return *this;}
	DataGrid<C,nbdim>& toblur_dir(double aperture,unsigned int d);

	// { NOT WORKING
	DataGrid<C,nbdim>& toblur_zeroborder(double aperture){for(unsigned int i=0;i<nbdim;i++) this->toblur_zeroborder_dir(aperture,i); return *this;}
	DataGrid<C,nbdim>& toblur_zeroborder_dir(double aperture,unsigned int d);
	// } NOT WORKING

	DataGrid<C,nbdim>& toresize_crude(unsigned int* newdims);
	DataGrid<C,nbdim> resize_crude(unsigned int* newdims) const;
	DataGrid<C,nbdim> resize_dir_crude(unsigned int dims, unsigned int size) const;

	DataGrid<C,nbdim>& toScaledWithConstantFrequencies(unsigned int dir, unsigned int factor);
	DataGrid<C,nbdim>& toInterleaved(unsigned int dir, unsigned int factor, unsigned int remainder);


	DataGrid<double ,nbdim>  ExpectedDistanceToOther(double exit_prob) const;


	template<class ARR_double> void ExpectedDistanceToOther(const DataGrid< ARR_double ,nbdim> &source, const DataGrid<double, 2>& transition);
	template<class ARR_double> void ExpectedDistanceToOther_instate(const DataGrid< ARR_double ,nbdim> &source,const DataGrid<double, 2>& transition, int state);

	void ExpectedDistanceToOther_I(const DataGrid< double ,nbdim> &source, const DataGrid<double, 2>& transition);
	void ExpectedDistanceToOther_instate_I(const DataGrid< double ,nbdim> &source,const DataGrid<double, 2>& transition, int state);

	void ExpectedDistanceToOther(const DataGrid< double ,nbdim+1> &source, const DataGrid<double, 2>& transition, bool state_zero_special = false);
	void ExpectedDistanceToOther_instate(const DataGrid< double ,nbdim+1> &source,const DataGrid<double, 2>& transition, int state);


	Vector< Tuple<unsigned int, nbdim> > get_ordered_Data_Coors() const;
	Vector< KeyElem<C, Tuple<unsigned int, nbdim> > > get_ordered_Data_Coors_Tagged() const;

	DataGrid< Tuple<unsigned int, nbdim> ,nbdim> LocalMaxMap() const;
	DataGrid< Tuple<unsigned int, nbdim> ,nbdim> LocalMinMap() const;

	DataGrid< unsigned int ,nbdim> SegementClimbToMax(unsigned int &nbsegs, bool dist_penal = false, bool knight = false , bool eq_join = true, bool mark_da_maxima = false, const DataGrid< bool ,nbdim>* ignore_filter = NULL) const;
	DataGrid< unsigned int ,nbdim> SegementClimbToMaxMax(unsigned int &nbsegs, bool eq_join = true) const;
	DataGrid< unsigned int ,nbdim> SegementClimbToMin(unsigned int &nbsegs, bool dist_penal = false, bool knight = false , bool eq_join = true, bool mark_da_maxima = false, const DataGrid< bool ,nbdim>* ignore_filter = NULL) const;

	DataGrid<C,nbdim> FFtransform_old() const;
	DataGrid<C,nbdim> invFFtransform_old() const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform() const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform() const;

	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_2dir(unsigned int dir, Bluestein& blue, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_2dir(unsigned int dir, Bluestein& blue, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_2(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_2(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

    // emulate a transform on a twice bigger array, where F[x] = F[-x] so that F'[w] = F'[-w], hence there is no need to have the
    DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_even(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_even(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_zeroborder(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_zeroborder(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_zeroborder_dir(unsigned int dir, Bluestein &blue, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_zeroborder_dir(unsigned int dir, Bluestein &blue, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

	void applycircleFilter(double minx, double maxx);
	DataGrid<C,nbdim> circleFilter(double minx, double maxx);

	DataGrid<unsigned int, nbdim> LabelConnected( bool (*isToLabel)(const C&), unsigned int *out_nblabels = NULL, Vector< Tuple<unsigned int, nbdim*2> >* out_boundingrect= NULL);

	template<class A_1, class A_2, class C_2, unsigned int size_2> void operator() (Oper2<A_1,A_2> const & op, DataGrid<C_2,size_2> const & _in ); // not a match
	template<class C_2, unsigned int size_2> void operator() (Oper2<C,C_2> const & op, DataGrid<C_2,size_2> const & _in); // match
	template<class A_1, class A_2, class C_2, unsigned int size_2> void operator() (Oper2<A_1,A_2> const & op, DataGrid<C_2,size_2> & _in ); // not a match
	template<class C_2, unsigned int size_2> void operator() (Oper2<C,C_2> const & op, DataGrid<C_2,size_2> & _in); // match

	template<class A_1, class A_2, class A_3, class C_2, unsigned int nbdim_2, class C_3, unsigned int nbdim_3> void operator() (Oper3<A_1,A_2,A_3> const & op, DataGrid<C_2, nbdim_2>  &, DataGrid<C_3, nbdim_3>  &); // not a match
	template<class C_2, unsigned int nbdim_2, class C_3, unsigned int nbdim_3> void operator() (Oper3<C,C_2,C_3> const & op, DataGrid<C_2, nbdim_2>  &, DataGrid<C_3, nbdim_3> &); // match

	class FiniteDifference:  LFHDECL_OPER2(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<C, nbdim>));
	class FiniteSum:  LFHDECL_OPER2(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<C, nbdim>));

	class Convolution : LFHDECL_OPER3(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<C, nbdim>), LFHCONCAT2(DataGrid<double, nbdim>));
	template<class D> class RankMap:  LFHDECL_OPER2(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<D, nbdim>));

	DataGrid<C, nbdim>& toFractalTexture(uint32_t mag = 8u, double blur_aperture = 0.0f, bool isSigned =false);
	DataGrid<C, nbdim>& toFractalTexture(uint32_t mag, double blur_aperture, const C minVal, const C maxVal);

	#ifdef Rcpp_hpp
    template<class RCL> void rdMatrix(const arma::Mat<RCL> &where);
    template<class RCL> void wrMatrix(arma::Mat<RCL> &where) const;
    #endif

}; // end of class DataGrid

// algorithm needs (input and model)
// step1 span parameter-space
// step2 learn
// step3 enjoy

// data can be stored as dense, ordered or hashed


// your friendly class which loves parallel iterators

// 1D tensor... Hashmap or Tuple

class OrderedTupleIterator{
public:
    uint8_t nbpages;
    uint8_t nbdims;
    uint8_t** data;
    uint8_t* curindex;
    bool start(int whichpage);
    bool next(int whichpage);
};


template<class C = Anything>
class SizedTensor{
public:
    Vector<C> data;
    Sparsity mapping;
    auto operator()(const uint32_t *coor, const Tuple<uint32_t> &dims)->decltype(data[0]);
    auto operator()(const uint32_t *coor, const Tuple<uint32_t> &dims)const ->decltype(data[0]);
};

template<class C, int D, class I = uint32_t, Tuple_flag Cflag =TUPLE_FLAG_NULL>
class TTensor{
public:
    Tuple<C, Cflag> data;
    I* coors;
    Tuple<uint32_t, D> dims;
};

class SparseTensor{ // polymorphic indexing structure
public:
    // if X is positive, dense with X dimentions
    // if X is 0, 1 dim ordered with offsets
    // if X is negative, -X dimensions are stored as sparse
	myHashmap<uint32_t, myHashmap<uint32_t> > indexes;
    Tuple<uint32_t> dims;
	class Iterator{
    public:
        SparseTensor& target;
        Iterator(SparseTensor& trg): target(trg){}
        /*Tuple<uint32_t> coor;
        C* cur;
        C* maxval;


        operator bool (){cur = target.data; maxval = cur + target.totsize(); coor.toZero(); return (target.totsize()!=0u);}
        bool operator++(int){uint32_t i; cur++; for(i=0;i<coor.getSize();i++) {if ((++coor[i]) < target.dims[i]) break; coor[i]=0;} return (i < coor.getSize());}
        Tuple<uint32_t> operator()()const{return coor;}
        C* operator->(){return cur;}
        C& operator*(){return *cur;}*/
	};

    SparseTensor::Iterator getIterator(){return Iterator(*this);}
    Anything operator()(const Tuple<uint32_t> &coor);
    AnythingConst operator()(const Tuple<uint32_t> &coor) const;
    void show(FILE*f =stdout, int level =0)const;
}; // end of class Tensor
/*
template< >
class SparseTensor<Anything>{
public:

	class KeyIterator{
    public:
        Tuple<uint32_t> coor;
        SparseTensor<Anything>* target;
        KeyIterator(){}
//        operator bool (){cur = target.data; maxval = cur + target.totsize(); coor.toZero(); return (target.totsize()!=0u);}
//        bool operator++(int){uint32_t i; cur++; for(i=0;i<coor.getSize();i++) {if ((++coor[i]) < target.dims[i]) break; coor[i]=0;} return (i < coor.getSize());}
        Tuple<uint32_t> operator()()const{return coor;}
        void initPartition(uint32_t id, uint32_t nbpartition);
	};
    SparseTensor<Anything>& toZero();
    Tuple< SparseTensor<Anything>::KeyIterator > getIterators(uint32_t nbiterators){Tuple<SparseTensor<Anything>::KeyIterator > fout; fout.toSize(nbiterators); for(uint32_t i=0;i<nbiterators;i++) fout[i].target = this; return fout;}
	void show(FILE*f =stdout, int level =0)const;
};
*/
// tensor has 1-4 dims, can be indexes or doubles



/*
class SymbolicType{
public:
    uint32_t anytype; // flag for attributes
    SymbolicType();
    ~SymbolicType();
    //uint32_t getBlockByteSize()const;
    SymbolicType& toZero();
    ERRCODE load(FILE *f, unsigned int size=0);
	ERRCODE save(FILE *f) const;
};*/


enum LFH_FUNCTION_enum{
    LFH_FUNCTION_ADDITION,
    LFH_FUNCTION_SUBSTRACTION,
    LFH_FUNCTION_MULTIPLICATION,
    LFH_FUNCTION_DIVISION,
    LFH_FUNCTION_LL_NORMAL,  // x, mu, sigma
    LFH_FUNCTION_LL_NEG_BINOMIAL // x, r, m
};
class FunctionStruct{
public:
    LFH_FUNCTION_enum what;

};

class ParallelTaskNode{
public:
    typedef uint32_t INDEX_TYPE;
    uint32_t alias;
    std::function<void(uint32_t)> fnc;
    Tuple<uint32_t> input_aliases;
    Tuple<uint32_t> thread_input;
	uint32_t& getIndex(){return alias;}
	uint32_t getIndex() const {return alias;}
};

// Optimization... gradient assent, greedy search and constrain programming
// it's all about defining a parameter-wise differentiable loss-function
// a lot is tensor operator, which uses iterators to populate
class SymbolicScope{
    template<class ...A> Tuple<uint32_t> mkAliases_routine(Vector<uint32_t> &fout, const char* val, A...);
    template<class ...A> Tuple<uint32_t> mkAliases_routine(Vector<uint32_t> &fout, string val, A...);
    template<class ...A> Tuple<uint32_t> mkAliases_routine(Vector<uint32_t> &fout, uint32_t val, A...);
    Tuple<uint32_t> mkAliases_routine(Vector<uint32_t> &fout);
    uint32_t currentTarget;
    void clearData(uint32_t alias);
public:
    dualHashmap<string, uint32_t> alias;
    myHashmap<uint32_t, uint32_t> typeflags; // 1 is double, 2 is ptr, 4 has dimensions, 8 is sparse , 16+ special classes
    myHashmap<uint32_t, myHashmap<string, uint32_t> > members; // aliases for related vars
    myHashmap<uint32_t, Tuple<double> > double_data;
    myHashmap<uint32_t, Tuple<uint32_t> > int_data;
    myHashmap<uint32_t, double* > double_ptr;
    myHashmap<uint32_t, uint32_t* > int_ptr;
    myHashmap<uint32_t, Tuple<uint32_t> > dim_data;

    myHashmap<ParallelTaskNode, void, defaultHashFnc< uint32_t > > fncs;

//    myHashmap<uint32_t, Dictionary> name_data;
    myHashmap<uint32_t, myHashmap<uint32_t, uint32_t> > index_data; // compressed one day? (aliasmap)

    //myHashmap<uint32_t, uint32_t> index_attribute;
    //myHashmap<uint32_t, OrderedTupleIterator > indexes;
    //Tuple<uint32_t> getDims(uint32_t alias)const;

    ThreadBase &tb;
    SymbolicScope (ThreadBase &_tb): tb(_tb){}

    void deleteAlias(uint32_t aliasA);
    void swapVars(uint32_t aliasA, uint32_t aliasB);

    SymbolicScope& operator[](const char* str);

    SymbolicScope& setTarget(uint32_t _target);
    SymbolicScope& operator=(uint32_t alias);
    SymbolicScope& operator+=(uint32_t alias);
    SymbolicScope& operator-=(uint32_t alias);
    SymbolicScope& operator*=(uint32_t alias);
    SymbolicScope& operator/=(uint32_t alias);

    SymbolicScope& operator=(Tuple<double> &&value);
    SymbolicScope& operator=(Tuple<uint32_t> &&value);
    SymbolicScope& operator=(const SparseTuple<double> &value);
    SymbolicScope& operator=(const Tuple<double> &value);
    SymbolicScope& operator=(const Tuple<uint32_t> &value);
    SymbolicScope& operator=(const TMatrix<double> &data);

    void optimize(const char* str);

    uint32_t assignFunctor( std::function<void(uint32_t)> fnc, Tuple<uint32_t> &input, Tuple<uint32_t> &thread_arg);
   // template<class F> uint32_t assignFunctor(const F& funct);


    void writeTo( TMatrix<double>& targ)const;
    void writeTo( Tuple<double>& targ)const;
    void runFunctor(uint32_t alias);

    SymbolicScope& cmp_ll_normal(uint32_t aliasX,uint32_t aliasMu,uint32_t aliasSigma);

    // populate iterator
    uint32_t prepareIterators(uint32_t alias, uint32_t nbthreads);
    // prepare argument variables as proxy variables for a loop
    bool first(uint32_t threadAlias, uint32_t threadInfo);
    bool next(uint32_t threadAlias, uint32_t threadInfo);

    uint32_t assignAlias();
    uint32_t assignAlias(string newname);
//    SparseTensor& operator[](uint32_t alias){return variables[alias];}
//    SparseTensor& operator[](string name){return (*this)[alias[name]];}
    template<class ...A> Tuple<uint32_t> mkAliases(A... args);
	ERRCODE load(FILE *f, unsigned int size=0);
	ERRCODE save(FILE *f) const;
	void show(FILE*f, int level);
};

class IndexMap{
public:
    uint32_t symbAlias;
    SparseTensor evaluate(const SymbolicScope& scope, const Tuple<uint32_t> &coor);
};
class FunctionNode;
class Symbolic{
public:
    SymbolicScope &symscp;
    uint32_t alias;
    Symbolic(SymbolicScope &target): symscp(target), alias(symscp.assignAlias()){}
    ~Symbolic(){symscp.deleteAlias(alias);}
    FunctionNode operator+(const Symbolic &other) const;
    FunctionNode operator-(const Symbolic &other) const;
    FunctionNode operator*(const Symbolic &other) const;
    FunctionNode operator/(const Symbolic &other) const;
    Symbolic& operator=(FunctionNode &fnc);
    Symbolic& operator=(FunctionNode &&fnc);
    Symbolic& operator=(const Tuple<double>& data);
    Symbolic& operator=(Tuple<double>&& data);
    void show(FILE* f =stdout, int level=0)const;
};

class FunctionNode{
public:
    SymbolicScope &symscp;
    LFH_FUNCTION_enum func;
    uint32_t alias;
    Tuple<uint32_t> argaliases;
    FunctionNode(SymbolicScope& _symscp): symscp(_symscp),alias(symscp.assignAlias()){}
    ~FunctionNode(){symscp.deleteAlias(alias);}
    FunctionNode& setFunction(LFH_FUNCTION_enum what);
    void setFunction(string name, LFH_FUNCTION_enum what, Tuple<uint32_t> &&input_aliases);
    void eval(uint32_t target);
    void eval(){return this->eval(alias);}
    void backward(uint32_t error_alias);
};

// takes 2 vectors to select row and collumn in a tensor, to get a vector as parameter
class TensorInnerProduct{
public:
    DataGrid<double, 3u> tensor;
    SymbolicScope& sm;
    uint32_t evalalias;

    Accessor< uint32_t, myHashmap<uint32_t, double> > left;
    Accessor< uint32_t, myHashmap<uint32_t, double> > right;
    IteratorMaker< uint32_t, Tuple<double, 0u, TUPLE_FLAG_REMOTE_MEMORY> > out_imk;
    uint32_t compute_nbthread;
    class Evaluator{
    public:
        TensorInnerProduct& trg;

        Evaluator(TensorInnerProduct& _trg): trg(_trg){}
        void operator()(uint32_t exec_alias);
    };
    TensorInnerProduct& setDims(const Tuple<uint32_t, 3u> &_dims){tensor.setSizes(_dims); return *this;}
    const uint32_t* getDims(const Tuple<uint32_t, 3u> &_dims)const {return tensor.dims;}
    TensorInnerProduct& toRand(){tensor.toRand(); return *this;}

    TensorInnerProduct(SymbolicScope& _sm);
    void run(ThreadBase &tb);

};




class GaussFunctor{
public:
    SymbolicScope& sm;
    uint32_t evalalias;

    class Evaluator{
        public:
        GaussFunctor& trg;
    //    Iterator< KeyElem<uint32_t, Tuple<double, 0u, TUPLE_FLAG_NULL> > > daite;
    //    Assessor< uint32_t, Tuple<double, 0u, TUPLE_FLAG_NULL> > damu;
        Evaluator(GaussFunctor& _trg): trg(_trg){}
        void operator()(uint32_t exec_alias);
    };

    class Backward{
    public:
        GaussFunctor& trg;
        Backward(GaussFunctor& _trg): trg(_trg){}
        void operator()(uint32_t exec_alias);
    };

    Accessor< uint32_t, Tuple<double, 0u, TUPLE_FLAG_REMOTE_MEMORY> > d_mu_acs;
    Accessor< uint32_t, Trianglix<double> > d_precisqrt_acs;

    Accessor< uint32_t, Tuple<double, 0u, TUPLE_FLAG_REMOTE_MEMORY> > mu_acs;
    Accessor< uint32_t, Trianglix<double> > precisqrt_acs;
    Accessor< uint32_t, Tuple<double, 0u, TUPLE_FLAG_REMOTE_MEMORY> > data_acs;
    IteratorMaker< uint32_t, double > ll_imk;

    uint32_t compute_nbthread;

    GaussFunctor(SymbolicScope& sm);
    void operator()(uint32_t exec_alias);

    void eval();

    void run(ThreadBase &tb);

};



/*
   class for symmetric square matrices
*/

class EigenFixer{
public:
    typedef Tuple<double> INTYPE;
    void operator()(Tuple<double> & eigenvals) const;
};

// 3D Symmetric Tensor
template<class A, class B, class C>
class Tetraglix{
public:
    Tuple<A> diag;
    B* plane;
    C* volum;
};

template<int nbstate>
class RevertibleHMM{
public:
    Trianglix<double> transition;
};

template<class C, unsigned int SIZE, unsigned int ORDER>
class Pyramidix{
};

template<class C, unsigned int SIZE>
class Pyramidix<C,SIZE, 1u> : public Tuple<C,SIZE> {
};

template<class C, unsigned int SIZE>
class Pyramidix<C,SIZE, 2u> : public Trianglix<C,SIZE> {
};


template<class C>
class Matrix{ // same encoding, different interpretation
	template<class O, class B> void alloc_add(Matrix<O>& fout, const Matrix<B>& other)const; // assumes fout.data == NULL
	template<class O, class B> void mkMult(Matrix<O>& fout, const B& other, Tuple< char, 0> *)const;
	template<class O, class B> void mkMult(Matrix<O>& fout, const B& other, Tuple< char, 1> *)const;
	template<class O, class B> void mkMult(Matrix<O>& fout, const B& other, Tuple< char, 2> *)const;
    void InitFromTrianglix(const C* data);
public:
    typedef YESNO< false > IS_COMMUTATIVE;
	unsigned int sizex;
	unsigned int sizey;
	C* data;
	Matrix(): data(NULL){}
	Matrix(unsigned int _sizex, unsigned int _sizey): sizex(_sizex), sizey(_sizey),data(new C[_sizex*_sizey]){}
	Matrix(const DataGrid<C,2>& other):sizex(other.dims[0]), sizey(other.dims[1]) {data = new C[sizex*sizey];}

    Matrix(const Trianglix<C,0u>& other):sizex(other.t_size), sizey(other.t_size) {data = new C[sizex*sizey];InitFromTrianglix(other.data);}
    template <unsigned int SIZE> Matrix(const Trianglix<C,SIZE>& other):sizex(SIZE), sizey(SIZE) {data = new C[sizex*sizey];InitFromTrianglix(other.data);}
    Matrix<C>& operator=(const Matrix<C>& other) {delete[](data); sizex = other.sizex; sizey = other.sizey; data = new C[sizex*sizey];for(unsigned int i=0;i<sizex*sizey;i++) data[i] = other.data[i]; }
    Matrix<C>& operator=(const Trianglix<C,0u>& other) {delete[](data); sizex = other.size; sizey = other.size; data = new C[sizex*sizey];InitFromTrianglix(other.data);}
    template <unsigned int SIZE> Matrix<C>& operator=(const Trianglix<C,SIZE>& other) {delete[](data); sizex = SIZE; sizey = SIZE; data = new C[sizex*sizey];InitFromTrianglix(other.data);}

	~Matrix(){delete[](data);}
    void setSizes(unsigned int x, unsigned int y){delete[](data); sizex = x; sizey=y; data = new C[sizex*sizey];}

	Matrix<C> operator+(const Matrix<C>& other)const;
	Matrix<C> operator-(const Matrix<C>& other)const;
	Matrix<C> operator-()const;
	Matrix<C>& operator+=(const Matrix<C>& other);
	Matrix<C>& operator-=(const Matrix<C>& other);
    Matrix<C>& toZero(){for(unsigned int i=0;i<sizex*sizey;i++) ExOp::toZero(data[i]); return(*this); }
    Matrix<C>& toRand(){for(unsigned int i=0;i<sizex*sizey;i++) ExOp::toRand(data[i]); return(*this); }


    Matrix<C> operator*(const Matrix<C>& other)const;

	Matrix<C> inverse() const;// TODO!

    double pnorm() const;
    double norm() const{return(sqrt(pnorm));}

    void show(FILE *f=stdout, int level=0)const;
};




template<class C, class D, class FUNC, unsigned int DIM>
class FuncGrid: public ConstGrid<C,DIM>{
	unsigned int start;
	unsigned int rangesize;
	public:
	Vector<D> data;
	FUNC* cl; // must implement one of:
	/*
	//	C operator()(const D*, const Tuple<unsigned int, DIM>& );
	//	C operator()(const D,const D,const D,... );
		C operator()(const Tuple<D,DIM> );
	*/

	FuncGrid(): start(0xFFFFFFFF) {}
	C operator()(const Tuple<unsigned int, DIM>&coor)const {
		if (start == 0xFFFFFFFF) return (*cl)(data.darray, coor);
		Tuple<unsigned int, DIM> ncoor = coor;
		unsigned int i;
		for(i=0;i<DIM;i++) ncoor[i] += start;
		return (*cl)(data.darray, ncoor);
		}
	void specialize_range(int i_start, int i_rangesize){ start = i_start; rangesize = i_rangesize;}
	void specialize_range_clear(){ start = 0xFFFFFFFF;}

	void getDims(Tuple<unsigned int, DIM> &o_dims) const{
		unsigned int dim = (start == 0xFFFFFFFF) ? data.getSize() : rangesize;
		unsigned int i;for(i=0;i<DIM;i++) o_dims[i] = dim;
		}
	void show(FILE* f = stdout) const{ // ignore shifts
		int ds = data.getSize();
		if (DIM == 1){
			/*Tuple<unsigned int, 1> ite1;
			for(ite1[0]=0;ite1[0]<ds;ite1[0]++) {
				ExOp::show( (*cl)(data.darray, ite1) , f,1); fprintf(f,"%c", ite1[0] == ds-1 ? '\n' : '\t');
			}*/
		}else if (DIM == 2){
			Tuple<unsigned int, 2> ite2;
			for(ite2[1]=0;ite2[1]<ds;ite2[1]++)
			for(ite2[0]=0;ite2[0]<ds;ite2[0]++) {
				ExOp::show((*cl)(data.darray, ite2) , f,1); fprintf(f,"%c", ite2[0] == ds-1 ? '\n' : '\t');
			}


		}
		}


		void CholeskyDecom(C* &f_out){ // to do, maybe
			int k,l;
			Tuple<unsigned int, DIM> coor;
			if (DIM == 2){
				k = data.getSize();
				f_out = new C[(k* (k+1)) /2];
				k=0;
				for(coor[1]=0;coor[1]<data.getSize();coor[1]++){
					for(coor[0]=0;coor[0]<coor[1];coor[0]++){
						f_out[k] = (*this)(coor);
				//		for(l = coor[0]-1; l>=0;l--) f_out[k] -= f_out[k+l-coor[0]] * f_out[k+l-coor[0]] * f_out[(l * (l+1)/2)-1];
						f_out[k] /= f_out[(coor[0] * (coor[0]+3))/2];

						k++;
					}
					f_out[k] = (*this)(coor);
					for(l = coor[0]-1; l>=0;l--) f_out[k] -= f_out[k+l-coor[0]] * f_out[k+l-coor[0]] * f_out[(l*(l+3))/2];
					k++;
				}
			}
		}
	};




class TiffFile{
	template<class TYPE> void writeDescriptor(TYPE min, TYPE max,char * &p, int nbsample);
public:
	FILE* f;
	bool inv;
	Vector<unsigned char> curflaglist;
	TiffFile(const char*, bool writeonly = false); // erase file if present
	~TiffFile();
	unsigned int curfp_pos;
	unsigned int endfile_pos;
	bool gotoNext();
	template <class C, unsigned int channels, Tuple_flag Cflag> bool fetch(DataGrid<Tuple<C,channels, Cflag>,2> &out);
	template <class C> bool fetch(DataGrid<C,3u>& f_out, char * imageType = NULL);
	template <class C> bool fetch(DataGrid<C,2u>& f_out, const unsigned int channel = 0);
	template <class C> bool fetch_grayscale(DataGrid<C,2>& f_out);
	bool isEmpty() const; //


	template <class C> bool put(DataGrid<C,3u>& f_out, bool updateRange =false);
	template <class C, unsigned int SIZE, Tuple_flag Cflag> bool put(DataGrid<Tuple<C, SIZE, Cflag>,2u>& f_out, bool updateRange =false);

	template <class C, class TYPE> bool put(DataGrid<C,3u>& f_out, TYPE min, TYPE max, bool updateRange =false);
	template <class C, class TYPE> bool put(DataGrid<C,2u>& f_out, TYPE min, TYPE max, bool updateRange =false);
	template <class C, class TYPE, unsigned int SIZE, Tuple_flag Cflag> bool put(DataGrid<Tuple<C,SIZE, Cflag>,2u>& f_out, TYPE min, TYPE max, bool updateRange =false);

    void savePLYasTIFF(const char * path);

	void addopt_Description(Vector< char* > & opt, char const * const descript) const;
	void addopt_Xscale(Vector< char* > & opt, double min_x, double max_x) const;
	void addopt_Yscale(Vector< char* > & opt, double min_y, double max_y) const;
	//	template <class TYPE>void addopt_Xscale(Vector< char* > & opt, TYPE min_x = ExCo<TYPE>::zero(), TYPE max_x = ExCo<TYPE>::zero() ) const;

	//	template <class C> bool put(const DataGrid<C,3>& f_out, Vector< char* > &options);


	//	template <class C, int channels, int flag, Tuple_flag Cflag> void put(vector< DataGrid<Tuple<C,channels, Cflag>,1> > &out);
	int flagType(int flagindex);

    class WriteScope{
	public:
        unsigned int nbflags;
        //Vector<char> extra_data;
    };
    WriteScope* current;
    void startFrameWrite();
    void endFrameWrite();

	template <class C> C getValue(int flagindex);
};


template< > vector<long int> TiffFile::getValue< vector<long int> >(int flagindex);
template< > vector<unsigned int> TiffFile::getValue< vector<unsigned int> >(int flagindex);
template< > unsigned int TiffFile::getValue<unsigned int>(int flagindex);
template< > int TiffFile::getValue<int>(int flagindex);

template<class C, int nbc>
class imageIO{
	public:
	static void importBMP(char* path, DataGrid<Tuple<C, nbc>,2> &out);
	static void exportBMP(char* path, DataGrid<Tuple<C, nbc>,2> &out);
};



template<class O, class I, int _out_dim, int _in_dim>
class ContinuousFunction : public Oper2< Tuple<O, _out_dim> , Tuple<I, _in_dim > >{
public:
	virtual void operator()(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &) const =0;
	virtual void derivative(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct) const;
	virtual void derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &, Tuple<I, _in_dim > &) const;
	void derivativeMatrix_default(TMatrix<O, _in_dim, _out_dim> &, Tuple<I, _in_dim > &) const;
	void derivative_from_difference(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct) const;
	void newtonStep(Tuple<I, _in_dim > &);
};

template<class O, class I, class A>
class FunctionScope{
public:
    O (*f)(I,A);
    A arg;
    FunctionScope(O (*_f)(I,A), A _arg): f(_f),arg(_arg){}
    O operator()(I i)const{return (*f)(i,arg);}
};

/*
template<class C, Tuple_flag Cflag>
class Tuple<C,1u, Cflag>{
public:
    myHashmap<uint32_t, C> data;
    typedef C TUPLE_TYPE;
    static C accessTupleType(const C& val,uint32_t offset){return val;}
    static C& setTupleType(C& val,uint32_t offset){return val;}


    class Iterator{
        uint32_t iteval;
        public:
        Tuple<C,1u, Cflag>& target;
        Iterator(Tuple<C,1u, Cflag>& trg): target(trg){}
        operator bool (){if (data.getSize() == 0) return false; iteval=0; return true;}
        bool operator++(int){iteval++; return iteval < data.getSize();}
        uint32_t operator()()const{data.deref_key(iteval);}
        C* operator->(){return &(data.deref(iteval));}
        C& operator*(){return data.deref(iteval);}
	};
    class ConstIterator{
        uint32_t iteval;
        public:
        const Tuple<C,1u, Cflag>& target;
        ConstIterator(const Tuple<C,1u, Cflag>& trg): target(trg){}
        operator bool (){if (data.getSize() == 0) return false; iteval=0; return true;}
        bool operator++(int){iteval++; return iteval < data.getSize();}
        uint32_t operator()()const{data.deref_key(iteval);}
        const C* operator->(){return &(data.deref(iteval));}
        const C& operator*(){return data.deref(iteval);}
	};
    Tuple<C,1u, Cflag>::Iterator getIterator(){return Tuple<C,1u, Cflag>::Iterator(*this);}
    Tuple<C,1u, Cflag>::ConstIterator getIterator()const{return Tuple<C,1u, Cflag>::ConstIterator(*this);}



};*/

template<class C, int sizex, int sizey> struct isTypeEquivalent< C*, TMatrix<C,sizex,sizey> > {enum {ans = true }; };
template<class C, int sizex, int sizey> struct isTypeEquivalent< TMatrix<C,sizex,sizey> , C*> {enum {ans = true }; };
template<class C, int sizex, int sizey, int size> struct isTypeEquivalent< Tuple<C,size> ,  TMatrix<C,sizex,sizey> > {enum {ans = true }; };
template<class C, int sizex, int sizey, int size> struct isTypeEquivalent< TMatrix<C,sizex,sizey> , Tuple<C,size> > {enum {ans = true }; };




#ifdef GNU_SCIENTIFIC_LIBRARY
class LogPvalSum{
public:
	unsigned int nb_independent;
	double value;

	LogPvalSum& toZero(){nb_independent =0; value=0.0f; return(*this);}
	void operator+=(const double val){if (ExOp::isValid(val)) {nb_independent++; value+= val;} }
	double operator()() const{ return (nb_independent > 1) ?  LogPvalue_SumLogPvalues_Ptail(value, nb_independent ) : value; }
	void compress(){if (nb_independent > 1) {value = LogPvalue_SumLogPvalues_Ptail(value, nb_independent ); nb_independent =1;} }
};
#endif


template<class C, unsigned int flag>
class GaussElem{
public:
    GaussElem();
    // Tuple<C,SIZE> getMean() const; inherited
};

template<class C, unsigned int SIZE>
//class GaussElem< Tuple<C, SIZE>, (ExCo<C>::IS_COMMUTATIVE::ans == 0) ? 0u : 1u >{ // COMMUTATIVE!
class GaussElem< Tuple<C, SIZE>, 0u >{ // COMMUTATIVE!
	public:
    double w,w2;
    Tuple<C, SIZE> mean;
    Trianglix<C, SIZE> cov;
    mutable C determinant; // set when needed
	static const int TARG2 =0;
	    typedef GaussElem< Tuple<C, SIZE>, 0u > SAFETYPE;
	static const bool IsPOD = Tuple<C, SIZE>::IsPOD && Trianglix<C, SIZE>::IsPOD;
	typedef Weight WEIGHT_TYPE;

    GaussElem(){ExOp::toZero(determinant);}

    GaussElem(const Tuple<WeightElem<C,2>, SIZE> &o, double _w=1.0f): w(_w), w2(0.0f){unsigned int i,j,k; for(i=0,k=0;i< SIZE;i++) {mean[i] = o[i].getMean(); for(j=0;j< i;j++) {cov.data[k++] =mean[i] * mean[j];} cov.data[k++] = o[i].getVar() + mean[i] * mean[i];}
		ExOp::toZero(determinant);
		}
    GaussElem(const Tuple<C, SIZE> &o, double _w=1.0f): w(_w), w2(_w*_w), mean(o * w), cov(o * sqrt(w)){ExOp::toZero(determinant); }
    GaussElem(const Tuple<C, SIZE> &_mean, const Trianglix<C, SIZE>& _cov, double _w, double _w2): w(_w), w2(_w2), mean(_mean), cov(_cov){ExOp::toZero(determinant);}

    // ExOp
    double getWeight()const{return w;}
	double getN()const{return (w == 0) ? 0.0f : w*w/w2;}
	inline unsigned int getSize() const{return SIZE;}
	void setSize(unsigned int size){mean.setSize(size); cov.setSize(size); }
 //   void show(FILE* f = stdout, int level =0)const;

	GaussElem< Tuple<C, SIZE>, 0 >& operator+=(const Tuple<C, SIZE>& shift){unsigned int i,j,k; for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) cov.data[k] += (((ExOp::mkTrju(shift[i]) * shift[j]) * w) + ExOp::mkTrju(shift[i]) * mean[j]+ ExOp::mkTrju(mean[i]) * shift[j]); mean += shift * w; return(*this);}
	GaussElem< Tuple<C, SIZE>, 0 >& operator-=(const Tuple<C, SIZE>& shift){unsigned int i,j,k; for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) cov.data[k] += (((ExOp::mkTrju(shift[i]) * shift[j]) * w) - ExOp::mkTrju(shift[i]) * mean[j]- ExOp::mkTrju(mean[i]) * shift[j]); mean -= shift * w; return(*this);}

	Tuple<C,SIZE> getMean() const {	return( mean / w );}
    Tuple<C,SIZE> getVar() const {Tuple<C,SIZE> fout;  double fact = (1.0f / (w*w - w2)); unsigned int i,j; for(i=0,j=0;j<SIZE;i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}
	Tuple<C,SIZE> getVar_biaised() const {Tuple<C,SIZE> fout;  double fact = (1.0f / (w*w)); unsigned int i,j; for(i=0,j=0;j<SIZE;i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}

    GaussElem< Tuple<C, SIZE>, 0 > & toZero(){w=0.0f;w2=0.0f;ExOp::toZero(mean); ExOp::toZero(cov);ExOp::toZero(determinant);return(*this);}
    GaussElem< Tuple<C, SIZE>, 0 > & toRand(){w=1.0f;w2=1.0f;ExOp::toRand(mean); cov = Trianglix<C, SIZE>(mean);ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, SIZE>, 0 > & operator=(const GaussElem< Tuple<C, SIZE>,TARG2 > &other){mean = other.mean;cov = other.cov;w = other.w;w2 = other.w2; determinant = other.determinant;return(*this);}


	GaussElem< Tuple<C, SIZE>, TARG2 >& operator+=(const GaussElem< Tuple<C, SIZE>,TARG2 > &other){w+=other.w; w2+=other.w2; mean += other.mean; cov += other.cov; ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, SIZE>, TARG2 > operator+(const GaussElem< Tuple<C, SIZE>, TARG2 > &other)const {return (GaussElem< Tuple<C, SIZE>, TARG2 >(mean +other.mean,  cov + other.cov, w+other.w, w2+other.w2));}

	// free mean addition
	GaussElem< Tuple<C, SIZE>, 0 >& addGaussElem_free_mean(const GaussElem< Tuple<C, SIZE>, TARG2 > &other);

	// REVERT SCOPE!
	GaussElem< Tuple<C, SIZE>, TARG2 >& operator-=(const GaussElem< Tuple<C, SIZE>,TARG2 > &other){w-=other.w; w2-=other.w2; mean -= other.mean; cov -= other.cov; ExOp::toZero(determinant); return(*this);}

	GaussElem< Tuple<C, SIZE>, 0u >& operator*=(const C& value){ mean *= value; cov *= ExOp::mkPowInt(value,2); return (*this);}
	GaussElem< Tuple<C, SIZE>, 0u >& operator*=(const Tuple<C,SIZE>& value){ mean *= value; unsigned int k=0; for(unsigned int i=0;i<SIZE;i++) for(unsigned int j=0;j<=i;j++) cov.data[k++] *= value[i] * value[j]; return (*this);}

	GaussElem< Tuple<C, SIZE>, 0u >& operator*=(const Weight& inc_w){w *= inc_w();w2 *= inc_w()*inc_w(); mean *= inc_w(); cov *= inc_w(); return (*this);}
	GaussElem< Tuple<C, SIZE>, 0u > operator*(const Weight& inc_w) const {return GaussElem< Tuple<C, SIZE>, 0u >(mean * inc_w(), cov * inc_w(), w * inc_w(),w2 * inc_w()*inc_w());}
	GaussElem< Tuple<C, SIZE>, 0u >& operator/=(const Weight& inc_w){w /= inc_w();w2 /= inc_w()*inc_w(); mean /= inc_w(); cov /= inc_w(); return (*this);}
	GaussElem< Tuple<C, SIZE>, 0u > operator/(const Weight& inc_w) const {return GaussElem< Tuple<C, SIZE>, 0u >(mean / inc_w(), cov / inc_w(), w / inc_w(),w2 / (inc_w()*inc_w()));}

	GaussElem< Tuple<C, SIZE>, 0u >& toHouseHolderMultiply(const Tuple<C, SIZE> &vec, const C &denum2);
	GaussElem< Tuple<C, SIZE>, 0u > mkHouseHolderMultiply(const Tuple<C, SIZE> &vec, const C &denum2) const{GaussElem< Tuple<C, SIZE>, 0u > fout(*this); fout.toHouseHolderMultiply(vec,denum2); return(fout);}


LFH_GOLD double likelihoodratio_dist(const GaussElem< Tuple<C, SIZE>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const;
    double bhattacharryya_dist(const GaussElem< Tuple<C, SIZE>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const; // weighted does not work


LFH_GOLD double LLikelihood(bool has_covariance =true) const{return (w == 0.0f)? 0.0f : Entropy(has_covariance) * -w;}
LFH_GOLD double Entropy(bool has_covariance =true ) const;

    double cmpLLikelihoodOf(const Tuple<C, SIZE> &what, bool has_covariance = true)const;
    double cmpLLikelihoodOf(const GaussElem<Tuple<C, SIZE>, 0u > &what, bool has_covariance = true) const;

    void setWeight(double _w){ mean *= (_w / w); cov *= (_w / w); w2 *= (_w / w); w = _w;ExOp::toZero(determinant);}

    Trianglix<C, SIZE> getCovariance()const;
	Trianglix<C, SIZE> getCovariance_biased()const;
	Trianglix<C, SIZE> getCorrelation()const;
	template<class D> void applyPartialCounts(const Trianglix<D, SIZE> &partialCounts, const TMatrix<C,SIZE,SIZE> &partialSums);
	template<class D, unsigned int O_SIZE> GaussElem< Tuple<C, O_SIZE>, 0u > operator*(const TMatrix<D,SIZE , O_SIZE> &) const;
	template<class D> GaussElem< Tuple<C, SIZE>, 0u > operator*(const Matrix<D> &);

    void setCovariance(const Trianglix<C, SIZE>& covar);
    void setCovariance(const GaussElem< Tuple<C, SIZE>, 0 >& other){setCovariance(other.getCovariance());}

//	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
    void show(FILE *f = stdout, int level=0)const{fprintf(f,"weight %e (%e) mean:\t", w ,w2); ExOp::show(mean / w,f,1); fprintf(f,"\nCovar\n"); ExOp::show(getCovariance(),f,0);  }

	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeSubGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;
	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeResidualGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;

#ifdef GNU_SCIENTIFIC_LIBRARY
	double PvalueOf(const GaussElem< Tuple<C, SIZE>, 0u > &x)const;
	LFH_GOLD double Pvalue_Hotelling(const GaussElem< Tuple<C, SIZE>, 0u > &x, bool LogPvalue = false) const;
	template<unsigned int S_SIZE> double PvalueOfResidual(const GaussElem< Tuple<C, SIZE>, 0u > &x , const Tuple<unsigned int, S_SIZE>& selected_dims)const; // TODO!
#endif
	void setMeanVaronDimention(double mean, double var, unsigned int dim);

	ERRCODE load(FILE *f, unsigned int size=0);
	ERRCODE save(FILE *f) const;
};

class UnknownScope{
   double outlier_LL;
public:
   double scp[8];
   double relerr;
   double antioss;
   double oldLL;


   UnknownScope(double aim);

   bool EMprocess(double* LLs, double &total_LL, uint32_t length, Tuple<double,2u> *outlier_prior = NULL); // replaces LLs by weights

   double EMregist(double LL); // return non-outlier weight! (0.0-1.0)
   bool EMfinit(double &LL); // return success/failure, and output total LL
};

template<class C>
class GaussElem< Tuple<C, 0u>, 0u >{ // COMMUTATIVE!
public:
    double w,w2;
    Tuple<C, 0u> mean;
    Trianglix<C, 0u> cov;
	typedef Weight WEIGHT_TYPE;
    mutable C determinant; // set when needed
	static const int TARG2 =0;
	typedef GaussElem< Tuple<C, 0u>, 0u > SAFETYPE;
    GaussElem(){ExOp::toZero(determinant);}

    GaussElem(const Tuple<WeightElem<C,2>, 0u> &o, double _w=1.0f): w(_w), w2(0.0f){
        this->setSize(o.getSize());
        unsigned int i,j,k; for(i=0,k=0;i< mean.getSize();i++) {mean[i] = o[i].getMean(); for(j=0;j< i;j++) {cov.data[k++] =mean[i] * mean[j];} cov.data[k++] = o[i].getVar() + mean[i] * mean[i];}
		ExOp::toZero(determinant);
	}
    GaussElem(const Tuple<C, 0u> &o, double _w=1.0f): w(_w), w2(_w*_w), mean(o * w), cov(o * sqrt(w)){ExOp::toZero(determinant); }
    GaussElem(const Tuple<C, 0u> &_mean, const Trianglix<C, 0u>& _cov, double _w, double _w2): w(_w), w2(_w2), mean(_mean), cov(_cov){ExOp::toZero(determinant);}

    // ExOp
    double getWeight()const{return w;}
	double getN()const{return (w == 0) ? 0.0f : w*w/w2;}

	void setSize(unsigned int size){mean.setSize(size); cov.setSize(size); }
	inline unsigned int getSize() const{return mean.getSize();}
	//   void show(FILE* f = stdout, int level =0)const;
	Tuple<C,0u> getMean() const {return( mean / w );}
    Tuple<C,0u> getVar() const {Tuple<C,0u> fout;  fout.setSize(mean.getSize()); double fact = (1.0f / (w*w - w2)); unsigned int i,j; for(i=0,j=0;j<mean.getSize();i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}
	Tuple<C,0u> getVar_biaised() const {Tuple<C,0u> fout; fout.setSize(mean.getSize()); double fact = (1.0f / (w*w)); unsigned int i,j; for(i=0,j=0;j<mean.getSize();i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}

    GaussElem< Tuple<C, 0u>, 0 > & toZero(){w=0.0f;w2=0.0f;ExOp::toZero(mean); ExOp::toZero(cov);ExOp::toZero(determinant);return(*this);}
    GaussElem< Tuple<C, 0u>, 0 > & toRand(){w=1.0f;w2=1.0f;ExOp::toRand(mean); cov = Trianglix<C, 0u>(mean);ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, 0u>, 0 > & operator=(const GaussElem< Tuple<C, 0u>,TARG2 > &other){mean = other.mean;cov = other.cov;w = other.w;w2 = other.w2; determinant = other.determinant;return(*this);}
    GaussElem< Tuple<C, 0u>, 0 > & toMemmove(GaussElem< Tuple<C, 0u>, 0 > &other){w = other.w; w2= other.w2; mean.toMemmove(other.mean); cov.toMemmove(other.cov); determinant = other.determinant; return *this;}

	GaussElem< Tuple<C, 0u>, TARG2 >& operator+=(const GaussElem< Tuple<C, 0u>,TARG2 > &other){
//		printf("time to crash %i and %i!\n", mean.getSize(), other.mean.getSize()); fflush(stdout);
	w+=other.w; w2+=other.w2; mean += other.mean; cov += other.cov; ExOp::toZero(determinant);
//		printf("time to crash done!\n"); fflush(stdout);
		 return(*this);
	}
    GaussElem< Tuple<C, 0u>, TARG2 > operator+(const GaussElem< Tuple<C, 0u>, TARG2 > &other)const {return (GaussElem< Tuple<C, 0u>, TARG2 >(mean +other.mean,  cov + other.cov, w+other.w, w2+other.w2));}
	// free mean addition
	GaussElem< Tuple<C, 0u>, 0 >& addGaussElem_free_mean(const GaussElem< Tuple<C, 0u>, TARG2 > &other);

	GaussElem< Tuple<C, 0u>, 0 >& operator+=(const Tuple<C, 0u>& shift){unsigned int i,j,k; for(j=0,k=0;j<mean.getSize();j++) for(i=0;i<=j;i++,k++) cov.data[k] += (((ExOp::mkTrju(shift[i]) * shift[j]) * w) + ExOp::mkTrju(shift[i]) * mean[j]+ ExOp::mkTrju(mean[i]) * shift[j]); mean += shift * w; return(*this);}
	GaussElem< Tuple<C, 0u>, 0 >& operator-=(const Tuple<C, 0u>& shift){unsigned int i,j,k; for(j=0,k=0;j<mean.getSize();j++) for(i=0;i<=j;i++,k++) cov.data[k] += (((ExOp::mkTrju(shift[i]) * shift[j]) * w) - ExOp::mkTrju(shift[i]) * mean[j]- ExOp::mkTrju(mean[i]) * shift[j]); mean -= shift * w; return(*this);}

	// REVERT SCOPE!
	GaussElem< Tuple<C, 0u>, TARG2 >& operator-=(const GaussElem< Tuple<C, 0u>,TARG2 > &other){w-=other.w; w2-=other.w2; mean -= other.mean; cov -= other.cov; ExOp::toZero(determinant); return(*this);}

	GaussElem< Tuple<C, 0u>, 0u >& operator*=(const C& value){ mean *= value; cov *= ExOp::mkPowInt(value,2); return (*this);}

	GaussElem< Tuple<C, 0u>, 0u >& operator*=(const Tuple<C, 0u>& value){ mean *= value; unsigned int k=0; for(unsigned int i=0;i<mean.tup_size;i++) for(unsigned int j=0;j<=i;j++) cov.data[k++] *= value[i] * value[j]; return (*this);}


	GaussElem< Tuple<C, 0u>, 0u >& operator*=(const Weight& inc_w){w *= inc_w();w2 *= inc_w()*inc_w(); mean *= inc_w(); cov *= inc_w(); return (*this);}
	GaussElem< Tuple<C, 0u>, 0u > operator*(const Weight& inc_w) const {return GaussElem< Tuple<C, 0u>, 0u >(mean * inc_w(), cov * inc_w(), w * inc_w(),w2 * inc_w()*inc_w());}
	GaussElem< Tuple<C, 0u>, 0u >& operator/=(const Weight& inc_w){w /= inc_w();w2 /= inc_w()*inc_w(); mean /= inc_w(); cov /= inc_w(); return (*this);}
	GaussElem< Tuple<C, 0u>, 0u > operator/(const Weight& inc_w) const {return GaussElem< Tuple<C, 0u>, 0u >(mean / inc_w(), cov / inc_w(), w / inc_w(),w2 /(inc_w()*inc_w()));}

	GaussElem< Tuple<C, 0u>, 0u >& toHouseHolderMultiply(const Tuple<C, 0u> &vec, const C &denum2);
	GaussElem< Tuple<C, 0u>, 0u > mkHouseHolderMultiply(const Tuple<C, 0u> &vec, const C &denum2) const{GaussElem< Tuple<C, 0u>, 0u > fout(*this); fout.toHouseHolderMultiply(vec,denum2); return(fout);}

	GaussElem< Tuple<C, 0u>, 0u >& operator+=(const KeyElem<double, SparseTuple<C> > &vec);

    double likelihoodratio_dist_underbound(const GaussElem< Tuple<C, 0u>, TARG2 > &other, unsigned int nb_dims, bool has_weight = true, bool has_covariance = true, C const * const this_det = NULL, C const * const  other_det = NULL) const;

	LFH_GOLD double likelihoodratio_dist(const GaussElem< Tuple<C, 0u>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const;
    double bhattacharryya_dist(const GaussElem< Tuple<C, 0u>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const;


	LFH_GOLD double LLikelihood(bool has_covariance =true ) const{return (w == 0.0f)? 0.0f : Entropy(has_covariance) * -w;}
	double cmpLLikelihoodOf(const Tuple<C, 0u> &what, bool has_covariance = true)const;
    double cmpLLikelihoodOf(const GaussElem<Tuple<C, 0u>, 0u > &what, bool has_covariance = true) const;

	LFH_GOLD double Entropy(bool has_covariance =true ) const;

    void setWeight(double _w){ mean *= (_w / w); cov *= (_w / w); w2 *= (_w / w); w = _w;ExOp::toZero(determinant);}

    Trianglix<C, 0u> getCovariance()const;
	Trianglix<C, 0u> getCovariance_biased()const;
	Trianglix<C, 0u> getCorrelation()const;
	template<class D> void applyPartialCounts(const Trianglix<D, 0u> &partialCounts, const TMatrix<C> &partialSums);

	template<class D> GaussElem< Tuple<C, 0u>, 0u > operator*(const Matrix<D> &);

    void setCovariance(const Trianglix<C, 0u>& covar);
    void setCovariance(const GaussElem< Tuple<C, 0u>, 0 >& other){setCovariance(other.getCovariance());}

	//	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
    void show(FILE *f = stdout, int level=0)const{fprintf(f,"mean:\t"); ExOp::show(mean / w,f,1); fprintf(f,"\nCovar\n"); ExOp::show(getCovariance(),f,0);  }

    GaussElem< Tuple<C, 0u>, 0u > mkSubGauss_first(unsigned int nbdim)const;

	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeSubGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;
	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeResidualGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;

	double PvalueOf(const GaussElem< Tuple<C, 0u>, 0u > &x)const;
	template<unsigned int S_SIZE> double PvalueOfResidual(const GaussElem< Tuple<C, 0u>, 0u > &x , const Tuple<unsigned int, S_SIZE>& selected_dims)const; // TODO!



	double getLogFoldinVariance_inMeanMajorAxis(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	double getLogFoldinVariance_Determinant(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	double getLogFoldinVariance_Diagonal(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	double getDifferenceinMean_inVarianceMajorAxis(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	#ifdef GNU_SCIENTIFIC_LIBRARY
	LFH_GOLD double Pvalue_Hotelling(const GaussElem< Tuple<C, 0u>, 0u > &x, bool LogPvalue = false) const;
	double Pvalue_LRTest_UnequalMean(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL,  bool LogPvalue = false) const; // Ho = mu_a = mu_b V_a = V_b   H1 = mu_a != mu_b V_a = V_b
	double Pvalue_LRTest_UnequalMeanandVariance(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL, bool LogPvalue = false) const; // Ho = mu_a = mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b
	double Pvalue_LRTest_UnequalVariance(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL, bool LogPvalue = false) const; // Ho = mu_a != mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b

	double Pvalue_for_sum(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL, bool LogPvalue = false);

	#endif
//	double Pvalue_LRTest_UnequalVariance_Masked(const GaussElem< Tuple<C, 0u>, 0u > &x, const Trianglix<bool> &mask, bool LogPvalue = false) const; // Ho = mu_a != mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b

	ERRCODE load(FILE *f, unsigned int size=0);
	ERRCODE save(FILE *f) const;

	Tuple<C> doImputation(const SparseTuple<C> &input, const Tuple<uint32_t> &missingindexes)const;
	Tuple<C> doImputationCentered(const SparseTuple<C> &input, const Tuple<uint32_t> &missingindexes)const; // relative to mean
	void addImputatedSampleCovar(Trianglix<double> &target, const SparseTuple<C> &input, const Tuple<uint32_t> &missingindexes)const; // relative to mean

	GaussElem<Tuple<C, 0u> > mkSubGaussElem(const Tuple<bool> filter)const;
	void ldFrom(const DataGrid<double, 2u> &data);

    // maximize likelihood of observation, where missing values are scatered around, uses EM to alternatively impute and maximize
	class TaskPartialObservationMaximize : public Event<uint32_t>{
        public:
        Tuple< Tuple<double,0u> > stats;
        //Tuple< GaussElem< Tuple<C, 0u>, 0u > > suff_stats;
        const SparseMatrix<double>& data;
        Tuple< Tuple<double> > imputed;
        Tuple<Tuple<uint32_t> > missIndexes;
        Tuple<double> missWeight;
        double baseLL;
        Tuple<double> ll_or_prob;


        Trianglix<double> precision_matrix;

        const GaussElem< Tuple<C, 0u>, 0u > &cur;
        const double* weight;
        ProgressBarPrint progb;
        bool include_uncertainty;

        Tuple<uint32_t> colrange;
        uint32_t taskID;
        uint32_t aprogr;
        double aprogr_factor;
        uint32_t thrprint;
        TaskPartialObservationMaximize(const SparseMatrix<double>& data, const GaussElem< Tuple<C, 0u>, 0u > &cur, double outlier_frequency, uint32_t nbthread = 4u, bool include_uncertainty = true);
        uint32_t operator()(uint32_t threadID);
	};
    template<class TB> void runTaskPartialObservationMaximize(TB &tb, const SparseMatrix<double>& data, uint32_t nbrow, uint32_t itemax, double outlier_frequency, Tuple<double,2u> *outlier_prior = NULL, bool include_uncertainty = true, Vector<double>* weight=NULL);

    class runTaskPartialObservationMaximizeMK2_task;
    void runTaskPartialObservationMaximizeMK2(const SparseMatrix<double>& data, uint32_t nbrow, uint32_t itemax, double outlier_frequency, Tuple<double,2u> *outlier_prior = NULL, bool include_uncertainty = true, Vector<double>* weight=NULL);



	class TaskStochasticMaximize : public Event<uint32_t>{
        public:
        ThreadBase &tb;
        const GaussElem< Tuple<C, 0u>, 0u > &cur;

        Tuple< Tuple<double,0u> > stats;
        //Tuple< GaussElem< Tuple<C, 0u>, 0u > > suff_stats;
        const SparseMatrix<double>& data;
        Tuple< Tuple<double> > imputed;
        Tuple<Tuple<uint32_t> > missIndexes;
        Tuple<double> missWeight;
        double baseLL;
        Tuple<double> ll_or_prob;

        Trianglix<double> precision_matrix;


        const double* weight;
        bool include_uncertainty;

        Tuple<uint32_t> colrange;
        uint32_t taskID;
        template<class F> TaskStochasticMaximize(ThreadBase &_tb, GaussElem< Tuple<C, 0u>, 0u >& host, const SparseMatrix<double>& data, F f, uint32_t nbrow, uint32_t maxite, double outlier_frequency,Tuple<double,2u> *outlier_prior, bool include_uncertainty = true, Vector<double>* weight = NULL);
        uint32_t operator()(uint32_t threadID);
	};
    template<class F> void runTaskStochasticMaximize(ThreadBase &tb, const SparseMatrix<double>& data, F f, uint32_t nbrow, uint32_t itemax, double outlier_frequency, Tuple<double,2u> *outlier_prior = NULL, bool include_uncertainty = true, Vector<double>* weight=NULL,
       FUNCTOREQUIRES_DECL(F,(Tuple<double, 2u>), (const Tuple<uint32_t,2u>&)));
    //                                                 double (F::*_f)(const Tuple<uint32_t, 2u>&) const = F::operator());
//    template<class F> double willitwork(const F &f  );
  //  double willitwork2(const function<double(const Tuple<unsigned int,2u>&)> &f);
 };



template<class C, unsigned int SIZE>
//class GaussElem< Tuple<C, SIZE>, (ExCo<C>::IS_COMMUTATIVE::ans == 0) ? 1u : 0u >{ // NON-COMMUTATIVE!
class GaussElem< Tuple<C, SIZE>, 1u >{ // NON-COMMUTATIVE!
public:
    double w,w2;
    Tuple<C, SIZE> mean;
	static const int TARG2 =1;
	typedef Weight WEIGHT_TYPE;
	typedef GaussElem< Tuple<C, SIZE>, 0u > SAFETYPE;
   // void show(FILE* f = stdout, int level =0)const;
	Tuple<C,SIZE> getMean() const {	return( (C)(mean / w) );}
//	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
    void show(FILE *f = stdout, int level=0)const{fprintf(f,"dont commute!\n");}
};
// Assumes C are commutative!

// Gradient Search scope stores the last given derivative of a Function on which gradient Ascent or Descent is Performed in order to check an bound on the non-linearity of the function.
class GradientSearchScope{
	int size;
	double* lastderiv;
	double lastvalue;
	double expdiff;
	double log_alpha;
public:
	GradientSearchScope(): lastderiv(NULL){}
	GradientSearchScope(const GradientSearchScope&);
	~GradientSearchScope() {delete[](lastderiv);}
	GradientSearchScope& operator=(const GradientSearchScope&);

	double operator()() const;

	void init(double initalalpha, int in_size);
	void registAscent(const double &value, double const * const deriv, double mod = 1.0f);
	void registDescent(const double &value, double const * const deriv, double mod = 1.0f);

	double updateAscent(const double &value, double * guess,  double const * const deriv); // returns a safe bound on the parameter changes;
	double updateDescent(const double &value, double * guess , double const * const deriv); // returns a safe bound on the parameter changes;

	void punish(double log_mag);

};

// allows entries with missing data
template<class C>
class PartialGaussElem{
    public:
    Trianglix< Tuple<C, 4u> > buffer;

    void setSize(int s);
    void getSize() const{return buffer.getSize();}
    void addEntry(const myHashmap<uint32_t, C>&, double weight = 1.0f);
    void addNonzero(const Tuple<C>& data, double weight = 1.0f);

    PartialGaussElem<C>& toAddFiltCol(const myHashmap<uint32_t,C>& data, const uint32_t* rows, double weight = 1.0f);

    PartialGaussElem<C> operator+(const PartialGaussElem&other) const;
    PartialGaussElem<C>& operator+=(const PartialGaussElem&other);
    void wrMean(Tuple<C> &fout) const;
    void wrVar(Tuple<C>  &fout) const;
    void wrVar_unbiaised(Tuple<C>  &fout) const;
    void wrVar(Trianglix<C> &fout) const;
};

class CurvatureSearchScope{
	unsigned int size;
	double* lastderiv; // (lastderiv, lastguess, tmpbuffer)
	double lastvalue;
	double expdiff;
	Tuple<double> curvature;
	Trianglix<double> curv;
	uint32_t alpha;
    unsigned int shrinkcount;
public:
    double epsilon;

	CurvatureSearchScope(): lastderiv(NULL){}
	CurvatureSearchScope(const CurvatureSearchScope&);
	~CurvatureSearchScope() {delete[](lastderiv);}
	CurvatureSearchScope& operator=(const CurvatureSearchScope&);

	void init(unsigned int in_size, double initalalpha, double _epsilon = pow(0.5, 53));
	void init(unsigned int in_size, const double * dderiv, double _epsilon = pow(0.5, 53)); // uses d2F/dx2 only!
	void init(unsigned int in_size, const Trianglix<double> &dderiv, double _epsilon = pow(0.5, 53));

	void initMK2(unsigned int in_size, const double * dderiv); // uses d2F/dx2 only!
    unsigned int  getSize()const{return size;}

    void updateRandom(double * guess); // needs at least 1 normal step done previously
	void updateGotNAN(double * guess); // returns a safe bound on the parameter changes;
	bool updateAscent(const double &value, double * guess,  double const * const deriv); // returns a safe bound on the parameter changes;
	double updateAscentVerbose(const double &value, double * guess,  double const * const deriv, bool debug); // returns a safe bound on the parameter changes;
	double updateDescent(const double &value, double * guess , double const * const deriv); // returns a safe bound on the parameter changes;

	double getRelErrorBound() const;
//	double updateAscentMK2(const double &value, double * guess,  double const * const deriv,  double const * const dbl_deriv); // returns a safe bound on the parameter changes;
	//double updateDescent(const double &value, double * guess , double const * const deriv,  double const * const dbl_deriv); // returns a safe bound on the parameter changes;

	double checkDerivative(const double &value, double * guess , double const * const deriv); // test if input derivative are correct

	void wrFinalGuess(double* guess)const;
    double getLastValue()const{return lastvalue;}
    void show(FILE* f = stdout, int level=0) const{fprintf(f,"Curvature scope in %i dimentions, last value %e\n",size,lastvalue);}
};

class FuzzyLineSearch{
public:
    double guess[3];
    double value[3];
    double initGuesses(double cur_best, double newguess);
    double solveArgmax(double& cur_best, double &nextguess) const; // find the best of the 3 guess, and return next guess
};

template<class C>
class FLineSearch{
public:
    C guess[3];
    C nextguess;
    double value[3];
    double initGuesses(C cur_best, double newguess);
    unsigned int solveArgmax(C& cur_best, double &nextguess) const; // find the best of the 3 guess, and return next guess
};

// when the number of parameters is large, but can be projected into a small number of directions
template<unsigned int DIMS>
class FuzzyCurvatureSearchScope{
public:
    Tuple<double, DIMS> alpha;
    double best;
    double operator[](int i){return alpha[i];}
    bool needsUndo(double val, bool is_first_time) const;
    void registUndo(double scale);
    bool regist(double val, double* derivmag, bool isInit = false);
};

enum LFHDOMAINMAP_enum{
    LFHDOMAINMAP_POSITIVE_REAL,
    LFHDOMAINMAP_ZERO_TO_ONE
};
template<class C, unsigned int DIMS, LFHDOMAINMAP_enum DOM>
class DomainMappedTuple : public Tuple<C, DIMS> {
    public:
    Tuple<C, DIMS> &target;
    DomainMappedTuple(Tuple<C, DIMS>& _target): target(_target){this->setSize(target.getSize());}
    Tuple<C, DIMS> eval()const{}
    void applyUpdate(Tuple<C, DIMS>& deriv, double alpha){}
};

// TODO!?
template<class C, unsigned int DIMS>
class DomainMappedTuple<C,DIMS,LFHDOMAINMAP_POSITIVE_REAL> : public Tuple<C, DIMS> {
    public:
    Tuple<C, DIMS> &target;
    DomainMappedTuple(Tuple<C, DIMS>& _target);
    Tuple<C, DIMS> eval()const;
    void applyUpdate(Tuple<C, DIMS>& deriv, double alpha);
};

template<class C, unsigned int DIMS>
class DomainMappedTuple<C,DIMS,LFHDOMAINMAP_ZERO_TO_ONE> : public Tuple<C, DIMS> {
    public:
    Tuple<C, DIMS> &target;
    DomainMappedTuple(Tuple<C, DIMS>& _target);
    Tuple<C, DIMS> eval()const;
    void applyUpdate(Tuple<C, DIMS>& deriv, double alpha);
};

/*
template<class C, class D>
class MaskMapping{
public:
    Tuple<uint32_t> head; //
    uint32_t keylength;
    Tuple<Vector<uint32_t> > stript_next;
    myHashmap<C, D> data;
    void insert(const C& key, const D& value){}
    void getIntersection(Vector<D> &res, C value, C mask) const{} // return all entries with  ((value ^ key) & mask) == 0
};*/

/*
template<class Key, class Data, class HashFnc>
class AsyncHashmap_handdle{
public:
    myHashmap<Key, Data, HashFnc>* target;
    myHashmap<Key, Data, HashFnc>* claimHashmap();
    void releaseHashmap();

};

template<class Key, class Data, class HashFnc>
class AsyncHashmap{
    AsyncStructure<myHashmap<Key, Data, HashFnc> > hashmap;
public:


};*/

template<class N, class P>
class SetPartition{
public:
    Tuple<uint32_t,2u> partitionless; // First Node iterator, Size of partition
    myHashmap<P, Tuple<uint32_t,2u> > partitions; // Partition, First Node iterator, Size of partition
    myHashmap<N, Tuple<uint32_t,3u> > nodes; // Node, partition iterator, prev/next elem iterator
    void addElem(const N&); // mk partitionless
    void addElem(const N&, const P&);
//    template<class COMP> initFromComparator(const COMP& comp, uint32_t nbelem);
};

class TableReader{ // Table read *must* have a header row, and '+' '-' '*' '/' '(' ')' are treated as special characters, no column name should start with a digit too
	unsigned int* read_info;
	unsigned int nbread;
	Vector<unsigned int> oper;
public:
	FILE *f;
	unsigned int nbcols;
	Vector<double> constants;

	string* current_row;
	string* col_name;
	TableReader(const char*, const Vector<const char*> &columns); // columns that are integer representation (only digits) are used as column offsets
	~TableReader();
	const char* operator[](const unsigned int i ) const{return current_row[i].c_str();}
	bool nextRow();
};

// store serial elements, assumes that store elements
template<class IDTYPE>
class SerialStore{
	void compress(char* buffer);
public:
	myHashmap<IDTYPE, KeyElem<unsigned int, unsigned int> > hash_scope;
	FILE *f;
	char *path;
	unsigned int wasted;
	unsigned int size;
    //unsigned int nb_readonly;
    FILE* readonly_handdle;

    SerialStore();
	SerialStore(const char * const f_path, bool write_only = false);
	~SerialStore();

    void setTarget(const char * const f_path, bool write_only = false); // , bool isInitialized = true

	template<class T> ERRCODE load(IDTYPE, T &);
	template<class T> ERRCODE save(IDTYPE, const T &);


    unsigned int itemSize(IDTYPE) const;
	template<class T> ERRCODE load_arrayitem(IDTYPE, T &, unsigned int posis, unsigned int unitsize = sizeof(T)) const;
	template<class T> ERRCODE save_arrayitem(IDTYPE, const T &, unsigned int posis, unsigned int unitsize = sizeof(T));
	void save_reservearray(IDTYPE, unsigned int size);

    FILE* getItemHandle(IDTYPE) const; // chunk size *must* be defined on advance for writing (using save_reservearray)
    FILE* getItemHandle(IDTYPE, unsigned int channel) const;
    void getItemHandle(IDTYPE, FILE*& filehanddle) const;


	void flush(char* newpath = NULL); //
	bool has(IDTYPE) const;

    void reload(); // reloads file *without saving modifications, if any*

    ERRCODE loadFrom(const char* const path);
    ERRCODE saveTo(const char* const path) const;



    void setNbReadOnly(unsigned int nb);
	void show(FILE* f = stdout) const;
};

template<class IDTYPE>
class SerialStore_readHanddle{
    public:
    SerialStore<IDTYPE>& target;
    FILE* f;
    SerialStore_readHanddle(SerialStore<IDTYPE>&);
    ~SerialStore_readHanddle();
    bool isValid()const{return f != NULL;}
};


// Scope for batch Fourier Transforms, on arrays of ANY size
class Bluestein{
public:
	unsigned char pre_pow2;
	unsigned int mult;
	unsigned char post_pow2;
	Vector< Complex<double> > bluewindow;
	unsigned int tsize;

	Bluestein(){}
	Bluestein(unsigned int size){setSize(size);}
	void setSize(unsigned int n_size);

	unsigned int getBufferSize()const;

	template<class C> void toFT(C* data) const;
	template<class C> void toIFT(C* data) const;

	template<class C> void toFT_even(C* data) const;
	template<class C> void toIFT_even(C* data) const;

	template<class C, class D> void toFT(D* o_data,const C* i_data) const;
	template<class C, class D> void toIFT(D* o_data,const C* i_data) const;

	template<class C> void alloc_resample_buffer(C* &iobuf, Bluestein* invblue) const;
	template<class C> void resample(C* resample_buffer, Bluestein* invblue) const;

	void show(FILE* f = stdout, int level =0)const;
};

enum CHART_CACHEC_FRAMES_enum{
	CHART_CACHEC_FRAMES_POINT_DENS=0
};


template<class C = unsigned char,unsigned int NBCHAN= 3u>
class CharArea{
	public:
	DataGrid<C ,3> image;
	Tuple<unsigned int , 4> inner_rect;
	Tuple<C, NBCHAN> color_bg;
	Tuple<C, NBCHAN> color_axe;
	Tuple<C, NBCHAN> color_axetext;

	Tuple<C, NBCHAN> color_min;
	Tuple<C, NBCHAN> color_max;

	myHashmap<uint32_t , void*> stored; // key is a CHART_CACHEC_FRAMES_enum

	void* getCached();

	void setDefaultStyle();
	void initialize(const Tuple<unsigned int,2> dims);
	void setAxes(double Xmin,double Ymin, double Xmax, double Ymax, bool Xlog= false, bool Ylog = false);

	template<class D, class E> void drawFrame(const DataGrid<D, 2> &fr,const E &min,const E &max);

	void overwriteFrame(const DataGrid<C, 3> &fr);

//	template<class D> DataGrid<C ,3> makeChart(const DataGrid<D,2> &image)const;

};


// structure store a list of strings.
// structure maintain a scope for substring searches
class Dictionary{
	uint32_t findWord(const char* word, int word_lenght) const;
	uint32_t insertWord(const char* word, int word_lenght);
	Tuple<uint32_t, 4u> findPrefix(const char* prefix, int word_lenght) const;

	uint32_t flags;
	public:
	Vector<char*> entries;
	Vector< Vector<uint32_t> > word_links;

	//myHashmap<uint32_t, uint32_t > mainhash;
	Vector< myHashmap<uint32_t, uint32_t> > hashes;
	Dictionary();
	~Dictionary();

	void setIgnoreCaps(bool ignore){ flags = (flags & 0xFFFFFFFE) | (ignore? 1: 0);}
	void addEntry(const char* str, bool allow_duplicate = true);

	void registerSubstrings(int maxlen); // maxlenght a multiple of 4 is expected, rounds up otherwise
	void findSubs(Vector<const char*> &fout, const char* query, uint32_t max_requested = 16u);
	void findIndexes(Vector<uint32_t> &fout, const char* query, uint32_t max_requested = 16u);
//	template<class C> findElements(Vector<C> &fout, const C* heap, const char* query, int max_requested = 16);
    int getSize() const{return entries.getSize();}

    uint32_t findEntry(const char* query);
    Vector<uint32_t> findEntries(const vector<string> &queries);

    uint32_t find(const char* what) const;

    Dictionary& toMemfree();
    Dictionary& toMemmove(Dictionary& source);


    void show(FILE*f =stdout, int level=0) const;
    ERRCODE save(FILE*f) const;
    ERRCODE load(FILE*f, unsigned int s=0);
};





template<class C>
class DicoElem{
	public:
	Dictionary base;
	Vector<C> data;
	DicoElem();
//	~dictionary();

	void setIgnoreCaps(bool ignore){base.setIgnoreCaps(ignore);}
	void addEntry(const char* str, const C &data);
	void findSubs(Vector<C> &fout, const char* query, int max_requested = 16);
    void show(FILE*f =stdout, int level=0) const;
};

class Annotation;

template<int HASANNOT>
class AnnotMaintenance{
public:
    void toMemmove(void*, void*) const;
    void save(void*, FILE* f) const;
    void load(void*, FILE* f) const;
    const Tuple<uint32_t>* hasAnnot(void*) const;
    Annotation* canAnnot() const;
};

template< >
class AnnotMaintenance<0>{
public:
    inline void toMemmove(void*, void*) const{}
    inline void save(void*, FILE* f) const{}
    inline void load(void*, FILE* f) const{}
    inline const Tuple<uint32_t>* hasAnnot(void*) const{return NULL;}
    inline Annotation* canAnnot() const{return NULL;}
};

// stores annotation for structure, which are assumed to reside in a static location in the memory
class Annotation{
public:
    Vector<Dictionary> dicos;
    myHashmap<void*, Tuple<uint32_t> > axes;
    Vector<void*> serialized_buffer;

    template<class C, unsigned int DIMS> Tuple<uint32_t, DIMS> getCoord(const C& target, const Tuple<string, DIMS> &names)const;
    template<class C, unsigned int DIMS> Tuple<string, DIMS> getNames(const C& target, const Tuple<uint32_t, DIMS> &names)const;

    template<class C, unsigned int DIMS> typename C::INNER_TYPE& operator()(C& ,const Tuple<string, DIMS> &names);
    template<class C, unsigned int DIMS> typename C::INNER_TYPE operator()(const C& ,const Tuple<string, DIMS> &names);

    ERRCODE read_Rtable(const char* path, TMatrix<double, 0u,0u>& target);
    ERRCODE write_Rtable(const char* path, const TMatrix<double, 0u,0u>& target) const;
    ERRCODE read_Rtable(const char* path, DataGrid<double, 2u>& target);
    ERRCODE write_Rtable(const char* path, const DataGrid<double, 2u>& target) const;

    bool hasDico(uint32_t dicoID) const{return((dicoID < dicos.getSize())&&(dicos[dicoID].entries.getSize() !=0 ));}
    uint32_t getDicoSize(uint32_t dicoID) const{return (dicoID < dicos.getSize()) ? dicos[dicoID].entries.getSize() : 0;}
    void genDicos(Tuple<uint32_t> &fout);
    void allocDicos(uint32_t nb){dicos.upSize(nb);}

    void useAnnotation()const; // makes annotation available (uses static scope)

    #ifdef Rcpp_hpp
    template<class T> void getAxeNamesAl(Rcpp::CharacterVector &object, T alias) const;
    template<class T> void setAxeNamesAl(const Rcpp::CharacterVector &object, T alias);
    template<class T> void rdAxesNamesAl(const Rcpp::NumericMatrix &object, T rowalias,T colalias);
    template<class T> void wrAxesNamesAl(Rcpp::NumericMatrix &object, T rowalias, T colalias) const;
    template<class T> void rdAxesNamesAl(const Rcpp::S4 &object, T rowalias, T  colalias);
    template<class T> void wrAxesNamesAl(Rcpp::S4 &object, T rowalias, T colalias) const;

    template<class T> void rdAxesNames(const Rcpp::NumericMatrix &object, const T &target);
    template<class T> void wrAxesNames(Rcpp::NumericMatrix &object, const T &target) const;
    template<class T> void rdAxesNames(const Rcpp::S4 &object, const T &target);
    template<class T> void wrAxesNames(Rcpp::S4 &object, const T &target) const;
    #endif

    template<class C> bool hasAnnotation(const C& target)const;
    template<class C> void cloneAnnotation(const C& target, const C& source);

    template<class C> unsigned int findDico(const C& target, uint32_t offset)const;
    template<class C> unsigned int findDico(const C* target, uint32_t offset)const;

    ERRCODE save(FILE*f) const;
    ERRCODE load(FILE*f, unsigned int s =0u);
};

class Metable{
public:
    SparseMatrix<uint32_t> data;
    Vector<uint32_t> annotAlias;
    uint32_t rowAlias;
    uint32_t colAlias;

    #ifdef Rcpp_hpp
    Rcpp::IntegerVector genFactor(uint32_t col) const;
    Rcpp::IntegerVector genFactor(string col) const;
    #endif
};


template<class NODE, int nbrel>
class Forest : public Vector< pair<Tuple<unsigned int, nbrel>, NODE > >{
public:
	unsigned int nbroots; // number of roots
	void cluster(Vector<NODE> &data, double (*metric)(const NODE&, const NODE&), NODE (*merge)(const NODE&, const NODE&) = NULL);
	void cluster_singlelink(Vector<NODE> &data, double (*metric)(const NODE&, const NODE&));
	void cluster_completelink(Vector<NODE> &data, double (*metric)(const NODE&, const NODE&));

    //class TaskClustering;


    template<class DATA, unsigned int TSIZE > void cluster_bhattacharryya( const Vector< Tuple<DATA, TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
    template<class DATA, unsigned int TSIZE > void cluster_likelihood_ratio( const Vector< Tuple< WeightElem<DATA, 2> , TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));

    template<class DATA, unsigned int TSIZE > void cluster_likelihood_ratio( const Vector< GaussElem< Tuple<DATA, TSIZE > > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));

    template<class DATA> void cluster_likelihood_ratio( const Vector< GaussElem< Tuple<DATA, 0u > > > &data);

    template<class DATA, class TB> void cluster_likelihood_ratio( const Vector< GaussElem< Tuple<DATA, 0u > > > &data, TB &threadbase);

	template<class DATA, unsigned int TSIZE > void cluster_bhattacharryya_simple( const Vector< Tuple< WeightElem<DATA, 2> , TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
    template<class DATA, unsigned int TSIZE > void cluster_bhattacharryya_complex( const Vector< GaussElem< Tuple<DATA, TSIZE > > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
	template<class DATA, unsigned int TSIZE > void cluster_likelihood( const Vector< Tuple< WeightElem<DATA, 2> , TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
	class Task_cluster_likelihood;
	ERRCODE cluster_likelihood(const SparseMatrix< double > &data);
	template<class DATA> void cluster(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&), DATA (*merge)(const DATA&, const DATA&) = NULL);
	template<class DATA> void cluster_singlelink(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));
	template<class DATA> void cluster_completelink(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));
template<class DATA> void cluster_majoritylink(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));


	template<class DATA> void partition_squarreroot(const Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));
	template<class DATA, class COMP> void partition_squarreroot_nbcluster(const Vector<DATA> &data, COMP (*metric)(const DATA&, const DATA&), int nbgroups);

	void saveGTRfile(const char* const path, const char* const  outer, double (*metric)(const NODE&, const NODE&)) const;
	void saveGTRfile(const char* const path, const char* const  outer, double (*xform)(const NODE&) = NULL) const; // simple convertion otherwise
    template<class DATA> void saveCDTfile(FILE* f, Vector<char*> &header , Vector<DATA> &data) const;
	template<class DATA, unsigned int SIZE, unsigned int ORDER> void saveCDTfile_W(FILE* f, Vector<char*> &header , Vector< Tuple<WeightElem<DATA, ORDER>, SIZE >  > &data) const;
	template<class DATA, unsigned int SIZE> void saveCDTfile_W(FILE* f, Vector<char*> &header , Vector< GaussElem< Tuple<DATA, SIZE > > > &data, const Tuple<unsigned int,SIZE> * permute = NULL) const;

	void makeLeftof(unsigned int left, unsigned int par);
	void makeRightof(unsigned int left, unsigned int par);
	bool hasParent(unsigned int node);
	bool hasLeft(unsigned int node);
	bool hasRight(unsigned int node);
	unsigned int getParent(unsigned int node)const;
	unsigned int getLeft(unsigned int node)const;
	unsigned int getRight(unsigned int node)const;
	void clear(unsigned int node);
	bool isclear(unsigned int node);

	void moveNode(unsigned int Node, unsigned int Target);

	Vector<uint32_t>& wrAsPermutation(Vector<uint32_t> &perm)const;

	void show(FILE* f =stdout, int level =0) const;
};

class IntegerSampler{
public:
    Tuple<uint32_t> prob;
    uint32_t imag;
    void init(double* prob, uint32_t length);
    uint32_t sample()const;
};

class CMPWindowDistibution{
public:
    Tuple<double> log_start_prior;
    TMatrix<double> transition; // diagonal contain sum of rates
    Tuple<double> condfactor;
    TMatrix<double> cached_coefs;

    Tuple<SparseTuple<double> > mkSamples(int32_t nbsamples) const;
    SparseTuple<double> mkSample() const{return std::move(this->mkSamples(1)[0]);}
//    SparseTuple<double> mkSample() const;
    CMPWindowDistibution& setRatesAndPrior(const Tuple<double> &prior, const TMatrix<double> &rates, uint32_t maxlength);
    // return LL or log probability for tuple

    class LearnScope{
    public:
        const CMPWindowDistibution& target;
        Tuple<double> d_log_start_prior;
        TMatrix<double> d_transition;
        double LL;
        Tuple<uint32_t> lenghtdist;

        LearnScope(const CMPWindowDistibution& _target):target(_target){}
        void init();
        double EMregist(const SparseTuple<double> &); // outputs LL
    };
    LearnScope mkEMscope()const;

    double operator()(const SparseTuple<double> &)const;
    void show(FILE*f=stdout, int level=0)const;
};

// simpler class with directed graph, the transition probability on the
// Related to : Discrete phase-type distribution
class MPReachDistribution{
public:
    Tuple<double> log_start_prior;
    Trianglix<double> log_transition; // The probability on the diagonal is probability to the end state (end of path)
    double underbound_for_exit_probability; // domain is constrained to have Expected path length below a threshold.
    // states are onlt allowed to decrease.
    class LearnScope{
    public:
        const MPReachDistribution& target;
        Tuple<uint32_t> start_count;
        Trianglix<uint32_t> transition_count;
        double LL;
        double prior_count;

        LearnScope(const MPReachDistribution& _target);
        void init();
        double EMregist(const SparseTuple<double> &); // outputs LL
    };

    Tuple<uint32_t, 2u> findCandidateStep(SparseTuple<double> &fout, const SparseTuple<double> &start_io, uint32_t& outdirhit, uint32_t nbdim, const double* derivative, double LL_prior_factor)const;
    uint32_t getSize()const{return log_start_prior.getSize();}
    double operator()(const SparseTuple<double> &)const;
    LearnScope mkEMscope()const{ MPReachDistribution::LearnScope fout(*this); return fout;}
    double learn(LearnScope& scp, double prior_count, bool do_show=false); // return exact LL *after* learning
    Tuple<SparseTuple<double> > mkSamples(int32_t nbsamples) const;
    SparseTuple<double> mkSample() const{return std::move(this->mkSamples(1)[0]);}


    template<class C> class Task_learnPermutation;
    template<class C> Tuple<uint32_t> learnPermutation(const SparseMatrix<C>& data, double prior_count, double &LL_to_change); // return permutation of data



    MPReachDistribution& setDefaultRatesAndPrior(uint32_t nbdims, double pathlength);
    MPReachDistribution& setRatesAndPrior(const Tuple<double> &prior, const Trianglix<double> &rates, double pathlength);
    void show(FILE*f=stdout, int level=0)const;
};


/*
template<class NODE, int nbrel>
class Forest<NODE,nbrel>::TaskClustering : public Event<uint32_t>{
    public:
    TaskClustering(uint32_t nbthread);
    void operator()(uint32_t threadID);

    template<class TB> void run(TB& thrbase);
};*/

// GSL functions

#ifdef GNU_SCIENTIFIC_LIBRARY
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
#endif


#include "DataGrid.hpp"
// #include "ExOp.hpp"

//#if __cplusplus >= 201103L
//  #include "require11.h"
//#endif


#include "stats.h"
#include "primitive.hpp"

} // end of namespace LFHPrimitive



#endif


