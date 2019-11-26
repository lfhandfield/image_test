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



#ifndef _defined_PrimitiveDeclare
#define _defined_PrimitiveDeclare

#define GNU_SCIENTIFIC_LIBRARY
#ifndef LFH_HAS_RUNNINGTIME_STATISTICS
#define LFH_HAS_RUNNINGTIME_STATISTICS true
#endif

//#include <armadillo>


#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#define __STDC_FORMAT_MACROS 1
#include <cinttypes>
#include <cstdint>
#include <chrono>
//#include <stdarg.h>
//#include <stddef.h>

#include <string.h>

#include <functional>
#include <thread>
#include <mutex>              // std::mutex, std::unique_lock
#include <condition_variable> // std::condition_variable
#include <future>
#define STDBINDFILTER


#include <time.h>
#include <cmath>
#include <map>

#include <stdint.h>
#include <list>
#include <queue>
#include <stack>
#include <vector>
#include <iostream>
//#include <mutex>
//#include <type_traits>

#include <cfloat>

using namespace std;

//define to include <SDL/SDL.h> instead fo "SDL.h", for example:
//#define LFH_NESTED_LIBRARIES


#define LFH_STR(tok) #tok

#define STRINGIFY(A) #A
#define MAKE_COMMA ,

typedef unsigned char ERRCODE;
/* No such file or directory */
#define ERRCODE_ILLEGAL 1
/* No such file or directory */
#define ERRCODE_NOFILE 2
/* 11: Try again */
#define ERRCODE_BUSY 11
/* 12: Out of memory */
#define ERRCODE_OUTOFMEM 12

// Target Table
//                     One client                 Many clients             Dependent Client (needs simutaneous destory )
//             Immovible        Ammovible  Immovible      Ammovible        Immovible      Ammovible
// Permanent   pointer        double-ptr     ptr          alias            ptr            alias/listener
// Temporary   double-ptr     double-ptr   alias          alias         alias/listener    alias/listener
//
#define LFH_DEFAULT_PORT "27015"

#define LFH_GOLD
#define LFH_FAULTY

#define LFH_OPT_CHECK(CMD) CMD
// #define LFH_OPT_CHECK(CMD)

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#define LFHDebug_RETURN_void(Test) if (!(Test)) return
#define LFHDebug_RETURN(Test, out) if (!(Test)) return(out)

#define LFHDebug(Test, String) if (Test) fprintf(stderr, String)
#define LFHDebug_LOOP(Test) while(true)

#define LFHPRIMITIVE_ANGCONVERT 0.000095873799242852576857380474343247f

//#define myexit(text) exit(MessageBox(NULL,text,"ERROR",MB_OK|MB_ICONSTOP))

#define myalive() fprintf(stdout, "Reached %s[%d] at %d tics\n", __FILE__, __LINE__, (uint32_t)clock());fflush(stdout)
#define myexit(MyExItTeXt) exit(1 | (0 * fprintf(stderr,"Exiting from %s[%d]:%s\n", __FILE__, __LINE__, "MyExItTeXt")))

#ifdef Rcpp_hpp
#define MYIOCHECK(CMD, VALUE) if ((my_tmpvariable_for_error_msgs = CMD) != VALUE) exit(1 | (0 * fprintf(stderr,"Exiting from %s[%d]: IO return value is %i, but %i was expected\n", __FILE__, __LINE__, my_tmpvariable_for_error_msgs, VALUE)))
#else // Rcpp_hpp
#define MYIOCHECK(CMD, VALUE) CMD
#endif // Rcpp_hpp

#define myexit1() exit(1 | fprintf(stderr, "Myexit from file %s line#%d\n", __FILE__, __LINE__))

#define mybreak(text) MessageBox(NULL,text,"ERROR",MB_OK|MB_ICONSTOP)

#ifndef M_PI_2
#define M_PI_2 ((double)1.5707963267948966192313216916398)
#endif
#ifndef M_LN10
#define M_LN10		2.30258509299404568402
#endif
#ifndef M_LN2
# define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif
#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880
#endif
#ifndef M_E
#define M_E		2.7182818284590452354
#endif

#ifndef M_PI
#define M_PI ((double)3.1415926535897932384626433832795)
#endif
#ifndef M_SQRT_1_2
#define M_SQRT_1_2 ((double)0.70710678118654752440084436210485)
#endif
#ifndef M_FLOAT_NAN
#define M_FLOAT_NAN ((float)0.0f/ 0.0f)
#endif
#ifndef M_2PI
#define M_2PI (6.283185307179586476925286766559)
#endif
#define AngConvert 0.000095873799242852576857380474343247

#define LFHCLONE(PoInTeR, ClAsSnAmE) ((PoInTeR == NULL) ? NULL : (ClAsSnAmE)(PoInTeR->clone()) )

#define LFHDECL_DESCTRUCTOR(ClAsSnAmE) virtual ~ClAsSnAmE(); ClAsSnAmE(const ClAsSnAmE & other); ClAsSnAmE& operator=(const ClAsSnAmE & other);
#define LFHDECL_VIRTUALCLONE(ClAsSnAmE) virtual void* clone()=0;
#define LFHDECL_CLONE(ClAsSnAmE) virtual void* clone() {return (void*) new ClAsSnAmE(*this) ;}

#define LFHDECL_COPYCON_BODY(nAmEsPaCe , ClAsSnAmE) nAmEsPaCe::ClAsSnAmE ::ClAsSnAmE(const nAmEsPaCe::ClAsSnAmE & other)
#define LFHDECL_COPYCON_BODY_AUTO(nAmEsPaCe , ClAsSnAmE) nAmEsPaCe::ClAsSnAmE ::ClAsSnAmE(const nAmEsPaCe::ClAsSnAmE & other)
#define LFHDECL_ASSIGN_BODY(nAmEsPaCe , ClAsSnAmE , BoDy) nAmEsPaCe::ClAsSnAmE & nAmEsPaCe::ClAsSnAmE ::operator=(const nAmEsPaCe::ClAsSnAmE & other)
#define LFHDECL_CLONEINIT(MeMbEr , TyPe) MeMbEr = (TyPe)other.MeMbEr->clone()
#define LFHDECL_CLONEASSIGN(MeMbEr , TyPe) delete(MeMbEr); MeMbEr = (TyPe)other.MeMbEr->clone()

#define LFHDECL_OPER1(ClAsS_A) public Oper1<ClAsS_A> { public: void operator()(ClAsS_A &) const;}
#define LFHDECL_OPER2(ClAsS_A, ClAsS_B) public Oper2<ClAsS_A, ClAsS_B> { public: void operator()(ClAsS_A &, ClAsS_B &) const;}
#define LFHDECL_OPER3(ClAsS_A, ClAsS_B, ClAsS_C) public Oper3<ClAsS_A, ClAsS_B, ClAsS_C> { public: void operator()(ClAsS_A &, ClAsS_B &, ClAsS_C &) const;}

// more of a reminder:
#define LFHCONCAT2E(a,b) a##b

#define LFHCONCAT2(a,b) a,b
#define LFHCONCAT3(a,b,c) a,b,c
#define LFHDECL_TRIVIAL_OPERATOR(oPeRaToR, oUtPuT) template<class TrIvIaLaRgUmEnT> oUtPuT& operator oPeRaToR(TrIvIaLaRgUmEnT const &other) {return( ((oUtPuT(*this)) oPeRaToR##= other) );}

#define LFHALIVE1(text) fprintf(stdout, "Debug flag: %s on %s::line#%d  (time=%i)\n", text, __FILE__, __LINE__, clock());fflush(stdout)
#define LFH_STOPPOPUP(a) printf("%c", printf(a) == -1 ? '\t' : '\n')
#define LFH_ENDPOPUP(a) exit(printf("%c", printf(a) == -1 ? '\t' : '\n'))

#define LFH_address unsigned int
#define myoffsetof(A,B,C) (int(&(((A*)0)->B)))

#define LFH_ALIVE printf("%i th line in %s reached within func %s\n", __LINE__ , __FILE__ ,  __func__) | fflush(stdout)
#define LFH_FUNCTION_MISSING printf("function  %s in %s is not done yet!\n", __func__ , __FILE__);fflush(stdout);exit(1)

#define LFH_ASSERT(A,B) if (!(A)) exit( 1 | fprintf(stderr, B ) | fprintf(stderr, "Assert Error at %s in %s\n", __func__ , __FILE__) )
#define LFH_FUNCTION_MISSING printf("function  %s in %s is not done yet!\n", __func__ , __FILE__);fflush(stdout);exit(1)


#define LFH_VERB(V,T) if (LFHPrimitive::verbose >= V) {printf("%s.%s[Line%i]: %s\n", __FILE__ ,  __func__ , __LINE__, T );fflush(stdout);}

#define ExCoMeMdEcLaRe(ClAsS_A) ERRCODE save(FILE* f) const; ERRCODE load(FILE* f, unsigned int arr_size = 0u); void show(FILE* f = stdout, int level = 0) const; string type_tostring()const; void zero();  void one(); void random();

#define LFH_NICE_ALLOCERROR(PRT,...)
//#define LFH_NICE_ALLOCERROR(PRT,...) if (PRT == nullptr) {thrbase.print("Memory allocation error in function %s (%s[Line:%d])\n", __func__, __FILE__, __LINE__);  thrbase.print(__VA_ARGS__); exit(12);}

//#define LFH_NICE_ALLOCERROR(PRT,CMD)

#define LFH_VALGRIND_MUTE(CMD)   CMD
//#define LFH_VALGRIND_MUTE(CMD)


// atomic operation
inline int fetch_and_add( int * variable, int value ) {
     asm volatile("lock; xaddl %%eax, %2;"
                  :"=a" (value)                  //Output
                  :"a" (value), "m" (*variable)  //Input
                  :"memory");
     return value;
 }

#ifdef GNU_SCIENTIFIC_LIBRARY
struct gsl_sf_result_struct {
	double val;
	double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
	double val;
	double err;
	int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;
#endif


namespace LFHPrimitive{



class ScopeVerbose{
    const char* text;
    public:
    ScopeVerbose(const char* what): text(what) {printf("%s <-\n", text);fflush(stdout);}
    ~ScopeVerbose(){printf("%s ->\n", text);fflush(stdout);}
};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Meta-programmation section
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	template <bool A> class YESNO{	public:	enum {ans = A};	};

	enum FLAGFORTYPE{
		FLAGFORTYPE_NULL =0,
		FLAGFORTYPE_CONST =1,
		FLAGFORTYPE_REF =2,
		FLAGFORTYPE_CONSTREF =3,
		};

	template< class A, FLAGFORTYPE isconst>
	class TYPEMANIP{public:typedef A TYPE;};


    template<int I> class classInt{public: const static int val = I;};
    template<unsigned int I> class classUInt{public: const static unsigned int val = I;};


	template<class A> class TYPEMANIP<A, FLAGFORTYPE_CONST>{public: typedef const A TYPE;};
	template<class A> class TYPEMANIP<A, FLAGFORTYPE_REF>{public:typedef A& TYPE;};
	template<class A> class TYPEMANIP<A, FLAGFORTYPE_CONSTREF>{	public:	typedef const A& TYPE;	};

    template<class A> class TYPEMANIP<const A, FLAGFORTYPE_NULL>{public: typedef A TYPE;};
	template<class A> class TYPEMANIP<const A, FLAGFORTYPE_REF>{public:typedef A& TYPE;};
	template<class A> class TYPEMANIP<const A, FLAGFORTYPE_CONSTREF>{	public:	typedef const A& TYPE;	};

    template<class A> class TYPEMANIP<A&, FLAGFORTYPE_NULL>{public: typedef A TYPE;};
	template<class A> class TYPEMANIP<A&, FLAGFORTYPE_CONST>{public: typedef const A TYPE;};
	template<class A> class TYPEMANIP<A&, FLAGFORTYPE_CONSTREF>{	public:	typedef const A& TYPE;	};

    template<class A> class TYPEMANIP<const A&, FLAGFORTYPE_NULL>{public: typedef A TYPE;};
	template<class A> class TYPEMANIP<const A&, FLAGFORTYPE_CONST>{public: typedef const A TYPE;};
	template<class A> class TYPEMANIP<const A&, FLAGFORTYPE_REF>{public:typedef A& TYPE;};

	extern unsigned int def_font[];
    enum metaop_type{
        METAOP_ZERO,
        METAOP_AND,
        METAOP_OR,
        METAOP_XOR,
        METAOP_PLUS,
        METAOP_MINU,
        METAOP_PROD,
        METAOP_DIVI,
        METAOP_MOD,
    };

enum ORTHOREF3_enum{
	ORTHOREF3_NULL=0,
	ORTHOREF3_FLIP_X=1,
	ORTHOREF3_FLIP_Y=2,
	ORTHOREF3_ROT_Z=3,
	ORTHOREF3_FLIP_Z=4,
	ORTHOREF3_ROT_Y=5,
	ORTHOREF3_ROT_X=6,
	ORTHOREF3_FLIP=7,
	ORTHOREF3_ROT_XP=10,
	ORTHOREF3_ROT_YP=44,
	ORTHOREF3_ROT_ZP=25,
	ORTHOREF3_ROT_XN=12,
	ORTHOREF3_ROT_YN=41,
	ORTHOREF3_ROT_ZN=26
};

    template<metaop_type METAOP, class A, class B, class TYPE = int>
    class MetaOperation{
    public:
        static const TYPE val=0;
    };

    template<class A, class B, class TYPE> class MetaOperation<METAOP_AND,A,B,TYPE>{ public: static const TYPE val=(A::val) & (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_OR,A,B,TYPE>{ public: static const TYPE val=(A::val) | (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_XOR,A,B,TYPE>{ public: static const TYPE val=(A::val) ^ (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_PLUS,A,B,TYPE>{ public: static const TYPE val=(A::val) + (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_MINU,A,B,TYPE>{ public: static const TYPE val=(A::val) - (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_PROD,A,B,TYPE>{ public: static const TYPE val=(A::val) * (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_DIVI,A,B,TYPE>{ public: static const TYPE val=(A::val) / (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_MOD,A,B,TYPE>{ public: static const TYPE val=(A::val) % (B::val); };

	// meant for small fatorial only! 12! at most
	template <int base, int nbfact = (base< 13) ? base : 0> class TEMPLATE_FACTORIAL{
public:
	enum {ans = base * TEMPLATE_FACTORIAL<base-1, nbfact-1>::ans};
};

template<unsigned int which> class TEMPLATE_PRIME{
public:
    static const unsigned int val = 0x7FFFFFFF;
};
template< >  class TEMPLATE_PRIME<1u>{public:	static const unsigned int val = 2;};
template< >  class TEMPLATE_PRIME<2u>{public:	static const unsigned int val = 3;};
template< >  class TEMPLATE_PRIME<3u>{public:	static const unsigned int val = 5;};
template< >  class TEMPLATE_PRIME<4u>{public:	static const unsigned int val = 7;};
template< >  class TEMPLATE_PRIME<5u>{public:	static const unsigned int val = 11;};
template< >  class TEMPLATE_PRIME<6u>{public:	static const unsigned int val = 13;};
template< >  class TEMPLATE_PRIME<7u>{public:	static const unsigned int val = 17;};
template< >  class TEMPLATE_PRIME<8u>{public:	static const unsigned int val = 19;};
template< >  class TEMPLATE_PRIME<9u>{public:	static const unsigned int val = 23;};
template< >  class TEMPLATE_PRIME<10u>{public:	static const unsigned int val = 29;};
template< >  class TEMPLATE_PRIME<11u>{public:	static const unsigned int val = 31;};
template< >  class TEMPLATE_PRIME<12u>{public:	static const unsigned int val = 37;};

template <int base> class TEMPLATE_FACTORIAL<base,0>{
public:
	enum {ans = 1};
};

template <int value> class TEMPLATE_CONSTANTS{
public:
    static constexpr double ONE_OVER_X = 1.0 / value;
};



template <class C, int S> class TEMPLATE_SIZECHECK{public:	static const bool isLT = sizeof(C) < S; static const bool isLE = sizeof(C) <= S; static const bool isEQ = sizeof(C) == S; static const bool isGT = sizeof(C) > S; static const bool isGE = sizeof(C) >= S; static const bool isNQ = sizeof(C) != S;};

template<unsigned int which>
class StConstList{};
template< > class StConstList<0u>{public: static const unsigned int prime = 2; static const unsigned int twopower = 1;};
template< > class StConstList<1u>{public: static const unsigned int prime = 3; static const unsigned int twopower = 2;};
template< > class StConstList<2u>{public: static const unsigned int prime = 5; static const unsigned int twopower = 4;};
template< > class StConstList<3u>{public: static const unsigned int prime = 7; static const unsigned int twopower = 8;};
template< > class StConstList<4u>{public: static const unsigned int prime = 11; static const unsigned int twopower = 16;};
template< > class StConstList<5u>{public: static const unsigned int prime = 13; static const unsigned int twopower = 32;};
template< > class StConstList<6u>{public: static const unsigned int prime = 17; static const unsigned int twopower = 64;};
template< > class StConstList<7u>{public: static const unsigned int prime = 19; static const unsigned int twopower = 128;};
template< > class StConstList<8u>{public: static const unsigned int prime = 23;	static const unsigned int twopower = 256; typedef char SINT;	typedef unsigned char UINT;};
template< > class StConstList<9u>{public: static const unsigned int prime = 29; static const unsigned int twopower = 512;};
template< > class StConstList<10u>{public: static const unsigned int prime = 31; static const unsigned int twopower = 1024;};
template< > class StConstList<11u>{public: static const unsigned int prime = 37; static const unsigned int twopower = 2048;};
template< > class StConstList<12u>{public: static const unsigned int prime = 43; static const unsigned int twopower = 4096;};
template< > class StConstList<16u>{public: static const unsigned int prime = 2;	typedef short SINT;	typedef unsigned short UINT;};
template< > class StConstList<32u>{public: static const unsigned int prime = 2;	typedef int SINT;	typedef unsigned int UINT;};
template< > class StConstList<64>{public: static const unsigned int prime = 2;typedef long int SINT; typedef unsigned long int UINT;};

// devides number by it prime decomposition, IF where each division by a prime power is comstrained to be smaller or equal to "base"
template <unsigned int val, unsigned int base, unsigned int ite = 1u, unsigned int div = 1u> class TEMPLATE_DIVIDE_BY_MAX_DIVISOR{
public:
static const unsigned int ans = ( TEMPLATE_PRIME<ite>::val > base ? val :
			 ((((val % TEMPLATE_PRIME<ite>::val) == 0u)&&(TEMPLATE_PRIME<ite>::val <= base / div)) ?
			  TEMPLATE_DIVIDE_BY_MAX_DIVISOR<val / TEMPLATE_PRIME<ite>::val, base, ite, div * TEMPLATE_PRIME<ite>::val >::ans   :
			  TEMPLATE_DIVIDE_BY_MAX_DIVISOR<val, base, ite+1, 1u>::ans)  );
};

template <unsigned int val, unsigned int base, unsigned int div> class TEMPLATE_DIVIDE_BY_MAX_DIVISOR<val, base,11u,div>{
public:	static const unsigned int ans = val;};

template <unsigned int base, unsigned int ite, unsigned int div> class TEMPLATE_DIVIDE_BY_MAX_DIVISOR<0u, base,ite,div>{
public:	static const unsigned int ans = 0;};

// (n+k)! / (n! * k!) safe to overflow if either n or k <= 6
template <unsigned int lenght,unsigned int nbdims> class TEMPLATE_TRIANGLE_NUMBER{ // not safe to overflows...
public:	enum {ans = (nbdims > lenght) ? TEMPLATE_TRIANGLE_NUMBER<nbdims, lenght>::ans : ((lenght + nbdims -1) * TEMPLATE_TRIANGLE_NUMBER<lenght, nbdims- 1 >::ans) / nbdims};};

// template <unsigned int nbdims> class TEMPLATE_TRIANGLE_NUMBER<0u, nbdims>{	public:	static const unsigned int ans = 0u;	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 0u>{	public:	static const unsigned int ans = 1u;	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 1u>{	public:	static const unsigned int ans = lenght;	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 2u>{	public: static const unsigned int ans = (TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,2u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,2u>::ans);};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 3u>{	public: static const unsigned int ans = (((lenght % 2) == 0 ? 2 : 1) * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,3u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,3u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+2,3u>::ans);	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 4u>{	public: static const unsigned int ans = (((lenght % 3) == 0 ? 3 : 1) * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,4u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,4u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+2,4u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+3,4u>::ans);	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 5u>{	public: static const unsigned int ans = ((lenght % 2) == 0 ? ((lenght % 4) == 0 ? 4 : 2) : 1) * ((lenght % 3) != 1 ? 3 : 1) * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,5u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,5u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+2,5u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+3,5u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+4,5u>::ans;	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 6u>{	public: static const unsigned int ans = (((lenght+1) & 2) ? 1 : 2) * ((lenght % 5) == 0 ? 5 : 1) * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,6u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,6u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+2,6u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+3,6u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+4,6u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+5,6u>::ans;	};


template <int base, int exp> class TEMPLATE_INT_POWER{
public: enum {ans = base * TEMPLATE_INT_POWER<base, exp - 1>::ans};};

template <int base> class TEMPLATE_INT_POWER<base,0>{
public:	enum {ans = 1};};

template <int val, class C, class D> class MT_IFTYPE{public:typedef C TYPE;};
template <class C, class D> class MT_IFTYPE<0,C,D>{public: typedef D TYPE;};

template <class C> class TEMP_IS_CONST{ public: static const bool ans = false;  };
template <class C> class TEMP_IS_CONST<const C>{ public: static const bool ans = true;  };

template <class C> class TEMP_IS_REF{ public: static const bool ans = false; };
template <class C> class TEMP_IS_REF<C&>{ public: static const bool ans = true;  };

class ExOp;

// macro to constrain functor to have member function, need to be last in the list of argument (sets default value), differs for declaration and definition
template<typename T> struct argument_type;
template<typename T, typename U> struct argument_type<T(U)> { typedef U type; };
template<typename T> struct argument_type<T()> { typedef uint32_t type; };


#define FUNCTOREQUIRES_DECL(F , O , A ) argument_type<void(O)>::type (argument_type<void(F)>::type::*_##F##var)(argument_type<void(A)>::type) const = argument_type<void(F)>::type::operator()
#define FUNCTOREQUIRES_DEF(F , O , A ) argument_type<void(O)>::type (argument_type<void(F)>::type::*_##F##var)(argument_type<void(A)>::type) const

#define FUNCTOREQUIRES_DECL1(F , O , A ) argument_type<void(O)>::type (argument_type<void(F)>::type::*_##F##var1)(A) const  = argument_type<void(F)>::type::operator()
#define FUNCTOREQUIRES_DEF1(F , O , A ) argument_type<void(O)>::type (argument_type<void(F)>::type::*_##F##var1)(A) const

#define FUNCTOREQUIRES_DECL2(F , O , A, B ) argument_type<void(O)>::type (argument_type<void(F)>::type::*_##F##var2)(A,B) const = argument_type<void(F)>::type::operator()
#define FUNCTOREQUIRES_DEF2(F , O , A, B ) argument_type<void(O)>::type (argument_type<void(F)>::type::*_##F##var2)(A,B) const

#define ACCEPTOR_DECL(F , A) void (argument_type<void(F)>::type::*_##F##acceptvar)(A) = argument_type<void(F)>::type::operator()
#define ACCEPTOR_DEF(F , A) void (argument_type<void(F)>::type::*_##F##acceptvar)(A)

#define ITERABLE_DECL(Z) typename argument_type<void(Z)>::type::ConstIterator (argument_type<void(Z)>::type::*_Fiterable_var)() const = &argument_type<void(Z)>::type::getIterator
#define ITERABLE_DEF(Z) typename argument_type<void(Z)>::type::ConstIterator (argument_type<void(Z)>::type::*_Fiterable_var)() const



// convention
template<char C> class LTRType{typedef void TYPE;};
template< > class LTRType<'c'>{typedef char TYPE;};
template< > class LTRType<'C'>{typedef unsigned char TYPE;};
template< > class LTRType<'s'>{typedef StConstList<16>::SINT TYPE;};
template< > class LTRType<'S'>{typedef StConstList<16>::UINT TYPE;};
template< > class LTRType<'l'>{typedef StConstList<64>::SINT TYPE;};
template< > class LTRType<'L'>{typedef StConstList<64>::UINT TYPE;};
template< > class LTRType<'i'>{typedef StConstList<32>::SINT TYPE;};
template< > class LTRType<'I'>{typedef StConstList<32>::UINT TYPE;};
template< > class LTRType<'f'>{typedef float TYPE;};
template< > class LTRType<'d'>{typedef double TYPE;};

template<class A> class IsLFHPrimitive{public: enum{ans = 0};};

template <class A, class B> class isTypeEquivalent {enum {ans = false };};

template<class lhs, class rhs> struct isTypeEqual {enum {ans = false }; };
template<class lhs> struct isTypeEqual<lhs, lhs> {enum { ans = true };  };

template<class A, class B> class STDRETTYPE2;

	template<class A, class B>
	class STDRETTYPE2{
		public:
		typedef A DEF_TYPE;
		typedef A PLUS_TYPE; // operator+
		typedef A MINU_TYPE; // operator-
		typedef A PROD_TYPE; // operator*
		typedef A DIVI_TYPE; // operator/
		typedef A OP2_TYPE; // operator()

        typedef A RMUL_TYPE; // mult_right  a *bH
		typedef A LMUL_TYPE; // mult_left   bH * a
        typedef A RDIV_TYPE; // mult_right  a / bH
		typedef A LDIV_TYPE; // mult_left   bH \ a
	};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// End of Meta-programattion
// -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
// Class Forward declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	class Function;
	template<class ARG = void> class Event;
	template<class KEY, class DATA, bool isMax = true> class ExtremumScope;
    class UnknownScope;
	template<typename A> class EnumBox;
	class WarH; // Handdles warnings!
	template<bool> class LinkAssert; // if false, gets and error
	class Weight; // simple double, may have special meaning
	template<class C, unsigned int nbstates> class Classifier;
	template<int buffersize> class SuperString;
	class Bluestein;
	template<class C> class TiffFileTypeInterpretation;
	//class PriorityQueue;
	// basic arrays

	#define TUPLE_FLAG_CONSTRUCT_MASK(bAse) (Tuple_flag)(0 & bAse)
	enum Tuple_flag{
		TUPLE_FLAG_NULL =0, // default 2^32 capacity
        TUPLE_FLAG_REMOTE_MEMORY=1024,
        TUPLE_FLAG_SYMBOLIC_HOSTED=1,
	};

	// defined in bastructs.h
	template<class C> class Complex;
	template<class C> class Quaternion;
	template<class C, unsigned int size=0u, Tuple_flag flag= TUPLE_FLAG_NULL> class Tuple;
    template<class C, unsigned int sizex=0u, unsigned int sizey=sizex, Tuple_flag TF = TUPLE_FLAG_NULL> class TMatrix;
	template<class C> class Matrix;
	template<class C> class SparseTuple;
	template<class C> class SparseMatrix;
	// defined in bastructs.h
    // defined in primitive.h
	class TiffFile;
    // defined in primitive.h


	template<class A> class Oper1;
	template<class A, class B> class Oper2;
	template<class A, class B, FLAGFORTYPE AT, FLAGFORTYPE BT> class OperMK2;
	template<class A, class B, class C> class Oper3;

	template<class C, bool hasone> class mantissa;

	//template<class C> class angle;

	template<class I, int ISIZE, class O, int OSIZE> class Continuous;


	template<class C, unsigned int order> class WeightElem_baseclass;
	template<class C, unsigned int order> class WeightElem;
	template<class C, unsigned int flag = 0u> class GaussElem; // order 2 with "covariances"
    template<class C> class PartialGaussElem;

	template<class C, class V> class GaussScope;
    template<class C, unsigned int SIZE = 0u>  class Trianglix;

    class StaticAttribs;
    class Annotation;

	//template<class C, int size> class Matrianglix; // squarre Matrix stored in triangular format!

enum LFHVectorflag{
	LFHVECTOR_NORMAL=0,
	LFHVECTOR_REMOTE=1, // data is not owned! cant change size or deallocate!
	LFHVECTOR_OWNED_LINKMEMORY=2, // data is owned and the Vector must maintain a scope access for its elements, scope is owned so it must be deleted uppon destruction
	LFHVECTOR_LINKMEMORY=3 // data is owned and the Vector must maintain a scope access for its elements
};
    template<class A, int flag = std::is_enum<A>::value > class ExCo; // (std::is_enum<A>::value) ? 1 : 0
	template<class C> class Vector;
	//template<class C, int size> class VectorArray;

    template <class C, int ASYNC_MAG> class AsyncBuffer;
	template <class C, unsigned int ASYNCBUFFER_MAG =0> class HeapTree; // ASYNCBUFFER_MAG =0 for no async buffer

	// abstract class, read_only
	template<class C, unsigned int nbdim> class ConstGrid;
	template<class C, unsigned int nbdim = 2u> class DataGrid;
	template<class C, class D, class FUNC, unsigned int nbdim> class FuncGrid;

	class SetComparison;
	template<class C> class IntervalSet;
	template<class C> class Interval;

	template<unsigned int nbdim> class GaussianDistribution;
	template<class C, int nbrel> class Forest; // 1 = parent, 2 = LR, 3= LRP, 4 = LRBA, 5 = LRBAP
    template<class C> class AppendixPtr;
    template<class C> class RessourcePtr;
	template <class K, class D> class KeyElem;

    template <class Key> class defaultHashFnc;
//    template <class Key, class Data, class HashFnc = defaultHashFnc< ExCo<Key>::INDEX_TYPE > > class myHashmap;
    template <class Key, class Data = void, class HashFnc = defaultHashFnc< typename MT_IFTYPE<isTypeEqual<Data,void>::ans, typename ExCo<Key>::INDEX_TYPE , Key >::TYPE   > > class myHashmap;
    template <class Key, class BackKey, class Data = void, class HashFnc = defaultHashFnc<Key> , class BackHashFnc = defaultHashFnc<BackKey> > class dualHashmap;
    template <class Key, class Data, class Category, class HashFnc = defaultHashFnc<Key> , class Category_HashFnc = defaultHashFnc<Category>  > class CategoryHashmap;
    template <class Node, class Partition> class SetPartition;
    template<class C, class D = void, class HF = defaultHashFnc< typename MT_IFTYPE<isTypeEqual<D,void>::ans, typename ExCo<C>::INDEX_TYPE , C >::TYPE  > > class AsyncHashmap;


    // end of list of classes that are definitly not enums:

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// End of Class Forward declarations
// -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
// Static Classes declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

enum LFHSTATIC_RESSOURCE_ENUM{
    LFHSTATIC_RESSOURCE_ANNOTATIONS
};

extern myHashmap<uint32_t, void*, defaultHashFnc<uint32_t> > lfhstatic_ressources;
extern int verbose;
extern int my_tmpvariable_for_error_msgs;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// End Static Classes declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

enum SETCMP_enum{
	SETCMP_EQUAL=0,
	SETCMP_MASKED_RANGE_ORDERED = 1, // true if (min(X) >= max(Y)) or (max(X) <= min(Y)) (which depends on wether 16 | 64 or 32 | 128 is set)
	// SETCMP_MASKED_RANGE_ORDERED = 256, // true if (max(X) >= min(Y)) or (min(X) <= max(Y)) (which depends on wether 16 | 64 or 32 | 128 is set)
	SETCMP_DISJOINT_BIT=2, // true if A x A y, x != y
	SETCMP_IS_NOT_SUBSET=4, // true if E x st A y, x != y
	SETCMP_IS_NOT_SUPERSET=8, // true if E y st A x, x != y
	SETCMP_GE= 13 |16 |64, // is there is 1 elem,
	SETCMP_LE= 13 |32 |128,
	SETCMP_GT= 15 |16 |64 , // true if A x A y, x != y
	SETCMP_LT= 15 |32 |128, // true if A x A y, x != y
	SETCMP_CMP_E_MASK= 13 |16 |32 |64  |128, // for >= <= operators  X & MASK == _GE or _LE
	SETCMP_CMP_T_MASK= 15 |16 |32 |64  |128, // for > < operators  X & MASK == _GT or _LT
    SETCMP_NOT_DISJOINT_OR_BOUNDED=12, //
    SETCMP_DISJOINT=14,
	SETCMP_MAX_GT_BIT=16, // true if max(X) > max(Y)
	SETCMP_MAX_LT_BIT=32, // true if max(X) < max(Y)
	SETCMP_MIN_GT_BIT=64, // true if min(X) > min(Y)
	SETCMP_MIN_LT_BIT=128, // true if min(X) < min(Y)
	SETCMP_MAX_BITS=16 | 32, // true if min(X) > min(Y)
	SETCMP_MIN_BITS=64 | 128, // true if min(X) < min(Y)
	SETCMP_MAX_GT=16 | 4,
	SETCMP_MAX_LT=32 | 8,
	SETCMP_MIN_GT=64 | 8,
	SETCMP_MIN_LT=128 | 4,
	SETCMP_PARTIAL_GREATER= 16 | 64,
	SETCMP_PARTIAL_SMALLER= 32 | 128,
	SETCMP_SUP_BOUNDED= 32 | 64 | 8,
	SETCMP_SUB_BOUNDED= 16 | 128 | 4,
	SETCMP_EXT_BITS = 16 | 32 | 64 | 128
};
class SetComparison{
	public:
		int comp;
		SetComparison();
		SetComparison(int v);

		template<class O> SetComparison(const O& , const O& , SETCMP_enum mask = SETCMP_CMP_E_MASK); // MASK FOR DESIRED INFO

		bool areMonotonicDisjoint() {return((comp & 3) != 0);}
		bool areDisjoint() {return((comp & SETCMP_DISJOINT_BIT) == SETCMP_DISJOINT_BIT);}
		bool doesIntersect() {return((comp & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT);}

        bool orderEqual() {return((comp & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT);}
        bool orderLE() {return((comp & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT);}
        bool orderEQ() {return((comp & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT);}


        bool partialGT()const{return ((comp & SETCMP_PARTIAL_GREATER) != 0) && ((comp & SETCMP_PARTIAL_SMALLER) == 0);}
        bool partialLT()const{return ((comp & SETCMP_PARTIAL_SMALLER) != 0) && ((comp & SETCMP_PARTIAL_GREATER) == 0);}
        bool partialEQ()const{return (comp & SETCMP_EXT_BITS) == 0;}

        bool minLT()const{return (comp & SETCMP_MIN_LT_BIT) == SETCMP_MIN_LT_BIT;}
        bool minLE()const{return (comp & SETCMP_MIN_GT_BIT) == 0;}
        bool minGT()const{return (comp & SETCMP_MIN_GT_BIT) == SETCMP_MIN_GT_BIT;}
        bool minGE()const{return (comp & SETCMP_MIN_LT_BIT) == 0;}
        bool minNQ()const{return (comp & SETCMP_MIN_BITS) != 0;}
        bool minEQ()const{return (comp & SETCMP_MIN_BITS) == 0;}

        bool maxLT()const{return (comp & SETCMP_MAX_LT_BIT) == SETCMP_MAX_LT_BIT;}
        bool maxLE()const{return (comp & SETCMP_MAX_GT_BIT) == 0;}
        bool maxGT()const{return (comp & SETCMP_MAX_GT_BIT) == SETCMP_MAX_GT_BIT;}
        bool maxGE()const{return (comp & SETCMP_MAX_LT_BIT) == 0;}
        bool maxNQ()const{return (comp & SETCMP_MAX_BITS) != 0;}
        bool maxEQ()const{return (comp & SETCMP_MAX_BITS) == 0;}


		bool doesContain() {return((comp & SETCMP_IS_NOT_SUBSET) == 0);}
		bool isContained() {return((comp & SETCMP_IS_NOT_SUPERSET) == 0);}



		int operator&(const int& mask){ comp &= mask; return(comp);}
		int operator&(const SETCMP_enum& mask){ comp &= (int)mask;return(comp);}

		void show(FILE* f= stdout) const{
			switch(comp){
				case 0: fprintf(f,"Set A = B\n"); break;
				case 1: fprintf(f,"Set A contains B\n"); break;
				case 2: fprintf(f,"Set B contains A\n"); break;
				case 3: fprintf(f,"Set are disjoint\n"); break;
				case 7: fprintf(f,"Set are disjoint and A>=B\n"); break;
				case 11: fprintf(f,"Set are disjoint and B>=A\n"); break;
			}
		};
	};
enum exo{
    EXOP_zero,
};





#define HAS_MEMBER_FUNC(func, name)                                        \
template<typename T,typename Sign>                                 \
struct name {                                                       \
typedef char yes[1];                                            \
typedef char no [2];                                            \
template <typename U, U> struct type_check;                     \
template <typename _1> static yes &chk(type_check<Sign, &_1::func> *); \
template <typename   > static no  &chk(...);                    \
static bool const ans = sizeof(chk<T>(0)) == sizeof(yes);     \
}

/* Makes a iterator structure with the following usage:
 * if (auto ite = ExOp::getIterator(object)) do{
 *   Key key      = ite();
 *   Data& data   = *ite;
 *   DataMember d = ite->member;
 * } while(ite++);
 */
HAS_MEMBER_FUNC(getIterator, Exlisten_getIterator);

HAS_MEMBER_FUNC(toMemmove, Exlisten_toMemmove);
HAS_MEMBER_FUNC(toMemswap, Exlisten_toMemswap);
HAS_MEMBER_FUNC(toMemfree, Exlisten_toMemfree);
HAS_MEMBER_FUNC(getWeight, Exlisten_getWeight);
HAS_MEMBER_FUNC(show, Exlisten_show);


HAS_MEMBER_FUNC(load, Exlisten_load);
HAS_MEMBER_FUNC(save, Exlisten_save);

HAS_MEMBER_FUNC(mkZero, Exlisten_mkZero);
HAS_MEMBER_FUNC(mkRand, Exlisten_mkRand);
HAS_MEMBER_FUNC(mkOne, Exlisten_mkOne);
HAS_MEMBER_FUNC(toZero, Exlisten_toZero);
HAS_MEMBER_FUNC(toUndefined, Exlisten_toUndefined); // uses NaN if double, or zero weight if weight is there

HAS_MEMBER_FUNC(toRand, Exlisten_toRand);
HAS_MEMBER_FUNC(toOne, Exlisten_toOne);
HAS_MEMBER_FUNC(toMin, Exlisten_toMin);
HAS_MEMBER_FUNC(toMax, Exlisten_toMax);



HAS_MEMBER_FUNC(toRightFlood, Exlisten_toRightFlood);
HAS_MEMBER_FUNC(toLeftFlood, Exlisten_toLeftFlood);
HAS_MEMBER_FUNC(toBitInverse, Exlisten_toBitInverse);

HAS_MEMBER_FUNC(isZero, Exlisten_isZero);
HAS_MEMBER_FUNC(isOne, Exlisten_isOne);

HAS_MEMBER_FUNC(lognorm, Exlisten_lognorm);
HAS_MEMBER_FUNC(pnorm, Exlisten_pnorm);
HAS_MEMBER_FUNC(pdist, Exlisten_pdist);

HAS_MEMBER_FUNC(mkMagnitude, Exlisten_mkMagnitude);
HAS_MEMBER_FUNC(mkOrdering, Exlisten_mkOrdering);
HAS_MEMBER_FUNC(tointpow, Exlisten_tointpow);
HAS_MEMBER_FUNC(toSquare, Exlisten_toSquare);
HAS_MEMBER_FUNC(toSqrt, Exlisten_toSqrt);
HAS_MEMBER_FUNC(toInverse, Exlisten_toInverse);
HAS_MEMBER_FUNC(toNegative, Exlisten_toNegative);
HAS_MEMBER_FUNC(toAbs, Exlisten_toAbs);
HAS_MEMBER_FUNC(mkAbs, Exlisten_mkAbs);
HAS_MEMBER_FUNC(toAbsoft, Exlisten_toAbsoft); // abs(x), but enforce a given a minimum value "K" using F(x,K) = K*(1 + ln(cosh(ln(2)*x/K) / ln(2)) (equivalent to K*(log2(2^(x/K)+2^(-x/K)))), its derivative is tanh(ln(2)*x/k)
HAS_MEMBER_FUNC(mkAbsoft, Exlisten_mkAbsoft); // abs(x), but enforce a given a minimum value "K" using F(x,K) = K*(1 + ln(cosh(ln(2)*x/K) / ln(2)) (equivalent to K*(log2(2^(x/K)+2^(-x/K)))), its derivative is tanh(ln(2)*x/k)
HAS_MEMBER_FUNC(toAbhard, Exlisten_toAbhard); // max of abs(x) and K (hence >= K)
HAS_MEMBER_FUNC(mkAbhard, Exlisten_mkAbhard); // max of abs(x) and K (hence >= K)
HAS_MEMBER_FUNC(mkAbsInverse, Exlisten_mkAbsInverse);
HAS_MEMBER_FUNC(mkAbsoftInverse, Exlisten_mkAbsoftInverse);
HAS_MEMBER_FUNC(mkAbhardInverse, Exlisten_mkAbhardInverse);

HAS_MEMBER_FUNC(operator=, Exlisten_toClone); // changes using implicit data type conversion
HAS_MEMBER_FUNC(toConvert, Exlisten_toConvert); // changes using explicit data type conversion

HAS_MEMBER_FUNC(toAdd, Exlisten_toAdd); // *this =  A + B
HAS_MEMBER_FUNC(operator+, Exlisten_mkadd);
HAS_MEMBER_FUNC(operator+=, Exlisten_toadd);

HAS_MEMBER_FUNC(toSubt, Exlisten_toSubt); // *this  =  A - B
HAS_MEMBER_FUNC(operator-, Exlisten_mksubt);
HAS_MEMBER_FUNC(operator-=, Exlisten_tosubt);

HAS_MEMBER_FUNC(toMult, Exlisten_toMult); // =  A * B
HAS_MEMBER_FUNC(operator*, Exlisten_mkmult); // =  A * B
HAS_MEMBER_FUNC(operator*=, Exlisten_tomult); // =  A * B

HAS_MEMBER_FUNC(toDivi, Exlisten_toDivi); // =  A * (1 / B)
HAS_MEMBER_FUNC(operator/=, Exlisten_todivi); // =  A * (1 / B)
HAS_MEMBER_FUNC(operator/, Exlisten_mkdivi); // =  A * (1 / B)

HAS_MEMBER_FUNC(addMult, Exlisten_addMult); // this +=  A * B
HAS_MEMBER_FUNC(subtMult, Exlisten_subtMult); // this -=  A * B




HAS_MEMBER_FUNC(operator<, Exlisten_lt);
HAS_MEMBER_FUNC(operator<=, Exlisten_le);
HAS_MEMBER_FUNC(operator>, Exlisten_gt);
HAS_MEMBER_FUNC(operator>=, Exlisten_ge);
HAS_MEMBER_FUNC(operator==, Exlisten_eq);
HAS_MEMBER_FUNC(operator!=, Exlisten_nq);
HAS_MEMBER_FUNC(operator(), Exlisten_function);

HAS_MEMBER_FUNC(mkrealproj, Exlisten_mkrealproj);
HAS_MEMBER_FUNC(mkimmaproj, Exlisten_mkimmaproj);
HAS_MEMBER_FUNC(mkjmmaproj, Exlisten_mkjmmaproj);
HAS_MEMBER_FUNC(mkkmmaproj, Exlisten_mkkmmaproj);


HAS_MEMBER_FUNC(mkSquare, Exlisten_mkSquare); // is A * A
HAS_MEMBER_FUNC(mkSqrt, Exlisten_mkSqrt); // is A * A
HAS_MEMBER_FUNC(mkTrjuProd, Exlisten_mkTrjuProd); // is A * trju(A)
HAS_MEMBER_FUNC(mkInverse, Exlisten_mkInverse);
HAS_MEMBER_FUNC(mkPosInverse, Exlisten_mkPosInverse); // abs(A^-1)
HAS_MEMBER_FUNC(mkNegative, Exlisten_mkNegative); // -A, aka operator-()
HAS_MEMBER_FUNC(mkNegInverse, Exlisten_mkNegInverse); // -abs(A^-1)
HAS_MEMBER_FUNC(mkPow, Exlisten_mkPow); // x ^ a where a is a double
HAS_MEMBER_FUNC(mkPowInt, Exlisten_mkPowInt); // x ^ a where a is an integer
HAS_MEMBER_FUNC(mkPowInvInt, Exlisten_mkPowInvInt); // x ^ (1/a) where a is an integer (available if mkInverse exists)
HAS_MEMBER_FUNC(toTrju, Exlisten_toTrju); // conjugate transpose (aka transjugate)
HAS_MEMBER_FUNC(mkTrju, Exlisten_mkTrju); // conjugate transpose (aka transjugate)
HAS_MEMBER_FUNC(wrDeterminant, Exlisten_wrDeterminant);
HAS_MEMBER_FUNC(wrTrace, Exlisten_wrTrace);
HAS_MEMBER_FUNC(wrBackMult, Exlisten_wrBackMult); // =  B * A   (or trju(trju(A) * trju(B)), but inner type is I_B*I_A)
HAS_MEMBER_FUNC(toBackMult, Exlisten_toBackMult); // =  B * A
HAS_MEMBER_FUNC(mkBackMult, Exlisten_mkBackMult); // =  B * A
HAS_MEMBER_FUNC(wrBackDivi, Exlisten_wrBackDivi); // =  B * (1/A)
HAS_MEMBER_FUNC(toBackDivi, Exlisten_toBackDivi); // =  B * (1/A)
HAS_MEMBER_FUNC(mkBackDivi, Exlisten_mkBackDivi); // =  B * (1/A)
HAS_MEMBER_FUNC(toInvDivi, Exlisten_toInvDivi); // =  (1/A) * B
HAS_MEMBER_FUNC(mkInvDivi, Exlisten_mkInvDivi); // =  (1/A) * B


HAS_MEMBER_FUNC(toInnerProd, Exlisten_toInnerProd); // =  trju(A) * B
HAS_MEMBER_FUNC(mkInnerProd, Exlisten_mkInnerProd); // =  trju(A) * B
HAS_MEMBER_FUNC(toOuterProd, Exlisten_toOuterProd); // =  A * trju(B)
HAS_MEMBER_FUNC(mkOuterProd, Exlisten_mkOuterProd); // =  A * trju(B)

HAS_MEMBER_FUNC(toBackInnerProd, Exlisten_toBackInnerProd); // =  trju(B) * A   (inner type: A*B)
HAS_MEMBER_FUNC(mkBackInnerProd, Exlisten_mkBackInnerProd); // =  trju(B) * A   (inner type: A*B)
HAS_MEMBER_FUNC(toBackOuterProd, Exlisten_toBackOuterProd); // =  B * trju(A)   (inner type: A*B)
HAS_MEMBER_FUNC(mkBackOuterProd, Exlisten_mkBackOuterProd); // =  B * trju(A)   (inner type: A*B)

HAS_MEMBER_FUNC(setcmp, Exlisten_setcmp);
HAS_MEMBER_FUNC(operator(), Exlisten_mainoper);
HAS_MEMBER_FUNC(isValid, Exlisten_isValid);
HAS_MEMBER_FUNC(toInvalid, Exlisten_toInvalid);

// HAS_MEMBER_FUNC(operator int, Exlisten_int_convert);

HAS_MEMBER_FUNC(getIndex, Exlisten_getIndex);

HAS_MEMBER_FUNC(wrFirstIterator, Exlisten_wrFirstIterator);
HAS_MEMBER_FUNC(wrLastIterator, Exlisten_wrLastIterator);
HAS_MEMBER_FUNC(wrNextIterator, Exlisten_wrNextIterator);
HAS_MEMBER_FUNC(wrPrevIterator, Exlisten_wrPrevIterator);
HAS_MEMBER_FUNC(wrEndIterator, Exlisten_wrEndIterator);

HAS_MEMBER_FUNC(mkouterprod, Exlisten_mkouterprod); // conjugate transpose
HAS_MEMBER_FUNC(mkvectorization, Exlisten_mkvectorization); // conjugate transpose

HAS_MEMBER_FUNC(type_dimentions, Exlisten_type_dimentions);

HAS_MEMBER_FUNC(mkgaussstat, Exlisten_mkgaussstat);

HAS_MEMBER_FUNC(type_tostring, Exlisten_type_tostring);


template<bool hasfunc>
class ExFn{// ExFn<true>
	public:
    double dummy();

    template<class A> inline static bool isZero(const A& what) {return what.isZero();}
    template<class A> inline static bool isOne(const A& what) {return what.isOne();}
    template<class A> inline static A& toZero(A& what) {return what.toZero();}
    template<class A> inline static A& toOne(A& what) {return what.toOne();}
    template<class A> inline static A& toRand(A& what) {return what.toRand();}

    template<class A> inline static A mkZero(const A& what) {return what.mkZero();}
    template<class A> inline static A mkOne(const A& what) {return what.mkOne();}
    template<class A> inline static A mkRand(const A& what) {return what.mkRand();}

    template<class A> inline static A& toUndefined(A& what) {return what.toUndefined();}
    template<class A> inline static A& toMin(A& what) {return what.toMin();}
    template<class A> inline static A& toMax(A& what) {return what.toMax();}
    template<class A> inline static A& toMin(A& what, const A& other) {return what.toMin(other);}
    template<class A> inline static A& toMax(A& what, const A& other) {return what.toMax(other);}


    // builders

    template<class A> inline static A mknegative(const A& a) {return -a;}
    template<class A> inline static A mkInverse(const A& a) {return a.mkInverse();}
    template<class A> inline static A mkPosInverse(const A& a) {return a.mkPosInverse();}
    template<class A> inline static A mkAbs(const A& a) {return a.mkAbs();}
    template<class A> inline static A mkAbsoft(const A& a, const double& b) {return a.mkAbsoft(b);}
    template<class A> inline static A mkAbhard(const A& a, const double& b) {return a.mkAbhard(b);}
    template<class A> inline static A mkAbsInverse(const A& a) {return a.mkAbs();}
    template<class A> inline static A mkAbsoftInverse(const A& a, const double& b) {return a.mkAbsoft(b);}
    template<class A> inline static A mkAbhardInverse(const A& a, const double& b) {return a.mkAbhard(b);}


	template<class A> inline static Tuple<typename ExCo<A>::INDEX_TYPE> mkOrdering(const A& a){return a.mkOrdering();}
	template<class A> inline static int32_t mkMagnitude(const A& a){return a.mkMagnitude();}
    template<class A> inline static A mkPow(const A& a, const double pow) {return a.mkPow(pow);}
    template<class A> inline static A mkPowInt(const A& a, const int pow) {return a.mkPowInt(pow);}
    template<class A> inline static A mkPowInvInt(const A& a, const int pow) {return a.mkPowInvInt(pow);}

    template<class A> inline static A mkSquare(const A& a) {return a.mkSquare();}
    template<class A> inline static typename ExCo<A>::TRJU_TYPE mkTrju(const A& a){return a.mkTrju();}
    template<class A> inline static A mkTrjuProd(const A& a){return a.mkTrjuProd();}
    template<class A> inline static A& toTrju(A& a){return a.toTrju();}
    template<class A, class B> static void wrDeterminant(const A& a, B& b){a.wrDeterminant(b);}
    template<class A, class B> static void wrTrace(const A& a, B& b){a.wrTrace(b);}



    // assigns

    template<class A> inline static A& toNegative(A& a) {return a.toNegative();}
    template<class A> inline static A& toAbs(A& a) {return a.toAbs();}
    template<class A> inline static A& toAbsoft(A& a, const double& b) {return a.toAbsoft(b);}
    template<class A> inline static A& toAbhard(A& a, const double& b) {return a.toAbhard(b);}
    template<class A> inline static A& toInverse(A& a) {return a.toInverse();}
    template<class A, class B> inline static A& toConvert(A& a, const B& b){return a.toConvert(b);}


    template<class A> inline static void tointpow(A& what, const int pow) {what.tointpow(pow);}
    template<class A> inline static A& toSquare(A& a) {return a.toSquare();}
    template<class A> inline static A& toSqrt(A& a) {return a.toSqrt();}

    template<class A> inline static	A& toClone(A& a, const A& b) {return (a = b);}
	template<class A> inline static	A& toAdd(A& a, const A& b) {return (a += b);}
	template<class A> inline static A& toSubt(A& a, const A& b) {return (a -= b);}
	template<class A> inline static A& toMult(A& a, const A& b) {return (a *= b);}
	template<class A> inline static A& toDivi(A& a, const A& b) {return (a /= b);}

    template<class A, class B> inline static A& toClone(A& a, const B& b) {return (a = b);}
	template<class A, class B> inline static A& toAdd(A& a, const B& b) {return (a += b);}
	template<class A, class B> inline static A& toSubt(A& a, const B& b) {return (a -= b);}
	template<class A, class B> inline static A& toMult(A& a, const B& b) {return (a *= b);}
	template<class A, class B> inline static A& toDivi(A& a, const B& b) {return (a /= b);}

	template<class A> inline static	A mkAdd(const A& a, const A& b) {return a + b;}
	template<class A> inline static A mkSubt(const A& a, const A& b) {return a - b;}
	template<class A> inline static A mkMult(const A& a, const A& b) {return a * b;}
	template<class A> inline static A mkDivi(const A& a, const A& b) {return a / b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PLUS_TYPE mkAdd(const A& a, const B& b) {return a + b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::MINU_TYPE mkSubt(const A& a, const B& b) {return a - b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkMult(const A& a, const B& b) {return a * b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::DIVI_TYPE mkDivi(const A& a, const B& b) {return a / b;}

    template<class A> inline static	A& toAdd(A& a, const A& b, const A& c) {return a.toAdd(b,c);}
	template<class A> inline static A& toSubt(A& a, const A& b, const A& c) {return a.toSubt(b,c);}
	template<class A> inline static A& toMult(A& a, const A& b, const A& c) {return a.toMult(b,c);}
	template<class A> inline static A& toDivi(A& a, const A& b, const A& c) {return a.toDivi(b,c);}
	template<class A,class B> inline static	A& toAdd(A& a, const B& b, const B& c){return a.toAdd(b,c);}
	template<class A,class B> inline static A& toSubt(A& a, const B& b, const B& c){return a.toSubt(b,c);}
	template<class A,class B> inline static A& toMult(A& a, const B& b, const B& c){return a.toMult(b,c);}
	template<class A,class B> inline static A& toDivi(A& a, const B& b, const B& c){return a.toDivi(b,c);}
    template<class A,class B,class C> inline static	A& toAdd(A& a, const B& b, const C& c){return a.toAdd(b,c);}
	template<class A,class B,class C> inline static A& toSubt(A& a, const B& b, const C& c){return a.toSubt(b,c);}
	template<class A,class B,class C> inline static A& toMult(A& a, const B& b, const C& c){return a.toMult(b,c);}
	template<class A,class B,class C> inline static A& toDivi(A& a, const B& b, const C& c){return a.toDivi(b,c);}

	template<class A> inline static A& addMult(A& a, const A& b, const A& c) {return a.addMult(b,c);}
	template<class A> inline static A& subtMult(A& a, const A& b, const A& c) {return a.subtMult(b,c);}
	template<class A,class B> inline static A& addMult(A& a, const B& b, const B& c){return a.addMult(b,c);}
	template<class A,class B> inline static A& subtMult(A& a, const B& b, const B& c){return a.subtMult(b,c);}
    template<class A,class B,class C> inline static	A& addMult(A& a, const B& b, const C& c){return a.addMult(b,c);}
	template<class A,class B,class C> inline static A& subtMult(A& a, const B& b, const C& c){return a.subtMult(b,c);}

	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkmult_trju(const A& a, const B& b) {return b * a;}

	template<class A> inline static A mkadd_alt(const A& a, const A& b) {A ca=a; return ca += b;}
	template<class A> inline static A mksubt_alt(const A& a, const A& b) {A ca=a; return ca -= b;}
	template<class A> inline static A mkmult_alt(const A& a, const A& b) {A ca=a; return ca *= b;}
	template<class A> inline static A mkdivi_alt(const A& a, const A& b) {A ca=a; return ca /= b;}

	template<class A, class B> inline static A mkadd_alt(const A& a, const B& b) {A ca=a; return ca += b;}
	template<class A, class B> inline static A mksubt_alt(const A& a, const B& b) {A ca=a; return ca -= b;}
	template<class A, class B> inline static A mkmult_alt(const A& a, const B& b) {A ca=a; return ca *= b;}
	template<class A, class B> inline static A mkdivi_alt(const A& a, const B& b) {A ca=a; return ca /= b;}


	template<class A> inline static typename ExCo<A>::REAL_TYPE mkrealproj(const A& a){return a.mkrealproj();}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkimmaproj(const A& a){return a.mkimmaproj();}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkjmmaproj(const A& a){return a.mkjmmaproj();}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkkmmaproj(const A& a){return a.mkkmmaproj();}

    template<class A> inline static typename ExCo<A>::GAUS_TYPE mkgaussstat(const A& a,double & weight){return a.mkgaussstat( weight);}

    template<class A> inline static typename ExCo<A>::INDEX_TYPE& getIndex(A& a){return a.getIndex();}
    template<class A> inline static typename ExCo<A>::INDEX_TYPE getIndex(const A& a){return a.getIndex();}
    template<class A> inline static typename ExCo<A>::INDEX_TYPE& getIndex(A*& a){return a->getIndex();}
    template<class A> inline static typename ExCo<A>::INDEX_TYPE getIndex(const A*& a){return a->getIndex();}

	template<class A, class B> inline static A& toBackMult(A& a, const B& b) {return (a.toBackMult(b));}
	template<class A, class B> inline static A mkBackMult(const A& a, const B& b) {return (a.mkBackMult(b));}


	template<class A> inline static typename ExCo<A>::OUTER_TYPE mkouterprod(const A& a, const A& b){return a.mkouterprod(b);}


    template<class A, class B> inline static bool isLT(const A& a, const B& b) {return a < b;}
    template<class A, class B> inline static bool isLE(const A& a, const B& b) {return a <= b;}
    template<class A, class B> inline static bool isGT(const A& a, const B& b) {return a > b;}
    template<class A, class B> inline static bool isGE(const A& a, const B& b) {return a >= b;}
    template<class A, class B> inline static bool isEQ(const A& a, const B& b) {return a == b;}
    template<class A, class B> inline static bool isNQ(const A& a, const B& b) {return a != b;}

    template<class A, class B> inline static bool comp_lti(const A& a, const B& b) {return b > a;}
    template<class A, class B> inline static bool comp_lei(const A& a, const B& b) {return b >= a;}
    template<class A, class B> inline static bool comp_gti(const A& a, const B& b) {return b < a;}
    template<class A, class B> inline static bool comp_gei(const A& a, const B& b) {return b <= a;}
    template<class A, class B> inline static bool comp_eqi(const A& a, const B& b) {return b == a;}
    template<class A, class B> inline static bool comp_nqi(const A& a, const B& b) {return b != a;}

    template<class A, class B> inline static bool comp_ltj(const A& a, const B& b) {return ExCo<A>::isLT(a,b);}
    template<class A, class B> inline static bool comp_lej(const A& a, const B& b) {return ExCo<A>::isLE(a,b);}
    template<class A, class B> inline static bool comp_gtj(const A& a, const B& b) {return ExCo<A>::isGT(a,b);}
    template<class A, class B> inline static bool comp_gej(const A& a, const B& b) {return ExCo<A>::isGE(a,b);}
    template<class A, class B> inline static bool comp_eqj(const A& a, const B& b) {return ExCo<A>::isEQ(a,b);}
    template<class A, class B> inline static bool comp_nqj(const A& a, const B& b) {return ExCo<A>::isNQ(a,b);}

	template<class F, class I> static typename ExCo<F>::template RETT<I>::TYPE comp3(F f, I f_in) {return f(f_in);}


	template<class A> inline static SETCMP_enum setcmp(const A &a, const A &b) {return a.setcmp(b);}

	template<class A> inline static ERRCODE save(const A& what, FILE *f) {return what.save(f);}
	template<class A> inline static ERRCODE load(A& what, FILE *f, unsigned int lenght) {return what.load(f, lenght);}
	template<class A> inline static ERRCODE save_ISA_pod(const A& what, FILE *f) {return (fwrite(&what,sizeof(A),1,f) == 1) ? 0 : 1;}
	template<class A> inline static ERRCODE load_ISA_pod(A& what, FILE *f, unsigned int lenght) {return (fread(&what,sizeof(A),1,f) == 1) ? 0 : 1 ;}

	template<class A> inline static A& toMemmove(A& a, A& o){return a.toMemmove(o);}
	template<class A> inline static A& toMemswap(A& a, A& o){return a.toMemswap(o);}
	template<class A> inline static A& toMemswap_hasmove(A& a, A& o){A tmp; tmp.toMemmove(o); o.toMemmove(a);return a.toMemmove(tmp);}

	template<class A> inline static A& toMemfree(A& a){return a.toMemfree();return a;}
	template<class A> inline static A& toMemmove_ISA_pod(A& a, A& o){return (a = o);}
	template<class A> inline static void toMemfree_ISA_pod(A& a){}
	template<class A> inline static double getWeight(const A& a){return a.getWeight; }


	template<class A> inline static double pdist(const A& a, const A& b){return a.pdist(b);}
	template<class A> inline static double norm(const A& a){return a.norm(); }
	template<class A> inline static double pnorm(const A& a){return a.pnorm(); }
	template<class A> inline static double lognorm(const A& a){return a.lognorm(); }

	template<class A> inline static A& toRightFlood(A& a){return a.toRightFlood();}
	template<class A> inline static A& toLeftFlood(A& a){return a.toLeftFlood();}
	template<class A> inline static A& toBitInverse(A& a){return a.toBitInverse();}


	template<class A> inline static void show(const A& a, char*&b , int level){return a.show(b,level);}
	template<class A> inline static void show(const A& a, FILE*f , int level){return a.show(f,level);}
	template<class A> inline static string type_tostring(const A& a){return  a.type_tostring();}


	template<class F, class A> inline static void compose(F &func, A &a){return func(a);}
	template<class F, class A> inline static void compose(F &func, const A &a){return func(a);}
	template<class F, class A> inline static void compose(const F &func, A &a){return func(a);}

	template<class F, class A, class B> inline static void compose(F &func, A &a, B &b){return func(a,b);}
	template<class F, class A, class B> inline static void compose(F &func, A &a, const B &b){return func(a,b);}
	template<class F, class A, class B> inline static void compose(F &func, const A &a, const B &b){return func(a,b);}
	template<class F, class A, class B> inline static void compose(const F &func, A &a, const B &b){return func(a,b);}
	template<class F, class A, class B> inline static void compose(const F &func, A &a, B &b){return func(a,b);}

	template<class A> inline static bool isValid(const A &a){return a.isValid();}
	template<class A> inline static bool isValid_cmpchk(const A &a){return ExCo<A>::isValid(a);}
	template<class A> inline static A& toInvalid(A &a){return a.toInvalid();}


    template<class A, class B> inline static A& mkInnerProd(A& a, const B& b, const B& c){return a.mkInnerProd(b,c);}
    template<class A, class B> inline static A& mkOuterProd(A& a, const B& b, const B& c){return a.mkOuterProd(b,c);}
    template<class A, class B, class C> inline static A& mkInnerProd(A& a, const B& b, const C& c){return a.mkInnerProd(b,c);}
    template<class A, class B, class C> inline static A& mkOuterProd(A& a, const B& b, const C& c){return a.mkInnerProd(b,c);}
    template<class A, class B> inline static A& toInnerProd(A& a, const B& b){return a.toInnerProd(b);}
    template<class A, class B> inline static A& toOuterProd(A& a, const B& b){return a.toOuterProd(b);}
    template<class A, class B> inline static A& toInnerProd(A& a, const B& b, const B& c){return a.toOuterProd(b,c);}
    template<class A, class B> inline static A& toOuterProd(A& a, const B& b, const B& c){return a.toOuterProd(b,c);}
    template<class A, class B, class C> inline static A& toInnerProd(A& a, const B& b, const C& c){return a.toOuterProd(b,c);}
    template<class A, class B, class C> inline static A& toOuterProd(A& a, const B& b, const C& c){return a.toOuterProd(b,c);}


    template<class A, class B> inline static A mkBackInnerProd(const A& a, const B& b){return a.mkBackInnerProd(b);}
    template<class A, class B> inline static A mkBackOuterProd(const A& a, const B& b){return a.mkBackOuterProd(b);}
    template<class A, class B> inline static A& toBackInnerProd(A& a, const B& b, const B& c){return a.toBackOuterProd(b,c);}
    template<class A, class B> inline static A& toBackOuterProd(A& a, const B& b, const B& c){return a.toBackOuterProd(b,c);}
    template<class A, class B, class C> inline static A& toBackInnerProd(A& a, const B& b, const C& c){return a.toBackOuterProd(b,c);}
    template<class A, class B, class C> inline static A& toBackOuterProd(A& a, const B& b, const C& c){return a.toBackOuterProd(b,c);}
};

template< >
class ExFn<false>{
public:

	template<class A> inline static bool isZero(const A& what) {return ExCo<A>::isZero(what);}
	template<class A> inline static bool isOne(const A& what) {return ExCo<A>::isOne(what);}
	template<class A> inline static A& toZero(A& what) {return ExCo<A>::toZero(what);}
	template<class A> inline static A& toOne(A& what) {return ExCo<A>::toOne(what);}
	template<class A> inline static A& toRand(A& what) {ExCo<A>::toRand(what); return(what);}
    template<class A, class B> inline static A& toConvert(A& a, const B& b){return ExCo<A>::toConvert(a, b);}

	template<class A> inline static A mkZero(const A& what);
    template<class A> inline static A mkOne(const A& what);
    template<class A> inline static A mkRand(const A& what);


    template<class A> inline static A& toUndefined(A& what) {return ExCo<A>::toUndefined(what);}
	template<class A> inline static A& toMin(A& what) {return ExCo<A>::toMinimum(what);}
	template<class A> inline static A& toMax(A& what) {return ExCo<A>::toMaximum(what);}
	template<class A> inline static A& toMin(A& what, const A& other) {ExCo<A>::mkMinimum(what,other);return what;}
	template<class A> inline static A& toMax(A& what, const A& other) {ExCo<A>::mkMaximum(what,other);return what;}
	// fundamental builders

	// builders

	template<class A> inline static A mknegative(const A& a) {return ExCo<A>::mknegative(a);}
	template<class A> inline static A mkInverse(const A& a) {return ExCo<A>::mkInverse(a);}

    template<class A> inline static A mkAbs(const A& a) {return ExCo<A>::mkAbs(a);}
    template<class A> inline static A mkAbsoft(const A& a, const double& b) {return ExCo<A>::mkAbsoft(a,b);}
    template<class A> inline static A mkAbhard(const A& a, const double& b) {return ExCo<A>::mkAbhard(a,b);}
    template<class A> inline static A mkAbsInverse(const A& a) {return ExCo<A>::mkAbsInverse(a);}
    template<class A> inline static A mkAbsoftInverse(const A& a, const double& b) {return ExCo<A>::mkAbsoftInverse(a,b);}
    template<class A> inline static A mkAbhardInverse(const A& a, const double& b) {return ExCo<A>::mkAbhardInverse(a,b);}


   	template<class A> inline static Tuple<typename ExCo<A>::INDEX_TYPE> mkOrdering(const A& a){return ExCo<A>::mkOrdering(a);}
	template<class A> inline static int32_t mkMagnitude(const A& a){return ExCo<A>::mkMagnitude(a);}
	template<class A> inline static A mkSquare(const A& a) {return ExCo<A>::mkSquare(a);}
	template<class A> inline static A mkPow(const A& what, const double pow) {return ExCo<A>::mkPow(what,pow);}
	template<class A> inline static A mkPowInt(const A& what, const int pow) {return ExCo<A>::mkPowInt(what,pow);}
	template<class A> inline static A mkPowInvInt(const A& what, const int pow) {return ExCo<A>::mkPowInvInt(what,pow);}
	template<class A> inline static typename ExCo<A>::TRJU_TYPE mkTrju(const A& a){return A(a);}
	template<class A> inline static A mkTrjuProd(const A& a){return ExCo<A>::mkTrjuProd(a);}
    template<class A> inline static A& toTrju(A& a){return ExCo<A>::mkTrjuProd(a);}
    template<class A, class B> static void wrDeterminant(const A& a, B& b){ExCo<A>::wrDeterminant(a,b);}
    template<class A, class B> static void wrTrace(const A& a, B& b){ExCo<A>::wrTrace(a,b);}

	// assigns

	template<class A> inline static A& toNegative(A& a) {return ExCo<A>::toNegative(a);}
	template<class A> inline static A& toAbs(A& a) {return ExCo<A>::toAbs(a);}
	template<class A> inline static A& toAbsoft(A& a, const double& b) {return ExCo<A>::toAbsoft(a,b);}
	template<class A> inline static A& toAbhard(A& a, const double& b) {return ExCo<A>::toAbhard(a,b);}
	template<class A> inline static A& toInverse(A& a) {return ExCo<A>::toInverse(a);}
	template<class A> inline static A& toSquare(A& a) {return ExCo<A>::toSquare(a);}
	template<class A> inline static A& toSqrt(A& a) {return ExCo<A>::toSqrt(a);}
	template<class A> inline static A& tointpow(A& what, const int pow) {return ExCo<A>::tointpow(what,pow);}

	template<class A> inline static	A& toClone(A& a, const A& b) {return ExCo<A>::toClone(a,b);}
	template<class A> inline static	A& toAdd(A& a, const A& b) {return ExCo<A>::toAdd(a,b);}
	template<class A> inline static A& toSubt(A& a, const A& b) {return ExCo<A>::toSubt(a,b);}
	template<class A> inline static A& toMult(A& a, const A& b) {return ExCo<A>::toMult(a,b);}
	template<class A> inline static A& toDivi(A& a, const A& b) {return ExCo<A>::toDivi(a,b);}
	template<class A,class B> inline static	A& toClone(A& a, const B& b) {return ExCo<A>::toClone(a,b);}

	template<class A,class B> inline static	A& toAdd(A& a, const B& b) {return ExCo<A>::toAdd(a,b);}
	template<class A,class B> inline static A& toSubt(A& a, const B& b) {return ExCo<A>::toSubt(a,b);}
	template<class A,class B> inline static A& toMult(A& a, const B& b) {return ExCo<A>::toMult(a,b);}
	template<class A,class B> inline static A& toDivi(A& a, const B& b) {return ExCo<A>::toDivi(a,b);}

	template<class A> inline static	A& toAdd(A& a, const A& b, const A& c);
	template<class A> inline static A& toSubt(A& a, const A& b, const A& c);
	template<class A> inline static A& toMult(A& a, const A& b, const A& c);
	template<class A> inline static A& toDivi(A& a, const A& b, const A& c);
	template<class A,class B> inline static	A& toAdd(A& a, const B& b, const B& c);
	template<class A,class B> inline static A& toSubt(A& a, const B& b, const B& c);
	template<class A,class B> inline static A& toMult(A& a, const B& b, const B& c);
	template<class A,class B> inline static A& toDivi(A& a, const B& b, const B& c);
    template<class A,class B,class C> inline static	A& toAdd(A& a, const B& b, const C& c);
	template<class A,class B,class C> inline static A& toSubt(A& a, const B& b, const C& c);
	template<class A,class B,class C> inline static A& toMult(A& a, const B& b, const C& c);
	template<class A,class B,class C> inline static A& toDivi(A& a, const B& b, const C& c);

    template<class A> inline static A& addMult(A& a, const A& b, const A& c);
	template<class A> inline static A& subtMult(A& a, const A& b, const A& c);
	template<class A,class B> inline static A& addMult(A& a, const B& b, const B& c);
	template<class A,class B> inline static A& subtMult(A& a, const B& b, const B& c);
    template<class A,class B,class C> inline static	A& addMult(A& a, const B& b, const C& c);
	template<class A,class B,class C> inline static A& subtMult(A& a, const B& b, const C& c);




//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PLUS_TYPE mkAdd(const A& a, const B& b) {return ExCo<A>::mkAdd(a,b);}
//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::MINU_TYPE mkSubt(const A& a, const B& b) {return ExCo<A>::mkSubt(a,b);}
//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkMult(const A& a, const B& b) {return ExCo<A>::mkMult(a,b); }
//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::DIVI_TYPE mkDivi(const A& a, const B& b) {return ExCo<A>::mkDivi(a,b);}
	//template<class A> inline static	A mkAdd(const A& a, const A& b)  {return ExCo<A>::mkAdd(a,b);}
	//template<class A> inline static A mkSubt(const A& a, const A& b)  {return ExCo<A>::mkSubt(a,b);}
	template<class A, class B> inline static A mkadd_alt(const A& a, const B& b) {return ExCo<A>::mkAdd(a, b);}
	template<class A, class B> inline static A mksubt_alt(const A& a, const B& b) {return ExCo<A>::mkSubt(a, b);}
	template<class A, class B> inline static A mkmult_alt(const A& a, const B& b) {return ExCo<A>::mkMult(a, b);}
	template<class A, class B> inline static A mkdivi_alt(const A& a, const B& b) {return ExCo<A>::mkDivi(a, b);}



//		template<class A> inline static A mkMult(const A& a, const A& b) {return ExFn< Exlisten_mkTrJu<decltype(b * a), decltype(b * a)& (ExCo<decltype(b * a)>:  >::mkmult_trju(a,b);}

//	template<class A, class B> inline static STDRETTYPE2<B,A>::PROD_TYPE mkMult(const A& a, const B& b) {return ExFn< Exlisten_mkmult<B, typename STDRETTYPE2<B,A>::PROD_TYPE (ExCo<B>::SAFETYPE::*)(const A& a)const>  >::mkmult_mkTrJu(a,b);}
	template<class A> inline static A mkAdd(const A& a, const A& b) {return ExFn< Exlisten_toadd< A , A& (ExCo<A>::SAFETYPE::*)(const A&) >::ans >::mkadd_alt(a,b);}
	template<class A> inline static A mkSubt(const A& a, const A& b) {return ExFn< Exlisten_tosubt< A , A& (ExCo<A>::SAFETYPE::*)(const A&) >::ans >::mksubt_alt(a,b);}
	template<class A> inline static A mkMult(const A& a, const A& b) {return ExFn< Exlisten_tomult< A , A& (ExCo<A>::SAFETYPE::*)(const A&) >::ans >::mkmult_alt(a,b);}
	template<class A> inline static A mkDivi(const A& a, const A& b) {return ExFn< Exlisten_todivi< A , A& (ExCo<A>::SAFETYPE::*)(const A&) >::ans >::mkdivi_alt(a,b);}

//  template<class A> inline static A mkDivi(const A& a, const A& b) {return ExFn< Exlisten_toDivi< A , A& (ExCo<A>::SAFETYPE::*)(const A&) >::ans >::mkdivi_alt(a,b);}

	template<class A, class B> inline static A mkAdd(const A& a, const B& b)  {return ExFn< Exlisten_toadd< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mkadd_alt(a,b);}
	template<class A, class B> inline static A mkSubt(const A& a, const B& b) {return ExFn< Exlisten_tosubt< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mksubt_alt(a,b);}
	template<class A, class B> inline static A mkMult(const A& a, const B& b) {return ExFn< Exlisten_tomult< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mkmult_alt(a,b);}
	template<class A, class B> inline static A mkDivi(const A& a, const B& b) {return ExFn< Exlisten_todivi< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mkdivi_alt(a,b);}

	template<class A> inline static typename ExCo<A>::OUTER_TYPE mkouterprod(const A& a, const A& b){return ExFn< Exlisten_mkouterprod< A , typename ExCo<A>::OUTER_TYPE (ExCo<A>::SAFETYPE::*)(const A&) const >::ans  >::mkouterprod(a,b);}


	template<class A, class B> inline static A& toBackMult(A& a, const B& b);
	template<class A, class B> inline static A mkBackMult(const A& a, const B& b);

	template<class A> inline static SETCMP_enum setcmp(const A &a, const A &b) {return ExCo<A>::setcmp(a,b);}

	template<class A> inline static ERRCODE save(const A& what, FILE *f) {return ExFn< ExCo<A>::IsPOD >::save_ISA_pod(what,f);}
	template<class A> inline static ERRCODE load(A& what, FILE *f, unsigned int lenght) { return ExFn< ExCo<A>::IsPOD >::load_ISA_pod(what,f,lenght);}
	template<class A> inline static ERRCODE save_ISA_pod(const A& what, FILE *f) {return ExCo<A>::save(what,f);}
	template<class A> inline static ERRCODE load_ISA_pod(A& what, FILE *f, unsigned int lenght) { return ExCo<A>::load(what,f,lenght);}

	template<class A> inline static A& toMemswap(A& a, A& o){return ExFn< Exlisten_toMemmove<A, A& (ExCo<A>::SAFETYPE::*)( A &) >::ans >::toMemswap_hasmove(a,o);}
	template<class A> inline static A& toMemswap_hasmove(A& a, A& o){return ExCo<A>::toMemswap(a,o);}

	template<class A> inline static A& toMemmove(A& a, A& o){return ExFn< ExCo<A>::IsPOD >::toMemmove_ISA_pod(a,o);}
	template<class A> inline static A& toMemfree(A& a){return ExFn< ExCo<A>::IsPOD >::toMemfree_ISA_pod(a);}
	template<class A> inline static A& toMemmove_ISA_pod(A& a, A& o){return ExCo<A>::toMemmove(a,o);}
	template<class A> inline static A& toMemfree_ISA_pod(A& a){return ExCo<A>::toMemfree(a);}

	template<class A> inline static double getWeight(const A& a){return 1.0f; }

	template<class A> inline static double pdist(const A& a, const A& b){return ExCo<A>::pdist(a,b);}
	template<class A> inline static double norm(const A& a){return ExCo<A>::norm(a);}
	template<class A> inline static double pnorm(const A& a){return ExCo<A>::pnorm(a);}
	template<class A> inline static double lognorm(const A& a){return ExCo<A>::lognorm(a);}

	template<class A> inline static A& toRightFlood(A& a){return ExCo<A>::toRightFlood(a);}
	template<class A> inline static A& toLeftFlood(A& a){return ExCo<A>::toLeftFlood(a);}
	template<class A> inline static A& toBitInverse(A& a){return ExCo<A>::toBitInverse(a);}


	template<class A> inline static void show(const A& a, char* &b , int level){return ExCo<A>::show(a,b,level);}
	template<class A> inline static void show(const A& a, FILE*f , int level){return ExCo<A>::show(a,f,level);}
	template<class A> inline static string type_tostring(const A& a){return ExCo<A>::type_tostring(a);}
	template<class F, class A, class B, unsigned int dims> inline static void compose(F &func, DataGrid<A,dims> &a, const DataGrid<B,dims> &b);

	template<class F, class A, class B, int DA, int DB> inline static void compose(F &func, DataGrid<A,DA> &a, const DataGrid<B,DB> &b);


	template<class A, class B> inline static bool isLT(const A& a, const B& b) {return( ExFn< Exlisten_gt<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_lti(a,b));}
	template<class A, class B> inline static bool isLE(const A& a, const B& b) {return( ExFn< Exlisten_ge<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_lei(a,b));}
	template<class A, class B> inline static bool isGT(const A& a, const B& b) {return( ExFn< Exlisten_lt<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_gti(a,b));}
	template<class A, class B> inline static bool isGE(const A& a, const B& b) {return( ExFn< Exlisten_le<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_gei(a,b));}
	template<class A, class B> inline static bool isEQ(const A& a, const B& b) {return( ExFn< Exlisten_eq<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_eqi(a,b));}
	template<class A, class B> inline static bool isNQ(const A& a, const B& b) {return( ExFn< Exlisten_nq<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_nqi(a,b));}

	template<class A, class B> inline static bool comp_lti(const A& a, const B& b) {return( ExFn< Exlisten_lt<ExCo<A>, bool (ExCo<B>::SAFETYPE::*)(const A&, const B&)const >::ans >::comp_ltj(a,b));}
	template<class A, class B> inline static bool comp_lei(const A& a, const B& b) {return( ExFn< Exlisten_le<ExCo<A>, bool (ExCo<B>::SAFETYPE::*)(const A&, const B&)const >::ans >::comp_lej(a,b));}
	template<class A, class B> inline static bool comp_gti(const A& a, const B& b) {return( ExFn< Exlisten_gt<ExCo<A>, bool (ExCo<B>::SAFETYPE::*)(const A&, const B&)const >::ans >::comp_gtj(a,b));}
	template<class A, class B> inline static bool comp_gei(const A& a, const B& b) {return( ExFn< Exlisten_ge<ExCo<A>, bool (ExCo<B>::SAFETYPE::*)(const A&, const B&)const >::ans >::comp_gej(a,b));}
	template<class A, class B> inline static bool comp_eqi(const A& a, const B& b) {return( ExFn< Exlisten_eq<ExCo<A>, bool (ExCo<B>::SAFETYPE::*)(const A&, const B&)const >::ans >::comp_eqj(a,b));}
	template<class A, class B> inline static bool comp_nqi(const A& a, const B& b) {return( ExFn< Exlisten_nq<ExCo<A>, bool (ExCo<B>::SAFETYPE::*)(const A&, const B&)const >::ans >::comp_nqj(a,b));}


	template<class A, class B> inline static bool comp_ltj(const A& a, const B& b) {return a < b;}
	template<class A, class B> inline static bool comp_lej(const A& a, const B& b) {return a <= b;}
	template<class A, class B> inline static bool comp_gtj(const A& a, const B& b) {return a >b;}
	template<class A, class B> inline static bool comp_gej(const A& a, const B& b) {return a >=b;}
	template<class A, class B> inline static bool comp_eqj(const A& a, const B& b) {return a ==b;}
	template<class A, class B> inline static bool comp_nqj(const A& a, const B& b) {return a !=b;}

	template<class A> inline static bool isValid(const A &a){return ExFn< Exlisten_isValid< ExCo<A> , bool (ExCo<A>::SAFETYPE::*)(const A&) >::ans >::isValid_cmpchk(a);}
	template<class A> inline static A& toInvalid(A &a){return ExFn< Exlisten_toInvalid< ExCo<A> , A& (ExCo<A>::SAFETYPE::*)() >::ans >::toInvalid(a);}
	template<class A> inline static bool isValid_cmpchk(const A &a){return true;}

	template<class F, class I, unsigned int SIZE> static Tuple< typename ExCo<F>::template RETT<I>::TYPE , SIZE > comp3(F f, Tuple<I, SIZE> f_in);
	template<class A> inline static typename ExCo<A>::GAUS_TYPE mkgaussstat(const A& a,double & weight){return ExCo<A>::mkgaussstat(a,weight);}

	template<class A> inline static typename ExCo<A>::INDEX_TYPE& getIndex(A& a){return ExCo<A>::getIndex(a);}
	template<class A> inline static typename ExCo<A>::INDEX_TYPE getIndex(const A& a){return ExCo<A>::getIndex(a);}
	template<class A> inline static typename ExCo<A>::INDEX_TYPE& getIndex(A*& a){return a;}
	template<class A> inline static typename ExCo<A>::INDEX_TYPE getIndex(const A*& a){return a;}

	template<class A> inline static typename ExCo<A>::REAL_TYPE mkrealproj(const A& a){return ExCo<A>::mkrealproj(a);}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkimmaproj(const A& a){return ExCo<A>::mkimmaproj(a);}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkjmmaproj(const A& a){return ExCo<A>::mkjmmaproj(a);}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkkmmaproj(const A& a){return ExCo<A>::mkkmmaproj(a);}

    template<class A> inline static A mkInnerProd(const A& a){return a * a;}
    template<class A> inline static A mkOuterProd(const A& a){return a * a;}
    template<class A, class B> inline static A mkInnerProd(const A& a, const B& b){return a * b;}
    template<class A, class B> inline static A mkOuterProd(const A& a, const B& b){return a * b;}
    template<class A, class B> inline static A& mkInnerProd(A& a, const B& b, const B& c){return a = b*c;}
    template<class A, class B> inline static A& mkOuterProd(A& a, const B& b, const B& c){return a.mkOuterProd(b,c);}
    template<class A, class B, class C> inline static A& mkInnerProd(A& a, const B& b, const C& c){return a.mkInnerProd(b,c);}
    template<class A, class B, class C> inline static A& mkOuterProd(A& a, const B& b, const C& c){return a.mkInnerProd(b,c);}
    template<class A, class B> inline static A& toInnerProd(A& a, const B& b, const B& c){return a = b *c;}
    template<class A, class B> inline static A& toOuterProd(A& a, const B& b, const B& c){return a = b *c;}
    template<class A, class B, class C> inline static A& toInnerProd(A& a, const B& b, const C& c){return a = b *c;}
    template<class A, class B, class C> inline static A& toOuterProd(A& a, const B& b, const C& c){return a = b *c;}

    template<class A, class B> inline static A mkBackInnerProd(const A& a, const B& b);
    template<class A, class B> inline static A mkBackOuterProd(const A& a, const B& b);
    template<class A, class B> inline static A& mkBackInnerProd(A& a, const B& b, const B& c){return a.mkBackInnerProd(b,c);}
    template<class A, class B> inline static A& mkBackOuterProd(A& a, const B& b, const B& c){return a.mkBackOuterProd(b,c);}
    template<class A, class B, class C> inline static A& mkBackInnerProd(A& a, const B& b, const C& c){return a.mkBackInnerProd(b,c);}
    template<class A, class B, class C> inline static A& mkBackOuterProd(A& a, const B& b, const C& c){return a.mkBackInnerProd(b,c);}
    template<class A, class B> inline static A& toBackInnerProd(A& a, const B& b, const B& c){return a.toBackOuterProd(b,c);}
    template<class A, class B> inline static A& toBackOuterProd(A& a, const B& b, const B& c){return a.toBackOuterProd(b,c);}
    template<class A, class B, class C> inline static A& toBackInnerProd(A& a, const B& b, const C& c){return a.toBackOuterProd(b,c);}
    template<class A, class B, class C> inline static A& toBackOuterProd(A& a, const B& b, const C& c){return a.toBackOuterProd(b,c);}

};

template<class A, class B>
class ExType{
    typedef A INNER_TYPE;
    typedef A OUTER_TYPE;
    typedef A MULT_TYPE;
    typedef A ADD_TYPE;
};

class ThreadLog_instance{
public:
    FILE* f;
    char* target;
    char* buffer;
    ThreadLog_instance(FILE* _f, char* _target, char* _buffer):f(_f),target(_target),buffer(_buffer){}
    operator char*(){return target;}
    ~ThreadLog_instance(){fprintf(f, "%s", buffer);fflush(f);}
};
class ThreadLog{
public:
    char buffer[1024];
    FILE* f;
    uint32_t nbscope;
    ThreadLog_instance operator()(){return ThreadLog_instance(f, buffer + nbscope, buffer);}
};



extern thread_local LFHPrimitive::ThreadLog thread_log;



	class ExOp{
	public:
		// class Properties
		template<class A> static void bitreverse(A& a){ExCo<A>::bitreverse(a);}
		template<class A> static void bytereverse(A& a){ExCo<A>::bytereverse(a);}
		//template<class A> static void next(A& a){ExOp_body<A>::next(a);}
		template<class A> static unsigned char upperbound_pow_of_2(const A& a){return ExCo<A>::upperbound_pow_of_2(a); }

		//template<class F> static typename ExCo<F>::DEFRETT comp2(F f, typename ExCo<F>:: f_in) {return f(f_in);}


		//		template<class F, class I> static auto comp3(F f, I f_in) -> int {return ExFn< Exlisten_mainoper< typename ExCo<F>::SAFETYPE, typename ExCo<F>::template RETT<I>::TYPE (ExCo<F>::SAFETYPE::*)(I b) >::ans  >::comp3(f,f_in);}
				template<class F, class I> static typename ExCo<F>::template RETT<I>::TYPE comp3(F f, I f_in) {return ExFn< Exlisten_mainoper< typename ExCo<F>::SAFETYPE, typename ExCo<F>::template RETT<I>::TYPE (ExCo<F>::SAFETYPE::*)(I b) >::ans  >::comp3(f,f_in);}

		//		template<class F> static int comp3(F f, double f_in) {return ExFn< true  >::comp3(f,f_in);}

		template<class F, class I> static typename ExCo<F>::template RETT<I>::TYPE comp2(F f, I f_in) {return f(f_in);}

	//	template<class F, class A, class B> static typename void apply(F& f, A& a, B& b){ }

		template<class O, class I> static void comp(O& o, void (*fun)(O&, const I&), const I & in){fun(o, in);}
		template<class O, class I> static void comp(O& o, O (*fun)(const I&), const I& in){o = fun(in);}
		template<class O, class I> static void comp(O& o, void (I::*fun)(O&), const I& in){in.fun(o);}
		template<class O, class I> static void comp(O& o, O (I::*fun)(), const I& in){o = in.fun();}
		template<class O, class I, class F> static void comp(O& o, void (F::*fun)(O&, const I&), const I& in){fun(o, in);}
		template<class O, class I, class F> static void comp(O& o, O (F::*fun)(const I&), const I& in){o = fun(in);}

        template<class O> static O& execMemfnc(O& object, typename ExCo<O>::FUNCTION_TYPE fnc){return ExCo<O>::execMemfnc(object, fnc);}

		template<class O, class I, class C, class D> static void comp(C& c, void (*fun)(O&, const I&), const D& in){
			I itmp = (I)in;
			O tmp;
			fun(tmp, itmp);
			c = (C) tmp;
			}
		template<class O, class I, class F, class C, class D> static void comp(C& c, void (F::*fun)(O&, const I&), const D& in){
			I itmp = (I)in;
			O tmp;
			fun(tmp, itmp);
			c = (C) tmp;
			}
		template<class O, class I, int size> static void comp(Tuple<O,size> &c, void (*fun)(O&, const I&), const Tuple<I, size>& in){
			for(unsigned int i =0; i < size;i++) fun(c[i], in[i]);
			}
		template<class O, class I, class F, int size> static void comp(Tuple<O,size> &c, void (F::*fun)(O&, const I&), const Tuple<I, size>& in){
			for(unsigned int i =0; i < size;i++) fun(c[i], in[i]);
			}
		template<class O, class I, unsigned int nbdim> static void comp(DataGrid<O,nbdim> &c, void (*fun)(O&, const I&), const DataGrid<I, nbdim>& in){
			c.setSizes(in.dims);
			for(unsigned int i = in.totsize(); i != 0xFFFFFFFF;i--) fun(c.data[i], in.data[i]);
			}
		template<class O, class I, class F, unsigned int nbdim> static void comp(DataGrid<O,nbdim> &c, void (F::*fun)(O&, const I&), const DataGrid<I,nbdim>& in){
			c.setSizes(in.dims);
			for(unsigned int i = in.totsize(); i != 0xFFFFFFFF;i--) fun(c.data[i], in.data[i]);
			}

	//	template<class A> static ExCo<A>::TYPE_norm norm_cmp_set(const A&);
	//	template<class A> static double norm_cmp_extract(const A&);

	/*	template<unsigned int query,class A> static bool hasProperties(const A& a){
			if ((query & STRUCTPROP_IS_VALID)&&(!ExOp_flagbody<STRUCTPROP_IS_VALID,A>::hasProperties(a)) )return(false);
			if ((query & STRUCTPROP_HAS_VALID_INVERSE)&&(!ExOp_flagbody<STRUCTPROP_IS_VALID,A>::hasProperties(a)) )return(false);
			if ((query & STRUCTPROP_CAN_COMMUTE)&&(!ExOp_flagbody<STRUCTPROP_CAN_COMMUTE,A>::hasProperties(a)) )return(false);
			return(true);
		}*/


	//	template<class A> static bool isValid(const A& a){return( ExCo<A>::isValid(a));}

		template<class A> static bool isValid(const A& a){return ExFn< Exlisten_isValid<A, bool (ExCo<A>::SAFETYPE::*)() const >::ans >::isValid(a);}

		template<class A> static double getWeight(const A& a){return(ExCo<A>::getWeight(a));}

		template<class A> static void move(const A& a,const A& b){
			if (ExCo<A>::isMobile::ans) a =b;
			else a.moveFrom(b);
		}

		template<class A> static A* new_class(const A & a) {return (new A(a));} // sementic goody!


		//template<class A> static void showonline(const A& val,FILE* out){ExCo<A>::showonline(val,out);}
		//template<class A> static void show(const A& val){ExCo<A>::show(val,stdout);}
		//template<class A> static void showonline(const A& val){ExCo<A>::showonline(val,stdout);}


		template<class A> static SetComparison compare(const A& a, const A& b) {return( SetComparison( ((a < b) ?  30 | 64 | 256 :  0 ) | ((a > b) ? 29 | 32 | 128 :  0) )); }

		template<class A> inline static double dist(const A& a, const A& b){return sqrt(ExFn< Exlisten_pdist<A, double (ExCo<A>::SAFETYPE::*)(const A&)const >::ans >::pdist(a,b));}
		template<class A> inline static double pnorm(const A& a){return ExFn< Exlisten_pnorm<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::pnorm(a);}
		template<class A> inline static double norm(const A& a){return ExFn< Exlisten_pnorm<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::norm(a);}
		template<class A> inline static double lognorm(const A& a){return ExFn< Exlisten_lognorm<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::lognorm(a);}

		template<class A> inline static bool isZero(const A& a){return ExFn< Exlisten_isZero<A, bool (ExCo<A>::SAFETYPE::*)()const >::ans >::isZero(a);}
		template<class A> inline static bool isOne(const A& a){return ExFn< Exlisten_isOne<A, bool (ExCo<A>::SAFETYPE::*)()const >::ans >::isOne(a);}
		template<class A> inline static A& toZero(A& a){return ExFn< Exlisten_toZero<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toZero(a);}
		template<class A> inline static A& toRand(A& a){return ExFn< Exlisten_toRand<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toRand(a);}
		template<class A> inline static A& toOne(A& a){return ExFn< Exlisten_toOne<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toOne(a);}
        template<class A> inline static A mkZero(const A& a) {return ExFn< Exlisten_mkZero<A, A& (ExCo<A>::SAFETYPE::*)()const >::ans >::mkZero(a);}
        template<class A> inline static A mkOne(const A& a) {return ExFn< Exlisten_mkOne<A, A& (ExCo<A>::SAFETYPE::*)()const >::ans >::mkOne(a);}
        template<class A> inline static A mkRand(const A& a) {return ExFn< Exlisten_mkRand<A, A& (ExCo<A>::SAFETYPE::*)()const >::ans >::mkRand(a);}

        template<class A> inline static A& toConvert(A& a, const A& b){return ExFn< Exlisten_toClone<A, A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toClone(a,b);} // same class, just clone
        template<class A, class B> inline static A& toConvert(A& a, const B& b){return ExFn< Exlisten_toConvert<A, A& (ExCo<A>::SAFETYPE::*)(const B& )const >::ans >::toConvert(a,b);}

		template<class A> inline static A& toUndefined(A& a){return ExFn< Exlisten_toUndefined<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toUndefined(a);}
		template<class A> inline static A& toMin(A& a){return ExFn< Exlisten_toMin<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toMin(a);}
		template<class A> inline static A& toMax(A& a){return ExFn< Exlisten_toMax<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toMax(a);}
		template<class A> inline static A& toMin(A& a, const A& other){return ExFn< Exlisten_toMin<A, A& (ExCo<A>::SAFETYPE::*)(const A& ) >::ans >::toMin(a,other);}
		template<class A> inline static A& toMax(A& a, const A& other){return ExFn< Exlisten_toMax<A, A& (ExCo<A>::SAFETYPE::*)(const A& ) >::ans >::toMax(a,other);}

		template<class A> inline static A mkone(){A fout; ExOp::toOne(fout); return fout;}
		template<class A> inline static A mkzero(){A fout; ExOp::toZero(fout); return fout;}

		template<class A> static A& toNegative(A& a){return ExFn< Exlisten_toNegative<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toNegative(a);}
		template<class A> static A& toAbs(A& a){return ExFn< Exlisten_toAbs<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toAbs(a);}
		template<class A> static A& toAbsoft(A& a,const double& b){return ExFn< Exlisten_toAbsoft<A, A& (ExCo<A>::SAFETYPE::*)(const double&) >::ans >::toAbsoft(a,b);}
		template<class A> static A& toAbhard(A& a,const double& b){return ExFn< Exlisten_toAbhard<A, A& (ExCo<A>::SAFETYPE::*)(const double&) >::ans >::toAbhard(a,b);}
		template<class A> static A& toInverse(A& a){return ExFn< Exlisten_toInverse<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toInverse(a);}
		template<class A> static A& toSquare(A& a){return ExFn< Exlisten_toSquare<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toSquare(a);}
		template<class A> static A& toSqrt(A& a){return ExFn< Exlisten_toSqrt<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toSqrt(a);}
		template<class A> static A& tointpow(A& a, const int pow){return ExFn< Exlisten_tointpow<A, A& (ExCo<A>::SAFETYPE::*)(const int) >::ans >::tointpow(a,pow);}

		template<class A> inline static A& toSubt(A& a); // synonymous with toNegative
		template<class A> inline static A& toDivi(A& a); // synonymous with toInverse

		template<class A> inline static A& toClone(A& a, const A& b){return ExFn< Exlisten_toClone<A, A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toClone(a,b);}
		template<class A> inline static A& toAdd(A& a, const A& b){return ExFn< Exlisten_toAdd<A, A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toAdd(a,b);}
		template<class A> inline static A& toSubt(A& a, const A& b){return ExFn< Exlisten_toSubt<A, A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toSubt(a,b);}
		template<class A> inline static A& toMult(A& a, const A& b){return ExFn< Exlisten_toMult<A, A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toMult(a,b);}
		template<class A> inline static A& toDivi(A& a, const A& b){return ExFn< Exlisten_toDivi<A, A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toDivi(a,b);}
		template<class A, class B> inline static A& toClone(A& a, const B& b){return ExFn< Exlisten_toClone<A, A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toClone(a,b);}
		template<class A, class B> inline static A& toAdd(A& a, const B& b){return ExFn< Exlisten_toAdd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toAdd(a,b);}
		template<class A, class B> inline static A& toSubt(A& a, const B& b){return ExFn< Exlisten_toSubt<A, A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toSubt(a,b);}
		template<class A, class B> inline static A& toMult(A& a, const B& b){return ExFn< Exlisten_toMult<A, A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toMult(a,b);}
		template<class A, class B> inline static A& toDivi(A& a, const B& b){return ExFn< Exlisten_toDivi<A, A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toDivi(a,b);}

		template<class A> inline static A& toAdd(A& a, const A& b, const A& c){return ExFn< Exlisten_toAdd<A, A& (ExCo<A>::SAFETYPE::*)(const A& b, const A& c) >::ans >::toAdd(a,b,c);}
		template<class A> inline static A& toSubt(A& a, const A& b, const A& c){return ExFn< Exlisten_toSubt<A, A& (ExCo<A>::SAFETYPE::*)(const A& b, const A& c) >::ans >::toSubt(a,b,c);}
		template<class A> inline static A& toMult(A& a, const A& b, const A& c){return ExFn< Exlisten_toMult<A, A& (ExCo<A>::SAFETYPE::*)(const A& b, const A& c) >::ans >::toMult(a,b,c);}
		template<class A> inline static A& toDivi(A& a, const A& b, const A& c){return ExFn< Exlisten_toDivi<A, A& (ExCo<A>::SAFETYPE::*)(const A& b, const A& c) >::ans >::toDivi(a,b,c);}
		template<class A, class B> inline static A& toAdd(A& a, const B& b, const B& c){return ExFn< Exlisten_toAdd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toAdd(a,b,c);}
		template<class A, class B> inline static A& toSubt(A& a, const B& b, const B& c){return ExFn< Exlisten_toSubt<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toSubt(a,b,c);}
		template<class A, class B> inline static A& toMult(A& a, const B& b, const B& c){return ExFn< Exlisten_toMult<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toMult(a,b,c);}
		template<class A, class B> inline static A& toDivi(A& a, const B& b, const B& c){return ExFn< Exlisten_toDivi<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toDivi(a,b,c);}
		template<class A, class B, class C> inline static A& toAdd(A& a, const B& b, const C& c){return ExFn< Exlisten_toAdd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toAdd(a,b,c);}
		template<class A, class B, class C> inline static A& toSubt(A& a, const B& b, const C& c){return ExFn< Exlisten_toSubt<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toSubt(a,b,c);}
		template<class A, class B, class C> inline static A& toMult(A& a, const B& b, const C& c){return ExFn< Exlisten_toMult<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toMult(a,b,c);}
		template<class A, class B, class C> inline static A& toDivi(A& a, const B& b, const C& c){return ExFn< Exlisten_toDivi<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toDivi(a,b,c);}

		template<class A> inline static A mkInnerProd(const A& a){return ExFn< Exlisten_mkInnerProd<A, A& (ExCo<A>::SAFETYPE::*)()const >::ans >::mkInnerProd(a);}
		template<class A> inline static A mkOuterProd(const A& a){return ExFn< Exlisten_mkOuterProd<A, A& (ExCo<A>::SAFETYPE::*)()const >::ans >::mkOuterProd(a);}
		template<class A, class B> inline static A mkInnerProd(const A& a, const B& b){return ExFn< Exlisten_mkInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkInnerProd(a,b);}
		template<class A, class B> inline static A mkOuterProd(const A& a, const B& b){return ExFn< Exlisten_mkOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkOuterProd(a,b);}
		template<class A, class B, class C> inline static A& mkInnerProd(A& a, const B& b, const C& c){return ExFn< Exlisten_toInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::mkInnerProd(a,b,c);}
		template<class A, class B, class C> inline static A& mkOuterProd(A& a, const B& b, const C& c){return ExFn< Exlisten_toOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::mkOuterProd(a,b,c);}
		template<class A, class B> inline static A& toInnerProd(A& a, const B& b){return ExFn< Exlisten_toInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const A& b, const B& c) >::ans >::toInnerProd(a,b);}
		template<class A, class B> inline static A& toOuterProd(A& a, const B& b){return ExFn< Exlisten_toOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const A& b, const B& c) >::ans >::toOuterProd(a,b);}
		template<class A, class B> inline static A& toInnerProd(A& a, const B& b, const B& c){return ExFn< Exlisten_toInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toInnerProd(a,b,c);}
		template<class A, class B> inline static A& toOuterProd(A& a, const B& b, const B& c){return ExFn< Exlisten_toOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toOuterProd(a,b,c);}
		template<class A, class B, class C> inline static A& toInnerProd(A& a, const B& b, const C& c){return ExFn< Exlisten_toInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toInnerProd(a,b,c);}
		template<class A, class B, class C> inline static A& toOuterProd(A& a, const B& b, const C& c){return ExFn< Exlisten_toOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toOuterProd(a,b,c);}

		template<class A, class B> inline static A mkBackInnerProd(const A& a, const B& b){return ExFn< Exlisten_toBackInnerProd<A, A (ExCo<A>::SAFETYPE::*)(const A& b, const B& c) >::ans >::mkBackInnerProd(a,b);}
		template<class A, class B> inline static A mkBackOuterProd(const A& a, const B& b){return ExFn< Exlisten_toBackOuterProd<A, A (ExCo<A>::SAFETYPE::*)(const A& b, const B& c) >::ans >::mkBackOuterProd(a,b);}

		template<class A, class B> inline static A& toBackInnerProd(A& a, const B& b, const B& c){return ExFn< Exlisten_toBackInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toBackInnerProd(a,b,c);}
		template<class A, class B> inline static A& toBackOuterProd(A& a, const B& b, const B& c){return ExFn< Exlisten_toBackOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const B& c) >::ans >::toBackOuterProd(a,b,c);}
		template<class A, class B, class C> inline static A& toBackInnerProd(A& a, const B& b, const C& c){return ExFn< Exlisten_toBackInnerProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toBackInnerProd(a,b,c);}
		template<class A, class B, class C> inline static A& toBackOuterProd(A& a, const B& b, const C& c){return ExFn< Exlisten_toBackOuterProd<A, A& (ExCo<A>::SAFETYPE::*)(const B& b, const C& c) >::ans >::toBackOuterProd(a,b,c);}



//		// special check for multiplication by weights
//		template<class A> inline static typename STDRETTYPE2<A,Weight>::PROD_TYPE mkMult(const A& a, const Weight& b){return ExFn< Exlisten_mkmult<A, typename STDRETTYPE2<A,B>::PROD_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkMult(a,b);}
//		template<class A> inline static A& toMult(A& a, const Weight& b){return ExFn< Exlisten_toMult<A, const A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toMult(a,b);}

		template<class A, class B> inline static A& toBackMult(A& a, const B& b){return ExFn< Exlisten_toBackMult<A, A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toBackMult(a,b);}
//		template<class A, class B> inline static const A& toRmult(A& a, const B& b){return ExFn< Exlisten_toRmult<A, void (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toRmult(a,b);}

		template<class A, class B> inline static A mkBackMult(const A& a, const B& b){return ExFn< Exlisten_mkBackMult<A, A& (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkBackMult(a,b);}
//		template<class A> inline static typename ExCo<A>::INNER_TYPE mkRmult(const A& a, const A& b){return ExFn< Exlisten_mkrmul<A, void (ExCo<A>::SAFETYPE::*)(const A& b)const >::ans >::toRmult(a,b);}
        template<class A> inline static A mkTrjuProd(const A& a){return ExFn< Exlisten_mkTrjuProd<A, A (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkTrjuProd(a);}
        template<class A> inline static typename ExCo<A>::TRJU_TYPE mkTrju(const A& a){return ExFn< Exlisten_mkTrju<A, typename ExCo<A>::TRJU_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkTrju(a);}
        template<class A> inline static A& toTrju(A& a){return ExFn< Exlisten_toTrju<A, A& (ExCo<A>::SAFETYPE::*)()  >::ans >::toTrju(a);}
        template<class A, class B> static void wrDeterminant(const A& a, B& b){return ExFn< Exlisten_wrDeterminant<A, void (ExCo<A>::SAFETYPE::*)(B&)const  >::ans >::wrDeterminant(a,b);}
        template<class A, class B> static void wrTrace(const A& a, B& b){return ExFn< Exlisten_wrTrace<A, void (ExCo<A>::SAFETYPE::*)(B&)const  >::ans >::wrTrace(a,b);}


		template<class A> inline static A mknegative(const A& a){return ExFn< Exlisten_mksubt<A, A (ExCo<A>::SAFETYPE::*)()const >::ans >::mknegative(a);}
		template<class A> inline static A mkInverse(const A& a){return ExFn< Exlisten_mkInverse<A, A (ExCo<A>::SAFETYPE::*)()const >::ans >::mkInverse(a);}

		template<class A> inline static A mkAbs(const A& a){return ExFn< Exlisten_mkAbs<A, A (ExCo<A>::SAFETYPE::*)()const >::ans >::mkAbs(a);}
		template<class A> inline static A mkAbsoft(const A& a,double b){return ExFn< Exlisten_mkAbsoft<A, A (ExCo<A>::SAFETYPE::*)(const double&)const >::ans >::mkAbsoft(a,b);}
		template<class A> inline static A mkAbhard(const A& a,double b){return ExFn< Exlisten_mkAbhard<A, A (ExCo<A>::SAFETYPE::*)(const double&)const >::ans >::mkAbhard(a,b);}
		template<class A> inline static A mkAbsInverse(const A& a){return ExFn< Exlisten_mkAbsInverse<A, A (ExCo<A>::SAFETYPE::*)()const >::ans >::mkAbsInverse(a);}
		template<class A> inline static A mkAbsoftInverse(const A& a,double b){return ExFn< Exlisten_mkAbsoftInverse<A, A (ExCo<A>::SAFETYPE::*)(const double&)const >::ans >::mkAbsoftInverse(a,b);}
		template<class A> inline static A mkAbhardInverse(const A& a,double b){return ExFn< Exlisten_mkAbhardInverse<A, A (ExCo<A>::SAFETYPE::*)(const double&)const >::ans >::mkAbhardInverse(a,b);}


		template<class A> inline static A mkSquare(const A& a){return ExFn< Exlisten_mkSquare<A, A (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkSquare(a);}
		template<class A> inline static A mkSqrt(const A& a){return ExFn< Exlisten_mkSqrt<A, A (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkSqrt(a);}

		template<class A> inline static typename ExCo<A>::REAL_TYPE mkrealproj(const A& a){return ExFn< Exlisten_mkrealproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkrealproj(a);}
		template<class A> inline static typename ExCo<A>::REAL_TYPE mkimmaproj(const A& a){return ExFn< Exlisten_mkimmaproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkimmaproj(a);}
		template<class A> inline static typename ExCo<A>::REAL_TYPE mkjmmaproj(const A& a){return ExFn< Exlisten_mkjmmaproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkjmmaproj(a);}
		template<class A> inline static typename ExCo<A>::REAL_TYPE mkkmmaproj(const A& a){return ExFn< Exlisten_mkkmmaproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkkmmaproj(a);}
        template<class A> inline static typename ExCo<A>::GAUS_TYPE mkgaussstat(const A& a,double weight = 1.0f){return ExFn< Exlisten_mkgaussstat<A, typename ExCo<A>::GAUS_TYPE (ExCo<A>::SAFETYPE::*)(double&)const  >::ans >::mkgaussstat(a,weight);}

        template<class A> inline static typename ExCo<A>::INDEX_TYPE& getIndex(A& a){return ExFn< Exlisten_getIndex<A, typename ExCo<A>::INDEX_TYPE& (ExCo<A>::SAFETYPE::*)() >::ans >::getIndex(a);}
        template<class A> inline static typename ExCo<A>::INDEX_TYPE getIndex(const A& a){return ExFn< Exlisten_getIndex<A, typename ExCo<A>::INDEX_TYPE (ExCo<A>::SAFETYPE::*)() const >::ans>::getIndex(a);}

        //template<class A> inline static typename ExCo<A>::INDEX_TYPE& getIndex(A*& a){return ExFn< Exlisten_getIndex<A, typename ExCo<A>::INDEX_TYPE& (ExCo<A>::SAFETYPE::*)() >::ans >::getIndex(a);}
        //template<class A> inline static typename ExCo<A>::INDEX_TYPE getIndex(const A*& a){return ExFn< Exlisten_getIndex<A, typename ExCo<A>::INDEX_TYPE (ExCo<A>::SAFETYPE::*)() const >::ans>::getIndex(a);}



		template<class A> inline static A mkAdd(const A& a, const A& b){return ExFn< Exlisten_mkadd<A, A (ExCo<A>::SAFETYPE::*)(const A& b)const >::ans >::mkAdd(a,b);}
		template<class A> inline static A mkSubt(const A& a, const A& b){return ExFn< Exlisten_mksubt<A, A (ExCo<A>::SAFETYPE::*)(const A& b)const >::ans >::mkSubt(a,b);}
		template<class A> inline static A mkDivi(const A& a, const A& b){return ExFn< Exlisten_mkdivi<A, A (ExCo<A>::SAFETYPE::*)(const A& b)const >::ans >::mkDivi(a,b);}

		template<class A> inline static A mkMult(const A& a, const A& b) {return ExFn< Exlisten_mkmult<A, A (ExCo<A>::SAFETYPE::*)(const A& b)const >::ans >::mkMult(a,b);}

		template<class A, class B> inline static typename STDRETTYPE2<A,B>::PLUS_TYPE mkAdd(const A& a, const B& b){return ExFn< Exlisten_mkadd<A, typename STDRETTYPE2<A,B>::PLUS_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkAdd(a,b);}
		template<class A, class B> inline static typename STDRETTYPE2<A,B>::MINU_TYPE mkSubt(const A& a, const B& b){return ExFn< Exlisten_mksubt<A, typename STDRETTYPE2<A,B>::MINU_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkSubt(a,b);}
		template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkMult(const A& a, const B& b){return ExFn< Exlisten_mkmult<A, typename STDRETTYPE2<A,B>::PROD_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkMult(a,b);}
		template<class A, class B> inline static typename STDRETTYPE2<A,B>::DIVI_TYPE mkDivi(const A& a, const B& b){return ExFn< Exlisten_mkdivi<A, typename STDRETTYPE2<A,B>::DIVI_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkDivi(a,b);}

		template<class A> inline static Tuple<typename ExCo<A>::INDEX_TYPE> mkOrdering(const A& a){return ExFn< Exlisten_mkOrdering<A, Tuple<typename ExCo<A>::INDEX_TYPE> (ExCo<A>::SAFETYPE::*)()const >::ans >::mkOrdering(a);}
		template<class A> inline static int32_t mkMagnitude(const A& a){return ExFn< Exlisten_mkMagnitude<A, int32_t (ExCo<A>::SAFETYPE::*)()const >::ans >::mkMagnitude(a);}
		template<class A> inline static A mkPow(const A& a, const double pow){return ExFn< Exlisten_mkPow<A, A (ExCo<A>::SAFETYPE::*)(const double)const >::ans >::mkPow(a,pow);}
		template<class A> inline static A mkPowInt(const A& a, const int pow){return ExFn< Exlisten_mkPowInt<A, A (ExCo<A>::SAFETYPE::*)(const int)const >::ans >::mkPowInt(a,pow);}
		template<class A> inline static A mkPowInvInt(const A& a, const int pow){return ExFn< Exlisten_mkPowInvInt<A, A (ExCo<A>::SAFETYPE::*)(const int)const >::ans >::mkPowInvInt(a,pow);}
		template<class A> inline static SETCMP_enum setcmp(const A& a, const A& b) {return ExFn< Exlisten_setcmp<A, void (ExCo<A>::SAFETYPE::*)(const A&)const >::ans >::setcmp(a,b);}

// comparisons

        template<class A, class B> inline static bool isLT(const A &a, const B &b);
        template<class A, class B> inline static bool isGT(const A &a, const B &b);
        template<class A, class B> inline static bool isLE(const A &a, const B &b);
        template<class A, class B> inline static bool isGE(const A &a, const B &b);
        template<class A, class B> inline static bool isEQ(const A &a, const B &b);
        template<class A, class B> inline static bool isNQ(const A &a, const B &b);



        // free any associated memory to an object, but does not free the object itself. If a pointer is given, identical to delete() ,then sets to NULL
        template<class A> inline static A& toMemfree(A& a){return ExFn< Exlisten_toMemfree<A, A& (ExCo<A>::SAFETYPE::*)( A &) >::ans >::toMemfree(a);}
		// same as operator =, but the source looses the ownership of its associated memory
		template<class A> inline static A& toMemmove(A& a, A& o){return ExFn< Exlisten_toMemmove<A, A& (ExCo<A>::SAFETYPE::*)( A &) >::ans >::toMemmove(a,o);}
		// same as operator =, but the source looses the ownership of its associated memory
		template<class A> inline static A& toMemswap(A& a, A& o){return ExFn< Exlisten_toMemswap<A, A& (ExCo<A>::SAFETYPE::*)( A &) >::ans >::toMemswap(a,o);}

		template<class A> inline static double getWeight(A& a){return  ExFn< Exlisten_getWeight<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::getWeight(a);}
		template<class A> inline static ERRCODE save(const A& w, FILE *f) {return ExFn< Exlisten_save<A, ERRCODE ( ExCo<A>::SAFETYPE::*)(FILE*) const >::ans >::save(w,f);}
		template<class A> inline static ERRCODE load(A& w, FILE *f, unsigned int lenght = 0xFFFFFFFF) {return ExFn< Exlisten_load<A, ERRCODE ( ExCo<A>::SAFETYPE::*)(FILE*, unsigned int)>::ans >::load(w,f,lenght);}


		template<class A> inline static A& toRightFlood(A& a){return ExFn< Exlisten_toRightFlood<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toRightFlood(a);}
		template<class A> inline static A& toLeftFlood(A& a){return ExFn< Exlisten_toLeftFlood<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toLeftFlood(a);}
		template<class A> inline static A& toBitInverse(A& a){return ExFn< Exlisten_toBitInverse<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toBitInverse(a);}

		// THE SHOW COMMAND! level 1: '\n' forbidden, level 2: '\t' forbidden, level 3: '[;]' forbidden, level 4: '(,)' forbidden,
		template<class A> inline static void show11(const A& a, char* b_out = thread_log(), int level = 0){return  ExFn< Exlisten_show<A, void (ExCo<A>::SAFETYPE::*)(char* &f, int level)const >::ans >::show(a,b_out,level);}
		template<class A> inline static void show(const A& a, char* &b_out, int level = 0){return  ExFn< Exlisten_show<A, void (ExCo<A>::SAFETYPE::*)(char* &f, int level)const >::ans >::show(a,b_out,level);}
		template<class A> inline static void show(const A& a, FILE* f_out = stdout, int level = 0){return  ExFn< Exlisten_show<A, void (ExCo<A>::SAFETYPE::*)(FILE* f, int level)const >::ans >::show(a,f_out,level);}
		template<class A> inline static void showf(const A& a, FILE* f_out = stdout, int level = 0){ExFn< Exlisten_show<A, void (ExCo<A>::SAFETYPE::*)(FILE* f, int level)const >::ans >::show(a,f_out,level); fflush(f_out);}
        template<class A> inline static string type_tostring(const A& a){return  ExFn< Exlisten_type_tostring<A, string (ExCo<A>::SAFETYPE::*)() const >::ans >::type_tostring(a);}

		// 1 arg operator
		template<class F, class A> inline static void comp(F &func, A &a){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A&) const>::ans)
			 >::compose(func,a);}
		template<class F, class A> inline static void comp(const F &func, A &a){ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A&) const>::ans)>::compose(func,a);}
		template<class F, class A> inline static void comp(F &func, const A &a){ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A&)>::ans)>::compose(func,a);}

		// 2  arg operator
		template<class F, class A, class B> inline static void comp(F &func, A &a, B &b){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , B&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A& , const B&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , B&) const>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const>::ans)
			 >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(const F &func, A &a, B &b){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , B&) const >::ans)
			    ||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const>::ans) >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(F &func, A &a, const B &b){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&)>::ans)
			    ||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const>::ans)
			    ||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A& , const B&) >::ans) >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(const F &func, A &a, const B &b){
			ExFn< Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const >::ans >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(F &func, const A &a, const B &b){
			ExFn< Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A& , const B&)>::ans >::compose(func,a,b);}

						// a * b^H
	//	template<class A> inline static void mkcomult(const A& a,const A& b)
};




template<class C>
class Magscaled{ // this is value * exp(exponent)
public:
    C value;
    double exponent;
    C operator()()const;
    Magscaled<C>()=default;
    Magscaled<C>(const C& value, double exponent);
    Magscaled<C>& operator*=(const Magscaled<C>&);
    Magscaled<C>& operator/=(const Magscaled<C>&);
    Magscaled<C> operator*(const Magscaled<C>&)const;
    Magscaled<C> operator/(const Magscaled<C>&)const;

    Magscaled<C>& operator*=(const C&);
    Magscaled<C>& operator/=(const C&);
    Magscaled<C> operator*(const C&)const;
    Magscaled<C> operator/(const C&)const;

    Magscaled<C>& monotonicAdd(const Magscaled<C>& smaller); // this += smaller, assumes this > smaller > 0 or this < smaller < 0
    Magscaled<C>& operator+=(const Magscaled<C>& other);

/*
    Magscaled<C>& operator+=(const Magscaled<C>&);
    Magscaled<C>& operator-=(const Magscaled<C>&);
    Magscaled<C> operator+(const Magscaled<C>&)const;
    Magscaled<C> operator-(const Magscaled<C>&)const;*/
    void show(FILE*f =stdout, int level=0)const;
};


LFH_GOLD	double gamma(double x); // robust even when x -> 0
LFH_GOLD	double lngamma(double x); // robust even when x -> 0
LFH_GOLD    void gamma(Magscaled<double> &fout, double x);


LFH_GOLD	double lngammaexp(double x); // log(gamma(exp(x)))
LFH_GOLD	double lngammadiff(double x, double y); // log(gamma(x+y)) - log(gamma(x)) robust even when y/x -> 0 or x -> 0
LFH_GOLD	double lngammasum(double x, double pivot); // log(gamma(pivot-x)) + log(gamma(pivot+x))


LFH_GOLD	double lngammachoose(double x, double y); // log(gamma(y+1)) - log(gamma(x+1)) - log(gamma(y-x))
LFH_GOLD	double lngammadiffexp(double x, double y); // log(gamma(exp(x)+ exp(y))) - log(gamma(exp(x)))
double erf_inverse(const double x);

double lnbeta(double a, double b); // robust even when x -> 0
double lnbetaexp(double a, double b); // log(beta(exp(a),exp(b)))





double log_exp_m1(double); // computes log(e^x - 1) , where x >= 0
double log_exp_m1(double); // computes log(e^x - 1) , where x >= 0
double d_log_exp_m1_dx(double); // computes d log(e^x - 1) dx , where x >= 0
//double log_1_mexp(double x){return(log_exp_m1(x) - x);} // computes log(1 - e^x) = log(e^x - 1) - x

// double tanh_scaled(); // = (tanh(x) + 1) / 2

double d_tanh(double);
double d2_tanh(double);


double log_e_x_p1(const double x); // Log(1 + exp(x)) (aka -log(logistic(-x)))    [[inverse of log_e_x_m1, its derivative is the logistic function]]
double log_e_x_m1(const double x); // (domain x > 0) Log(-1 + exp(x))             [[inverse of log_e_x_p1]]
double log_1_me_x(const double x); // (domain x < 0) Log(1 - exp(x))
double log_1_me_x_v2(const double x);

double log_1_me_nx(const double x); // (domain x > 0) Log(1 - exp(-x))
double d_log_1_me_nx_dx(const double x); // (domain x < 0) 1.0 / (1 - exp(-x))

// aliases
double logistic(double); // computes 1.0 / (1.0 + e^(-x)), works with infinity
double d_logistic_dx(double); // computes e^(-x) / (1.0 + e^(-x))^2
double logit(const double x);// log( x) - log(1 - x) , works with infinity
// double d_logit_dx(const double x);// (1-2x) / (x*(1 - x)) , works with infinity
double logit(const Magscaled<double> x);// log( x) - log(1 - x)
double logit_signif(const double x, const double e);// log( x * exp(e)) - log(1 - x * exp(e)) = -log(exp(-e - ln(x)) - 1.0) = -log_e_x_m1(-e - ln(x))

double powerlogistic(double x, double y); // domain:y>0 computes (1.0 + e^(-x)) ^ -y, y>0
double powerlogit(double x, double y); // domain:y>0  computes -log(x^(1/-y) - 1.0) = -log( e^(log(x)/-y) - 1.0)

// log((exp(log_x) + exp(log_y))/2)
double log_average(const double log_x, const double log_y);



double log_avg_exp(const double x,const double y); // log(0.5 *exp(x) + 0.5 *exp(y)) = log(1.0  + exp(y-x)) + x + log(0.5)



double logtransformed_sum(const double x,const double y); // log(exp(x) + exp(y))

double log_x_p1(const double x); // log(x+1)
double log_x_1plus_overx(const double x); // log(1 + x) / x
double log_e_x_1plus_overe_x(const double x); // log(1 + e^x) / e^x
double e_x_1minus_overx(const double x); // (exp(x) - 1) / x
double e_x_1minus(const double x); // exp(x) - 1
double one_over_e_x_1minus(const double x); // (domain x != 0), 1 / (exp(x) - 1)

LFH_GOLD	double Prob_NBdistrib(uint32_t x, double r, double p);

//log (2 - sum a_i) + log (sum exp(logerf(v) + log(w) - log(f))) + log(f)

LFH_GOLD	double randExponential(); // sample exponential distribution mean mean 1
LFH_GOLD	double randUniform();
uint32_t randSample(double* prob, uint32_t length); // sample integer with given probabilities

double sampleGaussian();
double sampleGaussian(double mu, double sigma);

double cumNormal(double x);
double lncumNormal(double x); double lncumNormalmmm(double x);
void cumNormal(Magscaled<double> &fout, double x);


double logitcumNormal(double x);
double d_lncumNormal_dx(double x);
double d2_lncumNormal_dx2(double x);

double logitPval_to_Zscore(double x);

double logP_to_Normal(const double mu, const double sigmasquare, const double logP = -randExponential());
double logP_to_Exponential(const double scale, const double logP = -randExponential());
double logP_to_Gamma(const double scale, const double shape, const double logP =-randExponential());

myHashmap<uint32_t, double, defaultHashFnc<uint32_t> > getHypergeomericMasses(uint32_t n, uint32_t K, uint32_t N, double probfraction = 2.0); // gets the probability masses that accounts for >= probfraction (or all if no threshold provided)
void getHypergeomericCover(Tuple<uint32_t> &fout, uint32_t n, uint32_t K, uint32_t N); // sample hypergeometric distribution at equidistant point



//double HypergeomericLogitPvalBound(uint32_t k, uint32_t n, uint32_t K, uint32_t N); // gets a bound negative logitpval, and upperbound for positive logit val
double HypergeomericLogitPval(double k, uint32_t n, uint32_t K, uint32_t N); // gets the probability masses that accounts for >= probfraction (or all if no threshold provided)
double HypergeomericLogitPval(uint32_t k, uint32_t n, uint32_t K, uint32_t N, double frac = 0.5); // gets the probability masses that accounts for >= probfraction (or all if no threshold provided)
double HypergeomericNormalLogitPval(double dev,double sigma, uint32_t n, uint32_t K, uint32_t N, bool show=false); // logit[ 0.5 + 0.5*\sum_k=0^K P(k) * erf( (dev -k) / (sigma *sqrt(2)) ) ]
double HypergeomericNormalPval_slow(double dev,double sigma, uint32_t n, uint32_t K, uint32_t N);

// just a 32bit extension to the exponent...
class ExtremeDouble{
public:
    double value;
    uint32_t ext_exp;
    ExtremeDouble(double value);
    ExtremeDouble& operator*=(const double&);
};
// 0, 1, 1, 1, 17, 31, 691, 10922, 929569, 3202291, 221930581, 9444233042, 56963745931, 29435334228302, 2093660879252671, 344502690252804724, 129848163681107301953, 868320396104950823611, 209390615747646519456961, 28259319101491102261334882
// 1, 2, 12, 45, 2520, 14175, 935550, 42567525, 10216206000, 97692469875, 18561569276250, 2143861251406875, 34806217964017500, 48076088562799171875, 9086380738369043484375, 3952575621190533915703125, 3920955016221009644377500000, 68739242628124575327993046875

#include "core.hpp"

} // end of namespace LFHPrimitive


#endif


