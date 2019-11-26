/*
 * core.hpp
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


 // namespace LFHPrimitive{

#include <inttypes.h>
#define LFH_typename_check(TyPeNaMe) \
template <class C>  int doexist_##TyPeNaMe( typename C::TyPeNaMe const * ); \
template <typename T>   char doexist_##TyPeNaMe( ... ); \
template<int hasfunc, class C> class isdef_##TyPeNaMe{public:typedef typename C::TyPeNaMe TYPE;}; \
template<class C> class isdef_##TyPeNaMe<sizeof(char), C>{public:typedef C TYPE; }

#define LFH_typename_check_dep(TyPeNaMe, WoUlDbEtYpE) \
template <class C>  int doexist_##TyPeNaMe( typename C::TyPeNaMe const * ); \
template <typename T>   char doexist_##TyPeNaMe( ... ); \
template<int hasfunc, class C> class isdef_##TyPeNaMe{public:typedef typename C::TyPeNaMe TYPE;}; \
template<class C> class isdef_##TyPeNaMe<sizeof(char), C>{public:typedef WoUlDbEtYpE TYPE; }

LFH_typename_check(DEF_TYPE);
LFH_typename_check(NEG_TYPE);
LFH_typename_check(TRJU_TYPE);
LFH_typename_check(LMUL_TYPE);
LFH_typename_check(INNER_TYPE);
LFH_typename_check(OUTER_TYPE);
LFH_typename_check(VECTOR_TYPE);
LFH_typename_check(COMPARATOR_TYPE);

LFH_typename_check(DERIVATIVE_TYPE);

LFH_typename_check_dep(SAFETYPE, ExCo<C>);
LFH_typename_check_dep(WEIGHT_TYPE, double);
LFH_typename_check(REAL_TYPE);
LFH_typename_check_dep(COMPLEX_TYPE, Complex<C>);
LFH_typename_check_dep(QUATERNION_TYPE, Quaternion<C>);
LFH_typename_check_dep(ITERATOR_TYPE, void);
LFH_typename_check_dep(GAUS_TYPE, LFHCONCAT2(WeightElem<C,2>) );

LFH_typename_check(INDEX_TYPE);

LFH_typename_check_dep(IS_ENUM, YESNO<false> );

LFH_typename_check_dep(IS_COMMUTATIVE, YESNO<false> );

LFH_typename_check_dep(IS_OWNED, YESNO<false> );
LFH_typename_check_dep(IS_POD, YESNO<false> );


template<class A, int flag>
class ExCo{
public:
    typedef A TYPE;
    typedef typename isdef_DEF_TYPE<sizeof(doexist_DEF_TYPE<A>(0)),A>::TYPE DEF_TYPE;
    typedef typename isdef_NEG_TYPE<sizeof(doexist_NEG_TYPE<A>(0)),A>::TYPE NEG_TYPE;
    typedef typename isdef_TRJU_TYPE<sizeof(doexist_TRJU_TYPE<A>(0)),A>::TYPE TRJU_TYPE;
    typedef typename isdef_LMUL_TYPE<sizeof(doexist_LMUL_TYPE<A>(0)),A>::TYPE LMUL_TYPE;

    typedef typename isdef_INNER_TYPE<sizeof(doexist_INNER_TYPE<A>(0)),A>::TYPE INNER_TYPE;
    typedef typename isdef_VECTOR_TYPE<sizeof(doexist_VECTOR_TYPE<A>(0)),A>::TYPE VECTOR_TYPE;
    typedef typename isdef_DERIVATIVE_TYPE<sizeof(doexist_DERIVATIVE_TYPE<A>(0)),A>::TYPE DERIVATIVE_TYPE;

    typedef typename isdef_OUTER_TYPE<sizeof(doexist_OUTER_TYPE<A>(0)),A>::TYPE OUTER_TYPE; // outer product result, (increase in the number of dimentions!)


    typedef typename isdef_REAL_TYPE<sizeof(doexist_REAL_TYPE<A>(0)),A>::TYPE REAL_TYPE;
    typedef typename isdef_COMPARATOR_TYPE<sizeof(doexist_COMPARATOR_TYPE<A>(0)),A>::TYPE COMPARATOR_TYPE;


    typedef typename isdef_COMPLEX_TYPE<sizeof(doexist_COMPLEX_TYPE<A>(0)),A>::TYPE COMPLEX_TYPE;
    typedef typename isdef_QUATERNION_TYPE<sizeof(doexist_QUATERNION_TYPE<A>(0)),A>::TYPE QUATERNION_TYPE;
    typedef typename isdef_GAUS_TYPE<sizeof(doexist_GAUS_TYPE<A>(0)),A>::TYPE GAUS_TYPE;

    typedef typename isdef_IS_COMMUTATIVE<sizeof(doexist_IS_COMMUTATIVE<A>(0)),A>::TYPE IS_COMMUTATIVE;
    typedef typename isdef_IS_OWNED<sizeof(doexist_IS_OWNED<A>(0)),A>::TYPE IS_OWNED;
    typedef typename isdef_IS_POD<sizeof(doexist_IS_POD<A>(0)),A>::TYPE IS_POD;

    typedef typename isdef_INDEX_TYPE<sizeof(doexist_INDEX_TYPE<A>(0)),A>::TYPE INDEX_TYPE;
    typedef A& (A::*FUNCTION_TYPE)();
    inline static A& execMemfnc(A& object, A& (A::*fnc)()){return (object.*fnc)();}


    //typedef typename isdef_SAFETYPE<sizeof(doexist_SAFETYPE<A>(0)),A>::TYPE SAFETYPE;
    typedef A SAFETYPE;


    static const bool IsPOD = IS_POD::ans ;
    static const bool NeedsAddLink = IS_OWNED::ans; // containers needs to update addresses in link registers


    template<class B> static inline A& toConversion(A& a, const B& b){return a = (A)b;} // tries implicit conversion if no explicit exists...



    static A mkOne(){A a; return ExOp::toOne(a);}
    static A mkZero(){A a; return ExOp::toZero(a);}
    static A random(){A a; a.toRand(); return(a);}
    static A intPow(const A& val, int power){
        switch(power){
            case -1: return(invert(val));
            case 0: return(mkOne());
            case 1: return(val);
            case 2: return(val*val);
            case 3: return(val*val*val);

        }
    }

    static A invert(const A& val){ return(val);} // needed! if possible
    static A floPow(const A& val, float power){ return(floPow(val,(double)power));} // needed! if possible
    static A floPow(const A& val, double power){ return(val);} // needed! if possible
    static double sign(const A& val);
    template<int v> static A intPow(const  A& what) {return(intPow(v));}
    static double arctan(const A& r,const A& i);

    static void show(const A& val, FILE* out = stdout, int level = 0){val.show(out,level);}

    static bool hasTransferFunction(){return(false);}// has function which move ressources, and leaves the origin ready for destruction

    static double pdist(const A& a, const A& b);
    static double norm(A& a){return ExCo<A>::norm(a);}
    static double norm_squarred(A& a){return ExCo<A>::norm_squarred(a);}

    static ERRCODE save(const A& what, FILE *f) {return (fwrite(&what,sizeof(A),1,f) == 1) ? 0 : 1;}
    static ERRCODE load(A& what, FILE *f, unsigned int lenght){return (1 == fread(&what,sizeof(A),1,f)) ? 0 : 1;}


    // build missing operations from primitives!

    // builds
    static A& toMemmove(A& a, A &b){return a = b;} // unknown, copy
    static A& toMemfree(A& a){return a;} // unknown, ignore
    inline static A mkAbs(const A& a);
    inline static A mknegative(const A& a);
    inline static A mkInverse(const A& a);
    inline static A mkSquare(const A& a);
    inline static A mkPowInt(const A& what, const int pow);
    inline static A mkPowInvInt(const A& what, const int pow); // uses mkInverse and mkPowInt

    inline static A mkTrjuProd(const A& what);
    inline static A mkTrju(const A& what);
    inline static A& toTrju(A& what);
    inline static A mkPosInverse(const A& a);
    inline static A mkNegInverse(const A& a);


    inline static A& getIndex(A& a){return a;}



    // compare

    inline static bool isGT(const A &,const A &);
    inline static bool isGE(const A &,const A &);
    inline static bool isLT(const A &,const A &);
    inline static bool isLE(const A &,const A &);
    inline static bool isEQ(const A &,const A &);
    inline static bool isNQ(const A &,const A &);


    static void tointpow(A& what, const int pow);
    //	static double pdist(const A& a, const A& b){return ExOp::pnorm(a-b);}
};

template<class A>
class ExCo<A,1u>{
public:
    typedef A TYPE;
    typedef ExCo<A> SAFETYPE;
    typedef void COMPLEX_TYPE;
    typedef void QUATERNION_TYPE;
    static const bool IsPOD = true;
    static void show(const A& what, FILE*f, int level){fprintf(f, "enum%i\n", (int32_t)what); }
    static A& toZero(A& what){memset(&what, '\0',sizeof(A)); return what;}
};
/*
template<class A>
	class ExCo<EnumBox<A> >{ // is enum, or unlucky class which needs an SAFETYPE defined
	public:

		typedef A TYPE;
        typedef typename isdef_DEF_TYPE<sizeof(doexist_DEF_TYPE<A>(0)),A>::TYPE DEF_TYPE;
		typedef typename isdef_NEG_TYPE<sizeof(doexist_NEG_TYPE<A>(0)),A>::TYPE NEG_TYPE;
		typedef typename isdef_TRJU_TYPE<sizeof(doexist_TRJU_TYPE<A>(0)),A>::TYPE TRJU_TYPE;
        typedef typename isdef_LMUL_TYPE<sizeof(doexist_LMUL_TYPE<A>(0)),A>::TYPE LMUL_TYPE;

		typedef typename isdef_INNER_TYPE<sizeof(doexist_INNER_TYPE<A>(0)),A>::TYPE INNER_TYPE;
        typedef typename isdef_VECTOR_TYPE<sizeof(doexist_VECTOR_TYPE<A>(0)),A>::TYPE VECTOR_TYPE;
		typedef typename isdef_DERIVATIVE_TYPE<sizeof(doexist_DERIVATIVE_TYPE<A>(0)),A>::TYPE DERIVATIVE_TYPE;

		typedef typename isdef_OUTER_TYPE<sizeof(doexist_OUTER_TYPE<A>(0)),A>::TYPE OUTER_TYPE; // outer product result, (increase in the number of dimentions!)


		typedef typename isdef_REAL_TYPE<sizeof(doexist_REAL_TYPE<A>(0)),A>::TYPE REAL_TYPE;
		typedef typename isdef_COMPLEX_TYPE<sizeof(doexist_COMPLEX_TYPE<A>(0)),A>::TYPE COMPLEX_TYPE;
		typedef typename isdef_QUATERNION_TYPE<sizeof(doexist_QUATERNION_TYPE<A>(0)),A>::TYPE QUATERNION_TYPE;
        typedef typename isdef_GAUS_TYPE<sizeof(doexist_GAUS_TYPE<A>(0)),A>::TYPE GAUS_TYPE;

        typedef typename isdef_IS_COMMUTATIVE<sizeof(doexist_IS_COMMUTATIVE<A>(0)),A>::TYPE IS_COMMUTATIVE;
        typedef typename isdef_IS_OWNED<sizeof(doexist_IS_OWNED<A>(0)),A>::TYPE IS_OWNED;
        typedef typename isdef_IS_POD<sizeof(doexist_IS_POD<A>(0)),A>::TYPE IS_POD;

        typedef ExCo<A> SAFETYPE;



		static const bool IsPOD = IS_POD::ans ;
		static const bool NeedsAddLink = IS_OWNED::ans; // containers needs to update addresses in link registers





		static A one(){A a; a.toOne(); return(a);}
		static A zero(){A a; a.toZero(); return(a);}
		static A random(){A a; a.toRand(); return(a);}
		static A intPow(const A& val, int power){
			switch(power){
				case -1: return(invert(val));
				case 0: return(one());
				case 1: return(val);
				case 2: return(val*val);
				case 3: return(val*val*val);

			}
		}

		static A invert(const A& val){ return(val);} // needed! if possible
		static A floPow(const A& val, float power){ return(floPow(val,(double)power));} // needed! if possible
		static A floPow(const A& val, double power){ return(val);} // needed! if possible
		static double sign(const A& val);
		template<int v> static A intPow(const  A& what) {return(intPow(v));}
		static double arctan(const A& r,const A& i);

		static void show(const A& val, FILE* out = stdout, int level = 0){val.show(out,level);}

		static void hasTransferFunction(){return(false);}// has function which move ressources, and leaves the origin ready for destruction

		static double pdist(const A& a, const A& b);
		static double norm(A& a){return ExCo<A>::norm(a);}
		static double norm_squarred(A& a){return ExCo<A>::norm_squarred(a);}

		static void save(const A& what, FILE *f) {
			//	if (hasfunc)what.save(f);
			//	else
			fwrite(&what,sizeof(A),1,f);
		} //
		static void load(A& what, FILE *f, unsigned int lenght){
			//		if (hasfunc) what.load(f, lenght);
			//	else
			fread(&what,sizeof(A),1,f);
		}

		// build missing operations from primitives!

		// buildstoNegative
		inline static A mkPowInt(const A& what, const int pow);
		inline static A mkPowInvInt(const A& what, const int pow); // uses mkInverse and mkPowInt


		// assignments

		inline static A& toinverse(A& what);
		inline static A& tonegative(A& what);
		inline static void tomult(A& a, const A& b);
		inline static A& toSquare(A& a);

		inline static A mkAdd(const A& a, const A&b);
		inline static A mksub(const A& a, const A&b);
		inline static A mkMult(const A& a, const A&b);
		inline static A mkDivi(const A& a, const A&b);

		// compare

		inline static bool isGT(const A &,const A &);
		inline static bool isGE(const A &,const A &);
		inline static bool isLT(const A &,const A &);
		inline static bool isLE(const A &,const A &);
		inline static bool isEQ(const A &,const A &);
		inline static bool isNQ(const A &,const A &);


		static void tointpow(A& what, const int pow);
		//	static double pdist(const A& a, const A& b){return ExOp::pnorm(a-b);}

	};

*/

// dead-end for object manipulation, cant tell the array size at this point
template<class A>
class ExCo<A*, 0>{
public:
	//typedef YESNO<false> IsPOD;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;
	typedef A* TYPE;
    typedef ExCo<A*> SAFETYPE;

	typedef ExCo<A*> COMPLEX_TYPE;
	typedef ExCo<A*> QUATERNION_TYPE;

	typedef typename isdef_INDEX_TYPE<sizeof(doexist_INDEX_TYPE<A>(0)),A>::TYPE INDEX_TYPE;




	static A * & toRand(A * & a){return a=NULL;}//{if (a != NULL) ExOp::toRand(*a); return a;}
	static A * & toOne(A * & a){return a=NULL;}//{if (a != NULL) ExOp::toOne(*a); return a;}
	static A * & toZero(A * & a){return a=NULL;}//{if (a != NULL) ExOp::toZero(*a); return a;}

    static A* & toMemfree(A* & a){delete(a); return (a =NULL);}
    static A* & toMemmove(A* & a, A* & b){a = b; b =NULL; return a;}
    static A* & toMemswap(A* & a, A* & b){A* tmp = b; b = a; return a =tmp;}
	static A* mkOne(){return( NULL);}
	static A* mkZero(){return( NULL);}
	static A* random(){return(NULL);}

	//static INDEX_TYPE getIndex(A * & a);

	static ERRCODE save(A * & what, FILE *f);
	static ERRCODE load(A * const & what, FILE *f, unsigned int lenght =0);
	static void show(A * const & what, FILE *f,int level);
};

/*
template<class A>
class ExCo<const A*, 0u>{
public:
	//typedef YESNO<false> IsPOD;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;


	static void show(A * const & what, FILE *f,int level);
};*/

// dead-end for object manipulation, cant tell the array size at this point
template< >
class ExCo<char*>{
public:
	//typedef YESNO<false> IsPOD;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;
	typedef char* TYPE;
	typedef ExCo<char*> SAFETYPE;

	typedef char* COMPLEX_TYPE;
	typedef char* QUATERNION_TYPE;
	typedef void* DERIVATIVE_TYPE;

    static char* & toMemfree(char *& a){delete(a); a =NULL; return a;}
    static char* & toMemmove(char*& a, char* & b){a = b; b =NULL; return a;}
    static char* & toMemswap(char*& a, char* & b){char* tmp = b; b = a; return a = tmp;}
	static char* & toRand(char *& a){a = NULL; return a;}
	static char* & toOne(char *& a) {a = NULL; return a;}
	static char* & toZero(char *& a){return a = NULL;}
	static char* mkOne(){return( NULL);}
	static char* mkZero(){return( NULL);}
	static char* random(){return(NULL);}
    static char* toUndefined(char* & a){delete[](a); a = NULL; return a;}

	static ERRCODE save(char* const & what, FILE *f) { fwrite(what,1,strlen(what)+1,f);return 0;}
	static ERRCODE load(char* & what, FILE *f, unsigned int lenght =0) {ERRCODE fout=0; char buffer[65536]; char* tmp = buffer; do {fread(tmp,sizeof(unsigned char), 1,f); } while(*(tmp++) != '\0'); what = new char[(unsigned int)(tmp - buffer)]; memcpy(what, buffer, sizeof(char) * (unsigned int)(tmp - buffer)); return fout;}

	static void show(char* const & what, FILE *f,int level) { if (level == 0) fprintf(f,"%s\n",what); else fprintf(f,"%s",what);}

};


template<class A, unsigned int SIZE>
class ExCo<A[SIZE], 0 >{
	public:
//        typedef typename ExCo<A>::SAFETYPE[SIZE] TYPE;
		typedef typename ExCo<A>::SAFETYPE SAFETYPE;


        typedef typename isdef_DEF_TYPE<sizeof(doexist_DEF_TYPE<A>(0)),A>::TYPE DEF_TYPE[SIZE];
		typedef typename isdef_NEG_TYPE<sizeof(doexist_NEG_TYPE<A>(0)),A>::TYPE NEG_TYPE[SIZE];
        typedef DataGrid<A,2> TRJU_TYPE;
        typedef DataGrid<A,2> LMUL_TYPE;
	//	typedef typename ExCo<A>::DERIVATIVE_TYPE DERIVATIVE_TYPE[SIZE];
        typedef A INNER_TYPE;

		typedef typename ExCo<A>::COMPLEX_TYPE* COMPLEX_TYPE;
		typedef typename ExCo<A>::QUATERNION_TYPE* QUATERNION_TYPE;

        template <class S> class SUBS_INNER{public: typedef Tuple<S, SIZE> TYPE;};

		static const bool IsPOD = ExCo<A>::IsPOD;
		static const bool NeedsAddLink = ExCo<A>::NeedsAddLink;


		static A (&toZero(A (&a)[SIZE]))[SIZE];
		static A (&toOne(A (&a)[SIZE]))[SIZE];
        static A (&toRand(A (&a)[SIZE]))[SIZE];
		static void minimum(A (&a)[SIZE]);
        static void maximum(A (&a)[SIZE]);
        static A (&toMemfree(A (&a)[SIZE]))[SIZE];
        static A (&toMemmove(A (&a)[SIZE], A (&b)[SIZE])) [SIZE];

		static void show(const A (&a)[SIZE], FILE* out=stdout, int level=0);

		static bool isValid(A (&a)[SIZE]);


		inline static A   (&toAdd(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE];
		inline static A  (&toSubt(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE];
		inline static A  (&toMult(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE];
		inline static A  (&toDivi(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE];

        inline static A (&toClone(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE];
        template <class B> inline static A (&toClone(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE];
        template <class B> inline static A (&toConvert(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE];

		template <class B> inline static A  (&toAdd(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE];
		template <class B> inline static A  (&toSubt(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE];
		template <class B> inline static A  (&toMult(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE];
		template <class B> inline static A  (&toDivi(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE];

		inline static double norm(const A (&a)[SIZE]);
		inline static double pnorm(const A (&a)[SIZE]);

		inline static ERRCODE save(const A (&a)[SIZE], FILE *f);
		inline static ERRCODE load(A (&a)[SIZE], FILE *f, unsigned int lenght =0);



		/*
		 static const bool IsPOD = true;
		 static const bool hasSpetoMemmove = false;


		 static double zero(){return(0.0f);}
		 static void zero(double &a){a =0.0f;}
		 static double one(){return(1.0f);}
		 static void one(double &a){a =1.0f;}
		 static double random(){return((((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32));  };
		 static double& toRand(double &a){a = (((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32); return a; };

		 static bool isValid(const double& a){
		 return(a + (-a) == 0);
		 }
		 static bool isnegative(const double& a){
		 return(a< 0.0f);
		 }
		 static double max(){return(DBL_MAX);}
		 static double min(){return(DBL_MIN);}

		 static double invert(const double& val){return(1.0f / val);}
		 static double intPow(const double& val, int power){return(pow(val, (double) power));}
		 static double floPow(const double& val, double power){return(pow(val, power));}
		 static double sign(const double& val){return(val >= 0 ? 1.0f : -1.0f);}
		 static double invintPow(const double& val, int power){return(pow(val, 1.0f / ((double)power)));}

		 static double arctan(const double& r,const double& i){return(atan2(r,i));}

		 static double doubleratio(const double& num,const double& den){return(num/den);}



		 static void show(const double &val, FILE* out=stdout, int level=0){fprintf(out, "%e", val); if (level == 0) fprintf(out, "\n"); }

		 static void save(const double& what, FILE *f) { fwrite(&what,sizeof(double), 1,f);}
		 static void load(double & what, FILE *f) { fread(&what,sizeof(double), 1,f);}
		 static void toMemmove(double& a,double& o){a=o;}

		 static double norm(const double& a){return fabs(a);}
		 static double pnorm(const double& a){return a*a ;}
		 */
	};


template< class R, class A>
class ExCo< R (*)(A), 0 >{ // function that returns A with arg B
    public:
    static const bool IsPOD = false;
    static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

    typedef void COMPLEX_TYPE;
    typedef void QUATERNION_TYPE;

    typedef ExCo< R (*)(A) > SAFETYPE;

    template< class C > class RETT{ // is a simple function: the return type is always the same for any input
    public:
        typedef R TYPE;
        C junkdata;
    };

    R operator()(A a); // fake declaration! (used to tell if input arguments are legal for function


//	R operator()(A&){} // fake declaration!
//	R operator()(const A&){} // fake declaration!
//	R operator()(A) const{} // fake declaration!
//	R operator()(A&) const{} // fake declaration!
//	R operator()(const A&) const{} // fake declaration!

};

template<class A, class B>
class ExCo<pair<A,B> , 0>{
public:
	typedef pair<typename ExCo<A>::TYPE,typename ExCo<B>::TYPE> TYPE;
	typedef ExCo<pair<A,B> >  SAFETYPE;
//	typedef pair< typename ExCo<A>::MAGNITUDE_TYPE, typename ExCo<B>::MAGNITUDE_TYPE> MAGNITUDE_TYPE;
	static const bool NeedsAddLink = ExCo<A>::NeedsAddLink ||  ExCo<B>::NeedsAddLink ; // containers needs to update addresses in link registers
	typedef pair< typename ExCo<A>::COMPLEX_TYPE, typename ExCo<B>::COMPLEX_TYPE>  COMPLEX_TYPE;
	typedef pair< typename ExCo<A>::QUATERNION_TYPE, typename ExCo<B>::QUATERNION_TYPE> QUATERNION_TYPE;
	static const bool IsPOD = ExCo<A>::IsPOD && ExCo<B>::IsPOD;
	inline static bool isZero(const pair<A,B> &a);
	inline static bool isOne(const pair<A,B> &a);
	static pair<A,B>& toZero(pair<A,B> &a);
	static pair<A,B>& toOne(pair<A,B> &a);
	static pair<A,B>& toRand(pair<A,B> &a);
	static bool isValid(const pair<A,B> &a);
	static pair<A,B>& toClone(pair<A,B> &a, const pair<A,B> &b);
	template<class A2, class B2> static pair<A,B>& toClone(pair<A,B> &a, const pair<A2,B2> &b);
	template<class A2, class B2> static pair<A,B>& toConvert(pair<A,B> &a, const pair<A2,B2> &b);
	static pair<A,B>& toMemmove(pair<A,B> &a, pair<A,B> &b);
	static pair<A,B>& toMemfree(pair<A,B> &a);
	static pair<A,B>& toAdd(pair<A,B> &a, const pair<A,B> &b);
	static pair<A,B>& toSubt(pair<A,B> &a, const pair<A,B> &b);
	static pair<A,B>& toMult(pair<A,B> &a, const pair<A,B> &b);
	static pair<A,B>& toDivi(pair<A,B> &a, const pair<A,B> &b);
	static void show(const pair<A,B> &a, FILE* out=stdout, int level=0);
	static string type_tostring(const pair<A,B> &a);
	static ERRCODE save(const pair<A,B> &a, FILE*f);
	static ERRCODE load(pair<A,B> &a, FILE*f, unsigned int l);
};

template<class A>
class ExCo< vector<A>, 0 >{
public:
    typedef ExCo<vector<A> > SAFETYPE;
//	typedef vector< typename ExCo<A>::MAGNITUDE_TYPE> MAGNITUDE_TYPE;
	static const bool NeedsAddLink = ExCo<A>::NeedsAddLink ; // containers needs to update addresses in link registers
	typedef vector< typename ExCo<A>::COMPLEX_TYPE>  COMPLEX_TYPE;
	typedef vector< typename ExCo<A>::QUATERNION_TYPE> QUATERNION_TYPE;

	static const bool IsPOD = false;

	inline static bool isZero(const vector<A> &a);
	inline static bool isOne(const vector<A> &a);
	inline static bool isValid(const vector<A> &a);
	inline static vector<A>& toMemfree(vector<A> &a);

    template<class B> static inline vector<A>& toConvert(vector<A> & a, const B& b){return b.wrStdVector(a);} // tries implicit conversion if no explicit exists...

/*	static void zero(vector<A> &a);
	static void one(vector<A> &a);
	static void toRand(vector<A> &a);

	static void show(const vector<A> &a, FILE* out=stdout, int level=0);
	static string type_tostring(const vector<A> &a);*/



};

template<class A>
class ExCo<const A, 0>{
public:
    typedef const A TYPE;
	typedef typename ExCo<A>::SAFETYPE SAFETYPE;
//	typedef typename ExCo<A>::MAGNITUDE_TYPE MAGNITUDE_TYPE;
	typedef typename ExCo<A>::COMPLEX_TYPE  COMPLEX_TYPE;
	typedef typename ExCo<A>::QUATERNION_TYPE QUATERNION_TYPE;
	static const bool NeedsAddLink = ExCo<A>::NeedsAddLink; // containers needs to update addresses in link registers
	static const bool IsPOD = true;
	static inline void show(const A &a, FILE* out=stdout, int level=0);
	static inline bool isValid(const A& a);
};

// static inline GaussScope<double, double> mkgaussstat(const cLaSs &a, double weight = 1.0f){return (GaussScope<double, double>((double)a, ((double)a)*((double)a),weight) );}


#define PRIMTYPE_CONSTRUCT_FUNCTION(cLaSs) \
static inline cLaSs mkMaximum(){cLaSs a; return toMaximum(a);} \
static inline cLaSs mkMinimum(){cLaSs a; return toMinimum(a);} \
static inline cLaSs mkZero(){cLaSs a; return(toZero(a));} \
static inline cLaSs mkOne(){cLaSs a; return(toOne(a));} \
static inline cLaSs random(){cLaSs a;  return toRand(a);} \
static inline cLaSs epsilon(){cLaSs a; epsilon(a); return(a);} /* minimum value that changes the representation of 1 when added */ \
static inline cLaSs delta(){cLaSs a; delta(a); return(a);} /* minimum value that is larger than 0 */ \
typedef cLaSs DEF_TYPE; \
typedef cLaSs NEG_TYPE; \
typedef cLaSs TRJU_TYPE; \
typedef cLaSs OUTER_TYPE; \
typedef cLaSs LMUL_TYPE; \
typedef cLaSs REAL_TYPE; \
typedef Complex<cLaSs>  COMPLEX_TYPE; \
typedef Quaternion<cLaSs> QUATERNION_TYPE; \
typedef GaussScope<double, double> GAUS_TYPE; \
typedef cLaSs INDEX_TYPE; \
typedef YESNO<true> IS_COMMUTATIVE; \
typedef YESNO<false> IS_OWNED; \
typedef YESNO<true> IS_POD; \
typedef cLaSs VECTOR_TYPE;\
typedef cLaSs (*FUNCTION_TYPE)(cLaSs);\
static inline cLaSs mkrealproj(const cLaSs& a){return a;}\
static inline cLaSs mkTrjuProd(const cLaSs& a){return a * a;}\
static inline cLaSs mkTrju(const cLaSs& a){return a;}\
static inline cLaSs& toTrju(cLaSs& a){return a;}\
static inline cLaSs mkimmaproj(const cLaSs& a){return mkZero();}\
static inline cLaSs mkjmmaproj(const cLaSs& a){return mkZero();}\
static inline cLaSs mkkmmaproj(const cLaSs& a){return mkZero();}\
static inline const cLaSs& getIndex(const cLaSs& a){return a;}\
static inline cLaSs& getIndex(cLaSs& a){return a;}\
static inline cLaSs& toMemfree(cLaSs& a){return a;}\
static inline cLaSs& toMemmove(cLaSs& a, cLaSs& b){return (a=b);}\
static inline cLaSs& toClone(cLaSs& a, const cLaSs& b){return (a=b);}

#define SIMPLE_FUNCTIONS(cLaSs)                                        \
static inline cLaSs mkAdd(const cLaSs &a, const cLaSs &b){return a+b;} \
static inline cLaSs mkSubt(const cLaSs &a, const cLaSs &b){return a-b;} \
static inline cLaSs mkMult(const cLaSs &a, const cLaSs &b){return a*b;} \
static inline cLaSs mkDivi(const cLaSs &a, const cLaSs &b){return a/b;} \
static inline cLaSs& toAdd(cLaSs &a, const cLaSs &b){return (a +=b);} \
static inline cLaSs& toSubt(cLaSs &a, const cLaSs &b){return (a -=b);} \
static inline cLaSs& toMult(cLaSs &a, const cLaSs &b){return (a *=b);} \
static inline cLaSs& toDivi(cLaSs &a, const cLaSs &b){return (a /=b);} \
static inline cLaSs mkSquare(const cLaSs &a){return a*a;} \
static inline cLaSs& toSquare(cLaSs &a){return a *= a;} \
static inline cLaSs mknegative(const cLaSs &a) {return -a;}\
static inline bool isGT(const cLaSs &a, const cLaSs &b){return a > b;} \
static inline bool isGE(const cLaSs &a, const cLaSs &b){return a >= b;}\
static inline bool isLT(const cLaSs &a, const cLaSs &b){return a < b;} \
static inline bool isLE(const cLaSs &a, const cLaSs &b){return a <= b;}\
static inline bool isEQ(const cLaSs &a, const cLaSs &b){return a == b;}\
static inline bool isNQ(const cLaSs &a, const cLaSs &b){return a != b;}\
static inline SETCMP_enum setcmp(const cLaSs &a, const cLaSs &b) {return ((a > b) ? SETCMP_GT : ((a == b) ? SETCMP_EQUAL : SETCMP_LT));}

#define SIMPLE_FUNCTIONS_FLOAT(cLaSs)                                   \
static inline void wrDeterminant(const cLaSs &a, float &b) {b = (float)a;}\
static inline void wrDeterminant(const cLaSs &a, double &b) {b = (double)a;}\
static inline void wrTrace(const cLaSs &a, float &b) {b = (float)a;}\
static inline void wrTrace(const cLaSs &a, double &b) {b = (double)a;}\
static inline cLaSs& toSqrt(cLaSs &a) {return a = mkSqrt(a);}\
static inline cLaSs& toMemswap(cLaSs&a, cLaSs&b){cLaSs tmp = b; b=a;return a=tmp;}

template< >
class ExCo<double>{
public:
	typedef double TYPE;
//	typedef short MAGNITUDE_TYPE;
	typedef void DERIVATIVE_TYPE;
	static const bool IsPOD = true;
	static const bool hasSpetoMemmove = false;

	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
	typedef ExOp SAFETYPE;

	static bool isZero(const double &a){return a == 0.0f;}
	static bool isOne(const double &a){return a == 1.0f;}
	static double& toInfinity(double &a){return a = INFINITY;}
	static double& toNegInfinity(double &a){return a = -INFINITY;}
	static bool isInfinite(double &a){return std::isinf(a);}
	static double& toZero(double &a){return a =0.0f;}
	static double& toOne(double &a){return a =1.0f;}
	static double& toRand(double &a){a = (((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32); return a;  };
	static void epsilon(double &a){a =DBL_EPSILON;} // minimum value that changes the representation of 1 when added
	static void delta(double &a){a = DBL_MIN;}  // minimum value that is larger than 0

	static double& toMaximum(double &a){a = DBL_MAX; return a;}
	static double& toMinimum(double &a){a = -DBL_MAX; return a;}
    static double& toUndefined(double &a){a = std::nan(""); return a;}



    inline static double& execMemfnc(double& object, double (*fnc)(double)){return object = (*fnc)(object);}

	PRIMTYPE_CONSTRUCT_FUNCTION(double)

	inline static double mkInversePlain(double a){return 1.0f / a;}
	inline static double& toinverse(double &a){return (a = 1.0f / a);}
	inline static double& toNegative(double &a){return (a = -a);}
	inline static double& toAbs(double &a){return a = fabs(a);}
	inline static double& tosquarre(double &a){ return(a *= a); }
	//inline static float toAbsoft(const float &a, const double& k){a = fabs(a); if (a < k) a = k; return a;}
	inline static double& toAbhard(double &a, const double& k){a = fabs(a); if (a < k) a = k; return a;}
	inline static double mkAbhard(const double &a, const double& k){double fout = fabs(a); return (fout < k) ? k : fout;}
	inline static double mkAbhardInverse(const double &a, const double& k){double fout = fabs(a); return 1.0 / ((fout < k) ? k : fout);}
    inline static double mkAbsoft(const double &a, const double& k){
        double fout;
        if ((fout = fabs(a)) > k*16) return fout;
        fout /= k;
        return log2(exp2(fout) + exp2(-fout)) * k;
        }
    inline static double& toAbsoft(double &a, double k){
        if ((a = fabs(a)) > k*16) return a;
        a /= k;
        return a = log2(exp2(a) + exp2(-a)) * k;
        }
    inline static double mkAbsoftInverse(const double &a, const double& k){
        float fout;
        if ((fout = fabs(a)) > k*16) return 1.0 / fout;
        fout /= k;
        return 1.0 / (log2(exp2(fout) + exp2(-fout)) * k);
        }
    inline static double& toAbsoftInverse(double &a, double k){
        if ((a = fabs(a)) > k*16) return a = (1.0 / a);
        a /= k;
        return a = 1.0 / (log2(exp2(a) + exp2(-a)) * k);
        }
	inline static double mkSqrt(const double &a){return sqrt(a);}
	inline static double mkInverse(const double &a){return(1.0f / a);}
	inline static double mkPowInt(const double& what, const int power){ return( pow(what, (double) power)); }
	inline static double mkPowInvInt(const double& what, const int power){ return( pow(what, 1.0f / ((double) power))); }


	SIMPLE_FUNCTIONS(double)
    SIMPLE_FUNCTIONS_FLOAT(double)

	static bool isValid(const double& a){
		return(a + (-a) == 0);
	}
	static bool isnegative(const double& a){
		return(a< 0.0f);
	}
	static double max(){return(DBL_MAX);}
	static double min(){return(DBL_MIN);}

	static double invert(const double& val){return(1.0f / val);}
	static double intPow(const double& val, int power){return(pow(val, (double) power));}
	static double floPow(const double& val, double power){return(pow(val, power));}
	static double sign(const double& val){return(val >= 0 ? 1.0f : -1.0f);}
	static double invintPow(const double& val, int power){return(pow(val, 1.0f / ((double)power)));}

	static double arctan(const double& r,const double& i){return(atan2(r,i));}

	static double doubleratio(const double& num,const double& den){return(num/den);}



	static void show(const double &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%e", val); if (level == 0) *(ptr++) = '\n';}
	static void show(const double &val, FILE* out=stdout, int level=0){fprintf(out, "%e", val); if (level == 0) fprintf(out, "\n"); }
	static string type_tostring(const double& a){return string("double");}

	static ERRCODE save(const double& what, FILE *f) {return (fwrite(&what,sizeof(double), 1,f) == 1) ? 0 : 1;}
	static ERRCODE load(double & what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(double), 1,f) == 1 ? 0 : 1;}

	static double lognorm(const double& a){return log(fabs(a));}
	static double norm(const double& a){return fabs(a);}
	static double pnorm(const double& a){return a*a;}


//	static double mkMult(const double& val ,const Weight& w){ return val * w.w;}
//	static double& tomult(double& val ,const Weight& w){ return val *= w.w;}


	template<unsigned int SIZE> static unsigned int getHyperDirection(const double (&a)[SIZE], const double (&b)[SIZE]){ // direction for hyper-ordering of arrays
    unsigned int best =0;

    typename StConstList<16>::UINT expb, expc;

    typename StConstList<64>::UINT xorb = (*((typename StConstList<64>::UINT *) &a)) ^ (*((typename StConstList<64>::UINT *) &b)); // xor and make unsigned
        // if exponent are different, the max eponent wins
        // if exponent are the same


		if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0xFFF0){
			if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0x8000) return(0);
			// most significant is max exponent

			((typename StConstList<16>::UINT *)&xorb)[3] = 0x0010; // most-significant turned on!
			expb = (((typename StConstList<16>::UINT *)&a)[3] > ((typename StConstList<16>::UINT *)&b)[3] ? ((typename StConstList<16>::UINT *)&a)[3] : ((typename StConstList<16>::UINT *)&b)[3]) >> 2;
			// gets greater exponent

		}else{
			if (xorb != 0) expb = (((typename StConstList<16>::UINT *)&a)[3]) >> 2;
			else  expb = 0;
			// gets exponent (same for both numbers)
		}




		for(unsigned int cur=1;cur<SIZE;cur++) {
		typename StConstList<64>::UINT xorc = (*((typename StConstList<64>::UINT *) &(a + cur))) ^ (*((typename StConstList<64>::UINT *) &(b + cur))); // xor and make unsigned
			if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0xFFF0){
				if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0x8000) return(cur);

				expc = (((typename StConstList<16>::UINT *)&(a+cur))[3] > ((typename StConstList<16>::UINT *)&(b+cur))[3] ? ((typename StConstList<16>::UINT *)&(a+cur))[3] : ((typename StConstList<16>::UINT *)&(b+cur))[3]) >> 2;
				// gets greater exponent
				((typename StConstList<16>::UINT *)&xorc)[3] = 0x0010;


			}else expc = (((typename StConstList<16>::UINT *)& (a +cur) )[3]) >> 2;

			if (expc > expb){
				if (expc > expb + 52){
					best = cur; expb = expc; xorb = xorc;
				}else{
					xorb >>=  expc - expb;
					if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
					expb = expc;
				}
			}else{
				if (expb <= expc + 52){  // guarrantied fail otherwise
					xorc >>=  expb - expc;
					if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
				}
			}

		}

		return(best);
	}



};

template< >
class ExCo<float>{
public:
//	typedef short MAGNITUDE_TYPE;
	typedef void DERIVATIVE_TYPE;
    typedef float TYPE;

	static const bool IsPOD = true;
	static const bool hasSpetoMemmove = false;

	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
	typedef ExOp SAFETYPE;

	static bool isZero(const float &a){return a == 0.0f;}
	static bool isOne(const float &a){return a == 1.0f;}
	static bool isInfinite(const float &a){return std::isinf(a);}
	static float& toZero(float &a){return a =0.0f;}
	static float& toOne(float &a){return a =1.0f;}
	static float& toRand(float &a){a = (((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32); return a;};
	static void epsilon(float &a){a =FLT_EPSILON;}
	static void delta(float &a){a = FLT_MIN;}
	static float& toMaximum(float &a){a = FLT_MAX; return a;}
	static float& toMinimum(float &a){a = -FLT_MAX; return a;}
    static float& toUndefined(float &a){a = nan(""); return a;}
    static float& toAbs(float &a){return a = fabs(a);}
    static float mkSqrt(const float &a){return sqrt(a);}
   // static float toAbsoft(const float &a, const double& k){float fout = fabs(a); return (fout > k * 16) ? fout : a / k  }
    inline static float& toAbhard(float &a, const double& k){a = fabs(a); if (a < k) a = k; return a;}
    inline static float mkAbhard(const float &a, const double& k){float fout = fabs(a); return (fout < k) ? k : fout;}
    inline static float mkAbhardInverse(const float &a, const double& k){float fout = fabs(a); return 1.0f / ((fout < k) ? k : fout);}
    inline static float mkAbsoft(const float &a, const double& k){
        float fout;
        if ((fout = fabs(a)) > k*16) return fout;
        fout /= k;
        return log2(exp2(fout) + exp2(-fout)) * k;
        }
    inline static float mkAbsoftInverse(const float &a, const double& k){
        float fout;
        if ((fout = fabs(a)) > k*16) return 1.0f / fout;
        fout /= k;
        return 1.0f  / (log2(exp2(fout) + exp2(-fout)) * k);
        }
	PRIMTYPE_CONSTRUCT_FUNCTION(float)


	static float mkInverse(const float &a){return(1.0f / a);}
	static float mkPowInt(const float& what, const int power){ return( pow(what, (double) power)); }
	static float mkPowInvInt(const float& what, const int power){ return( pow(what, 1.0f / ((double) power))); }

	inline static float& toinverse(float &a){return (a = 1.0f / a);}
	inline static float& toNegative(float &a){return (a = -a);}
	inline static float& tosquarre(float &a){ return(a *= a); }

	SIMPLE_FUNCTIONS(float)
    SIMPLE_FUNCTIONS_FLOAT(float)


	static bool isValid(const float& a){
		return(a + (-a) == 0);
	}
	static bool isnegative(const float& a){
		return(a< 0.0f);
	}
	static float max(){return(DBL_MAX);}
	static float min(){return(DBL_MIN);}

	static float invert(const float& val){return(1.0f / val);}
	static float intPow(const float& val, int power){return(pow(val, (float) power));}
	static float floPow(const float& val, float power){return(pow(val, power));}
	static float sign(const float& val){return(val >= 0 ? 1.0f : -1.0f);}
	static float invintPow(const float& val, int power){return(pow(val, 1.0f / ((float)power)));}
	static float arctan(const float& r,const float& i){return(atan2(r,i));}
	static float floatratio(const float& num,const float& den){return(num/den);}
    static void show(const double &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%e", val); if (level == 0) *(ptr++) = '\n';}
	static void show(const float &val, FILE* out=stdout, int level=0){fprintf(out, "%e", val); if (level == 0) fprintf(out, "\n"); }
    static string type_tostring(const float& a){return string("float");}
	static ERRCODE save(const float& what, FILE *f) { return (fwrite(&what,sizeof(float), 1,f) == 1) ? 0 : 1;}
	static ERRCODE load(float & what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(float), 1,f) == 1 ? 0 : 1;}

//	static float mkMult(const float& val ,const Weight& w){ return val * w();}
//	static float& tomult(float& val ,const Weight& w){ return val *= w();}

	static double lognorm(const float& a){return log(fabs(a));}
	static double norm(const float& a){return fabs(a);}
	static double pnorm(const float& a){return a*a ;}


	template<unsigned int SIZE> static unsigned int getHyperDirection(const float (&a)[SIZE], const float (&b)[SIZE]){ // direction for hyper-ordering of arrays
		unsigned int best =0;

		typename StConstList<16>::UINT expb, expc;

		typename StConstList<64>::UINT xorb = (*((typename StConstList<64>::UINT *) &a)) ^ (*((typename StConstList<64>::UINT *) &b)); // xor and make unsigned
		// if exponent are different, the max eponent wins
		// if exponent are the same


		if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0xFFF0){
			if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0x8000) return(0);
			// most significant is max exponent

			((typename StConstList<16>::UINT *)&xorb)[3] = 0x0010; // most-significant turned on!
			expb = (((typename StConstList<16>::UINT *)&a)[3] > ((typename StConstList<16>::UINT *)&b)[3] ? ((typename StConstList<16>::UINT *)&a)[3] : ((typename StConstList<16>::UINT *)&b)[3]) >> 2;
			// gets greater exponent

		}else{
			if (xorb != 0) expb = (((typename StConstList<16>::UINT *)&a)[3]) >> 2;
			else  expb = 0;
			// gets exponent (same for both numbers)
		}




		for(unsigned int cur=1;cur<SIZE;cur++) {
			typename StConstList<64>::UINT xorc = (*((typename StConstList<64>::UINT *) &(a + cur))) ^ (*((typename StConstList<64>::UINT *) &(b + cur))); // xor and make unsigned
			if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0xFFF0){
				if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0x8000) return(cur);

				expc = (((typename StConstList<16>::UINT *)&(a+cur))[3] > ((typename StConstList<16>::UINT *)&(b+cur))[3] ? ((typename StConstList<16>::UINT *)&(a+cur))[3] : ((typename StConstList<16>::UINT *)&(b+cur))[3]) >> 2;
				// gets greater exponent
				((typename StConstList<16>::UINT *)&xorc)[3] = 0x0010;


			}else expc = (((typename StConstList<16>::UINT *)& (a +cur) )[3]) >> 2;

			if (expc > expb){
				if (expc > expb + 52){
					best = cur; expb = expc; xorb = xorc;
				}else{
					xorb >>=  expc - expb;
					if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
					expb = expc;
				}
			}else{
				if (expb <= expc + 52){  // guarrantied fail otherwise
					xorc >>=  expb - expc;
					if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
				}
			}

		}

		return(best);
	}



};

template< >
class ExCo<bool >{
public:
//	typedef unsigned char MAGNITUDE_TYPE;
	typedef void DERIVATIVE_TYPE;
    typedef bool TYPE;
	static bool const IsPOD = true;
	static bool const hasSpetoMemmove = false;
	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
	typedef ExOp SAFETYPE;


	static bool isZero(const bool &a){return !a;}
	static bool isOne(const bool &a){return a;}

	static bool& toZero(bool & a){return a = false;}
	static bool& toOne(bool & a){return a = true;}
	static bool& toRand(bool & a){ a = (rand() & 255) ; return a;}
	static void epsilon(bool & a){ a = true;}
	static void delta(bool & a){ a = true;}
	static bool& toMinimum(bool & a){ a = false; return a;}
	static bool& toMaximum(bool & a){ a = true; return a;}


	PRIMTYPE_CONSTRUCT_FUNCTION(bool)


	static double sign(const int& val){return(val >= 0 ? 1.0f : -1.0f);}

    static void show(const double &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%c", val ? 'T' : 'F'); if (level == 0) *(ptr++) = '\n';}
	static void show(const bool &val, FILE* out=stdout, int level=0){fprintf(out, "%c", val ? 'T' : 'F'); if (level == 0) fprintf(out, "\n"); }
	static string type_tostring(const bool& a){return string("bool");}

	static ERRCODE save(const bool& what, FILE *f) {return (fwrite(&what,sizeof(char), 1,f) == 1) ? 0 : 1;}
	static ERRCODE load(bool& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(char), 1,f) == 1 ? 0 : 1;}

};


	template< >
	class ExCo<char>{
	public:
//		typedef unsigned char MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;

		static bool const IsPOD = true;
    typedef char TYPE;
		static bool const hasSpetoMemmove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;


        static bool isZero(const char&a){return a == '\0';}
        static bool isOne(const char&a){return a == 1;}

		static char& toZero(char& a){return a = '\0' ;}
		static char& toOne(char& a){return a = 1 ;}
		static char& toRand(char& a){ a = (rand() & 255) ; return a;}
		static void epsilon(char& a){ a = 1;}
		static void delta(char& a){ a = 1;}
		static char& toMinimum(char& a){ a = -128; return a;}
		static char& toMaximum(char& a){ a = 127; return a;}


		PRIMTYPE_CONSTRUCT_FUNCTION(char)
        SIMPLE_FUNCTIONS(char)

		static double sign(const int& val){return(val >= 0 ? 1.0f : -1.0f);}
        static char& toAbs(char & a){return a = (a & 0x80) ? -a : a;}
		static char mkAbs(const char &a){return (a & 0x80) ? -a : a;}

	    static void show(const char &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%c", val); if (level == 0) *(ptr++) = '\n';}
		static void show(const char& val, FILE* out=stdout, int level=0){fprintf(out, "%c", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const char& a){return string("char");}

		static ERRCODE save(const char& what, FILE *f) { return (fwrite(&what,sizeof(char), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load(char& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(char), 1,f) == 1 ? 0 : 1;}
		static void bitreverse(char& a){
			a = ((a << 4 ) & 0xF0) | ((a >> 4)& 0x0F);
			a = ((a << 2 ) & 0xCC) | ((a >> 2)& 0x33);
			a = ((a << 1 ) & 0xAA) | ((a >> 1)& 0x55);
		}
		static void bytereverse(char& a){}

		static int32_t mkMagnitude(const char& a){
				int32_t f,r;
				if (a & 0xF0) {r = a >> 4;f = 5;}
				else{r = a;f = 1;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const char (&a)[SIZE], const char (&b)[SIZE]){ // direction for hyper-ordering of arrays
			uint8_t xorb,xorc;
			*(char*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(char*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

		inline static char& toRightFlood(char& a){a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static char& toLeftFlood(char& a){a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static char& toBitInverse(char& a){return a ^= 0xFF;}

	};


	template<>
	class ExCo<uint8_t>{
	public:
		typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
    typedef uint8_t TYPE;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpetoMemmove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

        static bool isZero(const uint8_t &a){return a == '\0';}
        static bool isOne(const uint8_t &a){return a == 1;}

		static uint8_t& toZero(uint8_t & a){return a = 0 ;};
		static uint8_t& toOne(uint8_t & a){return a = 1 ;};
		static uint8_t& toRand(uint8_t & a){ a =  (rand() & 255) ;return a;};
		static void epsilon(uint8_t & a){ a = 1;}
		static void delta(uint8_t & a){ a = 1;}
		static uint8_t& toMinimum(uint8_t & a){ a = 0; return a;}
		static uint8_t& toMaximum(uint8_t & a){ a = 255; return a;}
        static uint8_t& toAbs(uint8_t & a){return a;}
		PRIMTYPE_CONSTRUCT_FUNCTION(uint8_t)
        SIMPLE_FUNCTIONS(uint8_t)
        static void show(const uint8_t &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%i", val); if (level == 0) *(ptr++) = '\n';}
        static void show(const uint8_t &val, FILE* out=stdout, int level=0){fprintf(out, "%i", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const uint8_t & a){return string("uint8_t");}
		static ERRCODE save(const uint8_t& what, FILE *f) {return (fwrite(&what,sizeof(uint8_t), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load( uint8_t& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(uint8_t), 1,f) == 1 ? 0 : 1;}

		static uint8_t upperbound_pow_of_2(const uint8_t &i){
			if (i <= 16){
				if (i <= 4){
					if (i > 2) return(2);
					else if (i < 2) return(0);
					else return(1);
				}else return( (i <=8) ? 3 : 4 );
			}else return ((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
		}
		static void bitreverse(uint8_t& a){
			a = ((a << 4 ) & 0xF0) | ((a >> 4)& 0x0F);
			a = ((a << 2 ) & 0xCC) | ((a >> 2)& 0x33);
			a = ((a << 1 ) & 0xAA) | ((a >> 1)& 0x55);
		}
		static void bytereverse(uint8_t& a){}
		static int32_t mkMagnitude(const uint8_t& a){
				int32_t f,r;
				if (a & 0xF0) {r = a >> 4;f = 5;}
				else{r = a;f = 1;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}

		template<unsigned int SIZE> static unsigned int getHyperDirection(const uint8_t (&a)[SIZE], const uint8_t (&b)[SIZE]){ // direction for hyper-ordering of arrays
			uint8_t xorb = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				uint8_t xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
		inline static uint8_t& toRightFlood(uint8_t& a){a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static uint8_t& toLeftFlood(uint8_t& a){a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static uint8_t& toBitInverse(uint8_t& a){return a ^= 0xFF;}
	};


	template< >
	class ExCo<short>{
	public:
//		typedef unsigned char MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
    typedef short TYPE;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpetoMemmove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;


        PRIMTYPE_CONSTRUCT_FUNCTION(int16_t)
        SIMPLE_FUNCTIONS(int16_t)

        static bool isZero(const short &a){return a == 0;}
        static bool isOne(const short &a){return a == 1;}
		static short& toZero(short & a){return a = 0;}
		static short& toOne(short & a){return a = 1;}
		static short& toRand(short & a ){ a = rand() ^ (rand() << 12);return a;};

		static short& toMinimum(short& a) {a=-32768;return a;}
		static short& toMaximum(short& a) {a=32767;return a;}
        static void epsilon(short & a){a = 1;}
		static void delta(short & a){a = 1;}

        static double sign(const short& val){return(0.0f);}

        static int16_t& toAbs(int16_t & a){return a = (a & 0x8000) ? -a : a;}
		static int16_t mkAbs(const int16_t &a){return (a & 0x8000) ? -a : a;}
        static void show(const int16_t &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%i", val); if (level == 0) *(ptr++) = '\n';}
		static void show(const int16_t &val, FILE* out=stdout, int level=0){fprintf(out, "%i", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const short & a){return string("short");}


		static ERRCODE save(const short& what, FILE *f) { return (fwrite(&what,sizeof(short), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load(short& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(short), 1,f) == 1 ? 0 : 1;}

		static void bitreverse(short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
			a = ((a << 4 ) & 0xF0F0) | ((a >> 4)& 0x0F0F);
			a = ((a << 2 ) & 0xCCCC) | ((a >> 2)& 0x3333);
			a = ((a << 1 ) & 0xAAAA) | ((a >> 1)& 0x5555);
		}
		static void bytereverse(short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
		}


		static int32_t mkMagnitude(const short& a){
				int32_t f;
				short r;
				if (a & 0xFF00) {r = a >> 8;f = 9;}
				else{r = a;f = 1;}
				if (r & 0xF0) {r >>= 4; f+=4;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const short (&a)[SIZE], const short (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned short xorb,xorc;
			*(unsigned short*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(unsigned short*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

		inline static short& toRightFlood(short& a){a |= a << 8;a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static short& toLeftFlood(short& a){a |= a >> 8;a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static short& toBitInverse(short& a){return a ^= 0xFFFF;}
	};


	template< >
	class ExCo<unsigned short>{
	public:
		static const bool IsPOD = true;
		typedef void DERIVATIVE_TYPE;
    typedef unsigned short TYPE;
		//	typedef YESNO<true> IsPOD;

		static bool const hasSpetoMemmove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

        PRIMTYPE_CONSTRUCT_FUNCTION(unsigned short)
        SIMPLE_FUNCTIONS(uint16_t)

        static bool isZero(const unsigned short &a){return a == 0;}
        static bool isOne(const unsigned short &a){return a == 1;}

		static unsigned short& toZero(unsigned short & a){return a = 0;}
		static unsigned short& toOne(unsigned short & a){return a = 1;}
		static unsigned short& toRand(unsigned short & a ){ a = rand() ^ (rand() << 12);return a;}
        static void epsilon(unsigned short & a){a = 1;}
		static void delta(unsigned short & a){a = 1;}
		static unsigned short& toMinimum(unsigned short& a) {a=0;return a;}
		static unsigned short& toMaximum(unsigned short& a) {a=65535;return a;}
        static unsigned short& toAbs(unsigned short & a){return a;}

        static void show(const uint16_t &val, char* &ptr, int level=0){ptr += sprintf(ptr, "%i", val); if (level == 0) *(ptr++) = '\n';}
		static void show(const uint16_t &val, FILE* out=stdout, int level=0){fprintf(out, "%u", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const unsigned short & a){return string("unsigned short");}

		static ERRCODE save(const unsigned short& what, FILE *f) {return (fwrite(&what,sizeof(unsigned short), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load( unsigned short& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(unsigned short), 1,f) == 1 ? 0 : 1;}

		static unsigned char upperbound_pow_of_2(const unsigned short &i){
			if (i <= 256){
				if (i <= 16){
					if (i <= 4){
						if (i > 2) return(2);
						else if (i < 2) return(0);
						else return(0);
					}else return( (i <=8) ? 3 : 4 );
				}else return((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
			}else return( (i<= 0x1000) ? ((i <= 0x0400) ? ((i <=0x200) ? 9 : 10) : ((i <=0x0800) ? 11 : 12 )) : ((i <= 0x4000) ? ((i <=0x2000) ? 13 : 14) : ((i <=0x8000) ? 15 : 16 )) );
		}
		static void bitreverse(unsigned short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
			a = ((a << 4 ) & 0xF0F0) | ((a >> 4)& 0x0F0F);
			a = ((a << 2 ) & 0xCCCC) | ((a >> 2)& 0x3333);
			a = ((a << 1 ) & 0xAAAA) | ((a >> 1)& 0x5555);
		}
		static void bytereverse(unsigned short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
		}
		static int32_t mkMagnitude(const unsigned short& a){
				int32_t f;
				unsigned short r;
				if (a & 0xFF00) {r = a >> 8;f = 9;}
				else{r = a;f = 1;}
				if (r & 0xF0) {r >>= 4; f+=4;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const unsigned short (&a)[SIZE], const unsigned short (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned short xorb = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				unsigned short xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
		inline static unsigned short& toRightFlood(unsigned short& a){a |= a << 8;a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static unsigned short& toLeftFlood(unsigned short& a){a |= a >> 8;a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static unsigned short& toBitInverse(unsigned short& a){return a ^= 0xFFFF;}
	};

	template< >
	class ExCo<int32_t>{
	public:
		typedef void DERIVATIVE_TYPE;
        typedef int32_t TYPE;
		static const bool IsPOD = true;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		//	typedef YESNO<true> IsPOD;
		typedef ExOp SAFETYPE;
		static bool const hasSpetoMemmove = false;


		PRIMTYPE_CONSTRUCT_FUNCTION(int32_t)
        SIMPLE_FUNCTIONS(int32_t)
        static bool isZero(const int32_t &a){return a == 0;}
        static bool isOne(const int32_t &a){return a == 1;}

		static int32_t& toZero(int32_t & a){return a = 0;}
		static int32_t& toOne(int32_t & a){return a = 1;}
		static int32_t& toRand(int32_t & a){ a = rand() ^ (rand() << 12) ^ (rand() << 24); return a;};

		static double sign(const int32_t& val){return(val >= 0 ? 1.0f : -1.0f);}
		static void show(const int32_t &val, FILE* out=stdout, int32_t level=0){fprintf(out, "%i", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const int32_t & a){return string("int");}

		static int32_t& toMinimum(int32_t& a) {a=0x80000000;return a;}
		static int32_t& toMaximum(int32_t& a) {a=0x7FFFFFFF;return a;}
		static int32_t& toAbs(int32_t & a){return a = (a & 0x80000000) ? -a : a;}
		static int32_t mkAbs(const int32_t &a){return (a & 0x80000000) ? -a : a;}
		static void epsilon(int32_t& a) {a=1;}
		static void delta(int32_t& a) {a=1;}

		static ERRCODE save(const int32_t& what, FILE *f) {return (fwrite(&what,sizeof(int32_t), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load(int32_t& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(int32_t), 1,f) == 1 ? 0 : 1;}

		static void bitreverse(int32_t& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
			a = ((a << 4 ) & 0xF0F0F0F0) | ((a >> 4)& 0x0F0F0F0F);
			a = ((a << 2 ) & 0xCCCCCCCC) | ((a >> 2)& 0x33333333);
			a = ((a << 1 ) & 0xAAAAAAAA) | ((a >> 1)& 0x55555555);
		}
		static void bytereverse(int32_t& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
		}
		static int32_t mkMagnitude(const int32_t& a){
				char f;
				int32_t r;
				if (a & 0xFFFF0000) {r = a >> 16;f = 17;}
				else{r = a;f = 1;}
				if (r & 0xFF00) {r >>= 8; f+=8;}
				if (r & 0xF0) {r >>= 4; f+=4;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const int32_t (&a)[SIZE], const int32_t (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned int xorb,xorc;
			*(unsigned int*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(unsigned int*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
		inline static int32_t& toRightFlood(int32_t& a){a |= a << 16;a |= a << 8;a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static int32_t& toLeftFlood(int32_t& a){a |= a >> 16;a |= a >> 8;a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static int32_t& toBitInverse(int32_t& a){return a ^= 0xFFFFFFFF;}
	};


	template< >
	class ExCo<uint32_t >{
	public:
typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
            typedef uint32_t TYPE;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpetoMemmove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

		PRIMTYPE_CONSTRUCT_FUNCTION(uint32_t)
        SIMPLE_FUNCTIONS(uint32_t)

        static bool isZero(const uint32_t &a){return a == 0;}
        static bool isOne(const uint32_t &a){return a == 1;}
        static uint32_t & toMemswap(uint32_t &a, uint32_t &b){uint32_t tmp = b; b=a;return a=tmp;}

		static uint32_t& toZero(uint32_t& a){return a = 0;}
		static uint32_t& toOne(uint32_t& a){return a = 1;}
		static uint32_t& toRand(uint32_t & a ){ a = rand() ^ (rand() << 12) ^ (rand() << 24);return a;}

		static uint32_t& toMinimum(uint32_t& a) {a=0;return a;}
		static uint32_t& toMaximum(uint32_t& a) {a=0xFFFFFFFF;return a;}
		static uint32_t& toAbs(uint32_t & a){return a;}
		static void epsilon(uint32_t& a) {a=1;}
		static void delta(uint32_t& a) {a=1;}

		static void show(const uint32_t&val, FILE* out=stdout, int level=0){fprintf(out, "%u", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const uint32_t & a){return string("uint32_t");}

		static ERRCODE save(const uint32_t& what, FILE *f) {return (fwrite(&what,sizeof(uint32_t), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load( uint32_t& what, FILE *f, uint32_t lenght=0) {return fread(&what,sizeof(uint32_t), 1,f) == 1 ? 0 : 1;}


		static unsigned char upperbound_pow_of_2(const uint32_t &i){
			if (i <= 65536){
				if (i <= 256){
					if (i <= 16){
						if (i <= 4){
							if (i > 2) return(2);
							else if (i < 2) return(0);
							else return(0);
						}else return( (i <=8) ? 3 : 4 );
					}else return((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
				}else return( (i<= 0x1000) ? ((i <= 0x0400) ? ((i <=0x200) ? 9 : 10) : ((i <=0x0800) ? 11 : 12 )) : ((i <= 0x4000) ? ((i <=0x2000) ? 13 : 14) : ((i <=0x8000) ? 15 : 16 )) );
			}else return( (i<= 0x01000000) ? ((i<= 0x00100000) ? ((i <= 0x00040000) ? ((i <=0x00020000) ? 17 : 18) : ((i <=0x00080000) ? 19 : 20 )) : ((i <= 0x00400000) ? ((i <=0x00200000) ? 21 : 22) : ((i <=0x00800000) ? 23 : 24 )) ) : ((i<= 0x10000000) ? ((i <= 0x04000000) ? ((i <=0x02000000) ? 25 : 26) : ((i <=0x08000000) ? 27 : 28 )) : ((i <= 0x40000000) ? ((i <=0x20000000) ? 29 : 30) : ((i <=0x80000000) ? 31 : 32 )) ));
		}


		static void bitreverse(uint32_t& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
			a = ((a << 4 ) & 0xF0F0F0F0) | ((a >> 4)& 0x0F0F0F0F);
			a = ((a << 2 ) & 0xCCCCCCCC) | ((a >> 2)& 0x33333333);
			a = ((a << 1 ) & 0xAAAAAAAA) | ((a >> 1)& 0x55555555);
		}
		static void bytereverse(uint32_t& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
		}

		static int32_t mkMagnitude(const uint32_t& a){
				int32_t f;
				uint32_t r;
				if (a & 0xFFFF0000) {r = a >> 16;f = 17;}
				else{r = a;f = 1;}
				if (r & 0xFF00) {r >>= 8; f+=8;}
				if (r & 0xF0) {r >>= 4; f+=4;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}

		template<unsigned int SIZE> static uint32_t getHyperDirection(const uint32_t (&a)[SIZE], const uint32_t (&b)[SIZE]){ // direction for hyper-ordering of arrays
			uint32_t xorb = a[0] ^ b[0]; // xor and make unsigned
			uint32_t best =0;
			for(uint32_t cur=1;cur<SIZE;cur++) {
				uint32_t xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
        inline static uint32_t readBits(uint32_t* data, uint32_t offset, uint32_t nbbit){
            uint32_t bstart = offset * nbbit;
            data += (bstart>>5); bstart &= 31;
            uint32_t bend = (bstart + nbbit-1) & 31;
            if (bend >= bstart) return (data[0] >> (bstart & 31)) & (0xFFFFFFFF >> (32 - nbbit));
            return (data[0] >> bstart) | ((data[1] & (0xFFFFFFFF >> (31 - bend))) << (32 - bstart));
        }
        inline static void writeBits(uint32_t* data, uint32_t offset, uint32_t nbbit, uint32_t value){
            uint32_t bstart = offset * nbbit;
            data += (bstart>>5); bstart &= 31;
            uint32_t bend = (bstart + nbbit-1) & 31;
            if (bend >= bstart) data[0] = (data[0] & ((0xFFFFFFFE << (bend & 31)) |(0xFFFFFFFF >> 32- bstart))) | (value << bstart);
            else{
                data[0] = (data[0] & (0xFFFFFFFF >> (32- bstart))) | (value << bstart);
                data[1] = (data[1] & (0xFFFFFFFE << bend)) | (value >> (32 - bstart));
            }
        }
		inline static uint32_t& toRightFlood(uint32_t& a){a |= a << 16;a |= a << 8;a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static uint32_t& toLeftFlood(uint32_t& a){a |= a >> 16;a |= a >> 8;a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static uint32_t& toBitInverse(uint32_t& a){return a ^= 0xFFFFFFFF;}
	};
	template< >
	class ExCo<int64_t>{
	public:
		typedef void DERIVATIVE_TYPE;
        typedef int64_t TYPE;
		static const bool IsPOD = true;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		//	typedef YESNO<true> IsPOD;
		typedef ExOp SAFETYPE;
		static bool const hasSpetoMemmove = false;


		PRIMTYPE_CONSTRUCT_FUNCTION(int64_t)
        SIMPLE_FUNCTIONS(int64_t)
        static bool isZero(const int64_t &a){return a == 0;}
        static bool isOne(const int64_t &a){return a == 1;}

		static int64_t& toZero(int64_t & a){return a = 0;}
		static int64_t& toOne(int64_t & a){return a = 1;}
		static int64_t& toRand(int64_t & a ){ a = (((int64_t)(( ((int32_t)rand()) << 12) ^ rand())) << 36) ^ (((int64_t)(rand() << 24)) ^ (rand() ^ (rand() << 12)));return a;};

		static double sign(const int64_t& val){return(val >= 0 ? 1.0f : -1.0f);}

    #ifdef __MINGW32__
		static void show(const int64_t &val, FILE* out=stdout, int64_t level=0){fprintf(out, "%" PRId64 "i", val); if (level == 0) fprintf(out, "\n"); }
    #else
		static void show(const int64_t &val, FILE* out=stdout, int64_t level=0){fprintf(out, "%li", val); if (level == 0) fprintf(out, "\n"); }
    #endif

        static string type_tostring(const int64_t & a){return string("int");}

		static int64_t& toMinimum(int64_t& a) {a=0x8000000000000000;return a;}
		static int64_t& toMaximum(int64_t& a) {a=0x7FFFFFFFFFFFFFFF;return a;}
		static int64_t& toAbs(int64_t & a){return a = (a & 0x8000000000000000) ? -a : a;}
		static int64_t mkAbs(const int64_t &a){return (a & 0x8000000000000000) ? -a : a;}
		static void epsilon(int64_t& a) {a=1;}
		static void delta(int64_t& a) {a=1;}


		static ERRCODE save(const int64_t& what, FILE *f) {return (fwrite(&what,sizeof(int64_t), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load(int64_t& what, FILE *f, unsigned int lenght=0) {return fread(&what,sizeof(int64_t), 1,f) == 1 ? 0 : 1;}

		static void bitreverse(int64_t& a){
			a = ((a << 32 ) & 0xFFFFFFFF00000000) | ((a >> 32)& 0x00000000FFFFFFFF);
			a = ((a << 16 ) & 0xFFFF0000FFFF0000) | ((a >> 16)& 0x0000FFFF0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00FF00FF00) | ((a >> 8)& 0x00FF00FF00FF00FF);
			a = ((a << 4 ) & 0xF0F0F0F0F0F0F0F0) | ((a >> 4)& 0x0F0F0F0F0F0F0F0F);
			a = ((a << 2 ) & 0xCCCCCCCCCCCCCCCC) | ((a >> 2)& 0x3333333333333333);
			a = ((a << 1 ) & 0xAAAAAAAAAAAAAAAA) | ((a >> 1)& 0x5555555555555555);
		}
		static void bytereverse(int64_t& a){
			a = ((a << 32 ) & 0xFFFFFFFF00000000) | ((a >> 32)& 0x00000000FFFFFFFF);
			a = ((a << 16 ) & 0xFFFF0000FFFF0000) | ((a >> 16)& 0x0000FFFF0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00FF00FF00) | ((a >> 8)& 0x00FF00FF00FF00FF);
		}
		static int32_t mkMagnitude(const int64_t& a){
				int32_t f;
				int64_t r;
                if (a & 0xFFFFFFFF00000000) {r = a >> 32;f = 33;}
				else{r = a;f = 1;}
				if (r & 0xFFFF0000) {r >>= 16; f+=16;}
				if (r & 0xFF00) {r >>= 8; f+=8;}
				if (r & 0xF0) {r >>= 4; f+=4;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const int64_t (&a)[SIZE], const int64_t (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned int xorb,xorc;
			*(unsigned int*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(unsigned int*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
		inline static int64_t& toRightFlood(int64_t& a){a |= a << 32;a |= a << 16;a |= a << 8;a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static int64_t& toLeftFlood(int64_t& a){a |= a >> 32;a |= a >> 16;a |= a >> 8;a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static int64_t& toBitInverse(int64_t& a){return a ^= 0xFFFFFFFFFFFFFFFF;}
	};


	template< >
	class ExCo<uint64_t >{
	public:
typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
            typedef uint64_t TYPE;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpetoMemmove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

		PRIMTYPE_CONSTRUCT_FUNCTION(uint64_t)
        SIMPLE_FUNCTIONS(uint32_t)

        static bool isZero(const uint64_t &a){return a == 0;}
        static bool isOne(const uint64_t &a){return a == 1;}
        static uint64_t & toMemswap(uint64_t &a, uint64_t &b){uint64_t tmp = b; b=a;return a=tmp;}

		static uint64_t& toZero(uint64_t & a){return a = 0;}
		static uint64_t& toOne(uint64_t & a){return a = 1;}
		static uint64_t& toRand(uint64_t & a ){ a = (((uint64_t)(( ((uint32_t)rand()) << 12) ^ rand())) << 36) ^ (((uint64_t)(rand() << 24)) ^ (rand() ^ (rand() << 12)));return a;}

		static uint64_t& toMinimum(uint64_t& a) {a=0;return a;}
		static uint64_t& toMaximum(uint64_t& a) {a=0xFFFFFFFFFFFFFFFF;return a;}
		static uint64_t& toAbs(uint64_t & a){return a;}
		static void epsilon(uint64_t& a) {a=1;}
		static void delta(uint64_t& a) {a=1;}

    #ifdef __MINGW32__
		static void show(const uint64_t&val, FILE* out=stdout, int level=0){fprintf(out, "%" PRIu64 "", val); if (level == 0) fprintf(out, "\n"); }
    #else
		static void show(const int64_t &val, FILE* out=stdout, int64_t level=0){fprintf(out, "%lu", val); if (level == 0) fprintf(out, "\n"); }
    #endif
        static string type_tostring(const uint64_t & a){return string("uint64_t");}

		static ERRCODE save(const uint64_t& what, FILE *f) {return (fwrite(&what,sizeof(uint64_t), 1,f) == 1) ? 0 : 1;}
		static ERRCODE load( uint64_t& what, FILE *f, uint64_t lenght=0) {return fread(&what,sizeof(uint64_t), 1,f) == 1 ? 0 : 1;}


		static unsigned char upperbound_pow_of_2(const uint64_t &i){
			if (i <= 65536){
				if (i <= 256){
					if (i <= 16){
						if (i <= 4){
							if (i > 2) return(2);
							else if (i < 2) return(0);
							else return(0);
						}else return( (i <=8) ? 3 : 4 );
					}else return((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
				}else return( (i<= 0x1000) ? ((i <= 0x0400) ? ((i <=0x200) ? 9 : 10) : ((i <=0x0800) ? 11 : 12 )) : ((i <= 0x4000) ? ((i <=0x2000) ? 13 : 14) : ((i <=0x8000) ? 15 : 16 )) );
			}else return( (i<= 0x01000000) ? ((i<= 0x00100000) ? ((i <= 0x00040000) ? ((i <=0x00020000) ? 17 : 18) : ((i <=0x00080000) ? 19 : 20 )) : ((i <= 0x00400000) ? ((i <=0x00200000) ? 21 : 22) : ((i <=0x00800000) ? 23 : 24 )) ) : ((i<= 0x10000000) ? ((i <= 0x04000000) ? ((i <=0x02000000) ? 25 : 26) : ((i <=0x08000000) ? 27 : 28 )) : ((i <= 0x40000000) ? ((i <=0x20000000) ? 29 : 30) : ((i <=0x80000000) ? 31 : 32 )) ));
		}


		static void bitreverse(uint64_t& a){
			a = ((a << 32 ) & 0xFFFFFFFF00000000) | ((a >> 32)& 0x00000000FFFFFFFF);
			a = ((a << 16 ) & 0xFFFF0000FFFF0000) | ((a >> 16)& 0x0000FFFF0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00FF00FF00) | ((a >> 8)& 0x00FF00FF00FF00FF);
			a = ((a << 4 ) & 0xF0F0F0F0F0F0F0F0) | ((a >> 4)& 0x0F0F0F0F0F0F0F0F);
			a = ((a << 2 ) & 0xCCCCCCCCCCCCCCCC) | ((a >> 2)& 0x3333333333333333);
			a = ((a << 1 ) & 0xAAAAAAAAAAAAAAAA) | ((a >> 1)& 0x5555555555555555);
		}
		static void bytereverse(int64_t& a){
			a = ((a << 32 ) & 0xFFFFFFFF00000000) | ((a >> 32)& 0x00000000FFFFFFFF);
			a = ((a << 16 ) & 0xFFFF0000FFFF0000) | ((a >> 16)& 0x0000FFFF0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00FF00FF00) | ((a >> 8)& 0x00FF00FF00FF00FF);
		}

		static int32_t mkMagnitude(const uint64_t& a){
				int32_t f;
				uint64_t r;
				if (a & 0xFFFFFFFF00000000) {r = a >> 32;f = 33;}
				else{r = a;f = 1;}
				if (r & 0xFFFF0000) {r >>= 16; f+=16;}
				if (r & 0xFF00) {r >>= 8; f+=8;}
				if (r & 0xF0) {r >>= 4; f+=4;}
				if (r & 0xC) {r >>= 2; f+=2;}
				if (r & 0x2) f++;
				else if (r == 0) f--;
				return f;
			}

		template<uint64_t SIZE> static uint64_t getHyperDirection(const uint64_t (&a)[SIZE], const uint64_t (&b)[SIZE]){ // direction for hyper-ordering of arrays
			uint64_t xorb = a[0] ^ b[0]; // xor and make unsigned
			uint64_t best =0;
			for(uint64_t cur=1;cur<SIZE;cur++) {
				uint64_t xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

		inline static uint64_t& toRightFlood(uint64_t& a){a |= a << 32;a |= a << 16;a |= a << 8;a |= a << 4; a |= a << 2; return a |= a << 1;}
		inline static uint64_t& toLeftFlood(uint64_t& a){a |= a >> 32;a |= a >> 16;a |= a >> 8;a |= a >> 4; a |= a >> 2; return a |= a >> 1;}
		inline static uint64_t& toBitInverse(uint64_t& a){return a ^= 0xFFFFFFFFFFFFFFFF;}
	};

	template< >
	class ExCo<string >{
	public:
        typedef string TYPE;

		static const bool IsPOD = false;
		static const bool NeedsAddLink = false;
		typedef ExCo<string> SAFETYPE;
		typedef string INDEX_TYPE;
		typedef string REAL_TYPE;
		typedef string COMPLEX_TYPE;
		typedef string DEF_TYPE;
		typedef string NEG_TYPE;
		typedef string TRJU_TYPE;
        typedef string LMUL_TYPE;
		typedef string INNER_TYPE;
        typedef string VECTOR_TYPE;
		typedef string OUTER_TYPE; // outer product result, (increase in the number of dimentions!)

		typedef string QUATERNION_TYPE;
		typedef string GAUS_TYPE;

        typedef YESNO<false> IS_COMMUTATIVE;
        typedef YESNO<false> IS_OWNED;
        typedef YESNO<false> IS_POD;

		static string& toMemmove(string & tar, string & scr){tar.swap(scr); scr.clear(); return tar;}
		static string& toMemfree(string & tar){tar.clear(); return tar;}
		static string mkZero(){return string("");}
		static string mkOne(){return string("1");}
		static void toZero(string& a){a = string("");}
		static void toOne(string& a){a = string("1");}
		static string& toRand(string& a){ a = string("Random"); return a;}
		static void show(const string& what, FILE *f,int level) {
            if (what.length() == 0){
                if (level == 0) fprintf(f,"empty string\n"); else fprintf(f,"empty string");
            }else{
                if (level == 0) fprintf(f,"%s\n",what.c_str()); else fprintf(f,"%s",what.c_str());
            }
        }
        static ERRCODE save(const string& what, FILE *f) {unsigned int s = what.length(); if (s == 0) return (fwrite(&s,sizeof(char), 1,f) == 1) ? 0 : 1; else return (fwrite(what.c_str(),sizeof(char), s+1,f) == s+1) ? 0 : 1;}
		static ERRCODE load( string& what, FILE *f, unsigned int lenght=0) { char buffer[65536]; char* tmp = buffer; do {fread(tmp,sizeof(unsigned char), 1,f); } while(*(tmp++) != '\0'); what = string(buffer);  return 0;}
	};

	template< >
	class ExCo<SETCMP_enum>{
        public:
		static void show(const SETCMP_enum& what, FILE *f,int level) {
            switch(what){
                case SETCMP_EQUAL: fprintf(f,"Equal");
                case SETCMP_DISJOINT: fprintf(f,"Disjoint");
                case SETCMP_GE: fprintf(f,"GorEqual");
                case SETCMP_LE: fprintf(f,"LorEqual");
                case SETCMP_GT: fprintf(f,"GreaterThan");
                case SETCMP_LT: fprintf(f,"SmallerThan");
                case SETCMP_IS_NOT_SUBSET: fprintf(f,"Superset");
                case SETCMP_IS_NOT_SUPERSET: fprintf(f,"Subset");
                case SETCMP_SUP_BOUNDED: fprintf(f,"ExtremaInside");
                case SETCMP_SUB_BOUNDED: fprintf(f,"ExtremaOutside");
                case SETCMP_MAX_GT: fprintf(f,"GreaterMaximum");
                case SETCMP_MAX_LT: fprintf(f,"SmallerMaximum");
                case SETCMP_MIN_GT: fprintf(f,"GreaterMinimum");
                case SETCMP_MIN_LT: fprintf(f,"SmallerMinimum");
                default: fprintf(f,"SetCMP0x%X", (unsigned int)what);
            }
            fprintf(f,"%c", (level == 0) ? '\n' : '\t');
            }

    };


	template<class A>
	class ExCo<std::function<void(A)>, 0u >{
        public:
        typedef ExCo<std::function<void(A)>, 0u > SAFETYPE;
        static const bool IsPOD = false;

        static std::function<void(A)>& toZero(std::function<void(A)>& what){return what;}

        static std::function<void(A)>& toClone(std::function<void(A)>& what, const std::function<void(A)>& from){return what = from;}
        static std::function<void(A)>& toMemmove(std::function<void(A)>& what, std::function<void(A)>& from){return what = from;}
		static void show(const std::function<void(A)>& what, FILE *f,int level) {
                fprintf(f,"functor\n");
            }

    };


// ExOp functions
template<class A, class B> inline bool ExOp::isLT(const A &a, const B &b){ return( ExFn< Exlisten_lt<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::isLT(a,b));}/**< \brief True if A is smaller than B \param Object A \param Object B*/
template<class A, class B> inline bool ExOp::isGT(const A &a, const B &b){ return( ExFn< Exlisten_gt<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::isGT(a,b));}/**< \brief True if A is larger than B \param Object A \param Object B*/
template<class A, class B> inline bool ExOp::isLE(const A &a, const B &b){ return( ExFn< Exlisten_le<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::isLE(a,b));}/**< \brief True if A is smaller or equal to B \param Object A \param Object B*/
template<class A, class B> inline bool ExOp::isGE(const A &a, const B &b){ return( ExFn< Exlisten_ge<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::isGE(a,b));}/**< \brief True if A is larger or equal to B \param Object A \param Object B*/
template<class A, class B> inline bool ExOp::isEQ(const A &a, const B &b){ return( ExFn< Exlisten_eq<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::isEQ(a,b));}/**< \brief True if A is equal to B \param Object A \param Object B*/
template<class A, class B> inline bool ExOp::isNQ(const A &a, const B &b){ return( ExFn< Exlisten_nq<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::isNQ(a,b));}/**< \brief True if A is not equal to B \param Object A \param Object B*/

#undef LFHTEMP
#define LFHTEMP template<class A>

LFHTEMP bool ExCo<vector<A>,0 >::isValid(const vector<A> &a){unsigned int i; for(i=0;i< a.size();i++) if (! ExOp::isValid(a[i])) break; return i == a.size();}
LFHTEMP bool ExCo<vector<A>,0 >::isZero(const vector<A> &a){unsigned int i; for(i=0;i< a.size();i++) if (! ExOp::isZero(a[i])) break; return i == a.size();}
LFHTEMP bool ExCo<vector<A>,0 >::isOne(const vector<A> &a){unsigned int i; for(i=0;i< a.size();i++) if (! ExOp::isOne(a[i])) break; return i == a.size();}
LFHTEMP vector<A>& ExCo<vector<A>,0 >::toMemfree(vector<A> &a){a.clear(); return(a);}




#undef LFHTEMP
#define LFHTEMP template<class A, unsigned int SIZE>

LFHTEMP A (&ExCo<A[SIZE], 0>::toZero( A (&a)[SIZE]))[SIZE]{	for(unsigned int i=0;i<SIZE;i++) ExOp::toZero(a[i]); return a;}
LFHTEMP A (&ExCo<A[SIZE], 0>::toRand( A (&a)[SIZE]))[SIZE]{	for(unsigned int i=0;i<SIZE;i++) ExOp::toRand(a[i]); return a;}
LFHTEMP A (&ExCo<A[SIZE], 0>::toOne( A (&a)[SIZE]))[SIZE]{	for(unsigned int i=0;i<SIZE;i++) ExOp::toOne(a[i]); return a;}
LFHTEMP void ExCo<A[SIZE], 0>::minimum( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toMin(a[i]);}
LFHTEMP void ExCo<A[SIZE], 0>::maximum( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toMax(a[i]);}

// A should not be POD
LFHTEMP A (&ExCo<A[SIZE], 0>::toMemfree(A (&a)[SIZE]))[SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toMemfree(a[i]);
    return a;
    }
// A should not be POD
LFHTEMP A (&ExCo<A[SIZE],0>::toMemmove(A (&a)[SIZE], A (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toMemmove(a[i], b[i]);
    return a;
}
LFHTEMP A (&ExCo<A[SIZE],0>::toClone(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toClone(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toAdd(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toAdd(a[i], b[i]);
    return a;
}
LFHTEMP A (&ExCo<A[SIZE],0>::toSubt(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toSubt(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toMult(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toMult(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toDivi(A (&a)[SIZE], const A (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toDivi(a[i], b[i]);
return a;}


LFHTEMP double ExCo<A[SIZE],0>::norm(const A (&a)[SIZE]){
    double norm =0;
    for(unsigned int i=0;i<SIZE;i++) norm += ExOp::pnorm(a[i]);
    return sqrt(norm);
}
LFHTEMP double ExCo<A[SIZE],0>::pnorm(const A (&a)[SIZE]){
    double norm =0;
    for(unsigned int i=0;i<SIZE;i++) norm += ExOp::pnorm(a[i]);
    return norm;
}

LFHTEMP void ExCo<A[SIZE],0>::show(const A (&a)[SIZE], FILE* f_out, int l){
	switch(l){
		case 0:
			for(unsigned int i=0;i<SIZE;i++) {
				ExOp::show(a[i],f_out,1);
				fprintf(f_out,"\n");
			}
			fprintf(f_out,"\n");
			break;
		case 1:
			for(unsigned int i=0;i<SIZE;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(a[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
			break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<SIZE;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(a[i],f_out,3);
			}
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<SIZE;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(a[i],f_out,4);
			}
			fprintf(f_out,")");
			break;
	}
}

LFHTEMP ERRCODE ExCo<A[SIZE],0>::save(const A (&a)[SIZE], FILE *f){ ERRCODE fout;
    fout = ExOp::save(a[0],f);
	for(unsigned int i=1;i<SIZE;i++) fout |= ExOp::save(a[i],f);
	return fout;
}

LFHTEMP ERRCODE ExCo<A[SIZE],0>::load(A (&a)[SIZE], FILE *f, unsigned int lenght){ ERRCODE fout=0;
	for(unsigned int i=0;i<SIZE;i++) fout |= ExOp::load(a[i], f,lenght/SIZE);
	return fout;
}


#undef LFHTEMP
#define LFHTEMP template<class A, unsigned int SIZE> template<class B>

LFHTEMP A (&ExCo<A[SIZE],0>::toClone(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toClone(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toConvert(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toConvert(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toAdd(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toAdd(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toSubt(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toSubt(a[i], b[i]);
return a;}
LFHTEMP A (&ExCo<A[SIZE],0>::toMult(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toMult(a[i], b[i]);
    return a;
}
LFHTEMP A (&ExCo<A[SIZE],0>::toDivi(A (&a)[SIZE], const B (&b)[SIZE])) [SIZE]{
    for(unsigned int i=0;i<SIZE;i++) ExOp::toDivi(a[i], b[i]);
    return a;
}


#undef LFHTEMP
#define LFHTEMP template<class A>

// const mappings
/*
LFHTEMP typename ExCo<A*>::INDEX_TYPE ExCo<A*  >::getIndex(A * & a){
    return (isTypeEqual<typename ExCo<A*>::INDEX_TYPE, A*>::ans) ? a : ExOp::getIndex(*a);
}*/
LFHTEMP void ExCo<A*,0>::show(A * const & what, FILE *f,int level) {ExOp::show(*what,f, level);}
LFHTEMP ERRCODE ExCo<A*,0>::save(A * & what, FILE *f) {return ExOp::save(*what,f);}
LFHTEMP ERRCODE ExCo<A*,0>::load(A * const & what, FILE *f, unsigned int lenght) {return ExOp::load(*what,f, lenght);}



//LFHTEMP void ExCo<const A*, 0u>::show(const A * & what, FILE *f,int level) {ExOp::show(*what,f, level);}

#undef LFHTEMP
#define LFHTEMP template<class A, int flag>

// assignments


LFHTEMP double ExCo<A,flag >::pdist(const A& a, const A& b){return ExOp::pnorm(a-b);}
LFHTEMP void ExCo<A,flag>::tointpow(A& what, const int pow) {
		int modp;
		if (pow >1){
		 modp = pow;
		A cp = what;
		while(modp > 1){
			ExOp::toSquare(what);
			if (modp & 1) ExOp::toMult(what, cp);
			modp = modp >> 1;
			}
		}else if (pow < 0){
			ExOp::toInverse(what);
			if (pow < -1){
			 modp = -pow;
			 A cp = what;
			 while(modp > 1){
			ExOp::toSquare(what);
			if (modp & 1) ExOp::toMult(what, cp);
			modp = modp >> 1;
			}
			}
			}else if (pow == 0) ExOp::toOne(what);
		}
LFHTEMP A ExCo<A,flag>::mkPowInt(const A& what, const int pow){
    int modp;
    if (pow >1){
        modp = pow;
        A cp = what;
        A fout = what;
        while(modp > 1){
            ExOp::toSquare(fout);
            if (modp & 1) ExOp::toMult(fout, cp);
            modp = modp >> 1;
        }
        return(fout);
    }else if (pow < 0){
        A fout = mkInverse(what);

        if (pow < -1){
         modp = -pow;
         A cp = fout;
         while(modp > 1){
        ExOp::toSquare(fout);
        if (modp & 1) ExOp::toMult(fout, cp);
        modp = modp >> 1;
        }
        }
        return(fout);
    }else if (pow == 0) return ExOp::mkone<A>();
    else return what;
}

LFHTEMP A ExCo<A,flag>::mkTrjuProd(const A& what){
    A fout = ExOp::mkTrju(fout);
    return fout *= what;
}
LFHTEMP A ExCo<A,flag>::mkTrju(const A& what){return what;}
LFHTEMP A& ExCo<A,flag>::toTrju(A& what){return what = ExOp::mkTrju(what);}


// builds

LFHTEMP A ExCo<A,flag>::mkPowInvInt(const A& what, const int pow){return ExOp::mkPowInt(ExOp::mkInverse(what), pow);}
LFHTEMP A ExCo<A,flag>::mknegative(const A& a){return ExOp::mkSubt(ExOp::mkzero<A>(), a);}
LFHTEMP A ExCo<A,flag>::mkInverse(const A& a){return ExOp::mkDivi(ExOp::mkone<A>(), a);}
LFHTEMP A ExCo<A,flag>::mkSquare(const A& a){return ExOp::mkMult(a,a);}
LFHTEMP A ExCo<A,flag>::mkAbs(const A& a){return ExOp::toAbs(A(a));}


	// assignments

LFHTEMP bool ExCo<A,flag>::isGT(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
LFHTEMP bool ExCo<A,flag>::isGE(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_E_MASK) != SETCMP_LE);}
LFHTEMP bool ExCo<A,flag>::isLT(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_T_MASK) == SETCMP_LT);}
LFHTEMP bool ExCo<A,flag>::isLE(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_E_MASK) != SETCMP_GE);}
LFHTEMP bool ExCo<A,flag>::isEQ(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
LFHTEMP bool ExCo<A,flag>::isNQ(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}



	// const mappings
#undef LFHTEMP
#define LFHTEMP template<class A>

LFHTEMP void ExCo<const A,0>::show(const A &a, FILE* fout, int level){ExCo<A>::show(a,fout,level);}
LFHTEMP bool ExCo<const A,0>::isValid(const A& a){return(ExCo<A>::isValid(a));}



#undef LFHTEMP
#define LFHTEMP template<class A, class B>


LFHTEMP bool ExCo<pair<A,B>,0>::isValid(const pair<A,B> &a){ return ExOp::isValid(a.first) &&  ExOp::isValid(a.second) ;}

LFHTEMP bool ExCo<pair<A,B>,0>::isZero(const pair<A,B> &a){return ExOp::isZero(a.first) && ExOp::isZero(a.second) ;}
LFHTEMP bool ExCo<pair<A,B>,0>::isOne(const pair<A,B> &a){return ExOp::isOne(a.first) && ExOp::isOne(a.second) ;}

LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toZero(pair<A,B> &a){ExOp::toZero(a.first); ExOp::toZero(a.second); return a;}
LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toOne(pair<A,B> &a){ExOp::toOne(a.first); ExOp::toOne(a.second); return a;}
LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toRand(pair<A,B> &a){ExOp::toRand(a.first); ExOp::toRand(a.second); return a;}

LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toClone(pair<A,B> &a, const pair<A,B> &b){ExOp::toClone(a.first,b.first); ExOp::toClone(a.second,b.second); return a;}
LFHTEMP template<class A2, class B2> pair<A,B>& ExCo<pair<A,B>,0>::toClone(pair<A,B> &a, const pair<A2,B2> &b){ExOp::toClone(a.first,b.first); ExOp::toClone(a.second,b.second); return a;}
LFHTEMP template<class A2, class B2> pair<A,B>& ExCo<pair<A,B>,0>::toConvert(pair<A,B> &a, const pair<A2,B2> &b){ExOp::toConvert(a.first,b.first); ExOp::toConvert(a.second,b.second); return a;}

LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toAdd(pair<A,B> &a, const pair<A,B> &b){ExOp::toAdd(a.first,b.first); ExOp::toAdd(a.second,b.second); return a;}
LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toSubt(pair<A,B> &a, const pair<A,B> &b){ExOp::toSubt(a.first,b.first); ExOp::toSubt(a.second,b.second); return a;}
LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toMult(pair<A,B> &a, const pair<A,B> &b){ExOp::toMult(a.first,b.first); ExOp::toMult(a.second,b.second); return a;}
LFHTEMP pair<A,B>& ExCo<pair<A,B>,0>::toDivi(pair<A,B> &a, const pair<A,B> &b){ExOp::toDivi(a.first,b.first); ExOp::toDivi(a.second,b.second); return a;}

LFHTEMP pair<A,B> & ExCo<pair<A,B>,0>::toMemfree(pair<A,B> &a){ExOp::toMemfree(a.first); ExOp::toMemfree(a.second);return a;}
LFHTEMP pair<A,B> & ExCo<pair<A,B>,0>::toMemmove(pair<A,B> &a,pair<A,B> &b){ExOp::toMemmove(a.first, b.first); ExOp::toMemmove(a.second, b.second); return a;}


LFHTEMP void ExCo<pair<A,B>,0>::show(const pair<A,B> &a, FILE* f_out, int l){
		switch(l){
		case 0:
			ExOp::show(a.first,f_out,1);
			fprintf(f_out,"\n");
			ExOp::show(a.second,f_out,1);
			fprintf(f_out,"\n\n");
			break;
		case 1:
			ExOp::show(a.first,f_out,2);
			fprintf(f_out,"\t");
			ExOp::show(a.second,f_out,2);
			fprintf(f_out,"\n");
			break;
		case 2:
			fprintf(f_out,"[");
			ExOp::show(a.first,f_out,3);
			fprintf(f_out,";");
			ExOp::show(a.second,f_out,3);
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			ExOp::show(a.first,f_out,4);
			fprintf(f_out,",");
			ExOp::show(a.second,f_out,4);
			fprintf(f_out,")");
			break;
	}

	}

LFHTEMP string ExCo<pair<A,B>,0>::type_tostring(const pair<A,B> &a){return string("pair<") + ExOp::type_tostring(a.first) +string(",")+ ExOp::type_tostring(a.second) + string(">");}
/*
LFHTEMP pair<A,B>& ExCo<pair<A,B> >::toMemmove(pair<A,B> &a, pair<A,B> &b){ExOp::toMemmove(a.first,b.first);ExOp::toMemmove(a.second,b.second); return a;}
LFHTEMP pair<A,B>& ExCo<pair<A,B> >::toMemfree(pair<A,B> &a){ExOp::toMemfree(a.first);ExOp::toMemfree(a.second); return a;}*/

LFHTEMP ERRCODE ExCo<pair<A,B>,0>::save(const pair<A,B> &a, FILE*f){ERRCODE fout = ExOp::save(a.first,f); return fout | ExOp::save(a.second,f); }
LFHTEMP ERRCODE ExCo<pair<A,B>,0>::load(pair<A,B> &a, FILE*f, unsigned int l){return ExOp::load(a.first,f) | ExOp::load(a.second,f);}

	// ExFn specializations


template<class F, class A, class B, unsigned int dims> void ExFn<false>::compose(F &func, DataGrid<A,dims> &a, const DataGrid<B,dims> &b)	{
    a.setSizes(b.dims);
    typename DataGrid<A,dims>::KeyIterator ite = a.getKeyIterator();
    if (ite.first()) do{
        ExOp::comp(func,(A&) a(ite()),(const B&) b(ite()));
    } while(ite.next());
}
template<class A, class B> inline A& ExFn<false>::toBackMult(A& a, const B& b){
    A tmp = b * a; return ExOp::toClone(a,tmp);
}
template<class F, class A, class B, unsigned int DA, unsigned int DB> inline static void compose(F &func, DataGrid<A,DA> &a, const DataGrid<B,DB> &b){
     if (DA < DB){
         a.setSizes(b.dims + DB-DA);
                    typename DataGrid<A,DA>::KeyIterator ite = a.getKeyIterator();
        if (ite.first()) do{
            ExOp::comp(func,(A&)a(ite()),(B&)b(ite()));
            } while(ite.next());
         }else{
     }
}
template<class A> inline A ExFn<false>::mkZero(const A& what) {A fout; return ExOp::toZero(fout);}
template<class A> inline A ExFn<false>::mkOne(const A& what) {A fout; return ExOp::toOne(fout);}
template<class A> inline A ExFn<false>::mkRand(const A& what) {A fout; return ExOp::toRand(fout);}


template<class A> inline A& ExFn<false>::toAdd(A& a, const A& b, const A& c) {return a = ExOp::mkAdd(b,c);}
template<class A> inline A& ExFn<false>::toSubt(A& a, const A& b, const A& c) {return a = ExOp::mkSubt(b,c);}
template<class A> inline A& ExFn<false>::toMult(A& a, const A& b, const A& c) {return a = ExOp::mkMult(b,c);}
template<class A> inline A& ExFn<false>::toDivi(A& a, const A& b, const A& c) {return a = ExOp::mkDivi(b,c);}
template<class A,class B> inline A& ExFn<false>::toAdd(A& a, const B& b, const B& c) {return a = ExOp::mkAdd(b,c);}
template<class A,class B> inline A& ExFn<false>::toSubt(A& a, const B& b, const B& c) {return a = ExOp::mkSubt(b,c);}
template<class A,class B> inline A& ExFn<false>::toMult(A& a, const B& b, const B& c) {return a = ExOp::mkMult(b,c);}
template<class A,class B> inline A& ExFn<false>::toDivi(A& a, const B& b, const B& c) {return a = ExOp::mkDivi(b,c);}
template<class A,class B,class C> inline A& ExFn<false>::toAdd(A& a, const B& b, const C& c) {return a = ExOp::mkAdd(b,c);}
template<class A,class B,class C> inline A& ExFn<false>::toSubt(A& a, const B& b, const C& c) {return a = ExOp::mkSubt(b,c);}
template<class A,class B,class C> inline A& ExFn<false>::toMult(A& a, const B& b, const C& c) {return a = ExOp::mkMult(b,c);}
template<class A,class B,class C> inline A& ExFn<false>::toDivi(A& a, const B& b, const C& c) {return a = ExOp::mkDivi(b,c);}


template<class A, class B> inline static A mkBackInnerProd(const A& a, const B& b){return ExOp::mkBackMult(a,b);}
template<class A, class B> inline static A mkBackOuterProd(const A& a, const B& b){return ExOp::mkBackMult(a,b);}

template<class A> inline A& ExFn<false>::addMult(A& a, const A& b, const A& c){return a += ExOp::mkMult(b,c);}
template<class A> inline A& ExFn<false>::subtMult(A& a, const A& b, const A& c){return a -= ExOp::mkMult(b,c);}
template<class A,class B> inline A& ExFn<false>::addMult(A& a, const B& b, const B& c){return a += ExOp::mkMult(b,c);}
template<class A,class B> inline A& ExFn<false>::subtMult(A& a, const B& b, const B& c){return a -= ExOp::mkMult(b,c);}
template<class A,class B,class C> inline A& ExFn<false>::addMult(A& a, const B& b, const C& c){return a += ExOp::mkMult(b,c);}
template<class A,class B,class C> inline A& ExFn<false>::subtMult(A& a, const B& b, const C& c){return a -= ExOp::mkMult(b,c);}

template<class A> inline A& ExOp::toSubt(A& a){return ExOp::toNegative(a);}
template<class A> inline A& ExOp::toDivi(A& a){return ExOp::toInverse(a);}



#undef LFHTEMP
#define LFHTEMP template<class C>
LFHTEMP C Magscaled<C>::operator()()const{return value * exp(exponent);}
LFHTEMP Magscaled<C> Magscaled<C>::operator*(const Magscaled<C> &other)const{return Magscaled<C>(value * other.value, exponent + other.exponent);}
LFHTEMP Magscaled<C>& Magscaled<C>::operator*=(const Magscaled<C> &other){value *= other.value; exponent += other.exponent;return *this;}
LFHTEMP Magscaled<C> Magscaled<C>::operator/(const Magscaled<C> &other)const{return Magscaled<C>(value / other.value, exponent - other.exponent);}
LFHTEMP Magscaled<C>& Magscaled<C>::operator/=(const Magscaled<C> &other){value *= other.value; exponent += other.exponent;return *this;}
LFHTEMP Magscaled<C> Magscaled<C>::operator*(const C &other)const{return Magscaled<C>(value * other, exponent);}
LFHTEMP Magscaled<C>& Magscaled<C>::operator*=(const C &other){value *= other;return *this;}
LFHTEMP Magscaled<C> Magscaled<C>::operator/(const C &other)const{return Magscaled<C>(value / other, exponent);}
LFHTEMP Magscaled<C>& Magscaled<C>::operator/=(const C &other){value *= other;return *this;}
LFHTEMP Magscaled<C>& Magscaled<C>::monotonicAdd(const Magscaled<C> &other){value += other.value * exp(other.exponent - exponent); return *this;}
LFHTEMP Magscaled<C>& Magscaled<C>::operator+=(const Magscaled<C> &other){
    if (exponent > other.exponent) value += other.value * exp(other.exponent - exponent);
    else{
        value *= exp(exponent - other.exponent);
        value += other.value;
        exponent = other.exponent;
    }
return *this;}
LFHTEMP void Magscaled<C>::show(FILE*f, int level)const{fprintf(f,"%e * e^%e%c", value, exponent, (level ==0) ? '\n' : '\t');}

