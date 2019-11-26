/*
 * bastructs.h
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

#ifndef _defined_Cstruct
#define _defined_Cstruct

#define GNU_SCIENTIFIC_LIBRARY
#ifndef LFH_HAS_RUNNINGTIME_STATISTICS
#define LFH_HAS_RUNNINGTIME_STATISTICS true
#endif

#include "core.h"

using namespace std;
namespace LFHPrimitive{


class Weight{ // double with special meaning
public:
	double w;
	Weight(){}
	Weight(const double&v):w(v){}
	const double& operator()()const{return(w);}
	void show(FILE* f = stdout, int level=0)const{if (level) fprintf(f,"%f", w); else fprintf(f,"%f\n", w);}
};
class Zscore{
    public:
    double scaled_zscore; // ~ N(0, weight)
    double weight;
    Zscore()=default;
    Zscore(double val):scaled_zscore(val), weight(1.0){}
    Zscore(double val,double _weight):scaled_zscore(val), weight(_weight){}
    operator double(){ return scaled_zscore / sqrt(weight);}
    Zscore operator+(const Zscore& other )const{return Zscore(scaled_zscore + other.scaled_zscore, weight + other.weight);}
    Zscore operator-(const Zscore& other )const{return Zscore(scaled_zscore - other.scaled_zscore, weight + other.weight);}
    Zscore& operator+(const Zscore& other ){scaled_zscore += other.scaled_zscore; weight += other.weight; return *this;}
    Zscore& operator-(const Zscore& other ){scaled_zscore -= other.scaled_zscore; weight += other.weight; return *this;}
};

// Anything can be int double char*, and optionally have an integer index
// The type itself holds the size of an unit within the 0xFF byte
enum LFH_ANYTHING_enum{
    LFH_ANYTHING_CHAR= 0x1,
    LFH_ANYTHING_SHORT= 0x2,
    LFH_ANYTHING_INT= 0x4,
    LFH_ANYTHING_LONG= 0x8,
    LFH_ANYTHING_UNSIGNED_CHAR= 0x101,
    LFH_ANYTHING_UNSIGNED_SHORT= 0x102,
    LFH_ANYTHING_UNSIGNED_INT= 0x104,
    LFH_ANYTHING_UNSIGNED_LONG= 0x108,
    LFH_ANYTHING_CHARX= 0x200,
    LFH_ANYTHING_CHAR1= 0x201,
    LFH_ANYTHING_CHAR2= 0x202,
    LFH_ANYTHING_CHAR4= 0x203,
    LFH_ANYTHING_CHAR5= 0x204,
    LFH_ANYTHING_FLOAT= 0x304,
    LFH_ANYTHING_DOUBLE= 0x308
};

#define ANYTHING_MACRO_FOR_FUNCTION_DECL(InTyPe) explicit operator InTyPe()const;\
	InTyPe operator+(const InTyPe& val)const;\
	InTyPe operator-(const InTyPe& val)const;\
	InTyPe operator*(const InTyPe& val)const;\
	InTyPe operator/(const InTyPe& val)const;\
	Anything& operator=(const InTyPe& val); \
	Anything& operator+=(const InTyPe& val); \
	Anything& operator-=(const InTyPe& val); \
	Anything& operator*=(const InTyPe& val); \
	Anything& operator/=(const InTyPe& val);
#define ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(InTyPe) explicit operator InTyPe()const;\
	InTyPe operator+(const InTyPe& val)const;\
	InTyPe operator-(const InTyPe& val)const;\
	InTyPe operator*(const InTyPe& val)const;\
	InTyPe operator/(const InTyPe& val)const;


#define ANYTHING_MACRO_FOR_FUNCTION_DEF(InTyPe) Anything::operator InTyPe()const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return (InTyPe) *(int8_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return (InTyPe) *(uint8_t*)ptr;\
        case LFH_ANYTHING_SHORT: return (InTyPe) *(int16_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return (InTyPe) *(uint16_t*)ptr;\
        case LFH_ANYTHING_INT: return (InTyPe) *(int32_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_INT: return (InTyPe) *(uint32_t*)ptr;\
        case LFH_ANYTHING_LONG: return (InTyPe) *(int64_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_LONG: return (InTyPe) *(uint64_t*)ptr;\
        case LFH_ANYTHING_FLOAT: return (InTyPe) *(float*)ptr;\
        case LFH_ANYTHING_DOUBLE: return (InTyPe) *(double*)ptr;\
    }}; AnythingConst::operator InTyPe()const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return (InTyPe) *(const int8_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return (InTyPe) *(const uint8_t*)ptr;\
        case LFH_ANYTHING_SHORT: return (InTyPe) *(const int16_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return (InTyPe) *(const uint16_t*)ptr;\
        case LFH_ANYTHING_INT: return (InTyPe) *(const int32_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_INT: return (InTyPe) *(const uint32_t*)ptr;\
        case LFH_ANYTHING_LONG: return (InTyPe) *(const int64_t*)ptr;\
        case LFH_ANYTHING_UNSIGNED_LONG: return (InTyPe) *(const uint64_t*)ptr;\
        case LFH_ANYTHING_FLOAT: return (InTyPe) *(const float*)ptr;\
        case LFH_ANYTHING_DOUBLE: return (InTyPe) *(const double*)ptr;\
}} InTyPe Anything::operator+(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) + val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)+val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) +val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)+val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)+val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)+val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)+val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)+val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)+val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)+val;\
}return val;} InTyPe AnythingConst::operator+(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) + val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)+val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) +val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)+val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)+val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)+val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)+val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)+val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)+val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)+val;\
}return val;} InTyPe Anything::operator-(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) - val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)-val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) -val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)-val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)-val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)-val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)-val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)-val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)-val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)-val;\
}return val;} InTyPe AnythingConst::operator-(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) - val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)-val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) -val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)-val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)-val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)-val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)-val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)-val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)-val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)-val;\
}return val;} InTyPe Anything::operator*(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) * val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)*val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) *val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)*val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)*val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)*val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)*val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)*val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)*val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)*val;\
}return val;} InTyPe AnythingConst::operator*(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) * val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)*val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) *val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)*val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)*val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)*val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)*val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)*val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)*val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)*val;\
}return val;} InTyPe Anything::operator/(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) / val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)/val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) /val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)/val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)/val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)/val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)/val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)/val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)/val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)/val;\
}return val;} InTyPe AnythingConst::operator/(const InTyPe& val)const{switch(anytype){\
        case LFH_ANYTHING_CHAR: return ((InTyPe)*(int8_t*)ptr) / val;\
        case LFH_ANYTHING_UNSIGNED_CHAR: return ((InTyPe)*(uint8_t*)ptr)/val;\
        case LFH_ANYTHING_SHORT: return ((InTyPe)*(int16_t*)ptr) /val;\
        case LFH_ANYTHING_UNSIGNED_SHORT: return ((InTyPe)*(uint16_t*)ptr)/val;\
        case LFH_ANYTHING_INT: return ((InTyPe)*(int32_t*)ptr)/val;\
        case LFH_ANYTHING_UNSIGNED_INT: return ((InTyPe)*(uint32_t*)ptr)/val;\
        case LFH_ANYTHING_LONG: return ((InTyPe)*(int64_t*)ptr)/val;\
        case LFH_ANYTHING_UNSIGNED_LONG: ((InTyPe)*(uint64_t*)ptr)/val;\
        case LFH_ANYTHING_FLOAT: return ((InTyPe)*(float*)ptr)/val;\
        case LFH_ANYTHING_DOUBLE:return ((InTyPe)*(double*)ptr)/val;\
}return val;}Anything& Anything::operator=(const InTyPe& val){switch(anytype){\
        case LFH_ANYTHING_CHAR: *(int8_t*)ptr = val;  break;\
        case LFH_ANYTHING_UNSIGNED_CHAR: *(uint8_t*)ptr = val; break;\
        case LFH_ANYTHING_SHORT: *(int16_t*)ptr = val; break;\
        case LFH_ANYTHING_UNSIGNED_SHORT: *(uint16_t*)ptr = val; break;\
        case LFH_ANYTHING_INT: *(int32_t*)ptr = val; break;\
        case LFH_ANYTHING_UNSIGNED_INT: *(uint32_t*)ptr = val; break;\
        case LFH_ANYTHING_LONG: *(int64_t*)ptr = val; break;\
        case LFH_ANYTHING_UNSIGNED_LONG: *(uint64_t*)ptr = val; break;\
        case LFH_ANYTHING_FLOAT: *(float*)ptr = val; break;\
        case LFH_ANYTHING_DOUBLE: *(double*)ptr = val; break;\
}return *this;} Anything& Anything::operator+=(const InTyPe& val){switch(anytype){\
        case LFH_ANYTHING_CHAR: *(int8_t*)ptr += val; break;\
        case LFH_ANYTHING_UNSIGNED_CHAR: *(uint8_t*)ptr += val; break;\
        case LFH_ANYTHING_SHORT: *(int16_t*)ptr += val; break;\
        case LFH_ANYTHING_UNSIGNED_SHORT: *(uint16_t*)ptr += val; break;\
        case LFH_ANYTHING_INT: *(int32_t*)ptr += val; break;\
        case LFH_ANYTHING_UNSIGNED_INT: *(uint32_t*)ptr += val; break;\
        case LFH_ANYTHING_LONG: *(int64_t*)ptr += val; break;\
        case LFH_ANYTHING_UNSIGNED_LONG: *(uint64_t*)ptr += val; break;\
        case LFH_ANYTHING_FLOAT: *(float*)ptr += val; break;\
        case LFH_ANYTHING_DOUBLE: *(double*)ptr += val; break;\
}return *this;} Anything& Anything::operator-=(const InTyPe& val){switch(anytype){\
        case LFH_ANYTHING_CHAR: *(int8_t*)ptr -= val; break;\
        case LFH_ANYTHING_UNSIGNED_CHAR: *(uint8_t*)ptr -= val; break;\
        case LFH_ANYTHING_SHORT: *(int16_t*)ptr -= val; break;\
        case LFH_ANYTHING_UNSIGNED_SHORT: *(uint16_t*)ptr -= val; break;\
        case LFH_ANYTHING_INT: *(int32_t*)ptr -= val; break;\
        case LFH_ANYTHING_UNSIGNED_INT: *(uint32_t*)ptr -= val; break;\
        case LFH_ANYTHING_LONG: *(int64_t*)ptr -= val; break;\
        case LFH_ANYTHING_UNSIGNED_LONG: *(uint64_t*)ptr -= val; break;\
        case LFH_ANYTHING_FLOAT: *(float*)ptr -= val; break;\
        case LFH_ANYTHING_DOUBLE: *(double*)ptr -= val; break;\
}return *this;} Anything& Anything::operator*=(const InTyPe& val){switch(anytype){\
        case LFH_ANYTHING_CHAR: *(int8_t*)ptr *= val; break;\
        case LFH_ANYTHING_UNSIGNED_CHAR: *(uint8_t*)ptr *= val; break;\
        case LFH_ANYTHING_SHORT: *(int16_t*)ptr *= val; break;\
        case LFH_ANYTHING_UNSIGNED_SHORT: *(uint16_t*)ptr *= val; break;\
        case LFH_ANYTHING_INT: *(int32_t*)ptr *= val; break;\
        case LFH_ANYTHING_UNSIGNED_INT: *(uint32_t*)ptr *= val; break;\
        case LFH_ANYTHING_LONG: *(int64_t*)ptr *= val; break;\
        case LFH_ANYTHING_UNSIGNED_LONG: *(uint64_t*)ptr *= val; break;\
        case LFH_ANYTHING_FLOAT: *(float*)ptr *= val; break;\
        case LFH_ANYTHING_DOUBLE: *(double*)ptr *= val; break;\
}return *this;} Anything& Anything::operator/=(const InTyPe& val){switch(anytype){\
        case LFH_ANYTHING_CHAR: *(int8_t*)ptr /= val; break;\
        case LFH_ANYTHING_UNSIGNED_CHAR: *(uint8_t*)ptr /= val; break;\
        case LFH_ANYTHING_SHORT: *(int16_t*)ptr /= val; break;\
        case LFH_ANYTHING_UNSIGNED_SHORT: *(uint16_t*)ptr /= val; break;\
        case LFH_ANYTHING_INT: *(int32_t*)ptr /= val; break;\
        case LFH_ANYTHING_UNSIGNED_INT: *(uint32_t*)ptr /= val; break;\
        case LFH_ANYTHING_LONG: *(int64_t*)ptr /= val; break;\
        case LFH_ANYTHING_UNSIGNED_LONG: *(uint64_t*)ptr /= val; break;\
        case LFH_ANYTHING_FLOAT: *(float*)ptr /= val; break;\
        case LFH_ANYTHING_DOUBLE: *(double*)ptr /= val; break;\
}return *this;}

class Anything{ // special structure that only interprets a void* as requested, does not hold/own any data!
public:
	uint8_t* ptr;
	uint32_t anytype; // pP tT g are pointers and strings, ***memory is never owned***
	Anything(uint8_t* _ptr, uint32_t _anytype): ptr(_ptr), anytype(_anytype){}
	ANYTHING_MACRO_FOR_FUNCTION_DECL(char)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(int8_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(int16_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(int32_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(int64_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(uint8_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(uint16_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(uint32_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(uint64_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(float)
	ANYTHING_MACRO_FOR_FUNCTION_DECL(double)
	void operator=(const string& val);
	Anything& toZero();
	Anything& toOne();
	template<class C> void operator=(C* val);
	explicit operator string()const;
};
class AnythingConst{ // special structure that only interprets a void* as requested, does not hold/own any data!
public:
    const uint8_t* ptr;
    uint32_t anytype;
	AnythingConst(){}
	AnythingConst(const uint8_t* _ptr, uint32_t _anytype): ptr(_ptr), anytype(_anytype){}
 	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(char)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(int8_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(int16_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(int32_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(int64_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(uint8_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(uint16_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(uint32_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(uint64_t)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(float)
	ANYTHING_MACRO_FOR_FUNCTION_DECL_CONST(double)
    explicit operator string()const;
};

class ProgressBarPrint{
    public:
    unsigned int state;
    unsigned int length;
    unsigned int lasttime;
    ProgressBarPrint(uint32_t _length);
    void start(const char*);
    void update(double fraction);
    void update(int nbdone, int total){return this->update(((double)nbdone) / total);}
    void finish();
};

} // end of namespace LFHPrimitive

#endif


