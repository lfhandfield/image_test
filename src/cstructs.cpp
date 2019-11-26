/*
 * bastructs.cpp
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

#include "bastructs.h"

namespace LFHPrimitive{

ProgressBarPrint::ProgressBarPrint(uint32_t _length):length(_length),state(0){}
void ProgressBarPrint::start(const char* text){
    unsigned int i;
    for(i=0;i<length;i++) putc(' ', stdout);
    fprintf(stdout, ": %s", text);
    for(i=length + strlen(text)+2;i != 0;i--) putc('\b', stdout);fflush(stdout);
    lasttime = clock();
    state =0;
}
void ProgressBarPrint::update(double fraction){
    if (((unsigned int)(fraction * length * 16)) > state) {
    if (state & 15) putc('\b', stdout);
    state &= 0xFFFFFFF0;
        while (((unsigned int)(fraction * length * 16)) > (state | 15)) {state += 16; putc('#', stdout);}
        state |= ((unsigned int)(fraction * length * 16)) & 15;
        switch(state & 15){
            case 1: putc('.', stdout);break;
            case 2: putc(',', stdout);break;
            case 3: putc(':', stdout);break;
            case 4: putc(';', stdout);break;
            case 5: putc('i', stdout);break;
            case 6: putc('j', stdout);break;
            case 7: putc('l', stdout);break;
            case 8: putc('!', stdout);break;
            case 9: putc('|', stdout);break;
            case 10: putc('L', stdout);break;
            case 11: putc('C', stdout);break;
            case 12: putc('G', stdout);break;
            case 13: putc('O', stdout);break;
            case 14: putc('Q', stdout);break;
            case 15: putc('@', stdout);break;
        }
        if ((clock() - lasttime) & 0xFFFFFC00 ) {fflush(stdout); lasttime =clock();}
    }
}
void ProgressBarPrint::finish(){
    if (state & 15) putc('\b', stdout);
    while (length > (state >> 4)) {state+= 16;putc('#', stdout);}
    fprintf(stdout,"\n");fflush(stdout);
}



ANYTHING_MACRO_FOR_FUNCTION_DEF(char)
ANYTHING_MACRO_FOR_FUNCTION_DEF(int8_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(uint8_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(int16_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(uint16_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(int32_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(uint32_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(int64_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(uint64_t)
ANYTHING_MACRO_FOR_FUNCTION_DEF(float)
ANYTHING_MACRO_FOR_FUNCTION_DEF(double)



Sparsity::Sparsity(): asize(0){}
Sparsity::~Sparsity(){if ((asize != 0)&&(hash_encoding > 1)) delete[](data);}
Sparsity& Sparsity::toMemfree(){if ((asize != 0)&&(hash_encoding > 1)) delete[](data); asize = 0; return *this;}

uint32_t Sparsity::operator[](uint32_t index)const{uint32_t fout;
    switch(hash_encoding){
        case 0: return (index < asize) ? index : 0xFFFFFFFF;
        case 1: return (index != 0) ? 0xFFFFFFFF : dsize;
    }
    if ((index >> max_magn) != 0) return 0xFFFFFFFF;
    if ((fout = find(index)) == 0xFFFFFFFF) return (index < dsize) ? index : 0xFFFFFFFF;
    else return (index < dsize) ? 0 : 0xFFFFFFFF;
}
uint32_t Sparsity::genAlias(){

}
void Sparsity::rehash(uint32_t maxalias, uint32_t maxsize, uint8_t n_hash_magn){
    if (n_hash_magn==0){}
    // newhashsize_mag is at least 5 here
    uint8_t n_size_magn = ExOp::mkMagnitude(maxsize);
    uint8_t n_max_magn = ExOp::mkMagnitude(maxalias);
    uint8_t n_hash_bytes;


    // occ-bit  data (topkey  next) up to 32bits
    uint8_t tbit = (n_size_magn + n_max_magn + 1);
    if (tbit <= 16){
        n_hash_bytes = 0;
    }else{
        n_hash_bytes = (tbit <= 32)? 1:2;
    }
    uint16_t *n_data = new uint16_t[1<< (n_hash_magn + n_hash_bytes)]; LFH_NICE_ALLOCERROR(n_data,"")
/*
    switch(n_hash_bytes | (hash_bytes <<2) ){
        case 0:
            for(uint32_t i = 0; (i >> hash_magn) == 0;i++){
                uint32_t d = (i & (0xFFFFFFFF >> (32-hash_magn))) | ((1 << size_magn) - (1 << hash_magn)) & (data + (i << hash_bytes))
            }
        break;
    }*/
}

uint32_t Sparsity::find(uint32_t alias)const{
    if (asize == 0) return 0xFFFFFFFF;
    switch(hash_encoding){
    case 0: return (alias < asize) ? alias : 0xFFFFFFFF;
    case 1: return (alias == dsize) ? 0 : 0xFFFFFFFF;
    case 2:{//:  16 bit  1-data-topkey-next
        unsigned long mask = 0xFFFF >> (16 - hash_magn);
        unsigned long cur = mask & alias;
        unsigned long high_mask = mask ^ (0xFFFF >> max_magn);
        uint16_t ncur;
        if (((data16[ cur ] >> max_magn) & 1) == 0) return (alias < asize) ? alias : 0xFFFFFFFF; // not found default behavior
        while((alias ^ data16[ cur ]) & high_mask){
            if (data16[ cur] & mask == cur) return (alias < asize) ? alias : 0xFFFFFFFF; // not found default behavior
            cur = mask & data16[ cur];
        }
        return (uint32_t)(data16[ cur ] >> (max_magn + 1));
        } //:
    case 3:{//: 2x16 bit  data-1 // topkey-next
        unsigned long mask = 0xFFFF >> (16 - hash_magn);
        unsigned long cur = mask & (alias << 1);
        uint16_t ncur;
        if ((data16[cur | 1] & 1) == 0) return (alias < asize) ? alias : 0xFFFFFFFF; // not found default behavior
        while((alias ^ data16[cur << 1]) >> hash_magn){
            if ((data16[ cur] & mask) == cur) return (alias < asize) ? alias : 0xFFFFFFFF; // not found default behavior
            cur = mask & data16[cur << 1] ;
        }
        return (uint32_t)(data16[ (cur << 1) | 1] >> 1);
    } //:
    case 4:{//: 32 bit  data-1-topkey-next

    }//:
    case 5:{//: 2x32 bit  data-1// topkey-next

    }//:

    // hashsize < arraysize < maxkeysize

    case 6:{//: extreme size >= 0x80000000, there is no hashing at all... only in degenerate cases...

    }

    }
}



uint32_t Sparsity::addEntry(uint32_t alias){
    if (asize == 0) {
        if (alias == 0) {hash_encoding =0;asize++;}
        else {hash_encoding =1;asize =1; dsize = alias;}
        return 0;
    }
    switch(hash_encoding){
        case 0:
            if (alias == asize) {asize++; return alias;}
            // TODO
        break;
        case 1:
            ExCo<uint32_t>::mkMagnitude(dsize | alias);

        break;

    }
    if (asize == 0){
        if (alias == 0) {asize++; hash_encoding = 0;}
        else{
            hash_encoding = 2;
        }
        return 0;
    }
    switch(hash_magn){
    case 0:
        if (alias == asize) {asize= 1; dsize = 1; return dsize;}
        else{

        }
    break;
    case 1:
        // if (alias >> max_magn) rehash(alias,  );

    break;
    }
}

void Sparsity::show(FILE*f, int level)const{
    if (asize == 0) {fprintf(f,"(empty Sparsity)\n"); return;}
    fprintf(f,"Sparsity with %i elements\n", asize);
    if (auto ite = this->getIterator()) do {
        printf("%i:%i\n", ite(), *ite);
    }while(ite++);
}

Sparsity::ConstIterator::operator bool (){
    if (target.asize == 0) return false;
    offset=0;
    nvac = target.dsize;
    switch(target.hash_encoding){
        case 0: cur_value = 0; return true;
        case 1: cur_value = target.dsize; return true;
        case 2:
        mask = 0xFFFF >> (16 - target.hash_magn);
        high_mask = mask ^ (0xFFFF >> target.max_magn);
        break;
    }
    if (nvac != 0) cur_value=0;
    else{

    }
return true;}

bool Sparsity::ConstIterator::operator++(int){
    offset++;
    switch(target.hash_encoding){
        case 0: cur_value = offset; return offset < target.asize;
        case 1: return false;
        case 2:
            if (offset < nvac) {cur_value = offset; return true;}
            if (offset == target.asize) return false;
    }
return true;}



Anything& Anything::toZero(){memset(ptr,'\0',anytype & 0xFF);}
Anything& Anything::toOne(){
    switch(anytype ){
        case LFH_ANYTHING_CHAR: *(int8_t*)ptr = 1; break;
        case LFH_ANYTHING_UNSIGNED_CHAR: *(uint8_t*)ptr = 1; break;
        case LFH_ANYTHING_SHORT: *(int16_t*)ptr = 1; break;
        case LFH_ANYTHING_UNSIGNED_SHORT:  *(uint16_t*)ptr = 1; break;
        case LFH_ANYTHING_INT: *(int32_t*)ptr = 1; break;
        case LFH_ANYTHING_UNSIGNED_INT: *(uint32_t*)ptr = 1; break;
        case LFH_ANYTHING_LONG: *(int64_t*)ptr = 1; break;
        case LFH_ANYTHING_UNSIGNED_LONG: *(uint64_t*)ptr = 1; break;
        case LFH_ANYTHING_FLOAT: *(float*)ptr = 1.0f; break;
        case LFH_ANYTHING_DOUBLE: *(double*)ptr = 1.0; break;
    }
return *this;}


Anything::operator string()const{
    if ((anytype & 63) == 't') return string(*(char**)ptr);
    char buffer[32];
    switch(anytype ){
        case LFH_ANYTHING_CHAR: sprintf(buffer, "%i", *(int8_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_CHAR: sprintf(buffer, "%i", *(uint8_t*)ptr);
        case LFH_ANYTHING_SHORT: sprintf(buffer, "%i", *(int16_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_SHORT: sprintf(buffer, "%i", *(uint16_t*)ptr);
        case LFH_ANYTHING_INT: sprintf(buffer, "%i", *(int32_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_INT: sprintf(buffer, "%i", *(uint32_t*)ptr);
        case LFH_ANYTHING_LONG: sprintf(buffer, "%" PRId64, *(int64_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_LONG: sprintf(buffer, "%" PRIu64, *(uint64_t*)ptr);
        case LFH_ANYTHING_FLOAT: sprintf(buffer, "%f", *(float*)ptr);
        case LFH_ANYTHING_DOUBLE: sprintf(buffer, "%e", *(double*)ptr);
    }
return(string(buffer));}
AnythingConst::operator string()const{
    if ((anytype & 63) == 't') return string(*(const char**)ptr);
    char buffer[32];
    switch(anytype){
        case LFH_ANYTHING_CHAR: sprintf(buffer, "%i", *(const int8_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_CHAR: sprintf(buffer, "%i", *(const uint8_t*)ptr);
        case LFH_ANYTHING_SHORT: sprintf(buffer, "%i", *(const int16_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_SHORT: sprintf(buffer, "%i", *(const uint16_t*)ptr);
        case LFH_ANYTHING_INT: sprintf(buffer, "%i", *(const int32_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_INT: sprintf(buffer, "%i", *(const uint32_t*)ptr);
        case LFH_ANYTHING_LONG: sprintf(buffer, "%" PRId64, *(const int64_t*)ptr);
        case LFH_ANYTHING_UNSIGNED_LONG: sprintf(buffer, "%" PRIu64, *(const uint64_t*)ptr);
        case LFH_ANYTHING_FLOAT: sprintf(buffer, "%f", *(const float*)ptr);
        case LFH_ANYTHING_DOUBLE: sprintf(buffer, "%e", *(const double*)ptr);
    }
return(string(buffer));}


uint8_t* Vector<Anything>::push_back_routine(){
	uint32_t rsize = (asize & 0x7FFFFFFF);
	uint32_t i;
	if (((asize + 1)^(asize))> rsize){
		if (asize & 0x80000000) asize &= 0x7FFFFFFf;
		else{
            uint8_t* swp = new uint8_t[((rsize+1) << 1) & (anytype & 0xFF) ]; LFH_NICE_ALLOCERROR(swp,"")
            if (rsize> 0){
                memcpy(swp, darray, rsize & (anytype & 0xFF) );
                delete[](darray);
            }
            darray = swp;
		}
	}
	return(darray + (anytype & 0xFF) * ((asize++) & 0x7FFFFFFF));
}
Vector<Anything>& Vector<Anything>::setSize_init(uint32_t nsize) {
	asize = nsize;
	if (nsize == 0) return *this;
	uint32_t i =1;
	while(i<= nsize) i = (i<< 1);
	darray = new uint8_t[i * (anytype & 0xFF)]; LFH_NICE_ALLOCERROR(darray,"")
return *this;}
Vector<Anything>& Vector<Anything>::setSize(uint32_t nsize) {
	if (nsize == (asize &0x7FFFFFFF)) return *this;
	if (asize != 0) delete[](darray);
	asize = nsize;
	if (nsize == 0) return *this;
	uint32_t i =1;
	while(i<= nsize) i = (i<< 1);
	darray = new uint8_t[i * (anytype & 0xFF)];LFH_NICE_ALLOCERROR(darray,"")
return *this;}
Vector<Anything>& Vector<Anything>::toMemswap(Vector<Anything> & v){
	union{uint32_t swasize;uint8_t swanytype;uint8_t* swdarray;};
	swdarray = darray; darray = v.darray; v.darray = swdarray;
	swanytype = anytype; anytype = v.anytype; v.anytype = swanytype;
	swasize = asize; asize = v.asize; v.asize = swasize;
return(*this);}
ERRCODE Vector<Anything>::load(FILE* f, unsigned int length){
    uint32_t rsize;
    if (fread(&rsize, 1, sizeof(uint32_t),f) != 1) return 1;
    if (rsize == 0) return 0;
    if (fread(&anytype, 1, sizeof(uint32_t),f) != 1) return 1;
    this->setSize(rsize);
return fread(darray, anytype & 0xFF,rsize, f) == rsize ? 0 : 1;}
ERRCODE Vector<Anything>::save(FILE* f) const{
    uint32_t rsize = asize & 0x7FFFFFFF;
    if (fwrite(&rsize, 1, sizeof(uint32_t),f) != 1) return 1;
    if (rsize == 0) return 0;
    if (fwrite(&anytype, 1, sizeof(uint32_t),f) != 1) return 1;;
return fwrite(darray, (anytype & 0xFF),rsize, f) ==rsize ? 0 : 1;}

void Vector<Anything>::show(FILE* f, int level) const{
    if (asize == 0) {fprintf(f, "Empty AnyVector\n"); return;}
    fprintf(f, "AnyVector with  %i %i-byte elems bytestream:\n", asize & 0x7FFFFFFF, anytype & 0xFF);
    ExOp::show(darray);
}

OptionStruct::OptionStruct(){}
OptionStruct::OptionStruct(string value){data.addEntry(value);}
OptionStruct::OptionStruct(const char* value){if (value != NULL); data.addEntry(value);}
OptionStruct::OptionStruct(const std::nullptr_t value){}
bool OptionStruct::has(const char* query){return(data.find(query) != 0xFFFFFFFF);}
void OptionStruct::populateFlags(const char* str, char sep){
    // TODO
}

} // end of namespace

