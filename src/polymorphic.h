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



#ifndef _defined_Polymorphic
#define _defined_Polymorphic

#define GNU_SCIENTIFIC_LIBRARY
#ifndef LFH_HAS_RUNNINGTIME_STATISTICS
#define LFH_HAS_RUNNINGTIME_STATISTICS true
#endif

#include "core.h"


namespace LFHPrimitive{
// Very nice idea, but destructor gets called too early with explicit conversions...
template<class K, class D = void>
class Iterator{
    public:
    std::function<bool(K*&,D*&)> ite;
    K* key;
    D* data;
    Iterator(std::function<bool(K*&,D*&)> _ite): ite(_ite){}
    Iterator(std::function<bool(K*&,D*&)> _ite, K* _key): ite(_ite), key(_key){}
    Iterator(std::function<bool(K*&,D*&)> _ite, K* _key, D* _data): ite(_ite), key(_key), data(_data){}
    ~Iterator(){}
    //template<class F> Iterator(F &f);
    operator bool(){return data != NULL;}
    bool operator++(int){return ite(key, data);}
    K operator()(){return *key;}
    D* operator->(){return data;}
    D& operator*(){return *data;}
    template<class F> Iterator<K,D>& operator=(F &func);
};

/*
template<class K, class D = void>
class Iterator{
    public:
    std::function<bool(K&,D&)> ite;
    KeyElem<K,D> buffer;
    Iterator(std::function<bool(K&,D&)> _ite): ite(_ite){}
    template<class F> Iterator(F &f);
    operator bool(){
        return ite(buffer.k, buffer.d);
        }
    bool operator++(int){return ite(buffer.k, buffer.d);}
    std::remove_reference<K>::type& operator()(){return buffer.k;}
    D* operator->(){return &(buffer.d);}
    typename std::remove_reference<D>::type&  operator*(){return buffer.d;}
    template<class F> Iterator<K,D >& operator=(F &func);
};

template<class K, class D = void>
class Iterator<K, D&>{
    public:
    std::function<bool(K&,D&)> ite;
    KeyElem<K,D> buffer;
    Iterator(std::function<bool(K&,D&)> _ite): ite(_ite){}
    template<class F> Iterator(F &f);
    operator bool(){
        return ite(buffer.k, buffer.d);
        }
    bool operator++(int){return ite(buffer.k, buffer.d);}
    std::remove_reference<K>::type& operator()(){return buffer.k;}
    D* operator->(){return &(buffer.d);}
    typename std::remove_reference<D>::type&  operator*(){return buffer.d;}
    template<class F> Iterator<K,D >& operator=(F &func);
};*/


template<class C>
class Iterator<C, void>{
    public:
    typedef typename ExCo<C>::TYPE KeyType;
    std::function<bool(C&)> ite;
    C data;
    template<class F> Iterator(F &f);
    operator bool(){return ite(data);}
    bool operator++(int){return ite(data);}
    KeyType operator()(){return ExOp::getIndex(data);}
    C* operator->(){return &data;}
    C& operator*(){return data;}
};
template<class K, class D = void>
class IteratorMaker{
public:
    std::function<Iterator<K,D>(uint32_t i,uint32_t nb)> mkite;
    IteratorMaker(){}
    IteratorMaker(std::function<Iterator<K,D>(uint32_t i,uint32_t nb)> _mkite):mkite(_mkite){}
    IteratorMaker<K,D>& operator=(const IteratorMaker<K,D>& _mkitor){mkite = _mkitor.mkite; return *this;}
    Iterator<K,D> operator()(uint32_t i, uint32_t nb){return mkite(i,nb);}
};
template<class K, class D>
class Accessor{
    public:
    std::function<bool(const typename std::remove_reference<K>::type&,D*&)> accss;
    Accessor(){}
    Accessor(std::function<bool(const typename std::remove_reference<K>::type&,D*&)> _accss):accss(_accss){}
    Accessor<K,D>& operator=(const Accessor<K,D>& _accssor){accss = _accssor.accss; return *this;}
//    Accessor<K,D>& operator=(std::function<bool(const typename std::remove_reference<K>::type&,D*&)> _accss){accss = _accss; return *this;}
    bool operator()(const typename std::remove_reference<K>::type& key){D* data; return accss(key,*data);}
    D& operator[](const typename std::remove_reference<K>::type& key){D* data; accss(key,data); return *data;} // warning, assumes entry exists! use with care...
};


#include "polymorphic.hpp"
} // end of namespace LFHPrimitive



#endif


