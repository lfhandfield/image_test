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

// contains KeyElem, Tuple, TMatrix, matriangle..., Complex, Quaternion, WeightElem

template<class KEY, class DATA,bool isMax> ExtremumScope<KEY,DATA,isMax>& ExtremumScope<KEY,DATA,isMax>::toZero(){ExOp::toMin(best); return(*this);}
template<class KEY, class DATA,bool isMax> void ExtremumScope<KEY,DATA,isMax>::init(const KEY& key, const DATA& data){if (ExOp::isValid(data)) {best  =data; best_key = key;}else ExOp::toMin(best);}
template<class KEY, class DATA,bool isMax> void ExtremumScope<KEY,DATA,isMax>::regist(const KEY& key, const DATA& data){if ((ExOp::isValid(data))&&(best < data)) {best  =data; best_key = key;}}
template<class KEY, class DATA> ExtremumScope<KEY,DATA,false>& ExtremumScope<KEY,DATA,false>::toZero(){ExOp::toMax(best);return(*this);}
template<class KEY, class DATA> void ExtremumScope<KEY,DATA,false>::init(const KEY& key, const DATA& data){if (ExOp::isValid(data)) {best  =data; best_key = key;}else ExOp::toMax(best);}
template<class KEY, class DATA> void ExtremumScope<KEY,DATA,false>::regist(const KEY& key, const DATA& data){if ((ExOp::isValid(data))&&(best > data)) {best  =data; best_key = key;}}

#undef LFHTEMP
#define LFHTEMP template <class C>

LFHTEMP AsyncStructure<C>::ReadAccess::operator bool(){
    int preval = target.sema.fetch_add(1);
    if (preval < 256) {state =1; return true;}
    else{target.sema.fetch_add(-1); state =0; return false;}
}
LFHTEMP AsyncStructure<C>::WriteAccess::operator bool(){
    int preval = target.sema.fetch_add(256);
    if (preval == 256) {state =256; return true;}
    else{target.sema.fetch_add(-256); state =0; return false;}
}


#undef LFHTEMP
#define LFHTEMP template <class C, int S>
LFHTEMP AsyncBuffer<C,S>::AsyncBuffer(): sema_wr(0), sema_rd(0), pos_wr(0u), pos_rd(0u), sema(0) {}

LFHTEMP AsyncBuffer<C,S>::ReadAccess::ReadAccess(AsyncBuffer& _target) :target(_target), mlock(_target.mtx){
    //target.cv.wait(mlock, [this]{return true;});
}

LFHTEMP bool AsyncBuffer<C,S>::ReadAccess::needsUpdate(){
    if (target.pos_rd != target.pos_wr) return(true);
    target.overmtx.lock();
    if (target.overflow.empty()) {target.overmtx.unlock(); return false;}
    do{ // tries to empty the overflow stack!
        uint32_t pos = target.sema.fetch_add(1);
        if (((pos - target.pos_rd) >> S) == 0){
            target.buffer[pos & ((1 << S) - 1)] = target.overflow.top();
            target.overflow.pop();
            target.pos_wr++;
        }else{
            target.sema.fetch_sub(1);
            target.overmtx.unlock();
            return true;
        }
    }while(!target.overflow.empty());
    target.overmtx.unlock();
    return true;
}
LFHTEMP C AsyncBuffer<C,S>::ReadAccess::update(){return target.buffer[(target.pos_rd++) & ((1 << S) - 1)];}

LFHTEMP void AsyncBuffer<C,S>::insert(const C &what){
    std::lock_guard<std::mutex> dalck(overmtx);
    uint32_t pos = sema.fetch_add(1);
    if (((pos - pos_rd) >> S) == 0){
        buffer[pos & ((1 << S) - 1)] = what;
        pos_wr++;
    }else{
        sema.fetch_sub(1);
        overflow.push(what);
    }
}

LFHTEMP ERRCODE AsyncBuffer<C,S>::rdSync(C& where){ // only use if no other thread can read;
    if (pos_rd == pos_wr) return 1; // buffer is empty
    ExOp::toMemmove(where, buffer[pos_rd]);
    pos_rd = (pos_rd + 1) & ((1 << S) - 1);
    return 0;
}
LFHTEMP ERRCODE AsyncBuffer<C,S>::rdAsync(C& where){
    if (pos_rd == pos_wr) return 1; // buffer is empty
    int pos = fetch_and_add(&sema_rd,1);
    if (pos){sema_rd--; return ERRCODE_BUSY;} // concurrent read
    ExOp::toMemmove(where, buffer[pos_rd]);
    pos_rd = (pos_rd + 1) & ((1 << S) - 1);
    sema_rd--;
    return 0;
}
LFHTEMP ERRCODE AsyncBuffer<C,S>::wrSync(const C& what){ // only use if no other thread can write;
    if ((((pos_wr + 1) ^ pos_rd) & ((1 << S) - 1)) == 0) return 1; // buffer is full
    buffer[pos_wr] = what;
    pos_wr = (pos_wr + 1) & ((1 << S) - 1);
    return 0;
}
LFHTEMP ERRCODE AsyncBuffer<C,S>::wrAsync(const C& what){
    if ((((pos_wr + 1) ^ pos_rd) & ((1 << S) - 1)) == 0) return 1; // buffer is full
    int pos = fetch_and_add(&sema_wr,1);
    if (pos){sema_wr--; return ERRCODE_BUSY;}// concurrent write
    buffer[pos_wr] = what;
    pos_wr = (pos_wr + 1) & ((1 << S) - 1);
    sema_wr--;
    return 0;
}
LFHTEMP ERRCODE AsyncBuffer<C,S>::mvSync(C& what){ // only use if no other thread can write;
    if ((((pos_wr + 1) ^ pos_rd) & ((1 << S) - 1)) == 0) return 1; // buffer is full
    ExOp::toMemmove(buffer[pos_wr],what);
    pos_wr = (pos_wr + 1) & ((1 << S) - 1);
    return 0;
}
LFHTEMP ERRCODE AsyncBuffer<C,S>::mvAsync(C& what){
    if ((((pos_wr + 1) ^ pos_rd) & ((1 << S) - 1)) == 0) return 1; // buffer is full
    int pos = fetch_and_add(&sema_wr,1);
    if (pos){sema_wr--; return ERRCODE_BUSY;}// concurrent write
    ExOp::toMemmove(buffer[pos_wr],what);
    pos_wr = (pos_wr + 1) & ((1 << S) - 1);
    sema_wr--;
    return 0;
}

	// class KeyElem
#undef LFHTEMP
#define LFHTEMP template <class C, class B>

LFHTEMP KeyElem<C,B>::KeyElem(const C& _k, const B& _n): k(_k),d(_n){}
/** \brief Checks that both Key and Data fields are 'VALID' using ExOp::isValid();
 */
LFHTEMP bool KeyElem<C,B>::isValid() const{return ( ExOp::isValid(k) && ExOp::isValid(d) );}
/** \brief Outputs text description of object
 *
 * \param File stream
 * \param Integer selecting embedded objects character separators
 *
 */
LFHTEMP void KeyElem<C,B>::show(FILE* out, int level) const{
		switch(level){
		case 0:
			fprintf(out,"Key: "); ExOp::show(k,out, 1);
			fprintf(out,"\nData: "); ExOp::show(d,out, 1);
			fprintf(out,"\n");
			break;
		case 1:
			fprintf(out,"Key: ["); ExOp::show(k,out, 2);
			fprintf(out,"]\tData: ["); ExOp::show(d,out, 2);
			fprintf(out,"]");
			break;
		case 2:
			fprintf(out,"K{"); ExOp::show(k,out, 3);
			fprintf(out,"}#D{"); ExOp::show(d,out, 3);
			fprintf(out,"}");
			break;
		case 3:
			fprintf(out,"("); ExOp::show(k,out, 4);
			fprintf(out,");("); ExOp::show(d,out, 4);
			fprintf(out,")");
			break;
		}
	}

#undef LFHTEMP
#define LFHTEMP template <class C,unsigned int TSIZE,Tuple_flag Cflag>
LFHTEMP double Tuple<C,TSIZE,Cflag>::weight(){return(data[TSIZE]);}


LFHTEMP void Tuple<C,TSIZE,Cflag>::zero(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toZero(data[i]);}
LFHTEMP void Tuple<C,TSIZE,Cflag>::random(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toRand(data[i]);}

LFHTEMP template<class D> Tuple<C,TSIZE,Cflag>::operator Tuple<D,TSIZE> () const{
		Tuple<D,TSIZE> out;
		for(unsigned int i=0;i<TSIZE;i++) out.data[i] = (D)data[i];
		return(out);
		}

LFHTEMP double Tuple<C,TSIZE,Cflag>::pnorm() const{
		double sum = ExOp::pnorm(data[0]);
		for(unsigned int i =1;i<TSIZE;i++) sum += ExOp::pnorm(data[i]);
		return(sum);
	}

LFHTEMP double Tuple<C,TSIZE,Cflag>::norm() const{
		double sum = ExOp::pnorm(data[0]);
		for(unsigned int i =1;i<TSIZE;i++) sum += ExOp::pnorm(data[i]);
		return(sqrt(sum));
	}


LFHTEMP template<class O> std::vector<O>& Tuple<C,TSIZE,Cflag>::wrStdVector(std::vector<O>& fout)const{
    fout.clear();
    fout.reserve(TSIZE);
    for(uint32_t i =0;i<TSIZE;i++) fout.push_back(data[i]);
return fout;}

LFHTEMP ERRCODE Tuple<C,TSIZE,Cflag>::save( FILE *f) const {
	if (ExCo<C>::IsPOD) return (fwrite(data,sizeof(C),TSIZE,f)== TSIZE) ? 0 : 1;
	ERRCODE fout = ExOp::save(data[0],f);
	for(unsigned int i=1;i<TSIZE;i++) fout |= ExOp::save(data[i],f);
}
LFHTEMP ERRCODE Tuple<C,TSIZE,Cflag>::load( FILE *f, unsigned int ch_TSIZE) { ERRCODE fout =0;
	if (ExCo<C>::IsPOD) fout |= fread(data,sizeof(C),TSIZE,f) == TSIZE ? 0 : 1;
	else for(unsigned int i=0;i<TSIZE;i++) fout |= ExOp::load(data[i],f, ch_TSIZE/TSIZE);
	return fout;
}

LFHTEMP bool Tuple<C,TSIZE,Cflag>::isValid() const{
unsigned int i;
for(i=0;i<TSIZE;i++) if (!ExOp::isValid(data[i])) break;
return(i == TSIZE);
}


LFHTEMP template<class O> auto Tuple<C,TSIZE,Cflag>::operator+(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkAdd((*this)[0], other[0])),TSIZE,Cflag>{
    Tuple< decltype( ExOp::mkAdd((*this)[0], other[0])),TSIZE,Cflag> f_out;
    for(unsigned int i=0;i<TSIZE;i++) f_out[i] = ExOp::mkAdd((*this)[i],other[i]);
return( f_out );}
LFHTEMP template<class O> auto Tuple<C,TSIZE,Cflag>::operator-(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkSubt((*this)[0], other[0])),TSIZE,Cflag>{
    Tuple< decltype( ExOp::mkSubt((*this)[0], other[0])),TSIZE,Cflag> f_out;
    for(unsigned int i=0;i<TSIZE;i++) f_out[i] = ExOp::mkSubt((*this)[i], other[i]);
return( f_out );}
LFHTEMP template<class O> auto Tuple<C,TSIZE,Cflag>::operator*(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkMult((*this)[0], other[0])),TSIZE,Cflag>{
    Tuple< decltype( ExOp::mkMult((*this)[0], other[0])),TSIZE,Cflag> f_out;
    for(unsigned int i=0;i<TSIZE;i++) f_out[i] = ExOp::mkMult((*this)[i], other[i]);
return( f_out );}
LFHTEMP template<class O> auto Tuple<C,TSIZE,Cflag>::operator/(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkDivi((*this)[0], other[0])),TSIZE,Cflag>{
    Tuple< decltype( ExOp::mkDivi((*this)[0], other[0])),TSIZE,Cflag> f_out;
    for(unsigned int i=0;i<TSIZE;i++) f_out[i] = ExOp::mkDivi((*this)[i], other[i]);
return( f_out );}
LFHTEMP void Tuple<C,TSIZE,Cflag>::show( FILE *f_out, int l) const{
	switch(l){
		case 1:
		case 0:
			for(unsigned int i=0;i<TSIZE;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(data[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
		break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<TSIZE;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(data[i],f_out,3);
			}
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<TSIZE;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(data[i],f_out,4);
			}
			fprintf(f_out,")");
			break;
	}
}

LFHTEMP string Tuple<C,TSIZE,Cflag>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",TSIZE); return string("Tuple<") + ExOp::type_tostring(data[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}

LFHTEMP GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >  Tuple<C,TSIZE,Cflag>::mkgaussstat(double &w) const{
    typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE outer;
    unsigned int i,j,k;
    if (ExCo<C>::IS_COMMUTATIVE::ans){
        for(i=0,k=0; i< TSIZE;i++) for(j=i;j<TSIZE;j++,k++) outer[k] = data[i] * data[j];
    }else{
        for(i=0; i< TSIZE;i++) for(j=0;j<TSIZE;j++) outer[i + j * TSIZE] = data[i] * data[j];
    }
    return GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >(*this, outer ,w);
    }

LFHTEMP static double overlap(const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&a, const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&b){

    return 0;
}


//LFHTEMP Tuple<C,TSIZE,Cflag>::Tuple(C const & other){	int i;	for(i = TSIZE-1;i>=0;i--) data[i] = other;}

//LFHTEMP
//Tuple<C,TSIZE,Cflag>::Tuple(Tuple<C,TSIZE,Cflag> const & clonefrom){memcpy(data,clonefrom.data, TSIZEof(C)*(TSIZE));}

LFHTEMP C& Tuple<C,TSIZE,Cflag>::operator[](int const pos) {return(data[pos]);}
LFHTEMP const C& Tuple<C,TSIZE,Cflag>::operator[](int const pos) const {return(data[pos]);}

LFHTEMP template<class O,unsigned int oTSIZE, Tuple_flag Oflag>
char Tuple<C,TSIZE,Cflag>::compare(const Tuple<O,oTSIZE, Oflag>& other) const{
	unsigned int i=0;
	unsigned int m = (TSIZE < oTSIZE) ? TSIZE : oTSIZE;
	for(;i<m;i++) if ((*this)[i] != other[i]) break;
	if (i == m) return( (TSIZE == oTSIZE) ? 0 : ((TSIZE < oTSIZE)  ? 2 : 1) );
	else return( ((*this)[i] < other[i]) ? 2 : 1 );
}

// Super Operators!
LFHTEMP
template<class I> void Tuple<C,TSIZE,Cflag>::operator() (Oper1<I> const & op){ // not a match
for(unsigned int i = 0;i < TSIZE;i++) data[i](op);
}
LFHTEMP
void Tuple<C,TSIZE,Cflag>::operator() (Oper1< C> const & op){ // match
for(unsigned int i = 0;i < TSIZE;i++)  op(data[i], data[i]);
}
LFHTEMP
template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2 ){ // not a match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::mkMaximum();i--) data[i](op,a_2.data[i]);
}

LFHTEMP
template<class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2){ // match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::mkMaximum();i--) op(data[i], a_2.data[i]);
}
LFHTEMP
template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2 ){ // not a match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::mkMaximum();i--) data[i](op,a_2.data[i]);
}
LFHTEMP
template<class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2){ // match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::mkMaximum();i--) op(data[i], a_2.data[i]);
}

	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3 ){ // not a match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) data[i](op,a_2.data[i],a_3.data[i]);
	}
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3 ){ // not a match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) data[i](op,a_2.data[i],a_3.data[i]);
	}
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> & a_3 ){ // not a match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) data[i](op,a_2.data[i],a_3.data[i]);
	}
	LFHTEMP
	template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3){ // match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) op(data[i], a_2.data[i], a_3.data[i]);
	}
	LFHTEMP
	template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3){ // match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) op(data[i], a_2.data[i], a_3.data[i]);
	}

	LFHTEMP
	template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> & a_3){ // match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) op(data[i], a_2.data[i], a_3.data[i]);
	}



LFHTEMP
	void Tuple<C,TSIZE,Cflag>::fourierTransform_routine(){
		// assumes TSIZE is a power of two, and C can be multiplied by a mycomplex number
		int step;
		int cur;
		int icur;
		mycomplex multip[TSIZE];
		multip[0] = mycomplex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  ((double)icur * M_PI) / step;
				multip[icur] = mycomplex(cos(tmpdouble),sin(tmpdouble));
			//	ExOp::show(multip[icur]);
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
	//					printf("%f\t%f\t",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
						C tmp = data[cur | step | icur];
					//	ExOp::show(tmp);
					//	ExOp::show(multip[step | icur]);
						tmp *= multip[step | icur];
					//	printf("gives\n");
					//	ExOp::show(tmp);
			//			printf("%f\t%f\t",tmp[0],tmp[1]);
						tmp += data[cur | icur];
				//		printf("%f\t%f\t",tmp[0],tmp[1]);
						data[cur | step | icur] *= multip[icur];
						data[cur | icur] += data[cur | step | icur];
						data[cur | step | icur] = tmp;
//						printf("%f\t%f\n",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
				}
			}
		}
	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::invfourierTransform_routine(){
		// assumes TSIZE is a power of two, and C can be multiplied by a mycomplex number
		int step;
		int cur;
		int icur;
		mycomplex multip[TSIZE];

		multip[0] = mycomplex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  -((double)icur * M_PI) / step;
				multip[icur] = mycomplex(cos(tmpdouble),sin(tmpdouble));
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
	//					printf("%f\t%f\t",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
						C tmp = data[cur | step | icur];
						tmp *= multip[step | icur];
			//			printf("%f\t%f\t",tmp[0],tmp[1]);
			//			data[cur | icur] *= half;
						tmp += data[cur | icur];
				//		printf("%f\t%f\t",tmp[0],tmp[1]);
						data[cur | step | icur] *= multip[icur];
						data[cur | icur] += data[cur | step | icur];
						data[cur | step | icur] = tmp;
//						printf("%f\t%f\n",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
				}
			}
		}
		mycomplex half = mycomplex(1/((double)TSIZE),0.0f);
		for(cur =0;cur<TSIZE;cur++) data[cur] *= half;

	}

LFHTEMP Tuple<C,TSIZE,Cflag>&  Tuple<C,TSIZE,Cflag>::toIntPows(const C& v){
    data[0] = v;
    unsigned int mp =1;
    for(unsigned int i=1;i<TSIZE;i++){
        data[i] = data[mp -1] * data[i- mp];
        if (mp<<1 == i+1) mp <<=1;
    }
return(*this);}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::fourierTransform_routine2(){
		// assumes TSIZE is a power of two, and C can be multiplied by a mycomplex number
		int step;
		int icur;
        step = TSIZE >> 1;
            for(icur=0;icur<step;icur++){
                    C tmp = data[icur] - data[step | icur];
                    data[ icur] += data[step | icur];
                    data[step | icur] = tmp;
            }
		if (TSIZE <= 2) return;
		mycomplex damultip;
		double tmpdouble;
		for(step=step >> 1;step != 0 ;step=step >> 1){
            // run
			for(icur =0;icur<TSIZE;icur+= step){
                tmpdouble = (M_PI * icur) / TSIZE;
			    damultip = mycomplex(cos(tmpdouble),sin(tmpdouble));

				for(;(icur & step) == 0 ; icur++){
                        data[ step | icur] *= damultip;
                        C tmp = data[icur];
                        tmp -= data[ step | icur];
						data[icur] += data[step | icur];
						data[step | icur] = tmp;
				}
			}
		}
	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::invfourierTransform_routine2(){
		// assumes TSIZE is a power of two, and C can be multiplied by a mycomplex number
				int step;
		int icur;

		if (TSIZE > 2) {
		mycomplex damultip;
		double tmpdouble;
		for(step= 1 ; (step << 1) < TSIZE ;step=step << 1){
            // run
			for(icur =0;icur<TSIZE;icur+= step){
                tmpdouble = (-M_PI * icur) / TSIZE;
			    damultip = mycomplex(cos(tmpdouble),sin(tmpdouble));
				for(;(icur & step) == 0 ; icur++){
                        C tmp = data[icur];
                        tmp -= data[ step | icur];
						data[icur] += data[step | icur];
						data[step | icur] = tmp;
						data[ step | icur] *= damultip;
				}
			}
		}
		}
		step = TSIZE >> 1;
            for(icur=0;icur<step;icur++){
                    C tmp = data[icur] - data[step | icur];
                    data[ icur] += data[step | icur];
                    data[step | icur] = tmp;
            }
		mycomplex half = mycomplex(1/((double)TSIZE),0.0f);
		for(icur =0;icur<TSIZE;icur++) data[icur] *= half;
	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::invfourierTransform_routine(mycomplex* array){
		// assumes TSIZE is a power of two, and C can be multiplied by a mycomplex number
		int step;
		int cur;
		int icur;
		mycomplex multip[TSIZE];

		multip[0] = mycomplex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  -((double)icur * M_PI) / step;
				multip[icur] = mycomplex(cos(tmpdouble),sin(tmpdouble));
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
						C tmp = array[cur | step | icur];
						tmp *= multip[step | icur];
						tmp += array[cur | icur];
						array[cur | step | icur] *= multip[icur];
						array[cur | icur] += array[cur | step | icur];
						array[cur | step | icur] = tmp;
				}
			}
		}
		mycomplex half = mycomplex(1/((double)TSIZE),0.0f);
		for(cur =0;cur<TSIZE;cur++) array[cur] *= half;

	}



LFHTEMP
	void Tuple<C,TSIZE,Cflag>::fourierTransform_routine(mycomplex* array){
		// assumes TSIZE is a power of two, and C can be multiplied by a mycomplex number
		int step;
		int cur;
		int icur;
		mycomplex multip[TSIZE];
		multip[0] = mycomplex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  ((double)icur * M_PI) / step;
				multip[icur] = mycomplex(cos(tmpdouble),sin(tmpdouble));
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
						C tmp = array[cur | step | icur];
						tmp *= multip[step | icur];
						tmp += array[cur | icur];
						array[cur | step | icur] *= multip[icur];
						array[cur | icur] += array[cur | step | icur];
						array[cur | step | icur] = tmp;
				}
			}
		}
	}


LFHTEMP
Tuple<mycomplex, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::bluesteinWindow(int subTSIZE){ // the uutput should be cached...
	Tuple<C, TSIZE,Cflag> out;
	int i;

	int j,s;
	j=0;
			for(i=0;i<subTSIZE;i++) {
				double ang = i*i*M_PI/ subTSIZE;
				out[j] =  mycomplex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE-subTSIZE+1;i++) {
				out[j] =  ExCo<C>::zero();
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE;i++) {
				double ang = (TSIZE-i)*(TSIZE-i)*M_PI/ subTSIZE;
				out[j] =  mycomplex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

		out.fourierTransform_routine();

	return(out);
	}

LFHTEMP
void Tuple<C,TSIZE,Cflag>::bluesteinWindow(mycomplex *& fout, int subTSIZE){ // the uutput should be cached...
	int i;
    fout = new mycomplex[TSIZE]; LFH_NICE_ALLOCERROR(fout,"")
	int j,s;
	j=0;
			for(i=0;i<subTSIZE;i++) {
				double ang = i*i*M_PI/ subTSIZE;
				fout[j] =  mycomplex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE-subTSIZE+1;i++) {
				fout[j] =  ExCo<C>::zero();
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE;i++) {
				double ang = (TSIZE-i)*(TSIZE-i)*M_PI/ subTSIZE;
				fout[j] =  mycomplex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

		Tuple<C,TSIZE,Cflag>::fourierTransform_routine(fout);

	}

LFHTEMP
template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::fourierTransform(const Tuple<mycomplex, superTSIZE,Cflag>& bluewindow) const{

		// implement Bluestein algorithm!

		// B windows : 2 {1.0+1.0i, 1.0-1.0i }
		// B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}


		Tuple<C, superTSIZE,Cflag> tmp;
		int i;
		int j=0;
		int s;

		for(i=0;i<TSIZE;i++){
				double ang = -i * i * M_PI / TSIZE;
				mycomplex factor = mycomplex(cos(ang),sin(ang));
				tmp[j] = data[i];
				tmp[j] *= factor;
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
		for(;i<superTSIZE;i++){
				tmp[j] = ExCo<C>::zero();
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		tmp.fourierTransform_routine();
//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

		Tuple<C, superTSIZE,Cflag> tmp2;
//		Tuple<C, superTSIZE,Cflag> tmpb = bluesteinWindow<superTSIZE>();

//		for(i=0;i<superTSIZE;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

		for(i=0;i<superTSIZE;i++){
				tmp2[j]	= tmp[i] * bluewindow[i];
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		tmp2.invfourierTransform_routine();
//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

		Tuple<C, TSIZE,Cflag> out;
		for(i=0;i<TSIZE;i++) {
			 double ang = -i * i * M_PI / TSIZE;
			mycomplex factor = mycomplex(cos(ang),sin(ang));
			 out[i] = tmp2[i];
			 out[i] *= factor;

			 }
		return(out);
	}

LFHTEMP
template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::invfourierTransform(const Tuple<mycomplex, superTSIZE,Cflag>& bluewindow) const{

		// implement Bluestein algorithm!

		// B windows : 2 {1.0+1.0i, 1.0-1.0i }
		// B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}


		Tuple<C, superTSIZE,Cflag> tmp;
		int i;
		int j=0;
		int s;

		for(i=0;i<TSIZE;i++){
				double ang = i * i * M_PI / TSIZE;
				mycomplex factor = mycomplex(cos(ang),sin(ang));
				tmp[j] = data[i];
				tmp[j] *= factor;
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
		for(;i<superTSIZE;i++){
				tmp[j] = ExCo<C>::zero();
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		tmp.fourierTransform_routine();
//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

		Tuple<C, superTSIZE,Cflag> tmp2;
//		Tuple<C, superTSIZE,Cflag> tmpb = bluesteinWindow<superTSIZE>();

//		for(i=0;i<superTSIZE;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

		for(i=0;i<superTSIZE;i++){
				tmp2[j]	= tmp[i];
				mycomplex factor =bluewindow[i];
				factor[1] = -factor[1];
				tmp2[j]	*= factor;
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		tmp2.invfourierTransform_routine();
//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

		Tuple<C, TSIZE,Cflag> out;
		for(i=0;i<TSIZE;i++) {
			 double ang = i * i * M_PI / TSIZE;
			mycomplex factor = mycomplex(cos(ang)/TSIZE,sin(ang)/TSIZE);
			 out[i] = tmp2[i];
			 out[i] *= factor;

			 }
		return(out);
	}

	LFHTEMP
	C Tuple<C,TSIZE,Cflag>::max() const{
		C _out = data[0];
		int i;
		for(i=1;i<TSIZE;i++) if (data[i] > _out) _out = data[i];
		return(_out);
	}
	LFHTEMP
	C Tuple<C,TSIZE,Cflag>::min() const{
		C _out = data[0];
		int i;
		for(i=1;i<TSIZE;i++) if (data[i] < _out) _out = data[i];
		return(_out);
	}

	LFHTEMP
	Tuple<C, TSIZE, Cflag> Tuple<C,TSIZE,Cflag>::normalize() const{
		C sum = data[0];
		Tuple<C, TSIZE, Cflag> _out;
		int i;
		for(i=1;i<TSIZE;i++) sum += data[i];
//		C sumi = ExCo<C>::invert(sum);
		for(i=0;i<TSIZE;i++) if (!ExCo<C>::isValid(_out[i] = data[i] / sum)) break;
		if (i == TSIZE) return(_out);
		// sum or sumi is NAN, normalize by the max beforehand

		C cax = max();
		C cix = min();
		sum = (cax > -cix) ? cax : -cix;
		if (sum  == ExCo<C>::zero()) {
			// prior guess
			_out[i] =  ExCo<C>::one() * (1.0f / TSIZE);
			for(i=1;i<TSIZE;i++) _out[i] = _out[0];
			return(_out);// all entries are zeroes
		}
		for(i=0;i<TSIZE;i++) _out[i] = data[i] / sum;
		sum = _out[0];
		for(i=1;i<TSIZE;i++) sum += _out[1];
		for(i=0;i<TSIZE;i++) _out[i] /= sum;
		return(_out);
		}
	LFHTEMP
	Tuple<C, TSIZE, Cflag> Tuple<C,TSIZE,Cflag>::normalize(C const & norm) const{
		C sum = data[0];
		Tuple<C, TSIZE, Cflag> _out;
		int i;
		for(i=1;i<TSIZE;i++) sum += data[i];
		C sumi = (ExCo<C>::invert(sum)) * norm;
		for(i=0;i<TSIZE;i++) _out[i] = data[i] * sumi;
		return(_out);
		}


LFHTEMP
template<unsigned int superTSIZE> Tuple<Complex<C>, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::fourierTransformReal() const{
		//Tuple<Complex<C>, TSIZE,Cflag>& _out;

		//Tuple<Complex<C>, TSIZE,Cflag>& _out;

		//int s

	}





	/*
template <class C, int nbdim>
const SetComparison& HyperPositionProjectCube<C,nbdim>::compare(const HyperPosition<C,nbdim> & other)  const {

	}
	*/
/*
LFHTEMP
template<class I, class O, class MC> Tuple<MC,TSIZE> Tuple<C,TSIZE,Cflag>::operator()(Oper<I,O> const & op) const{
	Tuple<MC,TSIZE> _out;
	int i;
	for(i=0;i<TSIZE;i++) _out = data[i](op);
	return(_out);
}

LFHTEMP
template<class I, class O> Tuple<O,TSIZE> Tuple<C,TSIZE,Cflag>::operator()(Oper< Tuple<I,TSIZE> ,Tuple<O,TSIZE> > const & op){
	Tuple<O,TSIZE> _out;
	int i;
	for(i=0;i<TSIZE;i++) _out = op(data[i]);
	return(_out);
}

LFHTEMP
template<class O> O Tuple<C,TSIZE,Cflag>::operator()(Oper< Tuple<C,TSIZE,Cflag> ,O> const & op) const {return(op(*this));}
	*/

LFHTEMP
Tuple<C,TSIZE,Cflag>::Tuple(C const* const clonefrom){memcpy(data,clonefrom, sizeof(C)*(TSIZE));}

LFHTEMP
template<class O, unsigned int oTSIZE>
Tuple<C,TSIZE,Cflag>::Tuple(Tuple<O,oTSIZE,Cflag> const & other){(*this) = other;}

LFHTEMP
const Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-() const{
	Tuple<C,TSIZE,Cflag> tmp;
	int i;
	for(i = TSIZE-1;i>=0;i--) tmp[i] = -data[i];
	return(Tuple<C,TSIZE,Cflag>(tmp));
	}

LFHTEMP Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(const C (&other)[TSIZE]){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) data[i] = other[i];
	return(*this);
	}

LFHTEMP Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(Tuple<C,TSIZE,Cflag> const & other){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) data[i] = other[i];
	return(*this);
	}

LFHTEMP Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::toMemmove(Tuple<C,TSIZE,Cflag> & other){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) ExOp::toMemmove(data[i],other[i]);
	return(*this);
}

LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(O const & other){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) data[i] = other;
	return(*this);
	}

LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(const O (&other)[TSIZE]){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) data[i] = other;
	return(*this);
	}

//LFHTEMP Tuple<C,TSIZE,Cflag>::operator Vector<C,LFHVECTOR_REMOTE> (){return Vector<C,LFHVECTOR_REMOTE>(data,TSIZE);}


LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator+=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] += other;
	return(*this);
	}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] -= other;
	return(*this);
	}

LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] *= other;
	return(*this);
}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] /= other;
	return(*this);
	}

LFHTEMP template <class O, unsigned int OSIZE, Tuple_flag OFLAG> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(const Tuple<O,OSIZE,OFLAG> & other){
	uint32_t dmax = other.getSize();
	dmax = (dmax > TSIZE) ? TSIZE : dmax;
	for(uint32_t i=0 ;i<dmax;i++) data[i] = other[i];
return(*this);}

	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator+=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] += other.data[i];
		return(*this);
	}
	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] -= other.data[i];
		return(*this);
	}

	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] *= other.data[i];
		return(*this);
	}
	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] /= other.data[i];
		return(*this);
	}

	LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(TMatrix<O,TSIZE,TSIZE> const & other){
		Tuple<C,TSIZE,Cflag> clo(*this);
		unsigned int i,k;

		for(i=0;i<TSIZE;i++) data[i] = other.data[i * TSIZE] * clo[0];
		for(k=1;k<TSIZE;k++){
			for(i=0;i<TSIZE;i++) data[i] += other.data[k + i * TSIZE] * clo[k];
		}
		return(*this);
	}

	LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(TMatrix<O,TSIZE,TSIZE> const & other){
		return(*this);
	}



LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator+(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) += other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator-(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) -= other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator*(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) *= other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator/(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) /= other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator+=(KeyElem<unsigned int, O> const & other) {data[other.k] += other.d; return(*this);}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-=(KeyElem<unsigned int, O> const & other) {data[other.k] -= other.d; return(*this);}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(KeyElem<unsigned int, O> const & other) {data[other.k] *= other.d; return(*this);}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(KeyElem<unsigned int, O> const & other) {data[other.k] /= other.d; return(*this);}


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	Tuple<C,TSIZE,Cflag> const &  Tuple<C,TSIZE,Cflag>::genConvolution_Gaussian(double std){
		Tuple<C,TSIZE,Cflag> out;
		int i,j;
		double tmp;
		for(i=0;i<TSIZE;i++){
			tmp = (i -(TSIZE >>1)) /std;
			out.data[i] = (C)(exp(-tmp*tmp/2.0f));
			for(j=0;j<TSIZE*100000000;j++) tmp += j;
			}

		return(out);
		}

			// OPERATOR!!!


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int pos>
	void Tuple<C,TSIZE,Cflag>::Selector<pos>::operator()(C & _out, Tuple<C,TSIZE,Cflag> & input) const{
		_out = input[pos];
	}


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::MakeTuple::operator()(Tuple<C,TSIZE,Cflag> & output, C & _in) const{
		for(unsigned int i=0;i<TSIZE;i++)	output[i] = _in;
	}

template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length ){
    unsigned int i;
    if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2))) return;
    C s = data[0] * vec[0];
    for(i=1;i<length;i++) s+= data[i] * vec[i];
    s = ExOp::mkTrju(s) / (-denum2);
    for(i=0;i<length;i++) data[i] += s * vec[i];
}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int fTSIZE>
	void Tuple<C,TSIZE,Cflag>::Concatenate<fTSIZE>::operator()(Tuple<C,TSIZE,Cflag> & _out, Tuple<C,fTSIZE,Cflag> & input, Tuple<C,TSIZE-fTSIZE,Cflag> & input2) const{
		int i;
		for(i=0;i<fTSIZE;i++) _out.data[i] = input[i];
		for(;i<TSIZE;i++) _out.data[i] = input2[i-fTSIZE];

	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int fTSIZE>
	void Tuple<C,TSIZE,Cflag>::Deconcatenate<fTSIZE>::operator()(Tuple<C,TSIZE,Cflag> & _out, Tuple<C,fTSIZE,Cflag> & out_2, Tuple<C,TSIZE+fTSIZE,Cflag> & input) const{
		int i;
		for(i=0;i<TSIZE;i++) _out[i] = input.data[i];
		for(;i<TSIZE+ fTSIZE;i++) out_2[i-TSIZE] = input.data[i];

	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int pos> void Tuple<C,TSIZE,Cflag>::SelectorWrite<pos>::operator()(Tuple<C,TSIZE,Cflag> & _out, C & _input) const{
		_out[pos] = _input;
	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int TSIZEin, unsigned int pos_in, unsigned int pos_out> void Tuple<C,TSIZE,Cflag>::SelectSelectorWrite<TSIZEin,pos_in,pos_out>::operator()(Tuple<C,TSIZE,Cflag> & _out, Tuple<C,TSIZEin,Cflag> & _input) const{
		_out[pos_out] = _input[pos_in];
	}

	//template<class C> void Square<C>::operator()(C & _out, C & input) const{_out =input*input;}
	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::L2Norm::operator()(C & _out, Tuple<C,TSIZE,Cflag> & input) const{
		_out = input[0] *input[0];
		int i;
		for(i=1;i<TSIZE;i++) _out += input[1] * input[1];
	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::FiniteDifferenceAssign::operator()(Tuple<C,TSIZE,Cflag> &_out) const{
		C tmp = _out[0] - _out[TSIZE-1];
		int i;
		for(i=0;i<TSIZE-1;i++){
			_out[i] = _out[i+1] - _out[i];
			}
			_out[i] = tmp;
		}

template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::FiniteDifference::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C,TSIZE,Cflag> & input) const{
    _out[TSIZE-1] = input[0] - input[TSIZE-1];
    int i;
    for(i=0;i<TSIZE-1;i++){
        _out[i] = input[i+1] - input[i];
    }
}
template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::ArgMax::operator()(int &_out, Tuple<C,TSIZE,Cflag> & input) const{
    int b =0;
    C curb = input[0];
    int i;
    for(i=1;i<TSIZE;i++){
        if (input[i] > curb) {curb = input[i]; b= i; }
        }
    _out = b;
}
template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::MassCenter::operator()(double &_out, Tuple<C,TSIZE,Cflag> & input) const{
    double tm = 0.0f;
    double sum = 0.0f;
    double tmp;
    int i;
    for(i=0;i<TSIZE;i++){
        tmp = (double) input[i]*input[i];
        sum += pow(2.0f,0.5f*(i+1-TSIZE)) * tmp;
        tm += tmp;
        }
    _out = sum / tm;
}
template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::MassCenter2::operator()(Tuple<double,3,Cflag> &_out, Tuple<C,TSIZE,Cflag> & input) const{
    double tm = (double) input[0];
    double sum = 0.0f;
    double tmp;
    int i;
    for(i=1;i<TSIZE;i++){
        tmp = (double) input[i];
        sum += i * tmp;
        tm += tmp;
        }
    _out[0] = sum / tm;
    _out[1] = tm;
    _out[2] = sum / tm;
}
template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::PopWeight::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C, TSIZE+1,Cflag> & input) const{
    int i;
    for(i=0;i<TSIZE;i++) _out[i] = input[i] / input[TSIZE];
}
template<class C, unsigned int TSIZE, Tuple_flag Cflag>
void Tuple<C,TSIZE,Cflag>::PushWeight::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C, TSIZE-1,Cflag> & input, C & w) const{
    int i;
    _out[TSIZE-1] = w;
    for(i=0;i<TSIZE-1;i++) _out[i] = input[i] * w;
}
template<class C, unsigned int TSIZE, Tuple_flag Cflag>
template<unsigned int winTSIZE>
void Tuple<C,TSIZE,Cflag>::Convolution<winTSIZE>::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C,TSIZE,Cflag> & input, Tuple<C, winTSIZE,Cflag> & convwind) const{
    int i,j;
    for(i=0;i<TSIZE;i++){
        if (i < (winTSIZE >> 1)) j=(winTSIZE >> 1) - i;
        else j = 0;
        _out[i] = input[i+j-(winTSIZE >> 1)] * convwind[j];
        for(j++;(j<winTSIZE)&&(i+j-(winTSIZE >> 1) < TSIZE);j++) _out[i] += input[i+j-(winTSIZE >> 1)] * convwind[j];
    }
}

LFHTEMP template<class I> Tuple<I> Tuple<C,TSIZE,Cflag>::mkOrdering()const{ Tuple<I> fout; fout.setSize(TSIZE);
return *this;}
LFHTEMP Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::toInverse(){
    for(unsigned int i=0;i< TSIZE;i++) ExOp::toInverse(data[i]);
return *this;}
LFHTEMP Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::mkInverse()const{ Tuple<C,TSIZE,Cflag> fout;
    for(unsigned int i=0;i< TSIZE;i++) fout.data[i] = ExOp::mkInverse(data[i]);
return fout;}
LFHTEMP Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::mkPowInt(int a)const{ Tuple<C,TSIZE,Cflag> fout;
    for(unsigned int i=0;i< TSIZE;i++) fout.data[i] = ExOp::mkPowInt(data[i],a);
return fout;}
LFHTEMP Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::mkPowInvInt(int a)const{ Tuple<C,TSIZE,Cflag> fout;
    for(unsigned int i=0;i< TSIZE;i++) fout.data[i] = ExOp::mkPowInvInt(data[i],a);
return fout;}
LFHTEMP void Tuple<C,TSIZE,Cflag>::wrDeterminant(C& fout)const{fout = data[0];
    for(unsigned int i=1;i< TSIZE;i++) fout *= data[i];
}
LFHTEMP template<class O> void Tuple<C,TSIZE,Cflag>::wrDeterminant(O& fout)const{ExOp::wrDeterminant(data[0], fout);
    O tmp;
    for(unsigned int i=1;i< TSIZE;i++) {ExOp::wrDeterminant(data[i], tmp); fout *= tmp;}
}
LFHTEMP void Tuple<C,TSIZE,Cflag>::wrTrace(C& fout)const{fout = data[0];
    for(unsigned int i=1;i< TSIZE;i++) fout += data[i];
}
LFHTEMP template<class O> void Tuple<C,TSIZE,Cflag>::wrTrace(O& fout)const{ExOp::wrTrace(data[0], fout);
    O tmp;
    for(unsigned int i=1;i< TSIZE;i++) {ExOp::wrTrace(data[i], tmp); fout += tmp;}
}
#undef LFHTEMP
#define LFHTEMP template<class C, Tuple_flag Cflag>

LFHTEMP Tuple<C,0u,Cflag>::Tuple(const Tuple<C,0u, Cflag> &other){setSize_init(other.tup_size); for(unsigned int i=0;i<tup_size;i++) data[i] = other.data[i];}
LFHTEMP template<unsigned int OSIZE> Tuple<C,0u,Cflag>::Tuple(const Tuple<C,OSIZE, Cflag> &other){tup_size = OSIZE; data = new C[OSIZE]; LFH_NICE_ALLOCERROR(data,"") for(unsigned int i=0;i<OSIZE;i++) data[i] = other.data[i];}
LFHTEMP template<class O> Tuple<C,0u,Cflag>::Tuple(const std::vector<O> &other){setSize_init(other.size()); for(uint32_t i=0;i<tup_size;i++) data[i] = other[i];}
LFHTEMP template<class O> Tuple<C,0u,Cflag>::Tuple(const Vector<O> &other){setSize_init(other.getSize()); for(uint32_t i=0;i<tup_size;i++) data[i] = other[i];}
LFHTEMP Tuple<C,0u,Cflag>::Tuple(std::initializer_list<C> dalist){setSize_init(dalist.size()); C* tar = data; const C* sour = std::begin(dalist); for(uint32_t i=0;i<tup_size;i++) new(tar++) C(*(sour++));}
LFHTEMP template<class O> Tuple<C,0u,Cflag>::Tuple(std::initializer_list<O> dalist){setSize_init(dalist.size()); C* tar = data; const O* sour = std::begin(dalist); for(uint32_t i=0;i<tup_size;i++) *(tar++) = C(*(sour++));}
LFHTEMP const C& Tuple<C,0u,Cflag>::operator[](unsigned int index)const{return data[index];}
LFHTEMP C& Tuple<C,0u,Cflag>::operator[](unsigned int index){return data[index];}
LFHTEMP void Tuple<C,0u,Cflag>::setSize_init(unsigned int s){
    tup_size = s;
    if (s !=0) {
        data = new C[tup_size];
        LFH_NICE_ALLOCERROR(data,"Tuple construction failed");
    }
}
LFHTEMP Tuple<C,0u,Cflag>& Tuple<C,0u,Cflag>::setSize(unsigned int s) {
    if (s == 0){
        if (tup_size) delete[](data);
        data = nullptr;
    }else if (s != tup_size) {
        if (tup_size) delete[](data);
        data = new C[s];  LFH_NICE_ALLOCERROR(data,"%i bytes of memory could not be allocated\n", (uint32_t)(tup_size * sizeof(C)))
    }
    tup_size=s;
return *this;}
LFHTEMP Tuple<C,0u,Cflag>& Tuple<C,0u,Cflag>::toResize(unsigned int s) {
    if (s == 0){
        if (tup_size) delete[](data);
        data = nullptr;
    }else if (s != tup_size) {
        C* ndata = new C[s];LFH_NICE_ALLOCERROR(ndata,"")
        int m = (tup_size < s) ? tup_size : s;
        for(int i = 0; i<m;i++) ExOp::toMemmove(ndata[i],data[i]);
        if (tup_size) delete[](data);
        data = ndata;
        LFH_OPT_CHECK(if ((tup_size)&&(data == NULL)) {fprintf(stderr, "%i bytes of memory could not be allocated\n", (uint32_t)(tup_size * sizeof(C)));exit(1);})
    }
    tup_size=s;
return *this;}
LFHTEMP Tuple<typename Tuple<C,0u,Cflag>::INDEX_TYPE> Tuple<C,0u,Cflag>::mkOrdering()const{ Tuple<typename Tuple<C,0u,Cflag>::INDEX_TYPE> fout; fout.setSize(tup_size);
    uint32_t i,j;

    if (tup_size < 5){
        if (tup_size < 3){
            if (tup_size == 1) fout[0] = 0;
            else if (tup_size == 2){
                if (data[0] < data[1]) {fout[0]= 0; fout[1] =1;}
                else {fout[0]= 1; fout[1] =0;}
            }
        }else if (tup_size  == 3){
            i = (data[0] < data[1]) ? 0 :1;
            fout[0] = (data[2] < data[i]) ? 2 : i;
            i = i^1;
            fout[2] = (data[i] < data[2]) ? i : 2;
            fout[1] = 3 - fout[2] - fout[0];
        }else{
            i = (data[0] < data[1]) ? 0 :1;
            j = (data[2] < data[3]) ? 2 :3;
            if (data[i] < data[j]) {fout[0] = i; fout[1] = j;}
            else  {fout[0] = j; fout[1] = i;}
            if (data[i^1] > data[j^1]) {fout[3] = i^1; j = j^1;}
            else  {fout[3] = j^1; j = i^1;}
            if (data[1] > data[j]) { data[2] = data[1]; data[1] = j;}
            else data[2] = j;
        }
        return(fout);
    }

    for(i=1; i < tup_size;){
        if (data[i] > data[i^1]) {
                fout[i] = i;
                i++;
                fout[i] = i;
        }else{
                fout[i] = i^1;
                i++;
                fout[i] = i^1;
        }
        i++;
    }
    i ^= 1;
    if (i < tup_size) fout[i] = i;
    uint32_t magn = 2;
    while( (tup_size >> magn) != 0){



    }



return *this;}
LFHTEMP Tuple<C,0u,Cflag>& Tuple<C,0u,Cflag>::toInverse(){
    for(unsigned int i=0;i< tup_size;i++) ExOp::toInverse(data[i]);
return *this;}
LFHTEMP Tuple<C,0u,Cflag> Tuple<C,0u,Cflag>::mkInverse()const{ Tuple<C,0u,Cflag> fout;
    fout.setSize(this->getSize());
    for(unsigned int i=0;i< tup_size;i++) fout.data[i] = ExOp::mkInverse(data[i]);
return fout;}
LFHTEMP Tuple<C,0u,Cflag> Tuple<C,0u,Cflag>::mkPowInt(int a)const{ Tuple<C,0u,Cflag> fout;
    fout.setSize(this->getSize());
    for(unsigned int i=0;i< tup_size;i++) fout.data[i] = ExOp::mkPowInt(data[i],a);
return fout;}
LFHTEMP Tuple<C,0u,Cflag> Tuple<C,0u,Cflag>::mkPowInvInt(int a)const{ Tuple<C,0u,Cflag> fout;
    fout.setSize(this->getSize());
    for(unsigned int i=0;i< tup_size;i++) fout.data[i] = ExOp::mkPowInvInt(data[i],a);
return fout;}
LFHTEMP void Tuple<C,0u,Cflag>::wrDeterminant(C& fout)const{
    if (tup_size == 0){ExOp::toOne(fout); return;}
    fout = data[0];
    for(unsigned int i=1;i< tup_size;i++) fout *= data[i];
}
LFHTEMP template<class O> void Tuple<C,0u,Cflag>::wrDeterminant(O& fout)const{ if (tup_size == 0){ExOp::toOne(fout); return fout;}
    ExOp::wrDeterminant(data[0], fout);
    O tmp;
    for(unsigned int i=1;i<tup_size;i++) {ExOp::wrDeterminant(data[i], tmp); fout *= tmp;}
}
LFHTEMP void Tuple<C,0u,Cflag>::wrTrace(C& fout)const{
    if (tup_size == 0){ExOp::toZero(fout); return;}
    fout = data[0];
    for(unsigned int i=1;i< tup_size;i++) fout += data[i];
}
LFHTEMP template<class O> void Tuple<C,0u,Cflag>::wrTrace(O& fout)const{ if (tup_size == 0){ExOp::toZero(fout); return fout;}
    ExOp::wrTrace(data[0], fout);
    O tmp;
    for(unsigned int i=1;i<tup_size;i++) {ExOp::wrTrace(data[i], tmp); fout += tmp;}
}
LFHTEMP C& Tuple<C,0u,Cflag>::push_back(){
    unsigned int i;
    if (tup_size == 0){
        data = new C[1]; LFH_NICE_ALLOCERROR(data,"")
        tup_size = 1;
    }else{
        C* tmp = data;
        data = new C[tup_size+1]; LFH_NICE_ALLOCERROR(data,"")
        for(i=0;i<tup_size;i++) ExOp::toMemmove(data[i],tmp[i]);
        delete[](tmp);tup_size++;
    }
return data[tup_size - 1];}
LFHTEMP void Tuple<C,0u,Cflag>::pop_back(){
    unsigned int i;
    if (tup_size < 2){
        if (tup_size == 1) delete[](data);
        tup_size = 0;
    }else{
        C* tmp = data; tup_size--;
        data = new C[tup_size]; LFH_NICE_ALLOCERROR(data,"")
        for(i=0;i<tup_size;i++) ExOp::toMemmove(data[i],tmp[i]);
        delete[](tmp);
    }
}
LFHTEMP void Tuple<C,0u,Cflag>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length ){
    unsigned int i;
    if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2)))  return;
    C s = data[0] * vec[0];
    for(i=1;i<length;i++) s+= data[i] * vec[i];
    s = ExOp::mkTrju(s) / (-denum2);
    for(i=0;i<length;i++) data[i] += s * vec[i];
}
LFHTEMP template<class O, uint32_t TSIZE> auto Tuple<C,0u,Cflag>::operator+(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkAdd((*this)[0u], other[0u])),0u,Cflag>{
    Tuple< decltype( ExOp::mkAdd((*this)[0u], other[0u])),0u,Cflag> f_out;
    for(unsigned int i=0;i<other.getSize();i++) f_out[i] = (*this)[i] + other[i];
return( f_out );}
LFHTEMP template<class O, uint32_t TSIZE> auto Tuple<C,0u,Cflag>::operator-(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkSubt((*this)[0u], other[0u])),0u,Cflag>{
    Tuple< decltype( ExOp::mkSubt((*this)[0u], other[0u])),0u,Cflag> f_out;
    for(unsigned int i=0;i<other.getSize();i++) f_out[i] = (*this)[i] - other[i];
return( f_out );}
LFHTEMP template<class O, uint32_t TSIZE> auto Tuple<C,0u,Cflag>::operator*(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkMult((*this)[0u], other[0u])),0u,Cflag>{
    Tuple< decltype( ExOp::mkMult((*this)[0u], other[0u])),0u,Cflag> f_out;
    for(unsigned int i=0;i<other.getSize();i++) f_out[i] = (*this)[i] * other[i];
return( f_out );}
LFHTEMP template<class O, uint32_t TSIZE> auto Tuple<C,0u,Cflag>::operator/(Tuple<O,TSIZE,Cflag> const & other) const
 -> Tuple< decltype( ExOp::mkDivi((*this)[0u], other[0u])),0u,Cflag>{
    Tuple< decltype( ExOp::mkDivi((*this)[0u], other[0u])),0u,Cflag> f_out;
    for(unsigned int i=0;i<other.getSize();i++) f_out[i] = (*this)[i] / other[i];
return( f_out );}
LFHTEMP template<class I> auto Tuple<C,0u,Cflag>::mkInnerProd(const I& other, ITERABLE_DEF(I)) const
 -> decltype((*this)[0] * other[0]){
    decltype((*this)[0] * other[0]) fout;
    uint32_t ite;
    uint32_t ind;
    bool first = true;
    if (auto ite = other.getIterator) do{
        if (ite() < this->getSize()){
            if (first) fout = (*this)[ite()] * (*ite);
            else fout += (*this)[ite()] * (*ite);
            first = false;
        }
    }while(ite++);
    if (first) ExOp::toZero(fout);
return fout;}
LFHTEMP template<class O> auto Tuple<C,0u,Cflag>::mkInnerProd(const SparseTuple<O>& other) const
 -> decltype((*this)[0] * other[0]){
    decltype((*this)[0] * other[0]) fout;
    uint32_t ite;
    uint32_t ind;
    bool first = true;
    for(ite=0;ite<other.getSize();ite++){
        if (other.deref_key(ite) < this->getSize()){
            if (first) fout = (*this)[other.deref_key(ite)] * other.deref(ite);
            else fout += (*this)[other.deref_key(ite)] * other.deref(ite);
            first = false;
        }
    }
    if (first) ExOp::toZero(fout);
return fout;}
LFHTEMP void Tuple<C,0u,Cflag>::show(char* &f_out, int l) const{
	switch(l){
    case 1:
    case 0:
        for(unsigned int i=0;i<tup_size;i++) {
            if (i != 0) *(f_out++) = '\t';
            ExOp::show(data[i],f_out,2);
        }
        if (l==0) *(f_out++) = '\n';
    break;
    case 2:
        *(f_out++) = '[';
        for(unsigned int i=0;i<tup_size;i++) {
            if (i != 0) *(f_out++) = ';';
            ExOp::show(data[i],f_out,3);
        }
        *(f_out++) = ']';
    break;
    default:
        *(f_out++) = '(';
        for(unsigned int i=0;i<tup_size;i++) {
            if (i != 0) *(f_out++) = ',';
            ExOp::show(data[i],f_out,4);
        }
        *(f_out++) = ')';
    break;
	}
}
LFHTEMP void Tuple<C,0u,Cflag>::show( FILE *f_out, int l) const{
	switch(l){
		case 1:
		case 0:
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(data[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
		break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(data[i],f_out,3);
			}
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(data[i],f_out,4);
			}
			fprintf(f_out,")");
			break;
	}
}
LFHTEMP Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::toMemmove(Vector<C>& source){
    this->setSize(source.getSize());
    unsigned int i;
    for(i=0;i<source.getSize();i++) ExOp::toMemmove(data[i], source[i]);
    source.toMemfree();
return *this;}
LFHTEMP Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::operator=(const Vector<C>& source){
    this->setSize(source.getSize());
    unsigned int i;
    for(i=0;i<source.getSize();i++) data[i] = source[i];
return *this;}
LFHTEMP template<class O, unsigned int S, Tuple_flag CFLAG2> auto Tuple<C,0u,Cflag>::mkInnerProd(const Tuple<O,S,CFLAG2>&other)const
 -> decltype((*this)[0] * other[0]){
    decltype((*this)[0] * other[0]) fout;
    if ((this->getSize() | other.getSize()) == 0) ExOp::toZero(fout);
    else{
        fout = data[0] * other.data[0];
        uint32_t i = (this->getSize() < other.getSize()) ? this->getSize() -1: other.getSize()-1;
        if (this->getSize() != other.getSize()) printf("warning, inner product has sizes %i and %i\n", this->getSize(), other.getSize() );
        for(;i>0;i--) fout += data[i] * other.data[i];
    }
return fout;}
LFHTEMP template<class O, unsigned int S, Tuple_flag CFLAG2> TMatrix<C,0u,0u> Tuple<C,0u,Cflag>::mkOuterProd(const Tuple<O,S,CFLAG2>&other)const{TMatrix<C,0u,0u> fout;

return fout;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D, uint32_t L>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::toMult(const TMatrix<O,0u,0u> & a, const Tuple<D,L> & b) const{
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D, uint32_t L>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::toMult(const Tuple<D,L> & a, const TMatrix<O,0u,0u> & b) const{
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::toMult(const TMatrix<O,0u,0u> & a, const SparseTuple<D> & b) const{
    uint32_t ite,y;
    C* tmp = a.data;
    this->setSize(a.sizes[1]);
    if (b.getSize() == 0) {this->toZero();return *this;}
    for(y=0;y<a.sizes[1];y++){
        data[y] = tmp[b.deref_key(0)] * b.deref(0);
        for(ite=1;ite <b.getSize();ite++) data[y] += tmp[b.deref_key(ite)] * b.deref(ite);
        tmp += a.sizes[0];
    }
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::toMult(const SparseTuple<D> & a, const TMatrix<O,0u,0u> & b) const{
    uint32_t ite,y;
    C* tmp = a.data;
    this->setSize(b.sizes[0]);
    if (a.getSize() == 0) {this->toZero();return *this;}
    tmp = b.data + (tmp[a.deref_key(ite)] * b.sizes[0]);
    for(y=0;y<b.sizes[1];y++) data[y] = *(tmp++) * a.deref(0);

    for(ite=1;ite <b.getSize();ite++){
        tmp = a.data + (tmp[a.deref_key(ite)] * b.sizes[0]);
        for(y=0;y<b.sizes[1];y++) data[y] += *(tmp++) * a.deref(ite);
    }
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::addMult(const TMatrix<O,0u,0u> & a, const SparseTuple<D> & b) const{
    uint32_t ite,y;
    C* tmp = a.data;
    this->setSize(a.sizes[1]);
    if (b.getSize() == 0) {this->toZero();return *this;}
    for(y=0;y<a.sizes[1];y++){
        for(ite=0;ite <b.getSize();ite++) data[y] += tmp[b.deref_key(ite)] * b.deref(ite);
        tmp += a.sizes[0];
    }
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::addMult(const SparseTuple<D> & a, const TMatrix<O,0u,0u> & b) const{
    uint32_t ite,y;
    C* tmp = a.data;
    this->setSize(b.sizes[0]);
    if (a.getSize() == 0) {this->toZero();return *this;}
    for(ite=0;ite <b.getSize();ite++){
        tmp = a.data + (tmp[a.deref_key(ite)] * b.sizes[0]);
        for(y=0;y<b.sizes[1];y++) data[y] += *(tmp++) * a.deref(ite);
    }
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::subtMult(const TMatrix<O,0u,0u> & a, const SparseTuple<D> & b) const{
    uint32_t ite,y;
    C* tmp = a.data;
    this->setSize(a.sizes[1]);
    if (b.getSize() == 0) {this->toZero();return *this;}
    for(y=0;y<a.sizes[1];y++){
        for(ite=0;ite <b.getSize();ite++) data[y] -= tmp[b.deref_key(ite)] * b.deref(ite);
        tmp += a.sizes[0];
    }
return *this;}
LFHTEMP template<class O, uint32_t Tuple_flag, class D>
Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::subtMult(const SparseTuple<D> & a, const TMatrix<O,0u,0u> & b) const{
    uint32_t ite,y;
    C* tmp = a.data;
    this->setSize(b.sizes[0]);
    if (a.getSize() == 0) {this->toZero();return *this;}
    for(ite=0;ite <b.getSize();ite++){
        tmp = a.data + (tmp[a.deref_key(ite)] * b.sizes[0]);
        for(y=0;y<b.sizes[1];y++) data[y] -= *(tmp++) * a.deref(ite);
    }
return *this;}
LFHTEMP Tuple<C,0u, Cflag>& Tuple<C,0u,Cflag>::permute(const Tuple<uint32_t>& permutation){
    Tuple<C,0u,Cflag> fout; fout.setSize(permutation.getSize());
    int i;
    if (permutation.getSize() > this->getSize()) return *this;
    for(i=0;i< permutation.getSize();i++) fout.data[permutation[i]] = std::move(data[i]);
    for(i=0;i< permutation.getSize();i++) data[i] = std::move(fout.data[i]);
return *this;}
LFHTEMP double Tuple<C,0u,Cflag>::pnorm() const{
	if (tup_size){
	double sum = ExOp::pnorm(data[0]);
	for(unsigned int i =1;i<tup_size;i++) sum += ExOp::pnorm(data[i]);
	return(sum);
	}else return 0;
}

LFHTEMP double Tuple<C,0u,Cflag>::norm() const{
	if (tup_size){
	double sum = ExOp::pnorm(data[0]);
	for(unsigned int i =1;i<tup_size;i++) sum += ExOp::pnorm(data[i]);
	return(sqrt(sum));
	}else return 0;
}


LFHTEMP template<class O,unsigned int oTSIZE, Tuple_flag Oflag> char Tuple<C,0u,Cflag>::compare(const Tuple<O,oTSIZE, Oflag>& other) const{
    unsigned int mins = tup_size > other.getSize() ? other.getSize() : tup_size;

    unsigned int i;
    for(i=0;i< mins;i++) if (ExOp::isNQ((*this)[i], other[i])) break;
    if (i < mins) return (ExOp::isLT((*this)[i], other[i])) ? 2 : 1;
    return (tup_size == other.getSize()) ? 0 : (tup_size == mins) ? 2 : 1;
}

LFHTEMP template<class O> std::vector<O>& Tuple<C,0u,Cflag>::wrStdVector(std::vector<O>& fout)const{
    fout.clear();
    fout.reserve(getSize());
    for(uint32_t i =0;i<getSize();i++) fout.push_back(data[i]);
return fout;}

LFHTEMP ERRCODE Tuple<C,0u,Cflag>::save( FILE *f) const {
	if (fwrite(&tup_size,sizeof(unsigned int),1,f) != 1) return 1;
	if (tup_size == 0) return 0;
	if (ExCo<C>::IsPOD) return (fwrite(data,sizeof(C),tup_size,f) == tup_size) ? 0 : 1;
	ERRCODE fout = ExOp::save(data[0],f);
	for(unsigned int i=1;i<tup_size;i++) fout |= ExOp::save(data[i],f);
	return fout;
}
LFHTEMP ERRCODE Tuple<C,0u,Cflag>::load( FILE *f, unsigned int ch_TSIZE) { ERRCODE fout =0;
	unsigned int i;
	fout = (fread(&i,sizeof(unsigned int),1,f) == 1) ? 0 : 1;
	this->setSize(i);
	if (ExCo<C>::IsPOD) {if (tup_size) fout |= fread(data,sizeof(C),tup_size,f);
	}else for(unsigned int i=0;i<tup_size;i++) fout |= ExOp::load(data[i],f, ch_TSIZE/tup_size);
	return fout;
}



LFHTEMP string Tuple<C,0u,Cflag>::type_tostring() const{return string("Tuple<") + ExOp::type_tostring(data[0]) +string(",0u>");}

#undef LFHTEMP
#define LFHTEMP template<class C>


LFHTEMP Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::operator+=(LFHPrimitive::Iterator<uint32_t, const C> other){
    bool loop = other;
    for(uint32_t i = 0; (i< tup_size)&&(loop);i++, loop= other++) data[i] += *other;
return *this;}
LFHTEMP Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::operator-=(LFHPrimitive::Iterator<uint32_t, const C> other){
    bool loop = other;
    for(uint32_t i = 0; (i< tup_size)&&(loop);i++, loop= other++) data[i] -= *other;
return *this;}
LFHTEMP Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::operator*=(LFHPrimitive::Iterator<uint32_t, const C> other){
    bool loop = other;
    for(uint32_t i = 0; (i< tup_size)&&(loop);i++, loop= other++) data[i] *= *other;
return *this;}
LFHTEMP Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::operator/=(LFHPrimitive::Iterator<uint32_t, const C> other){
    bool loop = other;
    for(uint32_t i = 0; (i< tup_size)&&(loop);i++, loop= other++) data[i] /= *other;
return *this;}


LFHTEMP void Tuple<C,0u,TUPLE_FLAG_REMOTE_MEMORY>::show(FILE* f_out, int l)const{
    switch(l){
		case 1:
		case 0:
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(data[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
		break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(data[i],f_out,3);
			}
			fprintf(f_out,"]");
        break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(data[i],f_out,4);
			}
			fprintf(f_out,")");
        break;
	}
}

#undef LFHTEMP
#define LFHTEMP template<Tuple_flag Cflag>

LFHTEMP Tuple<Anything,0u,Cflag>& Tuple<Anything,0u,Cflag>::setSize(uint32_t nsize){

return *this;}
LFHTEMP ERRCODE Tuple<Anything,0u,Cflag>::save( FILE *f) const {
	if (fwrite(&tup_size,sizeof(uint32_t),1,f) != 1) return 1;
	if (fwrite(&anytype,sizeof(uint32_t),1,f) != 1) return 1;
	if (tup_size == 0) return 0;
	return (fwrite(data,anytype & 0xFF,tup_size,f) == tup_size) ? 0 : 1;
}
LFHTEMP ERRCODE Tuple<Anything,0u,Cflag>::load( FILE *f, unsigned int ch_TSIZE) { ERRCODE fout =0;
	unsigned int i;
	fout = (fread(&i,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
	fout |= (fread(&anytype,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
	this->setSize(i);
	if (tup_size) fout |= (fread(data,anytype & 0xFF,tup_size,f) != tup_size ? 1 : 0);
	return fout;
}

#undef LFHTEMP
#define LFHTEMP template<class C,unsigned int order>

	LFHTEMP void WeightElem<C,order>::zero(){for(unsigned int i=0;i<order;i++) {w[i] =0.0f; ExOp::toZero(e[i]);} }
	LFHTEMP void WeightElem<C,order>::random(){for(unsigned int i=0;i<order;i++) ExOp::toRand(e[i]); for(unsigned int i=0;i<order;i++) w[i] =1.0f;}
	LFHTEMP void WeightElem<C,order>::one(){ExOp::toOne(e); for(unsigned int i=0;i<order;i++) {w[i] =1.0f;} }
	LFHTEMP ERRCODE WeightElem<C,order>::save( FILE *f) const {ERRCODE fout = ExOp::save(e,f); return fout | ExOp::save(w,f);}
	LFHTEMP ERRCODE WeightElem<C,order>::load( FILE *f, unsigned int TSIZE) {return ExOp::load(e,f) | ExOp::load(w,f);}
	LFHTEMP void WeightElem<C,order>::show( FILE *f_out, int l) const{
		switch(l){
			case 0:
				fprintf(f_out,"average:"); ExOp::show(getMean(), f_out,1);
				if (order >1) {fprintf(f_out,"\nvariance:"); ExOp::show(getVar(), f_out,1);}
				if (order >2) {fprintf(f_out,"\nskew:"); ExOp::show(getSkew(), f_out,1);}
				if (order >3) {fprintf(f_out,"\nkurt:"); ExOp::show(getKurt(), f_out,1);}
				fprintf(f_out,"\n");
			break;
			case 1:
				fprintf(f_out,"average:"); ExOp::show(getMean(), f_out,2);
				if (order >1) {fprintf(f_out,"\tvariance:"); ExOp::show(getVar(), f_out,2);}
				if (order >2) {fprintf(f_out,"\tskew:"); ExOp::show(getSkew(), f_out,2);}
				if (order >3) {fprintf(f_out,"\tkurt:"); ExOp::show(getKurt(), f_out,2);}
			break;
		//	case 2:
		//		fprintf(f_out,"[");
		//		for(int i=0;i<TSIZE;i++) {
		//			if (i != 0) fprintf(f_out,";");
		//			ExOp::show(data[i],f_out,3);
		//		}
		//		fprintf(f_out,"]");
		//		break;
		//	default:
		//		fprintf(f_out,"(");
		//		for(int i=0;i<TSIZE;i++) {
		//			if (i != 0) fprintf(f_out,",");
		//			ExOp::show(data[i],f_out,4);
		//		}
		//		fprintf(f_out,")");
		}
	}
LFHTEMP string WeightElem<C,order>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",order); return string("WeightElem<") + ExOp::type_tostring(e[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}

LFHTEMP void WeightElem<C, order>::replaceNan(C mean, C var, bool force_var_minimum){
    C tmp;
    if (!ExOp::isValid(tmp = this->getVar())){
        if (!ExOp::isValid(tmp = this->getMean())) e[0] = mean;
        w[0] = 1.0f;
        w[1] = 0.0f;
        e[1] = var + e[0]*e[0];
    }else if ((force_var_minimum)&&(tmp < var)){
        e[1] = var * (w[0] - (w[1] / w[0])) + (e[0] * e[0] / w[0]);
    }
}

	LFHTEMP WeightElem<C,order> WeightElem<C,order>::operator-(const WeightElem<C,order>& s) const{ WeightElem<C,order> f_out;
		f_out.w = w + s.w;
		for(int i=0;i< order;i++) f_out.e[i] = (i & 1) ? e[i] + s.e[i] : e[i] - s.e[i];
		return(f_out);
	}




	LFHTEMP void WeightElem<C, order>::OpMean::operator()(C & out, WeightElem<C, order> & in) const {out = in.getMean();}



	LFHTEMP void WeightElem<C, order>::OpVar::operator()(C & out, WeightElem<C, 2> & in) const {
		LinkAssert<(order>1) > ass;
		out = in.getVar();
//	printf("%f,%f\n");
	}


	LFHTEMP void WeightElem<C, order>::OpStd::operator()(C & out, WeightElem<C, 2> & in) const { out = ExCo<C>::invintPow(in.getVar(),2);}




#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int sizex, unsigned int sizey, Tuple_flag TF>
// TMatrix::TMatrix

LFHTEMP TMatrix<C,sizex,sizey,TF>::TMatrix(TMatrix<C, sizex, sizey, TF> const & idata){
    for(unsigned int k=0;k<sizex*sizey;k++) data[k] = idata.data[k];
}

	LFHTEMP TMatrix<C,sizex,sizey,TF>::TMatrix(DataGrid<C, 2> const & idata){
		Tuple<unsigned int,2> coor;
		for(coor[1] = 0; coor[1] < sizey;coor[1]++)
			for(coor[0] = 0; coor[0] < sizex;coor[0]++) data[coor[0] + coor[1]*sizex] = idata(coor);
	}

	LFHTEMP
	TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::operator=(TMatrix<C, sizex, sizey,TF> const & idata){
		for(unsigned int k=0;k<sizex*sizey;k++) data[k] = idata.data[k];
		return(*this);
	}

	LFHTEMP
	void TMatrix<C,sizex,sizey,TF>::operator()(Tuple<C, sizey> &_out, Tuple<C, sizex> &_in) const{
		_out = (*this) * _in;
	}

	LFHTEMP
	void TMatrix<C,sizex,sizey,TF>::derivative(Tuple<C, sizey> &_out, Tuple<C, sizex > &_in, int in_direct) const{
		int i;
		for(i=0;i<sizey;i++) _out[i] = data[in_direct + i * sizex];

	}

LFHTEMP void TMatrix<C,sizex,sizey,TF>::derivativeMatrix(TMatrix<C, sizex, sizey,TF> & _out, Tuple<C, sizex > &_in) const{_out = (*this);}
LFHTEMP void TMatrix<C,sizex,sizey,TF>::show(FILE* o, int level) const{
    unsigned int i,j;
    for(j=0;j<sizey;j++){
        for(i=0;i<sizex;i++){
            ExOp::show(data[i + j*sizex],o,1);
            fprintf(o,"%c", i == sizex-1 ? '\n' : '\t');
            }
        }

}
LFHTEMP TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,sizex,sizey,TF>::mkInverse() const{
    TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> tmp = this->mkTranspose();
    if (sizex == sizey){
        C det[8];
        if (sizex == 2){
            det[0] = data[0] * data[3] - data[1] * data[2];
            tmp.data[0] = data[3] / det[0];
            tmp.data[1] = -data[1] / det[0];
            tmp.data[2] = -data[2] / det[0];
            tmp.data[3] = data[0] / det[0];
            return tmp;
        }else if (sizex == 3){
            det[1] = data[0] * data[4] - data[1] * data[3];
            det[2] = data[1] * data[5] - data[2] * data[4];
            det[3] = data[2] * data[3] - data[0] * data[5];
            det[0] = data[6] * det[2] + data[7] * det[3] + data[8] * det[1];
            tmp.data[8] = det[1] / det[0];
            tmp.data[2] = det[2] / det[0];
            tmp.data[5] = det[3] / det[0];
            tmp.data[0] = (data[4]*data[8]-data[5]*data[7]) / det[0];
            tmp.data[3] = (data[5]*data[6]-data[3]*data[8]) / det[0];
            tmp.data[6] = (data[3]*data[7]-data[4]*data[6]) / det[0];
            tmp.data[1] = -(data[1]*data[8]-data[2]*data[7]) / det[0];
            tmp.data[4] = -(data[2]*data[6]-data[0]*data[8]) / det[0];
            tmp.data[7] = -(data[0]*data[7]-data[1]*data[6]) / det[0];
            return tmp;
        }else if (sizex == 4){
            det[0] = data[10] * data[15] - data[11] * data[14]; // 00XX
            det[1] = data[11] * data[12] - data[8] * data[15]; // X00X
            det[2] = data[8] * data[13] - data[9] * data[12]; // XX00
            det[3] = data[9] * data[14] - data[10] * data[13]; // 0XX0
            det[4] = data[9] * data[15] - data[11] * data[13]; // 0X0X
            det[5] = data[8] * data[14] - data[10] * data[12]; // X0X0
            tmp.data[0] = data[5] * det[0] - data[6] * det[4] + data[7] * det[3];
            tmp.data[4] = -(data[4] * det[0] + data[6] * det[1] + data[7] * det[5]);
            tmp.data[8] = data[4] * det[4] + data[5] * det[1] + data[7] * det[2];
            tmp.data[12] = -(data[4] * det[3] - data[5] * det[5] + data[6] * det[2]);
            tmp.data[1] = -(data[1] * det[0] - data[2] * det[4] + data[3] * det[3]);
            tmp.data[5] = data[0] * det[0] + data[2] * det[1] + data[3] * det[5];
            tmp.data[9] = -(data[0] * det[4] + data[1] * det[1] + data[3] * det[2]);
            tmp.data[13] = data[0] * det[3] - data[1] * det[5] + data[2] * det[2];
            det[6] = data[1] * data[7] - data[3] * data[5]; det[7] = -det[6] * det[5]; // 0X0X
            det[5] = data[0] * data[6] - data[2] * data[4]; det[7] -= det[5] * det[4]; // X0X0
            det[4] = data[3] * data[4] - data[0] * data[7]; det[7] -= det[4] * det[3]; // X00X
            det[3] = data[2] * data[7] - data[3] * data[6]; det[7] += det[3] * det[2]; // 00XX
            det[2] = data[1] * data[6] - data[2] * data[5]; det[7] -= det[2] * det[1]; // 0XX0
            det[1] = data[0] * data[5] - data[1] * data[4]; det[7] += det[1] * det[0]; // XX00
            tmp.data[2] = data[13] * det[3] - data[14] * det[6] + data[15] * det[2];
            tmp.data[6] = -(data[12] * det[3] + data[14] * det[4] + data[15] * det[5]);
            tmp.data[10] = data[12] * det[6] + data[13] * det[4] + data[15] * det[1];
            tmp.data[14] = -(data[12] * det[2] - data[13] * det[5] + data[14] * det[1]);
            tmp.data[3] = -(data[9] * det[3] - data[10] * det[6] + data[11] * det[2]);
            tmp.data[7] = data[8] * det[3] + data[10] * det[4] + data[11] * det[5];
            tmp.data[11] = -(data[8] * det[6] + data[9] * det[4] + data[11] * det[1]);
            tmp.data[15] = data[8] * det[2] - data[9] * det[5] + data[10] * det[1];
            for(uint32_t i=0;i<16;i++) tmp.data[i] /= det[7];
            return tmp;
        }


    }


    C* buf;
    double* normbuf;
    unsigned int i,z;

    if (sizex >= sizey){
        buf = new C[sizey * sizex]; LFH_NICE_ALLOCERROR(buf,"")
        normbuf = new double[(sizey << 1) | 1]; LFH_NICE_ALLOCERROR(normbuf,"")
        for(i=0;;i++){
            double sum2 =0.0f;
            for(z=1;z<sizex-i;z++) {
                buf[z] = tmp.data[i+z * sizey];
                sum2 += ExOp::pnorm(buf[z]);
            }
            buf[0] = tmp.data[i * (sizey +1)];
            sum2 += ExOp::pnorm(buf[0]);
            double tmpsum = sqrt(sum2);
            buf[0] += ((buf[0] < 0) ? -tmpsum: tmpsum);
            normbuf[(i << 1)] = (sum2 + fabs(tmp.data[i * (sizey +1)] *tmpsum));
            tmp.leftHouseHolderMultiply(buf, sizex - i, normbuf[(i << 1)], true);
            break;
            buf += sizex - i;
            if (i == sizey -1) break;


            tmp.rightHouseHolderMultiply(buf, sizex - i, normbuf[(i << 1)| 1],true);
        };
        tmp.show();

    }else{
        exit(1); // TODO!

    }
    delete[](buf);
    delete[](normbuf);
    /*

    DataGrid<double, 2> dainv;
    Vector<double> eigen;
    dainv.makeDiagonalizer_ofinverse(*this,eigen);
    TMatrix<double, sizex, sizey> fun(dainv);
    TMatrix<double, sizex, sizey> funt = fun.transpose();
    for(unsigned int i=0;i<sizex;i++) for(unsigned int j=0;j<sizey;j++) fun.data[i+j*sizex] *= eigen[i];
    fun *= funt;
    return(fun);
    */
    return tmp;
}
LFHTEMP TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,sizex,sizey,TF>::mkTranspose() const{ TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> f_out;
    for(unsigned int i= 0 ; i< sizex;i++)
        for(unsigned int j= 0 ; j< sizey;j++) f_out.data[j+i*sizey] = data[i+j*sizex];
    return(f_out);
}
LFHTEMP TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,sizex,sizey,TF>::mkTrju() const{ TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> f_out;
    for(unsigned int i= 0 ; i< sizex;i++)
        for(unsigned int j= 0 ; j< sizey;j++) f_out.data[j+i*sizey] = ExOp::mkTrju(data[i+j*sizex]);
return(f_out);}

	/*LFHTEMP
	C TMatrix<C,sizex,sizey,TF>::determinant() const{
		DataGrid<double, 2> dainv;
		Vector<double> eigen;
		dainv.makeDiagonalizer_ofinverse(*this,eigen);
		C fout = (C) eigen[0];
		for(unsigned int i=1;i<sizex;i++) fout *= eigen[i];
		return(fout);
	}*/

	LFHTEMP TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::toZero(){
		unsigned int i = sizex*sizey;
		for(i--;i != ExCo<unsigned int>::mkMaximum() ; i--) ExOp::toZero(data[i]);
		return(*this);
	}
	LFHTEMP TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::toOne(){
		unsigned int x,y,z;
		for(y=0,z=0;y<sizey;y++)
			for(x=0;x<sizex;x++,z++)
				if (x == y) ExOp::toOne(data[z]);
				else ExOp::toZero(data[z]);
		return(*this);
	}
	LFHTEMP TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::toRand(){
		unsigned int i = sizex*sizey;
		for(i--;i != ExCo<unsigned int>::mkMaximum() ; i--) ExOp::toRand(data[i]);
		return(*this);
	}
	LFHTEMP TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::toRandSymetric(){
		unsigned int x,y,z;
		for(y=0,z=0;y<sizey;y++)
			for(x=0;x<sizex;x++,z++) if (y > x) data[z] = data[y + x * sizex]; else ExOp::toRand(data[z]);
		return(*this);
	}

	LFHTEMP
	template<class O,Tuple_flag FO>
	TMatrix<C,sizex,sizey,TF> TMatrix<C,sizex,sizey,TF>::scale_rows(Tuple<O, sizey, FO> const &fact) const{
		TMatrix<C,sizex,sizey,TF> out;

		unsigned int i,j;
		for(i=0;i<sizex;i++) for(j=0;j<sizey;j++) out.data[i + j * sizex] = data[i + j * sizex] * fact[j];
		return(out);
	}

	LFHTEMP
	template<class O,Tuple_flag FO>
	TMatrix<C,sizex,sizey,TF> TMatrix<C,sizex,sizey,TF>::scale_cols(Tuple<O, sizex, FO> const &fact) const{
		TMatrix<C,sizex,sizey,TF> out;

		unsigned int i,j;
		for(i=0;i<sizex;i++) for(j=0;j<sizey;j++) out.data[i + j * sizex] = data[i + j * sizex] * fact[i];
		return(out);
	}

	LFHTEMP Trianglix<C,sizey> TMatrix<C,sizex,sizey,TF>::operator*(const Trianglix<C,sizex> &input)const{Trianglix<C,sizey> fout;
		// Tr_OUT = (this) * tr_IN * (this)^T
		TMatrix<C,sizex,sizey,TF> proj;
		unsigned int x,y,k;
		for(y=0;y< sizey;y++)
			for(x=0;x< sizex;x++){
				proj.data[x + y * sizex] = ExOp::mkMult(this->data[x + y * sizex],input.cell(x,x));
				for(k=0;k<x;k++) proj.data[x + y * sizex] += ExOp::mkMult(this->data[k + y * sizex],input.cell(k,x));
				for(k++;k<sizex;k++) proj.data[x + y * sizex] += ExOp::mkMult(this->data[k + y * sizex],ExOp::mkTrju(input.cell(k,x)));
			}

		proj.show();

		for(y=0;y< sizey;y++){
			for(x=0;x<= y;x++){
				fout.cell(y,x) = proj.data[ y* sizex] * ExOp::mkTrju(this->data[x * sizex]);
				for(k=1;k<sizex;k++) fout.cell(y,x) += proj.data[k + y* sizex] * ExOp::mkTrju(this->data[k + x * sizex]);
			}
		}
		return fout;
	}

    LFHTEMP
    TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,sizex,sizey,TF>::inverse_2() const{ TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> fout = this->mkTranspose();
        bool dir = (sizex < sizey);
        TMatrix<C,sizex,sizex> HH_X; HH_X.toOne();
        TMatrix<C,sizex,sizex> HH_Y; HH_Y.toOne();
        unsigned int x=0;
        unsigned int y=0;
        unsigned int z;
     //   if (dir)

     //   }
        double tmp;
        C* hh = new C[ (sizex < sizey) ? sizey:sizex]; LFH_NICE_ALLOCERROR(hh,"")
            double sum2 =0.0f;
			for(z=x+1;z<sizex;z++) {
				hh[z] = fout.data[y+z * sizey];
				sum2 += ExOp::pnorm(hh[z]);
			}

			hh[x] = fout.data[y+x * sizey];
			if ((sum2 != 0.0f)||(fabs(ExOp::pnorm(hh[x]) / sum2) > 1.0f / ExCo<double>::epsilon())) {
			sum2 += ExOp::pnorm(hh[x]);
			tmp = sqrt(sum2);
			hh[x] += ((hh[x] < 0) ? -tmp: tmp);
			tmp = (sum2 + fabs(fout.data[y+x * sizey] *tmp));

            }
            y++;
            dir = !dir;
        delete[](hh);

   //     HH_X.show();
    //    HH_Y.show();

return fout;}
LFHTEMP	const TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::LeftHouseHolderMultiply(const double * const vec, const double& sqrt_den, const int lenght, C* inner, bool hint){
			int inl = (hint) ? sizey + lenght - sizex : sizex; // hint tells that the first columns are immune to transformation

			if ((sqrt_den == 0.0f)||(!ExCo<double>::isValid(sqrt_den))) return *this;

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;

			for(j=0;j< inl;j++) inner[j] = data[ j + i* sizex] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += data[ j + i * sizex] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) data[j + i * sizex] += inner[j] * (vec[i] * fact);
			}

            return(*this);
		}

LFHTEMP void TMatrix<C,sizex,sizey,TF>::leftHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){
	int inl = (hint) ? sizey + lenght - sizex : sizex; // hint tells that the first columns are immune to transformation

	if ((sqrt_den == 0.0f)||(!ExCo<double>::isValid(sqrt_den))) return;

	C* min = data + sizex - inl +(sizey - lenght) * sizex;

	double fact = -1.0f / sqrt_den;

	int i,j;
	i=0;
	C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner,"")

	for(j=0;j< inl;j++) inner[j] = min[ j + i* sizex] * vec[i];

	for(i++;i<lenght;i++){
		for(j=0;j< inl;j++) inner[j] += min[ j + i * sizex] * vec[i];
	}
	for(i=0;i<lenght;i++){
		for(j=0;j< inl;j++) min[j + i * sizex] += inner[j] * (vec[i] * fact);
	}

	delete[](inner);
}

		// hint implies that the outside the region of interes are zeroes only!
LFHTEMP void TMatrix<C,sizex,sizey,TF>::rightHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){
	if ((sqrt_den == 0.0f)||(!ExCo<double>::isValid(sqrt_den))) return;

	int inl = (hint) ? sizex + lenght - sizey : sizey;
	int minj = sizey - inl;

	C* min = data + sizex - lenght + minj * sizex;

	double fact = -1.0f / sqrt_den;

	int i,j;
	i=0;
	C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner,"")


	for(j=0;j< inl;j++) inner[j] = min[ i + j * sizex] * vec[i];

	for(i++;i<lenght;i++){
		for(j=0;j< inl;j++) inner[j] += min[ i +  j * sizex] * vec[i];
	}
	for(i=0;i<lenght;i++){
		for(j=0;j< inl;j++) min[i + j * sizex] += inner[j] * (vec[i] * fact);
	}

	delete[](inner);
}

/*

	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey,TF>::operator+(A const & other) const{
		TMatrix<O,sizex,sizey> _out; int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) _out.data[i] = data[i] + other;
		return(_out);
	}
	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey,TF>::operator-(A const & other) const{
		TMatrix<O,sizex,sizey> _out;int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) _out.data[i] = data[i] - other;
		return(_out);
	}
	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey,TF>::operator*(A const & other) const{
		TMatrix<O,sizex,sizey> _out;int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) _out.data[i] = data[i] * other;
		return(_out);
	}

	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey,TF>::operator/(A const & other) const{
		TMatrix<O,sizex,sizey> _out;int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) _out.data[i] = data[i] / other;
		return(_out);
	}*/




LFHTEMP	template<class O, Tuple_flag OTF> TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::operator*=(TMatrix<O,sizex,sizex,OTF> const & other){
    C tmpr[sizex];
    unsigned int i,j,k;

    for(j = 0;j<sizey;j++){
    for(i=0;i<sizex;i++) {
        tmpr[i] = data[i + j*sizex];
        data[i + j*sizex] =tmpr[0] * other.data[i];
    }
    for(k=1;k<sizex;k++){
    for(i=0;i<sizex;i++) data[i + j*sizex] +=tmpr[k] * other.data[i+ k*sizex];
    }
    }
return(*this);}
LFHTEMP template<class O, Tuple_flag OTF> TMatrix<C,sizex,sizey,TF>& TMatrix<C,sizex,sizey,TF>::operator/=(TMatrix<O,sizex,sizex,OTF> const & other) {
    printf("no available...!\n"); exit(1);
    return(*this);
}

LFHTEMP template<class O, unsigned int osize, Tuple_flag OTF>
TMatrix<C,osize,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,sizex,sizey,TF>::operator*(TMatrix<O,osize,sizex,OTF> const & other) const{
    TMatrix<C,osize,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> _out;
    unsigned int i,j,k;
    for(k=0;k<osize;k++){
        for(j=0;j<sizey;j++){
            _out.data[k + j * osize] = data[j * sizex] * other.data[k];
            for(i=1;i<sizex;i++){
                _out.data[k + j * osize] += data[i + j * sizex] * other.data[k + i * osize];
            }
        }
    }
    return(_out);
}
LFHTEMP template<class O> auto TMatrix<C,sizex,sizey,TF>::operator*(Tuple<O,sizex> const & other)const
 -> Tuple<decltype(data[0]*other[0]),sizey>{
    Tuple<decltype(data[0]*other[0]),sizey> _out;
    unsigned int i,j;
    for(j=0;j<sizey;j++){
        _out.data[j] = data[j * sizex] * other.data[0];
        for(i=1;i<sizex;i++){
            _out.data[j] += data[i + j * sizex] * other.data[i];
        }
    }
return(_out);}
LFHTEMP template<class O> auto TMatrix<C,sizex,sizey,TF>::mkBackMult(Tuple<O,sizey> const & other)const
 -> Tuple<decltype(data[0]*other[0]),sizex>{
    Tuple<decltype(data[0]*other[0]),sizex> _out;
    unsigned int i,j;
    for(j=0;j<sizex;j++){
        _out.data[j] = data[j] * other.data[0];
        for(i=1;i<sizey;i++){
            _out.data[j] += data[j + i * sizex] * other.data[i];
        }
    }
return(_out);}

	LFHTEMP
	Tuple<C,sizey> TMatrix<C,sizex,sizey,TF>::operator*(Tuple<C,sizex> const & other) const{
		Tuple<C,sizey> _out;
		unsigned int i,j;
		for(j=0;j<sizey;j++){
			_out.data[j] = data[j * sizex] * other.data[0];
			for(i=1;i<sizex;i++){
				_out.data[j] += data[i + j * sizex] * other.data[i];
			}
		}
		return(_out);
		}
	LFHTEMP
	Continuous<C,sizex,C,sizey>* TMatrix<C,sizex,sizey,TF>::derivative(Tuple<C,sizex> const & where) const{
		return(new TMatrix<C,sizex,sizey,TF>(*this));
		}


//	LFHTEMP
//	PolyThing<C,0> TMatrix<C,sizex,sizey>::makeCharacteristic() const{
		/*
		PolyThing<C,0> _out;
		int i;
		_out.setOrder(sizex);

		_out[0] = data[0];
		for(i=1;i<sizex;i++) _out[0] *= data[i * (sizex+1)];
	*/
//	}

#undef LFHTEMP
#define LFHTEMP template<class C, Tuple_flag TF>
// TMatrix::TMatrix

LFHTEMP TMatrix<C,0u,0u,TF>::TMatrix(){sizes[0] =0;sizes[1] =0;}
LFHTEMP TMatrix<C,0u,0u,TF>::~TMatrix(){}
LFHTEMP TMatrix<C,0u,0u,TF>::TMatrix(const TMatrix<C,0u,0u,TF>& _clone): data(_clone.data){
    sizes[0] =_clone.sizes[0];
    sizes[1] =_clone.sizes[1];
}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::operator=(TMatrix<C,0u,0u,TF> const & other){
    this->setSizes(other.sizes[0],other.sizes[1]);
    uint32_t i = sizes[0]*sizes[1];
    if (i == 0) return *this;
    if (ExCo<C>::IsPOD) memcpy(data(),other.data(), sizeof(C) * i);
    else{
        for(i--;i != ExCo<uint32_t>::mkMaximum() ; i--) data[i] = other.data[i];
    }
return *this;}
LFHTEMP template<class O, Tuple_flag  OTF> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::operator=(TMatrix<O, 0u, 0u,OTF> const & other){
    this->setSizes(other.sizes[0],other.sizes[1]);
    uint32_t i = sizes[0]*sizes[1];
    for(i--;i != ExCo<uint32_t>::mkMaximum() ; i--) data[i] = other.data[i];
return *this;}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toZero(){
    uint32_t i = sizes[0]*sizes[1];
    for(i--;i != ExCo<uint32_t>::mkMaximum() ; i--) ExOp::toZero(data[i]);
    return(*this);
}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toOne(){
    unsigned int x,y,z;
    for(y=0,z=0;y<sizes[1];y++)
        for(x=0;x<sizes[0];x++,z++)
            if (x == y) ExOp::toOne(data[z]);
            else ExOp::toZero(data[z]);
    return(*this);
}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toRand(){
    uint32_t i = sizes[0]*sizes[1];
    for(i--;i != ExCo<uint32_t>::mkMaximum() ; i--) ExOp::toRand(data[i]);
    return(*this);
}
LFHTEMP C& TMatrix<C,0u,0u,TF>::operator()(const Tuple<unsigned int, 2u> &coor){return data[coor[0] + coor[1]* sizes[0]];}
LFHTEMP C TMatrix<C,0u,0u,TF>::operator()(const Tuple<unsigned int, 2u> &coor) const{return data[coor[0] + coor[1]* sizes[0]];}
LFHTEMP C& TMatrix<C,0u,0u,TF>::operator()(unsigned int x, unsigned int y){return data[x + y* sizes[0]];}
LFHTEMP C TMatrix<C,0u,0u,TF>::operator()(unsigned int x, unsigned int y) const{return data[x + y* sizes[0]];}

LFHTEMP template<class O,Tuple_flag OTF> TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,0u,0u,TF>::operator*(TMatrix<O, 0u, 0u,OTF> const & other)const{
    TMatrix<C,0u,0u,TF> fout; fout.setSizes(other.sizes[0],sizes[1]);
    unsigned int i,j,k;
    for(k=0;k<other.sizes[0];k++){
        for(j=0;j<sizes[1];j++){
            fout.data[k + j * other.sizes[0]] = data[j * sizes[0]] * other.data[k];
            for(i=1;i<sizes[0];i++){
                fout.data[k + j * other.sizes[0]] += data[i + j * sizes[0]] * other.data[k + i * other.sizes[0]];
            }
        }
    }
return(fout);}
LFHTEMP template<class O> TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,0u,0u,TF>::operator*(O const & other)const{
    TMatrix<C,0u,0u,TF> fout; fout.setSizes(sizes[0],sizes[1]);
    unsigned int i,j,k;
    for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) fout.data[i] = data[i] * other;
return fout;}
LFHTEMP template<class O,class P> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toAddOuterProd(const SparseTuple<O> & a, const SparseTuple<P> & b){
    uint32_t i,j;
    C* currow;
    for(i=0;i<a.getSize();i++){
        currow = data() + a.deref_key(i) * sizes[0];
        for(j=0;j<b.getSize();j++){
            currow[b.deref_key(j)] += a.deref(i) * ExOp::mkTrju(b.deref(j));
        }
    }
return *this;}
LFHTEMP template<class O,class P> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toAddBackOuterProd(const SparseTuple<O> & a, const SparseTuple<P> & b){
    uint32_t i,j;
    C* currow;
    for(i=0;i<a.getSize();i++){
        currow() = data() + a.deref_key(i) * sizes[0];
        for(j=0;j<b.getSize();j++){
            currow[b.deref_key(j)] += ExOp::mkMult(ExOp::mkTrju(b.deref(j)), a.deref(i));
        }
    }
return *this;}
LFHTEMP template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toAddOuterProd(const Tuple<O,S> & a, const SparseTuple<P> & b){
    uint32_t i,j;
    C* currow;
    if (a.getSize() != sizes[0]) printf("warning, incoherent sizes %i and %i\n", a.getSize(), sizes[0]);
    for(j=0;j<b.getSize();j++){
        currow = data() + b.deref_key(j) * sizes[0];
        if (b.deref_key(j)>= sizes[1]) printf("%i > %i !!!\n", b.deref_key(j), sizes[1]);
        for(i=0;i<a.getSize();i++){
            currow[i] += ExOp::mkMult(a[i], ExOp::mkTrju(b.deref(j)));
        }
    }
return *this;}
LFHTEMP template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toAddBackOuterProd(const Tuple<O,S> & a, const SparseTuple<P> & b){
    uint32_t i,j;
    C* currow;
    if (a.getSize() != sizes[1]) printf("warning, incoherent sizes %i and %i\n", a.getSize(), sizes[1]);
    for(i=0;i<a.getSize();i++){
        currow = data() + i * sizes[0];
        for(j=0;j<b.getSize();j++){
            if ((i==0)&&(b.deref_key(j)>= sizes[0])) printf("%i > %i !!!\n", b.deref_key(j), sizes[0]);
            currow[b.deref_key(j)] += ExOp::mkMult(ExOp::mkTrju(a[i]),b.deref(j));
        }
    }
return *this;}
LFHTEMP template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toAddOuterProd(const SparseTuple<O> & a, const Tuple<P,S> & b){
    uint32_t i,j;
    C* currow;
    if (b.getSize() != sizes[1]) printf("warning, incoherent sizes %i and %i\n", b.getSize(), sizes[1]);
    for(j=0;j<b.getSize();j++){
        currow() = data() + j * sizes[0];
        for(i=0;i<a.getSize();i++){
            if ((j==0)&&(a.deref_key(i)>= sizes[0])) printf("%i > %i !!!\n", a.deref_key(i), sizes[0]);
            currow[a.deref_key(i)] += ExOp::mkMult(a.deref(i),ExOp::mkTrju(b[j]));
        }
    }
return *this;}
LFHTEMP template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toAddBackOuterProd(const SparseTuple<O> & a, const Tuple<P,S> & b){
    uint32_t i,j;
    C* currow;
    if (b.getSize() != sizes[0]) printf("warning, incoherent sizes %i and %i\n", b.getSize(), sizes[0]);
    for(i=0;i<a.getSize();i++){
        currow = data() + a.deref_key(i) * sizes[0];
        if (a.deref_key(i)>= sizes[1]) printf("%i > %i !!!\n", a.deref_key(i), sizes[1]);
        for(j=0;j<b.getSize();j++){
            currow[j] += ExOp::mkMult(ExOp::mkTrju(a.deref(i)), b[j]);
        }
    }
return *this;}
LFHTEMP template<class O, Tuple_flag OFLAG> auto TMatrix<C,0u,0u,TF>::operator*(const Tuple<O, 0u, OFLAG> & a)const
 -> Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL>{
    Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL> fout;
    fout.setSize(sizes[0]);
    uint32_t i;
    const C* tmp;
    if (sizes[1] != a.getSize()) printf("warning: size mismatch! %i %i\n", sizes[1],  a.getSize());
    if (auto ite = a.getIterator()) {
        tmp = data() + (ite() * sizes[0]);
        for(i=0;i<sizes[0];i++) fout[i] = tmp[i] * (*ite);
        while(ite++){
            tmp = data() + (ite() * sizes[0]);
            for(i=0;i<sizes[0];i++) fout[i] += tmp[i] * (*ite);
        }
    }else fout.toZero();
return fout;}
LFHTEMP template<class O, Tuple_flag OFLAG> auto TMatrix<C,0u,0u,TF>::mkBackMult(const Tuple<O, 0u, OFLAG> & a)const
 -> Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL>{
    Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL> fout;
    fout.setSize(sizes[1]);
    uint32_t i;
    if (sizes[0] != a.getSize()) printf("warning: size mismatch! %i %i\n", sizes[0],  a.getSize());
    const C* tmp;
    if (auto ite = a.getIterator()) {
        tmp = data() + ite();
        for(i=0;i<sizes[1];i++) fout[i] = tmp[i * sizes[0]] * (*ite);
        while(ite++){
            tmp = data() + ite();
            for(i=0;i<sizes[1];i++) fout[i] += tmp[i * sizes[0]] * (*ite);
        }
    }else fout.toZero();
return fout;}
LFHTEMP template<class O> auto TMatrix<C,0u,0u,TF>::operator*(const SparseTuple<O> & a)const
 -> Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL>{
    Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL> fout;
    fout.setSize(sizes[0]);
    uint32_t i;
    const C* tmp;
    if (auto ite = a.getIterator()) {
        tmp = data() + ite() * sizes[0];
        for(i=0;i<sizes[0];i++) fout[i] = tmp[i] * (*ite);
        while(ite++){
            tmp = data() + ite() * sizes[0];
            for(i=0;i<sizes[0];i++) fout[i] += tmp[i] * (*ite);
        }
    }else fout.toZero();
return fout;}
LFHTEMP template<class O> auto TMatrix<C,0u,0u,TF>::mkBackMult(const SparseTuple<O> & a)const
 -> Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL>{
    Tuple<decltype((*this)(0,0)* a[0] ),0u,TUPLE_FLAG_NULL> fout;
    fout.setSize(sizes[1]);
    uint32_t i;
    const C* tmp;
    if (auto ite = a.getIterator()) {
        tmp = data() + ite();
        for(i=0;i<sizes[1];i++) fout[i] = tmp[i * sizes[0]] * (*ite);
        while(ite++){
            tmp = data() + ite();
            for(i=0;i<sizes[1];i++) fout[i] += tmp[i * sizes[0]] * (*ite);
        }
    }else fout.toZero();
return fout;}

LFHTEMP Tuple<C,0u> TMatrix<C,0u,0u,TF>::getColumn(uint32_t colID)const{  Tuple<C,0u> fout;fout.setSize(sizes[0]);
    const C* tmp = data() + (colID * sizes[0]);
    for(uint32_t i=0;i<sizes[0];i++) fout[i] = tmp[i];
return fout;}
LFHTEMP void TMatrix<C,0u,0u,TF>::getDims(Tuple<unsigned int, 2u> &o_dims) const{o_dims[0] = sizes[0];o_dims[1] = sizes[1];}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::setSizes(uint32_t x, uint32_t y){
    if (x * y == data.getSize()) return *this;
    if ((x != 0)&&((0xFFFFFFFF / x) < y)) {fprintf(stderr,"Allocated Dense Matrix is too big! %i x %i\n", x, y); exit(1);}
    data.setSize(x*y);
    sizes[0] = x;
    sizes[1] = y;
return *this;}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::toMemmove(TMatrix<C,0u,0u,TF>& other){
    data = other.data;
    sizes[0] = other.sizes[0];
    sizes[1] = other.sizes[1];
    ExOp::toZero(other.sizes);
return *this;}
LFHTEMP TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,0u,0u,TF>::mkTranspose() const{TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> fout;
    fout.setSizes(sizes[1],sizes[0]);
    for(unsigned int i= 0 ; i< sizes[0];i++)
        for(unsigned int j= 0 ; j< sizes[1];j++) fout.data[j+i*sizes[1]] = data[i+j*sizes[0]];
return(fout);}
LFHTEMP TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,0u,0u,TF>::mkTrju() const{TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> fout;
    fout.setSizes(sizes[1],sizes[0]);
    for(unsigned int i= 0 ; i< sizes[0];i++)
        for(unsigned int j= 0 ; j< sizes[1];j++) fout.data[j+i*sizes[1]] = ExOp::mkTrju(data[i+j*sizes[0]]);
return(fout);}
LFHTEMP TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,0u,0u,TF>::mkInverse() const{TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> fout;
    fout.setSizes(sizes[1],sizes[0]);
    if (sizes[0] == sizes[1]){
        C det[8];
        switch(sizes[0]){
        case 2:
            det[0] = data[0] * data[3] - data[1] * data[2];
            fout.data[0] = data[3] / det[0];
            fout.data[1] = -data[1] / det[0];
            fout.data[2] = -data[2] / det[0];
            fout.data[3] = data[0] / det[0];
        return fout;
        case 3:
            det[1] = data[0] * data[4] - data[1] * data[3];
            det[2] = data[1] * data[5] - data[2] * data[4];
            det[3] = data[2] * data[3] - data[0] * data[5];
            det[0] = data[6] * det[2] + data[7] * det[3] + data[8] * det[1];
            fout.data[8] = det[1] / det[0];
            fout.data[2] = det[2] / det[0];
            fout.data[5] = det[3] / det[0];
            fout.data[0] = (data[4]*data[8]-data[5]*data[7]) / det[0];
            fout.data[3] = (data[5]*data[6]-data[3]*data[8]) / det[0];
            fout.data[6] = (data[3]*data[7]-data[4]*data[6]) / det[0];
            fout.data[1] = -(data[1]*data[8]-data[2]*data[7]) / det[0];
            fout.data[4] = -(data[2]*data[6]-data[0]*data[8]) / det[0];
            fout.data[7] = -(data[0]*data[7]-data[1]*data[6]) / det[0];
        return fout;
        case 4:
            det[0] = data[10] * data[15] - data[11] * data[14]; // 00XX
            det[1] = data[11] * data[12] - data[8] * data[15]; // X00X
            det[2] = data[8] * data[13] - data[9] * data[12]; // XX00
            det[3] = data[9] * data[14] - data[10] * data[13]; // 0XX0
            det[4] = data[9] * data[15] - data[11] * data[13]; // 0X0X
            det[5] = data[8] * data[14] - data[10] * data[12]; // X0X0
            fout.data[0] = data[5] * det[0] - data[6] * det[4] + data[7] * det[3];
            fout.data[4] = -(data[4] * det[0] + data[6] * det[1] + data[7] * det[5]);
            fout.data[8] = data[4] * det[4] + data[5] * det[1] + data[7] * det[2];
            fout.data[12] = -(data[4] * det[3] - data[5] * det[5] + data[6] * det[2]);
            fout.data[1] = -(data[1] * det[0] - data[2] * det[4] + data[3] * det[3]);
            fout.data[5] = data[0] * det[0] + data[2] * det[1] + data[3] * det[5];
            fout.data[9] = -(data[0] * det[4] + data[1] * det[1] + data[3] * det[2]);
            fout.data[13] = data[0] * det[3] - data[1] * det[5] + data[2] * det[2];
            det[6] = data[1] * data[7] - data[3] * data[5]; det[7] = -det[6] * det[5]; // 0X0X
            det[5] = data[0] * data[6] - data[2] * data[4]; det[7] -= det[5] * det[4]; // X0X0
            det[4] = data[3] * data[4] - data[0] * data[7]; det[7] -= det[4] * det[3]; // X00X
            det[3] = data[2] * data[7] - data[3] * data[6]; det[7] += det[3] * det[2]; // 00XX
            det[2] = data[1] * data[6] - data[2] * data[5]; det[7] -= det[2] * det[1]; // 0XX0
            det[1] = data[0] * data[5] - data[1] * data[4]; det[7] += det[1] * det[0]; // XX00
            fout.data[2] = data[13] * det[3] - data[14] * det[6] + data[15] * det[2];
            fout.data[6] = -(data[12] * det[3] + data[14] * det[4] + data[15] * det[5]);
            fout.data[10] = data[12] * det[6] + data[13] * det[4] + data[15] * det[1];
            fout.data[14] = -(data[12] * det[2] - data[13] * det[5] + data[14] * det[1]);
            fout.data[3] = -(data[9] * det[3] - data[10] * det[6] + data[11] * det[2]);
            fout.data[7] = data[8] * det[3] + data[10] * det[4] + data[11] * det[5];
            fout.data[11] = -(data[8] * det[6] + data[9] * det[4] + data[11] * det[1]);
            fout.data[15] = data[8] * det[2] - data[9] * det[5] + data[10] * det[1];
            for(uint32_t i=0;i<16;i++) fout.data[i] /= det[7];
            return fout;
        }
    }
return fout;}
LFHTEMP TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> TMatrix<C,0u,0u,TF>::toInverse(){TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> dainv = this->mkInverse(); return this->toMemmove(dainv);}

LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::permuteRows(const Tuple<uint32_t> &permutation){
    if (permutation.getSize() != sizes[0]) {printf("Size mismatch\n"); return *this;}
    Tuple<C> buffer; buffer.setSize(permutation.getSize());
    int i,j;
    for(i=0;i<sizes[1];i++){
        for(j=0;j<sizes[0];j++) buffer[permutation[j]] = std::move(data[ i * sizes[0] + j]);
        for(j=0;j<sizes[0];j++) data[ i * sizes[0] + j] = std::move(buffer[j]);
    }
return *this;}
LFHTEMP TMatrix<C,0u,0u,TF>& TMatrix<C,0u,0u,TF>::permuteCols(const Tuple<uint32_t> &permutation){
    if (permutation.getSize() != sizes[1]) {printf("Size mismatch\n"); return *this;}
    Tuple<C> buffer; buffer.setSize(permutation.getSize());
    int i,j;
    for(i=0;i<sizes[0];i++){
        for(j=0;j<sizes[1];j++) buffer[permutation[j]] = std::move(data[ j * sizes[0] + i]);
        for(j=0;j<sizes[1];j++) data[ j * sizes[0] + i] = std::move(buffer[j]);
    }
return *this;}

LFHTEMP void TMatrix<C,0u,0u,TF>::show(FILE* f, int level)const{
    fprintf(f, "Matrix with sizes %i x %i\n", sizes[0], sizes[1]);
    unsigned int i,j;
    for(i=0;i<sizes[0];i++){
        for(j=0;j<sizes[1];j++){
            ExOp::show(data[i + j*sizes[0]],f,level+2);
            fprintf(f,"%c", (j == (sizes[1]-1)) ? '\n' : '\t');
        }
    }
    fprintf(f, "Matrix with sizes %i x %i\n", sizes[0], sizes[1]);
}
LFHTEMP ERRCODE TMatrix<C,0u,0u,TF>::load(FILE *f, unsigned int size){
    ERRCODE fout = (fread(sizes,sizeof(uint32_t),2,f)== 2) ? 0 : 1;
    this->setSizes(sizes[0], sizes[1]);
    uint32_t totsize = sizes[0] * sizes[1];
    if (totsize != 0) {
        if (ExCo<C>::IsPOD) {
            return fout |= (fread(data(),sizeof(C),totsize,f) == totsize ? 0 : 1);
        }else{
            for(uint32_t i =0; i< totsize;i++) if (fout = ExOp::load(data[i],f,size/totsize)) break;
        }
    }
    return fout;
}
LFHTEMP ERRCODE TMatrix<C,0u,0u,TF>::save(FILE *f) const{
    ERRCODE fout = (fwrite(sizes,sizeof(uint32_t),2,f) != 2) ? 0 : 1;
    if ((sizes[0] != 0)&&(sizes[1] != 0)) {
        uint32_t totsize = sizes[0] * sizes[1];
        if (totsize == 0) return 0;
        if (ExCo<C>::IsPOD) {
            return fout |= ((fwrite(data(),sizeof(C),totsize,f) != totsize) ? 0 : 1);
        }else{
            for(uint32_t i =0; i< totsize;i++) fout |= ExOp::save(data[i],f);
        }
    }
return fout;}

#ifdef Rcpp_hpp
LFHTEMP void TMatrix<C,0u,0u,TF>::rdMatrix(const Rcpp::NumericMatrix &where){
    uint32_t i,j;
    this->setSizes(where.nrow(),where.ncol());
    C* cur = data();
    for(j=0;j<where.ncol();j++){
        for(i=0;i<where.nrow();i++) *(cur++) = where(i,j);
    }
}

LFHTEMP void TMatrix<C,0u,0u,TF>::wrMatrix(Rcpp::NumericMatrix &where) const{
    where = Rcpp::NumericMatrix(sizes[0],sizes[1]);
    const C* cur = data();
    uint32_t i,j;
    for(j=0;j<where.ncol();j++){
        for(i=0;i<where.nrow();i++) where(i,j) = *(cur++);
    }
}
#endif




#undef LFHTEMP
#define LFHTEMP template<class C>

LFHTEMP template<class O> SparseTuple<C>::SparseTuple(const std::vector<O> &other){
}

LFHTEMP template<class D> C SparseTuple<C>::mkInnerMult(const SparseTuple<D>& other)const{
    uint32_t ite;
    uint32_t ind;
    C fout; fout.toZero();
    if (other.getSize() < this->getSize()){
        for(ite=0;ite<other.getSize();ite++){
            if ((ind = this->find(other.deref_key(ite)))!=0xFFFFFFFF) fout += other.deref(ite) * this->deref(ind);
        }
    }else{
        for(ite=0;ite<this->getSize();ite++){
            if ((ind = other.find(this->deref_key(ite)))!=0xFFFFFFFF) fout += other.deref(ind) * this->deref(ite);
        }
    }
return fout;}
LFHTEMP template<class D> SparseMatrix<C> SparseTuple<C>::mkOuterMult(const SparseTuple<D>& other)const{
    uint32_t ite;
    uint32_t ind;
    SparseMatrix<C> fout;
    /*if (other.getSize() < this->getSize()){
        for(ite=0;ite<other.getSize();ite++){
            if ((ind = this->find(other.deref_key(ite)))!=0xFFFFFFFF) fout += other.deref(ite) * this->deref(ind);
        }
    }else{
        for(ite=0;ite<this->getSize();ite++){
            if ((ind = other.find(this->deref_key(ite)))!=0xFFFFFFFF) fout += other.deref(ind) * this->deref(ite);
        }
    }*/
return fout;}
LFHTEMP template<class O> std::vector<O>& SparseTuple<C>::wrStdVector(std::vector<O>& fout)const{
    uint32_t cursize =0;
    for(uint32_t i=0;i<this->getSize();i++){


    }
}

LFHTEMP bool SparseTuple<C>::makeMaximumAsFirst(){
    if (this->getSize() <= 1u) return false;
    uint32_t maxcur = 0u;
    for(uint32_t i=1;i< this->getSize();i++) if (this->deref(i) > this->deref(maxcur)) maxcur =i;
    if (maxcur == 0u) return false;
    this->swapEntries(0u,maxcur);
return true;}

LFHTEMP template<class O> bool SparseTuple<C>::operator==(const SparseTuple<O>& other)const{
    uint32_t iteb;
    if (other.getSize() != this->getSize()) return false;
    if (auto ite = this->getIterator()) do{
        if ((iteb = other.find(ite())) == 0xFFFFFFFF) return false;
        if (*ite != other.deref(iteb)) return false;
    }while(ite++);
return true;}
LFHTEMP bool SparseTuple<C>::hasSameSparsity(const SparseTuple<C> &other)const{
    if (other.getSize() != this->getSize()) return false;
    if (auto ite = this->getIterator()) do{
        if (other.find(ite()) == 0xFFFFFFFF) return false;
    }while(ite++);
return(true);}


LFHTEMP SparseTuple<C>& SparseTuple<C>::toMult(const Trianglix<C,0u>& a, const SparseTuple<C>&b){
    uint32_t i,j,k;
    if (b.getSize() == 0) return this->toZero();
    for(i=0;i<this->getSize();i++){
        j=0;
        k = (b.deref_key(j) > this->deref_key(i)) ? ((b.deref_key(j)*(b.deref_key(j)+1))>>1) + this->deref_key(i) : ((this->deref_key(i)*(this->deref_key(i)+1))>>1) + b.deref_key(j);
        this->deref(i) = a.data[k] * b.deref(j);
        for(j++;j<b.getSize();j++){
            k = (b.deref_key(j) > this->deref_key(i)) ? ((b.deref_key(j)*(b.deref_key(j)+1))>>1) + this->deref_key(i) : ((this->deref_key(i)*(this->deref_key(i)+1))>>1) + b.deref_key(j);
            this->deref(i) += a.data[k] * b.deref(j);
        }
    }
return *this;}



#undef LFHTEMP
#define LFHTEMP template<class C>
LFHTEMP class SparseMatrix<C>::Task_getPrincipalComponents{
public:
    const SparseMatrix<C> &input;
    TMatrix<double> fout;
    Task_getPrincipalComponents(const SparseMatrix<C> &_input, int nbdims): input(_input){
        Tuple<uint32_t, 2u> coor;
        fout.setSizes(input.computeNBrows(),nbdims);
        fout.toRand();
    }
};
LFHTEMP TMatrix<double> SparseMatrix<C>::getPrincipalComponents(uint32_t nbcomp) const{
    SparseMatrix<C>::Task_getPrincipalComponents task(*this);
return task.fout;}

LFHTEMP TMatrix<double> SparseMatrix<C>::KendallColumns()const{TMatrix<double> fout;
    // rank of entries within columns are compared, with "0" being the only entry allowed to tie
    // the tau statistic counts the number of pair of pairs whose rank are concordant and discordant
    // \sum(i<j) (rank(x_i) > rank(y_i))^(rank(x_j) > rank(y_j))   - \sum(i<j) !(rank(x_i) > rank(y_i))^(rank(x_j) > rank(y_j))
    //
    // however equal ranks should be ignored, which changes the background distribution
    // assuming to ties does not exists, the distribution is normal(0, n(n-1)(2n-5)/2)
    //
    // we can normalize the difference to account for ties:

    // 2*(C - D) / sqrt((n * (n-1) - nb_0_in_x * (nb_0_in_x -1) )(n * (n-1) - nb_0_in_y * (nb_0_in_y -1)))

    // doing so the distribution of the above is normal in the limit:
    // (sqtr / )

    //


    uint32_t i,j,k,l;
    uint32_t ind;
    int32_t dadiff;
    uint32_t eqbuf[8];
    fout.setSizes(this->getNBcols(),this->getNBcols());
    Vector<Tuple<C, 2u> > inters;
    for(i=0; i< this->getNBcols();i++){
        for(j=i+1; j< this->getNBcols();j++){
            for(l= 0;l < 8;l++) eqbuf[l] = 0;
            if (auto ite = data[i].getIterator()) do{
                if ((ind = data[j].find(ite())) == 0xFFFFFFFF) {
                    eqbuf[ (*ite > 0) ? 4:5]++;
                    eqbuf[ (data[j].deref(ind) > 0) ? 6:7]++;
                    inters.push_back();
                    inters.last()[0] = *ite;
                    inters.last()[1] = data[j].deref(ind);
                }else eqbuf[(*ite > 0) ? 0 : 1]++;
            }while(ite++);
            if (auto ite = data[j].getIterator()) do{
                if ((ind = data[i].find(ite())) != 0xFFFFFFFF) eqbuf[(*ite > 0) ? 2 : 3]++;
            }while(ite++);

            dadiff = (int32_t)(eqbuf[2] * eqbuf[1] + eqbuf[0] * eqbuf[3]) - (int32_t)(eqbuf[3] * eqbuf[1] + eqbuf[0] * eqbuf[2]);

            for(l= 0;l < inters.getSize();l++){
                for(k= l+1;k < inters.getSize();k++) if ((inters[l][0] > inters[k][0])^(inters[l][1] > inters[k][1])) dadiff++; else dadiff--;
            }
            inters.toMemfree();
        }
    }
return fout;}

LFHTEMP Tuple<double> SparseMatrix<C>::UTestZScore(const Tuple<uint32_t> &listX, const Tuple<uint32_t> &listY, Tuple<double> *opt_logit_auroc, Tuple<double> * opt_logit_seen_auroc, bool print_prog, bool return_logitpval, bool hypergeo_correction)const{ // compute rank test for each column, for specific rows only statistics
    myHashmap<uint32_t, double> damass;
    Tuple<double> fout; fout.setSize(this->getNBcols());
    if (opt_logit_auroc) opt_logit_auroc->setSize(this->getNBcols());
    if (opt_logit_seen_auroc) opt_logit_seen_auroc->setSize(this->getNBcols());

    Tuple<uint32_t, 2u> zeros;
    Tuple<int64_t, 2u> tmpcnt;
    uint32_t ind,ind2;
    int64_t sum, zerosum;
    double tmpd, mu, erffact, totn, tie;
    myHashmap<C, Tuple<uint32_t, 2u> > cnts;
	ProgressBarPrint progbar(20);
    if (print_prog) progbar.start("Performing Hypergeo-U test");

    myHashmap<uint32_t, uint32_t> partit;
    if (auto ite = listX.getIterator()) do{
        partit[*ite] = 0;
    }while(ite++);
    if (auto ite = listY.getIterator()) do{
        partit[*ite] = 1;
    }while(ite++);
    totn = listX.getSize() + listY.getSize();


    for(uint32_t i=0;i<this->getNBcols();i++){
        if (print_prog) progbar.update(((double)i)/ this->getNBcols());
        if (data[i].getSize() < totn){ // faster by observations
            zeros[0] = listX.getSize();
            zeros[1] = listY.getSize();
            if (auto ite = data[i].getIterator()) do{
                if ((ind = partit.find(ite())) != 0xFFFFFFFF){
                    zeros[partit.deref(ind)]--;
                    cnts[*ite][partit.deref(ind)]++;
                }
            }while(ite++);
        }else{
            zeros.toZero();
            if (auto ite = listX.getIterator()) do{
                if ((ind = data[i].find(*ite)) == 0xFFFFFFFF) zeros[0]++;
                else cnts[data[i].deref(ind)][0]++;
            }while(ite++);
            if (auto ite = listY.getIterator()) do{
                if ((ind = data[i].find(*ite)) == 0xFFFFFFFF) zeros[1]++;
                else cnts[data[i].deref(ind)][1]++;
            }while(ite++);
        }
        if (!hypergeo_correction){
            cnts[(C)0][0] = zeros[0]; zeros[0] =0;
            cnts[(C)0][1] = zeros[1]; zeros[1] =0;
        }

        cnts.keysort();

        // HyperGeometric parameters:
        // N = listX.getSize() + listY.getSize()
        // K = zeros[0] + zeros[1]
        // n = listX.getSize()
        // k = zeros[0]



        if ((hypergeo_correction)&&((cnts.getSize() < 2)||(zeros[0] == listX.getSize())||(zeros[1] == listY.getSize()))) {
            if (cnts.getSize() == 0) {
                fout[i] = 0.0;
                if (opt_logit_auroc) (*opt_logit_auroc)[i] = 0.0;
            }else{ // true hypergeometric! only 2 values or 1 sample has only missing data
                fout[i] = -HypergeomericLogitPval(zeros[0], zeros[0] + zeros[1], listX.getSize(),listX.getSize() + listY.getSize());
                if (opt_logit_auroc){
                    tmpd = ((1.0 / listY.getSize()) * zeros[1]) - ((1.0 / listX.getSize())* zeros[0]);
                    (*opt_logit_auroc)[i] = log(1.0 - tmpd) - log(1.0 + tmpd);
                }
            }
            if (opt_logit_seen_auroc) (*opt_logit_seen_auroc)[i] = 0.0;
        }else{ // only zeroes== logit pval = 0

            tmpcnt[0] = listX.getSize() - zeros[0];
            tmpcnt[1] = listY.getSize() - zeros[1];
            zerosum = (tmpcnt[0] * zeros[1]) - (tmpcnt[1] * zeros[0]);
            // apply wilcox normal approximation to non-zero entries, first compute AUROC there
            double tie =0.0;
// library(InferN0) ; scptr <- infern0LoadFile("frozenallctlrs125.inf.scp"); tmp <-readRDS("tmptmp.rds"); output <- infern0ZeroCorrectedWilcox(scptr, tmp  == "Microglia", tmp  != "Microglia")

            // log( sum ( x_i ) ) +  log( - sum ( x_i ) )

            //log(ab) = log(a) +log(b)

            for(sum =0, ind=0;ind<cnts.getSize();ind++){
                tmpcnt -= cnts.deref(ind);
                tmpd = cnts.deref(ind)[1] + cnts.deref(ind)[0];
                tie += tmpd * (tmpd*tmpd -1.0);
                sum +=  (tmpcnt[0] * cnts.deref(ind)[1]) - (tmpcnt[1] * cnts.deref(ind)[0]);
                //if (i < 10) printf("%i[%i]: %i %i  (%i %i %i %i)\n", i, cnts.deref_key(ind), (tmpcnt[0] * cnts.deref(ind)[1]), (tmpcnt[1] * cnts.deref(ind)[0]), tmpcnt[0], cnts.deref(ind)[1], tmpcnt[1], cnts.deref(ind)[0]);
            }
            ind = listX.getSize() + listY.getSize() - zeros[0] - zeros[1];
            erffact = (((double)(listX.getSize() - zeros[0]) ) * (listY.getSize() - zeros[1])) * (((-tie / ind) / (ind-1)) + (ind+1) ) / 3.0;
            //printf("%i\t%e\t", i,erf(erffact * sum));
        // printf("zsum %i and sum %i\n", (int32_t)zerosum, (int32_t)sum);
        // printf("auroc = %e\n", 0.5 + 0.5 * ((((double)(zerosum + sum)) / listX.getSize()) / listY.getSize()));

            mu =   (((double)zerosum + sum) - (zeros[0] + zeros[1]) * listX.getSize());
        // printf("zsum %i and sum %i\n", (int32_t)zerosum, (int32_t)sum);
        // printf("auroc = %e\n", 0.5 + 0.5 * ((((double)(zerosum + sum)) / listX.getSize()) / listY.getSize()));

            /*damass = getHypergeomericMasses(zeros[0] + zeros[1], listX.getSize(),totn,0.99999999);
            if ( fabs(erffact * (mu + totn * damass.deref_key(0))) < -3.0){

            }else{
                tmpd = 0.0;  tie = 0.0;
                if (typename myHashmap<uint32_t, double>::Iterator ite = damass.getIterator()) do{
                    tmpd += (*ite) * erf(erffact * (mu + totn * ite())  ); // (*ite) * erf(erffact * (mu/ totn + ite()) * totn  );
                    tie += (*ite);
                }while(ite++);
                fout[i] = tmpd / tie;
                fout[i] = log(1.0 + fout[i]) - log(1.0 - fout[i]);
            }*/
            fout[i] = HypergeomericNormalLogitPval(-mu / totn,sqrt(erffact) / totn,zeros[0] + zeros[1],listX.getSize(),totn);

            if (!ExOp::isValid(tmpd)){
                printf("%i\t%i\t%i\t%i\n", zeros[0],zeros[1], listX.getSize(), listY.getSize());
                printf("%i\t%i\n", sum, zerosum);
                cnts.show();
                LFH_ALIVE; exit(1);
            }
            if (opt_logit_auroc){
                tmpd = zerosum + sum;
                mu = ((double)listX.getSize()) * listY.getSize();
                (*opt_logit_auroc)[i] = log(mu - tmpd) - log(mu + tmpd);
            }
            if (opt_logit_seen_auroc){
                tmpd = sum;
                mu = ((double)(listX.getSize() - zeros[0])) * (listY.getSize() - zeros[1]);
                (*opt_logit_seen_auroc)[i] = log(mu - tmpd) - log(mu + tmpd);
            }
            //printf("%e\t%e\n", fout[i], (*opt_logit_auroc)[i]);
        }
        cnts.toMemfree();
    }
    if (print_prog) progbar.finish();
    if (!return_logitpval){
        for(uint32_t i=0; i<fout.getSize();i++) fout[i] = logitPval_to_Zscore(fout[i]);
    }
return fout;}

// for each column, sample |list_p| from list_P and |list_n| from list_N
LFHTEMP Tuple<double> SparseMatrix<C>::UTestSamplingPvalue(const Tuple<double> &logitPval, uint32_t nb_P, uint32_t nb_N, const Tuple<uint32_t> &list_P, const Tuple<uint32_t> &list_N, uint32_t nbsubsamples)const{
    // compares the provided Utest rank statistic (as logitPval) to simililar rank statistic obtained some random samples from list_P and list_N
    myHashmap<uint32_t, double> damass;
    Tuple<double> fout; fout.setSize(this->getNBcols());
    Tuple<uint32_t, 2u> omegazeros;

    Tuple<uint32_t> sampledzeros_P; sampledzeros_P.setSize(nbsubsamples);
    Tuple<uint32_t> sampledzeros_N; sampledzeros_N.setSize(nbsubsamples);
    Tuple<uint32_t, 2u> zeros;
    Tuple<int64_t, 2u> tmpcnt;
    uint32_t ind,ind2;
    int64_t sum, zerosum;
    double tmpd, mu, erffact, totn, tie;
    myHashmap<C, Tuple<uint32_t, 2u> > cnts;


    myHashmap<uint32_t, uint32_t> partit;
    if (auto ite = list_P.getIterator()) do{
        partit[*ite] = 0;
    }while(ite++);
    if (auto ite = list_N.getIterator()) do{
        partit[*ite] = 1;
    }while(ite++);
    totn = list_P.getSize() + list_N.getSize();

    for(uint32_t i=0;i<this->getNBcols();i++){
        // both listP and listN are sampled
        if (data[i].getSize() < totn){ // faster by observations
            omegazeros[0] = list_P.getSize();
            omegazeros[1] = list_N.getSize();
            if (auto ite = data[i].getIterator()) do{
                if ((ind = partit.find(ite())) != 0xFFFFFFFF){
                    omegazeros[partit.deref(ind)]--;
                    cnts[*ite][partit.deref(ind)]++;
                }
            }while(ite++);
        }else{
            omegazeros.toZero();
            if (auto ite = list_P.getIterator()) do{
                if ((ind = data[i].find(*ite)) == 0xFFFFFFFF) omegazeros[0]++;
                else cnts[data[i].deref(ind)][0]++;
            }while(ite++);
            if (auto ite = list_N.getIterator()) do{
                if ((ind = data[i].find(*ite)) == 0xFFFFFFFF) omegazeros[1]++;
                else cnts[data[i].deref(ind)][1]++;
            }while(ite++);
        }
        // first, sample the number of zeros in both distributions
        getHypergeomericCover(sampledzeros_P,nb_P,omegazeros[0], list_P.getSize());
        getHypergeomericCover(sampledzeros_N,nb_N,omegazeros[1], list_N.getSize());



    }
return fout;}
LFHTEMP TMatrix<Zscore>  SparseMatrix<C>::UTestZScoreSplit(const Tuple< Tuple<uint32_t> > &lists, const Tuple<uint32_t> &ordering, uint32_t nbthreads, uint32_t splitgroups, const Tuple<uint32_t> &partit, TMatrix<double> *opt_logit_auroc,TMatrix<double> *opt_mean, TMatrix<double> *opt_drop_enrich,  bool print_prog, bool hypergeo_correction, bool downsample)const{ // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported
    TMatrix<Tuple<uint32_t> > splitlists;
    Tuple<uint32_t, 2u> coor;
    TMatrix<Zscore> fout;
    uint32_t i,l,k;
    Vector<KeyElem<uint32_t, uint32_t> > sortt;
    myHashmap<uint32_t, uint32_t> dapart;
    if (partit.getSize() == ordering.getSize()){
        for(i=0;i< partit.getSize(); i++) dapart[partit[i]]++;

    }else dapart[0] = 1;
    ThreadBase tb; tb.setSize(nbthreads);
  //  if (opt_logit_auroc) tb.flags["logit_auroc"] = opt_logit_auroc;
  //  if (opt_mean) tb.flags["mean"] = opt_logit_auroc;
  //  if (opt_drop_enrich) tb.flags["drop_enrich"] = opt_drop_enrich;
  //  if (return_logitpval) tb.flags["return_logitpval"] = true;
    if (hypergeo_correction) tb.flags["hypergeo_correction"] = true;
    if (print_prog) tb.flags["print.progress"] = true;
    tb.flags["print.progress"] = true;

    TMatrix<double> dasampling;

    splitlists.setSizes(lists.getSize() , splitgroups * dapart.getSize() );
    dasampling.setSizes(lists.getSize() , splitgroups * dapart.getSize() ).toZero();
    if (dapart.getSize() == 1){
        for(coor[0]=0;coor[0]<lists.getSize();coor[0]++){
            sortt.setSize(lists[coor[0]].getSize());
            for(i=0;i<lists[coor[0]].getSize();i++) {sortt[i].k = ordering[lists[coor[0]][i]]; sortt[i].d = lists[coor[0]][i];}
            sortt.sort();
            uint32_t j = lists[coor[0]].getSize() / splitgroups;
            k =0u;
            for(coor[1] =0;coor[1] < (j+1) * splitgroups - lists[coor[0]].getSize();coor[1]++) {
                splitlists(coor).setSize(j);
                for(i=0;i<splitlists(coor).getSize();i++) {splitlists(coor)[i] = sortt[k].d; dasampling(coor) += sortt[k].k; k++;}
                dasampling(coor) /= j;
            }
            for(;coor[1] <splitgroups;coor[1]++){
                splitlists(coor).setSize(j+1);
                for(i=0;i<splitlists(coor).getSize();i++) {splitlists(coor)[i] = sortt[k].d; dasampling(coor) += sortt[k].k; k++;}
                dasampling(coor) /= j+1;
            }
            sortt.toMemfree();
        }
    }else{
        Tuple< Vector<uint32_t> > subpart;
        subpart.setSize(dapart.getSize());
        for(coor[0]=0;coor[0]<lists.getSize();coor[0]++){
            for(i=0;i<lists[coor[0]].getSize();i++) {subpart[dapart.find(lists[coor[0]][i])].push_back(lists[coor[0]][i]);}
            coor[1] =0;
            for(l=0;l<dapart.getSize();l++){
                sortt.setSize(subpart[l].getSize());
                for(i=0;i<subpart[l].getSize();i++) {sortt[i].k = ordering[subpart[l][i]]; sortt[i].d = subpart[l][i];}
                sortt.sort();
                uint32_t j = subpart[l].getSize() / splitgroups;
                k =0u;
                for(;coor[1] < (j+l+1) * splitgroups - subpart[l].getSize();coor[1]++) {
                    splitlists(coor).setSize(j);
                    for(i=0;i<splitlists(coor).getSize();i++) {splitlists(coor)[i] = sortt[k].d; dasampling(coor) += sortt[k].k; k++;}
                    dasampling(coor) /= j;
                }
                for(;coor[1] < (l +1) * splitgroups;coor[1]++){
                    splitlists(coor).setSize(j+1);
                    for(i=0;i<splitlists(coor).getSize();i++) {splitlists(coor)[i] = sortt[k].d; dasampling(coor) += sortt[k].k; k++;}
                    dasampling(coor) /= j+1;
                }
                sortt.toMemfree();

                coor[0] += lists.getSize();
                subpart[l].toMemfree();
            }
            coor[0] -= lists.getSize() * dapart.getSize();
        }
    }
    printf("Mean sequencing depth:\n");
    ExOp::show(dasampling);
    double tmp;

    if (downsample){
        for(coor[1] =0;coor[1] < splitgroups * dapart.getSize();coor[1]++) {
            coor[0]=0;
            tmp = dasampling(coor);
            for(coor[0]++;coor[0]<lists.getSize();coor[0]++){
                if (dasampling(coor) < tmp) tmp = dasampling(coor);
            }
            for(coor[0]=0;coor[0]<lists.getSize();coor[0]++) dasampling(coor) = tmp / dasampling(coor);
        }
//        printf("Downsample probs:\n"); ExOp::show(dasampling);
    }

    /*
    Vector<uint32_t> sumsum;
    for(coor[0]=0;coor[0]<this->getNBcols();coor[0]++){
        for(coor[1]=0;coor[1]<data[coor[0]].getSize();coor[1]++){
            while(sumsum.getSize() <= data[coor[0]].deref_key(coor[1])) sumsum.push_back(0);
            sumsum[data[coor[0]].deref_key(coor[1])] += data[coor[0]].deref(coor[1]);
            }
    }

    coor[0]=0;
    WeightElem<double,2u> dathing;

    TMatrix<double> tmpreport; tmpreport.setSizes(lists.getSize()* dapart.getSize(),splitgroups);

    if (auto mite = splitlists.getIterator()) do{
       dathing.toZero();
       for(uint32_t k =0;k<mite->getSize();k++) {
           // if (splitlists(coor)[k] >= sumsum.getSize()) {printf("got weird! %i\n ", splitlists(coor)[k]);}else

            dathing += WeightElem<double, 2u>(sumsum[(*mite)[k]], 1.0);
       }
       tmpreport(mite()) = dathing.getMean();
       //printf("%e mean%e std%e\n", dathing.w[0], dathing.getMean(),sqrt(dathing.getVar()));
    }while(mite++);*/

	coor[0]=0;

//    for(coor[1] =0;coor[1] < splitgroups;coor[1]++) {printf("show which... %i\n",splitlists(coor).getSize());splitlists(coor).show();	}
    Tuple<TMatrix<double> > daout2;
    Tuple<TMatrix<double> > daout3;
    Tuple<TMatrix<double> > daout4;
    Tuple<TMatrix<double> > daout = this->UTestZScoreMulti(tb, splitlists, (opt_logit_auroc != NULL)? &daout2: NULL,  (opt_drop_enrich != NULL)? &daout3: NULL, (opt_mean != NULL)? &daout4: NULL , (downsample) ? &dasampling : NULL);
    double extrem[2]; uint32_t countnan;
    extrem[0] = extrem[1] = 0; countnan=0;
    for(i=0;i< daout.getSize();i++){
        if (auto mite = daout[i].getIterator()) do{
            if (!ExOp::isValid(*mite)) {countnan++; continue;}
            if  (*mite < extrem[0] ) extrem[0] = *mite;
            if  (*mite > extrem[1] ) extrem[1] = *mite;
        }while(mite++);
    }
//    printf("pre extreme Zscores are : %e to %e and %i nans\n", extrem[0], extrem[1], countnan);


    TMatrix<double> Zweight; Zweight.setSizes(lists.getSize() , splitgroups);
    TMatrix<uint32_t> Znbcell; Znbcell.setSizes(lists.getSize() , splitgroups);

    if (dapart.getSize() != 1){ // merging partitions

        if (auto mite = Zweight.getIterator()) do{
            coor = mite();
            l = 0;
            k = splitlists(coor).getSize();
            for(l++, coor[1] += splitgroups;l < dapart.getSize();l++, coor[1] += splitgroups) k += splitlists(coor).getSize();
            Znbcell(mite()) = k;
            coor = mite();
            Zweight(mite()) = ((double)splitlists(coor).getSize()) * (k- splitlists(coor).getSize());
            l=0;
            for(l++, coor[1] += splitgroups;l < dapart.getSize();l++, coor[1] += splitgroups) {
                Zweight(mite()) += ((double)splitlists(coor).getSize()) * (k- splitlists(coor).getSize());
            }
        }while(mite++);

        for(i=0;i<daout.getSize();i++){
            TMatrix<double> swout[4];
            swout[0].setSizes(lists.getSize() , splitgroups);
            swout[1].setSizes(lists.getSize() , splitgroups);
            swout[2].setSizes(lists.getSize() , splitgroups);
            swout[3].setSizes(lists.getSize() , splitgroups);
            if (auto mite = swout[0].getIterator()) do{
                l=0;
                coor = mite();
                double tmp = ((double)splitlists(coor).getSize()) * (Znbcell(mite()) - splitlists(coor).getSize());
                double sum = tmp;
                swout[0](mite()) = daout[i](coor) * sqrt(tmp);
                swout[1](mite()) = daout2[i](coor) * tmp;
                swout[2](mite()) = daout3[i](coor) * tmp;
                swout[3](mite()) = daout4[i](coor) * tmp;
                for(l++, coor[1] += splitgroups;l < dapart.getSize();l++, coor[1] += splitgroups) {
                    tmp = ((double)splitlists(coor).getSize()) * (Znbcell(mite())- splitlists(coor).getSize());
                    sum += tmp;
                    swout[0](mite()) += daout[i](coor) * sqrt(tmp);
                    swout[1](mite()) += daout2[i](coor) * tmp;
                    swout[2](mite()) += daout3[i](coor) * tmp;
                    swout[3](mite()) += daout4[i](coor) * tmp;
                }
                swout[0](mite()) /= sqrt(sum);
                swout[1](mite()) /= sum;
                swout[2](mite()) /= sum;
                swout[3](mite()) /= sum;
            }while(mite++);
        }
    }else{
        if (auto mite = splitlists.getIterator()) do{
           // printf("%i,%i: \t%i *%i\n",mite()[0],mite()[1], mite->getSize(), lists[mite()[0]].getSize()- mite->getSize());
            Zweight(mite()) = 1.0;
        }while(mite++);
    }

    fout.setSizes(daout.getSize(), splitgroups * lists.getSize());
    for(k =0u,coor[1] =0; coor[1] < splitgroups * lists.getSize(); coor[1]++){
        for(coor[0] =0; coor[0] < daout.getSize(); coor[0]++, k++) {
            double tmpweight = Zweight.data[coor[1]];
            fout.data[k] = Zscore(daout[coor[0]].data[coor[1]] * sqrt(tmpweight), tmpweight);
        }
    }
    if (opt_logit_auroc != NULL){
        opt_logit_auroc->setSizes(daout.getSize(), splitgroups * lists.getSize());
        for(k =0u,coor[1] =0; coor[1] < splitgroups * lists.getSize(); coor[1]++){
            for(coor[0] =0; coor[0] < daout.getSize(); coor[0]++, k++) opt_logit_auroc->data[k] = daout2[coor[0]].data[coor[1]];
        }
    }
	if (opt_mean != NULL){
        opt_mean->setSizes(daout.getSize(), splitgroups * lists.getSize());
        for(k =0u,coor[1] =0; coor[1] < splitgroups * lists.getSize(); coor[1]++){
            for(coor[0] =0; coor[0] < daout.getSize(); coor[0]++, k++) opt_mean->data[k] = daout3[coor[0]].data[coor[1]];
        }
    }
    if (opt_drop_enrich != NULL){
        opt_drop_enrich->setSizes(daout.getSize(), splitgroups * lists.getSize());
        for(k =0u,coor[1] =0; coor[1] < splitgroups * lists.getSize(); coor[1]++){
            for(coor[0] =0; coor[0] < daout.getSize(); coor[0]++, k++) opt_drop_enrich->data[k] = daout4[coor[0]].data[coor[1]];
        }
    }

    if (auto mite = fout.getIterator()){
        extrem[0] = extrem[1] = (double)*mite;countnan=0;
        while(mite++){
            double tmp = (double)*mite;
            if (!ExOp::isValid(tmp)) {countnan++; continue;}
            if  (tmp < extrem[0] ) extrem[0] = tmp;
            if  (tmp > extrem[1] ) extrem[1] = tmp;
        }

  //      printf("post extreme Zscores are : %e to %e and %i nans\n", extrem[0], extrem[1], countnan);
    }
return(fout);}
LFHTEMP  Tuple<TMatrix<double> > SparseMatrix<C>::UTestZScoreMulti(ThreadBase &tb, const TMatrix<Tuple<uint32_t> > &list_grid, Tuple< TMatrix<double> > *opt_logit_auroc,Tuple<TMatrix<double> > *opt_mean, Tuple<TMatrix<double> > *opt_drop_enrich, const TMatrix<double> *downsample)const{ // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported
class Task : public Event<uint32_t>{
public:
    ThreadBase &tb;
    const SparseMatrix<C> &data;
    const TMatrix<Tuple<uint32_t> > & list_grid;
    Tuple<TMatrix<double> > *opt_logit_auroc;
    Tuple<TMatrix<double> > *opt_mean;
    Tuple<TMatrix<double> > *opt_drop_enrich;
    Tuple<uint32_t> colrange;
    Tuple<TMatrix<double> > fout;
    myHashmap<uint32_t, uint32_t> partit;
    Tuple<uint32_t> sizes;
    Tuple<uint32_t> row_sizes;
    uint32_t totallength;

    uint32_t* downsampling;
    ~Task(){if (downsampling != NULL) delete[](downsampling);}

    Task(ThreadBase &_tb,const SparseMatrix<C>& _data, const TMatrix<Tuple<uint32_t> > & _list_grid,Tuple<TMatrix<double> > *_opt_logit_auroc, Tuple<TMatrix<double> > *_opt_mean,Tuple<TMatrix<double> > *_opt_drop_enrich):tb(_tb),data(_data),list_grid(_list_grid),opt_logit_auroc(_opt_logit_auroc), opt_mean(_opt_mean),opt_drop_enrich(_opt_drop_enrich),totallength(0u),downsampling(NULL){}
    uint32_t operator()(uint32_t threadID){
    myHashmap<uint32_t, double> damass;

    Tuple<uint32_t> zeros; zeros.setSize(list_grid.sizes[0]*list_grid.sizes[1]);
    int32_t row_zeros;
    Tuple<int64_t> tmpcnt;
    int64_t tmpcntsum[2];
    uint32_t ind,ind2;
    Tuple<int64_t> sum, zerosum;
    sum.setSize(list_grid.sizes[0]);
    zerosum.setSize(list_grid.sizes[0]);
    tmpcnt.setSize(list_grid.sizes[0]);
    double tmpd, mu, erffact, totn, tie;
    Tuple< myHashmap<C, Tuple<uint32_t> > > cnts; cnts.setSize(list_grid.sizes[1]);
    Tuple<uint32_t, 3u> coor3;
    uint32_t newval;
    uint32_t darand;
    for(uint32_t i=colrange[threadID];i<colrange[threadID+1];i++,tb.updateProgress(threadID)){ // every gene
        fout[i].setSizes(list_grid.sizes[0],list_grid.sizes[1]);
        if (opt_logit_auroc) (*opt_logit_auroc)[i].setSizes(list_grid.sizes[0],list_grid.sizes[1]);
        if (opt_mean) (*opt_mean)[i].setSizes(list_grid.sizes[0],list_grid.sizes[1]);
        if (opt_drop_enrich) (*opt_drop_enrich)[i].setSizes(list_grid.sizes[0],list_grid.sizes[1]);
        zeros = sizes;
        // step 1, fill count hashmaps
        if (downsampling){
            if (auto ite = data.data[i].getIterator()) do{
                if ((ind = partit.find(ite())) != 0xFFFFFFFF){
                    ind = partit.deref(ind);
                    if (downsampling[ind << 4] == 0xFFFFFFFF) newval = (*ite);
                    else{
                        newval =0;
                        for(coor3[2]=0 ; coor3[2] < ((*ite) % 15); coor3[2]++){
                            if (ExOp::toRand(darand) <= downsampling[ind << 4]) newval++;
                        }
                        for(coor3[2]=0 ; coor3[2] < (*ite) / 15; coor3[2]++){
                            if (ExOp::toRand(darand) <= downsampling[(ind << 4) | 8]){
                                if (darand <= downsampling[(ind << 4) | 4]){
                                    if (darand <= downsampling[(ind << 4) | 2]) newval += (darand <= downsampling[(ind << 4) | 1])? 0 : 1;
                                    else newval += (darand <= downsampling[(ind << 4) | 3])? 2 : 3;
                                }else if (darand <= downsampling[(ind << 4) | 6]) newval += (darand <= downsampling[(ind << 4) | 5])? 4 : 5;
                                else newval += (darand <= downsampling[(ind << 4) | 7])? 6 : 7;
                            }else if (darand <= downsampling[(ind << 4) | 12]){
                                if (darand <= downsampling[(ind << 4) | 10]) newval += (darand <= downsampling[(ind << 4) | 9])? 8 : 9;
                                else newval += (darand <= downsampling[(ind << 4) | 11])? 10 : 11;
                            }else if (darand <= downsampling[(ind << 4) | 14]) newval += (darand <= downsampling[(ind << 4) | 13])? 12 : 13;
                            else newval += (darand <= downsampling[(ind << 4) | 15])? 14 : 15;
                        }
                    }
                    if (newval != 0){
                        zeros[ind]--;
                        if ((ind2 = cnts[ind / list_grid.sizes[0]].find(newval)) == 0xFFFFFFFF){
                            cnts[ind / list_grid.sizes[0]][newval].setSize(list_grid.sizes[0]).toZero();
                            cnts[ind / list_grid.sizes[0]][newval][ind % list_grid.sizes[0]]++;
                        }else cnts[ind / list_grid.sizes[0]].deref(ind2)[ind % list_grid.sizes[0]]++;
                    }
                }
            }while(ite++);

        }else{

            if (auto ite = data.data[i].getIterator()) do{
                if ((ind = partit.find(ite())) != 0xFFFFFFFF){
                    ind = partit.deref(ind);
                    zeros[ind]--;
                    if ((ind2 = cnts[ind / list_grid.sizes[0]].find(*ite)) == 0xFFFFFFFF){
                        cnts[ind / list_grid.sizes[0]][*ite].setSize(list_grid.sizes[0]).toZero();
                        cnts[ind / list_grid.sizes[0]][*ite][ind % list_grid.sizes[0]]++;
                    }else cnts[ind / list_grid.sizes[0]].deref(ind2)[ind % list_grid.sizes[0]]++;
                }
            }while(ite++);
        }

        for(coor3[2] =0;coor3[2] < list_grid.sizes[1];coor3[2]++){
            row_zeros = 0;
            if (tb.flags.find("hypergeo_correction") != 0xFFFFFFFF) {
                for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++) row_zeros += zeros[coor3[0]];
                tmpcntsum[0] = 0u;

                for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                    tmpcntsum[0] += tmpcnt[coor3[1]] = sizes[coor3[0]] - zeros[coor3[0]]; // non-zero for current N
                    zerosum[coor3[1]] = (((int64_t)tmpcnt[coor3[1]]) * (row_zeros - zeros[coor3[0]])) - ((int64_t)((row_sizes[coor3[2]] - sizes[coor3[0]]) + zeros[coor3[0]] - row_zeros) * zeros[coor3[0]]);

                }
            }else{
                 //               printf("%i vs %i vs %i vs %i\n",zeros.getSize(),coor3[2] * list_grid.sizes[0], tmpcnt.getSize() ); fflush(stdout);
                row_zeros = 0;
                tmpcntsum[0] = 0u;
                for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                        if (cnts[coor3[2]].find(0u) == 0xFFFFFFFF) cnts[coor3[2]][0].setSize(list_grid.sizes[0]).toZero();
                        cnts[coor3[2]][0][coor3[1]] = zeros[coor3[0]];
                        zeros[coor3[0]] =0u;
                        zerosum[coor3[1]] = 0u;
                        tmpcnt[coor3[1]] = sizes[coor3[0]];
                        tmpcntsum[0] += tmpcnt[coor3[1]];
                }
            }

            // apply wilcox normal approximation to non-zero entries, first compute AUROC there
            cnts[coor3[2]].keysort();
            double tie =0.0;
            sum.toZero();
            for(ind=0;ind<cnts[coor3[2]].getSize();ind++){
                tmpcntsum[ (ind&1) ^ 1 ] = tmpcntsum[ind&1];
                for(coor3[1] = 0; coor3[1] < list_grid.sizes[0]; coor3[1]++){
                    tmpcnt[coor3[1]] -= cnts[coor3[2]].deref(ind)[coor3[1]];
                    tmpcntsum[(ind&1) ^ 1] -= cnts[coor3[2]].deref(ind)[coor3[1]];
                }
                tmpcntsum[ ind&1 ] -= tmpcntsum[ (ind&1) ^ 1 ];
                for(coor3[1] = 0; coor3[1] < list_grid.sizes[0]; coor3[1]++){
    //                    if ((i < 10)&&(coor3[1]==0)) printf("%i[%i]: %i %i  (%i %i %i %i)\n", i, cnts[coor3[2]].deref_key(ind), ((int64_t)tmpcnt[coor3[1]] * (tmpcntsum[ ind&1 ] - cnts[coor3[2]].deref(ind)[coor3[1]]))  ,(((int64_t)tmpcntsum[ (ind&1) ^ 1 ] - tmpcnt[coor3[1]]) * cnts[coor3[2]].deref(ind)[coor3[1]] )                    , tmpcnt[coor3[1]], tmpcntsum[ ind&1 ] - cnts[coor3[2]].deref(ind)[coor3[1]]  ,tmpcntsum[ (ind&1) ^ 1 ] - tmpcnt[coor3[1]], cnts[coor3[2]].deref(ind)[coor3[1]]);
                    sum[coor3[1]] +=  ((int64_t)tmpcnt[coor3[1]] * (tmpcntsum[ ind&1 ] - cnts[coor3[2]].deref(ind)[coor3[1]])) - (((int64_t)tmpcntsum[ (ind&1) ^ 1 ] - tmpcnt[coor3[1]]) * cnts[coor3[2]].deref(ind)[coor3[1]] );
                }

                tmpd = (double)(tmpcntsum[ ind&1 ]); tie += tmpd * (tmpd*tmpd -1.0);
            }
            if (opt_mean != NULL){
                for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++) (*opt_mean)[i].data[coor3[0]] = 0.0;
                for(ind=0;ind<cnts[coor3[2]].getSize();ind++){
                    for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                       (*opt_mean)[i].data[coor3[0]] += ((double) cnts[coor3[2]].deref_key(ind)) * cnts[coor3[2]].deref(ind)[coor3[1]];
                    }
                }
                /*if (i == 99) {
                    printf("test average!\n");
                    cnts[coor3[2]].show();
                    printf("is the following:\n");
                    for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++) printf("%e= %e / (%i + %i)\n", (*opt_mean)[i].data[coor3[0]] / (sizes[coor3[0]] + zeros[coor3[0]]), (*opt_mean)[i].data[coor3[0]] , sizes[coor3[0]] , zeros[coor3[0]]);
                }*/
                for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++) (*opt_mean)[i].data[coor3[0]] /= (sizes[coor3[0]] + zeros[coor3[0]]);

            }
                cnts[coor3[2]].toMemfree();
                ind = row_sizes[coor3[2]] - row_zeros;
                if (ind == 0){ // row only has zeros...
                    for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++) fout[i].data[coor3[0]] = 0.0;
                }else{
                   // if (hypergeo_correction){
                        tie = (((-tie / ind) / (ind-1)) + (ind+1) ) / 3.0;
                        if (!ExOp::isValid(tie)) tie = 0.0; // uses hypergeometric!
                        //printf("%i -> %e tie... non-zero %i\n", row_sizes[coor3[2]] - row_zeros,tie, (int)ind);
                        for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                            erffact = (((double)(sizes[coor3[0]] - zeros[coor3[0]]) ) * (row_sizes[coor3[2]] + zeros[coor3[0]] - sizes[coor3[0]] - row_zeros)) * tie;
                            //printf("erffact %e %e %e\n", erffact, (double)(sizes[coor3[0]] - zeros[coor3[0]]) , (double)((row_sizes[coor3[2]] + zeros[coor3[0]] - sizes[coor3[0]] - row_zeros)));
                            mu =  ((double)zerosum[coor3[1]] + sum[coor3[1]]) - row_zeros * sizes[coor3[0]];
                    //        printf("%i and %i\n", zerosum[coor3[1]], sum[coor3[1]] );
                            //mu =   (((double)zerosum + sum) - (zeros[0] + zeros[1]) * listP.getSize());
                            fout[i].data[coor3[0]] = HypergeomericNormalLogitPval(-mu / row_sizes[coor3[2]],sqrt(erffact) / row_sizes[coor3[2]],row_zeros,sizes[coor3[0]],row_sizes[coor3[2]]);
                            if (!ExOp::isValid(fout[i].data[coor3[0]])) printf("NAN generated by %i %i %i mu %e %e\n", row_sizes[coor3[2]],row_zeros,sizes[coor3[0]], mu,  erffact);
                    //        printf("%i %i %i, %e %e %e (%i %i %i)\n", i, coor3[2], coor3[1], fout[i].data[coor3[0]], -mu / row_sizes[coor3[2]],sqrt(erffact) / row_sizes[coor3[2]],row_zeros,sizes[coor3[0]],row_sizes[coor3[2]]);
                        }
                  /*  }else{ // NOT WORKING!
                        tmpd = (double)row_zeros; tie += tmpd * (tmpd*tmpd -1.0);

                        tie = (((-tie / ind) / (ind-1)) + (ind+1) ) / 3.0;
                        for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                            erffact = (((double)(sizes[coor3[0]] - zeros[coor3[0]]) ) * (row_sizes[coor3[2]] + zeros[coor3[0]] - sizes[coor3[0]] - row_zeros)) * tie;
                            mu =  ((double)zerosum[coor3[1]] + sum[coor3[1]]) - row_zeros * sizes[coor3[0]];
                            fout[i].data[coor3[0]] = -mu / sqrt(erffact);
                        }

                    }*/
                }
                if (opt_logit_auroc){
                    for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                        tmpd = zerosum[coor3[1]] + sum[coor3[1]];
                        mu = ((double)sizes[coor3[0]]) * (row_sizes[coor3[2]] - sizes[coor3[0]]);
                        (*opt_logit_auroc)[i].data[coor3[0]] = log(mu + tmpd) - log(mu - tmpd);
                    }
                }
                if (opt_drop_enrich){
                    for(coor3[1] = 0, coor3[0] = coor3[2] * list_grid.sizes[0]; coor3[1] < list_grid.sizes[0]; coor3[1]++,coor3[0]++){
                        if ((sizes[coor3[0]] - zeros[coor3[0]]) != 0)(*opt_drop_enrich)[i].data[coor3[0]] = (((double)(sizes[coor3[0]] - zeros[coor3[0]])) * (row_sizes[coor3[2]])) / (((double)(row_sizes[coor3[2]] - row_zeros)) * (sizes[coor3[0]]));
                        else (*opt_drop_enrich)[i].data[coor3[0]] = ((sizes[coor3[0]] == 0)||((row_sizes[coor3[2]] - row_zeros) == 0)) ? 1.0 : 0.0f;
                }   }
        }  }
        tb.finishProgress(threadID);
    return 0;}
    };
    Task datast(tb,*this,list_grid,opt_logit_auroc,opt_mean,opt_drop_enrich);

    if (downsample != NULL){
        // preparing downsampling scope
        datast.downsampling = new uint32_t[16 * list_grid.sizes[0] * list_grid.sizes[1]]; LFH_NICE_ALLOCERROR(datast.downsampling,"")

        for(int i=0;i<list_grid.sizes[0] * list_grid.sizes[1];i++){
            datast.downsampling[(i << 4)] = (uint32_t)( ldexp((*downsample).data[i], 32)); // if rand < val => retained
            if ((datast.downsampling[(i << 4)] == 0)&&((*downsample).data[i]>0.5)){
                datast.downsampling[(i << 4)] = 0xFFFFFFFF;
            }else{
            datast.downsampling[(i << 4) | 1] = (uint32_t)( ldexp( pow(1.0-(*downsample).data[i], 15.0), 32));  // 15! / x! * (15-x)!
            datast.downsampling[(i << 4) | 2] = datast.downsampling[(i << 4) | 1] + (uint32_t)( ldexp( 15.0 * (*downsample).data[i]*pow(1.0-(*downsample).data[i], 14.0), 32));
            datast.downsampling[(i << 4) | 3] = datast.downsampling[(i << 4) | 2] + (uint32_t)( ldexp( 105.0 * pow((*downsample).data[i], 2.0)*pow(1.0-(*downsample).data[i], 13.0), 32));
            datast.downsampling[(i << 4) | 4] = datast.downsampling[(i << 4) | 3] + (uint32_t)( ldexp( 455.0 * pow((*downsample).data[i], 3.0)*pow(1.0-(*downsample).data[i], 12.0), 32));
            datast.downsampling[(i << 4) | 5] = datast.downsampling[(i << 4) | 4] + (uint32_t)( ldexp( 1365.0 * pow((*downsample).data[i], 4.0)*pow(1.0-(*downsample).data[i], 11.0), 32));
            datast.downsampling[(i << 4) | 6] = datast.downsampling[(i << 4) | 5] + (uint32_t)( ldexp( 3003.0 * pow((*downsample).data[i], 5.0)*pow(1.0-(*downsample).data[i], 10.0), 32));
            datast.downsampling[(i << 4) | 7] = datast.downsampling[(i << 4) | 6] + (uint32_t)( ldexp( 5005.0 * pow((*downsample).data[i], 6.0)*pow(1.0-(*downsample).data[i], 9.0), 32));
            datast.downsampling[(i << 4) | 8] = datast.downsampling[(i << 4) | 7] + (uint32_t)( ldexp( 6435.0 * pow((*downsample).data[i], 7.0)*pow(1.0-(*downsample).data[i], 8.0), 32));
            datast.downsampling[(i << 4) | 15] = (uint32_t)( ldexp(1.0 - pow((*downsample).data[i], 15.0), 32));
            datast.downsampling[(i << 4) | 14] =  datast.downsampling[(i << 4) | 15] -  (uint32_t)(  ldexp( 15.0 * (1.0-(*downsample).data[i]) *pow((*downsample).data[i], 14.0), 32));
            datast.downsampling[(i << 4) | 13] =  datast.downsampling[(i << 4) | 14] -  (uint32_t)(  ldexp( 105.0 * pow(1.0-(*downsample).data[i],2.0) *pow((*downsample).data[i], 13.0), 32));
            datast.downsampling[(i << 4) | 12] =  datast.downsampling[(i << 4) | 13] -  (uint32_t)(  ldexp( 455.0 * pow(1.0-(*downsample).data[i],3.0) *pow((*downsample).data[i], 12.0), 32));
            datast.downsampling[(i << 4) | 11] =  datast.downsampling[(i << 4) | 12] -  (uint32_t)(  ldexp( 1365.0 * pow(1.0-(*downsample).data[i],4.0) *pow((*downsample).data[i], 11.0), 32));
            datast.downsampling[(i << 4) | 10] =  datast.downsampling[(i << 4) | 11] -  (uint32_t)(  ldexp( 3003.0 * pow(1.0-(*downsample).data[i],5.0) *pow((*downsample).data[i], 10.0), 32));
            datast.downsampling[(i << 4) | 9 ] =  datast.downsampling[(i << 4) | 10] -  (uint32_t)(  ldexp( 5005.0 * pow(1.0-(*downsample).data[i],6.0) *pow((*downsample).data[i], 9.0), 32));
          //  printf("%f -> mid %i vs %i\n", (*downsample).data[i], datast.downsampling[(i << 4) | 8], datast.downsampling[(i << 4) | 9] -  (uint32_t)(  ldexp( 6435.0 * pow((*downsample).data[i],7.0) *pow(1.0-(*downsample).data[i], 8.0), 32)));
          //  for(int j=0;j<16;j++) printf("%i: %f (%i)\n", j, (pow(0.5, 32.0) * datast.downsampling[(i << 4)|j]), datast.downsampling[(i << 4)|j]);
            }
        }
    }

    tb.initEqualRanges(datast.colrange, data.getSize());
    datast.sizes.setSize(list_grid.sizes[0]*list_grid.sizes[1]);
    datast.row_sizes.setSize(list_grid.sizes[1]).toZero();
    datast.fout.setSize(data.getSize());
    if (opt_logit_auroc) opt_logit_auroc->setSize(data.getSize());
	if (opt_mean) opt_mean->setSize(data.getSize());
	if (opt_drop_enrich) opt_drop_enrich->setSize(data.getSize());
    Tuple<uint32_t, 2u> coor2;
    Tuple<uint32_t, 3u> coor3;
    for(coor2[1] = 0 ; coor2[1]<list_grid.sizes[1]; coor2[1]++){
        for(coor2[0] = 0 ; coor2[0]<list_grid.sizes[0]; coor2[0]++){
            datast.totallength += list_grid(coor2).getSize();
            coor3[1] = coor2[0] + coor2[1] * list_grid.sizes[0];
            datast.sizes[coor3[1]] = list_grid(coor2).getSize();
            datast.row_sizes[coor2[1]] += list_grid(coor2).getSize();
            for(coor3[0] =0;coor3[0] < list_grid(coor2).getSize(); coor3[0]++){
                datast.partit[list_grid(coor2)[coor3[0]]] = coor3[1];
            }
        }
    }
    if (tb.flags.find("print.progress") != 0xFFFFFFFF) tb.startProgress("Performing Hypergeo-U test", data.getSize());
    tb.runTask(&datast);

   // double tmp;
    if (tb.flags.find("return_logitpval") == 0xFFFFFFFF){
        datast.totallength = datast.fout[0].sizes[1] * datast.fout[0].sizes[0];
        for(uint32_t j=0;j<datast.fout.getSize();j++) for(uint32_t i=0;i < datast.totallength;i++) {
                //tmp = fout[j].data[i];
                datast.fout[j].data[i] = logitPval_to_Zscore(datast.fout[j].data[i]);
               // if ((ExOp::isValid(tmp))&&(!ExOp::isValid(fout[j].data[i]))) printf("failed to convert to Zscore! %e -> %e\n", tmp,fout[j].data[i]);
        }
    }
return datast.fout;}

LFHTEMP template<class F> SparseMatrix<C> SparseMatrix<C>::subsetRows(F& f, ITERABLE_DEF(F) ) const{SparseMatrix<C> fout;
    myHashmap<uint32_t> coorhash; coorhash.toConvert(f);
    fout.setNBcols(this->getNBcols());
    for(uint32_t j=0;j<this->getNBcols();j++){
        for(uint32_t i =0; i< data[j].getSize();i++){
            if (coorhash.find( data[j].deref_key(i)) != 0xFFFFFFFF) fout.data[j][data[j].deref_key(i)] = data[j].deref(i);
        }
    }
return(fout);}
LFHTEMP DataGrid<C,2u> SparseMatrix<C>::mkDataGrid() const{ DataGrid<C,2u> fout;
    Tuple<uint32_t, 2u> coor; coor[1] = data.getSize();
    Vector<uint32_t> rows = this->compute_rowsizes();
    coor[0] = rows.getSize();
    /*
    if (coor[0] == 0u){
        for(coor[1]=0;coor[1]< data.getSize();coor[1]++){
            for(uint32_t j=0;j< data[coor[1]].getSize();j++){
                if (coor[0]<= data[coor[1]].deref_key(j)) coor[0] = data[coor[1]].deref_key(j)+1;
            }
        }
    }*/
    fout.toSizes(coor).toZero();
    for(coor[1]=0;coor[1]< data.getSize();coor[1]++){
        for(uint32_t j=0;j< data[coor[1]].getSize();j++){
            coor[0] = data[coor[1]].deref_key(j);
            fout(coor) = data[coor[1]].deref(j);
        }
    }
return fout;}


LFHTEMP uint32_t SparseMatrix<C>::computeNBrows() const { uint32_t fout =0u;
    unsigned int i,j;
    for(i=0;i< data.getSize();i++){
        for(j=0;j< data[i].getSize();j++){
            if (data[i].deref_key(j) >= fout) fout = data[i].deref_key(j) + 1u;
        }
    }
return fout;}

LFHTEMP Vector<uint32_t> SparseMatrix<C>::compute_rowsizes(uint32_t nbrow_guess) const { Vector<uint32_t> fout;
    unsigned int i,j;
    fout.setSize(nbrow_guess);
    if (nbrow_guess) fout.toZero();
    for(i=0;i< data.getSize();i++){
        for(j=0;j< data[i].getSize();j++){
            while(data[i].deref_key(j) >= fout.getSize()) fout.push_back(0);
            fout[data[i].deref_key(j)]++;
        }
    }
return fout;}
LFHTEMP void SparseMatrix<C>::memmoveColumn(uint32_t col, SparseTuple<C>& newcol){data[col].toMemmove( (myHashmap<uint32_t,C>&) newcol);}
LFHTEMP SparseMatrix<C>& SparseMatrix<C>::toTranspose(){ SparseMatrix<C> tmp = this->mkTranspose(); return this->toMemmove(tmp);}
LFHTEMP SparseMatrix<C> SparseMatrix<C>::mkTranspose() const{ SparseMatrix<C> fout;
    uint32_t i,j;
    for(j=0;j<data.getSize();j++){
        for(i=0;i<data[j].getSize();i++){
            if (data[j].deref_key(i) >= fout.data.getSize()) fout.data.upSize(data[j].deref_key(i)+1);
            fout.data[data[j].deref_key(i)][j] = data[j].deref(i);
        }
    }
return fout;}
LFHTEMP SparseMatrix<C>& SparseMatrix<C>::toFilter(double fraq_of_row, double fraq_of_col){
    Vector<KeyElem<uint32_t, uint32_t> > row_freq_rank;
    Vector<KeyElem<uint32_t, uint32_t> > col_freq_rank;
    unsigned int i;
    Vector<uint32_t> row_sizes = this->compute_rowsizes();
    for(i=0;i<row_sizes.getSize();i++) row_freq_rank.push_back(KeyElem<uint32_t, uint32_t>(row_sizes[i],i));
    for(i=0;i<data.getSize();i++) col_freq_rank.push_back(KeyElem<uint32_t, uint32_t>(data[i].getSize(),i));
    row_freq_rank.sort();
    col_freq_rank.sort();
    Vector<unsigned int> row_torem; row_torem.setSize((uint32_t)((1.0 - fraq_of_row) * row_sizes.getSize()) );
    Vector<unsigned int> col_torem; col_torem.setSize((uint32_t)((1.0 - fraq_of_col) * data.getSize()) );
    for(i=0;i<row_torem.getSize();i++) row_torem[i] = row_freq_rank[i].d;
    for(i=0;i<col_torem.getSize();i++) col_torem[i] = col_freq_rank[i].d;
    row_torem.sort();
    col_torem.sort();
    for(i=col_torem.getSize()-1;i != ExCo<uint32_t>::mkMaximum() ;i--) data.pop_swap(col_torem[i]);
    uint32_t j,k,l;
    for(i=row_torem.getSize()-1;i != ExCo<uint32_t>::mkMaximum() ;i--){
        l = row_sizes.getSize() - row_torem.getSize() + i;
        for(j=0;j< data.getSize();j++){
            data[j].erase_from_iterator(data[j].find(row_torem[i]));
            k = data[j].find(l);
            if (k != 0xFFFFFFFF) data[j].changeKey(k, row_torem[i]);

        }
    }
return *this;}

LFHTEMP SparseMatrix<C>& SparseMatrix<C>::toNormalizedRows(){
    myHashmap<uint32_t, Tuple<C , 2> > normbuf;
    unsigned int i,j,k;
    Tuple<C , 2> input;
    for(i=0;i<data.getSize();i++){
        for(j=0;j< data[i].getSize();j++){
            input[0] = data[i].deref(j);
            input[1] = input[0] * input[0];
            if ((k = normbuf.find(data[i].deref_key(j))) == 0xFFFFFFFF) normbuf = input;
            else normbuf.deref(k) += input;
        }
    }
    Vector<uint32_t> row_sizes = compute_rowsizes();
    for(j=0;j< normbuf.getSize();j++){
        input[0] = normbuf.deref(j)[0] / row_sizes[normbuf.deref_key(j)];
        input[1] = normbuf.deref(j)[1] / row_sizes[normbuf.deref_key(j)];
        input[1] -= input[0] * input[0];
        input[1] = pow(input[1], -0.5f);
        normbuf.deref(j)[1] = (ExOp::isValid(input[1])) ? input[1] : 1.0f;
        normbuf.deref(j)[0] = (ExOp::isValid(input[0])) ? input[0] : 0.0f;
    }
    for(i=0;i<data.getSize();i++){
        for(j=0;j< data[i].getSize();j++){
            data[i].deref(j) -= normbuf[data[i].deref_key(j)][0];
            data[i].deref(j) /= normbuf[data[i].deref_key(j)][1];
        }
    }
    return *this;
}

LFHTEMP void SparseMatrix<C>::wrDotProductKernel(Trianglix<C>& fout)const{
    fout.setSize(data.getSize()); fout.toZero();
    uint32_t i,j,k;
    uint32_t l,m;
    C tmp;
    for(i=0,k=0;i<data.getSize();i++,k++){
        for(j=0;j<i;j++,k++){
            for(l=0;l< data[i].getSize();l++){
                if (m == data[j].find(data[i].deref_key(l))) fout.data[k] += data[j].deref[m] * data[i].deref[l];
            }
        }
        for(l=0;l< data[i].getSize();l++) fout.data[k] += data[i].deref[l] * data[i].deref[l];
    }
}

LFHTEMP void SparseMatrix<C>::wrRowDistance(Trianglix<double> &fout)const{
    fout.setSize(data.getSize()); fout.toZero();
    uint32_t i,j,k;
    uint32_t l,m;
    C tmp;
    for(i=0,k=0;i<data.getSize();i++,k++){
        for(j=0;j<i;j++,k++){
            for(l=0;l< data[i].getSize();l++){
                if (m == data[j].find(data[i].deref_key(l))){
                    tmp = data[j].deref[m] - data[i].deref[l];
                    fout.data[k] += ExOp::norm(tmp);
                }
            }
        }
    }
return;}
LFHTEMP C& SparseMatrix<C>::addEntry(uint32_t row, uint32_t col){return data[col].addEntry(row);}
LFHTEMP C& SparseMatrix<C>::addEntry(Tuple<uint32_t, 2u> coor){return data[coor[1]].addEntry(coor[0]);}
LFHTEMP C SparseMatrix<C>::operator()(unsigned int x, unsigned int y) const{
    uint32_t ite;
    if ((ite = data[y].find(x)) != 0xFFFFFFFF) return data[y].deref(ite);
    C fout;
return ExOp::toZero(fout);}
LFHTEMP C& SparseMatrix<C>::operator()(unsigned int x, unsigned int y){
    uint32_t ite;
    if ((ite = data[y].find(x)) != 0xFFFFFFFF) return data[y].deref(ite);
return data[y].addEntry(x);}

LFHTEMP SparseMatrix<C>& SparseMatrix<C>::permuteRows(const Tuple<uint32_t> &permutation){
    for(int i=0;i < data.getSize();i++){
        for(int j=0;j<data[i].getSize();j++) data[i].deref_key(j) = permutation[data[i].deref_key(j)];
        data[i].rehash(data[i].hash_mag);
    }
return *this;}
LFHTEMP SparseMatrix<C>& SparseMatrix<C>::operator=(const SparseMatrix<C>& other){data = other.data;return *this;}
LFHTEMP template<class D> SparseMatrix<C>& SparseMatrix<C>::operator=(const SparseMatrix<D>& other){data = other.data;return *this;}

LFHTEMP template<class D> SparseMatrix<C>& SparseMatrix<C>::operator+=(SparseMatrix<D> const & other){ // TODO, maybe
    uint32_t inds[2];
    for(uint32_t col;col<data.getSize();col++){
        if (other.data[col].getSize() == 0) continue;
        for(inds[1]= 0;inds[1] < other.data[col].getSize() ; inds[1]++){
            if ((inds[0] = data[col].find(other.data[col].deref_key(inds[1]))) != 0xFFFFFFFF)  data[col].deref(inds[0]) += other.data[col].deref(inds[1]);
            else this->addEntry(col,other.data[col].deref_key(inds[1])) = other.data[col].deref(inds[1]);
        }
    }
return *this;}
LFHTEMP template<class D> SparseMatrix<C>& SparseMatrix<C>::operator-=(SparseMatrix<D> const & other){ // TODO, maybe
    uint32_t inds[2];
    for(uint32_t col;col<data.getSize();col++){
        if (other.data[col].getSize() == 0) continue;
        for(inds[1]= 0;inds[1] < other.data[col].getSize() ; inds[1]++){
            if ((inds[0] = data[col].find(other.data[col].deref_key(inds[1]))) != 0xFFFFFFFF)  data[col].deref(inds[0]) -= other.data[col].deref(inds[1]);
            else this->addEntry(col,other.data[col].deref_key(inds[1])) = -other.data[col].deref(inds[1]);
        }
    }
return *this;}

LFHTEMP template<class D> auto SparseMatrix<C>::operator*(SparseTuple<D> const & other) const
 -> SparseTuple<decltype( (*this)(0u,0u) * other[0u] )>{
    SparseTuple<decltype( (*this)(0u,0u) * other[0u] )> fout;
return fout;}
LFHTEMP SparseMatrix<C>& SparseMatrix<C>::toColAppend(const SparseMatrix<C>& other){data.toAppend(other.data);return(*this);}
LFHTEMP SparseMatrix<C>& SparseMatrix<C>::toColMemappend(SparseMatrix<C>& other){data.toMemappend(other.data);return(*this);}
LFHTEMP SparseMatrix<C>& SparseMatrix<C>::toSparseFromNonEqual(const DataGrid<C, 2u>& data, C nan_value){
    this->setSizes(data.dims[0], data.dims[1]);
    Tuple<uint32_t, 2u> coor;
    for(coor[1] =0;coor[1]<data.dims[1];coor[1]++){
        for(coor[0] =0;coor[0]<data.dims[0];coor[0]++){
            if (data(coor) != nan_value) this->addEntry(coor[0],coor[1]) = data(coor);
        }
    }
return *this;}

LFHTEMP SparseMatrix<C>::operator TMatrix<C>()const{  TMatrix<C> fout;
    uint32_t nbrow = this->computeNBrows();
    fout.setSizes(nbrow, data.getSize());
    fout.toZero();
    for(uint32_t j=0; j < data.getSize();j++){
        for(uint32_t i=0;i< data[j].getSize();i++) fout.data[data[j].deref_key(i) + j * nbrow] = data[j].deref(i);
    }
return fout;}

LFHTEMP void SparseMatrix<C>::sortByRowNumber(){for(uint32_t i=0;i<data.getSize();i++) data[i].keysort();}
LFHTEMP void SparseMatrix<C>::sortByRowNumber(uint32_t colID){ data[colID].keysort();}
LFHTEMP bool SparseMatrix<C>::makeMaximumAsFirst(uint32_t colID){
    if (data[colID].getSize() <= 1) return false;
    uint32_t maxcur = 0;
    uint32_t i;
    for(i=1;i< data[colID].getSize();i++){
        if (data[colID].deref(i) > data[colID].deref(maxcur)) maxcur =i;
    }
    if (maxcur ==0) return false;
    data[colID].swapEntries(0,maxcur);
return true;}

LFHTEMP bool SparsityComparator<C>::operator<(const SparsityComparator<C>& other) const{
    uint32_t ti;
    for(ti=0;ti < target->getSize();ti++){
        if (ti >= other.target->getSize()) return false;
        if (target->deref_key(ti) < other.target->deref_key(ti)) return true;
        else if (target->deref_key(ti) > other.target->deref_key(ti)) return false;
        ti++;
    }
return (ti < other.target->getSize());}
LFHTEMP bool SparsityComparator<C>::operator>(const SparsityComparator<C>& other) const{
    uint32_t ti;
    for(ti=0;ti < target->getSize();ti++){
        if (ti >= other.target->getSize()) return true;
        if (target->deref_key(ti) < other.target->deref_key(ti)) return false;
        else if (target->deref_key(ti) > other.target->deref_key(ti)) return true;
        ti++;
    }
return false;}
LFHTEMP bool SparsityComparator<C>::operator==(const SparsityComparator<C>& other) const{
    uint32_t ti;
    if (target->getSize() != other.target->getSize()) return false;
    for(ti=0;ti < target->getSize();ti++){
        if (ti >= other.target->getSize()) return false;
        if (target->deref_key(ti) < other.target->deref_key(ti)) return false;
        else if (target->deref_key(ti) > other.target->deref_key(ti)) return false;
    }
return true;}
LFHTEMP bool SparsityComparator<C>::operator>=(const SparsityComparator<C>& other) const{return !((*this) < other);}
LFHTEMP bool SparsityComparator<C>::operator<=(const SparsityComparator<C>& other) const{return !((*this) > other);}
LFHTEMP bool SparsityComparator<C>::operator!=(const SparsityComparator<C>& other) const{return !((*this) == other);}

LFHTEMP Vector<SparsityComparator<C> > SparseMatrix<C>::getSparseOrdered()const{Vector<SparsityComparator<C> > fout;
    fout.setSize(data.getSize());
    for(uint32_t i=0;i< data.getSize();i++) fout[i].target = &(data[i]);
    fout.sort();
return fout;}

/*
LFHTEMP SetPartition<uint32_t,uint32_t> SparseMatrix<C>::groupSparsity() const{ SetPartition<uint32_t,uint32_t> fout;

return fout;}*/

LFHTEMP ERRCODE SparseMatrix<C>::readDenseTable(const char* path_input, Vector<string> &rownames, Vector<string> &colnames, string &tname, OptionStruct optstr, bool is_floatingpoint, char separator){
    FILE *f = fopen(path_input, "r+");
    if (f == NULL) {fprintf(stderr,"Could not open '%s'!\n", path_input); return ERRCODE_NOFILE;}
    char buffer[256];
    char argbuffer[256];
    char argbuffer2[256];
    char sep[4];
    union{
        double dval;
        uint32_t ival;
    };
    rownames.toMemfree();
    colnames.toMemfree();
    sprintf(argbuffer,"%\[^%c\n\r]\%\[%c\n\r]",separator,separator);
    int tmp;
    if (2 != fscanf(f, argbuffer, buffer, sep)) {fprintf(stderr,"file '%s' format is not recognized!\n", path_input); fclose(f); return 1;}
    tname = string(buffer);
    if ((buffer[0] == '"')&&(buffer[tmp = (strlen(buffer)-1)] == '"')) { buffer[tmp] = '\0'; tname = string(buffer+1);}
    else tname = string(buffer);

//    int toomuch=0;
    while((sep[0] != '\n')&&(sep[0] != '\r')){
        if (2 != fscanf(f, argbuffer, buffer, sep)) {fprintf(stderr,"file '%s' format is not recognized!\n", path_input); fclose(f); return 1;}
//        printf("%i: %i and %s\n", toomuch, (int)sep[0], buffer); fflush(stdout);
        if ((buffer[0] == '"')&&(buffer[tmp = (strlen(buffer)-1)] == '"')) { buffer[tmp] = '\0'; colnames.push_back(string(buffer+1));}
        else colnames.push_back(string(buffer));
//        if (toomuch++ > 17000) {fclose(f); return 1;}
    }
    this->setNBcols(colnames.getSize());
    uint32_t r,c;
    r=0;
    if (is_floatingpoint){
        sprintf(argbuffer2,"%%lf\%\[%c\n\r]",separator);
        while(!feof(f)){
            if (2 != fscanf(f, argbuffer, buffer, sep)) {fprintf(stderr,"file '%s' format is not recognized!\n", path_input); fclose(f); return 1;}
            if ((buffer[0] == '"')&&(buffer[tmp = (strlen(buffer)-1)] == '"')) { buffer[tmp] = '\0'; rownames.push_back(string(buffer+1));}
            else rownames.push_back(string(buffer));

            c=0;
            while((sep[0] != '\n')&&(sep[0] != '\r')){
                fscanf(f,argbuffer2, &dval, &sep);
                if (dval != 0.0) data[c][r] = dval;
                c++;
            }
            if (c != colnames.getSize()){fprintf(stderr,"file '%s' format corrupted (%i columns < %i expected)!\n", path_input, c+1, colnames.getSize()+1); fclose(f); return 1;}
            r++;
        }
    }else{
        sprintf(argbuffer2,"%%i\%\[%c\n\r]",separator);
        //printf("%s da parse\n", argbuffer2); printf("%s\n",argbuffer2); fflush(stdout);
        while(!feof(f)){
            if (2 != fscanf(f, argbuffer, buffer, sep)) {fprintf(stderr,"file '%s' format is not recognized!\n", path_input); fclose(f); return 1;}
            if ((buffer[0] == '"')&&(buffer[tmp = (strlen(buffer)-1)] == '"')) { buffer[tmp] = '\0'; rownames.push_back(string(buffer+1));}
            else rownames.push_back(string(buffer));
            c=0;
            while((sep[0] != '\n')&&(sep[0] != '\r')){
                fscanf(f,argbuffer2, &ival, sep);
                if (dval != 0.0) data[c][r] = ival;
                c++;
            }
            if (c != colnames.getSize()){fprintf(stderr,"file '%s' format corrupted (%i columns < %i expected)!\n", path_input, c+1, colnames.getSize()+1); fclose(f); return 1;}
            r++;
        }
    }
    fclose(f);
return 0;}
LFHTEMP ERRCODE SparseMatrix<C>::readMTXTable(const char* path_prefix, Vector<string> &rownames, Vector<string> &colnames, string &tname, bool is_floatingpoint, bool is_min_addr_one){
    char buffer[1024];
    uint32_t plen = strlen(path_prefix);
    memcpy(buffer, path_prefix, plen); memcpy(buffer + plen, "genes.tsv", 10);
    rownames.toMemfree();
    FILE *f = fopen(buffer, "r+");
    if (f == NULL) {fprintf(stderr,"Could not read file '%s'!\n", buffer); return ERRCODE_NOFILE;}
    while(1 == fscanf(f,"%[^\n]\n",buffer)) rownames.push_back(string(buffer));
    fclose(f);
    colnames.toMemfree();
    memcpy(buffer, path_prefix, plen); memcpy(buffer + plen, "barcodes.tsv", 13);
    f = fopen(buffer, "r+");
    if (f == NULL) {fprintf(stderr,"Could not read file '%s'!\n", buffer); return ERRCODE_NOFILE;}
    while(1 == fscanf(f,"%[^\n]\n",buffer)) colnames.push_back(string(buffer));
    fclose(f);

    memcpy(buffer, path_prefix, plen); memcpy(buffer + plen, "matrix.mtx", 11);
    f = fopen(buffer, "r+");
    if (f == NULL) {fprintf(stderr, "could not open '%s'!\n", buffer); return ERRCODE_NOFILE;}
    tname = string("");
    while(1 == fscanf(f," %1[%]", buffer)) {
        if (tname.length() == 0) {fscanf(f,"%[^\n]\n",buffer); tname = string(buffer);} // first comment line detected
        else fscanf(f,"%*[^\n]\n"); // comment line detected
    }
    uint32_t buf[3];
    double dbuf;
    if (3 != fscanf(f,"%i %i %i\n", buf,buf+1,buf+2)) {fprintf(stderr, "Failed to read Header of file '%s' ", buffer); return 1;}
    if (buf[0] != rownames.getSize()) {fprintf(stderr, "Incorrect number of rows %i != %i expected from genes.tsv", buf[0], rownames.getSize()); return 1;}
    if (buf[1] != colnames.getSize()) {fprintf(stderr, "Incorrect number of rows %i != %i expected from barcodes.tsv",buf[1], colnames.getSize()); return 1;}

    this->data.setSize(buf[1]);
    if (is_floatingpoint){
        if (is_min_addr_one) while(3 == fscanf(f,"%i %i %lf\n", buf,buf+1,&dbuf)) this->data[buf[1]-1][buf[0]-1] = dbuf;
        else while(3 == fscanf(f,"%i %i %lf\n", buf,buf+1,&dbuf)) this->data[buf[1]][buf[0]] = dbuf;
    }else{
        if (is_min_addr_one) while(3 == fscanf(f,"%i %i %i\n", buf,buf+1,buf+2)) this->data[buf[1]-1][buf[0]-1] = buf[2];
        else while(3 == fscanf(f,"%i %i %i\n", buf,buf+1,buf+2)) this->data[buf[1]][buf[0]] = buf[2];
    }
return 0;}


LFHTEMP ERRCODE SparseMatrix<C>::writeMTXTable(const char* path_prefix, const Vector<string> &rownames, const Vector<string> &colnames, const string &tname, OptionStruct optstr, bool is_floatingpoint, bool is_min_addr_one)const{
    char buffer[1024];
    uint32_t plen = strlen(path_prefix);
    memcpy(buffer, path_prefix, plen);
    memcpy(buffer + plen, "genes.tsv", 10);
    FILE *f = fopen(buffer, "w+");
    if (f == NULL) {fprintf(stderr,"Could not create file '%s'!\n", buffer); return ERRCODE_NOFILE;}
    for(uint32_t i =0;i<rownames.getSize();i++) fprintf(f, "%s\n", rownames[i].c_str());
    fclose(f);
    memcpy(buffer + plen, "barcodes.tsv", 13);
    f = fopen(buffer, "w+");
    if (f == NULL) {fprintf(stderr,"Could not create file '%s'!\n", buffer); return ERRCODE_NOFILE;}
    for(uint32_t i =0;i<colnames.getSize();i++) fprintf(f, "%s\n", colnames[i].c_str());
    fclose(f);
    memcpy(buffer + plen, "matrix.mtx", 11);
    f = fopen(buffer, "w+");
    if (f == NULL) {fprintf(stderr,"Could not create file '%s'!\n", buffer); return ERRCODE_NOFILE;}
    fprintf(f, "% Table %s\n", tname.c_str());
    uint32_t i,j;
    for(i=0,j=0; i< data.getSize();i++) j += data[i].getSize();
    fprintf(f,"%i %i %i\n", rownames.getSize(), colnames.getSize(), j);
    if (is_floatingpoint){
        if (is_min_addr_one){
            for(i=0;i< data.getSize();i++){
                for(j=0;j< data[i].getSize();j++) fprintf(f, "%i %i %lf\n",data[i].deref_key(j)+1u ,i+1u, (double)data[i].deref(j));
            }
        }else{
            for(i=0;i< data.getSize();i++){
                for(j=0;j< data[i].getSize();j++) fprintf(f, "%i %i %lf\n",data[i].deref_key(j) ,i, (double)data[i].deref(j));
            }
        }
    }else{
        if (is_min_addr_one){
            for(i=0;i< data.getSize();i++){
                for(j=0;j< data[i].getSize();j++) fprintf(f, "%i %i %i\n",data[i].deref_key(j)+1u,i+1u, (int)data[i].deref(j));
            }
        }else{
            for(i=0;i< data.getSize();i++){
                for(j=0;j< data[i].getSize();j++) fprintf(f, "%i %i %i\n",data[i].deref_key(j),i, (int)data[i].deref(j));
            }
        }
    }
    fclose(f);

return 0;}

LFHTEMP ERRCODE SparseMatrix<C>::save(FILE* f) const{return data.save(f);}
LFHTEMP ERRCODE SparseMatrix<C>::load(FILE* f, unsigned int length){return data.load(f);}
LFHTEMP void SparseMatrix<C>::show(FILE* f, int level) const{
    Tuple<uint32_t, 2u> coor;
    uint32_t ite;
    C zero; ExOp::toZero(zero);
    /*Annotation* daannot;
    if (((ite = lfhstatic_ressources.find((uint32_t)LFHSTATIC_RESSOURCE_ANNOTATIONS)) !=0xFFFFFFFF)&&(( daannot = (Annotation*)lfhstatic_ressources.deref(ite)) != NULL)&&((ite = daannot->axes.find((void*)this))!=0xFFFFFFFF)){
        uint32_t dico[2];
        dico[0] = daannot->axes.deref(ite)[0];
        dico[1] = daannot->axes.deref(ite)[1];
        fprintf(f, "Sparse matrix with size %i,%i", row_sizes.getSize(), data.getSize());
        if (level == 0){
            fprintf(f, "\n");
            for(coor[0] = 0; coor[0] < data.getSize(); coor[0]++) fprintf(f, "\t%s",  daannot->dicos[dico[0]].entries[coor[0]]);
            for(coor[1] = 0; coor[1] < row_sizes.getSize(); coor[1]++){
                fprintf(f,"\n%s", daannot->dicos[dico[1]].entries[coor[1]]);
                for(coor[0] = 0;coor[0] < data.getSize();coor[0]++){
                    fprintf(f,"\t");
                    ite = data[coor[0]].find(coor[1]);
                    ExOp::show((ite == 0xFFFFFFFF) ? zero : data[coor[0]].deref(ite),f,2);
                }
            }
            fprintf(f, "\n");
        }
    }else{*/
       // fprintf(f, "Sparse matrix with size %i,%i", row_sizes.getSize(), data.getSize());
        /*if (level == 0){
            fprintf(f, "\n");
            for(coor[1] = 0; coor[1] < row_sizes.getSize(); coor[1]++){
                coor[0] = 0;
                while(true){
                    ite = data[coor[0]].find(coor[1]);
                    ExOp::show((ite == 0xFFFFFFFF) ? zero : data[coor[0]].deref(ite),f,2);
                    coor[0]++;
                    if (coor[0] == data.getSize()) break;

                    fprintf(f,"\t");
                }
                fprintf(f,"\n");
            }
        }*/
    //}
}

/*
LFHTEMP Tuple<uint32_t, 2u>  SparseMatrix<C>::KeyIterator::operator()()const{Tuple<uint32_t, 2u> fout;
    fout[0] = getRow();
    fout[1] = colID;
return fout;}*/
/*LFHTEMP myHashmap<uint32_t, uint32_t> SparseMatrix<C>::makeSparseClassPartitions()const{ myHashmap<uint32_t, uint32_t> fout;
    uint32_t i,k;
    // initialize partitions based of total number of entries!
    myHashmap<uint32_t, Tuple<uint32_t, 2u> > startparts;
    Tuple<uint32_t, 2u> startparts_dainput;
    for(i=0;i<data.getSize();i++){
        if ((k=startparts.find(data[i].getSize())) == 0xFFFFFFFF){
            startparts_dainput[0] = i;
            startparts_dainput[1] = i;
            startparts[data[i].getSize()] = startparts_dainput;
        }else{fout[i] = startparts.deref(k)[0];  startparts.deref(k)[0] = i;}
    }
    if (auto ite = startparts.getIterator()) do{
        fout[startparts.deref(k)[1]] = startparts.deref(k)[0];
    }while(ite++);

return fout;}*/

#ifdef Rcpp_hpp

LFHTEMP void SparseMatrix<C>::rdRcpp(const Rcpp::NumericMatrix &object){
    /*AnnotMaintenance<A> maint;
    Annotation* ann = maint.canAnnot();
    if (ann) ann->rdAxesNames(object, *this);*/
    this->setNBcols(object.cols());
    uint32_t nbr = object.rows();

    uint32_t x,y;
    for(y=0;y< data.getSize();y++){
        for(x=0;x< nbr;x++){
            if (object(x,y) != 0.0f) data[y][x] = object(x,y);
        }
    }
}
LFHTEMP void SparseMatrix<C>::wrRcpp(Rcpp::NumericMatrix &object, uint32_t nbrowsuggest)const{
    /*AnnotMaintenance<A> maint;
    Annotation* ann = maint.canAnnot();

    if (ann) ann->wrAxesNames(object, *this);*/
    uint32_t nbrow = this->computeNBrows();
    if (nbrow < nbrowsuggest) nbrow = nbrowsuggest;

    object = Rcpp::NumericMatrix(nbrow, data.getSize());
    uint32_t x,y;
    for(y=0;y< data.getSize();y++){
        for(x=0;x< data[y].getSize();x++){
            object(data[y].deref_key(x),y) = data[y].deref(x);
        }
    }
}


LFHTEMP void SparseMatrix<C>::rdRcppdgCMatrix(const Rcpp::S4 &mat, bool transpose){
   /* AnnotMaintenance<A> maint;
    Annotation* ann = maint.canAnnot();
    if (ann) ann->rdAxesNames(mat, *this);*/
    const int RTYPE = Rcpp::traits::r_sexptype_traits<C>::rtype;

    Rcpp::IntegerVector dims = mat.slot("Dim");
    if (transpose) this->setNBcols(dims[0]);
    else this->setNBcols(dims[1]);
    Rcpp::IntegerVector i = mat.slot("i");
    Rcpp::IntegerVector p = mat.slot("p");
    Rcpp::Vector<RTYPE> x = mat.slot("x");


    uint32_t ite;
    uint32_t pite =1;
    if (transpose){
        for(ite=0;ite< i.size() ;ite++){
            while((pite < p.size())&&(p[pite] == ite)) pite++;
            this->addEntry(pite-1,i[ite]) = x[ite];
        }
    }else{
        for(ite=0;ite< i.size() ;ite++){
            while((pite < p.size())&&(p[pite] == ite)) pite++;
            this->addEntry(i[ite], pite-1) = x[ite];
        }
    }
}
LFHTEMP void SparseMatrix<C>::wrRcppdgCMatrix(Rcpp::S4 &mat, uint32_t nbrowsuggest)const{ // Alloc with S4 mat("dgCMatrix")
    /*AnnotMaintenance<A> maint;
    Annotation* ann = maint.canAnnot();
    if (ann) ann->wrAxesNames(mat, *this);*/

    uint32_t tsize =0;
    uint32_t ite;
    for(ite=0; ite < data.getSize();ite++) tsize += data[ite].getSize();
    Rcpp::IntegerVector i = Rcpp::IntegerVector(tsize);
    Rcpp::IntegerVector p = Rcpp::IntegerVector(data.getSize()+1);
    Rcpp::IntegerVector dim = Rcpp::IntegerVector(2);
    Rcpp::Vector<Rcpp::traits::r_sexptype_traits<C>::rtype> xx = Rcpp::Vector<Rcpp::traits::r_sexptype_traits<C>::rtype>(tsize);

    uint32_t x,y;
    dim[1] = data.getSize();
    dim[0] = nbrowsuggest;
    Vector<int> daorder;
    for(y=0,ite=0;y< data.getSize();y++){
        p[y] =ite;
        if (data[y].getSize() == 0) continue;
        daorder.setSize(data[y].getSize());
        daorder.toOrderOf(&(data[y].heap[0]));
        if (dim[0] <= data[y].deref_key(daorder.last())) dim[0] = data[y].deref_key(daorder.last()) + 1;
        for(x=0;x< data[y].getSize();x++,ite++){
            i[ite] =  data[y].deref_key(daorder[x]);
            xx[ite] =  data[y].deref(daorder[x]);
        }
    }
    p[y] = ite;
    mat = Rcpp::S4("dgCMatrix");
    mat.slot("i") = i;
    mat.slot("p") = p;
    mat.slot("x") = xx;
    mat.slot("Dim") = dim;
}

LFHTEMP void SparseMatrix<C>::wrRcppRow(uint32_t row, SEXP &object)const{
    Rcpp::Vector<Rcpp::traits::r_sexptype_traits<C>::rtype> xx = Rcpp::Vector<Rcpp::traits::r_sexptype_traits<C>::rtype>(data.getSize());
    uint32_t ite;
    for(uint32_t i=0;i< data.getSize() ;i++) xx[i] = ((ite = data[i].find(row)) == 0xFFFFFFFF) ? 0 : data[i].deref(ite);
    object = xx;
}
LFHTEMP void SparseMatrix<C>::wrRcppCol(uint32_t col, SEXP &object, uint32_t nbrowsuggest)const{
    if (nbrowsuggest==0) nbrowsuggest = this->computeNBrows();
    Rcpp::Vector<Rcpp::traits::r_sexptype_traits<C>::rtype> xx = Rcpp::Vector<Rcpp::traits::r_sexptype_traits<C>::rtype>(nbrowsuggest);
    uint32_t ite;
    for(uint32_t i=0;i< nbrowsuggest;i++) xx[i] = ((ite = data[col].find(i)) == 0xFFFFFFFF) ? 0 : data[col].deref(ite);
    object = xx;
}

#endif
/*
	template<class C, int size>
	C Matrianglix<C,size>::determinant(){
		C _out = data[0];
		for(unsigned int i=1;i<size;i++) _out *= data[i * (size+1)];
		return(_out);
	}

	template<class C, int size>
	Matrianglix<C,size> Matrianglix<C,size>::operator*(const Matrianglix<C, size> & other){
		return(Matrianglix<C,size>(( ((TMatrix<C,size,size>) (*this)) * ((TMatrix<C,size,size>) (other)) )  ));
	}

	template<class C, int size>
	TMatrix<C,size,size> Matrianglix<C,size>::inverse(){
		int i,j,k;

		Matrianglix<C,size> _tmp;

		for(int i=0;i<size;i++) {
			_tmp.data[i* (size+1)] = ExCo<C>::invert(data[i* (size+1)]); // inverse of D

		}

		// back stubstitution!
		// lower
		for(i=1;i<size;i++){
			for(j=0;j<size-i;j++){
				_tmp.data[i*size+j*(size+1)] = -data[i*size+j*(size+1)];
				_tmp.data[i+j*(size+1)] = -data[i+j*(size+1)];
				for(k=1;k<i;k++) {
					_tmp.data[i*size+j*(size+1)] -= data[k+i*size+j*(size+1)] * _tmp.data[(i-k)*size+j*(size+1)];
					_tmp.data[i+j*(size+1)] -= data[k*size+i+j*(size+1)] * _tmp.data[(i-k)+j*(size+1)];
				}
				_tmp.data[i+j*(size+1)] *= _tmp.data[(j+i)* (size+1)]*_tmp.data[j*(size+1)];
			}

		}
		_tmp.perm = -perm;

	//	_tmp.show(stdout);
		TMatrix<C,size,size> _out;


		for(i=0;i<size;i++){
			_out.data[i * (size+1)] = _tmp.data[i* (size+1)];
			for(k=i+1;k<size;k++) _out.data[i * (size+1)] += _tmp.data[k + i* size] * _tmp.data[i + k* size];
			for(j=i+1;j<size;j++){
				_out.data[i+j  * size] = _tmp.data[i + j* size] * _tmp.data[j + j* size];
				_out.data[j+i  * size] = _tmp.data[j + i* size];
				for(k=j+1;k<size;k++) _out.data[i+j  * size] += _tmp.data[i + k* size] * _tmp.data[k + j* size];
				for(k=j+1;k<size;k++) _out.data[j+i  * size] += _tmp.data[k + i* size] * _tmp.data[j + k* size];
			}
		}




		return(_out);
	}*/




#undef LFHTEMP
#define LFHTEMP template <class C>


LFHTEMP	void Complex<C>::zero() {ExOp::toZero((*this)[0]);ExOp::toZero((*this)[1]);}
LFHTEMP	void Complex<C>::one() {ExOp::toOne((*this)[0]);ExOp::toZero((*this)[1]);}
LFHTEMP	void Complex<C>::random() {ExOp::toRand((*this)[0]);ExOp::toRand((*this)[1]);}


	LFHTEMP	C& Complex<C>::operator[](unsigned int ind) {return data[ind];}
//	LFHTEMP	C Complex<C>::operator[](unsigned int ind) const {return data[ind];}
	LFHTEMP	const C& Complex<C>::operator[](unsigned int ind) const {return data[ind];}


	LFHTEMP	SETCMP_enum Complex<C>::setcmp(const Complex<C> &a) const{
		bool v = data[0] < ExOp::mkzero<C>();
		if (v ^ (a.data[0] < ExOp::mkzero<C>())) return v ? SETCMP_LT : SETCMP_GT;
		else {
			double norm = ExOp::pnorm(data[0]) + ExOp::pnorm(data[1]);
			double onorm = ExOp::pnorm(a.data[0]) + ExOp::pnorm(a.data[1]);
			if (norm == onorm){
				return ExOp::setcmp(data[1], a.data[1]);
			} else return ((v)^(norm > onorm)) ? SETCMP_LT : SETCMP_GT;
		}
	}

	LFHTEMP bool Complex<C>::operator>(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
	LFHTEMP bool Complex<C>::operator>=(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_E_MASK) == SETCMP_GE);}
	LFHTEMP bool Complex<C>::operator<(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_T_MASK) == SETCMP_LT);}
	LFHTEMP bool Complex<C>::operator<=(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_E_MASK) == SETCMP_LE);}
	LFHTEMP bool Complex<C>::operator==(const Complex<C> &o)const{return(data[0] == o.data[0]) && (data[1] == o.data[1]);}
	LFHTEMP bool Complex<C>::operator!=(const Complex<C> &o)const{return(data[0] != o.data[0]) || (data[1] != o.data[1]);}




LFHTEMP	template<class A>
Complex<C>& Complex<C>::operator+=(const Complex<A>& other){
	ExOp::toAdd((*this)[0],other[0]);
	ExOp::toAdd((*this)[1],other[1]);
	return(*this);
	}

LFHTEMP	template<class A>
	Complex<C>& Complex<C>::operator-=(const Complex<A>& other){
		ExOp::toSubt((*this)[0],other[0]);
		ExOp::toSubt((*this)[1],other[1]);
		return(*this);
	}

LFHTEMP	template<class A>
Complex<C>& Complex<C>::operator*=(const Complex<A>& other){
	C tmp = ExOp::mkMult((*this)[0], other[1]);
	ExOp::toMult((*this)[0], other[0]);
//	ExOp::show((*this)[0]);
	ExOp::toSubt((*this)[0], ExOp::mkMult((*this)[1],other[1]));
//	ExOp::show((*this)[0]);
	ExOp::toMult((*this)[1], other[0]);
//	ExOp::show((*this)[1]);
	ExOp::toAdd((*this)[1], tmp);
	return(*this);
	}

LFHTEMP	template<class A> Complex<C>& Complex<C>::operator+=(const A& a){(*this)[0] += a;return(*this);}
LFHTEMP	template<class A> Complex<C>& Complex<C>::operator-=(const A& a){(*this)[0] -= a; return(*this);}
LFHTEMP	template<class A> Complex<C>& Complex<C>::operator*=(const A& a){
	(*this)[0] *= a;
	(*this)[1] *= a;
	return(*this);
}
LFHTEMP	template<class A> Complex<C>& Complex<C>::operator/=(const A& a){
	(*this)[0] /= a;
	(*this)[1] /= a;
	return(*this);
}

LFHTEMP	template<class A>
Complex< typename STDRETTYPE2<C,A>::PLUS_TYPE > Complex<C>::operator+(const Complex<A>& other) const{
	return(Complex((*this)[0] + other[0], (*this)[1] + other[1]));
}

LFHTEMP	template<class A>
Complex< typename STDRETTYPE2<C,A>::MINU_TYPE > Complex<C>::operator-(const Complex<A>& other) const{
	return(Complex((*this)[0] - other[0], (*this)[1] - other[1]));
}

LFHTEMP	template<class A> Complex< typename STDRETTYPE2<C,A>::PROD_TYPE > Complex<C>::operator*(Complex<A> const & other) const{
    return Complex< typename STDRETTYPE2<C,A>::PROD_TYPE >( ((*this)[0] *  other[0])-  ((*this)[1] *  other[1]) , ((*this)[0] *  other[1])+  ((*this)[1] *  other[0]));
    }
LFHTEMP	template<class A> Complex< typename STDRETTYPE2<C,A>::DIVI_TYPE > Complex<C>::operator/(Complex<A> const & other) const{

}

LFHTEMP	Complex< C >& Complex<C>::operator*=(Complex<C> const & other){
    C tmp = ((*this)[0] *  other[0])-  ((*this)[1] *  other[1]);
    (*this)[1] = ((*this)[0] *  other[1])+  ((*this)[1] *  other[0]);
    (*this)[0] = tmp;
    return *this;
    }
LFHTEMP	Complex< C >& Complex<C>::operator/=(Complex<C> const & other){
    double norm = ExOp::pnorm(other);
    C tmp = ((*this)[0] *  other[0]) +  ((*this)[1] *  other[1]);
    (*this)[1] = (((*this)[1] *  other[0])- ((*this)[0] *  other[1]))/ norm;
    (*this)[0] = tmp / norm;
    return *this;
}

LFHTEMP	template<class A> Complex< typename STDRETTYPE2<C,A>::PROD_TYPE > Complex<C>::operator*(A const & other) const{
    return Complex< typename STDRETTYPE2<C,A>::PROD_TYPE >((*this)[0] * other, (*this)[1] * other);
    }
LFHTEMP	template<class A> Complex< typename STDRETTYPE2<C,A>::DIVI_TYPE > Complex<C>::operator/(A const & other) const{
    return Complex< typename STDRETTYPE2<C,A>::PROD_TYPE >((*this)[0] / other, (*this)[1] / other);
    }

LFHTEMP Complex< typename ExCo<C>::NEG_TYPE > Complex<C>::operator-() const{
	return(Complex(-(*this)[0], -(*this)[1]));
}

LFHTEMP double Complex<C>::pnorm() const{ return(ExOp::pnorm(((Tuple<C,2>*)this)->data[0]) + ExOp::pnorm(((Tuple<C,2>*)this)->data[1])); }

LFHTEMP Complex<C> Complex<C>::inverse() const{
	Complex<C> _out;
	C tmpA = ((*this)[0] * (*this)[0]) + ((*this)[1] * (*this)[1]);
	_out[0] =(*this)[0] / tmpA;
	_out[1] =(*this)[1] / tmpA;
	return(Complex<C>(_out));
	}

LFHTEMP void Complex<C>::show(FILE* f, int level) const{
	switch(level){
		case 1:
            ExOp::show((*this)[0],f,level+1);fprintf(f," r\t");ExOp::show((*this)[1],f,level+1);fprintf(f," i");
		break;
		case 0:
            ExOp::show((*this)[0],f,level+1);fprintf(f," r\t");ExOp::show((*this)[1],f,level+1);fprintf(f," i\n");
		break;
		case 2:
            fprintf(f,"[");ExOp::show((*this)[0],f,level+1);fprintf(f," r;");ExOp::show((*this)[1],f,level+1);fprintf(f," i]");
			break;
		default:
            fprintf(f,"(");ExOp::show((*this)[0],f,level+1);fprintf(f," r,");ExOp::show((*this)[1],f,level+1);fprintf(f," i)");
			break;
	}
}
LFHTEMP string Complex<C>::type_tostring() const{return string("Complex<") + ExOp::type_tostring(data[0]) + string(">");}

LFHTEMP	Quaternion<C>::Quaternion() {}

LFHTEMP	const Quaternion<C>& Quaternion<C>::to_normal(const C& _i,const C& _j,const C& _k){
    return(*this);
}

LFHTEMP	const Quaternion<C>& Quaternion<C>::to_normal_and_scale(const C& _i,const C& _j,const C& _k,const C& _s){
    data[0] = sqrt(_s);
    data[1] = 0.0f;
    data[2] = 0.0f;
    data[3] = 0.0f;
    return(*this);
}

LFHTEMP	double Quaternion<C>::pnorm() const{return ExOp::pnorm(data[0]) + ExOp::pnorm(data[1]) +ExOp::pnorm(data[2]) + ExOp::pnorm(data[3]);}

LFHTEMP	const Quaternion<C>& Quaternion<C>::rotateX(double ang){
    double s = sin(ang); double c = cos(ang);
    C tmp = data[0] * c + (data[1] * s);
    data[1] =  data[1] * c + (data[0] * (-s)); data[0] = tmp;
    tmp = data[2] * c + (data[3] * (-s));
    data[3] = data[3] * c + (data[2] * s); data[2] = tmp;
    return(*this);
}
LFHTEMP	const Quaternion<C>& Quaternion<C>::rotateY(double ang){
    double s = sin(ang); double c = cos(ang);
    return(*this);
}
LFHTEMP	const Quaternion<C>& Quaternion<C>::rotateZ(double ang){
    double s = sin(ang); double c = cos(ang);
    return(*this);
}

LFHTEMP	void Quaternion<C>::mk_proj_matrix(TMatrix<C,3,3> &fout)const{
    fout.data[0] = data[0] * data[0];
    C tmp = data[1] * data[1];
    fout.data[4] = fout.data[0] - tmp;
    fout.data[0] += tmp;
    tmp = data[2] * data[2];
    fout.data[8] = fout.data[4] - tmp;
    fout.data[4] += tmp;
    fout.data[0] -= tmp;
    tmp = data[3] * data[3];
    fout.data[8] += tmp;
    fout.data[4] -= tmp;
    fout.data[0] -= tmp;

    fout.data[1] = data[0] * data[3] *2.0f;
    fout.data[5] = data[0] * data[1]*2.0f;
    fout.data[6] = data[0] * data[2]*2.0f;
    tmp = data[1] * data[2]*2.0f;
    fout.data[3] = tmp - fout.data[1];
    fout.data[1] += tmp;


    tmp = data[3] * data[2]*2.0f;
    fout.data[7] = tmp - fout.data[5];
    fout.data[5] += tmp;

    tmp = data[1] * data[3]*2.0f;
    fout.data[2] = tmp - fout.data[6];
    fout.data[6] += tmp;

}
LFHTEMP	void Quaternion<C>::mk_proj_matrix(TMatrix<C,4,4> &fout)const{
    fout.data[0] = data[0] * data[0];
    C tmp = data[1] * data[1];
    fout.data[5] = fout.data[0] - tmp;
    fout.data[0] += tmp;
    tmp = data[2] * data[2];
    fout.data[10] = fout.data[5] - tmp;
    fout.data[5] += tmp;
    fout.data[0] -= tmp;
    tmp = data[3] * data[3];
    fout.data[10] += tmp;
    fout.data[5] -= tmp;
    fout.data[0] -= tmp;

    fout.data[1] = data[0] * data[3]*2.0f;
    fout.data[6] = data[0] * data[1]*2.0f;
    fout.data[8] = data[0] * data[2]*2.0f;
    tmp = data[1] * data[2]*2.0f;
    fout.data[4] = tmp - fout.data[1];
    fout.data[1] += tmp;

    tmp = data[3] * data[2]*2.0f;
    fout.data[9] = tmp - fout.data[6];
    fout.data[6] += tmp;

    tmp = data[1] * data[3]*2.0f;
    fout.data[2] = tmp - fout.data[8];
    fout.data[8] += tmp;

    ExOp::toZero(fout.data[3]);
    ExOp::toZero(fout.data[7]);
    ExOp::toZero(fout.data[11]);
    ExOp::toZero(fout.data[12]);
    ExOp::toZero(fout.data[13]);
    ExOp::toZero(fout.data[14]);
    ExOp::toOne(fout.data[15]);
}
LFHTEMP	Tuple<C,3> Quaternion<C>::mkXvector() const{Tuple<C,3> fout;

    fout[0] = data[0] * data[0] + data[1] * data[1] - (data[2] * data[2] + data[3] * data[3]);
    fout[1] = ( data[1] * data[2] - data[0] * data[3]) *2.0f;
    fout[2] = ( data[1] * data[3] + data[0] * data[2]) *2.0f;

    return(fout);
}
LFHTEMP	Tuple<C,3> Quaternion<C>::mkYvector() const{Tuple<C,3> fout;
    fout[1] = data[0] * data[0] + data[2] * data[2] - (data[1] * data[1] + data[3] * data[3]);
    fout[2] = ( data[2] * data[3] - data[0] * data[1]) *2.0f;
    fout[0] = ( data[2] * data[1] + data[0] * data[3]) *2.0f;
    return(fout);
}
LFHTEMP	Tuple<C,3> Quaternion<C>::mkZvector() const{Tuple<C,3> fout;
    fout[2] = data[0] * data[0] + data[3] * data[3] - (data[2] * data[2] + data[1] * data[1]);
    fout[0] = ( data[3] * data[1] - data[0] * data[2]) *2.0f;
    fout[1] = ( data[3] * data[2] + data[0] * data[1]) *2.0f;
    return(fout);
}
LFHTEMP	template<class O> void Quaternion<C>::wrXvector(O* fout) const{
    fout[0] = data[0] * data[0] + data[1] * data[1] - (data[2] * data[2] + data[3] * data[3]);
    fout[1] = ( data[1] * data[2] - data[0] * data[3]) *2.0f;
    fout[2] = ( data[1] * data[3] + data[0] * data[2]) *2.0f;
}
LFHTEMP	template<class O> void Quaternion<C>::wrYvector(O* fout) const{
    fout[1] = data[0] * data[0] + data[2] * data[2] - (data[1] * data[1] + data[3] * data[3]);
    fout[2] = ( data[2] * data[3] - data[0] * data[1]) *2.0f;
    fout[0] = ( data[2] * data[1] + data[0] * data[3]) *2.0f;
}
LFHTEMP	template<class O> void Quaternion<C>::wrZvector(O* fout) const{
    fout[2] = data[0] * data[0] + data[3] * data[3] - (data[2] * data[2] + data[1] * data[1]);
    fout[0] = ( data[3] * data[1] - data[0] * data[2]) *2.0f;
    fout[1] = ( data[3] * data[2] + data[0] * data[1]) *2.0f;
}

LFHTEMP    const Quaternion<C>& Quaternion<C>::toUnitQuaternion(const Tuple<double, 3> value){
    double norm = sqrt(value[0] * value[0] + value[1] *value[1] + value[2]*value[2]);
    ExOp::toOne(data[1]);
    ExOp::toOne(data[2]);
    ExOp::toOne(data[3]);

    double ang = norm * M_PI * 0.5f;
    ExOp::toOne(data[0]); data[0] *= cos(ang);
    norm = sin(ang) / norm;
    data[1] *= value[0] * norm;
    data[2] *= value[1] * norm;
    data[3] *= value[2] * norm;
    return(*this);
}

// W2 + X2 -Y2 -Z2   2WZ + 2XY     2WY + 2XZ
// 2WZ - 2XY     W2 + Y2 - X2 -Z2   2WX + 2YZ
// 2WY + 2XZ         2WX - 2YZ   W2 + Z2 - Y2 -X2

// COS2 + 2SIN2X2 - 1   2SIN2WZ + 2SIN2XY     2WY + 2XZ
// 2COSSINZ - 2SIN2XY     COS2 + 2SIN2Y2 - 1   2WX + 2YZ
// 2COSSINY + 2SIN2XZ         2WX - 2YZ   COS2 + 2SIN2Z2 - 1

// COS2 + 2SIN2X2 - 1   2SIN2WZ + 2SIN2XY     2WY + 2XZ
// 2SINZ + (-1 + COS)XY     COS2 + 2SIN2Y2 - 1   2WX + 2YZ
// 2SINY + 2SIN2XZ         2WX - 2YZ   COS2 + 2SIN2Z2 - 1

LFHTEMP	const C&  Quaternion<C>::operator[](unsigned int w) const{return(data[w]);}
LFHTEMP	C& Quaternion<C>::operator[](unsigned int w){return(data[w]);}

LFHTEMP	void Quaternion<C>::zero() {ExOp::toZero((*this)[0]);ExOp::toZero((*this)[1]);ExOp::toZero((*this)[2]);ExOp::toZero((*this)[3]);}
LFHTEMP	void Quaternion<C>::one() {ExOp::toOne((*this)[0]);ExOp::toZero((*this)[1]);ExOp::toZero((*this)[2]);ExOp::toZero((*this)[3]);}
LFHTEMP	void Quaternion<C>::random() {ExOp::toRand((*this)[0]);ExOp::toRand((*this)[1]);ExOp::toRand((*this)[2]);ExOp::toRand((*this)[3]);}

LFHTEMP Quaternion<C> Quaternion<C>::operator-()const{Quaternion<C> fout; fout[0] = -(*this)[0];fout[1] = -(*this)[1];fout[2] = -(*this)[2];fout[3] = -(*this)[3]; return fout; }

LFHTEMP const Quaternion<C>& Quaternion<C>::inverse() const {
	Quaternion<C> _out;
	C tmpA = ((*this)[0] * (*this)[0]) + ((*this)[1] * (*this)[1])+ ((*this)[2] * (*this)[2])+ ((*this)[3] * (*this)[3]);
	_out[0] =(*this)[0] / tmpA;
	tmpA = -tmpA;
	_out[1] =(*this)[1] / tmpA;
	_out[2] =(*this)[2] / tmpA;
	_out[3] =(*this)[3] / tmpA;
	return(Complex<C>(_out));
}

LFHTEMP	template<class A>
Quaternion<C>& Quaternion<C>::operator+=(Quaternion<A> const & o){
		ExOp::toAdd(data[0],o[0]);
		ExOp::toAdd(data[1],o[1]);
		ExOp::toAdd(data[2],o[2]);
		ExOp::toAdd(data[3],o[3]);
		return(*this);
	}

LFHTEMP	template<class A>
Quaternion<C>& Quaternion<C>::operator-=(Quaternion<A> const & o){
		ExOp::toSubt(data[0],o[0]);
		ExOp::toSubt(data[1],o[1]);
		ExOp::toSubt(data[2],o[2]);
		ExOp::toSubt(data[3],o[3]);
		return(*this);
	}


LFHTEMP	template<Tuple_flag CF> Tuple<C, 3,CF> Quaternion<C>::operator*(const Tuple<C, 3,CF> &val) const{ Tuple<C, 3> fout;
    C nox,noy,noz;

    nox = -data[1] * data[0];
    noy = -data[2] * data[0];
    noz = -data[3] * data[0];

	C nax = data[2] * data[3];
	C nay = data[1] * data[3];
	C naz = data[1] * data[2];
	C rx2 = data[1] * data[1];
	C ry2 = data[2] * data[2];
	C rz2 = data[3] * data[3];
	C rw2 = data[0] * data[0];

    fout[0] = val[0] * (rw2 + rx2 - ry2 -rz2) + (val[1] * (naz-noz) + val[2] * (noy+nay)) * 2.0f;
    fout[1] = val[1] * (rw2 + ry2 - rx2 -rz2) + (val[0] * (naz+noz) + val[2] * (nax-nox)) * 2.0f;
    fout[2] = val[2] * (rw2 + rz2 - rx2 -ry2) + (val[0] * (nay-noy) + val[1] * (nax+nox)) * 2.0f;

    return fout;
}

LFHTEMP	template<class A,Tuple_flag CF> Tuple<A, 3,CF> Quaternion<C>::operator*(const Tuple<A, 3,CF> &val) const{ Tuple<A, 3> fout;
    C nox,noy,noz;

    nox = -data[1] * data[0];
    noy = -data[2] * data[0];
    noz = -data[3] * data[0];

	C nax = data[2] * data[3];
	C nay = data[1] * data[3];
	C naz = data[1] * data[2];
	C rx2 = data[1] * data[1];
	C ry2 = data[2] * data[2];
	C rz2 = data[3] * data[3];
	C rw2 = data[0] * data[0];

    fout[0] = val[0] * (rw2 + rx2 - ry2 -rz2) + (val[1] * (naz -noz) + val[2] * (noy +nay)) * 2.0f;
    fout[1] = val[1] * (rw2 + ry2 - rx2 -rz2) + (val[0] * (noz +naz) + val[2] * (nax -nox)) * 2.0f;
    fout[2] = val[2] * (rw2 + rz2 - rx2 -ry2) + (val[0] * (nay -noy) + val[1] * (nax +nox)) * 2.0f;
    return fout;
}

LFHTEMP	Quaternion<C> Quaternion<C>::operator*(const Quaternion<C> &other)const{Quaternion<C> fout;
    fout.data[0] = data[0]*other.data[0] -data[1]*other.data[1] -data[2]*other.data[2] -data[3]*other.data[3];
	fout.data[1] = data[0]*other.data[1] +data[1]*other.data[0] -data[2]*other.data[3] +data[3]*other.data[2];
	fout.data[2] = data[0]*other.data[2] +data[1]*other.data[3] +data[2]*other.data[0] -data[3]*other.data[1];
	fout.data[3] = data[0]*other.data[3] -data[1]*other.data[2] +data[2]*other.data[1] +data[3]*other.data[0];
    return fout;
}

LFHTEMP	template<class OC> Quaternion<C> Quaternion<C>::operator*(const Quaternion<OC> &other)const{Quaternion<C> fout;
    fout.data[0] = data[0]*other.data[0] -data[1]*other.data[1] -data[2]*other.data[2] -data[3]*other.data[3];
	fout.data[1] = data[0]*other.data[1] +data[1]*other.data[0] -data[2]*other.data[3] +data[3]*other.data[2];
	fout.data[2] = data[0]*other.data[2] +data[1]*other.data[3] +data[2]*other.data[0] -data[3]*other.data[1];
	fout.data[3] = data[0]*other.data[3] -data[1]*other.data[2] +data[2]*other.data[1] +data[3]*other.data[0];
    return fout;
}

LFHTEMP	Quaternion<C> Quaternion<C>::operator/(const Quaternion<C> &other)const{Quaternion<C> fout;
    double c = pow(ExOp::pnorm(other), -0.5f);
    fout.data[0] = (data[0]*other.data[0] +data[1]*other.data[1] +data[2]*other.data[2] +data[3]*other.data[3])* c;
	fout.data[1] = (-data[0]*other.data[1] +data[1]*other.data[0] +data[2]*other.data[3] -data[3]*other.data[2])* c;
	fout.data[2] = (-data[0]*other.data[2] -data[1]*other.data[3] +data[2]*other.data[0] +data[3]*other.data[1])* c;
	fout.data[3] = (-data[0]*other.data[3] +data[1]*other.data[2] -data[2]*other.data[1] +data[3]*other.data[0])* c;
    return fout;
}

LFHTEMP	template<class OC> Quaternion<C> Quaternion<C>::operator/(const Quaternion<OC> &other)const{Quaternion<C> fout;
    double c = pow(ExOp::pnorm(other), -0.5f);
    fout.data[0] = (data[0]*other.data[0] +data[1]*other.data[1] +data[2]*other.data[2] +data[3]*other.data[3]) * c;
	fout.data[1] = (-data[0]*other.data[1] +data[1]*other.data[0] +data[2]*other.data[3] -data[3]*other.data[2]) * c;
	fout.data[2] = (-data[0]*other.data[2] -data[1]*other.data[3] +data[2]*other.data[0] +data[3]*other.data[1]) * c;
	fout.data[3] = (-data[0]*other.data[3] +data[1]*other.data[2] -data[2]*other.data[1] +data[3]*other.data[0]) * c;
    return fout;
}

LFHTEMP	Quaternion<C>& Quaternion<C>::operator*=(const Quaternion<C> &other){
    C tmp[3];
    ExOp::toMemmove(tmp[0],data[0]);
    ExOp::toMemmove(tmp[1],data[1]);
    ExOp::toMemmove(tmp[2],data[2]);
    data[0] = tmp[0]*other.data[0] -tmp[1]*other.data[1] -tmp[2]*other.data[2] -data[3]*other.data[3];
	data[1] = tmp[0]*other.data[1] +tmp[1]*other.data[0] -tmp[2]*other.data[3] +data[3]*other.data[2];
	data[2] = tmp[0]*other.data[2] +tmp[1]*other.data[3] +tmp[2]*other.data[0] -data[3]*other.data[1];
	data[3] = tmp[0]*other.data[3] -tmp[1]*other.data[2] +tmp[2]*other.data[1] +data[3]*other.data[0];
    return *this;
}

LFHTEMP	Quaternion<C>& Quaternion<C>::operator/=(const Quaternion<C> &other){
    double c = pow(ExOp::pnorm(other), -0.5f);
    C tmp[3];
    ExOp::toMemmove(tmp[0],data[0]);
    ExOp::toMemmove(tmp[1],data[1]);
    ExOp::toMemmove(tmp[2],data[2]);
    data[0] = (tmp[0]*other.data[0] +tmp[1]*other.data[1] +tmp[2]*other.data[2] +data[3]*other.data[3]) *c;
	data[1] = (-tmp[0]*other.data[1] +tmp[1]*other.data[0] +tmp[2]*other.data[3] -data[3]*other.data[2]) *c;
	data[2] = (-tmp[0]*other.data[2] -tmp[1]*other.data[3] +tmp[2]*other.data[0] +data[3]*other.data[1]) *c;
	data[3] = (-tmp[0]*other.data[3] +tmp[1]*other.data[2] -tmp[2]*other.data[1] +data[3]*other.data[0]) *c;
    return *this;
}

LFHTEMP	Quaternion<C>& Quaternion<C>::toBackMult(const Quaternion<C> &other){
    C tmp[3];
    ExOp::toMemmove(tmp[0],data[0]);
    ExOp::toMemmove(tmp[1],data[1]);
    ExOp::toMemmove(tmp[2],data[2]);
    data[0] = other.data[0]*tmp[0] -other.data[1]*tmp[1] -other.data[2]*tmp[2] -other.data[3]*data[3];
	data[1] = other.data[0]*tmp[1] +other.data[1]*tmp[0] -other.data[2]*data[3] +other.data[3]*tmp[2];
	data[2] = other.data[0]*tmp[2] +other.data[1]*data[3] +other.data[2]*tmp[0] -other.data[3]*tmp[1];
	data[3] = other.data[0]*data[3] -other.data[1]*tmp[2] +other.data[2]*tmp[1] +other.data[3]*tmp[0];
    return *this;
}

LFHTEMP	Quaternion<C>& Quaternion<C>::toBackDivi(const Quaternion<C> &other){
    double c = pow(ExOp::pnorm(other), -0.5f);
    C tmp[3];
    ExOp::toMemmove(tmp[0],data[0]);
    ExOp::toMemmove(tmp[1],data[1]);
    ExOp::toMemmove(tmp[2],data[2]);
    data[0] = (other.data[0]*tmp[0] +other.data[1]*tmp[1] +other.data[2]*tmp[2] +other.data[3]*data[3]) *c;
	data[1] = (other.data[0]*tmp[1] -other.data[1]*tmp[0] +other.data[2]*data[3] -other.data[3]*tmp[2])*c;
	data[2] = (other.data[0]*tmp[2] -other.data[1]*data[3] -other.data[2]*tmp[0] +other.data[3]*tmp[1])*c;
	data[3] = (other.data[0]*data[3] +other.data[1]*tmp[2] -other.data[2]*tmp[1] -other.data[3]*tmp[0])*c;
    return *this;
}


LFHTEMP	template<class OC, class OB> void Quaternion<C>::wrMatrix(TMatrix<OC, 3,3>& fout, const OB& scale, bool transpose) const{
    C nox,noy,noz;
    if (transpose){
        nox = data[1] * data[0];
        noy = data[2] * data[0];
        noz = data[3] * data[0];
    }else{
        nox = -data[1] * data[0];
        noy = -data[2] * data[0];
        noz = -data[3] * data[0];
    }
	C nax = data[2] * data[3];
	C nay = data[1] * data[3];
	C naz = data[1] * data[2];
	C rx2 = data[1] * data[1];
	C ry2 = data[2] * data[2];
	C rz2 = data[3] * data[3];
	C rw2 = data[0] * data[0];

    OB nscale = scale / (rx2 +ry2 +rz2 +rw2);
    OB lf =  nscale * 2.0f ;

	fout.data[0] = (rw2 + rx2 - ry2 -rz2) * nscale;
	fout.data[3] = (naz+noz)  * lf;
	fout.data[6] = (nay-noy)  * lf;

	fout.data[1] = (naz-noz)  * lf;
	fout.data[4] = (rw2 + ry2 - rx2 -rz2)  * nscale;
	fout.data[7] = (nox+nax)  * lf;

	fout.data[2] = (noy+nay) * lf ;
	fout.data[5] = (nax-nox) * lf ;
	fout.data[8] = (rw2 + rz2 - rx2 -ry2) * nscale;
}



LFHTEMP	template<class OC, class OB> void Quaternion<C>::wrMatrix(TMatrix<OC, 4,4>& fout, const OB& scale, bool transpose, bool is_final) const{
    C nox,noy,noz;
    if (transpose){
        nox = data[1] * data[0];
        noy = data[2] * data[0];
        noz = data[3] * data[0];
    }else{
        nox = -data[1] * data[0];
        noy = -data[2] * data[0];
        noz = -data[3] * data[0];
    }
	C nax = data[2] * data[3];
	C nay = data[1] * data[3];
	C naz = data[1] * data[2];
	C rx2 = data[1] * data[1];
	C ry2 = data[2] * data[2];
	C rz2 = data[3] * data[3];
	C rw2 = data[0] * data[0];

    OB nscale = scale / (rx2 +ry2 +rz2 +rw2);
    OB lf =  nscale * 2.0f ;

	fout.data[0] = (rw2 + rx2 - ry2 -rz2) * nscale;
	fout.data[4] = (naz+noz)  * lf;
	fout.data[8] = (nay-noy)  * lf;

	fout.data[1] = (naz-noz)  * lf;
	fout.data[5] = (rw2 + ry2 - rx2 -rz2)  * nscale;
	fout.data[9] = (nox+nax)  * lf;

	fout.data[2] = (noy+nay) * lf ;
	fout.data[6] = (nax-nox) * lf ;
	fout.data[10] = (rw2 + rz2 - rx2 -ry2) * nscale;

    if (!is_final) return;
    ExOp::toZero(fout.data[3]);
    ExOp::toZero(fout.data[7]);
    ExOp::toZero(fout.data[11]);
    ExOp::toZero(fout.data[12]);
    ExOp::toZero(fout.data[13]);
    ExOp::toZero(fout.data[14]);
    ExOp::toOne(fout.data[15]);
}

LFHTEMP Quaternion<C> Quaternion<C>::mkInterpolated(const Quaternion<C>& origin, float factor) const{ Quaternion<C> fout;
    float tmp = 1.0f - factor;
    fout.data[0] = (factor * data[0]) + tmp * origin.data[0];
    fout.data[1] = (factor * data[1]) + tmp * origin.data[1];
    fout.data[2] = (factor * data[2]) + tmp * origin.data[2];
    fout.data[3] = (factor * data[3]) + tmp * origin.data[3];
    fout.toUnitary();
    return fout;
}

LFHTEMP C Quaternion<C>::dotProduct(const Quaternion<C>& other)const{
	return (data[0] * other[0]) + (data[1] * other[1]) + (data[2] * other[2]) + (data[3] * other[3]);
}

LFHTEMP Quaternion<C> Quaternion<C>::mkCloseInterpolated(const Quaternion<C>& origin, float factor) const{ Quaternion<C> fout;
    float tmp = (data[0] * origin[0]) + (data[1] * origin[1]) + (data[2] * origin[2]) + (data[3] * origin[3]) > 0 ? 1.0f - factor : factor - 1.0f;
    fout.data[0] = (factor * data[0]) + tmp * origin.data[0];
    fout.data[1] = (factor * data[1]) + tmp * origin.data[1];
    fout.data[2] = (factor * data[2]) + tmp * origin.data[2];
    fout.data[3] = (factor * data[3]) + tmp * origin.data[3];
    fout.toUnitary();
    return fout;
}
LFHTEMP Quaternion<C>& Quaternion<C>::toUnitary(){
    double norm = ExOp::pnorm(data[0]) + ExOp::pnorm(data[1]) + ExOp::pnorm(data[2]) + ExOp::pnorm(data[3]);
    norm = pow(norm, -0.5f);
    data[0] *= norm;    data[1] *= norm;    data[2] *= norm;    data[3] *= norm;
    return *this;
}
LFHTEMP void Quaternion<C>::show(FILE* f, int level) const{
    	switch(level){
		case 1:
		case 0:
            ExOp::show((*this)[0],f,level+1);fprintf(f," r\t");ExOp::show((*this)[1],f,level+1);fprintf(f," i\t");ExOp::show((*this)[2],f,level+1);fprintf(f," j\t");ExOp::show((*this)[3],f,level+1);fprintf(f," k\n");
		break;
		case 2:
            fprintf(f,"[");ExOp::show((*this)[0],f,level+1);fprintf(f," r;");ExOp::show((*this)[1],f,level+1);fprintf(f," i;");ExOp::show((*this)[2],f,level+1);fprintf(f," j;");ExOp::show((*this)[3],f,level+1);fprintf(f," k]");
			break;
		default:
           fprintf(f,"(");ExOp::show((*this)[0],f,level+1);fprintf(f," r,");ExOp::show((*this)[1],f,level+1);fprintf(f," i,");ExOp::show((*this)[2],f,level+1);fprintf(f," j,");ExOp::show((*this)[3],f,level+1);fprintf(f," k)");
			break;
	}
}
LFHTEMP string Quaternion<C>::type_tostring() const{return string("Quaternion<") + string(">");}

LFHTEMP void Trianglix<C,0u>::given_rotation_routine(double& c, double& s, double& c2,uint32_t offset, uint32_t maxsize){
    uint32_t i;
    C* rows[2];
    C tmp;
    rows[0] = data + ((offset + (offset +1)) >> 1);
    rows[1] = rows[0] + offset + 1;
    for(i=0;i<offset;i++){
        tmp = rows[1][i] * c + rows[0][i] * s;
        rows[0][i] =  rows[0][i] * c - rows[1][i] * s;
        rows[1][i] = tmp;
    }
    i++;
    rows[1][offset] = rows[1][offset] * (c2* 2.0 -1.0) + (rows[1][i] - rows[0][offset]) * c * s;
    tmp = rows[1][i] * c2 + rows[0][offset] * (c2 - 1.0);
    rows[0][offset] = rows[0][offset] * c2 + rows[1][i] * (c2 - 1.0);
    rows[1][i] = tmp;
    for(i++;i<maxsize;i++){
        rows[0] = ((i + (i +1)) >> 1) + offset;
        tmp = rows[0][1] * c + rows[0][0] * s;
        rows[0][0] =  rows[0][0] * c - rows[0][1] * s;
        rows[0][1] = tmp;
    }
}
LFHTEMP template<unsigned int SIZE, Tuple_flag TF> Trianglix<C, 0u>::Trianglix(const Tuple<C, SIZE, TF>& other): t_size(other.getSize()) {
    data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(other[i]) * other[j];
    }
LFHTEMP template<unsigned int SIZE, Tuple_flag TF> Trianglix<C, 0u>::Trianglix(const Tuple<C, SIZE, TF>& u, const Tuple<C, SIZE>& v): t_size(SIZE) {
    data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(u[i]) * v[j] + ExOp::mkTrju(v[i]) * u[j];
}
LFHTEMP Trianglix<C, 0u>::Trianglix(const Tuple<C, 0u>& other): t_size(other.getSize()) {
    data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(other[i]) * other[j];
    }
LFHTEMP Trianglix<C, 0u>::Trianglix(const Tuple<C, 0u>& u, const Tuple<C, 0u>& v): t_size(u.getSize()) {
    data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(u[i]) * v[j] + ExOp::mkTrju(v[i]) * u[j];
}
LFHTEMP Trianglix<C, 0u>::Trianglix(const Trianglix<C, 0u>& other): t_size(other.getSize()) {
    if (t_size){
        data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    }else data = NULL;
    unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i];}
LFHTEMP template<unsigned int SIZE> Trianglix<C, 0u>::Trianglix(const Trianglix<C, SIZE>& other): t_size(other.getSize()) {
    if (t_size){
        data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    }else data = NULL;
    unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i];
    }
LFHTEMP template<unsigned int SIZE> Trianglix<C, 0u>& Trianglix<C, 0u>::operator=(const Tuple<C, SIZE>& other){
    setSize(SIZE);
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(other[i]) * other[j];
    return(*this);
}
LFHTEMP Trianglix<C, 0u>& Trianglix<C, 0u>::toMemmove(Trianglix<C, 0u>& other){
	if (t_size) delete[](data);
	data = other.data;
	t_size = other.t_size;
	other.t_size = 0;
    return(*this);
}
LFHTEMP Trianglix<C, 0u>::~Trianglix(){if (t_size) {delete[](data);}}
LFHTEMP Trianglix<C, 0u>& Trianglix<C, 0u>::setSize(unsigned int s){if (s == t_size) return *this; if (t_size) {delete[](data);} t_size =s;
    if (t_size){
        data = new C[totsize()]; LFH_NICE_ALLOCERROR(data,"")
    }else data = NULL;
return *this;}
LFHTEMP double Trianglix<C, 0u>::pnorm() const {
	double fout=0; unsigned int i,j;
	C* k = data;
	for(j=0;j<this->getSize();j++){
		for(i=0;i<j;i++,k++) fout += ExOp::pnorm(*k);
		fout += ExOp::pnorm(*k) * 0.5f;k++;
	}
	fout *= 2.0f;
	return fout;
}
LFHTEMP Matrix<C> Trianglix<C, 0u>::makeMatrix()const{
    Matrix<C> damat;damat.setSizes(t_size,t_size);
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) {
        for(i=0;i<j;i++,k++) {damat(i,j) = data[k];damat(j,i) = data[k];}
        damat(j,j) = data[k]; k++;
    }
    return(damat);
    }

LFHTEMP Trianglix<C, 0u>::operator Matrix<C>()const{Matrix<C> damat;
    damat.setSizes(t_size,t_size);
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) {
        for(i=0;i<j;i++,k++) {damat(i,j) = data[k];damat(j,i) = data[k];}
        damat(j,j) = data[k]; k++;
    }
    return(damat);
}
LFHTEMP Trianglix<C, 0u>::operator TMatrix<C,0u,0u>()const{TMatrix<C,0u,0u> damat;
    damat.setSizes(t_size,t_size);
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) {
        for(i=0;i<j;i++,k++) {damat(i,j) = data[k];damat(j,i) = data[k];}
        damat(j,j) = data[k++];
    }
    return(damat);
}



LFHTEMP Trianglix<C, 0u> Trianglix<C,0u>::mkEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc) const{Trianglix<C,0u> fout = *this; fout.toEigenTransform(fnc);return fout;}
LFHTEMP template<class F> Trianglix<C,0u> Trianglix<C,0u>::mkEigenTransform(F fnc) const{Trianglix<C,0u> fout = *this; fout.toEigenTransform(fnc);return fout;}



LFHTEMP void Trianglix<C, 0u>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length, C* buf, bool hint ){
    // house holder multiplication, done twice!
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;
    C s;
        if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2))) return;
        buf[0] = data[0] * vec[0];
        s = data[0] * ExOp::mkTrju(vec[0]) * vec[0];

	for(j=1,k=1;j<length;j++,k++) {
        buf[j] = data[k] * vec[0]; buf[0] += ExOp::mkTrju(data[k++]) * vec[j];
        for(i=1;i<j;i++) {buf[j] += data[k] * vec[i]; buf[i] += ExOp::mkTrju(data[k++]) * vec[j];}
        s += data[k] * ExOp::mkTrju(vec[j]) * vec[j] + (buf[j] *  ExOp::mkTrju(vec[j]) + ExOp::mkTrju(buf[j]) * vec[j]);
        buf[j] += data[k] * vec[j];
	}

    if (!hint){
        for(;j<t_size;j++) {
            buf[j] = data[k++] * vec[0];
            for(i=1;i<length;i++) buf[j] += data[k++] * vec[i];
            k+= j - length+1;
        }
        for(j=0;j<t_size;j++) buf[j] /= denum2;
    }else for(j=0;j<length;j++) buf[j] /= denum2;
    s /= denum2 * denum2;

    for(j=0,k=0;j<length;j++) {
        for(i=0;i<j;i++) data[k++] += s * vec[j] *ExOp::mkTrju(vec[i]) - buf[j]*ExOp::mkTrju(vec[i]) - vec[j] * ExOp::mkTrju(buf[i]) ;
        data[k++] += s * vec[j] *ExOp::mkTrju(vec[j]) - buf[j]*ExOp::mkTrju(vec[j]) - vec[j] * ExOp::mkTrju(buf[j]);
    }
    if (hint) return;
    for(;j<t_size;j++) {
        for(i=0;i<length;i++) data[k++] -= buf[j]*ExOp::mkTrju(vec[i]);
        k+= j - length+1;
    }
}
LFHTEMP void Trianglix<C, 0u>::offdiagelimination_down(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)"
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,k;

    k = (col *(col+1))/2;
    for(i=0;i<=col;i++) {buf[i] = fact * data[k+i];}
    for(;i<t_size;i++) {k += i; buf[i] = fact * data[k+col];}
    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) data[k+i] += buf[i];
    k++;

    data[k+col] += buf[i] + ExOp::mkTrju(buf[i]) + fact * data[k-2] * ExOp::mkTrju(fact);
    for(i++;i<t_size;i++) {k += i; data[k+col] += buf[i];}
}
LFHTEMP void Trianglix<C, 0u>::offdiagelimination_up(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)"
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,k;

    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) buf[i] = fact * data[k+i];
    buf[i] = fact * data[k+i];
    k++;
    for(i++;i<t_size;i++) {k += i; buf[i] = fact * data[k+col];}

    k = ((col) *(col+1))/2;
    for(i=0;i<col;i++) data[k+i] += buf[i];
    data[k+col] += ExOp::mkTrju(buf[col+1]) * fact + ExOp::mkTrju(buf[i]);
    for(i++;i<t_size;i++) {k += i; data[k+col] += buf[i];}
}
LFHTEMP void Trianglix<C, 0u>::offdiagelimination_up_backroutine(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)"
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,k;

    k = ((col+2) *(col+1))/2;
    i =col+1;
    buf[i] = fact * data[k+i];
    k++;
    for(i++;i<t_size;i++) {k += i; buf[i] = fact * data[k+col];}

    k = ((col) *(col+1))/2;

    data[k+col] += ExOp::mkTrju(buf[col+1]) * fact;
    for(i=col+1;i<t_size;i++) {k += i; data[k+col] = buf[i];}
}
LFHTEMP void Trianglix<C, 0u>::QR_back(const C &factC,const C &factS,unsigned int col){
    C tmp;
    unsigned int k,i;
    k = (col *(col+1))/2;
    for(i=0;i<col;i++){
    tmp = data[k+i] * factC + data[k+i+col+1] * factS;
    data[k+i+col+1] = data[k+i+col+1] * factC - data[k+i] * factS;
    data[k+i] = tmp;

 //   printf("bepair = %i %i \n",k+i, k+i+col+1);
    }
    // A

    k += i+col+1;
 //   tmp = data[k-col-1] * factC + data[k] * factS;
 //   C tmp2 = -data[k-col-1] * factS + data[k] * factC;
 //   C tmp3 = data[k] * factC + data[k+1] * factS;
 //   C tmp4 = -data[k] * factS + data[k+1] * factC;
 //   printf("dapox = %i %i %i \n",k-col-1,k, k+1);
//    data[k] = - tmp * factS + tmp3 * factC; // or tmp2 * factC + tmp4 * factS
    tmp = data[k-col-1] * factS* factS + data[k+1] * factC* factC - 2.0f * data[k] * factC * factS;
    data[k] = data[k] * (factC * factC - factS * factS) - (data[k-col-1] - data[k+1]) * factC * factS;

    data[k-col-1] += data[k+1] -tmp;
    data[k+1] = tmp;

    for(i=col+2;i<t_size;i++){k+=i;
    tmp = data[k] * factC + data[k+1] * factS;
    data[k+1] = data[k+1] * factC - data[k] * factS;
    data[k] = tmp;
 //   printf("afpair = %i %i \n",k, k+1);
    }

}
LFHTEMP Trianglix<C, 0u>& Trianglix<C, 0u>::operator*=( const Trianglix<C, 0u>& other){
    unsigned int i,j,k,l,m,n;
    uint32_t mins = (t_size < other.t_size) ? t_size : other.t_size;
	if (mins == 0) return *this;
	uint32_t patsize = (mins & 1) ? mins * ((mins>>1)+1) : (mins>>1) * (mins+1);
	C* buffer = new C[patsize]; LFH_NICE_ALLOCERROR(buffer,"")

	for(j=0,k=0;j < mins;j++,k++) {
        for(i=0;i<j;i++,k++){
            m = (i * (i+1)) >> 1;
            n = (j * (j+1)) >> 1;
           // printf("K is %i\n", k); printf("mn0 is %i %i\n", m,n);
            buffer[k] = data[n] * ExOp::mkTrju(other.data[m]);
            for(l=1;l<=i;l++){
                m++;n++; //printf("mn1 is %i %i\n", m,n);
                buffer[k] += data[n] * ExOp::mkTrju(other.data[m]);
            }
            for(;l<=j;l++){
                m+=l;n++; //printf("mn2 is %i %i\n", m,n);
                buffer[k] += data[n] * other.data[m];
            }
            for(;l<mins;l++){
                m+=l;n+=l; //printf("mn3 is %i %i\n", m,n);
                buffer[k] += ExOp::mkTrju(data[n]) * other.data[m];
            }
        }

        buffer[k] = other.data[k] * data[k];
        m = (j * (j+1)) >> 1;
        for(l=0; l < j;l++,m++) {buffer[k] += data[m] * ExOp::mkTrju(other.data[m]);}
        l=j+1;
        for(m += l;l<mins; m+= l){
            buffer[k] += ExOp::mkTrju(data[m]) * other.data[m];l++;
        }
	}
	for(i=0;i<patsize;i++) data[i] = buffer[i];
	delete[](buffer);
	return *this;
}

LFHTEMP template<unsigned int SIZE> Tuple<C,SIZE,TUPLE_FLAG_NULL> Trianglix<C, 0u>::operator*( const Tuple<C, SIZE>& other) const{ Tuple<C,SIZE> fout;
    unsigned int i,j,k;
	if (t_size == 0) return fout;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;((j<SIZE)&&(j<t_size));j++) {
		fout[j] = other[0] * data[k];
		fout[0] += other[j] * ExOp::mkTrju(data[k]);
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mkTrju(data[k]);}
        fout[j] += other[j] * data[k++];
	}
	return fout;
}
LFHTEMP Tuple<C,0u,TUPLE_FLAG_NULL> Trianglix<C, 0u>::operator*( const Tuple<C,0u>& other) const{ Tuple<C,0u> fout; fout.setSize(other.getSize());
    if (t_size == 0) return fout;
    unsigned int i,j,k;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;j<t_size;j++) {
		fout[j] = other[0] * data[k];
		fout[0] += other[j] * ExOp::mkTrju(data[k]);
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mkTrju(data[k]);}
        fout[j] += other[j] * data[k++];
	}
	return fout;
}
LFHTEMP template<class O> auto Trianglix<C, 0u>::operator*(const SparseTuple<O> & a)const
 -> Tuple<decltype(this->data[0]* a[0] ),0u,TUPLE_FLAG_NULL>{
    Tuple<decltype(this->data[0]* a[0] ),0u,TUPLE_FLAG_NULL> fout;
    fout.setSize(t_size).toZero();
    uint32_t i,j;
    if (auto ite = a.getIterator()) do{
        const C* currow = data + ((ite() * (ite()+1u)) >> 1);
        for(i=0u;i<a.getSize();i++) {
            if ((j = a.deref_key(i)) < ite()){
            fout[ite()] += currow[j] * a.deref(i);
            fout[j] += currow[j] * (*ite);
        }   }
        fout[ite()] += currow[ite()] * (*ite);
    }while(ite++);
return fout;}
LFHTEMP template<class O> auto Trianglix<C, 0u>::mkInnerMult(const Tuple<O> & a)const
 -> decltype(this->data[0]* a[0] ){
    decltype(this->data[0]* a[0] ) fout;
    if (a.getSize() == 0) ExOp::toZero(fout);
    else{
        fout = data[0] * (a[0] * a[0]);
        for(uint32_t ite = 1 ; ite <- a.getSize(); ite++) {
            const C* currow = data + ((ite * (ite+1u)) >> 1);
            fout += currow[ite] * a[ite]* a[ite];
            for(uint32_t i=0u;i<ite;i++) fout +=  currow[i] * ((a[i] * a[ite]) * 2.0);
        }
    }
return fout;}
LFHTEMP template<class O> auto Trianglix<C, 0u>::mkInnerMult(const SparseTuple<O> & a)const
 -> decltype(this->data[0]* a[0] ){
    decltype(this->data[0]* a[0] ) fout;
    uint32_t i;
    if (auto ite = a.getIterator()){
        const C* currow = data + ((ite() * (ite()+1u)) >> 1);
        fout = currow[ite()] * ((*ite) * (*ite));
        while(true){
            for(i=0u;i<a.getSize();i++) if (a.deref_key(i) < ite()) fout +=  currow[a.deref_key(i)] * ((a.deref(i) * (*ite)) * 2.0);
            if (!(ite++)) break;
            currow = data + ((ite() * (ite()+1u)) >> 1);
            fout += currow[ite()] * (*ite)* (*ite);
        }
    }else ExOp::toZero(fout);
return fout;}
LFHTEMP template<class O> auto Trianglix<C, 0u>::mkInnerMult(const Tuple<O> & a, const Tuple<O> & b)const
 -> decltype(this->data[0]* a[0] ){
    decltype(this->data[0]* a[0] ) fout;
    uint32_t i,j,k;
    uint32_t minsize = (a.getSize() < b.getSize()) ?a.getSize() : b.getSize();
    if (minsize > t_size) minsize = t_size;
    if (minsize == 0) ExOp::toZero(fout);
    else{
        fout = data[0] * (a[0] * b[0]);
        for(j=1,k=1;j<t_size;j++) {
            for(i=0,k++;i<j;i++) fout += data[k++] * (a[i]*b[j] + a[j]*b[i]);
            fout += data[k++] * (a[j] * b[j]);
        }
    }
return fout;}

LFHTEMP template<class O> auto Trianglix<C, 0u>::mkInnerMult(const SparseTuple<O> & a, const SparseTuple<O> & b)const
 -> decltype(this->data[0]* a[0] ){
    decltype(this->data[0]* a[0] ) fout;
    uint32_t i;
    if (auto ite = a.getIterator()) {
        if (auto iteb = b.getIterator()){
            fout = data[(ite() < iteb()) ? ite() + ((iteb() * (iteb() +1))>>1) : iteb() + ((ite() * (ite() +1))>>1) ] * ((*ite)* (*iteb));
            while(iteb++){
                fout += data[(ite() < iteb()) ? ite() + ((iteb() * (iteb() +1))>>1) : iteb() + ((ite() * (ite() +1))>>1) ] * ((*ite)* (*iteb));
            }
            while(ite++){
                if (auto itec = b.getIterator()) do{
                    fout += data[(ite() < itec()) ? ite() + ((itec() * (itec() +1))>>1) : itec() + ((ite() * (ite() +1))>>1) ] * ((*ite)* (*itec));
                }while(itec++);
            }
        }else ExOp::toZero(fout);
    }else ExOp::toZero(fout);
return fout;}
LFHTEMP template<unsigned int SIZE> Tuple<C,SIZE> Trianglix<C, 0u>::divisionof( const Tuple<C, SIZE>& other) const{ Tuple<C,SIZE> fout;
    // TO DO ! partial inversion!
    return(fout);
}
LFHTEMP Tuple<C,0u> Trianglix<C, 0u>::divisionof( const Tuple<C,0u>& dev) const{ Tuple<C,0u> fout; fout.setSize(dev.getSize());

    // dev (P+P2/2)-1 dev
    Trianglix<C, 0u> invsc = (*this);
    Tuple<C, 0u> xdev = dev;
	if (getSize() < xdev.getSize()) exit(1);
	unsigned int siz = xdev.getSize();

    if (siz >2){
		//printf("alloc to die %i\n",siz); fflush(stdout);

    C* buf = new C[siz*2];
    C* hv = new C[siz *(siz+1)/2];
    double* normbuf = new double[siz-2]; LFH_NICE_ALLOCERROR(normbuf,"")

    unsigned int i,j,k;
	//	printf("%i\t%t\n", invsc.getSize(), xdev.getSize()); fflush(stdout);
		for(j=siz-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }

        buf[siz] = ExOp::mkInverse(invsc.data[0]);

        for(j=1;j<siz;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[siz+j-1]*invsc.data[k-1];
            invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[siz+j-1]*invsc.data[k-1];
            buf[siz+j] = ExOp::mkInverse(invsc.data[k]);
            }
        fout[j-1] = buf[siz+j-1] * xdev[j-1];

    for(j=siz-1;j>0;j--) {k = (j*(j+3))/2;fout[j-1] = buf[siz+j-1] * (xdev[j-1] -  fout[j] *invsc.data[k-1]);}

    for(j=2;j<siz;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j);

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (siz == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout[0] = (data[2] * dev[0] - ExOp::mkTrju(data[1]) * dev[1]) * denum;
        fout[1] = (data[0] * dev[1] - data[1] * dev[0]) * denum;
    }  else fout[0] = ExOp::mkInverse(data[0]) * dev[0];


    return(fout);
}
LFHTEMP template<unsigned int SIZE,Tuple_flag TF> C Trianglix<C, 0u>::Xformed_inner_product( const Tuple<C, SIZE,TF>& other) const{
    C s;ExOp::toZero(s);
    C s2;
    unsigned int i,j,k;
        for(j=0,k=0;j<other.getSize();j++) {
        if (j>0) s2 = data[k++] * other[0];
        else ExOp::toZero(s2);
        for(i=1;i<j;i++) s2 += data[k++] * other[i];
        // s += other[j] * ((ExOp::mkTrju(s2) + s2) + data[k++] * other[j]); WORKS FOR REAL ONLY!
        s +=  data[k++] * ExOp::mkTrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
        }
return s;}
/*LFHTEMP template<unsigned int SIZE,Tuple_flag TF> C Trianglix<C, 0u>::Xformed_inner_product_of_squarre( const Tuple<C, SIZE,TF>& other) const{
    C s;ExOp::toZero(s);
    C s2;
    unsigned int i,j,k;
        for(j=0,k=0;j<other.getSize();j++) {
        if (j>0) s2 = data[k++] * other[0];
        else ExOp::toZero(s2);
        for(i=1;i<j;i++) s2 += data[k++] * other[i];
        // s += other[j] * ((ExOp::mkTrju(s2) + s2) + data[k++] * other[j]); WORKS FOR REAL ONLY!
        s +=  data[k++] * ExOp::mkTrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
        }
return s;}*/
LFHTEMP C Trianglix<C, 0u>::Xformed_inner_product( const Tuple<C, 0u>& other) const{
    C s;
    C s2;
    unsigned int i,j,k;
    if (t_size == 0) {ExOp::toZero(s); return s;}
	s = data[0] * other[0] * ExOp::mkTrju(other[0]);

    if (ExOp::mkTrju(other[0]) * other[0] < 0.0f) { fprintf(stderr, "Should not be negative!: "); ExOp::show(ExOp::mkTrju(other[0]) * other[0],stderr); exit(1); }
	for(j=1,k=1;j<t_size;j++) {
        s2 = data[k++] * other[0];
        for(i=1;i<j;i++) s2 += data[k++] * other[i];
        s += 2.0f * other[j] * s2 + data[k++] * other[j]; // WORKS FOR REAL ONLY!
        // s +=  data[k++] * ExOp::mkTrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
		if (ExOp::mkTrju(other[j]) * other[j] < 0.0f) { fprintf(stderr,"Should not be negative!: "); ExOp::show(ExOp::mkTrju(other[j]) * other[j],stderr); exit(1); }
	}
return s;}
LFHTEMP template<class O, unsigned int SIZE, Tuple_flag OF> auto Trianglix<C, 0u>::Xformed_inner_product_of_squarre(const Tuple<O,SIZE,OF> & a)const
 -> decltype(ExOp::mkTrjuProd(this->data[0]* a[0])){
    Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL> fout;
    uint32_t i,k;
    uint32_t damax = a.getSize() < t_size ?  a.getSize(): t_size;
    for(k=0;k<damax;k++){
        const C* currow = data + ((k * (k+1u)) >> 1);
        fout[k] = currow[k] * a[k];
        for(i=0;i<k;i++) {
            fout[k] += currow[i] * a[i];
            fout[i] += currow[i] * a[k];
        }
    }
    fout[0] = ExOp::mkTrjuProd(fout[0]);
    for(k=1;k<fout.getSize();k++) fout[0] += ExOp::mkTrjuProd(fout[k]);
return fout[0];}

	LFHTEMP C Trianglix<C, 0u>::trace_of_product(const Trianglix<C,0u> &other) const{C fout; ExOp::toZero(fout);
		unsigned int i,j,k;
		for(j=0,k=0;j<t_size;j++) {
			for(i=0;i<j;i++,k++) fout += data[k] * ExOp::mkTrju(other.data[k]) + ExOp::mkTrju(data[k]) * other.data[k];
			fout += data[k] * other.data[k];
			k++;
		}
		return fout;
	}
LFHTEMP C Trianglix<C, 0u>::trace_of_division(const Trianglix<C,0u> &divisor) const{C fout; ExOp::toZero(fout);
    // dev (P+P2/2)-1 dev
    Trianglix<C, 0u> invsc = divisor;
    Trianglix<C, 0u> xdev = *this;
	if (getSize() < xdev.getSize()) exit(1);
	unsigned int siz = xdev.getSize();

    if (siz >2){
		//printf("alloc to die %i\n",siz); fflush(stdout);

    C* buf = new C[siz*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[siz *(siz+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[siz-2]; LFH_NICE_ALLOCERROR(normbuf,"")

    unsigned int i,j,k;
		for(j=siz-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,false);
        }

        buf[siz] = ExOp::mkInverse(invsc.data[0]);
        fout = xdev.data[0] * buf[siz];
        for(j=1;j<siz;j++){
            k = (j*(j+3))/2;
            xdev.offdiagelimination_down(-buf[siz+j-1]*invsc.data[k-1],j-1,hv);
            invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[siz+j-1]*invsc.data[k-1];
            buf[siz+j] = ExOp::mkInverse(invsc.data[k]);
            fout += xdev.data[k] * buf[siz+j];
            }

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (siz == 2){

		fout = (data[0] * divisor.data[2] + data[2] * divisor.data[0] + data[1] * ExOp::mkTrju(divisor.data[1]) + ExOp::mkTrju(data[1]) * divisor.data[1]) / (divisor.data[0] * divisor.data[2] - divisor.data[1] * ExOp::mkTrju(divisor.data[1]));


    }  else fout = data[0] * ExOp::mkInverse(divisor.data[0]);
	return fout;
}
LFHTEMP C Trianglix<C, 0u>::trace_of_inverse() const{C fout;
    if (t_size == 0) {ExOp::toZero(fout); return fout;}
    Tuple<C> eigen = this->getEigenValues();
    fout = ExOp::mkInverse(eigen[0]);
    for(uint32_t i=1;i<t_size;i++) fout += ExOp::mkInverse(eigen[i]);
    return fout;
}
LFHTEMP template<class O> Trianglix<C, 0u>& Trianglix<C, 0u>::tosymmetricMult(const Trianglix<O, 0u> & other){ // (other^(1/2)-T) * this * other^(1/2)
    // o = LDL
    // fout = o'to = LDL'tLDL'
    return *this;
    }
LFHTEMP template<class O> Trianglix<C, 0u>& Trianglix<C, 0u>::tosymmetricDivi(const Trianglix<O, 0u> & other){ // (other^(-1/2)-T) * this * other^(-1/2)

    return *this;
    }

LFHTEMP template<class O,unsigned int SIZE, Tuple_flag Cflag> auto Trianglix<C, 0u>::operator*(const Tuple<O,SIZE, Cflag> & a)const
 -> Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL>{
    Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL> fout;
    uint32_t i,k;
    uint32_t damax = a.getSize() < t_size ?  a.getSize(): t_size;
    for(k=0;k<damax;k++){
        const C* currow = data + ((k * (k+1u)) >> 1);
        fout[k] = currow[k] * a[k];
        for(i=0;i<k;i++) {
            fout[k] += currow[i] * a[i];
            fout[i] += currow[i] * a[k];
        }
    }
return fout;}



LFHTEMP template<class D, unsigned int SX, unsigned int SY> auto Trianglix<C, 0u>::mkInnerMult(const TMatrix<D,SX,SY> &other) const
 -> Trianglix<decltype(other.data[0] * this->data[0]* other.data[0] ), 0u> {
    Trianglix<decltype(other.data[0] * this->data[0]* other.data[0] ), 0u> fout;
    int maxshared = (other.nbRows() > t_size) ? t_size : other.nbRows();
    fout.toSize(other.nbCols()).toZero();
    TMatrix<decltype(other.data[0] * this->data[0]), 0,0> half; half.setSizes(other.nbCols(), maxshared);
    // multiply with lower triangular
    uint32_t i,j,k,l;
    for(i=0,k=0;i < maxshared;i++,k++){
        C tmp = data[k+i] * 0.5;
        for(l=0;l < other.nbCols(); l++) half.data[l+ i * other.nbCols()] = tmp * other.data[l * other.nbRows() + i];
        for(j=0;j<i;j++,k++){
            for(l=0;l < other.nbCols(); l++) half.data[l + i * other.nbCols()] += data[k] * other.data[l * other.nbRows() + j];
        }
    }
    for(i=0,k=0;i < other.nbCols();i++){
        for(j=0;j<i;j++,k++){
            for(l=0;l < maxshared; l++) fout.data[k] += half.data[j+l * other.nbCols()] * other.data[l + i * other.nbRows()] + half.data[i+l *other.nbCols()] * other.data[l+j * other.nbRows()];
        }
        for(l=0;l < maxshared; l++) fout.data[k] += half.data[i+l * other.nbCols()] * other.data[l+i * other.nbRows()];
        fout.data[k++] *= 2.0;
    }
return fout;}

LFHTEMP template<class D> auto Trianglix<C, 0u>::mkInnerMult(const SparseMatrix<D> &other) const
 -> Trianglix<decltype(other.data[0] * this->data[0]* other.data[0] ), 0u> {
    Trianglix<decltype(other.data[0] * this->data[0]* other.data[0] ), 0u> fout;
    /*int maxshared = (other.nbRows() > t_size) ? t_size : other.nbRows();
    fout.toSize(other.nbCols()).toZero();
    SparseMatrix<decltype(other.data[0] * this->data[0]), 0,0> half; half.toSizes(other.nbCols(), maxshared);
    // multiply with lower triangular
    uint32_t i,j,k,l;
    for(i=0,k=0;i < maxshared;i++,k++){
        C tmp = data[k+i] * 0.5;
        for(l=0;l < other.nbCols(); l++) half.data[l+ i * other.nbCols()] = tmp * other.data[l * other.nbRows() + i];
        for(j=0;j<i;j++,k++){
            for(l=0;l < other.nbCols(); l++) half.data[l + i * other.nbCols()] += data[k] * other.data[l * other.nbRows() + j];
        }
    }
    half.show();
    for(i=0,k=0;i < other.nbCols();i++){
        for(j=0;j<i;j++,k++){
            for(l=0;l < maxshared; l++) fout.data[k] += half.data[j+l * other.nbCols()] * other.data[l + i * other.nbRows()] + half.data[i+l *other.nbCols()] * other.data[l+j * other.nbRows()];
        }
        for(l=0;l < maxshared; l++) fout.data[k] += half.data[i+l * other.nbCols()] * other.data[l+i * other.nbRows()];
        fout.data[k++] *= 2.0;
    }*/
return fout;}


LFHTEMP template<class D, unsigned int SX, unsigned int SY> Trianglix<C, 0u> Trianglix<C, 0u>::mkOuterMult(const TMatrix<D,SX,SY> &other)const{
    int maxshared = (other.nbCols() > t_size) ? t_size : other.nbCols();
    Trianglix<C, 0u> fout; fout.setSize(other.nbRows()).toZero();

return fout;}
LFHTEMP template<class D, unsigned int SX, unsigned int SY> Trianglix<C, 0u> Trianglix<C, 0u>::dualDivision(const TMatrix<D,SX,SY> &other)const{
    Tuple<uint32_t,2u> osizes = other.getSizes();
}
LFHTEMP template<class D> Trianglix<C, 0u> Trianglix<C, 0u>::dualDivision(const Tuple<D*> &other, uint32_t lenght)const{
    Trianglix<C, 0u> fout; fout.setSize(other.getSize());
    //TODO
    return fout;
}
LFHTEMP template<class D, Tuple_flag Cflag> Tuple<C, 0u> Trianglix<C, 0u>::mkBackMult(const Tuple<D, 0u, Cflag> &other)const{Tuple<C, 0u> fout;
    uint32_t tt = other.getSize() < t_size ? other.getSize(): t_size;
    fout.setSize(tt);
    unsigned int i,j,k;
	if (tt == 0) return fout;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;j<tt;j++) {
		fout[j] = other[0] * data[k];
		fout[0] += other[j] * ExOp::mkTrju(data[k]);
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mkTrju(data[k]);}
        fout[j] += other[j] * data[k++];
	}
return fout;}
LFHTEMP template<class D, Tuple_flag Cflag> Tuple<C, 0u> Trianglix<C, 0u>::mkBackDivi(const Tuple<D, 0u, Cflag> &input)const{Tuple<C, 0u> fout;
    if (input.getSize() < 3){
        if (input.getSize() == 1) {fout.setSize(1); fout[0] = input[0] / this->data[0];}
        else if (input.getSize() == 2) {fout.setSize(2);
            C denum = ExOp::mkInverse(data[0] * data[2] - ExOp::mkTrjuProd(data[1]));
            fout[0] = denum * (input[0] * this->data[2] - input[1] * this->data[1]);
            fout[1] = denum * (input[1] * this->data[0] - input[0] * ExOp::mkTrju(this->data[1]));
        }
    }else{
        fout = input;
        HHscope hh(*this);
        hh.runFindInverse(*this);
        uint32_t j,k;
        for(j=t_size-1;j>1;j--) fout.HouseHolderMultiply(hh.hv+((j * (j+1))/2), hh.normbuf[j-2], j+1);
        k=0;
        fout[0] *= hh.tribuf[k++];
        for(j=1;j<t_size;j++){
            fout[j] -= fout[j-1] * hh.tribuf[k++];
            fout[j] *= hh.tribuf[k++];
            }
        for(j--;j!=0;j--){
            fout[j-1] -= fout[j] * hh.tribuf[k++];
        }
        for(j=2;j<t_size;j++) fout.HouseHolderMultiply(hh.hv+((j * (j+1))/2), hh.normbuf[j-2], j);
    }
return fout;}
LFHTEMP template<class D, Tuple_flag Cflag> Tuple<C, 0u> Trianglix<C, 0u>::mkInvDivi(const Tuple<D, 0u, Cflag> &input)const{Tuple<C, 0u> fout;
    fout.setSize(this->getSize());
    HHscope hh(*this);
    hh.runFindInverse(*this);



return fout;}




/*
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::operator*(const Matrix<C>& other) const{ // does not work with com
    Trianglix<C, 0u> fout; if (other.sizey != t_size) exit(1); fout.setSize(other.sizex);fout.toZero();
    C* s2 = new C[other.sizex];
    unsigned int i,j,k,l;
        for(j=0,k=0;j<t_size;j++) {

        if (j>0) {for(l=0;l<other.sizex;l++) s2[l] = data[k] * other.data[l]; k++;}
        else for(l=0;l<other.sizex;l++) ExOp::toZero(s2[l]);
        for(i=1;i<j;i++,k++) for(l=0;l<other.sizex;l++) s2[l] += data[k] * other.data[l + i * other.sizex];
        for(l=0;l<other.sizex;l++) s2[l] += ExOp::mkTrju(s2[l]) +data[k] * other.data[l + j * other.sizex]; // WORKS FOR REAL ONLY!
    //    for(l=0;l<other.sizex;l++) s2[l] = ExOp::mkTrju(s2[l]) + s2[l] + data[k] * other.data[l + j * other.sizex];
        k++;
            fout.data[0] += s2[0] * other.data[0 + j * other.sizex];
            fout.data[1] += s2[1] * other.data[0 + j * other.sizex] + s2[0] *  other.data[1 + j * other.sizex]; // a*c - b*d
            fout.data[2] += s2[1] * other.data[1 + j * other.sizex];
            fout.data[5] += s2[2] * other.data[2 + j * other.sizex];
            fout.data[6] += s2[3] * other.data[0 + j * other.sizex] + s2[0] *  other.data[3 + j * other.sizex]; // a*c - b*d
            fout.data[7] += s2[3] * other.data[1 + j * other.sizex] + s2[1] *  other.data[3 + j * other.sizex]; // a*c - b*d
        }
    delete[](s2);
    return fout;
}*/




LFHTEMP TMatrix<C,0u,0u,(Tuple_flag)0u>  Trianglix<C, 0u>::getSimulDiagonalizer(const Trianglix<C, 0u>& other, unsigned int osize) const { TMatrix<C,0u,0u,(Tuple_flag)0u> fout;
    if (osize == 0) osize = t_size;
    fout.setSizes(osize, t_size);
    ExOp::toOne(fout);

    return fout;
}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::inverse_OLD() const{

    Trianglix<C, 0u> fout(*this);

    if (t_size >2){
    C* buf = new C[t_size*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[t_size *(t_size+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[t_size-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        fout.data[0] = ExOp::mkInverse(fout.data[0]);
        buf[t_size] = fout.data[0];
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkInverse(fout.data[k] - ExOp::mkTrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }


    for(j=t_size-2;j<t_size;j--) fout.offdiagelimination_up( -buf[t_size+j],j,buf);

    for(j=2;j<t_size;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkInverse(data[0]);
    return(fout);



}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::inverse_MK2() const{
    Trianglix<C, 0u> fout(*this);

     Trianglix<C, 0u> backup(*this);
    if (t_size >3){
    C* buf = new C[t_size*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[t_size *(t_size+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    C* bufRQ = new C[t_size*2]; LFH_NICE_ALLOCERROR(bufRQ,"")
    double normnorm;
    double* normbuf = new double[t_size-2];  LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }


        bufRQ[0] = fout.data[0];
        bufRQ[1] = fout.data[1];


        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;

            normnorm = pow(ExOp::pnorm(bufRQ[((j-1)<<1)]) + ExOp::pnorm(fout.data[k-1] ), -0.5f);

            if (j > 1) fout.data[k-j-2] = bufRQ[((j-2)<<1) | 1 ] * (1.0f / normnorm);
        //    printf("%e %e %e\n", bufRQ[((j-1)<<1)],fout.data[k-1] , normnorm );
            bufRQ[((j-1)<<1)] *= normnorm;

         //   bufRQ[(j<<1)] = fout.data[k+j+1] * bufRQ[((j-1)<<1)];
          //  bufRQ[(j<<1)] = fout.data[k]*bufRQ[((j-1)<<1)]  ;
            if (j+1 < t_size) bufRQ[(j<<1) | 1] =  bufRQ[((j-1)<<1)] *fout.data[k+j+1];
            bufRQ[(j<<1)] = bufRQ[((j-1)<<1)] * fout.data[k];
            bufRQ[(j<<1)] -= fout.data[k-1] *bufRQ[((j-1)<<1) | 1 ] * normnorm;

         //   bufRQ[(j<<1)] -= bufRQ[((j-1)<<1) | 1]  *bufRQ[((j-1)<<1) | 1];

            //printf("partial x=%e y=%e  %e\n", bufRQ[(j<<1)], bufRQ[(j<<1) | 1] , fout.data[k-1] *bufRQ[((j-1)<<1) | 1 ] * normnorm );
            bufRQ[((j-1)<<1)| 1] = fout.data[k-1] * normnorm;

            }
        fout.data[k-1] = bufRQ[((j-2)<<1) | 1 ] * bufRQ[((j-1)<<1)];
     //   fout.data[k] =bufRQ[((j-1)<<1)];
       /* backup.show();
        printf("first rotation! %e %e\n", bufRQ[0] ,bufRQ[1] );
        backup.QR_back(bufRQ[0], bufRQ[1], 0);backup.show();
        printf("second rotation! %e %e\n", bufRQ[2] ,bufRQ[3] );
        backup.QR_back(bufRQ[2], bufRQ[3], 1);backup.show();
        printf("third rotation! %e %e\n", bufRQ[4] ,bufRQ[5] );
        backup.QR_back(bufRQ[4], bufRQ[5], 2);backup.show();
        fout.show();*/

        fout.data[0] = ExOp::mkInverse(fout.data[0]);
        buf[t_size] = fout.data[0];
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkInverse(fout.data[k] - ExOp::mkTrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }

    for(j=t_size-2;j<t_size;j--) {fout.offdiagelimination_up_backroutine( -buf[t_size+j],j,buf); if (j == t_size-2) fout.show();}



    for(j=2;j<t_size;j++) {fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);} // pts55

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 3){
        C denum = ExOp::mkInverse(data[0] * (data[4] * ExOp::mkTrju(data[4])- data[2] * data[5]) + ExOp::mkTrju(data[1]) * (data[1]*data[5] - ExOp::mkTrju(data[4])* data[3]) + ExOp::mkTrju(data[3])* ( data[2]*data[3] - data[1]*data[4]) );
        fout.data[0] = (data[2] * data[5] - data[4] * ExOp::mkTrju(data[4])) * denum;
        fout.data[2] = (data[0] * data[5] - data[3] * ExOp::mkTrju(data[3])) * denum;
        fout.data[5] = (data[0] * data[2] - data[1] * ExOp::mkTrju(data[1])) * denum;
        fout.data[1] = (ExOp::mkTrju(data[1]) * data[5] - data[4] * ExOp::mkTrju(data[3])) * denum;
        fout.data[3] = (ExOp::mkTrju(data[3]) * data[2] - data[1] * ExOp::mkTrju(data[4])) * denum;
        fout.data[4] = (ExOp::mkTrju(data[4]) * data[0] - data[1] * ExOp::mkTrju(data[3])) * denum;
    }else if (t_size == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkInverse(data[0]);
    return(fout);
}
LFHTEMP void Trianglix<C, 0u>::wrInverse_block(Trianglix<C, 0u>& fout, uint32_t _size) const{
    if (_size ==0) _size = t_size;
    fout.setSize(_size);
    if (_size < 3){
        if (_size < 2){
            if (_size == 1) fout.data[0] = ExOp::mkInverse(data[0]);
        }else{
            C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
            fout.data[0] = data[2] * denum;
            fout.data[1] = -data[1] * denum;
            fout.data[2] = data[0] * denum;
        }
        return;
    }
    Trianglix<C, 0u> matA,matB;

    //wrInverse_block(matA, _size >> 1);
    Tuple<C*> matrix; matrix.setSize( (_size+1) >> 1);
    for(int i= (_size >> 1) ; i < _size;i++) matrix[i - (_size >> 1)] = data + ((i * (i*1)) >> 1);
    matB = this->dualDivision(matrix, _size >> 1);
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkInverse_older() const{
    Trianglix<C, 0u> fout(*this);
    Trianglix<C, 0u> resid;

    if (t_size >2){
    C* buf = new C[t_size*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[t_size *(t_size+1)/2];  LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[t_size-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    C tmp;

    unsigned int i,j,k,l;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        // interative method for the

        Tuple<C,0u> tridiago; tridiago.setSize(t_size*2-1);
        Tuple<C,0u> x_tridiago; x_tridiago.setSize(t_size-1);
        Tuple<C,0u> refine_buf; refine_buf.setSize(t_size);
        tridiago[0] = fout.data[0];
        x_tridiago[0] = fout.data[1] / fout.data[0];
        j=0;
        while(true){
            k = ((j+1)*(j+4))/2;
            tridiago[(j<<1)| 1] = fout.data[k-1];
            j++;
            tridiago[(j<<1)] = fout.data[k];
            if (j == t_size-1) break;
            x_tridiago[j] = fout.data[k-1] / (fout.data[k] - fout.data[k-1] * x_tridiago[j-1]);
        }


        fout.data[0] = ExOp::mkInverse(fout.data[0]);
        buf[t_size] = fout.data[0];
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkInverse(fout.data[k] - ExOp::mkTrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }


    for(j=t_size-2;j<t_size;j--) fout.offdiagelimination_up( -buf[t_size+j],j,buf);

    // iterative refinement for tri-diagonal matrices!

 //   k = ((t_size-1) *t_size)/2;
 //   ExOp::show(fout.data[k] * tridiago[0] + fout.data[k+1] * tridiago[1]);
 //   for(i=2,k++;(i>>1) < t_size-1;i+=2, k++){
 //       ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]);
 //       }
 //   ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i]);


//   printf("before\n");   (Matrix<double>(fout) * tmptmp_triangle).show();

    for(l=0;l<t_size;l++){ // makes this n^3 like the hols
    for(j = t_size-1; j > 0; j--){
    k = (j *(j+1))/2;

 //  ExOp::show(fout.data[k] * tridiago[0] + fout.data[k+1] * tridiago[1]);
    refine_buf[0] = -fout.data[k] - ((fout.data[k+1] * tridiago[1]) / tridiago[0]);
    for(i=2,k++;(i>>1) < j ;i+=2, k++){
    //  ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]);
        refine_buf[(i>>1)] = ( (fout.data[k-1] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]  ) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
	}
    // objective is 1
    tmp = -1.0f;
    if (j < t_size-1)  tmp += tridiago[i+1] * fout.data[k+j+1];

  //  ExOp::show( tmp  + fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i]);
    tmp = ( tmp  + (fout.data[k-1] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i]) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
    fout.data[k] += tmp;// printf("%e refine\n",tmp);
    for(i=(j << 1)-2,k--; i != 0xFFFFFFFE;i-=2, k--){
    tmp = refine_buf[(i>>1)] - tmp * x_tridiago[(i>>1)];
    fout.data[k] += tmp;// printf("%e refine\n",tmp);
    }
    }

    fout.data[0] = (1.0f -  fout.data[1] * tridiago[1]) / tridiago[0];
    }


/*
    k = ((t_size-1) *t_size)/2;
    ExOp::show(fout.data[0] * tridiago[0] + fout.data[1] * tridiago[1]);
   refine_buf[0] = -fout.data[k] - ((-1.0f + fout.data[k+1] * tridiago[1]) / tridiago[0]);
    j=1;
    for(i=2,k+=j;(i>>1) < t_size-1;i+=2, k+=j){
        ExOp::show(fout.data[k-j] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+j+1] * tridiago[i+1]);
        refine_buf[(i>>1)] = ( (fout.data[k-j] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+j+1] * tridiago[i+1]  ) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);

        j++;
        }
    // objective is 1
    ExOp::show(fout.data[k-j] * tridiago[i-1] + fout.data[k] * tridiago[i]);
    tmp = ((fout.data[k-j] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i]) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
    fout.data[k] += tmp; printf("%e refine\n",tmp);
    for(i-=2,k-=j; i != 0xFFFFFFFE;i-=2, k-=j){
    tmp = refine_buf[(i>>1)] - tmp * x_tridiago[(i>>1)];
    fout.data[k] += tmp; printf("%e refine\n",tmp);
    }*/

 //   printf("after\n");   (Matrix<double>(fout) * tmptmp_triangle).show();

    for(j=2;j<t_size;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);

    resid = (*this) * fout;

    double errnorm =0;
    C one; ExOp::toOne(one);
    for(j=0,k=0;j<t_size;j++){
        resid.data[k] -= one;
        errnorm += ExOp::pnorm(resid.data[k]);
    }
    printf("inside id err %e\n", errnorm);
    resid.show();
    while(errnorm > 0.0000001f){
        for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            resid.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }


        for(j=2;j<t_size;j++) resid.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);
        fout -= resid;
        resid = (*this) * fout;
        for(j=0,k=0;j<t_size;j++){
            resid.data[k] -= one;
            errnorm += ExOp::pnorm(resid.data[k]);
        }
        printf("inside id err %e\n", errnorm);
        resid.show();
        errnorm =0.0f;
    }


    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkInverse(data[0]);
    return(fout);
}

LFHTEMP Trianglix<C, 0u>& Trianglix<C, 0u>::toRotate(int firstrow, double sine, double cosine){ // K <- QT*K*Q , where Q rotates row f, into f+1
    uint32_t i,j,k;
    k = (firstrow * (firstrow+1)) >> 1;
    C tmp;
    for(i=0;i< firstrow;i++){
        tmp = (data[k] * cosine) + (data[k+firstrow+1] * sine);
        data[k+firstrow+1] = (data[k] * sine) + (data[k+firstrow+1] * cosine);
        data[k++] = tmp;
    }

    j = k + firstrow+1;
    tmp = data[j] + ExOp::mkTrju(data[j]);
    data[j] = data[j] * (cosine*cosine) - ExOp::mkTrju(data[j]) * (sine * sine);
    data[j] += (data[j+1] - data[k]) * (sine * cosine);

    C tmpB = (data[k] * cosine + tmp * sine) * cosine + data[j+1] * sine * sine;
    data[j+1] = (data[k] * sine + tmp * cosine) * sine + data[j+1] * cosine * cosine;
    data[k] = tmpB;

    k = j + firstrow + 2;
    i +=2;
    for(;i< t_size;i++){
        tmp = (data[k] * cosine) + (data[k+1] * sine);
        data[k+1] = (data[k] * sine) + (data[k+1] * cosine);
        data[k] = tmp;
        k+= i+1;
    }

    return *this;
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkOuterInverseProduct(const Tuple<C*> outer) const{ Trianglix<C, 0u> fout;
    fout.setSize(outer.setSize());
    int i,j,k;
    C tmp;
    switch(t_size){
        case 0: return fout;
        case 1:
            tmp = ExOp::mkInverse(data[0]);
            for(i=0,k=0;i<outer.getSize();i++){
                for(j=0;j<=i;j++,k++) fout.data[k] = (outer[i][0] * ExOp::mkTrju(outer[j][0])) * tmp;
           //     fout.data[k++] = (outer[i][0] * ExOp::mkTrju(outer[i][0])) * tmp;
            }
        return fout;
        case 2: // todo
            tmp = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
            for(i=0,k=0;i<outer.getSize();i++){
                for(j=0;j<=i;j++,k++) fout.data[k] = ((outer[i][0] * ExOp::mkTrju(outer[j][0])) * data[2] + (outer[i][1] * ExOp::mkTrju(outer[j][1])) * data[0] - (outer[i][0] * ExOp::mkTrju(outer[j][1]) + ExOp::mkTrju(outer[j][0]) * outer[i][1]) * data[1]) * tmp;
              //  fout.data[k++] = ((outer[i][0] * ExOp::mkTrju(outer[i][0])) * data[2] + (outer[i][1] * ExOp::mkTrju(outer[i][1])) * data[0] - (outer[i][0] * ExOp::mkTrju(outer[i][1]) + ExOp::mkTrju(outer[i][0]) * outer[i][1]) * data[1]) * tmp;
            }
        return fout;
        default: // todo
            break;
    }
    return fout;
}

LFHTEMP void Trianglix<C, 0u>::toInmemInverse_routine(int start, int lenght){
    C* cur = data + ((start * (start+3))>>1);

    if (lenght < 4){
        C tmp[5];
        if (lenght == 2){
            C denum = ExOp::mkInverse(cur[0] * cur[start + 2] - ExOp::mkTrjuProd(cur[start + 1]));
            tmp[0] = cur[0] * denum;
            cur[0] = cur[start + 2] * denum;
            cur[start + 1] *= -denum;
            cur[start + 2] = tmp[0];
        }else{
            /*tmp[4] = ExOp::mkTrjuProd(data[4]);
            tmp[3] = ExOp::mkTrjuProd(data[3]);
            tmp[1] = ExOp::mkTrjuProd(data[1]);
            C denum = ExOp::mkInverse(cur[0] * tmp[0] - cur[start + 2] * data[5]) + ExOp::mkTrju(cur[start + 1]) * (cur[start + 1]*data[5] - ExOp::mkTrju(data[4])* data[3]) + ExOp::mkTrju(data[3])* ( cur[start + 2]*data[3] - cur[start + 1]*data[4]) );

            fout.data[1] = (ExOp::mkTrju(data[1]) * data[5] - data[4] * ExOp::mkTrju(data[3])) * denum;
            fout.data[3] = (ExOp::mkTrju(data[3]) * data[2] - data[1] * ExOp::mkTrju(data[4])) * denum;
            fout.data[4] = (ExOp::mkTrju(data[4]) * data[0] - data[1] * ExOp::mkTrju(data[3])) * denum;
            tmp[0] = (data[2] * data[5] - tmp[4] * denum;
            tmp[2] = (data[0] * data[5] - tmp[3]) * denum;
            fout.data[5] = (data[0] * data[2] - tmp[1]) * denum;
            cut[0] = tmp[0];
            cur[start+2] = tmp[2];*/
        }
    }else{


    }
};

LFHTEMP void Trianglix<C, 0u>::combineBlocks(const Trianglix<C, 0u> &xformed_D_inverse){
    // assumes blocks are currently [A^-1, B, [D - B^tA^-1B]^-1], fixes this to actual inverse
    uint32_t i,j,k;
    // need to multiply B twice, but multiply by A first
    // B <- B * A^-1
//    for(i=lenght;i < t_size;i++){


  //  }
}


LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkInverse_m2block() const{
    /*if (t_size < 5) return mkInverse();
    // minus 2 inverse...
    Trianglix<C, 0u> fout; fout.setSize(t_size);
    Trianglix<C, 0u> matA; matA.setSize(2);
    matA.data[0] = data[0];
    matA.data[1] = data[1];
    matA.data[2] = data[2];
    Tuple<C*> matB; matB.setSize(t_size-2);
    for(i=2;i<t_size;i++) matB[i-2] = data + ((i * (i+1)) >> 1);

    Trianglix<C, 0u> inverse1 = matA.mkOuterInverseProduct(matB);
    int i,j,k,pk;
    for(i=0,k=0;i< inverse1.getSize();i++){
        pk = 2 + (((i +2)*(i +3)) >> 1);
        for(j=0;j<=i;j++,k++) inverse1.data[k] =  data[pk + j] - inverse1.data[k];
    }
    Trianglix<C, 0u> inverse2 = inverse1.mkInverse();

    for(i=0,k=0;i< inverse1.getSize();i++){
        pk = 2 + (((i +2)*(i +3)) >> 1);
        for(j=0;j<=i;j++,k++) fout.data[pk+j] =  inverse2.data[k];
    }
*/

    /*int i = ((t_size -1) * (t_size +2))>> 1;
    matD.data[2] = data[i];
    matD.data[1] = data[i-1];
    matD.data[1] = data[((t_size -2) * (t_size +1))>> 1];*/

}


LFHTEMP C Trianglix<C, 0u>::maxEigenValue()const{ C daans; return daans;} // TODO!

/*LFHTEMP void Trianglix<C, 0u>::setSize(uint32_t _siz){
    tribuf(new C[_siz*3-2]),hv(new C[_siz *(_siz+1)/2]), normbuf(new double[_siz-2]);
}*/
LFHTEMP Trianglix<C, 0u>::HHscope::HHscope(uint32_t _siz): tribuf(new C[_siz*3-2]),hv(new C[_siz *(_siz+1)/2]), normbuf(new double[_siz-2]){}
LFHTEMP Trianglix<C, 0u>::HHscope::HHscope(const Trianglix<C, 0u> &target): tribuf(new C[target.t_size*3-2]),hv(new C[target.t_size *(target.t_size+1)/2]), normbuf(new double[target.t_size-2]){}


LFHTEMP Trianglix<C, 0u>::HHscope::~HHscope(){ delete[](tribuf); delete[](hv); delete[](normbuf);}
LFHTEMP void Trianglix<C, 0u>::HHscope::runFindInverse(const Trianglix<C, 0u> &target){
    Trianglix<C, 0u> invsc = target;
    unsigned int i,j,k;
    if (target.t_size == 1) tribuf[0] = ExOp::mkInverse(invsc.data[0]);
    else{
        for(j=target.t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, tribuf,true);
        }
        i=0;
        unsigned int l=target.t_size*3-3;
        tribuf[i] = ExOp::mkInverse(invsc.data[0]);
        for(j=1,k=0;;j++){
            tribuf[l] = invsc.data[k+j] * tribuf[i++];
            k = (j*(j+3))/2;
            tribuf[i++] = invsc.data[k-1];
            tribuf[i] = ExOp::mkInverse(invsc.data[k] - invsc.data[k-1]* tribuf[l--] );
            if (j+1 == target.t_size) break;
        }



    }
}
LFHTEMP void Trianglix<C, 0u>::HHscope::runFindEigen(const Trianglix<C, 0u> &target){
    Trianglix<C, 0u> invsc = target;
    unsigned int i,j,k;
    if (target.t_size <= 2) myexit("dont call this on small matrices");

    for(j=target.t_size-1;j>1;j--){
        k = (j * (j+1))/2;
        normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
        for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
        hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
        hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
        normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
        invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, tribuf,true);
    }
    i=0;
    unsigned int l=target.t_size*3-3;
    tribuf[i] = invsc.data[0];
    for(j=1,k=0;j < target.t_size;j++){
        k = (j*(j+3))/2;
        tribuf[i++] = invsc.data[k-1];
        tribuf[i++] = invsc.data[k];
    }
    this->show(target.t_size);
}

LFHTEMP double Trianglix<C, 0u>::HHscope::cptInvProjection(C* _out, const C* input, bool do_add_instead)const{
    printf("Tridiago:\n");

    return 0.0f;
}

LFHTEMP void Trianglix<C, 0u>::HHscope::show(int _size){
    printf("Tridiago:\n");
    uint32_t i;
    for(i = 0 ; i < _size*2-2;){
        ExOp::show(tribuf[i++],stdout,1);printf("\t");
        ExOp::show(tribuf[i++],stdout,1);printf("\n");
    }
    ExOp::show(tribuf[i++],stdout,1);printf("\n");
    return;
}

LFHTEMP void Trianglix<C, 0u>::wrSubTrianglix(Trianglix<C> &where, const Tuple<uint32_t> &dims)const{
    where.setSize(dims.getSize());
    uint32_t i,j,k,l;
    for(i=0,l=0;i<dims.getSize();i++){
        for(j=0;j<i;j++) {
            k = (dims[i] < dims[j]) ?  ((dims[j]* (dims[j]+1))>>1) + dims[i] : ((dims[i]* (dims[i]+1))>>1) + dims[j];
            where.data[l++] = data[k];
        }
        where.data[l++] = data[((dims[i]* (dims[i]+3))>>1)];
    }
}

LFHTEMP template<class F> TMatrix<C, 0u> Trianglix<C, 0u>::makeQDecomposition(F fnc, bool transpose) const{
    Tuple<C, 0u> lambdas;
    TMatrix<C, 0u> fout = this->makeQLambdaDecomposition(lambdas, transpose);
    uint32_t i,j,k;
    if (transpose){
        Tuple< decltype(fnc(lambdas[i])) , 0u> resbuf; resbuf.setSize(t_size);
        for(i = 0; i <  t_size; i++) resbuf[i] = fnc(lambdas[i]);
        for(i = 0,k=0; i <  t_size; i++){
            for(j=0;j<t_size;j++) fout.data[k++] *= resbuf[j];
        }
    }else{
        for(i = 0,k=0; i <  t_size; i++){
            auto result = fnc(lambdas[i]);
            for(j=0;j<t_size;j++) fout.data[k++] *= result;
        }
    }
return fout;}


#ifdef Rcpp_hpp
LFHTEMP template<class RCL> void Trianglix<C, 0u>::rdMatrix(const arma::Mat<RCL> &where){
    int i,j,k;
    this->setSize(where.n_rows);
    for(i=0,k=0;i<this->getSize();i++){
        for(j=0;j<=i;j++) data[k++] = where.at(i,j);
    }
}
LFHTEMP void Trianglix<C, 0u>::rdMatrix(const Rcpp::NumericMatrix &where){
    int i,j,k;
    this->setSize(where.nrow());
    for(i=0,k=0;i<this->getSize();i++){
        for(j=0;j<=i;j++) data[k++] = where(i,j);
    }
}
LFHTEMP template<class RCL> void Trianglix<C, 0u>::wrMatrix(arma::Mat<RCL> &where)const{
    where = arma::Mat<RCL>(this->getSize(),this->getSize());
    int i,j,k;
    for(i=0,k=0;i<this->getSize();i++){
        for(j=0;j<i;j++) {
            where.at(i,j) = data[k];
            where.at(j,i) = data[k++];
        }
        where.at(i,i) = data[k++];
    }
}

LFHTEMP template<class RCL> void Trianglix<C, 0u>::wrSubMatrix(arma::Mat<RCL> &where, Tuple<uint32_t> dims)const{
    where = arma::Mat<RCL>(dims.getSize(),dims.getSize());
    int i,j,k;
    for(i=0;i<dims.getSize();i++){
        for(j=0;j<i;j++) {
            k = (dims[i] < dims[j]) ?  ((dims[j]* (dims[j]+1))>>1) + dims[i] : ((dims[i]* (dims[i]+1))>>1) + dims[j];
            where.at(i,j) = data[k];
            where.at(j,i) = data[k];
        }
        where.at(i,i) = data[((dims[i]* (dims[i]+3))>>1)];
    }
}

LFHTEMP void Trianglix<C, 0u>::wrMatrix(Rcpp::NumericMatrix &where) const{
    where = Rcpp::NumericMatrix(this->getSize(),this->getSize());
    int i,j,k;
    for(i=0,k=0;i<this->getSize();i++){
        for(j=0;j<i;j++) {
            where(i,j) = data[k];
            where(j,i) = data[k++];
        }
        where(i,i) = data[k++];
    }
}

LFHTEMP template<class RCL> Trianglix<C, 0u>& Trianglix<C, 0u>::toMatrixRecomb(const arma::Mat<RCL> &ortho, const arma::Col<RCL> &eigen){
    this->setSize(ortho.n_rows);
    int i,j,k,l;
    for(i=0,k=0;i<ortho.n_rows;i++,k++){
        for(j=0;j<i;j++,k++) {
            data[k] = ortho.at(j,0) * ortho.at(i,0) * eigen.at(0);
            for(l=1;l< ortho.n_rows;l++) data[k] += ortho.at(j,l) * ortho.at(i,l) * eigen.at(l);
        }
        data[k] = ortho.at(i,0) * ortho.at(i,0) * eigen.at(0);
        for(l=1;l< ortho.n_rows;l++) data[k] += ortho.at(i,l) * ortho.at(i,l) * eigen.at(l);
    }
    return *this;
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkInverse() const{
    Trianglix<C, 0u> fout;
    arma::Mat<C> tmp;
    this->wrMatrix(tmp);
    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = 1.0 / eigval.at(i);
    return fout.toMatrixRecomb(eigvec, eigval);
    }

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbs() const{Trianglix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);
    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = fabs(eigval.at(i));

    return fout.toMatrixRecomb(eigvec, eigval);
}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbsInverse() const{Trianglix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);

    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = 1.0 / fabs(eigval.at(i));

    return fout.toMatrixRecomb(eigvec, eigval);
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbsoft(const double &k) const{Trianglix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);
    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = ExOp::mkAbsoft(eigval.at(i),k);

    return fout.toMatrixRecomb(eigvec, eigval);
}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbsoftInverse(const double &k) const{Trianglix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);

    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = ExOp::mkAbsoftInverse(eigval.at(i),k);

    return fout.toMatrixRecomb(eigvec, eigval);
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbhard(const double &k) const{Trianglix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);
    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = ExOp::mkAbhard(eigval.at(i),k);

    return fout.toMatrixRecomb(eigvec, eigval);
}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbhardInverse(const double &k) const{Trianglix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);

    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    for(int i =0;i< this->getSize();i++) eigval.at(i) = ExOp::mkAbhardInverse(eigval.at(i),k);

    return fout.toMatrixRecomb(eigvec, eigval);
}

LFHTEMP Trianglix<C,0u>& Trianglix<C,0u>::toEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc){
    C tmp[4];
    if (t_size < 3){
        if (t_size != 2) {
            if (t_size == 1) ExOp::execMemfnc(data[0],fnc);
        }else{
            //tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            ExOp::execMemfnc(tmp[3],fnc);
            tmp[2] = tmp[2] - tmp[1] ;
            ExOp::execMemfnc(tmp[2],fnc);
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        arma::Mat<C> tmp; this->wrMatrix(tmp);
        arma::Col<C> eigval;
        arma::Mat<C> eigvec;
        arma::eig_sym( eigval, eigvec, tmp );
        for(int i =0;i< this->getSize();i++) ExOp::execMemfnc(eigval.at(i),fnc);
        return this->toMatrixRecomb(eigvec, eigval);
    }
return *this;}
LFHTEMP template<class F> Trianglix<C,0u>& Trianglix<C,0u>::toEigenTransform(F fnc){
    C tmp[4];
    Tuple<double> eigbuf;
    uint32_t i;
    if (t_size < 3){
        if (t_size != 2) {
            if (t_size == 1) fnc(data[0]);
        }else{
            //tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            fnc(tmp[3]);
            tmp[2] = tmp[2] - tmp[1] ;
            fnc(tmp[2]);
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        arma::Mat<C> tmp; this->wrMatrix(tmp);
        arma::Col<C> eigval;
        arma::Mat<C> eigvec;
        arma::eig_sym( eigval, eigvec, tmp );
        for(int i =0;i< this->getSize();i++) fnc(eigval.at(i));
        return this->toMatrixRecomb(eigvec, eigval);
    }
return *this;}

LFHTEMP template<class F> Trianglix<C,0u>& Trianglix<C,0u>::toEigenTransformVec(F fnc){
    C tmp[4];
    Tuple<double> eigbuf;
    uint32_t i;
    if (t_size < 3){
        if (t_size != 2) {
            if (t_size == 1) {eigbuf.setSize(1); eigbuf[0] = data[0]; fnc(eigbuf); data[0] = eigbuf[0];}
        }else{
            //tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            eigbuf.setSize(2);
            eigbuf[0] = tmp[3];
            eigbuf[1] = tmp[2] - tmp[1];
            fnc(eigbuf);
            tmp[3] = eigbuf[0];
            tmp[2] = eigbuf[1];
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        arma::Mat<C> tmp; this->wrMatrix(tmp);
        arma::Col<C> eigval;
        arma::Mat<C> eigvec;
        arma::eig_sym( eigval, eigvec, tmp );
        eigbuf.setSize(t_size);
        for(i=0;i<t_size;i++) eigbuf[i] = eigval.at(i);
        fnc(eigbuf);
        for(i=0;i<t_size;i++) eigval.at(i) = eigbuf[i];
        return this->toMatrixRecomb(eigvec, eigval);
    }
return *this;}


/*
LFHTEMP void Trianglix<C,0u>::showthateigen(){
   arma::Mat<C> tmp; this->wrMatrix(tmp);
   arma::Col<C> eigval;
   arma::Mat<C> eigvec;
   arma::eig_sym( eigval, eigvec, tmp );
   for(unsigned int i=0;i< eigvec.n_rows;i++){
       printf("\n%e\t", eigval.at(i));
       for(uint32_t j=0;j< eigvec.n_rows;j++) {
           printf("\t%e", eigvec.at(i,j));
       }
   }
}*/

LFHTEMP Tuple<C, 0u> Trianglix<C, 0u>::getEigenValues() const{
Trianglix<C, 0u> fout(*this);
    Tuple<C, 0u> real_fout;
    real_fout.setSize(t_size);
    if (t_size >2){
        arma::Mat<C> tmp;
        this->wrMatrix(tmp);
        arma::Col<C> eigval;
        arma::Mat<C> eigvec;
        arma::eig_sym( eigval, eigvec, tmp );
        for(uint32_t i=0;i<real_fout.getSize();i++) real_fout[i] = eigval.at(i);
    }else if (t_size == 2) {
        real_fout[0] = fout.data[0] - fout.data[2];
        real_fout[1] = sqrt(real_fout[0] * real_fout[0] + (ExOp::mkTrju(fout.data[1]) * fout.data[1] * 4.0f) ) ;
        real_fout[0] = (fout.data[0] + fout.data[2] - real_fout[1]) *0.5f;
        real_fout[1] += real_fout[0];
    }else if (t_size == 1) real_fout[0] = data[0];
    return(real_fout);
}

LFHTEMP template<class D> SparseTuple<C> Trianglix<C, 0u>::mkInvDivi(const SparseTuple<D> &input) const{SparseTuple<C> fout;
    uint32_t k[3];
    if (input.getSize() < 3){
        if (input.getSize() == 1) fout[input.deref_key(0)] = ExOp::mkInverse(data[((input.deref_key(0) *(input.deref_key(0)+3))>>1)]) * input.deref(0);
        else if (input.getSize() == 2) {
            k[0] = (input.deref_key(0) * (input.deref_key(0)+ 3))>>1;
            k[1] = (input.deref_key(0) < input.deref_key(1)) ? ((input.deref_key(1) * (input.deref_key(1)+ 1))>>1) + input.deref_key(0) : ((input.deref_key(0) * (input.deref_key(0)+ 1))>>1) + input.deref_key(1);
            k[2] = (input.deref_key(1) * (input.deref_key(1)+ 3))>>1;
            C denum = ExOp::mkInverse(data[k[2]] * data[k[0]] - data[k[1]] * ExOp::mkTrju(data[k[1]]));
            fout[input.deref_key(0)] = (data[k[2]] * input.deref(0) + data[k[1]] * input.deref(1))/denum;
            fout[input.deref_key(1)] = (data[k[0]] * input.deref(1) + data[k[1]] * input.deref(0))/denum;
        }
    }else{
        Tuple<uint32_t> dadims; dadims.setSize(input.getSize());

        arma::Col<C> inputar = arma::Col<C>(input.getSize());
        for(k[0]=0;k[0]< input.getSize();k[0]++) {dadims[k[0]] = input.deref_key(k[0]); inputar.at(k[0]) = input.deref(k[0]);}
        arma::Mat<C> tmp;
        this->wrSubMatrix(tmp, dadims);
        arma::Col<C> armafout = arma::solve(tmp, inputar);
        for(k[0]=0;k[0]< input.getSize();k[0]++) fout[input.deref_key(k[0])] = armafout.at(k[0]);
    }
return fout;}


LFHTEMP TMatrix<C, 0u> Trianglix<C, 0u>::makeQLambdaDecomposition(Tuple<C, 0u> &out_l, bool transpose) const{ TMatrix<C, 0u> fout;
    arma::Mat<C> tmp; this->wrMatrix(tmp);
    arma::Col<C> eigval;
    arma::Mat<C> eigvec;
    arma::eig_sym( eigval, eigvec, tmp );
    out_l.setSize(this->getSize());
    for(int i =0;i< this->getSize();i++) out_l[i] = eigval.at(i);
    fout.setSizes(this->getSize(),this->getSize());
    uint32_t i,j,k;
    if (transpose){
        for(j =0,k=0;j< this->getSize();j++){
            for(int i =0;i< this->getSize() ;i++) fout.data[k++] = eigvec.at(j,i);
        }
    }else{
        for(j =0,k=0;j< this->getSize();j++){
            for(int i =0;i< this->getSize() ;i++) fout.data[k++] = eigvec.at(i,j);
        }
    }
return(fout);}

#else
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkInverse() const{
    Trianglix<C, 0u> fout(*this);
    Trianglix<C, 0u> resid;

    if (t_size >3){
    C* buf = new C[t_size*2];  LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[t_size *(t_size+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[t_size-2]; LFH_NICE_ALLOCERROR(normbuf,"")

    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        // interative method for the

        Tuple<C,0u> tridiago; tridiago.setSize(t_size*2-1);
        Tuple<C,0u> x_tridiago; x_tridiago.setSize(t_size-1);
        Tuple<C,0u> refine_buf; refine_buf.setSize(t_size);
        Vector<double*> rot_buffer;

/*
        C g,u,p;
        C s[3];
        int pr;

        pr =1;
        ExOp::toZero(s[0]);
        j = 0;
        // FORWARD!

        g = real_fout[0];
        p = g * g;
        s[2] = buf[0] / (p + buf[0]);
        s[pr] = p / (p + buf[0]);
        u = s[2] * (g + real_fout[1]);
        real_fout[0] = g + u;
        for(j=1;;j++){
            g = real_fout[j] - u;
            p = g * g / s[pr]; // (s[pr] == 0.0f) ? s[pr^1] * buf[j-1] :
            if (!ExOp::isValid(p)) p = s[pr^1] * buf[j-1];
            if (j + 1 == t_size) break;
            buf[j-1] = s[2] * (buf[j] + p);
            pr ^= 1;
            s[2] = buf[j] / (buf[j] + p);
            s[pr] = p / (buf[j] + p);
            u = s[2] * (g + real_fout[j+1]);
            real_fout[j] = g + u;
        }
        buf[j-1] = s[2] * p;
        real_fout[j] = g;*/


        /*tridiago[0] = fout.data[0];
        x_tridiago[0] = fout.data[1] / fout.data[0];
        j=0;
        while(true){
            k = ((j+1)*(j+4))/2;
            tridiago[(j<<1)| 1] = fout.data[k-1];
            j++;
            tridiago[(j<<1)] = fout.data[k];
            if (j == t_size-1) break;
            x_tridiago[j] = fout.data[k-1] / (fout.data[k] - fout.data[k-1] * x_tridiago[j-1]);
        }*/


        // DQH...HFH...HQT = I
        // QH...HFH...HQT = QH D^-1QT

        // D QH...HFH...HQ = QH...HRH...HQ

      //  rot_buffer.push_back(new double[2 * t_size -2];







        fout.data[0] = ExOp::mkInverse(fout.data[0]);
        buf[t_size] = fout.data[0];
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkInverse(fout.data[k] - ExOp::mkTrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }

    for(j=t_size-2;j<t_size;j--) fout.offdiagelimination_up( -buf[t_size+j],j,buf);







    // iterative refinement for tri-diagonal matrices!

 //   k = ((t_size-1) *t_size)/2;
 //   ExOp::show(fout.data[k] * tridiago[0] + fout.data[k+1] * tridiago[1]);
 //   for(i=2,k++;(i>>1) < t_size-1;i+=2, k++){
 //       ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]);
 //       }
 //   ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i]);


//   printf("before\n");   (Matrix<double>(fout) * tmptmp_triangle).show();
/*
    for(l=0;l<t_size;l++){ // makes this n^3 like the hols
    for(j = t_size-1; j > 0; j--){
    k = (j *(j+1))/2;

 //  ExOp::show(fout.data[k] * tridiago[0] + fout.data[k+1] * tridiago[1]);
    refine_buf[0] = -fout.data[k] - ((fout.data[k+1] * tridiago[1]) / tridiago[0]);
    for(i=2,k++;(i>>1) < j ;i+=2, k++){
    //  ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]);
        refine_buf[(i>>1)] = ( (fout.data[k-1] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]  ) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
	}
    // objective is 1
    tmp = -1.0f;
    if (j < t_size-1)  tmp += tridiago[i+1] * fout.data[k+j+1];

  //  ExOp::show( tmp  + fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i]);
    tmp = ( tmp  + (fout.data[k-1] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i]) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
    fout.data[k] += tmp;// printf("%e refine\n",tmp);
    for(i=(j << 1)-2,k--; i != 0xFFFFFFFE;i-=2, k--){
    tmp = refine_buf[(i>>1)] - tmp * x_tridiago[(i>>1)];
    fout.data[k] += tmp;// printf("%e refine\n",tmp);
    }
    }

    fout.data[0] = (1.0f -  fout.data[1] * tridiago[1]) / tridiago[0];
    }*/


/*
    k = ((t_size-1) *t_size)/2;
    ExOp::show(fout.data[0] * tridiago[0] + fout.data[1] * tridiago[1]);
   refine_buf[0] = -fout.data[k] - ((-1.0f + fout.data[k+1] * tridiago[1]) / tridiago[0]);
    j=1;
    for(i=2,k+=j;(i>>1) < t_size-1;i+=2, k+=j){
        ExOp::show(fout.data[k-j] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+j+1] * tridiago[i+1]);
        refine_buf[(i>>1)] = ( (fout.data[k-j] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+j+1] * tridiago[i+1]  ) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);

        j++;
        }
    // objective is 1
    ExOp::show(fout.data[k-j] * tridiago[i-1] + fout.data[k] * tridiago[i]);
    tmp = ((fout.data[k-j] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i]) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
    fout.data[k] += tmp; printf("%e refine\n",tmp);
    for(i-=2,k-=j; i != 0xFFFFFFFE;i-=2, k-=j){
    tmp = refine_buf[(i>>1)] - tmp * x_tridiago[(i>>1)];
    fout.data[k] += tmp; printf("%e refine\n",tmp);
    }*/

 //   printf("after\n");   (Matrix<double>(fout) * tmptmp_triangle).show();

    for(j=2;j<t_size;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);

    resid = (*this) * fout;

    double errnorm =0;
    C one; ExOp::toOne(one);
    for(j=0,k=0;j<t_size;j++){
        resid.data[k] -= one;
        errnorm += ExOp::pnorm(resid.data[k]);
        k += j+2;
    }
    while(errnorm > 0.0000001f){
        printf("inside id err %e\n", errnorm);
        Trianglix<C, 0u> target = this->inverse_MK2();
        (target * (*this)).show();

        resid.show();
        target = resid;
        for(j=t_size-1;j>1;j--){
            resid.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);
        }
        resid.show();


        resid.data[0] = ExOp::mkInverse(resid.data[0]);
        buf[t_size] = resid.data[0];
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= resid.data[k-1];
            resid.data[k] = ExOp::mkInverse(resid.data[k] - ExOp::mkTrju(resid.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(resid.data[k-1]);
            buf[t_size+j] = resid.data[k];
            }
        for(j=t_size-2;j<t_size;j--) resid.offdiagelimination_up( -buf[t_size+j],j,buf);

        for(j=2;j<t_size;j++) resid.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);
        resid.show();

        fout -= resid;

        printf("Objective\n");
        (target - (*this) * resid).show();
        resid = (*this) * fout;
        for(j=0,k=0;j<t_size;j++){
            resid.data[k] -= one;
            errnorm += ExOp::pnorm(resid.data[k]);
            k+= j+2;
        }
        printf("inside id err %e\n", errnorm);
        this->show();
        fout.show();
        errnorm =0.0f;
    }
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 3){ // 5.932993
        C denum = ExOp::mkInverse(data[0] * (data[4] * ExOp::mkTrju(data[4]) - data[2] * data[5]) + ExOp::mkTrju(data[1]) * (data[1]*data[5] - ExOp::mkTrju(data[4])* data[3]) + ExOp::mkTrju(data[3])* ( data[2]*data[3] - data[1]*data[4]) );
        fout.data[0] = (data[4] * ExOp::mkTrju(data[4]) - data[2] * data[5]) * denum;
        fout.data[2] = (data[3] * ExOp::mkTrju(data[3]) - data[0] * data[5]) * denum;
        fout.data[5] = (data[1] * ExOp::mkTrju(data[1]) - data[0] * data[2]) * denum;
        fout.data[1] = (ExOp::mkTrju(data[1]) * data[5] - data[4] * ExOp::mkTrju(data[3])) * denum;
        fout.data[3] = (ExOp::mkTrju(data[3]) * data[2] - data[1] * ExOp::mkTrju(data[4])) * denum;
        fout.data[4] = (ExOp::mkTrju(data[4]) * data[0] - data[1] * ExOp::mkTrju(data[3])) * denum;
    }else if (t_size == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkInverse(data[0]);
    return(fout);
}

LFHTEMP TMatrix<C, 0u> Trianglix<C, 0u>::makeQLambdaDecomposition(Tuple<C, 0u> &out_l, bool transpose) const{ TMatrix<C, 0u> fout;
    fout.setSizes(t_size, t_size);
    out_l.setSize(t_size);
    C tmp[3];
    if (t_size < 4){
        if (t_size < 2){
            if (t_size == 1) ExOp::toOne(fout.data[0]); out_l[0] = this->data[0];
        }else if (t_size == 2){
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            tmp[2] = (data[0] + data[2]) * 0.5;
            out_l[0] = tmp[2] + tmp[1];
            out_l[1] = tmp[2] - tmp[1];
            tmp[0] /= tmp[1] * 4.0;
            fout.data[0] = sqrt(tmp[0] - 0.5);
            fout.data[1] = sqrt(tmp[0] + 0.5);
            if (transpose) fout.data[1] = -fout.data[1];
            fout.data[2] = -fout.data[1];
            fout.data[3] = fout.data[0];
        }else{

        }
    }else{
        // TODO
    }

return fout;}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbsInverse() const{ Trianglix<C, 0u> fout; fout.setSize(this->getSize());
    if (this->getSize() < 4){
        if (this->getSize() < 2){
            if (this->getSize() == 1) fout.data[0] = ExOp::mkAbsInverse(data[0]);
        }else if (this->getSize() == 2) {
            // charastic : l^2 - l(a+b) + AB -C^2
            // abs(a+b/2 +- sqrt((a-b)^2 + 4C^2)
            // since a and b must be real, distriminant is non-negative
            C discr = ExOp::mkTrjuProd(data[2] - data[0]) + ExOp::mkTrjuProd(data[1] * 2.0);
            if (discr > ExOp::mkTrjuProd(data[2] + data[0])){
                ExOp::toInverse(discr);
                if (data[2] < -data[0]) ExOp::toNegative(discr);
                fout.data[0] = data[2] * discr;
                fout.data[1] = -data[1] * discr;
                fout.data[2] = data[0] * discr;
            }else{ // saddle!
                // TODO

            }
        }else{
            // TODO

        }
    }else {}// TODO
    return fout;
}

// log2(1 + 4^-x)
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbsoft(const double &k) const{ Trianglix<C, 0u> fout; fout.setSize(this->getSize()); return fout;}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbsoftInverse(const double & k) const{ Trianglix<C, 0u> fout; fout.setSize(this->getSize()); return fout;}
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbhard(const double &k) const{ Trianglix<C, 0u> fout; fout.setSize(this->getSize()); return fout;}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkAbhardInverse(const double & k) const{ Trianglix<C, 0u> fout; fout.setSize(this->getSize());
    if (this->getSize() < 4){
        if (this->getSize() < 2){
            if (this->getSize() == 1) fout.data[0] = ExOp::mkAbhardInverse(data[0],k);
        }else if (this->getSize() == 2) {
            // charastic : l^2 - l(a+b) + AB -C^2
            // abs(a+b/2 +- sqrt((a-b)^2 + 4C^2)
            // since a and b must be real, distriminant is non-negative
            C discr = ExOp::mkTrjuProd(data[2] - data[0]) + ExOp::mkTrjuProd(data[1] * 2.0);
            if (discr > ExOp::mkTrjuProd(data[2] + data[0])){
                ExOp::toInverse(discr);
                if (data[2] < -data[0]) ExOp::toNegative(discr);
                fout.data[0] = data[2] * discr;
                fout.data[1] = -data[1] * discr;
                fout.data[2] = data[0] * discr;
            }else{ // saddle!
                // TODO

            }
        }else{
            // TODO

        }
    }else {}// TODO
    return fout;
}
LFHTEMP Trianglix<C,0u>& Trianglix<C,0u>::toEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc){
    C tmp[5];
    if (t_size < 3){
        if (t_size != 2) {
            if (t_size == 1) ExOp::execMemfnc(data[0],fnc);
        }else{
            //tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            ExOp::execMemfnc(tmp[3],fnc);
            tmp[2] -= tmp[1];
            ExOp::execMemfnc(tmp[2],fnc);
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] *= tmp[3];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        // tri-digonalize, find eigenvalues, find QR rotations, transform eignevalues, recompose...
        uint32_t i,j,k,o;
        Trianglix<C, 0u> hh_scope = *this;
        double* normbuf = new double[hh_scope.t_size-2]; LFH_NICE_ALLOCERROR(normbuf,"")
        Trianglix<double> sinRot; sinRot.setSize(hh_scope.t_size-1);
        C* tridiago = new C[(hh_scope.t_size << 1)-1]; LFH_NICE_ALLOCERROR(tridiago,"")
        C* buf = new C[hh_scope.t_size]; LFH_NICE_ALLOCERROR(buf,"")
        C sum;
        C* vec;
        for(o=hh_scope.t_size-1;o>1;o--){
            vec = data+((o * (o+1))>>1);
            tridiago[(o<<1)] = vec[o];
            tridiago[(o<<1)-1] = vec[o-1];

            normbuf[o-2] = ExOp::pnorm(vec[0]);
            for(i=1;i<o-1;i++) {normbuf[o-2] += ExOp::pnorm(vec[i]);}
            vec[i] *= 1.0f - sqrt(1.0f + normbuf[o-2] / ExOp::pnorm(vec[i]));
            normbuf[o-2] += ExOp::pnorm(vec[i]);
            normbuf[o-2] *=0.5f;


            buf[0] = data[0] * vec[0];
            sum = data[0] * ExOp::mkTrju(vec[0]) * vec[0];
            for(j=1,k=1;j<o-1;j++) {
                buf[j] = data[k] * vec[0]; buf[0] += ExOp::mkTrju(data[k++]) * vec[j];
                for(i=1;i<j;i++) {buf[j] += data[k] * vec[i]; buf[i] += ExOp::mkTrju(data[k++]) * vec[j];}
                sum += data[k] * ExOp::mkTrju(vec[j]) * vec[j] + (buf[j] *  ExOp::mkTrju(vec[j]) + ExOp::mkTrju(buf[j]) * vec[j]);
                buf[j] += data[k] * vec[j];
                k++;
            }
            for(j=0;j<o-1;j++) buf[j] /= normbuf[o-2];
            sum /= normbuf[o-2] * normbuf[o-2];

            for(j=0,k=0;j<o-1;j++) {
                for(i=0;i<j;i++) data[k++] += sum * vec[j] *ExOp::mkTrju(vec[i]) - buf[j]*ExOp::mkTrju(vec[i]) - vec[j] * ExOp::mkTrju(buf[i]) ;
                data[k++] += sum * vec[j] *ExOp::mkTrju(vec[j]) - buf[j]*ExOp::mkTrju(vec[j]) - vec[j] * ExOp::mkTrju(buf[j]);
            }
        }
        for(o=0;o<3;o++) tridiago[o] = hh_scope.data[o];
        printf("HH scope:"); hh_scope.show();
        printf("tridiago\n");
        for(o=0;o<hh_scope.t_size-1;o++) printf("%e\t%e\n", tridiago[(o<<1)], tridiago[(o<<1) | 1] );
        printf("%e\n", tridiago[(o<<1)]);
        myalive();
        // finding eigen...
        for(o=hh_scope.t_size-1;o>1;o--){
            // init guess
            printf("%i-%i\n", (o<<1) -2, (o<<1));fflush(stdout);
            tmp[0] = (tridiago[o<<1] - tridiago[(o<<1) -2]) * 0.5;
            tmp[1] = tmp[0]*tmp[0] - tridiago[(o<<1)-1] * tridiago[(o<<1)-1];
            tmp[1] = sqrt(tmp[1]);
            printf("tmp %e\n", tmp[1]); fflush(stdout);

            // offdiag 1.931788  +- 0
            // -13.2780  -12.27894
            // 2.434205  0.8532858 5.7125089

            /*
            buf[o] = (tridiago[o<<1] + tridiago[(o<<1) -2]) + ((tmp[0] > 0) ? tmp[1] : -tmp[1]);
            do{
                printf("cur eigen: %e\n", buf[o]);
                for(i=0;i<o;i++){
                    sinRot.data[]

                }

            }while(tmp[0] != 0);*/


            // commit for eigen value... save rotations and update tridiago
            buf[o] = 15.7125089;
            k = (o * (o-1) >> 1);
            tridiago[0] -= buf[o];
            tmp[0] = tridiago[0];
            tmp[2] = tridiago[1];
            // angles : -0.9878495
            // 0.8263512
            for(i=0;i<o;i++){
                tmp[4] = tridiago[(i << 1) +1] * tridiago[(i << 1) +1] + tmp[i&1] * tmp[i&1];
                printf("here%i n2 = %e\n", i, tmp[4]); fflush(stdout);
                tmp[4] = sqrt(tmp[4]);
                printf("cosine %e\n", tmp[i&1] / tmp[4]);
                sinRot.data[k+i] = tridiago[(i << 1) +1] / tmp[4];
                tridiago[(i+1) << 1] -= buf[o];
                printf("ndcalc %e*%e %e %e\n", tridiago[(i+1) << 1], tmp[i&1], tmp[2|(i&1)],tmp[4]  );

                tmp[1-(i&1)] = (tridiago[(i+1) << 1] * fabs(tmp[i&1]) + tridiago[(i << 1)+1] * tmp[2|(i&1)] )/tmp[4] ;
                if (i+1 < o) tmp[3-(i&1)] = (tridiago[(i << 1)+3] * fabs(tmp[i&1]) ) / tmp[4] ;
                tmp[4] = fabs(tmp[i&1]) / tmp[4];
                if (i > 0){ // update prior off diago


                }
                printf("%e %e %e %e %e %e\n", tmp[i&1], tmp[4]*tmp[4], tridiago[(i+1)<<1] , sinRot.data[k+i] * sinRot.data[k+i], tmp[2 |(i&1)], 2.0 * tmp[4] * sinRot.data[k+i] );
                printf("%e %e %e\n", tmp[i&1]* tmp[4]*tmp[4], tridiago[(i+1)<<1] * sinRot.data[k+i] * sinRot.data[k+i], 2.0 * tmp[4] * sinRot.data[k+i] * tmp[2 |(i&1)] );
                tridiago[i << 1] = tmp[i&1] * tmp[4]*tmp[4] + ((i&1) ? -tridiago[(i+1)<<1]:tridiago[(i+1)<<1]) * sinRot.data[k+i] * sinRot.data[k+i] - 2.0 * tmp[4] * sinRot.data[k+i] * tmp[2 |(i&1)] + buf[o];


                printf("new dia %e %e\n", tmp[1-(i&1)], tridiago[i << 1] );
            }
            printf("tridiago\n");
            for(o=0;o<hh_scope.t_size-1;o++) printf("%e\t%e\n", tridiago[(o<<1)], tridiago[(o<<1) | 1] );
            printf("%e\n", tridiago[(o<<1)] );
        }
        this->toZero();

        tmp[0] = tridiago[0] - tridiago[2];
        tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + tridiago[1] * tridiago[1]);
        data[1] = tridiago[1] / tmp[1];
        tmp[2] = (tridiago[0] + tridiago[2]) * 0.5;
        buf[0] = tmp[2] + tmp[1];
        buf[1] = tmp[2] - tmp[1];

        printf("All eigenvals:");
        for(o=0;o<hh_scope.t_size;o++) {
            printf("%e\n",buf[o]);
            ExOp::execMemfnc(buf[o],fnc);
            data[((o * (o+3)) >>1)] = buf[o];
        }

        tmp[1] = (tmp[0] / tmp[1]) * 0.25;
        data[2] = (buf[0] + buf[1]) * 0.5;
        buf[0] -= buf[1];
        tmp[1] *= buf[0];
        data[0] = data[2] + tmp[1];
        data[2] -= tmp[1];
        data[1] *= buf[0] * 0.5;

        // apply rotations...


        // apply house holder...

        delete[](normbuf);
        delete[](tridiago);
        delete[](buf);
    }
return *this;}
LFHTEMP template<class F> Trianglix<C,0u>& Trianglix<C,0u>::toEigenTransform(F fnc){
    C tmp[4];
    if (t_size < 3){
        if (t_size != 2) {
            if (t_size == 1) fnc(data[0]);
        }else{
            tmp[0] = (data[0] - data[2]) * 0.5;
            tmp[1] = sqrt(tmp[0] * tmp[0] + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            fnc(tmp[3]);
            tmp[2] -= tmp[1] ;
            fnc(tmp[2]);
            tmp[1] = (tmp[0] / tmp[1]) * 0.5;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        HHscope hhscp(t_size);
        hhscp.runFindEigen(*this);

    }
return *this;}

LFHTEMP Tuple<C, 0u> Trianglix<C, 0u>::getEigenValues() const{
Trianglix<C, 0u> fout(*this);
    Tuple<C, 0u> real_fout; real_fout.setSize(t_size);
    if (t_size >2){
        C* buf = new C[t_size*2]; LFH_NICE_ALLOCERROR(buf,"")
        C* hv = new C[t_size *(t_size+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
        double* normbuf = new double[t_size-2];
        unsigned int i,j,k;
        for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }
        /*
        real_fout[0] = fout.data[0];
        buf[t_size] = ExOp::mkInverse(fout.data[0]);
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mkTrju(fout.data[k-1]) * buf[t_size+j-1];
            real_fout[j] = fout.data[k];
            buf[t_size+j] = ExOp::mkInverse(fout.data[k]);
        }*/
/*
        real_fout[0] = fout.data[0];
        for(j=0,k=1;j<t_size-1;){
            buf[j++] = fout.data[k] * ExOp::mkTrju(fout.data[k]);
            k++;
            real_fout[j] = fout.data[k];
            k += j+1;
        }*/


        real_fout[t_size-1] = fout.data[0];
        for(j=0,k=1;j<t_size-1;){
            buf[t_size-2 - (j++)] = fout.data[k] * ExOp::mkTrju(fout.data[k]);
            k++;
            real_fout[t_size-1 - j] = fout.data[k];
            k += j+1;
        }


        C g,u,p;
        C s[3];
        int pr;
        int loopy=0;
        double qual[2];
        do{
            for(loopy++;(loopy & 15) != 0;loopy++){
            //    printf("l%i\t%e\n",loopy , real_fout[0]);
            //    for(j=1;j<t_size;j++) printf("%e\t%e\n", buf[j-1], real_fout[j]);
                pr =1;
                ExOp::toZero(s[0]);
                j = 0;
                // FORWARD!

                g = real_fout[0];
                p = g * g;
                s[2] = buf[0] / (p + buf[0]);
                s[pr] = p / (p + buf[0]);
                u = s[2] * (g + real_fout[1]);
                real_fout[0] = g + u;
                for(j=1;;j++){
                    g = real_fout[j] - u;
                    p = g * g / s[pr]; // (s[pr] == 0.0f) ? s[pr^1] * buf[j-1] :
                    if (!ExOp::isValid(p)) p = s[pr^1] * buf[j-1];
                    if (j + 1 == t_size) break;
                    buf[j-1] = s[2] * (buf[j] + p);
                    pr ^= 1;
                    s[2] = buf[j] / (buf[j] + p);
                    s[pr] = p / (buf[j] + p);
                    u = s[2] * (g + real_fout[j+1]);
                    real_fout[j] = g + u;
                }
                buf[j-1] = s[2] * p;
                real_fout[j] = g;
                // BACKWARD, revert changes it seems...
    /*
                printf("l%i\t%e\n",loopy , real_fout[0]);
                for(j=1;j<t_size;j++) printf("%e\t%e\n", buf[j-1], real_fout[j]);
                j = t_size-1;
                p = g * g;
                ExOp::toZero(s[pr^1]);
                s[2] = buf[j-1] / (p + buf[j-1]);
                s[pr] = -p / (p + buf[j-1]);
                u = s[2] * (g + real_fout[j-1]);
                real_fout[j] = g + u;
                while(true){
                    j--;
                    g = real_fout[j] - u;
                    p = g * g / s[pr]; // (s[pr] == 0.0f) ? s[pr^1] * buf[j-1] :
                    if (!ExOp::isValid(p)) p = s[pr^1] * buf[j];
                    if (j == 0) break;
                    buf[j] = s[2] * (buf[j-1] - p);
                    pr ^= 1;
                    s[2] = buf[j-1] / (buf[j-1] - p);
                    s[pr] = p / (buf[j-1] - p);
                    u = s[2] * (g + real_fout[j-1]);
                    real_fout[j] = g + u;
                }
                buf[0] = s[2] * (-p);
                real_fout[0] = g;*/
            }
            qual[0] = log(ExOp::pnorm(real_fout[t_size-1]));
            qual[1] = 0.0f;
            for(j=0;j<t_size-1;j++){
                qual[0] += log(ExOp::pnorm(real_fout[j]));
                qual[1] += log(ExOp::pnorm(buf[j]));
            }
            qual[0] /= t_size;
            qual[1] /= t_size -1;
            //printf("%i: %e vs %e\n",loopy, qual[0], qual[1]);
        }while((loopy < 65536)&&(qual[0] < 32 + qual[1]));
        if (loopy >= 65536) fprintf(stderr, "Warning, could not achieve robust eigenvalue estimation!\n");

        delete[](hv);
        delete[](buf);
        delete[](normbuf);
    }/*else if (t_size == 3){ // TODO, maybe
        C cubic[3];
        cubic[2] = data[0] + data[2] + data[5];
        cubic[1] = (data[4] * ExOp::mkTrju(data[4])- data[2] * data[5]) + (data[3] * ExOp::mkTrju(data[3])- data[0] * data[5]) + (data[1] * ExOp::mkTrju(data[1])- data[0] * data[2]);
        cubic[0] = data[0] * (data[4] * ExOp::mkTrju(data[4])- data[2] * data[5]) + ExOp::mkTrju(data[1]) * (data[1]*data[5] - ExOp::mkTrju(data[4])* data[3]) + ExOp::mkTrju(data[3])* ( data[2]*data[3] - data[1]*data[4]);
        // need to find 3 roots for l^3 + cubic[2]l^2+ cubic[1]l+ cubic[0] =0
    }*/else if (t_size == 2) {
        real_fout[0] = fout.data[0] - fout.data[2];
        real_fout[1] = sqrt(real_fout[0] * real_fout[0] + (ExOp::mkTrju(fout.data[1]) * fout.data[1] * 4.0f) ) ;
        real_fout[0] = (fout.data[0] + fout.data[2] - real_fout[1]) *0.5f;
        real_fout[1] += real_fout[0];
    }else if (t_size == 1) real_fout[0] = data[0];
    return(real_fout);
}

LFHTEMP template<class D> SparseTuple<C> Trianglix<C, 0u>::mkInvDivi(const SparseTuple<D> &input) const{SparseTuple<C> fout;

return fout;}

LFHTEMP template<class F> Trianglix<C,0u>& Trianglix<C,0u>::toEigenTransformVec(F fnc){
return *this;}

#endif






/*
LFHTEMP C Trianglix<C, 0u>::log_determinant() const{
Trianglix<C, 0u> fout(*this);
   C log_prod;
    if (t_size >2){
    C* buf = new C[t_size*2];
    C* hv = new C[t_size *(t_size+1)/2];
    double* normbuf = new double[t_size-2];
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        log_prod = log(fout.data[0]);
        ExOp::toZero(hv[0]);
        if (fout.data[0] > hv[0]) hv[0] = fout.data[0];
        buf[t_size] = ExOp::mkInverse(fout.data[0]);
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mkTrju(fout.data[k-1]) * buf[t_size+j-1];
            log_prod  += log(fout.data[k]);
            if (fout.data[k] > hv[0]) hv[0] = fout.data[k];
            buf[t_size+j] = ExOp::mkInverse(fout.data[k]);
            }
        k =0;
        log_prod = (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        for(j=1;j<t_size;j++){k = (j*(j+3))/2;
        log_prod += (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        }

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2) return(log(data[0] *data[2] - data[1] * ExOp::mkTrju(data[1])));
    else if (t_size == 1) return(log(data[0]));
    else ExOp::toZero(log_prod);
    return(log_prod);
}*/
LFHTEMP double Trianglix<C, 0u>::bhattacharryya_partial(const Tuple<C, 0u> &dev ,const Trianglix<C, 0u>& other)const{
    double fout;
    // dev (P+P2/2)-1 dev
    C sum;
    C prod;
    Trianglix<C, 0u> invsc = (*this) + other;
    Tuple<C, 0u> xdev = dev;
    if (xdev.tup_size >2){
    C* buf = new C[xdev.tup_size*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[xdev.tup_size *(xdev.tup_size+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[xdev.tup_size-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=xdev.tup_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }

        buf[xdev.tup_size] = ExOp::mkInverse(invsc.data[0]);
        prod =  invsc.data[0];


        ExOp::toZero(hv[0]);
        if (invsc.data[0] > hv[0]) hv[0] = invsc.data[0];
        for(j=1;j<xdev.tup_size;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[xdev.tup_size+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[xdev.tup_size];
            else sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[xdev.tup_size+j-1];
            invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[xdev.tup_size+j-1]*invsc.data[k-1];
            if (invsc.data[k] > hv[0]) hv[0] = invsc.data[k];
            buf[xdev.tup_size+j] = ExOp::mkInverse(invsc.data[k]);
       //     sum += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[xdev.tup_size+j-1];
        prod = xdev.tup_size * log(0.5f);
        for(j=0;j<xdev.tup_size;j++) {k = (j*(j+3))/2;prod += (invsc.data[k] > 0.0f) ? log(invsc.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());}
        fout = 0.5f *( ExOp::norm(sum) + prod);


    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (xdev.tup_size == 2){
		fout = 0.5f; // todo

    }  else fout = 0.5f *(ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));


    return(fout);
    }
LFHTEMP C Trianglix<C, 0u>::Xformed_inner_product_of_inverse(const Tuple<C, 0u> &dev) const{
    double fout;


    // dev (P+P2/2)-1 dev
    C sum;
    Trianglix<C, 0u> invsc = (*this);
    Tuple<C, 0u> xdev = dev;
	if (getSize() < xdev.getSize()) exit(1);
	unsigned int siz = xdev.getSize();

    if (siz >2){
		//printf("alloc to die %i\n",siz); fflush(stdout);

    C* buf = new C[siz*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[siz *(siz+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[siz-2]; LFH_NICE_ALLOCERROR(normbuf,"")

    unsigned int i,j,k;
	//	printf("%i\t%t\n", invsc.getSize(), xdev.getSize()); fflush(stdout);
		for(j=siz-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }


        buf[xdev.tup_size] = ExOp::mkInverse(invsc.data[0]);


        for(j=1;j<siz;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[siz+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[siz];
            else sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[siz+j-1];
            invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[siz+j-1]*invsc.data[k-1];
            buf[siz+j] = ExOp::mkInverse(invsc.data[k]);
       //     sum += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[xdev.tup_size+j-1];
        fout = sum;


    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (siz == 2){
       fout = ExOp::mkTrju(dev[1]) * dev[0] * data[1];
       fout = (ExOp::mkTrju(dev[0]) * dev[0]* data[2] + ExOp::mkTrju(dev[1]) * dev[1]* data[0] - fout - ExOp::mkTrju(fout) ) * ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
    }  else fout = ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0]);


    return(fout);
    }
/*LFHTEMP void Trianglix<C, 0u>::eigen() const{
    Trianglix<C, 0u> fout(*this);
    if (getSize() >2){
    C* buf = new C[getSize()]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[getSize() *(getSize()+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[getSize()-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=getSize()-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf);
        }
        // tri-diago form!

        buf[0] = fout.data[0]; printf("eigen vals! %f", ExOp::mkrealproj(buf[0]) );
        for(j=1;j<getSize();j++){
            k = (j*(j+3))/2;
            buf[j] = data[k] - data[k-1] * ExOp::mkTrju(data[k-1]) / buf[j-1];printf("\t%f", ExOp::mkrealproj(buf[j]));
            }
        printf("\n");

        fout.offdiagelimination_down(-fout.data[1] / fout.data[0], 0, buf);


//     for(j=2;j<getSize();j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf);

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (getSize() == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (getSize() == 1) fout.data[0] = ExOp::mkInverse(data[0]);
    return(fout);
}*/
LFHTEMP void Trianglix<C, 0u>::show(FILE* f, int level)const{
    unsigned int i,ts,k,j;
    for(i=0,ts=totsize(),j=0,k=2;i<ts;i++) {ExOp::show(data[i],f,level+1); if (i==j) {fprintf(f,"\n"); j+= k++;} else fprintf(f,"\t");}
}
LFHTEMP string  Trianglix<C, 0u>::type_tostring() const{return string("Trianglix<") + ExOp::type_tostring(data[0]) +string(",0>");}
LFHTEMP C Trianglix<C, 0u>::log_determinant(uint32_t *fout_nb_negative_eigenvalues) const{

	Trianglix<C, 0u> fout(*this);
    C deter[2];
    C tmp;
    int safe_mag = 0;
    if (getSize() >2){
        C* buf = new C[getSize()*2]; LFH_NICE_ALLOCERROR(buf,"")
        C* hv = new C[getSize() *(getSize()+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
        double* normbuf = new double[getSize()-2]; LFH_NICE_ALLOCERROR(normbuf,"")
        unsigned int i,j,k;
        for(j=getSize()-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }
        delete[](hv);
        delete[](normbuf);

        buf[t_size*2-1] = fout.data[0];
        for(j=0,k=1;j<t_size-1;){
            buf[t_size-2 - (j++)] = fout.data[k] * ExOp::mkTrju(fout.data[k]);
            k++;
            buf[t_size*2-1 - j] = fout.data[k];
            k += j+1;
        }


        C g,u,p;
        C s[3];
        int pr;
        int loopy=0;
        double qual[2];
        do{
            for(loopy++;(loopy & 15) != 0;loopy++){
                //printf("l%i\t%e\n",loopy , real_fout[0]);
                //for(j=1;j<t_size;j++) printf("%e\t%e\n", buf[j-1], real_fout[j]);
                pr =1;
                ExOp::toZero(s[0]);
                j = 0;
                // FORWARD!

                g = buf[t_size];
                p = g * g;
                s[2] = buf[0] / (p + buf[0]);
                s[pr] = p / (p + buf[0]);
                u = s[2] * (g + buf[t_size+1]);
                buf[t_size] = g + u;
                for(j=1;;j++){
                    g = buf[t_size + j] - u;
                    p = g * g / s[pr]; // (s[pr] == 0.0f) ? s[pr^1] * buf[j-1] :
                    if (!ExOp::isValid(p)) p = s[pr^1] * buf[j-1];
                    if (j + 1 == t_size) break;
                    buf[j-1] = s[2] * (buf[j] + p);
                    pr ^= 1;
                    s[2] = buf[j] / (buf[j] + p);
                    s[pr] = p / (buf[j] + p);
                    u = s[2] * (g + buf[t_size+j+1]);
                    buf[t_size + j] = g + u;
                }
                buf[j-1] = s[2] * p;
                buf[t_size + j] = g;
                // BACKWARD, revert changes it seems...
    /*
                printf("l%i\t%e\n",loopy , real_fout[0]);
                for(j=1;j<t_size;j++) printf("%e\t%e\n", buf[j-1], real_fout[j]);
                j = t_size-1;
                p = g * g;
                ExOp::toZero(s[pr^1]);
                s[2] = buf[j-1] / (p + buf[j-1]);
                s[pr] = -p / (p + buf[j-1]);
                u = s[2] * (g + real_fout[j-1]);
                real_fout[j] = g + u;
                while(true){
                    j--;
                    g = real_fout[j] - u;
                    p = g * g / s[pr]; // (s[pr] == 0.0f) ? s[pr^1] * buf[j-1] :
                    if (!ExOp::isValid(p)) p = s[pr^1] * buf[j];
                    if (j == 0) break;
                    buf[j] = s[2] * (buf[j-1] - p);
                    pr ^= 1;
                    s[2] = buf[j-1] / (buf[j-1] - p);
                    s[pr] = p / (buf[j-1] - p);
                    u = s[2] * (g + real_fout[j-1]);
                    real_fout[j] = g + u;
                }
                buf[0] = s[2] * (-p);
                real_fout[0] = g;*/
            }
            qual[0] = log(ExOp::pnorm(buf[t_size*2-1]));
            qual[1] = 0.0f;
            for(j=0;j<t_size-1;j++){
                qual[0] += log(ExOp::pnorm(buf[t_size + j]));
                qual[1] += log(ExOp::pnorm(buf[j]));
            }
            qual[0] /= t_size;
            qual[1] /= t_size -1;
            //printf("%i: %e vs %e\n",loopy, qual[0], qual[1]);
        }while((loopy < 65536)&&(qual[0] < 32 + qual[1]));
        if (loopy >= 65536) fprintf(stderr, "Warning, could not achieve robust determinant estimation!\n");

        deter[0] = buf[t_size];
        ExOp::toOne(deter[1]);

        pr = 0;
        buf[getSize()] = ExOp::mkInverse(buf[t_size]);
        if (buf[t_size] < 0.0f) pr++;
        for(j=1;j<getSize();j++){
            if (buf[t_size+j] < 0.0f) pr++;
            tmp = deter[((j & 1)^1)] * buf[t_size+j];
            tmp -= deter[(j & 1)] * buf[j];
            if (!(ExOp::isValid(tmp))){ // overflow!
               deter[0] *= pow(0.5f, 300.0f);
               deter[1] *= pow(0.5f, 300.0f);
                tmp = deter[((j & 1)^1)] * buf[t_size+j];
                tmp -= deter[(j & 1)] * buf[j];
                safe_mag++;
            }else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
               deter[0] *= pow(2.0f, 300.0f);
               deter[1] *= pow(2.0f, 300.0f);
               tmp = deter[((j & 1)^1)] * buf[t_size+j];
               tmp -= deter[(j & 1)] * buf[j];
               safe_mag--;
            }
            deter[(j & 1)]  = tmp;

            }
        deter[0] = log(ExOp::norm(deter[(getSize() & 1)^1]));
        delete[](buf);
        if (fout_nb_negative_eigenvalues) fout_nb_negative_eigenvalues[0] = pr;
       return deter[0] + 300.0f * log(2.0f) * safe_mag;
    }else if (getSize() == 2) return(log(data[0] *data[2] - data[1] * ExOp::mkTrju(data[1])));
    else if (getSize() == 1) return(log(data[0]));
	else return(0.0f);
}

#undef LFHTEMP
#define LFHTEMP template <class C, unsigned int SIZE>

LFHTEMP void Trianglix<C, SIZE>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length, C* buf, bool hint ){
    // house holder multiplication, done twice!
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
	if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2))) return;
    unsigned int i,j,k;
    C s;

    buf[0] = data[0] * vec[0];
    s = data[0] * ExOp::mkTrju(vec[0]) * vec[0];
    for(j=1,k=1;j<length;j++) {
    buf[j] = data[k] * vec[0]; buf[0] += ExOp::mkTrju(data[k++]) * vec[j];
    for(i=1;i<j;i++) {buf[j] += data[k] * vec[i]; buf[i] += ExOp::mkTrju(data[k++]) * vec[j];}
    s += data[k] * ExOp::mkTrju(vec[j]) * vec[j] + (buf[j] *  ExOp::mkTrju(vec[j]) + ExOp::mkTrju(buf[j]) * vec[j]);
    buf[j] += data[k] * vec[j];
    k++;
    }
    if (!hint){
    for(;j<SIZE;j++) {
    buf[j] = data[k++] * vec[0];
    for(i=1;i<length;i++) buf[j] += data[k++] * vec[i];
        k+= j - length+1;
    }
    for(j=0;j<SIZE;j++) buf[j] /= denum2;
    }else for(j=0;j<length;j++) buf[j] /= denum2;
    s /= denum2 * denum2;

    for(j=0,k=0;j<length;j++) {
    for(i=0;i<j;i++) data[k++] += s * vec[j] *ExOp::mkTrju(vec[i]) - buf[j]*ExOp::mkTrju(vec[i]) - vec[j] * ExOp::mkTrju(buf[i]) ;
    data[k++] += s * vec[j] *ExOp::mkTrju(vec[j]) - buf[j]*ExOp::mkTrju(vec[j]) - vec[j] * ExOp::mkTrju(buf[j]);
    }
    if (hint) return;
    for(;j<SIZE;j++) {
        for(i=0;i<length;i++) data[k++] -= buf[j]*ExOp::mkTrju(vec[i]);
        k+= j - length+1;
    }
}
	LFHTEMP void Trianglix<C, SIZE>::CholeskyStep_up(const C * const vec, unsigned int length, C* buf){
		// = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
		unsigned int i,j,k,l;


		// computes the ( partial vec * T)

		buf[0] = data[0] * vec[0];
		for(i=1,k=1;i<SIZE;i++){
			buf[i] = data[k] * vec[0];
			if (i < length-1)	buf[0] += data[k] *  vec[i];
			else if (i == length -1) buf[0] += data[k] * (vec[length -1] - 1.0f);
			for(j=1,k++;j<i;j++,k++){
				if (j < length-1)	buf[i] += data[k] * vec[j];
				else if (j == length-1) buf[i] += data[k] * (vec[length -1] -1.0f);
				if (i < length-1)	buf[j] += data[k] * vec[i];
				else if (i == length-1) buf[j] += data[k] * (vec[length -1] -1.0f);
			}
			if (i < length-1)	buf[i] += data[k] * vec[i];
			else if (i == length-1) buf[i] += data[k] * (vec[length -1] -1.0f);
			k++;
		}

		//for(i=0;i<SIZE;i++) printf("%f\t", buf[i]);
		//printf("\n", buf[i]);



		k = (length *(length-1))/2;
		for(i=0;i<length-1;i++) data[k++] += 2.0f * buf[i];
		for(l=0;l< length-1;l++) data[k] += buf[l] * vec[l];
		data[k] += 2.0f * buf[i++] + buf[l] * (vec[l] -1.0f);
		for(k+=i+1;i<SIZE;k+=i+1) { data[k] += 2.0f * buf[i++];}
	}
LFHTEMP void Trianglix<C, SIZE>::offdiagelimination_down(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)"
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = (col *(col+1))/2;
    for(i=0;i<=col;i++) {buf[i] = fact * data[k+i];}
    for(;i<SIZE;i++) {k += i; buf[i] = fact * data[k+col];}
    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) data[k+i] += buf[i];
    k++;
    data[k+col] += buf[i] + ExOp::mkTrju(buf[i]) + fact * data[k-2] * ExOp::mkTrju(fact);
    for(i++;i<SIZE;i++) {k += i; data[k+col] += buf[i];}
}
LFHTEMP void Trianglix<C, SIZE>::offdiagelimination_up(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)"
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) buf[i] = fact * data[k+i];
    buf[i] = fact * data[k+i];
    k++;
    for(i++;i<SIZE;i++) {k += i; buf[i] = fact * data[k+col];}

    k = ((col) *(col+1))/2;
    for(i=0;i<col;i++) data[k+i] += buf[i];
    data[k+col] += ExOp::mkTrju(buf[col+1]) * fact + ExOp::mkTrju(buf[i]);
    for(i++;i<SIZE;i++) {k += i; data[k+col] += buf[i];}
}

LFHTEMP Trianglix<C, SIZE>::Trianglix(const Tuple<C, SIZE>& other){
    unsigned int i,j,k;
    for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(other[i]) * other[j];
}
LFHTEMP Trianglix<C, SIZE>::Trianglix(const Tuple<C, SIZE>& u, const Tuple<C, SIZE>& v){
    unsigned int i,j,k;
    for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(u[i]) * v[j] + ExOp::mkTrju(v[i]) * u[j];
}
LFHTEMP Trianglix<C, SIZE>& Trianglix<C, SIZE>::operator=(const Tuple<C, SIZE>& other){
    unsigned int i,j,k;
    for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mkTrju(other[i]) * other[j];
    return(*this);
}
LFHTEMP Tuple<C,SIZE> Trianglix<C, SIZE>::getEigenValues() const{
Trianglix<C, SIZE> fout(*this);
   Tuple<C,SIZE> real_fout;
    if (SIZE >2){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        real_fout[0] = fout.data[0];
        buf[SIZE] = ExOp::mkInverse(fout.data[0]);
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mkTrju(fout.data[k-1]) * buf[SIZE+j-1];
            real_fout[j] = fout.data[k];
            buf[SIZE+j] = ExOp::mkInverse(fout.data[k]);
            }

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2) {
    // P(x) = x2 - data[0] - data[2]x data[0] * data[2] - mktr[1] * data[1];
    // discr = (data[0] - data[1])^2 + mktr[1] * data[1]; hence, POSITIVE and real!
    real_fout[0] = fout.data[0] - fout.data[2];
    real_fout[1] = sqrt(real_fout[0] * real_fout[0] + (ExOp::mkTrju(fout.data[1]) * fout.data[1] * 4.0f) ) ;
    real_fout[0] = (fout.data[0] + fout.data[2] - real_fout[1]) *0.5f;
    real_fout[1] += real_fout[0];
    }else if (SIZE == 1) real_fout[0] = data[0];
    return(real_fout);

}


LFHTEMP C Trianglix<C, SIZE>::determinant() const{
Trianglix<C, SIZE> fout(*this);
    C deter[2];
    if (SIZE >2){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        deter[0] = fout.data[0];
        ExOp::toOne(deter[1]);

        buf[SIZE] = ExOp::mkInverse(fout.data[0]);
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            deter[(j & 1)] *= -(ExOp::mkTrju(fout.data[k-1]) * fout.data[k-1]);
            deter[(j & 1)] += deter[((j & 1)^1)] * fout.data[k];
            }

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
       return (SIZE & 1) ? -deter[0] : deter[1];
    }else if (SIZE == 2) return(data[0] *data[2] - data[1] * ExOp::mkTrju(data[1]));
    else if (SIZE == 1) return(data[0]);

}
LFHTEMP C Trianglix<C, SIZE>::log_determinant() const{
    Trianglix<C, SIZE> fout(*this);
    C deter[2];
    C tmp;
    int safe_mag = 0;
    if (SIZE >2){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }
    delete[](hv);
    delete[](buf);
    delete[](normbuf);

        deter[0] = fout.data[0];
        ExOp::toOne(deter[1]);

        buf[SIZE] = ExOp::mkInverse(fout.data[0]);
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            tmp = deter[((j & 1)^1)] * fout.data[k];
            tmp -= deter[(j & 1)] * (ExOp::mkTrju(fout.data[k-1]) * fout.data[k-1]);
            if (!(ExOp::isValid(tmp))){ // overflow!
               deter[0] *= pow(0.5f, 300.0f);
               deter[1] *= pow(0.5f, 300.0f);
                tmp = deter[((j & 1)^1)] * fout.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mkTrju(fout.data[k-1]) * fout.data[k-1]);
                safe_mag++;
            }else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
               deter[0] *= pow(2.0f, 300.0f);
               deter[1] *= pow(2.0f, 300.0f);
               tmp = deter[((j & 1)^1)] * fout.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mkTrju(fout.data[k-1]) * fout.data[k-1]);
               safe_mag--;
            }
            deter[(j & 1)]  = tmp;

            }
        deter[0] = log(ExOp::norm(deter[(SIZE & 1)^1]));


       return deter[0] + log(pow(2.0,300.0f)) * safe_mag;
    }else if (SIZE == 2) return(log(data[0] *data[2] - data[1] * ExOp::mkTrju(data[1])));
    else if (SIZE == 1) return(log(data[0]));
}
	LFHTEMP C Trianglix<C, SIZE>::maxEigenValue()const{
		Trianglix<C, SIZE> fout(*this);
		Tuple<C, SIZE> da_vec;
		C da_tmp[2];
		C da_norm;
		C da_dev;
		C da_old_norm;
		double tmptmp;
		if (SIZE >2){
			C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
			C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
			double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
			unsigned int i,j,k;
//			fout.show();
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
//				printf("da pre scale : %e\n", normbuf[j-2]);

				hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
				tmptmp = 1.0f -sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				if (ExOp::isValid(tmptmp)) {
					hv[k+i] *=  tmptmp;
					normbuf[j-2] += ExOp::pnorm(hv[k+i]);
					normbuf[j-2] *=0.5f;
				}else{ // hv[k+i] is too small, assume its zero!
//					printf("hello! has zero!\n");
					hv[k+i] = -sqrt(normbuf[j-2]);
				}

				ExOp::toZero(hv[k+j]);
//				printf("da scale : %e\n", normbuf[j-2]);
				fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);//fout.show();
			}
//			fout.show();
			delete[](normbuf);

			for(j=2;j<SIZE;j++){
				k = (j*(j+3))/2;
				fout.data[(j<<1)-1] = fout.data[k-1];
				fout.data[(j<<1)] = fout.data[k];
			}

			ExOp::toRand(da_vec);
//			fout.show();
			for(k=0;;k++){

			da_tmp[0] = da_vec[0];
			da_vec[0] = da_vec[0] * fout.data[0] + da_vec[1] * ExOp::mkTrju(fout.data[1]);
			da_tmp[1] = da_vec[0] * da_old_norm - da_tmp[0];// printf(" %e\t%e\n",da_tmp[0] , da_vec[0] );
			da_dev = da_tmp[1] * da_tmp[1];
 		    da_norm = da_vec[0]* da_vec[0];
			for(j=1;j<SIZE-1;j++){
				da_tmp[(j & 1)] = da_vec[j];
				da_vec[j] = da_tmp[((j & 1 ) ^ 1)]  * fout.data[(((j-1) << 1) | 1)]  + da_vec[j] * fout.data[(j<<1)] + da_vec[j+1] * ExOp::mkTrju(fout.data[((j<<1) | 1)]);
				da_tmp[((j & 1 ) ^ 1)] = da_tmp[(j & 1)] - da_vec[j]* da_old_norm; // printf(" %e\t%e\n",da_tmp[(j & 1)] , da_vec[j] );
				da_dev += da_tmp[((j & 1 ) ^ 1)] * da_tmp[((j & 1 ) ^ 1)];
				da_norm += da_vec[j]* da_vec[j];
			}
				da_tmp[(j & 1)] = da_vec[j];
				da_vec[j] = da_tmp[((j & 1 ) ^ 1)] * fout.data[(((j-1) << 1) | 1)] + da_vec[j] * fout.data[(j<<1)];
			//	printf(" %e\t%e\n",da_tmp[(j & 1)] , da_vec[j] );
				da_tmp[(j & 1)] -= da_vec[j] * da_old_norm;
				da_dev += da_tmp[(j & 1)] * da_tmp[(j & 1)];
				da_norm += da_vec[j]* da_vec[j];
			//	printf("devdev %e\t%e\n",da_norm * ExCo<double>::epsilon(), da_dev );
				if ((da_norm * ExCo<double>::epsilon()  > da_dev)||(k > SIZE*SIZE)) break;
				da_old_norm = ExOp::mkPowInvInt(da_norm, -2);
				da_vec *= da_old_norm; // ^ -0.5f

			}

			delete[](hv);
			delete[](buf);
			return ExOp::mkPowInvInt(da_norm, 2); // sqrt
		}else if (SIZE == 2){
			Tuple<C , 2> eigen = getEigenValues();
			return eigen[0] > eigen[1] ? eigen[0] : eigen[1];
		}  else return(fout.data[0]);
		return(da_norm);

	}
LFHTEMP  double Trianglix<C, SIZE>::pnorm() const {
	double fout=0; unsigned int i,j;
	C* k = data;
	for(j=0;j<SIZE;j++){
		for(i=0;i<j;i++,k++) fout += ExOp::pnorm(*k);
		fout += ExOp::pnorm(*k) * 0.5f;k++;
	}
	fout *= 2.0f;
	return fout;
}
	LFHTEMP template<unsigned int SIZ2> C Trianglix<C, SIZE>::Xformed_inner_product( const Tuple<C, SIZ2>& other) const{
		C s;ExOp::toZero(s);
		C s2;
		unsigned int i,j,k;

        for(j=0,k=0;(j<SIZ2) &&  (j<SIZE);j++) {
			if (j>0) s2 = data[k++] * other[0];
			else ExOp::toZero(s2);
			for(i=1;i<j;i++) s2 += data[k++] * other[i];
			// s += other[j] * ((ExOp::mkTrju(s2) + s2) + data[k++] * other[j]); WORKS FOR REAL ONLY!
			s +=  data[k++] * ExOp::mkTrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;

        }

		return s;
    }
	LFHTEMP C Trianglix<C, SIZE>::Xformed_inner_product( const Tuple<C, 0u>& other) const{
		C s;
		C s2;
		unsigned int i,j,k;
		if (other.size() == 0) {ExOp::toZero(s); return s;}
		s = data[0] * ExOp::mkTrju(other[0]) * other[0];
		if (ExOp::mkTrju(other[0]) * other[0] < 0.0f) { fprintf(stderr,"Should not be negative!: "); ExOp::show(ExOp::mkTrju(other[0]) * other[0],stderr); exit(1); }
		for(j=1,k=1;j<other.size();j++) {
			s2 = data[k++] * other[0];
			for(i=1;i<j;i++) s2 += data[k++] * other[i];
			s += 2.0f * other[j] * s2 + data[k++] * other[j]; // WORKS FOR REAL ONLY!
			// s +=  data[k++] * ExOp::mkTrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
			if (ExOp::mkTrju(other[j]) * other[j] < 0.0f) { fprintf(stderr,"Should not be negative!: "); ExOp::show(ExOp::mkTrju(other[j]) * other[j],stderr); exit(1); }
		}
		return s;
	}
/*
LFHTEMP C Trianglix<C, SIZE>::log_determinant() const{
Trianglix<C, SIZE> fout(*this);
   C log_prod;
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        ExOp::toZero(hv[0]);
        if (fout.data[0] > hv[0]) hv[0] = fout.data[0];
        buf[SIZE] = ExOp::mkInverse(fout.data[0]);
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mkTrju(fout.data[k-1]) * buf[SIZE+j-1];
            log_prod  += log(fout.data[k]);
            if (fout.data[k] > hv[0]) hv[0] = fout.data[k];
            buf[SIZE+j] = ExOp::mkInverse(fout.data[k]);
            }
        k =0;
        log_prod = (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        for(j=1;j<SIZE;j++){k = (j*(j+3))/2;
        log_prod += (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        }

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2) return(log(data[0] *data[2] - data[1] * ExOp::mkTrju(data[1])));
    else if (SIZE == 1) return(log(data[0]));
    return(log_prod);

} */

LFHTEMP Trianglix<C, SIZE> Trianglix<C, SIZE>::inverse() const{
    Trianglix<C, SIZE> fout(*this);

    if (SIZE >2){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        fout.data[0] = ExOp::mkInverse(fout.data[0]);
        buf[SIZE] = fout.data[0];

        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkInverse(fout.data[k] - ExOp::mkTrju(fout.data[k-1]) * buf[SIZE+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[SIZE+j] = fout.data[k];

            }


        for(j=SIZE-2;j<SIZE;j--){
            fout.offdiagelimination_up( -buf[SIZE+j],j,buf);
            }


    for(j=2;j<SIZE;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);

    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2){
        C denum = ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;

    }  else fout.data[0] =  ExOp::mkInverse(data[0]);
    return(fout);
}
// fout * fout' = inverse
LFHTEMP TMatrix<C,SIZE,SIZE> Trianglix<C, SIZE>::diagonalizer_of_inverse() const{ TMatrix<C,SIZE,SIZE> f_out;
    ExOp::toOne(f_out);

    Trianglix<C, SIZE> fout(*this);
    if (SIZE >1){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mkTrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mkTrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            f_out.LeftHouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        fout.data[0] = ExOp::mkInverse(fout.data[0]);
        buf[SIZE] = fout.data[0];

        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
             buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkInverse(fout.data[k] - ExOp::mkTrju(fout.data[k-1]) * buf[SIZE+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[SIZE+j] = fout.data[k];
            }


        for(j=SIZE-2;j<SIZE;j--){
            fout.offdiagelimination_up( -buf[SIZE+j],j,buf);
            }


    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else f_out.data[0] = pow(data[0], -0.5f);
    return(f_out);
    }
LFHTEMP double Trianglix<C, SIZE>::bhattacharryya_partial(const Tuple<C,SIZE> &dev ,const Trianglix<C,SIZE>& other)const{
    double fout;
    // dev (P+P2/2)-1 dev
    C sum,tmp;
    C deter[2];
    int safe_mag =0;
    Trianglix<C, SIZE> invsc = (*this) + other;
    Tuple<C,SIZE> xdev = dev;
    if (SIZE >1){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }



        deter[0] = invsc.data[0];
        ExOp::toOne(deter[1]);

         buf[SIZE] = ExOp::mkInverse(invsc.data[0]);
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;

            tmp = deter[((j & 1)^1)] * invsc.data[k];
            tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
            if (!(ExOp::isValid(tmp))){ // overflow!
               deter[0] *= pow(0.5f, 300.0f);
               deter[1] *= pow(0.5f, 300.0f);
                tmp = deter[((j & 1)^1)] * invsc.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
                safe_mag++;
            }else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
               deter[0] *= pow(2.0f, 300.0f);
               deter[1] *= pow(2.0f, 300.0f);
               tmp = deter[((j & 1)^1)] * invsc.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
               safe_mag--;
            }
            deter[(j & 1)]  = tmp;





            xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[SIZE];
            else sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
            invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            buf[SIZE+j] = ExOp::mkInverse(invsc.data[k]);
       //     sum += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];

        fout = 0.5f *( 0.5f * ExOp::norm(sum) + log(ExOp::norm(deter[((SIZE & 1)^1)])) + (300.0f * safe_mag - SIZE) *log(2.0f) );



    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else fout = 0.5f *(ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));


    return(fout);
    }
LFHTEMP C Trianglix<C, SIZE>::Xformed_inner_product_of_inverse(const Tuple<C, SIZE> &dev) const{
    double fout;
    // dev (P+P2/2)-1 dev
    C sum;
    C prod;
    Trianglix<C, SIZE> invsc = (*this);
    Tuple<C, SIZE> xdev = dev;
    if (SIZE >2){
    C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
    C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
    double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }

        buf[SIZE] = ExOp::mkInverse(invsc.data[0]);
        prod =  invsc.data[0];


        ExOp::toZero(hv[0]);
        if (invsc.data[0] > hv[0]) hv[0] = invsc.data[0];
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[SIZE];
            else sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
            invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            if (invsc.data[k] > hv[0]) hv[0] = invsc.data[k];
            buf[SIZE+j] = ExOp::mkInverse(invsc.data[k]);
       //     sum += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
        fout = sum;


    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2){
       fout = ExOp::mkTrju(dev[1]) * dev[0] * data[1];
       fout = (ExOp::mkTrju(dev[0]) * dev[0]* data[2] + ExOp::mkTrju(dev[1]) * dev[1]* data[0] - fout - ExOp::mkTrju(fout) ) * ExOp::mkInverse(data[0] * data[2] - data[1] * ExOp::mkTrju(data[1]));
    }  else fout = ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0]);


    return(fout);
    }
LFHTEMP C Trianglix<C, SIZE>::inv_Xformed_inner_product_singularguard(const Tuple<C,SIZE> &dev, double guard_fraction, double *ln_det) const{
        C fout;
		// dev (P+P2/2)-1 dev
        unsigned int count;
		C sum,tmp;
		C deter[2];
        C maxeigen; ExOp::toZero(maxeigen);
		int safe_mag =0;
		Trianglix<C, SIZE> invsc = (*this);
		Tuple<C,SIZE> xdev = dev;
		if (SIZE >1){
			C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
			C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
			double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
			unsigned int i,j,k;
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
				hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
				hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
				invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
				xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
			}



			deter[0] = invsc.data[0];
			ExOp::toOne(deter[1]);

			buf[SIZE] = invsc.data[0];
            maxeigen = fabs(buf[SIZE]);
			for(j=1;j<SIZE;j++){
				k = (j*(j+3))/2;

				/*
				tmp = deter[((j & 1)^1)] * invsc.data[k];
				tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
				if (!(ExOp::isValid(tmp))){ // overflow!
					deter[0] *= pow(0.5f, 300.0f);
					deter[1] *= pow(0.5f, 300.0f);
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
					safe_mag++;
				}else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
					deter[0] *= pow(2.0f, 300.0f);
					deter[1] *= pow(2.0f, 300.0f);
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
					safe_mag--;
				}
				deter[(j & 1)]  = tmp;*/

				xdev[j] -= ExOp::mkTrju(xdev[j-1]) * ExOp::mkInverse(buf[SIZE+j-1])*invsc.data[k-1];
				invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				buf[SIZE+j] = invsc.data[k];
                if (fabs(buf[SIZE+j]) > maxeigen) maxeigen = fabs(buf[SIZE+j]);
            }

            deter[0] = 0.0f;
            for(j=0;j<SIZE;j++){
                if (fabs(buf[SIZE+j]) > maxeigen * ExCo<C>::epsilon()){
                    deter[0] += log(fabs(buf[SIZE+j]));
                    count++;
                }
            }
            if (count == 0) count =SIZE;
            if ((guard_fraction == 0.0f) || (count == SIZE)){
				for(j=0;j<SIZE;j++) ExOp::toInverse(buf[SIZE+j]);
				if (ln_det != NULL) (*ln_det) = deter[0] * (((double)SIZE) / ((double)count));

			}else{

			for(j=0;j<SIZE;j++){
                if (fabs(buf[SIZE+j]) > maxeigen * ExCo<C>::epsilon()) ExOp::toInverse(buf[SIZE+j]);
                else buf[SIZE+j] = exp( -deter[0] / count) / guard_fraction;
            }
            if (ln_det != NULL) (*ln_det) = deter[0] * (((double)SIZE) / ((double)count));

			}
            fout = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[SIZE];
            for(j=1;j<SIZE;j++) fout += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];

		//	fout = (this_weight + other_weight) * 0.5f *( ExOp::norm(sum) + log(ExOp::norm(deter[((SIZE & 1)^1)])) + log(pow(2.0f,300.0f)) * safe_mag);



			delete[](hv);
			delete[](buf);
			delete[](normbuf);
		}else {
            fout = ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0]);
            if (ln_det != NULL) *ln_det = log(data[0]);
        }


		return(fout);
    }
LFHTEMP double Trianglix<C, SIZE>::weighted_bhattacharryya_partial(const Tuple<C,SIZE> &dev ,const Trianglix<C,SIZE>& other, double this_weight, double other_weight)const{
		double fout;
		// dev (P+P2/2)-1 dev
		C sum,tmp;
		C deter[2];
		int safe_mag =0;
		Trianglix<C, SIZE> invsc = ((*this) * (this_weight / (this_weight + other_weight))) + (other* (other_weight / (this_weight + other_weight))) ;
		Tuple<C,SIZE> xdev = dev;
		if (SIZE >1){
			C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
			C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
			double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
			unsigned int i,j,k;
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
				hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
				hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
				invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
				xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
			}



			deter[0] = invsc.data[0];
			ExOp::toOne(deter[1]);

			buf[SIZE] = ExOp::mkInverse(invsc.data[0]);
			for(j=1;j<SIZE;j++){
				k = (j*(j+3))/2;

				tmp = deter[((j & 1)^1)] * invsc.data[k];
				tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
				if (!(ExOp::isValid(tmp))){ // overflow!
					deter[0] *= pow(0.5f, 300.0f);
					deter[1] *= pow(0.5f, 300.0f);
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
					safe_mag++;
				}else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
					deter[0] *= pow(2.0f, 300.0f);
					deter[1] *= pow(2.0f, 300.0f);
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
					safe_mag--;
				}
				deter[(j & 1)]  = tmp;





				xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				if (j==1) sum = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[SIZE];
				else sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
				invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				buf[SIZE+j] = ExOp::mkInverse(invsc.data[k]);
				//     sum += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
			sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];

			fout = (this_weight + other_weight) * 0.5f *( ExOp::norm(sum) + log(ExOp::norm(deter[((SIZE & 1)^1)])) + log(pow(2.0f,300.0f)) * safe_mag);



			delete[](hv);
			delete[](buf);
			delete[](normbuf);
		}else fout = (this_weight + other_weight) * 0.5f *(ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));


		return(fout);
    }
	LFHTEMP double Trianglix<C, SIZE>::weighted_bhattacharryya_partial_MK2(const Tuple<C,SIZE> &dev ,const Trianglix<C,SIZE>& other, double this_weight, double other_weight)const{
		double fout;
		// dev (P+P2/2)-1 dev
		C sum,tmp;
		C deter[2];
		int safe_mag =0;
		Trianglix<C, SIZE> invsc = ((*this) * (this_weight / (this_weight + other_weight))) + (other* (other_weight / (this_weight + other_weight))) ;
		Tuple<C,SIZE> xdev = dev;
		if (SIZE >1){
			C* buf = new C[SIZE*2]; LFH_NICE_ALLOCERROR(buf,"")
			C* hv = new C[SIZE *(SIZE+1)/2]; LFH_NICE_ALLOCERROR(hv,"")
			double* normbuf = new double[SIZE-2]; LFH_NICE_ALLOCERROR(normbuf,"")
			unsigned int i,j,k;
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mkTrju(invsc.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);}
				hv[k+i] = ExOp::mkTrju(invsc.data[k+i]);
				hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
				invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
				xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
			}



			deter[0] = invsc.data[0];
			ExOp::toOne(deter[1]);

			buf[SIZE] = ExOp::mkInverse(invsc.data[0]);
			for(j=1;j<SIZE;j++){
				k = (j*(j+3))/2;

				tmp = deter[((j & 1)^1)] * invsc.data[k];
				tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
				if (!(ExOp::isValid(tmp))){ // overflow!
					deter[0] *= pow(0.5f, 300.0f);
					deter[1] *= pow(0.5f, 300.0f);
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
					safe_mag++;
				}else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
					deter[0] *= pow(2.0f, 300.0f);
					deter[1] *= pow(2.0f, 300.0f);
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mkTrju(invsc.data[k-1]) * invsc.data[k-1]);
					safe_mag--;
				}
				deter[(j & 1)]  = tmp;





				xdev[j] -= ExOp::mkTrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				if (j==1) sum = ExOp::mkTrju(xdev[0]) * xdev[0] * buf[SIZE];
				else sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
				invsc.data[k] -= ExOp::mkTrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				buf[SIZE+j] = ExOp::mkInverse(invsc.data[k]);
				//     sum += ExOp::mkTrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
			sum += ExOp::mkTrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];

			fout = (this_weight + other_weight) * 0.125f * ExOp::norm(sum) + 0.5f * log(ExOp::norm(deter[((SIZE & 1)^1)])) + 150.0f * log(2.0f) * safe_mag;



			delete[](hv);
			delete[](buf);
			delete[](normbuf);
		}else fout = (this_weight + other_weight) * 0.5f *(ExOp::mkTrju(dev[0]) * dev[0] * ExOp::mkInverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));


		return(fout);
    }

LFHTEMP template<class O> auto Trianglix<C, SIZE>::mkInvDivi(const Tuple<O,SIZE>& dev)
 -> Tuple<decltype(data[0]*dev[0]), SIZE>{
    Tuple<decltype(data[0]*dev[0]), SIZE> fout;

return fout;}
LFHTEMP template<class O> auto Trianglix<C, SIZE>::mkPosInvDivi(const Tuple<O,SIZE>& dev)
 -> Tuple<decltype(data[0]*dev[0]), SIZE>{
    Tuple<decltype(data[0]*dev[0]), SIZE> fout;

return fout;}


LFHTEMP Trianglix<C, SIZE> Trianglix<C, SIZE>::mkEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc) const{Trianglix<C, SIZE> fout = *this; fout.toEigenTransform(fnc);return fout;}
LFHTEMP template<class F> Trianglix<C, SIZE> Trianglix<C, SIZE>::mkEigenTransform(F fnc) const{Trianglix<C, SIZE> fout = *this; fout.toEigenTransform(fnc);return fout;}

#ifdef Rcpp_hpp

LFHTEMP Trianglix<C, SIZE>& Trianglix<C, SIZE>::toEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc){
    C tmp[4];
    if (SIZE < 3){
        if (SIZE == 1) ExOp::execMemfnc(data[0],fnc);
        else{
            tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            ExOp::execMemfnc(tmp[3],fnc);
            tmp[2] = tmp[2] - tmp[1] ;
            ExOp::execMemfnc(tmp[2],fnc);
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        arma::Mat<C> tmp; this->wrMatrix(tmp);
        arma::Col<C> eigval;
        arma::Mat<C> eigvec;
        arma::eig_sym( eigval, eigvec, tmp );
        for(int i =0;i< this->getSize();i++) ExOp::execMemfnc(eigval.at(i),fnc);
        return this->toMatrixRecomb(eigvec, eigval);
    }
return *this;}
LFHTEMP template<class F> Trianglix<C, SIZE>& Trianglix<C, SIZE>::toEigenTransform(F fnc){
    C tmp[4];
    if (SIZE < 3){
        if (SIZE == 1) fnc(data[0]);
        else{
            //tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            fnc(tmp[3]);
            tmp[2] = tmp[2] - tmp[1] ;
            fnc(tmp[2]);
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        arma::Mat<C> tmp; this->wrMatrix(tmp);
        arma::Col<C> eigval;
        arma::Mat<C> eigvec;
        arma::eig_sym( eigval, eigvec, tmp );
        for(int i =0;i< this->getSize();i++) fnc(eigval.at(i));
        return this->toMatrixRecomb(eigvec, eigval);
    }
return *this;}
LFHTEMP template<class RCL> void Trianglix<C,SIZE>::rdMatrix(const arma::Mat<RCL> &where){
    int i,j,k;
    if (where.nrow() < SIZE) myexit("input matrix is too small");
    for(i=0,k=0;i<SIZE;i++){
        for(j=0;j<=i;j++) data[k++] = where.at(i,j);
    }
}
LFHTEMP void Trianglix<C,SIZE>::rdMatrix(const Rcpp::NumericMatrix &where){
    int i,j,k;
    if (where.nrow() < SIZE) myexit("input matrix is too small");
    for(i=0,k=0;i<SIZE;i++){
        for(j=0;j<=i;j++) data[k++] = where(i,j);
    }
}
LFHTEMP template<class RCL> void Trianglix<C,SIZE>::wrMatrix(arma::Mat<RCL> &where)const{
    where = arma::Mat<RCL>(SIZE,SIZE);
    int i,j,k;
    for(i=0,k=0;i<SIZE;i++){
        for(j=0;j<i;j++) {
            where.at(i,j) = data[k];
            where.at(j,i) = data[k++];
        }
        where.at(i,i) = data[k++];
    }
}

LFHTEMP void Trianglix<C,SIZE>::wrMatrix(Rcpp::NumericMatrix &where) const{
    where = Rcpp::NumericMatrix(SIZE,SIZE);
    int i,j,k;
    for(i=0,k=0;i<SIZE;i++){
        for(j=0;j<i;j++) {
            where(i,j) = data[k];
            where(j,i) = data[k++];
        }
        where(i,i) = data[k++];
    }
}

LFHTEMP template<class RCL> Trianglix<C,SIZE>& Trianglix<C,SIZE>::toMatrixRecomb(const arma::Mat<RCL> &ortho, const arma::Col<RCL> &eigen){
    int i,j,k,l;
    for(i=0,k=0;i<SIZE;i++,k++){
        for(j=0;j<i;j++,k++) {
            data[k] = ortho.at(j,0) * ortho.at(i,0) * eigen.at(0);
            for(l=1;l< SIZE;l++) data[k] += ortho.at(j,l) * ortho.at(i,l) * eigen.at(l);
        }
        data[k] = ortho.at(i,0) * ortho.at(i,0) * eigen.at(0);
        for(l=1;l< SIZE;l++) data[k] += ortho.at(i,l) * ortho.at(i,l) * eigen.at(l);
    }
    return *this;
}

#else

LFHTEMP Trianglix<C, SIZE>& Trianglix<C, SIZE>::toEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc){
    C tmp[10];
    if (SIZE < 3){
        if (SIZE == 1) ExOp::execMemfnc(data[0],fnc);
        else{
            tmp[0] = (data[0] - data[2]) * 0.5;
            tmp[1] = sqrt(tmp[0] * tmp[0] + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            ExOp::execMemfnc(tmp[3],fnc);
            tmp[2] -= tmp[1] ;
            ExOp::execMemfnc(tmp[2],fnc);
            tmp[1] = (tmp[0] / tmp[1]) * 0.5;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] = tmp[3] * tmp[1];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{
        if (SIZE == 3){
            tmp[0] = data[0] * (data[2] * data[5] - data[4] * data[4]) + data[1] * (data[4] * data[3] - data[1] * data[5]) + data[3] * (data[1] * data[4] - data[2] * data[3]);
            tmp[1] = data[0] * (data[2] + data[5]) - data[1] * data[1] - data[3] * data[3] - data[4] * data[4] + data[2] * data[5];
            tmp[2] = data[0] + data[2] + data[5];
            double q = (tmp[1] - (tmp[2]*tmp[2] / 3.0) ) / 3.0;
            double r = (-1.5f*tmp[0] - ( ((tmp[2]*tmp[2]*tmp[2] / tmp) - 1.5f*tmp[2]*tmp[1] )/3.0) )/3.0;
            tmp = q*q*q + r*r;



        }




    }

return *this;}
LFHTEMP template<class F> Trianglix<C, SIZE>& Trianglix<C, SIZE>::toEigenTransform(F fnc){
    C tmp[4];
    if (SIZE < 3){
        if (SIZE == 1) fnc(data[0]);
        else{
            //tmp[0] = data[0] * data[2] -  data[1] * data[1];
            tmp[0] = data[0] - data[2];
            tmp[1] = sqrt(tmp[0] * tmp[0] * 0.25 + data[1] * data[1]);
            data[1] /= tmp[1];
            tmp[2] = (data[0] + data[2]) * 0.5;
            tmp[3] = tmp[2] + tmp[1];
            fnc(tmp[3]);
            tmp[2] = tmp[2] - tmp[1] ;
            fnc(tmp[2]);
            tmp[1] = (tmp[0] / tmp[1]) * 0.25;
            data[2] = (tmp[3] + tmp[2]) * 0.5;
            tmp[3] -= tmp[2];
            tmp[1] *= tmp[3];
            data[0] = data[2] + tmp[1];
            data[2] -= tmp[1];
            data[1] *= tmp[3] * 0.5;
        }
    }else{





    }

return *this;}
#endif // Rcpp_hpp


LFHTEMP template<class O, unsigned int OSIZE, Tuple_flag Cflag> auto Trianglix<C, SIZE>::operator*(const Tuple<O,OSIZE, Cflag> & a)const
 -> Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL>{
    Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL> fout;
    uint32_t i,k;
    uint32_t dmax = a.getSize();
    dmax = (dmax > SIZE) ? SIZE : dmax;
    for(k=0;k<dmax;k++){
        const C* currow = data + ((k * (k+1u)) >> 1);
        fout[k] = currow[k] * a[k];
        for(i=0;i<k;i++) {
            fout[k] += currow[i] * a[i];
            fout[i] += currow[i] * a[k];
        }
    }
return fout;}


LFHTEMP void Trianglix<C, SIZE>::show(FILE* f, int level)const{
    unsigned int i,k,j; // fprintf(f,"%u\t%u\n",totsize, sizeof(data) / sizeof(C) );
    for(i=0,j=0,k=2;i<totsize;i++) {ExOp::show(data[i],f,level+1); if (i==j) {fprintf(f,"\n"); j+= k++;} else fprintf(f,"\t");}
}
LFHTEMP string Trianglix<C, SIZE>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",SIZE); return string("Trianglix<") + ExOp::type_tostring(data[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}
LFHTEMP double Trianglix<C, SIZE>::WishartLogDensity(const Trianglix<C, SIZE>& var, unsigned int nbsamples) const{
	double fout;

	unsigned int i;
	fout = ((double)SIZE * nbsamples) * log(0.5f);
	fout -= ((double)nbsamples) * this->log_determinant();
	fout += ((double)(nbsamples - SIZE - 1)) * var.log_determinant();
	fout -= var.trace_of_division(*this);

	fout *=0.5f;
	for(i=0;i<SIZE;i++) fout -= lngamma( ((double)(nbsamples -i)) *0.5f);
	return fout;
}

template<int L> OptionStruct::OptionStruct(const string values[L]){for(int i=0;i <L ;i++) data.addEntry(values[i]);}
template<typename... VARTIC> OptionStruct::OptionStruct(string value, VARTIC... vartic) : OptionStruct(vartic...){data.addEntry(value);}

template<class S> uint32_t ThreadBase::callThatEvent(KeyElem<Event<S>* , typename Event<S>::TEMPLATE_ARG> ev){return (*ev.k)(ev.d);}

template<class S> ERRCODE ThreadBase::startEvent(Event<S>* ev, typename Event<S>::TEMPLATE_ARG  s){
    if (nbactive == nbthreads ) return 1;
    thrds[nbactive] = new std::thread(ThreadBase::callThatEvent<S>, KeyElem<Event<S>* , typename Event<S>::TEMPLATE_ARG>(ev, s) ); LFH_NICE_ALLOCERROR(thrds[nbactive] ,"")
    nbactive++;
    return 0;
}
template<class S> void ThreadBase::startEvent_ThenWait(Event<S>* ev, typename Event<S>::TEMPLATE_ARG s){(*ev)(s);this->waitForAllThreads();}


template<class F, typename... ARGS> void ThreadBase::submit(F & func, ARGS... args){
    {
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(std::bind(&F::operator(), std::ref(func), args...));
        #endif // STDBINDFILTER
    }
    condvar.notify_one();
}

template<class F, typename... ARGS> void ThreadBase::submit_ThenWait(F & func, ARGS... args){
   /* {
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(std::bind(&F::operator(), std::ref(func), args...));
        #endif // STDBINDFILTER
    }*/

    func(args...);
    {
        std::unique_lock<std::mutex> lock(endmut);
        endcondvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
    }
    if (async_progress_maintain != 0xFFFFFFFF){
        progb.finish();
        async_progress_maintain = 0xFFFFFFFF;
    }
    flushMsgs();
}
template<class F> void ThreadBase::submit_Array(F & func, uint32_t nbthr){
    {
        std::lock_guard<std::mutex> lock(mut);
        for(int i = 1;i<nbthr;i++){
            #ifdef STDBINDFILTER
            todolist.push_front(std::bind(&F::operator(), std::ref(func), i));
            #endif // STDBINDFILTER
        }
    }
    condvar.notify_all();
    if (nbthr) std::bind(&F::operator(), std::ref(func), 0)(); //
    std::unique_lock<std::mutex> lock(endmut);
    endcondvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
}

template<class F, class M, typename... ARGS> void ThreadBase::submitFunction(F & func, M memfunc, ARGS... args){
    {
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(std::bind(memfunc, std::ref(func), args...));
        #endif // STDBINDFILTER
    }
    condvar.notify_one();
}
template<class F, class M, typename... ARGS> void ThreadBase::submitFunction_ThenWait(F & func, M memfunc, ARGS... args){
    {
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(std::bind(memfunc, std::ref(func), args...));
        #endif // STDBINDFILTER

    }
    std::unique_lock<std::mutex> lock(endmut);
    endcondvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
    if (async_progress_maintain != 0xFFFFFFFF){
        progb.finish();
        async_progress_maintain = 0xFFFFFFFF;
    }
    flushMsgs();
}
template<class F, class M> void ThreadBase::submitFunction_Array(F & func, M memfunc, uint32_t nbthr){
    {
        std::lock_guard<std::mutex> lock(mut);
        for(int i = 1;i<nbthr;i++){
            #ifdef STDBINDFILTER
            todolist.push_front(std::bind(memfunc, std::ref(func), i));
            #endif // STDBINDFILTER
        }
    }
    condvar.notify_all();
    if (nbthr) std::bind(memfunc, std::ref(func), 0)(); //
    std::unique_lock<std::mutex> lock(endmut);
    endcondvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
}
template<class C> void ThreadBase::show(const C& what){
    char buffer[65536];
    char* ptr = buffer;
    ExOp::show(what,ptr,0);
    if (ptr[-1] == '\n') ptr[-1] = '\0';
    else *ptr = '\0';
    msgs.insert(string(buffer));
}


template<class C> void Anything::operator=(C* val){
    if ((anytype & 63) == 'p'){
        if (anytype & 64) {
            *(void**)ptr = (void*) val;
        }else{ // got to make a copy
            *(void**)ptr = (void*) new C(*val); LFH_NICE_ALLOCERROR(*(void**)ptr,"")
        }
    }
}



/*
LFHTEMP	template<class A>
const Quaternion<C>& Quaternion<C>::operator*=(Quaternion<A> const & other){
		C tmp[3];
		tmp[0] = (*this)[0];
		tmp[1] = (*this)[1];
		tmp[2] = (*this)[2];
		(*this)[0] = (*this)[0]*other[0] - (*this)[1]*other[1] -(*this)[2]*other[2] -(*this)[3]*other[3];
		(*this)[1] = (*this)[1]*other[0] + tmp[0]*other[1] -(*this)[3]*other[2] +(*this)[2]*other[3];
		(*this)[2] = (*this)[2]*other[0] + (*this)[3]*other[1] +tmp[0]*other[2] -tmp[1]*other[3];
		(*this)[3] = (*this)[3]*other[0] - tmp[2]*other[1] +tmp[1]*other[2] +tmp[0]*other[3];
		return(*this);
	}
	*/

 // end of namespace







