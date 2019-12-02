/*
 * DataGrid.hpp
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

#undef LFHTEMP
#define LFHTEMP template <class C, unsigned int nbdim>

	LFHTEMP DataGrid<C, nbdim>::DataGrid(void* owner): data(NULL) {}
	LFHTEMP DataGrid<C, nbdim>::DataGrid(unsigned int* d, void* owner): data(NULL){
		setSizes(d);
	}
	LFHTEMP DataGrid<C, nbdim>::DataGrid(const DataGrid<C,nbdim>& other){

		memcpy(dims, other.dims, sizeof(unsigned int) * nbdim);
		unsigned int ts = dims[0];
		for(unsigned int i=1;i<nbdim;i++) ts *= dims[i];
		data = new C[ts]; LFH_NICE_ALLOCERROR(data, "Could not allocate for DataGrid::DataGrid\n(size= %i)\n", (int)(ts*sizeof(C)))
        for(ts--;ts != 0xFFFFFFFF;ts--) data[ts] = other.data[ts];
	}


//	LFHTEMP	DataGrid<C, nbdim>::operator TypedFunctor<DataGrid<C,nbdim>, C, const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE& >()const{return TypedFunctor<DataGrid<C,nbdim>, C, const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE& >(*this);}
//    LFHTEMP DataGrid<C, nbdim>::operator Tuple<uint16_t,27u>() const{Tuple<uint16_t,27u> fout; fout.toRand();return fout;}
//    LFHTEMP DataGrid<C, nbdim>::operator std::function<C(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE& )> () const{return std::bind(&operator(), *this, );}

		LFHTEMP	template <class D> DataGrid<C, nbdim>::DataGrid(DataGrid<D, nbdim> const & other): data(NULL){
			//		printf("ata assign!\n");
			setSizes(other.getSizes());
			unsigned int i = totsize();
			for(i--;i!=0xFFFFFFFF;i--) data[i] = other.data[i];
		}

LFHTEMP	template<unsigned int fsize> DataGrid<C, nbdim>::DataGrid(const DataGrid<Tuple<C,fsize> ,nbdim-1> &o){
    unsigned int ts = fsize;
    dims[0] = ts;
    unsigned int i;
    for(i=1;i<nbdim;i++) ts *= dims[i]= o.dims[i-1];
    data = new C[ts]; LFH_NICE_ALLOCERROR(data, "Could not allocate for DataGrid\n(size= %i)\n",(int)(ts*sizeof(C)))
    typename  DataGrid< Tuple<C,fsize> ,nbdim-1>::KeyIterator ite = o.getKeyIterator();

    i=0;
    int j;
    if (ite.first()) do{
        for(j=0;j<fsize;j++) data[i++] = o(ite())[j];
    } while (ite.next());

}

LFHTEMP	template<unsigned int fsize> DataGrid<C, nbdim>::DataGrid(const DataGrid< C[fsize]  ,nbdim-1>&o){
    unsigned int ts = fsize;
    dims[0] = ts;
    unsigned int i;
    for(i=1;i<nbdim;i++) ts *= dims[i]= o.dims[i-1];
    data = new C[ts]; LFH_NICE_ALLOCERROR(data,"Could not allocate for DataGrid\n(size= %i)\n",(int)(ts*sizeof(C)))
    //typedef class DataGrid< C[fsize] ,nbdim-1>::KeyIterator itetype;
    //itetype ite = o.getKeyIterator();

    typename  DataGrid< C[fsize] ,nbdim-1>::KeyIterator ite = o.getKeyIterator();
    i=0;
    int j;
    if (ite.first()) do{
        for(j=0;j<fsize;j++) data[i++] = o(ite())[j];
    } while (ite.next());
}

LFHTEMP	template<unsigned int fsize> DataGrid<C, nbdim>::DataGrid(const DataGrid< DataGrid<C,fsize> , nbdim-fsize >&o){
    //todo
}
LFHTEMP DataGrid<C, nbdim>::~DataGrid(){if (data != NULL) delete[](data);}
LFHTEMP unsigned int DataGrid<C, nbdim>::totsize() const{ // total size of array
    unsigned int fout= dims[0];
    for(unsigned int i=1;i<nbdim;i++) fout *= dims[i];
return fout;}
LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::mkconcatenate(const DataGrid<C,nbdim> & other , unsigned int direction) const{ DataGrid<C,nbdim> fout;
    unsigned int i,j,k,chunksize, loopy;
    Tuple<unsigned int, nbdim> coor;
    loopy = chunksize =1;
    for(i=0;i<direction;i++) {coor[i] = dims[i]; chunksize *= dims[i];if (dims[i] != other.dims[i]) printf("invalid concatenate\n"); }
    coor[i] = dims[i] + other.dims[i];
    for(i++;i<nbdim;i++) {coor[i] = dims[i]; loopy *= dims[i]; if (dims[i] != other.dims[i]) printf("invalid concatenate\n");}
    fout.setSizes(coor);
    for(j=0,k=0;loopy != 0;loopy--){
        for(i=0; i < chunksize * dims[direction];i++,j++) fout.data[j+k] = this->data[j];
        for(i=0; i < chunksize * other.dims[direction];i++,k++) fout.data[j+k] = other.data[k];
    }
return fout;}
LFHTEMP DataGrid<C,nbdim-1> DataGrid<C, nbdim>::makeSlice(uint32_t which_dim, int coor)const{
    DataGrid<C,nbdim-1> fout;
    Tuple<unsigned int, nbdim-1> dasize;
    uint32_t i,j,k;
    unsigned int basechunk =1;
    unsigned int nbrep =1;
    for(i=0;i<which_dim;i++) {dasize[i] = dims[i]; basechunk *= dims[i];}
    for(;i<nbdim-1;i++) {dasize[i] = dims[i+1]; nbrep *= dims[i+1];}
    fout.setSizes(dasize);
    C* base = data + basechunk * coor;
    C* targ = fout.data;
    for(k=0; k< nbrep;k++){
        for(j=0; j< basechunk;j++) (*targ++) = (*base++);
        base += basechunk * (dims[which_dim] -1);
    }
return(fout);}
LFHTEMP DataGrid<C,nbdim+1> DataGrid<C, nbdim>::fromSlice(uint32_t which_dim) const{
    DataGrid<C,nbdim+1> fout;
    Tuple<unsigned int, nbdim+1> dasize;
    unsigned int i;
    unsigned int totsize =1;
    for(i=0;i<which_dim;i++) {dasize[i] = dims[i]; totsize*=dims[i];}
    for(;i<nbdim;i++) {dasize[i+1] = dims[i];totsize*=dims[i];}
    dasize[which_dim] =1;
    fout.setSizes(dasize);
    C* base = data;
    C* targ = fout.data;
    for(i=0; i< totsize;i++) (*targ++) = (*base++);
return fout;}

	LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::selectedSlices(unsigned int which_dim, Vector<unsigned int> coor){
        DataGrid<C,nbdim> fout;
        Tuple<unsigned int, nbdim> dasize;unsigned int i;
        for(i=0;i<which_dim;i++) {dasize[i] = dims[i]; totsize*=dims[i];}
        dasize[which_dim] = coor.getSize();
        for(i++;i<nbdim;i++) {dasize[i] = dims[i];totsize*=dims[i];}
        fout.setSizes(dasize);
		C* base = data;
		C* targ = fout.data;

        return fout;
	}


    LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toMemmove(DataGrid<C,nbdim>& source){
		delete[](data);
        data = source.data;
        source.data = NULL;
        memcpy(dims, source.dims, sizeof(unsigned int) * nbdim);
        return(*this);
    }

    LFHTEMP template<unsigned int fsize> DataGrid<Tuple<C,fsize>,nbdim-1> DataGrid<C, nbdim>::mkTupleGrid(const Tuple<unsigned int,fsize> &which) const{
    DataGrid<Tuple<C,fsize>,nbdim-1> fout;fout.setSizes(dims+1);

    Tuple<C,fsize>* cur= fout.data;
    C* icur = data;
    unsigned int i,ts,j;
    ts = fout.totsize();
    for(i=0;i<ts;i++,cur++, icur+= dims[0]){
        for(j=0;j<fsize;j++) (*cur)[j] = icur[which[j]];
        }

    return(fout);

    }


    LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::mkOrthoFlip(unsigned int w_s, unsigned int w_t) const{ DataGrid<C, nbdim> fout;
    	fout.setSizes(dims);
    	unsigned int coor[nbdim];
    	if (w_s >= w_t){
	    	if (w_s == w_t) {fout = (*this); return fout;}
	    	coor[0] = w_s; w_s = w_t; w_t = coor[0];
    	}
    	fout.dims[w_s] = dims[w_t];
    	fout.dims[w_t] = dims[w_s];
		unsigned int dir;
		ExOp::toZero(coor);
		C* cur = fout.data;
		while(true){
			*(cur++) = (*this)(coor);
			for(dir=0;dir<w_s;dir++) if (coor[dir] == dims[dir]-1) {coor[dir] =0;} else {coor[dir]++; break;}
			if (dir<w_s) continue;
			if (coor[w_t] == dims[w_t]-1) {coor[w_t] =0; dir++; } else {coor[w_t]++; continue;}
			for(;dir<w_t;dir++) if (coor[dir] == dims[dir]-1) {coor[dir] =0;} else {coor[dir]++; break;}
			if (dir<w_t) continue;
			if (coor[w_s] == dims[w_s]-1) {coor[w_s] =0; dir++; } else {coor[w_s]++; continue;}
			for(;dir<nbdim;dir++) if (coor[dir] == dims[dir]-1) {coor[dir] =0;} else {coor[dir]++; break;}
			if (dir == nbdim) break;
		}

	    return fout;
	    }
	LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::mkDimFlip(unsigned int which_dim) const{ DataGrid<C, nbdim> fout;
		fout.setSizes(dims);
		unsigned int coor[nbdim];
		unsigned int cooralt[nbdim];
		unsigned int altdir;
		unsigned int i;
		ExOp::toZero(coor); ExOp::toZero(cooralt);
		while(true){
			//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
			// prepare initial for backward!
			coor[which_dim] = 0;
			cooralt[which_dim] = getDimSize(which_dim)-1;
			while(coor[which_dim] < cooralt[which_dim]){
				fout(cooralt) = (*this)(coor);
				fout(coor) = (*this)(cooralt);
				coor[which_dim]++;cooralt[which_dim]--;
				}
				if (coor[which_dim] == cooralt[which_dim]) fout(coor) = (*this)(coor);
				// iterate througth other dimentions!
			for(altdir=1;altdir<nbdim;altdir++) if (coor[(which_dim + altdir) % nbdim] == getDimSize((which_dim + altdir) % nbdim)-1) {coor[(which_dim + altdir) % nbdim] =0; cooralt[(which_dim + altdir) % nbdim] =0;} else {coor[(which_dim + altdir) % nbdim]++;cooralt[(which_dim + altdir) % nbdim]++; break;}
			if (altdir == nbdim) break;
		}
		return(fout);
		}
	LFHTEMP const DataGrid<C,nbdim>& DataGrid<C, nbdim>::toDimFlip(unsigned int which_dim){
		unsigned int coor[nbdim];
		unsigned int cooralt[nbdim];
		unsigned int altdir;
		unsigned int i;
		C swap;
		ExOp::toZero(coor); ExOp::toZero(cooralt);
		while(true){
			//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
			// prepare initial for backward!
			coor[which_dim] = 0;
			cooralt[which_dim] = getDimSize(which_dim)-1;
			while(coor[which_dim] < cooralt[which_dim]){
				swap = (*this)(coor);
				(*this)(coor) = (*this)(cooralt);
				(*this)(cooralt) = swap;
				coor[which_dim]++;cooralt[which_dim]--;
				}
				// iterate througth other dimentions!
			for(altdir=1;altdir<nbdim;altdir++) if (coor[(which_dim + altdir) % nbdim] == getDimSize((which_dim + altdir) % nbdim)-1) {coor[(which_dim + altdir) % nbdim] =0; cooralt[(which_dim + altdir) % nbdim] =0;} else {coor[(which_dim + altdir) % nbdim]++;cooralt[(which_dim + altdir) % nbdim]++; break;}
			if (altdir == nbdim) break;
		}
		return (*this);
		}


	LFHTEMP
	template<class D> void DataGrid<C, nbdim>::matrixleftmultiply(const ConstGrid<D,2>& R){
		// R is squarre matrix!
		unsigned int coor[2];

		C* xformed= new C[dims[1]]; LFH_NICE_ALLOCERROR(xformed,"Could not allocate for DataGrid::matrixleftmultiply\n(size= %i)\n",(int)(dims[1]*sizeof(C)))
        unsigned int i,j;

		for(j=0;j<dims[0];j++){
			coor[0] =0;
			for(coor[1]=0;coor[1]<dims[1];coor[1]++) xformed[coor[1]] = this->data[j] * R(coor);
			for(coor[0]++;coor[0]<dims[1];coor[0]++){
				for(coor[1]=0;coor[1]<dims[1];coor[1]++) xformed[coor[1]] += this->data[j + coor[0] * dims[0]] * R(coor);
				}
			for(i=0;i<dims[1];i++) this->data[j + i * dims[0]] = xformed[i];
			}
		delete[](xformed);
		}

	LFHTEMP
	void DataGrid<C, nbdim>::settomatrixMultiply(const DataGrid<C,nbdim>& L, const DataGrid<C,nbdim>& R){
		// if dimentions does not fix, gets extended by identity matrix

		unsigned int coor[2];
		coor[0] = R.dims[0];
		coor[1] = L.dims[1];
		int min;
		bool w = (R.dims[0] < L.dims[1]);
		min = (w) ? R.dims[1]: L.dims[0];
		setSizes(coor);
		unsigned int k;

		for(coor[1] = 0;coor[1]<L.dims[1];coor[1]++){
			for(coor[0]=0;coor[0]<R.dims[0];coor[0]++){
			(*this).data[coor[0] + coor[1]*R.dims[0]] = L.data[coor[1]*L.dims[0]] * R.data[coor[0]];

			for(k=1;k<min;k++){
				(*this).data[coor[0] + coor[1]*R.dims[0]] += L.data[k + coor[1]*L.dims[0]] * R.data[coor[0]+ k*R.dims[0]];
			}
			if (w){
				if ((coor[0] >= k)&&(coor[0] < L.dims[0])) (*this).data[coor[0] + coor[1]*R.dims[0]] += L.data[coor[0]+ coor[1]*L.dims[0]];
			}else{
				if ((coor[1] >= k)&&(coor[1] < R.dims[1])) (*this).data[coor[0] + coor[1]*R.dims[0]] += R.data[coor[0]+ coor[1]*R.dims[0]];
			}

			}
		}


	}

		LFHTEMP
		void DataGrid<C, nbdim>::settomatrixMultiplyofTransposed(const DataGrid<C,nbdim>& L, const DataGrid<C,nbdim>& R){
			// if dimentions does not fix, gets extended by identity matrix

			unsigned int coor[2];
			coor[0] = R.dims[1];
			coor[1] = L.dims[1];
			int min;
			bool w = (R.dims[1] < L.dims[1]);
			min = (w) ? R.dims[0]: L.dims[0];
			setSizes(coor);
			unsigned int k;

			for(coor[1] = 0;coor[1]<L.dims[1];coor[1]++){
				for(coor[0]=0;coor[0]<R.dims[1];coor[0]++){
					(*this).data[coor[0] + coor[1]*R.dims[1]] = L.data[coor[1]*L.dims[0]] * R.data[coor[0]*R.dims[0]];

					for(k=1;k<min;k++){
						(*this).data[coor[0] + coor[1]*R.dims[1]] += L.data[k + coor[1]*L.dims[0]] * R.data[k + coor[0]*R.dims[0]];
					}
					if (w){
						if ((coor[0] >= k)&&(coor[0] < L.dims[0])) (*this).data[coor[0] + coor[1]*R.dims[1]] += L.data[coor[0]+ coor[1]*L.dims[0]];
					}else{
						if ((coor[1] >= k)&&(coor[1] < R.dims[0])) (*this).data[coor[0] + coor[1]*R.dims[1]] += R.data[coor[1]+ coor[0]*R.dims[0]];
					}

				}
			}


		}


	LFHTEMP
	void DataGrid<C, nbdim>::leftHouseHolderMultiply_lame(const double * const vec, const int lenght,bool hint){

		// assumes dims[0] == dims[1];
		int min = dims[0] < dims[1] ? dims[0] : dims[1];
		int i,j,k;
		C* tmp = new C[lenght]; LFH_NICE_ALLOCERROR(tmp,"")

		min -= lenght;

		for(i = (hint) ? min : 0;i< dims[0];i++){
			for(j=0;j< lenght;j++) {
				tmp[j] = this->data[i + (min + j) * dims[0]] * vec[j];
				this->data[i + (min + j) * dims[0]] *= (1 - vec[j] * vec[j]);
			}

			for(k = 0;k< lenght;k++){
				for(j = 0;j< lenght;j++) if (j != k) this->data[i + (min + j) * dims[0]] -= tmp[k] * vec[j];
			}

		}

		delete[](tmp);


	}
	LFHTEMP
	void DataGrid<C, nbdim>::rightHouseHolderMultiply_lame(const double * const vec, const int lenght,bool hint){

		// assumes dims[0] == dims[1];
		int min = dims[0] < dims[1] ? dims[0] : dims[1];
		int i,j,k;
		C* tmp = new C[lenght]; LFH_NICE_ALLOCERROR(tmp,"")
		min -= lenght;
		for(i = (hint) ? min : 0;i< dims[0];i++){
			for(j=0;j< lenght;j++) {
				tmp[j] = this->data[(min + j) + i * dims[0]] * vec[j];
				this->data[(min + j) + i * dims[0]] *= (1 - vec[j] * vec[j]);
			}

			for(k = 0;k< lenght;k++){
				for(j = 0;j< lenght;j++) if (j != k) this->data[(min + j) + i * dims[0]] -= tmp[k] * vec[j];
			}

		}

		delete[](tmp);


	}


		LFHTEMP
		void DataGrid<C, nbdim>::leftHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){

			int inl = (hint) ? lenght : dims[0];

			C* min = data + dims[0] - inl +(dims[1] - lenght) * dims[0];

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner ,"")

			for(j=0;j< inl;j++) inner[j] = min[ j + i* dims[0]] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += min[ j + i * dims[0]] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) min[j + i * dims[0]] += inner[j] * (vec[i] * fact);
			}

			delete[](inner);
		}

		// hint implies that the outside the region of interes are zeroes only!
		LFHTEMP
		void DataGrid<C, nbdim>::rightHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){

			int inl = (hint) ? lenght : dims[1];
			int minj = dims[1] - inl;

			C* min = data + dims[0] - lenght + minj * dims[0];

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner ,"")

			for(j=0;j< inl;j++) inner[j] = min[ i + j * dims[0]] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += min[ i +  j * dims[0]] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) min[i + j * dims[0]] += inner[j] * (vec[i] * fact);
			}

			delete[](inner);

		}

		LFHTEMP  template<class D>
		void DataGrid<C, nbdim>::leftHouseHolderMultiply(const D * const vec, const int lenght, const D& sqrt_den, bool hint){

			int inl = (hint) ? lenght : dims[0];

			C* min = data + dims[0] - inl +(dims[1] - lenght) * dims[0];

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner ,"")

			for(j=0;j< inl;j++) inner[j] = min[ j + i* dims[0]] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += min[ j + i * dims[0]] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) min[j + i * dims[0]] += inner[j] * (vec[i] * fact);
			}

			delete[](inner);
		}

		// hint implies that the outside the region of interes are zeroes only!
		LFHTEMP template<class D>
		void DataGrid<C, nbdim>::rightHouseHolderMultiply(const D * const vec, const int lenght, const D& sqrt_den, bool hint){

			int inl = (hint) ? lenght : dims[1];
			int minj = dims[1] - inl;

			C* min = data + dims[0] - lenght + minj * dims[0];

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner ,"")

			for(j=0;j< inl;j++) inner[j] = min[ i + j * dims[0]] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += min[ i +  j * dims[0]] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) min[i + j * dims[0]] += inner[j] * (vec[i] * fact);
			}

			delete[](inner);

		}


		// multiply on both sides, assume symetric squarre matrix!
		LFHTEMP
		void DataGrid<C, nbdim>::dualHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){

			if (true) {
				if (dims[1] != dims[0]) static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT;
			}

			int inl = (hint) ? lenght : dims[0];

			C* min = data + dims[0] - inl +(dims[1] - lenght) * dims[0];

			double fact = -2.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl]; LFH_NICE_ALLOCERROR(inner ,"")
			j = inl -1;
			i = lenght-1;
			inner[j] = min[j + i* dims[0]] * vec[i];
			C squarre =  inner[j]* vec[i];
			C tmp;
			int k = inl -1;
			for(j--;j>=0;j--) {inner[j] = min[ j + i* dims[0]] * vec[i]; squarre += inner[j] * vec[i]; inner[k] += inner[j];}

			for(i--,k--;i>=0;i--,k--){
				j=k;
				inner[j] += tmp = min[ j + i* dims[0]] * vec[i];
				squarre += tmp * vec[i];
				for(j--;j>=dims[0]-inl;j--) {
					inner[j] += tmp = min[ j + i* dims[0]] * vec[i];
					squarre += tmp * vec[i];
					inner[k] += tmp;
				}
				for(;j>=0;j--) {
					inner[j] += tmp = min[ j + i* dims[0]] * vec[i];
					squarre += tmp * vec[i];
				}
			}

			squarre *= pow(sqrt_den, -2.0f);

			for(i=0;i<lenght;i++){
				// diagonal!
				min[dims[0] - inl + i * (dims[0]+1)] += inner[dims[0] - inl + i] * (vec[i] * fact) + squarre * vec[i] * vec[i];
				for(j=0;j< inl;j++) min[j + i * dims[0]] += inner[j] * (vec[i] * fact) + squarre * vec[i] * vec[j];
			}

			delete[](inner);
			// M' = M - 2 * M * I + V * S * V
		}


	LFHTEMP
	void DataGrid<C, nbdim>::leftHouseHolderMultiply_lame(const double * const vec, const int lenght, const double& sqrt_den,bool hint){

		// assumes dims[0] == dims[1];
		int min = dims[0] < dims[1] ? dims[0] : dims[1];
		int i,j,k;
		C* tmp = new C[lenght]; LFH_NICE_ALLOCERROR(tmp ,"")
		min -= lenght;
		double fact = 1.0f / sqrt_den;
		if (!ExCo<double>::isValid(fact)) return; // either NAN, or den =0, dont do a thing in case it is den = 0
		for(i = (hint) ? min : 0;i< dims[0];i++){
			for(j=0;j< lenght;j++) {
				tmp[j] = this->data[i + (min + j) * dims[0]] * vec[j]* fact;
				this->data[i + (min + j) * dims[0]] *= (1 - vec[j] * vec[j]* fact);
			}

			for(k = 0;k< lenght;k++){
				for(j = 0;j< lenght;j++) if (j != k) this->data[i + (min + j) * dims[0]] -= tmp[k] * vec[j];
			}

		}

		delete[](tmp);


	}
	LFHTEMP
	void DataGrid<C, nbdim>::rightHouseHolderMultiply_lame(const double * const vec, const int lenght, const double& sqrt_den,bool hint){

		// assumes dims[0] == dims[1];
		int min = dims[0] < dims[1] ? dims[0] : dims[1];
		int i,j,k;
		C* tmp = new C[lenght]; LFH_NICE_ALLOCERROR(tmp ,"")
		min -= lenght;
		double fact = 1.0f / sqrt_den;
		if (!ExCo<double>::isValid(fact)) return; // either NAN, or den =0, dont do a thing in case it is den = 0
		for(i = (hint) ? min : 0;i< dims[0];i++){
			for(j=0;j< lenght;j++) {
				tmp[j] = this->data[(min + j) + i * dims[0]] * vec[j] * fact;
				this->data[(min + j) + i * dims[0]] *= (1 - vec[j] * vec[j]* fact);
			}

			for(k = 0;k< lenght;k++){
				for(j = 0;j< lenght;j++) if (j != k) this->data[(min + j) + i * dims[0]] -= tmp[k] * vec[j];
			}

		}

		delete[](tmp);


	}


	LFHTEMP
	void DataGrid<C, nbdim>::rightoffdiagelimination(double factor, int from, int to){
		unsigned int i;
		for(i = 0;i< dims[1];i++){
			this->data[to + i * dims[0]] += this->data[from + i * dims[0]] * factor;
		}
	}

	LFHTEMP
	void DataGrid<C, nbdim>::leftoffdiagelimination(double factor, int from, int to){
		unsigned int i;
		for(i = 0;i< dims[0];i++){
			this->data[i + to * dims[0]] += this->data[i + from * dims[0]] * factor;
		}
	}


	LFHTEMP
		void DataGrid<C, nbdim>::LinearTransform(const C* const target, C * _out){
			LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
			int i,k;
			for(k=0;k<dims[1];k++){
				_out[k] = target[0] * this->data[k* dims[0]];
				for(i=1;i<dims[0];i++){
					_out[k] += target[i] * this->data[i + k * dims[0]];
				}
			}
		}



	LFHTEMP
	void DataGrid<C, nbdim>::solveLinearSystem(const C* const target, C * _out){
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!

		DataGrid<double, nbdim> R;
		DataGrid<double, nbdim> L;


		int i,j,k;



		Tuple<unsigned int,2> coor;
		coor[0] = dims[1];
		//	coor[1] = (dims[0] < dims[1]) ? dims[0] : dims[1];
		coor[1] = dims[1];
		L.setSizes(coor);
		coor[0] = dims[0];
		//	coor[1] = (dims[0] < dims[1]) ? dims[0] : dims[1];
		coor[1] = dims[0];
		R.setSizes(coor);


		{//:
		KeyIterator ite = KeyIterator(R);
		if (ite.first()) do{
			coor = ite();
			R(coor) = (coor[0] == coor[1]) ? 1.0f : 0.0f;
		} while (ite.next());
		}//:

		{//:
			KeyIterator ite = KeyIterator(L);
			if (ite.first()) do{
				coor = ite();
				L(coor) = (coor[0] == coor[1]) ? 1.0f : 0.0f;
			} while (ite.next());
		}//:


		double sum2,tmp;



		double* hh = new double[dims[0]]; LFH_NICE_ALLOCERROR(hh ,"")
		for(i=0,j=0; i < dims[0]-1; j++){

		//for(i=0;i<dims[0];i++){
			sum2 = this->data[i*(dims[0]+1)]*this->data[i*(dims[0]+1)];
			for(k=j+1;k<dims[1];k++) {
				hh[k-j] = this->data[i + k * dims[0]];
				sum2 += hh[k-j] * hh[k-j];
			}



		//}
		tmp = sqrt(sum2);
		hh[0] = (this->data[i + j * dims[0]] + ((this->data[i + j * dims[0]] < 0) ? -tmp: tmp));
		L.leftHouseHolderMultiply(hh,dims[0]-j,(sum2 + fabs(this->data[i + j * dims[0]] *tmp)),false);
		this->leftHouseHolderMultiply(hh,dims[0]-j,(sum2 + fabs(this->data[i + j * dims[0]] *tmp)),false);
//		for(i=0;i<dims[1];i++)

		//	this->show(stdout);printf("\n");
		//	L.show(stdout);printf("\n");

	//	for(j=i+1;j<dims[0];j++) sum2 += this->data[i + j * dims[0]] * this->data[i + j * dims[0]];
			i++;

			sum2 = this->data[1+j*(dims[0]+1)]*this->data[1+j*(dims[0]+1)];
			for(k=i+1;k<dims[0];k++) {
				hh[k-i] = this->data[k + j * dims[0]];
				sum2 += hh[k-i] * hh[k-i];
			}
			tmp = sqrt(sum2);
			hh[0] = (this->data[i + j * dims[0]] + ((this->data[i + j * dims[0]] < 0) ? -tmp: tmp));

			R.rightHouseHolderMultiply(hh,dims[0]-i,(sum2 + fabs(this->data[i + j * dims[0]] *tmp)),false);
			this->rightHouseHolderMultiply(hh,dims[0]-i,(sum2 + fabs(this->data[i + j * dims[0]] *tmp)),false);

		//
		//	R.show(stdout);printf("\n");
		}
//		this->show(stdout);printf("\n");
		C* Xform_targ = new C[dims[1]];

		for(i=0;i<dims[1];i++){
			Xform_targ[i] = target[0] * L.data[i * dims[1]];
			for(j=1;j<dims[1];j++) Xform_targ[i] += target[j] * L.data[j + i * dims[1]];
//			printf("%f\n",Xform_targ[i]);
		}

		// Solving bidiagonal linear system;
		i = dims[0]-1;

		sum2 = fabs( this->data[dims[0] * (1+dims[0])] );
		for(i=0;i<dims[0]-1;i++){
			if (sum2 < fabs( this->data[i * (1+dims[0])] )) sum2 = fabs( this->data[i * (1+dims[0])] );
			if (sum2 < fabs( this->data[1+i * (1+dims[0])] )) sum2 = fabs( this->data[1+i * (1+dims[0])] );
		}


		while(fabs(this->data[i * (1+dims[0])] / sum2) < 1.0f / 65536.0f){
			Xform_targ[i] = ExCo<C>::zero();i--;
		}


		if ((i == dims[0]-1)||(fabs(this->data[1+i * (1+dims[0])] / sum2) < 1.0f / 65536.0f)){ // fully determined system!
			Xform_targ[i] = Xform_targ[i] / this->data[i * (1+dims[0])];
		for(i--;i>=0;i--){
			Xform_targ[i] = (Xform_targ[i] - ( this->data[1+i * (1+dims[0])] * Xform_targ[i+1] )) / this->data[i * (1+dims[0])];
		}
		}else{ // since bi-diagonal, 1 freedom is left
			j=i;
			C tmp = Xform_targ[i] / this->data[i * (1+dims[0])];
			C den = this->data[i * (1+dims[0])];
			C num = this->data[1+ i * (1+dims[0])];
			C sum0 = tmp;
			C sumK = (den/num) + (num / den);
			for(i--;i>=0;i--){
				if (fabs(this->data[1+i * (1+dims[0])] / sum2) < 1.0f / 65536.0f) break; // hit a off diagobal zero, the remaining is determined
				tmp = (Xform_targ[i] - ( this->data[1+i * (1+dims[0])] * tmp )) / this->data[i * (1+dims[0])];
				den *= this->data[i * (1+dims[0])];
				num *= this->data[1+ i * (1+dims[0])];
				sum0 = tmp - (this->data[i * (1+dims[0])] * sum0 / this->data[1+ i * (1+dims[0])]);
				sumK = (num  / den) + (this->data[i * (1+dims[0])] * sumK / this->data[1+ i * (1+dims[0])]);
			}
			Xform_targ[j+1] = -sum0 / sumK;
			for(i=j;i>=0;i--){
				Xform_targ[i] = (Xform_targ[i] - ( this->data[1+i * (1+dims[0])] * Xform_targ[i+1] )) / this->data[i * (1+dims[0])];
			}
		}
		sum2 =Xform_targ[5] * Xform_targ[5] ;
		for(i=0;i<5;i++){
//			printf("? %f\n",Xform_targ[i] *this->data[i * (1+dims[0])]   + Xform_targ[i+1] * this->data[1+i * (1+dims[0])]);
			sum2 +=  Xform_targ[i] * Xform_targ[i];
		}
//		printf("? %f\n",Xform_targ[i] *this->data[i * (1+dims[0])]);
//		printf("rerot norm %e\n", sqrt(sum2));
		for(i=0;i<dims[0];i++){
			_out[i] = Xform_targ[0] * R.data[i * dims[0]];
			for(j=1;j<dims[0];j++) _out[i] += Xform_targ[j] * R.data[j + i * dims[0]];
		}


		delete[](hh);
//		delete[](tmp_targ[0]);
//		delete[](tmp_targ[1]);




	}

	LFHTEMP
	void DataGrid<C, nbdim>::makeDiagonalizer(const ConstGrid<double,2> &symmetric_matrix, Vector<double> &eigenval){
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
		Tuple<unsigned int, 2> dims;
		symmetric_matrix.getDims(dims);
		setSizes(dims);
		eigenval.setSize(dims[0]);
		initeye();

		double* hh = new double[dims[0]]; LFH_NICE_ALLOCERROR(hh ,"")
		int i;

		unsigned int coor[2];
		double sum2;
		double tmp;

		DataGrid<double, 2> damat;

		damat.setSizes(dims);
		for(coor[0]=0;coor[0]<dims[0];coor[0]++) {
			for(coor[1]=0;coor[1]<dims[0];coor[1]++) {
				damat(coor) = symmetric_matrix(coor);
			}
		}







for(coor[0]=0;coor[0]+2<dims[0];coor[0]++){
			sum2 =0.0f;
			for(coor[1]=coor[0]+2;coor[1]<dims[0];coor[1]++) {
				hh[coor[1]] = damat(coor);
				sum2 += hh[coor[1]] * hh[coor[1]];
			}

			coor[1]=coor[0]+1;
			if ((sum2 != 0.0f)||(fabs(damat(coor) / sum2) > 1.0f / ExCo<double>::epsilon())) {
			sum2 += damat(coor)*damat(coor);
			tmp = sqrt(sum2);
			hh[coor[1]] = (damat(coor)+ ((damat(coor) < 0) ? -tmp: tmp));
			tmp = (sum2 + fabs(damat(coor) *tmp));
			rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
			damat.leftHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
			damat.rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
			}
		}



		for(coor[1]=0,coor[0]=0;;) {
			eigenval[coor[0]] = damat(coor);
			coor[0]++;
			if (coor[0] == dims[0]) break;
			tmp = damat(coor);

			sum2 = tmp / eigenval[coor[1]];

			rightoffdiagelimination( sum2,coor[0],coor[1]);
			coor[1]++;
			 damat(coor) -= tmp * sum2;
		}
		delete[](hh);
	}


		LFHTEMP
		void DataGrid<C, nbdim>::makeDiagonalizer_ofinverse(const ConstGrid<double,2> &symmetric_matrix, Vector<double> &eigenval){
			LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
			Tuple<unsigned int, 2> dims;
			symmetric_matrix.getDims(dims);
			setSizes(dims);
			eigenval.setSize(dims[0]);
			initeye();

            printf("hallo!\n");

			double* hh = new double[dims[0]]; LFH_NICE_ALLOCERROR(hh ,"")
			int i;

			unsigned int coor[2];
			double sum2;
			double tmp;

			DataGrid<double, 2> damat;

			damat.setSizes(dims);
			for(coor[0]=0;coor[0]<dims[0];coor[0]++) {
				for(coor[1]=0;coor[1]<dims[0];coor[1]++) {
					damat(coor) = symmetric_matrix(coor);
				}
			}







			for(coor[0]=0;coor[0]+2<dims[0];coor[0]++){

				sum2 =0.0f;
				for(coor[1]=coor[0]+2;coor[1]<dims[0];coor[1]++) {
					hh[coor[1]] = damat(coor);
					sum2 += hh[coor[1]] * hh[coor[1]];
				}

				coor[1]=coor[0]+1;
				if ((sum2 != 0.0f)||(fabs(damat(coor) / sum2) > 1.0f / ExCo<double>::epsilon())) {
				sum2 += damat(coor)*damat(coor);
				tmp = sqrt(sum2);
				hh[coor[1]] = (damat(coor)+ ((damat(coor) < 0) ? -tmp: tmp));
				tmp = (sum2 + fabs(damat(coor) *tmp));
				rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				damat.leftHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				damat.rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				}
			}



			for(coor[1]=0,coor[0]=0;;) {
				eigenval[coor[0]] = 1.0f / damat(coor);
				coor[0]++;
				if (coor[0] == dims[0]) break;
				tmp = damat(coor);

				sum2 = tmp * eigenval[coor[1]];

				rightoffdiagelimination(-sum2,coor[1],coor[0]);
				coor[1]++;
				damat(coor) -= tmp * sum2;
			}
			delete[](hh);
		}

		LFHTEMP
		template<unsigned int size> void DataGrid<C, nbdim>::makeDiagonalizer(const TMatrix<double,size,size> &symmetric_matrix, Vector<double> &eigenval){
			LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
			Tuple<unsigned int, 2> dims;
			dims[0] = size;
			dims[1] = size;
			setSizes(dims);
			eigenval.setSize(dims[0]);
			initeye();

			double* hh = new double[dims[0]]; LFH_NICE_ALLOCERROR(hh ,"")
			int i;

			unsigned int coor[2];
			double sum2;
			double tmp;

			DataGrid<double, 2> damat;

			damat.setSizes(dims);
			for(coor[0]=0;coor[0]<dims[0];coor[0]++) {
				for(coor[1]=0;coor[1]<dims[0];coor[1]++) {
					damat(coor) = symmetric_matrix(coor);
				}
			}







			for(coor[0]=0;coor[0]+2<dims[0];coor[0]++){


				sum2 =0.0f;
				for(coor[1]=coor[0]+2;coor[1]<dims[0];coor[1]++) {
					hh[coor[1]] = damat(coor);
					sum2 += hh[coor[1]] * hh[coor[1]];
				}

				coor[1]=coor[0]+1;
				if ((sum2 != 0.0f)||(fabs(damat(coor) / sum2) > 1.0f / ExCo<double>::epsilon())) {
				sum2 += damat(coor)*damat(coor);
				tmp = sqrt(sum2);
				hh[coor[1]] = (damat(coor)+ ((damat(coor) < 0) ? -tmp: tmp));
				tmp = (sum2 + fabs(damat(coor) *tmp));
				rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				damat.leftHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				damat.rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				}
			}


			for(coor[1]=0,coor[0]=0;;) {
				eigenval[coor[0]] = damat(coor);
				coor[0]++;
				if (coor[0] == dims[0]) break;
				tmp = damat(coor);

				sum2 = tmp / eigenval[coor[1]];

				rightoffdiagelimination( sum2,coor[0],coor[1]);
				coor[1]++;
				damat(coor) -= tmp * sum2;
			}
			delete[](hh);
		}


		LFHTEMP
		template<unsigned int size> void DataGrid<C, nbdim>::makeDiagonalizer_ofinverse(const TMatrix<double,size,size> &symmetric_matrix, Vector<double> &eigenval){
			LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
			Tuple<unsigned int, 2> tdims;
			tdims[0] = size;
			tdims[1] = size;
			setSizes(tdims);
			eigenval.setSize(tdims[0]);
			initeye();

			double* hh = new double[tdims[0]]; LFH_NICE_ALLOCERROR(hh ,"")

			unsigned int coor[2];
			double sum2;
			double tmp;

			DataGrid<double, 2> damat;

			damat.setSizes(tdims);
			for(coor[0]=0;coor[0]<tdims[0];coor[0]++) {
				for(coor[1]=0;coor[1]<tdims[1];coor[1]++) {
					damat(coor) = symmetric_matrix.data[coor[0] + coor[1] * tdims[0]];
				}
			}






			for(coor[0]=0;coor[0]+2<tdims[0];coor[0]++){

				sum2 =0.0f;
				for(coor[1]=coor[0]+2;coor[1]<tdims[0];coor[1]++) {
					hh[coor[1]] = damat(coor);
					sum2 += hh[coor[1]] * hh[coor[1]];
				}
				coor[1]=coor[0]+1;
				if ((sum2 != 0.0f)||(fabs(damat(coor) / sum2) > 1.0f / ExCo<double>::epsilon())) {
				sum2 += damat(coor)*damat(coor);
				tmp = sqrt(sum2);
				hh[coor[1]] = (damat(coor)+ ((damat(coor) < 0) ? -tmp: tmp));
				tmp = (sum2 + fabs(damat(coor) *tmp));
				rightHouseHolderMultiply(hh + coor[1],tdims[0]-coor[1],tmp,false);
				damat.leftHouseHolderMultiply(hh + coor[1],tdims[0]-coor[1],tmp,false);
				damat.rightHouseHolderMultiply(hh + coor[1],tdims[0]-coor[1],tmp,false);
				}
			}
// printf("trigiago:\n");			damat.show();            show();	          printf("then:\n");

			for(coor[1]=0,coor[0]=0;;) {
				eigenval[coor[0]] = 1.0f / damat(coor);
				coor[0]++;
				if (coor[0] == tdims[0]) break;
				tmp = damat(coor);

				sum2 = tmp * eigenval[coor[1]];

				rightoffdiagelimination(-sum2,coor[1],coor[0]);
				coor[1]++;
				damat(coor) -= tmp * sum2;
			}
      //      damat.show();            show();
			delete[](hh);
		}

		LFHTEMP
		void DataGrid<C, nbdim>::makeInverse(const ConstGrid<double,2> &symmetric_matrix){
			LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
			Tuple<unsigned int, 2> dims;
			symmetric_matrix.getDims(dims);
			setSizes(dims);
			unsigned int coor[2];
			double sum2;
			double tmp;
			switch(dims[0]){
				case 1:
					coor[0] = coor[1] = 0;
					this->data[0] = 1.0f / symmetric_matrix(coor);
				break;
				default:


			double* hh = new double[dims[0]]; LFH_NICE_ALLOCERROR(hh ,"")
			unsigned int i;


			DataGrid<double, 2> damat;


			damat.setSizes(dims);
			damat.initeye();
			for(coor[0]=0;coor[0]<dims[0];coor[0]++) {
				for(coor[1]=0;coor[1]<dims[0];coor[1]++) {
					(*this)(coor) = symmetric_matrix(coor);
				}
			}





			for(i=0;i<dims[0]-2;i++){

				coor[0] =i;
				sum2 =0.0f;
				for(coor[1]=coor[0]+2;coor[1]<dims[0];coor[1]++) {
					hh[coor[1]] = (*this)(coor);
					sum2 += hh[coor[1]] * hh[coor[1]];
				}

				coor[1]=coor[0]+1;
				sum2 += (*this)(coor)*(*this)(coor);
				tmp = sqrt(sum2);
				hh[coor[1]] = ((*this)(coor)+ (((*this)(coor) < 0) ? -tmp: tmp));
				tmp = (sum2 + fabs((*this)(coor) *tmp));
				damat.rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				leftHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
				rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
			}



			for(coor[1]=0,coor[0]=0;;) {
				hh[coor[0]] = 1.0f / (*this)(coor);
				coor[0]++;
				if (coor[0] == dims[0]) break;
				tmp = (*this)(coor);

				sum2 = tmp * hh[coor[1]];

				damat.rightoffdiagelimination(-sum2,coor[1],coor[0]);
				coor[1]++;
				(*this)(coor) -= tmp * sum2;
			}

			// This * Diag * This transposed!


			for(coor[1] =0;coor[1] <  dims[1]; coor[1]++){
				for(coor[0] =0;coor[0] < coor[1]; coor[0]++){
					(*this)(coor) = this->data[coor[1] + coor[0] *dims[0]];
				}

				for(;coor[0] <  dims[0]; coor[0]++){
					(*this)(coor) = damat.data[ coor[0] * dims[0]]* hh[0] * damat.data[coor[1] * dims[0]];
					for(i=1;i<dims[0];i++) (*this)(coor) += damat.data[i + coor[0] * dims[0]]* hh[i] * damat.data[i + coor[1] * dims[0]];
				}
			}


			delete[](hh);
			}
		}
LFHTEMP	void DataGrid<C, nbdim>::solve_sym_and_project(DataGrid<double, 2> &project, const ConstGrid<double,2> &symmetric_Matrix){

	}

	LFHTEMP
	void DataGrid<C, nbdim>::LinearSystem_fixpoint_solve(const DataGrid<C, nbdim>&target, double (*weights)(int a, int b)){




		}

	LFHTEMP
	void DataGrid<C, nbdim>::LinearSystem_residual(const DataGrid<C, nbdim>&guess, const DataGrid<C, nbdim>&target, const ConstGrid<double,2> &trans){
		int s_s = guess.dims[1];
		int s_d = guess.dims[0];
		Tuple<unsigned int, 2> coorA;
		Tuple<unsigned int, 2> coorB;
		double dist;
		int tmp;
		for(coorA[1]=0;coorA[1] < s_s; coorA[1]++){
			coorA[0] = coorA[1];
			dist = trans(coorA);
			for(coorA[0]=0; coorA[0] < s_d; coorA[0]++){
				(*this)(coorA) = target(coorA) - (guess(coorA) * dist);
			}
		}





			for(coorA[1]=0;coorA[1] < s_s; coorA[1]++)
				for(coorB[1]=0;coorB[1] < s_s; coorB[1]++){
					if (coorB[1] == coorA[1]) continue;
					coorB[0] = coorA[1];
					dist = trans(coorB);
					for(coorA[0]=0; coorA[0] < s_d; coorA[0]++){
						coorB[0] = coorA[0];
						(*this)(coorA) -= guess(coorB) * dist;
					}
			}

		}



	LFHTEMP
	void DataGrid<C, nbdim>::positivedefinite_inverse_and_involution(DataGrid<C, nbdim> &f_out, const DataGrid<C, nbdim> &invol){
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
		if ((dims[0] != dims[1])&&(dims[1] != invol.dims[0])) exit(1234);

		unsigned int coor[2];
		coor[0] = coor[1] = invol.dims[1];
		f_out.setSizes(coor);


	}


	LFHTEMP
	void DataGrid<C, nbdim>::setRow(const C* const _to, const Tuple<unsigned int, nbdim>& coor, const int direction){

		Tuple<unsigned int, nbdim> m_coor = coor;
		m_coor[direction] = dims[direction] -1;
		(*this)(m_coor) = _to[m_coor[direction]];
		for(m_coor[direction]--;m_coor[direction] != 0xFFFFFFFF;m_coor[direction]--) (*this)(m_coor) = _to[m_coor[direction]];
	}

LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator=(const DataGrid<C,nbdim>& other){
//	printf("datagrid assignement\n");
	setSizes(other.getSizes());
	unsigned int i = totsize();
	for(i--;i!=0xFFFFFFFF;i--) data[i] = other.data[i];
	return(*this);
}

LFHTEMP template<unsigned int S, Tuple_flag flag> DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator=(const DataGrid<Tuple<C, S, flag>,nbdim-1>& other){
	Tuple<unsigned int, nbdim> coor;
	unsigned int i;
	for(i=1;i<nbdim;i++) coor[i] = other.dims[i-1];
	dims[0] = S;
	this->setSizes(coor);
	C* d = this->data;
	C* o = other.data;
	i = totsize();
	for(i--;i!=0xFFFFFFFF;i--) *(d++) = *(o++);
	return *this;
}
LFHTEMP template<class O, unsigned int S, Tuple_flag flag> DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator=(const DataGrid<Tuple<O, S, flag>,nbdim-1>& other){
	Tuple<unsigned int, nbdim> coor;
	unsigned int i;
	for(i=1;i<nbdim;i++) coor[i] = other.dims[i-1];
	dims[0] = S;
	this->setSizes(coor);
	C* d = this->data;
	C* o = other.data;
	i = totsize();
	for(i--;i!=0xFFFFFFFF;i--) *(d++) = (C)*(o++);
	return *this;
}/*
LFHTEMP template<Tuple_flag flag> DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator=(const DataGrid<Tuple<C, 0u, flag>,nbdim-1>& other){
	LFH_FUNCTION_MISSING;
	return *this;
}
LFHTEMP template<class O, Tuple_flag flag> DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator=(const DataGrid<Tuple<O, 0u, flag>,nbdim-1>& other){
	LFH_FUNCTION_MISSING;
	return *this;
}*/

LFHTEMP	template<class O>
DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator=(const DataGrid<O,nbdim>& other){
//	printf("datagrid assignement\n");
	setSizes(other.getSizes());
	unsigned int i = totsize();
	for(i--;i!=0xFFFFFFFF;i--) data[i] = (C) other.data[i];
	return(*this);
}

LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator+=(const DataGrid<D, nbdim> & v){if (data == NULL) {this->setSizes(v.dims);this->toZero();} for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] += v.data[i]; return(*this);}

LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator-=(const DataGrid<D, nbdim> & v){if (data == NULL) {this->setSizes(v.dims);this->toZero();} for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] -= v.data[i]; return(*this);}
LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator*=(const DataGrid<D, nbdim> & v){if (data == NULL) {this->setSizes(v.dims);this->toZero();} for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] *= v.data[i]; return(*this);}
LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator/=(const DataGrid<D, nbdim> & v){if (data == NULL) {this->setSizes(v.dims);this->toZero();} for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] /= v.data[i]; return(*this);}
LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator+=(const D & v){for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] += v; return(*this);}
LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator-=(const D & v){for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] -= v; return(*this);}
LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator*=(const D & v){for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] *= v; return(*this);}
LFHTEMP	template<class D> DataGrid<C,nbdim>& DataGrid<C,nbdim>::operator/=(const D & v){for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] /= v; return(*this);}

LFHTEMP	template<class A, class B> DataGrid<C,nbdim>& DataGrid<C,nbdim>::toAddMult(const DataGrid<A, nbdim> & a, const B &b){if (data == NULL) {this->setSizes(a.dims);this->toZero();} for(unsigned int i = totsize()-1;i!=0xFFFFFFFF;i--) data[i] += a.data[i] * b; return(*this);}


LFHTEMP
template<class A_1, class A_2, class C_2, unsigned int size_2> void DataGrid<C, nbdim>::operator() (Oper2<A_1,A_2> const & op, const DataGrid<C_2,size_2> & a_2 ){ // not a match

for(unsigned int i = (*this).totsize()-1; i != 0xFFFFFFFF;i--) (*this)[i](op,a_2[i]);
}

LFHTEMP
template<class C_2, unsigned int size_2> void DataGrid<C, nbdim>::operator() (Oper2<C,C_2> const & op, const DataGrid<C_2,size_2> & a_2){ // match
//setSize(_in.getSize());

for(unsigned int i = (*this).totsize()-1; i != 0xFFFFFFFF;i--) op((*this)[i], a_2[i]);
}


LFHTEMP
template<class A_1, class A_2, class C_2, unsigned int size_2> void DataGrid<C, nbdim>::operator() (Oper2<A_1,A_2> const & op,  DataGrid<C_2,size_2> & a_2 ){ // not a match
for(unsigned int i = (*this).totsize()-1; i != 0xFFFFFFFF;i--) (*this)[i](op,a_2[i]);
}

LFHTEMP
template<class C_2, unsigned int size_2> void DataGrid<C, nbdim>::operator() (Oper2<C,C_2> const & op,  DataGrid<C_2,size_2> & a_2){ // match
//setSize(_in.getSize());
for(unsigned int i = (*this).totsize()-1; i != 0xFFFFFFFF;i--) op(data[i], a_2.data[i]);
}

		LFHTEMP
		template<class A_1, class A_2, class A_3, class C_2, unsigned int nbdim_2, class C_3, unsigned int nbdim_3>  void DataGrid<C, nbdim>::operator() (Oper3<A_1,A_2,A_3> const & op, DataGrid<C_2, nbdim_2>  &a_2, DataGrid<C_3, nbdim_3>  &a_3 ){ // not a match
			//setSize(_in.getSize());
			for(unsigned int i = (*this).totsize()-1; i != 0xFFFFFFFF;i--) (*this)[i](op,a_2[i],a_3[i]);
		}

		LFHTEMP
		template<class C_2, unsigned int nbdim_2, class C_3, unsigned int nbdim_3> void DataGrid<C, nbdim>::operator() (Oper3<C,C_2,C_3> const & op, DataGrid<C_2, nbdim_2>  &a_2, DataGrid<C_3, nbdim_3> &a_3){ // match
			//setSize(_in.getSize());
			for(unsigned int i = (*this).totsize()-1; i != 0xFFFFFFFF;i--) op((*this)[i], a_2[i], a_3[i]);
		}

LFHTEMP	unsigned int const * const DataGrid<C, nbdim>::getSizes() const {return(dims);}

    LFHTEMP void DataGrid<C, nbdim>::setSizes(unsigned int const * const _dims){
    memcpy(dims, _dims, sizeof(unsigned int) * nbdim);
    unsigned int ts = dims[0];
    for(unsigned int i=1;i<nbdim;i++) ts *= dims[i];
//		if (ExCo<C>::NeedsAddLink){
//			}
    if (data != NULL) {
        delete[](data);
    }
    data = new C[ts]; LFH_NICE_ALLOCERROR(data,"Could not allocate for DataGrid::setSizes\n(size= %i)\n",(int)(ts *sizeof(C)))
}
LFHTEMP DataGrid<C, nbdim>& DataGrid<C, nbdim>::setSizes(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE &in_dims){
    unsigned int i;
    unsigned int ts = dims[0] = in_dims[0];
    for(i=1;i<nbdim;i++) ts *= dims[i] = in_dims[i];
//		if (ExCo<C>::NeedsAddLink){
//		}
    if (data != NULL) {
    delete[](data);
    }
    data = new C[ts]; LFH_NICE_ALLOCERROR(data,"Could not allocate for DataGrid::setSizes\n(size= %i)\n", (int)(ts *sizeof(C)))
return *this;}

		LFHTEMP DataGrid<C, nbdim>& DataGrid<C, nbdim>::toZero(){
            if (data != NULL){
			C* cur = &(this->data[0]);
			for(unsigned int i = totsize();i>0;i--) ExOp::toZero(*(cur++)) ;
            }
            return(*this);
        }

		LFHTEMP DataGrid<C, nbdim>& DataGrid<C, nbdim>::toRand(){
			if (data != NULL){
			C* cur = &(this->data[0]);
			for(unsigned int i = totsize();i>0;i--)  ExOp::toRand(*(cur++)) ;
		}	return(*this);
    }

		LFHTEMP DataGrid<C, nbdim>& DataGrid<C, nbdim>::toOne(){
			if (data != NULL){
			C* cur = &((*this)[0]);
			for(unsigned int i = totsize();i>0;i--)  ExOp::toOne(*(cur++));
            }
        return(*this);
    }

LFHTEMP template<class O> DataGrid<O,nbdim> DataGrid<C, nbdim>::operator()(O (*fnct)(const C&)) const{
    DataGrid<O,nbdim> fout;
    fout.setSizes(this->dims);
    C* ptr = data + totsize();
    O* optr = fout.data;
    for(ptr--;data != ptr;ptr--) *(optr++) = fnct(*ptr);
return fout;}
LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator()(C& (*fnct)(C&)){
    C* ptr = data + totsize();
    for(ptr--;data != ptr;ptr--) fnct(*ptr);
    fnct(*ptr);
return *this;}


	LFHTEMP
	void DataGrid<C, nbdim>::initeye(){
		unsigned int modulo=1;
		unsigned int max = dims[nbdim-1];
		unsigned int tsize = max;
		unsigned int i;
		for(i=0;i<nbdim-1;i++){
			modulo = dims[i] * modulo + 1;
			tsize *= dims[i];
			if (max > dims[i]) max = dims[i];
			}
		max = modulo * (max-1);
		C* cur = this->data;
		for(i=0;i<=max;i++) *(cur++) = ((i % modulo) == 0) ? ExCo<C>::mkOne() : ExCo<C>::mkZero() ;
		for(;i<tsize;i++) *(cur++) = ExCo<C>::mkZero() ;
		}
	/*
	LFHTEMP
	Tuple<unsigned int,nbdim> DataGrid<C, nbdim>::getSizes() const{
		int i;
		Tuple<unsigned int,nbdim> &_out;
		for(i=0;i<nbdim;i++) _out[i] = dims[i];
		return(_out);
	}*/


//	LFHTEMP C& DataGrid<C, nbdim>::operator[](unsigned int p ){return(data[p]);}
//	LFHTEMP	C& DataGrid<C, nbdim>::operator[](int p){return(data[p]);}
//	LFHTEMP	C DataGrid<C, nbdim>::operator[](unsigned int p ) const{return(data[p]);}
//	LFHTEMP	C DataGrid<C, nbdim>::operator[](int p) const{return(data[p]);}



LFHTEMP C& DataGrid<C, nbdim>::operator()(unsigned int* p ){
    int pos=0;
    if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + p[i]) * dims[i-1];
    if (nbdim>=8) pos = (pos + p[7]) * dims[6];
    if (nbdim>=7) pos = (pos + p[6]) * dims[5];
    if (nbdim>=6) pos = (pos + p[5]) * dims[4];
    if (nbdim>=5) pos = (pos + p[4]) * dims[3];
    if (nbdim>=4) pos = (pos + p[3]) * dims[2];
    if (nbdim>=3) pos = (pos + p[2]) * dims[1];
    if (nbdim>=2) pos = (pos + p[1]) * dims[0];
return(this->data[pos + p[0]]);}
LFHTEMP C& DataGrid<C, nbdim>::operator()(int* p){
    int pos=0;
    if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + p[i]) * dims[i-1];
    if (nbdim>=8) pos = (pos + p[7]) * dims[6];
    if (nbdim>=7) pos = (pos + p[6]) * dims[5];
    if (nbdim>=6) pos = (pos + p[5]) * dims[4];
    if (nbdim>=5) pos = (pos + p[4]) * dims[3];
    if (nbdim>=4) pos = (pos + p[3]) * dims[2];
    if (nbdim>=3) pos = (pos + p[2]) * dims[1];
    if (nbdim>=2) pos = (pos + p[1]) * dims[0];
return(this->data[pos + p[0]]);}
LFHTEMP C DataGrid<C, nbdim>::operator()(unsigned int* p ) const{
    int pos=0;
    if (nbdim>=9) for(int i=nbdim-1;i>8;i--) pos = (pos + p[i]) * dims[i-1];
    if (nbdim>=8) pos = (pos + p[7]) * dims[6];
    if (nbdim>=7) pos = (pos + p[6]) * dims[5];
    if (nbdim>=6) pos = (pos + p[5]) * dims[4];
    if (nbdim>=5) pos = (pos + p[4]) * dims[3];
    if (nbdim>=4) pos = (pos + p[3]) * dims[2];
    if (nbdim>=3) pos = (pos + p[2]) * dims[1];
    if (nbdim>=2) pos = (pos + p[1]) * dims[0];
return(this->data[pos + p[0]]);}
LFHTEMP C DataGrid<C, nbdim>::operator()(int* p) const{
    int pos=0;
    if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + p[i]) * dims[i-1];
    if (nbdim>=8) pos = (pos + p[7]) * dims[6];
    if (nbdim>=7) pos = (pos + p[6]) * dims[5];
    if (nbdim>=6) pos = (pos + p[5]) * dims[4];
    if (nbdim>=5) pos = (pos + p[4]) * dims[3];
    if (nbdim>=4) pos = (pos + p[3]) * dims[2];
    if (nbdim>=3) pos = (pos + p[2]) * dims[1];
    if (nbdim>=2) pos = (pos + p[1]) * dims[0];
return(this->data[pos + p[0]]);}
LFHTEMP C& DataGrid<C, nbdim>::operator()(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE& p ){
    int pos=0;
    if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,i)) * dims[i-1];
    if (nbdim>=8) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,7)) * dims[6];
    if (nbdim>=7) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,6)) * dims[5];
    if (nbdim>=6) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,5)) * dims[4];
    if (nbdim>=5) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,4)) * dims[3];
    if (nbdim>=4) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,3)) * dims[2];
    if (nbdim>=3) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,2)) * dims[1];
    if (nbdim>=2) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,1)) * dims[0];
return(this->data[pos + Tuple<unsigned int,nbdim>::accessTupleType(p,0)]);}

LFHTEMP
C& DataGrid<C, nbdim>::operator()(const typename Tuple<int,nbdim>::TUPLE_TYPE& p){
    int pos=0;
    if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,i)) * dims[i-1];
    if (nbdim>=8) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,7)) * dims[6];
    if (nbdim>=7) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,6)) * dims[5];
    if (nbdim>=6) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,5)) * dims[4];
    if (nbdim>=5) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,4)) * dims[3];
    if (nbdim>=4) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,3)) * dims[2];
    if (nbdim>=3) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,2)) * dims[1];
    if (nbdim>=2) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,1)) * dims[0];
    return(this->data[pos + Tuple<int,nbdim>::accessTupleType(p,0)]);
}
LFHTEMP
C DataGrid<C, nbdim>::operator()(const typename Tuple<unsigned int,nbdim>::TUPLE_TYPE& p ) const{
    int pos=0;
    if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,i)) * dims[i-1];
    if (nbdim>=8) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,7)) * dims[6];
    if (nbdim>=7) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,6)) * dims[5];
    if (nbdim>=6) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,5)) * dims[4];
    if (nbdim>=5) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,4)) * dims[3];
    if (nbdim>=4) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,3)) * dims[2];
    if (nbdim>=3) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,2)) * dims[1];
    if (nbdim>=2) pos = (pos + Tuple<unsigned int,nbdim>::accessTupleType(p,1)) * dims[0];
return(this->data[pos + Tuple<unsigned int,nbdim>::accessTupleType(p,0)]);}
	LFHTEMP
	C DataGrid<C, nbdim>::operator()(const typename Tuple<int,nbdim>::TUPLE_TYPE& p) const{
        int pos=0;
        if (nbdim>=9) for(int i = nbdim-1;i>8;i--) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,i)) * dims[i-1];
        if (nbdim>=8) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,7)) * dims[6];
        if (nbdim>=7) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,6)) * dims[5];
        if (nbdim>=6) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,5)) * dims[4];
        if (nbdim>=5) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,4)) * dims[3];
        if (nbdim>=4) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,3)) * dims[2];
        if (nbdim>=3) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,2)) * dims[1];
        if (nbdim>=2) pos = (pos + Tuple<int,nbdim>::accessTupleType(p,1)) * dims[0];
        return(this->data[pos + Tuple<int,nbdim>::accessTupleType(p,0)]);
	}

LFHTEMP C DataGrid<C, nbdim>::linearInterpolation(const Tuple<double,nbdim> & coor) const{
    C* mini;
    double buf[nbdim];
    unsigned int offset[nbdim];
    unsigned int ttt;
    C partial[nbdim];

    ttt = (unsigned int)floor(coor[0]);
    if (ttt >= dims[1]-1) {mini = &(data[dims[0]-2]); buf[0] =1.0f;}
    else if (ttt < 0) {mini = data;  buf[0] =0.0f;}
    else {mini = &(data[ttt]);buf[0] = coor[0] - ttt;}
    if (nbdim > 1) {
        offset[1] = dims[0];
        ttt = (unsigned int) floor(coor[1]);
        if (ttt >= dims[1]-1) {ttt = dims[1]-2;  buf[1] =1.0f;}
        else if (ttt < 0){ttt =0;buf[1] = 0.0f;}
        else buf[1] = coor[1] - ttt;
        mini += offset[1] * ttt;
        for(uint32_t i=2;i<nbdim;i++) {
            offset[i] = dims[i-1] * offset[i-1];
            ttt = (unsigned int)floor(coor[i]);
            if (ttt >= dims[i]-1) {ttt = dims[i]-2; buf[i] =1.0f;}
            else if (ttt < 0){ttt =0;buf[i] = 0.0f;}
            else buf[i] = coor[i] - ttt;
            mini += offset[i] * ttt;
        }
        }
    uint32_t d=0;
    for(uint32_t i=0; ;i++){
        if ((i & 1u) == 0u) partial[0] = *(mini++);
        else {partial[0] = (1.0f - buf[0]) * partial[0] + buf[0] *(*(mini--));
            for(d=1u; (d<nbdim)&&((i >> d) & 1);d++) {partial[d] = (1.0f - buf[d]) * partial[d] + buf[d] * partial[d-1]; mini -= offset[d];}
            if (d == nbdim) return(partial[nbdim-1]);
            partial[d] = partial[d-1]; mini += offset[d];
        }
    }
}

	LFHTEMP
	typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1>,C>::TYPE DataGrid<C, nbdim>::operator()(unsigned int a) const{
		typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1>,C>::TYPE f_out;
 		if (nbdim == 1) f_out = *((typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1>,C>::TYPE * ) (data + a));
		else{
			DataGrid<C,nbdim-1> &conv = *((DataGrid<C,nbdim-1>*) (&f_out)) ;
			conv.setSizes(dims);
			unsigned int ts = dims[0];
			unsigned int i;
			for(i =1; i<nbdim-1;i++) ts *= dims[i];
			a = ts * a;
			for(i=0;i<ts;i++) conv.data[i] = data[a+i];
			}
		return f_out;
		}

	LFHTEMP
	typename MT_IFTYPE<nbdim-2, DataGrid<C,nbdim-2>,C>::TYPE DataGrid<C, nbdim>::operator()(unsigned int a,unsigned int b) const{
		typename MT_IFTYPE<nbdim-2, DataGrid<C,nbdim-2>,C>::TYPE f_out;
 		if (nbdim == 2) f_out = *((typename MT_IFTYPE<nbdim-2, DataGrid<C,nbdim-2>,C>::TYPE * ) (data + (a*dims[0] + b)));
		else{
			DataGrid<C,nbdim-2> &conv = *((DataGrid<C,nbdim-2>*) (&f_out)) ;
			conv.setSizes(dims);
			unsigned int ts = dims[0];
			unsigned int i;
			for(i =1; i<nbdim-2;i++) ts *= dims[i];
			a = ts * (a*dims[0] + b);
			for(i=0;i<ts;i++) conv.data[i] = data[a+i];
			}
		return f_out;
		}
	LFHTEMP
	typename MT_IFTYPE<nbdim-3, DataGrid<C,nbdim-3>,C>::TYPE DataGrid<C, nbdim>::operator()(unsigned int a,unsigned int b,unsigned int c) const{
		typename MT_IFTYPE<nbdim-3, DataGrid<C,nbdim-3>,C>::TYPE f_out;
 		if (nbdim == 3) f_out = *((typename MT_IFTYPE<nbdim-3, DataGrid<C,nbdim-3>,C>::TYPE * ) (data + ((a*dims[0] + b)*dims[1] + c)));
		else{
			DataGrid<C,nbdim-3> &conv = *((DataGrid<C,nbdim-3>*) (&f_out)) ;
			conv.setSizes(dims);
			unsigned int ts = dims[0];
			unsigned int i;
			for(i =1; i<nbdim-3;i++) ts *= dims[i];
			a = ts * ((a*dims[0] + b)*dims[1] + c);
			for(i=0;i<ts;i++) conv.data[i] = data[a+i];
			}
		return f_out;
		}
//	LFHTEMP	typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1,LFHVECTOR_REMOTE>,C&>::type DataGrid<C, nbdim>::operator()(unsigned int a){
//		typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1,LFHVECTOR_REMOTE>,C&>::type f_out;
//		return(f_out);
	//	if (nbdim == 1) return data[a]; // *((typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1,LFHVECTOR_REMOTE>,C&>::type * ) (data + a));
	//	else{
	/*		typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1,LFHVECTOR_REMOTE>,C>::type f_out;
			DataGrid<C,nbdim-1,LFHVECTOR_REMOTE> &conv = *((DataGrid<C,nbdim-1,LFHVECTOR_REMOTE>*) (&f_out)) ;
			conv.setSizes(dims);
			unsigned int ts = dims[0];
			unsigned int i;
			for(i =1; i<nbdim-1;i++) ts *= dims[i];
			a = ts * a;
			conv.data = data + a;
			return f_out;*/

			//		}
//		}

	LFHTEMP
	Tuple<unsigned int, nbdim> DataGrid<C, nbdim>::addressof(C const * const item) const{
		Tuple<unsigned int, nbdim> _out;
		unsigned int off = (unsigned int)(item - data);
		_out[0] = off % dims[0];
		if (nbdim>=2) {off /= dims[0]; _out[1] = off % dims[1];}
		if (nbdim>=3) {off /= dims[1]; _out[2] = off % dims[2];}
		if (nbdim>=4) {off /= dims[2]; _out[3] = off % dims[3];}
		if (nbdim>=5) {off /= dims[3]; _out[4] = off % dims[4];}
		if (nbdim>=6) {off /= dims[4]; _out[5] = off % dims[5];}
		if (nbdim>=7) {off /= dims[5]; _out[6] = off % dims[6];}
		if (nbdim>=8) {off /= dims[6]; _out[7] = off % dims[7];}
		return(_out);
	}
/*
	LFHTEMP template<int ssize>
	DataGrid<C, nbdim-ssize,LFHVECTOR_REMOTE> DataGrid<C, nbdim>::operator()(const Tuple<unsigned int,ssize>&p){
		DataGrid<C, nbdim-ssize,LFHVECTOR_REMOTE> _out;
		int i;
		for(i=0;i<nbdim-ssize;i++) _out.dims[0] = dims[0];
		_out.darray = (*this).darray;

		unsigned int pos=0;

		if (nbdim>=8) pos = (pos + p[7]) * dims[6];
		if (nbdim>=7) pos = (pos + p[6]) * dims[5];
		if (nbdim>=6) pos = (pos + p[5]) * dims[4];
		if (nbdim>=5) pos = (pos + p[4]) * dims[3];
		if (nbdim>=4) pos = (pos + p[3]) * dims[2];
		if (nbdim>=3) pos = (pos + p[2]) * dims[1];
		if (nbdim>=2) pos = (pos + p[1]) * dims[0];


		return(_out);
	}*/



LFHTEMP	unsigned int DataGrid<C, nbdim>::getDimSize(unsigned int d) {return(dims[d]);}
LFHTEMP void DataGrid<C, nbdim>::getDims(unsigned int* out_dim){memcpy(out_dim,dims,sizeof(unsigned int)*nbdim);}
LFHTEMP	typename Tuple<unsigned int, nbdim>::TUPLE_TYPE DataGrid<C, nbdim>::getDims() const{ typename Tuple<unsigned int, nbdim>::TUPLE_TYPE fout;
    for(unsigned int i=0;i<nbdim;i++) Tuple<unsigned int, nbdim>::setTupleType(fout,i) = dims[i];
return fout;}
LFHTEMP	void DataGrid<C, nbdim>::makeFreqConvolutionMap(int size, double low, double high){
    unsigned int coors[nbdim];
    size |= 1;
    int i;
    for(i=0;i<nbdim;i++) coors[i] = size;

    setSizes(coors);

    switch(nbdim){
        case 1:

            break;
        case 2:
            for( coors[0] =0;coors[0] < size;coors[0]++){
                for( coors[1] =0;coors[1] < size;coors[1]++){


                }
            }
            break;
    }
}

	LFHTEMP	DataGrid<C,nbdim> DataGrid<C, nbdim>::FFtransform_old() const{
		// assumes C can be multiplied by complex number!
		DataGrid<C,nbdim> _out;
		_out.setSizes(getSizes());
		int curdir,altdim;
		unsigned int coors[nbdim];
		for(curdir =0;curdir<nbdim;curdir++){
			Vector<mycomplex> bluewindow = Vector<mycomplex>::bluesteinWindow(dims[curdir]);
			Vector<C> f_in;

			f_in.setSize(dims[curdir]);
			memset(coors,'\0',sizeof(unsigned int)*nbdim);
			while(true){

				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					if (curdir ==0) f_in[coors[curdir]] = (*this)(coors);
					else f_in[coors[curdir]] = _out(coors);
				}
				Vector<C> f_out = f_in.fourierTransform(bluewindow);

				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					_out(coors) = f_out[coors[curdir]];
				}

				for(altdim=1;altdim<nbdim;altdim++) {
					if (coors[(altdim+curdir) % nbdim] == dims[(altdim+curdir)% nbdim]-1) coors[(altdim+curdir) % nbdim] =0;
					else {coors[(altdim+curdir)% nbdim]++; break;}
				}
				if (altdim == nbdim) break;
			}
		}
		return(_out);
	}

	LFHTEMP	DataGrid<C,nbdim> DataGrid<C, nbdim>::invFFtransform_old() const{
		// assumes C can be multiplied by complex number!
		DataGrid<C,nbdim> _out;
		_out.setSizes(getSizes());
		int curdir,altdim;
		unsigned int coors[nbdim];
		for(curdir =0;curdir<nbdim;curdir++){
			Vector<mycomplex> bluewindow = Vector<mycomplex>::bluesteinWindow(dims[curdir]);
			Vector<C> f_in;
			f_in.setSize(dims[curdir]);
			memset(coors,'\0',sizeof(unsigned int)*nbdim);
			while(true){

				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					if (curdir ==0) f_in[coors[curdir]] = (*this)(coors);
					else f_in[coors[curdir]] = _out(coors);
				}
				Vector<C> f_out = f_in.invfourierTransform(bluewindow);

				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					_out(coors) = f_out[coors[curdir]];
				}

				for(altdim=1;altdim<nbdim;altdim++) {
					if (coors[(altdim+curdir) % nbdim] == dims[(altdim+curdir)% nbdim]-1) coors[(altdim+curdir) % nbdim] =0;
					else {coors[(altdim+curdir)% nbdim]++; break;}
				}
				if (altdim == nbdim) break;
			}
		}

		return(_out);

	}

	LFHTEMP	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> DataGrid<C, nbdim>::FFtransform() const{
		Bluestein bl[nbdim];
		Tuple<unsigned int, nbdim> coor;

		Tuple<double, nbdim> vec;

		unsigned int maxb, tmp;
		for(unsigned int d=0;d<nbdim;d++){
			bl[d].setSize(dims[d]);
			tmp = bl[d].getBufferSize();
			if ((d == 0)||(tmp > maxb)) maxb = tmp;
			}
		typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb]; LFH_NICE_ALLOCERROR(buffer ,"")

		DataGrid< typename ExCo<double>::COMPLEX_TYPE ,nbdim> fout = this->FFtransform_2(bl, (typename  ExCo<double>::COMPLEX_TYPE * )buffer);
        delete[](buffer);

        return fout;
    }
	LFHTEMP	DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> DataGrid<C, nbdim>::invFFtransform() const{
		Bluestein bl[nbdim];
		Tuple<unsigned int, nbdim> coor;
		Tuple<double, nbdim> vec;
		unsigned int maxb, tmp;
		for(unsigned int d=0;d<nbdim;d++){
			bl[d].setSize(dims[d]);
			tmp = bl[d].getBufferSize();
			if ((d == 0)||(tmp > maxb)) maxb = tmp;
			}
		typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb]; LFH_NICE_ALLOCERROR(buffer ,"")

		DataGrid< typename ExCo<double>::COMPLEX_TYPE ,nbdim> fout = this->invFFtransform_2(bl, (typename  ExCo<double>::COMPLEX_TYPE * )buffer);
        delete[](buffer);
        return fout;
    }

	LFHTEMP	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> DataGrid<C, nbdim>::FFtransform_2(Bluestein* blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
		DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> fout[2];
		fout[0] = this->FFtransform_2dir(0, blue[0], buffer);
		for(unsigned int curdir =1;curdir<nbdim;curdir++) fout[(curdir & 1)] = fout[(1 ^ (curdir & 1))].FFtransform_2dir(curdir, blue[curdir], buffer);
		return fout[(1 ^ (nbdim & 1))];
	}

	LFHTEMP	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> DataGrid<C, nbdim>::FFtransform_2dir(unsigned int curdir, Bluestein& blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> _out;

			_out.setSizes(dims);
			unsigned int altdim;
			unsigned int coors[nbdim];
			ExOp::toZero(coors);
				while(true){
					for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){

						buffer[coors[curdir]] = typename ExCo<C>::COMPLEX_TYPE( (*this)(coors) );
					}
					blue.toFT(buffer);

					for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
						_out(coors) = buffer[coors[curdir]];
					}

					for(altdim=1;altdim<nbdim;altdim++) {
						if (coors[(altdim+curdir) % nbdim] == dims[(altdim+curdir)% nbdim]-1) coors[(altdim+curdir) % nbdim] =0;
						else {coors[(altdim+curdir)% nbdim]++; break;}
					}
					if (altdim == nbdim) break;
			}
			return(_out);
		}

LFHTEMP	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> DataGrid<C, nbdim>::invFFtransform_2(Bluestein* blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> fout[2];
	fout[0] = this->invFFtransform_2dir(0, blue[0], buffer);
	for(unsigned int curdir =1;curdir<nbdim;curdir++) fout[(curdir & 1)] = fout[(1 ^ (curdir & 1))].invFFtransform_2dir(curdir, blue[curdir], buffer);

return fout[(1 ^ (nbdim & 1))];
}

	LFHTEMP	DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> DataGrid<C, nbdim>::invFFtransform_2dir(unsigned int curdir, Bluestein& blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
			// assumes C can be multiplied by complex number!
			DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> _out;

			_out.setSizes(dims);
			unsigned int altdim;
			unsigned int coors[nbdim];


			memset(coors,'\0',sizeof(unsigned int)*nbdim);
			while(true){
					for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
						buffer[coors[curdir]] = (*this)(coors);
					}

				 blue.toIFT(buffer);

				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					_out(coors) = buffer[coors[curdir]];
				}

				for(altdim=1;altdim<nbdim;altdim++) {
					if (coors[(altdim+curdir) % nbdim] == dims[(altdim+curdir)% nbdim]-1  ) coors[(altdim+curdir) % nbdim] =0;
					else {coors[(altdim+curdir)% nbdim]++; break;}
				}
				if (altdim == nbdim) break;
		}
		return(_out);
	}


LFHTEMP	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> DataGrid<C, nbdim>::FFtransform_zeroborder(Bluestein* blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> _out;

	int curdir,altdim;
	unsigned int coors[nbdim];
	for(curdir =0;curdir<nbdim;curdir++) coors[curdir] = dims[curdir] *2;

	_out.setSizes(coors);

	for(curdir =0;curdir<nbdim;curdir++){


		memset(coors,'\0',sizeof(unsigned int)*nbdim);
		while(true){
			if (curdir == 0)
				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					buffer[coors[curdir]] = (*this)(coors);
				}
			else{
				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					buffer[coors[curdir]] = _out(coors);
				}
			}

			for(;coors[curdir]<dims[curdir]*2;coors[curdir]++){
				ExOp::toZero(buffer[coors[curdir]]);
			}

			blue[curdir].toFT(buffer);

			for(coors[curdir] =0;coors[curdir]<dims[curdir]*2;coors[curdir]++){
				_out(coors) = buffer[coors[curdir]];
			}

			for(altdim=1;altdim<nbdim;altdim++) {
				if (coors[(altdim+curdir) % nbdim] == ((altdim+curdir >= nbdim) ? _out.dims[(altdim+curdir)% nbdim] :  dims[(altdim+curdir)% nbdim])-1) coors[(altdim+curdir) % nbdim] =0;
				else {coors[(altdim+curdir)% nbdim]++; break;}
			}
			if (altdim == nbdim) break;
		}	}
	return(_out);
}

LFHTEMP	DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> DataGrid<C, nbdim>::invFFtransform_zeroborder(Bluestein* blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
	// assumes C can be multiplied by complex number!
	DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> _out;
	int curdir,altdim;
	unsigned int coors[nbdim];
	_out.setSizes(dims);



	for(curdir =0;curdir<nbdim;curdir++){

		memset(coors,'\0',sizeof(unsigned int)*nbdim);
		while(true){
			if (curdir == 0)
				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					buffer[coors[curdir]] = (*this)(coors);
				}
			else
				for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
					buffer[coors[curdir]] = _out(coors);
				}
			blue[curdir].toIFT(buffer);

			for(coors[curdir] =0;coors[curdir]<_out.dims[curdir];coors[curdir]++){
				_out(coors) = buffer[coors[curdir]];
			}

			for(altdim=1;altdim<nbdim;altdim++) {
				if (coors[(altdim+curdir) % nbdim] == ((altdim+curdir >= nbdim) ? dims[(altdim+curdir)% nbdim]/2 :  dims[(altdim+curdir)% nbdim]) -1) coors[(altdim+curdir) % nbdim] =0;
				else {coors[(altdim+curdir)% nbdim]++; break;}
			}
			if (altdim == nbdim) break;
		}	}

	DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> _out2;
	for(curdir =0;curdir<nbdim;curdir++) coors[curdir] = dims[curdir]/2;
	_out2.setSizes(coors);

	typename DataGrid<C, nbdim>::KeyIterator ite(_out2);

	if (ite.first()) do{
		_out2(ite()) = _out(ite());
	}while(ite.next());

	return(_out2);
}


		LFHTEMP	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> DataGrid<C, nbdim>::FFtransform_zeroborder_dir(unsigned int curdir, Bluestein &blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> _out;

			int altdim;
			unsigned int coors[nbdim];
			for(altdim=0;altdim<nbdim;altdim++) coors[altdim] = (altdim == curdir) ? dims[altdim] *2: dims[altdim] ;

			_out.setSizes(coors);

				memset(coors,'\0',sizeof(unsigned int)*nbdim);
				while(true){
					for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
						buffer[coors[curdir]] = (*this)(coors);
					}

					for(;coors[curdir]<_out.dims[curdir];coors[curdir]++){
						ExOp::toZero(buffer[coors[curdir]]);
					}
					blue.toFT(buffer);

					for(coors[curdir] =0;coors[curdir]<_out.dims[curdir];coors[curdir]++){
						_out(coors) = buffer[coors[curdir]];
					}

					for(altdim=1;altdim<nbdim;altdim++) {
						if (coors[(altdim+curdir) % nbdim] == ((altdim+curdir >= nbdim) ? _out.dims[(altdim+curdir)% nbdim] :  dims[(altdim+curdir)% nbdim])-1) coors[(altdim+curdir) % nbdim] =0;
						else {coors[(altdim+curdir)% nbdim]++; break;}
					}
					if (altdim == nbdim) break;
				}
			return(_out);
		}

		LFHTEMP	DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> DataGrid<C, nbdim>::invFFtransform_zeroborder_dir(unsigned int curdir, Bluestein &blue, typename ExCo<C>::COMPLEX_TYPE * buffer) const{
			// assumes C can be multiplied by complex number!
			DataGrid<typename ExCo<C>::COMPLEX_TYPE,nbdim> _out;
			int altdim;
			unsigned int coors[nbdim];
			coors[curdir] = dims[curdir]/2;
			_out.setSizes(dims);
				memset(coors,'\0',sizeof(unsigned int)*nbdim);
				while(true){
					for(coors[curdir] =0;coors[curdir]<dims[curdir];coors[curdir]++){
						buffer[coors[curdir]] = (*this)(coors);
					}
					blue.toIFT(buffer);

					for(coors[curdir] =0;coors[curdir]<_out.dims[curdir];coors[curdir]++){
						_out(coors) = buffer[coors[curdir]];
					}

					for(altdim=1;altdim<nbdim;altdim++) {
						if (coors[(altdim+curdir) % nbdim] == ((altdim+curdir >= nbdim) ? dims[(altdim+curdir)% nbdim]/2 :  dims[(altdim+curdir)% nbdim]) -1) coors[(altdim+curdir) % nbdim] =0;
						else {coors[(altdim+curdir)% nbdim]++; break;}
					}
					if (altdim == nbdim) break;
				}
			return(_out);
		}

	LFHTEMP
	void DataGrid<C, nbdim>::applycircleFilter(double min, double max){

		unsigned int coors[nbdim];
		memset(coors,'\0',sizeof(unsigned int)*nbdim);
		double dists2[nbdim];
		memset(dists2,'\0',sizeof(double)*nbdim);
	//	double dimnorm=0;
		double tmp=0;
		int altdim;
	//	for(altdim=0;altdim<nbdim;altdim++) dimnorm += dims[altdim]*dims[altdim];
		double trmax = 1.0f / (min * min);
		double trmin = 1.0f / (max * max);
	//	printf("%f,%f\n",trmin,trmax);

		C tmpzero = ExCo<C>::zero();

		while(true){
	//		printf("%f\n",dists2[0]);
			if (((dists2[0] < trmin)&&(dists2[0] != 0.0f))||(dists2[0]> trmax)) {
				(*this)(coors) = tmpzero;//
			}



			for(altdim=0;altdim<nbdim;altdim++) {
				if (coors[altdim] == dims[altdim]-1) coors[altdim] =0;
				else {coors[altdim]++;
					tmp =  ((double)((coors[altdim] < (dims[altdim] >> 1) ? coors[altdim] : dims[altdim] - coors[altdim]))) / dims[altdim];

					dists2[altdim] = tmp*tmp + ((altdim == nbdim-1) ? 0.0f: dists2[altdim+1]);
					for(altdim--;altdim<nbdim;altdim--) dists2[altdim] = dists2[altdim+1];
				break;
				}
			}
			if (altdim == nbdim) break;
		}
	}

	LFHTEMP
	DataGrid<C, nbdim> DataGrid<C, nbdim>::circleFilter(double min, double max){
		DataGrid<C,nbdim> _out;
		_out.setSizes(getSizes());

		unsigned int coors[nbdim];
		memset(coors,'\0',sizeof(unsigned int)*nbdim);
		double dists2[nbdim];
		memset(dists2,'\0',sizeof(double)*nbdim);
		//	double dimnorm=0;
		double tmp=0;
		unsigned int altdim;
		//	for(altdim=0;altdim<nbdim;altdim++) dimnorm += dims[altdim]*dims[altdim];
		double trmax = 1.0f / (min * min);
		double trmin = 1.0f / (max * max);
		//	printf("%f,%f\n",trmin,trmax);

		C tmpzero = ExCo<C>::zero();

		while(true){
			//		printf("%f\n",dists2[0]);
			if (((dists2[0] < trmin)&&(dists2[0] != 0.0f))||(dists2[0]> trmax)) {
				_out(coors) = tmpzero;//
			}else{
				_out(coors) = (*this)(coors);
			}



			for(altdim=0;altdim<nbdim;altdim++) {
				if (coors[altdim] == dims[altdim]-1) coors[altdim] =0;
				else {coors[altdim]++;
					tmp =  ((double)((coors[altdim] < (dims[altdim] >> 1) ? coors[altdim] : dims[altdim] - coors[altdim]))) / dims[altdim];

					dists2[altdim] = tmp*tmp + ((altdim == nbdim-1) ? 0.0f: dists2[altdim+1]);
					for(altdim--;altdim<nbdim;altdim--) dists2[altdim] = dists2[altdim+1];
					break;
				}
			}
			if (altdim == nbdim) break;
		}

		return(_out);
	}
#ifdef Rcpp_hpp
LFHTEMP template<class RCL> void DataGrid<C, nbdim>::rdMatrix(const arma::Mat<RCL> &where){
    if (nbdim != 2) return;
    Tuple<uint32_t,nbdim> coor;
    int i,j;
    coor[0] = where.n_rows;
    coor[1] = where.n_cols;
    this->setSizes(coor);
    C* cur = data;
    for(coor[1]=0;coor[1]<dims[1];coor[1]++){
        for(coor[0]=0;coor[0]<dims[0];coor[0]++) *(cur++) = where.at(coor[0],coor[1]);
    }
}
LFHTEMP template<class RCL> void DataGrid<C, nbdim>::wrMatrix(arma::Mat<RCL> &where)const{
    if (nbdim != 2) return;
    where = arma::Mat<RCL>(dims[0],dims[1]);
    const C* cur = data;
    Tuple<uint32_t,nbdim> coor;
    for(coor[1]=0;coor[1]<dims[1];coor[1]++){
        for(coor[0]=0;coor[0]<dims[0];coor[0]++) where.at(coor[0],coor[1]) = *(cur++);
    }
}
#endif

LFHTEMP void DataGrid<C, nbdim>::show(FILE* f, int level) const{
	unsigned int x,y;
	if (nbdim == 1){
		for(x=0;x<dims[0];x++) {
			ExOp::show(this->data[x] , f,1); fprintf(f,"%c", x == dims[0]-1 ? '\n' : '\t');
		}
	}else if (nbdim == 2){
		for(y=0;y<dims[1];y++) {
			for(x=0;x<dims[0];x++) {
			ExOp::show(this->data[x + y * dims[0]] , f,1); fprintf(f,"%c", x == dims[0]-1 ? '\n' : '\t');
			}
		}
	}
}

LFHTEMP ERRCODE DataGrid<C, nbdim>::save(FILE* f) const{
    if (fwrite(dims,sizeof(unsigned int), nbdim,f) != nbdim) return 1;
	unsigned int tsize = totsize();
	if (ExCo<C>::IsPOD) return (fwrite(this->data,sizeof(C), tsize,f) != tsize) ? 1 : 0;
    ERRCODE fout = ExOp::save(this->data[0],f);
    for(unsigned int i=1;i<tsize;i++) fout |= ExOp::save(this->data[i],f);
    return fout;
}

LFHTEMP	ERRCODE DataGrid<C, nbdim>::load(FILE* f, unsigned int lenght){
	unsigned int dbuf[nbdim];
	if (fread(dbuf,sizeof(unsigned int), nbdim,f) != nbdim) return 1;
	setSizes(dbuf);
	unsigned int tsize = totsize();
	if (ExCo<C>::IsPOD)	fread(this->data,sizeof(C), tsize,f);
	else{
		for(unsigned int i=0;i<tsize;i++) ExOp::load(this->data[i],f);
	}
	return 0;
}

LFHTEMP DataGrid< double ,nbdim> DataGrid<C, nbdim>::ExpectedDistanceToOther(double exit_prob) const{
	DataGrid< double ,nbdim> fout;
	fout.setSizes(dims);
	// no HMM, fast calculation
	ExOp::toZero(fout);

	unsigned int mincontact =0;
	unsigned int i,j,k;
	for(i =0;i<nbdim;i++) if (dims[i] > 1) mincontact++;

	Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
	Tuple<Tuple<unsigned int, nbdim>, nbdim * (nbdim -1) * 4> ncoor_knight;


	HeapTree< KeyElem< double, Tuple<unsigned int, nbdim> > > list_coor;

	KeyElem< double, Tuple<unsigned int, nbdim> > cur_in;
	KeyElem< double, Tuple<unsigned int, nbdim> > cur_danew;

	typename DataGrid<C, nbdim>::KeyIterator ite(*this);

	if (ite.first()) do{
		int nbind = (*this).get_indirectNeightbor(ite(), ncoor);
		int nbkni = (*this).get_knightNeightbor(ite(), ncoor_knight);
		C curv = (*this)(ite());
		cur_in.k = 1.0f /  exit_prob;
		for(i=0;i<nbind;i++) if (curv != (*this)(ncoor[i])) {
			k=0;
			for(j=0;j<nbdim;j++) {if (ite()[j] != ncoor[i][j]) k++;}
			if (cur_in.k > sqrt((double)k)) cur_in.k = sqrt((double)k);
		}
		for(i=0;i<nbkni;i++) if (curv != (*this)(ncoor_knight[i])) {
			if (cur_in.k > sqrt(5.0f)) cur_in.k = sqrt(5.0f);
		}
		fout(ite()) = cur_in.k;

		if (cur_in.k != (1.0f /  exit_prob)){
			cur_in.d = ite();
			list_coor.insert(cur_in);

		}

	} while(ite.next());

	while(!list_coor.isEmpty()){
		cur_in = list_coor.pop();
		if (cur_in.k != fout(cur_in.d)) continue; // got processed already

		int nbind = (*this).get_indirectNeightbor(cur_in.d, ncoor);
		int nbkni = (*this).get_knightNeightbor(cur_in.d, ncoor_knight);


		for(i=0;i<nbind;i++) {
			k=0;
			for(j=0;j<nbdim;j++) {if (ite()[j] != ncoor[i][j]) k++;}
			if (fout(ncoor[i]) > (cur_danew.k = cur_in.k + sqrt((double)k)) ) {cur_danew.d = ncoor[i]; fout(ncoor[i]) = cur_danew.k; list_coor.insert(cur_danew);}
		}
		for(i=0;i<nbkni;i++) {
			if (fout(ncoor_knight[i]) > (cur_danew.k = cur_in.k + sqrt(5.0f)) ) {cur_danew.d = ncoor_knight[i]; fout(ncoor_knight[i]) = cur_danew.k; list_coor.insert(cur_danew);}
		}
	}
	return(fout);
}


LFHTEMP	template<class ARR_double>
void DataGrid<C, nbdim>::ExpectedDistanceToOther(const DataGrid< ARR_double ,nbdim> &source, const DataGrid<double,2>& transition){
	setSizes(source.dims);
	int k;
	int nbstates = transition.dims[0];

	ExOp::toZero(*this);
	for(k=0;k<nbstates;k++){
		DataGrid< double,nbdim> partial;
		partial.ExpectedDistanceToOther_instate(source, transition, k);
		KeyIterator ite = KeyIterator(*this);
		if (ite.first()) do{
			(*this)(ite()) += source(ite())[k] * partial(ite());
		} while (ite.next());
	}
}
LFHTEMP	template<class ARR_double>
void DataGrid<C, nbdim>::ExpectedDistanceToOther_instate(const DataGrid< ARR_double ,nbdim> &source,const DataGrid<double,2>&  transition, int state){
	setSizes(source.dims);
	double max =0.0f;
	double tmp;
	int i,j,k,l;
	int nbstates = transition.dims[0];

	for(i=0;i<nbdim;i++) max += ((double)dims[i])* dims[i];
	max = sqrt(max);

	Vector< KeyElem< double, Tuple<unsigned int, nbdim> > > list;
	KeyIterator ite = KeyIterator(*this);
	if (ite.first()) do{
		list.push_back(KeyElem<double, Tuple<unsigned int, 2> >( ((double*)source(ite()))[state], ite()));
		(*this)(ite()) = max;
	} while (ite.next());

	list.sort();

	double transit[nbstates][1+ nbdim];

	for(i=0;i<nbdim;i++) {
		transit[i][0] = transition.data[state + j * nbstates]; // == 1.0f - exp(alpha)
		transit[i][nbdim] = log(1.0f - transition.data[state + j * nbstates]); // == 1.0f - exp(alpha)
		for(j=1;j<nbdim;j++) {
			transit[i][j] = 1.0f - exp(sqrt(1.0f + j) * transit[i][nbdim]); // exp(- alpha)
		}
		transit[i][nbdim] = 1.0f - exp(sqrt(5.0f) * transit[i][nbdim]);
	}

	for(j=0;j<nbstates;j++) {
		transit[state][j] = 1.0f;
		for(i=0;i<nbdim;i++) {
			if (i != state) transit[state][j] -= transit[i][j];
		}
	}


	Tuple<Tuple<unsigned int, 2>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
	Tuple<Tuple<unsigned int, 2>, nbdim * (nbdim -1) * 4> ncoor_knight;

	for(i=0;i<list.getSize();i++){
		int nbind = (*this).get_indirectNeightbor(list[i].d, ncoor);
		int nbkni = (*this).get_knightNeightbor(list[i].d, ncoor_knight);

		for(j=0;j<nbind;j++){
			l =0;
			for(k<0;k<nbdim;k++) if (ncoor[j][k] != list[i].d[k]) l++;
			max = 0.0f;
			tmp = source(list[i].d)[state] * transit[state][l-1];
			if (tmp != 0.0f){
			for(k<0;k<nbstates;k++) if (k != state) max += source(ncoor[j])[k] * transit[k][l-1];
			max = (*this)(list[i].d) / ( 1.0f + (max / tmp));
			if ( (*this)(ncoor[j]) > max)  (*this)(ncoor[j]) = max;
			}
		}
		for(j=0;j<nbkni;j++){
			tmp = source(list[i].d)[state] * transit[state][nbdim];
			if (tmp != 0.0f){
				for(k<0;k<nbstates;k++) if (k != state) max += source(ncoor_knight[j])[k] * transit[k][nbdim];
				max = (*this)(list[i].d) / ( 1.0f + (max / tmp));
				if ( (*this)(ncoor_knight[j]) > max)  (*this)(ncoor[j]) = max;
			}
		}
	}
}
LFHTEMP void DataGrid<C, nbdim>::ExpectedDistanceToOther(const DataGrid< double ,nbdim+1> &source, const DataGrid<double,2>& transition, bool state_zero_special){
	Tuple<unsigned int, nbdim> d;
	int k,l;
	for(k= 0;k<nbdim;k++) d[k] = source.dims[k+1];
	setSizes(d);

	int nbstates = transition.dims[0];
	if ((transition.dims[0] != transition.dims[1])||(transition.dims[0] != source.dims[0])) {fprintf(stderr,"Illegal Sizes!\n"); exit(1);}
	Tuple<unsigned int,nbdim+1> O_o;
	ExOp::toZero(*this);
	KeyIterator ite = KeyIterator(*this);
	for(k= (state_zero_special) ? 1 : 0;k<nbstates;k++){
		DataGrid< double,nbdim> partial;
		partial.ExpectedDistanceToOther_instate(source, transition, k);
		O_o[0] = k;
		if (ite.first()) do{
			for(l=0;l<nbdim;l++) O_o[l+1] = ite()[l];
			if  (ExCo<double>::isValid(partial(ite()))) (*this)(ite()) += (source(O_o) * partial(ite()));
			else (*this)(ite()) += source(O_o) * 256.0f;
			} while (ite.next());

	}
}
		/*
	LFHTEMP
	void DataGrid<C, nbdim>::ExpectedDistanceToOther_instate(const DataGrid< double ,nbdim+1> &source,const DataGrid<double,2>&  transition, int state){
		Tuple<unsigned int, nbdim> d;
		int i,j,k,l,z;
		HeapTree< defaultHeapBackPointer<double>* > pixelheap;
		DataGrid< defaultHeapBackPointer<double>, nbdim > pixelbuffer;
		for(k= 0;k<nbdim;k++) d[k] = source.dims[k+1];
		setSizes(d);
		pixelbuffer.setSizes(d);

		double max =0.0f;
		double tmp;

		int nbstates = transition.dims[0];
		Tuple<unsigned int,nbdim+1> coor_o;
//		Tuple<unsigned int,nbdim+1> coor_i;

		for(i=0;i<nbdim;i++) max += ((double)dims[i])* dims[i];
		max = sqrt(max);

		KeyIterator ite = KeyIterator(*this);
		unsigned int coor_mat[2];

		coor_mat[0] = state; // never changes!

		unsigned int max_dist=0; //
		for(i=0;i<nbdim;i++)  max_dist += source.dims[i+1];

		double decay; // P(Oi+1 != s | Ow = s)
		double transit[nbstates][1+ nbdim];
		double denum;

		coor_o[0] = state;
		if (ite.first()) do{

			for(i=0;i<nbdim;i++) coor_o[i+1] = ite()[i];
			coor_o[0] = state;
			coor_mat[1] = state;
  		    decay = transition(coor_mat) * source(coor_o);
			if (decay != 0.0f){
				denum = 0.0f;
			for(coor_o[0]=0;coor_o[0]<source.dims[0];coor_o[0]++){
				if (coor_o[0] == state) continue;
				coor_mat[1] = coor_o[0];
				denum += transition(coor_mat) * source(coor_o);
			}
			decay = 16.0f * decay / denum; // expected number of transitions, multiplied by 16
			}
			if ((decay > max_dist)||(!ExCo<double>::isValid(decay))) decay = max_dist;
			coor_o[0] = state;
			(*this)(curcoor) = pixelbuffer(curcoor).data;
			pixelbuffer(ite()).data = decay * source(coor_o);
			//printf("b: %e\n",decay);

			pixelheap.insert(&(pixelbuffer(ite())));
		} while (ite.next());




		for(i=0;i<nbstates;i++) {
			transit[i][0] = transition[state + i * nbstates]; // == 1.0f - exp(alpha)
			transit[i][nbdim] = log(1.0f - transition[state + i * nbstates]); // == 1.0f - exp(alpha)
			for(j=1;j<nbdim;j++) {
				transit[i][j] = 1.0f - exp(sqrt(1.0f + j) * transit[i][nbdim]); // exp(- alpha)
			}
			transit[i][nbdim] = 1.0f - exp(sqrt(5.0f) * transit[i][nbdim]);
		}

		for(j=0;j<nbdim;j++) {
			transit[state][j] = 1.0f;
			for(i=0;i<nbstates;i++) {
				if (i != state) transit[state][j] -= transit[i][j];
			}
		}

	//	for(j=0;j<=nbdim;j++) printf("%f\t%f\t%f\n", transit[j][0], transit[j][1], transit[j][2]);
		Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
		Tuple<Tuple<unsigned int, nbdim>, nbdim * (nbdim -1) * 4> ncoor_knight;

		double cmpeddist[nbdim+1];
//		printf("sizes %i\t%i\n",pixelbuffer.dims[0], pixelbuffer.dims[1]);fflush(stdout);
//		printf("sizes %i\t%i\n",(*this).dims[0], (*this).dims[1]);fflush(stdout);

		while(!(pixelheap.isEmpty())){
			Tuple<unsigned int, nbdim> curcoor = pixelbuffer.addressof(&(pixelheap.pop()));
	//		printf("%i\t%i\n",curcoor[0], curcoor[1]);fflush(stdout);
			pixelbuffer(curcoor).backptr = -1;

			(*this)(curcoor) = pixelbuffer(curcoor).data;

			for(z=0;z<nbdim;z++) coor_o[z+1] = curcoor[z];
			for(l=0;l<= nbdim;l++) {
				coor_o[0] = state;
				tmp = source(coor_o) * transit[state][0];
				cmpeddist[l] =0.0f;
				if (tmp !=0.0f){
				for(coor_o[0]=0;coor_o[0]<nbstates;coor_o[0]++) if (coor_o[0] != state) cmpeddist[l] += source(coor_o) * transit[coor_o[0]][0];
					cmpeddist[l] = ((*this)(curcoor)) * tmp / (tmp + cmpeddist[l]);
//					printf("transit %e\n", tmp / (tmp + cmpeddist[l]));
					if (!ExCo<double>::isValid(cmpeddist[l])) cmpeddist[l] =((*this)(curcoor));
				}
				cmpeddist[l] += ((l ==nbdim) ? sqrt(5) : sqrt((double)(l+1)));
			}
		//	printf("dadists:  %e\t%e\t%e\n", cmpeddist[0],cmpeddist[1],cmpeddist[2]);

			int nbind = (*this).get_indirectNeightbor(curcoor, ncoor);
			int nbkni = (*this).get_knightNeightbor(curcoor, ncoor_knight);

			for(j=0;j<nbind;j++){
				if (pixelbuffer(ncoor[j]).backptr == -1) continue;
				l =0;
				for(k=0;k<nbdim;k++) if (ncoor[j][k] != curcoor[k]) l++;
			//	printf("%f <= %f, %i\n",pixelbuffer(ncoor[j]).data, cmpeddist[l-1], l);
				if (pixelbuffer(ncoor[j]).data > cmpeddist[l-1]) {pixelbuffer(ncoor[j]).data = cmpeddist[l-1]; pixelheap.update(pixelbuffer(ncoor[j]).backptr);}
			}
			for(j=0;j<nbkni;j++){
				if (pixelbuffer(ncoor_knight[j]).backptr == -1) continue;
				if (pixelbuffer(ncoor_knight[j]).data > cmpeddist[nbdim]) {pixelbuffer(ncoor_knight[j]).data = cmpeddist[nbdim]; pixelheap.update(pixelbuffer(ncoor_knight[j]).backptr);}
			}

		}




	}*/

		 LFHTEMP
		 void DataGrid<C, nbdim>::ExpectedDistanceToOther_instate(const DataGrid< double ,nbdim+1> &source,const DataGrid<double,2>&  transition, int state){
		 Tuple<unsigned int, nbdim> d;
		 int i,j,k,l,z;
		 HeapTree< defaultHeapBackPointer<double>* > pixelheap;
		 DataGrid< defaultHeapBackPointer<double>, nbdim > pixelbuffer;
		 for(k= 0;k<nbdim;k++) d[k] = source.dims[k+1];
		 setSizes(d);
		 pixelbuffer.setSizes(d);

		 double max =0.0f;
		 double tmp;

		 int nbstates = transition.dims[0];
		 Tuple<unsigned int,nbdim+1> coor_o;
		 //		Tuple<unsigned int,nbdim+1> coor_i;





		 for(i=0;i<nbdim;i++) max += ((double)dims[i])* dims[i];
		 max = sqrt(max);

		 KeyIterator ite = KeyIterator(*this);
		 unsigned int coor_mat[2];

		 coor_mat[0] = state; // never changes! (trying to leave that one state!)

		 unsigned int max_dist=0; //
		 for(i=0;i<nbdim;i++)  max_dist += source.dims[i+1];

		 double decay; // P(Oi+1 != s | Ow = s)
		 double transit[nbstates][1+ nbdim];
		 double denum;

		 coor_o[0] = state;
		 if (ite.first()) do{

		 for(i=0;i<nbdim;i++) coor_o[i+1] = ite()[i];
		 coor_o[0] = state;
		 coor_mat[1] = state;
		 decay = transition(coor_mat) * source(coor_o);
		 if (decay != 0.0f){ // there is some probability of being in the state
		 denum = 0.0f;
		 for(coor_o[0]=0;coor_o[0]<source.dims[0];coor_o[0]++){
		 if (coor_o[0] == state) continue;
		 coor_mat[1] = coor_o[0];
		 denum += transition(coor_mat) * source(coor_o);
		 }
		 decay = 1.0f + (decay / denum); // expected number of transitions, multiplied by 16
		 }
		 if ((decay > max_dist)||(!ExCo<double>::isValid(decay))) decay = max_dist;

		 pixelbuffer(ite()).data = decay;
		 //printf("b: %e\n",decay);

		 pixelheap.insert(&(pixelbuffer(ite())));

			 (*this)(ite()) = pixelbuffer(ite()).data;
		 } while (ite.next());




		 for(i=0;i<nbstates;i++) {
		 transit[i][0] = transition.data[state + i * nbstates]; // == 1.0f - exp(alpha)
		 transit[i][nbdim] = log(1.0f - transition.data[state + i * nbstates]); // == 1.0f - exp(alpha)
		 for(j=1;j<nbdim;j++) {
		 transit[i][j] = 1.0f - exp(sqrt(1.0f + j) * transit[i][nbdim]); // exp(- alpha)
		 }
		 transit[i][nbdim] = 1.0f - exp(sqrt(5.0f) * transit[i][nbdim]);
		 }

		 for(j=0;j<nbdim;j++) {
		 transit[state][j] = 1.0f;
		 for(i=0;i<nbstates;i++) {
		 if (i != state) transit[state][j] -= transit[i][j];
		 }
		 }

		 //	for(j=0;j<=nbdim;j++) printf("%f\t%f\t%f\n", transit[j][0], transit[j][1], transit[j][2]);
		 Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
		 Tuple<Tuple<unsigned int, nbdim>, nbdim * (nbdim -1) * 4> ncoor_knight;

		 double cmpeddist[nbdim+1];
		 //		printf("sizes %i\t%i\n",pixelbuffer.dims[0], pixelbuffer.dims[1]);fflush(stdout);
		 //		printf("sizes %i\t%i\n",(*this).dims[0], (*this).dims[1]);fflush(stdout);

		 while(!(pixelheap.isEmpty())){
		 Tuple<unsigned int, nbdim> curcoor = pixelbuffer.addressof(&(pixelheap.pop()));
		 //		printf("%i\t%i\n",curcoor[0], curcoor[1]);fflush(stdout);
		 pixelbuffer(curcoor).backptr = -1;

		 (*this)(curcoor) = pixelbuffer(curcoor).data;

		 for(z=0;z<nbdim;z++) coor_o[z+1] = curcoor[z];
		 for(l=0;l<= nbdim;l++) {
		 coor_o[0] = state;
		 tmp = source(coor_o) * transit[state][0];
		 cmpeddist[l] =0.0f;
		 if (tmp !=0.0f){
		 for(coor_o[0]=0;coor_o[0]<nbstates;coor_o[0]++) if (coor_o[0] != state) cmpeddist[l] += source(coor_o) * transit[coor_o[0]][0];
		 cmpeddist[l] = ((*this)(curcoor)) * tmp / (tmp + cmpeddist[l]);
		 //					printf("transit %e\n", tmp / (tmp + cmpeddist[l]));
		 if (!ExCo<double>::isValid(cmpeddist[l])) cmpeddist[l] =((*this)(curcoor));
		 }
		 cmpeddist[l] += ((l ==nbdim) ? sqrt(5) : sqrt((double)(l+1)));
		 }
		 //	printf("dadists:  %e\t%e\t%e\n", cmpeddist[0],cmpeddist[1],cmpeddist[2]);

		 int nbind = (*this).get_indirectNeightbor(curcoor, ncoor);
		 int nbkni = (*this).get_knightNeightbor(curcoor, ncoor_knight);

		 for(j=0;j<nbind;j++){
		 if (pixelbuffer(ncoor[j]).backptr == -1) continue;
		 l =0;
		 for(k=0;k<nbdim;k++) if (ncoor[j][k] != curcoor[k]) l++;
		 //	printf("%f <= %f, %i\n",pixelbuffer(ncoor[j]).data, cmpeddist[l-1], l);
		 if (pixelbuffer(ncoor[j]).data > cmpeddist[l-1]) {pixelbuffer(ncoor[j]).data = cmpeddist[l-1]; pixelheap.update(pixelbuffer(ncoor[j]).backptr);}
		 }
		 for(j=0;j<nbkni;j++){
		 if (pixelbuffer(ncoor_knight[j]).backptr == -1) continue;
		 if (pixelbuffer(ncoor_knight[j]).data > cmpeddist[nbdim]) {pixelbuffer(ncoor_knight[j]).data = cmpeddist[nbdim]; pixelheap.update(pixelbuffer(ncoor_knight[j]).backptr);}
		 }

		 }




		 }





	LFHTEMP
	void DataGrid<C, nbdim>::ExpectedDistanceToOther_I(const DataGrid< double ,nbdim> &source, const DataGrid<double,2>& transition){
		setSizes(source.dims);
		int k;
		int nbstates = transition.dims[0];

		ExOp::toZero(*this);
		for(k=0;k<nbstates;k++){
			DataGrid< double,nbdim> partial;
			partial.ExpectedDistanceToOther_instate(source, transition, k);
			KeyIterator ite = KeyIterator(*this);
			if (ite.first()) do{
				(*this)(ite()) += source(ite())[k] * partial(ite());
			} while (ite.next());
		}
	}
	LFHTEMP
	void DataGrid<C, nbdim>::ExpectedDistanceToOther_instate_I(const DataGrid< double ,nbdim> &source,const DataGrid<double,2>&  transition, int state){
		setSizes(source.dims);
		double max =0.0f;
		double tmp;
		int i,j,k,l;
		int nbstates = transition.dims[0];

		for(i=0;i<nbdim;i++) max += ((double)dims[i])* dims[i];
		max = sqrt(max);

		Vector< KeyElem< double, Tuple<unsigned int, nbdim> > > list;
		KeyIterator ite = KeyIterator(*this);
		if (ite.first()) do{
			list.push_back(KeyElem<double, Tuple<unsigned int, 2> >( ((double*)source(ite()))[state], ite()));
			(*this)(ite()) = max;
		} while (ite.next());

		list.sort();

		double transit[nbstates][1+ nbdim];

		for(i=0;i<nbdim;i++) {
			transit[i][0] = transition.data[state + j * nbstates]; // == 1.0f - exp(alpha)
			transit[i][nbdim] = log(1.0f - transition.data[state + j * nbstates]); // == 1.0f - exp(alpha)
			for(j=1;j<nbdim;j++) {
				transit[i][j] = 1.0f - exp(sqrt(1.0f + j) * transit[i][nbdim]); // exp(- alpha)
			}
			transit[i][nbdim] = 1.0f - exp(sqrt(5.0f) * transit[i][nbdim]);
		}

		for(j=0;j<nbstates;j++) {
			transit[state][j] = 1.0f;
			for(i=0;i<nbdim;i++) {
				if (i != state) transit[state][j] -= transit[i][j];
			}
		}


		Tuple<Tuple<unsigned int, 2>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
		Tuple<Tuple<unsigned int, 2>, nbdim * (nbdim -1) * 4> ncoor_knight;

		for(i=0;i<list.getSize();i++){
			int nbind = (*this).get_indirectNeightbor(list[i].d, ncoor);
			int nbkni = (*this).get_knightNeightbor(list[i].d, ncoor_knight);

			for(j=0;j<nbind;j++){
				l =0;
				for(k<0;k<nbdim;k++) if (ncoor[j][k] != list[i].d[k]) l++;
				max = 0.0f;
				tmp = source(list[i].d)[state] * transit[state][l-1];
				if (tmp != 0.0f){
					for(k<0;k<nbstates;k++) if (k != state) max += source(ncoor[j])[k] * transit[k][l-1];
					max = (*this)(list[i].d) / ( 1.0f + (max / tmp));
					if ( (*this)(ncoor[j]) > max)  (*this)(ncoor[j]) = max;
				}
			}

			for(j=0;j<nbkni;j++){
				tmp = source(list[i].d)[state] * transit[state][nbdim];
				if (tmp != 0.0f){
					for(k<0;k<nbstates;k++) if (k != state) max += source(ncoor_knight[j])[k] * transit[k][nbdim];
					max = (*this)(list[i].d) / ( 1.0f + (max / tmp));
					if ( (*this)(ncoor_knight[j]) > max)  (*this)(ncoor[j]) = max;
				}
			}



		}

	}

	LFHTEMP
	Vector< KeyElem<C, Tuple<unsigned int, nbdim> > > DataGrid<C, nbdim>::get_ordered_Data_Coors_Tagged() const{
		Vector< KeyElem< C , Tuple<unsigned int, nbdim> > > _out;
		_out.setSize(this->totsize());
		unsigned int i;

		typename DataGrid<C,nbdim>::KeyIterator ite = getKeyIterator();
		i=0;
		if (ite.first()) do{
			_out[i] = KeyElem< C , Tuple<unsigned int, nbdim> >( (*this)(ite()) , ite());
			i++;
		} while (ite.next());
		_out.sort();

		return(_out);
	}

	LFHTEMP
	Vector< Tuple<unsigned int, nbdim> > DataGrid<C, nbdim>::get_ordered_Data_Coors() const{
		Vector< Tuple<unsigned int, nbdim> > _out;
		Vector< KeyElem< C , Tuple<unsigned int, nbdim> > > _tmp = (*this).get_ordered_Data_Coors_Tagged();
		unsigned int ts = (*this).totsize();
		unsigned int i;

		for(i=0;i<ts;i++) _out[i] = _tmp[i].d;

		return(_out);
	}


	LFHTEMP
	DataGrid< Tuple<unsigned int, nbdim> ,nbdim> DataGrid<C, nbdim>::LocalMaxMap() const{
		// for each point, find its physical local maxima and minima, than compute the excess physical distance of moving through current point.

		DataGrid< Tuple<unsigned int, nbdim> , nbdim> fout;

		fout.setSizes(dims);

		Vector< KeyElem< C ,Tuple<unsigned int, nbdim> > > ordered_coor = (*this).get_ordered_Data_Coors_Tagged();

		unsigned int i,j,s;
		Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
		C extr;
		for(i=ordered_coor.getSize()-1;i != 0xFFFFFFFF ;i--){
			extr = ordered_coor[i].k;
			fout(ordered_coor[i].d) = ordered_coor[i].d;
			int nbind = (*this).get_indirectNeightbor(ordered_coor[i].d, ncoor);
			extr = ordered_coor[i].k;

			for(j=0;j<nbind;j++){
				if (extr < (*this)(ncoor[j])) {extr = (*this)(ncoor[j]); fout(ordered_coor[i].d) = fout(ncoor[j]);}
			}
		}

		return(fout);

	}

		LFHTEMP
		DataGrid< Tuple<unsigned int, nbdim> ,nbdim> DataGrid<C, nbdim>::LocalMinMap() const{
			// for each point, find its physical local maxima and minima, than compute the excess physical distance of moving through current point.

			DataGrid< Tuple<unsigned int, nbdim> , nbdim> fout;

			fout.setSizes(dims);

			Vector< KeyElem< C ,Tuple<unsigned int, nbdim> > > ordered_coor = (*this).get_ordered_Data_Coors_Tagged();

			unsigned int i,j,s;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
			C extr;
			for(i=0;i<ordered_coor.getSize();i++){
				extr = ordered_coor[i].k;
				fout(ordered_coor[i].d) = ordered_coor[i].d;
				int nbind = (*this).get_indirectNeightbor(ordered_coor[i].d, ncoor);
				extr = ordered_coor[i].k;

				for(j=0;j<nbind;j++){
					if (extr > (*this)(ncoor[j])) {extr = (*this)(ncoor[j]); fout(ordered_coor[i].d) = fout(ncoor[j]);}
				}
			}

			return(fout);

		}


		LFHTEMP
		DataGrid< unsigned int ,nbdim> DataGrid<C, nbdim>::SegementClimbToMax(unsigned int &nbsegs, bool dist_penal, bool knight, bool equal_join, bool mark_maxima, const DataGrid< bool ,nbdim>* ignore_filter) const{
			DataGrid< unsigned int , nbdim> fout;

			fout.setSizes(dims);

			Vector< KeyElem< C ,Tuple<unsigned int, nbdim> > > ordered_coor = (*this).get_ordered_Data_Coors_Tagged();

			unsigned int i,j,k,l,s;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > mcoor;
			Tuple< Tuple<unsigned int , nbdim>,  nbdim * (nbdim -1) * 4   > ncoork;
			Tuple< Tuple<unsigned int , nbdim>,  nbdim * (nbdim -1) * 4   > mcoork;
			C extr;
			C curv;

		//	printf("%c in on!\n", mark_maxima ? 'Y' : 'n' );

			nbsegs = 0;
			stack< Tuple<unsigned int, nbdim> > join;

			if (equal_join) ExOp::toZero(fout);

			for(i=ordered_coor.getSize()-1;i != 0xFFFFFFFF ;i--){

				if ((ignore_filter != NULL)&&(!((*ignore_filter)(ordered_coor[i].d)))){ // ILLEGAL PIXEL
					nbsegs++;
					if (mark_maxima) fout(ordered_coor[i].d) = 0x80000000 | nbsegs;
					else fout(ordered_coor[i].d) = nbsegs;
					continue;
				}

				s = 0;
				curv = (*this)(ordered_coor[i].d);
				int nbind = (*this).get_indirectNeightbor(ordered_coor[i].d, ncoor);

				extr = (dist_penal)  ? ordered_coor[i].k - curv : ordered_coor[i].k ;
				if (dist_penal){
					for(j=0;j<nbind;j++){
						l=0;
						for(k=0;k<nbdim;k++) if (ncoor[j][k] != ordered_coor[i].d[k]) l++;
						if ((extr < (((*this)(ncoor[j]) - curv) * pow((double) l,-0.5f)) )&&((ignore_filter==NULL)||((*ignore_filter)(ncoor[j])))) {s = fout(ncoor[j]);extr = (((*this)(ncoor[j]) - curv) * pow((double) l,-0.5f));}
					}
				}else{
					for(j=0;j<nbind;j++){
						if ((extr < (*this)(ncoor[j]))&&((ignore_filter==NULL)||((*ignore_filter)(ncoor[j])))) {s = fout(ncoor[j]);extr = (*this)(ncoor[j]); }
					}
				}
				if (knight){
					nbind = (*this).get_knightNeightbor(ordered_coor[i].d, ncoork);
					if (dist_penal){
					for(j=0;j<nbind;j++){
						if ((extr < (((*this)(ncoork[j]) - curv) * pow(5.0f,-0.5f)) )&&((ignore_filter==NULL)||((*ignore_filter)(ncoork[j])))) {s = fout(ncoork[j]);extr = (((*this)(ncoork[j]) - curv) * pow(5.0f,-0.5f)); }
					}
					}else{
						for(j=0;j<nbind;j++){
							if ((extr < (*this)(ncoork[j]))&&((ignore_filter==NULL)||((*ignore_filter)(ncoork[j])))) {s = fout(ncoork[j]);extr = (*this)(ncoork[j]); }
						}
					}
				}
				extr = ordered_coor[i].k;

				if (s == 0) {
					if (equal_join){
						if (fout(ordered_coor[i].d) == 0) {

						join.push(ordered_coor[i].d);
						nbsegs++;
						if (mark_maxima) fout(ordered_coor[i].d) = 0x80000000 | nbsegs;
						else fout(ordered_coor[i].d) = nbsegs;

						while(!join.empty()){
							int nbind = (*this).get_indirectNeightbor(join.top(), ncoor);
							for(j=0;j<nbind;j++){
								if ((extr == (*this)(ncoor[j]))&&(fout(ncoor[j]) == 0)) {

									int nbindm = (*this).get_indirectNeightbor(ncoor[j], mcoor);
									for(k=0;k<nbindm;k++){
										if (extr < (*this)(mcoor[k])) break;
									}

									if (k == nbindm) {fout(ncoor[j]) = nbsegs; join.push(ncoor[j]);}
								}
							}
							join.pop();
						}
						}
					}else{
						nbsegs++;
						if (mark_maxima) fout(ordered_coor[i].d) = 0x80000000 | nbsegs;
						else fout(ordered_coor[i].d) = nbsegs;
					}

				} else fout(ordered_coor[i].d) = s & 0x7FFFFFFF;
			}
			return(fout);
		}


		LFHTEMP
		DataGrid< unsigned int ,nbdim> DataGrid<C, nbdim>::SegementClimbToMin(unsigned int &nbsegs, bool dist_penal, bool knight, bool equal_join, bool mark_maxima, const DataGrid< bool ,nbdim>* ignore_filter) const{
			DataGrid< unsigned int , nbdim> fout;

			fout.setSizes(dims);

			Vector< KeyElem< C ,Tuple<unsigned int, nbdim> > > ordered_coor = (*this).get_ordered_Data_Coors_Tagged();

			unsigned int i,j,k,l,s;
			Tuple<Tuple<unsigned int, nbdim>,(unsigned int) (TEMPLATE_INT_POWER<3,nbdim>::ans) -1u > ncoor;
			Tuple<Tuple<unsigned int, nbdim>,(unsigned int) (TEMPLATE_INT_POWER<3,nbdim>::ans) -1u > mcoor;
			Tuple< Tuple<unsigned int , nbdim>,  nbdim * (nbdim -1) * 4   > ncoork;
			Tuple< Tuple<unsigned int , nbdim>,  nbdim * (nbdim -1) * 4   > mcoork;
			C extr;
			C curv;

			nbsegs = 0;
			stack< Tuple<unsigned int, nbdim> > join;

			if (equal_join) ExOp::toZero(fout);

			for(i=0 ;i < ordered_coor.getSize() ;i++){

				if ((ignore_filter != NULL)&&(!((*ignore_filter)(ordered_coor[i].d)))){ // ILLEGAL PIXEL
					nbsegs++;
					if (mark_maxima) fout(ordered_coor[i].d) = 0x80000000 | nbsegs;
					else fout(ordered_coor[i].d) = nbsegs;
					continue;
				}
				s = 0;
				curv = (*this)(ordered_coor[i].d);
				unsigned int nbind = this->get_indirectNeightbor(ordered_coor[i].d, ncoor);

				extr = (dist_penal)  ? ordered_coor[i].k - curv : ordered_coor[i].k ;
				if (dist_penal){
					for(j=0;j<nbind;j++){
						l=0;
						for(k=0;k<nbdim;k++) if (ncoor[j][k] != ordered_coor[i].d[k]) l++;
						if ((extr > (((*this)(ncoor[j]) - curv) * pow((double) l,-0.5f)))&&((ignore_filter==NULL)||((*ignore_filter)(ncoor[j])))) {s = fout(ncoor[j]);extr = (((*this)(ncoor[j]) - curv) * pow((double) l,-0.5f)); }
					}
				}else{
					for(j=0;j<nbind;j++){
						if ((extr > (*this)(ncoor[j]))&&((ignore_filter==NULL)||((*ignore_filter)(ncoor[j])))) {s = fout(ncoor[j]);extr = (*this)(ncoor[j]); }
					}
				}
				if (knight){
					nbind = (*this).get_knightNeightbor(ordered_coor[i].d, ncoork);
					if (dist_penal){
					for(j=0;j<nbind;j++){
						if ((extr > (((*this)(ncoork[j]) - curv) * pow(5.0f,-0.5f)) )&&((ignore_filter==NULL)||((*ignore_filter)(ncoork[j])))) {s = fout(ncoork[j]);extr = (((*this)(ncoork[j]) - curv) * pow(5.0f,-0.5f)); }
					}
					}else{
						for(j=0;j<nbind;j++){
							if ((extr > (*this)(ncoork[j]))&&((ignore_filter==NULL)||((*ignore_filter)(ncoork[j])))) {s = fout(ncoork[j]);extr = (*this)(ncoork[j]); }
						}
					}
				}
				extr = ordered_coor[i].k;

				if (s == 0) {
					if (equal_join){
						if (fout(ordered_coor[i].d) == 0) {

						join.push(ordered_coor[i].d);
						nbsegs++;
						if (mark_maxima) fout(ordered_coor[i].d) = 0x80000000 | nbsegs;
						else fout(ordered_coor[i].d) = nbsegs;

						while(!join.empty()){
							unsigned int nbind = (*this).get_indirectNeightbor(join.top(), ncoor);
							for(j=0;j<nbind;j++){
								if ((extr == (*this)(ncoor[j]))&&(fout(ncoor[j]) == 0)) {

									unsigned int nbindm = (*this).get_indirectNeightbor(ncoor[j], mcoor);
									for(k=0;k<nbindm;k++){
										if (extr > (*this)(mcoor[k])) break;
									}

									if (k == nbindm) {fout(ncoor[j]) = nbsegs; join.push(ncoor[j]);}
								}
							}
							join.pop();
						}
						}
					}else{
						nbsegs++;
						if (mark_maxima) fout(ordered_coor[i].d) = 0x80000000 | nbsegs;
						else fout(ordered_coor[i].d) = nbsegs;
					}

				} else fout(ordered_coor[i].d) = s & 0x7FFFFFFF;
			}
			return(fout);
		}

		LFHTEMP
		DataGrid< unsigned int ,nbdim> DataGrid<C, nbdim>::SegementClimbToMaxMax(unsigned int &nbsegs, bool equal_join) const{
			DataGrid< unsigned int , nbdim> fout;

			fout.setSizes(dims);
			Vector< C> maxes;

			Vector< KeyElem< C ,Tuple<unsigned int, nbdim> > > ordered_coor = (*this).get_ordered_Data_Coors_Tagged();

			unsigned int i,j,k,s;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > mcoor;
			C extr;
			nbsegs = 0;
			stack< Tuple<unsigned int, nbdim> > join;

			if (equal_join) ExOp::toZero(fout);

			for(i=ordered_coor.getSize()-1;i != 0xFFFFFFFF ;i--){

				s = 0;
				int nbind = (*this).get_indirectNeightbor(ordered_coor[i].d, ncoor);
				extr = ordered_coor[i].k;

				for(j=0;j<nbind;j++){
					if (extr < (*this)(ncoor[j])) {s = fout(ncoor[j]);extr = (*this)(ncoor[j]); }
				}

				if (s == 0) { // current point is a local maxima!


					if (equal_join){
						if (fout(ordered_coor[i].d) == 0) {

							join.push(ordered_coor[i].d);
							nbsegs++;
							maxes.push_back(extr);
							//areas.push_back(KeyElem<C, unsigned int> >(extr,nbsegs));

							fout(ordered_coor[i].d) = nbsegs;

							while(!join.empty()){
								unsigned int nbind = (*this).get_indirectNeightbor(join.top(), ncoor);
								for(j=0;j<nbind;j++){
									if ((extr == (*this)(ncoor[j]))&&(fout(ncoor[j]) == 0)) {

										int nbindm = (*this).get_indirectNeightbor(ncoor[j], mcoor);
										for(k=0;k<nbindm;k++){
											if (extr < (*this)(mcoor[k])) break;
										}
										if (k == nbindm) {fout(ncoor[j]) = nbsegs; join.push(ncoor[j]);}
									}
								}
								join.pop();
							}
						}
					}else{
						nbsegs++;

						fout(ordered_coor[i].d) = nbsegs;
						maxes.push_back(extr);
						//areas.push_back(KeyElem<C, unsigned int> >(extr,nbsegs));
					}

				} else fout(ordered_coor[i].d) = s;
			}


			KeyElem< C , unsigned int> *transit = new KeyElem< C , unsigned int>[nbsegs];
			for(j=0;j<nbsegs;j++) {ExOp::toMin(transit[j].k); transit[j].d = j;}


			for(i=ordered_coor.getSize()-1;i != 0xFFFFFFFF ;i--){
				s = 0;
				int nbind = (*this).get_indirectNeightbor(ordered_coor[i].d, ncoor);
				extr = ordered_coor[i].k;
				s = fout(ordered_coor[i].d)-1;
				if (s != 0xFFFFFFFF){
					for(j=0;j<nbind;j++){
						k = fout(ncoor[j])-1;
						if (k != 0xFFFFFFFF){
						if ((maxes[k] != maxes[s])&&( extr < (*this)(ncoor[j]) )) {
							if (maxes[s] > maxes[k]) {if (maxes[s] + extr > transit[k].k) {transit[k].k = maxes[s] + extr; transit[k].d = s;}
							}else{ if (maxes[k] + extr > transit[s].k) {transit[s].k = maxes[k] + extr; transit[s].d = k;}}
							}
						}
					}
				}
			}


			s = 0;
			for(j=0;j<nbsegs;j++) {
				while((transit[ transit[j].d ].d != transit[j].d)&&(transit[j].d >= s)){
					transit[j].d = transit[ transit[j].d ].d;
				}

				if (transit[ j ].d >= s){
					if (transit[j].d != j) transit[ transit[j].d ].d = s;
					transit[ j ].d = s;
					s++;
				}
			}


			typename DataGrid< unsigned int , nbdim>::KeyIterator ite = fout.getKeyIterator();

			if (ite.first()) do{
				s = fout(ite());
				if (s !=0) fout(ite()) = 1+ transit[s-1].d;
			} while (ite.next());



			return(fout);
		}



		/*
		LFHTEMP template<unsigned int subdim >
		const DataGrid< unsigned int,nbdim - subdim>& DataGrid<C, nbdim>::SegementClimbToMax(unsigned int &nbsegs, DataGrid< unsigned int ,nbdim - subdim> &fout , bool eq_join) const{


			fout.setSizes(dims + subdim);

			Vector< KeyElem< C ,Tuple<unsigned int, nbdim> > > ordered_coor = (*this).get_ordered_Data_Coors_Tagged();

			unsigned int i,j,k,s;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > ncoor;
			Tuple<Tuple<unsigned int, nbdim>,(TEMPLATE_INT_POWER<3,nbdim>::ans) -1 > mcoor;
			C extr;
			nbsegs = 0;
			stack< Tuple<unsigned int, nbdim> > join;

			if (equal_join) ExOp::toZero(fout);

			for(i=ordered_coor.getSize()-1;i != 0xFFFFFFFF ;i--){

				s = 0;
				int nbind = (*this).get_indirectNeightbor(ordered_coor[i].d, ncoor);
				extr = ordered_coor[i].k;

				for(j=0;j<nbind;j++){
					if (extr < (*this)(ncoor[j])) {s = fout(ncoor[j]);extr = (*this)(ncoor[j]); }
				}

				if (s == 0) {
					if (equal_join){
						if (fout(ordered_coor[i].d) == 0) {

							join.push(ordered_coor[i].d);
							nbsegs++;
							fout(ordered_coor[i].d) = nbsegs;

							while(!join.empty()){
								int nbind = (*this).get_indirectNeightbor(join.top(), ncoor);
								for(j=0;j<nbind;j++){
									if ((extr == (*this)(ncoor[j]))&&(fout(ncoor[j]) == 0)) {

										int nbindm = (*this).get_indirectNeightbor(ncoor[j], mcoor);
										for(k=0;k<nbindm;k++){
											if (extr < (*this)(mcoor[k])) break;
										}
										if (k == nbindm) {fout(ncoor[j]) = nbsegs; join.push(ncoor[j]);}
									}
								}
								join.pop();
							}
						}
					}else{
						nbsegs++;
						fout(ordered_coor[i].d) = nbsegs;
					}

				} else fout(ordered_coor[i].d) = s;
			}




			return(fout);
		}
		LFHTEMP template<unsigned int subdim >
		const DataGrid< unsigned int,nbdim - subdim>& DataGrid<C, nbdim>::SegementClimbToMin(unsigned int &nbsegs, DataGrid< unsigned int ,nbdim - subdim> &fout , bool eq_join) const{
			return(fout);
		}*/



/*
	LFHTEMP
	const DataGrid<C,nbdim> & DataGrid<C, nbdim>::FFtransform(int dim, unsigned int size, unsigned int min){
		DataGrid<C,nbdim> out;
		unsigned int tmpdim[nbdim]; getDims(tmpdim);
		tmpdim[dim] = size*2;
		out.setSizes(tmpdim);

		unsigned int coors[nbdim];
		unsigned int altcoor =0;
		memset(coors,'\0',sizeof(unsigned int)*nbdim);
		int FFsize =1;
		int HHsize;
		int IIcur;
		int JJcur;
		while(FFsize < size*2) FFsize <<=1;
		C* buffer = new C[FFsize];
		Tuple<C,8> readbuffer;
		PolyThing<C> polynomial;
		int altdim =0;
		double offset;
		unsigned int tmp;
		int i;
		while(true){

			memset(buffer,'\0',sizeof(C)*FFsize);
			FFsize;
			IIcur =0;
			for(coors[dim] =0;coors[dim]<tmpdim[dim];coors[dim]++){
				buffer[IIcur] = ((*this)(coors));
				for(JJcur = FFsize;(JJcur & IIcur)&& (JJcur>2); JJcur>>=1) IIcur ^= JJcur;
				IIcur |= JJcur;
			}

			for(HHsize =2; HHsize <= FFsize; HHsize <<=1){
				for(IIcur=0;IIcur < FFsize;IIcur += (HHsize << 1)){
				//	for(i = 0;i< HHsize;i++)
				//		tmp = buffer[IIcur | i] + buffer[IIcur | i | HHsize];
				//		buffer[IIcur | i] = buffer[IIcur | i] + buffer[IIcur | i | HHsize];

				}
			}



			for(coors[dim]=0;coors[dim]<size*2;coors[dim]++) out(coors) = buffer[coors[dim]];

			for(altdim=0;altdim<nbdim;altdim++) {
				if (altdim != dim) {
					if (coors[altdim] == tmpdim[altdim]-1) coors[altdim] =0;
					else {coors[altdim]++; break;}
				}
			}
			if (altdim == nbdim) break;
		}


		return(out);
	}*/

	LFHTEMP
	void DataGrid<C, nbdim>::exponentialBlur(double std_dev, bool circular){
		unsigned int coor[nbdim];
		int dir,altdir;
		int i;
		C* buffer;
		C current;

		int cap;
		double decay = 1.0f+1.0f/(std_dev*std_dev);
		decay = decay - sqrt(decay*decay - 1.0f);
		double peak_w = (1-decay)/(1+decay);

	//	printf("%f\t%f\n",decay,peak_w);
		// peak weight is: (1-decay)/(1+decay), so that total weight is 1
		for(dir =0;dir<nbdim;dir++){
			memset(coor, '\0',sizeof(int)*nbdim);
			buffer = new C[getDimSize(dir)];
			while(true){
				//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
				// prepare initial for backward!
				coor[dir] = getDimSize(dir)-1;
				if (circular){
					current = (*this)(coor);
					for((coor[dir])--;coor[dir] != 0xFFFFFFFF;(coor[dir])--) current = (current * decay) + (*this)(coor);
					coor[dir] = getDimSize(dir)-1;
					buffer[coor[dir]] = (current * decay) + (*this)(coor);
				} else buffer[coor[dir]] = (*this)(coor);
				//				printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);


				for(coor[dir]--;coor[dir]!= 0xFFFFFFFF;coor[dir]--) buffer[coor[dir]] = (buffer[coor[dir]+1]  * decay) + (*this)(coor);

				//			printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
				// prepared initial for forward!
				cap = getDimSize(dir);
				coor[dir] = 0;
				current = (*this)(coor);
				if (circular){
					for(coor[dir]++;coor[dir]<cap;coor[dir]++) current = (current * decay) + (*this)(coor);
					coor[dir] = 0;
					current = (current * decay) + (*this)(coor);
				}

				(*this)(coor) = (current  + buffer[coor[dir]] * decay)* peak_w;
				for(coor[dir]++;coor[dir]<cap;coor[dir]++) {
//					current = (current * decay) + (*this)(coor);
//					(*this)(coor) = (current + (buffer[coor[dir]] * decay)) * peak_w;
					current = (current * decay) + (*this)(coor);
					buffer[coor[dir]] *= decay;
					(*this)(coor) = (current + buffer[coor[dir]] ) * peak_w;
				}


				// iterate througth other dimentions!
				for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == getDimSize((dir + altdir) % nbdim)-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
				if (altdir == nbdim) break;
			}
			delete[](buffer);
		}
	}

	LFHTEMP
	template <unsigned int size, Tuple_flag Cflag> DataGrid<Tuple<C, size,Cflag>, nbdim> DataGrid<C, nbdim>::exponentialBlur(Tuple<double,size,Cflag> apertures, bool circular){
		unsigned int coor[nbdim];
		int dir,altdir;
		int i;
		Tuple<C,size,Cflag>* buffer;
		Tuple<C,size,Cflag> current;
		DataGrid<Tuple<C,size,Cflag>, nbdim> _out;
		_out.setSizes(dims);


		Tuple<double,size> decay;
		Tuple<double,size> peak;
		double tmp;
		for(i=0;i<size;i++){
			tmp = 1.0f+1.0f/(apertures[i]*apertures[i]);
			decay[i] = tmp - sqrt(tmp*tmp - 1.0f);
			peak[i] = (1-decay[i])/(1+decay[i]);
		}
		int cap;
		// peak weight is: (1-decay)/(1+decay), so that total weight is 1
		for(dir =0;dir<nbdim;dir++){
			memset(coor, '\0',sizeof(int)*nbdim);
			buffer = new Tuple<C,size,Cflag>[getDimSize(dir)];
			while(true){
				//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
				// prepare initial for backward!
				coor[dir] = getDimSize(dir)-1;
				if (circular){
					current = ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					for((coor[dir])--;coor[dir] != 0xFFFFFFFF;(coor[dir])--) current = (current * decay) + ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					coor[dir] = getDimSize(dir)-1;
					buffer[coor[dir]] = (current * decay) + ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
				} else buffer[coor[dir]] = ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
				//				printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);


				for(coor[dir]--;coor[dir]!= 0xFFFFFFFF;coor[dir]--) buffer[coor[dir]] = (buffer[coor[dir]+1]  * decay) + ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));

				//			printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
				// prepared initial for forward!
				cap = getDimSize(dir);
				coor[dir] = 0;
				current = ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
				if (circular){
					for(coor[dir]++;coor[dir]<cap;coor[dir]++) current = (current * decay) + ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					coor[dir] = 0;
					current = (current * decay) + ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
				}

				_out(coor) = (current * peak) + (buffer[coor[dir]] * (peak * decay));
				for(coor[dir]++;coor[dir]<cap;coor[dir]++) {
					current = (current * decay) + ((dir == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					_out(coor) = (current + buffer[coor[dir]] * decay) * peak;
				}


				// iterate througth other dimentions!
				for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == getDimSize((dir + altdir) % nbdim)-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
				if (altdir == nbdim) break;
			}
			delete[](buffer);
		}
		return(_out);
	}

	LFHTEMP
	template <unsigned int size, Tuple_flag Cflag> DataGrid<Tuple<C, size,Cflag>, nbdim> DataGrid<C, nbdim>::crudeGaussianBlur(Tuple<double,size,Cflag> apertures, bool circular){
		unsigned int coor[nbdim];
		int dir,altdir;
		int i;
		Tuple<C,size,Cflag>* buffer;
		Tuple<C,size,Cflag> current;
		DataGrid<Tuple<C,size,Cflag>, nbdim> _out;
		_out.setSizes(dims);

		int loop;
		Tuple<double,size> decay;
		Tuple<double,size> peak;
		Tuple<double,size> boundm;
		double tmp;
		for(i=0;i<size;i++){
			tmp = 1.0f+2.0f/(apertures[i]*apertures[i]);
			decay[i] = tmp - sqrt(tmp*tmp - 1.0f);
//			decay[i] =0.0f;
			peak[i] = (1-decay[i])/(1+decay[i]);
			if (!circular) boundm[i] = 1.0f / (1 - decay[i]);
		}
		int cap;
		// peak weight is: (1-decay)/(1+decay), so that total weight is 1
		for(loop =0;loop < 2; loop++){
			for(dir =0;dir<nbdim;dir++){
				memset(coor, '\0',sizeof(int)*nbdim);
				buffer = new Tuple<C,size,Cflag>[getDimSize(dir)];
				while(true){
					//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
					// prepare initial for backward!
					coor[dir] = getDimSize(dir)-1;
					if (circular){
						current = ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
						for((coor[dir])--;coor[dir] != 0xFFFFFFFF;(coor[dir])--) current = (current * decay) + ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
						coor[dir] = getDimSize(dir)-1;
						buffer[coor[dir]] = (current * decay) + ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					} else buffer[coor[dir]] = ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor)) * boundm;  // / ((decay -1) * -1);
					//				printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);


					for(coor[dir]--;coor[dir]!= 0xFFFFFFFF;coor[dir]--) buffer[coor[dir]] = (buffer[coor[dir]+1]  * decay) + ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));

					//			printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
					// prepared initial for forward!
					cap = getDimSize(dir);
					coor[dir] = 0;
					current = ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					if (circular){
						for(coor[dir]++;coor[dir]<cap;coor[dir]++) current = (current * decay) + ((dir +loop == 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
						coor[dir] = 0;
						current = (current * decay) + ((dir + loop== 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
					}  else current *= boundm;

					buffer[coor[dir]] *= decay;
					buffer[coor[dir]] += current;
					_out(coor) = buffer[coor[dir]] * peak;
					for(coor[dir]++;coor[dir]<cap;coor[dir]++) {
						current *= decay;
						current += ((dir +loop== 0) ? Tuple<C,size,Cflag>((*this)(coor)) : _out(coor));
						(buffer[coor[dir]]) *= decay;
						(buffer[coor[dir]]) += current;
						(buffer[coor[dir]]) *= peak;
						_out(coor) = buffer[coor[dir]];
						//_out(coor) = Tuple<C,size,Cflag>((*this)(coor));
					}


					// iterate througth other dimentions!
					for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == getDimSize((dir + altdir) % nbdim)-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
					if (altdir == nbdim) break;
				}
				delete[](buffer);
			}
		}
		return(_out);
	}

		LFHTEMP
		 DataGrid<C, nbdim> DataGrid<C, nbdim>::crudeGaussianBlur(double apertures, bool circular){
		unsigned int coor[nbdim];
		unsigned int dir,altdir;

		C* buffer;
		C current;
		DataGrid<C, nbdim> _out;
		_out.setSizes(dims);

		unsigned int loop;
		double decay;
		double peak;
		double boundm;
		double tmp;
			tmp = 1.0f+2.0f/(apertures*apertures);
			decay = tmp - sqrt(tmp*tmp - 1.0f);
//			decay[i] =0.0f;
			peak = (1-decay)/(1+decay);
			if (!circular) boundm = 1.0f / (1 - decay);

		unsigned int cap;
		// peak weight is: (1-decay)/(1+decay), so that total weight is 1
		for(loop =0;loop < 2; loop++){
			for(dir =0;dir<nbdim;dir++){
				memset(coor, '\0',sizeof(int)*nbdim);
				buffer = new C[getDimSize(dir)];
				while(true){
					//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
					// prepare initial for backward!
					coor[dir] = getDimSize(dir)-1;
					if (circular){
						current = ((dir +loop == 0) ? (*this)(coor) : _out(coor));
						for((coor[dir])--;coor[dir] != 0xFFFFFFFF;(coor[dir])--) { current *= decay; current += ((dir +loop == 0) ? (*this)(coor) : _out(coor));}
						coor[dir] = getDimSize(dir)-1;
						buffer[coor[dir]] = (current * decay) + ((dir +loop == 0) ? (*this)(coor) : _out(coor));
					} else buffer[coor[dir]] = ((dir +loop == 0) ? (*this)(coor) : _out(coor)) * boundm;  // / ((decay -1) * -1);
					//				printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);


					for(coor[dir]--;coor[dir]!= 0xFFFFFFFF;coor[dir]--) buffer[coor[dir]] = (buffer[coor[dir]+1]  * decay) + ((dir +loop == 0) ? (*this)(coor) : _out(coor));

					//			printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
					// prepared initial for forward!
					cap = getDimSize(dir);
					coor[dir] = 0;
					current = ((dir +loop == 0) ? (*this)(coor) : _out(coor));
					if (circular){
						for(coor[dir]++;coor[dir]<cap;coor[dir]++) current = (current * decay) + ((dir +loop == 0) ? (*this)(coor) : _out(coor));
						coor[dir] = 0;
						current = (current * decay) + ((dir + loop== 0) ? (*this)(coor) : _out(coor));
					}  else current *= boundm;

					buffer[coor[dir]] *= decay;
					buffer[coor[dir]] += current;
					_out(coor) = buffer[coor[dir]] * peak;
					for(coor[dir]++;coor[dir]<cap;coor[dir]++) {
						current *= decay;
						current += ((dir +loop== 0) ? (*this)(coor) : _out(coor));
						(buffer[coor[dir]]) *= decay;
						(buffer[coor[dir]]) += current;
						(buffer[coor[dir]]) *= peak;
						_out(coor) = buffer[coor[dir]];
						//_out(coor) = Tuple<C,size,Cflag>((*this)(coor));
					}


					// iterate througth other dimentions!
					for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == getDimSize((dir + altdir) % nbdim)-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
					if (altdir == nbdim) break;
				}
				delete[](buffer);
			}
		}
		return(_out);
	}


	LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::Crop(const Tuple<unsigned int, nbdim> &min, const Tuple<unsigned int, nbdim> &max) const{
		DataGrid<C,nbdim> _out;
		int i,dir;
		Tuple<unsigned int, nbdim> coor;
		for(i=0;i<nbdim;i++) coor[i] = max[i] - min[i] + 1;
		_out.setSizes(&(coor[0]));

		Tuple<unsigned int, nbdim> altcoor;
		coor = min;
		ExOp::toZero(altcoor);
		do{
			_out(altcoor) = (*this)(coor);
			for(dir=0;dir<nbdim;dir++) if (coor[dir] == max[dir]) {altcoor[dir] =0;coor[dir] =min[dir];} else {coor[dir]++;altcoor[dir]++; break;}
		}while(dir < nbdim);
		return(_out);

	}


	LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toblur_dir(double aperture,unsigned int d){

		Bluestein bl;
		Tuple<unsigned int, nbdim> coor;


		unsigned int maxb;

		bl.setSize(this->dims[d]);
		maxb = bl.getBufferSize();
		typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb];

		DataGrid< double ,1> dawin;
		dawin.setSizes(&(this->dims[d]));

		double tmpd, sumd;
		sumd = 0.0f;

		double fact =  (1.0f/ (aperture * sqrt(2.0f)));
		for(coor[0] =0; coor[0] < dawin.dims[0]/2;coor[0]++){
			tmpd = fact * coor[0];
			tmpd = exp(-tmpd*tmpd);
			sumd += tmpd;
			dawin.data[coor[0]]  = tmpd;
		}
		for(; coor[0] < dawin.dims[0];coor[0]++){
			tmpd = fact * (dawin.dims[0] - coor[0]);
			tmpd = exp(-tmpd*tmpd);
			sumd += tmpd;
			dawin.data[coor[0]] = tmpd;
		}


		KeyIterator ite = KeyIterator(*this);
		DataGrid< typename ExCo<double>::COMPLEX_TYPE , 1> dacon = dawin.FFtransform_2dir(0, bl, (typename  ExCo<double>::COMPLEX_TYPE * ) buffer);
		DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ft = this->FFtransform_2dir(d, bl, (typename  ExCo<double>::COMPLEX_TYPE * ) buffer);

		//printf("%i is size!\n", maxb); fflush(stdout);
		//delete[](buffer); return *this;

		if (ite.first()) do{
			ft(ite()) *= dacon.data[ite()[d]];
		}while(ite.next());
		DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ift = ft.invFFtransform_2dir(d,bl,buffer);
		if (ite.first()) do{
			(*this)(ite()) = ( ExOp::mkrealproj(ift(ite())) ) /sumd;
		}while(ite.next());
		delete[](buffer);
		return *this;
	}

LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toblur_zeroborder_dir(double aperture,unsigned int d){
	//	printf("hasss!\n"); fflush(stdout);
	Bluestein bl;
	Tuple<unsigned int, nbdim> coor;
	//	printf("ha!\n"); fflush(stdout);

	unsigned int maxb;

	bl.setSize(this->dims[d]*2);
	maxb = bl.getBufferSize();
	typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb];
	//	printf("ha!\n"); fflush(stdout);
	DataGrid< double ,1> dawin;

    Tuple<unsigned int,1u> onecoor;
	onecoor[0] = this->dims[d]*2;
	dawin.setSizes(onecoor);

	/*unsigned int onecoor;
	onecoor = this->dims[d]*2;
	dawin.setSizes(onecoor);*/

	double tmpd, sumd;
	sumd = 0.0f;
	//	printf("ha!\n"); fflush(stdout);
	double fact =  (1.0f/ (aperture * sqrt(2.0f)));
	for(coor[0] =0; coor[0] < dawin.dims[0]/2;coor[0]++){
		tmpd = fact * coor[0];
		tmpd = exp(-tmpd*tmpd);
		sumd += tmpd;
		dawin.data[coor[0]]  = tmpd;
	}
	for(; coor[0] < dawin.dims[0];coor[0]++){
		tmpd = fact * (dawin.dims[0] - coor[0]);
		tmpd = exp(-tmpd*tmpd);
		sumd += tmpd;
		dawin.data[coor[0]] = tmpd;
	}

	//printf("ha!\n"); fflush(stdout);
	DataGrid< typename ExCo<double>::COMPLEX_TYPE , 1> dacon = dawin.FFtransform_2dir(0, bl, (typename  ExCo<double>::COMPLEX_TYPE * ) buffer);
	//printf("ha!\n"); fflush(stdout);


	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ft = this->FFtransform_zeroborder_dir(d, bl, (typename  ExCo<double>::COMPLEX_TYPE * ) buffer);

	//printf("ha!\n"); fflush(stdout);

	typename DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim>::KeyIterator ite = typename DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim>::KeyIterator(ft);
	if (ite.first()) do{
		ft(ite()) *= dacon.data[ite()[d]];
	}while(ite.next());
	//printf("ha!\n"); fflush(stdout);
	DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ift = ft.invFFtransform_zeroborder_dir(d,bl,buffer);
	//printf("ha!\n"); fflush(stdout);


	KeyIterator ite2 = KeyIterator(*this);
	if (ite2.first()) do{
		(*this)(ite2()) = ( ExOp::mkrealproj(ift(ite2())) ) /sumd;
	}while(ite2.next());
	delete[](buffer);
	//printf("ha!\n"); fflush(stdout);

	return *this;
}

LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toInterleaved(unsigned int dir, unsigned int factor, unsigned int remainder){
    Tuple<unsigned int, nbdim> coor;
    unsigned int adir;
    for(adir=0;adir< nbdim;adir++) coor[adir] =dims[adir];
    coor[dir] *= factor;
    DataGrid<C, nbdim> tfout;
    tfout.setSizes(coor);
    ExOp::toZero(coor);
    Tuple<unsigned int, nbdim> coora = coor;
    do{
        for(coor[dir] = 0; coor[dir] < tfout.dims[dir] ;coor[dir]++){
            if ((coor[dir] % factor) == remainder) {coora[dir] = coor[dir] / factor; tfout(coor) = (*this)(coora);}
            else ExOp::toZero(tfout(coor));
            }
        for(adir=0;adir<nbdim-1;adir++)
            if (coor[(dir + adir) % nbdim] == dims[(dir + adir) % nbdim]) {coora[(dir + adir) % nbdim] =0;coor[(dir + adir) % nbdim] =0;}
            else {coora[(dir + adir) % nbdim]++;coor[(dir + adir) % nbdim]++; break;}
    }while(adir<nbdim-1);

    this->toMemmove(tfout);
    return *this;
}

LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toScaledWithConstantFrequencies(unsigned int dir, unsigned int factor){
    Bluestein bl[2];

    bl[0].setSize(this->dims[dir]);
    bl[1].setSize(this->dims[dir] * factor);
    unsigned int maxb = bl[1].getBufferSize();
    typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb];

	DataGrid< typename ExCo<C>::COMPLEX_TYPE,nbdim > dacon = this->FFtransform_2dir(dir, bl[0], (typename  ExCo<C>::COMPLEX_TYPE * ) buffer);
    dacon.toInterleaved(dir, factor, 0);
    DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ift = dacon.invFFtransform_2dir(dir,bl[1],buffer);
    this->setSizes(ift.dims);
	typename DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim>::KeyIterator ite = typename DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim>::KeyIterator(ift);
	if (ite.first()) do{
		(*this)(ite()) = ExOp::mkrealproj(ift(ite()));
	}while(ite.next());
    return *this;
}

	LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toresize(unsigned int* newdims) {
		DataGrid<C,nbdim> fout[2];
		char c = 2;
		unsigned int i;
		for(i=0;i<nbdim-1;i++){
			if (dims[i] != newdims[i]){
				fout[c & 1] = ((c == 2) ? (*this) : fout[(c ^ 1)]).resize_dir(i, newdims[i]) ;
				c = (c ^ 1) & 1;
			}
		}
		if (dims[i] != newdims[i]){
		   if (c != 2) (*this) = fout[(c ^ 1)].resize_dir(i, newdims[i]);
		   else {fout[0] = (*this).resize_dir(i, newdims[i]); (*this) = fout[0];}
		}else if (c != 2) {
			(*this) = fout[(c ^ 1)];
		}
		return(*this);
	}

	LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::resize(unsigned int* newdims) const{
		DataGrid<C,nbdim> fout[2];
		char c = 2;
		unsigned int i;
		for(i=0;i<nbdim;i++){
			if (dims[i] != newdims[i]){
				fout[c & 1] = ((c == 2) ? (*this) : fout[(c ^ 1)]).resize_dir(i, newdims[i]) ;
				c = (c ^ 1) & 1;
			}
		}
		return ((c == 2) ? (*this) : fout[(c ^ 1)]);
	}

	LFHTEMP DataGrid<C,nbdim>& DataGrid<C, nbdim>::toresize_crude(unsigned int* newdims) {
		DataGrid<C,nbdim> fout[2];
		char c = 2;
		unsigned int i;
		for(i=0;i<nbdim-1;i++){
			if (dims[i] != newdims[i]){
				fout[c & 1] = ((c == 2) ? (*this) : fout[(c ^ 1)]).resize_dir_crude(i, newdims[i]) ;
				c = (c ^ 1) & 1;
			}
		}
		if (dims[i] != newdims[i]){
			if (c != 2) (*this) = fout[(c ^ 1)].resize_dir_crude(i, newdims[i]);
			else {fout[0] = (*this).resize_dir_crude(i, newdims[i]); (*this) = fout[0];}
		}else if (c != 2) {
			(*this) = fout[(c ^ 1)];
		}
		return(*this);
	}

	LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::resize_crude(unsigned int* newdims) const{
		DataGrid<C,nbdim> fout[2];
		char c = 2;
		unsigned int i;
		for(i=0;i<nbdim;i++){
		if (dims[i] != newdims[i]){
				fout[c & 1] = ((c == 2) ? (*this) : fout[(c ^ 1)]).resize_dir_crude(i, newdims[i]) ;
				c = (c ^ 1) & 1;
			}
		}
		return ((c == 2) ? (*this) : fout[(c ^ 1)]);
	}

	LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::resize_dir(unsigned int cdir, unsigned int n_size) const{
		if (n_size == dims[cdir]) return(*this);
	    DataGrid<C,nbdim> fout;
		Bluestein bl_f(dims[cdir]);
		Bluestein bl_i(n_size);
		unsigned int dir = bl_i.getBufferSize();
		Tuple<unsigned int, nbdim> coor;
		unsigned int bsize = bl_f.getBufferSize();
		if (bsize < dir) bsize = dir;
		typename ExCo<C>::COMPLEX_TYPE * FFTbuffer = new typename ExCo<C>::COMPLEX_TYPE[bsize];
		typename ExCo<C>::COMPLEX_TYPE dazero; ExOp::toZero(dazero);
        for(dir=0;dir<nbdim;dir++) coor[dir] = dims[dir];
		coor[cdir] = n_size;
		fout.setSizes(coor);
		unsigned int j,k;
		double factor = ((double)n_size) / dims[cdir];
		if (dims[cdir] < n_size){
		ExOp::toZero(coor);
			do{
				for(coor[cdir] =0;coor[cdir]<dims[cdir];coor[cdir]++) FFTbuffer[coor[cdir]] = (*this)(coor);
                for(;coor[cdir]<bsize;coor[cdir]++) FFTbuffer[coor[cdir]] = dazero;


				bl_f.toFT(FFTbuffer);

				k = n_size - dims[cdir];

				for(j= n_size -1; j > (n_size / 2)+1; j--) FFTbuffer[j+k] = FFTbuffer[j];
				if (k & 1) FFTbuffer[j] *= 0.5f;
				FFTbuffer[j+k] = FFTbuffer[j];

				for(j+= k-1; j > n_size/2; j--) ExOp::toZero(FFTbuffer[j]);

	//			for(;j>0;j--) tmp_out.push_back(tmp_out[dims[cdir] - j ]); // extending array
	//			for(j= dims[cdir]-1-k; j > dims[cdir] / 2 ; j--) tmp_out[j + k] = tmp_out[j]; // copying remains
	//			for(j= dims[cdir] / 2; j < (dims[cdir] / 2)+k ; j++) ExOp::toZero(tmp_out[j]); // writing zeros

				bl_i.toIFT(FFTbuffer);

			if (k & 1) 	for(coor[cdir] =0;coor[cdir] < n_size;coor[cdir]++) fout(coor) = ExOp::mkrealproj(FFTbuffer[n_size - coor[cdir]- 1]) *factor;
			else		for(coor[cdir] =0;coor[cdir] < n_size;coor[cdir]++) fout(coor) = ExOp::mkrealproj(FFTbuffer[coor[cdir]]) *factor;
				for(dir=1;dir<nbdim;dir++) if (coor[(dir+cdir) % nbdim] == dims[(dir+cdir) % nbdim]-1) coor[(dir+cdir) % nbdim] =0; else {coor[(dir+cdir) % nbdim]++; break;}
			} while (dir<nbdim);
        }else{

		ExOp::toZero(coor);
		do{
			for(coor[cdir] =0;coor[cdir]<dims[cdir];coor[cdir]++) FFTbuffer[coor[cdir]] = (*this)(coor);
			for(;coor[cdir]<bsize;coor[cdir]++) FFTbuffer[coor[cdir]] = dazero;


			bl_f.toFT(FFTbuffer);

			// removing highest frequencies
			k = dims[cdir] - n_size;

			// if (k & 1) todo

			for(j= n_size / 2; j < n_size; j++) FFTbuffer[j] = FFTbuffer[j + k];

			bl_i.toIFT(FFTbuffer);

			if (k & 1) 	for(coor[cdir] =0;coor[cdir] < n_size;coor[cdir]++) fout(coor) = ExOp::mkrealproj(FFTbuffer[n_size - coor[cdir]- 1]) *factor;
			else		for(coor[cdir] =0;coor[cdir] < n_size;coor[cdir]++) fout(coor) = ExOp::mkrealproj(FFTbuffer[coor[cdir]]) *factor;

			for(dir=1;dir<nbdim;dir++) if (coor[(dir+cdir) % nbdim] == dims[(dir+cdir) % nbdim]-1) coor[(dir+cdir) % nbdim] =0; else {coor[(dir+cdir) % nbdim]++; break;}
		} while (dir<nbdim);
		}
		delete[](FFTbuffer);
		return(fout);
	}

LFHTEMP DataGrid<C,nbdim> DataGrid<C, nbdim>::resize_dir_crude(unsigned int cdir, unsigned int n_size) const{
	if (n_size == dims[cdir]) return(*this);
	//printf("%i -> %i resize\n", dims[cdir], n_size);
	DataGrid<C,nbdim> fout;
	Tuple<unsigned int, nbdim> coor;
	unsigned int dir;
	for(dir=0;dir<nbdim;dir++) coor[dir] = dims[dir];
	coor[cdir] = n_size;
	fout.setSizes(coor);
	unsigned int j;
	double factor = ((double)(dims[cdir] - 1)) / (n_size -1);
	double offset;
	C buffer[8];
	C coefs[8];
	C tmp[3];
	unsigned int old_coor;
	Tuple<unsigned int, nbdim> cooro;

	ExOp::toZero(coor); ExOp::toZero(cooro);
		do{
			coor[cdir] =0;
			buffer[3] = (*this)(coor);
			buffer[4] = buffer[3];
			buffer[5] = buffer[3];
			buffer[6] = buffer[3];
			buffer[7] = buffer[3];
			old_coor =0;cooro[cdir] =0;
			for(; cooro[cdir] <n_size;cooro[cdir]++){
				offset = factor * (0.5f + cooro[cdir]) + 4.0f;
				while( offset - 0.5f  > coor[cdir]) {
					buffer[ (coor[cdir] & 7) ] = (coor[cdir] >=  dims[cdir]) ? buffer[ ( (coor[cdir]+7) & 7) ]: (*this)(coor);
					coor[cdir]++;
				}
				offset -= coor[cdir];
				if (coor[cdir] != old_coor) {old_coor = coor[cdir];
					tmp[0]   = (buffer[((coor[cdir]+7) & 7)] - buffer[((coor[cdir]+0) & 7)]) / (-7*92160);
					tmp[1]   = (buffer[((coor[cdir]+6) & 7)] - buffer[((coor[cdir]+1) & 7)]) / (5*18432);
					tmp[2]   = (buffer[((coor[cdir]+5) & 7)] - buffer[((coor[cdir]+2) & 7)]) / (-3*10240.0f);
					coefs[7] = (buffer[((coor[cdir]+4) & 7)] - buffer[((coor[cdir]+3) & 7)]) / (18432.0f);
					coefs[1] = 225*tmp[0]+441*tmp[1]+(tmp[2]+ coefs[7] * 9.0f) * 1225.0f;
					coefs[3] = -259*tmp[0]-499*tmp[1]-1299.0f*tmp[2]+ coefs[7] * -1891.0f;
					coefs[5] = 35*tmp[0]+59*tmp[1]+75.0f*tmp[2]+ coefs[7] * 83.0f;
					coefs[7] = -tmp[0]-tmp[1]-tmp[2]-coefs[7];
					tmp[0]   = (buffer[((coor[cdir]+7) & 7)] + buffer[((coor[cdir]+0) & 7)]) / (-92160);
					tmp[1]   = (buffer[((coor[cdir]+6) & 7)] + buffer[((coor[cdir]+1) & 7)]) / (18432);
					tmp[2]   = (buffer[((coor[cdir]+5) & 7)] + buffer[((coor[cdir]+2) & 7)]) / (-10240.0f);
					coefs[6] = (buffer[((coor[cdir]+4) & 7)] + buffer[((coor[cdir]+3) & 7)]) / (18432.0f);
					coefs[0] = 225*tmp[0]+441*tmp[1]+(tmp[2]+ coefs[6] * 9.0f) * 1225.0f;
					coefs[2] = -259*tmp[0]-499*tmp[1]-1299.0f*tmp[2] + coefs[6] * -1891.0f;
					coefs[4] = 35*tmp[0]+59*tmp[1]+75.0f*tmp[2]+ coefs[6] * 83.0f;
					coefs[6] = -tmp[0]-tmp[1]-tmp[2]-coefs[6];

				}
				tmp[0] = coefs[7];
				for( j= 6;  j != 0 ; j--) tmp[0] = coefs[j] + (tmp[0] * offset) ;

				fout(cooro) = coefs[0] + (tmp[0] * offset);
			}

			for(dir=1;dir<nbdim;dir++) if (coor[((dir+cdir) % nbdim)] == dims[((dir+cdir) % nbdim)]-1) {coor[((dir+cdir) % nbdim)] = 0 ; cooro[((dir+cdir) % nbdim)] = 0;} else {coor[((dir+cdir) % nbdim)]++; cooro[((dir+cdir) % nbdim)] = coor[((dir+cdir) % nbdim)]; break;}
		} while (dir<nbdim);
	return(fout);
}


	LFHTEMP const DataGrid<C,nbdim> & DataGrid<C, nbdim>::blur(const GaussianDistribution<nbdim> &gauss){

		DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim>  fout;
		fout.setSizes(dims);
		Bluestein bl[nbdim];
		Tuple<unsigned int, nbdim> coor;

		Tuple<double, nbdim> vec;

		unsigned int maxb, tmp;
		double tmpd, sumd;
		for(unsigned int d=0;d<nbdim;d++){
			bl[d].setSize(dims[d]);
			tmp = bl[d].getBufferSize();
			if ((d == 0)||(tmp > maxb)) maxb = tmp;
			}
		typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb];



		ExOp::toZero(coor);
		ExOp::toZero(vec);
		DataGrid< double ,nbdim> dawin;
		dawin.setSizes(dims);

		sumd = 0.0f;
	//	ExOp::toZero(dawin);//	dawin(coor) = 0.25f;

//		sumd =1.0f;
//		coor[0] = coor[0] +4;dawin(coor) = 1.00f;
//		coor[1] = coor[1] +1;dawin(coor) = 0.25f;
//		coor[1] = coor[1] +1;dawin(coor) = 0.25f;

		do{
			gauss(tmpd, vec);
			sumd += tmpd;
			dawin(coor) = tmpd;
			for(maxb=0;maxb<nbdim;maxb++) {
				if (coor[maxb] == dims[maxb]-1) {coor[maxb] =0; vec[maxb] = 0.0f;}
				else {
					coor[maxb]++;
					vec[maxb] = (coor[maxb] < dims[maxb] >> 1) ? coor[maxb] : ((double)coor[maxb]) - dims[maxb];
					break;
				}
			}
		} while (maxb<nbdim);


		KeyIterator ite = KeyIterator(*this);

		DataGrid< typename ExCo<double>::COMPLEX_TYPE ,nbdim> dacon = dawin.FFtransform_2(bl, (typename  ExCo<double>::COMPLEX_TYPE * )buffer);

		DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ft = this->FFtransform_2(bl,buffer);


		if (ite.first()) do{
			(ft(ite())) *= dacon(ite());
		}while(ite.next());

		DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ift = ft.invFFtransform_2(bl,buffer);
		if (ite.first()) do{
			(*this)(ite()) = ( ExOp::mkrealproj(ift(ite())) ) /sumd;
		}while(ite.next());
		delete[](buffer);
		return(*this);
	}


	LFHTEMP const DataGrid<C,nbdim> & DataGrid<C, nbdim>::blur_zeroborder(const GaussianDistribution<nbdim> &gauss){

			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim>  fout;
			fout.setSizes(dims);
			Bluestein bl[nbdim];
			Tuple<unsigned int, nbdim> coor;

			Tuple<double, nbdim> vec;

			unsigned int maxb, tmp;
			double tmpd, sumd;
			for(unsigned int d=0;d<nbdim;d++){
				bl[d].setSize(dims[d]*2);
				tmp = bl[d].getBufferSize();
				if ((d == 0)||(tmp > maxb)) maxb = tmp;
			}
			typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb];




			ExOp::toZero(vec);
			DataGrid< double ,nbdim> dawin;

			for(maxb=0;maxb<nbdim;maxb++) coor[maxb] = dims[maxb]*2;

			dawin.setSizes(coor);

			sumd = 0.0f;

			//	ExOp::toZero(dawin);//	dawin(coor) = 0.25f;

			//		sumd =1.0f;
			//		coor[0] = coor[0] +4;dawin(coor) = 1.00f;
			//		coor[1] = coor[1] +1;dawin(coor) = 0.25f;
			//		coor[1] = coor[1] +1;dawin(coor) = 0.25f;

			ExOp::toZero(coor); do{
				gauss(tmpd, vec);
				sumd += tmpd;
				dawin(coor) = tmpd;
				for(maxb=0;maxb<nbdim;maxb++) {
					if (coor[maxb] == ((dims[maxb]) << 1)-1) {coor[maxb] =0; vec[maxb] = 0.0f;}
					else {
						coor[maxb]++;
						vec[maxb] = (coor[maxb] < dims[maxb]) ? coor[maxb] : ((double)coor[maxb]) - (dims[maxb] << 1);
						break;
					}
				}
			} while (maxb<nbdim);



			DataGrid< typename ExCo<double>::COMPLEX_TYPE ,nbdim> dacon = dawin.FFtransform_2(bl, (typename  ExCo<double>::COMPLEX_TYPE * )buffer);
		/*	{
				TiffFile tftf("./check.tif");
				//	if (ite.first()) do{
				//		dawin(ite()) = dacon(ite())[0];
				//	}while(ite.next());
				tftf.put(dawin, 0.0f, 1.0f);
			}*/


			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ft = this->FFtransform_zeroborder(bl,buffer);
			KeyIterator iteb = KeyIterator(dawin);
			if (iteb.first()) do{
				ft(iteb()) *= dacon(iteb());
			}while(iteb.next());

				DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ift = ft.invFFtransform_zeroborder(bl,buffer);

			KeyIterator ite = KeyIterator(*this);

			if (ite.first()) do{
				(*this)(ite()) = ( ExOp::mkrealproj(ift(ite())) ) /sumd;
			}while(ite.next());
//				LFH_ALIVE;
			delete[](buffer);
			return(*this);
		}

		LFHTEMP DataGrid<C,nbdim> & DataGrid<C, nbdim>::convolvecircle_zeroborder(double radius){

			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim>  fout;
			fout.setSizes(dims);
			Bluestein bl[nbdim];
			Tuple<unsigned int, nbdim> coor;

			Tuple<double, nbdim> vec;

			unsigned int maxb, tmp;
			double tmpd, sumd;
			for(unsigned int d=0;d<nbdim;d++){
				bl[d].setSize(dims[d]*2);
				tmp = bl[d].getBufferSize();
				if ((d == 0)||(tmp > maxb)) maxb = tmp;
			}
			typename ExCo<C>::COMPLEX_TYPE * buffer = new typename ExCo<C>::COMPLEX_TYPE[maxb];

			ExOp::toZero(vec);
			DataGrid< double ,nbdim> dawin;

			for(maxb=0;maxb<nbdim;maxb++) coor[maxb] = dims[maxb]*2;

			dawin.setSizes(coor);

			sumd = 0.0f;

			//	ExOp::toZero(dawin);//	dawin(coor) = 0.25f;

			//		sumd =1.0f;
			//		coor[0] = coor[0] +4;dawin(coor) = 1.00f;
			//		coor[1] = coor[1] +1;dawin(coor) = 0.25f;
			//		coor[1] = coor[1] +1;dawin(coor) = 0.25f;

			ExOp::toZero(coor); do{
				tmpd = (coor[0] < (dawin.dims[0]>>1)) ? coor[0] * coor[0] : (dawin.dims[0] - coor[0]) * (dawin.dims[0] - coor[0]);
				for(maxb=1;maxb<nbdim;maxb++) tmpd +=  (coor[maxb] < (dawin.dims[maxb]>>1)) ? coor[maxb] * coor[maxb] : (dawin.dims[maxb] - coor[maxb]) * (dawin.dims[maxb] - coor[maxb]);
				if (sqrt(tmpd) < radius+0.5f){
					if (sqrt(tmpd) < radius- 0.5f) {sumd += 1.0f;dawin(coor) = 1.0f;}
					else {sumd += (radius+ 0.5f) - sqrt(tmpd);dawin(coor) = (radius+ 0.5f) - sqrt(tmpd);}}
				for(maxb=0;maxb<nbdim;maxb++) {
					if (coor[maxb] == ((dims[maxb]) << 1)-1) {coor[maxb] =0; vec[maxb] = 0.0f;}
					else {
						coor[maxb]++;
						vec[maxb] = (coor[maxb] < dims[maxb]) ? coor[maxb] : ((double)coor[maxb]) - (dims[maxb] << 1);
						break;
					}
				}
			} while (maxb<nbdim);

			DataGrid< typename ExCo<double>::COMPLEX_TYPE ,nbdim> dacon = dawin.FFtransform_2(bl, (typename  ExCo<double>::COMPLEX_TYPE * )buffer);
			/* {
			typename  DataGrid< double ,nbdim>::KeyIterator ite = dawin.getKeyIterator();
			 TiffFile tftf("./check.tif");
			 	if (ite.first()) do{
			 		dawin(ite()) = dawin(ite());
			 	}while(ite.next());
			 tftf.put(dawin, 0.0f, 1.0f);
			 }*/


			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ft = this->FFtransform_zeroborder(bl,buffer);
			KeyIterator iteb = KeyIterator(dawin);
			if (iteb.first()) do{
				ft(iteb()) *= dacon(iteb());
			}while(iteb.next());

			DataGrid< typename ExCo<C>::COMPLEX_TYPE ,nbdim> ift = ft.invFFtransform_zeroborder(bl,buffer);

			KeyIterator ite = KeyIterator(*this);

			if (ite.first()) do{
				(*this)(ite()) = (ift(ite())[0]) /sumd;
			}while(ite.next());
			//				LFH_ALIVE;
			delete[](buffer);
			return(*this);
		}

	LFHTEMP DataGrid<C,nbdim> & DataGrid<C, nbdim>::blur_crude(const GaussianDistribution<nbdim> &gauss){
		unsigned int coor[nbdim];
		unsigned int dir,altdir;

		C* buffer;
		C current;

		unsigned int loop;
		double decay[nbdim];
		double peak[nbdim];
		double boundm[nbdim];
		double tmp;


		for(loop =0;loop < nbdim; loop++){

		 tmp = 1.0f-2.0f*gauss.ihvar.data[loop * (nbdim+1)];
//		 tmp = 1.0f+2.0f;
		decay[loop] = tmp - sqrt(tmp*tmp - 1.0f);
		peak[loop] = 2.0f*(1-decay[loop])/(1+decay[loop]);
		boundm[loop] = 1.0f / (1 - decay[loop]);
		}


		unsigned int cap;
		// peak weight is: (1-decay)/(1+decay), so that total weight is 1
		for(loop =0;loop < 2; loop++){
			for(dir =0;dir<nbdim;dir++){
				memset(coor, '\0',sizeof(int)*nbdim);
				buffer = new C[getDimSize(dir)];
				while(true){
					//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
					// prepare initial for backward!
					coor[dir] = getDimSize(dir)-1;

					buffer[coor[dir]] = (*this)(coor) * boundm[dir];  // / ((decay -1) * -1);
					for(coor[dir]--;coor[dir]!= 0xFFFFFFFF;coor[dir]--) buffer[coor[dir]] = (buffer[coor[dir]+1]  * decay[dir]) + (*this)(coor);

					//			printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
					// prepared initial for forward!
					cap = getDimSize(dir);
					coor[dir] = 0;
					current = (*this)(coor);

					current *= boundm[dir];

					buffer[coor[dir]] += current;
					(*this)(coor) =  buffer[coor[dir]] * peak[dir];
					for(coor[dir]++;coor[dir]<cap;coor[dir]++) {
						current *= decay[dir];
						(buffer[coor[dir]]) += current;
						current += (*this)(coor);
						(buffer[coor[dir]]) *= peak[dir];
						(*this)(coor) = buffer[coor[dir]];

					}


					// iterate througth other dimentions!
					for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == getDimSize((dir + altdir) % nbdim)-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
					if (altdir == nbdim) break;
				}
				delete[](buffer);
			}
		}
	return(*this);
	}


LFHTEMP DataGrid<C,nbdim> & DataGrid<C, nbdim>::blur_crude_dir(unsigned int dir, double wind){
	unsigned int coor[nbdim];
	unsigned int altdir;

	C* buffer;
	C current;

	unsigned int loop;
	double decay;
	double peak;
	double boundm;
	double tmp;


	tmp = 1.0f-wind*wind;

	decay = tmp - sqrt(tmp*tmp - 1.0f);
	peak = 2.0f*(1-decay)/(1+decay);
	boundm = 1.0f / (1 - decay);


	unsigned int cap;
	// peak weight is: (1-decay)/(1+decay), so that total weight is 1
	for(loop =0;loop < 2; loop++){
			memset(coor, '\0',sizeof(int)*nbdim);
			buffer = new C[getDimSize(dir)];
			while(true){
				//		printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
				// prepare initial for backward!
				coor[dir] = getDimSize(dir)-1;

				buffer[coor[dir]] = (*this)(coor) * boundm;  // / ((decay -1) * -1);
				for(coor[dir]--;coor[dir]!= 0xFFFFFFFF;coor[dir]--) buffer[coor[dir]] = (buffer[coor[dir]+1]  * decay) + (*this)(coor);

				//			printf("%i,%i,%i,in here!\n",coor[0],coor[1],coor[2]);fflush(stdout);
				// prepared initial for forward!
				cap = getDimSize(dir);
				coor[dir] = 0;
				current = (*this)(coor);

				current *= boundm;

				buffer[coor[dir]] += current;
				(*this)(coor) =  buffer[coor[dir]] * peak;
				for(coor[dir]++;coor[dir]<cap;coor[dir]++) {
					current *= decay;
					(buffer[coor[dir]]) += current;
					current += (*this)(coor);
					(buffer[coor[dir]]) *= peak;
					(*this)(coor) = buffer[coor[dir]];
				}
				// iterate througth other dimentions!
				for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == getDimSize((dir + altdir) % nbdim)-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
				if (altdir == nbdim) break;
			}
			delete[](buffer);
	}
	return(*this);
}

/*	LFHTEMP
	DataGrid<C,nbdim> & DataGrid<C, nbdim>::resize(double* factors){
		 DataGrid<C,nbdim> out[2];
		 if (nbdim == 1){
			 if ((factors[0] != 1.0f)&&(factors[0] <= 0.0f)){
			 out[0] = this->resizedirection(factors[0],0);
			 (*this) = out[0];
			 }
		 }else{
		 int c=2;
		 int i;
		 for(i=0;i<nbdim;i++){
			 if ((factors[i] != 1.0f)&&(factors[i] <= 0.0f)){
				 out[c & 1] = ((c & 2) ? *this : out[(c ^ 1)] ).resizedirection(factors[i],i);
				 c = (c ^ 1) & 1;
				 }
			 }
			 if (c != 2) (*this) = out[(c ^ 1)];
		}
		return(*this);
		}*/

	LFHTEMP
	DataGrid<unsigned int, nbdim> DataGrid<C, nbdim>::LabelConnected( bool (*isToLabel)(const C&), unsigned int *out_nblabels, Vector< Tuple<unsigned int,nbdim*2> >* out_boundingrect){
		DataGrid<unsigned int, nbdim> _out;
		_out.setSizes(dims);
		ExOp::toZero(_out);
		unsigned int curlabel = 0;
		int i;

		if (out_boundingrect){
			out_boundingrect->push_back(Tuple<unsigned int, nbdim*2>());
			for(i=0;i<nbdim;i++) {(*out_boundingrect)[0][i] = 0; (*out_boundingrect)[0][i+nbdim] = dims[i];}
		}

	//	Tuple<unsigned int, nbdim> neigh[nbdim*2];
		Tuple<unsigned int, nbdim> cur;

		KeyIterator ite = KeyIterator(*this);

		if (ite.first()) do{
			if (_out(ite()) == 0){
				if (isToLabel((*this)(ite()))) {
					cur = ite();
					curlabel++;
					if (out_boundingrect){
						out_boundingrect->push_back(Tuple<unsigned int, nbdim*2>());
						for(i=0;i<nbdim;i++) (*out_boundingrect)[curlabel][i] = (*out_boundingrect)[curlabel][i+nbdim] = cur[i];

					}


					stack< Tuple<unsigned int, nbdim> > stuff;

					stuff.push(cur);

					while(!stuff.empty()){
						cur = stuff.top(); stuff.pop();
						if ((_out(cur) == 0)&&(isToLabel((*this)(cur)))){
							_out(cur) = curlabel;
							if (out_boundingrect){
								for(i=0;i<nbdim;i++){
									if ((*out_boundingrect)[curlabel][i] > cur[i]) (*out_boundingrect)[curlabel][i] = cur[i];
									if ((*out_boundingrect)[curlabel][i+nbdim] < cur[i]) (*out_boundingrect)[curlabel][i+nbdim] = cur[i];
								}
							}
							for(i=0;i<nbdim;i++){
								if (cur[i] >0) {cur[i]--; stuff.push(cur); cur[i]++;}
								if (cur[i] <dims[i]-1) {cur[i]++; stuff.push(cur); cur[i]--;}
							}
						}
					}

					if (out_boundingrect) for(i=0;i<nbdim;i++){
						(*out_boundingrect)[curlabel][i+nbdim] -= (*out_boundingrect)[curlabel][i];
						(*out_boundingrect)[curlabel][i+nbdim] += 1;
					}
				}
			}

		} while (ite.next());

		if (out_nblabels != NULL) (*out_nblabels) = curlabel;


		return(_out);
	}

	LFHTEMP
	DataGrid<int, nbdim> DataGrid<C, nbdim>::directNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&) , bool equalitycount){
		Tuple<unsigned int, nbdim> coor;
		Tuple<unsigned int, nbdim> altcoor;
		DataGrid<int, nbdim> _out;

		_out.setSizes(dims);
		ExOp::toZero(_out);
		unsigned int dir;
		ExOp::toZero(coor);
		SetComparison tmpcmp;
		bool doit = true;
		do{
			for(dir = 0; dir<nbdim;dir++){
				if (coor[dir] == dims[dir]-1){
					//if (flag & 1) {
					//	altcoor = coor * KeyElem<unsigned int, unsigned int>(dir,0);
					//} else
					 doit = false;
				}else altcoor = coor + KeyElem<unsigned int, unsigned int>(dir,1);
				if (doit){
					tmpcmp = query( (*this)(coor), (*this)(altcoor));
					switch(tmpcmp & SETCMP_CMP_E_MASK){
						case SETCMP_EQUAL:
							if (equalitycount){
								_out(coor)++;
								_out(altcoor)++;
							}
							break;
						case SETCMP_GE:
							_out(coor)++;
							break;
						case SETCMP_LE:
							_out(altcoor)++;
							break;


					}
				}else doit = true;

			}





		for(dir=0;dir<nbdim;dir++) if (coor[dir] == dims[dir]-1) coor[dir] =0; else {coor[dir]++; break;}
		} while(dir < nbdim);
		return(_out);
	}
	LFHTEMP
	DataGrid<int, nbdim> DataGrid<C, nbdim>::indirectNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&) , bool equalitycount){
		Tuple<unsigned int, nbdim> coor;
		Tuple<unsigned int, nbdim> altcoor;
		DataGrid<int, nbdim> _out;

		_out.setSizes(dims);
		ExOp::toZero(_out);
		unsigned int dir, altdir;
		ExOp::toZero(coor);
		SetComparison tmpcmp;
		Tuple<unsigned int, nbdim> start;
		bool flag = false;
		do{

				for(dir =0;dir<nbdim;dir++){
					if (coor[dir] != dims[dir]-1){
						altdir = dir;
						altcoor[altdir] = ((flag & 1)&&( coor[altdir] == dims[altdir]-1)) ? 0 : coor[altdir] + 1;
						if (flag & 1){
						for(altdir++;altdir<nbdim;altdir++) start[altdir] = altcoor[altdir] = (coor[altdir] == 0) ? dims[altdir] - 1: coor[altdir]-1;
						}else{
						for(altdir++;altdir<nbdim;altdir++) start[altdir] = altcoor[altdir] = (coor[altdir] == 0) ? 0: coor[altdir]-1;
						}

						do {

							tmpcmp = query( (*this)(coor), (*this)(altcoor));
							switch(tmpcmp & SETCMP_CMP_E_MASK){
								case SETCMP_EQUAL:
									if (equalitycount){
										_out(coor)++;
										_out(altcoor)++;
									}
									break;
								case SETCMP_GE:
									_out(coor)++;
									break;
								case SETCMP_LE:
									_out(altcoor)++;
									break;
								default:
									printf("no!\%in", (int)(tmpcmp & SETCMP_CMP_E_MASK));
							}

							if (flag & 1){
								for(altdir=dir+1;altdir<nbdim;altdir++) if ((altcoor[altdir] == coor[altdir]+1)||((altcoor[altdir] == 0)&&(coor[altdir] == dims[altdir]-1))) altcoor[altdir] =start[altdir]; else {altcoor[altdir] -= (altcoor[altdir] == dims[altdir]-1) ? altcoor[altdir] : 1; break;}
							}else{
								for(altdir=dir+1;altdir<nbdim;altdir++) if ((altcoor[altdir] == coor[altdir]+1)||(altcoor[altdir] == dims[altdir]-1)) altcoor[altdir] =start[altdir]; else {altcoor[altdir]++; break;}
							}

						} while(altdir < nbdim);

					}
					altcoor[dir] = coor[dir];
				}

			for(dir=0;dir<nbdim;dir++) if (coor[dir] == dims[dir]-1) coor[dir] =0; else {coor[dir]++; break;}
		} while(dir < nbdim);
		return(_out);
	}

	LFHTEMP
	DataGrid<int, nbdim> DataGrid<C, nbdim>::knightNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&) , bool equalitycount){
		Tuple<unsigned int, nbdim> coor;
		Tuple<unsigned int, nbdim> altcoor;
		DataGrid<int, nbdim> _out;

		_out.setSizes(dims);
		ExOp::toZero(_out);
		unsigned int dir;
		ExOp::toZero(coor);
		SetComparison tmpcmp;
		bool doit = true;
		do{





			for(dir=0;dir<nbdim;dir++) if (coor[dir] == dims[dir]-1) coor[dir] =0; else {coor[dir]++; break;}
		} while(dir < nbdim);
		return(_out);
	}

	LFHTEMP unsigned int DataGrid<C, nbdim>::get_directNeightbor(const Tuple<unsigned int , nbdim> & coor, Tuple< Tuple<unsigned int , nbdim>, nbdim * 2> &_out )const{
		unsigned int nb_out=0;
		unsigned int dir;
        for(dir =0;dir<nbdim;dir++) {
                if (coor[dir] != 0) {_out[nb_out] = coor; _out[nb_out++][dir]--;}
                if (coor[dir] < dims[dir]-1) {_out[nb_out] = coor; _out[nb_out++][dir]++;}
            }
		return(nb_out);
	}

	LFHTEMP unsigned int DataGrid<C, nbdim>::get_indirectNeightbor(const Tuple<unsigned int , nbdim> & coor, Tuple< Tuple<unsigned int , nbdim>, (unsigned int)(TEMPLATE_INT_POWER<3,nbdim>::ans) -1u > &_out )const{
		unsigned int nb_out=0;
		unsigned int dir;

		Tuple<unsigned int , nbdim> cur;

		for(dir =0;dir<nbdim;dir++) cur[dir] = (coor[dir] == 0) ? 0 : coor[dir]-1;

		do{
			if (cur != coor) _out[nb_out++] = cur;
			for(dir=0;dir<nbdim;dir++) if ((coor[dir]+1 == cur[dir])||(cur[dir] == dims[dir]-1)) cur[dir] = (coor[dir] != 0) ? coor[dir]-1 : 0; else {cur[dir]++; break;}

		} while(dir < nbdim);

		return(nb_out);
	}


	LFHTEMP unsigned int DataGrid<C, nbdim>::get_knightNeightbor(const Tuple<unsigned int , nbdim> & coor, Tuple< Tuple<unsigned int , nbdim>,  nbdim * (nbdim -1) * 4   > &_out )const{
		unsigned int nb_out=0;
		unsigned int dir;
		unsigned int adir;
		unsigned int ccc;
		for(dir=0;dir<nbdim;dir++){
			for(adir=dir+1;adir<nbdim;adir++){

				ccc = (coor[adir] > 1 ? 3 : (coor[adir] == 1 ? 1 : 0)) |  (coor[adir] < dims[adir]-2 ? 12 : (coor[adir] == dims[adir]-2 ? 4 : 0));
				if (coor[dir] > 1) {
					if (ccc & 1) {_out[nb_out] = coor;_out[nb_out][dir] -= 2; _out[nb_out][adir] -= 1;nb_out++;}
					if (ccc & 2) {_out[nb_out] = coor;_out[nb_out][dir] -= 1; _out[nb_out][adir] -= 2;nb_out++;}
					if (ccc & 4) {_out[nb_out] = coor;_out[nb_out][dir] -= 2; _out[nb_out][adir] += 1;nb_out++;}
					if (ccc & 8) {_out[nb_out] = coor;_out[nb_out][dir] -= 1; _out[nb_out][adir] += 2;nb_out++;}
					} else if (coor[dir] == 1){
					if (ccc & 2) {_out[nb_out] = coor;_out[nb_out][dir] -= 1; _out[nb_out][adir] -= 2;nb_out++;}
					if (ccc & 8) {_out[nb_out] = coor;_out[nb_out][dir] -= 1; _out[nb_out][adir] += 2;nb_out++;}
				}
				if (coor[dir] < dims[dir]-2) {
					if (ccc & 1) {_out[nb_out] = coor;_out[nb_out][dir] += 2; _out[nb_out][adir] -= 1;nb_out++;}
					if (ccc & 2) {_out[nb_out] = coor;_out[nb_out][dir] += 1; _out[nb_out][adir] -= 2;nb_out++;}
					if (ccc & 4) {_out[nb_out] = coor;_out[nb_out][dir] += 2; _out[nb_out][adir] += 1;nb_out++;}
					if (ccc & 8) {_out[nb_out] = coor;_out[nb_out][dir] += 1; _out[nb_out][adir] += 2;nb_out++;}
				} else if (coor[dir] == dims[dir]-2) {
					if (ccc & 2) {_out[nb_out] = coor;_out[nb_out][dir] += 1; _out[nb_out][adir] -= 2;nb_out++;}
					if (ccc & 8) {_out[nb_out] = coor;_out[nb_out][dir] += 1; _out[nb_out][adir] += 2;nb_out++;}
				}
			}
		}
		return(nb_out);
	}


	LFHTEMP
	void DataGrid<C, nbdim>::FiniteDifference::operator()(DataGrid<C, nbdim>  &_out , DataGrid<C, nbdim> &_in) const{

		int coor[nbdim];
		memset(coor,'\0',sizeof(int)*nbdim);
		int dir;
		int dims[nbdim];
		C sum;
		for(dir=0;dir<nbdim;dir++) {dims[dir] = _in.getDimSize(dir);}
		while(true){

			C tmp = _in(coor);
			for(dir=0;dir<nbdim;dir++) {
				coor[dir] = coor[dir] == dims[dir]-1 ? 0 : coor[dir]+1;
				C tmp2 = _in(coor) -  tmp;
				if (dir == 0) sum = tmp2 * tmp2;
				else sum += tmp2 * tmp2;
				coor[dir] = coor[dir] == 0 ? dims[dir]-1 : coor[dir]-1 ;
			}
			_out(coor) = sum;
			for(dir=0;dir<nbdim;dir++) if (coor[dir] == dims[dir]-1) coor[dir] =0; else {coor[dir]++; break;}
			if (dir == nbdim) break;
		}


	}
	LFHTEMP
	void DataGrid<C, nbdim>::FiniteSum::operator()(DataGrid<C, nbdim>  &_out , DataGrid<C, nbdim> &_in) const{

		int coor[nbdim];
		memset(coor,'\0',sizeof(int)*nbdim);
		int dir;
		int dims[nbdim];
		C sum;
		for(dir=0;dir<nbdim;dir++) {dims[dir] = _in.getDimSize(dir);}
		while(true){

			C tmp = _in(coor);
			for(dir=0;dir<nbdim;dir++) {
				coor[dir] = coor[dir] == dims[dir]-1 ? 0 : coor[dir]+1;
				C tmp2 = _in(coor) + tmp;
				if (dir == 0) sum = tmp2 * tmp2;
				else sum += tmp2 * tmp2;
				coor[dir] = coor[dir] == 0 ? dims[dir]-1 : coor[dir]-1 ;
			}
			_out(coor) = sum;
			for(dir=0;dir<nbdim;dir++) if (coor[dir] == dims[dir]-1) coor[dir] =0; else {coor[dir]++; break;}
			if (dir == nbdim) break;
		}


	}

	LFHTEMP
	template <class D>
	void DataGrid<C, nbdim>::RankMap<D>::operator()(DataGrid<C, nbdim>  &_out , DataGrid<D, nbdim> &_in) const{

		vector< KeyElem<C, int> > pixellist;
		int totsize =1;
		int i;
		for(i=0;i<nbdim;i++) totsize *= _in.dims[i];
		for(i=0;i<totsize;i++) pixellist.push_back( KeyElem<C, int>(_in[i], i));

//		std::sort(pixellist.begin(), pixellist.end());

		_out[pixellist[0].d] = 0;
		C j = 0;
		for(i=1;i<totsize;i++) {
			if (pixellist[i-1] != pixellist[i]) j = ((C) i) / totsize;
			_out[pixellist[i].d] = j;
		}


	}

LFHTEMP
		DataGrid<C,2> DataGrid<C, nbdim>::pseudoInverse(DataGrid<C,nbdim>*L, DataGrid<C,nbdim>*R) const{
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!

			unsigned int coor[2];
		DataGrid<C,2> f_out;

		DataGrid<C,nbdim> mat[2];

		bool swap = dims[1] > dims[0];



		int na = swap ? 1 : 0;

		C* hh = new C[dims[na]];

		C* tmparr;
		if (L == NULL){
			coor[1] = coor[0] = dims[1];
			mat[na].setSizes(coor);
			mat[na].initeye();
		} else mat[na] = (*L);

			na = 1 - na;
		if (R == NULL){
			coor[0] = coor[1] = dims[0];
			mat[na].setSizes(coor);
			mat[na].initeye();
		} else mat[na] = (*R);

			if (swap){
				coor[0] = dims[1];
				coor[1] = dims[0];
				f_out.setSizes(coor);
				tmparr = ((Vector<C>&)f_out).darray;
				for(coor[0]=0;coor[0]<dims[0];coor[0]++)
					for(coor[1]=0;coor[1]<dims[1];coor[1]++) *(tmparr++) = (*((Vector<C>*)this))[coor[0] + coor[1] *dims[0]];
			}else{
				coor[0] = dims[0];
				coor[1] = dims[1];
				f_out.setSizes(coor);
				tmparr = ((Vector<C>&)f_out).darray;
				for(coor[1]=0;coor[1]<dims[1];coor[1]++)
					for(coor[0]=0;coor[0]<dims[0];coor[0]++) *(tmparr++) = (*((Vector<C>*)this))[coor[0] + coor[1] *dims[0]];
			}





		C sum2,tmp;

		int nb = 1-na;
		// bi-diagonalization!
		int i,j,k;
			unsigned int maxd = f_out.dims[0];
		for(i=0,j=0; i < f_out.dims[1]; j++){


			sum2 = ExCo<C>::intPow(f_out.data[j + i*maxd],2);
			for(k=j+1;k<maxd;k++) {
				hh[k-j] = f_out.data[k + i*maxd];
				sum2 += ExCo<C>::intPow(hh[k-j],2);
			}

			tmp += ExCo<C>::invintPow(sum2,2);


			if  (ExCo<C>::isnegative(f_out.data[j + i*maxd])){
			hh[0] =  f_out.data[j + i*maxd] - tmp;
			sum2 -= f_out.data[j + i * maxd] *tmp;
			}else{
				hh[0] =  f_out.data[j + i*maxd] + tmp;
				sum2 += f_out.data[j + i* maxd] *tmp;
			}


			f_out.rightHouseHolderMultiply(hh,maxd-j,sum2,false);

			mat[1].leftHouseHolderMultiply(hh,maxd-j,sum2,false);


			if (i+2 == f_out.dims[1]) break;
			i++;


			sum2 = ExCo<C>::intPow(f_out.data[j + i*maxd],2);
			for(k=i+1;k<f_out.dims[1];k++) {
				hh[k-i] = f_out.data[j + k*maxd];
				sum2 += ExCo<C>::intPow(hh[k-i],2);
			}

			tmp += ExCo<C>::invintPow(sum2,2);
			if  (ExCo<C>::isnegative(f_out.data[j + i*maxd])){
				hh[0] =  f_out.data[j + i*maxd] - tmp;
				sum2 -= f_out.data[j + i * maxd] *tmp;
			}else{
				hh[0] =  f_out.data[j + i*maxd] + tmp;
				sum2 += f_out.data[j + i * maxd] *tmp;
			}
			f_out.leftHouseHolderMultiply(hh,f_out.dims[1]-i,sum2,false);
			mat[0].rightHouseHolderMultiply(hh,f_out.dims[1]-i,sum2,false);


		}


			DataGrid<C,nbdim> tmptmp;
			printf("\n");
			f_out.show(stdout);

			tmptmp.settomatrixMultiply(mat[0],f_out);
			printf("\n");
			tmptmp.show(stdout);

			mat[0].settomatrixMultiply(tmptmp, mat[1]);

			printf("is initmatrix?\n");

			mat[0].show(stdout);
			fprintf(stderr, "old code used!\n");
			exit(1);

			// find zeroes!

			for(i=0,j=0; i < f_out.dims[1]; j++){
				tmp = ExCo<C>::intPow(f_out.data[j + i*maxd],2);
				if ((i ==0)||(sum2 > tmp)) tmp = sum2;
				i++;
				if (i == f_out.dims[1]) break;
				tmp = ExCo<C>::intPow(f_out.data[j + i*maxd],2);
				if ((i ==0)||(sum2 > tmp)) tmp = sum2;
			}
			sum2 *= ExCo<double>::epsilon();

			// if an element on the diagonal is 0, we must make its left neighbor 0 as well, rotation of rows bottom up!

			for(i=0; i < f_out.dims[1]; i++){
				tmp = ExCo<C>::intPow(f_out.data[j + i*maxd],2);
				hh[i] = (tmp < sum2) ? ExCo<C>::zero() : ExCo<C>::intPow(f_out.data[j + i*maxd],-1);
			}

			// multiply diagonal to left

			for(i=0; i < mat[0].dims[0]; i++){
				for(j=0; j < mat[0].dims[1]; j++){
					mat[0].data[i + j*mat[0].dims[0]] = hh[i] * mat[0].data[i + j*mat[0].dims[0]];
				}
			}


			if (swap){
				coor[0] = mat[1].dims[0];
				coor[1] = mat[0].dims[1];
				f_out.setSizes(coor);

				for(i=0; i < f_out.dims[0]; i++){
					for(j=0; j < f_out.dims[1]; j++){
						k = dims[0]-1;
						f_out.data[i + j * f_out.dims[0]] = mat[0].data[i + k * mat[0].dims[0]] * mat[1].data[k + j *mat[1].dims[0]];
						for(k--;k>=0;k--) f_out.data[i + j * f_out.dims[0]] += mat[0].data[i + k * mat[0].dims[0]] * mat[1].data[k + j *mat[1].dims[0]];
					}
				}
			}else{
				coor[0] = mat[0].dims[0];
				coor[1] = mat[1].dims[1];
				f_out.setSizes(coor);

				for(i=0; i < f_out.dims[0]; i++){
					for(j=0; j < f_out.dims[1]; j++){
						k = dims[1]-1;
						f_out.data[i + j * f_out.dims[0]] = mat[1].data[i + k * mat[1].dims[0]] * mat[0].data[k + j *mat[0].dims[0]];
						for(k--;k>=0;k--) f_out.data[i + j * f_out.dims[0]] += mat[0].data[i + k * mat[1].dims[0]] * mat[1].data[k + j *mat[0].dims[0]];
					}
				}
			}


			// combine Left and Right








		delete[](hh);
		return(f_out);
	}

	LFHTEMP
	DataGrid<C,nbdim> DataGrid<C, nbdim>::Inverse() const{
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!


	}

LFHTEMP template<class O> DataGrid<C,nbdim>& DataGrid<C, nbdim>::operator+=(KeyElem< Tuple<double,nbdim>, O > const & other){
	unsigned int bcoor =0;
	unsigned int k = nbdim-1;
	bcoor = (unsigned int) other.k[k];

	for(k--;k!= 0xFFFFFFFF;k--) bcoor = ((unsigned int)other.k[k]) + bcoor * ((unsigned int)dims[k]);
	data[bcoor] += other.d;

	return(*this);
	}


LFHTEMP DataGrid< Tuple<C, nbdim> ,nbdim> DataGrid<C, nbdim>::makeGradient() const{ DataGrid< Tuple<C, nbdim> ,nbdim> fout; fout.setSizes(dims);
        unsigned int dir,altdir;
        Tuple<unsigned int, nbdim> coor,altcoor;
        for(dir=0;dir<nbdim;dir++){
            ExOp::toZero(coor);
            while(true){
                altcoor = coor;
                coor[dir]=0;

                altcoor[dir] = coor[dir] +1;
                fout(coor)[dir] = (*this)(altcoor) *2.0f;
                altcoor[dir] = coor[dir];
                fout(coor)[dir] -= (*this)(altcoor) *2.0f;
                for(coor[dir]++;coor[dir]< dims[dir]-1;coor[dir]++){
                    altcoor[dir] = coor[dir] +1;
                    fout(coor)[dir] = (*this)(altcoor);
                    altcoor[dir] = coor[dir] -1;
                    fout(coor)[dir] -= (*this)(altcoor);
                }
                altcoor[dir] = coor[dir] ;
                fout(coor)[dir] = (*this)(altcoor) *2.0f;
                altcoor[dir] = coor[dir]-1;
                fout(coor)[dir] -= (*this)(altcoor) *2.0f;



                for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] == dims[(dir + altdir) % nbdim]-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
                if (altdir == nbdim) break;
            }
        }
        return(fout);
    }

LFHTEMP DataGrid<double,nbdim> DataGrid<C, nbdim>::makeGradientNorm() const{DataGrid<double,nbdim> fout; fout.setSizes(dims);
    DataGrid< Tuple<C, nbdim> ,nbdim> grgr = this->makeGradient();
    typename DataGrid<C, nbdim>::KeyIterator ite = this->getKeyIterator();
    Tuple<C, nbdim> tmp; unsigned int i; double tmpn;
    if (ite.first()) do{
        tmp = grgr(ite());
        tmpn = ExOp::pnorm(tmp[0]);
        for(i=1;i<nbdim;i++) tmpn += ExOp::pnorm(tmp[i]);
        fout(ite()) = sqrt(tmpn);
    } while(ite.next());
    return(fout);
    }

LFHTEMP DataGrid<unsigned int,nbdim> DataGrid<C, nbdim>::makeWatershedSortIndex(unsigned int &nbpix,  bool min_gradient_bassin, bool tree_skip_save) const{ // no_cheese inclues lone surrounded pixels


	//double imratio = ((double)totsize()) /  ((double)min_area_fold);
	unsigned int tsize;

  //  bool has_cheese_up = false;

	KeyIterator ite = getKeyIterator();

    Tuple< unsigned int, 3>* stat;
    Tuple<unsigned int, nbdim>* stat_coor;
    GaussElem<C>* stat_gaus;
    HeapTree< KeyElem<double, unsigned int> >* loc_merges;
    DataGrid<unsigned int,nbdim> fout;
    DataGrid<unsigned int,nbdim>* metapix;

    DataGrid< Tuple<unsigned int,nbdim>, nbdim> groups;
    GaussElem<C> prior,tmpprob; ExOp::toZero(prior);
    unsigned int i=0;
    unsigned int j,k,l;
	KeyElem<double,  unsigned int > toins;
    Tuple<unsigned int, nbdim> cursor;
	Tuple<Tuple<unsigned int, nbdim>, 2 * nbdim > ncoor;
	unsigned int nbind;
	myHashmap<unsigned int, unsigned int> exist;
    if (min_gradient_bassin){
        printf("Performing Watershed transform:"); fflush(stdout);
        DataGrid<double , nbdim> gradient = this->makeGradientNorm();
        fout = gradient.SegementClimbToMin(tsize);

        printf(" DONE\nPreparing Probabilistic clustering:"); fflush(stdout);
        groups.setSizes(dims);
        stat_gaus = new GaussElem<C>[tsize];
        stat = new Tuple< unsigned int, 3>[tsize << 1];
        stat_coor = new Tuple<unsigned int, nbdim>[tsize];
        loc_merges = new HeapTree< KeyElem<double, unsigned int> >[tsize];
        for(i=0;i<tsize;i++) {
            stat[i][0] = 0;

            stat[i][2] = i;
        }

        if (ite.first()) do{
            i = fout(ite());
            if (i == tsize) {i=0; fout(ite()) = 0;} // uses the 0 label
            if (stat[i][0] == 0){ // first pix!
                stat[i][0] = 1;
                stat_gaus[i] = GaussElem<C>((*this)(ite()));
                prior += stat_gaus[i];
                stat_coor[i] = ite();
                groups(ite()) = ite();
            }else{
                stat[i][0] +=1;
                tmpprob = GaussElem<C>((*this)(ite()));
                prior += tmpprob;
                stat_gaus[i] += tmpprob;
                groups(ite()) = groups(stat_coor[i]);
                groups(stat_coor[i]) = ite();
            }
        } while(ite.next());
         prior.setWeight(0.01f);

        if (ite.first()) do{
            stat_gaus[fout(ite())] += prior;
        } while(ite.next());

        printf(" DONE\nInitializing Merge Queue:"); fflush(stdout);
        for(i=0;i<tsize;i++){
            if (stat[i][0] == 0) continue;
            cursor = stat_coor[i];
            do{

                nbind = (*this).get_directNeightbor(cursor, ncoor);
                for(j=0;j<nbind;j++) exist[fout(ncoor[j])] =1;

            cursor = groups(cursor);
            }while(cursor != stat_coor[i]);

            l=0;
            for(j=0;j<exist.heap.getSize();j++){
                toins.d = exist.heap[j].k.k; // printf("%i\n", toins.d.second);
                if (toins.d < i) {
                    toins.k = stat_gaus[ i ].bhattacharryya_dist(stat_gaus[ toins.d ]);
                    if (ExCo<double>::isValid(toins.k)) loc_merges[ i ].insert(toins);
                }
                if (toins.d != i) l++;
            }
            stat[i][1] = l;
            exist.toMemfree();
        }
        metapix = new DataGrid<unsigned int , nbdim>(fout);
    } else {


        tsize = (*this).totsize();fout.setSizes(dims);

        groups.setSizes(dims);
        stat_gaus = new GaussElem<C>[tsize];
        stat = new Tuple< unsigned int, 3>[tsize << 1];
        stat_coor = new Tuple<unsigned int, nbdim>[tsize];
        loc_merges = new HeapTree< KeyElem<double, unsigned int> >[tsize];


        i=0;
        if (ite.first()) do{
            stat_gaus[i] = GaussElem<C>((*this)(ite()));
            prior += stat_gaus[i];
            fout(ite()) = i;
            groups(ite()) = ite();

            stat[i][0] = 0;
            stat[i][1] = 10;
            stat[i][2] = i; // = OxFFFFFFFF => done!

            stat_coor[i++] = ite();
        }while(ite.next());

        printf("computed prior! %i\n",(int)i); fflush(stdout);

        prior.setWeight(0.01f);

        for(i=0;i<tsize;i++) stat_gaus[i] += prior;

        printf("Initalizing for merging\n"); fflush(stdout);
        for(i=0;i<tsize;i++){
		nbind = (*this).get_directNeightbor(stat_coor[i], ncoor);
		for(j=0;j<nbind;j++) if (i > (toins.d = fout(ncoor[j]))) {
			toins.k = stat_gaus[i].bhattacharryya_dist(stat_gaus[toins.d]);
                if (ExCo<double>::isValid(toins.k)) {
                    loc_merges[i].insert(toins);
                }
			}
            stat[i][1] = nbind;
        }

    }




    HeapTree< KeyElem<double,  unsigned int> > merges;
    HeapTree< KeyElem<double,  unsigned int> > areas;

	for(i=0;i<tsize;i++){
        if (!(loc_merges[i].isEmpty())) {
            toins = loc_merges[i].top();
            toins.d = i;
            merges.insert(toins);

        }else  toins.d = i;
        toins.k = stat_gaus[i].w;

    }

	i =  tsize;


    printf("Start merging\n"); fflush(stdout);



	while(!merges.isEmpty()){

        if (!areas.isEmpty()) {
            // area too small now, forced merge!
            toins = areas.pop();
            if (stat[toins.d][2] == 0xFFFFFFFF) continue; // was merged, all good!
            stat[i][0] = toins.d;
            loc_merges[ stat[stat[i][0]][2] ].toMemfree();


            cursor = stat_coor[stat[stat[i][0]][2]];
            nbind = (*this).get_directNeightbor(cursor, ncoor);
            for(j=0;j<nbind;j++) {exist[fout(ncoor[j])] = 1;}
            exist[stat[i][0]] = 1;
            exist[0xFFFFFFFF] = 1;
            while (groups(cursor) != stat_coor[stat[stat[i][0]][2]] ){ // remove points from group is no *new* neighbors are close to elem
                nbind = (*this).get_directNeightbor(groups(cursor), ncoor);
                k=exist.getSize();
                for(j=0;j<nbind;j++) {
                    l = fout(ncoor[j]);
                    exist[l] = 1;
                }
                if (k == exist.getSize()) { fout(groups(cursor)) = 0xFFFFFFFF; groups(cursor) = groups(groups(cursor));}
                else cursor = groups(cursor);
            }
            l=0;
            for(j=0;j<exist.heap.getSize();j++){
                toins.d = exist.heap[j].k.k; // printf("%i\n", toins.d.second);
                if  (toins.d == stat[i][0]) continue;
                if  (toins.d == 0xFFFFFFFF) continue;
                l++;
                if (stat[ toins.d ][2] == 0xFFFFFFFF) {fprintf(stderr,"illegal? %i, %i\n", (int)toins.d , (int)tsize ); exit(1);}
                toins.k = stat_gaus[ stat[stat[i][0]][2]  ].bhattacharryya_dist(stat_gaus[ stat[ toins.d ][2] ]);
                if (ExCo<double>::isValid(toins.k)) toins.k = ExCo<double>::mkMaximum();
                loc_merges[ stat[stat[i][0]][2]  ].insert(toins);
            }
        //    printf("%i has %i neigh\n", stat[i][0],l);
            if (loc_merges[ stat[stat[i][0]][2]  ].isEmpty()) {
                    printf("try %i?\n",(int)stat[stat[toins.d][2]][0]);
                    for(j=0;j<exist.heap.getSize();j++) printf("%i ssible!\n", (int)exist.heap[j].k.k);
                    printf("%i got impossible!\n", (int)stat[i][0]); exit(1);
                    }
            exist.toMemfree();
            toins = loc_merges[ stat[stat[i][0]][2]  ].pop();
            stat[i][1] = toins.d;

        }else{
	//	printf("%i\n", merges.getSize());
		toins = merges.pop();

        l = stat[toins.d][2]; // printf("alive R! %i %i\n", toins.d ,l ); fflush(stdout);
        if (l == 0xFFFFFFFF) continue; // no longer exists
        stat[i][0] = toins.d;
        if (loc_merges[l].isEmpty()) thrbase.terminate("empty?\n");
        toins = loc_merges[l].pop();
     //   printf("alive Y! %i\n", toins.d); fflush(stdout);
   		if (stat[toins.d][2] == 0xFFFFFFFF) {
            //already done... NEXT!
      //     printf("alive T!\n"); fflush(stdout);
            while (!(loc_merges[l].isEmpty())) {
                toins = loc_merges[l].top();
                if (stat[toins.d][2] != 0xFFFFFFFF) {
                    toins.d = stat[i][0];
                    merges.insert(toins);
                    break;
                }
                loc_merges[l].pop();
            }
            continue;
        }
        stat[i][1] = toins.d;
        }



		if ((i % 500)==0) {
        printf("%i out of %i,   (merge dist = %e) \n",(int)i-tsize, (int)tsize -1, toins.k); fflush(stdout);
        }


        if (stat_gaus[stat[stat[i][0]][2]].w > stat_gaus[stat[stat[i][1]][2]].w){
            // swap, bigger area last!
            l = stat[i][0]; stat[i][0] = stat[i][1]; stat[i][1] = l;
        }
  //      printf("alive Z!\n"); fflush(stdout);
        stat[i][2] = stat[stat[i][0]][2]; // steals ID

      //  printf("alive A! %i\n",  stat[i][2]); fflush(stdout);

        stat_gaus[stat[i][2]] += stat_gaus[stat[stat[i][1]][2]];
  //      printf("aliveb! %i\n",stat[stat[i][1]][2]); fflush(stdout);


        loc_merges[stat[i][2]].toMemfree();
        loc_merges[stat[stat[i][1]][2]].toMemfree();


		// MARK AS DONE!
        if (stat[i][0]  == stat[i][1]) thrbase.terminate("got equal...\n");

        if ((stat[stat[i][0]][2] == 0xFFFFFFFF)||(stat[stat[i][1]][2]) == 0xFFFFFFFF) printf("WTH!\n");


        toins.k = stat_gaus[stat[i][2]].w;
        toins.d = i;



		if (exist.heap.getSize() != 0) printf("OGM! %i \n",(int)exist.heap.getSize());
		exist[stat[i][0]] = 7;
        exist[stat[i][1]] = 7;
        exist[i] = 7;
        exist[0xFFFFFFFF] = 7;
		// inserting new!
		cursor = stat_coor[stat[stat[i][0]][2]];
        fout(cursor) =i;
		nbind = (*this).get_directNeightbor(cursor, ncoor);
       for(j=0;j<nbind;j++) {l = fout(ncoor[j]);if (exist.find(l) == 0xFFFFFFFF) exist[l] = 1; else exist[l] |= 1;}


		while (groups(cursor) != stat_coor[stat[stat[i][0]][2]] ){ // remove points from group is no *new* neighbors are close to elem
            nbind = (*this).get_directNeightbor(groups(cursor), ncoor);
            k=0;
            for(j=0;j<nbind;j++) {
                l = fout(ncoor[j]);
                if (exist.find(l) == 0xFFFFFFFF) {exist[l] = 1; k=1;}
                else exist[l] |= 1;
                }
                if (k == 0) { fout(groups(cursor)) = 0xFFFFFFFF; groups(cursor) = groups(groups(cursor));}
                else {
                fout(groups(cursor)) = i;cursor = groups(cursor);
                }
		}

   //     printf("alive B!\n"); fflush(stdout);

		cursor = stat_coor[stat[stat[i][1]][2]];
        fout(cursor) =i;
		nbind = (*this).get_directNeightbor(cursor, ncoor);
       for(j=0;j<nbind;j++) {l = fout(ncoor[j]);if (exist.find(l) == 0xFFFFFFFF) exist[l] = 2;else exist[l] |= 2;}

		while (groups(cursor) != stat_coor[stat[stat[i][1]][2]] ){ // remove points from group is no *new* neighbors are close to elem
            nbind = (*this).get_directNeightbor(groups(cursor), ncoor);
            k=0;
            for(j=0;j<nbind;j++) {
                l = fout(ncoor[j]);
                if (exist.find(l) == 0xFFFFFFFF) {exist[l] = 2; k=1;}
                else exist[l] |= 2;
                }
                if (k == 0) {
                fout(groups(cursor)) = 0xFFFFFFFF; groups(cursor) = groups(groups(cursor));}
                else {
                fout(groups(cursor)) = i;cursor = groups(cursor);
                }
		}

		ncoor[0] = groups(stat_coor[stat[stat[i][0]][2]]); groups(stat_coor[stat[stat[i][0]][2]]) = groups(stat_coor[stat[stat[i][1]][2]]); groups(stat_coor[stat[stat[i][1]][2]]) = ncoor[0];
		stat[stat[i][0]][2] = 0xFFFFFFFF;
        stat[stat[i][1]][2] = 0xFFFFFFFF;

        l=0;
        for(j=0;j<exist.heap.getSize();j++){
            if ((exist.heap[j].k.d & 4) == 0) {
                toins.d = exist.heap[j].k.k; // printf("%i\n", toins.d.second);
                if (stat[ toins.d ][2] == 0xFFFFFFFF) {fprintf(stderr,"illegal? %i, %i\n", (int)toins.d , (int)tsize ); exit(1);}
                toins.k = stat_gaus[ stat[i][2] ].bhattacharryya_dist(stat_gaus[ stat[ toins.d ][2] ]);
                l++;
                if (ExCo<double>::isValid(toins.k)) loc_merges[ stat[i][2] ].insert(toins);

                if ((exist.heap[j].k.d & 3) ==3){
                    stat[stat[toins.d][2]][1] = stat[stat[toins.d][2]][1] -1;
                    if (stat[stat[toins.d][2]][1] == 1){
                    //    printf("%i ask for %i!\n", i, toins.d);
                        stat[stat[toins.d][2]][0] = i;
                        toins.k = stat_gaus[ stat[toins.d][2]  ].w;
                        areas.insert(toins);
                    }
                }
                /*
                if ((exist.heap[j].first.d & 3) ==3){ // check if surrournded!
                    cursor = stat_coor[stat[toins.d][2]];
                    nbind = (*this).get_directNeightbor(cursor, ncoor);
                    for(l=0;l<nbind;l++) if (i != fout(ncoor[l])) break;
                    if (l == nbind){

                       	for(cursor = groups(cursor);cursor != stat_coor[stat[toins.d][2]];cursor = groups(cursor)){ // remove points from group is no *new* neighbors are close to elem
                            nbind = (*this).get_directNeightbor(cursor, ncoor);
                            for(l=0;l<nbind;l++) if (i != fout(ncoor[l])) break;
                            if (l != nbind) break;

                        }
                        if (cursor == stat_coor[stat[toins.d][2]]){ // got surrounded by new area!
                            toins.k = stat_gaus[ stat[toins.d][2]  ].w;
                            areas.insert(toins);
        }
        }
        }*/

        }
        }

    //    printf("alive C!\n"); fflush(stdout);
        if (l == 1){ // has 1 neighbor only, merge now!
            toins.k = stat_gaus[ stat[i][2] ].w;
            toins.d = i;
            areas.insert(toins);
        } else stat[stat[i][2]][1] = l;
        exist.toMemfree();



        if (!(loc_merges[ stat[i][2] ].isEmpty())) {
            toins = loc_merges[ stat[i][2] ].top();
            toins.d = i;
            merges.insert(toins);
        }
        if (!(stat[i][2] < tsize)) exit(1);
        i++;


      //  if (i> 600000) break;
    }

    printf("merges finihed! %i inner, %i expected\n",(int)i, (int)tsize+tsize-1);fflush(stdout);


	for(j=0; j < i ;j++){
        if (stat[j][2] != 0xFFFFFFFF) break;
    }
    k = j +1;
   while (k < i){
        if (stat[k][2] != 0xFFFFFFFF) {
                stat[i][0] = j;
                stat[i][1] = k;
                stat[i][2] = stat[j][2];
                stat[j][2] = 0xFFFFFFFF;
                stat[k][2] = 0xFFFFFFFF;
                j = i++;
                if (i  == tsize+tsize-1) break;
            }
        k++;
    }

    printf("merges patching ! finihed! %i inner, %i expected\n",(int)i, (int)tsize+tsize-1);fflush(stdout);

	stack<unsigned int> s;

	j=0;
	s.push(i-1);
	while(s.size() > 0){
	i = s.top(); s.pop();

	if (i < tsize){ // isleaf!
        stat[i][0] =j;
		fout(stat_coor[i]) = j++;

	}else{
		s.push(stat[i][1]);
        s.push(stat[i][0]); if (tree_skip_save) j++;
	}

	}

    nbpix = j;
    if (min_gradient_bassin){

        if (ite.first()) do{
            fout(ite()) = stat[(*metapix)(ite())][0];
        }while(ite.next());

        delete(metapix);
    }

	printf("rendered segmentation done!\n");fflush(stdout);

    delete[](stat);
    delete[](stat_coor);
    delete[](loc_merges);
	delete[](stat_gaus);

	return fout;
}



LFHTEMP DataGrid<unsigned int, nbdim> DataGrid<C, nbdim>::makeSlicesIndexes(const Vector<double> &fractions)const{	DataGrid<unsigned int, nbdim> fout; fout.setSizes(dims);
	Vector< KeyElem< C, Tuple< unsigned int, nbdim > > > pixelList;
	typename DataGrid<C, nbdim>::KeyIterator ite = this->getKeyIterator();
	C sum, sum2; ExOp::toZero(sum); ExOp::toZero(sum2);
	if (ite.first()) do{
		sum += (*this)(ite());
		pixelList.push_back(KeyElem< C, Tuple< unsigned int, nbdim > >((*this)(ite()),ite()));
	} while(ite.next());
	pixelList.sort();
	unsigned int i=0;
	unsigned int j=0;
	for(; j< pixelList.getSize();j++){
		sum2 += pixelList[j].k;
		if (sum2 > sum * fractions[i]){
			i++;
			if (i == fractions.getSize()) break;
			}
			fout(pixelList[j].d) = i;
		}
	for(; j< pixelList.getSize();j++) fout(pixelList[j].d) = i;

	return fout;
	}
	 /*
	LFHTEMP
	DataGrid<C,nbdim> DataGrid<C, nbdim>::inverse_ofPosDef() const{ // inverse, requires positive definite!
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
		Tuple<unsigned int, 2> dims;
		getDims(dims);
		setSizes(dims);
		C* eigenval = new C[dims[0]];


		C* hh = new double[dims[0]];
		int i;

		unsigned int coor[2];
		C sum2;
		C tmp;

		DataGrid<C, 2> f_out;

		f_out.setSizes(dims);
		for(coor[0]=0;coor[0]<dims[0];coor[0]++) {
			for(coor[1]=0;coor[1]<dims[0];coor[1]++) {
				f_out(coor) = symmetric_matrix(coor);
			}
		}







		for(i=0;i<dims[0]-2;i++){

			coor[0] =i;
			sum2 = ExCo<C>::zero();
			for(coor[1]=coor[0]+2;coor[1]<dims[0];coor[1]++) {
				hh[coor[1]] = damat(coor);
				sum2 += hh[coor[1]] * hh[coor[1]];
			}

			coor[1]=coor[0]+1;
			sum2 += damat(coor)*damat(coor);
			tmp = ExCo<C>::invintPow(sum2, 2);
			hh[coor[1]] = (damat(coor)+ ((damat(coor) < 0) ? -tmp: tmp));
			tmp = (sum2 + fabs(damat(coor) *tmp));
			rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
			damat.leftHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
			damat.rightHouseHolderMultiply(hh + coor[1],dims[0]-coor[1],tmp,false);
		}



		for(coor[1]=0,coor[0]=0;;) {
			eigenval[coor[0]] = damat(coor);
			coor[0]++;
			if (coor[0] == dims[0]) break;
			tmp = damat(coor);

			sum2 = tmp / eigenval[coor[1]];

			rightoffdiagelimination( sum2,coor[0],coor[1]);
			coor[1]++;
			damat(coor) -= tmp * sum2;
		}
		delete[](hh);



		return(f_out);

	}*/

	LFHTEMP
	void DataGrid<C, nbdim>::update_pseudoInverse(const DataGrid<C,nbdim>& target){
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
		if ((dims[0] != target.dims[1])||(dims[1] != target.dims[0])) static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT;

	}

	LFHTEMP
	void DataGrid<C, nbdim>::update_Inverse(const DataGrid<C,nbdim>& target){
		LinkAssert<(nbdim ==2)> ass;// nbdims should be 2!
		if ((dims[0] != target.dims[1])||(dims[1] != target.dims[0])||(target.dims[0] != target.dims[1])) static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT;

	}

LFHTEMP void DataGrid<C, nbdim>::drawLine(const Tuple<unsigned int,nbdim> &a, const Tuple<unsigned int,nbdim> &b, const C& value){
	unsigned int i,be,t;
	be = a[0] >= b[0] ? a[0] - b[0] : b[0] - a[0];
	for(i=1;i<nbdim;i++) {t = a[i] >= b[i] ? a[i] - b[i] : b[i] - a[i]; if (t > be) be=t;}
	Tuple<unsigned int, nbdim> coor = a;
	Tuple<unsigned int, nbdim> res;
    for(i=0;i<nbdim;i++) res[i] = be >> 1;
	(*this)(coor) = value;
	for(t=0;t<be;t++){
		for(i=0;i<nbdim;i++) if (a[i] < b[i]) {
			if ( (res[i] += (b[i] - a[i])) >= be ) {res[i] -=be; coor[i]++;}
			}else if ( (res[i] += (a[i] - b[i])) >= be ) {res[i] -=be; coor[i]--;}
		(*this)(coor) = value;
}   }

LFHTEMP template<unsigned int size> void DataGrid<C, nbdim>::drawLine(const Tuple<unsigned int,nbdim-1> &a, const Tuple<unsigned int,nbdim-1> &b, const Tuple<C,size> &value){
	unsigned int i,be,t;
	be = a[0] - b[0];
	be = a[0] >= b[0] ? a[0] - b[0] : b[0] - a[0];
	for(i=1;i<nbdim-1;i++) {t = a[i] >= b[i] ? a[i] - b[i] : b[i] - a[i]; if (t > be) be=t;}

	Tuple<unsigned int, nbdim> coor;
	for(i=0;i<nbdim-1;i++) coor[i+1] = a[i];
	Tuple<unsigned int, nbdim-1> res; ExOp::toZero(res);
	for(coor[0]=0;coor[0]<dims[0];coor[0]++) (*this)(coor) = value[coor[0]];
	for(t=0;t<be;t++){
		for(i=0;i<nbdim-1;i++) if (a[i] < b[i]) {
			if ( (res[i] += (b[i] - a[i])) >= be ) {res[i] -=be; coor[i+1]++;}
		}else if ( (res[i] += (a[i] - b[i])) >= be ) {res[i] -=be; coor[i+1]--;}
		for(coor[0]=0;coor[0]<dims[0];coor[0]++) (*this)(coor) = value[coor[0]];
	}
	}

LFHTEMP void DataGrid<C, nbdim>::drawChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what){
	Tuple<unsigned int,nbdim> tmptmp = where;
	for(tmptmp[1] = where[1]; (tmptmp[1] != 0xFFFFFFFF) &&(tmptmp[1]+12 > where[1]);tmptmp[1]--){
		int daindex = (((int)which)-32) * 3  +((where[1] - tmptmp[1]) >> 2);
		int base = ((where[1] - tmptmp[1]) & 3) <<3;
		for(tmptmp[0] = where[0]; (tmptmp[0] < dims[0]) &&(tmptmp[0] < where[0]+8);tmptmp[0]++){
			if (((def_font[daindex]) >> (base + tmptmp[0] - where[0])) & 1) (*this)(tmptmp) = what;
		}
	}
}

LFHTEMP void DataGrid<C, nbdim>::drawCharFlip(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what){
	Tuple<unsigned int,nbdim> tmptmp = where;
	for(tmptmp[1] = where[1]; (tmptmp[1] < dims[1]) &&(tmptmp[1]< where[1]+ 12);tmptmp[1]++){
		int daindex = (((int)which)-32) * 3  +((tmptmp[1] - where[1]) >> 2);
		int base = ((tmptmp[1] - where[1]) & 3) <<3;
		for(tmptmp[0] = where[0]; (tmptmp[0] < dims[0]) &&(tmptmp[0] < where[0]+8);tmptmp[0]++){
			if (((def_font[daindex]) >> (base + tmptmp[0] - where[0])) & 1) (*this)(tmptmp) = what;
		}
	}
}

LFHTEMP void DataGrid<C, nbdim>::drawUpChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what){
	Tuple<unsigned int,nbdim> tmptmp = where;

	for(tmptmp[0] = where[0]; (tmptmp[0] != 0xFFFFFFFF) &&(tmptmp[0]+12 > where[0]);tmptmp[0]--){
		int daindex = (((int)which)-32) * 3  +((where[0] - tmptmp[0]) >> 2);
		int base = 7| (((where[0] - tmptmp[0]) & 3) <<3);
		for(tmptmp[1] = where[1]; (tmptmp[1] < dims[1]) &&(tmptmp[1] < where[1]+8);tmptmp[1]++){
			if (((def_font[daindex]) >> (base + where[1] - tmptmp[1])) & 1) (*this)(tmptmp) = what;
		}
	}
}

LFHTEMP void DataGrid<C, nbdim>::drawDownChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what){
	Tuple<unsigned int,nbdim> tmptmp = where;
	for(tmptmp[0] = where[0]; (tmptmp[0] < dims[2]) &&(tmptmp[0] < where[0]+12);tmptmp[0]++){
		int daindex = (((int)which)-32) * 3  +((tmptmp[0]- where[0]) >> 2);
		int base = ((( tmptmp[0] - where[0]) & 3) <<3);
		for(tmptmp[1] = where[1]; (tmptmp[1] < dims[1]) &&(tmptmp[1] < where[1]+8);tmptmp[1]++){
			if (((def_font[daindex]) >> (base + tmptmp[1] - where[1])) & 1) (*this)(tmptmp) = what;
		}
	}
}


LFHTEMP template<unsigned int size> void DataGrid<C, nbdim>::drawChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what){
	Tuple<unsigned int,size> tmptmp;
	for(tmptmp[0]=0;tmptmp[0]<nbdim-1;tmptmp[0]++) tmptmp[tmptmp[0]+1] = where[tmptmp[0]];

	for(tmptmp[2] = where[1]; (tmptmp[2] != 0xFFFFFFFF) &&(tmptmp[2]+12 > where[1]);tmptmp[2]--){
		int daindex = (((int)which)-32) * 3  +((where[1] - tmptmp[2]) >> 2);
		int base = ((where[1] - tmptmp[2]) & 3) <<3;
		for(tmptmp[1] = where[0]; (tmptmp[1] < dims[1]) &&(tmptmp[1] < where[0]+8);tmptmp[1]++){
			if (((def_font[daindex]) >> (base + tmptmp[1] - where[0])) & 1) {for(tmptmp[0]=0;tmptmp[0]<dims[0];tmptmp[0]++) (*this)(tmptmp) = what[tmptmp[0]];}
		}
	}
}

LFHTEMP template<unsigned int size> void DataGrid<C, nbdim>::drawCharFlip(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what){
	Tuple<unsigned int,size> tmptmp;
	for(tmptmp[0]=0;tmptmp[0]<nbdim-1;tmptmp[0]++) tmptmp[tmptmp[0]+1] = where[tmptmp[0]];
	for(tmptmp[2] = where[1]; (tmptmp[2] < dims[1]) &&(tmptmp[2]< where[1]+ 12);tmptmp[2]++){
		int daindex = (((int)which)-32) * 3  +((tmptmp[2]- where[1]) >> 2);
		int base = ((tmptmp[2] - where[1]) & 3) <<3;
		for(tmptmp[1] = where[0]; (tmptmp[1] < dims[1]) &&(tmptmp[1] < where[0]+8);tmptmp[1]++){
			if (((def_font[daindex]) >> (base + tmptmp[1] - where[0])) & 1) {for(tmptmp[0]=0;tmptmp[0]<dims[0];tmptmp[0]++) (*this)(tmptmp) = what[tmptmp[0]];}
		}
	}
}

LFHTEMP template<unsigned int size> void DataGrid<C, nbdim>::drawUpChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what){
	Tuple<unsigned int,size> tmptmp;
	for(tmptmp[0]=0;tmptmp[0]<nbdim-1;tmptmp[0]++) tmptmp[tmptmp[0]+1] = where[tmptmp[0]];

	for(tmptmp[1] = where[0]; (tmptmp[1] != 0xFFFFFFFF) &&(tmptmp[1]+12 > where[0]);tmptmp[1]--){
		int daindex = (((int)which)-32) * 3  +((where[0] - tmptmp[1]) >> 2);
		int base = 7| (((where[0] - tmptmp[1]) & 3) <<3);
		for(tmptmp[2] = where[1]; (tmptmp[2] < dims[2]) &&(tmptmp[2] < where[1]+8);tmptmp[2]++){
			if (((def_font[daindex]) >> (base + where[1] - tmptmp[2])) & 1) {for(tmptmp[0]=0;tmptmp[0]<dims[0];tmptmp[0]++) (*this)(tmptmp) = what[tmptmp[0]];}
		}
	}
}

LFHTEMP template<unsigned int size> void DataGrid<C, nbdim>::drawDownChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what){
	Tuple<unsigned int,size> tmptmp;
	for(tmptmp[1] = where[0]; (tmptmp[1] < dims[2]) &&(tmptmp[1] < where[0]+12);tmptmp[1]++){
		int daindex = (((int)which)-32) * 3  +((tmptmp[1]- where[0]) >> 2);
		int base = ((( tmptmp[1] - where[0]) & 3) <<3);
		for(tmptmp[2] = where[1]; (tmptmp[2] < dims[2]) &&(tmptmp[2] < where[1]+8);tmptmp[2]++){
			if (((def_font[daindex]) >> (base + tmptmp[2] - where[1])) & 1) {for(tmptmp[0]=0;tmptmp[0]<dims[0];tmptmp[0]++) (*this)(tmptmp) = what[tmptmp[0]];}
		}
	}
}



LFHTEMP DataGrid<C, nbdim>& DataGrid<C, nbdim>::toFractalTexture(uint32_t im_mag, double blur_aperture, bool isSigned){

	DataGrid<double,2> im;
	Tuple<unsigned int,2 >coor, acoor,tmp;
	coor[0] = 1<<im_mag; coor[1] = 1<<im_mag;
	im.setSizes(coor);	this->setSizes(coor);
	int sum,tq;

	unsigned int coormask = (1 << im_mag) -1;

	unsigned int rmag = 1;
	short rrr;
	for(unsigned int mag = 1u << (im_mag-1);mag >0; mag >>= 1){

		if (mag == (1u << (im_mag-1))){
			ExOp::toRand(rrr); rrr /= 256;
			ExOp::toZero(coor); im(coor) = rrr;
			coor[0]= mag;  coor[1]= mag; im(coor) = 1-rrr;
		}else{
			for(coor[1]= mag; (coor[1] >> im_mag) == 0; coor[1] += mag*2){
				acoor[1] = coor[1] - mag;
				acoor[0] = 0;
				tq = im(acoor);
				for(coor[0]= mag; (coor[0] >> im_mag) == 0; coor[0] += mag*2){
					sum = tq;
					acoor[1] = (coor[1] + mag) & coormask;
					sum += im(acoor);
					acoor[0] = (coor[0] + mag) & coormask;
					sum += im(acoor);
					acoor[1] = coor[1] -mag;
					sum += (tq = im(acoor));
					ExOp::toRand(rrr);
					im(coor) = ( (rrr >> 4) + (rrr >> rmag)) + (sum >> 2);
				}
			}
		}
		rmag++;

		for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] += mag){

			tmp[0] = (coor[1] > mag) ? coor[1] -mag : coor[1] + ((1 << im_mag) - mag) ;
			acoor[0] = 0;
			acoor[1] = coor[1];tq = im(acoor);
			for(coor[0] = mag ;(coor[0] >> im_mag) == 0; coor[0] += mag*2){
				sum = tq;
				acoor[0] = coor[0]; acoor[1] = tmp[0]; sum += im(acoor); acoor[1] =  coor[1] + mag; sum += im(acoor);
				acoor[1] = coor[1];
				acoor[0] = (coor[0] + mag) & coormask;
				sum += (tq = im(acoor));
				 ExOp::toRand(rrr);
				im(coor) = ( (rrr >> 4) + (rrr >> rmag)) + (sum >> 2);
			}
			coor[1] += mag;

			tmp[0] = (coor[1] + mag) & coormask;

			acoor[0] = (1 << im_mag) - mag;
			acoor[1] = coor[1]; tq = im(acoor);
			for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] += mag*2){
				sum = tq;
				acoor[0] = coor[0]; acoor[1] = coor[1] -mag; sum += im(acoor); acoor[1] = tmp[0]; sum += im(acoor);
				acoor[1] = coor[1]; acoor[0] += mag;
				sum += (tq = im(acoor));
				ExOp::toRand(rrr);
				im(coor) = ( (rrr >> 4) + (rrr >> rmag)) + (sum >> 2);
			}
		}
		rmag++;
	}

	if (blur_aperture != 0.0f) im.toblur(blur_aperture);

   WeightElem<double, 2> norm; ExOp::toZero(norm);
	for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] ++){
		for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] ++){
			norm += WeightElem<double, 2>(im(coor));
		}
	}
	double mmm[2]; mmm[0] = norm.getMean(); mmm[1] = pow(2.0f * norm.getVar(),-0.5f);
	if (isSigned){
		for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] ++){
			for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] ++){
				(*this)(coor) = ExCo<C>::mkMaximum() * erf((im(coor)  - mmm[0]) * mmm[1]);
			}
		}
	}else{
		for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] ++){
			for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] ++){
				(*this)(coor) = ExCo<C>::mkMaximum() * (erf((im(coor)  - mmm[0]) * mmm[1]) * 0.5f + 0.5f);
			}
		}
	}
	return *this;
}

LFHTEMP DataGrid<C, nbdim>& DataGrid<C, nbdim>::toFractalTexture(uint32_t im_mag, double blur_aperture, const C minV, const C maxV){

	DataGrid<double,2> im;
	Tuple<unsigned int,2 >coor, acoor,tmp;
	coor[0] = 1<<im_mag; coor[1] = 1<<im_mag;
	im.setSizes(coor);	this->setSizes(coor);
	int sum,tq;

	unsigned int coormask = (1 << im_mag) -1;

	unsigned int rmag = 1;
	short rrr;
	for(unsigned int mag = 1u << (im_mag-1);mag >0; mag >>= 1){

		if (mag == (1u << (im_mag-1))){
			ExOp::toRand(rrr); rrr /= 256;
			ExOp::toZero(coor); im(coor) = rrr;
			coor[0]= mag;  coor[1]= mag; im(coor) = 1-rrr;
		}else{
			for(coor[1]= mag; (coor[1] >> im_mag) == 0; coor[1] += mag*2){
				acoor[1] = coor[1] - mag;
				acoor[0] = 0;
				tq = im(acoor);
				for(coor[0]= mag; (coor[0] >> im_mag) == 0; coor[0] += mag*2){
					sum = tq;
					acoor[1] = (coor[1] + mag) & coormask;
					sum += im(acoor);
					acoor[0] = (coor[0] + mag) & coormask;
					sum += im(acoor);
					acoor[1] = coor[1] -mag;
					sum += (tq = im(acoor));
					ExOp::toRand(rrr);
					im(coor) = ( (rrr >> 4) + (rrr >> rmag)) + (sum >> 2);
				}
			}
		}
		rmag++;

		for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] += mag){

			tmp[0] = (coor[1] > mag) ? coor[1] -mag : coor[1] + ((1 << im_mag) - mag) ;
			acoor[0] = 0;
			acoor[1] = coor[1];tq = im(acoor);
			for(coor[0] = mag ;(coor[0] >> im_mag) == 0; coor[0] += mag*2){
				sum = tq;
				acoor[0] = coor[0]; acoor[1] = tmp[0]; sum += im(acoor); acoor[1] =  coor[1] + mag; sum += im(acoor);
				acoor[1] = coor[1];
				acoor[0] = (coor[0] + mag) & coormask;
				sum += (tq = im(acoor));
				 ExOp::toRand(rrr);
				im(coor) = ( (rrr >> 4) + (rrr >> rmag)) + (sum >> 2);
			}
			coor[1] += mag;

			tmp[0] = (coor[1] + mag) & coormask;

			acoor[0] = (1 << im_mag) - mag;
			acoor[1] = coor[1]; tq = im(acoor);
			for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] += mag*2){
				sum = tq;
				acoor[0] = coor[0]; acoor[1] = coor[1] -mag; sum += im(acoor); acoor[1] = tmp[0]; sum += im(acoor);
				acoor[1] = coor[1]; acoor[0] += mag;
				sum += (tq = im(acoor));
				ExOp::toRand(rrr);
				im(coor) = ( (rrr >> 4) + (rrr >> rmag)) + (sum >> 2);
			}
		}
		rmag++;
	}

	if (blur_aperture != 0.0f) im.toblur(blur_aperture);

   WeightElem<double, 2> norm; ExOp::toZero(norm);
	for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] ++){
		for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] ++){
			norm += WeightElem<double, 2>(im(coor));
		}
	}
	double mmm[2]; mmm[0] = norm.getMean(); mmm[1] = pow(2.0f * norm.getVar(),-0.5f);

    for(coor[1]= 0; (coor[1] >> im_mag) == 0; coor[1] ++){
        for(coor[0]= 0; (coor[0] >> im_mag) == 0; coor[0] ++){
            (*this)(coor) = minV + (maxV-minV) * (0.5f * erfc((im(coor)  - mmm[0]) * mmm[1]));
        }
    }

	return *this;
}

 // end of namescape
