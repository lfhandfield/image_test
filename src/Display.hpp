/*
 * display.hpp
 *
 * Copyright (C) 2019 Louis-Francois Handfield
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

namespace LFHDisplay{


#undef LFHTEMP
#define LFHTEMP template <class C>


LFHTEMP void EventQueue<C>::runTo(C arg, uint32_t n_time){
    if (getSync()){
        proc_async_routine();
        while( (p_queue.isEmpty() == false)&&((time = p_queue.top().time) - n_time > 0xF0000000)){
            Event<C>* nev = (Event<C>*)p_queue.pop().ev;
            uint32_t o = (*nev)(arg);
            if (o){
                if (o != 0xFFFFFFFF) this->insert_sync(time + o,nev);
            }else delete(nev);
            proc_async_routine();
        }
        time = n_time;
        freeSync();
    }
}

LFHTEMP void EventQueue<C>::startQueueThread(){}


LFHTEMP void EventQueue<C>::proc_async_routine(){
    while(async_write != async_read){
        p_queue.insert(buffer[async_read]);
        async_read = (async_read + 1) & 255;
    }
}

LFHTEMP void EventQueue<C>::insert_async(uint32_t time, Event<C>* ev){
    buffer[async_write] = EventUnit(time,ev);
    while(((async_write + 1) & 255) == async_read) {
        if (getSync()){
            this->proc_async_routine();
        }else thrbase.terminate("ASYNC SPIKE B!\n");
        freeSync();
    }
    async_write = (async_write + 1) & 255;
}

LFHTEMP void EventQueue<C>::insert_asap(Event<C>* ev){
    buffer[async_write] = EventUnit(time,ev);
    while(((async_write + 1) & 255) == async_read) {
        if (getSync()){
            this->proc_async_routine();
        }else thrbase.terminate("ASYNC SPIKE A!\n");
        freeSync();
    }
    async_write = (async_write + 1) & 255;
}


LFHTEMP bool EventQueue<C>::pop(Event<C>*& fout){ // always pop!
    this->proc_async_routine();
    if (p_queue.isEmpty()) return false;
    fout = (Event<C>*) p_queue.pop().ev;
    return true;
    }

LFHTEMP bool EventQueue<C>::pop_exec(C arg){ // always pop!
    this->proc_async_routine();
    if (p_queue.isEmpty()) return false;
    Event<C>* tmp = (Event<C>*)p_queue.pop().ev;
    uint32_t i = (*tmp)(arg);
    if (i+1 > 1) p_queue.insert(EventUnit(i,tmp));
    else if (i == 0) delete(tmp);
    return true;
    }
	/*
	template<class C> BitmapRessource::BitmapRessource(const LFHPrimitive::DataGrid<C,3> &input, const C* range){

			glGenBuffers(1 , &slot);
    		glBindBuffer(GL_PIXEL_UNPACK_BUFFER, slot);
			w = input.dims[1];
			h = input.dims[2];

    		unsigned char* buffer = new unsigned char[4 * w * h ];
    		unsigned char* cur= buffer;
    		unsigned int coor[3];

    		C fa,fm;
    		if (range){
	    		fa = range[0];
	    		fm = 255.0f / (range[1] - range[0]);

	    		}
    		switch(input.dims[0]){
	    		case 3:
    		for(coor[2]=0;coor[2]<input.dims[2];coor[2]++)
    		for(coor[1]=0;coor[1]<input.dims[1];coor[1]++) {
	    		cur[0] = 255;
	    		if (range){
	    		coor[0] =2; cur[1] = (input(coor)-fa) * fm;
	    		coor[0] =1; cur[2] = (input(coor)-fa) * fm;
	    		coor[0] =0; cur[3] = (input(coor)-fa) * fm;
	    		}else{
	    		coor[0] =2; cur[1] = input(coor);
	    		coor[0] =1; cur[2] = input(coor);
	    		coor[0] =0; cur[3] = input(coor);
    			}
	    		cur+= 4;
	    		}
	    		break;
	    		}


    		glBufferData(GL_PIXEL_UNPACK_BUFFER,sizeof(char)*4*input.dims[1]*input.dims[2],buffer,GL_STATIC_COPY);


			delete[](buffer);
	}

*/


template<class C> void GUITextArea::setDico(DicoElem<C>& whay){this->setDico(&whay.base);}

#undef LFHTEMP
#define LFHTEMP template <class TARG>
/*
 LFHTEMP AliasPtr<TARG>::AliasPtr(const AliasPtr<TARG>& input): alias(input.alias){
    if (alias) {
        unsigned int ite = AliasBank.find(alias);
        if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");exit(1);}
        AliasBank.deref(ite).second++;
    }
    }

 LFHTEMP AliasPtr<TARG>::AliasPtr(const TARG* t){
    unsigned int ite = AliasOf.find((void *)t);
    if (ite != 0xFFFFFFFF){
        alias = AliasOf.deref(ite);
        ite = AliasBank.find(alias);
        AliasBank.deref(ite).second++;
    }else{ alias = 0; printf("warning! target has no alias! cant link!\n");}
    }

 LFHTEMP  AliasPtr<TARG>& AliasPtr<TARG>::operator=(const AliasPtr<TARG>& input){
    clear();
    alias = input.alias;
    if (alias) {
        unsigned int ite = AliasBank.find(alias);
        if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");exit(1);}
        AliasBank.deref(ite).second++;
    }
    return (*this);
}

 LFHTEMP  AliasPtr<TARG>& AliasPtr<TARG>::operator=(const TARG* t){
    clear();
    printf("supper assign!\n");
    unsigned int ite = AliasOf.find((void *)t);
    if (ite != 0xFFFFFFFF){
        alias = AliasOf.deref(ite);
        printf("got %i assign!\n",alias );
        ite = AliasBank.find(alias);
        if (ite == 0xFFFFFFFF)  printf("got null ite!\n" );
        else AliasBank.deref(ite).second++;
    }else{ alias = 0; printf("warning! target has no alias! cant link!\n");}
    return(*this);
}

 LFHTEMP  void AliasPtr<TARG>::clear(){
    if (alias){
        unsigned int ite = AliasBank.find(alias);
        if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");exit(1);}
        if ((-- (AliasBank.deref(ite).second)) == 0) AliasBank.erase_from_iterator(ite);
        alias = 0;
    }
    }

 LFHTEMP  bool AliasPtr<TARG>::isValid()const {
    if (alias){
        unsigned int ite = AliasBank.find(alias);
        if (ite == 0xFFFFFFFF) {fprintf(stderr,"Alias entry deleted too early!\n");exit(1);}
        return AliasBank.deref(ite).first != NULL;
    }else return false;
}

 LFHTEMP const TARG* AliasPtr<TARG>::operator->()const{
   unsigned int ite =  AliasBank.find(alias);
    if (ite == 0xFFFFFFFF) return(NULL);
    else{
        return (TARG*) AliasBank.deref(ite).first;
    }
}

 LFHTEMP const TARG* AliasPtr<TARG>::operator*()const{
   unsigned int ite =  AliasBank.find(alias);
    if (ite == 0xFFFFFFFF) return(NULL);
    else{
        return (TARG*) AliasBank.deref(ite).first;
    }
}

 LFHTEMP AliasPtr<TARG>& AliasPtr<TARG>::memmove(AliasPtr<TARG>& from ){
    clear();
    alias = from.alias;
    from.alias =0u;
    return *this;
}
*/

template<class A> void GUIArea::insertGUI(A*a){ subs.push_back(a->GUI_alias); ctrl_state.gui_objects_ptr[a->GUI_alias].parent_alias = this->GUI_alias; }
template<class A> void GUIArea::removeGUI(A*a){ this->removeGUI(a->GUI_alias); }


#undef LFHTEMP
#define LFHTEMP template <class Q, class A, class I>
LFHTEMP ThreadQueue<Q,A,I>::ThreadQueue(A _scope): scope(_scope){}
LFHTEMP void ThreadQueue<Q,A,I>::kill(){int fout; keep_running = false; if (thr.joinable()) thr.join();}
LFHTEMP ERRCODE ThreadQueue<Q,A,I>::startThread(){
    thr = std::thread(std::bind(&ThreadQueue<Q,A,I>::operator(), std::ref(*this)));
    return thr.joinable() ? 0 : 1;
}
LFHTEMP uint32_t ThreadQueue<Q,A,I>::operator()(){
    keep_running = true;
    while(keep_running){
        while(this->run(scope));
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }
return 0;}
/*
#undef LFHTEMP
#define LFHTEMP template <class Q, class I>
LFHTEMP void ThreadQueue<Q,void,I>::kill(){int fout; keep_running = false; SDL_WaitThread(thread, &fout);}
LFHTEMP ERRCODE ThreadQueue<Q,void,I>::startThread(){
    thread = Controlstate::create_thread(this);
    return (thread == NULL) ? 1 : 0;
}
LFHTEMP uint32_t ThreadQueue<Q,void,I>::operator()(){
    keep_running = true;
    while(keep_running){
        while(this->run());
        SDL_Delay(1024);
        }
    return 0;
}
#undef LFHTEMP
#define LFHTEMP template <class Q, class A>
LFHTEMP ThreadQueue<Q,A,void>::ThreadQueue(A _scope): scope(_scope){}
LFHTEMP void ThreadQueue<Q,A,void>::kill(){int fout; keep_running = false; SDL_WaitThread(thread, &fout);}
LFHTEMP ERRCODE ThreadQueue<Q,A,void>::startThread(){
    thread = Controlstate::create_thread(this);
    return (thread == NULL) ? 1 : 0;
}
LFHTEMP uint32_t ThreadQueue<Q,A,void>::operator()(){
    keep_running = true;
    while(keep_running){
        while(this->run(scope));
        SDL_Delay(1024);
        }
    return 0;
}

*/

} // end of namespace
