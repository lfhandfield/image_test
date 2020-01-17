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

thread_local ThreadLog thread_log;
LFHPrimitive::ThreadBase thrbase;


ThreadBase::ThreadBase():nbactive(0), progb(20),async_progress_maintain(0xFFFFFFFF),nbthreads(0),running(true){
}

ThreadBase::~ThreadBase(){this->toMemfree(); stopThreadArray(); joinThreads();}
ThreadBase& ThreadBase::toMemfree(){while(nbactive) {thrds[--nbactive]->join();delete(thrds[nbactive]);}
    if (nbthreads) {delete[](thrds); nbthreads = 0;}
}

ERRCODE ThreadBase::startEvent(Event<void>* ev){
    if (nbactive == nbthreads) return 1;
    //thrds[nbactive++] = new std::thread(Event<void>::operator(), ev);
    thrds[nbactive] = new std::thread(ThreadBase::callThatEventVoid, ev); LFH_NICE_ALLOCERROR(thrds[nbactive],"")
    nbactive++;
    return 0;
}
void ThreadBase::waitForAllThreads(){
    if (async_progress_maintain != 0xFFFFFFFF){
        while(nbactive>0){async_progress_maintain = nbactive--; thrds[nbactive]->join(); delete(thrds[nbactive]);}
        progb.finish();
        async_progress_maintain = 0xFFFFFFFF;
    }else while(nbactive>0){async_progress_maintain = nbactive--; thrds[nbactive]->join(); delete(thrds[nbactive]);}
}
uint32_t ThreadBase::callThatEventVoid(Event<void>* ev) {return (*ev)();}
void ThreadBase::submit(std::function<void()> what){
    {
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(what);
        #endif // STDBINDFILTER
    }
    condvar.notify_one();
}
void ThreadBase::submit_ThenWait(std::function<void()> what){
    /*{
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(what);
        #endif // STDBINDFILTERwaitForAllThreads
    }

    std::unique_lock<std::mutex> innerlock(mut);
    while(true){
        condvar.wait(innerlock,[this]{return true;});
        if (todolist.empty()) break;
        std::function<void()> fout = std::move(todolist.back());
        todolist.pop_back();
        nbactivethr++;
        innerlock.unlock();
    }
    innerlock.unlock();*/
    what();
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

int ThreadBase::printLog(const char *__format, ...) {
    std::unique_lock<std::mutex> lck(thrbase.accessFinalMutex());
    va_list __va; va_start(__va, __format);
    // std::fprintf(log, "thr%i:", std::this_thread::get_id());
    int __res = std::vfprintf(log, __format, __va);
    va_end(__va);
    fflush(log);
return __res;}
int ThreadBase::printLogF(const char *__format, ...) {
    va_list __va; va_start(__va, __format);
    std::fprintf(log, "thr%X:", std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(log, __format, __va);
    va_end(__va);
    fflush(log);
return __res;}
int ThreadBase::printf(const char *__format, ...) {
    std::unique_lock<std::mutex> lck(thrbase.accessFinalMutex());
    va_list __va; va_start(__va, __format);
    std::fprintf(stdout,"thr%X:", std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(stdout, __format, __va);
    va_end(__va);
return __res;}
int ThreadBase::printf_l(const char *__format, ...) {
    va_list __va; va_start(__va, __format);
    std::fprintf(stdout,"thr%X:", std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(stdout,__format, __va);
    va_end(__va);
return __res;}
int ThreadBase::fprintf(FILE* f, const char *__format, ...) {
    if (f == NULL) f = log;
    std::unique_lock<std::mutex> lck(thrbase.accessFinalMutex());
    va_list __va; va_start(__va, __format);
    std::fprintf(f, "thr%X:", std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(stdout, __format, __va);
    va_end(__va);
return __res;}
int ThreadBase::fprintf_l(FILE* f,const char *__format, ...) {
    if (f == NULL) f = log;
    va_list __va; va_start(__va, __format);
    std::fprintf(f, "thr%X:", std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(f, __format, __va);
    va_end(__va);
    fflush(f);
return __res;}

int ThreadBase::print(const char *__format, ...) {
    va_list __va; va_start(__va, __format);
    char buffer[65536];
    int __res = std::vsprintf(buffer, __format, __va);
    va_end(__va);
    msgs.insert(string(buffer));
return __res;}


void ThreadBase::terminate(const char* __format, ...){
    running = false;
    va_list __va; va_start(__va, __format);
    char buffer[65536];
    int __res = std::vsprintf(buffer, __format, __va);
    va_end(__va);
    msgs.insert(string(buffer));
exit(1);}


std::function<void()> ThreadBase::getFunc(){
    std::unique_lock<std::mutex> lock(mut);
    condvar.wait(lock,[this]{return !todolist.empty(); });
    std::function<void()> fout = std::move(todolist.back());
    todolist.pop_back();
return(fout);}

void ThreadBase::operator()(){
    while(true){
        std::unique_lock<std::mutex> lock(mut); // ,std::defer_lock
        condvar.wait(lock,[this]{return !todolist.empty()||(!running); });
        if (!running) break;
        std::function<void()> fout = std::move(todolist.back());
        todolist.pop_back();
        nbactivethr++;
        lock.unlock();
        fout();
        lock.lock();
        nbactivethr--;
        endcondvar.notify_one();
    }
    condvar.notify_one();
}
void ThreadBase::joinAll(){
    std::unique_lock<std::mutex> lock(endmut);
    endcondvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
    /*while(true){
        mut.lock();
        if (!running) {mut.unlock(); break;}
        if (todolist.empty()){
            mut.unlock();
            {
                std::unique_lock<std::mutex> lock(endmut);
                endcondvar.wait(lock,[this]{return (nbactivethr==0)||(!running)||(!todolist.empty()); });
                lock.unlock();
            }
            if ((!running)||((nbactivethr==0)&&(todolist.empty()))) break;
        }else{
            std::function<void()> fout = std::move(todolist.back());
            todolist.pop_back();
            mut.unlock();
            fout();
        }
    }*/
}
ERRCODE ThreadBase::startThread(uint32_t thrID, std::function< int(uint32_t) > fnc){
    if (dedicated.find(thrID) != 0xFFFFFFFF) return 1;
    dedicated[thrID] = new std::thread(std::bind(fnc, thrID));
return 0;}
void ThreadBase::joinThread(uint32_t thrID){dedicated[thrID]->join(); delete(dedicated[thrID]); dedicated.erase(thrID);}
void ThreadBase::startThreadArray(uint32_t _nbthreads){
    running = true;
    nbactivethr =0;
    nbthreads = _nbthreads;
    mainthreadid = std::this_thread::get_id();
    #ifdef STDBINDFILTER
    for(int i=0;i<_nbthreads;i++) futures.emplace_back(std::bind(&ThreadBase::operator(), std::ref(*this)));
    #endif // STDBINDFILTER
}
void ThreadBase::stopThreadArray(){
    running = false;
    condvar.notify_all();
    fflush(stdout);
    for(int i=0;i<futures.size();i++) futures[i].join();
    futures.clear();
    flushMsgs();
}
void ThreadBase::joinThreads(){
    if (auto ite = dedicated.getIterator()) do{
        (*ite)->join();
        delete(*ite);
    }while(ite++);
    dedicated.toMemfree();
    flushMsgs();
}
void ThreadBase::flushMsgs(FILE *f){
    int cnt =0;
    string msgstar,old;
    if (msgs.pop(old)){
        while (msgs.pop(msgstar)){
        if (old == msgstar) cnt++;
        else {
        if (cnt) {std::fprintf(f,"(x%i) %s\n", cnt+1, old.c_str()); cnt=0;}
        else std::fprintf(f,"%s\n", old.c_str());
        old = std::move(msgstar);
        }
        }
        if (cnt) std::fprintf(f,"(x%i) %s\n", cnt+1, old.c_str());
        else std::fprintf(f,"%s\n", old.c_str());
    }
}
void ThreadBase::startEvent_ThenWait(Event<void>* ev){(*ev)();this->waitForAllThreads();}
void ThreadBase::test(){
    Tuple<int,4u> hoho;
    Tuple<double,4u> haha;
    /*hoho.show();
    std::printf("hello!\n");
    thrds[0] = CRT_THREAD(hoho, toZero);
    thrds[1] = CRT_THREAD(haha, toOne);
    thrds[2] = std::thread(static_call_test, 2);
    thrds[3] = std::thread(static_call_test, 5);
    std::printf("sent!\n");
    thrds[0].join();
    thrds[1].join();
    thrds[2].join();
    thrds[3].join();
    hoho.show();*/
    //thrds[0] = CRT_THREAD(hoho, operator+=, 17);

    //std::function<void()> _p = std::bind(Tuple<int,4u>::addition<double>, hoho,haha)

    //static_cast<void(*)(decltype(hoho), decltype(haha))>(operator==<char, std::string::traits_type, std::string::allocator_type>),

    // that lambda line
    //thrds[0] = std::thread( [&hoho](decltype(haha) rhs){return hoho += rhs * 7;}, haha);
    // that lambda line

   // thrds[0].join();
    std::printf("end comfirmed!\n");
}
void ThreadBase::startProgress(const char* text, uint32_t totalsteps){
    async_progress_maintain = 0u;
    async_progress = 0;
    async_progress_max = totalsteps;
    progb.start(text);
}
void ThreadBase::updateProgress(uint32_t threadID){
    if (threadID != async_progress_maintain) async_progress++;
    else progb.update(((double)async_progress++) / async_progress_max);
}
void ThreadBase::finishProgress(uint32_t threadID){
    if (threadID == async_progress_maintain) {async_progress_maintain = 0xFFFFFFFF; progb.finish();}
}

void ThreadBase::initEqualRanges(Tuple<uint32_t> &fout, uint32_t nbelems, uint32_t nb_threads) const{
    if (nb_threads == 0) nb_threads = nbthreads;
    fout.setSize(nb_threads+1);
    for(uint32_t i =0 ; i <= nb_threads;i++) fout[i] = (i *  nbelems) / nb_threads;
}

} // end of namespace

