/*
 * display.h
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

#ifndef _defined_LFHDisplay
#define _defined_LFHDisplay


#include "primitive.h"
//#include "Advanced.h"
//#include "Threedim.h"


//#include <GL/glew.h>
//#define GLFW_DLL
//#include <GLFW/glfw3.h>

//#include "./glad.h"


#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#define SHADER_STRINGIFY(V,A) "#version " #V "\n" #A
//#define AL_AL_H

//#ifdef AL_AL_H
//#include <AL/al.h>
//#include <AL/alc.h>
// #include <AL/alut.h>
//#endif


/*
#ifdef LFH_NESTED_LIBRARIES
#include <SDL/SDL.h>
#include <SDL/SDL_net.h>
#include <SDL/SDL_mixer.h>
#else${COMPILER_LIB}
*/

//#ifdef Rcpp_hpp
//#include "SDL2/SDL.h"
//#else
// #include "SDL2/SDL.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_net.h>
//#include "SDL_mixer.h"
//#endif

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

/*
	animation capture script

super 1liner script:

set scope:

*****************************
***  toggle VH addon     ***
*** select armature only ***
*** ROOT bone roll must be 0 ***
import os
basepath = "C:/prog/meshbuild/XXXX
*****************************

super do-it-all 1-liner script:
*****************************
for it in range(0,120): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()
for it in range(0,60): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()
for it in range(0,2): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()
for it in range(0,20): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()

for it in range(0,96): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()
for it in range(0,6): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()

ypos
for it in range(0,280): frstr = str(it + 10000); frstr = frstr[1:];bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()

super 1liner script:

*****************************
bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True);os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo();frstr = '0001'
*****************************

super 2-liner script for 10 frames:

*****************************
frstr = '0000'
frstr = frstr if frstr[3:] == '0' else frstr[:-2] + str(1+int(frstr[2:]));frstr;
for it in '0123456789': frstr = frstr[:-1] + it;bpy.context.scene.frame_set(int(frstr));bpy.ops.object.mode_set(mode='POSE');bpy.ops.ed.undo_push();bpy.ops.pose.armature_apply();bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+frstr+".obj",use_selection=True, axis_forward='-Y', axis_up='Z');os.remove(basepath+frstr+".mtl");os.remove(basepath+frstr+".obj");bpy.ops.ed.undo()
*****************************

mesh save from command line: (select mesh & bones first)
*****
bpy.ops.export_scene.obj('EXEC_DEFAULT', filepath=basepath+".obj",use_selection=True, axis_forward='Y', axis_up='Z');
****



windowed mode:

bpy.ops.export_scene.obj('INVOKE_REGION_WIN', filepath="d:/junk.obj",use_selection=True)

 use_triangles=False,
         use_edges=True,
         use_normals=False,
         use_uvs=True,
         use_materials=True,
         use_apply_modifiers=True,
         use_blen_objects=True,
         group_by_object=False,
         group_by_material=False,
         keep_vertex_order=False,
         use_vertex_groups=False,
         use_nurbs=True,
         use_selection=True,
         use_animation=False,
         global_matrix=None,
         path_mode='AUTO',
         use_bones=True


context, filepath="./fun.obj",
         False,
         True,
         False,
         True,
         True,
         True,
         True,
         False,
         False,
         False,
         False,
         True,
         True,
         False,
         None,
         'AUTO',
         True

*/


// server vs emulator

//#include <utility>
//#define __NO_STD_VECTOR // Use cl::vector instead of STL version
//#include <CL/cl.hpp>

using namespace LFHPrimitive;

/*
   Texture   vs   buffer
   power of 2     any
   1 stencil      many stencil, but separate
   3D             2.1D only
                  accumulate stuff
*/
namespace LFHDisplay{
class Controlstate;
class MyWindow;
class renderMode;
class GUIrenderMode;
class GUIObject;
class GUIArea;
class GUITextArea;
class GUIGrid;
#ifdef AL_AL_H
    class SoundRessource;
#endif


class InputState;
class Physical;

#define WAVS_uint32_t uint32_t
#define TEX2_uint32_t uint32_t
typedef uint32_t Guialias;

enum GUIMSG_Enum : uint32_t{ // internal and external events
    GUIMSG_NULL=0,
    LFHGUI_DRAW=1,
    LFHGUI_STATIC_DRAW=2,
    LFHGUI_ALIAS_DRAW=3,
    LFHGUI_WINDOW_OPENCLOSE=4,
    LFHGUI_SET_BG_TEXTURE=5,
    LFHGUI_SET_DATA=6, // data_ interpreted by target class
    LFHGUI_MOUSE_BUTTON, // fired for every push and click!
    LFHGUI_MOUSE_PUSH,
    LFHGUI_MOUSE_STARE,
    LFHGUI_QUERY_RECT_POS,
    LFHGUI_REQUEST_POINTER,
    LFHGUI_TEXTEDIT_DONE,
    LFHGUI_TEXTEDIT_CHANGE,
    LFHGUI_KEY_TYPE,
    GUIMSG_MOUSE_ENTER,
    GUIMSG_MOUSE_EXIT,
    GUIMSG_MOUSE_MOVE, // ongoing event
    GUIMSG_MOUSE_DRAG, // ongoing event
    GUIMSG_TAGGED_TEXT_CLICK,
    LFHGUI_SCOPE_FOCUS,
    LFHGUI_MOUSE_STENCILID_CHANGE,
    LFHGUI_QUITPROCESS,



    GUIMSG_TEXT_CHANGE,
	GUIMSG_VALUE_CHANGE,
	LFHGUIEVENT_NULL,
    LFHGUIEVENT_GOT_VISIBLE,
    LFHGUIEVENT_GOT_INVISIBLE,
    LFHGUIEVENT_FOCUS_GAIN,
    LFHGUIEVENT_FOCUS_LOSS,
    LFHGUIEVENT_MOUSE_STARE,

    GUIMSG_MOUSE_DOWN_LBUTTON, // fast response before drag or Click DClick
    GUIMSG_MOUSE_DRAG_LBUTTON, // resolved as a click
    GUIMSG_MOUSE_DROP_LBUTTON, // resolved as a click
    GUIMSG_MOUSE_CLICK_LBUTTON,
    GUIMSG_MOUSE_CONTEXT_CLICK_LBUTTON, // click while dragging
    GUIMSG_MOUSE_DCLICK_LBUTTON,// resolved as a long press click
    GUIMSG_MOUSE_DOWN_RBUTTON,
    GUIMSG_MOUSE_DRAG_RBUTTON, // resolved as a click
    GUIMSG_MOUSE_DROP_RBUTTON, // resolved as a click
    GUIMSG_MOUSE_CLICK_RBUTTON,
    GUIMSG_MOUSE_CONTEXT_CLICK_RBUTTON, // click while dragging
    GUIMSG_MOUSE_DCLICK_RBUTTON,// resolved as a long press click
    GUIMSG_MOUSE_DOWN_MBUTTON,
    GUIMSG_MOUSE_DRAG_MBUTTON, // resolved as a click
    GUIMSG_MOUSE_DROP_MBUTTON, // resolved as a click
    GUIMSG_MOUSE_CLICK_MBUTTON,
    GUIMSG_MOUSE_CONTEXT_CLICK_MBUTTON,
    GUIMSG_MOUSE_DCLICK_MBUTTON,// resolved as a long press click
    LFHGUIEVENT_MOUSE_WHEELUP,
    LFHGUIEVENT_MOUSE_WHEELDOWN
};

enum RELPOS_enum : uint32_t{ // window positioning
    RELPOS_CENTERED =0,
    RELPOS_LEFT =1,
    RELPOS_RIGHT =2,
    RELPOS_IDENTICAL =4,
    RELPOS_BEFORE =5,
    RELPOS_AFTER =6,
    RELPOS_TOP =1<<3,
    RELPOS_BOTTOM =2<<3,
    RELPOS_TL_CORNER =0x9,
    RELPOS_TR_CORNER =0xA,
    RELPOS_BL_CORNER =0x11,
    RELPOS_BR_CORNER =0x12,
    RELPOS_YIDENTICAL =4<<3,
    RELPOS_YBEFORE =5<<3,
    RELPOS_YAFTER =6<<3,
    };
enum mycustomattributes : uint32_t{ // shader_attribute offsets
	ATTRIBUTE_POSITION,            //0
	ATTRIBUTE_NORMAL,
	ATTRIBUTE_TEXTURE_COORDINATES,  //1
	ATTRIBURE_BONEID,
	ATTRIBURE_BONEWEIGHT,
	ATTRIBURE_CHARID,
	ATTRIBUTE_POSITION_ALT,
	ATTRIBUTE_NORMAL_ALT
};
enum Guitask_Enum{ // obsolete
	GUITASK_DRAW_STATIC=0,
};
enum LFHGUI_SHADERS_enum{
    LFHGUI_SHADERS_TEXT=0,
    LFHGUI_SHADERS_FRAME=2,
    LFHGUI_SHADERS_FRAME_ALIAS=3
};
enum GUITEXTURES_enum : uint32_t{
    GUITEXTURES_TEXT=0,
    GUITEXTURES_TEXTOFFSETS=1,
    GUITEXTURES_FRAME=2,
    GUITEXTURES_BUTTON=3,
    GUITEXTURES_SCROLL=4,
    GUITEXTURES_DEFAULT_PALETTE=5,
    GUITEXTURES_DEFAULT_PALTEX=6,
    GUITEXTURES_BUTTONGRID=7,
    GUITEXTURES_BUTTONGRID_MASK=8,
	GUITEXTURES_GUIWIDGETS=9
};
enum GUISOUND_enum : uint32_t{
    GUISOUND_MOUSEOVER=0,
    GUISOUND_TEXTOFFSETS=1,
    GUISOUND_FRAME=2,
    GUISOUND_BUTTON=3,
    GUISOUND_SCROLL=4,
    GUISOUND_DEFAULT_PALETTE=5
};
enum GUIATTRIB_enum{
	GUIATTRIB_ARROW_TL_COLOR1,
	GUIATTRIB_ARROW_BR_COLOR1

};
bool waitingForSuccess(ERRCODE output);

class GUImessage : public Event<void>{
    public:
    typedef unsigned int MSG_KEY;
    typedef unsigned int MSG_FILTER;
    unsigned int msg_key;
    GUIMSG_Enum type; // msg_filter;
    union{
        unsigned int parameter;
        void* target;
        unsigned short coor[2];
    };
    unsigned int parameter2;
    bool validFilter(const unsigned int &);
    unsigned int operator()();
    };
class ProcessState{
	public:
	virtual ~ProcessState(){}
	virtual int OnKeyDown(const SDL_KeyboardEvent& Event)=0; //OnKeyDown(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
	virtual int OnKeyUp(const SDL_KeyboardEvent& Event)=0; //OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
	virtual int OnMaintain()=0;//OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
	virtual int listen(GUImessage&)=0;
};
class StringcaptureProcess : public ProcessState{
public:
	~StringcaptureProcess();
	int OnKeyDown(const SDL_KeyboardEvent& Event); //OnKeyDown(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
	int OnKeyUp(const SDL_KeyboardEvent& Event); //OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
	int OnMaintain();//OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
	int listen(GUImessage&){return 1;}
};
class TextOptions{
public:
    unsigned int font;
    int shift[2];
    bool flip_ydir;
    char centered;
	TextOptions(){shift[0] =0;shift[1] =0; font = 10; flip_ydir = false;centered = 0;}
};
class TextDisplayData{
	public:
	char* text;
	uint16_t position[2];
	uint32_t color;
};
enum GUITEXT_LAYOUT_enum{
    GUITEXT_LAYOUT_TRUNCATE=0, //default
    GUITEXT_LAYOUT_WRAP, // requires dimention[0] !=0
    GUITEXT_LAYOUT_CENTERED,
    GUITEXT_LAYOUT_CENTERED_ITEMS,  // rows are superposed
    GUITEXT_LAYOUT_WRAP_CENTERED
};
class GuiTextAttribute;
class GUIStyle{
public:
    unsigned int text_borders[4];
//            LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> BG_tex;
//           LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> Border_tex;
    GUITEXT_LAYOUT_enum textlayout;
    unsigned int border_width;
    unsigned int scroll_width;

	TextOptions font_op;
	GLfloat colors[16][9];
	uint16_t colmap_index[12];
    uint16_t dimentions[2];

	myHashmap<unsigned int, void*> gui_attributes;
	GUIStyle();
	GUIStyle& toZero(){return *this;}

	void setDefDimentions(uint16_t width=0,uint16_t height=0);
	void setToAreaDefault();
	void setToMenuDefault();

	void setBorderSizes(int _border_size, int _scroll_size, bool update_text_borders = true);

	void setColorHIS(int color_index, double hue, double bright, double sat);
	void setColorHISpair(int color_index, double hue, double bright, double sat, double hue2, double bright2, double sat2);

	void setStateColormap(int channel, int color_index, int flag);

	bool setFrameUniforms(GLint shaderID, int32_t* rect, const int32_t* par_rect, bool isForeground = false);
	void setColorUniforms(GLint shaderID, int32_t state);

    int makeTextMesh(GuiTextAttribute &texta, const char* utf8string, int cursor_pos = -1);
    // offsets marks the newline starts
    int makeTextMesh(GuiTextAttribute &texta, const Vector<const char*> &utf8strings, Vector<uint32_t> *offsets = NULL);
    void makeTextMesh(GuiTextAttribute &texta, const Vector<TextDisplayData> &utf8strings);

    void drawTextMesh(const GuiTextAttribute &texta, const int32_t* rect, const int32_t* parrect);
    void drawTextMeshAlias(const GuiTextAttribute &texta, const int32_t* rect, const int32_t* parrect);
};
class RessourceLoader{
public:
    virtual void useTexture(const uint32_t)=0;
    virtual void useSound(const uint32_t)=0;
    virtual void useTextureExt(uint32_t ID)=0;
    virtual void allocTextureExt(uint32_t ID)=0;
    virtual void deallocTextureExt(uint32_t ID)=0;
};
template<class QUEUE, class ARGUMENT = typename QUEUE::INNER_TYPE, class INDEX = typename QUEUE::INDEX_TYPE>
class ThreadQueue : public QUEUE, public Event< >{  // got to get rid to this
public:
    bool keep_running;
    ARGUMENT scope;
    std::thread thr;
    ThreadQueue(ARGUMENT _scope);
    ERRCODE startThread();
    uint32_t operator()();
    void kill();
};

// has a custom compare method
class EventUnit{
public:
	static const bool IsPOD =true;
	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
    unsigned int time;
    void* ev;
    EventUnit(){}
    template<class C> EventUnit(unsigned int _time, Event<C>* _ev):time(_time), ev((void*)_ev){}
    	bool operator>(const EventUnit &o)const { return (time - o.time-1) < 0x80000000;}
		bool operator>=(const EventUnit &o)const{ return (time - o.time) <= 0x80000000;}
		bool operator<(const EventUnit &o)const{ return (o.time - time-1) < 0x80000000;}
		bool operator<=(const EventUnit &o)const{ return (o.time - time) <= 0x80000000;}
		bool operator==(const EventUnit &o)const{ return ev != o.ev; }
		bool operator!=(const EventUnit &o)const{ return ev == o.ev; }
};

template<class C>
class EventQueue{
    void proc_async_routine();
public:
    HeapTree< EventUnit > p_queue; // many thread, ordered in regions
    int semaphore;

    unsigned int time;
    EventUnit buffer[256];

    EventUnit  async_buffer[256];
    unsigned char async_read;
    unsigned char async_write;

    void ThreadProgram();
    void startQueueThread();

    EventQueue():semaphore(0),async_read(0),async_write(0){}
    void insert_async(unsigned int time, Event<C>* ev);
    void insert_asap(Event<C>* ev); /**<  */
    void runTo(C arg,unsigned int n_time);

    unsigned int getSize()const{return p_queue.getSize() + ((async_write - async_read) & 255);}

    bool getSync() {if (fetch_and_add(&semaphore, 1) == 0) return true; semaphore--; return false;}
    void freeSync() {semaphore--;}

    //{ requires to be Synched
    bool pop(Event<C>*& fout); // always pop!
    bool pop_exec(C arg); // always pop!
	void insert_sync(unsigned int time, Event<C>* ev){p_queue.insert(EventUnit(time, ev));} // can only be performed by the queue thread!
    void insert_sync_offset(unsigned int time_offset, Event<C>* ev){p_queue.insert(EventUnit(time_offset + time, ev));}
    //} requires to be Synched

};

template< >
class EventQueue<void>{
    void proc_async_routine();
public:
    HeapTree< EventUnit > p_queue; // many thread, ordered in regions
    int semaphore;

    unsigned int time;
    EventUnit buffer[256];

    EventUnit  async_buffer[256];
    unsigned char async_read;
    unsigned char async_write;

    void ThreadProgram();
    void startQueueThread();

    EventQueue():semaphore(0),async_read(0),async_write(0){}
    void insert_async(unsigned int time, Event< >* ev);
    void insert_asap(Event< >* ev); /**<  */
    void runTo(unsigned int n_time);

    unsigned int getSize()const{return p_queue.getSize() + ((async_write - async_read) & 255);}

    bool getSync() {if (fetch_and_add(&semaphore, 1) == 0) return true; semaphore--; return false;}
    void freeSync() {semaphore--;}

    //{ requires to be Synched
    bool pop(Event< >*& fout); // always pop!
    bool pop_exec(); // always pop!
	void insert_sync(unsigned int time, Event< >* ev){p_queue.insert(EventUnit(time, ev));} // can only be performed by the queue thread!
    void insert_sync_offset(unsigned int time_offset, Event< >* ev){p_queue.insert(EventUnit(time_offset + time, ev));}
    //} requires to be Synched

};

class Latent : public Event< >{
public:
    EventQueue<void > queue; // processes that can wait, and may take much time
    unsigned int state;
    unsigned int internal_state;
    virtual unsigned int operator()();
    };
enum INPUTSTATE_enum{
	INPUTSTATE_PRESSED_SHIFT=0,
	INPUTSTATE_PRESSED_MOUSELEFT=1,
	INPUTSTATE_PRESSED_MOUSERIGHT=2,
	INPUTSTATE_PRESSED_MOUSECENTER=3
};
// stores variables needed by all GUIobject that *also* are shared in a array of objects
class GUIHead{
public:
    GUIObject* target;
    Guialias parent_alias;
    GUIObject* operator->(){return target;}
    GUIHead& toZero(){target = NULL;return *this;}
    bool isValid(){return target != NULL;}
};
class Controlstate{
    friend class MyWindow;
    friend class GUImessage;
    friend class StringcaptureProcess;
    void compiledefaultshaders();
    void notifyStencilChange(int n_mouse_stencil, GLint *color_and_stencilID, bool hasextra); // only called from mywindow->render();
    Vector<ProcessState*> states;
    public:

    char def_path[256];
    int def_path_start;

    unsigned char volume_music;
    unsigned char volume_sounds;

    static bool init_openGL();
    static void check_openGL();
    static bool init_SDL(const char* const name, const char* const prod);
    static void clean_openGL();
#ifdef _SDL_NET_H
    static bool init_SDL_NET();
#endif

#ifdef _SDL_MIXER_H
    static void clean_SDL_mixer();
    static bool init_SDL_mixer();
#endif
#ifdef AL_AL_H
    static bool init_openAL();
    static void clean_openAL();
    ALCdevice *aldev;
    ALCcontext *alcon;
#endif
    static GLuint compileshader_routine(const char * const vcode, const char * const code, const char* shadername);
    static GLuint compileshader_advroutine(const char * const vcode, const char * const code, unsigned int flag, const char* shadername);

    static void close_procstates(unsigned int to_level);

    bool alias_storm;

    void create_latent_thread();

    void main_control_loop(RessourceLoader* res_loader);

    //uint32_t max_reserved_guialias;
    RessourceLoader* ressourceHanddle;

    GLuint datext_shader,sstext_shader;
    GLuint datext2_shader,sstext2_shader;
    GLuint daframe_shader;
    GLuint daframe_alias_shader;
    GLuint daicon_shader;
    GLuint daicon_alias_shader;
    GLuint daheat_shader;

   // GLuint quadindex_buffer; // TODO!
   // static void textChooseSize(Sint16* text, unsigned int f, unsigned short &x,unsigned short &y);

    char* usePrefPath(const char* const filename);

    unsigned short findLine(Uint16* &text, unsigned int s, unsigned int w, unsigned short& nb_space);
    //void drawtext(Uint16* text, unsigned int f,unsigned int x,unsigned int y,unsigned int &w,unsigned int &h);

    unsigned int wordwidth(Uint16* &text);
    void convert(Uint16*, const char*);
    void setForgroundObject(Guialias alias) {foreground_alias = alias;}
	void clearForgroundObject() {foreground_alias = 0u;}

    void enter_stringcapture_mode(GUITextArea* guitext = NULL);
    void exit_stringcapture_mode();

    MyWindow* curwin;
    unsigned short mouse_click_milli;
    unsigned int last_mouse_button_event_time;
    unsigned int mouse_stencil;
    uint32_t mouse_stencilID;
    uint32_t mouse_array_offset; // set if text or gui in within an array
    unsigned int mouse_coor[2];
    unsigned int mouse_push_coor[2];
    char stringCapture[1024];
    unsigned int stringCapture_cur;
    unsigned int stringCapture_endptr; // *not* string length (it is unicode!)
    uint32_t foreground_alias;
	uint32_t tooltipalias;
    uint32_t button_state;
    bool isLeftButtonDown()const {return ((button_state & 7) == 1) ;}
    bool isRightButtonDown()const {return ((button_state & 7) == 2) ;}
    bool isMiddleButtonDown()const {return ((button_state & 7) == 3) ;}
    bool isMouseDragging()const {return ((button_state & 8) != 0) ;}
    bool isLeftShiftDown()const {return ((button_state & 16) != 0) ;}
    bool isLeftCtrlDown()const {return ((button_state & 32) != 0) ;}
    bool isLeftAltDown()const {return ((button_state & 64) != 0) ;}


    LFHPrimitive::myHashmap<unsigned int, Event<GUImessage>*> gui_events;
    LFHPrimitive::myHashmap<unsigned int, GUIHead> gui_objects_ptr;
    LFHPrimitive::myHashmap<unsigned int, char> ThreadID_mask;
    LFHPrimitive::myHashmap<unsigned int, GUIStyle> gui_styles;

    EventQueue<void > timequeue;
    EventQueue<void > latentqueue; // processes ran by other thread

    Latent passive; // processes that can wait, and may take much time...

    bool fetchMessage(GUImessage&);
	bool query(INPUTSTATE_enum ) const;

	void manifestAsTooltip(uint32_t alias); // (r <=max_reserved_guialias)||
    uint32_t mkGuiAlias(){uint32_t r; do{ExOp::toRand(r);} while((gui_objects_ptr.find(r) != 0xFFFFFFFF)); return r;}

    void operator<<(ProcessState* proc){states.push_back(proc);}
    };
extern Controlstate ctrl_state;
// need to be created/deleted by openGL synced process;


   // special tags
   // substitution:
   // \\ -> \     \< -> <
   // \sX -> set space value
   // \SX -> set char spacer value
   // \fX -> set font size
   // \xX  \yX -> move cursor

class FormatedText{
public:
    unsigned char* txt; // special characters: \ and <
    unsigned int txtlen;
    unsigned int linkTag;
    unsigned int color;
    unsigned int BGcolor;
    unsigned char fontsize;
    unsigned char spacesize;
    unsigned char paragraphtag; // has newline/tab/doublespace...
    unsigned int w,h;
    FormatedText();
    ~FormatedText();

    void layout();
    void draw();
};
enum GUIFLAG_FLAGENUM{
    GUIFLAG_FLAGENUM_TEXT_MANIFEST= 0x08000000,
    GUIFLAG_FLAGENUM_CROP_BOUNDS= 0x10000000,
    GUIFLAG_FLAGENUM_DIRTYTEXT= 0x20000000,
    GUIFLAG_FLAGENUM_DISABLED= 0x40000000,
    GUIFLAG_FLAGENUM_INVISIBLE= 0x80000000
};
class GUIObject{
    friend class GUITextArea;
    friend class GUIButton;
    friend class GUIArray;
    friend class GUIArea;
    GUIObject(){} // used by allowed array constructors only
    void initialize_routine(uint32_t alias, uint32_t style_alias);
    virtual int32_t* accessRect()=0;
    virtual void wrSubPosition(int32_t *fout, const GUIObject* subptr) const{}
    public:
    static const bool IsPOD = false;
    static const bool IsComplex = false;
    static const bool NeedsAddLink = false;
    Guialias GUI_alias;
	uint32_t gui_flags; /* reserved: 0x80000000 [invisible flag] 0x40000000 [disabled flag] 0x20000000 [dirty-text] 0x10000000 [crop-boundaries]*/
    Guialias stare_alias; Uint32 stare_time;
    uint32_t styleID;

    GUIObject(unsigned int alias, unsigned int style_alias);
    GUIObject(unsigned int alias, unsigned int style_alias, const char* text);

    void wrPosition(int32_t *fout) const;

    virtual ~GUIObject(){}
 //   virtual void wrRect(unsigned int * target) const {memcpy(target,rect,sizeof(unsigned int)<<2); return;}
    virtual int hasScrollBar() const{return false;}
    virtual void setText(const char* _newtext){}
    virtual char* getText(){return(NULL);}
    virtual int nbTextLines(){return(0);}
    virtual void setPositionRelativeTo(const Guialias& other_alias, RELPOS_enum, unsigned int pad_x =0, unsigned int pad_y =0);
    inline void setPositionRelativeTo(const Guialias& other_alias, int value){setPositionRelativeTo(other_alias, (RELPOS_enum)value);}

    virtual void update()=0; // if the object is selected, it may receive updates
    void setDimentions(unsigned int width, unsigned int height);

    // { openGL synced

/** \brief Draw Object, needs to be synced with openGL
 *  \param par_rect, if absent, current object 'rect' is a absolute position, otherwise, 'rect' is relative to par_rect and may be cropped
 */
    virtual void draw(bool mouse_over, const int32_t* par_rect = NULL)=0;
/** \brief Draw Object, needs to render a color matching its GUI_alias on non-transparent pixels
 *  \param par_rect, if absent, current object 'rect' is a absolute position, otherwise, 'rect' is relative to par_rect and may be cropped
 */
    virtual void drawAlias(bool is_text,const int32_t* par_rect = NULL)=0;

	void setVisible(bool isVisible);
	bool isVisible()const{return (gui_flags & 0x80000000) == 0;}
	void setEnabled(bool isEnabled); // if disabled, cannot be selected clicked, only drawn
	bool isEnabled()const{return (gui_flags & 0x40000000) == 0;}


    virtual ERRCODE manifest()=0;
    virtual ERRCODE vanish()=0;
    // } openGL synced

//			virtual void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const =0;
//            virtual void drawstatic()=0; // OFFSCREEN RENDERING!
//            virtual bool loadstatic()=0; // return (isunchanged) // OFFSCREEN RENDERING!

    //virtual int GUImouseClick(unsigned int m_x, unsigned int m_y, unsigned int button_state_change){return 1;}
    virtual GUIMSG_Enum processGUIevent(const GUIMSG_Enum event){return event;}
//     virtual void reposition(const Tuple<unsigned int, 4> &rect);
    virtual void onResize(){}
    virtual int getIntValue()const{return (int)this->getValue();}
    virtual double getValue()const=0;

    virtual const int32_t* getRect()const=0;
    void setPos(const int32_t* pos);
    void setRect(const int32_t* rect);


    virtual GUIObject& deref(int index){return this[index];}
    virtual uint32_t getArrayIDfromMouse() const;
};
class GuiTextAttribute{
public:
    unsigned int glbuffer_indexes[2]; // indexes to render bg rect
    char* text;
//    uint16_t text_pos[2];


    /*
    struct TextBlock{
        int16_t pos[2];
        int16_t lenght;
        uint32_t link_alias;
    };
    Vector< TextBlock > draw_data;*/
    /*
    class GUITextLink : public GUIObject{
    public:
        char* target;
        GUITextLink();

        void draw(bool mouse_over, const int32_t* par_rect){}
        void drawAlias(bool is_text,const int32_t* par_rect){}
        ERRCODE manifest(){return 0;}
        ERRCODE vanish(){return 0;}
        double getValue()const{ return atof(target);}
    };*/

     GuiTextAttribute(): text(NULL){glbuffer_indexes[0] =0;}
     GuiTextAttribute(const char* _inatt){glbuffer_indexes[0] =0; unsigned int i = strlen(_inatt) +1;text = new char[i]; memcpy(text, _inatt,i);}
     GuiTextAttribute(const GuiTextAttribute&)=delete;
     GuiTextAttribute& operator=(const GuiTextAttribute& src)=delete;

     bool hasText() const{return glbuffer_indexes[0] != 0;}
     ~GuiTextAttribute(){delete[](text); if (glbuffer_indexes[0] !=0) glDeleteBuffers(1, glbuffer_indexes);}
     GuiTextAttribute& toMemmove(GuiTextAttribute& src){text = src.text; src.text =0; glbuffer_indexes[0] = src.glbuffer_indexes[0];glbuffer_indexes[1] = src.glbuffer_indexes[1]; src.glbuffer_indexes[0]=0; /*text_pos[0] =src.text_pos[0];text_pos[1] = src.text_pos[1];*/ return *this;}

};
enum Guiobject_State_Enum{
    GUIOBJECT_STATE_NULL=0,
    GUIOBJECT_STATE_VISIBLE =1,
    GUIOBJECT_STATE_ACTIVE =2,
    GUIOBJECT_STATE_SCROLLBAR =4,
    GUIOBJECT_STATE_CLOSE_BUTTON =8,
    GUIOBJECT_STATE_RESIZE_MOVE_BUTTON =16,
    GUIOBJECT_STATE_HAS_BORDERS =32
};
class GUIArea : public GUIObject{
    int32_t* accessRect(){return rect;}
    void wrSubPosition(int32_t *fout, const GUIObject* subptr) const;
    public:
//    LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> wind_bit;
//    LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> inner_wind_bit;
    const int32_t* getRect()const{return rect;}
    unsigned int bg_texture;
    unsigned int last_focus_time;
    Guiobject_State_Enum area_state;

    LFHPrimitive::Vector<Uint32> subs;
    int32_t rect[8];

    GUIArea(unsigned int alias, unsigned int style_alias) : GUIObject(alias,style_alias) {}

    virtual ~GUIArea(){}
    void update(); // if the object is selected, it may receive updates
    void draw(bool mouse_over, const int32_t* par_rect = NULL);
    void drawAlias(bool is_text,const int32_t* par_rect = NULL); // draw outline, no textures and
//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;
    ERRCODE manifest();
    ERRCODE vanish();

    void setText(char*);
    template<class A> void insertGUI(A*);
    template<class A> void removeGUI(A*);
    void insertGUI(Uint32, RELPOS_enum, unsigned int pad_x =0, unsigned int pad_y =0);
    void removeGUI(Uint32);
//    void insertGUI(const GUIObject& ob){this->insertGUI(ob.GUI_alias);}
    void removeGUI(const GUIObject& ob){this->removeGUI(ob.GUI_alias);}

	//void setForgroundSub(int alias);
    GUIMSG_Enum  processGUIevent(const GUIMSG_Enum event);

	void drawSubs();
	void drawAliasSubs(bool is_text);



	void setAreaDimentions(unsigned int width, unsigned int height, bool is_equal_to_real);

    double getValue()const{return 0;}

    uint32_t getArrayIDfromMouse()const;
};

// holds text with foreign aliases
class GUITextCloud : public GUIObject{
    int32_t* accessRect(){return rect;}
public:
    int nbdim;
    int32_t rect[4];
    Tuple< Tuple<double, 2u> > project;
    class TextContext {
    public:
        uint32_t tileID;
        char* text;
        TextContext();
        TextContext(const TextContext& other)=delete;
        TextContext(const char* _text);
        ~TextContext();
        TextContext& toMemfree();
        TextContext& toZero();
        TextContext& toMemmove(TextContext& other);
        TextContext& operator=(const TextContext& other)=delete;
        void show(FILE* f = stdout, int lvl =0)const;
    };

    class TileData{
    public:
        Tuple<uint32_t, 16u> index;
        Tuple<uint32_t, 17u> offsets;
        uint32_t dirty;
        GuiTextAttribute texta;
        TileData()=default;
        TileData(const TileData& other)=delete;
        TileData& operator=(const TileData& other)=delete;
        TileData& toMemmove(TileData& other);
        TileData& toZero(){return *this;}
        void show(FILE* f = stdout, int lvl =0)const;
    };
    // group of text sharing glbuffer (and font details)
    myHashmap<uint32_t, TileData> texta; // guialiases  and foreign aliases
    myHashmap<uint32_t, TextContext> textmap;
    std::function< Tuple<float, 3u>(uint32_t) > position_holder;

    void insert(uint32_t alias, char* text);
    GUITextCloud(unsigned int alias, unsigned int style_alias);
    ~GUITextCloud();
    void setText(Vector<string> datext); // uses memmove!
    void draw(bool mouse_over, const int32_t* par_rect);
    void drawAlias(bool is_text,const int32_t* par_rect);
    double getValue()const{return 0;}
    ERRCODE manifest(){return 0;}
    ERRCODE vanish(){return 0;}
    void update(){} // if the object is selected, it may receive updates
    const int32_t* getRect()const{return rect;}
};
// assumes subs are arrays
class GUIArray : public GUIObject{
    LFHPrimitive::Vector<Uint32> subs;
    int32_t* accessRect(){return rect;}
    uint32_t cursize;
    void wrSubPosition(int32_t *fout, const GUIObject* subptr) const;
    public:
    const int32_t* getRect()const{return rect;}
    int32_t rect[8];
    const int maxsize;

    int nb_columns;
    int selected;
    uint32_t offset[2];

    GUIArray(unsigned int maxnbitem, unsigned int alias, unsigned int style_alias);
    ERRCODE manifest(){return 0;}
    ERRCODE vanish(){return 0;}
    void setArrayFormat(int nb_columns, int elem_width, int elem_height);
    void insertGUI(Uint32, RELPOS_enum, unsigned int pad_x =0, unsigned int pad_y =0);
    void removeGUI(Uint32);
    void update();
    GUIMSG_Enum processGUIevent(const GUIMSG_Enum event);
    //int	GUImouseClick(unsigned int m_x, unsigned int m_y, unsigned int button_state_change);
    void draw(bool mouse_over, const int32_t* par_rect = NULL);
    void drawAlias(bool is_text,const int32_t* par_rect = NULL);
	double getValue()const{return selected;}
	void setCursize(int32_t _size);	int32_t getCursize()const{return cursize;}
	// void setPos();
	uint32_t getArrayIDfromMouse()const;
};
class GUIScroll : public GUIObject{
    int32_t* accessRect(){return rect;}

    public:
    int32_t rect[4];
    double  min;
    double  max;
    double  value;
    double  minstep; // or interval size
    uint32_t flags;

    GUIScroll(unsigned int alias,unsigned int style_alias);
	GUIScroll&	operator=(const GUIScroll&){exit(1);return *this;}
	const int32_t* getRect()const{return rect;}
    void    draw(bool mouse_over, const int32_t* par_rect = NULL);
    void    drawAlias(bool is_text,const int32_t* par_rect = NULL);
    void	update();
    ERRCODE	manifest();
    ERRCODE	vanish();
    double  getValue()const{return value;}
    void	onResize();
	void	setScrollingAndRange(double val, double minstep, double min, double max, bool has_range = false);

    //int	GUImouseClick(unsigned int m_x, unsigned int m_y, unsigned int button_state);
    GUIMSG_Enum 	processGUIevent(const GUIMSG_Enum event);
};
class GUIValueBar : public GUIObject{
    int32_t* accessRect(){return rect;}

    public:
    int32_t rect[4];
    double  value;
    double  pivot_value;
    double  max_value;
    char* cur_str;
    unsigned int state;
    GUIValueBar(unsigned int alias,unsigned int style_alias);
	GUIValueBar&	operator=(const GUIScroll&);
	~GUIValueBar();
	const int32_t* getRect()const{return rect;}
    void setText(const char*);
    void    draw(bool mouse_over, const int32_t* par_rect = NULL);
    void    drawAlias(bool is_text,const int32_t* par_rect = NULL);
    void	update();
    ERRCODE	manifest();
    ERRCODE	vanish();
    double  getValue()const{return value;}
    void	onResize();
    //int	GUImouseClick(unsigned int m_x, unsigned int m_y, unsigned int button_state);
    GUIMSG_Enum 	processGUIevent(const GUIMSG_Enum event);
};
class GUIButton : public GUIObject{
    int32_t* accessRect(){return rect;}

    GUIButton(): GUIObject(){}
    public:
//    LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> wind_bit;
//    LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> inner_wind_bit;
    GuiTextAttribute texta;
    int32_t rect[4];
    unsigned int	bg_texture;
    unsigned int	innerect[4];
    unsigned int	last_focus_time;
    unsigned int	state; //(Need New Text Mesh)(Has Text Mesh)()(is Disabled)()
    Guiobject_State_Enum area_state;
    GUIButton(unsigned int alias, unsigned int style_alias);
    GUIButton(unsigned int alias, unsigned int style_alias, const char* text);
    const int32_t* getRect()const{return rect;}
	static GUIButton* mkArray(uint32_t array_size, unsigned int alias, unsigned int style_alias);
    GUIObject& deref(int index){return this[index];}

    void	update(); // if the object is selected, it may receive updates
    void	draw(bool mouse_over, const int32_t* par_rect = NULL);
    void	drawAlias(bool is_text,const int32_t* par_rect = NULL);

    void	onResize();
//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;
    ERRCODE	manifest();
    ERRCODE	vanish();

    void	setText(const char*);
    double	getValue()const{return 0;}
    GUIMSG_Enum 	    processGUIevent(const GUIMSG_Enum event);
};

class MyWindow : public GUIArea{
    public:
    bool isfullscr;
    SDL_Window* Surf_Display;
    SDL_Renderer* renderer;
    SDL_GLContext glcontext;
    bool display_text_capture;
 //   LFHPrimitive::DataGrid<Uint32, 2> rgba_layer; // console layer
 //   LFHPrimitive::DataGrid<Uint32, 2> alias_layer; // console layer
 //   Uint32 last_mouse_alias;
//        LFHPrimitive::RessourcePtr<LFHDisplay::BitmapRessource> wind_bit;

    GLuint framebufferID[3];
//    GLuint renderbufferID;

    GLuint last_mouse_pos[2];
    float depth[9];
    GLint stencilval[4];
    float last_mouse_depth;
    float last_mouse_ddepth[2]; // estimates of derivative in depth

    LFHPrimitive::Vector<renderMode*> render_list;
    MyWindow() : GUIArea(0,0) {LFH_ALIVE; exit(1);}
    MyWindow(unsigned int alias,unsigned int style_alias,int sizex,int sizey,RELPOS_enum posis, bool full, int winid, const char* win_name = NULL);
    ~MyWindow();

    void setWindowName(const char* newname);


    LFHPrimitive::myHashmap<uint32_t, GUIHead> textlinks;


    void operator<<(renderMode *);
    void operator>>(renderMode *);


    GLvoid glPrint(const char *fmt, ...);
    void resize(void);
//	void setupDIB();
//	void setupPixelFormat();
//	void resizeDIB();
//	void setupPalette();

    void render();
    void renderAlias();
    void render_static();
    virtual void idleFunc();

    GUIMSG_Enum  processGUIevent(const GUIMSG_Enum event){return event;}
    void (*getIdleFunc())();
    void setProjection();
};
class renderEvent : public LFHPrimitive::Event< >{
    MyWindow* win;
    public:
    renderEvent(MyWindow* _win):win(_win){}
    unsigned int operator()(){
	win->render(); return 0;
	}
};
class renderMode{
    public:
    virtual void draw(MyWindow*)=0;
    virtual void drawPostGUI(MyWindow*){}
    virtual void drawAlias(MyWindow*)=0;
    virtual ~renderMode(){}
};
#ifndef LFHINTERRUPTS
#define LFHINTERRUPTS LFH_interrupts
enum LFH_interrupts{
    LFHI_CLOSE=0,
    LFHI_OPEN=1,
    LFHI_SAVE=2
};
#endif

enum lfhgui_type{
    LFHGUI_TYPE_NULL=0,
    LFHGUI_TYPE_FRAME=1, // may has menu and scroll bars
    LFHGUI_TYPE_TEXT=2
};



		/*
#define LFHGUI_STATE_INVISIBLE 1
#define LFHGUI_STATE_SELECTED 2
#define LFHGUI_STATE_DISABLED 3
		class GUIObject{
			public:
			unsigned int rect[4];
			lfhgui_type guitype;
			unsigned int id; // for messages;
			unsigned int state;
			void* data;
			GUIObject(){}
			GUIObject(lfhgui_type i_type, void* i_data, unsigned int* i_rect);
			void draw();
			void draw_ID();
			};
*/
			// text, button, glidebar, menum choices tabs


template<class C, int order>	class TimedValueFunction;
template<class C, int order=1>	class TimedValueFunction{
	public:
	LFHPrimitive::Tuple<C, order+1> data;
	int time;
	double scale;
	TimedValueFunction(){}
	TimedValueFunction(const C& i_data):scale(0.0f){data[0] = i_data;}
	void initScale(double i_scale){data[1] = LFHPrimitive::ExCo<C>::mkZero();scale = i_scale;}
	C operator()(int q_time = clock())const{
		if (scale == 0.0f) return(data[0]);
		else return(data[0] + (data[1] * (scale * (q_time - time)  )));
		}
	void setDeriv(const C& i_data, int i_time = clock()){
		data[0] += data[1] * (scale * (i_time - time));time = i_time;
		data[1] = i_data;
		}
	void addDeriv(const C& i_data, int i_time = clock()){
		data[0] += data[1] * (scale * (i_time - time));time = i_time;
		data[1] += i_data;
		}
	void setTimedChange(const C& i_target, int i_delay){ // todo
		data[0] += data[1] * (scale * (i_delay - time));time += i_delay;
		data[1] = i_target;
		}
    TimedValueFunction<C, 1>& toZero(){ExOp::toZero(data); return *this;}
	};

class InputState{
    public:
    unsigned int mouse[2];
    unsigned int mouse_click[2];
    unsigned int mouse_state;
    bool MouseLDown(){return(mouse_state & 1);}
};
extern InputState istate;
class MouseStare : public LFHPrimitive::Event< >{
public:
	Uint32 x,y;
	Uint32 oldtime;
	Uint32 oldalias;
	Uint32 windowalias;
	MouseStare(Uint32 _oldtime,Uint32 _oldalias,Uint32 _windowalias, Uint32 _x, Uint32 _y) :  x(_x), y(_y), oldtime(_oldtime), oldalias(_oldalias), windowalias(_windowalias){}
	unsigned int operator()();
};
class GUITaskEvent : public LFHPrimitive::Event< >{
public:
	Guitask_Enum task;
	Uint32 guiID;
	GUITaskEvent(Uint32 _guiID,Guitask_Enum _task) : task(_task),guiID(_guiID){}
	unsigned int operator()();
};
class GUITextArea : public GUIObject{ // can be modified!
	bool isStatic() const {return(state & 1);}
	bool isSelected() const {return(state & 2);}
	bool isInModification() const {return(state & 4);}
    int32_t* accessRect(){return rect;}

	GUITextArea(): GUIObject(){}
	public:
    int32_t rect[4];
	unsigned int state;
	char* cur_str;
	Dictionary* dico;
	unsigned int dico_index;
    GuiTextAttribute texta;

	unsigned int modstamp;

	GUITextArea(unsigned int alias, unsigned int style_alias);
	static GUITextArea* mkArray(uint32_t array_size, unsigned int alias, unsigned int style_alias);
    GUIObject& deref(int index){return this[index];}
    const int32_t* getRect()const{return rect;}
	void (*onUpdate)();
	void update(){if (onUpdate) onUpdate();}
	void setText(const char*);
	void draw(bool isover, const int32_t* par_rect = NULL);
	void drawAlias(bool is_text,const int32_t* par_rect = NULL);

	void startTextCapture();
	void onResize();
	void setDico(Dictionary* whay);
	template<class C> void setDico(DicoElem<C>& whay);


//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;
	ERRCODE manifest();
	ERRCODE vanish();

	double getValue()const;
	int getIntValue()const;

    //int GUImouseClick(unsigned int m_x, unsigned int m_y, unsigned int button_state);
    GUIMSG_Enum  processGUIevent(const GUIMSG_Enum event);
};
class GUIProgBar : public GUIObject{ // can be modified!
	bool isStatic() const {return(state & 1);}
	bool isSelected() const {return(state & 2);}
    int32_t* accessRect(){return rect;}

	public:
    int32_t rect[4];
	LFHPrimitive::Event<void> *task;
	unsigned int state;
    GuiTextAttribute texta;
	unsigned int styleID;
	unsigned int modstamp;

	GUIProgBar(unsigned int alias,unsigned int style_alias);
	virtual ~GUIProgBar();
	void (*onUpdate)();
	void update(){if (onUpdate) onUpdate();}
	void setText(const char*);
	void draw(bool isover, const int32_t* par_rect = NULL);
	void drawAlias(bool is_text,const int32_t* par_rect = NULL);
    const int32_t* getRect()const{return rect;}
	void onResize();
	void attachTask(LFHPrimitive::Event<void> *_task);

//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;
	ERRCODE manifest();
	ERRCODE vanish();

	double getValue()const;
};
    /*
    class GUIImageArea : public GUIObject{
        bool isStatic() const {return(state & 1);}
        bool isSelected() const {return(state & 2);}
        bool isInModification() const {return(state & 4);}
        public:
        unsigned int state;
        DataGrid<char, 3> customImage;
        GUIImageArea(unsigned int given_alias) : GUIObject(given_alias){}
        const unsigned int* getRect() const{return(rect);}
        void (*onUpdate)();
        void update(){if (onUpdate) onUpdate();}
        void draw();
//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;

        void drawstatic();
        bool loadstatic();
    };*/

class DisplayIconData{
    int texture;
	public:
	double colors[4][8];
	double UVdata[3];
	int animPeriod;
	char* text;

	DisplayIconData();
	~DisplayIconData();
	DisplayIconData& toMemmove(DisplayIconData& other);
	void setColor(int index, double hue, double bright, double sat);
	void setColorPair(int index, double hue, double bright, double sat, double hue2, double bright2, double sat2);
    void setTexture(int id);
    int getTexture()const;

	void setUVsquarre(double size, double x, double y, double border_size);
	void setUVquad(int index);
	void setUVanim16(int milliPeriod, int icon_width, int downsize_mag = 0);
	void setText(const char* text);
	void setUVuniform(GLint location)const;
};

class GUIIcon : public GUIObject{
    int32_t* accessRect(){return rect;}

	public:
    int32_t rect[4];
	unsigned int bg_texture;
    unsigned int innerect[4];
    unsigned int last_focus_time;
    unsigned int state; //(Need New Text Mesh)(Has Text Mesh)()(is Disabled)()
	unsigned int styleID;


    Guiobject_State_Enum area_state;
    GuiTextAttribute texta;
    GUIIcon(unsigned int given_alias,unsigned int style_alias);
    GUIIcon(unsigned int given_alias,unsigned int style_alias, const char* text);
    void	update(); // if the object is selected, it may receive updates
    void	draw(bool isover, const int32_t* par_rect = NULL);
    void	drawAlias(bool is_text,const int32_t* par_rect = NULL);
    void	onResize();
//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;
    ERRCODE	manifest();
    ERRCODE	vanish();
    const int32_t* getRect()const{return rect;}
    void	setEnabled(bool isEnabled);
    bool	isEnabled() const{return ((state & 8) == 0);}
    void	setText(const char*);
    double	getValue()const{return 0;}
};
class GUIGrid : public GUIObject{
    int32_t* accessRect(){return rect;}
public:
    int32_t rect[4];
	Vector<DisplayIconData> grid;
	unsigned int offset[2];
	unsigned int grid_size[2];
	unsigned int state;
    GuiTextAttribute texta;
	GUIGrid(unsigned int alias,unsigned int style_alias);

    void	update(); // if the object is selected, it may receive updates
    void	draw(bool isover, const int32_t* par_rect = NULL);
    void	drawAlias(bool is_text,const int32_t* par_rect = NULL);

	DisplayIconData& getIconData(int offset);
	void setText(int offset, const char*);
    const int32_t* getRect()const{return rect;}

	void setGridSize(int x, int y);

	void wrSlotPosize(float* fout, int slot)const;

	int getGridOffset(int mouse_x, int mouse_y)const;

    virtual ERRCODE manifest(){return 0;}
    virtual ERRCODE vanish(){return 0;}

    virtual double getValue()const{return 0.0f;}
};
class GUIAttribList : public GUIObject{
	bool isStatic() const {return(state & 1);}
	bool isSelected() const {return(state & 2);}
	bool isInModification() const {return(state & 4);}
    int32_t* accessRect(){return rect;}

	public:
    int32_t rect[4];
	unsigned int state;
	DataGrid< char*, 2 > cur_strs;
	GUIAttribList(unsigned int alias, unsigned int style_alias) : GUIObject(alias,style_alias){}
	void (*onUpdate)();
	void update(){if (onUpdate) onUpdate();}
	void draw(bool isover, const int32_t* par_rect = NULL);
//			void draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const;
	void enterValue(unsigned int col, unsigned int row, double val);
    const int32_t* getRect()const{return rect;}
	ERRCODE manifest();
	ERRCODE vanish();
	double getValue()const{return 0;}
};
class GUILog : public GUIArea{
    int32_t* accessRect(){return rect;}
public:
    int32_t rect[4];
	KeyElem< pair<unsigned char, unsigned short> , char*> log[65536];
	unsigned short indexes[256][2];
    GuiTextAttribute texta;

	//unsigned short sizes[256];
	unsigned short main_index;
	Tuple<unsigned short> channels;

	int state; // 0x80:needs mesh_update  0x100:got 65536>= messages
	uint32_t lasttime;
	GUILog(unsigned int given_alias, unsigned int style_alias);
	uint32_t colorof(unsigned int type) const;
	void insertText(unsigned char type, const char *str){this->insertText_withlength(type,str,strlen(str));}
	void insertText_withlength(unsigned char type, const char *str, unsigned int len);

    const int32_t* getRect()const{return rect;}
	void draw(bool isover, const int32_t* par_rect = NULL);
	void drawAlias(bool is_text,const int32_t* par_rect = NULL);

	void startTextCapture();
	void filterIncr();
	ERRCODE manifest();
	ERRCODE vanish();
	void onResize();
//    bool loadstatic();
    double getValue()const{return 0;}
};
class GUIDropList : public GUIObject{
    void expand();
    void select(int which);
    int32_t* accessRect(){return rect;}

    public:
    int32_t rect[4];
    GuiTextAttribute texta;
    Vector<char*> options;
    int state;
    int current_selection;
//    int tmpsize;
    GUIDropList(unsigned int given_alias, unsigned int style_alias);
    ~GUIDropList();
    GUIDropList(const GUIDropList& other)=delete;
    GUIDropList(GUIDropList&& other)=delete;
    GUIDropList& operator=(const GUIDropList& other)=delete;
    GUIDropList& operator=(GUIDropList&& other)=delete;

    void leachOptions(Vector<char*> &_options);
    const int32_t* getRect()const{return rect;}
    void draw(bool draw, const int32_t* par_rect = NULL);
    void drawAlias(bool is_text,const int32_t* par_rect = NULL);
    void onResize();
    ERRCODE	manifest();
    ERRCODE vanish();
    virtual GUIMSG_Enum  processGUIevent(const GUIMSG_Enum event);

    int getIntValue()const{return current_selection;}
    double getValue()const{return (double) current_selection;}
    void update();
};
class GUIHeatMap : public GUIObject{
    int32_t* accessRect(){return rect;}

    public:
    int32_t rect[4];

    DataGrid<float, 2u> grid;
    GLuint glbuffer_index;

    uint32_t state;
    uint16_t state_lastrow_modified;
    uint16_t selection[4];
    float selection_value;

    double scroll_rect[4];
    double color_range[2];

    float value_range[2]; // min max

    GUIHeatMap(unsigned int given_alias,unsigned int style_alias);
    ERRCODE manifest(){return 0;}
    ERRCODE vanish(){return 0;}
    double getValue()const{return 0.0f;}
    void update(){} // if the object is selected, it may receive updates

    uint32_t getDataWidth()const {return grid.dims[0];}

    void setMapDims(int _row, int _col);
    void setDataRow(int slot, const float* data);

    void setColorRange(float _min, float _max){value_range[0] = _min; value_range[1] = _max;}
    void setHorizontalRange(float _min, float _max){scroll_rect[0] = _min; scroll_rect[2] = _max - _min;}
    const int32_t* getRect()const{return rect;}
    void draw(bool isover, const int32_t* par_rect = NULL);
    void drawAlias(bool is_text,const int32_t* par_rect = NULL);

    GUIHeatMap& toMemmove(DataGrid<float, 2u> &newdata);

    Tuple<unsigned int, 2u> getMousePos() const;
    GUIMSG_Enum  processGUIevent(const GUIMSG_Enum event);
};
class GUICumulus : public GUIObject{
    public:
    GUICumulus(unsigned int given_alias,unsigned int style_alias);
    void draw(bool isover, const int32_t* par_rect = NULL);
    void drawAlias(bool is_text,const int32_t* par_rect = NULL);
    GUIMSG_Enum  processGUIevent(const GUIMSG_Enum event);
};

}; // namespace end

#ifdef _SDL_NET_H
// #include "./NetEvents.h"
#endif


#include "Display.hpp"

#endif



