#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>
#include <thread>
#include <mutex>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_net.h>
#include "bastructs.h"
#include "display.h"

using namespace LFHPrimitive;

class DefaultRessourceLoader : public LFHDisplay::RessourceLoader{
    public:
    GLuint textures[10];

    DefaultRessourceLoader();
    void loadTexture(char* path, GLuint& where, int flag);
    void useTexture(const GUITEXTURES_enum);
    void useSound(const GUISOUND_enum);
    void useTextureExt(unsigned int ID);
    void allocTextureExt(unsigned int ID);
    void deallocTextureExt(unsigned int ID);
};


class Task : public LFHPrimitive::ArgumentParser, public LFHDisplay::ProcessState, public LFHDisplay::renderMode{
	public:
	LFHDisplay::MyWindow* dawin;
	Tuple<uint32_t, 5> data;
	TASK_MEMBER_DEFINITIONS
        int OnKeyDown(const SDL_KeyboardEvent& Event); //OnKeyDown(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
        int OnKeyUp(const SDL_KeyboardEvent& Event); //OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
        int OnMaintain();//OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
        int listen(LFHDisplay::GUImessage&);
	void draw(MyWindow*);
	void drawAlias(MyWindow*);
};

