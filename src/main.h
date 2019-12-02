#include "Display.h"

using namespace LFHPrimitive;

class DefaultRessourceLoader : public LFHDisplay::RessourceLoader{
    public:
    GLuint textures[10];

    DefaultRessourceLoader();
    void loadTexture(char* path, GLuint& where, int flag);
    void useTexture(const uint32_t);
    void useSound(const uint32_t);
    void useTextureExt(uint32_t);
    void allocTextureExt(uint32_t);
    void deallocTextureExt(uint32_t);
};

#define TASK_MEMBER_DEFINITIONS Taskscope(); int defstore(char* const * token, int nbtoken); void store(char* const * token, int nbtoken); void nbaddtoken(char const * const token, int& min, int& max); void help();

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

