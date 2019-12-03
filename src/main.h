#include "Display.h"

using namespace LFHPrimitive;

class DefaultRessourceLoader : public LFHDisplay::RessourceLoader{
    public:
    GLuint textures[10];

    DefaultRessourceLoader();
    void loadTexture(const char* path, GLuint& where, int flag);
    void useTexture(const uint32_t);
    void useSound(const uint32_t);
    void useTextureExt(uint32_t);
    void allocTextureExt(uint32_t);
    void deallocTextureExt(uint32_t);
};

#define TASK_MEMBER_DEFINITIONS int defstore(char* const * token, int nbtoken); void store(char* const * token, int nbtoken); void nbaddtoken(char const * const token, int& min, int& max); void help();

enum STYLS_enum{
	STYLS_DEFAULT
};

enum GUIID_enum:{
	GUIID_WINDOW=0,
	GUIID_MAIN_DROP
};

class Task : public LFHPrimitive::ArgumentParser, public LFHDisplay::ProcessState, public LFHDisplay::renderMode{
	public:
	LFHDisplay::MyWindow* dawin;
	Tuple<uint32_t, 5> data;
	TASK_MEMBER_DEFINITIONS
		
	LFHDisplay::GUIDropList dd_menu;
	
	Task(); 
        int OnKeyDown(const SDL_KeyboardEvent& Event); //OnKeyDown(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
        int OnKeyUp(const SDL_KeyboardEvent& Event); //OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
        int OnMaintain();//OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
        int listen(LFHDisplay::GUImessage&);
	void draw(LFHDisplay::MyWindow*);
	void drawAlias(LFHDisplay::MyWindow*);
};

