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


class Task{
	public:
	LFHDisplay::MyWindow* dawin;
	Tuple<uint32_t, 5> data;
};

