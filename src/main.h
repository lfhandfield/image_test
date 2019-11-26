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

using namespace LFHPrimitive;
class Task{
	public:
	Tuple<uint32_t, 5> data;
};

