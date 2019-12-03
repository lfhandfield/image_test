/*
 * display.cpp
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



#include "./Display.h"

namespace LFHDisplay{
// LFHPrimitive::ResManager<LFHDisplay::BitmapRessource> LFHDisplay::bitmap_mana;


bool colorIndexMode = false;
bool doubleBuffered = false;
bool renderToDIB = false;

/*
LRESULT APIENTRY
WndProc(
    HWND hWnd,
    UINT message,
    WPARAM wParam,
    LPARAM lParam)
{
    switch (message) {
    case WM_CREATE:
	return 0;
    case WM_DESTROY:
	PostQuitMessage(0);
	return 0;
    case WM_SIZE:
    case WM_PALETTECHANGED:
	// Update palette mapping if this *is not* the active window.

	break;
    case WM_QUERYNEWPALETTE:
	// Update palette mapping if this *is* the active window.

	break;
    case WM_PAINT:
	// Update the window.  Don't use the device context returned by
	// BeginPaint as it won't have the right palette selected into it.
	break;
    case WM_CHAR:
	switch ((int)wParam) {
	case VK_ESCAPE:
	    DestroyWindow(hWnd);
	    return 0;
	case VK_SPACE:

	//    if (idleFunc) {
	//	idleFunc = NULL;
	//    } else {
	//	idleFunc = doRedraw;
	//    }
	default:
	    break;
	}
	break;
    default:
	break;
    }
    // Deal with any unprocessed messages
    return DefWindowProc(hWnd, message, wParam, lParam);
}

*/


InputState istate;
Controlstate ctrl_state;



unsigned int Latent::operator()(){
    state =1;
    while(state == 1) {
        internal_state = state;
        if (queue.getSync()){
            while (queue.pop_exec());
            queue.freeSync();
        }
        SDL_Delay(1024);
        }
    state = 1;
    return 0;
    }

bool Controlstate::init_openGL(){
    /*if (!glfwInit()) {
        fprintf(stderr, "ERROR: could not start GLFW3\n");
        return 1;
      }
    glewExperimental = GL_TRUE;*/
    #ifndef __MINGW32__
    /*glewExperimental=GL_TRUE;
    GLenum err = glewInit();

    if (GLEW_OK != err)
    {
        // Problem: glewInit failed, something is seriously wrong.
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));

    }*/
    #endif
    return true;
}

void Controlstate::check_openGL(){
    const GLubyte* entry;
    entry = glGetString(GL_VENDOR); printf("Vendor: %s\n", (entry) ? (const char*) entry : "Unknown");
    entry = glGetString(GL_RENDERER); printf("Renderer: %s\n", (entry) ? (const char*) entry : "Unknown");
    entry = glGetString(GL_VERSION); printf("OpenGL Version: %s\n", (entry) ? (const char*) entry : "Unknown");
    entry = glGetString(GL_SHADING_LANGUAGE_VERSION); printf("Shader Version: %s\n", (entry) ? (const char*) entry : "Unknown");
    entry = glGetString(GL_EXTENSIONS); printf("Extensions: %s\n", (entry) ? (const char*) entry : "Unknown");
}


/*
inline void
checkErr(cl_int err, const char * name){
if (err != CL_SUCCESS) {
std::cerr << "ERROR: " << name
<< " (" << err << ")" << std::endl;
exit(EXIT_FAILURE);
}
}


bool Controlstate::init_openCL(){
    cl_int err;
    cl::vector< cl::Platform > platformList;
    cl::Platform::get(&platformList);
    checkErr(platformList.size()!=0 ? CL_SUCCESS : -1, "cl::Platform::get");
    std::cerr << "Platform number is: " << platformList.size() << std::endl;std::string platformVendor;
    platformList[0].getInfo((cl_platform_info)CL_PLATFORM_VENDOR, &platformVendor);
    std::cerr << "Platform is by: " << platformVendor << "\n";
    cl_context_properties cprops[3] =
    {CL_CONTEXT_PLATFORM, (cl_context_properties)(platformList[0])(), 0};cl::Context context(
    CL_DEVICE_TYPE_CPU,
    cprops,
    NULL,
    NULL,
    &err);
    checkErr(err, "Conext::Context()");
    return true;
}*/

bool Controlstate::init_SDL(const char* const name, const char* const prod){
    ctrl_state.mouse_click_milli = 200;
    ctrl_state.ThreadID_mask[SDL_ThreadID()] = 1;
    int return_happy_valgrind=0;
    if ((return_happy_valgrind = SDL_Init(SDL_INIT_EVERYTHING)) != 0) return false;
    if ((name != NULL)&&(prod != NULL)){
        char* path = SDL_GetPrefPath(prod,name);
        ctrl_state.def_path_start = strlen(path);
        if (ctrl_state.def_path_start > 240) thrbase.terminate("path to file is too long!\n");
        memcpy(ctrl_state.def_path, path, ctrl_state.def_path_start);
        SDL_free(path);
    }
    atexit(SDL_Quit);
    return true;
}

char* Controlstate::usePrefPath(const char* const filename){
    strcpy(def_path + def_path_start, filename);
    return def_path;
}


#ifdef _SDL_NET_H

bool Controlstate::init_SDL_NET(){
    int return_happy_valgrind=0;
    if ((return_happy_valgrind =SDLNet_Init()) != 0) {
        std::cerr << "SDLNet_Init: " << SDLNet_GetError() << std::endl;
        return false;
        }
    atexit(SDLNet_Quit);
    return true;
}

#endif

void Controlstate::clean_openGL(){
}


#ifdef _SDL_MIXER_H
bool Controlstate::init_SDL_mixer(){
	printf("trying to open Audio!\n");
	//Initialize SDL_mixer with our chosen audio settings
	if(Mix_OpenAudio(44100, AUDIO_S16SYS, 2, 4096) != 0) printf("Unable to initialize audio: %s\n", Mix_GetError());
    else return true;
    return false;}
void Controlstate::clean_SDL_mixer(){
	//Need to make sure that SDL_mixer and SDL have a chance to clean up
	Mix_CloseAudio();
}
#endif

#ifdef AL_AL_H
bool Controlstate::init_openAL(){
        ctrl_state.aldev = alcOpenDevice(NULL); // open default device
        if (ctrl_state.aldev != NULL) {
            ctrl_state.alcon=alcCreateContext(ctrl_state.aldev,NULL); // create context
            if (ctrl_state.alcon != NULL) {
                alcMakeContextCurrent(ctrl_state.alcon); // set active context
            }
        }

        ALfloat listenerPos[]={0.0f,0.0f,0.0f};
        ALfloat listenerVel[]={0.0f,0.0f,0.0f};
        ALfloat listenerOri[]={0.0f,0.0f,1.0f, 0.0f,1.0f,0.0f};
    	alListenerfv(AL_POSITION,listenerPos);
    	alListenerfv(AL_VELOCITY,listenerVel);
    	alListenerfv(AL_ORIENTATION,listenerOri);
    return true;
}
void Controlstate::clean_openAL(){

}
#endif

void printInfoLog(GLhandleARB obj) {
    int infologLength = 0;
    int charsWritten  = 0;
    char *infoLog;
    /*glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infologLength);
    if (infologLength > 0){
	infoLog = (char *)malloc(infologLength);
	//glGetInfoLogARB(obj, infologLength, &charsWritten, infoLog);
	//printf("%s\n",infoLog);
	free(infoLog);
    }*/
    printf("\n");
}

GLuint Controlstate::compileshader_routine(const char * const vcode, const char * const  code, const char * const name){
    GLuint dashade = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint davshade = glCreateShader(GL_VERTEX_SHADER);
    GLchar* tmp_code; unsigned int lenght;
    lenght= strlen(vcode)+1; tmp_code = new GLchar[lenght]; memcpy(tmp_code,vcode,lenght);
    glShaderSource(davshade,1,(const GLchar**)&tmp_code,NULL);delete[](tmp_code);
	glCompileShader(davshade);
    char *infoLog;
    GLint suclength=0; glGetShaderiv(davshade,GL_COMPILE_STATUS,&suclength);
    if (suclength == GL_FALSE){
        glGetShaderiv(davshade, GL_INFO_LOG_LENGTH, &suclength);
        if (suclength){
            infoLog = new char[suclength];
            glGetShaderInfoLog(davshade,suclength,&suclength,infoLog);
            if (suclength) printf("%s vshader error: %s\n",name, infoLog);
            else printf("%s vshader error with no log\n",name);
            delete[](infoLog);
        }else printf("%s vshader error with no log alloc\n",name);
        //glDeleteShader(davshade); davshade =0;
    }

//	if (glGetError()) {printf("error for %s:",name);printInfoLog(davshade);}
    lenght = strlen(code)+1; tmp_code = new GLchar[lenght]; memcpy(tmp_code ,code,lenght);
    glShaderSource(dashade,1,(const GLchar**)&tmp_code,NULL);delete[](tmp_code);
	glCompileShader(dashade);
	    suclength=0; glGetShaderiv(dashade,GL_COMPILE_STATUS,&suclength);
    if (suclength == GL_FALSE){
        glGetShaderiv(dashade, GL_INFO_LOG_LENGTH, &suclength);
        if (suclength){
            infoLog = new char[suclength];
            glGetShaderInfoLog(dashade,suclength,&suclength,infoLog);
            if (suclength) printf("%s sshader error: %s\n",name, infoLog);
            else printf("%s sshader error with no log\n",name);
            delete[](infoLog);
        }else printf("%s sshader error with no log  alloc\n",name);
        //glDeleteShader(dashade); dashade =0;
    }
//	if (glGetError()) {printf("error for %s:",name);printInfoLog(dashade);}
    GLuint fout;
    LFH_VALGRIND_MUTE(fout = 0;)
    fout = glCreateProgram();
    glAttachShader(fout , davshade);
	if (glGetError()) {printf("Vsh error for %s:",name);printInfoLog(fout);}
    glAttachShader(fout , dashade);
    if (glGetError()) {printf("Fsh error for %s:",name);printInfoLog(fout);}
    glLinkProgram(fout);

     if (glGetError()) {printf("Lnk error for %s:",name);printInfoLog(fout);}
    else {
        printf("shader success! %i:%s", fout, name);
        glUseProgram(fout);
        if (glGetError()){
            printf("\t cannot be used! isprogram %c\n", glIsProgram(fout) ? 'Y' : 'N');

        }else printf("\n");
    }

    return(fout);
}
GLuint Controlstate::compileshader_advroutine(const char * const vcode, const char * const  code, unsigned int flag, const char * const name){
	if  (glCreateShader == NULL) {printf("glCreateShader is null!\n!!!!"); exit(1);}
	
GLuint dashade = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint davshade = glCreateShader(GL_VERTEX_SHADER);
    GLchar* tmp_code; unsigned int lenght;
    lenght= strlen(vcode)+1; tmp_code = new GLchar[lenght]; memcpy(tmp_code,vcode,lenght);
    glShaderSource(davshade,1,(const GLchar**)&tmp_code,NULL);delete[](tmp_code);glCompileShader(davshade);
    char *infoLog;
    GLint suclength=0; glGetShaderiv(davshade,GL_COMPILE_STATUS,&suclength);
    if (suclength == GL_FALSE){
        glGetShaderiv(davshade, GL_INFO_LOG_LENGTH, &suclength);
        if (suclength){
            infoLog = new char[suclength];
            glGetShaderInfoLog(davshade,suclength,&suclength,infoLog);
            if (suclength) printf("%s vshader error: %s\n",name, infoLog);
            else printf("%s vshader error with no log\n",name);
            delete[](infoLog);
        }else printf("%s vshader error with no log alloc\n",name);
        //glDeleteShader(davshade); davshade =0;
    }
    lenght = strlen(code)+1; tmp_code = new GLchar[lenght]; memcpy(tmp_code ,code,lenght);
    glShaderSource(dashade,1,(const GLchar**)&tmp_code,NULL);delete[](tmp_code);glCompileShader(dashade);
    suclength=0; glGetShaderiv(dashade,GL_COMPILE_STATUS,&suclength);
    if (suclength == GL_FALSE){
        glGetShaderiv(dashade, GL_INFO_LOG_LENGTH, &suclength);
        if (suclength){
            infoLog = new char[suclength];
            glGetShaderInfoLog(dashade,suclength,&suclength,infoLog);
            if (suclength) printf("%s sshader error: %s\n",name, infoLog);
            else printf("%s sshader error with no log\n",name);
            delete[](infoLog);
        }else printf("%s sshader error with no log  alloc\n",name);
        //glDeleteShader(dashade); dashade =0;
    }
    //if ((dashade == 0)||(davshade == 0)) return 0;
    GLuint fout;
    LFH_VALGRIND_MUTE(fout = 0;)
    fout = glCreateProgram();
    glAttachShader(fout, davshade);
//	if (glGetError()) {printf("error for %s:",name);printInfoLog(fout);}
    glAttachShader(fout, dashade);
	if  (glBindAttribLocation == NULL) {printf("glBindAttribLocation is null!\n!!!!"); exit(1);}
 //   if (glGetError()) {printf("error for %s:",name);printInfoLog(fout);}
    glBindAttribLocation(fout,ATTRIBUTE_POSITION, "glcVertex");//printInfoLog(fout);
    if (flag & 1){
        glBindAttribLocation(fout,ATTRIBUTE_NORMAL, "glcNormal");//printInfoLog(fout);
        glBindAttribLocation(fout,ATTRIBUTE_TEXTURE_COORDINATES, "glcMultiTexCoord0");//printInfoLog(fout);
    }
    if (flag & 2){
        glBindAttribLocation(fout,ATTRIBURE_BONEID, "boneId");//printInfoLog(fout);
        glBindAttribLocation(fout,ATTRIBURE_BONEWEIGHT, "boneWeight");//printInfoLog(fout);
    }
    if (flag & 4){
        glBindAttribLocation(fout,ATTRIBURE_CHARID, "chr_id");//printInfoLog(fout);
    }
    if (flag & 8){
        glBindAttribLocation(fout,ATTRIBURE_BONEID, "boneData");//printInfoLog(fout);
    }
    if (flag & 16){
        glBindAttribLocation(fout,ATTRIBUTE_POSITION_ALT, "glxVertex");
        glBindAttribLocation(fout,ATTRIBUTE_NORMAL_ALT, "glxNormal");
    }

    glLinkProgram(fout);
//	if (glGetError()) {printf("error for %s:",name);printInfoLog(fout);}
    /*if (flag & 1) printf("got location: %i\t%i (for n and tx)\n",glGetAttribLocation(fout, "glcNormal")
                                      ,glGetAttribLocation(fout, "glcMultiTexCoord0"));
    if (flag & 2) printf("got location: %i\t%i (for bones)\n",glGetAttribLocation(fout, "boneId")
                                      ,glGetAttribLocation(fout, "boneWeight"));
    if (flag & 4) printf("got location: %i\n",glGetAttribLocation(fout, "chr_id"));
    printf("got location: %i (vertexes)\n",glGetAttribLocation(fout, "glcVertex"));*/
    //printInfoLog(fout);
    GLint status[2];
    if (glGetError()) {printf("L error for %s:",name);printInfoLog(fout);}
    else {
        glUseProgram(fout);
        if (glGetError()){
            printf("Advanced shader %s compiled but cannot be used! isprogram %c\n",name, glIsProgram(fout) ? 'Y' : 'N');
        }else printf("Advanced shader %s compiled!\n",name);
    }

    return(fout);
}

void Controlstate::compiledefaultshaders(){
   LFH_ALIVE;

     ctrl_state.datext_shader = ctrl_state.compileshader_advroutine( SHADER_STRINGIFY(130,
    in vec4 glcVertex;
    in vec4 chr_id;

 //   flat out vec4 font_color;
    uniform vec2 pixscale;
    uniform vec2 vPos;
	//uniform vec2 tex_offset;
    uniform int tex_height;
    uniform sampler2D tFontcoor;
    uniform vec4 def_color;
    uniform vec4 bgcolor_base;
    uniform float bgcolor_factor;
    out vec2 texcoor;
    out float channel;
    out vec4 fcolor;
    out vec4 bcolor;
    out vec4 bgcolor;
    void main(void){
        //vec2 tmpp = // ivec2( int(chr_id.s * 255),int(chr_id.t * 255));
        vec4 fcoor = texture(tFontcoor, chr_id.xy * (255.0f/256.0f));
        //texcoor = tmp.xy + ivec2(int(chr_id.z)* tmp.z -(tmp.z>>1), tex_height*int(chr_id.z));
       // texcoor = vec2(tmpfff.x +  chr_id.z* (tmpfff.z -1/ 256.0f), chr_id.w) ;

        vec2 vdir = vec2( ((gl_VertexID & 2)!= 0) ? 0.0f : 1.0f ,  ((gl_VertexID ^ (gl_VertexID >> 1)) & 1) == 1 ? 0.0f : 1.0f);
        texcoor = vec2(fcoor.x + (vdir.x * 10 - 1.0) /512.0f ,(20*floor(fcoor.y * 31.875)+ vdir.y *  19.0 +1.0f)/512.0f   );
        channel = exp2(15.0f - fract(fcoor.y * 31.875)*8.0f) / 255.0f;

        //texcoor = vec2((float(tmpfff.x) + chr_id.z * 16.0f ) /256.0f ,(float(tmpfff.y & 31)+(1.0f-chr_id.w)) * tex_height );
        // texcoor = vec2(float(tmpfff.x)/256.0f + (chr_id.z * (float(tmpfff.z & 31)/256.0f) ) , float((tmpfff.y & 31) *tex_height) /256.0f +(1.0f-chr_id.w) * tex_height) ; //* tmpfff.x
        gl_Position = vec4((vPos.x + ((glcVertex.x * 0.5 + vdir.x) * (tex_height>>1) )) * pixscale.x -1.0f,
                           1.0f - ( (vPos.y + glcVertex.y * (tex_height >>1) + vdir.y *float(tex_height-1))) * pixscale.y,-1,1);
//        gl_Position = vec4(vdir.x,vdir.y,-1,1);

        fcolor = (chr_id.w < 0.1f) ? def_color : vec4( fract(chr_id.z * 15.9375) * 1.06666f ,chr_id.z* 1.06666f,fract(chr_id.w * 15.9375)* 1.06666f,chr_id.w * 1.06666f);
        bcolor = bgcolor_base + (fcolor * bgcolor_factor);
        bgcolor = ((gl_VertexID & 4)==0) ? vec4(0.5f,0.5f,0.5f,0.5f) : vec4(0.0f,0.0f,0.0f,0.0f);
        } // (1-z) * f = p-b
    ), SHADER_STRINGIFY(130,
  //  \#version 130
    uniform sampler2D tFontdata;
    uniform sampler2D tPalette;
    uniform vec4 boundrect;
    in float channel;
    in vec4 fcolor;
    in vec4 bcolor;
    in vec4 bgcolor;
  //  flat in vec4 font_color;
    in vec2 texcoor;
        void main(){
            vec4 tmp = texture(tFontdata, texcoor);
            if ((gl_FragCoord.x < boundrect.x)
               ||(gl_FragCoord.y < boundrect.y)
               ||(gl_FragCoord.x > boundrect.z)
               ||(gl_FragCoord.y > boundrect.w)) discard;
            gl_FragColor = (fract(tmp.x * channel) < 0.5f) ? bgcolor : ((fract(tmp.y * channel) < 0.5f) ? fcolor: bcolor);
            }
    ),5, "datext_shader");
   LFH_ALIVE;
    glUseProgram(datext_shader);
    glUniform1i(glGetUniformLocation(datext_shader, "tPalette"),2);
    glUniform1i(glGetUniformLocation(datext_shader, "tFontdata"),0);
    glUniform1i(glGetUniformLocation(datext_shader, "tFontcoor"),1);
   LFH_ALIVE;
ctrl_state.datext2_shader = ctrl_state.compileshader_advroutine( SHADER_STRINGIFY(130,
    in vec4 glcVertex;
    in vec4 chr_id;

 //   flat out vec4 font_color;
    uniform vec2 pixscale;
	//uniform vec2 tex_offset;
    uniform int tex_height;
    uniform sampler2D tFontcoor;
    uniform vec4 def_color;
    uniform vec4 bgcolor_base;
    uniform vec4 uPosition;
    uniform float bgcolor_factor;
    out vec2 texcoor;
    out float channel;
    out vec4 fcolor;
    out vec4 bcolor;
    out vec4 bgcolor;
    void main(void){
        vec4 fcoor = texture(tFontcoor, chr_id.xy * (255.0f/256.0f));
        vec2 vdir = vec2( ((gl_VertexID & 2)!= 0) ? 0.0f : 1.0f ,  ((gl_VertexID ^ (gl_VertexID >> 1)) & 1) == 1 ? 0.0f : 1.0f);
        texcoor = vec2(fcoor.x + (vdir.x * 10 - 1.0) /512.0f ,(20*floor(fcoor.y * 31.875)+ vdir.y * 19.0f +1.0f)/512.0f   );
        channel = exp2(15.0f - fract(fcoor.y * 31.875)*8.0f) / 255.0f;
        vec4 vRelative = gl_ProjectionMatrix *gl_ModelViewMatrix * uPosition;
        gl_Position = vec4((vRelative.x / vRelative.w) + ((glcVertex.x * 0.5 + vdir.x) * 10 * uPosition.z ) * pixscale.x,
                           (vRelative.y / vRelative.w) - ( glcVertex.y * 10 + vdir.y * 19  ) * pixscale.y * uPosition.z,-1,1);


        fcolor = (chr_id.w < 0.1f) ? def_color : vec4( fract(chr_id.z * 15.9375) * 1.06666f ,chr_id.z* 1.06666f,fract(chr_id.w * 15.9375)* 1.06666f,chr_id.w * 1.06666f);
        bcolor = bgcolor_base + (fcolor * bgcolor_factor);
        bgcolor = ((gl_VertexID & 4)==0) ? vec4(0.5f,0.5f,0.5f,0.5f) : vec4(0.0f,0.0f,0.0f,0.0f);
        } // (1-z) * f = p-b
    ), SHADER_STRINGIFY(130,
  //  \#version 130
    uniform sampler2D tFontdata;
    uniform sampler2D tPalette;
    uniform vec4 boundrect;
    in float channel;
    in vec4 fcolor;
    in vec4 bcolor;
    in vec4 bgcolor;
  //  flat in vec4 font_color;
    in vec2 texcoor;
        void main(){
            vec4 tmp = texture(tFontdata, texcoor);
        /*    if ((gl_FragCoord.x < boundrect.x)
               ||(gl_FragCoord.y < boundrect.y)
               ||(gl_FragCoord.x > boundrect.z)
               ||(gl_FragCoord.y > boundrect.w)) discard;*/
            gl_FragColor = (fract(tmp.x * channel) < 0.5f) ? bgcolor : ((fract(tmp.y * channel) < 0.5f) ? fcolor: bcolor);
            }
    ),5, "datext2_shader");
   LFH_ALIVE;
    glUseProgram(datext2_shader);
    glUniform1i(glGetUniformLocation(datext2_shader, "tPalette"),2);
    glUniform1i(glGetUniformLocation(datext2_shader, "tFontdata"),0);
    glUniform1i(glGetUniformLocation(datext2_shader, "tFontcoor"),1);
   LFH_ALIVE;



    ctrl_state.sstext_shader = ctrl_state.compileshader_advroutine( SHADER_STRINGIFY(130,
    in vec4 glcVertex;
    in vec4 chr_id;
    //varying vec2 texcoor;
    //varying float channel;
    uniform vec2 pixscale;
    uniform vec2 vPos;
    uniform vec3 id_color;
    uniform int tex_height;
    uniform sampler2D tFontcoor;
    out vec4 icolor;
    out vec3 acolor;
    out vec2 texcoor;
    out float channel;
    void main(void){
        vec4 fcoor = texture(tFontcoor, chr_id.xy * (255.0f/256.0f));
        vec2 vdir = vec2( ((gl_VertexID & 2)!= 0) ? 0.0f : 1.0f ,  ((gl_VertexID ^ (gl_VertexID >> 1)) & 1) == 1 ? 0.0f : 1.0f);
        texcoor = vec2(fcoor.x + (vdir.x * 10 - 1.0) /512.0f ,(20*floor(fcoor.y * 31.875)+ vdir.y * 19.0f +1.0f)/512.0f   );
        channel = exp2(15.0f - fract(fcoor.y * 31.875)*8.0f) / 255.0f;
        gl_Position = vec4((vPos.x + ((glcVertex.x * 0.5 + vdir.x) * (tex_height>>1) )) * pixscale.x -1.0f,
                           1.0f - ( (vPos.y + glcVertex.y * (tex_height >>1) + vdir.y *float(tex_height-1))) * pixscale.y,-1,1);
        icolor = vec4(id_color.x,id_color.y,id_color.z, (1.0 / 255) * (gl_VertexID >> 2));
        acolor.y = 1.0;
        acolor.z = 1.0;
        acolor.x = ((gl_VertexID & 4)==0) ? 0.5f : 0.0f;
        } // (1-z) * f = p-b
    ), SHADER_STRINGIFY(130,
  //  \#version 130
    uniform sampler2D tFontdata;
    uniform sampler2D tPalette;
    uniform vec4 boundrect;
    in float channel;
    in vec4 icolor;
    in vec3 acolor;
    in vec2 texcoor;
        void main(){
            vec4 tmp = texture(tFontdata, texcoor);
            if ((gl_FragCoord.x < boundrect.x)
               ||(gl_FragCoord.y < boundrect.y)
               ||(gl_FragCoord.x > boundrect.z)
               ||(gl_FragCoord.y > boundrect.w)) discard;
            if (((fract(tmp.x * channel) < 0.5f) ? acolor.x : ((fract(tmp.y * channel) < 0.5f) ? acolor.y : acolor.z)) == 0) discard;
            gl_FragColor = icolor;
            }
    ),5, "sstext_shader");
   LFH_ALIVE;
    glUseProgram(sstext_shader);
    glUniform1i(glGetUniformLocation(sstext_shader, "tPalette"),2);
    glUniform1i(glGetUniformLocation(sstext_shader, "tFontdata"),0);
    glUniform1i(glGetUniformLocation(sstext_shader, "tFontcoor"),1);
   LFH_ALIVE;

    ctrl_state.sstext2_shader = ctrl_state.compileshader_advroutine( SHADER_STRINGIFY(130,
    in vec4 glcVertex;
    in vec4 chr_id;
    uniform vec2 pixscale;
    uniform int tex_height;
    uniform sampler2D tFontcoor;
    uniform vec4 def_color;
    uniform vec4 bgcolor_base;
    uniform vec4 uPosition;
    uniform float bgcolor_factor;
    uniform vec3 id_color;
    out vec2 texcoor;
    out float channel;
    out vec4 icolor;
    out vec3 acolor;
    void main(void){
        vec4 fcoor = texture(tFontcoor, chr_id.xy * (255.0f/256.0f));
        vec2 vdir = vec2( ((gl_VertexID & 2)!= 0) ? 0.0f : 1.0f ,  ((gl_VertexID ^ (gl_VertexID >> 1)) & 1) == 1 ? 0.0f : 1.0f);
        texcoor = vec2(fcoor.x + (vdir.x * 10 - 1.0) /512.0f ,(20*floor(fcoor.y * 31.875)+ vdir.y * 19.0f +1.0f)/512.0f   );
        channel = exp2(15.0f - fract(fcoor.y * 31.875)*8.0f) / 255.0f;
        vec4 vRelative = gl_ProjectionMatrix *gl_ModelViewMatrix * uPosition;
        gl_Position = vec4((vRelative.x / vRelative.w) + ((glcVertex.x * 0.5 + vdir.x) * 10 * uPosition.z *99 ) * pixscale.x,
                           (vRelative.y / vRelative.w) - ( glcVertex.y * 10 + vdir.y * 19  ) * pixscale.y * uPosition.z*99,-1,1);
        icolor = vec4(id_color.x,id_color.y,id_color.z, (1.0 / 255) * (gl_VertexID >> 2));
        acolor.y = 1.0;
        acolor.z = 1.0;
        acolor.x = ((gl_VertexID & 4)==0) ? 0.5f : 0.0f;
        } // (1-z) * f = p-b
    ), SHADER_STRINGIFY(130,
  //  \#version 130
    uniform sampler2D tFontdata;
    uniform sampler2D tPalette;
    uniform vec4 boundrect;
    in float channel;
    in vec4 icolor;
    in vec3 acolor;
    in vec2 texcoor;
        void main(){
            vec4 tmp = texture(tFontdata, texcoor);
            //if (((fract(tmp.x * channel) < 0.5f) ? acolor.x : ((fract(tmp.y * channel) < 0.5f) ? acolor.y : acolor.z)) == 0) discard;
            gl_FragColor = icolor;
            }
    ),5, "sstext2_shader");
   LFH_ALIVE;
    glUseProgram(sstext2_shader);
    glUniform1i(glGetUniformLocation(sstext2_shader, "tPalette"),2);
    glUniform1i(glGetUniformLocation(sstext2_shader, "tFontdata"),0);
    glUniform1i(glGetUniformLocation(sstext2_shader, "tFontcoor"),1);
   LFH_ALIVE;


/*
 ctrl_state.datext_shader = ctrl_state.compileshader_routine( STRINGIFY(
    \#version 130
//	varying vec2 texcoor;
  //  varying out float channel;
  //  uniform vec4 pixscale;
  //  uniform int tex_height;
  //  uniform sampler2D tex_coors;
    void main(void){
        //vec2 tmpp = // ivec2( int(chr_id.s * 255),int(chr_id.t * 255));
     //   vec4 fcoor = texture2D(tex_coors, gl_MultiTexCoord0.xy * (511.0f/512.0f));
        //texcoor = tmp.xy + ivec2(int(chr_id.z)* tmp.z -(tmp.z>>1), tex_height*int(chr_id.z));
       // texcoor = vec2(tmpfff.x +  chr_id.z* (tmpfff.z -1/ 256.0f), chr_id.w) ;

       // vec2 vdir = vec2( ((gl_VertexID & 2)!= 0) ? 0.0f : 1.0f ,  ((gl_VertexID ^ (gl_VertexID >> 1)) & 1) == 1 ? 0.0f : 1.0f);
    //    vec2 vdir = gl_Vertex.zw;

  //      texcoor = vec2(fcoor.x + (vdir.x * (tex_height >>1) ) /512.0f ,(tex_height*floor(fcoor.y * 31.875)+ vdir.y * float(tex_height-1) +1.0f)/512.0f   );

        //channel = 1 << ((tmpfff>> 24) & 31);
     //   channel = chansel * 256.0f / 255.0f;
//        channel = exp2(15.0f - fract(fcoor.y * 31.875)*8.0f) / 255.0f;

        //texcoor = vec2((float(tmpfff.x) + chr_id.z * 16.0f ) /256.0f ,(float(tmpfff.y & 31)+(1.0f-chr_id.w)) * tex_height );
        // texcoor = vec2(float(tmpfff.x)/256.0f + (chr_id.z * (float(tmpfff.z & 31)/256.0f) ) , float((tmpfff.y & 31) *tex_height) /256.0f +(1.0f-chr_id.w) * tex_height) ; /*//* tmpfff.x
        //gl_Position = vec4((pixscale.z + gl_Vertex.x * (tex_height >>1)+ vdir.x* (tex_height >>1) ) * pixscale.x -1.0f,1.0f - ( pixscale.w + gl_Vertex.y * (tex_height >>1) + vdir.y *float(tex_height-1)) * pixscale.y,-1,1);
        gl_Position = vec4(0.0,0.0,-1.0,1.0);

     //   color = vec4( fract(gl_MultiTexCoord0.z * 15.9375) * 1.06666f ,gl_MultiTexCoord0.z* 1.06666f,fract(gl_MultiTexCoord0.w * 15.9375)* 1.06666f,gl_MultiTexCoord0.w * 1.06666f);
        } // (1-z) * f = p-b
    ), STRINGIFY(
    \#version 130
 //   uniform sampler2D tex_font;
  //  uniform sampler2D tex_color;
 //   uniform vec4 boundrect;
 //   in float channel;
 //   in vec2 texcoor;
        void main(){
          //  vec4 tmp = texture(tex_font, texcoor);
        //    if ((gl_FragCoord.x < boundrect.x)
        //       ||(gl_FragCoord.y < boundrect.y)
        //       ||(gl_FragCoord.x > boundrect.z)
        //       ||(gl_FragCoord.y > boundrect.w)) discard;
          //  if (fract(tmp.x * channel) < 0.5f) discard;

           // vec4 col = texture(tex_color
         //   gl_FragColor = (fract(tmp.y * channel) < 0.5f) ? vec4(1.0,0.7,0.8,1.0): vec4(1.0,0.0,0.2,1.0);
            //gl_FragColor = vec4(float(tmp.x)/4,float(tmp.y)/4,float(tmp.z)/4,1);
            gl_FragColor = vec4(1.0,1.0,1.0,1.0);
            }
    ),"datext_shader");*/

    // palette rules
    // if y != z and w != 1.0 : palette color lookup

	// uses 2 textures, foreground and backgroupd,
	// foreground manages transparency and overwrites/deforms background
    glGetError();
    ctrl_state.daframe_shader = ctrl_state.compileshader_routine( SHADER_STRINGIFY(130,
    uniform vec4 boundrect; // min_x min_y max_x max_y
	uniform vec4 boundtexcoor; // min_u min_v max_u max_v
    uniform vec2 windowsize;
    out vec2 texcoor;
    void main(void){
        texcoor = vec2( (gl_Vertex.x < 0.5f) ? boundtexcoor.x : boundtexcoor.z ,  (gl_Vertex.y < 0.5f) ? boundtexcoor.y : boundtexcoor.w);
        //vec2 vdir = vec2( (gl_Vertex.x < 0.5f) ? 0.0f : 1.0f ,  (gl_Vertex.y < 0.5f) ? 0.0f : 1.0f);
        gl_Position = vec4( -1.0f + 2.0f * ( (gl_Vertex.x < 0.5f) ? boundrect.x : boundrect.z) / windowsize.x
                          , -1.0f + 2.0f *( (gl_Vertex.y < 0.5f) ? boundrect.w : boundrect.y ) / windowsize.y ,-1,1);
        }
    ), SHADER_STRINGIFY(130,
    uniform sampler2D tBackground;
    uniform sampler2D tPalette;
    uniform sampler2D tBorder;

    uniform vec4 boundrect;
    uniform vec2 borderUV;
    uniform float timeloop;
    uniform vec4 maskRBG;

    uniform mat3 color; // HBGx2, Trans, Spec,Emit,
    uniform mat3 color2;
    uniform mat3 color3;
    uniform mat3 color4;
    in vec2 texcoor;
        void main(){
            vec2 lt;
            lt.x = (texcoor.x < borderUV.x) ? 0.25 * texcoor.x / borderUV.x : ((texcoor.x > 1.0f - borderUV.x) ? 1.0 - 0.25 * (1.0f - texcoor.x) / borderUV.x : 0.25 + 0.5 * fract(0.015625 * gl_FragCoord.x));
            lt.y = (texcoor.y < borderUV.y) ? 0.25 * texcoor.y / borderUV.y : ((texcoor.y > 1.0f - borderUV.y) ? 1.0 - 0.25 * (1.0f - texcoor.y) / borderUV.y : 0.25 + 0.5 * fract(0.015625 * gl_FragCoord.y));

            vec4 colselect = texture(tBorder,lt);
			vec4 colmask = texture(tBackground,lt);
			if ((abs(lt.x -0.5) > 0.375f)||(abs(lt.y -0.5) > 0.375f)) colmask = vec4(0,0,0,0);
			float grayfact;
			vec4 colA;
			vec4 colB;
			float modul = colselect.x * 4.047619047619047619047619047619f;
            vec4 uni_col;
            colA.w = 1.0f;
            colB.w = 1.0f;
			if (colselect.x < 255.0f / 510.0f){
				if (colselect.x < 127.0f / 510.0f){
					colA.xyz = color[0];
					colB.xyz = color[0];
					colA.w = 0.0f;
				}else{
					colA.xyz = color2[0];
					colB.xyz = color2[1];
                    modul -= 64.0f / 63.0f;
				}
			}else{
				if (colselect.x < 383.0f / 510.0f){
					colA.xyz = color3[0];
					colB.xyz = color3[1];
					modul -= 128.0f / 63.0f;
				}else{
					colA.xyz = color4[0];
					colB.xyz = color4[1];
					modul -= 192.0f / 63.0f;
				}
			}
			uni_col = mix(colA, colB, 0.5f - 0.5f*mix(cos(6.283185307179586476925286766559f * colselect.z), cos(6.283185307179586476925286766559f *(colselect.z+timeloop)) ,modul) );
			uni_col.y = ((colselect.y < 0.5f) ? uni_col.y * colselect.y: 0.5 - (1.0 - uni_col.y) * (1-colselect.y)) * 2.0;
			if (colselect.w >= 0.75) grayfact = 16 * (1.0 - colselect.w) * (1.0 - colselect.w);
			else grayfact = ((colselect.w < 0.5f) ? uni_col.z * colselect.w : 0.5 - (1.0 - uni_col.z) * (0.75-colselect.w) * 2.0) * 2.0; // * ( (colselect.w >= 0.75) ?  uni_col.z * 20 * (0.8 - colselect.w) : 1.0f);
			vec4 palcol = texture(tPalette, vec2(uni_col.x, (30.0 * uni_col.y + 1.0) / 128.0 ));
			palcol.xyz = mix(vec3(uni_col.y,uni_col.y,uni_col.y) ,colselect.w < 0.75 ? palcol.xyz : colselect.xyz ,(1.0-grayfact) * uni_col.z);
			palcol.xyz = mix(colmask.xyz, palcol.xyz ,uni_col.w);
			palcol.w = uni_col.w + (1.0 - uni_col.w) * colmask.w;
			palcol.xyz = mix(palcol.xyz, maskRBG.xyz, maskRBG.w);
			if (palcol.w < 0.125f) discard;
			gl_FragColor = palcol;
			/*gl_FragColor = vec4(0,0,0,0);*/
            }
    ),"daframe_shader");

    glUseProgram(daframe_shader);
    glUniform1i(glGetUniformLocation(daframe_shader, "tBackground"),2);
    glUniform1i(glGetUniformLocation(daframe_shader, "tBorder"),0);
    glUniform1i(glGetUniformLocation(daframe_shader, "tPalette"),1);
   LFH_ALIVE;

//            lt.x = (gl_FragCoord.x < boundrect.x + 32* bordersize) ? 0.5 + 0.0078125 * (gl_FragCoord.x - boundrect.x) / bordersize : ((gl_FragCoord.x + 32* bordersize> boundrect.z) ? 1.0 - 0.0078125 * (boundrect.z - gl_FragCoord.x) / bordersize:  0.5 * fract(0.015625 * gl_FragCoord.x));

    ctrl_state.daframe_alias_shader = ctrl_state.compileshader_routine( SHADER_STRINGIFY(130,
//\#version 130
    uniform vec4 boundrect;
	uniform vec4 boundtexcoor; // min_u min_v max_u max_v
    uniform vec2 windowsize;
    out vec2 texcoor;
    void main(void){

        texcoor = vec2( (gl_Vertex.x < 0.5f) ? boundtexcoor.x : boundtexcoor.z ,  (gl_Vertex.y < 0.5f) ? boundtexcoor.y : boundtexcoor.w);
        vec2 vdir = vec2( (gl_Vertex.x < 0.5f) ? 0.0f : 1.0f ,  (gl_Vertex.y < 0.5f) ? 0.0f : 1.0f);
        gl_Position = vec4( -1.0f + 2.0f * ( (gl_Vertex.x < 0.5f) ? boundrect.x : boundrect.z) / windowsize.x
                          , -1.0f + 2.0f * ( (gl_Vertex.y < 0.5f) ? boundrect.w : boundrect.y) / windowsize.y ,-1,1);
        } // (1-z) * f = p-b
    ), SHADER_STRINGIFY(130,
//\#version 130
    //uniform isampler2D tex_font;
    //uniform sampler2D tex_color;
    uniform sampler2D tBorder;
    uniform vec4 boundrect;
    uniform vec2 borderUV;

    uniform vec2 color;
    uniform vec2 color2;
    uniform vec2 color3;
    uniform vec2 color4;


    uniform vec4 id_color;
    in vec2 texcoor;
        void main(){
            //int tmp = texture(tex_font, texcoor).r;
            vec2 lt;

            lt.x = (texcoor.x < borderUV.x) ? 0.25 * texcoor.x / borderUV.x : ((texcoor.x > 1.0f - borderUV.x) ? 1.0 - 0.25 * (1.0f - texcoor.x) / borderUV.x : 0.25 + 0.5 * fract(0.015625 * gl_FragCoord.x));
            lt.y = (texcoor.y < borderUV.y) ? 0.25 * texcoor.y / borderUV.y : ((texcoor.y > 1.0f - borderUV.y) ? 1.0 - 0.25 * (1.0f - texcoor.y) / borderUV.y : 0.25 + 0.5 * fract(0.015625 * gl_FragCoord.y));

            vec4 colselect = texture(tBorder,lt);
            vec2 trans;
            float modul = colselect.x * 4.047619047619047619047619047619f;

            if (colselect.x < 255.0f / 510.0f){
				if (colselect.x < 127.0f / 510.0f){
					trans = color;
				}else{
					trans = color2;
                    modul -= 64.0f / 63.0f;
				}
			}else{
				if (colselect.x < 383.0f / 510.0f){
					trans = color3;
					modul -= 128.0f / 63.0f;
				}else{
					trans = color4;
					modul -= 192.0f / 63.0f;
				}
			}
            float seltr = mix(trans.x, trans.y, 0.5f - 0.5f*mix(cos(6.283185307179586476925286766559f * colselect.z), cos(6.283185307179586476925286766559f *(colselect.z)) ,modul) );

            if (seltr < 0.125) discard;
            gl_FragColor = id_color;
            }
    ),"daframe_alias_shader");


    glUseProgram(daframe_alias_shader);
    glUniform1i(glGetUniformLocation(daframe_alias_shader, "tBorder"),0);


	// (color-channel,bright-grayness,time-stamp,aura)
	// if color-channel < 0.25, uniform color lookup if performed.

	// static colors, dynamic colors


    ctrl_state.daicon_shader = ctrl_state.compileshader_routine( SHADER_STRINGIFY(130,
//\#version 130
    uniform vec4 boundrect;
    uniform vec4 tex_rectUV;
    uniform vec4 mask_rectUV;
    uniform vec2 windowsize;
    out vec2 texUV;
    out vec2 maskUV;
    void main(void){
        vec2 vdir = vec2( (gl_Vertex.x < 0.5f) ? 0.0f : 1.0f ,  (gl_Vertex.y < 0.5f) ? 0.0f : 1.0f);
        gl_Position = vec4( -1.0f + 2.0f * (0.5f + ((gl_Vertex.x < 0.5f) ? boundrect.x : boundrect.z)) / windowsize.x
                          , -1.0f + 2.0f * (0.5f + ((gl_Vertex.y < 0.5f) ? boundrect.w : boundrect.y)) / windowsize.y ,1,1);
        texUV = vec2( (gl_Vertex.x < 0.5f) ? tex_rectUV.x : tex_rectUV.z, (gl_Vertex.y < 0.5f) ? tex_rectUV.w : tex_rectUV.y);
        maskUV = vec2( (gl_Vertex.x < 0.5f) ? mask_rectUV.x : mask_rectUV.z, (gl_Vertex.y < 0.5f) ? mask_rectUV.w : mask_rectUV.y);
        }
    ), SHADER_STRINGIFY(130,
//\#version 130
    uniform vec4 boundrect;
    uniform sampler2D tMotif;
    uniform sampler2D tMask;
    uniform sampler2D tPalette;
    uniform float timeloop;
    uniform vec4 maskRBG;
    uniform vec4 color; // hue, bright, gray
    uniform vec4 color2;
    uniform vec4 color3;
    uniform vec4 color4;
    uniform vec4 color5;
    uniform vec4 color6;
    uniform vec4 color7;
    in vec2 texUV;
    in vec2 maskUV;
        void main(){
            vec4 colselect = texture(tMotif, texUV);
			vec4 colmask = texture(tMask, maskUV);
			float grayfact;
			vec4 colA;
			vec4 colB;
			float modul = colselect.x * 4.047619047619047619047619047619f;
			if (colselect.x < 255.0f / 510.0f){
				if (colselect.x < 127.0f / 510.0f){
					colA = vec4(color.x,color.y,color.z,0.0f);
					colB = color;
				}else{
					colA = color2;
					colB = color3;
					modul -= 64.0f / 63.0f;
				}
			}else{
				if (colselect.x < 383.0f / 510.0f){
					colA = color4;
					colB = color5;
					modul -= 128.0f / 63.0f;
				}else{
					colA = color6;
					colB = color7;
					modul -= 192.0f / 63.0f;
				}
			}
			vec4 uni_col = mix(colA, colB, 0.5f - 0.5f*mix(cos(6.283185307179586476925286766559f * colselect.z), cos(6.283185307179586476925286766559f *(colselect.z+timeloop)) ,modul) );

			//uni_col = backcolor;
			//if ((color.x > 0.5)&&(colselect.x < 63.0f / 510.0f)) uni_col = backcolor.xyz;
			//else{
				uni_col.y = ((colselect.y < 0.5f) ? uni_col.y * colselect.y: 0.5 - (1.0 - uni_col.y) * (1-colselect.y)) * 2.0;
			//	uni_col.z = ((colselect.w < 0.5f) ? uni_col.z * colselect.w: 0.5 - (1.0 - uni_col.z) * (1-colselect.z)) * 2.0;
			//}
			if (colselect.w >= 0.75) grayfact = 16 * (1.0 - colselect.w) * (1.0 - colselect.w);
			else grayfact = ((colselect.w < 0.5f) ? uni_col.z * colselect.w : 0.5 - (1.0 - uni_col.z) * (0.75-colselect.w) * 2.0) * 2.0; // * ( (colselect.w >= 0.75) ?  uni_col.z * 20 * (0.8 - colselect.w) : 1.0f);

			vec4 palcol = texture(tPalette, vec2(uni_col.x, (30 * uni_col.y + 1) / 128 ));
			palcol.xyz = mix(vec3(uni_col.y,uni_col.y,uni_col.y) ,colselect.w < 0.75 ? palcol.xyz : colselect.xyz ,(1.0-grayfact) * uni_col.z);
			palcol.w = uni_col.w * colmask.w;
			palcol.xyz = mix(palcol.xyz, maskRBG.xyz, maskRBG.w);
			//if (colselect.w <= 0.95)

			//palcol = vec4(colselect.w,colselect.w,colselect.w,1.0);
			gl_FragColor = palcol;
            }
    ),"daicon_shader");


    glUseProgram(daicon_shader);
    glUniform1i(glGetUniformLocation(daicon_shader, "tMotif"),0);
    glUniform1i(glGetUniformLocation(daicon_shader, "tMask"),1);
    glUniform1i(glGetUniformLocation(daicon_shader, "tPalette"),2);

	daicon_alias_shader = ctrl_state.compileshader_routine( SHADER_STRINGIFY(130,
//\#version 130
    uniform vec4 boundrect;
    uniform vec2 windowsize;
    uniform vec4 tex_rectUV;
    uniform vec4 mask_rectUV;
    out vec2 texUV;
    out vec2 maskUV;

    void main(void){
        vec2 vdir = vec2( (gl_Vertex.x < 0.5f) ? 0.0f : 1.0f ,  (gl_Vertex.y < 0.5f) ? 0.0f : 1.0f);
        gl_Position = vec4( -1.0f + 2.0f * ((gl_Vertex.x < 0.5f) ? boundrect.x : boundrect.z) / windowsize.x
                          , -1.0f + 2.0f * ((gl_Vertex.y < 0.5f) ? boundrect.w : boundrect.y ) / windowsize.y ,-1,1);
        texUV = vec2( (gl_Vertex.x < 0.5f) ? tex_rectUV.x : tex_rectUV.z, (gl_Vertex.y < 0.5f) ? tex_rectUV.w : tex_rectUV.y);
        maskUV = vec2( (gl_Vertex.x < 0.5f)? mask_rectUV.x : mask_rectUV.z, (gl_Vertex.y < 0.5f) ? mask_rectUV.w : mask_rectUV.y);
        }
    ), SHADER_STRINGIFY(130,
//\#version 130
    uniform sampler2D tMotif;
    uniform sampler2D tMask;
    uniform float vHeight;
    in vec2 texUV;
    in vec2 maskUV;
    uniform vec4 id_color;
        void main(){
            vec4 colselect= texture(tMotif, texUV);
			vec4 colmask = texture(tMask, maskUV);
            if (colmask.w <= 0.25) discard;
			if ((colselect.x < 127.0f / 510.0f)&&( abs(colselect.z - 0.5) > 0.25)) discard;
            gl_FragColor = id_color;
            }
    ),"daicon_alias_shader");
    glUseProgram(daicon_alias_shader);
    glUniform1i(glGetUniformLocation(daicon_alias_shader, "tMotif"),0);
    glUniform1i(glGetUniformLocation(daicon_alias_shader, "tMask"),1);
   LFH_ALIVE;
    ctrl_state.daheat_shader = ctrl_state.compileshader_routine( SHADER_STRINGIFY(130,
//\#version 130
    uniform vec4 boundrect;
    uniform vec2 windowsize;
    uniform vec4 projection;
    uniform vec4 selection;
    uniform float timephase;
    uniform vec2 val_range;
    out vec3 palUV;

    void main(void){
        float cval = gl_Vertex.z < val_range.x ? 0.0f : (gl_Vertex.z > val_range.y ? 1.0f : (gl_Vertex.z - val_range.x) / (val_range.y - val_range.x));
        vec2 tmppos;
        if (gl_Vertex.x < projection.x) tmppos.x = 0.0f;
        else if (gl_Vertex.x > projection.x + projection.z) tmppos.x = projection.z;
        else tmppos.x =gl_Vertex.x - projection.x;
        if (gl_Vertex.y < projection.y) tmppos.y = 0.0f;
        else if (gl_Vertex.y > projection.y + projection.w)  tmppos.y = projection.w;
        else tmppos.y =gl_Vertex.y - projection.y;
        gl_Position = vec4( -1.0f + 2.0f * (0.5f + boundrect.x + (tmppos.x * (boundrect.z - boundrect.x)) / projection.z ) / windowsize.x
                          , -1.0f + 2.0f * (0.5f + boundrect.w + (tmppos.y * (boundrect.y - boundrect.w)) / projection.w ) / windowsize.y
                          , 0.0,1);

        vec2 cpos = gl_Vertex.xy;
        cpos.x += (gl_Vertex.w < 0.5) ? 0.5 : -0.5;
        cpos.y += (abs(gl_Vertex.w-0.5) > 0.25) ? 0.5 : -0.5;
        if ((selection.x > cpos.x)||(selection.z < cpos.x)||(selection.y > cpos.y)||(selection.w < cpos.y)){
            palUV = vec3(0.8f - 0.7f * cval,0.0078125f+ cval *0.234375f, 1.0);
        }else{
            palUV = vec3(1.3f - 0.7f * cval, 0.0078125f+ (timephase + cval) *0.1171875f, 0.2);
        }
        }
    ), SHADER_STRINGIFY(130,
//\#version 130
    uniform sampler2D tPalette;
    uniform vec4 boundrect;
    in vec3 palUV;
        void main(){
            gl_FragDepth = 0.0f;
			gl_FragColor = mix(vec4(palUV.y * 4,palUV.y * 4,palUV.y * 4,1.0f), texture(tPalette, palUV.xy), palUV.z);
            }
    ),"daheat_shader");
   LFH_ALIVE;
    glUseProgram(daheat_shader);
    glUniform1i(glGetUniformLocation(daheat_shader, "tPalette"),0);
   LFH_ALIVE;
}
void Controlstate::create_latent_thread(){
}

bool Controlstate::query(INPUTSTATE_enum what) const{
	switch(what){
		case INPUTSTATE_PRESSED_MOUSELEFT: return (button_state & 1) != 0;
		case INPUTSTATE_PRESSED_MOUSERIGHT: return (button_state & 2) != 0;
		case INPUTSTATE_PRESSED_MOUSECENTER: return (button_state & 4) != 0;
		default: return (button_state & 8) != 0;
	}
}
void Controlstate::main_control_loop(RessourceLoader* res_loader){
    foreground_alias = 0;
    ressourceHanddle = res_loader;
    SDL_Event event;
    unsigned int d;
    //LFHPrimitive::Event<GUImessage>* tmp_event;
    GUImessage tmp_message;
	tooltipalias = 0;
    stringCapture_cur = 0x1000;
    stringCapture_endptr = 0x1000;
LFH_ALIVE;
    ctrl_state.compiledefaultshaders();
LFH_ALIVE;
	ctrl_state.create_latent_thread();
	int i,gui_out;
	ctrl_state.mouse_coor[0] =0;ctrl_state.mouse_coor[1] =0;
LFH_ALIVE;
    mouse_click_milli=0;
    last_mouse_button_event_time =0;
    ctrl_state.mouse_stencil =0;
    ctrl_state.mouse_stencilID =0;
    ctrl_state.mouse_array_offset=0;
    ctrl_state.mouse_push_coor[0] =0; ctrl_state.mouse_push_coor[1] =0;
    ctrl_state.foreground_alias=0;
    ctrl_state.button_state=0;
LFH_ALIVE;



	// makes
	/*glGenBuffers(1, &(quadindex_buffer));
	uint16_t* quad_buffer = new uint16_t[98304];
	for(*/

LFH_ALIVE;
    while(ctrl_state.states.size() != 0) {
        if (!thrbase.isRunning()) {
            thrbase.stopThreads();
            ctrl_state.close_procstates(0);
            break;
        }
        timequeue.runTo(SDL_GetTicks());
LFH_ALIVE;
        d = ctrl_state.states.size();
        while(SDL_PollEvent(&event)) {
        switch(event.type) {
        case SDL_MOUSEMOTION:
            mouse_coor[0] = event.motion.x;
            mouse_coor[1] = event.motion.y;
            tmp_message.coor[0] = event.motion.x;
            tmp_message.coor[1] = event.motion.y;
            tmp_message.msg_key = ((mouse_stencil & 0xFE) == 2) ? mouse_stencilID : curwin->GUI_alias;
            if (button_state & 7){
                if (button_state & 8) tmp_message.type = GUIMSG_MOUSE_DRAG;
                else{
                    if (((mouse_push_coor[0] - mouse_coor[0]+1) < 3)||(((mouse_push_coor[1] - mouse_coor[1])+1) < 3)){ // is equal to +- 1 pixel
                        button_state |= 8;// start dragging!
                        switch(button_state & 7){
                            case 1:tmp_message.type = GUIMSG_MOUSE_DRAG_LBUTTON; break;
                            case 2:tmp_message.type = GUIMSG_MOUSE_DRAG_RBUTTON; break;
                            case 3:tmp_message.type = GUIMSG_MOUSE_DRAG_MBUTTON; break;
                        }
                    } else break; // nothing to report!
                }
            }else tmp_message.type = GUIMSG_MOUSE_MOVE;
            tmp_message.msg_key = ((mouse_stencil & 0xFE) == 2) ? mouse_stencilID : curwin->GUI_alias;
            tmp_message.type = gui_objects_ptr[tmp_message.msg_key]->processGUIevent(tmp_message.type);
            if (tmp_message.type != GUIMSG_NULL) tmp_message();
        break;
        case SDL_MOUSEBUTTONDOWN:
            //printf("is in a gui down(ID:%i)!\n", mouse_stencilID);
            tmp_message.msg_key = ((mouse_stencil & 0xFE) == 2) ? mouse_stencilID : curwin->GUI_alias;
            last_mouse_button_event_time = SDL_GetTicks();
            timequeue.insert_sync(last_mouse_button_event_time + mouse_click_milli, new MouseStare(last_mouse_button_event_time, tmp_message.msg_key, curwin->GUI_alias,event.button.x, event.button.y));
            tmp_message.coor[0] = event.button.x;
            tmp_message.coor[1] = event.button.y;
            mouse_coor[0] = event.button.x;
            mouse_coor[1] = event.button.y;
            if ((button_state & 7) == 0){
                mouse_push_coor[0] = mouse_coor[0];
                mouse_push_coor[1] = mouse_coor[1];
                switch(event.button.button) {
                    case SDL_BUTTON_LEFT: button_state |= 1; tmp_message.type = GUIMSG_MOUSE_DOWN_LBUTTON; break;
                    case SDL_BUTTON_RIGHT: button_state |= 2; tmp_message.type = GUIMSG_MOUSE_DOWN_RBUTTON; break;
                    case SDL_BUTTON_MIDDLE: button_state |= 3; tmp_message.type = GUIMSG_MOUSE_DOWN_MBUTTON; break;
                }
                tmp_message.parameter2 = button_state;
            //    printf("%i is gui id, %i is found\n", tmp_message.msg_key, gui_objects_ptr.find(tmp_message.msg_key) ); fflush(stdout);
                // render to get subID!
                curwin->stencilval[0] = mouse_stencil; curwin->renderAlias();


                tmp_message.type = gui_objects_ptr[tmp_message.msg_key]->processGUIevent(tmp_message.type);
                if (tmp_message.type != GUIMSG_NULL) tmp_message();
            } // currently dragging, ignore!



            break;
        case SDL_MOUSEBUTTONUP: {
                //printf("is in a gui up (ID:%i)!\n", mouse_stencilID);
            tmp_message.msg_key = ((mouse_stencil & 0xFE) == 2) ? mouse_stencilID : curwin->GUI_alias;
            tmp_message.coor[0] = event.button.x;
            tmp_message.coor[1] = event.button.y;
            mouse_coor[0] = event.button.x;
            mouse_coor[1] = event.button.y;
            if (!gui_objects_ptr[tmp_message.msg_key].isValid()) exit(1);
            last_mouse_button_event_time = SDL_GetTicks();
            switch(event.button.button) {
                case SDL_BUTTON_LEFT:
                    if ((button_state & 7) == 1){
                        tmp_message.type = (button_state & 8) ? GUIMSG_MOUSE_DROP_LBUTTON : GUIMSG_MOUSE_CLICK_LBUTTON; button_state &= 0xFFFFFFF0;
                    }else{
                        tmp_message.type = GUIMSG_MOUSE_CONTEXT_CLICK_LBUTTON;
                    }
                     break;
                case SDL_BUTTON_RIGHT:
                    if ((button_state & 7) == 2){
                        tmp_message.type = (button_state & 8) ? GUIMSG_MOUSE_DROP_RBUTTON : GUIMSG_MOUSE_CLICK_RBUTTON; button_state &= 0xFFFFFFF0;
                    }else{
                        tmp_message.type = GUIMSG_MOUSE_CONTEXT_CLICK_RBUTTON;
                    }
                    break;
                case SDL_BUTTON_MIDDLE:
                    if ((button_state & 7) == 3){
                        tmp_message.type = (button_state & 8) ? GUIMSG_MOUSE_DROP_MBUTTON : GUIMSG_MOUSE_CLICK_MBUTTON; button_state &= 0xFFFFFFF0;
                    }else{
                        tmp_message.type = GUIMSG_MOUSE_CONTEXT_CLICK_MBUTTON;
                    }
                    break;
                default:printf("fired up %i\n", event.button.button); //case SDL_BUTTON_MIDDLE: button_state &= 0xFFFFFFFF ^ 4; break;
            }
            tmp_message.type =gui_objects_ptr[tmp_message.msg_key]->processGUIevent(tmp_message.type);
            if (tmp_message.type != GUIMSG_NULL) tmp_message();
        }break;
        case SDL_MOUSEWHEEL:{
            event.key.keysym.sym =  (event.wheel.y < 0) ? SDLK_SLEEP +1 : SDLK_SLEEP+2;
            while (d--){
                if ((ctrl_state.states[d]->OnKeyUp(event.key)) == 0) break;
            }
        }break;
        case SDL_KEYDOWN: {
            switch(event.key.keysym.sym){
                case SDLK_LALT: button_state |= 64; break;
                case SDLK_LCTRL: button_state |= 32; break;
                case SDLK_LSHIFT: button_state |= 16; break;
                default: break;
            }
            if ((stringCapture_cur & 0x1000) == 0x1000){
                while (d--){
                    if ((ctrl_state.states[d]->OnKeyDown(event.key)) == 0) break;
                }
            }
            break;
        }
        case SDL_KEYUP: {
            int leng;
            switch(event.key.keysym.sym){
                case SDLK_LALT: button_state &= 0xFFFFFFBF; break;
                case SDLK_LCTRL: button_state &= 0xFFFFFFDF; break;
                case SDLK_LSHIFT: button_state &= 0xFFFFFFEF; break;
                default: break;
            }
            if ((stringCapture_cur & 0x1000) == 0x1000){
                while (d--){
                    if ((ctrl_state.states[d]->OnKeyUp(event.key)) == 0) break;
                }
            }else{

                if (button_state & 96){
                    if ((button_state & 112) == 32){
                        switch(event.key.keysym.sym) {
                        case SDLK_LEFT:
                        if (ctrl_state.stringCapture_cur != 0){
                            do ctrl_state.stringCapture_cur--;
                            while ((ctrl_state.stringCapture_cur != 0)&&(stringCapture[ctrl_state.stringCapture_cur] != 32));
                            ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                        }
                        break;
                        case SDLK_RIGHT:
                        if (stringCapture[ctrl_state.stringCapture_cur] != '\0'){
                            do ctrl_state.stringCapture_cur++;
                            while ((stringCapture[ctrl_state.stringCapture_cur] & 223) != 0);
                            ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                        }
                        break;
                        case 'v':{
                            char* clip = SDL_GetClipboardText();
                            if (clip){
                                gui_out = strlen(stringCapture + ctrl_state.stringCapture_endptr);
                                i =strlen(clip);
                                memcpy(stringCapture+ ctrl_state.stringCapture_endptr + i, stringCapture+ ctrl_state.stringCapture_endptr, gui_out+1);
                                memcpy(stringCapture+ ctrl_state.stringCapture_endptr, clip, i);
                                ctrl_state.stringCapture_endptr += i;
                                ctrl_state.stringCapture_cur += i;
                                ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                                SDL_free(clip);
                            }
                        }break;
                        case 'c': SDL_SetClipboardText(stringCapture); break;
                        case 'x': SDL_SetClipboardText(stringCapture); ctrl_state.stringCapture_cur = ctrl_state.stringCapture_endptr =0; ctrl_state.stringCapture[ctrl_state.stringCapture_cur] = '\0'; ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;break;
                        }
                    }
                }else{
                    switch(event.key.keysym.sym) {
                    case SDLK_UP: printf("got up arrow!\n"); break;
                    case SDLK_DOWN: printf("got down arrow!\n"); break;
                    case SDLK_LEFT:
                        if (ctrl_state.stringCapture_cur != 0){
                            do ctrl_state.stringCapture_cur--;
                            while ((stringCapture[ctrl_state.stringCapture_cur] & 192) == 128);
                            ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                        }
                        break;
                    case SDLK_RIGHT:
                        if (stringCapture[ctrl_state.stringCapture_cur] != '\0'){
                            do ctrl_state.stringCapture_cur++;
                            while ((stringCapture[ctrl_state.stringCapture_cur] & 192) == 128);
                            ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                        }
                        break;
                    case SDLK_DELETE:
                        if (ctrl_state.stringCapture[ctrl_state.stringCapture_cur] != '\0'){
                            leng=1;
                            while((ctrl_state.stringCapture[ctrl_state.stringCapture_cur +leng] & 192) == 128) leng++;
                            for(d= ctrl_state.stringCapture_cur+leng; ctrl_state.stringCapture[d] != '\0';d++) ctrl_state.stringCapture[d-leng] = ctrl_state.stringCapture[d];
                            ctrl_state.stringCapture[d-leng] = ctrl_state.stringCapture[d];
                            ctrl_state.stringCapture_endptr -= leng;
                            ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                        }
                        break;
                    case '\r': exit_stringcapture_mode(); break;
                    case '\b':
                            if (ctrl_state.stringCapture_cur){
                            leng=1;
                            while((ctrl_state.stringCapture[ctrl_state.stringCapture_cur -leng] & 192) == 128) leng--;
                            for(d= ctrl_state.stringCapture_cur; ctrl_state.stringCapture[d] != '\0';d++) ctrl_state.stringCapture[d-leng] = ctrl_state.stringCapture[d];
                            ctrl_state.stringCapture[d-leng] = ctrl_state.stringCapture[d];
                            ctrl_state.stringCapture_endptr -= leng;
                            ctrl_state.stringCapture_cur -= leng;
                            ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
                            } break;
                    }
                    d=0;
                }
            }
            break;
        }
        case SDL_TEXTEDITING:{
            if ((stringCapture_cur & 0x1000) == 0){
            stringCapture_endptr = stringCapture_cur;
            for(d=0;event.edit.text[d] != '\0';d++) ctrl_state.stringCapture[d + ctrl_state.stringCapture_endptr] = event.edit.text[d];
            stringCapture_endptr += d;
            stringCapture[ctrl_state.stringCapture_endptr] = '\0';
			// ((GUITextArea*)gui_objects_ptr[stringCapture_guiID])->state |= 5;
            }
            d=0;
        }break;
        case SDL_TEXTINPUT:{
            int leng;
            if ((stringCapture_cur & 0x1000) == 0){
            leng=strlen(event.text.text);
            if (leng){
                stringCapture_endptr +=leng;
                for(d=stringCapture_endptr; d >= stringCapture_cur + leng ;d--) {
                    ctrl_state.stringCapture[d] = ctrl_state.stringCapture[d - leng];
                }
                for(d=0;d < leng;d++) ctrl_state.stringCapture[stringCapture_cur + d] = event.text.text[d];
                stringCapture_cur +=leng;
                ((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state |= 5;
            }
            }
            d=0;
        }break;
        case SDL_JOYAXISMOTION: {
       //        OnJoyAxis(event.jaxis.which,event.jaxis.axis,event.jaxis.value);

        }break;
        case SDL_JOYBALLMOTION: {
      //         OnJoyBall(event.jball.which,event.jball.ball,event.jball.xrel,event.jball.yrel);
            break;
        }
        case SDL_JOYHATMOTION: {
       //       OnJoyHat(event.jhat.which,event.jhat.hat,event.jhat.value);
            break;
        }
        case SDL_JOYBUTTONDOWN: {
     //          OnJoyButtonDown(event.jbutton.which,event.jbutton.button);
            break;
        }
        case SDL_JOYBUTTONUP: {
      //         OnJoyButtonUp(event.jbutton.which,event.jbutton.button);
            break;
        }
        case SDL_QUIT: // alt-f4!
            passive.state = 0;
            close_procstates(0);
            return;

      //         OnExit();
        case SDL_SYSWMEVENT: {
            //Ignore
            break;
        }

//        case SDL_VIDEORESIZE: {
      //         OnResize(event.resize.w,event.resize.h);
 //           break;
  //      }

  /*      case SDL_VIDEOEXPOSE: {
      //         OnExpose();
            break;
        }*/
        case SDL_WINDOWEVENT:{

            switch (event.window.event) {
            case SDL_WINDOWEVENT_RESIZED:
                SDL_Log("Window %d resized to %dx%d",
                        event.window.windowID, event.window.data1,
                        event.window.data2);
            break;
           /* case SDL_WINDOWEVENT_SHOWN: SDL_Log("Window %d shown", event.window.windowID); break;
            case SDL_WINDOWEVENT_HIDDEN: SDL_Log("Window %d hidden", event.window.windowID); break;
            case SDL_WINDOWEVENT_EXPOSED: SDL_Log("Window %d exposed", event.window.windowID); break;
            case SDL_WINDOWEVENT_MOVED:
                SDL_Log("Window %d moved to %d,%d",
                        event.window.windowID, event.window.data1,
                        event.window.data2);
                break;

            case SDL_WINDOWEVENT_MINIMIZED:
                SDL_Log("Window %d minimized", event.window.windowID);
                break;
            case SDL_WINDOWEVENT_MAXIMIZED:
                SDL_Log("Window %d maximized", event.window.windowID);
                break;
            case SDL_WINDOWEVENT_RESTORED:
                SDL_Log("Window %d restored", event.window.windowID);
                break;
            case SDL_WINDOWEVENT_ENTER:
                SDL_Log("Mouse entered window %d",
                        event.window.windowID);
                break;
            case SDL_WINDOWEVENT_LEAVE:
                SDL_Log("Mouse left window %d", event.window.windowID);
                break;
            case SDL_WINDOWEVENT_FOCUS_GAINED:
                SDL_Log("Window %d gained keyboard focus",
                        event.window.windowID);
                break;
            case SDL_WINDOWEVENT_FOCUS_LOST:
                SDL_Log("Window %d lost keyboard focus",
                        event.window.windowID);
                break;
            case SDL_WINDOWEVENT_CLOSE:

                SDL_Log("Window %d closed", event.window.windowID);
                            tmp_message.type = LFHGUI_QUITPROCESS;
            tmp_message();
                break;*/
            default:
               // SDL_Log("Window %d got unknown event %d", event.window.windowID, event.window.event);
                break;
            }
        }break;
        default: {
            printf("%i obtained!\n", event.type);
       //     OnUser(event.user.type,event.user.code,event.user.data1,event.user.data2);
            break;
        }
    }
        d = ctrl_state.states.size();
        }

        if (d > ctrl_state.states.size()){ // some scope popped, notify current
            if (ctrl_state.states.size() == 0) break;
            tmp_message.type = LFHGUI_SCOPE_FOCUS;
            printf("will %i %i", d, ctrl_state.states.size());
            tmp_message();
        }
/*
        if (fetchMessage(tmp_message)){ //
            for(d = ctrl_state.states.size(); (d--);) if (ctrl_state.states[d]->listen(tmp_message) == 0) break;
            if (d == 0){
                if (tmp_message.type == LFHGUI_QUITPROCESS) break;
            }
        }
*/
         for(d = ctrl_state.states.size(); (d--) != 0 ;) ctrl_state.states[d]->OnMaintain();

        latentqueue.pop_exec();
        latentqueue.pop_exec();
        latentqueue.pop_exec();
        latentqueue.pop_exec();

    }
    passive.state = 0;
    // exiting! delete elems in stack;
    close_procstates(0);
}

void Controlstate::notifyStencilChange(int n_mouse_stencil, GLint *color_and_stencilID, bool hasextra){
  //  printf("stencil %i->%i  (ID=%X->%X)\n", mouse_stencil, n_mouse_stencil, mouse_stencilID, n_mouse_stencilID);
	if (tooltipalias !=0){
		curwin->removeGUI(tooltipalias);
		tooltipalias = 0;
	}
    GUImessage tmp_message;
    if ((mouse_stencil & 0xFE) == 2){ // GUI is losing focus
        tmp_message.msg_key = ctrl_state.mouse_stencilID;
        tmp_message.type = gui_objects_ptr[mouse_stencilID]->processGUIevent(GUIMSG_MOUSE_EXIT);
        if (tmp_message.type != GUIMSG_NULL) tmp_message();
    }

    mouse_stencil = n_mouse_stencil;
    if (hasextra) {mouse_stencilID = (color_and_stencilID[0] & 0xFFFFFF00) | color_and_stencilID[1]; mouse_array_offset = color_and_stencilID[0] &0xFF;}
    else mouse_stencilID = color_and_stencilID[0];

    if ((mouse_stencil & 0xFE) == 2){ // GUI is losing focus
		if (!gui_objects_ptr[mouse_stencilID].isValid()){
			printf("Warning, illegal mouseID %X!!!\n", mouse_stencilID);
			mouse_stencilID = curwin->GUI_alias;
            mouse_array_offset = 0xFFFFFFFF;
		}else{
            mouse_array_offset = ctrl_state.gui_objects_ptr[gui_objects_ptr[mouse_stencilID].parent_alias]->getArrayIDfromMouse();
		}
        tmp_message.type = gui_objects_ptr[mouse_stencilID]->processGUIevent(GUIMSG_MOUSE_ENTER);
        tmp_message.msg_key = ctrl_state.mouse_stencilID;
        if (tmp_message.type != GUIMSG_NULL) tmp_message();
    }else{
        tmp_message.msg_key = ctrl_state.mouse_stencil;
        tmp_message.type = LFHGUI_MOUSE_STENCILID_CHANGE;
        tmp_message();
    }
}

void Controlstate::close_procstates(unsigned int to_level){
    unsigned int d;
    for(d= ctrl_state.states.size(); (d--) > to_level ;){
        ProcessState* curproc = ctrl_state.states[d];
        ctrl_state.states.pop_back();
        delete(curproc);
    }
}
void Controlstate::enter_stringcapture_mode(GUITextArea* guitext){
	if (guitext == NULL) {
		stringCapture_endptr =stringCapture_cur =0;foreground_alias =0;
	}else{
		ctrl_state.foreground_alias = guitext->GUI_alias;
		if (guitext->cur_str == NULL){
			stringCapture_endptr =stringCapture_cur =0;
		}else{
			strcpy(ctrl_state.stringCapture, guitext->cur_str);
			stringCapture_endptr =stringCapture_cur = strlen(ctrl_state.stringCapture);
		}
		guitext->state |= 4;
	}
	ctrl_state.stringCapture[stringCapture_cur] = '\0';
	SDL_StartTextInput();

}
void Controlstate::exit_stringcapture_mode(){
    GUImessage tmp_message;
	stringCapture_cur = 0x1000;  stringCapture_endptr = 0x1000;
	tmp_message.msg_key = foreground_alias;
	tmp_message.type = LFHGUI_TEXTEDIT_DONE;
	((GUITextArea*)gui_objects_ptr[foreground_alias].target)->state &= 0xFFFFFFFB;
	((GUITextArea*)gui_objects_ptr[foreground_alias].target)->setText(stringCapture);
	if (gui_objects_ptr[foreground_alias].parent_alias != 0) this->clearForgroundObject();
	int d = ctrl_state.states.size();
	while (d--){
		if ((ctrl_state.states[d]->listen(tmp_message)) == 0) break;
	}
}
/*
void Controlstate::drawtext(Uint16* text,unsigned int f, unsigned int x,unsigned int y,unsigned int &w,unsigned int &h){
        unsigned int s = gui_fonts[f]->w +1;
        Uint16* pcur = text;
//        Uint16* lcur = text;
        Uint16* wcur = text;
        unsigned int j = wordwidth(pcur);

        unsigned int vs = (gui_fonts[f]->h / 96) + 1;
        unsigned int pos[2]; pos[0] = x; pos[1] =y;
        if (w != 0) { // width is given!
        unsigned short nb_space;
        Uint32 k =0;
        while(wcur[0] != '\0'){
        //    printf("pist %i\t%i\t%i\n", wcur - lcur, text - lcur,wcur[0] );
            unsigned int l = findLine(wcur,s,w,nb_space);
        //    printf("dist %i\t%i\n", wcur - lcur, text - lcur);
            for(;text != wcur;text++){
            switch((*text) >> 8){
                case 0:
                    if ((*text) == '\0') break;
                    if ((*text) == ' ') {pos[0] +=nb_space;break;}
                    j+=4;
                    if ((!((*text) & 0x80 )) && ((*text) & 0x60 )) {gui_fonts[f]->draw_Vslice(pos[0],pos[1], (((*text) - 32) * (vs-1) ),vs-1); pos[0] += s;}
                break;
                }
            }
            pos[1] -= vs; pos[0] = x;
            if (wcur[0] == ' ') {wcur++; text++;}
            if (k++ > 1000) exit(1);
        }
        h = y - pos[1] +vs -2;
        }
    }*/
unsigned short Controlstate::findLine(Uint16* &gtext, unsigned int s, unsigned int w, unsigned short& nbsp){unsigned short fout=0;
    bool loopy = true;
    unsigned short t=0;
    Uint16* text = gtext;nbsp =0;
    for(;loopy;text++){

    switch((*text) >> 8){
        case 0:
            if ((*text) == '\0') {loopy = false;break;}
            if ((*text) == ' ') {
                t +=s;
                if (t > w) {loopy = false;break;}
                else {gtext = text; fout = t;nbsp++;}
            }
            if ((!((*text) & 0x80 )) && ((*text) & 0x60 )) t += s;
        break;
        }
    }

    if (fout == 0) {gtext = text; fout = t; nbsp =10;}
    else nbsp = (w - fout + (nbsp+1) * s) / nbsp;
    return fout;
}
unsigned int Controlstate::wordwidth(Uint16* &text){
    unsigned int j;
    bool loopy = true;
    for(j=0;loopy;text++){
        switch((*text) >> 8){
            case 0:
                if (((*text) == '\0') || ((*text) == ' '))  {loopy = false;break;}
                j+=4;
            break;
            }
        }
    return j;
}
void Controlstate::convert(Uint16* targ, const char* sss){
        const char* s = sss;
        while((*s) !='\0') *(targ++) = (unsigned short) *(s++);
        *targ = 0;
    }

void Controlstate::manifestAsTooltip(uint32_t alias){
    GUIObject* obj =  gui_objects_ptr[alias].target;
    int32_t pos[2];
	if (tooltipalias != 0) ctrl_state.curwin->removeGUI(tooltipalias);
	if (obj == NULL) {fprintf(stderr,"GUIObject for tooltip does not exist!\n"); exit(1);}
	if ((int32_t)mouse_coor[0] < (obj->getRect()[2]>>1)) pos[0] = 0;
	else if ((int32_t)mouse_coor[0] > curwin->rect[2] - (obj->getRect()[2]>>1)) pos[0] = curwin->rect[2] - obj->getRect()[2];
	else pos[0] = mouse_coor[0] - (obj->getRect()[2]>>1);

	if ((int32_t)mouse_coor[1] < (obj->getRect()[3]>>1))  pos[1] = 0;
	else if ((int32_t)mouse_coor[1] > curwin->rect[3] - (obj->getRect()[3]>>1))  pos[1] = curwin->rect[3] - obj->getRect()[3];
	else pos[1] = mouse_coor[1] - (obj->getRect()[3]>>1);

	obj->setPos(pos);
	tooltipalias = alias;
    ctrl_state.curwin->insertGUI(obj);
}

unsigned int GUImessage::operator()(){
    unsigned int d = ctrl_state.states.getSize();
    while(d>0) {
		d--;
		if (ctrl_state.states[d]->listen(*this) == 0) break;
	}
    return(0);
}
bool GUImessage::validFilter(const unsigned int &){return(true);}
    StringcaptureProcess::~StringcaptureProcess(){
    //    GUImessage damsg;
    //    if (ctrl_state.stringC_cur){
    //        damsg
    //    }
    }
int StringcaptureProcess::OnKeyDown(const SDL_KeyboardEvent& event){ //OnKeyDown(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
// printf("%i\n",  event.keysym.sym);
switch(event.keysym.sym){
    case '\b':
        if (ctrl_state.stringCapture_cur){
            while ((ctrl_state.stringCapture[--ctrl_state.stringCapture_cur] & 0xC0) == 0x80) ;
            ctrl_state.stringCapture[ctrl_state.stringCapture_cur] = '\0';
            return(0);
        }
    case SDLK_RETURN:
        ctrl_state.states.pop_back(); delete(this); return(0);
    break;
    default:
//            if (event.keysym.sym >= SDLK_NUMLOCK) return 0;
        //  printf("%i\t%i\n", (int)event.keysym.sym, (int) event.keysym.unicode);
        ctrl_state.stringCapture[ctrl_state.stringCapture_cur++] = event.keysym.sym; ctrl_state.stringCapture[ctrl_state.stringCapture_cur] = '\0';

    }
return(0);}
int StringcaptureProcess::OnKeyUp(const SDL_KeyboardEvent& event){ //OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
    return(0);
}
int StringcaptureProcess::OnMaintain(){//OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
return 0;
}

GUIStyle::GUIStyle(): textlayout(GUITEXT_LAYOUT_TRUNCATE){
    setBorderSizes(16,16);
	dimentions[0] = 0;
	for(int i =0;i<8;i++) colmap_index[i] = i % 4;
}

/** \brief Sets default dimentions for the gui object, will overwrite dimension whenever object is moved (unless width =0)
 *
 * \param uint16_t width
 * \param uint16_t height
 * \return
 *
 */
void GUIStyle::setDefDimentions(uint16_t width,uint16_t height){
    dimentions[0]= width;dimentions[1]= height;
    }
/** \brief Sets default dimentions for the gui object, will overwrite dimension whenever object is moved (unless width =0)
 *
 * \param uint16_t border_size no-interaction band width
 * \param uint16_t scroll_size band width if object interaction on boundary is enabled
 * \return
 *
 */
void GUIStyle::setBorderSizes(int _border_size, int _scroll_size, bool update_text_borders){
    border_width = _border_size;
    scroll_width = _scroll_size;
    if (update_text_borders){
        text_borders[0] = border_width >> 1;
        text_borders[1] = border_width >> 1;
        text_borders[2] = border_width >> 1;
        text_borders[3] = border_width >> 1;
    }
}

/** \brief Sets color data at a given color index, no gradient in color
 *
 * \param color_index, 0-16
 * \param hue, 0.0f-1.0f
 * \param bright, 0.0f-1.0f
 * \param sat, 0.0f-1.0f
 * \return
 *
 */
void GUIStyle::setColorHIS(int color_index, double hue, double bright, double sat){
	colors[color_index][0] = hue;
	colors[color_index][1] = bright;
	colors[color_index][2] = sat;
	colors[color_index][3] = hue;
	colors[color_index][4] = bright;
	colors[color_index][5] = sat;
    colors[color_index][6] = 1.0f;
	colors[color_index][7] = 1.0f;
	colors[color_index][8] = 1.0f;
}

/** \brief Sets color data at a given color index
 *
 * \param color_index, 0-16
 * \param hue, 0.0f-1.0f
 * \param bright, 0.0f-1.0f
 * \param sat, 0.0f-1.0f
 * \param hue2, 0.0f-1.0f
 * \param bright2, 0.0f-1.0f
 * \param sat2, 0.0f-1.0f
 * \return
 *
 */
void GUIStyle::setColorHISpair(int color_index, double hue, double bright, double sat, double hue2, double bright2, double sat2){
	colors[color_index][0] = hue;
	colors[color_index][1] = bright;
	colors[color_index][2] = sat;
	colors[color_index][3] = hue2;
	colors[color_index][4] = bright2;
	colors[color_index][5] = sat2;
	colors[color_index][6] = 1.0f;
    colors[color_index][7] = 1.0f;
    colors[color_index][8] = 1.0f;
}

/** \brief Links a state and channel to a color index
 *
 * \param channel, 0-3
 * \param color_index, 0-16
 * \param flag, 0-1
 * \return
 *
 */
void GUIStyle::setStateColormap(int channel, int color_index, int flag){
    if (flag == 1){
        colmap_index[channel+ 4] = color_index;
    }else{
        colmap_index[channel] = color_index;
    }
}

/** \brief Sets the uniforms {"boundrect", "boundtexcoor"}
 *
 * \param shaderID
 * \param rect,
 * \param bound
 * \return
 *
 */
bool GUIStyle::setFrameUniforms(GLint shaderID, int32_t* rect, const int32_t* par_rect, bool isForeground){
    GLint locs[2];
    locs[0] = glGetUniformLocation(shaderID, "boundrect");
    locs[1] = glGetUniformLocation(shaderID, "boundtexcoor");
    if (par_rect != NULL){
        GLfloat coor[4];
        GLfloat UV[4];
        if ((rect[1] < par_rect[5])&&(!isForeground)){
            if (rect[1] + rect[3] < par_rect[5]) return false;
            UV[1] = (((GLfloat)(rect[3] - (par_rect[5] - rect[1]))) / rect[3]);
            coor[3] = ctrl_state.curwin->rect[3] - par_rect[1];
        }else{
            UV[1] = 1.0f;
            coor[3] = ctrl_state.curwin->rect[3] - par_rect[1] - (rect[1] - par_rect[5]);
        }
        if ((rect[1] + rect[3] > par_rect[5] + par_rect[3])&&(!isForeground)){
            if (rect[1] > par_rect[5] + par_rect[3]) return false;
            UV[3] = ((GLfloat)(rect[1] + rect[3] - par_rect[5] - par_rect[3])) / rect[3];
            coor[1] = ctrl_state.curwin->rect[3] - par_rect[1] - par_rect[3];
        }else{
            UV[3] = 0.0f;
            coor[1] = ctrl_state.curwin->rect[3] - par_rect[1] - (rect[1] - par_rect[5]) - rect[3];
        }
        glUniform4f(locs[0], rect[0] + par_rect[0] - par_rect[4], coor[1], rect[0] + par_rect[0] + rect[2]- par_rect[4], coor[3]);
        glUniform4f(locs[1], 0.0f,UV[1],1.0f,UV[3]);
    }else{
        glUniform4f(locs[0], rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
        glUniform4f(locs[1], 0.0f,0.0f,1.0f,1.0f);
    }

    glUniform2f(glGetUniformLocation(shaderID, "borderUV"), ((GLfloat)border_width) / rect[2] , ((GLfloat)border_width) / rect[3]);

    return true;
}

void GUIStyle::setColorUniforms(GLint shaderID, int32_t state){
    if (state){
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color"), 1, false, colors[colmap_index[4]]);
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color2"), 1, false, colors[colmap_index[5]]);
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color3"), 1, false, colors[colmap_index[6]]);
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color4"), 1, false, colors[colmap_index[7]]);
    }else{
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color"), 1, false, colors[colmap_index[0]]);
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color2"), 1, false, colors[colmap_index[1]]);
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color3"), 1, false, colors[colmap_index[2]]);
        glUniformMatrix3fv(glGetUniformLocation(shaderID, "color4"), 1, false, colors[colmap_index[3]]);
    }
}


int GUIStyle::makeTextMesh(GuiTextAttribute &texta, const char* utf8string, int cursor_pos){
    unsigned int i, cursor;
    Vector< uint32_t> unicode;
    uint32_t color = 0x00000000;
    texta.glbuffer_indexes[1] = 0;
    if (texta.glbuffer_indexes[0] == 0) glGenBuffers(1, texta.glbuffer_indexes);
    glBindBuffer(GL_ARRAY_BUFFER, texta.glbuffer_indexes[0]);
    //printf("single makeing textmesh%i with %s\n", texta.glbuffer_indexes[0], utf8string);

    int wrap_width;
    if (textlayout != GUITEXT_LAYOUT_TRUNCATE){
        if (dimentions[0] == 0) {fprintf(stderr, "GUIText layout error! no default width\n"); exit(1);}
        wrap_width = (dimentions[0] - text_borders[0] - text_borders[2]) / font_op.font;
    } else wrap_width = 0;

    for(const char* cur = utf8string; ;cur++){
        if ((int)(cur -  utf8string) == cursor_pos) {unicode.push_back(0x10000); texta.glbuffer_indexes[1]++;}
        if ((*cur) == '\0') break;
        if ((*cur) & 128){
            if ((*cur) & 32){
                if ((*cur) & 16){
                    if ((*cur) & 8){
                        if (((*cur) & 255) ==  255){
                                // control character, selects color!
                                unicode.push_back( ((((unsigned int)cur[1]) << 8 ) & 0xFF00) | (((unsigned int)cur[2]) & 0xFF) | 0x30000 );
                                cur+=2;
                            }// else... realllllly? 5byte and more?
                    }else{
                        unicode.push_back((((unsigned int)cur[3]) & 0x3F) | ((((unsigned int)cur[2])<< 6) & 0xFC0)| ((((unsigned int)cur[1])<< 12) & 0x3F000) | ((((unsigned int)cur[1])<< 12) & 0x1C0000) );
                        cur+=3;texta.glbuffer_indexes[1]++;
                    }
                 //   printf("the hell?\n"); exit(1);
                }else{ // 3 byte
                    unicode.push_back((((unsigned int)cur[2]) & 0x3F) | ((((unsigned int)cur[1])<< 6) & 0xFC0)| ((((unsigned int)cur[0])<< 12) & 0xF000) );
                    cur+=2;texta.glbuffer_indexes[1]++;
                }
            }else{ // 2 byte
                unicode.push_back((((unsigned int)cur[1]) & 0x3F) | ((((unsigned int)cur[0])<< 6) & 0x7C0) );
                texta.glbuffer_indexes[1]++;
                cur++;
            }
        }else switch(*cur){
            case '\n':
            case ' ':
            case '\t':
            case '\r':
                unicode.push_back(((unsigned int)*cur) | 0x10000);
            break;
            case '\a':
                color =  0x30000;
                color |= (cur[1] & 64) ? ((cur[1] & 15) + 9) : (cur[1] & 15);
                color |= (((cur[2] & 64) ? ((cur[2] & 15) + 9) : (cur[2] & 15)) << 4);
                color |= (((cur[3] & 64) ? ((cur[3] & 15) + 9) : (cur[3] & 15)) << 8);
                color |= (((cur[4] & 64) ? ((cur[4] & 15) + 9) : (cur[4] & 15)) << 12);
                unicode.push_back( color );
                cur+=4;
            break;
            default:
            texta.glbuffer_indexes[1]++;
                unicode.push_back(((unsigned int)*cur));
            break;
        }
    }
    float* vbuf = new float[texta.glbuffer_indexes[1]<<4];
    int off[2]; ExOp::toZero(off);
    color = 0x00000000;
    for(i=0,cursor=0;cursor<unicode.getSize();cursor++){
        if (unicode[cursor] & 0x10000){
            if (unicode[cursor] & 0x20000){
                color = unicode[cursor] << 16;
            }else{
                switch(unicode[cursor]){
                case '\0' | 0x10000:
                    vbuf[i] = off[0] - 1;
                    vbuf[i|1] = (font_op.flip_ydir) ? -off[1] : off[1];
                    vbuf[i|2] = 0.0f;
                    vbuf[i | 14] = 1.0f;
                    vbuf[i | 10] = 0.625f;
                    vbuf[i | 6] = 0.375f;
                    ((uint32_t*)vbuf)[i|3] =  124;
                    vbuf[i | 12] = vbuf[i | 8] = vbuf[i | 4] =vbuf[i];
                    vbuf[i | 13] = vbuf[i | 9] = vbuf[i | 5] =vbuf[i|1];

                    ((uint32_t*)vbuf)[i | 7 ] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
                    ((uint32_t*)vbuf)[i | 11] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
                    ((uint32_t*)vbuf)[i | 15] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
                    i+= 16;

                    break; // cursor
                case ' '|0x10000: off[0]+=2; break;
                case '\t'|0x10000: off[0] += 8; break;
                case '\n'|0x10000: off[0] = 0; off[1] += 2;  break;
                }
            }
        }else{
        vbuf[i] = off[0];
        vbuf[i|1] = (font_op.flip_ydir) ? -off[1] : off[1];
        vbuf[i|2] = 0.0f;
        vbuf[i | 14] = 1.0f;
        vbuf[i | 10] = 0.625f;
        vbuf[i | 6] = 0.375f;
        if ( ((unicode[cursor] >= 0x3040)&& (unicode[cursor] <= 0x30FF))
           ||((unicode[cursor] >= 0x4E00)&& (unicode[cursor] <= 0x9FFF))
           ||((unicode[cursor] >= 0xFF01)&& (unicode[cursor] <= 0xFF60)) ) off[0]+=4;
        else off[0]+=2;

        if ((off[0]/2 >= wrap_width)&&(wrap_width != 0)){
            off[0] = 0;
            off[1] += 2;
        }
        ((uint32_t*)vbuf)[i|3] =  color | unicode[cursor];
        vbuf[i | 12] = vbuf[i | 8] = vbuf[i | 4] =vbuf[i];
        vbuf[i | 13] = vbuf[i | 9] = vbuf[i | 5] =vbuf[i|1];

        ((uint32_t*)vbuf)[i | 7 ] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
        ((uint32_t*)vbuf)[i | 11] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
        ((uint32_t*)vbuf)[i | 15] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
        i+= 16;
    }
    }
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (texta.glbuffer_indexes[1]<<4), vbuf, GL_STATIC_DRAW);
    delete[](vbuf);
    texta.glbuffer_indexes[1] <<=2;
    return off[1]+2;
}
int GUIStyle::makeTextMesh(GuiTextAttribute &texta, const Vector<const char*> &utf8strings, Vector<uint32_t> *offsets){
    unsigned int i, cursor,l;
    Vector< unsigned int > unicode;
    texta.glbuffer_indexes[1] = 0;
    if (texta.glbuffer_indexes[0] == 0) {
        glGenBuffers(1, texta.glbuffer_indexes);
    }
    glBindBuffer(GL_ARRAY_BUFFER, texta.glbuffer_indexes[0]);
    if (utf8strings.getSize() == 0){
        texta.glbuffer_indexes[1] =0; return 0;
    }
    //printf("vec makeing textmesh%i with %s\n", texta.glbuffer_indexes[0], utf8strings[0]);

    glEnableVertexAttribArray(ATTRIBUTE_POSITION);
    glEnableVertexAttribArray(ATTRIBURE_CHARID);
    glVertexAttribPointer(ATTRIBUTE_POSITION,3, GL_FLOAT, GL_FALSE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(0));
    glVertexAttribPointer(ATTRIBURE_CHARID,4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(sizeof(float)*3));

    uint32_t color;

    if (offsets != NULL) {offsets->setSize(1); (*offsets)[0] = 0;}

    int wrap_width;
    float xbase;
    if (textlayout == GUITEXT_LAYOUT_WRAP){
        if (dimentions[0] == 0) {fprintf(stderr, "GUIText layout error! no default width\n"); exit(1);}
        wrap_width = (dimentions[0] - text_borders[0] - text_borders[2]) / font_op.font;
        xbase = 0.5f * wrap_width;
    } else {wrap_width = 0; xbase =0;}

    int xcount =0;
    int lastbrpt =0;
    l=0;if (utf8strings.getSize() != 0) do{
		if (utf8strings[l] != NULL)
        for(const char* cur = utf8strings[l]; *cur != 0;cur++){
            //    printf("%c: %i\n", *cur, (int)*cur);
            if ((*cur) & 128){
                if ((*cur) & 32){
                    if ((*cur) & 16){
                        if ((*cur) & 8){
                            if (((*cur) & 255) ==  255){
                                // control character, selects color!
                                unicode.push_back( ((((unsigned int)cur[1]) << 8 ) & 0xFF00) | (((unsigned int)cur[2]) & 0xFF) | 0x30000 );
                                cur+=2;
                            }// else... realllllly? 5byte and more?
                        }else{ // 4 byte not really supported
                            unicode.push_back((((unsigned int)cur[3]) & 0x3F) | ((((unsigned int)cur[2])<< 6) & 0xFC0)| ((((unsigned int)cur[1])<< 12) & 0x3F000) | ((((unsigned int)cur[1])<< 12) & 0x1C0000) );
                            cur+=3;texta.glbuffer_indexes[1]++;
                        }

                    }else{ // 3 byte
                        unicode.push_back((((unsigned int)cur[2]) & 0x3F) | ((((unsigned int)cur[1])<< 6) & 0xFC0)| ((((unsigned int)cur[0])<< 12) & 0xF000) );
                        cur+=2;texta.glbuffer_indexes[1]++;
                    }
                }else{ // 2 byte
                    unicode.push_back((((unsigned int)cur[1]) & 0x3F) | ((((unsigned int)cur[0])<< 6) & 0x7C0) );
                    texta.glbuffer_indexes[1]++;
                    cur++;
                }
            }else switch(*cur){
                case ' ':
                    lastbrpt = xcount;
                case '\n':
                case '\t':
                case '\r':
                    unicode.push_back(((unsigned int)*cur) | 0x10000);
                break;
                case '\a':
                    color =  0x30000;
                    color |= (cur[1] & 64) ? ((cur[1] & 15) + 9) : (cur[1] & 15);
                    color |= (((cur[2] & 64) ? ((cur[2] & 15) + 9) : (cur[2] & 15)) << 4);
                    color |= (((cur[3] & 64) ? ((cur[3] & 15) + 9) : (cur[3] & 15)) << 8);
                    color |= (((cur[4] & 64) ? ((cur[4] & 15) + 9) : (cur[4] & 15)) << 12);
                    unicode.push_back( color );
                    cur+=4;
                break;
                default:
                texta.glbuffer_indexes[1]++;
                    unicode.push_back(((unsigned int)*cur));
                break;
            }
            xcount++;
            if (xcount == wrap_width){
                if (unicode[unicode.getSize()- wrap_width + lastbrpt] == 0x10020) unicode[unicode.getSize()- wrap_width + lastbrpt] = 0x10000 | '\n';
                else unicode.push_back(0x10000 | '\n');
                xcount = wrap_width;
            }
        }
        l++;
        if (l == utf8strings.getSize()) break;
        unicode.push_back(((unsigned int) '\n') | 0x10000);
	} while(true);

	if (unicode.last() != (0x10000 | '\n')) unicode.push_back(0x10000 | '\n');

	//printf("predicted length: %i\n", texta.glbuffer_indexes[1]);
    float* vbuf = new float[texta.glbuffer_indexes[1]<<4];
    int off[2]; ExOp::toZero(off);
    color = 0x00000000;
    int wordind =0;
    float cenoff;
    for(i=0,cursor=0;cursor<unicode.getSize();cursor++){

        if (unicode[cursor] & 0x10000){
            if (unicode[cursor] & 0x20000){
                color = unicode[cursor] << 16;
            }else{
                switch(unicode[cursor]){
                case ' '|0x10000: off[0]+=2; break;
                case '\t'|0x10000: off[0] += 8; break;
                case '\n'|0x10000:
                    if (textlayout == GUITEXT_LAYOUT_CENTERED){
                        cenoff = 0.5 * off[0];
                        while(wordind < i) {vbuf[wordind] -= cenoff; wordind += 4;}
                        off[1] += 2;
                    }else if (textlayout == GUITEXT_LAYOUT_CENTERED_ITEMS){
                        cenoff = 0.5 * off[0];
                        while(wordind < i) {vbuf[wordind] -= cenoff; wordind += 4;}
                    }else off[1] += 2;
                    off[0] = 0;
                    if (offsets != NULL) offsets->push_back(i >> 4);
                break;
                }
            }
        }else{
            vbuf[i] = xbase + off[0];
            vbuf[i|1] = (font_op.flip_ydir) ? -off[1] : off[1];
            vbuf[i|2] = 0.0f;
            vbuf[i | 14] = 1.0f;
            vbuf[i | 10] = 0.625f;
            vbuf[i | 6] = 0.375f;
            if ( ((unicode[cursor] >= 0x3040)&& (unicode[cursor] <= 0x30FF))
               ||((unicode[cursor] >= 0x4E00)&& (unicode[cursor] <= 0x9FFF))
               ||((unicode[cursor] >= 0xFF01)&& (unicode[cursor] <= 0xFF60)) ) off[0]+=4;
            else off[0]+=2;
            ((uint32_t*)vbuf)[i|3] =   color | unicode[cursor];
            vbuf[i | 12] = vbuf[i | 8] = vbuf[i | 4] =vbuf[i];
            vbuf[i | 13] = vbuf[i | 9] = vbuf[i | 5] =vbuf[i|1];
            vbuf[i | 15] = vbuf[i |11] = vbuf[i | 7] =vbuf[i|3];
         //   ((uint32_t*)vbuf)[i | 7 ] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
         //   ((uint32_t*)vbuf)[i | 11] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
         //   ((uint32_t*)vbuf)[i | 15] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
            i+= 16;
        }
    }
    if ((offsets != NULL)&&(offsets->getSize() <= utf8strings.getSize())) offsets->push_back(i >> 4);

    if ((textlayout == GUITEXT_LAYOUT_CENTERED)||(textlayout == GUITEXT_LAYOUT_CENTERED_ITEMS)) {
        cenoff = 0.5 * off[0];
        while(wordind < i) {vbuf[wordind] -= cenoff; wordind += 4;}
        if (textlayout == GUITEXT_LAYOUT_CENTERED_ITEMS) cenoff = 0.5;
        else cenoff = 0.5 * off[1];
        wordind = 1;
        while(wordind < i) {vbuf[wordind] -= cenoff; wordind += 4;}
    }


    //printf("actual length: %i\n", i>>4);

    /*
    if ((font_op.centered & 1)&&(i >= 4)){
        l = ((wrap_width << 1) - (int)vbuf[i-4]) >> 1;
        if (l != 0){
            j = i-4;
            vbuf[j] += l;
            while(j>0){
                if (vbuf[j|1] != vbuf[j-3]) break;
                j-= 4;
                vbuf[j] += l;
            }
        }
    }*/
    //printf("actual length: %i\n", i>>4);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (texta.glbuffer_indexes[1]<<4), vbuf, GL_STATIC_DRAW);
    delete[](vbuf);
    texta.glbuffer_indexes[1] <<=2;
    return off[1];
}

void GUIStyle::makeTextMesh(GuiTextAttribute &texta, const Vector<TextDisplayData> &utf8strings){
    unsigned int i, cursor,l;
    Vector< unsigned int > unicode;
    texta.glbuffer_indexes[1] = 0;
    if (texta.glbuffer_indexes[0] == 0) glGenBuffers(1, texta.glbuffer_indexes);
    glBindBuffer(GL_ARRAY_BUFFER, texta.glbuffer_indexes[0]);
    //printf("tdd makeing textmesh%i with %s\n", texta.glbuffer_indexes[0], utf8strings[0].text);

    uint32_t color;
    l=0;if (utf8strings.getSize() != 0) do{
		if (utf8strings[l].text != NULL)
        for(const char* cur = utf8strings[l].text; *cur != 0;cur++){
            //    printf("%c: %i\n", *cur, (int)*cur);
            if ((*cur) & 128){
                if ((*cur) & 32){
                    if ((*cur) & 16){
                        if ((*cur) & 8){
                            if (((*cur) & 255) ==  255){
                                // control character, selects color!
                                unicode.push_back( ((((unsigned int)cur[1]) << 8 ) & 0xFF00) | (((unsigned int)cur[2]) & 0xFF) | 0x30000 );
                                cur+=2;
                            }// else... realllllly? 5byte and more?
                        }else{ // 4 byte not really supported
                            unicode.push_back((((unsigned int)cur[3]) & 0x3F) | ((((unsigned int)cur[2])<< 6) & 0xFC0)| ((((unsigned int)cur[1])<< 12) & 0x3F000) | ((((unsigned int)cur[1])<< 12) & 0x1C0000) );
                            cur+=3;texta.glbuffer_indexes[1]++;
                        }

                    }else{ // 3 byte
                        unicode.push_back((((unsigned int)cur[2]) & 0x3F) | ((((unsigned int)cur[1])<< 6) & 0xFC0)| ((((unsigned int)cur[0])<< 12) & 0xF000) );
                        cur+=2;texta.glbuffer_indexes[1]++;
                    }
                }else{ // 2 byte
                    unicode.push_back((((unsigned int)cur[1]) & 0x3F) | ((((unsigned int)cur[0])<< 6) & 0x7C0) );
                    texta.glbuffer_indexes[1]++;
                    cur++;
                }
            }else switch(*cur){
                case '\n':
                case ' ':
                case '\t':
                case '\r':
                    unicode.push_back(((unsigned int)*cur) | 0x10000);
                break;
                case '\a':
                    color =  0x30000;
                    color |= (cur[1] & 64) ? ((cur[1] & 15) + 9) : (cur[1] & 15);
                    color |= (((cur[2] & 64) ? ((cur[2] & 15) + 9) : (cur[2] & 15)) << 4);
                    color |= (((cur[3] & 64) ? ((cur[3] & 15) + 9) : (cur[3] & 15)) << 8);
                    color |= (((cur[4] & 64) ? ((cur[4] & 15) + 9) : (cur[4] & 15)) << 12);
                    unicode.push_back( color );
                    cur+=4;
                break;
                default:
                texta.glbuffer_indexes[1]++;
                    unicode.push_back(((unsigned int)*cur));
                break;
            }
        }
        l++;
        if (l == utf8strings.getSize()) break;
        unicode.push_back(((unsigned int) '\n') | 0x10000);
	} while(true);
    float* vbuf = new float[texta.glbuffer_indexes[1]<<4];
    float off[2]; ExOp::toZero(off);
    color = 0x00000000;

	int tent = 0;
	off[0] = ((float)utf8strings[tent].position[0]) / (font_op.font >> 1);
	off[1] = ((float)utf8strings[tent].position[1]) / font_op.font;
    for(i=0,cursor=0;cursor<unicode.getSize();cursor++){
        if (unicode[cursor] & 0x10000){
            if (unicode[cursor] & 0x20000){
                color = unicode[cursor] << 16;
            }else{
                switch(unicode[cursor]){
                case ' '|0x10000: off[0]+=2; break;
                case '\t'|0x10000: off[0] += 8; break;
                case '\n'|0x10000:
					tent++;
					off[0] = ((float)utf8strings[tent].position[0]) / (font_op.font >> 1);
					off[1] = ((float)utf8strings[tent].position[1]) / font_op.font;
				break;
                }
            }
        }else{
            vbuf[i] = off[0];
            vbuf[i|1] = off[1];
            vbuf[i|2] = 1.0f;
            if ( ((unicode[cursor] >= 0x3040)&& (unicode[cursor] <= 0x30FF))
               ||((unicode[cursor] >= 0x4E00)&& (unicode[cursor] <= 0x9FFF))
               ||((unicode[cursor] >= 0xFF01)&& (unicode[cursor] <= 0xFF60)) ) off[0]+=4;
            else off[0]+=2;
            ((uint32_t*)vbuf)[i|3] =   color | unicode[cursor];
            vbuf[i | 12] = vbuf[i | 8] = vbuf[i | 4] =vbuf[i];
            vbuf[i | 13] = vbuf[i | 9] = vbuf[i | 5] =vbuf[i|1];
         //   vbuf[i | 14] = vbuf[i |10] = vbuf[i | 6] =vbuf[i|2];
            vbuf[i | 15] = vbuf[i |11] = vbuf[i | 7] =vbuf[i|3];
        //   ((uint32_t*)vbuf)[i | 7 ] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
        //   ((uint32_t*)vbuf)[i | 11] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
        //   ((uint32_t*)vbuf)[i | 15] = 0x00000000 | ((uint32_t*)vbuf)[i | 3];
            i+= 16;
        }
    }
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (texta.glbuffer_indexes[1]<<4), vbuf, GL_STATIC_DRAW);
    delete[](vbuf);
    texta.glbuffer_indexes[1] <<=2;
}

/** \brief Render text
 *
 * \param shaderID
 * \param rect,
 * \param bound
 * \return
 *
 */
void GUIStyle::drawTextMesh(const GuiTextAttribute &texta, const int32_t* obrect, const int32_t* parrect){
	if (texta.glbuffer_indexes[0] == 0) return;
	glUseProgram(ctrl_state.datext_shader);
  //  glUniform1i(glGetUniformLocation(ctrl_state.datext_shader, "tex_color"), 0);
    glUniform1i(glGetUniformLocation(ctrl_state.datext_shader, "tex_height"), 20); // 16);
    glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "def_color"), 1.0f,0.0f,0.0f,1.0f); // 16);
    glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "bgcolor_base"), 0.0f,0.0f,0.0f,0.5f); // 16);
    glUniform1f(glGetUniformLocation(ctrl_state.datext_shader, "bgcolor_factor"), 0.5f); // 16);
    glUniform2f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3]);
    float tpos[2];
    if (textlayout == GUITEXT_LAYOUT_CENTERED){
        tpos[0] = obrect[0] + (obrect[2] >> 1) + font_op.shift[0];
        tpos[1] = obrect[1] + (obrect[3] >> 1) + font_op.shift[1];
    }else{
        tpos[0] = obrect[0] + text_borders[0] + font_op.shift[0];
        tpos[1] = obrect[1] + text_borders[0] + font_op.shift[1];
    }
    if (parrect == NULL){
        glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "boundrect"), ((float)obrect[0] + text_borders[0]), (float)(ctrl_state.curwin->rect[3]-1  - obrect[1] - obrect[3] + text_borders[3]), (float)(obrect[0] + obrect[2] - text_borders[2]), (float)(ctrl_state.curwin->rect[3]-1 - obrect[1] -  text_borders[1]) );
    }else{
        glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "boundrect"), ((float) parrect[0]), (float)(ctrl_state.curwin->rect[3]-1  - parrect[3] - parrect[1]), (float)(parrect[0]  + parrect[2]), (float)(ctrl_state.curwin->rect[3]-1 - parrect[1]) );
        tpos[0] += parrect[0] - parrect[4];
        tpos[1] += parrect[1] - parrect[5];
    }
    glUniform2f(glGetUniformLocation(ctrl_state.datext_shader, "vPos"), tpos[0] , tpos[1]);
   // glUniform2f(glGetUniformLocation(ctrl_state.datext_shader, "tex_offset"), obrect[0] + t_opt.hborder[0], 0.0f);

    glActiveTexture(GL_TEXTURE2); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXTOFFSETS);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXT);
    glBindBuffer(GL_ARRAY_BUFFER, texta.glbuffer_indexes[0]);
    glEnableVertexAttribArray(ATTRIBUTE_POSITION);
    glEnableVertexAttribArray(ATTRIBURE_CHARID);
    glVertexAttribPointer(ATTRIBUTE_POSITION,3, GL_FLOAT, GL_FALSE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(0));
    glVertexAttribPointer(ATTRIBURE_CHARID,4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(sizeof(float)*3));
    glDrawArrays(GL_QUADS, 0, texta.glbuffer_indexes[1]);
    }
/** \brief Render text
 *
 * \param shaderID
 * \param rect,
 * \param bound
 * \return
 *
 */
void GUIStyle::drawTextMeshAlias(const GuiTextAttribute &texta, const int32_t* obrect, const int32_t* parrect){
    if (texta.glbuffer_indexes[1] == 0) return;
	glUseProgram(ctrl_state.sstext_shader);
  //  glUniform1i(glGetUniformLocation(ctrl_state.datext_shader, "tex_color"), 0);
    glUniform1i(glGetUniformLocation(ctrl_state.sstext_shader, "tex_height"), 20); // 16);
    glUniform4f(glGetUniformLocation(ctrl_state.sstext_shader, "def_color"), 1.0f,0.0f,0.0f,1.0f); // 16);
    glUniform4f(glGetUniformLocation(ctrl_state.sstext_shader, "bgcolor_base"), 0.0f,0.0f,0.0f,0.5f); // 16);
    glUniform1f(glGetUniformLocation(ctrl_state.sstext_shader, "bgcolor_factor"), 0.5f); // 16);
    glUniform2f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3]);

    float tpos[2];
    if (textlayout == GUITEXT_LAYOUT_CENTERED){
        tpos[0] = obrect[0] + (obrect[2] >> 1) + font_op.shift[0];
        tpos[1] = obrect[1] + (obrect[3] >> 1) + font_op.shift[1];
    }else{
        tpos[0] = obrect[0] + text_borders[0] + font_op.shift[0];
        tpos[1] = obrect[1] + text_borders[0] + font_op.shift[1];
    }
    if (parrect == NULL){
        glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "boundrect"), ((float)obrect[0] + text_borders[0]), (float)(ctrl_state.curwin->rect[3]-1  - obrect[1] - obrect[3] + text_borders[3]), (float)(obrect[0] + obrect[2] - text_borders[2]), (float)(ctrl_state.curwin->rect[3]-1 - obrect[1] -  text_borders[1]) );
    }else{
        glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "boundrect"), ((float) parrect[0]), (float)(ctrl_state.curwin->rect[3]-1  - parrect[3] - parrect[1]), (float)(parrect[0]  + parrect[2]), (float)(ctrl_state.curwin->rect[3]-1 - parrect[1]) );
        tpos[0] += parrect[0] - parrect[4];
        tpos[1] += parrect[1] - parrect[5];
    }
    glUniform2f(glGetUniformLocation(ctrl_state.datext_shader, "vPos"), tpos[0] , tpos[1]);
   // glUniform2f(glGetUniformLocation(ctrl_state.datext_shader, "tex_offset"), obrect[0] + t_opt.hborder[0], 0.0f);

    glActiveTexture(GL_TEXTURE2);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXTOFFSETS);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXT);
    glBindBuffer(GL_ARRAY_BUFFER, texta.glbuffer_indexes[0]);
    glEnableVertexAttribArray(ATTRIBUTE_POSITION);
    glEnableVertexAttribArray(ATTRIBURE_CHARID);
    glVertexAttribPointer(ATTRIBUTE_POSITION,3, GL_FLOAT, GL_FALSE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(0));
    glVertexAttribPointer(ATTRIBURE_CHARID,4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(sizeof(float)*3));
    glDrawArrays(GL_QUADS, 0, texta.glbuffer_indexes[1]);
    }
void GUIStyle::setToAreaDefault(){
	setBorderSizes(12,12);
}

void GUIStyle::setToMenuDefault(){
	setBorderSizes(64,64);
}

	MyWindow::MyWindow(unsigned int alias, unsigned int style_alias, int sizex,int sizey, RELPOS_enum posis, bool full, int nCmdShow,const char* winname): GUIArea(alias,style_alias),isfullscr(full),last_mouse_depth(0.0){

    SDL_GL_SetAttribute(SDL_GL_RED_SIZE,        8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,      8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,       8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,      8);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,      24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE,    8);
    SDL_GL_SetAttribute(SDL_GL_BUFFER_SIZE,        32);
   /* SDL_GL_SetAttribute(SDL_GL_ACCUM_RED_SIZE,    8);
    SDL_GL_SetAttribute(SDL_GL_ACCUM_GREEN_SIZE,    8);
    SDL_GL_SetAttribute(SDL_GL_ACCUM_BLUE_SIZE,    8);
    SDL_GL_SetAttribute(SDL_GL_ACCUM_ALPHA_SIZE,    8);*/
  //  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS,  1);
  //  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,  2);
    rect[0] = 0;
    rect[1] = 0;
    rect[4] = 0;
    rect[5] = 0;

//SDL_WINDOW_FULLSCREEN |

    SDL_DisplayMode current;

    int i,px,py;
    for(i = 0; i < SDL_GetNumVideoDisplays(); ++i){   // Get current display mode of all displays.
        if (SDL_GetCurrentDisplayMode(i, &current) == 0){
            switch( (((int)posis) & 7) ){
                case 1: px = 0; break;
                case 2: px = current.w - sizex; break;
                default: px = (current.w - sizex) >> 1;
            }
            switch( ((((int)posis) >> 3 )& 7)){
                case 1: py = 0; break;
                case 2: py = current.h - sizey; break;
                default: py = (current.h - sizey) >> 1;
            }
            break;
        }
    }
    if (i >= SDL_GetNumVideoDisplays()){px = SDL_WINDOWPOS_UNDEFINED; py = SDL_WINDOWPOS_UNDEFINED;}

    char bufname[256];
    if (winname) strcpy(bufname, winname);
    else strcpy(bufname, "untitled");


    if((Surf_Display = SDL_CreateWindow(bufname,px,py, sizex, sizey, SDL_WINDOW_OPENGL )) == NULL) {
        fprintf(stderr, "Failed to create Window!\n");
        fprintf(stderr, "SDL error: %s\n", SDL_GetError());

        Surf_Display = SDL_CreateWindow("gamewin",SDL_WINDOWPOS_UNDEFINED,SDL_WINDOWPOS_UNDEFINED, sizex, sizey, 0 );
        printf("Can open a non-opengl window? %c\n", Surf_Display == NULL ? 'N' : 'Y');
        exit(1);
    }

    int tmpdim[2];
    SDL_GetWindowSize(Surf_Display,tmpdim,tmpdim+1);
    rect[6] = rect[2] = tmpdim[0];
    rect[7] = rect[3] = tmpdim[1];

   // if (NULL == (renderer =SDL_CreateRenderer(Surf_Display, -1, 0))) exit(1);
    if ((glcontext = SDL_GL_CreateContext(Surf_Display)) == NULL) myexit(SDL_GetError(void));

    SDL_GL_SetSwapInterval(1);

    /*if(!gladLoadGL()) {
        printf("Something went wrong!\n");
        exit(-1);
    }*/


	if (glClearColor == NULL) printf("function is still not linked...!!!\n");
  //  SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
  //  SDL_RenderClear(renderer);
  //  SDL_RenderPresent(renderer);
 // SDL_HWSURFACE | SDL_DOUBLEBUF | SDL_OPENGL
 //   SDL_EnableUNICODE(1);
    GLfloat light0Pos[4] = { 0.0F, 0.0F, 1.0F, 1.0F };
    glClearColor(0.5f, 0.0f, 0.75f, 1.0F);
    glClearIndex((GLfloat) 0.0f);
    setProjection();
    glTranslatef(0.0F, 0.0F, -0.5F);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.00f);
	glCullFace(GL_BACK);
    glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
    glEnable(GL_LIGHT0);
    glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.4f);
    if (!colorIndexMode) {
	glEnable(GL_COLOR_MATERIAL);
    }
    glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glEnable(GL_TEXTURE_2D);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);       // Blending Function For Translucency Based On Source Alpha
	glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glColor4f(1,1,1,1);
    // shader!
    //last_mouse_alias=0;
//  wind_bit->update(sizex,sizey,rgba_layer.data);

/*
void glFramebufferParameteri?(GLenum target?, GL_FRAMEBUFFER_DEFAULT_WIDTH, GLint param?);
target? is the location where the framebuffer object is bound. To set the width, set pname? to ; to set the height, use GL_FRAMEBUFFER_DEFAULT_HEIGHT.

*/
    framebufferID[0] = 0;
    framebufferID[1] = 0;
    framebufferID[2] = 0;

    if (glGenFramebuffers != NULL) {
    if (glFramebufferTexture2D  == NULL) myexit("No extern frame buffer texture!");
        glGenFramebuffers(3, framebufferID);
        printf("Cur Error Code %i\n", (int)glGetError());
        for(int i =0 ; i < 3; i++){
            if (framebufferID[i] == 0) {printf("Could not allocate framebuffer no%i\n", i); ExOp::show(framebufferID); exit(1);}
        }
    }

//    glGenTextures(1,&renderbufferID);

    last_mouse_ddepth[0] =0.0; // happy valgrind!
    last_mouse_ddepth[1] =0.0; // happy valgrind!
    memset(depth, '\0', 9 * sizeof(float)); // happy valgrind!
    memset(stencilval, '\0', 3 * sizeof(GLint)); // happy valgrind!


   // glGenRenderbuffers(1,&renderbufferID);

	} //         	if (pcol.w &lt; 0.5) discard;
/*
        	col *= (1.0 - fact);
        	col += pcol * fact;

        	if (col.w &lt; 0.0625f) discard;

        	col.xyz *= (1.0 - fact);
        	col.xyz += pcol.xyz * fact;
*/
MyWindow::~MyWindow(){
    SDL_GL_DeleteContext(glcontext);
    SDL_DestroyWindow(Surf_Display);
//    glDeleteTextures(1,&renderbufferID);
    if (glGenFramebuffers != NULL) glDeleteFramebuffers(3, framebufferID);
}
void MyWindow::setWindowName(const char* newname){
    // SDL_WM_SetCaption(newname,newname); sdl1.2
}
void MyWindow::resize(void){
   setProjection();
    glViewport(0, 0, rect[2], rect[3]);
}
void MyWindow::setProjection(void){
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glMatrixMode(GL_MODELVIEW);glLoadIdentity();
}
void MyWindow::operator<<(renderMode* _nrend){
	render_list.push_back(_nrend);
}
void MyWindow::operator>>(renderMode* _nrend){
    unsigned int i;
    for(i=0;i<render_list.getSize();i++) if (render_list[i] == _nrend ) break;
    if (render_list.getSize()) render_list.pop_swap(i);
}
void MyWindow::render(){

//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
    unsigned int i;




    // reserved stancils
    // 0: uninitialized
    // 1: selected
    // 2: gui
    //

//
/*     glBindFramebuffer(GL_FRAMEBUFFER, framebufferID[2]);
glFramebufferTexture2D(
GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, renderbufferID, 0);
glDrawBuffer(GL_NONE);*/
/*     glPixelStorei( GL_PACK_ROW_LENGTH, rect[2] );
    glPixelStorei( GL_PACK_SKIP_ROWS, 0 );
    glPixelStorei( GL_PACK_SKIP_PIXELS, 0 );
    glPixelStorei( GL_PACK_ALIGNMENT, 1 );*/

//      glBindRenderbuffer(GL_RENDERBUFFER, renderbufferID);
//       glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, rect[2], rect[3]);
//      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, renderbufferID);
//    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );


    //for(i=0;i<render_list.size();i++) render_list[i]->draw(this);
    for(i=0;i<render_list.size();i++) {
        if (render_list[i] == NULL) printf("Trying to render a sub that does not exist...\n");
        else {
            LFHBreakpoint();
            render_list[i]->draw(this);
        }
    }

    glMatrixMode(GL_PROJECTION);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    glEnable(GL_STENCIL_TEST);
    glLoadIdentity();
    glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
    glTranslatef(-1.0f,1.0f,0.0f);
    glScalef(2.0f / (rect[2]),-2.0f / (rect[3]),1.0f);
    glMatrixMode(GL_MODELVIEW);glLoadIdentity();

	//GLint vmaj,vmin;
	//glGetIntegerv(GL_MAJOR_VERSION,&vmaj);
	//glGetIntegerv(GL_MINOR_VERSION,&vmin);

//-printf("OpenGl version %i.%i    %s\n", glGetString(GL_VERSION));

    glDepthFunc(GL_ALWAYS);glDisable(GL_CULL_FACE);

    //ctrl_state.alias_storm = clock() & 0x0400;
    ctrl_state.alias_storm = false;
	if (ctrl_state.alias_storm){
        glDisable( GL_MULTISAMPLE );glBlendFunc(GL_ONE,GL_ZERO);
        this->drawAliasSubs(false);
        if (ctrl_state.foreground_alias != 0){
            ctrl_state.gui_objects_ptr[ctrl_state.foreground_alias]->drawAlias(ctrl_state.gui_objects_ptr[ ctrl_state.gui_objects_ptr[ctrl_state.foreground_alias].parent_alias]->getRect() );
        }
        glEnable( GL_MULTISAMPLE );glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	}else{
        this->drawSubs();
        if (ctrl_state.foreground_alias != 0){
            ctrl_state.gui_objects_ptr[ctrl_state.foreground_alias]->draw( ((ctrl_state.mouse_stencil& 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == ctrl_state.foreground_alias) , ctrl_state.gui_objects_ptr[ ctrl_state.gui_objects_ptr[ctrl_state.foreground_alias].parent_alias]->getRect() );
        }
	}




	glDepthFunc(GL_LEQUAL);glEnable(GL_CULL_FACE);
    for(i=0;i<render_list.size();i++) render_list[i]->drawPostGUI(this);



//    if (wind_bit) wind_bit->draw(0,0);
     while(true){
        GLenum curerr;
        curerr = GL_NO_ERROR;
        LFH_VALGRIND_MUTE(curerr = glGetError();)
        if (curerr == GL_NO_ERROR) break;
        switch(curerr){
            case GL_INVALID_ENUM:printf("An unacceptable value is specified for an enumerated argument.\n"); break;
            case GL_INVALID_VALUE: printf("A numeric argument is out of range.\n");break;
            case GL_INVALID_OPERATION: printf("The specified operation is not allowed in the current state.\n");break;
            case GL_STACK_OVERFLOW: printf("This command would cause a stack overflow.\n");break;
            case GL_STACK_UNDERFLOW: printf("This command would cause a stack underflow.\n");break;
            case GL_OUT_OF_MEMORY: printf("There is not enough memory left to execute the command. The state of the GL is undefined, except for the state of the error flags, after this error is recorded.\n");break;
        //    case GL_TABLE_TOO_LARGE: printf("The specified table exceeds the implementation's maximum supported table size.\n");break;
            default:
                printf("GL error with code %i.\n", (int)curerr);break;
        }
     }


//    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);


    last_mouse_depth = depth[4];
    for(i=0;i<9;i++) depth[i] = i;

    unsigned int j = rect[3] - ctrl_state.mouse_coor[1]-1;
    i = ctrl_state.mouse_coor[0];
    if (i == (uint32_t)(rect[2]-1)) i-=2;
    else if (i !=0u ) i--;
    if (j == (uint32_t)(rect[3]-1))j-=2;
    else if (j !=0u ) j--;


    //ExOp::show(depth);
    glReadPixels(i++, j++, 3, 3, GL_DEPTH_COMPONENT, GL_FLOAT, depth);
    //printf("mmm %e\n", last_mouse_depth);
    //ExOp::show(depth);

    glReadPixels(i, j, 1, 1, GL_STENCIL_INDEX, GL_INT, stencilval);


    last_mouse_ddepth[0] = (2.0f *(depth[5] - depth[3]) + depth[2] + depth[8] - depth[0] - depth[6]) /8.0f;
    last_mouse_ddepth[1] = (2.0f *(depth[7] - depth[1]) + depth[6] + depth[8] - depth[0] - depth[2]) /8.0f;
    bool tmp_storm;
    if ((ctrl_state.button_state & 7) == 0) {
        this->renderAlias();
    }else{ // some mouse button is down, anchor to current object

        if (((ctrl_state.mouse_stencil & 0xFE) == 2)&&(stencilval[0] == 3)){ // mouse is on some text, get offset
            glEnable(GL_SCISSOR_TEST); glDisable( GL_MULTISAMPLE );glBlendFunc(GL_ONE,GL_ZERO);
            glScissor(i,j,1,1);
            glClear(GL_DEPTH_BUFFER_BIT);


            glDisable(GL_CULL_FACE);glDepthFunc(GL_ALWAYS);
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glScalef(1.0f,-1.0f,1.0f);
            glTranslatef(-1.0f,1.0f,0.0f);
            glScalef(2.0f / (rect[2]),2.0f / (rect[3]),0.0f);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();


            ctrl_state.gui_objects_ptr[ctrl_state.mouse_stencilID]->drawAlias(ctrl_state.gui_objects_ptr[ ctrl_state.gui_objects_ptr[ctrl_state.mouse_stencilID].parent_alias]->getRect() );

            glDepthFunc(GL_LEQUAL);glEnable(GL_CULL_FACE);
            //glReadBuffer(GL_COLOR_ATTACHMENT0);
            glReadPixels(i, j, 1, 1,  GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, &ctrl_state.mouse_array_offset);
            glEnable( GL_MULTISAMPLE );glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);    glDisable(GL_SCISSOR_TEST);
        }
    }
    SDL_GL_SwapWindow(Surf_Display);
}
void MyWindow::renderAlias(){
    uint32_t posy = rect[3] - ctrl_state.mouse_coor[1]-1;
    if (stencilval[0] > 1u){
        bool tmp_storm = ctrl_state.alias_storm; ctrl_state.alias_storm  = false;
        glEnable(GL_SCISSOR_TEST); glDisable( GL_MULTISAMPLE );glBlendFunc(GL_ONE,GL_ZERO);
        glScissor(ctrl_state.mouse_coor[0],posy,1,1);
        glClear(GL_DEPTH_BUFFER_BIT);

        if ((stencilval[0] & 0xFE) == 2u){
            glDisable(GL_CULL_FACE);glDepthFunc(GL_ALWAYS);
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glScalef(1.0f,-1.0f,1.0f);
            glTranslatef(-1.0f,1.0f,0.0f);
            glScalef(2.0f / (rect[2]),2.0f / (rect[3]),0.0f);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            this->drawAliasSubs(stencilval[0] == 3);

            if (ctrl_state.foreground_alias != 0){
                ctrl_state.gui_objects_ptr[ctrl_state.foreground_alias]->drawAlias(stencilval[0] == 3, ctrl_state.gui_objects_ptr[ ctrl_state.gui_objects_ptr[ctrl_state.foreground_alias].parent_alias]->getRect() );
            }
            glDepthFunc(GL_LEQUAL);glEnable(GL_CULL_FACE);
        }else{
            for(int k=0;k<render_list.size();k++) render_list[k]->drawAlias(this);
        }
        glReadPixels(ctrl_state.mouse_coor[0], posy, 1, 1,  GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, stencilval+2);
        glReadPixels(ctrl_state.mouse_coor[0], posy, 1, 1, GL_STENCIL_INDEX, GL_INT, stencilval+3);
        ctrl_state.notifyStencilChange(stencilval[0], stencilval+2,true);
        ctrl_state.alias_storm = tmp_storm;
        glEnable( GL_MULTISAMPLE );glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);    glDisable(GL_SCISSOR_TEST);
    }else{
        stencilval[2] = 0;
        if ((stencilval[0] == 0)&&(ctrl_state.mouse_stencil != 0)) ctrl_state.notifyStencilChange(0, stencilval+2, false);
    }
}



void MyWindow::render_static(){/*
    LFHPrimitive::ExOp::toZero(rgba_layer);
    Tuple<unsigned int, 2> pos;
    Uint32 color = 0xFF00FFFF;
    pos[0] =30;
    pos[1] = 40;
    rgba_layer.drawChar('A', pos,color);
    pos[0] =60;
    pos[1] = 40;
    rgba_layer.drawCharFlip('A', pos,color);
    GUIObject* tmp_subs;
    for(pos[0] = 0;pos[0]< subs.size();pos[0]++){
        tmp_subs = ctrl_state.gui_objects_ptr[subs[pos[0]]];
        if (tmp_subs == NULL) exit(1);
        tmp_subs->draw_static(rgba_layer,alias_layer);
        }

    wind_bit->update(rect[2],rect[3],rgba_layer.data);*/
/*
    LFHPrimitive::ExOp::toZero(alias_layer);
    msg.type = LFHGUI_ALIAS_DRAW;
    msg.target = &alias_layer;
	for(unsigned int i=0;i<gui_alias_todraw.size();i++) {msg.key = gui_alias_todraw[i];static_guiscope.GUI_msg_base.dispatch(msg);}
*/
    }
/*
void MyWindow::unloadGUI(GUID_Id what){
    for(unsigned int i=0;i<subs.size();i++) if (subs[i] == (Uint32)what) subs.pop_swap(i);
    GUIObject* cur = ctrl_state.gui_objects_ptr[(Uint32)what];
    delete(cur);
    }*/
GLvoid MyWindow::glPrint(const char *fmt, ...)    // Custom GL "Print" routine
{
	char text[256];	                // Holds our string
	va_list ap;	        	// Pointer to list of arguments
	if (fmt == NULL)		// If there's no text
		return;			// Do nothing
	va_start(ap, fmt);		// Parses the string for variables
	vsprintf(text, fmt, ap);	// And converts symbols to actual numbers
	va_end(ap);			// Results are stored in text
	glPushAttrib(GL_LIST_BIT);	// Pushes the display list bits
//	glListBase(base - 32);		// Sets the base character to 32
	glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);	// Draws the display list text
	glPopAttrib();						// Pops the display list bits
}
void MyWindow::idleFunc(){
	render();
	}
/*
	string MyWindow::openfilename(char *filter, HWND owner) {
  OPENFILENAME ofn;
  char fileName[MAX_PATH] = "";
  ZeroMemory(&ofn, sizeof(ofn));
  ofn.lStructSize = sizeof(OPENFILENAME);
  ofn.hwndOwner = owner;
  ofn.lpstrFilter = filter;
  ofn.lpstrFile = fileName;
  ofn.nMaxFile = MAX_PATH;
  ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
  ofn.lpstrDefExt = "";
  string fileNameStr;
  if ( GetOpenFileName(&ofn) )
     fileNameStr = fileName;
  return fileNameStr;
}

string MyWindow::savefilename(char *filter, HWND owner) {
  OPENFILENAME ofn;
  char fileName[MAX_PATH] = "";
  ZeroMemory(&ofn, sizeof(ofn));
  ofn.lStructSize = sizeof(OPENFILENAME);
  ofn.hwndOwner = owner;
  ofn.lpstrFilter = filter;
  ofn.lpstrFile = fileName;
  ofn.nMaxFile = MAX_PATH;
  ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
  ofn.lpstrDefExt = "";
  string fileNameStr;
  if ( GetSaveFileName(&ofn) )
     fileNameStr = fileName;
  return fileNameStr;
}
*/

/*
GUIArea::GUIArea(Vector<GUIObject> &where,unsigned int* p_drect, unsigned int *p_trect){
	memcpy(trect,p_trect,sizeof(trect));
	//glGenBuffers(1,&buf_id);
	where.push_back(GUIObject(LFHGUI_TYPE_FRAME, (void*)this,p_drect));
	}
*/

void GUIArea::setAreaDimentions(unsigned int width, unsigned int height, bool is_equal_to_real){
    rect[6] = width;
    rect[7] = height;
    if (is_equal_to_real){
        rect[2] = rect[6];
        rect[3] = rect[7];
        rect[4] = 0;
        rect[5] = 0;
    }else{
        if (rect[4] >= rect[6] - rect[2]) rect[4] = rect[6] - rect[2]-1;
        if (rect[5] >= rect[7] - rect[3]) rect[5] = rect[7] - rect[3]-1;
    }
}

void GUIArea::insertGUI(Uint32 a, RELPOS_enum relpos_type, unsigned int pad_x, unsigned int pad_y){
    GUIObject* obj = ctrl_state.gui_objects_ptr[a].target;
    subs.push_back(a);
    ctrl_state.gui_objects_ptr[a].parent_alias = this->GUI_alias;
    GUIStyle& istyle = ctrl_state.gui_styles[obj->styleID];
    if (istyle.dimentions[0] != 0) obj->setDimentions(istyle.dimentions[0],istyle.dimentions[1]);

    int32_t* crect = obj->accessRect();
    switch(((int)relpos_type) & 7){
        case 0: // centered
            crect[0] = ((rect[2] - crect[2]) >> 1) + pad_x;
        break;
        case 1: // left
            crect[0] = pad_x;
        break;
        case 2: // right
            crect[0] = (rect[2] - crect[2]-1) - pad_x;
        break;
    }
    switch(((int)relpos_type) & 56){
        case 0: // centered
            crect[1] = ((rect[3] - crect[3]) >> 1) + pad_y;
        break;
        case 8: // left
            crect[1] = pad_y;
        break;
        case 16: // right
            crect[1] = (rect[3] - crect[3]-1)-pad_y;
        break;
    }
    obj->onResize();
}
void GUIArea::removeGUI(Uint32 a){ for(unsigned int i=0;i< subs.size();i++) if (subs[i] == a) subs.pop_swap(i); ctrl_state.gui_objects_ptr[a].parent_alias = a;}




GUIObject::GUIObject(unsigned int alias, unsigned int style_alias): GUI_alias(alias),gui_flags(0), stare_alias(0), styleID(style_alias)  {ctrl_state.gui_objects_ptr[alias].target = this; ctrl_state.gui_objects_ptr[alias].parent_alias = alias;}
GUIObject::GUIObject(unsigned int alias, unsigned int style_alias,const char* _text): GUI_alias(alias),gui_flags(0),stare_alias(0), styleID(style_alias) {ctrl_state.gui_objects_ptr[alias].target = this;this->setText(_text); ctrl_state.gui_objects_ptr[alias].parent_alias = alias;}
void GUIObject::initialize_routine(uint32_t alias, uint32_t style_alias){
    GUI_alias = alias;
    styleID = style_alias;
    gui_flags = 0;
    stare_alias = 0;
}
void GUIObject::setDimentions(unsigned int width, unsigned int height){
    int32_t* rect = this->accessRect();
    rect[2] = width;rect[3] = height;
}
void GUIObject::setPos(const int32_t* pos){
	int32_t* rect = this->accessRect();
	rect[0] = pos[0];rect[1] = pos[1];
}
void GUIObject::setRect(const int32_t* _rect){
	int32_t* rect = this->accessRect();
	rect[0] = _rect[0];rect[1] = _rect[1];rect[2] = _rect[2];rect[3] = _rect[3];
}
void GUIObject::wrPosition(int32_t *fout) const{
    const int32_t * rect;
    uint32_t par = ctrl_state.gui_objects_ptr[GUI_alias].parent_alias;
    if (par == GUI_alias){
        rect = this->getRect();
        fout[0] = rect[0];
        fout[1] = rect[1];
    }else{
        ctrl_state.gui_objects_ptr[par]->wrSubPosition(fout, this);
    }
}
void GUIObject::setPositionRelativeTo(const Guialias& other_alias, RELPOS_enum relpos_type, unsigned int pad_x, unsigned int pad_y){
	int32_t* rect = this->accessRect();
	const int32_t* otherrect = ctrl_state.gui_objects_ptr[other_alias]->getRect();
	switch(((int)relpos_type) & 7){
		case 0: // centered
			rect[0] = ((otherrect[2] - rect[2]) >> 1) + pad_x;
		break;
		case 1: // left
			rect[0] = pad_x;
		break;
		case 2: // right
			rect[0] = (otherrect[2] - rect[2]-1) - pad_x;
		break;
		case 4: // centered
			rect[0] = otherrect[0] + pad_x;
		break;
		case 5: // before
			rect[0] = otherrect[0] - rect[2] - pad_x;
		break;
		case 6: // after
			rect[0] = otherrect[0] + otherrect[2] + pad_x;
		break;
	}
	switch(((int)relpos_type) & 56){
		case 0: // centered
			rect[1] = ((otherrect[3] - rect[3]) >> 1) + pad_y;
		break;
		case 8: // left
			rect[1] = pad_y;
		break;
		case 16: // right
			rect[1] = (otherrect[3] - rect[3]-1)-pad_y;
		break;
		case 32: // centered
			rect[1] = otherrect[1] + pad_y;
		break;
		case 40: // before
			rect[1] = otherrect[1] - rect[3]-pad_y;
		break;
		case 48: // after
			rect[1] = otherrect[1] + otherrect[3]+pad_y;
		break;
	}
	rect[0] += otherrect[0];
	rect[1] += otherrect[1];
	this->onResize();
}
void GUIObject::setVisible(bool is_visible){
	if ((is_visible)^((gui_flags & 0x80000000)==0)){
		gui_flags ^= 0x80000000;
		if (is_visible) this->manifest();
		else this->vanish();
	}
}
void GUIObject::setEnabled(bool isEnabled){
    if (isEnabled) gui_flags &= 0xBFFFFFFF;
    else gui_flags |= 0x40000000;
}
uint32_t GUIObject::getArrayIDfromMouse() const{return 0;}

void GUIArea::wrSubPosition(int32_t *fout, const GUIObject* subptr) const{
    fout[0] = rect[0] + subptr->getRect()[0];
    fout[1] = rect[1] + subptr->getRect()[1];
}
void GUIArea::update(){}
void GUIArea::draw(bool mouse_over, const int32_t* par_rect){
	if (this->gui_flags & 0x80000000) return;


	GLuint dashader = ctrl_state.daframe_shader;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];

	glUseProgram(dashader);

	bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
	if (is_over){
		glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
		gui_style.setColorUniforms(dashader,1);
	}else{
	    gui_style.setColorUniforms(dashader,0);
	}


//	glUniform1f(glGetUniformLocation(ctrl_state.daframe_shader, "bordersize"),gui_style.border_width);
//      glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color1"), 0.24f,1.0,0.5,0.0f); // outside



   // glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color5"), 0.62f,1.0,0.5,1.0f);
   // glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color6"), 0.74f,1.0,0.5,1.0);
   // glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color7"), 0.86f,1.0,0.5,1.0f);


	glActiveTexture(GL_TEXTURE2);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
	glActiveTexture(GL_TEXTURE1);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
	glActiveTexture(GL_TEXTURE0);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
	glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
	glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
	glUniform4f(glGetUniformLocation(dashader, "boundtexcoor"), 0.0f,0.0f,1.0f,1.0f);
    glUniform2f(glGetUniformLocation(dashader, "borderUV"), ((GLfloat)gui_style.border_width) / rect[2] , ((GLfloat)gui_style.border_width) / rect[3]);

	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
	if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
	//glUseProgram(daicon_shader);
	return this->drawSubs();
}
void GUIArea::drawSubs(){
	unsigned int i;
	GUIObject* tmpob;

	for(i=0;i<subs.size();i++) {

     //  if ((ctrl_state.mouse_stencil & 0xFE) == 2u) printf("%X, ",  subs[i]);
		if (ctrl_state.foreground_alias == subs[i]) continue;
		if (ctrl_state.tooltipalias == subs[i]) glDisable(GL_STENCIL_TEST);
		tmpob = ctrl_state.gui_objects_ptr[subs[i]].target;
		if (!tmpob) exit(1);

		tmpob->draw( ((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == subs[i]), rect);
		if (ctrl_state.tooltipalias == subs[i]) glEnable(GL_STENCIL_TEST);
	}
	//if ((ctrl_state.mouse_stencil & 0xFE) == 2u) printf("== %X?\n", ctrl_state.mouse_stencilID);
}
void GUIArea::drawAlias(bool is_text, const int32_t* par_rect){
	if (this->gui_flags & 0x80000000) return;


    GLuint dashader = ctrl_state.daframe_alias_shader;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    glUseProgram(dashader);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else {
        glUniform4f(glGetUniformLocation(dashader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
        glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    }

    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);


    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),gui_style.border_width);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
    glUniform4f(glGetUniformLocation(dashader, "boundtexcoor"), 0.0f,0.0f,1.0f,1.0f);
    glUniform2f(glGetUniformLocation(dashader, "borderUV"), ((GLfloat)gui_style.border_width) / rect[2] , ((GLfloat)gui_style.border_width) / rect[3]);

	glBegin(GL_QUADS); glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();

	return this->drawAliasSubs(is_text);
}

void GUIArea::drawAliasSubs(bool is_text){
    unsigned int i;
    GUIObject* tmpob;
	for(i=0;i<subs.size();i++) {
		if (ctrl_state.foreground_alias == subs[i]) continue;
		if (ctrl_state.tooltipalias == subs[i]) continue;
		tmpob = ctrl_state.gui_objects_ptr[subs[i]].target;
		if (!tmpob) exit(1);
		tmpob->drawAlias(is_text, rect);
	}
}

        //void GUIArea::drawstatic(){
        //  printf("stdr %i\t%i\n",innerect[2], innerect[3] ); fflush(stdout);
        //    if (!(wind_bit)) wind_bit = new BitmapRessource();
        //    glBindFramebuffer(GL_FRAMEBUFFER, wo->framebufferID[0]);
        //    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
       /*   if (text){
                char* p;
                unsigned int pos[2]; pos[0] = rect[0]; pos[1] =rect[1];
                if (text){
                for(p = cur_str; *p != '\0';p++){
                    if ((!((*p) & 0x80 )) && ((*p) & 0x60 )) {fontres->draw_Vslice(pos[0],pos[1], (((*p) - 32) << 4),16); pos[0] += 13;}
                    }
                }
            }*/
        //   printf("done? %c\n", (wind_bit) ? 'Y' : 'N'); fflush(stdout);
        //    wind_bit->grab(innerect[2], innerect[3]);
        //    glBindFramebuffer(GL_FRAMEBUFFER, 0);
        //}
ERRCODE GUIArea::manifest(){ERRCODE fout=0;
    unsigned int i;
    GUIObject* tmpob;
    for(i=0;i<subs.size();i++) {tmpob = ctrl_state.gui_objects_ptr[subs[i]].target; if (!tmpob) exit(1); fout |= tmpob->manifest();}
    // ctrl_state.latentqueue.insert_async(0, new GUITaskEvent(GUI_alias, GUITASK_DRAW_STATIC));
    return fout;
}
ERRCODE GUIArea::vanish(){ERRCODE fout=0;
	unsigned int i;
	GUIObject* tmpob;
	for(i=0;i<subs.size();i++) {tmpob = ctrl_state.gui_objects_ptr[subs[i]].target; if (!tmpob) exit(1); fout |= tmpob->vanish();}
	return fout;
}

GUIMSG_Enum GUIArea::processGUIevent(const GUIMSG_Enum event){
    return GUIMSG_NULL;
}


            /*
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    unsigned int bsize=8;
    LFHPrimitive::RessourcePtr<TextureRessource> text_gui; text_gui = TEX2_GU_SBARS;
    text_gui->use();
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f,0.75f); glVertex2i(bsize,0);
    glTexCoord2f(0.0f,1.0f);glVertex2i(bsize,bsize-1);
    glTexCoord2f(0.0625f * ((rect[2] >> 1) - bsize),1.0f);glVertex2i(rect[2]-1-bsize,bsize-1);
    glTexCoord2f(0.0625f * ((rect[2] >> 1) - bsize),0.75f);glVertex2i(rect[2]-1-bsize,0);
	glTexCoord2f(0.0f,0.75f); glVertex2i(bsize,rect[3]-bsize);
    glTexCoord2f(0.0f,1.0f);glVertex2i(bsize,rect[3]-1);
    glTexCoord2f(0.0625f * ((rect[2] >> 1) - bsize),1.0f);glVertex2i(rect[2]-1-bsize,rect[3]-1);
    glTexCoord2f(0.0625f * ((rect[2] >> 1) - bsize),0.75f);glVertex2i(rect[2]-1-bsize,rect[3]-bsize);
	glTexCoord2f(0.0f,0.75f); glVertex2i(0,bsize);
    glTexCoord2f(0.0625f * ((rect[3] >> 1) - bsize),0.75f);glVertex2i(0,rect[3]-1-bsize);
    glTexCoord2f(0.0625f * ((rect[3] >> 1) - bsize),1.0f);glVertex2i(bsize-1,rect[3]-1-bsize);
    glTexCoord2f(0.0f,1.0f);glVertex2i(bsize-1,bsize);
	glTexCoord2f(0.0f,0.75f); glVertex2i(rect[2]-bsize,bsize);
    glTexCoord2f(0.0625f * ((rect[3] >> 1) - bsize),0.75f);glVertex2i(rect[2]-bsize,rect[3]-1-bsize);
    glTexCoord2f(0.0625f * ((rect[3] >> 1) - bsize),1.0f);glVertex2i(rect[2]-1,rect[3]-1-bsize);
    glTexCoord2f(0.0f,1.0f);glVertex2i(rect[2]-1,bsize);
    glEnd();
    glBindTexture(GL_TEXTURE_2D,0);
	glBegin(GL_QUADS);
//	glColor3f(0.0f,0.0f,1.0f);glVertex2i(0,0);
 //   glColor3f(1.0f,1.0f,0.0f);glVertex2i(0,rect[3]-1);
 //   glColor3f(1.0f,0.0f,0.0f);glVertex2i(rect[2]-1,rect[3]-1);
//    glColor3f(0.0f,1.0f,1.0f);glVertex2i(rect[2]-1,0);
	glColor3f(0.0f,0.0f,1.0f);glVertex2i(bsize,bsize);
    glColor3f(1.0f,1.0f,0.0f);glVertex2i(bsize,rect[3]-1-bsize);
    glColor3f(1.0f,0.0f,0.0f);glVertex2i(rect[2]-1-bsize,rect[3]-1-bsize);
    glColor3f(0.0f,1.0f,1.0f);glVertex2i(rect[2]-1-bsize,bsize);
    glEnd();
            for(i=0;i<subs.size();i++) {tmpob = ctrl_state.gui_objects_ptr[subs[i]]; tmpob->drawstatic();}

            if (!(wind_bit)) wind_bit = new BitmapRessource();
            wind_bit->grab(innerect[2], innerect[3]);

		void GUIArea::draw_static(LFHPrimitive::DataGrid<unsigned char, 3> &buff) const{

                coor[0] =4;
                coor[1] = innerect[2];
                coor[2] = innerect[3];
                buff.setSizes(coor);
                if (bg_texture == bg_texture+10){
                }else {
                    coor[0] =4;
                    coor[1] = 4;
                    coor[2] = 4;
                    bgtex.setSizes(coor);
                    LFHPrimitive::ExOp::toRand(bgtex);
                }
                for(coor[2]=0;coor[2]<buff.dims[2];coor[2]++){
                    altcoor[2] = coor[2] % bgtex.dims[2];
                for(coor[1]=0;coor[1]<buff.dims[1];coor[1]++){
                    altcoor[1] = coor[1] % bgtex.dims[1];
                    for(coor[0]=0;coor[0]<4;coor[0]++) {altcoor[0] = coor[0];buff(coor) = bgtex(altcoor);}
                }
                }

                // draw Borders

        }
        void GUIArea::draw_alias(LFHPrimitive::DataGrid<pair<unsigned int,unsigned int>, 2> &buff) const{
            LFHPrimitive::Tuple<unsigned int, 2> coor;
            LFHPrimitive::Tuple<unsigned int, 3> altcoor;
                coor[0] = innerect[2];
                coor[1] = innerect[3];
            buff.setSizes(coor);
            pair<unsigned int,unsigned int> input;
            input.first = GUI_alias;
            input.second = 0;
            for(coor[1]=0;coor[1]< buff.dims[1];coor[1]++)
            for(coor[0]=0;coor[0]< buff.dims[0];coor[0]++){
                buff(coor) = input;
            }
        }*/
/*
void GUIArea::draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map)const{
    LFHPrimitive::DataGrid<Uint32, 2> rgba_map_tmp;
    LFHPrimitive::DataGrid<Uint32, 2> alias_map_tmp;
    Tuple<unsigned int, 2> dims; dims[0] = innerect[2]; dims[1] = innerect[3]; rgba_map_tmp.setSizes(dims); alias_map_tmp.setSizes(dims);
    ExOp::toZero(rgba_map_tmp);
    Tuple<unsigned int, 2> coor;
    ExOp::toRand(rgba_map_tmp); ExOp::toRand(rgba_map_tmp);
    for(dims[1] = 0; dims[1] < innerect[3];dims[1]++)
        for(dims[0] = 0; dims[0] < innerect[2];dims[0]++){
        alias_map_tmp(dims) = GUI_alias;
        }
    unsigned int i;
    GUIObject* tmpob;
    for(i=0;i<subs.size();i++) {tmpob = ctrl_state.gui_objects_ptr[subs[i]]; if (!tmpob) exit(1); tmpob->draw_static(rgba_map_tmp,alias_map_tmp);}
    LFHPrimitive::Tuple<unsigned int, 3> altcoor;
    Uint32 dapix;
   // if (area_state & GUIOBJECT_STATE_HAS_BORDERS){
                LFHPrimitive::TiffFile tf("d:/prog/Display/Images/guisbars2.tif");
// area_state | GUIOBJECT_STATE_SCROLLBAR
// area_state | GUIOBJECT_STATE_CLOSE_BUTTON
                LFHPrimitive::DataGrid<unsigned char, 3> bgtex;
    for(dims[1] = innerect[1]; dims[1] < innerect[1] + rect[3]-8;dims[1]++){coor[1] = rect[1]+8 + dims[1] - innerect[1];
        for(dims[0] = innerect[0]; dims[0] < innerect[0] + rect[2]-8;dims[0]++){ coor[0] = rect[0]+8 + dims[0] - innerect[0];
        rgba_map(coor) = rgba_map_tmp(dims);
        alias_map(coor) = alias_map_tmp(dims);
    }
    }
                tf.fetch(bgtex);
                for(altcoor[2]=0;altcoor[2]<8;altcoor[2]++){
                for(altcoor[1]=0;altcoor[1]<8;altcoor[1]++){
                    coor[0] = rect[0] + 7 - altcoor[1];
                    coor[1] = rect[1] + 7 - altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    coor[0] = rect[0] + rect[2] - 8 + altcoor[1];
                    coor[1] = rect[1] + 7 - altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    coor[0] = rect[0] + rect[2] - 8 + altcoor[1];
                    coor[1] = rect[1] + rect[3] - 8 + altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    coor[0] = rect[0] + 7 - altcoor[1];
                    coor[1] = rect[1] + rect[3] - 8 + altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    }
                }

                for(coor[1]=8+rect[1];coor[1]<rect[1] +rect[3]-8;coor[1]++){
                for(altcoor[2]=88;altcoor[2]<96;altcoor[2]++){
                    altcoor[1] = coor[1] & 127;
                    coor[0] = rect[0] + 95 - altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    coor[0] = rect[0] + rect[2] - 96 + altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    }
                }
                for(coor[0]=8+rect[0];coor[0]<rect[0] +rect[2]-8;coor[0]++){
                for(altcoor[2]=88;altcoor[2]<96;altcoor[2]++){
                    altcoor[1] = coor[0] & 127;
                    coor[1] = rect[1] + 95 - altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    coor[1] = rect[1] + rect[3] - 96 + altcoor[2];
                    altcoor[0]=3;if (bgtex(altcoor) != 0){
                    for(altcoor[0]=0;altcoor[0]<4;altcoor[0]++)  ((unsigned char*)&dapix)[altcoor[0]] = bgtex(altcoor);
                    rgba_map(coor) = dapix; alias_map(coor) = GUI_alias;}
                    }
                }
   //        }
    }
*/
	void GUIArea::setText(char*){
	}
uint32_t GUIArea::getArrayIDfromMouse()const{
    return 0;
}

GUIArray::GUIArray(unsigned int maxnbitem, unsigned int alias, unsigned int style_alias) : GUIObject(alias,style_alias), cursize(0), maxsize(maxnbitem), nb_columns(1), selected(0){
    // happy valgrind!
    memset(rect, '\0', sizeof(rect));
    memset(offset, '\0', sizeof(offset));
}


void GUIArray::wrSubPosition(int32_t *fout, const GUIObject* subptr)const{
    fout[0] = rect[0] + subptr->getRect()[0];
    fout[1] = rect[1] + subptr->getRect()[1];
}
void GUIArray::setArrayFormat(int _nb_columns, int elem_width, int elem_height){
    nb_columns =_nb_columns;
    offset[0]=elem_width;
    offset[1]=elem_height;
}
void GUIArray::setCursize(int32_t _size){
    cursize = _size;
    rect[4] =0;
    rect[5] =0;
    rect[6] = offset[0] * nb_columns;
    rect[7] = offset[1] * (_size / nb_columns);
    if (rect[7] < rect[3]) rect[7] = rect[3];
}
void GUIArray::insertGUI(Uint32 a, RELPOS_enum relpos_type, unsigned int pad_x, unsigned int pad_y){
    GUIObject* obj = ctrl_state.gui_objects_ptr[a].target;
    subs.push_back(a);
    int i;
    int32_t* crect;
    ctrl_state.gui_objects_ptr[a].parent_alias = this->GUI_alias;
    for(i=0;i< maxsize;i++){
        crect = obj->deref(i).accessRect();
        switch(((int)relpos_type) & 7){
            case 0: // centered
                crect[0] = ((offset[0] - crect[2]) >> 1) + pad_x;
            break;
            case 1: // left
                crect[0] = pad_x;
            break;
            case 2: // right
                crect[0] = (offset[0] - crect[2]-1) - pad_x;
            break;
        }
        switch(((int)relpos_type) & 56){
            case 0: // centered
                crect[1] = ((offset[1] - crect[3]) >> 1) + pad_y;
            break;
            case 8: // left
                crect[1] = pad_y;
            break;
            case 16: // right
                crect[1] = (offset[1] - crect[3]-1)-pad_y;
            break;
        }
        crect[0] += offset[0] * (i % nb_columns);
        crect[1] += offset[1] * (i / nb_columns);
        this->onResize();
    }

    }


void GUIArray::removeGUI(Uint32 a){for(unsigned int i=0;i< subs.size();i++) if (subs[i] == a) subs.pop_swap(i); ctrl_state.gui_objects_ptr[a].parent_alias = a;}

void GUIArray::update(){

}

GUIMSG_Enum GUIArray::processGUIevent(const GUIMSG_Enum event){
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];

    if ((rect[7] != rect[3])&&(ctrl_state.mouse_coor[0] > rect[0]+ rect[2] - gui_style.scroll_width)) { // in scroll area!
        if (event == GUIMSG_MOUSE_CLICK_LBUTTON){
            if (ctrl_state.mouse_coor[1] > rect[1]+ rect[3] - gui_style.scroll_width) rect[5] = (rect[5] + 30 < rect[7] - rect[3]) ? rect[5] + 30 : rect[7] - rect[3] ;
            else if  (ctrl_state.mouse_coor[1] < rect[1]+ gui_style.scroll_width) rect[5] = (rect[5] > 30) ? rect[5] - 30 : 0;
            else rect[5] = ((ctrl_state.mouse_coor[1] - rect[1]- gui_style.scroll_width) * (rect[7] - rect[3])) / (rect[3] - gui_style.scroll_width*2);
        }else if ((event == GUIMSG_MOUSE_DRAG)||(event == GUIMSG_MOUSE_DROP_LBUTTON)||(event == GUIMSG_MOUSE_DRAG_LBUTTON)){
            if (ctrl_state.mouse_coor[1] > rect[1]+ rect[3] - gui_style.scroll_width) rect[5] = rect[7] - rect[3];
            else if  (ctrl_state.mouse_coor[1] < rect[1]+ gui_style.scroll_width) rect[5] = 0;
            else rect[5] = ((ctrl_state.mouse_coor[1] - rect[1]- gui_style.scroll_width) * (rect[7] - rect[3])) / (rect[3] - gui_style.scroll_width*2);
            return GUIMSG_NULL;
        }
	}


	return event;
}


void GUIArray::draw(bool is_over, const int32_t* par_rect){
    uint32_t i,j;
    int32_t subshift[4];
    subshift[0] = 0;
    GUIObject* tmpob;



    //if (this->gui_flags & 0x80000000) return;


	GLuint dashader = ctrl_state.daframe_shader;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];

	glUseProgram(dashader);

	bool is_over2 = ((ctrl_state.mouse_stencil & 0xFE)== 2u);
	if ((!is_over)&&(is_over2)){
        for(i=0;i<subs.size();i++) if (subs[i] == ctrl_state.mouse_stencilID) break;
        is_over2 = (i < subs.size());

	}



	if (is_over) glStencilFunc(GL_ALWAYS, 1, 0xFFFF);


    gui_style.setColorUniforms(dashader,0);
	glUniform1i(glGetUniformLocation(dashader, "border_index"),0);
	glUniform1i(glGetUniformLocation(dashader, "tPalette"),1);
	glUniform1i(glGetUniformLocation(dashader, "tex_color"),2);

	//glUniform1f(glGetUniformLocation(ctrl_state.daframe_shader, "bordersize"),gui_style.border_width);
	glActiveTexture(GL_TEXTURE2);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
	glActiveTexture(GL_TEXTURE1);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
	glActiveTexture(GL_TEXTURE0);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
	glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );


	if (rect[7] == rect[3]) { // no vertical scroll
        for(j=0;j<cursize;j++) {
            subshift[1] = ctrl_state.curwin->rect[3] - rect[1] - offset[1] * (j / nb_columns);
            subshift[0] = rect[0] + offset[0] * (j % nb_columns);
            glUniform4f(glGetUniformLocation(dashader, "boundrect"), subshift[0], subshift[1] - offset[1], subshift[0] + offset[0], subshift[1]);
            if ((is_over2)&&(ctrl_state.curwin->rect[3] - ctrl_state.mouse_coor[1] < subshift[1])&&(ctrl_state.curwin->rect[3] - ctrl_state.mouse_coor[1] >= subshift[1] - offset[1])){
                gui_style.setColorUniforms(dashader,1);
                glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
                gui_style.setColorUniforms(dashader,0);
            }else{
                glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
            }
        }
	}else{
		if (ctrl_state.mouse_coor[0] > rect[0] + rect[2] - gui_style.scroll_width) is_over2 = false;
	//    rect[5] = ((clock() & 4095) * (rect[7] - rect[3])) >> 12;
        subshift[2] = offset[0];
        subshift[3] = offset[1];
        for(j=0;j<cursize;j++) {
            subshift[1] = offset[1] * (j / nb_columns);
            subshift[0] = offset[0] * (j % nb_columns);
            if (gui_style.setFrameUniforms(dashader,subshift,rect)){
                if ((is_over2)&&(ctrl_state.mouse_coor[1] > rect[1] + subshift[1] - rect[5])&&( ctrl_state.mouse_coor[1] <= rect[1] + subshift[1] - rect[5] + offset[1])){
                    gui_style.setColorUniforms(dashader,1);
                    glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
                    gui_style.setColorUniforms(dashader,0);
                }else{
                    glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
                }
            }
        }



	dashader = ctrl_state.daicon_shader;
	glUseProgram(dashader);
//    glUniform1f(glGetUniformLocation(dashader, "bordersize"),1);

	glUniform4f(glGetUniformLocation(dashader, "color"), ((double)((clock() & 0x7FFF) >>7))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color2"), gui_style.colors[2][0], gui_style.colors[2][1], gui_style.colors[2][2], 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color3"), gui_style.colors[3][0], gui_style.colors[3][1], gui_style.colors[3][2], 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color4"), 0.4f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color5"), 0.6f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color6"), 0.7f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color7"), 0.9f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
    glUniform1f(glGetUniformLocation(dashader, "timeloop"), ((double)((clock() & 0xFFF) >>4))/ 256.0f );
	glUniform4f(glGetUniformLocation(dashader, "maskRBG"), 0.0f,0.0f,0.0f,0.0f);
    glActiveTexture(GL_TEXTURE2);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
            ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_SCROLL);
            glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
            glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f, 0.75f, 1.0f, 1.0f);

        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1]- rect[3] +gui_style.scroll_width , rect[0]+rect[2]-1,ctrl_state.curwin->rect[3] - rect[1] - gui_style.scroll_width);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
            glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - rect[3] , rect[0]+rect[2]-1, ctrl_state.curwin->rect[3] - rect[1]- rect[3] +gui_style.scroll_width);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - gui_style.scroll_width , rect[0]+rect[2]-1 , ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();

        j = gui_style.scroll_width + (((rect[3] - (gui_style.scroll_width * 3)) * rect[5]) / (rect[7] - rect[3]));
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - j - gui_style.scroll_width , rect[0]+rect[2]-1 , ctrl_state.curwin->rect[3] - rect[1] - j );

        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();


	}

	if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);

    for(j=0;j<subs.size();j++) {
        tmpob = ctrl_state.gui_objects_ptr[subs[j]].target; if (!tmpob) exit(1);
        if (tmpob == NULL) {fprintf(stderr, "GUIArray has an unexisting sub!\n"); exit(1);}
        if (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == subs[j])){
            for(i=0;i<cursize;i++){
                tmpob->deref(i).draw( ctrl_state.mouse_array_offset == i , rect);
            }
        }else{
            for(i=0;i<cursize;i++){
                tmpob->deref(i).draw(false, rect);
            }
        }
	}
}

void GUIArray::drawAlias(bool is_text, const int32_t* par_rect){

    uint32_t i,j;
    int32_t subshift[4];
    GUIObject* tmpob;

    //if (this->gui_flags & 0x80000000) return;


    GLuint dashader = ctrl_state.daframe_alias_shader;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    glUseProgram(dashader);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else {
        glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * ((GUI_alias>>24) & 255),(1.0f / 255) * ((GUI_alias>>16) & 255), (1.0f / 255) * ((GUI_alias>>8) & 255),(1.0f / 255) * (GUI_alias & 255));
        glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    }
    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);

    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );

    subshift[2] = offset[0];
    subshift[3] = offset[1];

	if (rect[7] == rect[3]) { // no vertical scroll
        for(j=0;j<cursize;j++) {
            subshift[1] = ctrl_state.curwin->rect[3] - rect[1] - offset[1] * (j / nb_columns);
            subshift[0] = rect[0] + offset[0] * (j % nb_columns);
            glUniform4f(glGetUniformLocation(dashader, "boundrect"), subshift[0], subshift[1] - offset[1], subshift[0] + offset[0], subshift[1]);
            glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
        }
	}else{
        for(j=0;j<cursize;j++) {
            subshift[1] = offset[1] * (j / nb_columns);
            subshift[0] = offset[0] * (j % nb_columns);
            if (gui_style.setFrameUniforms(dashader,subshift,rect)){
                glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
            }
        }

        dashader = ctrl_state.daicon_alias_shader;
        glUseProgram(dashader);
        glActiveTexture(GL_TEXTURE1);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
        //glBindTexture(GL_TEXTURE_2D,0);
        glActiveTexture(GL_TEXTURE0);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
        glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
        glUniform4f(glGetUniformLocation(dashader, "id_color"), (1.0f / 255) * ((GUI_alias>>24) & 255),(1.0f / 255) * ((GUI_alias>>16) & 255), (1.0f / 255) * ((GUI_alias>>8) & 255),(1.0f / 255) * (GUI_alias & 255));
        glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
        glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f, 0.75f, 1.0f, 1.0f);

        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1]- rect[3] +gui_style.scroll_width , rect[0]+rect[2]-1,ctrl_state.curwin->rect[3] - rect[1] - gui_style.scroll_width);
        glBegin(GL_QUADS); glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();
        glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);

        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - rect[3] , rect[0]+rect[2]-1, ctrl_state.curwin->rect[3] - rect[1]- rect[3] +gui_style.scroll_width);
        glBegin(GL_QUADS); glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - gui_style.scroll_width , rect[0]+rect[2]-1 , ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS); glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();

        j = gui_style.scroll_width + (((rect[3] - (gui_style.scroll_width * 3)) * rect[5]) / (rect[7] - rect[3]));
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0]+rect[2] -gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - j - gui_style.scroll_width , rect[0]+rect[2]-1 , ctrl_state.curwin->rect[3] - rect[1] - j );
        glBegin(GL_QUADS); glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();
	}

	for(j=0;j<subs.size();j++) {
        tmpob = ctrl_state.gui_objects_ptr[subs[j]].target; if (!tmpob) exit(1);
        if (tmpob == NULL) {fprintf(stderr, "GUIArray has an unexisting sub!\n"); exit(1);}
        for(i=0;i<cursize;i++){
            tmpob->deref(i).drawAlias(is_text,rect);
        }
	}
}

uint32_t GUIArray::getArrayIDfromMouse()const{ // TODO, account for columns
    //double subshift[2];
    if (rect[7] == rect[3]) {
        return (ctrl_state.mouse_coor[1] -rect[1]) / offset[1];
    }else{
        //printf(" %o / %i\n", offset[1]);
        return (ctrl_state.mouse_coor[1] - rect[1] + rect[5]) / offset[1];
    }
    return 0;
}

GUIScroll::GUIScroll(unsigned int alias, unsigned int style_alias): GUIObject(alias,style_alias), flags(0){

}
void	GUIScroll::update(){
}
ERRCODE	GUIScroll::manifest(){ return 0;
}
ERRCODE	GUIScroll::vanish(){ return 0;
}


void GUIScroll::draw(bool mouse_over, const int32_t* par_rect){
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
	bool is_over = (((ctrl_state.mouse_stencil & 0xFE)== 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
	if (is_over){
		glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
		//glUniform4f(glGetUniformLocation(dashader, "color1"), 0.12f,1.0,0.7,1.0f); // stk bar
		//glUniform4f(glGetUniformLocation(dashader, "color3"), 0.36f,0.2,0.7,1.0f); // round
	}else{
		//glUniform4f(glGetUniformLocation(dashader, "color1"), 0.12f,1.0,0.5,1.0f); // stk bar
		//glUniform4f(glGetUniformLocation(dashader, "color3"), 0.36f,0.2,0.7,1.0f); // round
	}

	GLuint dashader = ctrl_state.daicon_shader;
	glUseProgram(dashader);
//    glUniform1f(glGetUniformLocation(dashader, "bordersize"),1);

    //gui_style.setColorUniforms(dashader,1);

	glUniform4f(glGetUniformLocation(dashader, "color"), ((double)((clock() & 0x7FFF) >>7))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color2"), gui_style.colors[2][0], gui_style.colors[2][1], gui_style.colors[2][2], 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color3"), gui_style.colors[3][0], gui_style.colors[3][1], gui_style.colors[3][2], 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color4"), 0.4f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color5"), 0.6f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color6"), 0.7f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color7"), 0.9f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
    glUniform1f(glGetUniformLocation(dashader, "timeloop"), ((double)((clock() & 0xFFF) >>4))/ 256.0f );
	glUniform4f(glGetUniformLocation(dashader, "maskRBG"), 0.0f,0.0f,0.0f,0.0f);
    glActiveTexture(GL_TEXTURE2); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0);


    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);


 //   if (this->button_size != 50) {
//		printf("Got to fix size... %i,%i -> %i\n", this->rect[2], this->rect[3], this->button_size);
//		this->button_size = 50;
//	}



	if (gui_style.scroll_width > 0){
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_SCROLL);
        glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.30f, 0.30f, 0.45f, 0.45f);
        glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.625f, 0.75f, 0.625f, 1.0f);

        // draw scroll background!
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + gui_style.scroll_width, ctrl_state.curwin->rect[3] - rect[1]- rect[3] , rect[0]+rect[2]-1 - gui_style.scroll_width,ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS); glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();

        // draw left button
        if (rect[2] > rect[3]){
            glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + gui_style.scroll_width, ctrl_state.curwin->rect[3] - rect[1]);
		}else{
            glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - gui_style.scroll_width, rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
		}
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
		glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();

        // draw right button
        if (rect[2] > rect[3]){
            glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + rect[2] - gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2]-1, ctrl_state.curwin->rect[3] - rect[1]);
        }else{
            glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1] - rect[3] + gui_style.scroll_width);
        }
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1); glEnd();

        dashader = ctrl_state.daframe_shader;
        glUseProgram(dashader);
        glActiveTexture(GL_TEXTURE2); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
        glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
        glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTON);
        //glUniform1f(glGetUniformLocation(dashader, "bordersize"),gui_style.border_width);
        glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
        if (!gui_style.setFrameUniforms(dashader, rect, par_rect)) return;
        if (rect[2] > rect[3]){
            if (flags & 1){
                int offset = gui_style.scroll_width + (int)(((value - min) / (max - min)) * (rect[2] - gui_style.scroll_width * 2));
                glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offset, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + offset + (int)(((rect[2] - gui_style.scroll_width * 2) * minstep) / (max - min)), ctrl_state.curwin->rect[3] - rect[1]);
            }else{
                int offset = gui_style.scroll_width + (int)(((value - min) / (max - min)) * (rect[2] - gui_style.scroll_width * 3));
                glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offset, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + offset + gui_style.scroll_width, ctrl_state.curwin->rect[3] - rect[1]);
            }
        }else{
            if (flags & 1){
                int offset = gui_style.scroll_width + (int)(((value - min) / (max - min)) * (rect[3] - gui_style.scroll_width * 2));
                glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - offset - (int)(((rect[3] - gui_style.scroll_width * 2) * minstep) / (max - min)), rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1] - offset);
            }else{
                int offset = gui_style.scroll_width + (int)(((value - min) / (max - min)) * (rect[3] - gui_style.scroll_width * 3));
                glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - offset - gui_style.scroll_width, rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1] - offset);
            }
        }
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
	}
    if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}

void	GUIScroll::drawAlias(bool is_text, const int32_t* par_rect){
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
	GLuint dashader = ctrl_state.daicon_alias_shader;
	glUseProgram(dashader);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(dashader, "id_color"), (1.0f / 255) * ((GUI_alias>>24) & 255),(1.0f / 255) * ((GUI_alias>>16) & 255),(1.0f / 255) * ((GUI_alias>>8) & 255), (1.0f / 255) * (GUI_alias & 255));
	if (gui_style.scroll_width > 0){
        glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.30f, 0.30f, 0.45f, 0.45f);
        glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.30f, 0.30f, 0.45f, 0.45f);
        glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + gui_style.scroll_width, ctrl_state.curwin->rect[3] - rect[1]- rect[3] , rect[0]+rect[2]-1 - gui_style.scroll_width,ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
        glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + gui_style.scroll_width, ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + rect[2] - gui_style.scroll_width-1, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2]-1, ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
	}
}
void	GUIScroll::onResize(){
	/*
	if (this->rect[2] > this->rect[3]){
		this->button_size = (this->rect[2] * 4 >= this->rect[3]) ?  this->rect[3] : ((double)this->rect[3] * this->rect[3]) / this->rect[2];
	}else{
		this->button_size = (this->rect[3] * 4 >= this->rect[2]) ? -this->rect[2] : ((double)-this->rect[2] * this->rect[2]) / this->rect[3];
	}*/
}
void	GUIScroll::setScrollingAndRange(double _val, double _minstep, double _min, double _max, bool has_range){
	min = _min;
    max = _max;
    value = _val;
    if (has_range){
        minstep = _minstep;
        flags |= 1;
    }else{
        flags &= 0xFFFFFFFE;
        minstep = _minstep;
    }
}

GUIMSG_Enum	GUIScroll::processGUIevent(const GUIMSG_Enum event){
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    GUImessage tmp_message;
    double effmax;
    GUIMSG_Enum ev = (event != GUIMSG_MOUSE_DRAG) ? event : (flags & 1) && (ctrl_state.isRightButtonDown()) ? GUIMSG_MOUSE_DROP_RBUTTON : GUIMSG_MOUSE_DROP_LBUTTON;
	switch(ev){
	case GUIMSG_MOUSE_DRAG_LBUTTON:
    case GUIMSG_MOUSE_DROP_LBUTTON:
    case GUIMSG_MOUSE_CLICK_LBUTTON:
        effmax = (flags & 1) ? max - minstep : max;
        if (rect[2] > rect[3]){
            if (ctrl_state.mouse_coor[0] < rect[0] + gui_style.scroll_width){
                value -= minstep;
                if (value < min) value = min;
            }else if (ctrl_state.mouse_coor[0] >= rect[0] + rect[2] - gui_style.scroll_width){
                value += minstep;
                if (value > effmax ) value = effmax;
            }else{
                value = min + ((effmax - min) * (ctrl_state.mouse_coor[0] - rect[0]-gui_style.scroll_width)) / (rect[2]- gui_style.scroll_width*2);
            }
        }else{
            if (ctrl_state.mouse_coor[1] < rect[1] + gui_style.scroll_width){
                value -= minstep;
                if (value < min) value = min;
            }else if (ctrl_state.mouse_coor[1] >= rect[1] + rect[3] - gui_style.scroll_width){
                value += minstep;
                if (value > effmax ) value = effmax;
            }else{
                value = min + ((effmax - min) * (ctrl_state.mouse_coor[1] - rect[1]-gui_style.scroll_width)) / (rect[3]- gui_style.scroll_width*2);
            }
        }
    return GUIMSG_VALUE_CHANGE;
    case GUIMSG_MOUSE_DRAG_RBUTTON:
    case GUIMSG_MOUSE_CLICK_RBUTTON:
        if (flags & 1){
            if (rect[2] > rect[3]){
                if (ctrl_state.mouse_coor[0] < rect[0] + gui_style.scroll_width){
                    minstep += value - min;
                    value = min;
                    flags &= 0xFFFFFFFD;
                }else if (ctrl_state.mouse_coor[0] >= rect[0] + rect[2] - gui_style.scroll_width){
                    minstep = max - value;
                    flags |= 2;
                }else{
                    effmax = min + ((max - min) * (ctrl_state.mouse_coor[0] - rect[0]-gui_style.scroll_width)) / (rect[2]- gui_style.scroll_width*2);
                    if (effmax < value + 0.5f * minstep){
                        minstep += value - effmax;
                        value = effmax;
                        flags &= 0xFFFFFFFD;
                    }else{
                        minstep = effmax - value;
                        flags |= 2;
                    }
                }
            }else{
                if (ctrl_state.mouse_coor[1] < rect[1] + gui_style.scroll_width){
                    minstep += value - min;
                    value = min;
                    flags &= 0xFFFFFFFD;
                }else if (ctrl_state.mouse_coor[1] >= rect[1] + rect[3] - gui_style.scroll_width){
                    minstep = max - value;
                    flags |= 2;
                }else{
                    effmax = min + ((max - min) * (ctrl_state.mouse_coor[1] - rect[1]-gui_style.scroll_width)) / (rect[3]- gui_style.scroll_width*2);
                    if (effmax < value + 0.5f * minstep){
                        minstep += value - effmax;
                        value = effmax;
                        flags &= 0xFFFFFFFD;
                    }else{
                        minstep = effmax - value;
                        flags |= 2;
                    }
                }
            }
        }
        return GUIMSG_VALUE_CHANGE;
    case GUIMSG_MOUSE_DROP_RBUTTON:
        if (rect[2] > rect[3]){
            effmax =  min + ((ctrl_state.mouse_coor[0] < rect[0] + gui_style.scroll_width) ? 0.0f : ((max - min) * (ctrl_state.mouse_coor[0] - rect[0]-gui_style.scroll_width)) / (rect[2]- gui_style.scroll_width*2));
        }else{
            effmax =  min + ((ctrl_state.mouse_coor[1] < rect[1] + gui_style.scroll_width) ? 0.0f : ((max - min) * (ctrl_state.mouse_coor[1] - rect[1]-gui_style.scroll_width)) / (rect[3]- gui_style.scroll_width*2));
        }

        if (flags & 2){
            if (effmax < value){
                value = effmax;
                minstep = 0;
            }else{
                minstep = ((effmax  < max) ? effmax : max) - value;
            }
        }else{
            if (value + minstep < effmax){
                minstep = 0;
                value = (effmax > max) ? max : effmax;
            }else{
                minstep += value - effmax;
                value = effmax;
            }
        }
        return GUIMSG_VALUE_CHANGE;
        default: return GUIMSG_NULL;
	}
	return GUIMSG_NULL;
}

GUIValueBar::GUIValueBar(unsigned int alias, unsigned int style_alias): GUIObject(alias,style_alias), cur_str(NULL),state(0){

}
GUIValueBar::~GUIValueBar(){delete[](cur_str);}
void	GUIValueBar::update(){
}
ERRCODE	GUIValueBar::manifest(){ return 0;
}
ERRCODE	GUIValueBar::vanish(){ return 0;
}

GUIValueBar&	GUIValueBar::operator=(const GUIScroll&){exit(1);return *this;}

void GUIValueBar::setText(const char* newtext){
    if (cur_str) delete[](cur_str);
    unsigned int i = strlen(newtext) +1;
    cur_str = new char[i]; memcpy(cur_str, newtext,i);
	state |= 1;
    }
void GUIValueBar::draw(bool mouse_over, const int32_t* par_rect){ return;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
	bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
	if (is_over){
		glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
		//glUniform4f(glGetUniformLocation(dashader, "color1"), 0.12f,1.0,0.7,1.0f); // stk bar
		//glUniform4f(glGetUniformLocation(dashader, "color3"), 0.36f,0.2,0.7,1.0f); // round
	}else{
		//glUniform4f(glGetUniformLocation(dashader, "color1"), 0.12f,1.0,0.5,1.0f); // stk bar
		//glUniform4f(glGetUniformLocation(dashader, "color3"), 0.36f,0.2,0.7,1.0f); // round
	}

	GLuint dashader = ctrl_state.daicon_shader;
	glUseProgram(dashader);
//    glUniform1f(glGetUniformLocation(dashader, "bordersize"),1);

	glUniform4f(glGetUniformLocation(dashader, "color"), ((double)((clock() & 0x7FFF) >>7))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color2"), gui_style.colors[2][0], gui_style.colors[2][1], gui_style.colors[2][2], 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color3"), gui_style.colors[3][0], gui_style.colors[3][1], gui_style.colors[3][2], 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color4"), 0.4f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color5"), 0.6f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color6"), 0.7f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color7"), 0.9f + ((double)((clock() & 0x1FFF) >>5))/ 256.0f, 0.5f, 0.5f, 1.0f);
    glUniform1f(glGetUniformLocation(dashader, "timeloop"), ((double)((clock() & 0xFFF) >>4))/ 256.0f );
	glUniform4f(glGetUniformLocation(dashader, "maskRBG"), 0.0f,0.0f,0.0f,0.0f);
    glActiveTexture(GL_TEXTURE2);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0);

    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);


 //   if (this->button_size != 50) {
//		printf("Got to fix size... %i,%i -> %i\n", this->rect[2], this->rect[3], this->button_size);
//		this->button_size = 50;
//	}


	if (gui_style.border_width > 0){
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + gui_style.border_width, ctrl_state.curwin->rect[3] - rect[1]);
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();

		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.875f, 0.875f, 1.0f, 1.0f);
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_SCROLL);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();

		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + rect[2] - gui_style.border_width, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.875f, 0.875f, 1.0f, 1.0f);
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_SCROLL);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
		int offset = gui_style.border_width + (int)((value / max_value) * (rect[2] - gui_style.border_width * 3));
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offset, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + offset + gui_style.border_width, ctrl_state.curwin->rect[3] - rect[1]);
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.875f, 0.875f, 1.0f, 1.0f);
		ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_SCROLL);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
	}
    if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}
void	GUIValueBar::drawAlias(bool is_text, const int32_t* par_rect){ return;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
	GLuint dashader = ctrl_state.daicon_alias_shader;
	glUseProgram(dashader);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(dashader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
	glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);
	glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f, 0.25f, 0.5f, 0.5f);
	if (gui_style.border_width > 0){
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + gui_style.border_width, ctrl_state.curwin->rect[3] - rect[1]);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + rect[2] - gui_style.border_width, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
		int offset = gui_style.border_width + (int)((value / max_value) * (rect[2] - gui_style.border_width * 3));
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.75f, 0.75f, 1.0f, 1.0f);
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offset, ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + offset + gui_style.border_width, ctrl_state.curwin->rect[3] - rect[1]);
	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();
	}
}
void GUIValueBar::onResize(){}

GUIMSG_Enum	GUIValueBar::processGUIevent(const GUIMSG_Enum event){
	/*switch(event){
	case

	}*/
	return GUIMSG_NULL;
}
GUIButton::GUIButton(unsigned int alias, unsigned int style_alias) : GUIObject(alias,style_alias), state(0){
}
GUIButton::GUIButton(unsigned int alias, unsigned int style_alias, const char* _text) : GUIObject(alias, style_alias), texta(_text),state(1) {
}
GUIButton* GUIButton::mkArray(uint32_t array_size, unsigned int alias, unsigned int style_alias){
    GUIButton* fout = new GUIButton[array_size];
    for(int i=0;i<array_size;i++){
        fout[i].initialize_routine(alias, style_alias);
        fout[i].state = 0;
    }
    ctrl_state.gui_objects_ptr[alias].target = fout;
    ctrl_state.gui_objects_ptr[alias].parent_alias = alias;
    return fout;
}
void GUIButton::onResize(){
}
void GUIButton::update(){
} // if the object is selected, it may receive updates
void GUIButton::draw(bool mouse_over, const int32_t* par_rect){

    GLuint dashader = ctrl_state.daframe_shader;
    glUseProgram(dashader);
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];

	if (mouse_over){
		glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
		gui_style.setColorUniforms(dashader,1);
	}else{
	    gui_style.setColorUniforms(dashader,0);

	}

    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),gui_style.border_width);
    glActiveTexture(GL_TEXTURE2);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
	glActiveTexture(GL_TEXTURE1);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
	glActiveTexture(GL_TEXTURE0);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTON);
	/*
	glActiveTexture(GL_TEXTURE2);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture();*/
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);

    if (!gui_style.setFrameUniforms(dashader, rect, par_rect)) return;


	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();

    Vector<const char* > text_input;
    if (state & 1){
        text_input.push_back(texta.text);


        gui_style.makeTextMesh(texta, text_input);

        state |= 2;
        state &= 0xFFFFFFFE;
    }

    if (state & 2) {
        glDisable(GL_STENCIL_TEST);
        gui_style.drawTextMesh(texta, this->rect, par_rect);
        /*if (par_rect == NULL) ctrl_state.curwin->drawTextMesh(glbuffer_indexes,this->texop, this->rect);
        else{
            int32_t a_rect[4];
            a_rect[0] = rect[0] + par_rect[0];
            a_rect[1] = rect[1] + par_rect[1];
            a_rect[2] = rect[2];
            a_rect[3] = rect[3];
            ctrl_state.curwin->drawTextMesh(glbuffer_indexes,this->texop, a_rect);
        }*/
        glEnable(GL_STENCIL_TEST);
    }
    if (mouse_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}
void GUIButton::drawAlias(bool is_text, const int32_t* par_rect){
    GLuint dashader = ctrl_state.daframe_alias_shader;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    glUseProgram(dashader);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(ctrl_state.daframe_alias_shader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);

    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTON);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
    if (!gui_style.setFrameUniforms(dashader, rect, par_rect)) return;
	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();
}
ERRCODE GUIButton::manifest(){return 0;}
ERRCODE GUIButton::vanish(){return 0;}
void GUIButton::setText(const char* newtext){
    if (texta.text) delete[](texta.text);
    unsigned int i = strlen(newtext) +1;
    texta.text = new char[i]; memcpy(texta.text, newtext,i);
    state |= 1;
}

GUIMSG_Enum GUIButton::processGUIevent(const GUIMSG_Enum event){
	switch(event){
		case GUIMSG_MOUSE_ENTER:
			ctrl_state.ressourceHanddle->useSound(GUISOUND_MOUSEOVER);
		return  GUIMSG_NULL;
		case GUIMSG_MOUSE_MOVE: return GUIMSG_NULL;
		case GUIMSG_MOUSE_DRAG: return GUIMSG_NULL;
		default: return event;
	}
}


DisplayIconData::DisplayIconData(): texture(0), animPeriod(0), text(NULL){
}

DisplayIconData::~DisplayIconData(){
	if (text) delete[](text);
}

DisplayIconData& DisplayIconData::toMemmove(DisplayIconData& other){
	ExOp::toMemmove(text, other.text);
	int i,j;
	texture = other.texture;
	for(i=0;i<3;i++) UVdata[i] = other.UVdata[i];
	for(i=0;i<4;i++) for(j=0;j<8;j++) colors[i][j] = other.colors[i][j];
	animPeriod = other.animPeriod;
	return *this;
}
void DisplayIconData::setTexture(int id){
//    printf("Set texture to %X\n", id);
    if (id == texture) return;
    if (texture != 0) ctrl_state.ressourceHanddle->deallocTextureExt(texture);
    texture = id;
//    printf("query texture %X:%i\n",(intptr_t)&texture, texture);
    if (texture != 0) ctrl_state.ressourceHanddle->allocTextureExt(texture);
//    printf("query texture %X:%i\n",(intptr_t)&texture, texture);
}

int DisplayIconData::getTexture()const{
//    printf("query texture %X:%i\n",(intptr_t)&texture, texture);
    return(texture);
    }


void DisplayIconData::setText(const char* newtext){
    if (text) delete[](text);
    unsigned int i = strlen(newtext) +1;
    if (i == 1) text = NULL;
    else{text = new char[i]; memcpy(text, newtext,i);}
}

void DisplayIconData::setColor(int index, double hue, double bright, double sat){
	colors[index][0] = hue;
	colors[index][1] = bright;
	colors[index][2] = sat;
	colors[index][3] = 1.0f;
	colors[index][4] = hue;
	colors[index][5] = bright;
	colors[index][6] = sat;
	colors[index][7] = 1.0f;

}
void DisplayIconData::setColorPair(int index, double hue, double bright, double sat, double hue2, double bright2, double sat2){
	colors[index][0] = hue;
	colors[index][1] = bright;
	colors[index][2] = sat;
	colors[index][3] = 1.0f;
	colors[index][4] = hue2;
	colors[index][5] = bright2;
	colors[index][6] = sat2;
	colors[index][7] = 1.0f;
}
void DisplayIconData::setUVsquarre(double _size, double x, double y, double border_size){
	animPeriod = 0;
	UVdata[0] = x - border_size;
	UVdata[1] = y - border_size;
	UVdata[2] = _size + border_size * 2;
}
void DisplayIconData::setUVquad(int index){
	animPeriod = 0;
	UVdata[0] = 0.25f*(index % 4);
	UVdata[1] = 0.25f*(index / 4);
	UVdata[2] = 0.25f;
}
void DisplayIconData::setUVanim16(int milliPeriod, int width, int downsize_mag){
    animPeriod = milliPeriod;
	UVdata[0] = 1.0f / 1024.0f; //(width << 3) ;
	UVdata[1] = 63.0f/64.0f;
}


void DisplayIconData::setUVuniform( GLint location)const{
    if (animPeriod){
        int tst = clock() << 4;
        float frac = fmod(((double)tst) / (animPeriod >> 1), 2.0f); // 0-2.0f
        switch(( tst /  animPeriod) & 15 ){
            case 0: glUniform4f(location, UVdata[0], UVdata[0] * (1.0f + frac) + UVdata[1], UVdata[0] + UVdata[1], UVdata[0] * (1.0f + frac)); break;
            case 1: glUniform4f(location, UVdata[0], UVdata[0] * (3.0f + frac) + UVdata[1], UVdata[0] + UVdata[1], UVdata[0] * (3.0f + frac)); break;
            case 2: glUniform4f(location, UVdata[0], UVdata[0] * (5.0f + frac) + UVdata[1], UVdata[0] + UVdata[1], UVdata[0] * (5.0f + frac)); break;
            case 3: glUniform4f(location, UVdata[0] * (1.0f + frac), UVdata[0] * 7.0f + UVdata[1], UVdata[0] * (1.0f + frac) + UVdata[1], UVdata[0]* 7.0f ); break;
            case 4: glUniform4f(location, UVdata[0] * 3.0f, UVdata[0] * (7.0f - frac) + UVdata[1], UVdata[0] * 3.0f + UVdata[1], UVdata[0] * (7.0f - frac)); break;
            case 5: glUniform4f(location, UVdata[0] * (3.0f + frac), UVdata[0] * 5.0f + UVdata[1], UVdata[0] * (3.0f + frac) + UVdata[1], UVdata[0] * 5.0f); break;
            case 6: glUniform4f(location, UVdata[0] * 5.0f, UVdata[0] * (5.0f + frac) + UVdata[1], UVdata[0] * 5.0f + UVdata[1], UVdata[0] * (5.0f + frac)); break;
            case 7: glUniform4f(location, UVdata[0] * (5.0f + frac), UVdata[0] * 7.0f + UVdata[1], UVdata[0] * (5.0f + frac) + UVdata[1], UVdata[0] * 7.0f); break;
            case 8: glUniform4f(location, UVdata[0] * 7.0f, UVdata[0] * (7.0f - frac) + UVdata[1], UVdata[0] * 7.0f + UVdata[1], UVdata[0] * (7.0f - frac)); break;
            case 9: glUniform4f(location, UVdata[0] * 7.0f, UVdata[0] * (5.0f - frac) + UVdata[1], UVdata[0] * 7.0f + UVdata[1], UVdata[0] * (5.0f - frac) ); break;
            case 10: glUniform4f(location, UVdata[0] * 7.0f, UVdata[0] * (3.0f - frac) + UVdata[1], UVdata[0] * 7.0f + UVdata[1], UVdata[0] * (3.0f - frac) ); break;
            case 11: glUniform4f(location, UVdata[0] * (7.0f - frac), UVdata[0] + UVdata[1], UVdata[0] * (7.0f - frac) +UVdata[1], UVdata[0]); break;
            case 12: glUniform4f(location, UVdata[0] * 5.0f, UVdata[0] * (1.0f + frac) + UVdata[1], UVdata[0] * 5.0f  +UVdata[1], UVdata[0]* (1.0f + frac) ); break;
            case 13: glUniform4f(location, UVdata[0] * (5.0f - frac), UVdata[0] * 3.0f + UVdata[1], UVdata[0] * (5.0f - frac) +UVdata[1], UVdata[0] * 3.0f ); break;
            case 14: glUniform4f(location, UVdata[0] * 3.0f, UVdata[0] * (3.0f - frac) + UVdata[1], UVdata[0] * 3.0f +UVdata[1], UVdata[0]* (3.0f - frac) ); break;
            case 15: glUniform4f(location, UVdata[0] * (3.0f - frac), UVdata[0] + UVdata[1], UVdata[0] * (3.0f - frac) +UVdata[1], UVdata[0]); break;

            default:
                glUniform4f(location, UVdata[0], UVdata[0]+ UVdata[1], UVdata[0]+UVdata[1], UVdata[0]);
            break;
        }

    }else{
		glUniform4f(location, UVdata[0], UVdata[1]+ UVdata[2], UVdata[0]+UVdata[2], UVdata[1]);
    }
}

GUIIcon::GUIIcon(unsigned int alias, unsigned int style_alias) : GUIObject(alias,style_alias), state(0){

	styleID = 666;
}
GUIIcon::GUIIcon(unsigned int alias, unsigned int style_alias, const char* _text) : GUIObject(alias, style_alias, _text),state(1) {

}
void GUIIcon::onResize(){
}
void GUIIcon::update(){
} // if the object is selected, it may receive updates
void GUIIcon::draw(bool mouse_over, const int32_t* par_rect){
    GLuint dashader = ctrl_state.daframe_shader;
    glUseProgram(dashader);
    bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    if (is_over){
        glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
        gui_style.setColorUniforms(dashader,1);
    }else{
        gui_style.setColorUniforms(dashader,0);
    }

    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),32);

	glActiveTexture(GL_TEXTURE2);
	ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTON);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();

    Vector<const char* > text_input;
    if (state & 1){
        text_input.push_back(texta.text);

        gui_style.makeTextMesh(texta, text_input);

        state |= 2;
        state &= 0xFFFFFFFE;
    }

    if (state & 2) {
        glDisable(GL_STENCIL_TEST);
        //ctrl_state.curwin->drawTextMesh(glbuffer_indexes,gui_style.font_op, this->rect);
        gui_style.drawTextMesh(texta, this->rect, par_rect);
        glEnable(GL_STENCIL_TEST);
    }
    if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}
void GUIIcon::drawAlias(bool is_text, const int32_t* par_rect){
    GLuint dashader = ctrl_state.daframe_alias_shader;
    glUseProgram(dashader);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(ctrl_state.daframe_alias_shader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);
    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),32.0);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTON);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();
}
ERRCODE GUIIcon::manifest(){return 0;}
ERRCODE GUIIcon::vanish(){return 0;}
void GUIIcon::setText(const char* newtext){
    if (texta.text) delete[](texta.text);
    unsigned int i = strlen(newtext) +1;
    texta.text = new char[i]; memcpy(texta.text, newtext,i);
    state |= 1;
}


GUIGrid::GUIGrid(unsigned int alias, unsigned int style_alias): GUIObject(alias,style_alias), state(0){

}

void GUIGrid::setGridSize(int x, int y){
	grid_size[0] =x;
	grid_size[1] =y;
	grid.setSize(x*y);
	// DisplayIconData
}

int GUIGrid::getGridOffset(int x, int y) const{
	return (((x - rect[0]) * grid_size[0]) / (rect[2])) + grid_size[0] * (((y - rect[1]) * grid_size[1]) / (rect[3]));
}


void GUIGrid::wrSlotPosize(float* fout, int slot)const{
	fout[0] = rect[0] + (0.5f + (slot % grid_size[0])) * rect[2] / grid_size[0];
    fout[1] = rect[1] + (0.5f + (slot / grid_size[0])) * rect[3] / grid_size[1];
    fout[2] = rect[2] / grid_size[0];
}
DisplayIconData& GUIGrid::getIconData(int x){return grid[x];}
void GUIGrid::setText(int offset, const char* newtext){
	grid[offset].setText(newtext);
	state |= 1;
}
void GUIGrid::update(){

}
void GUIGrid::draw(bool mouse_over, const int32_t* par_rect){
    if (this->gui_flags & 0x80000000) return;
	GUIStyle& gui_style = ctrl_state.gui_styles[styleID];
    if (state & 1){
		{ //:
        //	printf("grid uses %i\n", glbuffer_indexes[0] );
            Vector<TextDisplayData> dadata;
            int i,x,y;
            dadata.setSize(grid_size[0]*grid_size[1]);
            for(i=0,y=0;y<grid_size[1];y++)	for(x = 0; x< grid_size[0];x++,i++){
                dadata[i].text = grid[i].text;
                dadata[i].color = -1;
                dadata[i].position[0] = (rect[2] * (1 + (x <<1))) / (grid_size[0] <<1);
                dadata[i].position[1] = (rect[3] * (1 + (y <<1))) / (grid_size[1] <<1);
            }
            gui_style.makeTextMesh(texta, dadata);
		} //:
        state &= 0xFFFFFFFE;
		state |= 2;
    }

    GLuint dashader = ctrl_state.daicon_shader;
	bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
	if (is_over){
		glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
	}
    glUseProgram(dashader);
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
//    glUniform1f(glGetUniformLocation(dashader, "bordersize"),1);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);

	glUniform4f(glGetUniformLocation(dashader, "color"), 0.0f, 0.5f, 0.0f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color2"), 0.0f, 0.5f, 0.0f, 1.0f);
	glUniform4f(glGetUniformLocation(dashader, "color3"), 0.0f, 0.6f, 0.0f, 1.0f);
    glUniform1f(glGetUniformLocation(dashader, "timeloop"), ((double)((clock() & 0xFFF) >>4))/ 256.0f );
	glUniform4f(glGetUniformLocation(dashader, "maskRBG"), 0.0f,0.0f,0.0f,0.0f);
    glActiveTexture(GL_TEXTURE2);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
	int x,y;
	int offsets_grid[2];
	offsets_grid[0] = rect[2] / grid_size[0];
	offsets_grid[1] = rect[3] / grid_size[1];

	for(int i =0 ; i < grid.getSize();i++){
		x = i % grid_size[0]; y = i / grid_size[0];
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offsets_grid[0] * x, ctrl_state.curwin->rect[3] - rect[1] - offsets_grid[1] * (y+1), rect[0] + offsets_grid[0] * (x+1), ctrl_state.curwin->rect[3] - rect[1] - offsets_grid[1] * y);


		if (x == 0) x = (grid_size[0] == 1) ? 3:0;
		else x = (x == grid_size[0]-1) ? 2:1;
		if (y == 0) y = (grid_size[1] == 1) ? 3:2;
		else y = (y == grid_size[1]-1) ? 0:1;
//		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f*x, 0.25f*y, 0.25f*(x+1), 0.25f*(y+1));
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f*x, 0.25f*y, 0.25f*(x+1), 0.25f*(y+1));
        glBegin(GL_QUADS);
            glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
        glEnd();
	}



    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);

    double grid_tick_offset = 0.05f;

	for(uint32_t i =0 ; i < grid.getSize();i++){
        //printf("drawing slot %i in %X\n", i, this);
        if (grid[i].getTexture() == 0) continue;
		x = i % grid_size[0]; y = i / grid_size[0];
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offsets_grid[0] * x, ctrl_state.curwin->rect[3] - rect[1] - offsets_grid[1] * (y+1), rect[0] + offsets_grid[0] * (x+1), ctrl_state.curwin->rect[3] - rect[1] - offsets_grid[1] * y);


		if (x == 0) x = (grid_size[0] == 1) ? 3:0;
		else x = (x == grid_size[0]-1) ? 2:1;
		if (y == 0) y = (grid_size[1] == 1) ? 3:2;
		else y = (y == grid_size[1]-1) ? 0:1;

		//glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);

		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f*x, 0.25f*y, 0.25f*(x+1), 0.25f*(y+1));
        grid[i].setUVuniform(glGetUniformLocation(dashader, "tex_rectUV"));
		ctrl_state.ressourceHanddle->useTextureExt(grid[i].getTexture());
		//else ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);

		glUniform4f(glGetUniformLocation(dashader, "color"), grid[i].colors[0][0], grid[i].colors[0][1], grid[i].colors[0][2], grid[i].colors[0][3]);
		glUniform4f(glGetUniformLocation(dashader, "color2"), grid[i].colors[1][0], grid[i].colors[1][1], grid[i].colors[1][2], grid[i].colors[1][3]);
		glUniform4f(glGetUniformLocation(dashader, "color3"), grid[i].colors[1][4], grid[i].colors[1][5], grid[i].colors[1][6], grid[i].colors[1][7]);
		glUniform4f(glGetUniformLocation(dashader, "color4"), grid[i].colors[2][0], grid[i].colors[2][1], grid[i].colors[2][2], grid[i].colors[2][3]);
		glUniform4f(glGetUniformLocation(dashader, "color5"), grid[i].colors[2][4], grid[i].colors[2][5], grid[i].colors[2][6], grid[i].colors[2][7]);
		glUniform4f(glGetUniformLocation(dashader, "color6"), grid[i].colors[3][0], grid[i].colors[3][1], grid[i].colors[3][2], grid[i].colors[3][3]);
		glUniform4f(glGetUniformLocation(dashader, "color7"), grid[i].colors[3][4], grid[i].colors[3][5], grid[i].colors[3][6], grid[i].colors[3][7]);

	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();
	}

	if (state & 2) gui_style.drawTextMesh(texta, this->rect, par_rect);
	//ctrl_state.curwin->drawTextMesh(glbuffer_indexes,dastyle.font_op, this->rect);

	if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);

}

void GUIGrid::drawAlias(bool is_text, const int32_t* par_rect){
	if (this->gui_flags & 0x80000000) return;
    GLuint dashader = ctrl_state.daicon_alias_shader;
    glUseProgram(dashader);
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(dashader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
//    glUniform1f(glGetUniformLocation(dashader, "bordersize"),1);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
	//glBindTexture(GL_TEXTURE_2D,0);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);


	int x,y;
	int offsets_grid[2];
	offsets_grid[0] = rect[2] / grid_size[0];
	offsets_grid[1] = rect[3] / grid_size[1];

	for(uint32_t i =0 ; i < grid.getSize();i++){
		x = i % grid_size[0]; y = i / grid_size[0];
		glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0] + offsets_grid[0] * x, ctrl_state.curwin->rect[3] - rect[1] - offsets_grid[1] * (y+1), rect[0] + offsets_grid[0] * (x+1), ctrl_state.curwin->rect[3] - rect[1] - offsets_grid[1] * y);


		if (x == 0) x = (grid_size[0] == 1) ? 3:0;
		else x = (x == grid_size[0]-1) ? 2:1;
		if (y == 0) y = (grid_size[1] == 1) ? 3:2;
		else y = (y == grid_size[1]-1) ? 0:1;
		glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.35f, 0.35f, 0.4f, 0.4f);
		//glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.25f*x, 0.25f*y, 0.25f*(x+1), 0.25f*(y+1));
		glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.25f*x, 0.25f*y, 0.25f*(x+1), 0.25f*(y+1));
        glBegin(GL_QUADS);
            glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
        glEnd();
	}

}
GUITextArea::GUITextArea(unsigned int alias, unsigned int style_alias): GUIObject(alias,style_alias), state(0), cur_str(NULL), dico(NULL){}
GUITextArea* GUITextArea::mkArray(uint32_t array_size, unsigned int alias, unsigned int style_alias){
    GUITextArea* fout = new GUITextArea[array_size];
    for(int i=0;i<array_size;i++){
        fout[i].initialize_routine(alias, style_alias);
        fout[i].state = 0;
        fout[i].cur_str = NULL;
        fout[i].dico = NULL;
    }
    ctrl_state.gui_objects_ptr[alias].target = fout;
    ctrl_state.gui_objects_ptr[alias].parent_alias = alias;
return fout;}
void GUITextArea::onResize(){}
void GUITextArea::setDico(Dictionary* whay){
	dico = whay;
	dico_index = 0xFFFFFFFF;
}
void GUITextArea::startTextCapture(){
    ctrl_state.enter_stringcapture_mode(this);
    ctrl_state.setForgroundObject(GUI_alias);
}
void GUITextArea::setText(const char* newtext){
    if (cur_str) delete[](cur_str);
    unsigned int i = strlen(newtext) +1;
    cur_str = new char[i]; memcpy(cur_str, newtext,i);
	state |= 1;
	Vector<uint32_t> options_text;
    if (dico != NULL){
        dico->findIndexes(options_text,ctrl_state.stringCapture,1);
        dico_index =(options_text.size()==0) ?  0xFFFFFFFF : options_text[0];
    }
}
void GUITextArea::draw(bool mouse_over, const int32_t* par_rect){
	if (gui_flags & 0x80000000) return;
	Vector<const char*> options_text;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
	GLuint dashader = ctrl_state.daframe_shader;
    bool is_over = (( (ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
    if (is_over) glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
    int tsize;
    if (gui_style.border_width != 0) {
        glUseProgram(dashader);
        glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
        if (!gui_style.setFrameUniforms(dashader, rect, par_rect, GUI_alias == ctrl_state.foreground_alias)) return;
    //glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "boundrect"), ((float)texop.left-6), (float)(ctrl_state.curwin->rect[3] - texop.top - texop.height-1-6), (float)(texop.left + texop.width+6), (float)(ctrl_state.curwin->rect[3] - texop.top-1+6) );

        gui_style.setColorUniforms(dashader,0);


        glActiveTexture(GL_TEXTURE1);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
        glActiveTexture(GL_TEXTURE0);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);

        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
	}

    if (state & 1){


		if (state & 4){
			if (dico != NULL){
				//printf("finding subs of '%s'\n", ctrl_state.stringCapture);
				options_text.push_back(ctrl_state.stringCapture);
				dico->findSubs(options_text,ctrl_state.stringCapture,5);
				tsize = gui_style.makeTextMesh(texta, options_text);
			}else{
                tsize = gui_style.makeTextMesh(texta, ctrl_state.stringCapture, ctrl_state.stringCapture_cur);
			}
		}else {
		    tsize = gui_style.makeTextMesh(texta, cur_str);
		    //options_text.push_back(cur_str);
		}
	//	printf("text area uses %i\n", glbuffer_indexes[0]);
/*
        if (rect[3] != 12 + 8 * tsize){
            rect[1] -= (12 + 8 * tsize - rect[3]);

        }*/
        rect[3] = 12 + gui_style.font_op.font * tsize;

        state |= 2;
        state &= 0xFFFFFFFE;
    }
    if (state & 2) {
        if ((is_over)&&(gui_flags & GUIFLAG_FLAGENUM_TEXT_MANIFEST)){
            glStencilFunc(GL_ALWAYS, 3, 0xFFFF);
            gui_style.drawTextMesh(texta, this->rect, par_rect);
        }else gui_style.drawTextMesh(texta, this->rect, par_rect);
    }
    if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}
void GUITextArea::drawAlias(bool is_text, const int32_t* par_rect){
    if (gui_flags & 0x80000000) return;
	GLuint dashader = ctrl_state.daframe_alias_shader;
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    if (gui_style.border_width != 0) {
        glUseProgram(dashader);
        if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
        else glUniform4f(glGetUniformLocation(dashader, "id_color"),
                        (1.0f / 255) * ((GUI_alias>>24) & 255),
                        (1.0f / 255) * ((GUI_alias>>16) & 255),
                        (1.0f / 255) * ((GUI_alias>>8) & 255),
                        (1.0f / 255) * (GUI_alias & 255));
        glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
        /*
        glUniform1f(glGetUniformLocation(dashader, "color"), 1.0f); // inside
        glUniform1f(glGetUniformLocation(dashader, "color2"), 1.0f); // stk bar

        glUniform1f(glGetUniformLocation(dashader, "color3"), 1.0f); // round
        glUniform1f(glGetUniformLocation(dashader, "color4"), 1.0); // around
        glActiveTexture(GL_TEXTURE0);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
        glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );

        gui_style.setFrameUniforms(dashader, rect, par_rect, GUI_alias == ctrl_state.foreground_alias);
        glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
        */


    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);

        //glUniform1f(glGetUniformLocation(dashader, "bordersize"),gui_style.border_width);
        glActiveTexture(GL_TEXTURE0);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
        glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );

        if (!gui_style.setFrameUniforms(dashader, rect, par_rect, GUI_alias == ctrl_state.foreground_alias)) return;

 //       glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
        glBegin(GL_QUADS);
            glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
        glEnd();

	}

    if (gui_flags & GUIFLAG_FLAGENUM_TEXT_MANIFEST){
        gui_style.drawTextMeshAlias(texta, this->rect, par_rect);
    }

}

ERRCODE GUITextArea::manifest(){return true;}
ERRCODE GUITextArea::vanish(){return true;}

GUIMSG_Enum GUITextArea::processGUIevent(const GUIMSG_Enum event){ // DOES NOT WORK IF INBEDDED!
    switch(event){
    case GUIMSG_MOUSE_DOWN_LBUTTON:
        gui_flags |= (uint32_t) GUIFLAG_FLAGENUM_TEXT_MANIFEST;
        return GUIMSG_NULL;
    case GUIMSG_MOUSE_DCLICK_LBUTTON:
        gui_flags &= (0xFFFFFFFF ^ (uint32_t) GUIFLAG_FLAGENUM_TEXT_MANIFEST);
        break;
    case GUIMSG_MOUSE_DROP_LBUTTON:
        gui_flags &= (0xFFFFFFFF ^ (uint32_t) GUIFLAG_FLAGENUM_TEXT_MANIFEST);
        break;
    case GUIMSG_MOUSE_CLICK_LBUTTON:
        gui_flags &= (0xFFFFFFFF ^ (uint32_t) GUIFLAG_FLAGENUM_TEXT_MANIFEST);
        if (state & 4){
            if (dico != NULL){
                GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
                if (ctrl_state.mouse_coor[1] < rect[1]+6+ 2*gui_style.font_op.font) return GUIMSG_NULL;
                Vector<uint32_t> options_text;
                dico->findIndexes(options_text,ctrl_state.stringCapture,5);
                uint32_t index = (ctrl_state.mouse_coor[1] - rect[1] - (6+ 2*gui_style.font_op.font))/(2*gui_style.font_op.font);
                if (index >= options_text.getSize()) index =  options_text.getSize()-1u;
                if (options_text.getSize() == 0) return GUIMSG_NULL;
                dico_index = options_text[index];
                int len = strlen(dico->entries[dico_index]);
                memcpy(ctrl_state.stringCapture,dico->entries[dico_index],len+1);
            }
            state |= 1;
            ctrl_state.exit_stringcapture_mode();
            return GUIMSG_NULL;
        }
    default: break;
	}
	return event;
}

double GUITextArea::getValue()const{ return (dico == NULL) ? ((cur_str == NULL) ? 0.0f : atof(cur_str)) : dico_index ;}
int GUITextArea::getIntValue()const{ return (dico == NULL) ? ((cur_str == NULL) ? 0.0f : atoi(cur_str)) : dico_index;}
GUIProgBar::GUIProgBar(unsigned int alias, unsigned int style_alias):  GUIObject(alias,style_alias), task(NULL) {}
GUIProgBar::~GUIProgBar(){}
void GUIProgBar::attachTask(LFHPrimitive::Event<void> *_task){
	char buffer[256];
	_task->wrName(buffer);
	if (texta.text) delete[](texta.text);
    unsigned int i = strlen(buffer) +1;
    texta.text = new char[i]; memcpy(texta.text, buffer,i);
	state |= 1;
	task = _task;
}
void GUIProgBar::setText(const char* _ntext){
    if (texta.text) delete[](texta.text);
    int i = strlen(_ntext) +1;
    texta.text = new char[i]; memcpy(texta.text, _ntext,i);
    state |= 1;
}
ERRCODE GUIProgBar::manifest(){return 0;}
ERRCODE GUIProgBar::vanish(){return 0;}
void GUIProgBar::onResize(){}

void GUIProgBar::draw(bool mouse_over, const int32_t* par_rect){
	Vector<const char*> options_text;
	GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
    if (is_over) glStencilFunc(GL_ALWAYS, 1, 0xFFFF);

	char buffer[256];
	if (task != NULL){
		task->wrName(buffer);
		if (strcmp(texta.text, buffer) != 0) this->setText(buffer);
	}


    if (state & 1){
		options_text.push_back(texta.text);

        gui_style.makeTextMesh(texta, options_text);
        rect[3] = 12 + gui_style.font_op.font * 2;
        state |= 2;
        state &= 0xFFFFFFFE;
    }
    //
    glUseProgram(ctrl_state.daframe_shader);
    glUniform2f(glGetUniformLocation(ctrl_state.daframe_shader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
    glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);

//glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "boundrect"), ((float)texop.left-6), (float)(ctrl_state.curwin->rect[3] - texop.top - texop.height-1-6), (float)(texop.left + texop.width+6), (float)(ctrl_state.curwin->rect[3] - texop.top-1+6) );

    //glUniform1f(glGetUniformLocation(ctrl_state.daframe_shader, "bordersize"),8.0f);

	double phase = (task != NULL) ? task->getProgress() : 0.0f;
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color"), phase,1.0,0.1,1.0f); // inside
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color2"), 0.12f,1.0,0.5,1.0f); // stk bar
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color3"), 0.36f,0.2,0.4,1.0f); // round
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color4"), phase,1.0,0.6,1.0); // around
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color5"), phase,1.0,0.5,1.0f);
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color6"), 0.74f,1.0,0.5,1.0);
	glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "color7"), 0.86f,1.0,0.5,1.0f);


    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);

	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();


    if (state & 2) gui_style.drawTextMesh(texta, this->rect, par_rect);
    //ctrl_state.curwin->drawTextMesh(glbuffer_indexes,gui_style.font_op, this->rect);

    if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}
void GUIProgBar::drawAlias(bool is_text, const int32_t* par_rect){
	GLuint dashader = ctrl_state.daframe_alias_shader;
    glUseProgram(dashader);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(ctrl_state.daframe_alias_shader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);
    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),8.0f);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
}
double GUIProgBar::getValue()const{return task == NULL ? 0.0f : task->getProgress();}



/*
void GUITextArea::draw_static(LFHPrimitive::DataGrid<Uint32, 2>& rgba_map , LFHPrimitive::DataGrid<Uint32, 2>& alias_map) const{
    if (cur_str == NULL) return;
    Uint32 color;
    Tuple<unsigned int,2> position;
    position[0] = rect[0];
    position[1] = rect[1];
    color = 0xFFFFFFFF;
    for(char* ptr = cur_str;(*ptr) !='\0'; ptr++) rgba_map.drawCharFlip((*ptr),position,color); position[0] += 10;
    color = GUI_alias;
    position[0] = rect[0];
    position[1] = rect[1];
    for(char* ptr = cur_str;(*ptr) !='\0'; ptr++) alias_map.drawCharFlip((*ptr),position,color); position[0] += 10;
	}
*/


void GUIAttribList::draw(bool mouse_over, const int32_t* par_rect){
	}
ERRCODE GUIAttribList::manifest(){return false;}
ERRCODE GUIAttribList::vanish(){return false;}
/*
void GUIAttribList::listen(const GUImessage &msg){
        switch(msg.type){
            case LFHGUI_DRAW:
            break;
            case LFHGUI_WINDOW_OPENCLOSE:
            break;
            case LFHGUI_ALIAS_DRAW:break;
            case LFHGUI_SET_DATA:{//:
                // expects an Vector<char*>
                DataGrid<char*, 2 > &list = *(DataGrid<char*, 2 >*)msg.target;
                cur_strs.setSizes(list.dims);
                DataGrid<char*, 2 >::KeyIterator ite = list.getKeyIterator();
                if (ite.first()) do{
                    cur_strs(ite()) = cloneString(list(ite()));
                }while(ite.next());
            } break;//:
            case LFHGUI_QUERY_RECT_POS:{
                ((int*)msg.target)[0] = rect[0];
                ((int*)msg.target)[1] = rect[1];
                ((int*)msg.target)[2] = rect[2];
                ((int*)msg.target)[3] = rect[3];
                ((int*)msg.target)[4] = rect[0];
            } break;
            case LFHGUI_REQUEST_POINTER: *(GUIAttribList**)msg.target = this; return;
            default:printf("UnHanddled Message!\n");
        }
    }*/

void GUIAttribList::enterValue(unsigned int col, unsigned int row, double val) {
    char buffer[256];unsigned int coor[2];
    coor[0] = col; coor[1]= row;
    delete[](cur_strs(coor));
    snprintf(buffer, sizeof(buffer), "%f", val);
    cur_strs(coor) = cloneString(buffer);
    ExOp::show(cur_strs.dims);
}
//void GUIImageArea::draw(bool mouse_over, const int32_t* par_rect){}

unsigned int MouseStare::operator()(){
    GUImessage tmp_message;
    int guiout;
   // MyWindow* win =(MyWindow*)(ctrl_state.gui_objects_ptr[windowalias].target);
    if ((ctrl_state.last_mouse_button_event_time == oldtime)&&(ctrl_state.mouse_stencilID == oldalias)&&((ctrl_state.mouse_stencil& 0xFE) == 2u)){
        if (ctrl_state.button_state & 7){
            //ctrl_state.button_state |= 8;
            guiout = ctrl_state.gui_objects_ptr[oldalias]->processGUIevent(GUIMSG_MOUSE_DOWN_LBUTTON);
            tmp_message.type = LFHGUI_MOUSE_PUSH;
        }else{
            guiout = ctrl_state.gui_objects_ptr[oldalias]->processGUIevent(LFHGUIEVENT_MOUSE_STARE);
            tmp_message.type = LFHGUI_MOUSE_STARE;
        }
        tmp_message.msg_key = oldalias;
        tmp_message.coor[0] = x;
        tmp_message.coor[1] = y;
        tmp_message.parameter2 = ctrl_state.button_state;
        if (guiout != 0) tmp_message();
    }
    return(0);
}
unsigned int GUITaskEvent::operator()(){
    switch(task){
        case GUITASK_DRAW_STATIC:
            glMatrixMode(GL_PROJECTION);glLoadIdentity();
            glMatrixMode(GL_MODELVIEW);glLoadIdentity();
            glTranslatef(-1.0f,-1.0f,0.0f);
            glScalef(2.0f / ctrl_state.curwin->rect[2],2.0f / ctrl_state.curwin->rect[3],1.0f);
            glBindTexture(GL_TEXTURE_2D,0);
          //  glRasterPos2i(0,0);
          //  ctrl_state.gui_objects_ptr[guiID]->loadstatic();
        break;
        }
    return(0);
    }



GUILog::GUILog(unsigned int given_alias, unsigned int style_alias):GUIArea(given_alias, style_alias),main_index(0),state(0){
    lasttime =SDL_GetTicks() - 5000;
    //ctrl_state.latentqueue.insert_async(0, new GUITaskEvent(GUI_alias, GUITASK_DRAW_STATIC));
    indexes[0][0] = 0xFFFE;
    indexes[0][1] = 0xFFFE;
    for(unsigned int i=1;i< 256; i++) {indexes[i][0] = 0xFFFF;indexes[i][1] = 0xFFFF;}
    log[0xFFFE].k.first = 1;
    log[0xFFFF].k.first = 0;
    }
void GUILog::onResize(){

}
Uint32 GUILog::colorof(unsigned int type)const{
    return 0xFFFFFFFF;
}
void GUILog::startTextCapture(){
    state |= 1;
    ctrl_state.enter_stringcapture_mode();
}
void GUILog::insertText_withlength(unsigned char type, const char *str, unsigned int len){
    lasttime =SDL_GetTicks();
    if (state & 0x100) {
        delete[](log[main_index].d);

        // delete those strings!
    } else if (main_index == 65535) state |= 0x100;
    log[main_index].d = new char[len+4];
    memcpy(log[main_index].d+3, str, len);
    log[main_index].d[0] = 0xFF;
    log[main_index].d[len+3] = '\0';
    if (type == 1) {
        *(unsigned short*)(log[main_index].d+1) = 0xF888 | (rand() & 0x777);
    } else *(unsigned short*)(log[main_index].d+1) = 0xFFFF;
    if (indexes[type][0] == indexes[type][1]){
        log[main_index].k.second =0;
        indexes[type][0] = main_index;
        if (log[indexes[type][1]].k.first != type){ // has 0 items
            indexes[type][1] = main_index;
        }
    }else{ // at least 2 items
        log[main_index].k.second = indexes[type][0] ^ indexes[type][1];
        log[indexes[type][0]].k.second ^= main_index ^ indexes[type][1];
        log[indexes[type][1]].k.second ^= main_index ^ indexes[type][0];
        indexes[type][0] = main_index;
    }
    log[main_index++].k.first = type;
    state |= 0x80;
}
void GUILog::filterIncr(){
    if (channels.getSize() == 0) channels.push_back(0);
    else if (channels[0] == 1) channels.pop_back();
    else channels[0]++;
    state |= 0x80;
}
ERRCODE GUILog::manifest(){
    return 0;
}
ERRCODE GUILog::vanish(){
    return 0;
}
void GUILog::draw(bool mouse_over, const int32_t* par_rect){
    Vector<const char*> todraw;
    unsigned int i;
    unsigned short tcur[2];
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];


    if (state & 0x80) {
        state ^= 0x80;
        switch(channels.getSize()) {
            case 0: // no filter
                if (state & 0x100){
                }else for(i=main_index-1;i+1 != 0; i--) todraw.push_back(log[i].d);
            break;
            case 1:
                if (indexes[channels[0]][0] == indexes[channels[0]][1]){
                    // 1 or 0 entries;
                    if (log[indexes[channels[0]][0]].k.first == channels[0]) todraw.push_back(log[indexes[channels[0]][0]].d);
                }else{
                    // >=2 entries
                    tcur[0] = indexes[channels[0]][0];
                    tcur[1] = log[tcur[0]].k.second ^ indexes[channels[0]][1];

                    todraw.push_back(log[indexes[channels[0]][0]].d);
                    while(tcur[1] != indexes[channels[0]][1]){
                        todraw.push_back(log[tcur[1]].d);
                        tcur[0] ^= log[tcur[1]].k.second;
                        if (tcur[0] == indexes[channels[0]][1]) break;
                        todraw.push_back(log[tcur[0]].d);
                        tcur[1] ^= log[tcur[0]].k.second;
                    }
                    todraw.push_back(log[indexes[channels[0]][1]].d);
                }
            break;
        }
        gui_style.makeTextMesh(texta, todraw);
        }
  //  glDisable(GL_CULL_FACE);
    glUseProgram(ctrl_state.daframe_shader);
    bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
    if (is_over){
        glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
        gui_style.setColorUniforms(ctrl_state.daframe_shader,1);
    }else{
        gui_style.setColorUniforms(ctrl_state.daframe_shader,0);
    }


    //glUniform1f(glGetUniformLocation(ctrl_state.daframe_shader, "bordersize"),16.0f);




    glUniform2f(glGetUniformLocation(ctrl_state.daframe_shader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
//    glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "boundrect"), ((float)texop.left-6), (float)(ctrl_state.curwin->rect[3] - texop.top - texop.height-1-6), (float)(texop.left + texop.width+6), (float)(ctrl_state.curwin->rect[3] - texop.top-1+6) );
    if (!gui_style.setFrameUniforms(ctrl_state.daframe_shader, rect, par_rect)) return;


    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();

	gui_style.drawTextMesh(texta, this->rect, par_rect);
    if (SDL_GetTicks() - lasttime > 5000) return;
}

void GUILog::drawAlias(bool is_text, const int32_t* par_rect){
GLuint dashader = ctrl_state.daframe_alias_shader;
GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    glUseProgram(dashader);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(ctrl_state.daframe_alias_shader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);
    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),16.0f);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
   // glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
	if (!gui_style.setFrameUniforms(dashader, rect, par_rect)) return;
	glBegin(GL_QUADS);
		glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);
	glEnd();

}

GUIDropList::GUIDropList(unsigned int given_alias,unsigned int style_alias) : GUIObject(given_alias, style_alias),state(0){}
GUIDropList::~GUIDropList(){for(uint32_t i=0;i<options.getSize();i++) delete[](options[i]);}

void GUIDropList::leachOptions(Vector<char*> &_options){
    if (options.getSize() != 0) {printf("old had %i at %p\n", options.getSize(), (void*) options.darray);}
    for(uint32_t i=0;i<options.getSize();i++) {printf("old addr %p\n", (void*) options[i]);delete[](options[i]);}
    options.toMemmove(_options); state |= 1; current_selection=0;
    printf("%p and now old had %i at %p %i\n",(void*) this, options.getSize(), (void*) options.darray, current_selection);
}


void GUIDropList::update(){}
void GUIDropList::onResize(){}
void GUIDropList::draw(bool mouse_over, const int32_t* par_rect){
    Vector<const char*> options_text;
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
 //   unsigned int i;
	int optite;
    bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));
    if (is_over) glStencilFunc(GL_ALWAYS, 1, 0xFFFF);
    if (state & 1){
        if (state & 4){ // expended
            printf("%p doom: %i, %p and %i\n", (void*) this,options.getSize(),options.darray , current_selection ); fflush(stdout);
            options_text.push_back(options[current_selection]);
            for(optite=0;optite< (int)options.getSize();optite++) if (optite != current_selection) options_text.push_back(options[optite]);
        }else if (current_selection < (int)options.getSize()){
            options_text.push_back(options[current_selection]);
        }
        gui_style.makeTextMesh(texta, options_text);
        rect[3] = 12 + 2 * gui_style.font_op.font * options_text.getSize();
        state |= 2;
        state &= 0xFFFFFFFE;
    }


    glUseProgram(ctrl_state.daframe_shader);
    glUniform2f(glGetUniformLocation(ctrl_state.daframe_shader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );


    if (!gui_style.setFrameUniforms(ctrl_state.daframe_shader, rect, par_rect, GUI_alias == ctrl_state.foreground_alias)) return;
  //  glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);

//glUniform4f(glGetUniformLocation(ctrl_state.daframe_shader, "boundrect"), ((float)texop.left-6), (float)(ctrl_state.curwin->rect[3] - texop.top - texop.height-1-6), (float)(texop.left + texop.width+6), (float)(ctrl_state.curwin->rect[3] - texop.top-1+6) );


    //glUniform1f(glGetUniformLocation(ctrl_state.daframe_shader, "bordersize"),8.0f);
    gui_style.setColorUniforms(ctrl_state.daframe_shader,0);

    glActiveTexture(GL_TEXTURE2);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALTEX);
    glActiveTexture(GL_TEXTURE1);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);

    glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();


    if (state & 2){
        gui_style.drawTextMesh(texta, this->rect, par_rect);
        //ctrl_state.curwin->drawTextMesh(glbuffer_indexes,this->texop, this->rect);
    }
    if (is_over) glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}
void GUIDropList::drawAlias(bool is_text, const int32_t* par_rect){
    GLuint dashader = ctrl_state.daframe_alias_shader;
    glUseProgram(dashader);
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(ctrl_state.daframe_alias_shader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255),
                    (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    if (!gui_style.setFrameUniforms(dashader, rect, par_rect, GUI_alias == ctrl_state.foreground_alias)) return;
    glUniform2f(glGetUniformLocation(dashader, "color"), 0.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color2"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color3"), 1.0,1.0);
    glUniform2f(glGetUniformLocation(dashader, "color4"), 1.0,1.0);
    //glUniform1f(glGetUniformLocation(dashader, "bordersize"),8.0f);
    glActiveTexture(GL_TEXTURE0);
    ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_FRAME);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );

	glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();

}
GUIMSG_Enum GUIDropList::processGUIevent(const GUIMSG_Enum event){
    switch(event){
        case GUIMSG_MOUSE_DOWN_LBUTTON:
            if (options.getSize() != 0) this->expand();
        return GUIMSG_NULL;
        case GUIMSG_MOUSE_DROP_LBUTTON:{
            if (options.getSize() != 0) {
                if ((ctrl_state.mouse_coor[0] < (uint32_t)rect[0])||(ctrl_state.mouse_coor[1] < (uint32_t)rect[1])) {this->select(current_selection); return GUIMSG_NULL;}
                GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
                // rect[3] = 12 + 16 * options_text.getSize();
                int nsel = (ctrl_state.mouse_coor[1]- rect[1]-6) / (gui_style.font_op.font * 2);
               // printf("dansel %i \n", nsel);
                if (nsel <= 0) {
                    this->select(current_selection);
                }else {
                    if (nsel >= (int)options.getSize()) nsel = (int)(options.getSize())-1;
                    this->select((nsel <= (int)current_selection) ? nsel-1 : nsel);
                }
            }
        }return GUIMSG_NULL;
        default: return GUIMSG_NULL;
    }
}
ERRCODE GUIDropList::manifest(){
//    state |= 2;
    return 0;
}
ERRCODE GUIDropList::vanish(){
    state &= 0xFFFFFFFD;
    return 0;
}
void GUIDropList::expand(){
    state |= 5;
	ctrl_state.setForgroundObject(GUI_alias);
}
void GUIDropList::select(int which){
    state |= 1;
    state &= 0xFFFFFFFB;
    GUImessage tmpmsg;
    tmpmsg.msg_key = this->GUI_alias;
    bool didchange =(current_selection != which);
    tmpmsg.type = GUIMSG_TEXT_CHANGE;
    current_selection = which;
    if (didchange) tmpmsg();
    ctrl_state.clearForgroundObject();
}
GUIHeatMap::GUIHeatMap(unsigned int given_alias,unsigned int style_alias) : GUIObject(given_alias, style_alias), state(0){glbuffer_index = 0;
value_range[0] = 0.0f;
value_range[1] = 0.0f;
selection[2] = 0;
}
void GUIHeatMap::setMapDims(int _row, int _col){
    Tuple<uint32_t,2> coor; coor[0] = _col; coor[1] = _row;
    grid.setSizes(coor);
    state |= 2;
    scroll_rect[2] = _col;
    scroll_rect[3] = _row;
    scroll_rect[0] = 0.0f;
    scroll_rect[1] = 0.0f;
    grid.toRand();
}
void GUIHeatMap::setDataRow(int slot, const float* data){
    state |= 2;
    //state_lastrow_modified = slot;
    uint32_t i;
    float* ptr = grid.data + grid.dims[0] * slot;
    for(i = 0 ; i< grid.dims[0];i++) ptr[i] = data[i];
}
GUIHeatMap& GUIHeatMap::toMemmove(DataGrid<float, 2u> &newdata){
    grid.toMemmove(newdata);
    state |= 2;
    scroll_rect[2] = grid.dims[0];
    scroll_rect[3] = grid.dims[1];
    scroll_rect[0] = 0.0f;
    scroll_rect[1] = 0.0f;
    //state_lastrow_modified = 0xFFFFFFFF;
return *this;}
void GUIHeatMap::draw(bool mouse_over, const int32_t* par_rect){
    if (state & 2){
        if (glbuffer_index == 0) glGenBuffers(1,&glbuffer_index);
        glBindBuffer(GL_ARRAY_BUFFER, glbuffer_index);
     //   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glbuffer_indexes[1]);

     //   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * 6 * grid.dims[0] * grid.dims[1], NULL, GL_STATIC_DRAW);

        float* heatdata = new float[16 * grid.dims[0] * grid.dims[1]];
     //   uint32_t* heatindex = new uint32_t[6 * grid.dims[0] * grid.dims[1]];
        for(uint32_t hite=0;hite<grid.dims[0] * grid.dims[1];hite++){
            heatdata[(hite<<4)] = ((float)(hite % grid.dims[0]));
            heatdata[(hite<<4)|1] = ((float)(hite / grid.dims[0]));
            hite <<= 4;
            heatdata[hite|4] = heatdata[hite];
            heatdata[hite|13] = heatdata[hite|1];
            heatdata[hite|12] = heatdata[hite|8] = heatdata[hite] + 1.0f;
            heatdata[hite|9] = heatdata[hite|5] = heatdata[hite|1] + 1.0f;
            heatdata[hite|14] = heatdata[hite|10] = heatdata[hite|6] =heatdata[hite|2] = fabs(grid.data[hite>>4]);
            heatdata[hite|3] = 0.125f;
            heatdata[hite|7] = 0.375f;
            heatdata[hite|11] = 0.625f;
            heatdata[hite|15] = 0.875f;
            hite >>= 4;
        }
     //   delete[](heatindex);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 16 * grid.dims[0] * grid.dims[1], heatdata, GL_STATIC_DRAW);

        delete[](heatdata);
        state &= 0xFFFFFFFD;
    }/*else if (state & 4){


    }*/

    if (glbuffer_index != 0){
        glUseProgram(ctrl_state.daheat_shader);

        glUniform4f(glGetUniformLocation(ctrl_state.daheat_shader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1] - rect[3], rect[0] + rect[2], ctrl_state.curwin->rect[3] - rect[1]);
        glUniform2f(glGetUniformLocation(ctrl_state.daheat_shader, "windowsize"), ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3] );
        glUniform4f(glGetUniformLocation(ctrl_state.daheat_shader, "projection"), scroll_rect[0], scroll_rect[1], scroll_rect[2], scroll_rect[3]);

        if (selection[2] == 0){
            glUniform4f(glGetUniformLocation(ctrl_state.daheat_shader, "selection"),0,0,0,0);
        }else{
            glUniform4f(glGetUniformLocation(ctrl_state.daheat_shader, "selection"),(float)selection[0],(float)selection[1],(float)(selection[0]+selection[2]),(float)(selection[1]+selection[3]));
            int clo = SDL_GetTicks();
            glUniform1f(glGetUniformLocation(ctrl_state.daheat_shader, "timephase"), 0.0009765625f * ((clo & 512) ? clo & 511 : 512 -(clo & 511) ));
        }

        glUniform2f(glGetUniformLocation(ctrl_state.daheat_shader, "val_range"), value_range[0], value_range[1]);
        glActiveTexture(GL_TEXTURE0);
        ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);

        glBindBuffer(GL_ARRAY_BUFFER, glbuffer_index);
        glVertexAttribPointer(ATTRIBUTE_POSITION,4, GL_FLOAT, GL_FALSE, sizeof(float) *4 , BUFFER_OFFSET(0));
        glDrawArrays(GL_QUADS,0, (grid.dims[0] * grid.dims[1])<<2);
    }
}
void GUIHeatMap::drawAlias(bool is_text, const int32_t* par_rect){
    GUIStyle& gui_style =  ctrl_state.gui_styles[styleID];
	GLuint dashader = ctrl_state.daicon_alias_shader;
	glUseProgram(dashader);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID_MASK);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_BUTTONGRID);
    glUniform2f(glGetUniformLocation(dashader, "windowsize"),  ctrl_state.curwin->rect[2], ctrl_state.curwin->rect[3]);
    if (ctrl_state.alias_storm) glUniform4f(glGetUniformLocation(dashader, "id_color"),(1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), (1.0f / 255) * (rand() & 0xFF), 1.0);
    else glUniform4f(glGetUniformLocation(dashader, "id_color"), (1.0f / 255) * ((GUI_alias>>24) & 255),(1.0f / 255) * ((GUI_alias>>16) & 255),(1.0f / 255) * ((GUI_alias>>8) & 255), (1.0f / 255) * (GUI_alias & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    glUniform4f(glGetUniformLocation(dashader, "mask_rectUV"), 0.30f, 0.30f, 0.45f, 0.45f);
    glUniform4f(glGetUniformLocation(dashader, "tex_rectUV"), 0.30f, 0.30f, 0.45f, 0.45f);
    glUniform4f(glGetUniformLocation(dashader, "boundrect"), rect[0], ctrl_state.curwin->rect[3] - rect[1]- rect[3] , rect[0]+rect[2]-1,ctrl_state.curwin->rect[3] - rect[1]);
    glBegin(GL_QUADS);glVertex2i(0,0);glVertex2i(1,0);glVertex2i(1,1);glVertex2i(0,1);glEnd();
}
Tuple<unsigned int, 2u> GUIHeatMap::getMousePos() const{ Tuple<unsigned int, 2u> fout;
    if ((int32_t)ctrl_state.mouse_coor[0] <= rect[0]) fout[0] = (int)(scroll_rect[0]);
    else if ((int32_t)ctrl_state.mouse_coor[0] >= rect[0] + rect[2]-1) fout[0] = (int)(scroll_rect[0] + scroll_rect[2]);
    else fout[0] = (int)(scroll_rect[0] + (scroll_rect[2] *(ctrl_state.mouse_coor[0] - rect[0]))/ rect[2]);
    if ((int32_t)ctrl_state.mouse_coor[1] <= rect[1]) fout[1] = (int)(scroll_rect[1]);
    else if ((int32_t)ctrl_state.mouse_coor[1] >= rect[1] + rect[3]-1) fout[1] = (int)(scroll_rect[1] + scroll_rect[3]);
    else fout[1] = (int)(scroll_rect[1] + (scroll_rect[3] *(ctrl_state.mouse_coor[1] - rect[1]))/ rect[3]);
    return fout;
}

GUIMSG_Enum GUIHeatMap::processGUIevent(const GUIMSG_Enum event){
    Tuple<unsigned int, 2u> coor;
    switch(event){
        case GUIMSG_MOUSE_CLICK_LBUTTON:
            coor = this->getMousePos();
            selection[0] = coor[0];
            selection[1] = coor[1];
            selection[2] = 1;
            selection[3] = 1;
            selection_value = grid(coor);
        return GUIMSG_VALUE_CHANGE;
        default: return event;
    }
}

GUITextCloud::TextContext::TextContext(): text(NULL){}
GUITextCloud::TextContext::TextContext(const char* _text){uint32_t l =strlen(_text)+1; text = new char[l]; memcpy(text,_text, l);}
GUITextCloud::TextContext::~TextContext(){this->toMemfree();}
GUITextCloud::TextContext& GUITextCloud::TextContext::toZero(){if (text) delete[](text); return *this;}
GUITextCloud::TextContext& GUITextCloud::TextContext::toMemfree(){if (text) delete[](text); return *this;}
GUITextCloud::TextContext& GUITextCloud::TextContext::toMemmove(GUITextCloud::TextContext& other){
    this->toMemfree();
    text = other.text; other.text=NULL;tileID = other.tileID;
return *this;}
void GUITextCloud::TextContext::show(FILE* f, int lvl)const{fprintf(f,"%X: \"%s\"%c", tileID, text, (lvl == 0) ? '\n' : '\t');}
GUITextCloud::TileData& GUITextCloud::TileData::toMemmove(GUITextCloud::TileData& other){
    ExOp::toMemmove(index, other.index);
    ExOp::toMemmove(offsets, other.offsets);
    ExOp::toMemmove(dirty, other.dirty);
    ExOp::toMemmove(texta, other.texta);
return *this;}
void GUITextCloud::TileData::show(FILE* f, int lvl)const{ExOp::show(index, f, lvl+1);ExOp::show(offsets, f, lvl+1); fprintf(f,"d %i\n", dirty);}

GUITextCloud::GUITextCloud(unsigned int alias, unsigned int style_alias) : GUIObject(alias, style_alias){}
GUITextCloud::~GUITextCloud(){}
void GUITextCloud::insert(uint32_t alias, char* text){
    GUITextCloud::TextContext newcontext(text);
    textmap.addEntry(alias).toMemmove(newcontext);
    uint32_t i,j;
    uint32_t newalias;
    for(i=0;i<texta.getSize();i++){
        for(j=0;j<16;j++) if (texta.deref(i).index[j] == 0) break;
        if (j < 16){
            texta.deref(i).index[j] = alias;
            texta.deref(i).dirty |= (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT;
            textmap[alias].tileID = texta.deref_key(i);
            break;
        }
    }
    if (i == texta.getSize()){
        newalias = ctrl_state.mkGuiAlias();
        // new tile needed...
        texta[newalias].index[0] = alias;
        for(j=1;j<16;j++) texta[newalias].index[j] = 0;
        texta[newalias].dirty = (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT;
        textmap[alias].tileID = newalias;
    }
}

void GUITextCloud::setText(Vector<string> datext){
    /*gui_flags |= (uint32_t) GUIFLAG_FLAGENUM_DIRTYTEXT;
    uint32_t i,j;
    if (texts.getSize() != 0){
        for(i=0; i < texts.getSize();i++) delete[](texts[i]);
    }
    texts.setSize(datext.getSize());
    for(i=0; i < texts.getSize();i++) {
        j = datext[i].length() +1;
        texts[i] = new char[j];
        memcpy(texts[i], datext[i].c_str(),j);
    }

    // fake coordinates:
    nbdim = 2;
    project.setSize(nbdim);
    project[0][0] = 1.0f;
    project[0][1] = 0.0f;
    project[1][0] = 0.0f;
    project[1][1] = 1.0f;

    pos.setSize(nbdim * texts.getSize());

    for(i=0; i < nbdim * texts.getSize();i++) pos[i] = ((double) rand())/ RAND_MAX;*/
}


void GUITextCloud::draw(bool mouse_over, const int32_t* par_rect){
    GUIStyle& istyle = ctrl_state.gui_styles[styleID];
    //bool is_over = (((ctrl_state.mouse_stencil & 0xFE) == 2u)&&(ctrl_state.mouse_stencilID == GUI_alias));

    glStencilFunc(GL_ALWAYS, 3, 0xFFFF);
    uint32_t i,j,k;
    Tuple<double, 2u> coor;

    GLint dashader = ctrl_state.datext2_shader;
    glUseProgram(dashader);
    glUniform1i(glGetUniformLocation(dashader, "tex_height"), 20);
    glUniform4f(glGetUniformLocation(dashader, "def_color"), 1.0f,0.0f,0.0f,1.0f);
    glUniform4f(glGetUniformLocation(dashader, "bgcolor_base"), 0.0f,0.0f,0.0f,0.5f);
    glUniform1f(glGetUniformLocation(dashader, "bgcolor_factor"), 0.5f);

    istyle.setFrameUniforms(dashader,rect,par_rect);

    glActiveTexture(GL_TEXTURE2); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXTOFFSETS);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXT);
    Tuple<float, 3u> pos;
  // glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "boundrect"), 0.0, 0.0, (float)(ctrl_state.curwin->rect[2]-1), (float)(ctrl_state.curwin->rect[3]-1) );

    for(i=0;i< texta.getSize();i++){
        if (texta.deref(i).dirty & (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT){
            Vector<const char*> dainput;
            Vector<uint32_t> offsets;
            for(j=0;j<16;j++) if (texta.deref(i).index[j] != 0) {
                dainput.push_back(textmap[texta.deref(i).index[j]].text);
            }
            istyle.makeTextMesh(texta.deref(i).texta,dainput,&offsets);
          //  printf("glinded %i\n", texta.deref(i).d.glbuffer_indexes[0]);

            // populate uint32_t offset[2];
            for(j=0,k=0;j<16;j++){
                texta.deref(i).offsets[j] = (offsets[k] << 2);
                if (texta.deref(i).index[j] != 0) k++;
            }
            texta.deref(i).offsets[16] = (offsets[k] << 2);
            texta.deref(i).dirty -= (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT;
        }

        glBindBuffer(GL_ARRAY_BUFFER, texta.deref(i).texta.glbuffer_indexes[0]);
        glEnableVertexAttribArray(ATTRIBUTE_POSITION);
        glEnableVertexAttribArray(ATTRIBURE_CHARID);
        glVertexAttribPointer(ATTRIBUTE_POSITION,3, GL_FLOAT, GL_FALSE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(0));
        glVertexAttribPointer(ATTRIBURE_CHARID,4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(sizeof(float)*3));
        glUniform2f(glGetUniformLocation(dashader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3]);

        for(j=0;j<16;j++){
            if (texta.deref(i).index[j] != 0) {
                //printf("got %i\n", texta.deref(i).index[j]);
                pos = position_holder(texta.deref(i).index[j]);
              //  glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], obrect[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , obrect[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);
              //  glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], pos[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , pos[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);

                glUniform4f(glGetUniformLocation(dashader, "uPosition"), pos[0], pos[1], pos[2] , 1.0);

                //glUniform4f(glGetUniformLocation(dashader, "uPosition"),  pos[0], pos[1], pos[2], 1.0 );
                //glUniform4f(glGetUniformLocation(dashader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], 0 , 0);
                glDrawArrays(GL_QUADS, texta.deref(i).offsets[j], texta.deref(i).offsets[1 + j] - texta.deref(i).offsets[j]);
            }

        }

    }
/*
LFHTEMP void GUITextCloud<F>::draw(bool mouse_over, const int32_t* par_rect){
    GUIStyle& istyle = ctrl_state.gui_styles[styleID];

    uint32_t i,j,k;
    Tuple<double, 2u> coor;

    GLint dashader = ctrl_state.datext_shader;
    glUseProgram(dashader);
    glUniform1i(glGetUniformLocation(dashader, "tex_height"), 20); // 16);
    glUniform4f(glGetUniformLocation(dashader, "def_color"), 1.0f,0.0f,0.0f,1.0f); // 16);
    glUniform4f(glGetUniformLocation(dashader, "bgcolor_base"), 0.0f,0.0f,0.0f,0.5f); // 16);
    glUniform1f(glGetUniformLocation(dashader, "bgcolor_factor"), 0.5f); // 16);

    glActiveTexture(GL_TEXTURE2); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXTOFFSETS);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXT);


   // glEnableVertexAttribArray(ATTRIBUTE_POSITION);
   // glEnableVertexAttribArray(ATTRIBURE_CHARID);
   // glVertexAttribPointer(ATTRIBUTE_POSITION,3, GL_FLOAT, GL_FALSE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(0));
   // glVertexAttribPointer(ATTRIBURE_CHARID,4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(sizeof(float)*3));
    uint32_t obrect[4];
    obrect[0] =200;
    obrect[1] =200;
    obrect[2] =400;
    obrect[3] =400;
    Tuple<float, 3u> pos;
glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "boundrect"), ((float)pos[0] + istyle.text_borders[0]), (float)(ctrl_state.curwin->rect[3]-1  - pos[1] - obrect[3] + istyle.text_borders[3]), (float)(pos[0] + obrect[2] - istyle.text_borders[2]), (float)(ctrl_state.curwin->rect[3]-1 - obrect[1] - istyle.text_borders[1]) );

    for(i=0;i< texta.getSize();i++){
        if (texta[i].k[16] & (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT){
            Vector<const char*> dainput;
            Vector<uint32_t> offsets;
            for(j=0;j<16;j++) if (texta[i].k[j] != 0) {
                dainput.push_back(textmap[texta[i].k[j]].text);
            }


            istyle.makeTextMesh(texta[i].d,dainput,&offsets);
            offsets.show();fflush(stdout);
            // populate uint32_t offset[2];
            for(j=0,k=0;j<16;j++){
                texta[i].k[j+17] = (offsets[k] << 2);
                if (texta[i].k[j] != 0) k++;
            }
            texta[i].k[33] = (offsets[k] << 2);
            texta[i].k[16] -= (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT;
        }
        glBindBuffer(GL_ARRAY_BUFFER, texta[i].d.glbuffer_indexes[0]);

        for(j=0;j<16;j++){
            if (texta[i].k[j] != 0) {
                pos = position_holder(texta[i].k[j], clock());
              //  glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], obrect[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , obrect[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);
              //  glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], pos[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , pos[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);

                glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], pos[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , pos[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);
                //glUniform4f(glGetUniformLocation(dashader, "uPosition"),  pos[0], pos[1], pos[2], 1.0 );
                //glUniform4f(glGetUniformLocation(dashader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], 0 , 0);
                glDrawArrays(GL_QUADS, texta[i].k[17 + j], texta[i].k[18 + j] - texta[i].k[17 + j]);
            }

        }
    }
*/
    glStencilFunc(GL_ALWAYS, 2, 0xFFFF);
}

void GUITextCloud::drawAlias(bool is_text, const int32_t* par_rect){
    if (!is_text) return;
    GUIStyle& istyle = ctrl_state.gui_styles[styleID];
    uint32_t i,j,k;
    Tuple<double, 2u> coor;

    GLint dashader = ctrl_state.datext2_shader;
    glUseProgram(dashader);
    glUniform1i(glGetUniformLocation(dashader, "tex_height"), 20); // 16);
    glUniform4f(glGetUniformLocation(dashader, "def_color"), 1.0f,0.0f,0.0f,1.0f); // 16);
    glUniform4f(glGetUniformLocation(dashader, "bgcolor_base"), 0.0f,0.0f,0.0f,0.5f); // 16);
    glUniform1f(glGetUniformLocation(dashader, "bgcolor_factor"), 0.5f); // 16);
    istyle.setFrameUniforms(dashader,rect,par_rect);

    glUniform3f(glGetUniformLocation(dashader, "id_color"),
                    (1.0f / 255) * ((GUI_alias>>24) & 255),
                    (1.0f / 255) * ((GUI_alias>>16) & 255),
                    (1.0f / 255) * ((GUI_alias>>8) & 255));
    glStencilFunc(GL_ALWAYS, GUI_alias & 0xFF, 0xFFFF);
    //printf("da alias... %X\n", GUI_alias);

    glActiveTexture(GL_TEXTURE2); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_DEFAULT_PALETTE);
    glActiveTexture(GL_TEXTURE1); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXTOFFSETS);
    glActiveTexture(GL_TEXTURE0); ctrl_state.ressourceHanddle->useTexture(GUITEXTURES_TEXT);
    Tuple<float, 3u> pos;
    glUniform2f(glGetUniformLocation(dashader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3]);
    for(i=0;i< texta.getSize();i++){break;
        if (texta.deref(i).dirty & (uint32_t)GUIFLAG_FLAGENUM_DIRTYTEXT) continue;

        glBindBuffer(GL_ARRAY_BUFFER, texta.deref(i).texta.glbuffer_indexes[0]);
        glEnableVertexAttribArray(ATTRIBUTE_POSITION);
        glEnableVertexAttribArray(ATTRIBURE_CHARID);
        glVertexAttribPointer(ATTRIBUTE_POSITION,3, GL_FLOAT, GL_FALSE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(0));
        glVertexAttribPointer(ATTRIBURE_CHARID,4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float)*3+sizeof(char) * 4 , BUFFER_OFFSET(sizeof(float)*3));

        for(j=0;j<16;j++){
            if (texta.deref(i).index[j] != 0) {
                pos = position_holder(texta.deref(i).index[j]);
              //  glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], obrect[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , obrect[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);
              //  glUniform4f(glGetUniformLocation(ctrl_state.datext_shader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], pos[0] + istyle.text_borders[0] + istyle.font_op.shift[0] , pos[1] + istyle.text_borders[0] + istyle.font_op.shift[1]);

                glUniform4f(glGetUniformLocation(dashader, "uPosition"), pos[0], pos[1], pos[2] , 1.0);
                //glUniform4f(glGetUniformLocation(dashader, "uPosition"),  pos[0], pos[1], pos[2], 1.0 );
                //glUniform4f(glGetUniformLocation(dashader, "pixscale"), 2.0f / ctrl_state.curwin->rect[2], 2.0f / ctrl_state.curwin->rect[3], 0 , 0);
                glDrawArrays(GL_QUADS, texta.deref(i).offsets[j], texta.deref(i).offsets[1 + j] - texta.deref(i).offsets[j]);
            }
        }
    }
}



GUICumulus::GUICumulus(unsigned int given_alias,unsigned int style_alias): GUIObject(given_alias,style_alias){
}
void GUICumulus::draw(bool mouse_over, const int32_t* par_rect){
}
void GUICumulus::drawAlias(bool is_text, const int32_t* par_rect){
}
GUIMSG_Enum  GUICumulus::processGUIevent(const GUIMSG_Enum event){
    return GUIMSG_NULL;
}

FormatedText::FormatedText():txt(NULL){}
FormatedText::~FormatedText(){delete[](txt);}
void FormatedText::draw(){
    unsigned char* cur = txt;
    char mode = 1; // mode =0 to exit;
    while(mode){
        switch(mode){
        case 1:
            while(mode == 1){
                switch(*(cur++)){
                    case '\0': mode =0; break;
                    case ' ': glTranslatef((float)spacesize,0,0); break;
                    case '\\':
                        switch(*(cur++)){
                            case 's': spacesize = *(cur++); break;
                            default:
                                printf("unrecognized text format flag '\\%c'\n", cur[-1]);
                        }
                    break;
                    case '<': mode =2; break;
                }
            }
        break;
        case 2:
            while(mode == 2){
                mode =0;
            }
        break;
        case 3:
            while(mode == 3){
                mode =0;
            }
        break;
        }
    }
}
bool waitingForSuccess(ERRCODE output){if (output == 0) return false; SDL_Delay(1); return true;}



void EventQueue<void>::runTo(uint32_t n_time){
    if (getSync()){
        proc_async_routine();
        while( (p_queue.isEmpty() == false)&&((time = p_queue.top().time) - n_time > 0xF0000000)){
                Event< >* nev = (Event< >*)p_queue.pop().ev;
                uint32_t o = (*nev)();
                if (o){
                    if (o != 0xFFFFFFFF) this->insert_sync(time + o,nev);
                }else delete(nev);
                proc_async_routine();
            }
        time = n_time;
        freeSync();
    }
}

void EventQueue<void>::startQueueThread(){}


void EventQueue<void>::proc_async_routine(){
    while(async_write != async_read){
        p_queue.insert(buffer[async_read]);
        async_read = (async_read + 1) & 255;
    }
}

void EventQueue<void>::insert_async(uint32_t time, Event< >* ev){
    buffer[async_write] = EventUnit(time,ev);
    while(((async_write + 1) & 255) == async_read) {
        if (getSync()){
            this->proc_async_routine();
        }else {printf("ASYNC SPIKE!\n"); exit(1);}
        freeSync();
    }
    async_write = (async_write + 1) & 255;
}

void EventQueue<void>::insert_asap(Event< >* ev){
    buffer[async_write] = EventUnit(time,ev);
    while(((async_write + 1) & 255) == async_read) {
        if (getSync()){
            this->proc_async_routine();
        }else {printf("ASYNC SPIKE!\n"); exit(1);}
        freeSync();
    }
    async_write = (async_write + 1) & 255;
}


bool EventQueue<void>::pop(Event< >*& fout){ // always pop!
    this->proc_async_routine();
    if (p_queue.isEmpty()) return false;
    fout = (Event< >*) p_queue.pop().ev;
    return true;
    }

bool EventQueue<void>::pop_exec(){ // always pop!
    this->proc_async_routine();
    if (p_queue.isEmpty()) return false;
    Event< >* tmp = (Event< >*)p_queue.pop().ev;
    uint32_t i = (*tmp)();
    if (i+1 > 1) p_queue.insert(EventUnit(i,tmp));
    else if (i == 0) delete(tmp);
    return true;
    }

} // namespace end
