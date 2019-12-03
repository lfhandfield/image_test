#include "main.h"

DefaultRessourceLoader::DefaultRessourceLoader(){
    this->loadTexture("/opt/display/Images/ft04data.tif", textures[0], 8);
    this->loadTexture("/opt/display/Images/ft04coor.tif", textures[1], 8);
    this->loadTexture("/opt/display/Images/guisbars.tif", textures[2], 0);
    this->loadTexture("/opt/display/Images/guibuttn.tif", textures[3], 0);
    this->loadTexture("/opt/display/Images/guiicons.tif", textures[4], 0);
    this->loadTexture("/opt/display/Images/palette0.tif", textures[5], 0);
    this->loadTexture("/opt/display/Images/guitexts.tif", textures[6], 0);
    this->loadTexture("/opt/display/Images/guibgrid.tif", textures[7], 0);
    this->loadTexture("/opt/display/Images/guigrmsk.tif", textures[8], 0);
    this->loadTexture("/opt/display/Images/gtoolbar.tif", textures[9], 0);
}

void DefaultRessourceLoader::loadTexture(const char* path, GLuint& slot, int flag){
    LFHPrimitive::TiffFile tf(path, 'r');
	LFHPrimitive::DataGrid<unsigned char, 3> im;
    Tuple<unsigned int,3> coooo; // tmp

	glGenTextures(1, &slot);
    glBindTexture(GL_TEXTURE_2D,slot);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    if (flag & 1){
        if (!tf.fetch(im)) {fprintf(stderr,"Found no frames in %s\n", path); exit(1);};
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, im.dims[1], im.dims[2], 0, GL_RED_INTEGER, GL_UNSIGNED_INT, im.data);
    }else{
        if (!tf.fetch(im)) {fprintf(stderr,"Found no frames in %s\n", path); exit(1);};
        if (im.dims[2] > 256){
            Tuple<unsigned int,3> ndim; GLint nsize;
            glGetIntegerv(GL_MAX_TEXTURE_SIZE, &nsize);
            if (im.dims[2] > nsize){
                 ndim[0] = im.dims[0]; ndim[2] = nsize; ndim[1] = ndim[2];
                im.toresize_crude(ndim);
            }
        }
	   

        switch(im.dims[0]){
        case 3:
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, im.dims[1], im.dims[2], 0, GL_RGB, GL_UNSIGNED_BYTE, im.data);
        break;
        case 4:
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, im.dims[1], im.dims[2], 0, GL_RGBA, GL_UNSIGNED_BYTE, im.data);
        break;
        }
    }


    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    if (flag & 8){
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }else{
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    }

}
void DefaultRessourceLoader::useTexture(const uint32_t which){glBindTexture(GL_TEXTURE_2D, textures[(int) which]);}
void DefaultRessourceLoader::useTextureExt(uint32_t ID){}
void DefaultRessourceLoader::allocTextureExt(uint32_t ID){}
void DefaultRessourceLoader::deallocTextureExt(uint32_t ID){}
void DefaultRessourceLoader::useSound(const uint32_t){}

int main(int argc, char** argv){
	Task task;
task(argc, argv);}

Task::Task(){}

void Task::nbaddtoken(char const * const token, int& min, int& max){
    switch(*token){
    case '\0':min =0; break;
    case 'o': min =1; break;
    }
}
void Task::store(char* const * token, int nbtoken){
    /*switch(*token){
    case 'o': min =1; break;
    }*/
}
void Task::help(){
        printf("This makes a interactive 3D overlay.\n");
        printf("Arguments: [2]\n");
        printf("        (in file) input tif image\n");
        printf("        (in file) output tif image\n");
        printf("        (out file) output tif image\n");

        printf("\n");
        printf("Flags\n\n");
        printf("\t-s:\tShow output in Preview (works on MacOS)\n");
        printf("\t-o (FILE f='output.tif') : path for output.\n");
        printf("Version 1.0\n");

}
int Task::OnMaintain(){//OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
    dawin->render();
return(0);}
int Task::OnKeyDown(const SDL_KeyboardEvent& event){ //OnKeyDown(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
    /*switch(event.keysym.sym){
    case '3':{ //:

    }break;//:
    default:
    ffr->makeDefault();
    }*/
return(0);}
int Task::OnKeyUp(const SDL_KeyboardEvent& event){ //OnKeyUp(Event.key.keysym.sym,Event.key.keysym.mod,Event.key.keysym.unicode);
return(0);}
int Task::listen(LFHDisplay::GUImessage& msg){
return(0);}
void Task::draw(LFHDisplay::MyWindow*){
}
void Task::drawAlias(LFHDisplay::MyWindow*){
}
int Task::defstore(char* const * token, int nbtoken){
    if (!LFHDisplay::Controlstate::init_SDL(NULL,NULL)) return 1;
    if (!LFHDisplay::Controlstate::init_openGL()) return 1;
    LFHDisplay::ctrl_state << this;
    dawin = new LFHDisplay::MyWindow(1u,0u,1024,768,LFHDisplay::RELPOS_RIGHT,false, 0);
    LFHDisplay::ctrl_state.curwin = dawin;
    (*LFHDisplay::ctrl_state.curwin) << this;
    DefaultRessourceLoader* resl = new DefaultRessourceLoader();
    LFH_ALIVE;
    LFHDisplay::ctrl_state.main_control_loop(resl);
    LFH_ALIVE;
    delete(resl);
    LFH_ALIVE;
return 0;}



