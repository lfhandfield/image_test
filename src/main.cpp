#include "main.h"

DefaultRessourceLoader::DefaultRessourceLoader(){
    this->loadTexture("/opt/display/Images/ft04data.tif", textures[0], 8);
    this->loadTexture("/opt/display/Images/ft04coor.tif", textures[1], 8);
    this->loadTexture("/opt/display/Images/guisbars.tif", textures[2], 0);
    this->loadTexture("/opt/display/Images/guibuttn.tif", textures[(int)GUITEXTURES_BUTTON], 0);
    this->loadTexture("/opt/display/Images/guiicons.tif", textures[4], 0);
    this->loadTexture("/opt/display/Images/palette0.tif", textures[5], 0);
    this->loadTexture("/opt/display/Images/guitexts.tif", textures[6], 0);
    this->loadTexture("/opt/display/Images/guibgrid.tif", textures[7], 0);
    this->loadTexture("/opt/display/Images/guigrmsk.tif", textures[8], 0);
    this->loadTexture("/opt/display/Images/gtoolbar.tif", textures[9], 0);
}

void DefaultRessourceLoader::loadTexture(char* path, GLuint& slot, int flag){

    LFHPrimitive::TiffFile tf(path);
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
void DefaultRessourceLoader::useTexture(const GUITEXTURES_enum which){glBindTexture(GL_TEXTURE_2D, textures[(int) which]);}
void DefaultRessourceLoader::useTextureExt(unsigned int ID){}
void DefaultRessourceLoader::allocTextureExt(unsigned int ID){}
void DefaultRessourceLoader::deallocTextureExt(unsigned int ID){}
void DefaultRessourceLoader::useSound(const GUISOUND_enum){}


int main(int argc, char** argv){
	Task task;
task(argc, argv);}
