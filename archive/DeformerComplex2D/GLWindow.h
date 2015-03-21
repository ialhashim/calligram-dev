#ifndef GL_WINDOW
#define GL_WINDOW

#include <stdlib.h>
#include <GL/glut.h>
#include "Grid.h"

class GLWindow
{

private:
   int width;
   int height;
   float scale;
   float xTranslate;
   float yTranslate;
   float oldX;
   float oldY;

   bool zooming;
   bool translating;
   bool movingPolygon;
   bool creatingPolygon;
   bool addingHandle;
   bool movingHandle;

   Grid grid;
public:  

   void Init();
   void SetGrid (int gridSize) {grid.Resize(gridSize); }
   void Key(unsigned char key, int x, int y);
   void Draw( void );
   void Resize( const int width, const int height );
   void Mouse( int button, int state, int x, int y );
   void MouseMotion( int x, int y );
   void DrawString(char *aString, bool append, float x, float y);
   bool LoadTexture(char * textureFile);

};


#endif
