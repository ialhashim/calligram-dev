#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <iostream>
#include "GLWindow.h"


using namespace std;


void Init();
void Key(unsigned char key, int x, int y);
void Draw( void );
void Resize( const int width, const int height );
void Mouse( int button, int state, int x, int y );
void MouseMotion(int x, int y );
void SpecialKey(int key, int x, int y);



/** Not object oriented at all.. **/
static GLWindow glutWindow;

int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutInitWindowSize(640,480);
  glutInitWindowPosition(200, 150); 
  glutCreateWindow("CS176 HW3");

  // Glut window set image
  // Clut window set grid resolution
  glutWindow.SetGrid(100);
  glutWindow.LoadTexture("frog1.png");
  Init();
  glutReshapeFunc(Resize);
  glutKeyboardFunc(Key);
  glutMouseFunc(Mouse);
  glutMotionFunc (MouseMotion);
  glutDisplayFunc(Draw);
  glutMainLoop();
  return 0;             
}


void Init() {
   glutWindow.Init();
}

void Key(unsigned char key, int x, int y) {
   glutWindow.Key(key, x, y);
}

void Draw( void ) {
   glutWindow.Draw();
}

void Resize( const int width, const int height ) {
   glutWindow.Resize(width, height);
}

void Mouse( int button, int state, int x, int y ) {
   glutWindow.Mouse(button, state, x, y);
}

void MouseMotion( int x, int y ) {
   glutWindow.MouseMotion(x, y);
}
