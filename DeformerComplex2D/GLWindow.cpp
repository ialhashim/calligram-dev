#include <iostream>
#include <math.h>
#include "GLWindow.h"
#include "pngLoad.h"

#include <QImage>

using namespace std;


void GLWindow::Init() {
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glClearColor( .8, .8, .8, 0 );
   glClearDepth( 1 );
   scale = 1;
   xTranslate = 0;
   yTranslate = 0;
   zooming = false;
   translating = false;
   movingPolygon = false;
   creatingPolygon = false;

   grid.isCauchy = true;
}

void GLWindow::Key(unsigned char key, int x, int y) {
   switch (key) {
      case 'w':
         grid.drawStyle = Grid::DRAW_WIRE;
         break;

      case 't':
         grid.drawStyle = Grid::DRAW_TEXTURE;
         break;

      case 'c':
         grid.ClearPlygon();
         break;

      case 'r':
         grid.ResetPolygon();
         break;

      case 'm':
         grid.toggleEditMode();
         break;

      case 'h':
         grid.hideControls = !grid.hideControls;
         break;

       case 's':
         grid.isCauchy = !grid.isCauchy;
         grid.ResetPolygon();
         break;
   }

   glutPostRedisplay();
}


void GLWindow::Mouse( int button, int state, int x, int y ) {
   if (state != GLUT_DOWN) {
      translating = false;
      zooming = false;
      movingPolygon = false;
      creatingPolygon = false;
      addingHandle = false;
      movingHandle = false;

      oldX = oldY = -1;
      return;
   }

   if(!(glutGetModifiers() & GLUT_ACTIVE_SHIFT)) { 
      if (grid.editMode == DEFORM_CAGE) {
          if (grid.isCauchy && button == GLUT_LEFT_BUTTON)
             movingPolygon = true;

          else if (button == GLUT_LEFT_BUTTON)
             movingHandle = true;
      }
      else  {
          if (button == GLUT_LEFT_BUTTON)
             creatingPolygon = true;

          if ((!grid.isCauchy) && (button == GLUT_RIGHT_BUTTON))
             addingHandle = true;
      }
   }

   if(glutGetModifiers() & GLUT_ACTIVE_SHIFT) { 
      if (button == GLUT_LEFT_BUTTON) {
         translating = true;
      }

      if (button == GLUT_RIGHT_BUTTON) {
         zooming = true;
      }
   }

   float myX, myY;

   myX =  (float)(x*2.0 - width)/(float)scale - xTranslate;
   myY = -(float)(y*2.0 - height)/(float)scale - yTranslate;

   if (creatingPolygon) {
      grid.AddPolygonPoint(myX, myY);
      glutPostRedisplay();
   }

   if (addingHandle) {
      grid.AddHandle(myX, myY);
      glutPostRedisplay();
   }

   oldX = myX;
   oldY = myY;
}

void GLWindow::MouseMotion( int x, int y ) {
   float myX, myY;
   myX =  (float)(x * 2.0 - width)/(float)scale - xTranslate;
   myY = -(float)(y * 2.0 - height)/(float)scale - yTranslate;

   if (oldX >= 0 && oldY >= 0) {

      if (movingPolygon) {
         grid.MovePolygonPoint(oldX, oldY, (myX - oldX), (myY-oldY));
         glutPostRedisplay();
      }

      else if (movingHandle) {
         grid.MoveHandle(oldX, oldY, (myX - oldX), (myY-oldY));
         glutPostRedisplay();
      }

      else if (translating) {
         xTranslate += (myX - oldX);
         yTranslate += (myY - oldY);
         glutPostRedisplay();
      }

      else if (zooming) {
         if ((myY-oldY) > 0)
            scale *= 1.05;
         else 
            scale /= 1.05;
         glutPostRedisplay();
      }
   }

   oldX =  (float)(x * 2.0 - width)/(float)scale - xTranslate;
   oldY = -(float)(y * 2.0 - height)/(float)scale - yTranslate;
}

void GLWindow::Resize(int _width, int _height) {
   if( _width <= 0 || _height <= 0 ) return;
   width = _width;
   height = _height;
   glViewport(0,0, _width, _height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   int vPort[4];

   glGetIntegerv(GL_VIEWPORT, vPort);

   glOrtho(-vPort[2], vPort[2], -vPort[3], vPort[3], -1000, 1000);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   scale = _height*1.75;
   xTranslate = -.5;
   yTranslate = -.5;
}

void GLWindow::Draw() {
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
   glLoadIdentity();

   glScalef(scale, scale, 1);
   glTranslatef(xTranslate, yTranslate, 0);

   glLineWidth(1);
   glColor3f(.3f, .3f, .3f);

   glBegin(GL_LINE_LOOP);
      glVertex2f (.0f,  .0f);
      glVertex2f (.0f, 1.0f);
      glVertex2f(1.0f, 1.0f);
      glVertex2f(1.0f, 0.0f);
   glEnd();

   grid.Draw();

   glColor3f(0.0f, 0.0f, 0.0f);

   if (grid.isCauchy)
       DrawString("Cauchy-Green Coordinates, ", false, 0.0,  -.03);
   else
       DrawString("P2P Cauchy-Green Coordinates, ", false, 0.0, -.03);

   if (grid.editMode == CONSTRUCT_CAGE)
       DrawString("Control Point Insert Mode", true, 0, 0);
   else
       DrawString("Deformation Mode", true, 0, 0);

   glFlush();
   glutSwapBuffers();
}

void GLWindow::DrawString(char *aString, bool append, float x, float y)    {
    int i = 0;
    void *font = GLUT_BITMAP_HELVETICA_12;

    if (!append)
        glRasterPos2f(x, y);

    while (char c = aString[i++]) {
       glutBitmapCharacter(font, c);
    }
}

bool GLWindow::LoadTexture(char * filename) {
   unsigned long int w, h;
   unsigned char *imageData;
   //pngLoad(filename, &w, &h, &imageData);



   GLfloat borderColour[] = { 1.0f, 1.0f, 1.0f, 1.0f };   
   glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColour);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

   //glTexImage2D(GL_TEXTURE_2D, 0, 4, tex.width(), tex.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.bits());
   // NOTE: DOESN'T WORK WITH ALPHA CHANNELS
   //glTexImage2D(GL_TEXTURE_2D, 0, 4, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, imageData);

   QImage image(filename);
   //glTexImage2D(GL_TEXTURE_2D, 0, 4, image.width(), image.height(), 0, GL_RGB, GL_UNSIGNED_BYTE, image.bits() );
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width(), image.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, image.bits());

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);

   grid.drawStyle = Grid::DRAW_TEXTURE;

   return true;
}

