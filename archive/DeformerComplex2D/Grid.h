#ifndef GRID
#define GRID
#include "Vector2.h"

#define LAMBDA .005

typedef enum {CONSTRUCT_CAGE = 0, DEFORM_CAGE = 1} EditMode;
class Grid
{
private:
   int size;
   Vector2 ** points;
   bool ** points_valid;

   Vector2 ** barycentricCoords;
   Vector2 *** P_Matrix;
   Vector2 *** C_Matrix;

   std::vector<Vector2> polygon;
   std::vector<Vector2> movedPolygon;
   std::vector<Vector2> handles;
   std::vector<Vector2> movedHandles;

public:
   static const int DRAW_WIRE = 0;
   static const int DRAW_TEXTURE = 1;
   EditMode editMode;
   bool hideControls;
   bool isCauchy;
   int drawStyle;

   Grid();
   void Resize(int n);

   void Draw();
   void toggleEditMode();
   void AddPolygonPoint(const float & x, const float & y);
   void MovePolygonPoint(const float & x, const float & y, const float & dx, const float & dy);

   void AddHandle(const float & x, const float & y);
   void MoveHandle(const float & x, const float & y, const float & dx, const float & dy);

   void ClearPlygon() {
       handles.clear();
       polygon.clear();
       movedPolygon.clear();
       movedHandles.clear();

       updatePointsValid();
   }

   void ResetPolygon() {
       for (unsigned int i = 0; i < polygon.size(); ++i)    {
           movedPolygon[i] = polygon[i];
       }

       for (unsigned int i = 0; i < handles.size(); i++)    {
           movedHandles[i] = handles[i];
       }
   }

   void updatePointsValid();
   void computeBarycentricCoordinates();
   Vector2 ComputeC_j(int j, Vector2 z) const;
   Vector2 ComputeD_j(int j, Vector2 z) const;

   bool isPointInPolygon(const Vector2 &p) const;
   int getNumWindsAroundPoint(const Vector2 &p) const;
   double distanceLeftFromLine(const Vector2 &thePoint, const Vector2 &v1, const Vector2 &v2) const;
   ~Grid(void);
};

#endif
