#include <stdlib.h>
#include <GL/glut.h>
#include "Grid.h"
#include "allocMultidim.h"
#include "MatrixStripped.hh"
#include "matrix.hh"

using namespace std;

Grid::Grid()

{
   points = NULL;
   P_Matrix = NULL;
   C_Matrix = NULL;
   points_valid = NULL;

   editMode = CONSTRUCT_CAGE;

   size = 0;
}

void Grid::toggleEditMode()
{
    if (editMode == CONSTRUCT_CAGE) {
        editMode = DEFORM_CAGE;

        computeBarycentricCoordinates();
        updatePointsValid();  
    }
    else
        editMode = CONSTRUCT_CAGE;

    ResetPolygon();
}

void Grid::AddPolygonPoint(const float & x, const float & y)
{
   if (x >= 0 && x <= 1 && y >= 0 && y <= 1) {
      polygon.push_back(Vector2(x, y));
      movedPolygon.push_back(Vector2(x, y));
   }
}

void Grid::AddHandle(const float &x, const float &y)
{
   if (x >= 0 && x <= 1 && y >= 0 && y <= 1) {
      handles.push_back(Vector2(x, y));
      movedHandles.push_back(Vector2(x, y));
   }
}

/* Recomputes what meshpoints are within the source polygon */
void Grid::updatePointsValid()
{
    int i, j;
    if (polygon.size() != 0)    {
       for(i = 0; i < size; i++)    {
           for (j = 0; j < size; j++)   {
               points_valid[i][j] = isPointInPolygon(points[i][j]);
           }
       }
    }

    else    { /* Draw all points when there is no source/target polygon */
       for(i = 0; i < size; i++)    {
           for (j = 0; j < size; j++)   {
               points_valid[i][j] = true;
           }
       }
    }
}

/* Computes the barycentric coordinates for the specified source polygon */
/* Uses regular Cauchy-Green coordinates if isCauchy == 1,
 * Uses P2P Cauchy-Green coordinates if isCauchy == 0. */
void Grid::computeBarycentricCoordinates()
{
    if (P_Matrix != NULL)
       free(P_Matrix);

    if (C_Matrix != NULL)
       free(C_Matrix);

    int n = polygon.size();
    P_Matrix = (Vector2 ***)allocMultidim(sizeof(Vector2), 3, size, size, n);
    C_Matrix = (Vector2 ***)allocMultidim(sizeof(Vector2), 3, size, size, n);

    int row, col;
    for (row = 0; row < size; row++) {
        for (col = 0; col < size; col++) {
            Vector2 z = points[row][col];

            int j;
            for (j = 0; j < n; j++) {
                C_Matrix[row][col][j] = ComputeC_j(j, z);            
            }
        }
    }
    
    if (!isCauchy)   {
        // Do the P2P coordinate computation
        //
        // Sample w values at the midpoints of cage edges
        // Construct augmented matrix A = [[C][lambda * D]
        //
        // A is (p + k) x n, where p is the number of control points
        //  and k is the number of samples. I am going to make one sample
        //  per cage edge, so k = n = polygon.size().
        //
        // C, in this case, is a pxn matrix evalulated at the original
        //  (untransformed control points): C_j(r_i)
        //
        // Pseudo-invert A
        int p = handles.size();
        int n = polygon.size();
        int k = n;

        MatrixStripped<Vector2> A(p + k, n);

        /* Insert the [C] part of the A matrix */
        int i, j;
        for (i = 0; i < p; i++) {
            Vector2 z = handles[i];
            for (j = 0; j < n; j++) {
                A(i, j) = ComputeC_j(j, z);
            }
        }

        /* Insert the [lambda * D] part of the A matrix */
        for (i = 0; i < k; i++) {
            /* Choose ws to be midpoints of the cage edges */
            /* w_i is the midpoint of edge (polygon[i], polygon[i + 1])*/
            Vector2 pt = polygon[i];
            Vector2 nextPt = polygon[(i + 1) % n];

            /* Average x coordinates and y coordinates to get midpoint */
            Vector2 w((pt.x() + nextPt.x()) / 2.0, (pt.y() + nextPt.y()) / 2.0);

            for (j = 0; j < n; j++) {
                // Compute d_j(w), where w is a midpoint 
                A(p + i, j) = LAMBDA * ComputeD_j(j, w);
            }
        }

        /* Compute Pseudoinverse of A:
         * A^+ = ((A ^ t * A) ^ -1) * (A ^ t) */
        /* Note, A^t is actually the complex conjugate of the "regular" transpose
         * that MatrixStripped performs */
        MatrixStripped<Vector2> A_transpose(A.transpose());
        
        for (i = 0; i < A_transpose.getrows(); i++) {
            for (j = 0; j < A_transpose.getcols(); j++)
                A_transpose(i, j) = A_transpose(i, j).conjugate();
        }
        
        MatrixStripped<Vector2> A_A_transpose(A_transpose * A);

        /* Now, import into a matrix object for inversion */
        matrix<Vector2> matrixForInversion(n, n);

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                matrixForInversion.setvalue(i, j, A_A_transpose(i, j));
        }

        /* Use the matrix class to invert the matrix */
        matrixForInversion.invert();

        /* Import back into my stripped matrix class that is easier to use */
        MatrixStripped<Vector2> N_computation(n, n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                bool success;
                Vector2 N_value;
                matrixForInversion.getvalue(i, j, N_value, success);
                N_computation(i, j) = N_value;
                assert(success);
            }
        }

        N_computation *= A_transpose;
            
        /* Now, N is the first p columns of N_computation
         * So P can be conputed:
         * P_j = Sum_k(C_k(z) * N_k,j) */
        for (row = 0; row < size; row++)    {
            for (col = 0; col < size; col++)    {
                for (j = 0; j < p; j++) {
                    int k_iter;
                    Vector2 value(0, 0);
                    for (k_iter = 0; k_iter < n; k_iter++)
                        value += C_Matrix[row][col][k_iter] * N_computation(k_iter, j);

                    P_Matrix[row][col][j] = value;
                }
            }
        }
    }
}

Vector2 Grid::ComputeC_j(int j, Vector2 z) const
{
    int n = polygon.size();
    int prevIndex = j - 1;
    if (prevIndex < 0)
        prevIndex += n;

    int nextIndex = j + 1;
    if (nextIndex >= n)
        nextIndex -= n;

    Vector2 z_j_m = polygon[prevIndex];
    Vector2 z_j   = polygon[j];
    Vector2 z_j_p = polygon[nextIndex];
    
    Vector2 A_j_p = z_j_p - z_j;
    Vector2 A_j   = z_j   - z_j_m;

    Vector2 B_j_p = z_j_p - z;
    Vector2 B_j   = z_j   - z;
    Vector2 B_j_m = z_j_m - z;

    return Vector2(0, - 1 / (2 * M_PI)) *
        ((B_j_p / A_j_p) * ((B_j_p / B_j).c_log()) -
        (B_j_m / A_j) * ((B_j / B_j_m).c_log()));
}

Vector2 Grid::ComputeD_j(int j, Vector2 z) const
{
    int n = polygon.size();
    int prevIndex = j - 1;
    if (prevIndex < 0)
        prevIndex += n;

    int nextIndex = j + 1;
    if (nextIndex >= n)
        nextIndex -= n;

    Vector2 z_j_m = polygon[prevIndex];
    Vector2 z_j   = polygon[j];
    Vector2 z_j_p = polygon[nextIndex];
    
    Vector2 B_j_p = z_j_p - z;
    Vector2 B_j   = z_j   - z;
    Vector2 B_j_m = z_j_m - z;

    Vector2 coeff = Vector2(0, - 1 / (2 * M_PI));
    return (coeff / (B_j_m * B_j)) - (coeff / (B_j * B_j_p));
}

/************************* isPointInPolygon *************************
* Uses the winding number to determine if a point is in the transformed
* polygon or not. If the winding number is greater than zero then it
* is in the polygon.
********************************************************************/
bool Grid::isPointInPolygon(const Vector2 &aPoint) const
{
    return (getNumWindsAroundPoint(aPoint) != 0);   
}
 
/********************** getNumWindsAroundPoint **********************
* Determines the "winding number," the number of times the polygon
* "winds around" the specified point. If this winding number is
* nonzero, it means thePoint is within theTrace.
********************************************************************/
int Grid::getNumWindsAroundPoint(const Vector2 &p) const
{
    int winds = 0;
    int i;

    int n = movedPolygon.size();
    for (i = 0; i < n; i++)    {
        Vector2 curr = polygon[i];
        int next_i = i + 1;

        if (next_i >  n)
            next_i -= n;

        Vector2 next = polygon[(i + 1) % movedPolygon.size()];
        
        if (curr.y() <= p.y())  {/* Edge starts below point */
            /* If edge crosses to above the point with the point at its left */
            /* (Counterclockwise around point) */
            winds += (next.y() > p.y() && distanceLeftFromLine(p, curr, next) > 0);
        }

        else    {/* Edge starts above point */
            /* If edge crosses to below the point with the point at its right */
            /* (Clockwise around point) */
            winds -= (next.y() <= p.y() && distanceLeftFromLine(p, curr, next) < 0);
        }
    }

    return winds;
}

/*********************** distanceLeftFromLine ***********************
* Function use to determine whether thePoint is to the left or the
* right of--or on--the line running through v1 and v2. Directions
* are relative to the orientation of the line (v1, v2).
* Case > 0:
*   Point is to the "left" of the line.
* Case = 0:
*   Point is on the line.
* Case < 0:
*   Point is to the "right" of the line.
********************************************************************/
double Grid::distanceLeftFromLine(const Vector2 &thePoint, const Vector2 &v1, const Vector2 &v2) const
{
    return (thePoint.x() - v1.x()) * (v2.y() - v1.y()) -
        (v2.x() - v1.x()) * (thePoint.y() - v1.y());
}

void Grid::MovePolygonPoint(const float & x, const float & y, const float & dx, const float & dy)
{
   if (movedPolygon.size() < 1) return;
   
   int closest = 0;
   Vector2 v(x, y);
   float dist = (movedPolygon[0] - v).dotItself();

   for (unsigned int i = 1; i < movedPolygon.size(); ++i) {
      if ((movedPolygon[i] - v).dotItself() < dist) {
         closest = i;
         dist = (movedPolygon[i] - v).dotItself();
      }
   }

   movedPolygon[closest].x() += dx;
   movedPolygon[closest].y() += dy;
}

void Grid::MoveHandle(const float & x, const float & y, const float & dx, const float & dy)
{
   if (movedHandles.size() < 1) return;
   
   int closest = 0;
   Vector2 v(x, y);
   float dist = (movedHandles[0] - v).dotItself();

   for (unsigned int i = 1; i < movedHandles.size(); ++i) {
      if ((movedHandles[i] - v).dotItself() < dist) {
         closest = i;
         dist = (movedHandles[i] - v).dotItself();
      }
   }

   movedHandles[closest].x() += dx;
   movedHandles[closest].y() += dy;
}


void Grid::Resize(int n)
{
   if (points) {

      for (int i = 0; i < size; ++i)    {
         delete [] points[i];
         delete [] barycentricCoords[i];
      }

      delete [] points;
      delete [] barycentricCoords;
      free(points_valid);
   }

   size = n;

   points = new Vector2 *[n];
   points_valid = (bool **)allocMultidim(sizeof(bool), 2, size, size);

   barycentricCoords = new Vector2 *[n];

   for (int i = 0; i < size; ++i)   {
      points[i] = new Vector2[size];
      barycentricCoords[i] = new Vector2[size];
   }

   float dx = 1.0/(float)(size-1);

   for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
         points[i][j].x() = i*dx;
         points[i][j].y() = j*dx;
      }
   }

   drawStyle = DRAW_WIRE;
   updatePointsValid();
}

void Grid::Draw()
{
   Vector2 **transformed_pts = (Vector2 **)allocMultidim(sizeof(Vector2), 2, size, size);

   int row, col;
   bool isCageConstruction = (editMode == CONSTRUCT_CAGE);

   int n = polygon.size();
   int p = handles.size();
   // Transform all the points
   for (row = 0; row < size; row++)    {
        for (col = 0; col < size; col++)    {
           if (polygon.size() > 0 && (!isCageConstruction)) {
               Vector2 value = Vector2(0, 0);

               int j;
               if (isCauchy)    {
                   for (j = 0; j < n; j++)
                       value += C_Matrix[row][col][j] * movedPolygon[j];
               }
               else {
                   for (j = 0; j < p; j++)
                       value += P_Matrix[row][col][j] * movedHandles[j];
               }

               transformed_pts[row][col] = value;
           }

           else
               transformed_pts[row][col] = points[row][col];
        }
    }

   glLineWidth(.5f);

   if (drawStyle == DRAW_WIRE) {
      glDisable(GL_TEXTURE_2D);
      glColor3f(.1, .4, .1);

      for (int i = 0; i < size - 1; ++i) {
         for (int j = 0; j < size - 1; ++j) {
            /* only draw polygons with points within the source polygon */
            if (isCageConstruction || (points_valid[i][j] && points_valid[i + 1][j] && points_valid[i][j + 1]))   {
                glBegin(GL_LINE_LOOP);

                Vector2 v1 = transformed_pts[i][j];
                Vector2 v2 = transformed_pts[i + 1][j];
                Vector2 v3 = transformed_pts[i][j + 1];

                glVertex2f(v1.x(), v1.y());
                glVertex2f(v2.x(), v2.y());
                glVertex2f(v3.x(), v3.y());
                glEnd();
             }

            /* only draw polygons with points within the source polygon */
            if (isCageConstruction || (points_valid[i + 1][j + 1] && points_valid[i + 1][j] && points_valid[i][j + 1]))   {
                glBegin(GL_LINE_LOOP);

                Vector2 v1 = transformed_pts[i + 1][j];
                Vector2 v2 = transformed_pts[i + 1][j + 1];
                Vector2 v3 = transformed_pts[i][j + 1];

                glVertex2f(v1.x(), v1.y());
                glVertex2f(v2.x(), v2.y());
                glVertex2f(v3.x(), v3.y());

                glEnd();
            }
         }
      }
   }

   else {
      glColor3f(1, 1, 1);

      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glEnable(GL_TEXTURE_2D);
      glMatrixMode(GL_TEXTURE);
      glLoadIdentity();	

      glMatrixMode(GL_MODELVIEW);

      for (int i = 0; i < size-1; ++i) {
         for (int j = 0; j < size-1; ++j) {

            /* only draw polygons with points within the source polygon */
            if (isCageConstruction || (points_valid[i][j] && points_valid[i + 1][j] && points_valid[i][j + 1]))   {
                glBegin(GL_POLYGON);

                Vector2 v1 = transformed_pts[i][j];
                Vector2 v2 = transformed_pts[i + 1][j];
                Vector2 v3 = transformed_pts[i][j + 1];

                glTexCoord2f(points[i][j].x(), points[i][j].y());
                glVertex2f(v1.x(), v1.y());

                glTexCoord2f(points[i+1][j].x(), points[i+1][j].y());
                glVertex2f(v2.x(), v2.y());

                glTexCoord2f(points[i][j+1].x(), points[i][j+1].y());
                glVertex2f(v3.x(), v3.y());

                glEnd();
            }

            /* only draw polygons with points within the source polygon */
            if (isCageConstruction || (points_valid[i + 1][j + 1] && points_valid[i + 1][j] && points_valid[i][j + 1]))   {
                glBegin(GL_POLYGON);

                Vector2 v1 = transformed_pts[i + 1][j];
                Vector2 v2 = transformed_pts[i + 1][j + 1];
                Vector2 v3 = transformed_pts[i][j + 1];

                glTexCoord2f(points[i+1][j].x(), points[i+1][j].y());
                glVertex2f(v1.x(), v1.y());

                glTexCoord2f(points[i+1][j+1].x(), points[i+1][j+1].y());
                glVertex2f(v2.x(), v2.y());

                glTexCoord2f(points[i][j+1].x(), points[i][j+1].y());
                glVertex2f(v3.x(), v3.y());

                glEnd();
            }
         }
      }
   }

   free(transformed_pts);

   glDisable(GL_TEXTURE_2D);

   glLineWidth(1.7);

   glPointSize(4.2);

   glColor3f(.7f, .7f, .8f);

   if (editMode == CONSTRUCT_CAGE && !hideControls)  {
       glBegin(GL_LINE_LOOP);

       /* Draw cage edges */
       for (unsigned int i = 0; i < polygon.size(); ++i) {
          glVertex2f(polygon[i].x(), polygon[i].y());
       }

       glEnd();

       glColor3f(.25f, .25f, .65f);

       /* Draw cage vertices */
       glBegin(GL_POINTS);

       for (unsigned int i = 0; i < polygon.size(); ++i)
          glVertex2f(polygon[i].x(), polygon[i].y());

       /* Draw handles */
       if (!isCauchy)   {
           glColor3f(0.0f, 0.0f, 1.0f);
           for (unsigned int i = 0; i < handles.size(); ++i)    {
              glVertex2f(handles[i].x(), handles[i].y());
           }
       }

       glEnd();
   }

   if (editMode == DEFORM_CAGE && !hideControls) {
       glColor3f(.8f, .7f, .8f);


       if (isCauchy)    {

           /* Draw transformed cage edges */
           glBegin(GL_LINE_LOOP);
           for (unsigned int i = 0; i < polygon.size(); ++i)
              glVertex2f(movedPolygon[i].x(), movedPolygon[i].y());

           glEnd();

           glColor3f(.85f, .25f, .65f);
           glBegin(GL_POINTS);

           /* Draw transformed cage vertices */
           for (unsigned int i = 0; i < polygon.size(); ++i) {
              glVertex2f(movedPolygon[i].x(), movedPolygon[i].y());
           }

           glEnd();
       }

       else {
           /* Draw transformed handles */
           glColor3f(1.0f, 0.0f, 0.0f);

           glBegin(GL_POINTS);

           for (unsigned int i = 0; i < handles.size(); ++i) {
              glVertex2f(movedHandles[i].x(), movedHandles[i].y());
           }
           glEnd();
       }

   }
}

Grid::~Grid(void)
{
   if (!points) return;

   for (int i = 0; i < size; ++i)   {
      delete [] points[i];
      delete [] barycentricCoords[i];
   }

   delete [] points;
   delete [] barycentricCoords;
   free(P_Matrix);
   free(C_Matrix);
   free(points_valid);
}

