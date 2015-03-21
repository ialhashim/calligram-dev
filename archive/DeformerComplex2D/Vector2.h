#ifndef VECTOR2_CM
#define VECTOR2_CM

#define NOMINMAX
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>


class Vector2
{

   static const int vDim = 2;

private:
   double vec[vDim];

   inline void copy (const Vector2 & o){
      for (int i = 0; i < vDim; i++)
         vec[i] = o.vec[i];
   }

   /* NOTE: REMOVING THIS WILL BREAK THE DISTANCE_TO_SEGMENT CODE
   inline const Vector2& operator*=(const Vector2& rhs){
      for (int i = 0; i < vDim; i++)
         vec[i] = vec[i]*rhs.vec[i];		
      return *this;
   }
   */

public:

   static Vector2 e1() {return Vector2(1, 0);}
   static Vector2 e2() {return Vector2(0, 1);}
   static Vector2 ones() {return Vector2(1, 1);}
   static Vector2 zero() {return Vector2(0, 0);}

   Vector2(void){
      for (int i = 0; i < vDim; i++)
         vec[i] = 0;
   }

   Vector2(double x, double y){
      vec[0] = x; vec[1] = y; 
   }

   Vector2(const Vector2 & o){
      copy(o);
   }

   // Comparison
   inline bool operator==(const Vector2& rhs) const{
      return (fabs(vec[0] - rhs.vec[0]) < 1.e-8 && 
         fabs(vec[1] - rhs.vec[1]) <  1.e-8);
   }

   inline bool operator!=(const Vector2& rhs) const{
      return (!(*this == rhs));
   }
   

   // Assignment
   inline const Vector2& operator=(const double rhs){
      vec[0] = rhs;
      vec[1] = 0;
      return *this;
   }

   inline const Vector2& operator=(const Vector2& rhs){
      if (&rhs != this)
         copy(rhs);
      return *this;
   }

   inline const Vector2& operator+=(const Vector2& rhs){
      for (int i = 0; i < vDim; i++)
         vec[i] = vec[i]+rhs.vec[i];		
      return *this;
   }

   inline const Vector2& operator-=(const Vector2& rhs){
      for (int i = 0; i < vDim; i++)
         vec[i] = vec[i]-rhs.vec[i];		
      return *this;
   }

   inline const Vector2& operator*=(const Vector2& rhs){
      double a = vec[0];
      double b = vec[1];
      double c = rhs.vec[0];
      double d = rhs.vec[1];

      // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
      vec[0] = a * c - b * d;
      vec[1] = a * d + b * c;

      return *this;
   }

   inline const Vector2& operator*=(const double& rhs){
      for (int i = 0; i < vDim; i++)
         vec[i] = vec[i] * rhs;		
      return *this;
   }

   inline const Vector2& operator/=(const Vector2& rhs){
      double a = vec[0];
      double b = vec[1];
      double c = rhs.vec[0];
      double d = rhs.vec[1];

      double denominator = c * c + d * d;

      vec[0] = (a * c + b * d) / denominator;
      vec[1] = (b * c - a * d) / denominator;

      return *this;
   }

   inline const Vector2& operator/=(const double& rhs){
      for (int i = 0; i < vDim; i++)
         vec[i] = vec[i]/rhs;		
      return *this;
   }

   // Arithmetic
   inline Vector2 operator+(const Vector2& rhs) const {
      Vector2 tmp(*this);
      tmp += rhs;
      return tmp;
   }


   inline Vector2 operator-(const Vector2& rhs) const {
      Vector2 tmp(*this);
      tmp -= rhs;
      return tmp;
   }	

   double dot(const Vector2& rhs) const {
      double sum = 0;
      for (int i = 0; i < vDim; i++)
         sum += rhs.vec[i]*vec[i];
      return sum;
   }

   double dotItself() const {
      double sum = 0;
      for (int i = 0; i < vDim; i++)
         sum += vec[i]*vec[i];
      return sum;
   }

   inline Vector2 operator*(const double& rhs) const {
      Vector2 tmp(*this);
      tmp *= rhs;
      return tmp;
   }
   
   inline Vector2 operator*(const Vector2& rhs) const {
      Vector2 tmp(*this);
      tmp *= rhs;
      return tmp;
   }

   inline Vector2 operator/(const double& rhs) const {
      Vector2 tmp(*this);
      tmp /= rhs;
      return tmp;
   }

   inline Vector2 operator/(const Vector2& rhs) const {
      Vector2 tmp(*this);
      tmp /= rhs;
      return tmp;
   }

   //Acess
   inline double & operator[] (int i){
      //if (i >= vDim || i < 0 ){
      //	std::cerr << "Out of bounds vector accessing: [" << i << "]" << std::endl;
      //	return vec[0];
      //}
      //else
      return vec[i];
   }

   double crossProductZ(const Vector2& rhs) const{
      return vec[0]*rhs[1] - vec[1]*rhs[0];
   }

   inline double operator[] (int i) const{
      //if (i >= vDim || i < 0){
      //	std::cerr << "Out of bounds vector accessing: [" << i << "]" << std::endl;
      //	return 0;
      //}
      //else
      return vec[i];
   }

   double x() const {return vec[0];}
   double & x() {return vec[0];}
   double y() const {return vec[1];}
   double & y() {return vec[1];}

   // Other
   double abs () const {
      double sum = 0;
      for (int i = 0; i < vDim; i++)
         sum += vec[i]*vec[i];
      return sqrt(sum);
   }

   // Compute complex logarithm
   Vector2 c_log() const   {
      double r = hypot(vec[0], vec[1]);
      double theta = atan2(vec[1], vec[0]);

      return Vector2(log(r), theta);
   }

   Vector2 conjugate() const    {
      return Vector2(vec[0], -vec[1]);
   }

   Vector2 norm() const {
      Vector2 norm;
      double length = abs();
      for (int i = 0; i < vDim; i++)
         norm[i] = vec[i]/length;
      return norm;
   }

   void normalize(){
      double length = abs();
      for (int i = 0; i < vDim; i++)
         vec[i] = vec[i]/length;
   }

   void bbox( Vector2& min, Vector2& max ){
      if( vec[0] < min.vec[0] ) min.vec[0] = vec[0];
      if( vec[0] > max.vec[0] ) max.vec[0] = vec[0];
      if( vec[1] < min.vec[1] ) min.vec[1] = vec[1];
      if( vec[1] > max.vec[1] ) max.vec[1] = vec[1];
   }

   void negate () {
      vec[0] = -vec[0];
      vec[1] = -vec[1];
   }

   // Input/Output
   void print(std::ostream& out, char * name) const{ 
      for (int i = 0; i < vDim; i++)
         out << name << "(" << (i+1) << ")= " << vec[i] << "; ";
      out << std::endl;
   } 

   friend std::ostream& operator<<(std::ostream& out,const Vector2 & m) { 
      out << m.vec[0];
      for (int i = 1; i < vDim; i++)
         out << "\t" << m.vec[i];
      return out; 
   } 


   friend std::istream& operator>>(std::istream& in, Vector2 & m) { 
      in >> std::ws;
      for (int i = 0; i < vDim; i++)
         in >> m.vec[i];
      return in; 
   } 


   friend Vector2 operator*(const double& rhs, const Vector2 & v){
      Vector2 res;
	  for (int i = 0; i < v.vDim; i++)
         res[i] = v[i]*rhs;		
      return res;
   }

   double triangleArea(const Vector2 & s1, const Vector2 & s2) {
      return 0.5*fabs((s1 - *this).crossProductZ(s2 - *this));
   }

   double distanceToSegment(const Vector2 & s1, const Vector2 & s2, Vector2 & proj, double & centerDist) {
      double dist;
      double edgeLen = (s1-s2).abs();
      dist = 2*triangleArea(s1, s2)/edgeLen;
      Vector2 nor(-s1.y() + s2.y(), s1.x() - s2.x());
      nor.normalize();
      nor *= dist;
      proj = *this;
      if (nor.dot(s1 - *this))
         proj -= nor;
      else
         proj += nor;
      centerDist = ((s1 + s2)*0.5 - proj).dotItself();

      if ((*this - s1).dot(s2 - s1) < 0 || (*this - s2).dot(s1 - s2) < 0) {
         dist = (s1 - *this).abs();
         double dist1 = (s2 - *this).abs();
         proj = s1;
         if (dist1 < dist) {
            dist = dist1;
            proj = s2;
         }
      }
      return dist;
   }


   bool isInVector(std::vector<Vector2> vecBig) {
      for (unsigned int i = 0; i < vecBig.size(); i++)
         //if (fabs(vecBig[i].x() - x()) < 1.e-6 && fabs(vecBig[i].y() - y()) < 1.e-6)
         if (vecBig[i] == *this)
            return true;
      return false;
   }

   // isLeft(): tests if a point is Left|On|Right of an infinite line.
   //    Return: >0 for P left of the line through P0 and P1
   //            =0 for P on the line < 0 for P2 right of the line
   bool isLeft(Vector2 P0, Vector2 P1) {
      double check = (P1.x() - P0.x())*(y() - P0.y()) - (x() - P0.x())*(P1.y ()- P0.y());
      if (check > 0)
         return true;
      //if (check > -1.e-18 && ((P0-P1).dotItself() < (P0-*this).dotItself()))
      //   return true;
      return false;

   }

   void rotate(const double & angleClockwise) {
      Vector2 tmp = *this;
      vec[0] = cos(angleClockwise)*tmp.x() - sin(angleClockwise)*tmp.y();
      vec[1] = sin(angleClockwise)*tmp.x() + cos(angleClockwise)*tmp.y();
   }

   static bool segmentIntersect(const Vector2 & x1, const Vector2 & x2, const Vector2 & x3, const Vector2 & x4,
      double &x, double &y) {
         bool intersect = true;
         double co = (x4.y() - x3.y())*(x2.x() - x1.x()) - (x4.x()-x3.x())*(x2.y()-x1.y());

         if (fabs(co) < 10e-12)
            intersect = false;
         else {
            double ua = ((x4.x() - x3.x())*(x1.y() - x3.y()) - (x4.y()-x3.y())*(x1.x()-x3.x()))/co;
            if (ua < 0 || ua > 1) 
               intersect = false;
            else {
               double ub = ((x2.x() - x1.x())*(x1.y() - x3.y()) - (x2.y()-x1.y())*(x1.x()-x3.x()))/co;
               if (ub < 0 || ub > 1)
                  intersect = false;
               else {
                  x = x1.x() + ua*(x2.x() - x1.x());
                  y = x1.y() + ua*(x2.y() - x1.y());
                  //cout << "Show[Graphics[{{PointSize[0.02], Point[{" << x1.x() << ", " << x1.y() << "}]}, {PointSize[0.02], Point[{" 
                  //   <<  x2.x() << "," << x2.y() << "}], Point[{"<< x3.x() << ", " << x3.y() << "}]}, {PointSize[0.02], Point[{"
                  //   << x4.x() << ", " << x4.y() << "}]}, {PointSize[0.02], RGBColor[1, 0, 0], Point[{"
                  //   << x << ", " << y << "}]}, Line[{{" << x1.x() << ", " << x1.y() << "}, {" << x2.x() << ", " << x2.y() << "}}], Line[{{"
                  //   << x3.x() << ", " << x3.y() << "}, {" << x4.x() << ", " << x4.y() << "}}]}]]" << endl;
               }
            }
         }
         return intersect;
   }


};



#endif

