/************************************************************************/
/**

   \file       pldist.c
   
   \version    V1.1
   \date       07.07.14
   \brief      Calculate distance from a point to a line
   
   \copyright  (c) University of Reading / Dr. Andrew C. R. Martin 1999-2014
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College London,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

-  V1.1  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Geometry
   #FUNCTION  blPointLineDistance()
   Calculates the shortest distance from a point P to a line between
   points P1 and P2. This value is returned and the point on the line
   can also be returned. See also blDistPtLine()

*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
/*#include "MathType.h"*/

/************************************************************************/
/* Defines and macros
*/
#define CROSSPRODUCT(p1,p2,p3) \
        (p3)[0] = (p1)[1]*(p2)[2] - (p1)[2]*(p2)[1]; \
        (p3)[1] = (p1)[2]*(p2)[0] - (p1)[0]*(p2)[2]; \
        (p3)[2] = (p1)[0]*(p2)[1] - (p1)[1]*(p2)[0]
#define DOTPRODUCT(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2])

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>float blPointLineDistance(float Px, float Py, float Pz,
                          float P1x, float P1y, float P1z,
                          float P2x, float P2y, float P2z,
                          float *Rx, float *Ry, float *Rz,
                          float *frac)
   ------------------------------------------------------
*//**
   \param[in]     Px          Point x coordinate
   \param[in]     Py          Point y coordinate
   \param[in]     Pz          Point z coordinate
   \param[in]     P1x         Line start x coordinate
   \param[in]     P1y         Line start y coordinate
   \param[in]     P1z         Line start z coordinate
   \param[in]     P2x         Line end x coordinate
   \param[in]     P2y         Line end y coordinate
   \param[in]     P2z         Line end z coordinate
   \param[out]    *Rx         Nearest point on line x coordinate
   \param[out]    *Ry         Nearest point on line y coordinate
   \param[out]    *Rz         Nearest point on line z coordinate
   \param[out]    *frac       Fraction along P1-P2 of R
   \return                       Distance from P to R

   Calculates the shortest distance from a point P to a line between
   points P1 and P2. This value is returned.

   If the Rx,Ry,Rz pointers are all non-NULL, then the point on the
   line nearest to P is output.

   If the frac pointer is non-NULL then the fraction of R along the
   P1-P2 vector is output. Thus:
         R==P1 ==> frac=0
         R==P2 ==> frac=1
   Thus if (0<=frac<=1) then the point R is within the line segment
   P1-P2

-  16.11.99 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
float blPointLineDistance(float Px, float Py, float Pz,
                         float P1x, float P1y, float P1z,
                         float P2x, float P2y, float P2z,
                         float *Rx, float *Ry, float *Rz,
                         float *frac)
{
   
   float (A)[3];
   float (u)[3];
   float (Q)[3];
   float (PQ)[3];
   float (PR)[3];
   float (QP)[3];
   float (QP2)[3];
   float  alen, len, f;
   
   
   /* Calculate vector from P1 to P2                                    */
   A[0] = P2x - P1x;
   A[1] = P2y - P1y;
   A[2] = P2z - P1z;
   print(A)
   /* Calculate length of this vector                                   */
   alen = sqrt(DOTPRODUCT(A,A));

   /* If the two ends of the line are coincident then return the distance
      from either of them
   */
   if(alen==(float)0.0)
   {
      len = sqrt((Px-P1x)*(Px-P1x) + 
                 (Py-P1y)*(Py-P1y) + 
                 (Pz-P1z)*(Pz-P1z));
      if(frac!=NULL)
         *frac = 0.0;
      if(Rx != NULL && Ry != NULL && Rz != NULL)
      {
         *Rx = P1x;
         *Ry = P1y;
         *Rz = P1z;
      }

      return(len);
   }

   /* Calculate the unit vector along A                                 */
   u[0] = A[0] / alen;
   u[1] = A[1] / alen;
   u[2] = A[2] / alen;
   
   /* Select Q as any point on A, we'll make it P1                      */
   Q[0] = P1x;
   Q[1] = P1y;
   Q[2] = P1z;
   
   /* Calculate vector PQ                                               */
   PQ[0] = Q[0] - Px;
   PQ[1] = Q[1] - Py;
   PQ[2] = Q[2] - Pz;

   /* Vector PR is the cross product of PQ and the unit vector
      along A (i.e. u)
   */
   CROSSPRODUCT(PQ, u, PR);
   
   /* And the length of that vector is the length we want               */
   len = sqrt(DOTPRODUCT(PR,PR));


   if(frac != NULL || Rx != NULL || Ry != NULL || Rz != NULL)
   {
      /*** OK we now know how far the point is from the line, so we   ***
       *** now want to calculate where the closest point (R) on the   ***
       *** line is to point P                                         ***/

      /* Find the projection of QP onto QP2                             */
      QP[0] = Px - Q[0];
      QP[1] = Py - Q[1];
      QP[2] = Pz - Q[2];
      
      QP2[0] = P2x - Q[0];
      QP2[1] = P2y - Q[1];
      QP2[2] = P2z - Q[2];
      
      f = DOTPRODUCT(QP, QP2) / sqrt(DOTPRODUCT(QP2, QP2));
      if(frac != NULL)
      {
         *frac = f/alen;
      }
      
      /* Find point R: this is the fraction f of the unit vector along 
         P1-P2 added onto Q 
      */
      if(Rx != NULL && Ry != NULL && Rz != NULL)
      {
         *Rx = Q[0] + f * u[0];
         *Ry = Q[1] + f * u[1];
         *Rz = Q[2] + f * u[2];
      }
   }
   
   return(len);
}


/************************************************************************/
#ifdef DEMO
int main(int argc, char **argv)
{
   float Px, Py, Pz,
        P1x, P1y, P1z,
        P2x, P2y, P2z,
        Rx, Ry, Rz, d, f;
 
   Px = 5;   Py = 2;   Pz = 0;
   P1x = 5;  P1y = 2;  P1z = 0;
   P2x = 10; P2y = 2;  P2z = 0;
   
   d = blPointLineDistance(Px, Py, Pz,
                           P1x, P1y, P1z,
                           P2x, P2y, P2z,
                           &Rx, &Ry, &Rz, &f);
   
   printf("*** Distance is %f; Point R is %f %f %f; f is %f\n",
          d,Rx,Ry,Rz,f);
   
   return(0);
}
#endif
