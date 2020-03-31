/* -------------------------------------------------------------------
   File "cap3d.geo"

   This file is the geometrical description used by GMSH to produce
   the file "cap3d.msh".
   ------------------------------------------------------------------- */

/* Definition of some parameters for geometrical dimensions 
 */
zc = 10.0 ; zp = 8.0; rc = 0.1 ; rp = 4.0 ; a = 1.0;
h = 0.1;

/* Definition of gemetrical points */

Point(1)  = { 0.0  , 0.0 ,  -zc , h    } ;
Point(2)  = { rc   , 0.0 ,  -zc , h    } ;
Point(3)  = { rc   , 0.0 ,  -a  , h    } ;
Point(4)  = { rc   , 0.0 ,  0.0 , h/15 } ;
Point(5)  = { rc+a , 0.0 ,  0.0 , h    } ;
Point(6)  = { rp   , 0.0 ,  0.0 , h*10 } ;
Point(7)  = { rp   , 0.0 ,  zp  , h*10 } ;
Point(8)  = { 0.0  , 0.0 ,  zp  , h*10 } ;
Point(9)  = { 0.0  , 0.0 ,  a   , h    } ;
Point(10) = { 0.0  , 0.0 ,  0.0 , h/15 } ;
Point(11) = { 0.0  , 0.0 ,  -a  , h    } ;

/* Definition of gemetrical lines */
Line(1)  = {1,2};   Line(2)  = {2,3};  Line(3) = {3,4};
Line(4)  = {4,5};   Line(5)  = {5,6};  Line(6) = {6,7};
Line(7)  = {7,8};   Line(8)  = {8,9};  Line(9) = {9,10}; 	
Line(10) = {10,11}; Line(11) = {11,1};

/* Definition of geometrical surfaces */

Line Loop(12) = {1,2,3,4,5,6,7,8,9,10,11};   
Plane Surface(13) = {12};

Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{13}; }
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{56}; }
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{99}; }
Extrude {{0,0,1}, {0,0,0}, Pi/2} { Surface{142}; }

Physical Volume(124) = {1,2, 3, 4};
/*Physical Surface(143) = {177,180,142,};*/
