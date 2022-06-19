#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
/*void give_angles(double u1, double v1, double f,
                 double u2, double v2,
                 double& phi, double& theta, double& psi)
{
double a1, b1, c1, d1, e1, f1;
double a2, b2, c2, d2, e2, f2;

a1 = u1*v1/f;
b1 = f + v1*v1/f;
c1 = -f - u1*u1/f;
d1 = -u1*v1/f;
e1 = v1;
f1 = -u1;

a2 = u2*v2/f;
b2 = f + v2*v2/f;
c2 = -f - u2*u2/f;
d2 = -u2*v2/f;
e2 = v2;
f2 = -u2;

double m[3][3];
m[0][0] = a1*a1 + b1*b1 + a2*a2 + b2*b2;
m[0][1] = m[1][0] = a1*c1 + b1*d1 + a2*c2 + b2*d2;
m[0][2] = m[2][0] = a1*e1 + b1*f1 + a2*e2 + b2*f2;
m[1][1] = c1*c1 + d1*d1 + c2*c2 + d2*d2;
m[1][2] = m[2][1] = c1*e1 + d1*f1 + c2*e2 + d2*f2;
m[2][2] = e1*e1 + f1*f1 + e2*e2 + f2*f2;

double determinant = 0;
double minv[3][3];

for(int i=0;i<3;i++)
    determinant = determinant + (m[0][i]*(m[1][(i+1)%3]*m[2][(i+2)%3] - m[1][(i+2)%3]*m[2][(i+1)%3]));

for(int i=0;i<3;i++)
{
    for(int j=0;j<3;j++)
        minv[i][j] = ((m[(i+1)%3][(j+1)%3] * m[(i+2)%3][(j+2)%3]) - (m[(i+1)%3][(j+2)%3]*m[(i+2)%3][(j+1)%3]))/ determinant;

}

phi = (minv[0][0]*a1 + minv[0][1]*c + minv[0][2]*e) * (u2-u1) + (minv[0][0]*b + minv[0][1]*d + minv[0][2]*f) * (v2-v1);
theta = (minv[1][0]*a1 + minv[1][1]*c + minv[1][2]*e) * (u2-u1) + (minv[1][0]*b + minv[1][1]*d + minv[1][2]*f) * (v2-v1);
psi = (minv[2][0]*a + minv[2][1]*c + minv[2][2]*e) * (u2-u1) + (minv[2][0]*b + minv[2][1]*d + minv[2][2]*f) * (v2-v1);

}
*/


//sm_TM_CP_F = focal length in cm
//pixel size = in cm
void sm_TM_CP (double sm_TM_CP_prevmat_1[][3],
               double sm_TM_CP_prevmat_2[][3],
               double sm_TM_CP_F,
               double sm_TM_CP_predmat[][3],
               int N,
               double pixel_size)
{
double phi = 0, theta = 0, psi = 0;
double u, v, f = sm_TM_CP_F;
double dun, dvn;
double m[2*N][3];

for (int k=0; k<N; k++)
{
    u = sm_TM_CP_prevmat_2[k][0] * pixel_size;
    v = sm_TM_CP_prevmat_2[k][1] * pixel_size;
    m[2*k][0] = u*v/f;
    m[2*k][1] = -f - u*u/f;
    m[2*k][2] = v;
    m[2*k+1][0] = f + v*v/f;
    m[2*k+1][1] = -u*v/f;
    m[2*k+1][2] = -u;
}

for (int k=0; k<2*N; k++)
cout<<m[k][0]/f<<' '<<m[k][1]/f<<' '<<m[k][2]/f<<endl;


double m2[3][3] = {0,0,0,0,0,0,0,0,0};

for (int k=0; k<2*N; k++)
{
   m2[0][0] += m[k][0]*m[k][0];
   m2[1][0] += m[k][1]*m[k][0];
   m2[2][0] += m[k][2]*m[k][0];
   m2[0][1] += m[k][0]*m[k][1];
   m2[1][1] += m[k][1]*m[k][1];
   m2[2][1] += m[k][2]*m[k][1];
   m2[0][2] += m[k][0]*m[k][2];
   m2[1][2] += m[k][1]*m[k][2];
   m2[2][2] += m[k][2]*m[k][2];
}

for (int k=0; k<3; k++)
cout<<m2[k][0]/f/f<<' '<<m2[k][1]/f/f<<' '<<m2[k][2]/f/f<<endl;
cout<<endl<<endl;

double determinant = 0;
double minv[3][3];

for(int i=0;i<3;i++)
    determinant = determinant + (m2[0][i]*(m2[1][(i+1)%3]*m2[2][(i+2)%3] - m2[1][(i+2)%3]*m2[2][(i+1)%3]));

for(int i=0;i<3;i++)
{
    for(int j=0;j<3;j++)
        minv[i][j] = ((m2[(i+1)%3][(j+1)%3] * m2[(i+2)%3][(j+2)%3]) - (m2[(i+1)%3][(j+2)%3]*m2[(i+2)%3][(j+1)%3]))/ determinant;
}

for (int k=0; k<3; k++)
cout<<minv[k][0]*f*f<<' '<<minv[k][1]*f*f<<' '<<minv[k][2]*f*f<<endl;

cout<<endl<<endl;

for (int k=0; k<2*N; k++)
cout<<m[k][0]/f<<' '<<m[k][1]/f<<' '<<m[k][2]/f<<endl;

cout<<endl<<endl;

//WRONG OUTPUT IS FROM THE FOLLOWING FOR LOOP ----------------------------------------------
double m3[3][2*N];
for (int k=0; k<2*N; k++)
{
    m3[0][k] = minv[0][0]*m[k][0] + minv[0][1]*m[k][1] + minv[0][2]*m[k][2];
    m3[1][k] = minv[1][0]*m[k][0] + minv[1][1]*m[k][1] + minv[1][2]*m[k][2];
    m3[2][k] = minv[2][0]*m[k][0] + minv[2][1]*m[k][1] + minv[2][2]*m[k][2];
}

for (int k=0; k<2*N; k++)
cout<<m3[k][0]*f<<' '<<m3[k][1]*f<<' '<<m3[k][2]*f<<endl;

cout<<endl<<endl;

double col[2*N];
for (int k=0; k<N; k++)
{
    col[2*k] = sm_TM_CP_prevmat_2[k][0] - sm_TM_CP_prevmat_1[k][0];
    col[2*k+1] = sm_TM_CP_prevmat_2[k][1] - sm_TM_CP_prevmat_1[k][1];
}


for (int k=0; k<2*N; k++)
{
    phi += m3[0][k]*col[k];
    theta += m3[1][k]*col[k];
    psi += m3[2][k]*col[k];
}

for (int k=0; k<N; k++)
{
    u = sm_TM_CP_prevmat_2[k][0] * pixel_size;
    v = sm_TM_CP_prevmat_2[k][1] * pixel_size;

    dun = (u*v/f)*phi + (-f - u*u/f)*theta + v*psi;
    dvn = (f + v*v/f)*phi + (-u*v/f)*theta + (-u)*psi;

    sm_TM_CP_predmat[k][0] = sm_TM_CP_prevmat_2[k][0] + dun/pixel_size;
    sm_TM_CP_predmat[k][1] = sm_TM_CP_prevmat_2[k][1] + dvn/pixel_size;
    //sm_TM_CP_predmat[k][0] = dun/pixel_size;
    //sm_TM_CP_predmat[k][1] = dvn/pixel_size;

    sm_TM_CP_predmat[k][2] = sm_TM_CP_prevmat_1[k][2];
}


}


