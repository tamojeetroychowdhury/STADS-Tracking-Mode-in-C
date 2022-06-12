#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void give_angles(double u, double v, double f,
                 double du, double dv,
                 double& phi, double& theta, double& psi)
{
double a, b, c, d, e, g;
a = u*v/f;
b = f + v*v/f;
c = -f - u*u/f;
d = -u*v/f;
e = v;
f = -u;
double m[3][3];
m[0][0] = a*a + b*b;
m[0][1] = m[1][0] = a*c + b*d;
m[0][2] = m[2][0] = a*e + b*f;
m[1][1] = c*c + d*d;
m[1][2] = m[2][1] = c*e + d*f;
m[2][2] = e*e + f*f;

double determinant;
double minv[3][3];

for(int i=0;i<3;i++)
    determinant = determinant + (m[0][i]*(m[1][(i+1)%3]*m[2][(i+2)%3] - m[1][(i+2)%3]*m[2][(i+1)%3]));

for(int i=0;i<3;i++)
{
    for(int j=0;j<3;j++)
        minv[i][j] = ((m[(i+1)%3][(j+1)%3] * m[(i+2)%3][(j+2)%3]) - (m[(i+1)%3][(j+2)%3]*m[(i+2)%3][(j+1)%3]))/ determinant;

}

phi = (minv[0][0]*a + minv[0][1]*c + minv[0][2]*e)*du + (minv[0][0]*b + minv[0][1]*d + minv[0][2]*f)*dv;
theta = (minv[1][0]*a + minv[1][1]*c + minv[1][2]*e)*du + (minv[1][0]*b + minv[1][1]*d + minv[1][2]*f)*dv;
psi = (minv[2][0]*a + minv[2][1]*c + minv[2][2]*e)*du + (minv[2][0]*b + minv[2][1]*d + minv[2][2]*f)*dv;

}



//sm_TM_CP_F = focal length in cm
//pixel size = in cm
void sm_TM_CP (double sm_TM_CP_prevmat_1[][3],
               double sm_TM_CP_prevmat_2[][3],
               double sm_TM_CP_F,
               double sm_TM_CP_predmat[][3],
               int N,
               double pixel_size)
{
double phi, theta, psi;
double u, v, du, dv, f = sm_TM_CP_F;
double dun, dvn;

for (int k=0; k<N; k++)
{
    u = sm_TM_CP_prevmat_1[k][0] * pixel_size;
    v = sm_TM_CP_prevmat_1[k][1] * pixel_size;
    du = (sm_TM_CP_prevmat_2[k][0] - sm_TM_CP_prevmat_1[k][0]) * pixel_size;
    dv = (sm_TM_CP_prevmat_2[k][1] - sm_TM_CP_prevmat_1[k][1]) * pixel_size;
    give_angles(u,v,f,du,dv,phi,theta,psi);

    dun = (u*v/f)*phi + (-f - u*u/f)*theta + v*psi;
    dvn = (f + v*v/f)*phi + (-u*v/f)*theta + (-u)*psi;

    sm_TM_CP_predmat[k][0] = sm_TM_CP_prevmat_2[k][0] + dun/pixel_size;
    sm_TM_CP_predmat[k][1] = sm_TM_CP_prevmat_2[k][1] + dvn/pixel_size;
    sm_TM_CP_predmat[k][2] = sm_TM_CP_prevmat_1[k][2];
}


}
