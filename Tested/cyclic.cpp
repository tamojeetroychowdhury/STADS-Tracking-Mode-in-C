#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sm_TM_consts.h"
#include "sm_TM_CP.h"
#include "sm_TM_RBM.h"
#include "sm_TM_CP_inputsort.h"
#include "sm_TM_SNT_matched.h"

#include <fstream>
#include <iostream>
using namespace std;

int main()
{
int len1, len2, len3;
cin>>len1>>len2;


double prev1[40][3];


double prev2[40][3];

double p1[40][3];
double p2[40][3];
int comm;
double pred[40][3];

double truemat[40][2];
double matchmat[40][3];

int snt[5060][11];

fstream newfile;
newfile.open("sm_TM_SNT.txt", ios::in);

int x, c=0;
while (newfile>>x)
{snt[c/11][c%11] = x;
c++;
}


fstream file;
file.open("Guide_Star_5060.txt", ios::in);

double gd_sc[5060][4];
c=0;
double y;
while (file>>y)
{gd_sc[c/4][c%4] = y;
c++;
}

int good = 0;
int num = 18;
double r = num*pixel_size;

double e = 0.000001;
double snt_out[5][3], snt_match=0;

for (int h=0; h<len1; h++) cin>>prev1[h][2]>>prev1[h][0]>>prev1[h][1];
for (int h=0; h<len2; h++) cin>>prev2[h][2]>>prev2[h][0]>>prev2[h][1];


while (true)
{
comm = sm_TM_CP_inputsort(prev1, len1, prev2, len2, p1, p2);

//cout<<"done4"<<endl<<endl;

sm_TM_CP(p1, p2, FOCAL_LENGTH, pred, comm, 1);
//cout<<"done5"<<endl<<endl;

cin>>len3;


for (int h=0; h<len3; h++) cin>>truemat[h][0]>>truemat[h][1];

good = sm_TM_RBM(pred, comm, truemat, len3, r, true, true, matchmat);

cout<<endl<<endl;
for (int i=0; i<good; i++)
    cout<<matchmat[i][2]<<' '<<matchmat[i][0]<<' '<<matchmat[i][1]<<endl;

cout<<endl<<endl<<good<<" out of "<<comm<<endl<<endl;


snt_match = sm_TM_SNT_match(matchmat, good, truemat, len3, snt, gd_sc, FOCAL_LENGTH, 1, 20, e, snt_out);

for (int i=0; i<snt_match; i++)
    cout<<snt_out[i][0]<<' '<<snt_out[i][1]<<' '<<snt_out[i][2]<<endl;

cout<<endl<<endl<<snt_match<<endl;

for (int i=0; i<len2; i++)
{
prev1[i][0] = prev2[i][0];
prev1[i][1] = prev2[i][1];
prev1[i][2] = prev2[i][2];}
//cout<<"done1"<<endl;
len1 = len2;

//for (int i=0; i<len1; i++) cout<<prev1[i][2]<<' '<<prev1[i][0]<<' '<<prev1[i][1]<<endl;
//cout<<endl<<endl;

for (int i=0; i<good; i++)
{
prev2[i][0] = matchmat[i][0];
prev2[i][1] = matchmat[i][1];
prev2[i][2] = matchmat[i][2];
}
//cout<<"done2"<<endl;

for (int i = 0; i<snt_match; i++)
{
prev2[i+good][0] = snt_out[i][0];
prev2[i+good][1] = snt_out[i][1];
prev2[i+good][2] = snt_out[i][2];
}
len2 = good + snt_match;
//cout<<"done3"<<endl;

//for (int i=0; i<len2; i++) cout<<prev2[i][2]<<' '<<prev2[i][0]<<' '<<prev2[i][1]<<endl;
//cout<<endl<<endl;

}
}
