#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sm_TM_consts.h"
#include "sm_TM_CP.h"
#include "sm_TM_RBM.h"
#include "sm_TM_CP_inputsort.h"
#include "sm_TM_SNT_matched.h"
#include "sm_GC.h"
#include "sm_TM_SNT.h"

#include "input1.h"
//using namespace std;

int main()
{
int tot_len3 = 0;

for (int z=0; z<7; z++)
{
comm = sm_TM_CP_inputsort(prev1, len1, prev2, len2, p1, p2);

sm_TM_CP(p1, p2, FOCAL_LENGTH, pred, comm, 1);

len3 = arr_len3[z];

for (int h=0; h<len3; h++) {truemat[h][0] = arr_truemat[h+tot_len3][0];  truemat[h][1] = arr_truemat[h+tot_len3][1];}
tot_len3 = tot_len3 + len3;


good = sm_TM_RBM(pred, comm, truemat, len3, r, true, true, matchmat);

printf("\n");
printf("\n");

for (int i=0; i<good; i++)
{
    printf("%f ", matchmat[i][2]);
    printf("%f ", matchmat[i][0]);
    printf("%f\n", matchmat[i][1]);
}

printf("\n\n");
printf("%d out of ", good);
printf("%d", comm);
printf("\n\n");


snt_match = sm_TM_SNT_match(matchmat, good, truemat, len3, snt, sm_GC, FOCAL_LENGTH, 1, 20, e, snt_out);

for (int i=0; i<snt_match; i++)
    {
    printf("%f ", snt_out[i][2]);
    printf("%f ", snt_out[i][0]);
    printf("%f\n", snt_out[i][1]);
    }


    //cout<<snt_out[i][0]<<' '<<snt_out[i][1]<<' '<<snt_out[i][2]<<endl;

printf("\n\n");
printf("%f", snt_match);
printf("\n\n");

for (int i=0; i<len2; i++)
{
prev1[i][0] = prev2[i][0];
prev1[i][1] = prev2[i][1];
prev1[i][2] = prev2[i][2];}
len1 = len2;


for (int i=0; i<good; i++)
{
prev2[i][0] = matchmat[i][0];
prev2[i][1] = matchmat[i][1];
prev2[i][2] = matchmat[i][2];
}

for (int i = 0; i<snt_match; i++)
{
prev2[i+good][0] = snt_out[i][0];
prev2[i+good][1] = snt_out[i][1];
prev2[i+good][2] = snt_out[i][2];
}
len2 = good + snt_match;

}
}
