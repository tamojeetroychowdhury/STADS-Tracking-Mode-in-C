#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "functions.h"
//#include "cortos.h"
#include "input1.h"

//double arr_truemat[2000][2];
void main()
{
int i = 0;


int z=0;
for (z=0; z<len3; z++)
{
truemat[z][0] = z;
truemat[z][1] = arr_truemat[z][0];
truemat[z][2] = arr_truemat[z][1];
}


int match = 0;
int tot_len3 = len3;


double new_matched[40][3];
match = tracking(prev1, len1, prev2, len2, truemat, len3, new_matched);

printf("%d\n", match);
int run_lism = 0, run_tm = 0;

for (z=1; z < 5; z++)
{
    len3 = arr_len3[z];
    for (i=0; i<len2; i++)
    {
        prev1[i][0] = prev2[i][0];
        prev1[i][1] = prev2[i][1];
        prev1[i][2] = prev2[i][2];
    }
    len1 = len2;

    for (i=0; i<match; i++)
    {
        prev2[i][0] = new_matched[i][0];
        prev2[i][1] = new_matched[i][1];
        prev2[i][2] = new_matched[i][2];
    }

    len2 = match;

    for (i=0; i<len3; i++)
    {
        truemat[i][0] = i;
        truemat[i][1] = arr_truemat[tot_len3 + i][0];
        truemat[i][2] = arr_truemat[tot_len3 + i][1];
    }

    tot_len3 = tot_len3 + len3;

    match = tracking(prev1, len1, prev2, len2, truemat, len3, new_matched);
    /*if (match < 4)
    {
        match = lism(truemat, len3, new_matched);
        run_lism++;
    }
    else run_tm++;*/

}

printf("\n\n%d\n%d\n\n",run_lism,run_tm);
//cortos_exit(0);
}