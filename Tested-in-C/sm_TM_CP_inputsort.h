#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "bubblesort.h"

int sm_TM_CP_inputsort(double prevmat_1_raw[][3],
                       int N1,
                       double prevmat_2_raw[][3],
                       int N2,
                       double prevmat_1[][3],
                       double prevmat_2[][3])
{
bubbleSort3(prevmat_1_raw, N1, 2);
bubbleSort3(prevmat_2_raw, N2, 2);

int i=0, j=0, k=0;
//cout<<N1<<' '<<N2<<endl;
while (i<N1 and j<N2)
{
    //cout<<"fu"<<endl;
    //while (j<N2)
        {
            //cout<<"f"<<endl;
            if (prevmat_1_raw[i][2] == prevmat_2_raw[j][2])
            {
               prevmat_1[k][0] = prevmat_1_raw[i][0];
               prevmat_1[k][1] = prevmat_1_raw[i][1];
               prevmat_1[k][2] = prevmat_1_raw[i][2];
               prevmat_2[k][0] = prevmat_2_raw[j][0];
               prevmat_2[k][1] = prevmat_2_raw[j][1];
               prevmat_2[k][2] = prevmat_2_raw[j][2];
               i++; j++; k++;
               //cout<<"success"<<endl;
            }

            else if (prevmat_1_raw[i][2] < prevmat_2_raw[j][2]) i++; //cout<<"increment"<<endl;}
            else if (prevmat_1_raw[i][2] > prevmat_2_raw[j][2]) j++; //cout<<"increment"<<endl;}
        }
}
return k;
}
