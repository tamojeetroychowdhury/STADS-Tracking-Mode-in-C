#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sm_TM_consts.h"
#include "sm_TM_CP.h"
#include "sm_TM_RBM.h" //Remove this to test only CP

#include <iostream>
using namespace std;

int main()
{//return 0;}
int N = 8;
double xi = -129.9038106, xf = -67.4012702;
double prev1[8][3] = {0 , 0 , 1005 ,
-15 , 28 , 3854 ,
11 , 30 , 3348 ,
28 , 35 , 1207 ,
77 , 28 , 1746 ,
-28 , 31 , 1882 ,
-15 , 32 , 2140 ,
-33 , 33 , 126
}; //COPY-PASTE FROM TESTCASE.PY


double prev2[8][3] = {0 , 0 , 1005 ,
-16 , 30 , 3854 ,
11 , 30 , 3348 ,
28 , 34 , 1207 ,
67 , 24 , 1746 ,
-30 , 34 , 1882 ,
-16 , 34 , 2140 ,
-36 , 36 , 126
}; //COPY-PASTE FROM TESTCASE.PY


double pred[8][3];
sm_TM_CP(prev1, prev2, FOCAL_LENGTH, pred, N, pixel_size);

//TO CHECK CP OUTPUT--------------------
/*for (int i=0; i<N; i++)
    cout<<pred[i][0]<<' '<<pred[i][1]<<' '<<pred[i][2]<<endl;
*/


//TO CHECK CP + RBM OUTPUT---------------
double truemat[8][2] = {0 , 0 ,
-17 , 32 ,
11 , 30 ,
27 , 33 ,
60 , 21 ,
-33 , 37 ,
-17 , 36 ,
-40 , 40
}; //COPY-PASTE FROM TESTCASE.PY
double matchmat[8][5];
sm_TM_RBM(pred, N, truemat, N, 15, true, true, matchmat);

for (int i=0; i<N; i++)
    cout<<matchmat[i][0]<<' '<<matchmat[i][1]<<' '<<matchmat[i][2]<<' '<<matchmat[i][3]<<' '<<matchmat[i][4]<<' '<<endl;


}
