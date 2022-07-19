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
{//return 0;}
int N = 23;

double prev1[26][3] = {-0.00122999501482241 , 0.00125580038164328 , 31 ,
-0.000531993442872951 , 0.000765151914021973 , 32 ,
-0.000862182080487385 , 0.000363145420354094 , 409 ,
-0.000531892252678178 , 0.000764741735390563 , 633 ,
-0.0014369987045634 , -0.00102570490443067 , 1005 ,
-0.000500556043443351 , 0.00127614111994934 , 1500 ,
-0.00144991031830279 , -0.00138917730500854 , 1746 ,
-0.00146615557960649 , -0.0013720442572969 , 1882 ,
-0.00139474157640222 , -0.00103822965881446 , 2140 ,
-0.00179258251323874 , 0.00102973128832262 , 2377 ,
-0.00165210648780167 , 0.00129941895300686 , 2436 ,
0.00171846758372918 , -0.00100012236733451 , 4345 ,
-0.00157318199900453 , 0.00109874771533396 , 4530 ,
-0.000761855403927383 , -0.000246950375642829 , 4915 ,
-0.00104419340708097 , -0.00102332907952122 , 6155 ,
-0.000572864684356113 , 0.000211476328867852 , 6310 ,
-0.00156448350396059 , 0.00137219984700099 , 6542 ,
-0.00143045950403053 , -0.000765255203497532 , 6672 ,
0.00172789510202329 , -0.00110889056736703 , 6817 ,
-0.000148534567096163 , 0.000958741982493309 , 7059 ,
0.000245190735877534 , -0.000725118108939005 , 7060 ,
0.000646110080428695 , -0.000630491524486825 , 7777 ,
-0.00144078963944998 , -0.00138998943830333 , 7863 ,
-0.00137014911482598 , -3.01049373070904e-05 , 7955 ,
-0.000409660403687119 , 0.00016097055343155 , 8040 ,
-0.0014062830286159 , -0.000728109656076586 , 8397
}; //COPY-PASTE FROM TESTCASE.PY


double prev2[26][3] = {-0.00130719346127088 , 0.00133631720768881 , 31 ,
-0.000610109485485412 , 0.000844018864966206 , 32 ,
-0.000941194813208317 , 0.000442695600028771 , 409 ,
-0.000610009160581121 , 0.000843608408556263 , 633 ,
-0.00151900242942944 , -0.000944966593855277 , 1005 ,
-0.000577574915043792 , 0.00135501880850296 , 1500 ,
-0.00153267270060044 , -0.00130839354240266 , 1746 ,
-0.00154888394617148 , -0.0012912277522024 , 1882 ,
-0.00147676643273265 , -0.000957580167698553 , 2140 ,
-0.00187039710779468 , 0.00111145834346702 , 2377 ,
-0.00172932328258253 , 0.00138089790790137 , 2436 ,
0.00163628961711149 , -0.00092605125446533 , 4345 ,
-0.00165079322038177 , 0.00117999944249153 , 4530 ,
-0.000842151449495101 , -0.000167656477597185 , 4915 ,
-0.00112615415239213 , -0.00094342028807233 , 6155 ,
-0.000652168047857722 , 0.000290379543027856 , 6310 ,
-0.00164152455231371 , 0.0014535001987358 , 6542 ,
-0.00151191962098205 , -0.00068453433302011 , 6672 ,
0.00164546744354573 , -0.00103481456485777 , 6817 ,
-0.000226196783436107 , 0.00103677958402096 , 7059 ,
0.000163872424393803 , -0.000647965976656813 , 7060 ,
0.000564963300564863 , -0.000554204180504797 , 7777 ,
-0.00152355279952771 , -0.00130922467055703 , 7863 ,
-0.0014500663659824 , 5.05183806280626e-05 , 7955 ,
-0.000489059377651189 , 0.000239514083725375 , 8040 ,
-0.00148766243031322 , -0.000647440136109568 , 8397
 }; //COPY-PASTE FROM TESTCASE.PY

double p1[26][3];
double p2[26][3];
int comm;
comm = sm_TM_CP_inputsort(prev1, 26, prev2, 26, p1, p2);

double pred[62][3];
sm_TM_CP(p1, p2, FOCAL_LENGTH, pred, comm, 1);

//TO CHECK CP OUTPUT--------------------
for (int i=0; i<comm; i++)
    cout<<pred[i][0]<<' '<<pred[i][1]<<' '<<pred[i][2]<<endl;
cout<<endl<<endl;


//TO CHECK CP + RBM OUTPUT---------------
double truemat[23][2] = {-0.00138432184075037 , 0.00141696467181087 ,
-0.000688146529540794 , 0.000923008099381291 ,
-0.00102012990313077 , 0.000522365822938691 ,
-0.000688047067040722 , 0.000922597360808969 ,
-0.00160092760789223 , -0.000864119569224323 ,
-0.000654516803285583 , 0.00143402341863436 ,
-0.00161535476468288 , -0.00122750484632309 ,
-0.00163153226397643 , -0.00121030605941861 ,
-0.00155871222177141 , -0.00087682228685355 ,
0.00155422222728845 , -0.000851883044359351 ,
-0.00172833712570783 , 0.00126138219755144 ,
-0.000922365855286326 , -8.82490113973506e-05 ,
-0.00120803215545988 , -0.000863404517194805 ,
-0.000731390162504474 , 0.000369399938524261 ,
-0.00159330245232515 , -0.000603702058581278 ,
0.0015631511777276 , -0.000960642490625731 ,
-0.000303777252911798 , 0.00111493941216153 ,
8.26485018791513e-05 , -0.000570709179339672 ,
0.000483914333535938 , -0.000477812769290401 ,
-0.00160623553714822 , -0.00122835501833311 ,
-0.00152990928272043 , 0.000131260351790489 ,
-0.000568375238622674 , 0.000318173552890409 ,
-0.00156896447573098 , -0.000566658939972863
}; //COPY-PASTE FROM TESTCASE.PY
double matchmat[23][5];
int good = 0;
int num = 18;
double r = num*pixel_size;
good = sm_TM_RBM(pred, comm, truemat, 23, r, true, true, matchmat);

for (int i=0; i<good; i++)
    cout<<matchmat[i][0]<<' '<<matchmat[i][1]<<' '<<matchmat[i][2]<<' '<<matchmat[i][3]<<' '<<matchmat[i][4]<<' '<<endl;

cout<<endl<<endl<<good<<" out of "<<comm<<endl<<endl;

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



double e = 0.0000001;
double snt_out[5][3], snt_match=0;
snt_match = sm_TM_SNT_match(matchmat, good, truemat, 23, snt, gd_sc, FOCAL_LENGTH, 1, 18, e, snt_out);

for (int i=0; i<snt_match; i++)
    cout<<snt_out[i][0]<<' '<<snt_out[i][1]<<' '<<snt_out[i][2]<<endl;

cout<<endl<<endl<<snt_match<<endl;


}
