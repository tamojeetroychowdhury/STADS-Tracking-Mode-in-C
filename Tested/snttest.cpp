#include <fstream>

#include <iostream>
using namespace std;

int main(){
int snt[5060][11];

fstream newfile;
newfile.open("sm_TM_SNT.txt", ios::in);

char tp[100];
int x, c=0;
while (newfile>>x)
{snt[c/11][c%11] = x;
c++;
}


}

