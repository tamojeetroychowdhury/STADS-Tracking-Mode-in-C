#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mainlism.h"

int lism(double unmatched[][2], int len_unmatched, double matched[][3])
{
    long double unmatched_input[40][3];
    int i=0;
    for (i=0; i<len_unmatched; i++)
    {
        unmatched_input[i][0] = i;
        unmatched_input[i][1] = unmatched[i][0];
        unmatched_input[i][2] = unmatched[i][1];
    }

    long double body_vecs[8876][4];
    int sm_IS[8876][2];
    int lism_matched = 0;

    lism_matched = sm(unmatched_input, len_unmatched, body_vecs, sm_IS);
    //double x=0, y=0, f=FOCAL_LENGTH;
    double mag = 0;
    for (i=0; i<lism_matched; i++)
    {        
        mag = FOCAL_LENGTH/body_vecs[i][3];
        matched[i][0] = body_vecs[i][1]*mag;
        matched[i][1] = body_vecs[i][2]*mag;
        matched[i][2] = sm_IS[i][1];
    }
    return lism_matched;
}