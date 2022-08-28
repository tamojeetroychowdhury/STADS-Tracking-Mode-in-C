#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "constants.h"
#include "sm_K_vec_arr.h"
#include "sm_GC.h"
#include "sm_TM_SNT.h"


//RGA functions

/*
void sort(double centroids_st[][3], unsigned long pixel_track[MAX_STARS], int tot_stars){
	int i, j, k;
    double temp_st[3];
	unsigned long temp_pt;
    i = 0;
	for (i = 0; i < tot_stars; i++){
        j = 0;
		for (j = 0; j < tot_stars-i-1; j++){
			if(pixel_track[j] < pixel_track[j + 1]){
				temp_pt = pixel_track[j];
				pixel_track[j] = pixel_track[j + 1];
				pixel_track[j + 1] = temp_pt;
                k = 0;
				for (k = 0; k < 3; k++){
					temp_st[k] = centroids_st[j][k];
					centroids_st[j][k] = centroids_st[j + 1][k];
					centroids_st[j + 1][k] = temp_st[k];
				}
			}
		}
	}
}

void getData(unsigned short p_i, unsigned short p_j, unsigned short* star_num, unsigned long x_sum[], unsigned long y_sum[], unsigned long pixel_sum[], unsigned short num_pixels[], short arr_out_img[BREADTH + 2][LENGTH + 2]){
    // base case
    if (arr_out_img[p_j][p_i] <= THRESHOLD)
        return;

    // keeping track of the sums
    x_sum[*star_num] += (unsigned long)arr_out_img[p_j][p_i] * p_i;
    y_sum[*star_num] += (unsigned long)arr_out_img[p_j][p_i] * p_j;
    pixel_sum[*star_num] += (unsigned long)arr_out_img[p_j][p_i];
    num_pixels[*star_num] += 1;
    arr_out_img[p_j][p_i] = 0;

    // recursive calls
    getData(p_i - 1, p_j, star_num, x_sum, y_sum, pixel_sum, num_pixels, arr_out_img);
    getData(p_i + 1, p_j, star_num, x_sum, y_sum, pixel_sum, num_pixels, arr_out_img);
    getData(p_i, p_j - 1, star_num, x_sum, y_sum, pixel_sum, num_pixels, arr_out_img);
    getData(p_i, p_j + 1, star_num, x_sum, y_sum, pixel_sum, num_pixels, arr_out_img);

    return;
}



void regionGrowth(short arr_out_img[BREADTH + 2][LENGTH + 2], double centroids_st[MAX_STARS][3], int* tot_stars)
{
    int valid_stars = 0;
    // clear out all arrays
    //unsigned short stars = 0;
    //for (stars = 0; stars < MAX_STARS; stars++)
    //{
        //x_sum[stars] = {0};
        //y_sum[stars] = {0};
        //pixel_sum[stars] = {0};
        //num_pixels[stars] = {0};
    //}

    //star_num = 0;
    //valid_stars = 0;
    //tot_stars = 0;

    
    unsigned short star_num = 0;
    unsigned long x_sum[MAX_STARS] = {0};
    unsigned long y_sum[MAX_STARS] = {0};
    unsigned long pixel_sum[MAX_STARS] = {0};
    unsigned short num_pixels[MAX_STARS] = {0};
    unsigned long pixel_track[MAX_STARS] = {0};


    int i = 0;
    int j = 0;
    int k = 0;

    for (j = 1; j < BREADTH + 1; j += SKIP_PIXELS)
        for (i = 1; i < LENGTH + 1; i += SKIP_PIXELS)
            if (arr_out_img[j][i] > THRESHOLD)
            {
                getData(i, j, &star_num, x_sum, y_sum, pixel_sum, num_pixels, arr_out_img);
                star_num++;
            }
    
    for (k = 0; k < star_num; k++)
    {   
        
        if ((num_pixels[k] <= STAR_MAX_PIXEL) & (num_pixels[k] > STAR_MIN_PIXEL))
        {
            valid_stars++;
            centroids_st[valid_stars-1][0] = valid_stars;
            centroids_st[valid_stars-1][1] = ((double)x_sum[k] / (double)pixel_sum[k] - ((double)(LENGTH / 2) + 0.5)) * PIXEL_WIDTH;
            centroids_st[valid_stars-1][2] = (-1 * ((double)y_sum[k] / (double)pixel_sum[k] - ((double)(BREADTH / 2) + 0.5))) * PIXEL_WIDTH;
            pixel_track[valid_stars-1] = pixel_sum[k];        
        }
    }

    *tot_stars += valid_stars;
    sort(centroids_st, pixel_track, *tot_stars);
	if (*tot_stars > NUM_MAX_STARS)
		*tot_stars = NUM_MAX_STARS;
    return;
}

//SM functions
void sm_4_star_circulate(double sm_3D_vecs[][4], int *N_circ, int N_i)
{   
    int v;
    int j;
    int k; 

    int last = 0;
    for (j = N_i - 1; j >= 0; j--)
    {
        if (sm_3D_vecs[j][0] != -1)
        {
            last = j;
            break;
        }
    }
    double curr[4] = {sm_3D_vecs[last][0], sm_3D_vecs[last][1], sm_3D_vecs[last][2], sm_3D_vecs[last][3]};
    double var;
    for (k = 0; k < N_i; k++)
    {
        if (sm_3D_vecs[k][0] != -1)
        {
            for (v = 0; v < 4; v++)
            {
                var = sm_3D_vecs[k][v];
                sm_3D_vecs[k][v] = curr[v];
                curr[v] = var;
            }
        }
    }
    *N_circ++;
}

int already_matched(int sm_IS[][2], int indx, int N_i)
{   
    int i;
    for (i = 0; i < N_i; i++)
    {
        if (sm_IS[i][0] == indx)
        {
            return 1;
        }
    }
    return 0;
}

void sm_4_star(double four_stars[][4], double sm_3D_vecs[][4], int sm_IS[][2], double body_vecs_IS[][4], int sm_K_vec_arr[][3], int *N_match, int N_i, double q, double m, int N_is, double matched[][3])
{   
    int i=0;
    int j=0;
    int k=0;

    int SIM[N_GC][6];
    memset(SIM, 0, N_GC * sizeof(SIM[0]));

    int SIM_indx_arr[1000];
    memset(SIM_indx_arr,0, 1000*sizeof(SIM_indx_arr[0]));

    int SIM_flags[N_GC];
    memset(SIM_flags, 0, N_GC*sizeof(SIM_flags[0]));

    int top_indx = 0;

    

    
    double p[6];
    int ct = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = i + 1; j < 4; j++)
        {
             double norm1 = sqrt(four_stars[i][1] * four_stars[i][1] + four_stars[i][2] * four_stars[i][2] + four_stars[i][3] * four_stars[i][3]);
             double norm2 = sqrt(four_stars[j][1] * four_stars[j][1] + four_stars[j][2] * four_stars[j][2] + four_stars[j][3] * four_stars[j][3]);
            p[ct] = fabs(four_stars[i][1] * four_stars[j][1] + four_stars[i][2] * four_stars[j][2] + four_stars[i][3] * four_stars[j][3]) / (norm1 * norm2);
            ct++;
        }
    }

    int checks[4][6] = {{1, 1, 1, 0, 0, 0},
                        {1, 0, 0, 1, 1, 0},
                        {0, 1, 0, 1, 0, 1},
                        {0, 0, 1, 0, 1, 1}};

    for (j = 0; j < 6; j++)
    {
        double sin_j = sqrt(1 - (p[j] * p[j]));
        int k_top = ceil((cos(DELTA) * p[j] + sin_j * sin(DELTA) - q) / m);
        int k_bot = floor((cos(DELTA) * p[j] - sin_j * sin(DELTA) - q) / m);
        if (k_top <= 0 || k_bot >= 224792)
        {
            printf("bad values : k_bot = %d, k_top = %d\n", k_bot, k_top);
            continue;
        }
        else
        {
            if (k_top > 224792)
            {
                k_top = 224792;
            }
            if (k_bot < 0)
            {
                k_bot = 1;
            }
            int k_start = sm_K_vec_arr[k_bot - 1][2] + 1;
            int k_end = sm_K_vec_arr[k_top - 1][2];
            if (k_start == k_end)
            {
                SIM[sm_K_vec_arr[k_end - 1][0]][j] = 1;
                SIM_indx_arr[top_indx]=sm_K_vec_arr[k_end-1][0];
                top_indx ++;
            }
            else
            {
                for (i = k_start; i <= k_end; i++)
                {
                    SIM[sm_K_vec_arr[i - 1][0]][j] = 1; 
                    SIM[sm_K_vec_arr[i - 1][1]][j] = 1;

                    if(SIM_flags[sm_K_vec_arr[i - 1][0]] == 0) {
                        SIM_flags[sm_K_vec_arr[i - 1][0]] = 1;
                        SIM_indx_arr[top_indx]=sm_K_vec_arr[i-1][0];
                        top_indx++;
                    }

                    if(SIM_flags[sm_K_vec_arr[i - 1][1]] == 0) {
                        SIM_flags[sm_K_vec_arr[i - 1][1]] = 1;
                        SIM_indx_arr[top_indx]=sm_K_vec_arr[i-1][1];
                        top_indx++;
                    }
                    
                }
            }
        }
    }


    for (j = 0; j < 4; j++)
    {
        int matched_rows = 0;
        int temp = 0;  

        for(i=0; i<top_indx; i++)
        {
            SIM_flags[SIM_indx_arr[i]] = 1;
        }

        for (i = 0; i <top_indx; i++)
        {
             k = SIM_indx_arr[i];

            if(SIM_flags[SIM_indx_arr[i]]==1)
            {

                if (SIM[k][0] == checks[j][0] && SIM[k][1] == checks[j][1] && SIM[k][2] == checks[j][2] && SIM[k][3] == checks[j][3] && SIM[k][4] == checks[j][4] && SIM[k][5] == checks[j][5])
                {
                    matched_rows++;
                    temp = k;
                }
                SIM_flags[k]=0;
            }

        }
        
        if (matched_rows == 1)
        {
            int flag = already_matched(sm_IS, (int)four_stars[j][0], N_is);
            if (flag == 0)
            {
                (*N_match)++;
                for (k = 0; k < N_i; k++)
                {
                    if ((int)sm_3D_vecs[k][0] == (int)four_stars[j][0])
                    {
                        sm_3D_vecs[k][0] = -1;
                        break;
                    }
                }
                k = 0;
                for (k = 0; k < N_GC; k++)
                {
                    if (sm_IS[k][0] == -1)
                    {
                        sm_IS[k][0] = (int)four_stars[j][0];
                        body_vecs_IS[k][0] = four_stars[j][0];
                        sm_IS[k][1] = temp;
                        i = 1;
                        for (i = 1; i < 4; i++){
                            body_vecs_IS[k][i] = four_stars[j][i];
                            // matched[k][i-1] = four_stars[j][i-1];
                        }
                        
                        break;
                    }
                }
                // k=0;
                // for (k = 0; k < 40; k++)
                // {
                //     matched[k][2] = body_vecs_IS[k][0];
                // }
                // k=0;
                // for (k = 0; k < 40; k++)
                // {
                //     i=0;
                //     for (i = 0; i < 2; i++)
                //     {
                //         matched[k][i] = body_vecs_IS[k][i+1];
                //     }
                // }
            }
        }
    }
    
}

void sm_gnrt_3D_vec(double sm_3D_vecs[][4], double sm_sorted_UIS[][3], int N_i)
{   int i;
    for (i = 0; i < N_i; i++)
    {
        double local = sqrt((sm_sorted_UIS[i][1] / FOCAL_LENGTH) * (sm_sorted_UIS[i][1] / FOCAL_LENGTH) + (sm_sorted_UIS[i][2] / FOCAL_LENGTH) * (sm_sorted_UIS[i][2] / FOCAL_LENGTH) + 1);
        sm_3D_vecs[i][0] = sm_sorted_UIS[i][0];
        sm_3D_vecs[i][1] = (sm_sorted_UIS[i][1] / (FOCAL_LENGTH * local));
        sm_3D_vecs[i][2] = (sm_sorted_UIS[i][2] / (FOCAL_LENGTH * local));
        sm_3D_vecs[i][3] = 1 / local;
        //printf("%.16f %.16f %.16f\n", sm_3D_vecs[i][3], sm_3D_vecs[i][1], sm_3D_vecs[i][2]);

    }
}

void sm_validate(double sm_3D_vecs[][4], int sm_IS[][2], double body_vecs_IS[][4], double sm_GC[][4], int *N_is, int N_i, double tol, double p_1, double p_2){
    int i, j;
    int N_new = *N_is;
    int votes[N_i];
    memset(votes, 0, N_i * sizeof(votes[0]));
    // for (i = 0; i<N_i; i++){
    //     for (j = 0; j<4; j++){
    //         printf("%Lf         ", sm_3D_vecs[i][j]);
    //     }
    //     printf("\n");
    // }

    for (i = 0; i < N_i; i++){
        if (sm_IS[i][0] == -1)
            continue;
        for (j = i+1; j< N_i-1; j++){
            if (sm_IS[j][0] != -1){
                double d_ij = fabs(body_vecs_IS[i][1]*body_vecs_IS[j][1] + body_vecs_IS[i][2]*body_vecs_IS[j][2] + body_vecs_IS[i][3]*body_vecs_IS[j][3]);
                // printf("i: %d, j: %d, norm1: %Lf, norm2: %Lf, d_ij: %Lf\n", i, j, norm1, norm2, d_ij);
                double d_ij_gc = fabs(sm_GC[sm_IS[i][1] - 1][1]*sm_GC[sm_IS[j][1] - 1][1] + sm_GC[sm_IS[i][1] - 1][2]*sm_GC[sm_IS[j][1] - 1][2] + sm_GC[sm_IS[i][1] - 1][3]*sm_GC[sm_IS[j][1] - 1][3]);
                // printf("i: %d, j: %d, norm_gc_1: %Lf, norm_gc_2: %Lf, d_ij_gc: %Lf\n", i, j, norm_gc_1, norm_gc_2, d_ij_gc);
                if (fabs(d_ij/d_ij_gc - 1) < tol/100){
                    votes[i]++;
                    votes[j]++;
                }
                if((votes[i] < 0) || (votes[j] < 0))
                    printf("Idx i: %d, j: %d\n", i, j);
            }
        }
    }
    int N_LB = p_1*(N_new)/100;
    // printf("\nN_LB: %d \n", N_LB);
    for (i = 0; i < N_i; i++){
        if (sm_IS[i][0] == -1)
            continue;
        if (votes[i] < N_LB){
            sm_IS[i][0] = -1;
            body_vecs_IS[i][0] = -1;
            // printf("%d discarded, %d votes \n", i, votes[i]);
            N_new--;
        }
    }
    int N_UB = p_2*(N_new)/100;
    // printf("N_UB: %d \n", N_UB);
    for (i = 0; i < N_i; i++){
        if (sm_IS[i][0] == -1)
            continue;
        if (votes[i] < N_UB){
            sm_IS[i][0] = -1;
            body_vecs_IS[i][0] = -1;
            // printf("%d discarded, %d votes \n\n\n", i, votes[i]);
            N_new--;
        }
    }
    *N_is = N_new;
}

void bubbleSort(double arr[][3], int n)
{
    int i, j;
    double temp[3];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j][1] * arr[j][1] + arr[j][2] * arr[j][2] > arr[j + 1][1] * arr[j + 1][1] + arr[j + 1][2] * arr[j + 1][2])
            {
                // swap the elements
                temp[0] = arr[j][0];
                temp[1] = arr[j][1];
                temp[2] = arr[j][2];
                arr[j][0] = arr[j + 1][0];
                arr[j][1] = arr[j + 1][1];
                arr[j][2] = arr[j + 1][2];
                arr[j + 1][0] = temp[0];
                arr[j + 1][1] = temp[1];
                arr[j + 1][2] = temp[2];
            }
        }
    }
}

void starMatching(double centroids_st[MAX_STARS][3], int tot_stars, double data[3][MAX_STARS], int input_ids[50], int star_ids[50], int* matched_stars, double matched[][3]){
    // Maximum number of iterations
    int N_max = tot_stars - 1;

    // N_i is specially used to pass onto other functions(remains unchanged throughout the code)
    int N_i = tot_stars;

    // N_uis is a SM variable for the total number of stars detected from the FE block(changes value in main)
    int N_uis = tot_stars;

    // Required maximum number of matched stars
    int N_th = 20;

    // Number of times control goes to sm_4_star_circulate
    int N_circ = 0;

    // Number of identified stars
    int N_is = 0;

    double tol = TOL;
    double p_1 = P1;
    double p_2 = P2;

    // Constants for using the k vector table (to be used in the 4 star matching)
    double m = (Y_MAX - Y_MIN + 2 * EPSILON) / (N_KVEC_PAIRS - 1);
    double q = Y_MIN - m - EPSILON;

    // Array for storing the Identified Stars
    int sm_IS[N_GC][2];
    // Array for storing body-frame vectors of Identified stars
    double body_vecs_IS[N_GC][4];

    //Sorting the UIS according to euclidean distance
    bubbleSort(centroids_st, N_i);

    // Initialize a block of memory to -1
    memset(sm_IS, -1, N_GC * sizeof(sm_IS[0]));
    memset(body_vecs_IS, -1, N_GC * sizeof(body_vecs_IS[0]));

    // Array to store body frame vectors of extracted stars from FE block
    double sm_3D_vecs[N_i][4];

    // Generate body frame vectors for extracted stars
    sm_gnrt_3D_vec(sm_3D_vecs, centroids_st, N_i);

    int i = 1;
    for (i = 1; i <= N_max; i++){
        
            if (N_uis >= 4 && N_is < N_th){
            // Variable for storing the number of stars matched in a particular iteration    
            int N_match = 0;

            // This will store the 4 stars used from the sm_3D_vecs table
            double four_stars[4][4];

            // Here the variable countt is used just to count whether 4 stars have been extracted
            int countt = 0, j = 0;

            // Pick body vectors of 4 stars from sm_3D_vecs
            for (countt = 0, j = 0; j < N_i && countt < 4; j++){
                if ((int)sm_3D_vecs[j][0] != -1){
                    int k = 0;
                    for (k = 0; k < 4; k++){
                        four_stars[countt][k] = sm_3D_vecs[j][k];
                    }
                    countt++;
                }
            }

            // Perform 4 star N star algorithm
            sm_4_star(four_stars, sm_3D_vecs, sm_IS, body_vecs_IS, sm_K_vec_arr, &N_match, N_i, q, m, N_is, matched);

            // Decrement matched stars from those detected from FE block for next iteration
            N_uis -= N_match;

            // Increment identified stars with matched stars
            N_is += N_match;

            // If no star is detected or we have not circulated the sm_3D_vecs exhaustively
            if (N_match == 0 && N_circ <= 2 * N_i){

                sm_4_star_circulate(sm_3D_vecs, &N_circ, N_i);
            
                if (N_circ >= 2 * N_i){
                    break;
                }
            }
        }

        // End the loop if the conditions "N_uis >= 4 && N_is < N_th" are not satisfied
        else{
            break;
        }
    }

    sm_validate(sm_3D_vecs, sm_IS, body_vecs_IS, sm_GC, &N_is, N_i, tol, p_1, p_2);

    // Index variable for organising SM output
    int data_index = 0;

    double x_inertial[N_is], y_inertial[N_is], z_inertial[N_is];

    // Arrays to store body frame components (sm_3D_vec)
    double x_body[N_is], y_body[N_is], z_body[N_is];

    // If "some" stars match
    if (N_is != 0){
        *matched_stars = N_is;
    
        i = 0;
        for (i = 0; i < N_GC; i++){
            if (sm_IS[i][0] != -1){

                // Storing input ids and SSP ids in respective arrays
                input_ids[data_index] = sm_IS[i][0];
                star_ids[data_index] = sm_IS[i][1];

                // Storing body frame components of matched stars in respective arrays
                x_body[data_index] = body_vecs_IS[i][3];
                y_body[data_index] = body_vecs_IS[i][1];
                z_body[data_index] = body_vecs_IS[i][2];

                // Storinng intertial frame components of matched stars in respective arrays
                x_inertial[data_index] = sm_GC[sm_IS[i][1] - 1][1];
                y_inertial[data_index] = sm_GC[sm_IS[i][1] - 1][2];
                z_inertial[data_index] = sm_GC[sm_IS[i][1] - 1][3];

                // Increment data_index for next matched star
                data_index++;
            }
        }

        // Arranging SM output vectors compatible with Estimation block
        int j = 0;
        for (j = 0; j < N_is; j++){
            data[0][j] = x_inertial[j];
            data[0][j + N_is] = x_body[j];
            data[1][j] = y_inertial[j];
            data[1][j + N_is] = y_body[j];
            data[2][j] = z_inertial[j];
            data[2][j + N_is] = z_body[j];

        }
    }

    // If no star is matched
    else{
        printf("\nCOULD NOT MATCH STARS\n");
        
    }
}

// Estimation structures
typedef struct Vec Vec; 
struct Vec{
    double x, y, z;
};

typedef struct Matrix_3 Matrix_3;
struct Matrix_3{
    double elements[3][3];
};

//Estimation functions
double absoluteValue(double x){
    if (x < 0){
        x = -x;
    }
    return x;
}

double square_root(double x){
    double guess = 1;
    while (absoluteValue((guess * guess) / x - 1.0) >= 0.00001)
    {
        guess = ((x / guess) + guess) / 2;
    }
    return guess;
}

void Vec_construct(Vec *temp, double x, double y, double z){
    temp->x = x;
    temp->y = y;
    temp->z = z;
}

void cross(Vec *temp, Vec *vec1, Vec *vec2){
    temp->x = vec1->y * vec2->z - vec1->z * vec2->y;
    temp->y = vec1->z * vec2->x - vec1->x * vec2->z;
    temp->z = vec1->x * vec2->y - vec1->y * vec2->x;
}

double dot(Vec *vec1, Vec *vec2){
    return (vec1->x * vec2->x + vec1->y * vec2->y + vec1->z * vec2->z);
}

void scale_vec(Vec *temp, double k){
    temp->x *= k;
    temp->y *= k;
    temp->z *= k;
}

void add_vec(Vec *temp, Vec *vec1, Vec *vec2){
    temp->x = vec1->x + vec2->x;
    temp->y = vec1->y + vec2->y;
    temp->z = vec1->z + vec2->z;
}

void matrix_construct(Matrix_3 *temp, unsigned char x){
    if (x == 0)
    { // Null Matrix
        unsigned char i = 0;
        unsigned char j = 0;
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                (temp->elements)[i][j] = 0;
            }
        }
    }
    if (x == 1)
    {
        unsigned char i = 0;
        unsigned char j = 0;
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                if (i == j)
                    (temp->elements)[i][j] = 1;
                else
                    (temp->elements)[i][j] = 0;
            }
        }
    }
}

void add_matrix(Matrix_3 *temp, Matrix_3 *m1, Matrix_3 *m2){
    unsigned char i = 0;
    unsigned char j = 0;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            (temp->elements)[i][j] = (m1->elements)[i][j] + (m2->elements)[i][j];
        }
    }
}

void scale_matrix(Matrix_3 *temp, Matrix_3 *M, double k){
    unsigned char i = 0;
    unsigned char j = 0;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            temp->elements[i][j] = (M->elements[i][j]) * k;
        }
    }
}

double det(Matrix_3 *temp){
    double x;
    x = temp->elements[0][0] * (temp->elements[1][1] * temp->elements[2][2] - temp->elements[1][2] * temp->elements[2][1]) - temp->elements[0][1] * (temp->elements[1][0] * temp->elements[2][2] - temp->elements[1][2] * temp->elements[2][0]) + temp->elements[0][2] * (temp->elements[1][0] * temp->elements[2][1] - temp->elements[1][1] * temp->elements[2][0]);
    return x;
}

double trace(Matrix_3 *temp){
    return (temp->elements[0][0] + temp->elements[1][1] + temp->elements[2][2]);
}

void T(Matrix_3 *temp, Matrix_3 *m1){
    unsigned char i = 0;
    unsigned char j = 0;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            temp->elements[i][j] = m1->elements[j][i];
        }
    }
}

void adjoint(Matrix_3 *temp, Matrix_3 *m){
    temp->elements[0][0] = m->elements[1][1] * m->elements[2][2] - m->elements[1][2] * m->elements[2][1];
    temp->elements[1][0] = m->elements[1][2] * m->elements[2][0] - m->elements[1][0] * m->elements[2][2];
    temp->elements[2][0] = m->elements[1][0] * m->elements[2][1] - m->elements[2][0] * m->elements[1][1];
    temp->elements[0][1] = m->elements[2][1] * m->elements[0][2] - m->elements[0][1] * m->elements[2][2];
    temp->elements[1][1] = m->elements[0][0] * m->elements[2][2] - m->elements[2][0] * m->elements[0][2];
    temp->elements[2][1] = m->elements[2][0] * m->elements[0][1] - m->elements[0][0] * m->elements[2][1];
    temp->elements[0][2] = m->elements[0][1] * m->elements[1][2] - m->elements[1][1] * m->elements[0][2];
    temp->elements[1][2] = m->elements[1][0] * m->elements[0][2] - m->elements[0][0] * m->elements[1][2];
    temp->elements[2][2] = m->elements[0][0] * m->elements[1][1] - m->elements[1][0] * m->elements[0][1];
}

void outer_product(Matrix_3 *temp, Vec *v1, Vec *v2){
    temp->elements[0][0] = v1->x * v2->x;
    temp->elements[0][1] = v1->x * v2->y;
    temp->elements[0][2] = v1->x * v2->z;
    temp->elements[1][0] = v1->y * v2->x;
    temp->elements[1][1] = v1->y * v2->y;
    temp->elements[1][2] = v1->y * v2->z;
    temp->elements[2][0] = v1->z * v2->x;
    temp->elements[2][1] = v1->z * v2->y;
    temp->elements[2][2] = v1->z * v2->z;
}

void matmul(Vec *temp, Matrix_3 *M, Vec *v){
    temp->x = M->elements[0][0] * v->x + M->elements[0][1] * v->y + M->elements[0][2] * v->z;
    temp->y = M->elements[1][0] * v->x + M->elements[1][1] * v->y + M->elements[1][2] * v->z;
    temp->z = M->elements[2][0] * v->x + M->elements[2][1] * v->y + M->elements[2][2] * v->z;
}

void quest(double Q[4], double data[3][MAX_STARS], int N){
    Vec b[N];
    Vec r[N];

    double epsilon_est = 0.001;
    double a[N];
    int c;
    for (c = 0; c < N; c++)
    {
        a[c] = 0;
        a[c] = 1.0 / N;
    }
    double guess = 1;
    Matrix_3 B, S;
    matrix_construct(&B, 0);
    matrix_construct(&S, 0);
    Vec z;
    Vec_construct(&z, 0, 0, 0);
    unsigned char counter = 0;
    double K;
    double norm_2_z;
    double det_S;
    double trace_B;
    double zTSz;
    double zTSSz;
    Matrix_3 adj_S;
    Matrix_3 check;
    int i = 0;
    for (i = 0; i < N * 2; i++)
    { 
        if (i < N)
        {
            // Construct r vectors according to the data read
            Vec_construct(&r[i], data[0][i], data[1][i], data[2][i]);
        }
        else
        {
             // Construct b vectors according to the data read
            Vec_construct(&b[i - N], data[0][i], data[1][i], data[2][i]);
        }
    }

    // Perform Quest
    double f(double x)
    { // This is the function whose roots are found using Newton Rapshon
        double first = (x * x - trace_B * trace_B + K) * (x * x - trace_B * trace_B - norm_2_z);
        double second = (x - trace_B) * (zTSz + det_S);
        return (first - second - zTSSz);
    }
    double f_bar(double x)
    { // This is the derivative of the function whose roots are found using Newton Rapshon
        double temp = 2 * x * (2 * x * x - 2 * trace_B * trace_B + K - norm_2_z);
        return (temp - zTSz - det_S);
    }

    do
    {
        if (counter == 1)
        { // Perform rotation abt x --- {x,-y,-z}
            for (i = 0; i < N; i++)
            {
                r[i].y = -1 * r[i].y;
                r[i].z = -1 * r[i].z;
            }
        }
        if (counter == 2)
        { // Perform rotation abt y --- {-x,y,-z}
            for (i = 0; i < N; i++)
            {
                r[i].x = -1 * r[i].x;
                r[i].z = -1 * r[i].z;
            }
        }
        if (counter == 3)
        { // Perform rotation abt z --- {-x,-y,z}
            for (i = 0; i < N; i++)
            {
                r[i].y = -1 * r[i].y;
                r[i].x = -1 * r[i].x;
            }
        }
        for (i = 0; i < N; i++)
        {
            Matrix_3 temp;
            outer_product(&temp, &(b[i]), &(r[i]));
            scale_matrix(&temp, &temp, a[i]);
            add_matrix(&B, &B, &temp);
        }
        for (i = 0; i < N; i++)
        {
            Vec temp;
            cross(&temp, &(b[i]), &(r[i]));
            scale_vec(&temp, a[i]);
            add_vec(&z, &z, &temp);
        }

        trace_B = trace(&B);
        T(&S, &B);              
        add_matrix(&S, &S, &B); 
        adjoint(&adj_S, &S);    
        K = trace(&adj_S);
        norm_2_z = z.x * z.x + z.y * z.y + z.z * z.z;
        det_S = det(&S);  
        Vec temp1, temp2; 
        matmul(&temp1, &S, &z);
        zTSz = dot(&z, &temp1);
        matmul(&temp2, &S, &temp1);
        zTSSz = dot(&z, &temp2);

        double rho = guess; 
        for (i = 0; i < 10; i++)
        {
            double x1 = f(rho);
            double x2 = f_bar(rho);
            rho = rho - x1 / x2;
        }
        rho += trace_B;
        Matrix_3 I;
        matrix_construct(&I, 1);
        scale_matrix(&I, &I, rho);
        scale_matrix(&check, &S, -1);
        add_matrix(&check, &check, &I);

        counter++;
    } while ((det(&check)) < epsilon_est);

    Vec q_3;
    double q_4 = det(&check);
    Matrix_3 adj_check;
    adjoint(&adj_check, &check);
    matmul(&q_3, &adj_check, &z);

    if (counter == 2)
    {
        Vec e;
        Vec_construct(&e, 1, 0, 0);
        double new_4 = -dot(&e, &q_3);
        cross(&q_3, &e, &q_3);
        scale_vec(&e, q_4);
        add_vec(&q_3, &q_3, &e);
        q_4 = new_4;
    }
    else if (counter == 3)
    {
        Vec e;
        Vec_construct(&e, 0, 1, 0);
        double new_4 = -dot(&e, &q_3);
        cross(&q_3, &e, &q_3);
        scale_vec(&e, q_4);
        add_vec(&q_3, &q_3, &e);
        q_4 = new_4;
    }
    else if (counter == 4)
    {
        Vec e;
        Vec_construct(&e, 0, 0, 1);
        double new_4 = -dot(&e, &q_3);
        cross(&q_3, &e, &q_3);
        scale_vec(&e, q_4);
        add_vec(&q_3, &q_3, &e);
        q_4 = new_4;
    }
    double alpha;
    double x1 = (q_3.x);
    double x2 = (q_3.y);
    double x3 = (q_3.z);
    double x4 = (q_4);
    alpha = square_root(x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4);
    Q[0] = x1 / alpha;
    Q[1] = x2 / alpha;
    Q[2] = x3 / alpha;
    Q[3] = x4 / alpha;
    return;
}

// void OILS(short arr_out_img[BREADTH + 2][LENGTH + 2]){
int OILS(double centroids_st[][3], int tot_stars, double matched[][3]){
    
    // // Store start time of FE and OILS
    //  uint64_t FE_START = cortos_get_clock_time();

    // int tot_stars = 0;
    // double centroids_st[MAX_STARS][3];

    // // Perform Feature Extraction
    // regionGrowth(arr_out_img, centroids_st, &tot_stars);
    // // Store stop time of FE
    //  uint64_t FE_END = cortos_get_clock_time();

    // printf("%d\n", tot_stars);

    // // Print FE output on the terminal
     int i = 0;
     int j = 0;
    
    // Variable passed by reference to starMatcing function to get number of mathched stars
    int matched_stars = 0;
    // Array to store input ids of matched stars
    int input_ids[MAX_STARS] = {0};
    // Array to store SSP ids of matched stars
    int star_ids[MAX_STARS] = {0};
   
    // Array to store the input of Estimation block
    double data[3][MAX_STARS];

    

    // Store start time of SM
    // uint64_t SM_START = cortos_get_clock_time();

    // Perform Star Matching
    starMatching(centroids_st, tot_stars, data, input_ids, star_ids, &matched_stars, matched);

    // Store stop time of SM
    // uint64_t SM_END = cortos_get_clock_time();

    // Prints the number of stars matched by the SM block on the terminal
    printf("\nTotal matched stars = %d\n\n", matched_stars);

    // Prints the inertial frame vectors and body frame vectors of matched stars on the terminal
    if(matched_stars != 0){
        for (j = 0; j < matched_stars; j++){
            // printf("%d %d ", input_ids[j], star_ids[j]);
            matched[j][2] = star_ids[j];

            // printf("%.16f %.16f %.16f ", data[0][j], data[1][j], data[2][j]);
             matched[j][0] = data[0][j+matched_stars];
             matched[j][1] = data[1][j+matched_stars];

            // printf("%.16f %.16f %.16f\n",data[0][j + matched_stars], data[1][j + matched_stars], data[2][j + matched_stars]);
        }




        // Store start time of Estimation
        // uint64_t EST_START = cortos_get_clock_time();
        
        // Array to store the quaternion output of Estimation block
        // double out[4];

        // // Perform QUEST
        // quest(out, data, matched_stars);

        // // Printing the quaternion values on the terminal
        // printf("\n");
        // printf("quat_3.x = %.15f\n", out[0]);
        // printf("quat_3.y = %.15f\n", out[1]);
        // printf("quat_3.z = %.15f\n", out[2]);
        // printf("quat_4   = %.15f\n", out[3]);
    


        //Store stop time of Estimation and OILS
        // uint64_t EST_END = cortos_get_clock_time();

        // CoRTOS operation
        // uint32_t t_FEs = FE_START & (0xffffffff);
        // uint32_t t_FEe = FE_END & (0xffffffff);
        // uint32_t t_SMs = SM_START & (0xffffffff);
        // uint32_t t_SMe = SM_END & (0xffffffff);
        // uint32_t t_ESTs = EST_START & (0xffffffff);
        // uint32_t t_ESTe = EST_END & (0xffffffff);
        
        // Print run times
        printf("\n");
        // printf("FE__run_time =  %u\n", t_FEe - t_FEs);
        // printf("SM__run_time =  %u\n", t_SMe - t_SMs);
        // printf("EST_run_time =  %u\n", t_ESTe - t_ESTs);
        // printf("Net_run_time =  %u\n", t_ESTe - t_FEs);
    }
    return matched_stars;
}
*/


// TM functions

//sorts 2 column array in ascending order of column indexed ind
void bubbleSort2(double arr[][2], int n, int ind)
{
    int i, j;
    long double temp[2];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j][ind] > arr[j + 1][ind])
            {
                // swap the elements
                temp[0] = arr[j][0];
                temp[1] = arr[j][1];

                arr[j][0] = arr[j + 1][0];
                arr[j][1] = arr[j + 1][1];

                arr[j + 1][0] = temp[0];
                arr[j + 1][1] = temp[1];
            }
        }
    }
}

//sorts 3 column array in ascending order of column indexed ind
void bubbleSort3(double arr[][3], int n, int ind)
{
    int i, j;
    double temp[3];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j][ind] > arr[j + 1][ind])
            {
                // swap the elements
                temp[0] = arr[j][0];
                temp[1] = arr[j][1];
                temp[2] = arr[j][2];

                arr[j][0] = arr[j + 1][0];
                arr[j][1] = arr[j + 1][1];
                arr[j][2] = arr[j + 1][2];

                arr[j + 1][0] = temp[0];
                arr[j + 1][1] = temp[1];
                arr[j + 1][2] = temp[2];
            }
        }
    }
}

//angular distance cos_theta between two stars in Guide Star Catalogue
double sm_TM_id_angdist(int id_1,
                        int id_2,
                        double sm_GD_SC[5060][4])
{
double dot, norm1, norm2, dist;
dot = sm_GD_SC[id_1][3]*sm_GD_SC[id_2][3] + sm_GD_SC[id_1][1]*sm_GD_SC[id_2][1] + sm_GD_SC[id_1][2]*sm_GD_SC[id_2][2];
norm1 = sqrt(sm_GD_SC[id_1][3]*sm_GD_SC[id_1][3] + sm_GD_SC[id_1][1]*sm_GD_SC[id_1][1] + sm_GD_SC[id_1][2]*sm_GD_SC[id_1][2]);
norm2 = sqrt(sm_GD_SC[id_2][3]*sm_GD_SC[id_2][3] + sm_GD_SC[id_2][1]*sm_GD_SC[id_2][1] + sm_GD_SC[id_2][2]*sm_GD_SC[id_2][2]);

dist = dot/(norm1*norm2);
return dist;

}

//angular distance cos_theta between two stars with pixel coordinates x_1,y_1 and x_2,y_2
double sm_TM_centroid_angdist(double x_1, double y_1,
                              double x_2, double y_2,
                              double sm_TM_F,
                              double pixel_size)
{
double dot, norm1, norm2, dist;
dot = (x_1*pixel_size)*(x_2*pixel_size) + (y_1*pixel_size)*(y_2*pixel_size) + (sm_TM_F)*(sm_TM_F);
norm1 = sqrt((x_1*pixel_size)*(x_1*pixel_size) + (y_1*pixel_size)*(y_1*pixel_size) + (sm_TM_F)*(sm_TM_F));
norm2 = sqrt((x_2*pixel_size)*(x_2*pixel_size) + (y_2*pixel_size)*(y_2*pixel_size) + (sm_TM_F)*(sm_TM_F));

dist = dot/(norm1*norm2);
return dist;
}

//picks out common star ID's from previous two frames prevmat_1_unsorted and prevmat_2_unsorted, and stores them 
//in frames prevmat_1 and prevmat_2
int sm_TM_CP_inputsort(double prevmat_1_raw[][3],
                       int N1,
                       double prevmat_2_raw[][3],
                       int N2,
                       double sm_TM_prevmat_1[][3],
                       double sm_TM_prevmat_2[][3])
{
bubbleSort3(prevmat_1_raw, N1, 2);

bubbleSort3(prevmat_2_raw, N2, 2);


int i=0, j=0, k=0;
// for(i = 0; i < 20; i++){
// printf("%f %f %d\n", prevmat_1_raw[i][0], prevmat_1_raw[i][1], (int)prevmat_1_raw[i][2]);
// }
// printf("\n");
// for(i = 0; i < 20; i++){
// printf("%f %f %d\n", prevmat_2_raw[i][0], prevmat_2_raw[i][1], (int)prevmat_2_raw[i][2]);
// }

while ((i<N1) && (j<N2))
{

        
            if (prevmat_1_raw[i][2] == prevmat_2_raw[j][2])
            {
               sm_TM_prevmat_1[k][0] = prevmat_1_raw[i][0];
               sm_TM_prevmat_1[k][1] = prevmat_1_raw[i][1];
               sm_TM_prevmat_1[k][2] = prevmat_1_raw[i][2];
               sm_TM_prevmat_2[k][0] = prevmat_2_raw[j][0];
               sm_TM_prevmat_2[k][1] = prevmat_2_raw[j][1];
               sm_TM_prevmat_2[k][2] = prevmat_2_raw[j][2];
               i++; j++; k++;
               //printf("%d %d %d\n", i, j, k);

            }

            else if (prevmat_1_raw[i][2] < prevmat_2_raw[j][2]) i++;
            else if (prevmat_1_raw[i][2] > prevmat_2_raw[j][2]) j++;
        
}
return k;
}

//takes in two previous frames of common stars, and outputs the predicted matrix sm_TM_CP_predmat
//containing predicted centroids with star ID's in next frame
//All frames have x coordinate in first column, y in second column, star ID in third
void sm_TM_CP (double sm_TM_CP_prevmat_1[][3],
               double sm_TM_CP_prevmat_2[][3],
               double sm_TM_CP_F,
               double sm_TM_CP_predmat[][3],
               int N,
               double pixel_size)
{
//phi, theta, psi are the roll, pitch, yaw obtained using previous two frames
double phi = 0, theta = 0, psi = 0;
double u, v, f = sm_TM_CP_F;
double dun, dvn;
double m[2*N][3];

//m is the (2N x 3) matrix made from the star coordinates. Refer CDR for actual terms in the matrix.
int k=0;
for (k=0; k<N; k++)
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

//m2 stores the matrix product (m transpose)m.
double m2[3][3] = {0,0,0,0,0,0,0,0,0};

for (k=0; k<2*N; k++)
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

//determinant of m2 matrix
double determinant = 0;
double minv[3][3];

int i=0;
for(i=0;i<3;i++)
    determinant = determinant + (m2[0][i]*(m2[1][(i+1)%3]*m2[2][(i+2)%3] - m2[1][(i+2)%3]*m2[2][(i+1)%3]));

//inverse of m2 matrix
int j=0;
for(i=0;i<3;i++)
{
    for(j=0;j<3;j++)
        minv[i][j] = ((m2[(i+1)%3][(j+1)%3] * m2[(i+2)%3][(j+2)%3]) - (m2[(i+1)%3][(j+2)%3]*m2[(i+2)%3][(j+1)%3]))/ determinant;
}


//to store roll, pitch, yaw: m3 holds the matrix product inverse[(m_transpose)m] x m
double m3[3][2*N];
for (k=0; k<2*N; k++)
{
    m3[0][k] = minv[0][0]*m[k][0] + minv[0][1]*m[k][1] + minv[0][2]*m[k][2];
    m3[1][k] = minv[1][0]*m[k][0] + minv[1][1]*m[k][1] + minv[1][2]*m[k][2];
    m3[2][k] = minv[2][0]*m[k][0] + minv[2][1]*m[k][1] + minv[2][2]*m[k][2];
}


//storing delta_u and delta_v in col array
double col[2*N];
for (k=0; k<N; k++)
{
    col[2*k] = (sm_TM_CP_prevmat_2[k][0] - sm_TM_CP_prevmat_1[k][0]) * pixel_size;
    col[2*k+1] = (sm_TM_CP_prevmat_2[k][1] - sm_TM_CP_prevmat_1[k][1]) * pixel_size;
}

//solve for roll pitch yaw
for (k=0; k<2*N; k++)
{
    phi += m3[0][k]*col[k];
    theta += m3[1][k]*col[k];
    psi += m3[2][k]*col[k];
}


//estimation of next centroids based on the roll pitch yaw obtained from previous two frames
for (k=0; k<N; k++)
{
    u = sm_TM_CP_prevmat_2[k][0] * pixel_size;
    v = sm_TM_CP_prevmat_2[k][1] * pixel_size;

    dun = (u*v/f)*phi + (-f - u*u/f)*theta + v*psi;
    dvn = (f + v*v/f)*phi + (-u*v/f)*theta + (-u)*psi;

    sm_TM_CP_predmat[k][0] = sm_TM_CP_prevmat_2[k][0] + dun/pixel_size;
    sm_TM_CP_predmat[k][1] = sm_TM_CP_prevmat_2[k][1] + dvn/pixel_size;
    sm_TM_CP_predmat[k][2] = sm_TM_CP_prevmat_1[k][2];
}
}

//checks for stars using radius based matching as outlined in the CDR.
int sm_TM_RBM (double sm_TM_RBM_predmat[][3], //predicted centroids from CP block
                int N, //length of predicted centroid array
                double sm_TM_RBM_truemat[][2], //true centroids from feature extraction
                int M, //length of true centroid array
                double sm_TM_RBM_radius, //radius of window to check for matches
                double sm_TM_RBM_matchmat[][3] //matched star array with true centroids stored
                )

{
int stars_in_range = 0;
int sole_matched = -1;
double r = sm_TM_RBM_radius;
int matched = -1;

int i=0;
int j=0;
double dx=0;
double dy=0;
for (i=0; i<N; i++)
{
    stars_in_range = 0;
    sole_matched = -1;
    /*for (int j = start_ind; j<M; j++)
    {
        double dx = abs(sm_TM_RBM_truemat[j][0] - sm_TM_RBM_predmat[i][0]);
        if (dx<r)
        {
            start_ind = j;
            double dy = abs(sm_TM_RBM_truemat[j][1] - sm_TM_RBM_predmat[i][1]);
            if (dy<r)
            {
                stars_in_range++;
                if (stars_in_range > 1) break;
                sole_matched = j;
            }
        }
    }*/
    for (j=0; j<N; j++)
    {
        dx = (sm_TM_RBM_truemat[j][0] - sm_TM_RBM_predmat[i][0]);
        dy = (sm_TM_RBM_truemat[j][1] - sm_TM_RBM_predmat[i][1]);
        if ((dx<r) && (dx>(-r)) && (dy<r) && dy>(-r))
        {
            sole_matched = j;
            stars_in_range++;
        }
    }
    if (stars_in_range == 1)
    {
        matched++;
        sm_TM_RBM_matchmat[matched][0] = sm_TM_RBM_truemat[sole_matched][0];
        sm_TM_RBM_matchmat[matched][1] = sm_TM_RBM_truemat[sole_matched][1];
        sm_TM_RBM_matchmat[matched][2] = sm_TM_RBM_predmat[i][2];
    }
    //cout<<stars_in_range<<endl;
}

return matched + 1;
}

int sm_TM_SNT_match(double sm_TM_RBM_matchmat[][3], //array of stars matched in RBM
                    int L, //length of above array
                    double fe_output[][2], //entire output of centroids (x,y) from FE
                    int N, //length of above array
                    int sm_TM_SNT[5060][11], //stores star neighbourhood table
                    double sm_GD_SC[5060][4], //stores guide star catalogue
                    double sm_TM_CP_F, //focal length
                    double pixel_size, 
                    int Nth, //matching threshold
                    double error, //permitted error in angular distance
                    double sm_TM_SNT_output[][3] //stores SNT matched stars in the array
                    ) 
{
double fe_unmatched[N][2];
int num_unmatched = 0;

//to remove the already matched stars of FE output from consideration
int j=0;
int i=0;
for (j=0; j < N; j++)
{
    int already_matched = 0;

    for (i=0; i < L; i++)
    {
        if ((sm_TM_RBM_matchmat[i][0] == fe_output[j][0]) && (sm_TM_RBM_matchmat[i][1] == fe_output[j][1]))
        {
            already_matched = 1;
            break;
        }
    }
    if (already_matched==0)
    {
        fe_unmatched[num_unmatched][0] = fe_output[j][0];
        fe_unmatched[num_unmatched][1] = fe_output[j][1];
        num_unmatched++;
    }
}

//check the unmatched centroids against the matched stars' neighbours
//running = index of unmatched star that we want to match
//new_matched = number of matched stars, 
//              also the index to put in a newly matched star into array

int new_matched = 0, running = 0;
for (running = 0; running < num_unmatched; running++)
{
    int done = 0; //if FE unmatched star gets matched, done becomes 1
    {
        if (new_matched == Nth - L) break;
        int curr_ref_star = sm_TM_RBM_matchmat[i][2];
        if (curr_ref_star > 5060) continue;
        j = 0;

        while (sm_TM_SNT[curr_ref_star-1][j] != 0)
        {
            double x1 = sm_TM_RBM_matchmat[i][0];
            double y1 = sm_TM_RBM_matchmat[i][1];
            double x2 = fe_unmatched[running][0];
            double y2 = fe_unmatched[running][1];

            double q1 = sm_TM_id_angdist((curr_ref_star-1), (sm_TM_SNT[curr_ref_star-1][j]-1), sm_GD_SC);
            double q2 = sm_TM_centroid_angdist(x1, y1, x2, y2, sm_TM_CP_F, pixel_size);
            double q = q1 - q2;

            if ((q<error) && (q>(-error)))
            {
                int new_id = sm_TM_SNT[curr_ref_star-1][j];
                
                for (i=0; i<new_matched; i++)
                    if (sm_TM_SNT_output[i][2] == new_id) continue;

                sm_TM_SNT_output[new_matched][2] = sm_TM_SNT[curr_ref_star-1][j];
                sm_TM_SNT_output[new_matched][0] = fe_unmatched[running][0];
                sm_TM_SNT_output[new_matched][1] = fe_unmatched[running][1];

                new_matched++;
                done = 1;
            }
            j++;
            if (done==1) break; //break 

        }
        if (done==1) break;
    }
}

return new_matched;

}


//verify if stars are correctly matched by checking pairwise angular distance
//of identified star ID's and image stars.
int verify(double sm_TM_unverified_output[][3], int len, double sm_TM_verified_output[][3], double error)
{
    int count = 0;
    int i=0, j=0;
    double d1 = 0, d2 = 0, d=0;
    int okay = 1;
    for (i=0; i<len; i++)
    {
        okay = 1;
        for (j=i+1; j<len; j++)
        {
            d1 = sm_TM_id_angdist(sm_TM_unverified_output[i][0], sm_TM_unverified_output[j][0], sm_GC);
            d2 = sm_TM_centroid_angdist(sm_TM_unverified_output[i][1], sm_TM_unverified_output[i][2], sm_TM_unverified_output[j][1], sm_TM_unverified_output[j][2], 0.036, 1);
            d = d1 - d2;
            if ((d>error) || (d<-error)) okay = 0;
        }

        if (okay == 1)
        {
            sm_TM_verified_output[count][0] = sm_TM_unverified_output[i][0];
            sm_TM_verified_output[count][1] = sm_TM_unverified_output[i][1];
            sm_TM_verified_output[count][2] = sm_TM_unverified_output[i][2];
            count++;
        }
    }
    return count;
}

/*int lism(double unmatched[][3], int len_unmatched, double matched[][3])
{
    // double unmatched_input[len_unmatched][3];
    int i=0;
    // for (i=0; i<len_unmatched; i++)
    // {
    //     unmatched_input[i][0] = i;
    //     unmatched_input[i][1] = unmatched[i][0];
    //     unmatched_input[i][2] = unmatched[i][1];
    // }
    

    double body_vecs[8876][4];
    int sm_IS[8876][2];
    int lism_matched = 0;

    lism_matched = OILS(unmatched, len_unmatched, matched);
    i = 0;
    // for (i=0; i<lism_matched; i++)
    // {
    //     matched[i][0] = body_vecs[i][1]*FOCAL_LENGTH;
    //     matched[i][1] = body_vecs[i][2]*FOCAL_LENGTH;
    //     matched[i][2] = sm_IS[i][1];
    // }
    
    return lism_matched;
}
*/

int tracking(double prev1[][3], //previous frame 1 with x,y and star ID
             int len1, //length of above array
             double prev2[][3], //previous frame 2 with x,y and star ID
             int len2, //length of above array
             double trumat[][3], //new image frame with (x,y) of centroids in column 2 and 3
             int len3, //length of above array
             double sm_TM_unverified_output[][3] //TM matched final output
             )
{
int i = 0, temp = 0;
// for(i = 0; i < len3; i++){
// printf("%f %f %f\n", trumat[i][0], trumat[i][1], trumat[i][2]);
// }

double truemat[40][2];

i = 0;

double prev1_new[40][3], prev2_new[40][3];

for (i = 0; i <len1; i++){
    temp = prev1[i][0];
    prev1_new[i][0] = prev1[i][1];
    prev1_new[i][1] = prev1[i][2];
    prev1_new[i][2] = temp;
}

for (i = 0; i <len2; i++){
    temp = prev2[i][0];
    prev2_new[i][0] = prev2[i][1];
    prev2_new[i][1] = prev2[i][2];
    prev2_new[i][2] = temp;
}

i=0; temp = 0;

for (i = 0; i <len3; i++){
    truemat[i][0] = trumat[i][1];
    truemat[i][1] = trumat[i][2];
}

//for(i = 0; i < len3; i++){
// printf("%f %f\n", truemat[i][0], truemat[i][1]);
// }

{
double p1[40][3];
double p2[40][3]; 
int common_in_both_prev;

 //for(i = 0; i < 20; i++){
 //printf("%f %f %d\n", prev1[i][0], prev1[i][1], (int)prev1[i][2]);
 //}
 //printf("%d %d\n", len1, len2);
common_in_both_prev = sm_TM_CP_inputsort(prev1_new, len1, prev2_new, len2, p1, p2);

// for(i = 0; i < 20; i++){
//printf("%f %f %d\n", p1[i][0], p1[i][1], (int)p1[i][2]);
//}

// for(i = 0; i < 20; i++){
// printf("%f %f %d\n", p2[i][0], p2[i][1], (int)p2[i][2]);
// }


double pred[40][3];
sm_TM_CP(p1, p2, FOCAL_LENGTH, pred, common_in_both_prev, 1);

 //for(i = 0; i < 20; i++){
 //printf("%f %f %d\n", pred[i][0], pred[i][1], (int)pred[i][2]);
 //}

int num_RBM_matched = 0;
double matchmat[40][3];
num_RBM_matched = sm_TM_RBM(pred, common_in_both_prev, truemat, len3, r, matchmat);

printf("\n");



int i=0;
for (i=0; i<num_RBM_matched; i++)
{
    printf("%d ", (int)matchmat[i][2]);
    printf("%f ", matchmat[i][0]);
    printf("%f\n", matchmat[i][1]);
}

printf("\n\n");
printf("%d out of ", num_RBM_matched);
printf("%d", common_in_both_prev);
printf("\n\n");

double snt_out[10][3];
int num_SNT_matched = 0;
if (num_RBM_matched<Nth) 
    num_SNT_matched = sm_TM_SNT_match(matchmat, num_RBM_matched, truemat, len3, snt, sm_GC, FOCAL_LENGTH, 1, 20, e, snt_out);

for (i=0; i<num_SNT_matched; i++)
    {
    printf("%d ", (int)snt_out[i][2]);
    printf("%f ", snt_out[i][0]);
    printf("%f\n", snt_out[i][1]);
    }


    //cout<<snt_out[i][0]<<' '<<snt_out[i][1]<<' '<<snt_out[i][2]<<endl;

printf("\n");
printf("%d", num_SNT_matched);
printf("\n");

/*for (i=0; i<len2; i++)
{
prev1[i][0] = prev2[i][0];
prev1[i][1] = prev2[i][1];
prev1[i][2] = prev2[i][2];}
len1 = len2;
*/
//double sm_TM_unverified_output[40][3];

for (i=0; i<num_RBM_matched; i++)
{
sm_TM_unverified_output[i][0] = matchmat[i][2];
sm_TM_unverified_output[i][1] = matchmat[i][0];
sm_TM_unverified_output[i][2] = matchmat[i][1];
}

for (i = 0; i<num_SNT_matched; i++)
{
sm_TM_unverified_output[i+num_RBM_matched][0] = snt_out[i][2];
sm_TM_unverified_output[i+num_RBM_matched][1] = snt_out[i][0];
sm_TM_unverified_output[i+num_RBM_matched][2] = snt_out[i][1];
}

//double sm_TM_verified_output[40][3];
//int num_TM_verified = verify(sm_TM_unverified_output, num_RBM_matched + num_SNT_matched, sm_TM_verified_output, e);

int total = num_RBM_matched + num_SNT_matched;
for (i=0; i<total; i++)
{
    printf("%d %f %f\n", (int)sm_TM_unverified_output[i][0], sm_TM_unverified_output[i][1], sm_TM_unverified_output[i][2]);
}


printf("%d\n\n", total);

return total;
}
}






