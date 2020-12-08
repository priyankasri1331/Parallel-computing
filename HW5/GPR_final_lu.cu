#include "iostream"
#include <omp.h>
#include<stdlib.h>
#include<iomanip>
#include<math.h>
#include <iostream>
#include <cstdlib>
#include <cuda.h>
#include <cmath> 
#include <vector>
#include <chrono>

using namespace std;

int num_threads;

//returns a zero matrix of requested size
void comp_inverse(vector<vector<double> > &B)
{
    
    int n = B.size();
    //vector<double> x(n);
    vector<vector<double> > A (n, vector<double>(2 * n,0));
    double ratio;
    int i,j,k;
    

    for(int i = 0; i < n; i++)
    {

        for(int j  = 0; j < n; j++)
        {
            A[i][j] = B[i][j];
        }
    }

    /* Augmenting Identity Matrix of Order n */

    for(i=0;i<n;i++)
    {

        for(j=0;j<n;j++)
        {
            if(i==j)
            {
                A[i][j+n] = 1;
            }
            else
            {
                A[i][j+n] = 0;
            }
        }
    }

    /* Applying Gauss Jordan Elimination */
    
    for(i=0;i<n;i++)
    {
        if(A[i][i] == 0.0)
        {
            cout<<"Mathematical Error!";
            exit(0);
        }
        for(j=0;j<n;j++)
        {
            if(i!=j)
            {
                ratio = A[j][i]/A[i][i];
                for(k=1;k<2*n;k++)
                {
                    A[j][k] = A[j][k] - ratio*A[i][k];
                }
            }
        }
    }
    /* Row Operation to Make Principal Diagonal to 1 */
    
    for(i=0;i<n;i++)
    {
        for(j=n;j<2*n;j++)
        {
        A[i][j] = A[i][j]/A[i][i];
        }
    }

    //copy from the augmented matrix to matrix A;

    for(int i = 0; i < n ; i++)
    {

        for(int j = n; j < n *2; j++)
        {
            B[i][j - n] = A[i][j];
        }
    }

    return;
}
//Matrix multiplication
vector<vector<double> > matrix_mul( vector<vector<double> > &A, vector<vector<double> > &B)
{
    int A_row = A.size();
    int A_col = A[0].size();
    int B_row = B.size();
    int B_col = B[0].size();

    vector<vector<double> > temp(A_row, std::vector<double>(B_col, 0));

    for(int i = 0; i < A_row; i++)
    {
        for(int j = 0; j < B_col; j++)
        {
            for(int k = 0; k < A_col; k++)
            {
                //cout << i << "____" << j << "___________" << k << endl;
                temp[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    //display(temp);
    //cout << "_________________";
    return temp;
}
//returns a matrix with 0s
vector<vector<double> > zeros(int n_rows, int n_cols)
{
    
    vector<vector<double> > zeros_matrix(n_rows, vector<double>(n_cols,0));
    
    for(int i = 0; i < n_rows; i++)
    {
        for(int j = 0; j < n_cols; j++)
        {
            zeros_matrix[i][j] = 0;
        }

    }
    return zeros_matrix;

}

void display_mat(vector<vector<double>> &A , int n)
{
    for(int i = 0; i < n; i++)
    {
	for(int j = 0; j < n;j++)
    	{
        cout << A[i][j] << "\t";
    	}
    cout << endl;
	}
return;
}
// Matrix with random values in it
void create_random_vec(vector<vector<double> > &A)
{
    int n_rows = A.size();
    int n_cols = A[0].size();

    for(int i = 0; i < n_rows; i++)
    {

        for(int j = 0; j < n_cols; j++)
        {
            A[i][j] = (double) rand()/RAND_MAX;
        }
    }
    return;
}
//initialize mxm grid of points
void init_mxm(int m, vector<vector<double> > *XY, double h)
{

    int idx = 0;
    for(int i = 1; i<= m; i++)
    {

	for(int j = 1; j <= m; j++)
        {           
            (*XY)[idx][0] = (double)i * h;
            (*XY)[idx][1] = (double)j * h;
            idx = idx + 1;
        }
        
    }

    return;
}

void mat_trans_pose(vector<vector<double> > XY, vector<vector<double> > *XYt)
{
    int m = XY.size();
    int n = XY[0].size();
    

    for(int i = 0; i < m; i++)
    {
        for(int j =0;j < n; j++)
        {
            (*XYt)[j][i] = XY[i][j];

        }

    }
    return;

}
//Initialize observed data vector f

void init_obs_data_f(vector<vector<double> > *XY, vector<vector<double> > *f, int n)
{
    
    vector<vector<double>> rand_matrix(n, vector<double>(1,0));
    
    create_random_vec(rand_matrix);

    int n_rows = (*XY).size();
    int n_cols = (*XY)[0].size();

    vector<vector<double>> temp(n_rows, vector<double>(n_cols,0));
    vector<vector<double>> temp_col(n, vector<double>(1,0));
   

    for(int i = 0; i < n; i++)
    {
        (*f)[i][0] = 0.1*(rand_matrix[i][0] - 0.5);

    }

 
    for(int i = 0; i < n_rows; i++)
    {

        for(int j = 0; j < n_cols; j++)
        {
            temp[i][j] = (*XY)[i][j] - 0.5;
        }

    }


    double sum;

    for(int i = 0; i < n_rows; i++)
    {
        sum = 0;
        for(int j = 0; j < n_cols; j++)
        {
            sum += (temp[i][j]) * (temp[i][j]);

        }
        temp_col[i][0] = sum;
    }




    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < 1; j++)
        {
            (*f)[i][j] = (*f)[i][j] + 1.0 - temp_col[i][j];
        }
    }


    return;
}

//initialize K matrix

void init_K(vector<vector<double> > *XY, vector<vector<double> > *K, int n)
{
    (*K) = zeros(n , n);

    vector<double> dist(1,(*XY)[0].size());
    double sum;

    for(int a = 0; a < n; a++)
    {

	for(int b = 0; b < n; b++)
        {
 
            for(int c = 0; c < ((*XY)[0].size()); c++)
            {
                dist[c] =  ((*XY)[a][c]) - ((*XY)[b][c]);

            }
        
            sum  = 0;
            for(int d = 0; d < ((*XY)[0].size()); d++)
            {
                sum += (dist[d]) * (dist[d]);
            }

            (*K)[a][b] = exp (-(sum));

        }
    }

    return;
}

//Create an identity matrix
void eye(vector<vector<double> > *ide_mat, int n)
{

    for(int i  = 0; i < n; i++)
    {
        
        for(int j = 0; j < n; j++)
        {
            if(i == j)
            {
                (*ide_mat)[i][j] = 1;
            }
            else
            {
                (*ide_mat)[i][j] = 0;
            }
            
        }
    }

    return;
}

__global__ void L_gpu(double * L, double * U, int n, int index){
    int thread_id = (blockDim.x*blockIdx.x)+threadIdx.x;
    if(thread_id<(n-(index+1))) {
        L[((index+1+thread_id)*n)+index] = U[((index+1+thread_id)*n)+index]/U[(index*n)+index];

    }
    //__syncthreads();
}

__global__ void U_gpu(double * L, double * U, int n, int index, int row_index){
    int thread_id = (blockDim.x*blockIdx.x)+threadIdx.x;
    if(thread_id<(n-index)) {
        //for(int col=index; col<dimension; col++){
        U[(row_index*n)+index+thread_id] -= L[(row_index*n)+index]*U[(index*n)+index+thread_id];
        //}
    }
    //__syncthreads();
}


//Performs LU factorization
void LUfac(vector<vector<double> > *L, vector<vector<double> > *U, int n)
{

    double *L_new = (double*) malloc(n*n*sizeof(double));
    double *U_new = (double*) malloc(n*n*sizeof(double));

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            L_new[(i*n) + j] = (*L)[i][j];
            U_new[(i*n) + j] = (*U)[i][j];
        }
    } 

    double *gpu_U;
    double *gpu_L;
    
    int xxx = cudaMalloc((void **)&gpu_U, n*n*sizeof(double));
    int yyy = cudaMalloc((void **)&gpu_L, n*n*sizeof(double));


    xxx = cudaMemcpy(gpu_U, U_new, n*n*sizeof(double), cudaMemcpyHostToDevice);
    yyy = cudaMemcpy(gpu_L, L_new, n*n*sizeof(double), cudaMemcpyHostToDevice);
    
    //cout << xxx << "   " << yyy << "memcpy" << endl;

    //dim3 threadsPerBlock(2, 2);
    //dim3 numBlocks(n / threadsPerBlock.x + 1, n / threadsPerBlock.y + 1);
    dim3 threadsPerBlock(32);
    dim3 numBlocks(n / threadsPerBlock.x + 1);

	int sum = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n; i++)
    {
       L_gpu<<<numBlocks,threadsPerBlock>>>(gpu_L, gpu_U, n, i);
        for(int row = (i + 1); row < n; row++)
            {
                U_gpu<<<numBlocks,threadsPerBlock>>>(gpu_L, gpu_U, n, i, row);
                //*******************************************************
                //Enable below code segment if GPU functions are disabled
                //*********************************************************
                //The factor that the above row elements are to be multiplied with
                /*double fac = (*U)[row][i] / (*U)[i][i];
                for(int col = 0; col < n; col++)
                {              
                    sum += 1;
                    (*U)[row][col] = (*U)[row][col] - ((fac)*((*U)[i][col]));
                }
                (*L)[row][i] = fac;*/
            }

    }
    auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur_ms = stop_time - start_time;
    std::cout << "LU exec time GPU: " << dur_ms.count() << "ms" << std::endl;



    double *outU;
    double *outL;
    
    outU = (double*)malloc(sizeof(double)*n*n);
    outL = (double*)malloc(sizeof(double)*n*n);
    xxx = cudaMemcpy( outU, gpu_U, n*n*sizeof(double),cudaMemcpyDeviceToHost);
    yyy = cudaMemcpy( outL, gpu_L, n*n*sizeof(double),cudaMemcpyDeviceToHost);

    //cout << xxx << "   memcpy2" << endl;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            (*U)[i][j] = outU[(i*n)+j];
            (*L)[i][j] = outL[(i*n)+j];
        }
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(i == j)
                (*L)[i][j] = 1;
        }

    }
    cudaFree(gpu_U);
    cudaFree(gpu_L);

    return;
}


void init_k(vector<double> *k, vector<double> *rstar,  vector<vector<double>> *XY, int n)
{
    for(int i =0; i < n; i++)
    {
        (*k)[i] = 0.0;

    }
    vector<double> d(n,1);

    for(int i = 0; i < n; i++)
    {
        
        for(int j = 0; j < (*rstar).size(); j++)
        {
            d[j] = (*rstar)[j] - (*XY)[i][j];


        }

        double sum = 0;
        

        for(int k = 0; k < (*rstar).size(); k++)
        {
            sum += (d[k])*(d[k]);
        }
        (*k)[i] = exp(-(sum));


    }

    return;

}

double compute_pred_val(vector<vector<double>> *U, vector<vector<double>> *L,  vector<vector<double>> *f,  vector<double> *k)
{

    int n = (*U).size();
    vector<vector<double> > temp1((*L).size(), vector<double>(1, 0));
    vector<vector<double> > temp2((*U).size(), vector<double>(1,0));
    vector<vector<double> > k_new(1, vector<double>((*U).size(),0));
    vector<vector<double> > fstar(1,vector<double>(1,0));

    //display(*U);
    comp_inverse(*U);
    //cout << "inv_U" << endl;
    //display(*U);
    
    comp_inverse(*L);
    //cout << "inv_L" << endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < (*U).size(); i++)
    {
        k_new[0][i] = (*k)[i];

    }

    //cout << "k_new" << endl;
 
    //Forward substitution
    temp1 = matrix_mul(*L, *f);

    //cout << "temp1" << endl;
    //Back substitution
    temp2 = matrix_mul(*U, temp1);
    //cout << "temp2" << endl;

    fstar = matrix_mul(k_new, temp2);
    //cout << "fstar" << endl;
    auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur_s = stop_time - start_time;
    std::cout << "Solver routine time GPU: " << dur_s.count() << "ms" << std::endl;
    cout << endl;
    return fstar[0][0];
}

__global__ void tuk(double *U, double *K, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
	//printf("+++++++++++++++++++++++++++");
    	if (i < n && j < n)
	{
		//printf("Test\n");
		//printf("%d %d\n", U[i][j], K[i][j]);
        	U[(i*n)+j]= 0.01 * U[(i*n)+j] + K[(i*n)+j];
	}

	//printf("%d %d\n", U[i][j], K[i][j]);
}



double GPR(int m, vector<double> *rstar)
{

    //initialize mxm grid of points

    int n = m*m;
    double h = (double) 1/(m + 1);
    //double t = 0.1;
    double fstar;

    //Create space for all matrices
    vector<vector<double>> XY(n, vector<double>(2,0));
    vector<vector<double>> f(n, vector<double>(1,0));
    vector<vector<double>> K(n, vector<double>(n,0));
    vector<vector<double>> L(n, vector<double>(n,0));
    vector<vector<double>> U(n, vector<double>(n,0));
    vector<vector<double>> LU(n, vector<double>(n,0));
    vector<vector<double>> k_check(n, vector<double>(n,0));
    vector<double> k(n,1);


    XY = zeros(n, 2);

    init_mxm(m, &XY, h);
    //cout << "init_mxm" << endl;



    init_obs_data_f(&XY, &f, n);
    //cout << "init_obs_data_f" << endl;

    init_K(&XY, &K, n);

    //cout << "init_K" << endl;
    eye(&U, n);
    //cout << "eye" << endl;

    //cout << "LU" << endl;
    double *gpu_U;
    double *gpu_K;
    
    int xxx = cudaMalloc((void **)&gpu_U, n*n*sizeof(double));
    int yyy = cudaMalloc((void **)&gpu_K, n*n*sizeof(double));

   // if(xxx == cudaErrorMemoryAllocation) cout << "Success 1" << endl;

    //cout << xxx << "   xxx" << endl;

    double *U_new = (double*) malloc(n*n*sizeof(double));
    double *K_new = (double*) malloc(n*n*sizeof(double));

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            U_new[(i*n) + j] = U[i][j];
            K_new[(i*n) + j] = K[i][j];
        }
    } 


    xxx = cudaMemcpy(gpu_U, U_new, n*n*sizeof(double), cudaMemcpyHostToDevice);
    yyy = cudaMemcpy(gpu_K, K_new, n*n*sizeof(double), cudaMemcpyHostToDevice);
    
    //cout << xxx << "   " << yyy << "memcpy" << endl;

    dim3 threadsPerBlock(16,16);
    dim3 numBlocks(n / threadsPerBlock.x + 1, n / threadsPerBlock.y + 1);

    //cout << "Started !!" << n / threadsPerBlock.x << endl;
    tuk<<<numBlocks, threadsPerBlock>>>(gpu_U, gpu_K, n);
    //cout << "Endl  !!" << endl;

    double *outU;
    outU = (double*)malloc(sizeof(double)*n*n);
    xxx = cudaMemcpy( outU, gpu_U, n*n*sizeof(double),cudaMemcpyDeviceToHost);

    //cout << xxx << "   memcpy2" << endl;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            U[i][j] = outU[(i*n)+j];
        }
    }

	//*******************************************************************
        //Enable below segment if GPU is disabled

	//**********************************************************************

	/*for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			U[i][j] = 0.01 * U[i][j] + K[i][j];
		}			
	}*/




    cudaFree(gpu_U);
    cudaFree(gpu_K);
	
	//display_mat(U,n);
	//cout << "___________________" << endl;
//display_mat(U, n);   

    LUfac(&L, &U, n);

    init_k(&k, &(*rstar),&XY,n);

    fstar = compute_pred_val(&U, &L, &f, &k);

    return fstar;

}



int main(int argc, char *argv[])
{
    vector<double> rstar(2,1);

    cudaSetDevice(1);

    int m = std::atoi(argv[1]);
    rstar[0] = std::atof(argv[2]);
    rstar[1] = std::atof(argv[3]);

    auto start_time = std::chrono::high_resolution_clock::now();
    double fstar = GPR(m, &rstar);
    auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur_s = stop_time - start_time;
    std::cout << "Total time GPU: " << dur_s.count()/1000 << "s" << std::endl;
    cout << endl;
    cout << "fstar=" << fstar;
    cout << endl;

    return 0;
}
