#include <math.h>
#include "iostream"
#include <omp.h>
#include<stdlib.h>
#include<iomanip>
#include<math.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;

#define PI 3.14

void eye(vector<vector<double> > *ide_mat, int n)
{
    #pragma omp parallel for
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

void comp_inverse(vector<vector<double> > &B)
{
    
    int n = B.size();
    //vector<double> x(n);
    vector<vector<double> > A (n, vector<double>(2 * n,0));
    double ratio;
    int i,j,k;
    
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++)
    {
	//#pragma omp parallel for
        for(int j  = 0; j < n; j++)
        {
            A[i][j] = B[i][j];
        }
    }

    /* Augmenting Identity Matrix of Order n */
    #pragma omp parallel for collapse(2)
    for(i=0;i<n;i++)
    {
        //#pragma omp parallel for
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
    
    //#pragma omp parallel for collapse(2)
    //#pragma parallel task
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
		//#pragma omp parallel for
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
    //#pragma omp parallel for
    for(int i = 0; i < n ; i++)
    {
        //#pragma omp parallel for
        for(int j = n; j < n *2; j++)
        {
            B[i][j - n] = A[i][j];
        }
    }

    return;
}
void display_vec(vector<double> &A)
{
    for(int i = 0; i < A.size(); i++)
    {
        cout << A[i] << "\t";
    }
    cout << endl;
}

void display_vec_int(vector<int> &A)
{
    for(int i = 0; i < A.size(); i++)
    {
        cout << A[i] << "\t";
    }
    cout << endl;
}


void display( vector<vector<double> > &A) 
{

    for (int i=0; i<A.size(); i++) 
    { 
        for (int j=0; j<A[0].size(); j++)
        {
            //cout << "I'm here";
            cout << A[i][j] << " "; 

        } 

        cout << endl; 
    } 
} 

vector<vector<double> > zeros(int n_rows, int n_cols)
{
    
    vector<vector<double> > zeros_matrix(n_rows, vector<double>(n_cols,0));
    
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n_rows; i++)
    {
        //#pragma omp parallel for
        for(int j = 0; j < n_cols; j++)
        {
            zeros_matrix[i][j] = 0;
        }

    }
    return zeros_matrix;

}

vector<vector<double> > ones(int n_rows, int n_cols)
{
    
    vector<vector<double> > ones_matrix(n_rows, vector<double>(n_cols,0));
    
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n_rows; i++)
    {
        //#pragma omp parallel for
        for(int j = 0; j < n_cols; j++)
        {
            ones_matrix[i][j] = 1;
        }

    }
    return ones_matrix;

}

vector<vector<double> > matrix_mul( vector<vector<double> > &A, vector<vector<double> > &B)
{
    int A_row = A.size();
    int A_col = A[0].size();
    int B_row = B.size();
    int B_col = B[0].size();

    vector<vector<double> > temp(A_row, std::vector<double>(B_col, 0));

    #pragma omp parallel for
    for(int i = 0; i < A_row; i++)
    {
        for(int j = 0; j < B_col; j++)
        {
            //#pragma omp parallel for
            for(int k = 0; k < A_col; k++)
            {
                //cout << i << "____" << j << "___________" << k << endl;
                temp[i][j] += A[i][k]*B[k][j];
            }
        }
    }

    return temp;
}

void mat_trans_pose(vector<vector<double> > XY, vector<vector<double> > *XYt)
{
    int m = XY.size();
    int n = XY[0].size();
    
    //#pragma omp parallel for
    for(int i = 0; i < m; i++)
    {
        for(int j =0;j < n; j++)
        {
            (*XYt)[j][i] = XY[i][j];

        }

    }
    return;

}

void init_mxm(int m, vector<vector<double> > *XY, double h)
{

    int idx = 0;
    for(int i = 1; i<= m; i++)
    {
        //#pragma omp parallel for
	for(int j = 1; j <= m; j++)
        {           
            (*XY)[idx][0] = (double)i * h;
            (*XY)[idx][1] = (double)j * h;
            idx = idx + 1;
        }
        
    }

    return;
}

void create_random_vec(vector<vector<double> > &A)
{
    int n_rows = A.size();
    int n_cols = A[0].size();
    //#pragma omp parallel for
    for(int i = 0; i < n_rows; i++)
    {
	//#pragma omp parallel for
        for(int j = 0; j < n_cols; j++)
        {
            A[i][j] = (double) rand()/RAND_MAX;
        }
    }
    return;
}


vector<vector<double> > kernel(vector<vector<double> > *X, vector<vector<double> > *Y, vector<vector<double> > *L)
{
    
    int X_size = (*X).size();
    int Y_size = (*Y).size();
    
    vector<vector<double> > X_first_col(X_size, vector<double>(1,0));
    vector<vector<double> > X_second_col(X_size, vector<double>(1,0));
    vector<vector<double> > Y_first_col(Y_size, vector<double>(1,0));
    vector<vector<double> > Y_second_col(Y_size, vector<double>(1,0));

    vector<vector<double> > temp_1(1, vector<double>(Y_size,0));
    vector<vector<double> > temp_3(X_size, vector<double>(1,0));
    vector<vector<double> > temp_2(X_size, vector<double>(Y_size,0));
    vector<vector<double> > temp_4(X_size, vector<double>(Y_size,0));
    vector<vector<double> > temp_6(X_size, vector<double>(Y_size,0));
    vector<vector<double> > temp_8(X_size, vector<double>(Y_size,0));
    vector<vector<double> > Y_transpose_1(Y_first_col[0].size(), vector<double>(Y_first_col.size(),0));
    vector<vector<double> > Y_transpose_2(Y_first_col[0].size(), vector<double>(Y_first_col.size(),0));
    vector<vector<double> > term_1(X_size, vector<double>(Y_size,0));
    vector<vector<double> > term_2(X_size, vector<double>(Y_size,0));



    #pragma omp parallel for
    for(int i = 0; i < X_size; i++)
    {
        X_first_col[i][0] = (*X)[i][0];
    }
    
    #pragma omp parallel for
    for(int i = 0; i < X_size; i++)
    {
        X_second_col[i][0] = (*X)[i][1];
    }

    #pragma omp parallel for
    for(int i = 0; i < Y_size; i++)
    {
        Y_first_col[i][0] = (*Y)[i][0];
    }
    
    #pragma omp parallel for
    for(int i = 0; i < Y_size; i++)
    {
        Y_second_col[i][0] = (*Y)[i][1];
    }

    mat_trans_pose(Y_first_col,&Y_transpose_1);

    mat_trans_pose(Y_second_col,&Y_transpose_2);
    

    temp_1 = ones(1,Y_size);
    
    temp_2 = matrix_mul(X_first_col, temp_1);

 
    temp_3 = ones(X_size,1);  

    temp_4 = matrix_mul(temp_3, Y_transpose_1);
    
    temp_6 = matrix_mul(X_second_col, temp_1);

    temp_8 = matrix_mul(temp_3, Y_transpose_2);


    #pragma omp parallel for collapse(2)
    for(int i = 0; i < term_1.size(); i++)
    {
        for(int j =0;j < term_1[0].size(); j++)
        {
            term_1[i][j] = (temp_2[i][j] - temp_4[i][j]);
            term_1[i][j] = term_1[i][j] * term_1[i][j] /((*L)[0][0] * (*L)[0][0]);

            term_2[i][j] = (temp_6[i][j] - temp_8[i][j]);
            term_2[i][j] = term_2[i][j] * term_2[i][j] /((*L)[1][0] * (*L)[1][0]);

            term_1[i][j] = 1/sqrt(2*PI) * exp(-(1.0/2.0)*(term_1[i][j] + term_2[i][j]));



        }

    }  

    return term_1;

}

void init_obs_data_f(vector<vector<double> > *XY, vector<vector<double> > *f, int n, int m)
{
    
    vector<vector<double>> rand_matrix(n, vector<double>(1,0));
    vector<vector<double>> kernel_sp(n, vector<double>(1,0));
    vector<vector<double>> A(1,vector<double>(2,0));;
    vector<vector<double>> B(2,vector<double>(1,0));

    vector<vector<double>> C(2,vector<double>(1,0));
    vector<vector<double>> D(n,vector<double>(1,0));

    C[0][0] = 0.2;
    C[1][0] = 0.1;

    D = matrix_mul(*XY, C);

    //cout << "_________" << endl;
    
    
    create_random_vec(rand_matrix);

    for(int i = 0; i < 2; i++)
    {
        A[0][i] = 0.25;
    }
    

    for(int i = 0; i < 2; i++)
    {
        B[i][0] = 2/m;
    }

    kernel_sp = kernel(XY, &A, &B);
   
    #pragma omp parallel for 
    for(int i = 0; i < n; i++)
    {
        (*f)[i][0] = (0.02*(rand_matrix[i][0] - 0.5)) + kernel_sp[i][0] + D[i][0];

    }

    return;
}

vector<double> LParameter(double start, double end, double stride, int m)
{

    std::vector<double> Lparam;
    double sum = start;

    while(sum < end)
    {
        Lparam.push_back(sum);
        sum += stride;

    }

    for(int i = 0; i < Lparam.size(); i++)
    {
        Lparam[i] = Lparam[i]/m;
    }

    return Lparam;

}

void LUfac(vector<vector<double> > *L, vector<vector<double> > *U, int n)
{

    unsigned long int sum = 0;
    //auto start_time = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n; i++)
    {
        #pragma omp parallel for 
	for(int row = (i + 1); row < n; row++)
        {
            //The factor that the above row elements are to be multiplied with
            double fac = (*U)[row][i] / (*U)[i][i];
            for(int col = 0; col < n; col++)
            {              
                sum += 1;
 		(*U)[row][col] = (*U)[row][col] - ((fac)*((*U)[i][col]));
                
            }
            (*L)[row][i] = fac;
        }

    }
    //auto stop_time = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> dur_ms = stop_time - start_time;
    //std::cout << "LU exec time: " << dur_ms.count() << "ms" << std::endl;
    //cout << "no of operations = " << sum << endl;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
	    sum += 1;
            if(i == j)
                (*L)[i][j] = 1;
        }

    }

    return;
}

vector<vector<double>> compute_pred_val(vector<vector<double>> *U, vector<vector<double>> *L,  vector<vector<double>> *f,  vector<vector<double> > *k, vector<int> *itrain, vector<int> *itest)
{

    int n = (*U).size();
    vector<vector<double> > temp1((*L).size(), vector<double>(1, 0));
    vector<vector<double> > temp2((*U).size(), vector<double>(1,0));
    vector<vector<double> > k_dash((*k)[0].size(), vector<double>((*k).size(),0));
    vector<vector<double> > ftest(sqrt(n),vector<double>(1,0));
    vector<vector<double> > f_train((*itrain).size(),vector<double>(1,0));

    comp_inverse(*U);
    
    comp_inverse(*L);

    mat_trans_pose(*k, &k_dash);

    #pragma omp parallel for

    for(int i = 0; i < (*itrain).size(); i++)
    {
        
        f_train[i][0] = (*f)[(*itrain)[i]][0];
    }


    //Forward substitution
    temp1 = matrix_mul(*L, f_train);


    temp2 = matrix_mul(*U, temp1);


    ftest = matrix_mul(k_dash, temp2);

    return ftest;
}

vector<vector<double> > GPR(vector<vector<double> > *XY, vector<vector<double> > *f,vector<int> *itest, vector<int> *itrain, double t, vector<vector<double> > *l)
{

    int n = (*XY).size();

    vector<vector<double>> K0(n, vector<double>(n,0));

    K0 = kernel(XY, XY, l);

    int ntest = (*itest).size();
    int ntrain = (*itrain).size();

    vector<vector<double>> K(ntrain, vector<double>(ntrain,0));
    
	
    #pragma omp parallel for
    for(int i=0; i<ntrain; i++)
    {
        for(int j=0; j<ntrain; j++)
        {
            K[i][j] = K0[(*itrain)[i]][(*itrain)[j]];
        }
    }


    vector<vector<double>> L(ntrain, vector<double>(ntrain,0));
    vector<vector<double>> U(ntrain, vector<double>(ntrain,0));


    eye(&U, ntrain);
    
 
    #pragma omp parallel for 
    for(int i = 0; i < ntrain; i++)
    {
        
        for(int j = 0; j < ntrain; j++)
        {
            U[i][j] = t * U[i][j] + K[i][j];
        }
    }

    LUfac(&L, &U, ntrain);


    vector<vector<double>> k(ntrain, vector<double>(ntest,0));
    #pragma omp parallel for
    for(int i=0; i<ntrain; i++)
    {
        for(int j=0; j<ntest; j++)
        {
            k[i][j] = K0[(*itrain)[i]][(*itest)[j]];
        }
    }


    vector<vector<double>> ftest(ntest, vector<double>(1,0));

    ftest = compute_pred_val(&U, &L, f, &k, itrain, itest);


    return ftest;
}

int main()
{

    int m = 32;
    int n = m*m;
    double h = (double) 1/(m + 1);
    int ntest;
    int ntrain;


    
    //double t = 0.01;
    double fstar;

    vector<vector<double>> XY(n, vector<double>(2,0));

    vector<vector<double>> f(n, vector<double>(1,0));

    vector<double> rperm;

    XY = zeros(n, 2);

    auto start_time = std::chrono::high_resolution_clock::now();

    init_mxm(m, &XY, h);

    init_obs_data_f(&XY, &f, n, m);


    if(n >= 10)
    {
        ntest = round(0.1*n);
    }
    else
    {
        ntest = 1;
    }
    
    ntrain = n - ntest;

    for(int i = 0; i < n; i++)
    {
        rperm.push_back(i);
    }
    //std::random_device rd;
    //std::mt19937 g(rd());
    srand(0);
    random_shuffle(rperm.begin(), rperm.end()); 

    vector<int> itest(ntest,0);
    vector<int> itrain(ntrain,0);

    for(int i = 0; i < ntest; i++)
    {
        itest[i] = rperm[i];
    }

    for(int i = 0; i < ntrain; i++)
    {
        itrain[i] = rperm[i + ntest];
    }

    vector<double> Lparam = LParameter(0.25, 10, 0.5, m);

    vector<vector<double>> L_param_vec(2, vector<double>(1,0));

    vector<vector<double> > ftest(itest.size(),vector<double>(1,0)); 


    vector<vector<double> > MSE(Lparam.size(),vector<double>(Lparam.size(),0));
    vector<vector<double> > err(itest.size(),vector<double>(1,0)); 
    vector<vector<double> > err_transpose(1,vector<double>(itest.size(),0));
    vector<vector<double> > temp_vec(1,vector<double>(1,0));

    //Lparam.size()


	double min = 100000000;
	double lparam_1 = 0;
	double lparam_2 = 0;

    int i = 0;
    int j = 0;
  
    //
    //pragma omp parallel for collapse(2) shared(min,ftest,MSE,err,L_param_vec,lparam_1,lparam_2) private(i,j)
    //#pragma omp parallel for  ordered schedule(static,16) 
    for(i = 0; i < Lparam.size(); i++)
    {
        //#pragma omp parallel for
        for(j = 0; j < Lparam.size(); j++)
        {
	    //#pragma omp critical
	    //#pragma omp ordered
            L_param_vec[0][0] = Lparam[i];
            //#pragma omp critical
            //#pragma omp ordered
            L_param_vec[1][0] = Lparam[j];
            //#pragma omp critical
            //#pragma omp task
            //#pragma omp critical
           //#pragma omp ordered
            ftest = GPR(&XY, &f, &itest, &itrain, 0.5, &L_param_vec);
            
            //#pragma omp simd            
            for(int k = 0; k < itest.size(); k++)
            {
                err[k][0] = f[itest[k]][0] - ftest[k][0];
            }

	    //#pragma omp ordered
            mat_trans_pose(err, &err_transpose);
            //#pragma omp ordered
            temp_vec = matrix_mul(err, err_transpose);
            MSE[i][j] = temp_vec[0][0] / (double) ntest;
	    //#pragma omp ordered
	    if(min > MSE[i][j])
	    {
		min = MSE[i][j];
		lparam_1 = Lparam[i];
		lparam_2 = Lparam[j];
            }
            	//cout << "Lparam(l1) = " << Lparam[i] << "\t" << "Lparam(l2) = " << Lparam[j] << "\t" << "MSE = " << MSE[i][j] << endl;

        }
        
    }
    
    cout << "Lparam(l1) = " << lparam_1 << "\t" << "Lparam(l2) = " << lparam_2 << "\t" << "MSE_min = " << min << endl;

    auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur_ms = stop_time - start_time;
    std::cout << "GPR exec time: " << dur_ms.count() << "ms" << std::endl;

    return 0;
}
