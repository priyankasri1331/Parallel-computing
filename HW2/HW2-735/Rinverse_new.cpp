#include <iostream>
#include <cstdlib>
#include <cmath> 

#define N 2

using namespace std;


void display(double **A) 
{ 
    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
} 


class create_triu
{
    public:
    int n;
    //double **A;
    create_triu(int _n)
    {
        n = _n;
    }

    double** create_random_matrix()
    {
        double **A = new double*[n];
        for(int i = 0; i < n; i++)
            A[i] = new double[n];

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                A[i][j] = (double) rand() / RAND_MAX;
            }
        }
        return (double**) A;
    }

    double** absol(double **A)
    {  
        double x;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                x = A[i][j];
                A[i][j] = abs(x);
            }
        }
        return (double**) A;
    }

    double* sum(double **A)
    {
        double *A_arr = new double[n];
        double temp_sum = 0.0;
        
        for(int i = 0; i < n; i++)
        {
            temp_sum = 0.0;
            for(int j = 0; j < n; j++)
            {
                temp_sum += A[j][i]; 
            }
            A_arr[i] = temp_sum;
        }
        return (double*) A_arr;
    }

    double** diag(double *B)
    {
        double **A_diag = new double *[n];
        for(int i = 0; i < n; i++)
            A_diag[i] = new double[n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(i == j)
                {
                    A_diag[i][j] = B[i];
                }
                else
                {
                    A_diag[i][j] = 0.0;
                }
                
            }
        }
        return (double**) A_diag;
    }

    double** add(double **A, double **B)
    {
        double **A_add = new double *[n];
        for(int i = 0; i < n; i++)
            A_add[i] = new double[n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                A_add[i][j] = A[i][j] + B[i][j];
            }
        }
        return (double**) A_add;
    }

    double** triu(double **A)
    {
        double **A_triu = new double *[n];
        for(int i = 0; i < n; i++)
            A_triu[i] = new double[n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if((i == j) || (i < j))
                {
                    A_triu[i][j] = A[i][j];
                }
                else
                {
                    A_triu[i][j] = 0.0;
                }
                
            }
        }
        return (double**) A_triu;
    }

    double** getCofactor(double **A, int p, int q, int n)
    { 
        double **A_cofac = new double *[n];
        for(int i = 0; i < n; i++)
            A_cofac[i] = new double[n];
        //p = 0, q = 0;
        int a = 0, b = 0;
        for(int i = 0; i < n; i++)
        {
            if(i != p)
            {
                cout << "i = " << i << p << endl;
                b = 0;
                for(int j = 0; j < n; j++)
                {
                    if(j != q)
                    {
                        cout << "j = " << j << q << endl;
                        A_cofac[a][b] = A[i][j];
                        cout << A_cofac[a][b] << A[i][j] << endl;
                        b++;
                        
                    }
                }
                a++;
            }
        }
        display(A_cofac);
        return A_cofac;
    }

    double determinant(double **A, int n)
    {
        double D = 0;
        double **temp; //= new double *[n];
        //for(int i = 0; i < n; i++)
        //    temp[i] = new double[n];

        if(n == 1)
            return A[0][0];

        int sign = 1;
        //for(int k = 0; k < n; k++)
        //{
        for(int j = 0; j < n; j++)
            {
                cout << " ## j:" <<j<< "\n";  
                temp = getCofactor(A, 0, j, n);
                display(temp);
                D += (double) (sign * A[0][j] * determinant(temp, n - 1));

                sign = -sign;
            }
        //}

        return D;
    }

    //double adjoint()

    //double inverse()

    











};

/*double** compute_inverse(double **R)
{
    int n = N;
    double **Ri = new double *[n];
    for(int i = 0; i < n; i++)
        Ri[i] = new double[n];

    if (n == 1)
    {
        Ri[0][0] = 1/R[0][0];
    }
    else if(n < 16)
    {
        Ri = inv(R);
    }
    {
        int n1 = round(n/2);
        if (n1 == 2)
        {

        }

        Ri = R;
        Ri();

    }
    return (double**) Ri;
}*/


int main()
{
    int n = N;

    if(n == 0)
    {
        cout << "n should be greater than 0";
        exit(0);
    }
    double **A;
    double *B;
    double **C;
    double **D;
    double **E;
    double **F;
    double G;

    create_triu obj1(n);
    A = obj1.create_random_matrix();
    A = obj1.absol(A);
    B = obj1.sum(A);
    C = obj1.diag(B);
    D = obj1.add(A,C);
    //E = obj1.triu(D);


    //E = obj1.getCofactor(A, 0, 2, 2);



    G = obj1.determinant(E,n); 
    //E has the triangular matrix

    display(A);
    cout << endl;
    display(E);
    cout << endl;
    cout << G << endl;

    free(A);
    free(B);
    free(C);
    free(D);
    free(E);
    free(F);
    

    //F = compute_inverse(E);
    return 0;
}
/*
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cout << A[i][j] << endl;
        }
    }

    for(int i = 0; i < n; i++)
    {
        cout << B[i] << endl;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << C[i][j] << endl;
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << D[i][j] << endl;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << E[i][j] << endl;
    }
*/