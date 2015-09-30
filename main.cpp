#include <iostream>
#include <cmath>
//#include "armadillo"
#include <lib.h>
#include <algorithm>

using namespace std;
//using namespace arma;
#include "time.h"



//void Jacobi_rotation();
void Jacobi_rotations(double **,int); //matrix double ** A, dimention of the matrix int n
void matr_rotation(double **,double ,double,int, int, int);
void print_out_matrix(double **,int);
void find_min_eigenvalues(double **,int);

void matrix_construct(double,double,int,double**);

int main()
{
    clock_t start, finish; // declaring start and final time of the whole program
    start =clock();
    int n = 3;
    double **A;
    A=new double*[n];

    for(int i=0;i<n;i++){
        A[i]=new double[n];
        for(int j=0;j<n;j++){
            A[i][j]=0;//n*i+j;

        }

    }
    A[0][0]=1;A[0][1]=2;A[0][2]=3;
    A[1][0]=2;A[1][1]=2;A[1][2]=4;
    A[2][0]=3;A[2][1]=4;A[2][2]=3;
    //printing out the matrix
    print_out_matrix(A,n);
    //a)finding the max value of matrix:
//    Jacobi_rotations(A,n);

    //constructing a matrix to solve Schrodinger with potential V:
    cout<< endl<<"Quantum eigenv"<<endl;
    double ro_min=0; //domain : from ro_min to ro_max
    double ro_max=3.5;
    n=800;     // the domain will be split up in n points

    double **B;
    B=new double*[n];
    matrix_construct(ro_min,ro_max,n,B);


 //   print_out_matrix(B,n);

    Jacobi_rotations(B,n);
    cout <<"very excite!: "<<endl;

    find_min_eigenvalues(B,n);
    finish=clock();
    cout<< endl<< "time that it took to run the program: "<< ((finish-start)/CLOCKS_PER_SEC)<<endl;
    cout <<"n_step = "<<n<<endl;
    return 0;
}

void Jacobi_rotations(double ** A_, int n_){
    cout<< "initializing Jacobi_rotations " << n_ <<" points"<<endl;
    int amount_of_trans=0;
    double epsilon = pow(10,-8); //epsilon is going to be our boundary for 0
    int i_max;int j_max;
    double max;
    while(1){
    // PART 1 :
    // a for loop that finds max non diag value:
    max=0;
    for (int i=0;i<n_;i++){
        for(int j=0;j<n_;j++){
            if(abs(*(*(A_+i)+j))>abs(max) && i!=j){
                max=*(*(A_+i)+j);j_max=j;i_max=i;
            }


        }
    }
    //this line terminates the loop when the matrix is diagonalized:
    if(epsilon>abs(max)){/*cout<<"max value is : "<<max;*/break;}
    int k;int l;
    if(i_max>j_max){l=i_max;k=j_max;}
    else           {l=j_max;k=i_max;}


    // PART 2
    //introducing tau as it is written in a), to find the angle:
    double tau;
    tau = (A_[l][l]-A_[k][k])/(2*A_[k][l]);

    //using the quadratic formula t=-tau+-sqrt(1+tau^2) and choosing the one with least absolute value...
    //... meaning least rotation.
    double t1=-tau+sqrt(1+pow(tau,2));double t2=-tau-sqrt(1+pow(tau,2));
    t2 =abs(t1);t2=abs(t2);
    double t;
    if(t1>t2){ t=t1;}
    else{t=t2;}

    double c = 1./sqrt(1+pow(t,2));
    double s = t*c;
    //  PART 3 updating the matrix:
    // updating the matrix A to A' by updating values with matrix multiplication A'=StAS
    matr_rotation(A_,c,s,k,l,n_);
    amount_of_trans+=1;


      }
     cout<< "Jacobi_rotations ended, amount of transformations: "<<amount_of_trans<<endl;
}

void matr_rotation(double **A, double c,double s,int k, int l,int n){//implementing the rotation by given formulas
    double b_i_k[n];double b_i_l[n];double b_k_k;double b_l_l;


    for(int i=0;i<n;i++){
        b_i_k[i]=A[i][k]*c-A[i][l]*s;
        b_i_l[i]=A[i][l]*c+A[i][k]*s;

    }
    b_k_k=A[k][k]*c*c-2*A[k][l]*c*s+A[l][l]*s*s;
    b_l_l=A[l][l]*c*c+2*A[k][l]*c*s+A[k][k]*s*s;
    for(int i=0;i<n;i++){
        //setting column values:
        A[i][k]=b_i_k[i];
        A[i][l]=b_i_l[i];
        //setting row values:
        A[k][i]=b_i_k[i];
        A[l][i]=b_i_l[i];
    }
    A[k][k]=b_k_k;
    A[l][l]=b_l_l;
    A[k][l]=0;
    A[l][k]=0;
    //checking if b_l_k and b_k_l is really zero

}
void print_out_matrix(double **A,int n){
  //  cout<<"printing out the matrix: "<<endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}
void matrix_construct(double ro_min, double ro_max,int n,double **A){
    cout <<"initializing a construct of "<<n<<" point matrix A"<<endl;
    //defining step length:
    double h=(ro_max-ro_min)/n;
    double ro=ro_min;
//    double **A;
//    A=new double*[n-1];

    double d1=2./(h*h);

    double e1=-1/(h*h);

    for(int i=0;i<n;i++){
        A[i]=new double[n];
        ro=ro_min+h*i;
      //  cout<< ro<<endl;

        for(int j=0;j<n;j++){

            if(i==j){A[i][j]=d1+ro*ro;}
            else if(j==(i-1) || j==(i+1)){A[i][j]=e1;}
            else{
                A[i][j]=0;
            }

        }

    }
    cout<< "ending matrix construct"<<endl;

}
void find_min_eigenvalues(double **A,int n){
    double *a;
    a=new double[n];
    double min=100;
    for(int i=0;i<n;i++){
        a[i]=A[i][i];
        if(a[i]<min){min=a[i];}
//        cout<<"very excite : "<< a[i]<<endl;
    }
    sort(a,a+n);
    for(int i=0;i<3;i++){cout<< "eig value: "<< a[i]<<endl;}
//    cout<< "finish with the min value: "<< min<<endl;
}
