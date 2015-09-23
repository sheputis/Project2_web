#include <iostream>
#include <cmath>
//#include "armadillo"


using namespace std;
//using namespace arma;

//void Jacobi_rotation();
void max_non_diag(double **,int); //matrix double ** A, dimention of the matrix int n

int main()
{

    int n = 2;
    double **A;
    A=new double*[n];

    for(int i=0;i<n;i++){
        A[i]=new double[n];
        for(int j=0;j<n;j++){
            A[i][j]=n*i+j;
            //cout<<A[i][j];
        }
    }


    max_non_diag(A,n);


    /*
    for(int i=0;i<n;i++){

        for(int j=0;i<n;i++){
       //     A[i][j]=1;
            cout<<"heyVV";
        }
    }
*/


   // cout << "Hello World!"<<endl;
    return 0;
}

void max_non_diag(double ** A_, int n_){
    double epsilon = pow(10,-1); //epsilon is going to be our boundary for 0
    int i_max;int j_max;
    double max=1;
    while(max>epsilon){
    // PART 1 :
    // a for loop that finds max non diag value:
    double max=0;
    for (int i=0;i<n_;i++){

        for(int j=0;j<n_;j++){
            if(*(*(A_+i)+j)>max && i!=j){
                max=*(*(A_+i)+j);j_max=j;i_max=i;
            }
            cout<<*(*(A_+i)+j)<<endl;

        }
    }
    cout <<"maximum value :"<<max<<endl<<"i and j for max values:"<<i_max<<" , "<<j_max<<endl;
    int k;int l;
    if(i_max>j_max){l=i_max;k=j_max;}
    else           {l=j_max;k=i_max;}


    // PART 2
    //introducing tau as it is written in a), to find the angle:
    double tau;
    tau = ((*(*(A_+l)+l))-(*(*(A_+k)+k)))/(2*(*(*(A_+k)+l))); // a very confusing way of pointing to variables
    //using the quadratic formula t=-tau+-sqrt(1+tau^2) and choosing the one with least absolute value...
    //... meaning least rotation.
    double t1=-tau+sqrt(1+pow(tau,2));double t2=-tau-sqrt(1+pow(tau,2));
    t2 =abs(t1);t2=abs(t2);
    double t;
    if(t1>t2){ t=t1;}
    else{t=t2;}
    cout <<"tangens , smaller, bigger: "<<t<<" "<<t1<<t2<<endl;
    //cosine c and sine s is :
    double c = 1./sqrt(1+pow(t,2));
    double s = t*c;
    cout << "cosine and sine: "<<c<<" "<<s;


    //  PART 3 updating the matrix:
    // updating the matrix A to A' by updating values with matrix multiplication A'=StAS
    break;
    }
}

