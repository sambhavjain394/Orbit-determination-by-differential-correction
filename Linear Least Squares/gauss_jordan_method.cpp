#include<bits/stdc++.h>
using namespace std;
# define ll long long

//Function for printing Matrix
void Print_Matrix(vector<vector<long double> > A)
{
    for(ll i=0;i<A.size();i++)
    {
        for(ll j=0;j<A[i].size();j++)
        {
            cout << A[i][j] << "  ";
        }
      printf("\n");
    }
    return ;
}

// Function for finding transpose of a matrix
vector<vector<long double> > transpose(vector<vector<long double> > A)
{
    vector<vector<long double> > trans(A[0].size(),vector<long double>(A.size()));

    for(ll i=0;i<A.size();i++)
        for(ll j=0;j<A[i].size();j++)
    {
        trans[j][i]=A[i][j];
    }
    return trans;
}

// Function for finding matrix multiplication
vector<vector<long double> > multiplication(vector<vector<long double> > A, vector<vector<long double> > B)
{
     if(A[0].size()!=B.size()){
        cout<<"Matrix dimensions are not compatible for multiplication (column1 != row2)\n";
       exit(0);
   }
    vector<vector<long double> > C(A.size(), vector<long double>(B[0].size()));
    for (ll i=0;i<A.size();i++)
        for (ll j=0; j< B[0].size();j++)
        { C[i][j]=0;
            for (ll k=0; k< B.size(); k++)
            {
                C[i][j]+=A[i][k]*B[k][j];
            }
        }
    return C;
}

//Function for creating Augmented Matrix
vector<vector<long double> > Augmented_Matrix(vector<vector<long double> > A)
{
    vector<vector<long double> > Augmented_Matrix(A.size(),vector<long double>(2*(A[0].size())));

    for(ll i=0;i<A.size();i++)
    {
         for(ll j=0;j<2*A[i].size();j++)
        {
            if (j == (i + A.size()))
                 Augmented_Matrix[i][j] = 1;
            else if (j >= A.size())
                Augmented_Matrix[i][j] = 0;
            else
                Augmented_Matrix[i][j] = A[i][j];
        }
    }

    long double temp = 0;
    for (ll i = Augmented_Matrix.size() - 1; i > 0; i--){
        // Swapping each and every element of the two rows
        if (Augmented_Matrix[i - 1][0] < Augmented_Matrix[i][0]){
         for (ll j = 0; j < Augmented_Matrix[0].size(); j++){
         temp = Augmented_Matrix[i][j];
         Augmented_Matrix[i][j] = Augmented_Matrix[i - 1][j];
         Augmented_Matrix[i - 1][j] = temp;
            }
        }
    }
    return Augmented_Matrix;
}


//Calculation of Inverse of Matrix
vector<vector<long double> > Inverse_Matrix(vector<vector<long double> > A)
{
    vector<vector<long double> > Inverse_Matrix(A.size(),vector<long double>((A.size())));

   if(A[0].size()!=A.size()){
       cout<<" Matrix Inverse is not possible (It is not Square Matrix) \n";
       exit(0);
   }

   vector<vector<long double> > A1;
   A1 = Augmented_Matrix(A);

  long double temp=0;
    for (ll i = 0; i < A1.size(); i++){
        for (ll j = 0; j < A1.size(); j++){
            if (j != i){
            temp = A1[j][i] / A1[i][i];

                for (ll k = 0; k < A1[0].size(); k++){
                        A1[j][k] -= A1[i][k] * temp;
                }
            }
        }
    }

    for (ll i = 0; i < A1.size(); i++){
        temp = A1[i][i];
        for (ll j = 0; j < A1[0].size(); j++){
            A1[i][j] = A1[i][j] / temp;
        }
    }

    for(ll i=0;i<A1.size();i++){
        for(ll j=0;j<A1.size();j++)
        {
            Inverse_Matrix[i][j]=A1[i][j+A1.size()];
            if (Inverse_Matrix[i][j]!=Inverse_Matrix[i][j]){
                cout<<" Matrix Inverse is not possible (Determinant is zero) \n";
                exit(0);
            }
        }
    }

   return Inverse_Matrix;
}


int main()
{
    pair<long double,long double> p[8];

    p[0]=make_pair(1,1);
    p[1]=make_pair(2,1);
    p[2]=make_pair(3,2);
    p[3]=make_pair(4,3);
    p[4]=make_pair(5,3);
    p[5]=make_pair(6,4);
    p[6]=make_pair(7,4);
    p[7]=make_pair(8,6);

    vector<vector<long double> > A;
    int n=8;
    for(ll i=0;i<n;i++)
    {
        vector<long double> temp;
        temp.push_back(1LL);
        temp.push_back(p[i].first);

        A.push_back(temp);
    }

    vector<vector<long double> > b;
    for(ll i=0; i<n; i++)
    {
        vector<long double> temp;
        temp.push_back(p[i].second);
        b.push_back(temp);
    }

    vector<vector<long double> > X;
    X=Inverse_Matrix(multiplication(transpose(A),A));
    Print_Matrix(X);

    vector<vector<long double> > Y;
    Y=multiplication(Inverse_Matrix(multiplication(transpose(A),A)),multiplication(transpose(A),b));
    Print_Matrix(Y);

    vector<vector<long double> > Z;

long double W1[] = {1, 4, 7};
long double W2[] = {2, 5, 8};
long double W3[] = {3, 6, 9};

ll n1 = 3;
    for(ll i=0; i<n1; i++){
        vector<long double> temp1;
        temp1.push_back(W1[i]);
        temp1.push_back(W2[i]);
        temp1.push_back(W3[i]);
        Z.push_back(temp1);
    }
    Print_Matrix(Z);
    Z=Inverse_Matrix(Z);
    Print_Matrix(Z);
}
