#include<bits/stdc++.h>
using namespace std;
# define ll long long
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
// Function for finding cofactors of matrix
vector<vector<long double> > cofactor(vector<vector<long double> >A,int p, int q)
{
    int i=0, j=0;
    int n;
    n=A.size();
    vector<vector<long double> > B(n-1,vector<long double>(n-1));
    for (int row =0; row<n; row++)
    {
        for (int col=0; col<n; col++)
        {
            if (row!= p && col != q)
            {
                B[i][j++] = A[row][col];
                if (j == n-1)
                {
                    j=0;
                    i++;
                }
            }
        }
    }
    return B;
}

//Function for finding determinant of matrix
long double determinant(vector<vector<long double> > A)
{
    long double D;
    D=0;
    if (A.size()==1)
        return A[0][0];
    int sign=1;
    // Finding the determinant of matrix using 1st row
    for (int i=0; i<A.size(); i++)
    {
        D+=sign*A[0][i]*determinant(cofactor(A,0,i));

        // Changing sign for each alternative element
        sign=-sign;
    }
    return D;
}

// Function for finding adjoint of matrix
vector<vector<long double> > adjoint(vector<vector<long double> > A)
{
    int N=A.size();
    vector<vector<long double> > B(N,vector<long double>(N));
    if (N==1)
    {
        B[0][0]=1;
        return B;
    }
    int sign=1;
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            sign=((i+j)%2==0)? 1 : -1;
            B[j][i]=sign*determinant(cofactor(A,i,j));
        }
    }
    return B;
}


// Function for finding matrix inverse
vector<vector<long double> > inverse(vector<vector<long double> > A)
{
    long double D;
    D=determinant(A);
    if (A.size()!=A[0].size())
    {
       cout<<"Matrix Inverse cannot be calculated-- Not a Square matrix";
       exit(0);
    }
    if (determinant(A)==0)
    {
       cout<<"Matrix Inverse cannot be calculated-- Singular Matrix";
       exit(0);
    }
    int N=A.size();
    vector<vector<long double> > adj(N,vector<long double>(N));
    adj=adjoint(A);

    vector<vector<long double> > C(N,vector<long double>(N));

    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
    {
     C[i][j]=adj[i][j]/D;
    }
    return C;
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
    for(ll i=0;i<n;i++){
        vector<long double> temp;
        temp.push_back(1LL);
        temp.push_back(p[i].first);
        A.push_back(temp);
    }

    vector<vector<long double> > b;
    for(ll i=0; i<n; i++){
        vector<long double> temp;
        temp.push_back(p[i].second);
        b.push_back(temp);
    }

    vector<vector<long double> > X;
    X=inverse(multiplication(transpose(A),A));
    for (int i=0; i<X.size();i++){
        for (int j=0; j<X[0].size();j++)
            cout<<X[i][j]<<" ";
        cout<< endl;
    }

    vector<vector<long double> > Y;
    Y=multiplication(inverse(multiplication(transpose(A),A)),multiplication(transpose(A),b));
    for (int i=0; i<Y.size();i++){
        for (int j=0; j<Y[0].size();j++)
            cout<<Y[i][j]<<" ";
        cout<< endl;
    }
}
