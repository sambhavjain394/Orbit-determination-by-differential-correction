#include<bits/stdc++.h>
#include <cerrno>
#include <cstdlib>
using namespace std;
# define ll long long

// Defining constant pi
long double pi = 4*atan(1);

// Earth Parameters
long double R_Earth=6378.1363*1000;      // (m)
long double omega_1=7.292115*pow(10,-5); // (rad/sec)
long double ecc_Earth=0.081819221456;
long double Arcs=3600*180/pi;
long double MJD_J2000=51544.5;
long double dtr=pi/180;

// Station Parameters
struct station
{
    string name;
    long double alt;
    long double lat;
    long double longi;
};

// Function for saving station data
// OMITTED for PRIVACY ISSUES of ISRO

// Remainder evaluation on dividing two numbers
long double find_mod(long double A, long double B)
{
    long double find_mod;

    if ((A/B)>0)
        find_mod = A/B - floor(A/B);
    else if ((A/B)<0)
        find_mod = A/B - ceil(A/B);

    return find_mod;
}

////////////////////////////////////////////////////
//Section - Printing Matrix
////////////////////////////////////////////////////

// Printing the 2D vector
void print_matrix_2d_ld (vector<vector<long double> > Data)
{
    for (int i=0; i<Data.size();i++)
    {
        for (int j=0; j<Data[0].size();j++)
        {
          cout<<fixed<<setprecision(10)<<Data[i][j]<<" ";
        }
        cout<<endl;
    }
}

// Print 2D String vector Data
void print_matrix_2d_str (vector<vector<string> > Data)
{
    for (int i=0; i<Data.size();i++)
    {
        for (int j=0; j<Data[0].size();j++)
        {
          cout<<Data[i][j]<<" ";
        }
        cout<<endl;
    }
}

// Printing the 1D vector Data
void print_matrix_1d_ld (vector<long double> Data)
{
    for (int i=0; i<Data.size();i++)
    {
        cout<<Data[i]<<" "<<endl;
    }
}

//////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Rotation Matrix x, y and z
/////////////////////////////////////////////////////

// Function for defining Rotation matrix along x axis
vector<vector<long double>> R1(long double theta)
{
    vector<vector<long double>> Rot_Mat;
    vector<long double> temp;
    temp.push_back(1);
    temp.push_back(0);
    temp.push_back(0);
    Rot_Mat.push_back(temp);
    temp.clear();
    temp.push_back(0);
    temp.push_back(cos(theta));
    temp.push_back(sin(theta));
    Rot_Mat.push_back(temp);
    temp.clear();
    temp.push_back(0);
    temp.push_back(-sin(theta));
    temp.push_back(cos(theta));
    Rot_Mat.push_back(temp);

    return Rot_Mat;
}

// Function for defining Rotation matrix along y axis
vector<vector<long double>> R2(long double theta)
{
    vector<vector<long double>> Rot_Mat;
    vector<long double> temp;
    temp.push_back(cos(theta));
    temp.push_back(0);
    temp.push_back(-sin(theta));
    Rot_Mat.push_back(temp);
    temp.clear();
    temp.push_back(0);
    temp.push_back(1);
    temp.push_back(0);
    Rot_Mat.push_back(temp);
    temp.clear();
    temp.push_back(sin(theta));
    temp.push_back(0);
    temp.push_back(cos(theta));
    Rot_Mat.push_back(temp);

    return Rot_Mat;
}

// Function for defining Rotation matrix along z axis
vector<vector<long double>> R3(long double theta)
{
    vector<vector<long double>> Rot_Mat;
    vector<long double> temp;
    temp.push_back(cos(theta));
    temp.push_back(sin(theta));
    temp.push_back(0);
    Rot_Mat.push_back(temp);
    temp.clear();
    temp.push_back(-sin(theta));
    temp.push_back(cos(theta));
    temp.push_back(0);
    Rot_Mat.push_back(temp);
    temp.clear();
    temp.push_back(0);
    temp.push_back(0);
    temp.push_back(1);
    Rot_Mat.push_back(temp);

    return Rot_Mat;
}
////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Matrix Transformation
/////////////////////////////////////////////////////

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
vector<vector<long double>> Inverse_Matrix(vector<vector<long double>> A)
{
    long double temp;

    int order=A.size();

    vector<vector<long double>> aug_mat(order,vector<long double>(2*order));

    for (int i=0; i<order; i++)
        for (int j=0; j<2*order; j++)
    {
        if(j<order)
            aug_mat[i][j]=A[i][j];
        else
            if (j == (i+order) )
                aug_mat[i][j]=1;
        else
            aug_mat[i][j]=0;
    }

    for (int i=order -1; i>0; i--)
    {
     if (aug_mat[i-1][0] < aug_mat[i][0])
        for (int j=0; j< 2*order; j++)
     {   temp = aug_mat[i][j];
         aug_mat[i][j] = aug_mat[i-1][j];
         aug_mat[i-1][j] = temp;
     }
    }

    for (int i=0; i< order; i++)
    {
        for (int j=0; j< order; j++)
        {
            if (j!=i)
            {
                temp=aug_mat[j][i]/aug_mat[i][i];
                for (int k=0; k<2*order; k++)
                    aug_mat[j][k]-=aug_mat[i][k]*temp;
            }
        }
    }
    for (int i=0; i< order; i++)
    {
        temp=aug_mat[i][i];
        for (int j=0; j<2*order; j++)
        {
            aug_mat[i][j]=aug_mat[i][j]/temp;

        }
    }
    for (int i=0; i< order; i++)
        for (int j=order; j<2*order; j++)
            A[i][j-order]=aug_mat[i][j];

    return A;
}

// Calculate Addition of two vector
vector<vector<long double> > Add(vector<vector<long double> > Data1, vector<vector<long double> > Data2)
{
    vector<vector<long double> > Data(Data1.size(), vector<long double>(Data1[0].size()));
    if (Data1.size()==Data2.size() && Data1[0].size()==Data2[0].size())
    {
        for(int i=0; i<Data1.size(); i++)
        {
            for(int j=0; j<Data1[0].size(); j++)
            {
                Data[i][j] = Data1[i][j] + Data2[i][j];
            }
        }
    }
    else
    {
        cout<<"Vector addition is not possible";
        exit(0);
    }
    return Data;
}

///////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Reading Input File
/////////////////////////////////////////////////////

// Function to read text files
vector<vector<string> > read_text_file(string file_name)
{
  vector<vector<string> > Data;
  string line;

  ifstream myfile (file_name.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
            vector<string> temp_row;
            string temp_elem="";
            int i=0;
            while (i< line.size())
            {
                if (line[i]=='\t' || line[i]==' ')
                {
                    temp_row.push_back(temp_elem);
                    temp_elem="";
                    while (line[i]=='\t'  || line[i]==' ' && i<line.size())
                        i++;
                }

                temp_elem+=line[i];
                i++;
            }

            if (temp_elem.size()>0)
            {
                temp_row.push_back(temp_elem);
            }
            Data.push_back(temp_row);

        line.clear();
    }
  }
  else cout << "Unable to open file";
  myfile.close();
  return Data;
}

/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Finding Epoch Seconds
/////////////////////////////////////////////////////

// Function to find epoch seconds
vector<vector<long double> > get_epoch_seconds(vector<vector<long double> > Data_measured, vector<vector<long double> > Data1)
{
    //Data1-Reference Time for which nominal data is provided
    //Data2-Time of Data measured epoch

  //Storing values of time from Data Given
  vector<vector<long double> > Data2(Data_measured.size(),vector<long double>(7));
    for (int i=0; i<Data2.size();i++)
    {
        for (int j=0; j<Data2[0].size();j++)
        {
           Data2[i][j] = Data_measured[i][j];
        }
    }
  //Storing Epoch of measured data in seconds from reference time
  vector<vector<long double> > Time_seconds(Data2.size(),vector<long double>(1));

  long double millisecond1, millisecond2;
  tm time1, time2;
  //set values for time1...
  time1.tm_year = Data1[0][0] - 1900; //years counted staring from 1900
  time1.tm_mon = Data1[0][1] - 1;     //Jan == 0, Feb == 1 ... Dec == 11
  time1.tm_mday = Data1[0][2];        // 1 - 31
  time1.tm_hour = Data1[0][3];        // 0 - 23
  time1.tm_min = Data1[0][4];         // 0 - 59
  time1.tm_sec = Data1[0][5];         // 0 - 59
  millisecond1 = Data1[0][6];

  for(int i=0; i<Data2.size(); i++)
  {
    //set values for time2...
    time2.tm_year = Data2[i][0] - 1900; //years counted staring from 1900
    time2.tm_mon = Data2[i][1] - 1;     //Jan == 0, Feb == 1 ... Dec == 11
    time2.tm_mday = Data2[i][2];        // 1 - 31
    time2.tm_hour = Data2[i][3];        // 0 - 23
    time2.tm_min = Data2[i][4];         // 0 - 59
    time2.tm_sec = Data2[i][5];         // 0 - 59
    millisecond2 = Data2[i][6];

    //convert into seconds numbers
    time_t secs1=mktime(&time1);
    time_t secs2=mktime(&time2);
    time_t difference_in_seconds = secs2 - secs1;
    //Calculating seconds
    long double seconds;
    seconds = difference_in_seconds;

    if (millisecond2 < millisecond1)
    {
        millisecond2 += 1000;
        seconds -= 1;
    }
    //Calculating time including millisecond
    Time_seconds[i][0]=seconds + 0.001*(millisecond2-millisecond1);
  }
  return Time_seconds;
}

/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Lagrange interpolation of 7th order
/////////////////////////////////////////////////////

// Function for 7th order Lagrange interpolation
vector<vector<long double> > lagrange_interpolate(vector< vector<long double> > Data1 , vector< vector<long double> > Time_seconds)
{
    //Data1 is Nominal Data

    //Calculating Interpolated values at given time
    vector<vector<long double> > Result(Time_seconds.size(), vector<long double>(7));

    long double term1_1;
    long double term2_2;
    long double term3_3;
    long double term4_4;
    long double term5_5;
    long double term6_6;

    long double term1;
    long double term2;
    long double term3;
    long double term4;
    long double term5;
    long double term6;

    int i1;
    int i2;

  for(int k=0; k<Time_seconds.size();k++)
  {
    term1_1 = 0;
    term2_2 = 0;
    term3_3 = 0;
    term4_4 = 0;
    term5_5 = 0;
    term6_6 = 0;

    i1=0;
    i2=0;
    while(Data1[i1][2]<=Time_seconds[k][0])
    {
        i1=i1+1;
    }

    if(i1>4 && i1<(Data1.size()-3))
    {
        i2=i1-4;
    }
    else if(i1<4)
    {
        i2=i1;
    }
    else
    {
        i2=i1-7;
    }

    for(int i=i2;i<(i2+7);i++)
    {
        term1 = Data1[i][3];
        term2 = Data1[i][4];
        term3 = Data1[i][5];
        term4 = Data1[i][6];
        term5 = Data1[i][7];
        term6 = Data1[i][8];

        for(int j=i2;j<(i2+7);j++)
        {
            if (j!=i)
            {
                term1 = term1*(Time_seconds[k][0] - Data1[j][2])/(Data1[i][2] - Data1[j][2]);
                term2 = term2*(Time_seconds[k][0] - Data1[j][2])/(Data1[i][2] - Data1[j][2]);
                term3 = term3*(Time_seconds[k][0] - Data1[j][2])/(Data1[i][2] - Data1[j][2]);
                term4 = term4*(Time_seconds[k][0] - Data1[j][2])/(Data1[i][2] - Data1[j][2]);
                term5 = term5*(Time_seconds[k][0] - Data1[j][2])/(Data1[i][2] - Data1[j][2]);
                term6 = term6*(Time_seconds[k][0] - Data1[j][2])/(Data1[i][2] - Data1[j][2]);
            }
        }
        // Add current term to result
        term1_1 += term1;
        term2_2 += term2;
        term3_3 += term3;
        term4_4 += term4;
        term5_5 += term5;
        term6_6 += term6;
    }
    Result[k][0] = Time_seconds[k][0];
    Result[k][1] = term1_1;
    Result[k][2] = term2_2;
    Result[k][3] = term3_3;
    Result[k][4] = term4_4;
    Result[k][5] = term5_5;
    Result[k][6] = term6_6;
  }
  return Result;
}
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Evaluating UT1, Julian Date
/////////////////////////////////////////////////////

// Calculating MJD
long double JD_to_MJD(long double JD)
{
    long double MJD;
    MJD = JD - 2400000.5;
    return MJD;
}

long double MJD_to_JD(long double MJD)
{
    long double JD;
    JD = MJD + 2400000.5;
    return JD;
}

// Finding Mean Obliquity
long double MeanObliquity(long double Mjd_TT)
{
	long double temp;

	const long double T = (Mjd_TT - MJD_J2000)/36525.0;

	temp = 23.439291 - 0.0130042*T - (1.64E-07)*T*T + (5.04E-07)*T*T*T; //[deg]

	temp = temp*dtr; //[rad]

    return (temp);
}

// Calculating Julian Date
vector<vector<long double> > Get_JD(vector<vector<long double> > Measured_Data)
{
    long double year[Measured_Data.size()], mon[Measured_Data.size()], date[Measured_Data.size()], hour[Measured_Data.size()], minute[Measured_Data.size()], sec[Measured_Data.size()], ms[Measured_Data.size()];

    for(int i=0; i<Measured_Data.size(); i++)
    {
        year[i] = Measured_Data[i][0];
        mon[i] = Measured_Data[i][1];
        date[i] = Measured_Data[i][2];
        hour[i] = Measured_Data[i][3];
        minute[i] = Measured_Data[i][4];
        sec[i] = Measured_Data[i][5];
        ms[i] = Measured_Data[i][6];
    }

    //Julian Date Calculation
    vector<vector<long double> > JD(Measured_Data.size(), vector<long double>(1));
    for (int i=0; i<Measured_Data.size(); i++)
    {
        JD[i][0] = 367*year[i] - int(7*(year[i] + int((mon[i]+9)/12))/4) + int(275*mon[i]/9) + date[i] + 1721013.5 + ((((sec[i]+0.001*ms[i])/60+minute[i])/60+ hour[i])/24);
        //JD[i][1] = hour[i];
        //JD[i][2] = minute[i];
        //JD[i][3] = sec[i] + 0.001*ms[i];
    }
    return JD;
}

// Calculating Julian Date
vector<vector<long double> > Get_JD_UT1 ( vector<vector<long double> > JD_UTC , vector<vector<long double> > UT1_UTC )
{
    //Julian Date Calculation
    long double Mjd_UTC, Mjd_UT1;
    vector<vector<long double> > JD_UT1(JD_UTC.size(), vector<long double>(1));
    for (int i=0; i<JD_UTC.size(); i++)
    {
        Mjd_UTC = JD_to_MJD(JD_UTC[i][0]);
        Mjd_UT1 = Mjd_UTC + UT1_UTC[i][0]/86400.00;
        JD_UT1[i][0] = MJD_to_JD(Mjd_UT1);

    }
    return JD_UT1;
}

// Calculating Julian Date
vector<vector<long double> >  Get_JD_TT   ( vector<vector<long double> > JD_UTC , long double delta_AT )
{

    //Julian Date Calculation
    long double Mjd_UTC, Mjd_TT;
    vector<vector<long double> > JD_TT(JD_UTC.size(), vector<long double>(1));
    for (int i=0; i<JD_UTC.size(); i++)
    {
        Mjd_UTC = JD_to_MJD(JD_UTC[i][0]);
        Mjd_TT = Mjd_UTC + (delta_AT+32.184)/86400.00;
        JD_TT[i][0] = MJD_to_JD(Mjd_TT);

    }
    return JD_TT;
}
///////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Evaluating Corrections for Transformation in Coordinate System
/////////////////////////////////////////////////////


vector<long double> NutAngles_IAU1980(long double Mjd_TT)
{
	int i, N_coeff = 106;

	// Variables
	long double l, lp, F, D, Om;  //Mean arguments of luni-solar motion
	long double arg;

	// Constants
	const long double T  = (Mjd_TT - MJD_J2000)/36525.0;
	const long double T2 = T*T;
	const long double T3 = T2*T;

	const long double rev = 360.0*3600.0;  // arcsec/revolution [1degree = 3600 arcsec]

	const long double C[106][9] =
	{
		//
		// l  l' F  D Om    dpsi    *T     deps     *T       #
		//
		{  0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 },   //   1
		{  0, 0, 0, 0, 2,   20620,    2,   -8950,    5 },   //   2
		{ -2, 0, 2, 0, 1,     460,    0,    -240,    0 },   //   3
		{  0,-2, 2,-2, 1,     -20,    0,      10,    0 },   //   7
		{  2, 0,-2, 0, 1,      10,    0,       0,    0 },   //   8
		{  0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 },   //   9
		{  0, 1, 0, 0, 0,   14260,  -34,     540,   -1 },   //  10
		{  0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 },   //  11
		{  2, 0,-2, 0, 0,     110,    0,       0,    0 },   //   4
		{ -2, 0, 2, 0, 2,     -30,    0,      10,    0 },   //   5
		{  1,-1, 0,-1, 0,     -30,    0,       0,    0 },   //   6
		{  0,-1, 2,-2, 2,    2170,   -5,    -950,    3 },   //  12
		{  0, 0, 2,-2, 1,    1290,    1,    -700,    0 },   //  13
		{  2, 0, 0,-2, 0,     480,    0,      10,    0 },   //  14
		{  0, 0, 2,-2, 0,    -220,    0,       0,    0 },   //  15
		{  0, 2, 0, 0, 0,     170,   -1,       0,    0 },   //  16
		{  0, 1, 0, 0, 1,    -150,    0,      90,    0 },   //  17
		{  0, 2, 2,-2, 2,    -160,    1,      70,    0 },   //  18
		{  0,-1, 0, 0, 1,    -120,    0,      60,    0 },   //  19
		{ -2, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  20
		{  0,-1, 2,-2, 1,     -50,    0,      30,    0 },   //  21
		{  2, 0, 0,-2, 1,      40,    0,     -20,    0 },   //  22
		{  0, 1, 2,-2, 1,      40,    0,     -20,    0 },   //  23
		{  1, 0, 0,-1, 0,     -40,    0,       0,    0 },   //  24
		{  2, 1, 0,-2, 0,      10,    0,       0,    0 },   //  25
		{  0, 0,-2, 2, 1,      10,    0,       0,    0 },   //  26
		{  0, 1,-2, 2, 0,     -10,    0,       0,    0 },   //  27
		{  0, 1, 0, 0, 2,      10,    0,       0,    0 },   //  28
		{ -1, 0, 0, 1, 1,      10,    0,       0,    0 },   //  29
		{  0, 1, 2,-2, 0,     -10,    0,       0,    0 },   //  30
		{  0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 },   //  31
		{  1, 0, 0, 0, 0,    7120,    1,     -70,    0 },   //  32
		{  0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 },   //  33
		{  1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 },   //  34
		{  1, 0, 0,-2, 0,   -1580,    0,     -10,    0 },   //  35
		{ -1, 0, 2, 0, 2,    1230,    0,    -530,    0 },   //  36
		{  0, 0, 0, 2, 0,     630,    0,     -20,    0 },   //  37
		{  1, 0, 0, 0, 1,     630,    1,    -330,    0 },   //  38
		{ -1, 0, 0, 0, 1,    -580,   -1,     320,    0 },   //  39
		{ -1, 0, 2, 2, 2,    -590,    0,     260,    0 },   //  40
		{  1, 0, 2, 0, 1,    -510,    0,     270,    0 },   //  41
		{  0, 0, 2, 2, 2,    -380,    0,     160,    0 },   //  42
		{  2, 0, 0, 0, 0,     290,    0,     -10,    0 },   //  43
		{  1, 0, 2,-2, 2,     290,    0,    -120,    0 },   //  44
		{  2, 0, 2, 0, 2,    -310,    0,     130,    0 },   //  45
		{  0, 0, 2, 0, 0,     260,    0,     -10,    0 },   //  46
		{ -1, 0, 2, 0, 1,     210,    0,    -100,    0 },   //  47
		{ -1, 0, 0, 2, 1,     160,    0,     -80,    0 },   //  48
		{  1, 0, 0,-2, 1,    -130,    0,      70,    0 },   //  49
		{ -1, 0, 2, 2, 1,    -100,    0,      50,    0 },   //  50
		{  1, 1, 0,-2, 0,     -70,    0,       0,    0 },   //  51
		{  0, 1, 2, 0, 2,      70,    0,     -30,    0 },   //  52
		{  0,-1, 2, 0, 2,     -70,    0,      30,    0 },   //  53
		{  1, 0, 2, 2, 2,     -80,    0,      30,    0 },   //  54
		{  1, 0, 0, 2, 0,      60,    0,       0,    0 },   //  55
		{  2, 0, 2,-2, 2,      60,    0,     -30,    0 },   //  56
		{  0, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  57
		{  0, 0, 2, 2, 1,     -70,    0,      30,    0 },   //  58
		{  1, 0, 2,-2, 1,      60,    0,     -30,    0 },   //  59
		{  0, 0, 0,-2, 1,     -50,    0,      30,    0 },   //  60
		{  1,-1, 0, 0, 0,      50,    0,       0,    0 },   //  61
		{  2, 0, 2, 0, 1,     -50,    0,      30,    0 },   //  62
		{  0, 1, 0,-2, 0,     -40,    0,       0,    0 },   //  63
		{  1, 0,-2, 0, 0,      40,    0,       0,    0 },   //  64
		{  0, 0, 0, 1, 0,     -40,    0,       0,    0 },   //  65
		{  1, 1, 0, 0, 0,     -30,    0,       0,    0 },   //  66
		{  1, 0, 2, 0, 0,      30,    0,       0,    0 },   //  67
		{  1,-1, 2, 0, 2,     -30,    0,      10,    0 },   //  68
		{ -1,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  69
		{ -2, 0, 0, 0, 1,     -20,    0,      10,    0 },   //  70
		{  3, 0, 2, 0, 2,     -30,    0,      10,    0 },   //  71
		{  0,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  72
		{  1, 1, 2, 0, 2,      20,    0,     -10,    0 },   //  73
		{ -1, 0, 2,-2, 1,     -20,    0,      10,    0 },   //  74
		{  2, 0, 0, 0, 1,      20,    0,     -10,    0 },   //  75
		{  1, 0, 0, 0, 2,     -20,    0,      10,    0 },   //  76
		{  3, 0, 0, 0, 0,      20,    0,       0,    0 },   //  77
		{  0, 0, 2, 1, 2,      20,    0,     -10,    0 },   //  78
		{ -1, 0, 0, 0, 2,      10,    0,     -10,    0 },   //  79
		{  1, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  80
		{ -2, 0, 2, 2, 2,      10,    0,     -10,    0 },   //  81
		{ -1, 0, 2, 4, 2,     -20,    0,      10,    0 },   //  82
		{  2, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  83
		{  1, 1, 2,-2, 2,      10,    0,     -10,    0 },   //  84
		{  1, 0, 2, 2, 1,     -10,    0,      10,    0 },   //  85
		{ -2, 0, 2, 4, 2,     -10,    0,      10,    0 },   //  86
		{ -1, 0, 4, 0, 2,      10,    0,       0,    0 },   //  87
		{  1,-1, 0,-2, 0,      10,    0,       0,    0 },   //  88
		{  2, 0, 2,-2, 1,      10,    0,     -10,    0 },   //  89
		{  2, 0, 2, 2, 2,     -10,    0,       0,    0 },   //  90
		{  1, 0, 0, 2, 1,     -10,    0,       0,    0 },   //  91
		{  0, 0, 4,-2, 2,      10,    0,       0,    0 },   //  92
		{  3, 0, 2,-2, 2,      10,    0,       0,    0 },   //  93
		{  1, 0, 2,-2, 0,     -10,    0,       0,    0 },   //  94
		{  0, 1, 2, 0, 1,      10,    0,       0,    0 },   //  95
		{ -1,-1, 0, 2, 1,      10,    0,       0,    0 },   //  96
		{  0, 0,-2, 0, 1,     -10,    0,       0,    0 },   //  97
		{  0, 0, 2,-1, 2,     -10,    0,       0,    0 },   //  98
		{  0, 1, 0, 2, 0,     -10,    0,       0,    0 },   //  99
		{  1, 0,-2,-2, 0,     -10,    0,       0,    0 },   // 100
		{  0,-1, 2, 0, 1,     -10,    0,       0,    0 },   // 101
		{  1, 1, 0,-2, 1,     -10,    0,       0,    0 },   // 102
		{  1, 0,-2, 2, 0,     -10,    0,       0,    0 },   // 103
		{  2, 0, 0, 2, 0,      10,    0,       0,    0 },   // 104
		{  0, 0, 2, 4, 2,     -10,    0,       0,    0 },   // 105
		{  0, 1, 0, 1, 0,      10,    0,       0,    0 }    // 106
	};

	// Mean arguments of luni-solar motion

	l  = find_mod( 485866.733 + (1325.0*rev +  715922.633)*T  + 31.310*T2 + 0.064*T3, rev );  //mean anomaly of the Moon

	lp = find_mod( 1287099.804 + (99.0*rev + 1292581.224)*T   -  0.577*T2 - 0.012*T3, rev );  //mean anomaly of the Sun

	F  = find_mod( 335778.877 + (1342.0*rev +  295263.137)*T  - 13.257*T2 + 0.011*T3, rev );  //mean argument of latitude

	D  = find_mod( 1072261.307 + (1236.0*rev + 1105601.328)*T -  6.891*T2 + 0.019*T3, rev );  //mean longitude elongation of the Moon from the Sun

	Om = find_mod(  450160.280 - ( 5.0*rev +  482890.539)*T   +  7.455*T2 + 0.008*T3, rev );  //mean longitude of the ascending node


	// Nutation in longitude and obliquity [rad]
	long double deps = 0.0;	//initializing to zero in order to avoid storing garbage value
	long double dpsi = 0.0;	//initializing to zero in order to avoid storing garbage value

	for (i=0; i<N_coeff; i++)
	{
		arg   =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*Om ) / Arcs;

		dpsi +=  ( C[i][5]+C[i][6]*T ) * sin(arg);//Nuation in longitude
		deps +=  ( C[i][7]+C[i][8]*T ) * cos(arg);//Nuation in obliquity
	}

	dpsi = (1.0E-5) * (dpsi)/Arcs; //[rad]
	deps = (1.0E-5) * (deps)/Arcs; //[rad]

	//printf("\n");
	//printf("Nutation in Longitude (dpsi) and Obliquity (deps) [rad]: %lf %lf \n", *dpsi, *deps);
	//printf("Nutation in Longitude (dpsi) and Obliquity (deps) [deg]: %lf %lf \n", (*dpsi)*rtd, (*deps)*rtd);
	//printf("\n");
	vector<long double> data;
	data.push_back(dpsi);
	data.push_back(deps);

	return data;
}


// Finding Equation of Equinox
long double EqnEquinox_IAU1980(long double Mjd_TT,long double dPsi, long double dEps)
{
	long double dpsi, deps;  //Nutation angles
	long double Cor_dpsi, Cor_deps;

	long double longAscNodeLunOrb;

	long double term1, term2;

	long double temp;

	const long double r  = 360.0; //[deg] (1 revolution)

	const long double T  = (Mjd_TT - MJD_J2000)/36525.0;

	vector<long double> angle_op;
	// Nutation in longitude and obliquity
	angle_op=NutAngles_IAU1980(Mjd_TT); //calling NutAngles function
	dpsi=angle_op[0];
	deps=angle_op[1];

	Cor_dpsi = dpsi + dPsi;	//Corrected Nutation angle in Longitude [rad]
	Cor_deps = deps + dEps;	//Corrected Nutation angle in ecliptic  [rad]

	longAscNodeLunOrb = 125.04455501 - (   5*r + 134.1361851)*T + 0.0020756*T*T + (2.139E-06)*T*T*T - (1.650E-08)*T*T*T*T;	//[deg]

	while(longAscNodeLunOrb < -360.0)  longAscNodeLunOrb += 360.0; //[deg]
	while(longAscNodeLunOrb > +360.0)  longAscNodeLunOrb -= 360.0; //[deg]

	longAscNodeLunOrb = longAscNodeLunOrb*dtr; //[rad]

	term1 = (0.002649*(1.0/3600.0)*dtr)*sin(longAscNodeLunOrb); //               [rad]
	//term2 = (-0.000013*(1.0/3600.0)*dtr)*cos(longAscNodeLunOrb);// Montenbruck [rad]
	term2 = (0.000063*(1.0/3600.0)*dtr)*sin(2.0*longAscNodeLunOrb);// Vallado    [rad]

	// Equation of the equinoxes
	//temp = dpsi * cos( MeanObliquity(Mjd_TT) ); //[rad]
	temp = Cor_dpsi * cos( MeanObliquity(Mjd_TT) ) + term1 + term2; //[rad]

	/*printf("\n");
	printf("Eqn. of Equinox [rad]: %lf\n", temp);
	printf("\n");*/

	return  (temp);
}


long double GMST(long double Mjd_UT1)
{
	// Variables
	double Mjd_0, UT1, T_0, T, gmst, temp;

	// Constants
	const double Secs = 86400.0;        // seconds per day

	// Mean Sidereal Time
	Mjd_0 = floor(Mjd_UT1);

	UT1   = (Mjd_UT1 - Mjd_0)*Secs;          // [s]

	T_0   = (Mjd_0  - MJD_J2000)/36525.0;

	T     = (Mjd_UT1 - MJD_J2000)/36525.0;

	gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1 + (0.093104 - 6.2E-06*T)*T*T; // [s]

	temp = find_mod(gmst,Secs)*2*pi;	//[rad]

	return  temp;
}


long double GAST_IAU1980(long double JD_UT1, long double JD_TT, long double dPsi, long double dEps)
{
    long double theta_AST, temp, Mjd_UT1, Mjd_TT;
    Mjd_UT1=JD_to_MJD(JD_UT1);
    Mjd_TT=JD_to_MJD(JD_TT);

	theta_AST = GMST(Mjd_UT1) + EqnEquinox_IAU1980(Mjd_TT, dPsi, dEps);	//calling functions GMST & EqnEquinox
    temp = find_mod(theta_AST,2*pi)*2*pi; //calling modulo function

	return temp;
}

vector<vector<long double> > GHAMatrix_IAU1980(long double JD_UT1, long double JD_TT, long double dPsi, long double dEps)
{
	double temp;

	temp = GAST_IAU1980(JD_UT1, JD_TT, dPsi, dEps);

    vector<vector<long double> > Siderial_Time_Matrix;
    Siderial_Time_Matrix=R3(temp);//calling function R_z about z-axis (rotation through GAST) and storing outut in matrix Theta

	return Siderial_Time_Matrix;
}


//Calculate Polar Motion Matrix
vector<vector<long double>> Polar_Motion_Matrix(long double x_p, long double y_p)
{

    vector<vector<long double>> Polar_Motion_Matrix(3, vector<long double>(3,0));
    Polar_Motion_Matrix=multiplication( R2(-x_p),R1(-y_p));

    return Polar_Motion_Matrix;
}

/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
// Section - Frame conversion from TOD to ECEF
/////////////////////////////////////////////////////

// Function for Calculating conversion from TOD to ECEF
vector<vector<long double> > Nominal_TOD_to_ECEF(vector<vector<long double> > Measured_Data, vector<vector<long double> > Data_Nominal_TOD, vector<vector<long double> > Polar_angle, vector<vector<long double> > LOD ,vector<vector<long double> > delta_psi, vector<vector<long double> > delta_epsilon, vector<vector<long double> > JD_UT1, vector<vector<long double> > JD_TT)
{
    vector<vector<long double>> Data_Nominal_ECEf(Measured_Data.size(), vector<long double>(7));

    vector<vector<long double>> PM(3, vector<long double>(3,0));
    vector<vector<long double>> Theta(3, vector<long double>(3,0));
    vector<vector<long double>> dTheta(3, vector<long double>(3,0));

    //long double T_TDB;
    long double AST;
    long double omega_Earth;

    vector<vector<long double>> r_ECEF(1, vector<long double>(3));
    vector<vector<long double>> v_ECEF(1, vector<long double>(3));

    vector<vector<long double>> r_TOD(1, vector<long double>(3));
    vector<vector<long double>> v_TOD(1, vector<long double>(3));


    for(int i=0;i<Measured_Data.size();i++)
    {
        r_TOD[0][0] = Data_Nominal_TOD[i][1];
        r_TOD[0][1] = Data_Nominal_TOD[i][2];
        r_TOD[0][2] = Data_Nominal_TOD[i][3];

        v_TOD[0][0] = Data_Nominal_TOD[i][4];
        v_TOD[0][1] = Data_Nominal_TOD[i][5];
        v_TOD[0][2] = Data_Nominal_TOD[i][6];

        //T_TDB = get_T_TDB(Measured_Data[i]);

        PM = Polar_Motion_Matrix(Polar_angle[i][0]/Arcs, Polar_angle[i][1]/Arcs);
        Theta = GHAMatrix_IAU1980(JD_UT1[i][0], JD_TT[i][0], delta_psi[i][0]/Arcs, delta_epsilon[i][0]/Arcs);

        omega_Earth = (72921151.467064 - 0.843994809*1000*LOD[i][0])*pow(10,-12); // [rad/s]; Vallado Book

        AST = GAST_IAU1980(JD_UT1[i][0], JD_TT[i][0], delta_psi[i][0]/Arcs, delta_epsilon[i][0]/Arcs);

        dTheta[0][0] = -omega_Earth*sin(AST); dTheta[0][1] = omega_Earth*cos(AST);  dTheta[0][2] = 0.0;
        dTheta[1][0] = -omega_Earth*cos(AST); dTheta[1][1] = -omega_Earth*sin(AST); dTheta[1][2] = 0.0;
        dTheta[2][0] = 0.0;                   dTheta[2][1] = 0.0;                   dTheta[2][2] = 0.0;

        r_ECEF = transpose(multiplication(PM,multiplication(Theta,transpose(r_TOD))));
        v_ECEF = transpose(Add(multiplication(PM,multiplication(dTheta,transpose(r_TOD))),multiplication(PM,multiplication(Theta,transpose(v_TOD)))));

        Data_Nominal_ECEf[i][0] = Data_Nominal_TOD[i][0];
        Data_Nominal_ECEf[i][1] = r_ECEF[0][0];
        Data_Nominal_ECEf[i][2] = r_ECEF[0][1];
        Data_Nominal_ECEf[i][3] = r_ECEF[0][2];

        Data_Nominal_ECEf[i][4] = v_ECEF[0][0];
        Data_Nominal_ECEf[i][5] = v_ECEF[0][1];
        Data_Nominal_ECEf[i][6] = v_ECEF[0][2];
    }
    return Data_Nominal_ECEf;
}

// Site Coordinate from SEZ to ECEf frame
vector<vector<long double>> siteSVECEF(vector<long double> altgd, vector<long double> latgd, vector<long double> lon)
{
    vector<vector<long double>> svECEF;

    for (int j=0;j<altgd.size();j++)
    {
	int i;
	vector<long double> temp(6);
	double rdel, rk, cearth;


    //------  find rdel and rk components of site vector  ---------//
    cearth = R_Earth / sqrt(1.0 - (ecc_Earth*ecc_Earth*sin(latgd[j])*sin(latgd[j]) ) );
    rdel   = (cearth + altgd[j])*cos(latgd[j]);
    rk     = ((1.0 - ecc_Earth*ecc_Earth)*cearth + altgd[j])*sin(latgd[j]);

   	//--------------- Site position vector  -----------------//
   	temp[0] = rdel * cos(lon[j]);
    temp[1] = rdel * sin(lon[j]);
    temp[2] = rk;
	//--------------------------------------------------------//
   	//---------------  find site velocity vector  -----------------//
   	for(i=0; i<3; i++)
	 temp[i+3] = 0.0; // Velocity is zero bcz the coordinate system is fixed to the Earth w.r.t. ECEF frame
   	//-------------------------------------------------------------//
   	svECEF.push_back(temp);
    }
   	return svECEF;
}

////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Evaluating Measured Data (rho, rho_dot) from given Nominal Data
/////////////////////////////////////////////////////

// Function to calculate measurement variables Razel
long double RAZEL(long double measurement_type, vector<long double> sv_ECEF, vector<long double> sv_site_ECEF,long double site_latgd,long double site_lon)
{
	int i;
    long double Measure=0;

	vector<vector<long double>> Rho_SEZ, DRho_SEZ;
	long double rSiteECEF[3] = {0.0, 0.0, 0.0}, vSiteECEF[3] = {0.0, 0.0, 0.0};
	rSiteECEF[0]=sv_site_ECEF[0];
	rSiteECEF[1]=sv_site_ECEF[1];
	rSiteECEF[2]=sv_site_ECEF[2];
    vSiteECEF[0]=sv_site_ECEF[3];
    vSiteECEF[1]=sv_site_ECEF[4];
    vSiteECEF[2]=sv_site_ECEF[5];

	vector<vector<long double>> Rho_ECEF, DRho_ECEF;
	vector<long double> temp;
	long double r_ECEF[3] = {0.0, 0.0, 0.0}, v_ECEF[3] = {0.0, 0.0, 0.0};

	long double Rho,DRho;

	//-------------------------------//
	for(i=0; i<3; i++)
	{
		r_ECEF[i] = sv_ECEF[i];   // Satellite position vector in ECEF frame
		v_ECEF[i] = sv_ECEF[i+3]; // Satellite velocity vector in ECEF frame
	}
	//------------------------------//

	//----- ECEF range vector from site to satellite -----//
	temp.push_back(r_ECEF[0]-rSiteECEF[0]);
	Rho_ECEF.push_back(temp);
	temp.clear();
	temp.push_back(r_ECEF[1]-rSiteECEF[1]);
	Rho_ECEF.push_back(temp);
	temp.clear();
    temp.push_back(r_ECEF[2]-rSiteECEF[2]);
	Rho_ECEF.push_back(temp);
	temp.clear();

    temp.push_back(v_ECEF[0]-vSiteECEF[0]);
    DRho_ECEF.push_back(temp);
	temp.clear();
	temp.push_back(v_ECEF[1]-vSiteECEF[1]);
    DRho_ECEF.push_back(temp);
	temp.clear();
	temp.push_back(v_ECEF[2]-vSiteECEF[2]);
    DRho_ECEF.push_back(temp);
	temp.clear();

	//----------- Convert to SEZ for calculations -------------//
	Rho_SEZ=multiplication(R2(0.5*pi-site_latgd),multiplication(R3(site_lon),Rho_ECEF));

	DRho_SEZ=multiplication(R2(0.5*pi-site_latgd),multiplication(R3(site_lon),DRho_ECEF));

	//-------------  Range, Azimuth and Elevation -----------------------------//
	Rho = sqrt(Rho_ECEF[0][0]*Rho_ECEF[0][0]+Rho_ECEF[1][0]*Rho_ECEF[1][0] + Rho_ECEF[2][0]*Rho_ECEF[2][0]); // Range

    DRho = (Rho_SEZ[0][0]*DRho_SEZ[0][0]+Rho_SEZ[1][0]*DRho_SEZ[1][0]+Rho_SEZ[2][0]*DRho_SEZ[2][0])/Rho;

    if (measurement_type==3)
    {
        Measure = Rho;
    }
    if (measurement_type==6)
    {
        Measure = DRho;
    }

    return Measure;
}

//////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//Section - Evaluating Final Matrix that include various outputs
/////////////////////////////////////////////////////

// Function to calculate Nominal values for Measurement Data
vector<vector<long double> > Final_Data_Matrix(vector<vector<long double> > Measured_Data, vector<vector<long double> > Data_nominal_ECEF ,vector<vector<long double> > svsite_ECEF, vector<long double> site_altgd, vector<long double> site_latgd, vector<long double> site_lon, vector<vector<long double> > JD_UT1)
{
    vector<vector<long double> > Final_Data;

    for (int i=0; i<Measured_Data.size();i++)
    {
        long double measure_type;
        measure_type=Measured_Data[i][7]; // From the file, 7th column consists of Measurement type

        vector<long double> sv_ECEF(6);
        for(int j=0;j<6;j++)
        {
            sv_ECEF[j] = Data_nominal_ECEF[i][j+1];
        }

        long double Nominal_measure;
        Nominal_measure=RAZEL(measure_type,sv_ECEF,svsite_ECEF[i],site_latgd[i],site_lon[i]);

        vector<long double> temp;
        temp.push_back(Data_nominal_ECEF[i][0]); // Epoch seconds
        temp.push_back(Data_nominal_ECEF[i][1]); // Nominal r_ijk in ECEF
        temp.push_back(Data_nominal_ECEF[i][2]); // Nominal r_ijk in ECEF
        temp.push_back(Data_nominal_ECEF[i][3]); // Nominal r_ijk in ECEF
        temp.push_back(Data_nominal_ECEF[i][4]); // Nominal v_ijk in ECEF
        temp.push_back(Data_nominal_ECEF[i][5]); // Nominal v_ijk in ECEF
        temp.push_back(Data_nominal_ECEF[i][6]); // Nominal v_ijk in ECEF
        temp.push_back(Nominal_measure);         // Nominal Measurement
        temp.push_back(Measured_Data[i][8]);      // Actual Measurement
        temp.push_back(svsite_ECEF[i][0]);   // Station r_ijk in ECEF
        temp.push_back(svsite_ECEF[i][1]);   // Station r_ijk in ECEF
        temp.push_back(svsite_ECEF[i][2]);   // Station r_ijk in ECEF
        temp.push_back(svsite_ECEF[i][3]);   // Station v_ijk in ECEF
        temp.push_back(svsite_ECEF[i][4]);   // Station v_ijk in ECEF
        temp.push_back(svsite_ECEF[i][5]);   // Station v_ijk in ECEF
        temp.push_back(measure_type);
        temp.push_back(JD_UT1[i][0]);
        temp.push_back(site_latgd[i]);       // Station latitude
        temp.push_back(site_lon[i]);         // Station longitude
        Final_Data.push_back(temp);
    }

    return Final_Data;
}

/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
// Section - Functions for Evaluating Delta_x
/////////////////////////////////////////////////////

// Evaluating Partial Derivative Matrix
vector<vector<long double> > PD_Matrix(vector<vector<long double> > Final_Data, long double Delta_i, long double no_of_obs)
{
    vector<vector<long double> > A(no_of_obs, vector<long double>(6));
    vector<long double> sv_ECEF(6);
    vector<long double> sv_site_ECEF(6);
    long double measure_type, site_latgd, site_lon;

    vector<long double> sv_ECEF1(6);

    for(int i=0; i<6; i++)
    {
        sv_ECEF[i] = Final_Data[0][i+1];
        sv_site_ECEF[i] = Final_Data[0][i+9];
    }

    site_latgd=Final_Data[0][17];
    site_lon=Final_Data[0][18];
    measure_type = Final_Data[0][15];

    for(int i=0; i<no_of_obs; i++)
    {
        for(int j=0; j<6; j++)
        {
            sv_ECEF1 = sv_ECEF;
            sv_ECEF1[j] =  Delta_i*sv_ECEF1[j];
            A[i][j] = (RAZEL(measure_type, sv_ECEF1, sv_site_ECEF, site_latgd, site_lon)-RAZEL(measure_type, sv_ECEF, sv_site_ECEF, site_latgd, site_lon))/(sv_ECEF1[j]-sv_ECEF[j]);
        }
    }

    return A;
}

// Evaluating At_W_A in the Algorithm
vector<vector<long double> > At_W_A(vector<vector<long double> > Final_Data, long double Delta_i, long double no_of_obs, vector<vector<long double> > W)
{
    vector<vector<long double> > A;
    A = PD_Matrix(Final_Data, Delta_i, no_of_obs);

    vector<vector<long double> > At_W_A;
    At_W_A = multiplication(transpose(A),multiplication(W,A));

    return At_W_A;
}

// Evaluating At_W_b in the Algorithm
vector<vector<long double> > At_W_b(vector<vector<long double> > Final_Data, long double Delta_i, long double no_of_obs, vector<vector<long double> > W)
{
    vector<vector<long double> > b(no_of_obs,vector<long double> (1));

    for(int i=0; i<no_of_obs; i++)
        b[i][0]=Final_Data[0][8]-Final_Data[0][7];

    vector<vector<long double> > A;
    A = PD_Matrix(Final_Data, Delta_i, no_of_obs);

    vector<vector<long double> > At_W_b(6,vector<long double> (1));
    At_W_b = multiplication(transpose(A),multiplication(W,b));

    return At_W_b;
}

// Evaluation of RMS value
long double RMS_eval(vector<vector<long double> > Final_Data, long double no_of_obs, vector<vector<long double> > W)
{
    long double N = Final_Data.size();
    long double RMS=0;

    vector<vector<long double> > b(no_of_obs,vector<long double> (N));

    for(int i=0; i<no_of_obs; i++)
    {
        for(int j=0; j<N; j++)
        {
            b[i][j]=Final_Data[j][8]-Final_Data[j][7];
            RMS += (b[i][j]*W[i][i]*b[i][j])/((N-1)*no_of_obs);
        }
    }
    RMS = sqrt(RMS);
    return RMS;
}

/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
// Section - int main
/////////////////////////////////////////////////////

int main ()
{
// The code takes all the data in SI units

// If we want to add any ground station in the data, then we have to add in initialize_stations() function

// Code assumes consistency in format of files for measurement and nominal vector generation file

// Reading Data from text file

// Initializing Station Data
vector<station> Data_station;
Data_station=initialize_stations();

vector<vector<string> > Data1;
//vector<vector<string> > Data1_temp;
vector<vector<string> > Data2;
vector<vector<string> > Data3;

Data1=read_text_file("Data_Measured_Table.txt");
Data2=read_text_file("Nominal_Data_Table.txt");
Data3=read_text_file("EOP_data.txt");
long double no_of_obs = 1;  // rho, beta, theta, etc Here it is only rho/rho_dot

// Taking only range values
//for (int i=0; i<Data1_temp.size();i++)
//{
//  if (Data1_temp[i][7]=="6")
//       Data1.push_back(Data1_temp[i]);
//}

vector<vector<long double> > Data_Measured;
    vector<long double> site_alti,site_lat,site_lon;
    long double temp1;
    // Checking the station coordinates first and storing them
    for (int i=0; i<Data1.size();i++)
        for (int j=0; j<Data_station.size();j++)
        if (Data1[i][10]==Data_station[j].name) // 11th column stores the station names in Measurement Data
        {
            vector<long double> temp;
            for (int j=0; j<Data1[0].size()-2;j++) // Removing last two columns from file as one is orbit no. and other is station name
            {
                char * e;
                temp1=strtold(Data1[i][j].c_str(),&e);
                temp.push_back(temp1);
            }
            site_alti.push_back(Data_station[j].alt);
            site_lat.push_back(Data_station[j].lat);
            site_lon.push_back(Data_station[j].longi);
            Data_Measured.push_back(temp);
            break;
        }

// Converting Vector String Data to Vector Long Double : Nominal Data Temp
vector<vector<long double> > Nominal_Data_temp;
long double temp_ld1;
    for (int i=0; i<Data2.size();i++)
    {
        vector<long double> temp1;
        for (int j=0; j<Data2[0].size();j++)
        {
            char * e;
            temp_ld1=strtold(Data2[i][j].c_str(),&e);
            temp1.push_back(temp_ld1);
        }
    Nominal_Data_temp.push_back(temp1);
    }

// Converting Vector String Data to Vector Long Double : Data_UT1_X_Y
vector<vector<long double> > Data_EOP;
long double temp_ld2;
    for (int i=0; i<Data3.size();i++)
    {
        vector<long double> temp2;
        for (int j=0; j<Data3[0].size();j++)
        {
            char * e;
            temp_ld2=strtold(Data3[i][j].c_str(),&e);
            temp2.push_back(temp_ld2);
        }
    Data_EOP.push_back(temp2);
    }
//print_matrix_2d_ld(Data_EOP);

// Storing Delta_UT1   unit sec
//         Polar_angle unit arcsec
//         LOD         unit sec
//         d_psi       unit arcsec
//         d_epsilon   unit arcsec
vector<vector<long double> > Polar_angle(Data_Measured.size(),vector<long double> (2));
vector<vector<long double> > Delta_UT1(Data_Measured.size(),vector<long double> (1));
vector<vector<long double> > LOD(Data_Measured.size(),vector<long double> (1));
vector<vector<long double> > d_psi(Data_Measured.size(),vector<long double> (1));
vector<vector<long double> > d_epsilon(Data_Measured.size(),vector<long double> (1));

long double k1;
long double k2;

for (int i=0; i<Data_Measured.size();i++)
{
    for(int j=0; j<Data_EOP.size();j++)
    {
        if (Data_Measured[i][0]==Data_EOP[j][0] && Data_Measured[i][1]==Data_EOP[j][1] && Data_Measured[i][2]==Data_EOP[j][2])
        {
            Polar_angle[i][0] = Data_EOP[j][3];
            Polar_angle[i][1] = Data_EOP[j][4];
            Delta_UT1[i][0] = Data_EOP[j][5];
            LOD[i][0] = Data_EOP[j][6];
            d_psi[i][0] = Data_EOP[j][7];
            d_epsilon[i][0] = Data_EOP[j][8];
            break;
        }
    }
}

// Defining reference date
vector<vector<long double> > Reference_date(1,vector<long double>(7));
long double p[7]={2020,6,8,0,0,0,0};
    for(int i=0; i<7; i++)
    {
        Reference_date[0][i] = p[i];
    }

// Finding Epoch seconds for Measurement Data
vector<vector<long double> > epoch_seconds(Data_Measured.size(),vector<long double>(1));
epoch_seconds = get_epoch_seconds(Data_Measured,Reference_date);

//Finding Nominal Data in TOD by interpolating at Measurement epochs
vector<vector<long double> > Data_nominal_TOD(Data_Measured.size(),vector<long double>(7));
Data_nominal_TOD = lagrange_interpolate(Nominal_Data_temp, epoch_seconds);
//print_matrix_2d_ld(Data_nominal_TOD);

// Converting the nominal data from TOD frame into SI units: from second column to last column, converting km into m
for (int i=0; i<Data_nominal_TOD.size();i++)
    for (int j=1; j<Data_nominal_TOD[0].size();j++)
        Data_nominal_TOD[i][j]*=1000;

//Finding Julian Date
vector<vector<long double> > JD_UTC;
JD_UTC = Get_JD(Data_Measured);
//print_matrix_2d_ld(JD_UTC);

//Finding Julian Date
vector<vector<long double> > JD_UT1;
JD_UT1 = Get_JD_UT1(JD_UTC,Delta_UT1);
//print_matrix_2d_ld(JD_UT1);

//Finding Julian Date
vector<vector<long double> > JD_TT;
long double delta_AT=37;
JD_TT = Get_JD_TT(JD_UTC,delta_AT);
//print_matrix_2d_ld(JD_TT);

//Finding Nominal Data in ECEF frame from TOD frame
vector<vector<long double> > Data_nominal_ECEF(Data_Measured.size(),vector<long double>(7));
Data_nominal_ECEF = Nominal_TOD_to_ECEF(Data_Measured, Data_nominal_TOD, Polar_angle, LOD, d_psi, d_epsilon, JD_UT1, JD_TT);
//print_matrix_2d_ld(Data_nominal_ECEF);

// Converting site coordinates into ECEF Frame
vector<vector<long double>> svsite_ECEF;
svsite_ECEF=siteSVECEF(site_alti,site_lat,site_lon);

vector<vector<long double> > Final_Data;
Final_Data=Final_Data_Matrix(Data_Measured, Data_nominal_ECEF, svsite_ECEF, site_alti, site_lat, site_lon, JD_UT1);
//print_matrix_2d_ld(Final_Data);

// Printing residuals
//for(int i=0; i< Final_Data.size(); i++)
//    cout<<fixed<<setprecision(10)<<abs(Final_Data[i][8]-Final_Data[i][7])<<endl;

//////////////////////////////////////////////////////////////////////////////////////////
//
// Calculating value of Weight matrix
vector<vector<long double> > X_nominal(1,vector<long double> (6));
    for(int i=0; i<X_nominal[0].size(); i++)
    {
        X_nominal[0][i] = Nominal_Data_temp[0][i+3]*1000;
    }

// Calculation of Weight Matrix
long double weight_of_obs[]={1};
vector<vector<long double> > W(no_of_obs,vector<long double>(no_of_obs));
    for(int i=0; i<W.size();i++)
    {
        for(int j=0; j<W[0].size();j++)
        {
            if(i==j)
                W[i][j]=weight_of_obs[i];
            else
                W[i][j]=0;
        }
    }

vector<vector<long double> > Final_Data_row(1,vector<long double> (Final_Data[0].size()));
vector<vector<long double> > At_W_A1(6,vector<long double> (6,0));
vector<vector<long double> > At_W_b1(6,vector<long double> (1,0));
long double Delta_i = 1.00001;

    for(int i=0; i<Final_Data.size(); i++)
    {
        for(int j=0; j<Final_Data[0].size();j++)
        {
            Final_Data_row[0][j] = Final_Data[i][j];
        }
        At_W_A1 = Add(At_W_A1, At_W_A(Final_Data_row, Delta_i, no_of_obs, W));
        At_W_b1 = Add(At_W_b1, At_W_b(Final_Data_row, Delta_i, no_of_obs, W));
    }


vector<vector<long double> > delta_x;
delta_x = multiplication(Inverse_Matrix(At_W_A1),At_W_b1);
cout<<"Delta X is given by:"<< endl;
print_matrix_2d_ld(delta_x);

long double RMS;
RMS = RMS_eval(Final_Data, no_of_obs, W);
cout<<"RMS Value is: "<<RMS<<endl;

vector<vector<long double> > X_nominal_new(1,vector<long double> (6));
X_nominal_new = Add(X_nominal,transpose(delta_x));

cout<<"Old nominal"<<endl;
print_matrix_2d_ld(X_nominal);
cout<<"New nominal"<<endl;
print_matrix_2d_ld(X_nominal_new);

return 0;
}

