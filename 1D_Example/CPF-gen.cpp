#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "shift.h"
#include "init.h"
using namespace std;


#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

 using namespace std; 
 struct Ran {
  /*static*/ float y2;
  /*static*/ int use_last/* = 0*/;

  /*static*/ unsigned long mt[N]; /* the array for the state vector  */
  /*static*/ int mti/*=N+1*/; /* mti==N+1 means mt[N] is not initialized */

  /* initializes mt[N] with a seed */
  void init_genrand(unsigned long s)
  {
      mt[0]= s & 0xffffffffUL;
      for (mti=1; mti<N; mti++) {
    mt[mti] = 
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
      }
  }

  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  /* slight change for C++, 2004/2/26 */
  void init_by_array(unsigned long init_key[], int key_length)
  {
      int i, j, k;
      init_genrand(19650218UL);
      i=1; j=0;
      k = (N>key_length ? N : key_length);
      for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
      + init_key[j] + j; /* non linear */
    mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    i++; j++;
    if (i>=N) { mt[0] = mt[N-1]; i=1; }
    if (j>=key_length) j=0;
      }
      for (k=N-1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
      - i; /* non linear */
    mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    i++;
    if (i>=N) { mt[0] = mt[N-1]; i=1; }
      }

      mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
  }

  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long genrand_int32(void)
  {
      unsigned long y;
      static unsigned long mag01[2]={0x0UL, MATRIX_A};
      /* mag01[x] = x * MATRIX_A  for x=0,1 */

      if (mti >= N) { /* generate N words at one time */
    int kk;

    if (mti == N+1)   /* if init_genrand() has not been called, */
        init_genrand(5489UL); /* a default initial seed is used */

    for (kk=0;kk<N-M;kk++) {
        y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<N-1;kk++) {
        y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mti = 0;
      }
    
      y = mt[mti++];

      /* Tempering */
      y ^= (y >> 11);
      y ^= (y << 7) & 0x9d2c5680UL;
      y ^= (y << 15) & 0xefc60000UL;
      y ^= (y >> 18);

      return y;
  }

  /* generates a random number on [0,0x7fffffff]-interval */
  long genrand_int31(void)
  {
      return (long)(genrand_int32()>>1);
  }

  /* generates a random number on [0,1]-real-interval */
  double genrand_real1(void)
  {
      return genrand_int32()*(1.0/4294967295.0); 
      /* divided by 2^32-1 */ 
  }

  /* generates a random number on [0,1)-real-interval */
  double genrand_real2(void)
  {
      return genrand_int32()*(1.0/4294967296.0); 
      /* divided by 2^32 */
  }

  /* generates a random number on (0,1)-real-interval */
  double genrand_real3(void)
  {
      return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
      /* divided by 2^32 */
  }

  /* generates a random number on [0,1) with 53-bit resolution*/
  double genrand_res53(void) 
  { 
      unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
      return(a*67108864.0+b)*(1.0/9007199254740992.0); 
  } 
  /* These real versions are due to Isaku Wada, 2002/01/09 added */
  
  //custom interface part
  Ran() {use_last=0;mti=N+1;/*dummy constructor*/}
  
  void seed(int ISEED) {init_genrand(ISEED);}
  
  Ran(int ISEED) {use_last=0;mti=N+1;seed(ISEED);}
  
  float flt(){return genrand_real3();}
  /* boxmuller.c           Implements the Polar form of the Box-Muller Transformation
   * (c) Copyright 1994, Everett F. Carter Jr.
   * Permission is granted by the author to use
   * this software for any application provided this
   * copyright notice is preserved.
   */
  float box_muller(float m, float s)  /* normal random variate generator */
  {               /* mean m, standard deviation s */
    float x1, x2, w, y1;
    if (use_last)           /* use value from previous call */
    {
      y1 = y2;
      use_last = 0;
    }
    else
    {
      do {
        x1 = 2.0 * flt() - 1.0;
        x2 = 2.0 * flt() - 1.0;
        w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0 );

      w = sqrt( (-2.0 * log( w ) ) / w );
      y1 = x1 * w;
      y2 = x2 * w;
      use_last = 1;
    }

    return( m + y1 * s );
  }
  inline float gauss_distr(){return box_muller(0.0,1.0);}
};

 // return a uniformly distributed random number
double uniformRandom()
{
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}
 // return a normally distributed random number
double normalRandom()
{
  double u1=uniformRandom();
  double u2=uniformRandom();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

/*
void simulation(double xs, double xa, double dt, Ran *eran, ofstream& file)
{
  

  double milestone = 1.0;
  double milestone_step = 0.05;

  bool abs = false;
  do {
    double xtshift = 8*(1 - x) - 7.493877551020408*cos(2.754*(-0.24332 + x))*(-y + 0.17006802721088435*sin(2.754*(-0.24332 + x))) - 21.240280675137495*cos(14*(0.23166 + x))*(-z + 0.1896453631708705*sin(14*(0.23166 + x))) - 
   48.81598667776852*cos(14.657*(0.47268 + x))*(-p + 0.4163197335553705*sin(14.657*(0.47268 + x))) - 23.112727272727273*cos(7.945*(0.47501 + x))*(-w + 0.36363636363636365*sin(7.945*(0.47501 + x)));
    double ytshift = 16*(-y + 0.17006802721088435*sin(2.754*(-0.24332 + x)));
    double ztshift = 8*(-z + 0.1896453631708705*sin(14*(0.23166 + x)));
    double wtshift = 8*(-w + 0.36363636363636365*sin(7.945*(0.47501 + x)));
    double ptshift = 8*(-p + 0.4163197335553705*sin(14.657*(0.47268 + x)));
    double dWx = sqrt(dt)* (eran->gauss_distr());
    double dWy = sqrt(dt)* (eran->gauss_distr());
    double dWz = sqrt(dt)* (eran->gauss_distr());
    double dWw = sqrt(dt)* (eran->gauss_distr());
    double dWp = sqrt(dt)* (eran->gauss_distr());
    x = x + xtshift * dt + sqrt(2.0)*dWx;  //diffusion and Wiener process
    y = y + ytshift * dt + sqrt(2.0)*dWy;
    z = z + ztshift * dt + sqrt(2.0)*dWz;
    w = w + wtshift * dt + sqrt(2.0)*dWw;
    p = p + ptshift * dt + sqrt(2.0)*dWp;

    //if (x <= xs){ //passage happened
      //x = xs; 

    //reflecting boundary
    if (x > xa) x = 2 * xa - x;

    //cout << "\n" << dt * ncount << "\t" << x << "\t" << milestone;
    ncount++;

    if (x<milestone){
     	//file << dt * ncount << "\t";
        cout << dt * ncount << "\t" << flush;
     	milestone -= milestone_step;
    }
  
  } while (x > xs);
  file << "\n";
  
  return;
}
*/

void timestep(int n, double* x, double* dx, double* dWt, double dt, Ran* eran){
    //Calculate Weiner process step dWt
    for(int i=0; i<n; i++){
        dWt[i]=sqrt(dt)* (eran->gauss_distr());
    }

    //update position based on Ito expression
    shift(dWt, dt, x, dx);

    //Increment x
    for(int i=0; i<n; i++){
        x[i] += dx[i];
    }
}

void fpt(int n, double xs, double xa, double dt, Ran *eran, ofstream& flux_file, ofstream& cpf_file){
    double* x = new double[n];
    double* dx = new double[n];
    double* dWt = new double[n];

    unsigned long long ncount = 0;
    unsigned long long i = 0;

    while (i<100000000){
        //Initialize position of particle (x[0]=xs)
//         double* p = new double[n-1];
//         for (int i=0; i<n-1; i++){
//             p[i] = (eran->gauss_distr());
//         }
        init(xs, x);
//         delete[] p;


        //Run simulation
        do{
            timestep(n, x, dx, dWt, dt, eran);
            i++;
            if (i%50==0) cpf_file << x[0] << "\n";
        } while (x[0] < xa);
        ncount++;
        flux_file << dt * i << "\t" << ncount << "\n";
    }

    delete[] x;
    delete[] dx;
    delete[] dWt;
}

int main(int narg, char** arg)
{
    ofstream flux_file, cpf_file;
    flux_file.open("flux.dat", ios::out);
    cpf_file.open("cpf.dat", ios::out);

    Ran eran(0);
    fpt(1, 0.0, 1.0, 0.0001, &eran, flux_file, cpf_file);

    flux_file.close();
    cpf_file.close();

    return 0;
}