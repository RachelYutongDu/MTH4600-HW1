
// This file includes some standard C library functions and
// declares the functions found in "Definitions.h".
#include "Declarations.h"

// MTUniform is the Mersenne twister, MWCUniform is the multiply with carry
//    generator that we discussed in class.


int main() {

   double MWCUniform (unsigned int);

   // Declare the variables
   int j, k, m, n;
   double U[5], Z, p, q, mu, sigma, M[5];
   int *X;

   // Allocate space for X[0], ... ,X[63999] and initialize each to 0
   X = (int *) calloc (64000,sizeof (int));

   //seed the random number generator
   n = GetInteger ("How many simulations should I do?...");

   // Start the clock
   Time();

   //Begin simulations
   for (k = 1; k <= n; k++) {

      //Generate a realization of (U1,U2,U3) and then(m1,m2,m3)
      //May use MWCUniform, MTUniform as the RNG
      for (j = 1; j <=3; j++) {

         U[j] = MTUniform(0);

	 // M[j] is integer and 0 <= M[j] <= 39
	 M[j] = int (U[j]*40);  
      }
      
      //Compute the value m for (m1,m2,m3), where 0 <= m <= 63999
      m = M[1]*40*40 + M[2]*40 + M[3];

      // Increase the appropriate counter Xm by 1
      X[m] ++;

}


   // Now each X[m] should be binomial(n,p), where p = 1/64000
   p = 1.0/64000.0;
   q = 1.0 - p;
   mu = n*p;
   sigma = sqrt(n*p*q);


   // If n is very large, the numbers Zm = (X[m] - mu)/sigma should be approximately Normal(0,1) by CLT
   for (m = 0; m <= 63999; m++) {
      Z = (X[m] -mu)/sigma;
      NormalHistogram(Z, 40,0);
   }

   //Create the Tax files for viewing
   NormalHistogram(0,40,1);

   printf("Computations took %8.2f seconds.\n", Time());

   Exit ();

}

#include "Definitions.h"



////////////////////////////////////////////////////////////////////////////////
// The following code implements the Multiply-With-Carry RNG.
////////////////////////////////////////////////////////////////////////////////
double MWCUniform (unsigned int seed) {

   // Static variables retain their values between function calls.
   static unsigned int a1, a0, n1, n0, c1, c0, initialized=0;

   unsigned int w1, w0, x1, x0, y1, y0, z1, z0, carry, N;

   // This function is called only be MWCUniform and follows below.
   void Split (unsigned int, unsigned int *, unsigned int *);

   // Seed the MWC function when it is passed a non-zero seed or it is not yet
   //    initialized.
   // Represent a, n, and c as with pairs of 16-bit numbers:
   //  a = a1 * 2^16 + a0 (the multiplier)
   //  n = n1 * 2^16 + n0 (initially the seed)
   //  c = c1 * 2^16 + c0 (the increment, initially 1)
   if (seed || !initialized) {
      seed = (seed ? seed : 1);
      Split (4294967118, &a1, &a0);
      Split (seed, &n1, &n0);
      c0 = 1;                                 // first 16 bits of "c"
      c1 = 0;                                 // second 16 bits
      initialized = 1;
   }

   // Generate the next number in the sequence.

   // Compute w1, w0, where a0 * n0 = w1 * 2^16 + w0.
   Split (a0 * n0, &w1, &w0);

   // Compute x1, x0, where a1 * n0 = x1 * 2^16 + x0.
   Split (a1 * n0, &x1, &x0);

   // Compute y1, y0, where a0 * n1 = y1 * 2^16 + y0.
   Split (a0 * n1, &y1, &y0);

   // Compute z1, z0, where a1 * n1 = z1 * 2^16 + z0.
   Split (a1 * n1, &z1, &z0);

   // Compute next n0, n1, c0, c1, where an + c = c1 * 2^48 + c0 * 2^32 + n1 * 2^16 + n0.

   // First compute n0 and whatever is carried.
   Split (w0 + c0, &carry, &n0);

   // Compute n1 and whatever is carried.
   Split (carry + w1 + x0 + y0 + c1, &carry, &n1);

   // Compute c0 and whatever is carried.
   Split (carry + x1 + y1 + z0, &carry, &c0);

   // Compute c1; nothing is carried.
   c1 = z1 + carry;

   // Re-assemble n1 and n0.
   N = (n1 << 16) + n0;

   return ((N + 0.5) / 4294967296.0);

}

// This function computes 16-bit numbers x1 and x0 such that x = x1 * 2^16 + x0.
// Pointers *x1 and *x0 are used because a function can return only one value.
void Split (unsigned int x, unsigned int *x1, unsigned int *x0) {

   // First 16 bits of x.
   *x0 = x & 0xffff;

   // Second 16 bits of x.
   *x1 = x >> 16;

   return;

}










