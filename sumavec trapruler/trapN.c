//
//  main.c
//  trap mpi
//
//  Created by Luis Alberto Sanchez Moreno on 9/03/18.
//  Copyright © 2018 Luis Alberto Sanchez Moreno. All rights reserved.
//

/*

 * Usage:    mpirun -np <number of processes> ./trap < number of trapezoids>
 *
 * Algorithm:
 *    1.  Each process calculates "its" interval of
 *        integration.
 *    2.  Each process estimates the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
 *
 * Note:  f(x) is all hardwired to x*x.
 *
 */
#include <stdio.h>
#include <math.h>

/* We'll be using MPI routines, definitions, etc. */
#include <mpi.h>

void Get_data(int p, int my_rank, double* a_p, double* b_p, int* n_p);

double Trap(double local_a, double local_b, int local_n,
           double h);    /* Calculate local area  */

double f(double x); /* function we're integrating */

double Get_max(double my_val, int my_rank, int p);


int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    double      a;         /* Left endpoint             */
    double      b;         /* Right endpoint            */
    int         n;         /* Number of trapezoids      */
    double      h;         /* Trapezoid base length     */
    double      local_a;   /* Left endpoint my process  */
    double      local_b;   /* Right endpoint my process */
    int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
    double      my_int;    /* Integral over my interval */
    double      r_int;     /* Received integral         */
    double      total;     /* Total integral            */
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;
    double      start, finish, my_elapsed, elapsed;

    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    Get_data(p, my_rank, &a, &b, &n);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */

    /* Length of each process' interval of
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    my_int = Trap(local_a, local_b, local_n, h);

    /* Add up the areas calculated by each process */
    if (my_rank == 0) {
        total = my_int;
        for (source = 1; source < p; source++) {
            MPI_Recv(&r_int, 1, MPI_DOUBLE, source, tag,
                MPI_COMM_WORLD, &status);
            total += r_int;
        }
    } else {  
        MPI_Send(&my_int, 1, MPI_DOUBLE, dest,
            tag, MPI_COMM_WORLD);
    }

    finish = MPI_Wtime();
    my_elapsed = finish-start;
    elapsed = Get_max(my_elapsed, my_rank, p);

    /* Print the result */
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n",
            n);
        printf("of the integral from %f to %f = %.15f\n",
            a, b, total);
        printf("Elapsed time = %e seconds\n", elapsed);
    }

    /* Shut down MPI */
    MPI_Finalize();

    return 0;
} /*  main  */

/*------------------------------------------------------------------
 * Function:     Get_data
 * Purpose:      Read in the data on process 0 and send to other
 *               processes
 * Input args:   p, my_rank
 * Output args:  a_p, b_p, n_p
 */
void Get_data(int p, int my_rank, double* a_p, double* b_p, int* n_p) {
   int        q;
   MPI_Status status;

   if (my_rank == 0) {
      printf("Enter a, b, and n\n");
      scanf("%lf %lf %d", a_p, b_p, n_p);

      for (q = 1; q < p; q++) {
         MPI_Send(a_p, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD);
         MPI_Send(b_p, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD);
         MPI_Send(n_p, 1, MPI_INT, q, 0, MPI_COMM_WORLD);
      }
   } else {
      MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
   }
}  /* Get_data */

/*------------------------------------------------------------------
 * Function:     Trap
 * Purpose:      Estimate a definite integral using the trapezoidal
 *               rule
 * Input args:   local_a (my left endpoint)
 *               local_b (my right endpoint)
 *               local_n (my number of trapezoids)
 *               h (stepsize = length of base of trapezoids)
 * Return val:   Trapezoidal rule estimate of integral from
 *               local_a to local_b
 */
double Trap(
          double  local_a   /* in */,
          double  local_b   /* in */,
          int     local_n   /* in */,
          double  h         /* in */) {
    double integral;   /* Store result in integral  */
    double x;
    int i;

    integral = (f(local_a) + f(local_b))/2.0;
    for (i = 1; i <= local_n-1; i++) {
        x = local_a + i*h;
        integral += f(x);
    }
    integral = integral*h;

    return integral;
} /*  Trap  */

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
    double return_val;

    return_val = sin(exp(x));

    return return_val;
} /* f */


/*------------------------------------------------------------------
 * Function:    Get_max
 * Purpose:     Compute the global max of a collection of doubles --
 *              one per process
 * Input args:  all
 * Return val:  Global max on process 0, 0 other processes.
 */
double Get_max(double my_val, int my_rank, int p) {
   int q;
   double global_max = 0.0, val;
   MPI_Status status;

   if (my_rank == 0) {
      global_max = my_val;
      for (q = 1; q < p; q++) {
         MPI_Recv(&val, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, &status);
         if (val > global_max) global_max = val;
      }
   } else {
      MPI_Send(&my_val, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   }

   return global_max;
}  /* Get_max */