//https://people.math.sc.edu/Burkardt/cpp_src/gsl_test/gsl_test.cpp

# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iostream>
# include <iomanip>

using namespace std;

# include <gsl/gsl_eigen.h>
# include <gsl/gsl_fft_complex.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_poly.h>
# include <gsl/gsl_qrng.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_sf_dawson.h>
# include <gsl/gsl_multiroots.h>
# include <gsl/gsl_vector.h>

int main ( );
void bessel_j0_values ( int *n_data, double *x, double *fx );
void bessel_j1_values ( int *n_data, double *x, double *fx );
void dawson_values ( int *n_data, double *x, double *fx );
void gsl_eigen_nonsymm_test ( );
void gsl_fft_complex_test ( );
double gsl_integration_f ( double x, void *params );
void gsl_integration_qng_test ( );
void gsl_linalg_LU_test ( );
double gsl_monte_f ( double *x, size_t dim, void *params );
void gsl_multiroot_fsolver_test ( );
void gsl_multiroot_print_state ( size_t iter, gsl_multiroot_fsolver *s );
void gsl_poly_eval_test ( );
void gsl_qrng_niederreiter_2_test ( );
void gsl_qrng_sobol_test ( );
void gsl_sf_bessel_J0_test ( );
void gsl_sf_bessel_J1_test ( );
void gsl_sf_coupling_3j_test ( );
void gsl_sf_coupling_6j_test ( );
void gsl_sf_coupling_9j_test ( );
void gsl_sf_dawson_test ( );
int rosenbrock_f ( const gsl_vector *x, void *params, gsl_vector *f );
void three_j_values ( int *n_data, double *j1, double *j2, double *j3, 
  double *m1, double *m2, double *m3, double *fx );
void six_j_values ( int *n_data, double *j1, double *j2, double *j3, 
  double *j4, double *j5, double *j6, double *fx );
void nine_j_values ( int *n_data, double *j1, double *j2, double *j3, 
  double *j4, double *j5, double *j6, double *j7, double *j8, double *j9,
  double *fx );
void timestamp ( );
//
//  Rosenbrock parameter structure.
//
struct rparams
{
  double a;
  double b;
};

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the GSL tests.
//
//  Modified:
//
//    09 March 2017
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "GSL_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test GSL, the GNU Scientific Library.\n";

  gsl_eigen_nonsymm_test ( );

  gsl_fft_complex_test ( );

  gsl_integration_qng_test ( );

  gsl_linalg_LU_test ( );

  gsl_multiroot_fsolver_test ( );

  gsl_poly_eval_test ( );

  gsl_qrng_niederreiter_2_test ( );

  gsl_qrng_sobol_test ( );

  gsl_sf_bessel_J0_test ( );
  gsl_sf_bessel_J1_test ( );

  gsl_sf_coupling_3j_test ( );
  gsl_sf_coupling_6j_test ( );
  gsl_sf_coupling_9j_test ( );

  gsl_sf_dawson_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "GSL_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void bessel_j0_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_VALUES returns some values of the J0 Bessel function.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double fx_vec[N_MAX] = {
    -0.1775968E+00, -0.3971498E+00, -0.2600520E+00,  0.2238908E+00, 
     0.7651976E+00,  1.0000000E+00,  0.7651977E+00,  0.2238908E+00, 
    -0.2600520E+00, -0.3971498E+00, -0.1775968E+00,  0.1506453E+00, 
     0.3000793E+00,  0.1716508E+00, -0.0903336E+00, -0.2459358E+00, 
    -0.1711903E+00,  0.0476893E+00,  0.2069261E+00,  0.1710735E+00, 
    -0.0142245E+00 };
  double x_vec[N_MAX] = {
    -5.0E+00, -4.0E+00, -3.0E+00, -2.0E+00, 
    -1.0E+00,  0.0E+00,  1.0E+00,  2.0E+00, 
     3.0E+00,  4.0E+00,  5.0E+00,  6.0E+00, 
     7.0E+00,  8.0E+00,  9.0E+00, 10.0E+00, 
    11.0E+00, 12.0E+00, 13.0E+00, 14.0E+00, 
    15.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80*

void bessel_j1_values ( int *n_data, double *x, double *fx )

//****************************************************************************80*
//
//  Purpose:
//
//    BESSEL_J1_VALUES returns some values of the J1 Bessel function.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double fx_vec[N_MAX] = {
    0.3275791E+00,  0.0660433E+00, -0.3390590E+00, -0.5767248E+00, 
   -0.4400506E+00,  0.0000000E+00,  0.4400506E+00,  0.5767248E+00, 
    0.3390590E+00, -0.0660433E+00, -0.3275791E+00, -0.2766839E+00, 
   -0.0046828E+00,  0.2346364E+00,  0.2453118E+00,  0.0434728E+00, 
   -0.1767853E+00, -0.2234471E+00, -0.0703181E+00,  0.1333752E+00, 
    0.2051040E+00 };
  double x_vec[N_MAX] = { 
    -5.0E+00, -4.0E+00, -3.0E+00, -2.0E+00, 
    -1.0E+00,  0.0E+00,  1.0E+00,  2.0E+00, 
     3.0E+00,  4.0E+00,  5.0E+00,  6.0E+00, 
     7.0E+00,  8.0E+00,  9.0E+00, 10.0E+00, 
    11.0E+00, 12.0E+00, 13.0E+00, 14.0E+00, 
    15.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80*

void dawson_values ( int *n_data, double *x, double *fx )

//****************************************************************************80*
//
//  Purpose: 
//
//    DAWSON_VALUES returns some values of Dawson's integral.
//
//  Discussion:
//
//    D(X) = exp ( -X^2 ) * Integral ( 0 <= Y <= X ) exp ( Y^2 ) dY
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 22

  double fx_vec[N_MAX] = { 
    0.0000000000E+00, 0.0993359924E+00, 0.1947510334E+00, 0.2826316650E+00, 
    0.3599434819E+00, 0.4244363835E+00, 0.4747632037E+00, 0.5105040576E+00, 
    0.5321017071E+00, 0.5407243187E+00, 0.5410442246E+00, 0.5380795069E+00, 
    0.5262066800E+00, 0.5072734964E+00, 0.4833975174E+00, 0.4565072375E+00, 
    0.4282490711E+00, 0.3999398943E+00, 0.3725593490E+00, 0.3467727691E+00, 
    0.3229743193E+00, 0.3013403889E+00 };
  double x_vec[N_MAX] = { 
    0.0E+00, 0.1E+00, 0.2E+00,          0.3E+00, 
    0.4E+00, 0.5E+00, 0.6E+00,          0.7E+00, 
    0.8E+00, 0.9E+00, 0.9241388730E+00, 1.0E+00, 
    1.1E+00, 1.2E+00, 1.3E+00,          1.4E+00, 
    1.5E+00, 1.6E+00, 1.7E+00,          1.8E+00, 
    1.9E+00, 2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0E+00;
    *fx = 0.0E+00;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }
  return;

# undef N_MAX
}
//****************************************************************************80

void gsl_eigen_nonsymm_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_EIGEN_NONSYMM_TEST tests GSL_EIGEN_NONSYMM.
//
//  Modified:
//
//    09 March 2017
//
//  Reference:
//
//    Mark Gelassi, Jim Davies, James Tyler, Bryan Gough,
//    Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi,
//    GNU Scientific Library Reference Manual.
//
{
  double data[] = { -1.0, 1.0, -1.0, 1.0,
                    -8.0, 4.0, -2.0, 1.0,
                    27.0, 9.0, 3.0, 1.0,
                    64.0, 16.0, 4.0, 1.0 };
//
//  Convert C array to GSL matrix.
//
  gsl_matrix_view m = gsl_matrix_view_array ( data, 4, 4 );
//
//  Allocate GSL matrix and vector space for eigenvectors and eigenvalues.
//
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc ( 4, 4 );
  gsl_vector_complex *eval = gsl_vector_complex_alloc ( 4 );
//
// Set up workspace W.
//
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc ( 4 );

  cout << "\n";
  cout << "GSL_EIGEN_NONSYMM_TEST:\n";
  cout << "  GSL_EIGEN_NONYSMM computes the eigenvalues and eigenvectors\n";
  cout << "  of a nonsymmetric matrix.\n";
//
//  Compute the eigenvalues and eigenvectors.
//
  gsl_eigen_nonsymmv ( &m.matrix, eval, evec, w );
//
//  Free the workspace memory.
//
  gsl_eigen_nonsymmv_free ( w );
//
//  Sort the results.
//
  gsl_eigen_nonsymmv_sort ( eval, evec, GSL_EIGEN_SORT_ABS_DESC );
//
//  Print the results.
//  
  for ( int i = 0; i < 4; i++ )
  {
    gsl_complex eval_i = gsl_vector_complex_get ( eval, i );
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column ( evec, i );

    printf ( "\n" );
    printf ( "  Eigenvalue(%d) = %g + %gi\n", 
      i, GSL_REAL ( eval_i ), GSL_IMAG ( eval_i ) );

    printf ( "\n" );
    printf ( "  Eigenvector(%d) = \n", i );
    for ( int j = 0; j < 4; j++ )
    {
      gsl_complex z = gsl_vector_complex_get ( &evec_i.vector, j );
      printf ( "  %g + %gi\n", GSL_REAL ( z ), GSL_IMAG ( z ) );
    }
  }
//
//  Free memory.
//
  gsl_vector_complex_free ( eval );
  gsl_matrix_complex_free ( evec );

  return;
}
//****************************************************************************80

void gsl_fft_complex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_FFT_COMPLEX_TEST tests GSL_FFT_COMPLEX.
//
//  Modified:
//
//    09 March 2017
//
//  Reference:
//
//    Mark Gelassi, Jim Davies, James Tyler, Bryan Gough,
//    Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi,
//    GNU Scientific Library Reference Manual.
//
{ 
  double data[2*128];
  int i;

  cout << "\n";
  cout << "GSL_FFT_COMPLEX_TEST:\n";
  cout << "  GSL_FFT_COMPLEX computes the fast Fourier transform\n";
  cout << "  of a complex vector of data.\n";

  for ( i = 0; i < 128; i++ )
  {
    data[0+i*2] = 0.0;
    data[1+i*2] = 0.0;
  }
  data[0+0*2] = 1.0;
//
//  Insert a pulse.
//  Make it symmetric so the transform will actually be real.
//
  for ( i = 1; i <= 10; i++ )
  {
    data[0+i*2] = 1.0;
    data[0+(128-i)*2] = 1.0;
  }
//
//  Print the input.
//
  printf ( "\n" );
  printf ( "  Input data (a symmetric pulse):\n" );
  printf ( "\n" );
  for ( i = 0; i < 128; i++ )
  {
    printf ( "%3d %e %e\n", i, data[0+i*2], data[1+i*2] );
  }
//
//  Compute the FFT.
//
  gsl_fft_complex_radix2_forward ( data, 1, 128 );
//
//  Print the scaled output.
//
  printf ( "\n" );
  printf ( "  Output data:\n" );
  printf ( "\n" );
  for ( i = 0; i < 128; i++ )
  {
    printf ( "%3d %e %e\n", 
      i, data[0+i*2] / sqrt ( 128.0 ), data[1+i*2] / sqrt ( 128.0 ) );
  }

  return;
}
//****************************************************************************80

double gsl_integration_f ( double x, void *params )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_INTEGRATION_F evaluates the integrand for quadrature test.
//
//  Discussion:
//
//    Integral ( 0 <= x <= 1 ) 2/(2+sin(10*pi*x)) = 2/sqrt(3).
//
//  Modified:
//
//    09 March 2017
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Gelassi, Jim Davies, James Tyler, Bryan Gough,
//    Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi,
//    GNU Scientific Library Reference Manual.
//
//  Parameters:
//
//    Input, double *X, the evaluation point.
//
//    Input, void *PARAMS, optional parameters.
//
//    Output, double GSL_INTEGRATION_F, the value of the integrand at X.
//
{
  double value;

  value = 2.0 / ( 2.0 + sin ( 10.0 * M_PI * x ) );
  
  return value;
}
//****************************************************************************80

void gsl_integration_qng_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_INTEGRATION_QNG_TEST tests GSL_INTEGRATION_QNG.
//
//  Modified:
//
//    09 March 2017
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Gelassi, Jim Davies, James Tyler, Bryan Gough,
//    Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi,
//    GNU Scientific Library Reference Manual.
//
{
  double a;
  double abserr;
  double abserr_exact;
  double b;
  double epsabs;
  double epsrel;
  int error_code;
  gsl_function F;
  size_t neval;
  double result;
  double result_exact;

  cout << "\n";
  cout << "GSL_INTEGRATION_QNG_TEST:\n";
  cout << "  GSL_INTEGRATION_QNG uses the Gauss-Kronrod rule to estimate\n";
  cout << "  the integral of f(x) over [a,b].\n";

  F.function = &gsl_integration_f;
  F.params = NULL;
  a = 0.0;
  b = 1.0;
//
//  If we reduce EPSABS much more, then the code returns an error.
//
  epsabs = 1.0E-01;
  epsrel = 0.0E+00;

  error_code = gsl_integration_qng ( &F, a, b, epsabs, epsrel, &result, 
    &abserr, &neval );

  if ( error_code != 0 )
  {
    cout << "\n";
    cout << "  GSL_INTEGRATION_QNG returned error code = " << error_code << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  Integral estimate = " << result 
         << ", error estimate = " << abserr << "\n";
    cout << "  Number of function evaluations was " << neval << "\n";
    result_exact = 2.0 / sqrt ( 3.0 );
    abserr_exact = fabs ( result - result_exact );
    cout << "  Exact integral = " << result_exact 
         << ", Exact error = " << abserr_exact << "\n";
  }

  return;
}
//****************************************************************************80

void gsl_linalg_LU_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_LINAL_LU_TEST tests GSL_LINALG_LU.
//
//  Modified:
//
//    09 March 2017
//
//  Reference:
//
//    Mark Gelassi, Jim Davies, James Tyler, Bryan Gough,
//    Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi,
//    GNU Scientific Library Reference Manual.
//
{
  double a_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };

  double b_data[] = { 1.0, 2.0, 3.0, 4.0 };
//
//  Create GSL_MATRIX M.
//
  gsl_matrix_view m = gsl_matrix_view_array ( a_data, 4, 4 );
//
//  Create GSL_VECtOR B.
//
  gsl_vector_view b = gsl_vector_view_array ( b_data, 4 );
//
//  Allocate GSL_VECTOR X (solution).
//
  gsl_vector *x = gsl_vector_alloc ( 4 );
  
  int s;
//
//  Allocate GSL_PERMUTATION P (pivot vector).
//
  gsl_permutation *p = gsl_permutation_alloc ( 4 );

  cout << "\n";
  cout << "GSL_LINALG_LU_TEST\n";
  cout << "  GSL_LINALG_LU_DECOMP computes LU factors of a matrix A.\n";
  cout << "  GSL_LINALG_LU_SOLVE solves linear systems A*x=b.\n";
//
//  LU factor the matrix M.
//
  gsl_linalg_LU_decomp ( &m.matrix, p, &s );
//
//  Solve M*x=b.
//
  gsl_linalg_LU_solve ( &m.matrix, p, &b.vector, x );
//
//  Display the solution.
//
  printf ( "\n" );
  printf ( "  Computed solution vector X: \n" );
  printf ( "\n" );

  gsl_vector_fprintf ( stdout, x, "%g" );
//
//  Free memory.
//
  gsl_permutation_free ( p );
  gsl_vector_free ( x );

  return;
}
//****************************************************************************80

double gsl_monte_f ( double *x, size_t dim_num, void *params )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_MONTE_F evaluates the integrand for a Monte Carlo test.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Gelassi, Jim Davies, James Tyler, Bryan Gough,
//    Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi,
//    GNU Scientific Library Reference Manual.
//
//  Parameters:
//
//    Input, double *X, the evaluation point.
//
//    Input, size_t DIM, the spatial dimension.
//
//    Input, void *PARAMS, optional parameters.
//
//    Output, double GSL_MONTE_F, the value of the integrand at X.
//
{
  double c;
  double den;
  int dim;
  double value;

  c = 1.0 / pow ( 2.0 * M_PI, dim_num );

  den = 1.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    den = den * cos ( x[dim] );
  }

  den = 1.0 - den;

  if ( den == 0.0 )
  {
    value = 100000.0;
  }
  else
  {
    value = c / den;
  }
  
  return value;
}
//****************************************************************************80

void gsl_multiroot_fsolver_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_MULTIROOT_FSOLVER_TEST tests the N-dimensional root finder.
//
//  Modified:
//
//    05 August 2005.
//
//  Author:
//
//    John Burkardt
//
{
  size_t i;
  size_t iter = 0;
  const size_t n = 2;
  struct rparams p = { 1.0, 10.0 };
  gsl_multiroot_fsolver *s;
  int status;
  const gsl_multiroot_fsolver_type *T;
  gsl_vector *x;
  double x_init[2] = { -10.0, -5.0 };

  cout << "\n";
  cout << "GSL_MULTIROOT_FSOLVER_TEST:\n";
  cout << "  Demonstrate the ability to find a root of a\n";
  cout << "  set of nonlinear equations.\n";
  cout << "\n";
  cout << "  In this case, we have two functions in two unknowns,\n";
  cout << "  and the only root is X = (1,1).\n";
  cout << "\n";

  gsl_multiroot_function f = { &rosenbrock_f, n, &p };
  x = gsl_vector_alloc ( n );

  gsl_vector_set ( x, 0, x_init[0] );
  gsl_vector_set ( x, 1, x_init[1] );

  T = gsl_multiroot_fsolver_hybrids;

  s = gsl_multiroot_fsolver_alloc ( T, 2 );

  gsl_multiroot_fsolver_set ( s, &f, x );

  gsl_multiroot_print_state ( iter, s );

  do 
  {
    iter++;

    status = gsl_multiroot_fsolver_iterate ( s );

    gsl_multiroot_print_state ( iter, s );

    if ( status )
    {
      break;
    }
    status = gsl_multiroot_test_residual ( s->f, 1.0E-07 );

  } while ( status == GSL_CONTINUE && iter < 1000 );

  printf ( "  status = %s\n", gsl_strerror ( status ) );
//
//  Free memory.
//
  gsl_multiroot_fsolver_free ( s );
  gsl_vector_free ( x );

  return;
}
//****************************************************************************80

void gsl_multiroot_print_state ( size_t iter, gsl_multiroot_fsolver *s )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_MULTIROOT_PRINT_STATE prints the state of the root finding iteration.
//
//  Modified:
//
//    05 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, SIZE_T ITER, the iteration counter.
//
//    Input, GSL_MULTIROOT_FSOLVER *S, the solver state.
//
{
  printf ( "  iter = %3lu x = %.3f %.3f f(x) = %.3e  %.3e\n",
    iter,
    gsl_vector_get ( s->x, 0 ),
    gsl_vector_get ( s->x, 1 ),
    gsl_vector_get ( s->f, 0 ),
    gsl_vector_get ( s->f, 1 ) );

  return;
}
//****************************************************************************80

void gsl_poly_eval_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_POLY_EVAL_TEST tests the GSL real polynomial evaluator.
//
//  Modified:
//
//    08 March 2017
//
//  Author:
//
//    John Burkardt
//
{
  double c[4] = {-18.0, 27.0, -10.0, 1.0 };
  int i;
  int len = 4;
  double p;
  double x;

  cout << "\n";
  cout << "GSL_POLY_EVAL_TEST:\n";
  cout << "  GSL_POLY_EVAL evaluates a real polynomial.\n";
  cout << "  p(x) = (x-1)*(x-3)*(x-6)\n";
  cout << "       = x^3 - 10x^2 + 27x - 18.\n";
  cout << "\n";
  for ( i = 0; i <= 20; i++ )
  {
    x = ( double ) i / 2.0;
    p = gsl_poly_eval ( c, len, x );
    cout << "  " << setw(8) << x
         << "  " << setw(8) << p << "\n";
  }

  return;
}
//****************************************************************************80

void gsl_qrng_niederreiter_2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_QRNG_NIEDERREITER_2_TEST tests the GSL Niederreiter sequence routine.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N 25

  int i;
  int j;
  double v[M];

  gsl_qrng * q = gsl_qrng_alloc ( gsl_qrng_niederreiter_2, M );

  cout << "\n";
  cout << "GSL_QRNG_NIEDERREITER_2_TEST:\n";
  cout << "  GSL_QRNG_ALLOC sets aside space for a sequence;\n";
  cout << "  GSL_QRNG_SOBOL requests the Niederreiter_2 sequence;\n";
  cout << "  GSL_QRNG_GET gets the next entry of the requested sequence;\n";
  cout << "\n";
  cout << "  Determine the first " << N << " points of the Niederreiter2\n";
  cout << "  quasi-random sequence in " << M << " dimensions.\n";
  cout << "\n";
  cout << "     I       X(I)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    gsl_qrng_get ( q, v );

    cout << setw(6)  << i << "  ";
    for ( j = 0; j < M; j++ )
    {
      cout << setw(12) << v[j] << "  ";
    }
    cout << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void gsl_qrng_sobol_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_QRNG_SOBOL_TEST tests the GSL Sobol sequence routine.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N 25

  int i;
  int j;
  double v[M];

  gsl_qrng * q = gsl_qrng_alloc ( gsl_qrng_sobol, M );

  cout << "\n";
  cout << "GSL_QRNG_SOBOL_TEST:\n";
  cout << "  GSL_QRNG_ALLOC sets aside space for a sequence;\n";
  cout << "  GSL_QRNG_SOBOL requests the Sobol sequence;\n";
  cout << "  GSL_QRNG_GET gets the next entry of the requested sequence;\n";
  cout << "\n";
  cout << "  Determine the first " << N << " points of the Sobol\n";
  cout << "  quasi-random sequence in " << M << " dimensions.\n";
  cout << "\n";
  cout << "     I       X(I)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    gsl_qrng_get ( q, v );

    cout << setw(6)  << i << "  ";
    for ( j = 0; j < M; j++ )
    {
      cout << setw(12) << v[j] << "  ";
    }
    cout << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void gsl_sf_bessel_J0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_SF_BESSEL_J0_TEST tests the GSL Bessel J0 evaluation routine.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "GSL_SF_BESSEL_J0_TEST:\n";
  cout << "  Evalute the J0 Bessel function.\n";
  cout << "\n";
  cout << "        X            Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( &n_data, &x, &fx1 );

    if ( n_data <= 0 )
    {
      break;
    }
    fx2 = gsl_sf_bessel_J0 ( x );
    cout                    << "  "
         << setw(12) << x   << "  "
         << setw(12) << fx1 << "  "
         << setw(12) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void gsl_sf_bessel_J1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_SF_BESSEL_J1_TEST tests the GSL Bessel J1 evaluation routine.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "GSL_SF_BESSEL_J1_TEST:\n";
  cout << "  Evalute the J1 Bessel function.\n";
  cout << "\n";
  cout << "        X            Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j1_values ( &n_data, &x, &fx1 );

    if ( n_data <= 0 )
    {
      break;
    }
    fx2 = gsl_sf_bessel_J1 ( x );
    cout                    << "  "
         << setw(12) << x   << "  "
         << setw(12) << fx1 << "  "
         << setw(12) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void gsl_sf_coupling_3j_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_SF_COUPLING_3J_TEST demonstrates GSL_SF_COUPLING_3J.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;
  int twoj1;
  int twoj2;
  int twoj3;
  int twom1;
  int twom2;
  int twom3;

  cout << "\n";
  cout << "GSL_SF_COUPLING_3J_TEST:\n";
  cout << "  GSL_SF_COUPLING_3J returns values of\n";
  cout << "  the Wigner 3J coefficient.\n";
  cout << "\n";
  cout << "      J1      J2      J3      M1      M2      M3        THREE_J\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    three_j_values ( &n_data, &j1, &j2, &j3, &m1, &m2, &m3, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "  " << setw(6) << j1
         << "  " << setw(6) << j2
         << "  " << setw(6) << j3
         << "  " << setw(6) << m1
         << "  " << setw(6) << m2
         << "  " << setw(6) << m3
         << "  " << setprecision(16) << setw(24) << fx << "\n";

    twoj1 = round ( 2.0 * j1 );
    twoj2 = round ( 2.0 * j2 );
    twoj3 = round ( 2.0 * j3 );
    twom1 = round ( 2.0 * m1 );
    twom2 = round ( 2.0 * m2 );
    twom3 = round ( 2.0 * m3 );

    fx2 = gsl_sf_coupling_3j ( twoj1, twoj2, twoj3, twom1, twom2, twom3 );

    cout << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << setprecision(16) << setw(24) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void gsl_sf_coupling_6j_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_SF_COUPLING_6J_TEST demonstrates GSL_SF_COUPLING_6J.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  int n_data;
  int twoj1;
  int twoj2;
  int twoj3;
  int twoj4;
  int twoj5;
  int twoj6;

  cout << "\n";
  cout << "GSL_SF_COUPLING_6J_TEST:\n";
  cout << "  GSL_SF_COUPLING_6J returns values of\n";
  cout << "  the Wigner 6J coefficient.\n";
  cout << "\n";
  cout << "      J1      J2      J3      J4      J5      J6        SIX_J\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    six_j_values ( &n_data, &j1, &j2, &j3, &j4, &j5, &j6, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "  " << setw(6) << j1
         << "  " << setw(6) << j2
         << "  " << setw(6) << j3
         << "  " << setw(6) << j4
         << "  " << setw(6) << j5
         << "  " << setw(6) << j6
         << "  " << setprecision(16) << setw(24) << fx << "\n";

    twoj1 = round ( 2.0 * j1 );
    twoj2 = round ( 2.0 * j2 );
    twoj3 = round ( 2.0 * j3 );
    twoj4 = round ( 2.0 * j4 );
    twoj5 = round ( 2.0 * j5 );
    twoj6 = round ( 2.0 * j6 );

    fx2 = gsl_sf_coupling_6j ( twoj1, twoj2, twoj3, twoj4, twoj5, twoj6 );

    cout << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << setprecision(16) << setw(24) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void gsl_sf_coupling_9j_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_SF_COUPLING_9J_TEST demonstrates GSL_SF_COUPLING_9J.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  double j7;
  double j8;
  double j9;
  int n_data;
  int twoj1;
  int twoj2;
  int twoj3;
  int twoj4;
  int twoj5;
  int twoj6;
  int twoj7;
  int twoj8;
  int twoj9;

  cout << "\n";
  cout << "GSL_SF_COUPLING_9J_TEST:\n";
  cout << "  GSL_SF_COUPLING_9J returns values of\n";
  cout << "  the Wigner 9J coefficient.\n";
  cout << "\n";
  cout << "      J1      J2      J3      J4      J5      J6"
       << "      J7      J8      J9        NINE_J\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    nine_j_values ( &n_data, &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    cout << "  " << setw(6) << j1
         << "  " << setw(6) << j2
         << "  " << setw(6) << j3
         << "  " << setw(6) << j4
         << "  " << setw(6) << j5
         << "  " << setw(6) << j6
         << "  " << setw(6) << j7
         << "  " << setw(6) << j8
         << "  " << setw(6) << j9
         << "  " << setprecision(16) << setw(24) << fx << "\n";

    twoj1 = round ( 2.0 * j1 );
    twoj2 = round ( 2.0 * j2 );
    twoj3 = round ( 2.0 * j3 );
    twoj4 = round ( 2.0 * j4 );
    twoj5 = round ( 2.0 * j5 );
    twoj6 = round ( 2.0 * j6 );
    twoj7 = round ( 2.0 * j7 );
    twoj8 = round ( 2.0 * j8 );
    twoj9 = round ( 2.0 * j9 );

    fx2 = gsl_sf_coupling_9j ( twoj1, twoj2, twoj3, twoj4, twoj5, twoj6,
      twoj7, twoj8, twoj9 );

    cout << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << "      "
         << "  " << setprecision(16) << setw(24) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void gsl_sf_dawson_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GSL_SF_DAWSON_TEST tests the GSL Dawson function routine.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "GSL_SF_DAWSON_TEST:\n";
  cout << "  Evalute Dawson's integral.\n";
  cout << "\n";
  cout << "        X            Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    dawson_values ( &n_data, &x, &fx1 );

    if ( n_data <= 0 )
    {
      break;
    }
    fx2 = gsl_sf_dawson ( x );
    cout                    << "  "
         << setw(12) << x   << "  "
         << setw(12) << fx1 << "  "
         << setw(12) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void nine_j_values ( int *n_data, double *j1, double *j2, double *j3, 
  double *j4, double *j5, double *j6, double *j7, double *j8, double *j9,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    NINE_J_VALUES returns some values of the Wigner 9J function.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *J1, *J2, *J3, *J4, *J5, *J6, *J7, *J8, *J9, 
//    the arguments of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 9

  double fx_vec[N_MAX] = { 
     0.0004270039294528318,
    -0.001228915451058514,
    -0.0001944260688400887,
     0.003338419923885592,
    -0.0007958936865080434,
    -0.004338208690251972,
     0.05379143536399187,
     0.006211299937499411,
     0.03042903097250921 };
  double j1_vec[N_MAX] = {
    1.0,
    1.5,
    2.0,
    1.0, 
    1.5,
    2.0,
    0.5,
    1.0,
    1.5 };
  double j2_vec[N_MAX] = {
    8.0,
    8.0,
    8.0,
    3.0,
    3.0, 
    3.0,
    0.5,
    0.5,
    0.5 };
  double j3_vec[N_MAX] = {
    7.0,
    7.0,
    7.0,
    2.0,
    2.0,
    2.0,
    1.0,
    1.0,
    1.0 };
  double j4_vec[N_MAX] = { 
    6.5,
    6.5,
    6.5,
    4.0,
    4.0,
    4.0,
    2.0,
    2.0,
    2.0 };
  double j5_vec[N_MAX] = {
    7.5,
    7.5,
    7.5,
    1.5,
    1.5,
    1.5,
    1.0,
    1.0,
    1.0 };
  double j6_vec[N_MAX] = { 
    7.5,
    7.5,
    7.5,
    3.0,
    3.0,
    3.0,
    1.5,
    1.5,
    1.5 };
  double j7_vec[N_MAX] = { 
    6.0,
    6.0,
    6.0,
    3.5,
    3.5,
    3.5,
    1.5,
    1.5,
    1.5 };
  double j8_vec[N_MAX] = { 
    10.0,
    10.0,
    10.0,
     2.0,
     2.0,
     2.0,
     0.5,
     0.5,
     0.5 };
  double j9_vec[N_MAX] = { 
    6.0,
    6.0,
    6.0,
    2.0,
    2.0,
    2.0,
    1.5,
    1.5,
    1.5 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *j1 = 0.0;
    *j2 = 0.0;
    *j3 = 0.0;
    *j4 = 0.0;
    *j5 = 0.0;
    *j6 = 0.0;
    *j7 = 0.0;
    *j8 = 0.0;
    *j9 = 0.0;
    *fx = 0.0;
  }
  else
  {
    *j1 = j1_vec[*n_data-1];
    *j2 = j2_vec[*n_data-1];
    *j3 = j3_vec[*n_data-1];
    *j4 = j4_vec[*n_data-1];
    *j5 = j5_vec[*n_data-1];
    *j6 = j6_vec[*n_data-1];
    *j7 = j7_vec[*n_data-1];
    *j8 = j8_vec[*n_data-1];
    *j9 = j9_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int rosenbrock_f ( const gsl_vector *x, void *params, gsl_vector *f )

//****************************************************************************80
//
//  Purpose:
//
//    ROSENBROCK_F evaluates the Rosenbrock function.
//
//  Modified:
//
//    05 August 2005.
//
//  Author:
//
//    John Burkardt
//
{
  double a = ( ( struct rparams *) params )-> a;
  double b = ( ( struct rparams *) params )-> b;

  const double x0 = gsl_vector_get ( x, 0 );
  const double x1 = gsl_vector_get ( x, 1 );
  const double y0 = a * ( 1.0 - x0 );
  const double y1 = b * ( x1 - x0 * x0 );
  gsl_vector_set ( f, 0, y0 );
  gsl_vector_set ( f, 1, y1 );

  return GSL_SUCCESS;
}
//****************************************************************************80

void six_j_values ( int *n_data, double *j1, double *j2, double *j3, 
  double *j4, double *j5, double *j6, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    SIX_J_VALUES returns some values of the Wigner 6J function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      SixJSymbol[{j1,j2,j3},{j4,j5,j6}]
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *J1, *J2, *J3, *J4, *J5, *J6, the arguments 
//    of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 15

  double fx_vec[N_MAX] = { 
     0.03490905138373300,  
    -0.03743025039659792,  
     0.01890866390959560,  
     0.007342448254928643, 
    -0.02358935185081794,  
     0.01913476955215437,  
     0.001288017397724172, 
    -0.01930018366290527,  
     0.01677305949382889,  
     0.005501147274850949, 
    -0.02135439790896831,  
     0.003460364451435387, 
     0.02520950054795585,  
     0.01483990561221713,  
     0.002708577680633186 };
  double j1_vec[N_MAX] = {
    1.0, 
    2.0, 
    3.0, 
    4.0, 
    5.0, 
    6.0, 
    7.0, 
    8.0, 
    9.0, 
   10.0, 
   11.0, 
   12.0, 
   13.0, 
   14.0, 
   15.0 };
  double j2_vec[N_MAX] = {
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0, 
    8.0 };
  double j3_vec[N_MAX] = {
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0, 
    7.0 };
  double j4_vec[N_MAX] = { 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5, 
    6.5 };
  double j5_vec[N_MAX] = {
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5 };
  double j6_vec[N_MAX] = { 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5, 
    7.5 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *j1 = 0.0;
    *j2 = 0.0;
    *j3 = 0.0;
    *j4 = 0.0;
    *j5 = 0.0;
    *j6 = 0.0;
    *fx = 0.0;
  }
  else
  {
    *j1 = j1_vec[*n_data-1];
    *j2 = j2_vec[*n_data-1];
    *j3 = j3_vec[*n_data-1];
    *j4 = j4_vec[*n_data-1];
    *j5 = j5_vec[*n_data-1];
    *j6 = j6_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void three_j_values ( int *n_data, double *j1, double *j2, double *j3, 
  double *m1, double *m2, double *m3, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    THREE_J_VALUES returns some values of the Wigner 3J function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}]
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *J1, *J2, *J3, *M1, *M2, *M3, the arguments 
//    of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 8

  double fx_vec[N_MAX] = { 
     0.2788866755113585, 
    -0.09534625892455923, 
    -0.06741998624632421, 
     0.1533110351679666, 
    -0.1564465546936860, 
     0.1099450412156551, 
    -0.05536235693131719, 
     0.01799835451137786 };
  double j1_vec[N_MAX] = { 
    1.0, 
    2.0, 
    3.0, 
    4.0, 
    5.0, 
    6.0, 
    7.0, 
    8.0 };
  double j2_vec[N_MAX] = { 
    4.5, 
    4.5, 
    4.5, 
    4.5, 
    4.5, 
    4.5, 
    4.5, 
    4.5 };
  double j3_vec[N_MAX] = { 
    3.5, 
    3.5, 
    3.5, 
    3.5, 
    3.5, 
    3.5, 
    3.5, 
    3.5 };
  double m1_vec[N_MAX] = { 
    1.0, 
    1.0, 
    1.0, 
    1.0, 
    1.0, 
    1.0, 
    1.0, 
    1.0 };
  double m2_vec[N_MAX] = { 
    -3.5, 
    -3.5, 
    -3.5, 
    -3.5, 
    -3.5, 
    -3.5, 
    -3.5, 
    -3.5 };
  double m3_vec[N_MAX] = {
    2.5, 
    2.5, 
    2.5, 
    2.5, 
    2.5, 
    2.5, 
    2.5, 
    2.5 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *j1 = 0.0;
    *j2 = 0.0;
    *j3 = 0.0;
    *m1 = 0.0;
    *m2 = 0.0;
    *m3 = 0.0;
    *fx = 0.0;
  }
  else
  {
    *j1 = j1_vec[*n_data-1];
    *j2 = j2_vec[*n_data-1];
    *j3 = j3_vec[*n_data-1];
    *m1 = m1_vec[*n_data-1];
    *m2 = m2_vec[*n_data-1];
    *m3 = m3_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
