#include <stdio.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

int
func (double t,const double y[], double f[],
      void *params)
{
  (void)(t); 
  double mu = *(double *)params;
  double k,p;
  k = y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
  p= k*sqrt(k);
       //y[0]=x , y[1]=y , y[2]=z  
  f[0] = y[3];  //Vx
  f[1] = y[4];  //Vy
  f[2] = y[5];  //Vz
  f[3] = -mu*y[0]/p;   // -mu*x/(x^2+y^2+z^2)^3/2
  f[4] = -mu*y[1]/p;   // -mu*y/(x^2+y^2+z^2)^3/2
  f[5] = -mu*y[2]/p;   // -mu*z/(x^2+y^2+z^2)^3/2
  return GSL_SUCCESS;
}

int
jac (double t,const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t); 
  double mu = *(double *)params;
    double k;
  k = y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
  double r=k*sqrt(k);
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 6, 6);

  // fill the Jacobian matrix 
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 1.0);
  gsl_matrix_set (m, 0, 4, 0.0);
  gsl_matrix_set (m, 0, 5, 0.0);
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, 0.0);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 1, 4, 1.0);
  gsl_matrix_set (m, 1, 5, 0.0);
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 3, 0.0);
  gsl_matrix_set (m, 2, 4, 0.0);
  gsl_matrix_set (m, 2, 5, 1.0);
  gsl_matrix_set (m, 3, 0, mu*(2*y[0]*y[0]-(y[1]*y[1]+y[2]*y[2]))/(r*k));
  gsl_matrix_set (m, 3, 1, 3*mu*y[0]*y[1]/(r*k));
  gsl_matrix_set (m, 3, 2, 3*mu*y[0]*y[2]/(r*k));
  gsl_matrix_set (m, 3, 3, 0.0);
  gsl_matrix_set (m, 3, 4, 0.0);
  gsl_matrix_set (m, 3, 5, 0.0);
  gsl_matrix_set (m, 4, 0, 3*mu*y[1]*y[0]/(r*k));
  gsl_matrix_set (m, 4, 1, mu*(2*y[1]*y[1]-(y[0]*y[0]+y[2]*y[2]))/(r*k));
  gsl_matrix_set (m, 4, 2, 3*mu*y[1]*y[2]/(r*k));
  gsl_matrix_set (m, 4, 3, 0.0);
  gsl_matrix_set (m, 4, 4, 0.0);
  gsl_matrix_set (m, 4, 5, 0.0);
  gsl_matrix_set (m, 5, 0, 3*mu*y[2]*y[0]/(r*k));
  gsl_matrix_set (m, 5, 1, 3*mu*y[2]*y[1]/(r*k));
  gsl_matrix_set (m, 5, 2, mu*(2*y[2]*y[2]-(y[0]*y[0]+y[1]*y[1]))/(r*k));
  gsl_matrix_set (m, 5, 3, 0.0);
  gsl_matrix_set (m, 5, 4, 0.0);
  gsl_matrix_set (m, 5, 5, 0.0);

// set explicit t dependence of f[i] 
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[3] = 0.0;
  dfdt[4] = 0.0;
  dfdt[5] = 0.0;  
  return GSL_SUCCESS;
}

int
main (void)
{
    double k3=0.01720209895, mu ;
    mu= k3*k3;
  gsl_odeiv2_system sys = {func, jac, 6, &mu};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
                                  1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 10000;
  double y[6] = { -6.685075643890280E-01, 1.706558779772845E+00, -1.123043436831628E+00, -1.207950278464835E-02,-3.384553931794463E-03,3.350958907738935E-03};
   //double y[6]= {1, 0, 0, 0, 0.1, 0};
  for (i = 0; i <= t1; i++)
    {
      double ti = i;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      printf ("%.5le %.5le %.5le %.5le %.5le %.5le %.5le\n", t, y[0], y[1],  y[2], y[3],  y[4], y[5]);
    }

  gsl_odeiv2_driver_free (d);
  return 0;
}



