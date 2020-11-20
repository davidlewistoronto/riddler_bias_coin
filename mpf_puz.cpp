//---------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#pragma hdrstop


#define my_assert(x, y) if (!(x)) {fprintf (stderr, "assert botched for %s at %d %s\n", y, __LINE__, __FILE__);barf("assertion failure",__LINE__);}

#define max_poly_order	50
#define max_bin_iters	20
#define max_newton_iters	300
#define root_converge_val		1e-14
#define mpf_converge_val		1e-140
#define max_inputs		20

#define n_combin_table 21

int combin_table [n_combin_table] [n_combin_table];

bool debug_solve = false;

#include <boost/multiprecision/gmp.hpp>

using namespace std;
using namespace boost::multiprecision;



int mpf_precision = 512;
#define max_coeff	10

void barf (char *msg, int lc) {
	fprintf (stderr, "in line %d we are displeased by your %s\n", lc, msg);
	exit (1);
}

double mpf_get_root (int n_coeff, double *coeff, double xstart)
{	mpf_t mp_coeff [max_coeff];
	mpf_t fval;
	mpf_t fval_abs;
	mpf_t dfval;
	mpf_t x;
	mpf_t dx;
	mpf_t xpow;
	mpf_t f_converge_tol;
	mpf_t t1;
	mpf_t mp_icoeff;
	mpf_t damp;
	int i;
	int icoeff;
	int fval_cvg_flag;
	int iiter;
	double ret;

	for (i = 0; i < n_coeff; i++)
	{	mpf_init2 (mp_coeff [i], mpf_precision);
		mpf_set_d (mp_coeff [i], coeff [i]);
	}
	mpf_init2 (fval, mpf_precision);
	mpf_init2 (fval_abs, mpf_precision);
	mpf_init2 (dfval, mpf_precision);
	mpf_init2 (x, mpf_precision);
	mpf_init2 (dx, mpf_precision);
	mpf_init2 (xpow, mpf_precision);
	mpf_init2 (f_converge_tol, mpf_precision);
	mpf_init2 (t1, mpf_precision);
	mpf_init2 (mp_icoeff, mpf_precision);
	mpf_init2 (damp, mpf_precision);

	mpf_set_d (x, xstart);
	mpf_set_d (f_converge_tol, mpf_converge_val);
	mpf_set_d (damp, 1.0);

	/* eval fval */
	mpf_set_d (fval, 0.0);
	mpf_set_d (xpow, 1.0);
	for (icoeff = 0; icoeff < n_coeff; icoeff++)
	{	mpf_mul (t1, xpow, mp_coeff [icoeff]);
		mpf_add (fval, t1, fval);
		mpf_mul (xpow, xpow, x);
	}
	mpf_abs (fval_abs, fval);
	fval_cvg_flag = mpf_cmp (fval_abs, f_converge_tol);
	iiter = 0;
	
	while (iiter < max_newton_iters && fval_cvg_flag > 0)
	{	/* eval dfval */
		mpf_set_d (dfval, 0);
		mpf_set_d (xpow, 1.0);
		for (icoeff = 1; icoeff < n_coeff; icoeff++)
		{	mpf_set_d (mp_icoeff, (double) icoeff);
			mpf_mul (t1, xpow, mp_coeff [icoeff]);
			mpf_mul (t1, t1, mp_icoeff);
			mpf_add (dfval, dfval, t1);
			mpf_mul (xpow, xpow, x);
		}

		/* set x = x - f / dfdx */
		mpf_div (dx, fval, dfval);
		mpf_mul (dx, dx, damp);
		mpf_sub (x, x, dx);

#ifdef debug_mp_root
		printf ("f: ");	
		mpf_out_str (stdout, 10, 100, fval);
		printf ("\n");
		printf ("df: ");	
		mpf_out_str (stdout, 10, 100, dfval);
		printf ("\n");
		printf ("dx: ");	
		mpf_out_str (stdout, 10, 100, dx);
		printf ("\n");
		printf ("x: ");	
		mpf_out_str (stdout, 10, 100, x);
		printf ("\n");
#endif

		/* eval fval */
		mpf_set_d (fval, 0.0);
		mpf_set_d (xpow, 1.0);
		for (icoeff = 0; icoeff < n_coeff; icoeff++)
		{	mpf_mul (t1, xpow, mp_coeff [icoeff]);
			mpf_add (fval, t1, fval);
			mpf_mul (xpow, xpow, x);
		}
		mpf_abs (fval_abs, fval);
		fval_cvg_flag = mpf_cmp (fval_abs, f_converge_tol);
		iiter++;

	}
#ifdef print_mpf
	printf ("mpf ");
	mpf_out_str (stdout, 10, 20, x);
	printf ("\n");
#endif

	ret = mpf_get_d (x);

	for (i = 0; i < n_coeff; i++)
	{	mpf_clear (mp_coeff [i]);
	}
	mpf_clear (fval);
	mpf_clear (fval_abs);
	mpf_clear (dfval);
	mpf_clear (x);
	mpf_clear (dx);
	mpf_clear (xpow);
	mpf_clear (f_converge_tol);
	mpf_clear (t1);
	mpf_clear (mp_icoeff);
	mpf_clear (damp);

	my_assert (iiter < max_newton_iters, "newton convergence");
	return ret;
}

#define my_max(x,y) ((x)>(y)?(x):(y))

//---------------------------------------------------------------------------
class my_poly
{	public:
		double coeff [max_poly_order];
		int n_coeff;
		my_poly ();
		void copy_to (my_poly *p);
		double eval_at (double x);
		double eval_integral_at (double x);
		double eval_deriv_at (double x);
		void get_integral (my_poly *p);
		void get_deriv (my_poly *p);
		void scale_x (double s);
		void translate_x (double x, my_poly *p);
		double get_root (double xstart);
		void get_all_roots (double *roots, int &nroots, bool precise);
		void print (void);
		void printlf (char *s);


};




my_poly::my_poly ()
{	int i;

	n_coeff = 1;
	for (i = 0; i < max_poly_order; i++)
	{	coeff [i] = 0.0;
	}
}

void my_poly::copy_to (my_poly *p){
	int i;

	p->n_coeff = n_coeff;
	for (i = 0; i < n_coeff; i++)
	{	p->coeff [i] = coeff [i];
	}
}

void my_poly::scale_x (double s)
{	int i;
	double spow;

	spow = 1.0;
	for (i = 0; i < n_coeff; i++)
	{	coeff [i] = spow * coeff [i];
		spow *= s;
	}
}

void my_poly::translate_x (double dx, my_poly *p)
{	int i;
	int j;
	double dxpow;

	p->n_coeff = n_coeff;
	for (i = 0; i < n_coeff; i++)
	{	p->coeff [i] = 0.0;
	}

	for (i = 0; i < n_coeff; i++)
	{	dxpow = 1.0;
		for (j = i; j >= 0; j--)
		{	p->coeff [j] += coeff [i] * dxpow * combin_table [i] [j];
			dxpow *= dx;
		}
	}
}


double my_poly::eval_at (double x)
{	int i;
	double xpow;
	double r;

	r = 0;
	xpow = 1.0;
	for (i = 0; i < n_coeff; i++)
	{	r += xpow * coeff [i];
		xpow *= x;
	}
	return r;
}

double my_poly::eval_integral_at (double x)
{	int i;
	double xpow;
	double r;

	r = 0;
	xpow = x;
	for (i = 0; i < n_coeff; i++)
	{	r += xpow * coeff [i] / (i + 1);
		xpow *= x;
	}
	return r;
}

void my_poly::get_integral (my_poly *p_int)
{	int i;

	my_assert (n_coeff < max_poly_order, " integrating poly with too many coeff");
	p_int->coeff [0] = 0;
	for (i = 0; i < n_coeff; i++)
	{	p_int->coeff [i + 1] = coeff [i] / (i + 1);
	}
	p_int->n_coeff = n_coeff + 1;
}

void my_poly::get_deriv (my_poly *p_int)
{	int i;

	for (i = 0; i < n_coeff - 1; i++)
	{	p_int->coeff [i] = coeff [i + 1] * (i + 1);
	}
	if (n_coeff == 1)
	{	p_int->n_coeff = 1;
		p_int->coeff [0] = 0.0;
	}
	else
	{	p_int->n_coeff = n_coeff - 1;
	}
}

double my_poly::eval_deriv_at (double x)
{	int i;
	double xpow;
	double r;

	r = 0;
	xpow = 1;
	for (i = 1; i < n_coeff; i++)
	{	r += xpow * coeff [i] * i;
		xpow *= x;
	}
	return r;
}

void my_poly::print (void)
{	int i;
	for (i = 0; i < n_coeff; i++)
	{	printf (" %d:%g", i, coeff [i]);
	}
}

void my_poly::printlf (char * s)
{	printf ("%s", s);
	print ();
	printf ("\n");
}

#define newton_get_root

#ifdef newton_get_root
double my_poly::get_root (double xstart)
{	int iiter;
	double p_val;
	double x;

	iiter = 0;
	x = xstart;
	p_val = eval_at (x);
	while (iiter < max_newton_iters && fabs (p_val) > root_converge_val)
	{	x = x - p_val / eval_deriv_at (x);
		p_val = eval_at (x);
		iiter++;
	}
	my_assert (iiter < max_newton_iters, "newton convergence");
	return x;
}
#else
double my_poly::get_root (double xstart)
{	int iiter;
	double p_val;
	double x;
	double damp;

	damp = 1;
	iiter = max_newton_iters;
	while (damp > .001 && iiter == max_newton_iters)
	{	iiter = 0;
		x = xstart;
		p_val = eval_at (x);
		while (iiter < max_newton_iters && fabs (p_val) > root_converge_val)
		{	x = x - damp * p_val / eval_deriv_at (x);
			p_val = eval_at (x);
			iiter++;
		}
		damp *= 0.9;
	}
	my_assert (iiter < max_newton_iters, "newton convergence");
	return x;
}
#endif


void init_combin_table (void)
{	int i, j;
	combin_table [0] [0] = 1;
	for (i = 1; i < n_combin_table; i++)
	{	combin_table [i] [0] = 1;
		for (j = 1; j <= i; j++)
		{	combin_table [i] [j] = combin_table [i - 1] [j - 1] + combin_table [i - 1] [j];
		}
	}
}

int minus_one_to_n (int n)
{	return (1 - 2 * (n & 1));
}

void my_poly::get_all_roots (double *roots, int &nroots, bool precise)
{   int n_deriv_roots;
	my_poly my_deriv;
	double deriv_roots [max_poly_order];
	int iroot;
	double lower_root_bound;
	double upper_root_bound;
	int icoeff;
	double root_test1;
	double root_test2;
	double root_next;
	double root_eval;
	int iiter;

	lower_root_bound = 0.0;
	upper_root_bound = 0.5 + root_converge_val;
	while (coeff [n_coeff - 1] == 0 && n_coeff > 1)
	{	n_coeff--;
	}
	if (n_coeff == 1)
	{	nroots = 0;
	}
	else if (n_coeff == 2)
	{	roots [0] = coeff [0] / coeff [1];
		if (roots [0] >= lower_root_bound && roots [0] <= upper_root_bound)
		{	nroots = 1;
		}
		else
		{	nroots = 0;
		}
	}
	else
	{  	get_deriv (&my_deriv);
		my_deriv.get_all_roots (deriv_roots, n_deriv_roots, false);

		/* find bounds on roots */
		nroots = 0;
		for (iroot = 0; iroot <= n_deriv_roots; iroot++)
		{	if (iroot == 0)
			{	root_test1 = lower_root_bound;
			}
			else
			{	root_test1 = deriv_roots [iroot - 1];
			}
			if (iroot == n_deriv_roots)
			{	root_test2 = upper_root_bound;
			}
			else
			{	root_test2 = deriv_roots [iroot];
			}

			if (eval_at (root_test1) * eval_at (root_test2) <= 0)
			{
				iiter = 0;
//				while (iiter < max_bin_iters && fabs (root_test1 - root_test2) >= root_converge_val &&
//						fabs (eval_at (root_test1) - eval_at (root_test2)) >= root_converge_val)
				while (iiter < max_bin_iters && fabs (root_test1 - root_test2) >= root_converge_val)
				{	root_next = 0.5 * (root_test1 + root_test2);
					root_eval = eval_at (root_next);
					if (eval_at (root_test1) >= 0 & eval_at (root_test2) <= 0)
					{	if (root_eval >= 0)
						{	root_test1 = root_next;
						}
						else
						{	root_test2 = root_next;
						}
					}
					else
					{   if (root_eval >= 0)
						{	root_test2 = root_next;
						}
						else
						{	root_test1 = root_next;
						}
					}
				}

#ifdef use_secant
				iiter = 0;
				while (iiter < max_newton_iters && fabs (root_test1 - root_test2) >= root_converge_val &&
						fabs (eval_at (root_test1) - eval_at (root_test2)) >= root_converge_val)
				{	root_next = root_test1 - eval_at (root_test1) * (root_test2 - root_test1) / (eval_at (root_test2) - eval_at (root_test1));
					root_test1 = root_test2;
					root_test2 = root_next;
					iiter++;
				}
#endif
				root_next = 0.5 * (root_test1 + root_test2);
                if (precise)
                {
#ifdef use_ieee_root
					root_next = get_root (root_next);
#else
					root_next = mpf_get_root (n_coeff, coeff, root_next);
#endif
                }
				if (root_next >= lower_root_bound && root_next <= upper_root_bound && (nroots == 0 || root_next - roots [nroots - 1] > root_converge_val))
				{	roots [nroots] = root_next;
					nroots++;
				}
			}
		}
	}
}

void solve_set (int n_one_set [max_inputs], int n_inputs)
{	int i_one;
	my_poly s_poly;
	int i_zero;
	int n_one;
	int i_term;
	int ip;
	double p;
	double psolve;
	double proots [max_poly_order];
	int nroots;
	int iroot;

	if (debug_solve)
	{	printf ("solve");
		for (i_one = 0; i_one <= n_inputs; i_one++)
		{	printf (" %d", n_one_set [i_one]);
		}
	}

	s_poly.n_coeff = n_inputs + 1;
	for (i_one = 0; i_one <= n_inputs; i_one++)
	{	n_one = n_one_set [i_one];
		i_zero = n_inputs - i_one;
		/* poly is n_one * (p ^ i_one * (1 - p) ^ i_zero) */

		/* coefficients of (1 - p) ^ i_zero are combin_table [i_zero] [i] * (-1 ^ (i) */

		for (i_term = 0; i_term <= i_zero; i_term++)
		{	s_poly.coeff [i_one + i_term] += n_one * combin_table [i_zero] [i_term] * minus_one_to_n (i_term);
		}
	}
	if (debug_solve)
	{	s_poly.printlf (" poly");
	}

	/* the set of functions has symmetry about both x = .5 and y = .5 so we can consider only
	 * functions with f(0) = 0, and search for a root in [0..0.5].
	 * That is, if there is a function f(x), there also exists f(1-x), 1-f(x), and 1-f(1-x).
	 */

	if (s_poly.eval_at (0.0) <= 1e-9)
	{	/* solve for poly = 0.5 */
	
		s_poly.coeff [0] -= 0.5;
		/* only search for a root in 0..0.5 */
#ifdef solve_incr
		for (ip = 0; ip < 500; ip++)
		{	p = ip * .001;
			if (s_poly.eval_at (p) * s_poly.eval_at (p + .001) <= 0)
			{
//			psolve = s_poly.get_root (p);
			psolve = mpf_get_root (s_poly.n_coeff, s_poly.coeff, p);
				printf ("root %1.17g\n", psolve);
			}
		}
#else
		s_poly.get_all_roots (proots, nroots, true);
		for (iroot = 0; iroot < nroots; iroot++)
		{	printf ("mr %1.17g\n", proots [iroot]);
		}
#endif
	}

}

void solve_puz_recur (int n_one_set [max_inputs], int ilevel, int n_inputs)
{	int i_one;

	if (ilevel == n_inputs + 1)
	{	solve_set (n_one_set, n_inputs);
	}
	else
	{	for (i_one = 0; i_one <= combin_table [n_inputs] [ilevel]; i_one++)
		{	n_one_set [ilevel] = i_one;
			solve_puz_recur (n_one_set, ilevel + 1, n_inputs);
		}
	}
}

void solve_puz (int n_inputs)
{	int n_one_set [max_inputs];

	solve_puz_recur (n_one_set, 0, n_inputs);
}

#pragma argsused
int main(int argc, char* argv[])
{   int n_inputs;

	if (argc != 2 || sscanf (argv [1], "%d", &n_inputs) != 1)
	{	fprintf (stderr, "bad args\n");
		exit (1);
	}

	init_combin_table ();

//	test_poly ();

	solve_puz (n_inputs);

	return 0;
}
