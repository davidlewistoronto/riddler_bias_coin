//---------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#pragma hdrstop


#define my_assert(x, y) if (!(x)) {fprintf (stderr, "assert botched for %s at %d %s\n", y, __LINE__, __FILE__);barf("assertion failure",__LINE__);}

#define max_poly_order	50
#define max_root_iters	50
#define root_converge_val		1e-14
#define max_inputs		20

#define n_combin_table 21

int combin_table [n_combin_table] [n_combin_table];
bool debug_solve = false;

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
		void scale_x (double s);
		void translate_x (double x, my_poly *p);
		double get_root (double xstart);
		void print (void);
		void printlf (char *s);


};



void barf (char *msg, int lc) {
	fprintf (stderr, "in line %d we are displeased by your %s\n", lc, msg);
	exit (1);
}

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
	while (iiter < max_root_iters && fabs (p_val) > root_converge_val)
	{	x = x - p_val / eval_deriv_at (x);
		p_val = eval_at (x);
		iiter++;
	}
	my_assert (iiter < max_root_iters, "newton convergence");
	return x;
}
#else
double my_poly::get_root (double xstart)
{	int iiter;
	double p_val;
	double x;
	double damp;

	damp = 1;
	iiter = max_root_iters;
	while (damp > .001 && iiter == max_root_iters)
	{	iiter = 0;
		x = xstart;
		p_val = eval_at (x);
		while (iiter < max_root_iters && fabs (p_val) > root_converge_val)
		{	x = x - damp * p_val / eval_deriv_at (x);
			p_val = eval_at (x);
			iiter++;
		}
		damp *= 0.9;
	}
	my_assert (iiter < max_root_iters, "newton convergence");
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

void solve_set (int n_one_set [max_inputs], int n_inputs)
{	int i_one;
	my_poly s_poly;
    int i_zero;
    int n_one;
    int i_term;
    int ip;
    double p;
    double psolve;

	printf ("solve");
    for (i_one = 0; i_one <= n_inputs; i_one++)
    {	printf (" %d", n_one_set [i_one]);
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
    s_poly.printlf (" poly");

    /* solve for poly = 0.5 */

    s_poly.coeff [0] -= 0.5;
    for (ip = 0; ip < 1000; ip++)
    {	p = ip * .001;
    	if (s_poly.eval_at (p) * s_poly.eval_at (p + .001) <= 0)
        {	psolve = s_poly.get_root (p);
        	printf ("root %g\n", psolve);
        }
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
//---------------------------------------------------------------------------

