//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stddef.h>
#include <ctime>
#include <climits>
#include <cstdio>
#include <float.h>

#if defined(WIN32) || defined(__WIN32__)
#	ifndef WIN32
#		define WIN32
#	endif
#	include <windows.h> /* for GetTickCount */
#	include <process.h> /* for getpid */
#	ifdef _MSC_VER  // for Visual C++
#		define getpid _getpid
#	endif
#else  // not Windows
#	include <sys/time.h> /* for getpid */
#	include <unistd.h>
#endif


#include "Rmathdefines.h"
#include "maths.h"
#include "mtrand/mtrand.h"

#if defined(_MSC_VER) || defined(__BORLANDC__) // Visual C++ or Borland
#	include <algorithm>
	using std::max;
	using std::min;
#	define fmax2 max
#	define fmin2 min
#else					// GCC
#	define fmax2 fmax
#	define fmin2 fmin
#endif



long long int get_timediff(bool start) {
#ifdef WIN32
	static  DWORD  prevtime;
	if (start) {
		prevtime = GetTickCount();
		return 0;
	} else {
		return (long long int) (GetTickCount() - prevtime);
	}
#else
	static long long int prevtime;
	timeval  tv;
	gettimeofday(&tv, NULL);
	if (start) {
		prevtime = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
		return 0;
	} else {
		long long int curtime = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
		return curtime - prevtime;
	}
#endif
}


long long int get_time_msec(void) {
#ifdef WIN32
	return GetTickCount();
#else
	timeval time;
	gettimeofday(&time, NULL);
	return (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif
}

unsigned long int random_seed(void) {
#ifdef DEBUG_RAND
	return 666;
#else
#ifdef WIN32
	return (GetTickCount() % LONG_MAX) + getpid();
#else
	timeval time;
	gettimeofday(&time, NULL);
	return (time.tv_usec % LONG_MAX)  + getpid();
#endif
#endif
}

// Source for MTRand:
// http://bedaux.net/mtrand/
//

double unif_rand(void) {
	static MTRand drand(random_seed());
	double ret = drand();
	return (ret < 1.0) ? ret : myfmod(ret, 1.0);
}


//unsigned long int random_seed(void) {
//	return (GetTimeUs64() % LONG_MAX) + getpid();
//}

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */
/*
int64_t GetTimeUs64(void) {
#ifdef WIN32
	 Windows
	FILETIME ft;
	LARGE_INTEGER li;
	 Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	 * to a LARGE_INTEGER structure.
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;
	int64_t ret = li.QuadPart;
	ret -= 116444736000000000LL;  Convert from file time to UNIX epoch time.
	ret /= 10;  From 100 nano seconds (10^-7) to 1 millisecond (10^-6) intervals
#else  Linux
	struct timeval tv;
	gettimeofday(&tv, NULL);
	int64_t ret = tv.tv_usec;
	 Adds the seconds (10^0) after converting them to microseconds (10^-6)
	ret += (tv.tv_sec * 1e+06);
#endif
	return ret;
}
*/


//double unif_rand(void) {
//	static double v[256];
//	static bool set = false;
//	static size_t idx = 0;
//	if(!set) {
//		for (int i = 0; i < 256; ++i) v[i] = __unif_rand();
//		set = true;
//	}
//	if(++idx > 255) idx = 0;
//	return v[idx];
//}


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<class T>
int compare_numeric(const void * a, const void * b) {
//	return (* (T*)(a) - *(T*)(b)); // c-style casting
	return sgn(*static_cast<const T*>(a) - *static_cast<const T*>(b));
}

#include <stdio.h>

void mean_sd(double* x, const size_t n, double * result) {
	double mean, sum = 0, sum2 = 0, dn;
	dn = (double) n;
	for (size_t i = 0; i < n; ++i) {
		sum += x[i];
		sum2 += pow2(x[i]);
	}
	mean = sum / dn;
	result[0] = mean;
	result[1] = sqrt(((sum2 / dn) - pow2(mean)) * dn / (dn - 1.0));
}

void quantile(double* x, const size_t n, const double * probs, const size_t np,
		double * qs) {
	double* index = new double[np];
	int* lo = new int[np];
	int* hi = new int[np];

	for (size_t i = 0; i < np; ++i) {
		index[i] = (n - 1) * probs[i];
		lo[i] = (int) floor(index[i]);
		hi[i] = (int) ceil(index[i]);
	}
	qsort(x, n, sizeof(double), compare_numeric<double>);
	for (size_t i = 0; i < np; ++i) {
		qs[i] = x[lo[i]];
	}
	double h;
	for (size_t i = 0; i < np; ++i) {
		if (index[i] > lo[i]) {
			h = index[i] - lo[i];
			qs[i] = (1 - h) * qs[i] + h * x[hi[i]];
		}
	}
	delete[] hi;
	delete[] lo;
	delete[] index;
}

double myfmod(double x1, double x2) {
	return x1 - floor(x1 / x2) * x2;
}

double myfmod(long int x1, long int x2) {
	return x1 - floor((double) x1 / (double) x2) * x2;
}

double myfmod(size_t x1, size_t x2) {
	return x1 - floor(double(x1) / double(x2)) * x2;
}

void rand_perm(int k, int *a) {
	int x = 0;
	while (x < k) {
		int n = (int) runif(0, k);
		int temp = a[x];
		a[x++] = a[n];
		a[n] = temp;
	}
}

void sincos(const double a, double & sinv, double & cosv) {
	double z(tan(a / 2)), zsq(z * z);
	sinv = 2 * z / (1 + zsq);
	cosv = (1 - zsq) / (1 + zsq);
}

double angle(double x, double y) {
	double ret;
	if ((x == 0) && (y == 0)) {
		ret = unif_rand() * M2_PI;
	} else if (x != 0) {
		ret = myatan(y / x);
		if (x < 0) {
			ret += M_PI;
		} else if (y < 0) {
			ret += M2_PI;
		}
	} else {
		if (y > 0) {
			ret = M_PI_2;
		} else {
			ret = M_PI + M_PI_2;
		}
	}

	return (ret);
}

void compute_summary(double* values, size_t size, numeric_summary output) {
	static int counter = 0;
	if (size != 0) {
		++counter;
		//printf("compute_summary (%i) of N = %i: \n", counter, (int) size);

		qsort(values, size, sizeof(double), compare_numeric<double>);
		double probs[NUMSUM_QTL_SIZE] = NUMSUM_QTL_PROBS;
		quantile(values, size, probs, NUMSUM_QTL_SIZE, output);
//		for (int i = 0; i < 4; ++i) {
//			if(output[i] > output[i + 1]) {
//				fprintf(stderr, "Error: in 'compute_summary', q[%g] = %g > q[%g] = %g \n", probs[i], output[i], probs[i + 1], output[i + 1]);
//				for (size_t i = 0; i < size; ++i)
//					printf("-- value[%d] %g \n", (int) i, values[i]);
//
//				for (int i = 0; i < 5; ++i) {
//					printf("-- q[%g] = %g \n", probs[i], output[i]);
//				}
//				exit(111);
//			}
//		}

		mean_sd(values, size, output + NUMSUM_QTL_SIZE);

		//printf("mean: %g, sd: %g \n", output[5], output[6]);

	} else {
		for (int i = 0; i < NUMSUM_QTL_SIZE + 2; ++i) output[i] = 0.0;
	}

	/*min, q1, median, q2, max, mean, sd*/
}

int sample1_with_prob(double *prob, const size_t nprob) {
	//double cum_prob[nprob], sum_prob = 0;
	//double* cum_prob = NULL;   // Pointer to int, initialize to nothing.
	double* cum_prob = new double[nprob];  // Allocate n ints and save ptr in a.
	double sum_prob = 0;
	for (size_t i = 0; i < nprob; i++) {
		sum_prob += prob[i];
		cum_prob[i] = sum_prob;
	}
	double x = unif_rand() * sum_prob;
	size_t j;
	for (j = 0; j < nprob; j++)
		if (x < cum_prob[j])
			break;

	delete[] cum_prob;  // When done, free memory pointed to by a.
	//cum_prob = NULL;     // Clear a to prevent using invalid memory reference.
	return (int) j;
}


void sample_with_prob(size_t * result, size_t n, double *prob, const size_t nprob) {
	double* cum_prob = new double[nprob];  // Allocate n ints and save ptr in a.
	double sum_prob = 0;
	for (size_t i = 0; i < nprob; i++) {
		sum_prob += prob[i];
		cum_prob[i] = sum_prob;
	}

	for (size_t i = 0; i < n; ++i) {
		double x = unif_rand() * sum_prob;
		size_t k;
		for (k = 0; k < nprob; k++)
			if (x < cum_prob[k])
				break;
		result[i] = k;
	}

	delete[] cum_prob;  // When done, free memory pointed to by a.
	//cum_prob = NULL;     // Clear a to prevent using invalid memory reference.
}


size_t sample1(size_t n) {
//	return rand() % n;
	return (size_t) floor(unif_rand() * n);
}

double fsign(double x, double y) {
	return ((y >= 0) ? fabs(x) : -fabs(x));
}

int imax2(int x, int y) {
	return (x < y) ? y : x;
}

int imin2(int x, int y) {
	return (x < y) ? x : y;
}

//#if defined(_MSC_VER) || defined(__BORLANDC__)
//// http://www.johndcook.com/cpp_expm1.html
//double expm1(double x) {
//	if (fabs(x) < 1e-5)
//	return x + 0.5 * x * x;
//	else
//	return exp(x) - 1.0;
//}
//#endif

// END R stuff

double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p) {
	double p_, q, r, val;

	R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

	if (sigma < 0)
		return ML_NAN;
	if (sigma == 0)
		return mu;

	p_ = R_DT_qIv(p);/* real lower_tail prob. p */
	q = p_ - 0.5;

	if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
		r = .180625 - q * q;
		val =
				q
						* (((((((r * 2509.0809287301226727
								+ 33430.575583588128105) * r
								+ 67265.770927008700853) * r
								+ 45921.953931549871457) * r
								+ 13731.693765509461125) * r
								+ 1971.5909503065514427) * r
								+ 133.14166789178437745) * r
								+ 3.387132872796366608)
						/ (((((((r * 5226.495278852854561
								+ 28729.085735721942674) * r
								+ 39307.89580009271061) * r
								+ 21213.794301586595867) * r
								+ 5394.1960214247511077) * r
								+ 687.1870074920579083) * r
								+ 42.313330701600911252) * r + 1.);
	} else { /* closer than 0.075 from {0,1} boundary */

		/* r = min(p, 1-p) < 0.075 */
		if (q > 0)
			r = R_DT_CIv(p);/* 1-p */
		else
			r = p_;/* = R_DT_Iv(p) ^=  p */

		r =
				sqrt(
						-((log_p
								&& ((lower_tail && q <= 0)
										|| (!lower_tail && q > 0))) ?
								p : /* else */log(r)));
		/* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

		if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
			r += -1.6;
			val =
					(((((((r * 7.7454501427834140764e-4
							+ .0227238449892691845833) * r
							+ .24178072517745061177) * r
							+ 1.27045825245236838258) * r
							+ 3.64784832476320460504) * r
							+ 5.7694972214606914055) * r + 4.6303378461565452959)
							* r + 1.42343711074968357734)
							/ (((((((r * 1.05075007164441684324e-9
									+ 5.475938084995344946e-4) * r
									+ .0151986665636164571966) * r
									+ .14810397642748007459) * r
									+ .68976733498510000455) * r
									+ 1.6763848301838038494) * r
									+ 2.05319162663775882187) * r + 1.);
		} else { /* very close to  0 or 1 */
			r += -5.;
			val =
					(((((((r * 2.01033439929228813265e-7
							+ 2.71155556874348757815e-5) * r
							+ .0012426609473880784386) * r
							+ .026532189526576123093) * r
							+ .29656057182850489123) * r + 1.7848265399172913358)
							* r + 5.4637849111641143699) * r
							+ 6.6579046435011037772)
							/ (((((((r * 2.04426310338993978564e-15
									+ 1.4215117583164458887e-7) * r
									+ 1.8463183175100546818e-5) * r
									+ 7.868691311456132591e-4) * r
									+ .0148753612908506148525) * r
									+ .13692988092273580531) * r
									+ .59983220655588793769) * r + 1.);
		}

		if (q < 0.0)
			val = -val;
		/* return (q >= 0.)? r : -r ;*/
	}
	return mu + sigma * val;
}

double exp_rand(void) {
	/* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
	/* The highest n (here 8) is determined by q[n-1] = 1.0 */
	/* within standard precision */
	const static double q[] = { 0.6931471805599453, 0.9333736875190459,
			0.9888777961838675, 0.9984959252914960, 0.9998292811061389,
			0.9999833164100727, 0.9999985691438767, 0.9999998906925558,
			0.9999999924734159, 0.9999999995283275, 0.9999999999728814,
			0.9999999999985598, 0.9999999999999289, 0.9999999999999968,
			0.9999999999999999, 1.0000000000000000 };
	double a, u, ustar, umin;
	int i;

	a = 0.;
	/* precaution if u = 0 is ever returned */
	u = unif_rand();
	while (u <= 0.0 || u >= 1.0)
		u = unif_rand();
	for (;;) {
		u += u;
		if (u > 1.0)
			break;
		a += q[0];
	}
	u -= 1.;

	if (u <= q[0])
		return a + u;

	i = 0;
	ustar = unif_rand();
	umin = ustar;
	do {
		ustar = unif_rand();
		if (ustar < umin)
			umin = ustar;
		i++;
	} while (u > q[i]);
	return a + umin * q[0];
}

double norm_rand(void) {
	double u1;
//    int i;
#define BIG 134217728 /* 2^27 */
	/* unif_rand() alone is not of high enough precision */
	u1 = unif_rand();
	u1 = (int) (BIG * u1) + unif_rand();
	return qnorm5(u1 / BIG, 0.0, 1.0, 1, 0);
}

//------------------------------------------------------------------------------

double runif(double a, double b) {
	return a + (unif_rand() * (b - a));
}

double rpois(double mu) {
	/* Factorial Table (0:9)! */
	const static double fact[10] = { 1., 1., 2., 6., 24., 120., 720., 5040.,
			40320., 362880. };

	/* These are static --- persistent between calls for same mu : */
	static int l, m;

	static double b1, b2, c, c0, c1, c2, c3;
	static double pp[36], p0, p, q, s, d, omega;
	static double big_l;/* integer "w/o overflow" */
	static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/

	/* Local Vars  [initialize some for -Wall]: */
	double del, difmuk = 0., E = 0., fk = 0., fx, fy, g, px, py, t, u = 0., v,
			x;
	double pois = -1.;
	int k, kflag, big_mu, new_big_mu = false;

	if (!R_FINITE(mu) || mu < 0)
		return ML_NAN;

	if (mu <= 0.)
		return 0.;

	big_mu = mu >= 10.;
	if (big_mu)
		new_big_mu = false;

	if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

		if (big_mu) {
			new_big_mu = true;
			/* Case A. (recalculation of s,d,l	because mu has changed):
			 * The poisson probabilities pk exceed the discrete normal
			 * probabilities fk whenever k >= m(mu).
			 */
			muprev = mu;
			s = sqrt(mu);
			d = 6. * mu * mu;
			big_l = floor(mu - 1.1484);
			/* = an upper bound to m(mu) for all mu >= 10.*/
		} else { /* Small mu ( < 10) -- not using normal approx. */

			/* Case B. (start new table and calculate p0 if necessary) */

			/*muprev = 0.;-* such that next time, mu != muprev ..*/
			if (mu != muprev) {
				muprev = mu;
				m = imax2(1, (int) mu);
				l = 0; /* pp[] is already ok up to pp[l] */
				q = p0 = p = exp(-mu);
			}

			REPEAT {
				/* Step U. uniform sample for inversion method */
				u = unif_rand();
				if (u <= p0)
					return 0.;

				/* Step T. table comparison until the end pp[l] of the
				 pp-table of cumulative poisson probabilities
				 (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
				if (l != 0) {
					for (k = (u <= 0.458) ? 1 : imin2(l, m); k <= l; k++)
						if (u <= pp[k])
							return (double) k;
					if (l == 35) /* u > pp[35] */
						continue;
				}
				/* Step C. creation of new poisson
				 probabilities p[l..] and their cumulatives q =: pp[k] */
				l++;
				for (k = l; k <= 35; k++) {
					p *= mu / k;
					q += p;
					pp[k] = q;
					if (u <= q) {
						l = k;
						return (double) k;
					}
				}
				l = 35;
			} /* end(REPEAT) */
		}/* mu < 10 */

	} /* end {initialize persistent vars} */

	/* Only if mu >= 10 : ----------------------- */

	/* Step N. normal sample */
	g = mu + s * norm_rand();/* norm_rand() ~ N(0,1), standard normal */

	if (g >= 0.) {
		pois = floor(g);
		/* Step I. immediate acceptance if pois is large enough */
		if (pois >= big_l)
			return pois;
		/* Step S. squeeze acceptance */
		fk = pois;
		difmuk = mu - fk;
		u = unif_rand(); /* ~ U(0,1) - sample */
		if (d * u >= difmuk * difmuk * difmuk)
			return pois;
	}

	/* Step P. preparations for steps Q and H.
	 (recalculations of parameters if necessary) */

	if (new_big_mu || mu != muprev2) {
		/* Careful! muprev2 is not always == muprev
		 because one might have exited in step I or S
		 */
		muprev2 = mu;
		omega = M_1_SQRT_2PI / s;
		/* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
		 * approximations to the discrete normal probabilities fk. */

		b1 = one_24 / mu;
		b2 = 0.3 * b1 * b1;
		c3 = one_7 * b1 * b2;
		c2 = b2 - 15. * c3;
		c1 = b1 - 6. * b2 + 45. * c3;
		c0 = 1. - b1 + 3. * b2 - 15. * c3;
		c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
	}

	if (g >= 0.) {
		/* 'Subroutine' F is called (kflag=0 for correct return) */
		kflag = 0;
		goto Step_F;
	}

	REPEAT {
		/* Step E. Exponential Sample */

		E = exp_rand(); /* ~ Exp(1) (standard exponential) */

		/*  sample t from the laplace 'hat'
		 (if t <= -0.6744 then pk < fk for all mu >= 10.) */
		u = 2 * unif_rand() - 1.;
		t = 1.8 + fsign(E, u);
		if (t > -0.6744) {
			pois = floor(mu + s * t);
			fk = pois;
			difmuk = mu - fk;

			/* 'subroutine' F is called (kflag=1 for correct return) */
			kflag = 1;

			Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

			if (pois < 10) { /* use factorials from table fact[] */
				px = -mu;
				py = pow(mu, pois) / fact[(int) pois];
			} else {
				/* Case pois >= 10 uses polynomial approximation
				 A_0-A_7 for accuracy when advisable */
				del = one_12 / fk;
				del = del * (1. - 4.8 * del * del);
				v = difmuk / fk;
				if (fabs(v) <= 0.25) {
					px = fk * v * v
							* (((((((A_7 * v + A_6) * v + A_5) * v + A_4) * v
									+ A_3) * v + A_2) * v + A_1) * v + A_0)
							- del;
				} else {/* |v| > 1/4 */
					px = fk * log(1. + v) - difmuk - del;
				}
				py = M_1_SQRT_2PI / sqrt(fk);
			}
			x = (0.5 - difmuk) / s;
			x *= x;/* x^2 */
			fx = -0.5 * x;
			fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
			if (kflag > 0) {
				/* Step H. Hat acceptance (E is repeated on rejection) */
				if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
					break;
			} else {
				/* Step Q. Quotient acceptance (rare case) */
				if (fy - u * fy <= py * exp(px - fx))
					break;
			}
		} /* t > -.67.. */
	}
	return pois;
}

double rgamma(double a, double scale) {
	/* Constants : */
	const static double sqrt32 = 5.656854;
	const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

	/* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
	 * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
	 * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
	 */
	const static double q1 = 0.04166669;
	const static double q2 = 0.02083148;
	const static double q3 = 0.00801191;
	const static double q4 = 0.00144121;
	const static double q5 = -7.388e-5;
	const static double q6 = 2.4511e-4;
	const static double q7 = 2.424e-4;

	/* State variables [FIXME for threading!] :*/
	static double aa = 0.;
	static double aaa = 0.;
	static double s, s2, d; /* no. 1 (step 1) */
	static double q0, b, si, c;/* no. 2 (step 4) */

	double e, p, q, r, t, u, v, w, x, ret_val;

	if (!R_FINITE(a) || !R_FINITE(scale) || a <= 0.0 || scale <= 0.0)
		return ML_NAN;

	if (a < 1.) { /* GS algorithm for parameters a < 1 */
		e = 1.0 + exp_m1 * a;
		REPEAT {
			p = e * unif_rand();
			if (p >= 1.0) {
				x = -log((e - p) / a);
				if (exp_rand() >= (1.0 - a) * log(x))
					break;
			} else {
				x = exp(log(p) / a);
				if (exp_rand() >= x)
					break;
			}
		}
		return scale * x;
	}

	/* --- a >= 1 : GD algorithm --- */

	/* Step 1: Recalculations of s2, s, d if a has changed */
	if (a != aa) {
		aa = a;
		s2 = a - 0.5;
		s = sqrt(s2);
		d = sqrt32 - s * 12.0;
	}
	/* Step 2: t = standard normal deviate,
	 x = (s,1/2) -normal deviate. */

	/* immediate acceptance (i) */
	t = norm_rand();
	x = s + 0.5 * t;
	ret_val = x * x;
	if (t >= 0.0)
		return scale * ret_val;

	/* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
	u = unif_rand();
	if (d * u <= t * t * t)
		return scale * ret_val;

	/* Step 4: recalculations of q0, b, si, c if necessary */

	if (a != aaa) {
		aaa = a;
		r = 1.0 / a;
		q0 =
				((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r + q2) * r
						+ q1) * r;

		/* Approximation depending on size of parameter a */
		/* The constants in the expressions for b, si and c */
		/* were established by numerical experiments */

		if (a <= 3.686) {
			b = 0.463 + s + 0.178 * s2;
			si = 1.235;
			c = 0.195 / s - 0.079 + 0.16 * s;
		} else if (a <= 13.022) {
			b = 1.654 + 0.0076 * s2;
			si = 1.68 / s + 0.275;
			c = 0.062 / s + 0.024;
		} else {
			b = 1.77;
			si = 0.75;
			c = 0.1515 / s;
		}
	}
	/* Step 5: no quotient test if x not positive */

	if (x > 0.0) {
		/* Step 6: calculation of v and quotient q */
		v = t / (s + s);
		if (fabs(v) <= 0.25)
			q = q0
					+ 0.5 * t * t
							* ((((((A_7 * v + A_6) * v + A_5) * v + A_4) * v
									+ A_3) * v + A_2) * v + A_1) * v;
		else
			q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);

		/* Step 7: quotient acceptance (q) */
		if (log(1.0 - u) <= q)
			return scale * ret_val;
	}

	REPEAT {
		/* Step 8: e = standard exponential deviate
		 *	u =  0,1 -uniform deviate
		 *	t = (b,si)-double exponential (laplace) sample */
		e = exp_rand();
		u = unif_rand();
		u = u + u - 1.0;
		if (u < 0.0)
			t = b - si * e;
		else
			t = b + si * e;
		/* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
		if (t >= -0.71874483771719) {
			/* Step 10:	 calculation of v and quotient q */
			v = t / (s + s);
			if (fabs(v) <= 0.25)
				q = q0
						+ 0.5 * t * t
								* ((((((A_7 * v + A_6) * v + A_5) * v + A_4) * v
										+ A_3) * v + A_2) * v + A_1) * v;
			else
				q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
			/* Step 11:	 hat acceptance (h) */
			/* (if q not positive go to step 8) */
			if (q > 0.0) {
				w = expm1(q);
				/*  ^^^^^ original code had approximation with rel.err < 2e-7 */
				/* if t is rejected sample again at step 8 */
				if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
					break;
			}
		}
	} /* REPEAT .. until  `t' is accepted */
	x = s + 0.5 * t;
	return scale * x * x;
}

double rnorm(double mu, double sigma) {
	if (!R_FINITE(mu) || !R_FINITE(sigma) || sigma < 0.0)
		return ML_NAN;

	if (sigma == 0.0)
		return mu;
	else
		return mu + sigma * norm_rand();
}

double rexp(double scale) { // scale == 1 / rate
	if (!R_FINITE(scale) || scale <= 0.0)
		return ML_NAN;
	return scale * exp_rand();
}

double rchisq(double df) {
	if (!R_FINITE(df) || df <= 0.0)
		return ML_NAN;
	return rgamma(df / 2.0, 2.0);
}

//------------------------------------------------------------------------------

double rinvgauss(double mu, double lambda) {
	if (mu <= 0)
		return ML_NAN; // stop("mu must be positive");
	if (lambda <= 0)
		return ML_NAN; // stop("lambda must be positive");
	double y2, u, r2, r1;

	y2 = rchisq(1);
	u = unif_rand();
	r2 = mu / (2 * lambda)
			* (2 * lambda + mu * y2
					+ sqrt(4 * lambda * mu * y2 + (mu * mu * y2 * y2)));
	r1 = (mu * mu) / r2;

	return u < (mu / (mu + r1)) ? r1 : r2;
}

/*
 *  DESCRIPTION
 *    Random variates from the negative binomial distribution.
 *  NOTES
 *    x = the number of failures before the n-th success
 *  METHOD
 *    Generate lambda as gamma with shape parameter n and scale
 *    parameter p/(1-p).  Return a Poisson deviate with mean lambda.
 */

double rnbinom(double size, double prob) {
	if (!R_FINITE(size) || !R_FINITE(prob) || size <= 0 || prob <= 0
			|| prob > 1.)
		/* prob = 1 is ok, PR#1218 */
		return ML_NAN;
	return (prob == 1.) ? 0. : rpois(rgamma(size, (1. - prob) / prob));
}

double rnbinom_mu(double size, double mu) {
	if (!R_FINITE(size) || !R_FINITE(mu) || size <= 0 || mu < 0)
		return ML_NAN;
	return (mu == 0) ? 0 : rpois(rgamma(size, mu / size));
}

double rcauchy(double location, double scale) {
	if (ISNAN(location) || !R_FINITE(scale) || scale < 0)
		ML_ERR_return_NAN;
	if (scale == 0. || !R_FINITE(location))
		return location;
	else
		return location + scale * tan(M_PI * unif_rand());
}

double rwrpcauchy(double location, double rho) {
	double result;

	if (rho < 0.0 || rho > 1.0)	ML_ERR_return_NAN;

	if(rho == 0)
		result = unif_rand() * M2_PI;
	else if(rho == 1)
		result = location;
	else {
		result = myfmod(rcauchy(location, -log(rho)), M2_PI);
	}
	return result;
}

double rwrpnormal (double mu, double rho) {
    if (rho == 0)
        return runif(0, M2_PI);
    else if (rho == 1)
        return mu;
    else {
        return myfmod(rnorm(mu, sqrt(-2 * log(rho))), M2_PI);
    }
}

//------------------------------------------------------------------------------

/* Reference:
 * R. C. H. Cheng (1978).
 * Generating beta variates with nonintegral shape parameters.
 * Communications of the ACM 21, 317-322.
 * (Algorithms BB and BC)
 */

#include <float.h>	/* DBL_MAX_EXP, LDBL_MAX_EXP */

#define expmax	(DBL_MAX_EXP * M_LN2)/* = log(DBL_MAX) */

double rbeta(double aa, double bb) {
	double a, b, alpha;
	double r, s, t, u1, u2, v, w, y, z;

	int qsame;
	/* FIXME:  Keep Globals (properly) for threading */
	/* Uses these GLOBALS to save time when many rv's are generated : */
	static double beta, gamma, delta, k1, k2;
	static double olda = -1.0;
	static double oldb = -1.0;

	if (aa <= 0. || bb <= 0. || (!R_FINITE(aa) && !R_FINITE(bb)))
		return ML_NAN;

	if (!R_FINITE(aa))
		return 1.0;

	if (!R_FINITE(bb))
		return 0.0;

	/* Test if we need new "initializing" */
	qsame = (olda == aa) && (oldb == bb);
	if (!qsame) {
		olda = aa;
		oldb = bb;
	}

	a = fmin2(aa, bb);
	b = fmax2(aa, bb); /* a <= b */
	alpha = a + b;

#define v_w_from__u1_bet(AA) 			\
	    v = beta * log(u1 / (1.0 - u1));	\
	    if (v <= expmax) {			\
	    	w = AA * exp(v);		\
	    	if(!R_FINITE(w)) w = DBL_MAX;	\
	    } else				\
			w = DBL_MAX

	if (a <= 1.0) { /* --- Algorithm BC --- */

		/* changed notation, now also a <= b (was reversed) */

		if (!qsame) { /* initialize */
			beta = 1.0 / a;
			delta = 1.0 + b - a;
			k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
			k2 = 0.25 + (0.5 + 0.25 / delta) * a;
		}
		/* FIXME: "do { } while()", but not trivially because of "continue"s:*/
		for (;;) {
			u1 = unif_rand();
			u2 = unif_rand();
			if (u1 < 0.5) {
				y = u1 * u2;
				z = u1 * y;
				if (0.25 * u2 + z - y >= k1)
					continue;
			} else {
				z = u1 * u1 * u2;
				if (z <= 0.25) {
					v_w_from__u1_bet(b);
					break;
				}
				if (z >= k2)
					continue;
			}

			v_w_from__u1_bet(b);

			if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
				break;
		}
		return (aa == a) ? a / (a + w) : w / (a + w);

	} else { /* Algorithm BB */

		if (!qsame) { /* initialize */
			beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
			gamma = a + 1.0 / beta;
		}
		do {
			u1 = unif_rand();
			u2 = unif_rand();

			v_w_from__u1_bet(a);

			z = u1 * u1 * u2;
			r = gamma * v - 1.3862944;
			s = a + r - w;
			if (s + 2.609438 >= 5.0 * z)
				break;
			t = log(z);
			if (s > t)
				break;
		} while (r + alpha * log(alpha / (b + w)) < t);

		return (aa != a) ? b / (b + w) : w / (b + w);
	}
#undef v_w_from__u1_bet
}

int chebyshev_init(double *dos, int nos, double eta) {
	int i, ii;
	double err;

	if (nos < 1)
		return 0;

	err = 0.0;
	i = 0; /* just to avoid compiler warnings */
	for (ii = 1; ii <= nos; ii++) {
		i = nos - ii;
		err += fabs(dos[i]);
		if (err > eta) {
			return i;
		}
	}
	return i;
}

double chebyshev_eval(double x, const double *a, const int n) {
	double b0, b1, b2, twox;
	int i;

	if (n < 1 || n > 1000)
		return ML_NAN;

	if (x < -1.1 || x > 1.1)
		return ML_NAN;

	twox = x * 2;
	b2 = b1 = 0;
	b0 = 0;
	for (i = 1; i <= n; i++) {
		b2 = b1;
		b1 = b0;
		b0 = twox * b1 - b2 + a[n - i];
	}
	return (b0 - b2) * 0.5;
}

// From R Mathlib (/src/nmath/log1p.c)
#ifndef HAVE_LOG1P
double log1p(double x)
{
	/* series for log1p on the interval -.375 to .375
	 *				     with weighted error   6.35e-32
	 *				      log weighted error  31.20
	 *			    significant figures required  30.93
	 *				 decimal places required  32.01
	 */
	const static double alnrcs[43] = {
		+.10378693562743769800686267719098e+1,
		-.13364301504908918098766041553133e+0,
		+.19408249135520563357926199374750e-1,
		-.30107551127535777690376537776592e-2,
		+.48694614797154850090456366509137e-3,
		-.81054881893175356066809943008622e-4,
		+.13778847799559524782938251496059e-4,
		-.23802210894358970251369992914935e-5,
		+.41640416213865183476391859901989e-6,
		-.73595828378075994984266837031998e-7,
		+.13117611876241674949152294345011e-7,
		-.23546709317742425136696092330175e-8,
		+.42522773276034997775638052962567e-9,
		-.77190894134840796826108107493300e-10,
		+.14075746481359069909215356472191e-10,
		-.25769072058024680627537078627584e-11,
		+.47342406666294421849154395005938e-12,
		-.87249012674742641745301263292675e-13,
		+.16124614902740551465739833119115e-13,
		-.29875652015665773006710792416815e-14,
		+.55480701209082887983041321697279e-15,
		-.10324619158271569595141333961932e-15,
		+.19250239203049851177878503244868e-16,
		-.35955073465265150011189707844266e-17,
		+.67264542537876857892194574226773e-18,
		-.12602624168735219252082425637546e-18,
		+.23644884408606210044916158955519e-19,
		-.44419377050807936898878389179733e-20,
		+.83546594464034259016241293994666e-21,
		-.15731559416479562574899253521066e-21,
		+.29653128740247422686154369706666e-22,
		-.55949583481815947292156013226666e-23,
		+.10566354268835681048187284138666e-23,
		-.19972483680670204548314999466666e-24,
		+.37782977818839361421049855999999e-25,
		-.71531586889081740345038165333333e-26,
		+.13552488463674213646502024533333e-26,
		-.25694673048487567430079829333333e-27,
		+.48747756066216949076459519999999e-28,
		-.92542112530849715321132373333333e-29,
		+.17578597841760239233269760000000e-29,
		-.33410026677731010351377066666666e-30,
		+.63533936180236187354180266666666e-31,
	};

#ifdef NOMORE_FOR_THREADS
	static int nlnrel = 0;
	static double xmin = 0.0;

	if (xmin == 0.0) xmin = -1 + sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)); */
	if (nlnrel == 0) /* initialize chebychev coefficients */
	nlnrel = chebyshev_init(alnrcs, 43, DBL_EPSILON/20);/*was .1*d1mach(3)*/
#else
# define nlnrel 22
	const static double xmin = -0.999999985;
	/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
#endif

	if (x == 0.) return 0.;/* speed */
	if (x == -1) return ML_NEGINF;
	if (x < -1) return ML_NAN;

	if (fabs(x) <= .375) {
		/* Improve on speed (only);
		 again give result accurate to IEEE double precision: */
		if(fabs(x) < .5 * DBL_EPSILON)
		return x;

		if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
		return x * (1 - .5 * x);
		/* else */
		return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
	}
	/* else */
	if (x < xmin) {
		/* answer less than half precision because x too near -1 */
		perror("full precision may not have been achieved in 'log1p'");
	}
	return log(1 + x);
}
#endif
#ifndef HAVE_EXPM1
double expm1(double x) {
	double y, a = fabs(x);

	if (a < DBL_EPSILON)
		return x;
	if (a > 0.697)
		return exp(x) - 1; /* negligible cancellation */

	if (a > 1e-8)
		y = exp(x) - 1;
	else
		/* Taylor expansion, more accurate in this range */
		y = (x / 2 + 1) * x;

	/* Newton step for solving   log(1 + y) = x   for y : */
	/* WARNING: does not work for y ~ -1: bug in 1.5.0 */
	y -= (1 + y) * (log1p(y) - x);
	return y;
}
#endif


//#include "toms708.c"
//
//double pbeta_raw(double x, double pin, double qin, int lower_tail, int log_p) {
//	double x1 = 0.5 - x + 0.5, w, wc;
//	int ierr;
//	bratio(pin, qin, x, x1, &w, &wc, &ierr, log_p); /* -> ./toms708.c */
//	/* ierr = 8 is about inaccuracy in extreme cases */
//	//if (ierr && (ierr != 8 || log_p))
////	MATHLIB_WARNING(_("pbeta_raw() -> bratio() gave error code %d"), ierr);
//	return lower_tail ? w : wc;
//} /* pbeta_raw() */
//
//double pbeta(double x, double pin, double qin, int lower_tail, int log_p) {
//#ifdef IEEE_754
//	if (ISNAN(x) || ISNAN(pin) || ISNAN(qin)) return x + pin + qin;
//#endif
//
//	if (pin <= 0 || qin <= 0)
//		ML_ERR_return_NAN;
//
//	if (x <= 0)
//		return R_DT_0;
//	if (x >= 1)
//		return R_DT_1;
//	return pbeta_raw(x, pin, qin, lower_tail, log_p);
//}

/*double pbinom(double x, double n, double p, int lower_tail, int log_p) {
 #ifdef IEEE_754
 if (ISNAN(x) || ISNAN(n) || ISNAN(p))
 return x + n + p;
 if (!R_FINITE(n) || !R_FINITE(p)) ML_ERR_return_NAN;

 #endif
 if (R_D_nonint(n))
 ML_ERR_return_NAN;
 n = R_D_forceint(n);
 PR#8560: n=0 is a valid value
 if (n < 0 || p < 0 || p > 1)
 ML_ERR_return_NAN;

 if (x < 0)
 return R_DT_0;
 x = floor(x + 1e-7);
 if (n <= x)
 return R_DT_1;
 return pbeta(p, x + 1, n - x, !lower_tail, log_p);
 }*/

/*static double do_search(double y, double *z, double p, double n, double pr,
 double incr) {
 if (*z >= p) {
 search to the left
 REPEAT {
 double newz;
 if (y == 0 || (newz = pbinom(y - incr, n, pr, l._t.TRUE, log_p FALSE)) < p)
 return y;
 y = fmax2(0, y - incr);
 *z = newz;
 }
 } else {  search to the right
 REPEAT {
 y = fmin2(y + incr, n);
 if (y == n || (*z = pbinom(y, n, pr, l._t.TRUE, log_pFALSE)) >= p)
 return y;
 }
 }
 return NAN;
 }*/

/*
 double qbinom(double p, double n, double pr, int lower_tail, int log_p) {
 double q, mu, sigma, gamma, z, y;

 #ifdef IEEE_754
 if (ISNAN(p) || ISNAN(n) || ISNAN(pr))
 return p + n + pr;
 #endif
 if (!R_FINITE(n) || !R_FINITE(pr))
 ML_ERR_return_NAN;
 if log_p is true, p = -Inf is a legitimate value
 if (!R_FINITE(p) && !log_p)
 ML_ERR_return_NAN;

 if (n != floor(n + 0.5))
 ML_ERR_return_NAN;
 if (pr < 0 || pr > 1 || n < 0)
 ML_ERR_return_NAN;

 R_Q_P01_boundaries(p, 0, n);

 if (pr == 0. || n == 0)
 return 0.;

 q = 1 - pr;
 if (q == 0.)
 return n;  covers the full range of the distribution
 mu = n * pr;
 sigma = sqrt(n * pr * q);
 gamma = (q - pr) / sigma;

 Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
 * FIXME: This is far from optimal [cancellation for p ~= 1, etc]:
 if (!lower_tail || log_p) {
 p = R_DT_qIv(p);  need check again (cancellation!):
 if (p == 0.)
 return 0.;
 if (p == 1.)
 return n;
 }
 temporary hack --- FIXME ---
 if (p + 1.01 * DBL_EPSILON >= 1.)
 return n;

 y := approx.value (Cornish-Fisher expansion) :
 z = qnorm5(p, 0., 1., lower_tailTRUE, log_pFALSE);
 y = floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);

 if (y > n)  way off
 y = n;

 #ifdef DEBUG_qbinom
 REprintf("  new (p,1-p)=(%7g,%7g), z=qnorm(..)=%7g, y=%5g\n", p, 1-p, z, y);
 #endif
 z = pbinom(y, n, pr, lower_tailTRUE, log_pFALSE);

 fuzz to ensure left continuity:
 p *= 1 - 64 * DBL_EPSILON;

 if (n < 1e5)
 return do_search(y, &z, p, n, pr, 1);
 Otherwise be a bit cleverer in the search
 {
 double incr = floor(n * 0.001), oldincr;
 do {
 oldincr = incr;
 y = do_search(y, &z, p, n, pr, incr);
 incr = fmax2(1., floor(incr / 100.));
 } while (oldincr > 1 && incr > n * 1e-15);
 return y;
 }
 }
 */
/*

 double rbinom(double nin, double pp) {
 FIXME: These should become THREAD_specific globals :

 static double c, fm, npq, p1, p2, p3, p4, qn;
 static double xl, xll, xlr, xm, xr;

 static double psave = -1.0;
 static int nsave = -1;
 static int m;

 double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
 double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
 int i, ix, k, n;

 if (!R_FINITE(nin))
 return ML_NAN;
 r = floor(nin + 0.5);
 if (r != nin)
 return ML_NAN;
 if (!R_FINITE(pp) ||
 n=0, p=0, p=1 are not errors <TSL>
 r < 0 || pp < 0. || pp > 1.)
 return ML_NAN;

 if (r == 0 || pp == 0.)
 return 0;
 if (pp == 1.)
 return r;

 if (r >= INT_MAX) evade integer overflow,
 and r == INT_MAX gave only even values
 return qbinom(unif_rand(), r, pp, lower_tail0, log_p0);
 else
 n = (int) r;

 p = fmin2(pp, 1. - pp);
 q = 1. - p;
 np = n * p;
 r = p / q;
 g = r * (n + 1);

 Setup, perform only when parameters change [using static (globals):

 FIXING: Want this thread safe
 -- use as little (thread globals) as possible

 if (pp != psave || n != nsave) {
 psave = pp;
 nsave = n;
 if (np < 30.0) {
 inverse cdf logic for mean less than 30
 qn = pow(q, (double) n);
 goto L_np_small;
 } else {
 ffm = np + p;
 m = (int)  ffm;
 fm = m;
 npq = np * q;
 p1 = (int) (2.195 * sqrt(npq) - 4.6 * q) + 0.5;
 xm = fm + 0.5;
 xl = xm - p1;
 xr = xm + p1;
 c = 0.134 + 20.5 / (15.3 + fm);
 al = (ffm - xl) / (ffm - xl * p);
 xll = al * (1.0 + 0.5 * al);
 al = (xr - ffm) / (xr * q);
 xlr = al * (1.0 + 0.5 * al);
 p2 = p1 * (1.0 + c + c);
 p3 = p2 + c / xll;
 p4 = p3 + c / xlr;
 }
 } else if (n == nsave) {
 if (np < 30.0)
 goto L_np_small;
 }

 -------------------------- np = n*p >= 30 : -------------------
 REPEAT {
 u = unif_rand() * p4;
 v = unif_rand();
 triangular region
 if (u <= p1) {
 ix = (int) (xm - p1 * v + u);
 goto finis;
 }
 parallelogram region
 if (u <= p2) {
 x = xl + (u - p1) / c;
 v = v * c + 1.0 - fabs(xm - x) / p1;
 if (v > 1.0 || v <= 0.)
 continue;
 ix = (int) x;
 } else {
 if (u > p3) {  right tail
 ix = (int) (xr - log(v) / xlr);
 if (ix > n)
 continue;
 v = v * (u - p3) * xlr;
 } else { left tail
 ix = (int) (xl + log(v) / xll);
 if (ix < 0)
 continue;
 v = v * (u - p2) * xll;
 }
 }
 determine appropriate way to perform accept/reject test
 k = abs(ix - m);
 if (k <= 20 || k >= npq / 2 - 1) {
 explicit evaluation
 f = 1.0;
 if (m < ix) {
 for (i = m + 1; i <= ix; i++)
 f *= (g / i - r);
 } else if (m != ix) {
 for (i = ix + 1; i <= m; i++)
 f /= (g / i - r);
 }
 if (v <= f)
 goto finis;
 } else {
 squeezing using upper and lower bounds on log(f(x))
 amaxp = (k / npq)
 * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
 ynorm = -k * k / (2.0 * npq);
 alv = log(v);
 if (alv < ynorm - amaxp)
 goto finis;
 if (alv <= ynorm + amaxp) {
 stirling's formula to machine accuracy
 for the final acceptance/rejection test
 x1 = ix + 1;
 f1 = fm + 1.0;
 z = n + 1 - fm;
 w = n - ix + 1.0;
 z2 = z * z;
 x2 = x1 * x1;
 f2 = f1 * f1;
 w2 = w * w;
 if (alv
 <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w)
 + (ix - m) * log(w * p / (x1 * q))
 + (13860.0
 - (462.0
 - (132.0
 - (99.0 - 140.0 / f2)
 / f2) / f2) / f2)
 / f1 / 166320.0
 + (13860.0
 - (462.0
 - (132.0
 - (99.0 - 140.0 / z2)
 / z2) / z2) / z2)
 / z / 166320.0
 + (13860.0
 - (462.0
 - (132.0
 - (99.0 - 140.0 / x2)
 / x2) / x2) / x2)
 / x1 / 166320.0
 + (13860.0
 - (462.0
 - (132.0
 - (99.0 - 140.0 / w2)
 / w2) / w2) / w2)
 / w / 166320.)
 goto finis;
 }
 }
 }

 L_np_small:
 ---------------------- np = n*p < 30 : -------------------------

 REPEAT {
 ix = 0;
 f = qn;
 u = unif_rand();
 REPEAT {
 if (u < f)
 goto finis;
 if (ix > 110)
 break;
 u -= f;
 ix++;
 f *= (g / ix - r);
 }
 }
 finis: if (psave > 0.5)
 ix = n - ix;
 return (double) ix;
 }

 */
