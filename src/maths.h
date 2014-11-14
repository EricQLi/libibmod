//---------------------------------------------------------------------------

#ifndef mathsH
#define mathsH
//---------------------------------------------------------------------------
//#include <stddef.h>
//#include <stdio.h>
#include <cstdio>
//#include "uint.h"
#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif


#ifdef _MSC_VER
#	ifndef _USE_MATH_DEFINES
#		define _USE_MATH_DEFINES
#	endif
#	define fmax std::max
#	define fmin std::min
#endif
#include <cmath>


#ifdef __BORLANDC__
#	include <limits>
#	ifndef INFINITY
#		define INFINITY std::numeric_limits<double>::infinity()
#	endif
#	ifndef NAN
#		define NAN std::numeric_limits<double>::quiet_NaN()
#	endif
#	define isnan(x) _isnan(x)
#	define ML_POSINF	INFINITY
#	define ML_NEGINF	((-1.0) / 0.0)
#	define ML_NAN		NAN
#elif defined(_MSC_VER)
#	include <float.h>
#	ifndef INFINITY
#		define INFINITY (DBL_MAX+DBL_MAX)
#	endif
#	ifndef NAN
#		define NAN (INFINITY-INFINITY)
#	endif
#	define ML_POSINF INFINITY
#	define ML_NEGINF	-INFINITY
#	define ML_NAN		NAN
#	define ISNAN _isnan
#else
#	define ML_POSINF	(1.0 / 0.0)
#	define ML_NEGINF	((-1.0) / 0.0)
#	define ML_NAN		(0.0 / 0.0)
#	define ISNAN std::isnan
#endif


//#ifdef _MSC_VER
//typedef __int32 int32_t;
//typedef unsigned __int32 uint32_t;
//typedef __int64 int64_t;
//typedef unsigned __int64 uint64_t;
//#else

#ifdef __GNUC__
#define HAVE_LOG1P
#define HAVE_EXPM1
#endif

#include <stdint.h>
//#endif

unsigned long int random_seed(void);

long long int get_timediff(bool start = false);
long long int get_time_msec(void);



// Extra math constants (some from R mathlib
#define M2_PI 6.283185307179586
#define M_2PI M2_PI
#define M_1_SQRT_2PI 0.39894228040143270
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#define M_2_E 0.367879441171442
#define M_SQRT_2PI 2.506628274631
#define M_SQRT_32  5.656854249492380195206754896838 /* sqrt(32) */
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667

#define mysqrt(x) sqrt(x)
#define mylog(x) log(x)
#define myatan(x) atan(x)
#define mypow(x, y) ((y == 2)? x * x : pow(x, y))

double unif_rand(void);
double runif(double a, double b);

/**
 * random integer min <= x < max
 * @param min
 * @param max
 * @return
 */
inline int irand(int min, int max) {
	return (int) floor(runif(min, max));
}
inline long int irand(long int min, long int max) {
	return (long int) floor(runif(min, max));
}
inline size_t irand(size_t min, size_t max) {
	return (size_t) floor(runif(min, max));
}

/**
 * Compute the integer remainder from the division of numerator by denominator.
 * Differs from '%' in that the result is always positive.
 * @param a
 * @param b
 * @return
 */
inline int imod(const int a, const int b) {
	int res = a % b;
	if (res < 0)
		return res + b;
	return res;
}



void rand_perm(int k, int *a);
int rand_pos(int n, int size, int* result, int start_at = -1);

//double dist(double x, double y);
double angle(double x, double y);
double myfmod(double x1, double x2)  __attribute__ ((pure));
int sample1_with_prob(double *prob, const size_t nprob);
void sample_with_prob(size_t * result, size_t n, double *prob, const size_t nprob);

void sincos(const double a, double & sinv, double & cosv);

size_t sample1(size_t n);

extern inline double squared(double x) {
	return x * x;
}


#ifdef _MSC_VER
#define dist(a, b) _hypot(a, b)
#else
#define dist(a, b) hypot(a, b)
#endif

/**
 * Convert floating-point [0,1] to boolean
 * @param x a floating-point number
 * @return true if x > 0.5
 */
extern inline bool to_bool(double x) {
	return x > .5;
}

/**
 * Probabilistic rounding. Round numeric x.y to x with probability 0.y, or to x+1 with probability 1-0.y
 * @param x floating point number
 * @return integer after rounding
 */
extern inline int round_prob(double x) {
	double flop, intp;
	flop = modf(x, &intp);
	return (int) intp + ((unif_rand() < flop) ? 1 : 0);
}

/**
 * Scale a number to [0,1] limits
 * @param x the number to be scaled
 * @param xmin upper bound
 * @param xmax lower bound
 * @return number within 0 and 1
 */
inline double to01(double x, double xmin, double xmax) {
	return (x - xmin) / (xmax - xmin);
}

inline double pow2(double x) {
	return x * x;
}

template<typename T>
T array_sum(T* arr, size_t n) {
	T sum = 0;
	for (size_t i = 0; i < n; i++)
		sum += arr[i];
	return sum;
}


template<class T>
int intersection(T* &result, const T *arr1, const T *arr2, int m, int n) {
	int i = 0, j = 0, len = 0;
	T* tmp = new T[(m > n) ? m : n];

	while (i < m && j < n) {
		if (arr1[i] < arr2[j])
			++i;
		else if (arr2[j] < arr1[i])
			++j;
		else { /* if arr1[i] == arr2[j] */
			tmp[len++] = arr2[j++];
			++i;
		}
	}
	result = new T[len];
	for (int i = 0; i < len; ++i) result[i] = tmp[i];
	delete[] tmp;
	return len;
}


template<class T>
int intersection_len(const T *arr1, const T *arr2, const int m, const int n) {
	int i = 0, j = 0, len = 0;
	while (i < m && j < n) {
		if (arr1[i] < arr2[j])
			++i;
		else if (arr2[j] < arr1[i])
			++j;
		else { /* if arr1[i] == arr2[j] */
			++len; ++j; ++i;
		}
	}
	return len;
}

void quantile(double* x, const size_t n, const double * probs, const size_t np,
		double * qs);

//void quantile(double* x, const size_t n, const double prob, double & qs);
/**
 * Calculate mean and std. dev. simultanously
 * @param x numeric array
 * @param n array length
 * @param result double[2] to store the result (first value is mean, second is std. dev.)
 */
void mean_sd(double* x, const size_t n, double * result);

/**
 * Finds minimum and maximum value in an array
 * @param x	numeric array
 * @param n	array length
 * @param min numeric to store the minimum
 * @param max numeric to store the maximum
 */
template<class T>
void minmax(T* x, size_t n, T & min, T & max) {
	max = x[0], min = x[0];
	for (size_t i = 0; i < n; i++) {
		if (max < x[i])
			max = x[i];
		if (min > x[i])
			min = x[i];
	}
}

template<class T>
inline T fit(T x, T xmin, T xmax) {
	if (x < xmin)
		return xmin;
	if (x > xmax)
		return xmax;
	return x;
	// return std::max(std::min(x, xmax), xmin);
}

template<class T>
inline void fiti(T & x, const T & xmin, const T & xmax) {
	if (x < xmin)
		x = xmin;
	else if (x > xmax)
		x = xmax;
}


template<class T>
inline void torus(T & x, const T & xmin, const T & xmax) {
	if (x < xmin)
		x = xmax;
	else if (x > xmax)
		x = xmin;
}

template<class T>
inline bool is_within(T x, T xmin, T xmax) {
	if (x < xmin)
		return false;
	if (x > xmax)
		return false;
	return true;
	// return std::max(std::min(x, xmax), xmin);
}


#define NUMSUM_QTL_SIZE 7
#define NUMSUM_QTL_PROBS { 0, .025, .25, .5, .75, .975, 1 }

typedef double numeric_summary[NUMSUM_QTL_SIZE + 2];

/**
 * Calculate numerical summary (similar to the one in R)
 * @param values
 * @param size
 * @param output calculates: min, q1, median, q2, max, mean, sd
 */
void compute_summary(double* values, size_t size, numeric_summary output);


template<class T>
void histogram1(T* x, size_t n, int* bins, size_t n_bins, T xmin, T xmax) {
//	double real_max, real_min;
//	minmax(x, n, real_min, real_max);
//	if(xmin > real_min || xmax < real_max) {
//		fprintf(stderr, "Warning: 'xmin' or 'xmax' do not span range of 'x'. \n");
//		xmin = real_min;
//		xmax = real_max;
//	}
	bool warning_given = false;
	for (size_t i = 0; i < n_bins; ++i)
		bins[i] = 0;
	double d = (xmax - xmin) / (double) (n_bins);
	for (size_t i = 0; i < n; ++i) {
		double z = (x[i] - xmin) / d;
		if (is_within(x[i], xmin, xmax)) {
			bins[(size_t) ceil(z) - ((z == 0.0) ? 0 : 1)]++;
		} else if (!warning_given) {
			fprintf(stderr,
					"Warning: some 'x' not counted. 'xmin' or 'xmax' do not span range of 'x'. \n");
			warning_given = true;
		}
	}
}

template<class T>
void histogram2(T* x, size_t n, int* bins, size_t n_bins, T xmin, T xmax) {
// TODO: log bins
	//double real_max, real_min;
	if (xmax == xmin)
		xmax += 1;
	for (size_t i = 0; i < n_bins; ++i)
		bins[i] = 0;
	double d;
	d = (xmax - xmin) / (double) (n_bins);
	for (size_t i = 0; i < n; ++i) {
		double z = (fit(x[i], xmin, xmax) - xmin) / d;
		bins[(size_t) ceil(z) - ((z == 0.) ? 0 : 1)]++;
	}
}

/* R mathlib functions */
int imax2(int x, int y);
int imin2(int x, int y);

#ifndef HAVE_LOG1P
double log1p(double x);
#endif
#ifndef HAVE_EXPM1
double expm1(double x);
#endif


/* R r* random numbers */

double runif(double a, double b);
double rpois(double mu);
double rgamma(double a, double scale);
double rnorm(double mu, double sigma);
double rexp(double scale);
double rchisq(double df);
double rinvgauss(double mu, double lambda);

double rnbinom(double size, double prob);
double rnbinom_mu(double size, double mu);

double rbinom(double nin, double pp);

double rbeta(double aa, double bb);

double rcauchy(double location, double scale);
double rwrpcauchy(double location, double rho);
double rwrpnormal (double mu, double rho);

#endif
