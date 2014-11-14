/*
 * Rmathdefines.h
 *
 *  Created on: 2 wrz 2013
 *      Author: s12kb2
 */

#ifndef RMATHDEFINES_H_
#define RMATHDEFINES_H_


/// R stuff:
#define REPEAT for(;;)



//#ifdef __BORLANDC__
//#	include <limits>
//#	ifndef INFINITY
//#		define INFINITY std::numeric_limits<double>::infinity()
//#	endif
//#	ifndef NAN
//#		define NAN std::numeric_limits<double>::quiet_NaN()
//#	endif
//#	define isnan(x) _isnan(x)
//#	define ML_POSINF	INFINITY
//#	define ML_NEGINF	((-1.0) / 0.0)
//#	define ML_NAN		NAN
//#elif defined(_MSC_VER)
//#	include <float.h>
//#	ifndef INFINITY
//#		define INFINITY (DBL_MAX+DBL_MAX)
//#	endif
//#	ifndef NAN
//#		define NAN (INFINITY-INFINITY)
//#	endif
//#	define ML_POSINF INFINITY
//#	define ML_NEGINF	-INFINITY
//#	define ML_NAN		NAN
//#else
//#	define ML_POSINF	(1.0 / 0.0)
//#	define ML_NEGINF	((-1.0) / 0.0)
//#	define ML_NAN		(0.0 / 0.0)
//#endif


#ifdef __USE_ISOC99
//#	ifdef _GLIBCXX_CMATH
//#		define R_FINITE(x)  std::isfinite(x)
//#	else
#		define R_FINITE(x)  std::isfinite(x)
//#	endif
#else
#	define R_FINITE(x)  _finite(x)
#endif

#define A_0	-0.5
#define A_1	 0.3333333
#define A_2	-0.2500068
#define A_3	 0.2000118
#define A_4	-0.1661269
#define A_5	 0.1421878
#define A_6	-0.1384794
#define A_7	 0.1250060

#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */


#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
				 : R_D_Cval(p))
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
				 : R_D_Lval(p))


#define  ML_ERR_return_NAN return NAN;

#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)






#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    return ML_NAN;				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    return ML_NAN;				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }


#endif /* RMATHDEFINES_H_ */
