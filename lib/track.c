#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "sync.h"
#include "track.h"
#include "base.h"

#ifndef M_PI
#define M_PI 3.1415926545
#endif

#ifndef M_PI_2
#define M_PI_2 (M_PI * 2.0)
#endif

// Modeled after the line y = x
double LinearInterpolation(double p) {
	return p;
}

// Modeled after the parabola y = x^2
double QuadraticEaseIn(double p) {
	return p * p;
}

// Modeled after the parabola y = -x^2 + 2x
double QuadraticEaseOut(double p) {
	return -(p * (p - 2));
}

// Modeled after the piecewise quadratic
// y = (1/2)((2x)^2)             ; [0, 0.5)
// y = -(1/2)((2x-1)*(2x-3) - 1) ; [0.5, 1]
double QuadraticEaseInOut(double p) {
	if(p < 0.5) {
		return 2 * p * p;
	} else {
		return (-2 * p * p) + (4 * p) - 1;
	}
}

// Modeled after the cubic y = x^3
double CubicEaseIn(double p) {
	return p * p * p;
}

// Modeled after the cubic y = (x - 1)^3 + 1
double CubicEaseOut(double p) {
	double f = (p - 1);
	return f * f * f + 1;
}

// Modeled after the piecewise cubic
// y = (1/2)((2x)^3)       ; [0, 0.5)
// y = (1/2)((2x-2)^3 + 2) ; [0.5, 1]
double CubicEaseInOut(double p) {
	if(p < 0.5) {
		return 4 * p * p * p;
	} else {
		double f = ((2 * p) - 2);
		return 0.5 * f * f * f + 1;
	}
}

// Modeled after the quartic x^4
double QuarticEaseIn(double p) {
	return p * p * p * p;
}

// Modeled after the quartic y = 1 - (x - 1)^4
double QuarticEaseOut(double p) {
	double f = (p - 1);
	return f * f * f * (1 - p) + 1;
}

// Modeled after the piecewise quartic
// y = (1/2)((2x)^4)        ; [0, 0.5)
// y = -(1/2)((2x-2)^4 - 2) ; [0.5, 1]
double QuarticEaseInOut(double p)  {
	if(p < 0.5) {
		return 8 * p * p * p * p;
	} else {
		double f = (p - 1);
		return -8 * f * f * f * f + 1;
	}
}

// Modeled after the quintic y = x^5
double QuinticEaseIn(double p)  {
	return p * p * p * p * p;
}

// Modeled after the quintic y = (x - 1)^5 + 1
double QuinticEaseOut(double p)  {
	double f = (p - 1);
	return f * f * f * f * f + 1;
}

// Modeled after the piecewise quintic
// y = (1/2)((2x)^5)       ; [0, 0.5)
// y = (1/2)((2x-2)^5 + 2) ; [0.5, 1]
double QuinticEaseInOut(double p)  {
	if(p < 0.5) {
		return 16 * p * p * p * p * p;
	} else {
		double f = ((2 * p) - 2);
		return  0.5 * f * f * f * f * f + 1;
	}
}

// Modeled after quarter-cycle of sine wave
double SineEaseIn(double p) {
	return sin((p - 1) * M_PI_2) + 1;
}

// Modeled after quarter-cycle of sine wave (different phase)
double SineEaseOut(double p) {
	return sin(p * M_PI_2);
}

// Modeled after half sine wave
double SineEaseInOut(double p) {
	return 0.5 * (1 - cos(p * M_PI));
}

// Modeled after shifted quadrant IV of unit circle
double CircularEaseIn(double p) {
	return 1 - sqrt(1 - (p * p));
}

// Modeled after shifted quadrant II of unit circle
double CircularEaseOut(double p) {
	return sqrt((2 - p) * p);
}

// Modeled after the piecewise circular function
// y = (1/2)(1 - sqrt(1 - 4x^2))           ; [0, 0.5)
// y = (1/2)(sqrt(-(2x - 3)*(2x - 1)) + 1) ; [0.5, 1]
double CircularEaseInOut(double p) {
	if(p < 0.5) {
		return 0.5 * (1 - sqrt(1 - 4 * (p * p)));
	} else {
		return 0.5 * (sqrt(-((2 * p) - 3) * ((2 * p) - 1)) + 1);
	}
}

// Modeled after the exponential function y = 2^(10(x - 1))
double ExponentialEaseIn(double p) {
	return (p == 0.0) ? p : pow(2, 10 * (p - 1));
}

// Modeled after the exponential function y = -2^(-10x) + 1
double ExponentialEaseOut(double p) {
	return (p == 1.0) ? p : 1 - pow(2, -10 * p);
}

// Modeled after the piecewise exponential
// y = (1/2)2^(10(2x - 1))         ; [0,0.5)
// y = -(1/2)*2^(-10(2x - 1))) + 1 ; [0.5,1]
double ExponentialEaseInOut(double p) {
	if(p == 0.0 || p == 1.0) return p;
	
	if(p < 0.5) {
		return 0.5 * pow(2, (20 * p) - 10);
	} else {
		return -0.5 * pow(2, (-20 * p) + 10) + 1;
	}
}

// Modeled after the damped sine wave y = sin(13pi/2*x)*pow(2, 10 * (x - 1))
double ElasticEaseIn(double p) {
	return sin(13 * M_PI_2 * p) * pow(2, 10 * (p - 1));
}

// Modeled after the damped sine wave y = sin(-13pi/2*(x + 1))*pow(2, -10x) + 1
double ElasticEaseOut(double p) {
	return sin(-13 * M_PI_2 * (p + 1)) * pow(2, -10 * p) + 1;
}

// Modeled after the piecewise exponentially-damped sine wave:
// y = (1/2)*sin(13pi/2*(2*x))*pow(2, 10 * ((2*x) - 1))      ; [0,0.5)
// y = (1/2)*(sin(-13pi/2*((2x-1)+1))*pow(2,-10(2*x-1)) + 2) ; [0.5, 1]
double ElasticEaseInOut(double p) {
	if(p < 0.5) {
		return 0.5 * sin(13 * M_PI_2 * (2 * p)) * pow(2, 10 * ((2 * p) - 1));
	} else {
		return 0.5 * (sin(-13 * M_PI_2 * ((2 * p - 1) + 1)) * pow(2, -10 * (2 * p - 1)) + 2);
	}
}

// Modeled after the overshooting cubic y = x^3-x*sin(x*pi)
double BackEaseIn(double p) {
	return p * p * p - p * sin(p * M_PI);
}

// Modeled after overshooting cubic y = 1-((1-x)^3-(1-x)*sin((1-x)*pi))
double BackEaseOut(double p) {
	double f = (1 - p);
	return 1 - (f * f * f - f * sin(f * M_PI));
}

// Modeled after the piecewise overshooting cubic function:
// y = (1/2)*((2x)^3-(2x)*sin(2*x*pi))           ; [0, 0.5)
// y = (1/2)*(1-((1-x)^3-(1-x)*sin((1-x)*pi))+1) ; [0.5, 1]
double BackEaseInOut(double p) {
	if(p < 0.5) {
		double f = 2 * p;
		return 0.5 * (f * f * f - f * sin(f * M_PI));
	} else {
		double f = (1 - (2*p - 1));
		return 0.5 * (1 - (f * f * f - f * sin(f * M_PI))) + 0.5;
	}
}

double BounceEaseOut(double p) {
	if(p < 4/11.0) {
		return (121 * p * p)/16.0;
	} else if(p < 8/11.0) {
		return (363/40.0 * p * p) - (99/10.0 * p) + 17/5.0;
	} else if(p < 9/10.0) {
		return (4356/361.0 * p * p) - (35442/1805.0 * p) + 16061/1805.0;
	} else {
		return (54/5.0 * p * p) - (513/25.0 * p) + 268/25.0;
	}
}

double BounceEaseIn(double p) {
	return 1 - BounceEaseOut(1 - p);
}

double BounceEaseInOut(double p) {
	if(p < 0.5) {
		return 0.5 * BounceEaseIn(p*2);
	} else {
		return 0.5 * BounceEaseOut(p * 2 - 1) + 0.5;
	}
}

static double interpolate(double t, enum key_type e) {
	switch (e) {
		case KEY_STEP:
			return 0.0;
		case KEY_LINEAR:
			return t;
		case KEY_SMOOTH:
			return t * t * (3.0 - 2.0 * t);
		case KEY_RAMP:
			return pow(t, 2.0);

		case KEY_IN_QUAD:
			return QuadraticEaseIn(t);
		case KEY_OUT_QUAD:
			return QuadraticEaseOut(t);
		case KEY_IN_OUT_QUAD:
			QuadraticEaseInOut(t);
		case KEY_IN_CUBIC:
			return CubicEaseIn(t);
		case KEY_OUT_CUBIC:
			return CubicEaseOut(t);
		case KEY_IN_OUT_CUBIC:
			return CubicEaseInOut(t);
		case KEY_IN_QUART:
			return QuarticEaseIn(t);
		case KEY_OUT_QUART:
			return QuarticEaseOut(t);
		case KEY_IN_OUT_QUART:
			return QuarticEaseInOut(t);
		case KEY_IN_QUINT:
			return QuinticEaseIn(t);
		case KEY_OUT_QUINT:
			return QuinticEaseOut(t);
		case KEY_IN_OUT_QUINT:
			return QuinticEaseInOut(t);
		case KEY_IN_SINE:
			return SineEaseIn(t);
		case KEY_OUT_SINE:
			return SineEaseOut(t);
		case KEY_IN_OUT_SINE:
			return SineEaseInOut(t);
		case KEY_IN_CIRC:
			return CircularEaseIn(t);
		case KEY_OUT_CIRC:
			return CircularEaseOut(t);
		case KEY_IN_OUT_CIRC:
			return CircularEaseInOut(t);
		case KEY_IN_EXPO:
			return ExponentialEaseIn(t);
		case KEY_OUT_EXPO:
			return ExponentialEaseOut(t);
		case KEY_IN_OUT_EXPO:
			return ExponentialEaseInOut(t);
		case KEY_IN_ELASTIC:
			return ElasticEaseIn(t);
		case KEY_OUT_ELASTIC:
			return ElasticEaseOut(t);
		case KEY_IN_OUT_ELASTIC:
			return ElasticEaseInOut(t);
		case KEY_IN_BACK:
			return BackEaseIn(t);
		case KEY_OUT_BACK:
			return BackEaseOut(t);
		case KEY_IN_OUT_BACK:
			return BackEaseInOut(t);
		case KEY_IN_BOUNCE:
			return BounceEaseIn(t);
		case KEY_OUT_BOUNCE:
			return BounceEaseOut(t);
		case KEY_IN_OUT_BOUNCE:
			return BounceEaseInOut(t);

		default:
			assert(0);
			return t;
	}
}

double sync_get_val(const struct sync_track *t, double row) {
	int idx, irow;

	/* If we have no keys at all, return a constant 0 */
	if (!t->num_keys)
		return 0.0f;

	irow = (int)floor(row);
	idx = key_idx_floor(t, irow);

	/* at the edges, return the first/last value */
	if (idx < 0)
		return t->keys[0].value;
	if (idx > (int)t->num_keys - 2)
		return t->keys[t->num_keys - 1].value;

	double start = t->keys[idx].value;
	double end = t->keys[idx + 1].value;
	double row_a = t->keys[idx].row;
	double row_b = t->keys[idx + 1].row;
	double p = (row - row_a) / (row_b - row_a);

	return start + (end - start) * interpolate(p, t->keys[idx].type);
}

int sync_find_key(const struct sync_track *t, int row) {
	int lo = 0, hi = t->num_keys;

	/* binary search, t->keys is sorted by row */
	while (lo < hi) {
		int mi = (lo + hi) / 2;
		assert(mi != hi);

		if (t->keys[mi].row < row)
			lo = mi + 1;
		else if (t->keys[mi].row > row)
			hi = mi;
		else
			return mi; /* exact hit */
	}
	assert(lo == hi);

	/* return first key after row, negated and biased (to allow -0) */
	return -lo - 1;
}

#ifndef SYNC_PLAYER
int sync_set_key(struct sync_track *t, const struct track_key *k) {
	int idx = sync_find_key(t, k->row);
	if (idx < 0) {
		/* no exact hit, we need to allocate a new key */
		void *tmp;
		idx = -idx - 1;
		tmp = realloc(t->keys, sizeof(struct track_key) *
			(t->num_keys + 1));
		if (!tmp)
			return -1;
		t->num_keys++;
		t->keys = tmp;
		memmove(t->keys + idx + 1, t->keys + idx,
			sizeof(struct track_key) * (t->num_keys - idx - 1));
	}
	t->keys[idx] = *k;
	return 0;
}

int sync_del_key(struct sync_track *t, int pos) {
	void *tmp;
	int idx = sync_find_key(t, pos);
	assert(idx >= 0);
	memmove(t->keys + idx, t->keys + idx + 1,
		sizeof(struct track_key) * (t->num_keys - idx - 1));
	assert(t->keys);
	tmp = realloc(t->keys, sizeof(struct track_key) *
		(t->num_keys - 1));
	if (t->num_keys != 1 && !tmp)
		return -1;
	t->num_keys--;
	t->keys = tmp;
	return 0;
}
#endif
