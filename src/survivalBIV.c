
#include <math.h>
#include <R.h>
#include <string.h>

/*
Author:
	Artur Araújo (adapted from 'rcmp')

Description:
	Compares two double values.

Parameters:
	x[in]		first double
	y[in]		second double

Return value:
	Returns -1 if x is lower than y.
	Returns 1 if x is greater than y.
	Returns 0 otherwise.
*/

static int cmp_doubles(
	const double x,
	const double y)
{
	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
} // cmp_doubles

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'x' by increasing order with 'indx' and 'indy' alongside.

Parameters:
	x[inout]		pointer to vector's first element
	indx[inout]		pointer to vector's first element
	indy[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_univ_index(
	double *const x,
	int *const indx,
	int *const indy,
	const int n)
{
	double v;
	int i, j, h, iv, iu;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			v = x[i];
			iv = indx[i];
			iu = indy[i];
			j = i;
			while (j >= h && cmp_doubles(x[j - h], v) > 0) {
				x[j] = x[j-h];
				indx[j] = indx[j-h];
				indy[j] = indy[j-h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			indy[j] = iu;
		}
	}
} // sort_univ_index

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'x' by increasing order with 'indx', 'y' and 'indy' alongside.

Parameters:
	x[inout]		pointer to vector's first element
	indx[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	indy[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_biv(
	double *const x,
	int *const indx,
	double *const y,
	int *const indy,
	const int n)
{
	double v, u;
	int i, j, h, iv, iu;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			v = x[i];
			iv = indx[i];
			u = y[i];
			iu = indy[i];
			j = i;
			while (j >= h && cmp_doubles(x[j - h], v) > 0) {
				x[j] = x[j-h];
				indx[j] = indx[j-h];
				y[j] = y[j-h];
				indy[j] = indy[j-h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = u;
			indy[j] = iu;
		}
	}
} // sort_biv

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'x' by increasing order with 'indx', 'y' and 'z' alongside.

Parameters:
	x[inout]		pointer to vector's first element
	indx[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	z[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_biv_stime(
	double *const x,
	int *const indx,
	double *const y,
	double *const z,
	const int n)
{
	double v, u, iu;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			v = x[i];
			iv = indx[i];
			u = y[i];
			iu = z[i];
			j = i;
			while (j >= h && cmp_doubles(x[j - h], v) > 0) {
				x[j] = x[j-h];
				indx[j] = indx[j-h];
				y[j] = y[j-h];
				z[j] = z[j-h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = u;
			z[j] = iu;
		}
	}
} // sort_biv_stime

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'x' by increasing order with 'indx' and 'y' alongside.

Parameters:
	x[inout]		pointer to vector's first element
	indx[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_biv_time(
	double *const x,
	int *const indx,
	double *const y,
	const int n)
{
	double v, u;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			v = x[i];
			iv = indx[i];
			u = y[i];
			j = i;
			while (j >= h && cmp_doubles(x[j - h], v) > 0) {
				x[j] = x[j-h];
				indx[j] = indx[j-h];
				y[j] = y[j-h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = u;
		}
	}
} // sort_biv_time

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'x' by increasing order with 'indx', 'y' and 'z' alongside.

Parameters:
	x[inout]		pointer to vector's first element
	indx[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	z[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_biv_data(
	double *const x,
	double *const indx,
	double *const y,
	double *const z,
	const int n)
{
	double v, iv, u, iu;
	int i, j, h;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			v = x[i];
			iv = indx[i];
			u = y[i];
			iu = z[i];
			j = i;
			while (j >= h && cmp_doubles(x[j - h], v) > 0) {
				x[j] = x[j-h];
				indx[j] = indx[j-h];
				y[j] = y[j-h];
				z[j] = z[j-h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = u;
			z[j] = iu;
		}
	}
} // sort_biv_data

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'indx' by decreasing order with x and y alongside.

Parameters:
	indx[inout]		pointer to vector's first element
	x[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_double_rev(
	double *const indx,
	double *const x,
	double *const y,
	const int n)
{
	double iv, v, u;
	int i, j, h;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			u = y[i];
			j = i;
			while (j >= h && cmp_doubles(indx[j - h], iv) < 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				y[j] = y[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
			y[j] = u;
		}
	}
} // sort_double_rev

/*
Author:
	Artur Araújo (adapted from 'icmp')

Description:
	Compares two int values.

Parameters:
	x[in]		first int
	y[in]		second int

Return value:
	Returns -1 if x is lower than y.
	Returns 1 if x is greater than y.
	Returns 0 otherwise.
*/

static int cmp_ints(
	const int x,
	const int y)
{
	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
} // cmp_ints

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'indx' by decreasing order with 'x' alongside.

Parameters:
	indx[inout]		pointer to vector's first element
	x[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_int_rev_index(
	int *const indx,
	int *const x,
	const int n)
{
	int i, j, h, v, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			j = i;
			while (j >= h && cmp_ints(indx[j - h], iv) < 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
		}
	}
} // sort_int_rev_index

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'indx' by decreasing order with x and y alongside.

Parameters:
	indx[inout]		pointer to vector's first element
	x[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_int_rev(
	int *const indx,
	double *const x,
	double *const y,
	const int n)
{
	double v, u;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			u = y[i];
			j = i;
			while (j >= h && cmp_ints(indx[j - h], iv) < 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				y[j] = y[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
			y[j] = u;
		}
	}
} // sort_int_rev

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Restores the elements of a sorted double vector back
		to their previous order.

Parameters:
	indx[inout]		pointer to index vector's first element
	x[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_back_double(
	int *const indx,
	double *const x,
	const int n)
{
	double v;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			j = i;
			while (j >= h && cmp_ints(indx[j - h], iv) > 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
		}
	}
} // sort_back_double

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'indx' by increasing order with 'x', 'indy' and 'y' alongside.

Parameters:
	indx[inout]		pointer to vector's first element
	x[inout]		pointer to vector's first element
	indy[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_back_univ_weights(
	int *const indx,
	double *const x,
	int *const indy,
	double *const y,
	const int n)
{
	double v, u;
	int i, j, h, iv, iu;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			iu = indy[i];
			u = y[i];
			j = i;
			while (j >= h && cmp_ints(indx[j - h], iv) > 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				indy[j] = indy[j-h];
				y[j] = y[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
			indy[j] = iu;
			y[j] = u;
		}
	}
} // sort_back_univ_weights

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'indx' by increasing order with 'x', 'indy' and 'y' alongside.

Parameters:
	indx[inout]		pointer to vector's first element
	x[inout]		pointer to vector's first element
	indy[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_back_univ_weights_double(
	int *const indx,
	double *const x,
	double *const indy,
	double *const y,
	const int n)
{
	double v, u, iu;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			iu = indy[i];
			u = y[i];
			j = i;
			while (j >= h && cmp_ints(indx[j - h], iv) > 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				indy[j] = indy[j-h];
				y[j] = y[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
			indy[j] = iu;
			y[j] = u;
		}
	}
} // sort_back_univ_weights_double

/*
Author:
	Artur Araújo

Description:
	Sorts time1 by increasing order with delta alongside.
	In the subsets of constant time1 observations,
		events are put first and censures last.

Parameters:
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	index[inout]	pointer to index first element
	t1 [in]			time1 value (defines last index)
	len[in]			length of time1, delta and index

Return value:
	Returns last index.

Remarks:
	Vectors time1, delta and index must have the same length.
*/

static int sort_univ_surv_index(
	double *const time1,
	int *const delta,
	int *const index,
	const double t1,
	const int len)
{
	register int i = 0;
	int j, k, end = len/2;
	sort_univ_index(time1, delta, index, len); // sort time1 and delta with index
	if (time1[end] > t1) end = 0;
	for (; end < len; end++) {
		if (time1[end] > t1) break; // determine last index
	}
	while (i < end) { // loop through the sample until last index is reached
		for (i++, j = 1; i < end && time1[i] == time1[i-1]; i++) { // loop through the sample until time1 changes or last index is reached
			j++; // count equal times
		}
		if (j > 1) { // if there are equal times
			k = i-j;
			sort_int_rev_index(&delta[k], &index[k], j); // put censored observations last
		}
	}
	return end;
} // sort_univ_surv_index

/*
Author:
	Artur Araújo

Description:
	Sorts time1 by increasing order with delta alongside.
	In the subsets of constant time1 observations,
		events are put first and censures last.

Parameters:
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	index[inout]	pointer to index first element
	t1 [in]			time1 value (defines last index)
	len[in]			length of time1, delta and index

Return value:
	Returns last index.

Remarks:
	Vectors time1, delta and index must have the same length.
*/

static int sort_univ_surv_index_double(
	double *const time1,
	double *const delta,
	int *const index,
	const double t1,
	const int len)
{
	register int i = 0;
	int j, k, end = len/2;
	sort_biv_time(time1, index, delta, len); // sort time1 and delta with index
	if (time1[end] > t1) end = 0;
	for (; end < len; end++) {
		if (time1[end] > t1) break; // determine last index
	}
	while (i < end) { // loop through the sample until last index is reached
		for (i++, j = 1; i < end && time1[i] == time1[i-1]; i++) { // loop through the sample until time1 changes or last index is reached
			j++; // count equal times
		}
		if (j > 1) { // if there are equal times
			k = i-j;
			revsort(&delta[k], &index[k], j); // put censored observations last
		}
	}
	return end;
} // sort_univ_surv_index_double

/*
Author:
	Artur Araújo

Description:
	Sorts Stime by increasing order with status alongside.
	In the subsets of constant Stime observations,
		events are put first and censures last.

Parameters:
	Stime[inout]	pointer to Stime first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	time2[inout]	pointer to time2 first element
	len[in]			length of Stime, status, time1 and time2

Return value:
	This function doesn't return a value.

Remarks:
	Vectors Stime, status, time1 and time2 must have the same length.
*/

static void sort_biv_surv_stime(
	double *const Stime,
	int *const status,
	double *const time1,
	double *const time2,
	const int len)
{
	register int i = 0;
	int j, k;
	sort_biv_stime(Stime, status, time1, time2, len); // sort data
	while (i < len) { // loop through the sample until last index is reached
		for (i++, j = 1; i < len && Stime[i] == Stime[i-1]; i++) { // loop through the sample until Stime changes or last index is reached
			j++; // count equal Stimes
		}
		if (j > 1) { // if there are equal Stimes
			k = i-j;
			sort_int_rev(&status[k], &time1[k], &time2[k], j); // put censored observations last
		}
	}
	return;
} // sort_biv_surv_stime

/*
Author:
	Artur Araújo

Description:
	Sorts Stime by increasing order with m alongside.
	In the subsets of constant Stime observations,
		m observations are sorted by decreasing order.

Parameters:
	Stime[inout]	pointer to Stime first element
	m[inout]		pointer to status first element
	time1[inout]	pointer to time1 first element
	time2[inout]	pointer to time2 first element
	len[in]			length of Stime, m, time1 and time2

Return value:
	This function doesn't return a value.

Remarks:
	Vectors Stime, m, time1 and time2 must have the same length.
*/

static void sort_biv_surv(
	double *const Stime,
	double *const m,
	double *const time1,
	double *const time2,
	const int len)
{
	register int i = 0;
	int j, k;
	sort_biv_data(Stime, m, time1, time2, len); // sort data
	while (i < len) { // loop through the sample until last index is reached
		for (i++, j = 1; i < len && Stime[i] == Stime[i-1]; i++) { // loop through the sample until Stime changes or last index is reached
			j++; // count equal Stimes
		}
		if (j > 1) { // if there are equal Stimes
			k = i-j;
			sort_double_rev(&m[k], &time1[k], &time2[k], j); // put censored observations last
		}
	}
	return;
} // sort_biv_surv

/*
Author:
	Artur Araújo

Description:
	Computes a single probability value at a specified time index.

Parameters:
	time[in]		pointer to time first element
	status[in]		pointer to status first element
	len[in]			pointer to length of time and status
	end[in]			pointer to index to compute the probability at
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time and status must be sorted first.
	Vectors time and status must have the same length.
*/

void KaplanMeierValue(
	const double *const time,
	const int *const status,
	const int *const len,
	const int *const end,
	double *const surv)
{
	register int i = 0;
	int n, d;
	*surv = 1;
	while (i < *end) { // loop through the sample until last index is reached
		n = *len-i; // count the living
		d = status[i]; // initialize dead count
		for (i++; i < *end && time[i] == time[i-1]; i++) { // loop until time changes or last index is reached
			d += status[i]; // count the dead
		}
		*surv *= 1-(double)d/n; // compute survival probability
	}
	return;
} // KaplanMeierValue

/*
Author:
	Artur Araújo

Description:
	Sorts time and status and then calls 'KaplanMeierValue'.

Parameters:
	time[inout]		pointer to time first element
	status[inout]	pointer to status first element
	len[in]			pointer to length of time and status
	t[in]			pointer to time to compute the probability at
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time and status must have the same length.
*/

void KaplanMeierValueSort(
	double *const time,
	int *const status,
	const int *const len,
	const double *const t,
	double *const surv)
{
	int end = *len/2;
	rsort_with_index(time, status, *len); // use internal R sorting to sort time and status
	if (time[end] > *t) end = 0;
	for (; end < *len; end++) {
		if (time[end] > *t) break; // determine last index
	}
	KaplanMeierValue(time, status, len, &end, surv); // compute survival probability
	return;
} // KaplanMeierValueSort

/*
Author:
	Artur Araújo

Description:
	Computes probabilities indexed at the unique times up to a specified
		time index.

Parameters:
	time[in]		pointer to time first element
	status[in]		pointer to status first element
	len[in]			pointer to length of time and status
	end[in]			pointer to index to compute the probability at
	unique[out]		pointer to unique first element (vector of unique times)
	surv[out]		pointer to survival probabilities vector
	lenu[out]		pointer to length of unique and surv vectors

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time and status must be sorted first.
	Vectors time and status must have the same length.
*/

void KaplanMeierVector(
	const double *const time,
	const int *const status,
	const int *const len,
	const int *const end,
	double *const unique,
	double *const surv,
	int *const lenu)
{
	register int i = 0, j = 0;
	int n, d;
	double p = 1;
	while (i < *end) { // loop through the sample until last index is reached
		n = *len-i; // count the living
		d = status[i]; // initialize dead count
		for (i++; i < *end && time[i] == time[i-1]; i++) { // loop until time changes or last index is reached
			d += status[i]; // count the dead
		}
		p *= 1-(double)d/n; // compute survival probability
		unique[j] = time[i-1];
		surv[j++] = p;
	}
	*lenu = j;
	return;
} // KaplanMeierVector

/*
Author:
	Artur Araújo

Description:
	Sorts time and status and then calls 'KaplanMeierVector'.

Parameters:
	time[inout]		pointer to time first element
	status[inout]	pointer to status first element
	len[in]			pointer to length of time and status
	t[in]			pointer to time to compute the probability at (defines last index)
	unique[out]		pointer to unique first element (vector of unique times)
	surv[out]		pointer to survival probabilities vector
	lenu[out]		pointer to length of unique and surv vectors

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time and status must have the same length.
*/

void KaplanMeierVectorSort(
	double *const time,
	int *const status,
	const int *const len,
	const double *const t,
	double *const unique,
	double *const surv,
	int *const lenu)
{
	int end = *len/2;
	rsort_with_index(time, status, *len); // use internal R sorting to sort time and status
	if (time[end] > *t) end = 0;
	for (; end < *len; end++) {
		if (time[end] > *t) break; // determine last index
	}
	KaplanMeierVector(time, status, len, &end, unique, surv, lenu); // compute survival probabilities vector
	return;
} // KaplanMeierVectorSort

/*
Author:
	Artur Araújo

Description:
	Computes a single probability value at a specified time index.
	This function implements the weights version of the Kaplan-Meier
		estimator.

Parameters:
	time2[in]		pointer to time2 first element
	status[in]		pointer to status first element
	weights[in]		pointer to weights first element
	delta[in]		pointer to delta first element
	len[in]			pointer to length of time2, status, weights, and delta
	end[in]			pointer to index to compute the probability at
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, weights and delta must be sorted first.
	Vectors time2, status, weights and delta must have the same length.
*/

void WeightedKaplanMeierValue(
	const double *const time2,
	const int *const status,
	const double *const weights,
	const int *const delta,
	const int *const len,
	const int *const end,
	double *const surv)
{
	register int i;
	double n, d;
	*surv = 1;
	for (i = *len-1, n = 0; i >= *end; i--) { // loop in reverse order until end index is reached
		n += delta[i]*weights[i]; // initialize living weighting
	}
	while (i >= 0) { // loop through the sample in reverse order until zero index is reached
		n += delta[i]*weights[i];
		d = status[i]*weights[i]; // initialize dead weighting
		for (i--; i >= 0 && time2[i] == time2[i+1]; i--) { // loop in reverse order until time changes or zero index is reached
			n += delta[i]*weights[i]; // weight the living
			d += status[i]*weights[i]; // weight the dead
		}
		if (n > 0) *surv *= 1-d/n; // compute survival probability
	}
	return;
} // WeightedKaplanMeierValue

/*
Author:
	Artur Araújo

Description:
	Sorts time2, status, weights and delta and then calls 'WeightedKaplanMeierValue'.

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	weights[inout]	pointer to weights first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, weights, and delta
	t[in]			pointer to time to compute the probability at (defines last index)
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, weights and delta must have the same length.
*/

void WeightedKaplanMeierValueSort(
	double *const time2,
	int *const status,
	double *const weights,
	int *const delta,
	const int *const len,
	const double *const t,
	double *const surv)
{
	int end = *len/2;
	sort_biv(time2, status, weights, delta, *len); // sort data
	if (time2[end] > *t) end = 0;
	for (; end < *len; end++) {
		if (time2[end] > *t) break; // determine last index
	}
	WeightedKaplanMeierValue(time2, status, weights, delta, len, &end, surv); // compute survival probability
	return;
} // WeightedKaplanMeierValueSort

/*
Author:
	Artur Araújo

Description:
	Computes probabilities indexed at the unique times up to a specified
		time index.
	This function implements the weights version of the Kaplan-Meier
		estimator.

Parameters:
	time2[in]		pointer to time2 first element
	status[in]		pointer to status first element
	weights[in]		pointer to weights first element
	delta[in]		pointer to delta first element
	len[in]			pointer to length of time2, status, weights, and delta
	end[in]			pointer to index to compute the probability at
	unique[out]		pointer to unique first element (vector of unique times)
	surv[out]		pointer to survival probabilities vector
	istart[out]		pointer to start index of unique and surv vectors

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, weights and delta must be sorted first.
	Vectors time2, status, weights and delta must have the same length.
	istart is the start index in R's notation (one more than C's notation)
*/

void WeightedKaplanMeierVector(
	const double *const time2,
	const int *const status,
	const double *const weights,
	const int *const delta,
	const int *const len,
	const int *const end,
	double *const unique,
	double *const surv,
	int *const istart)
{
	register int i, j = *len-1;
	double n, d;
	for (i = *len-1, n = 0; i >= *end; i--) { // loop in reverse order until end index is reached
		n += delta[i]*weights[i]; // initialize living weighting
	}
	while (i >= 0) { // loop through the sample in reverse order until zero index is reached
		n += delta[i] * weights[i];
		d = status[i] * weights[i]; // initialize dead weighting
		for (i--; i >= 0 && time2[i] == time2[i+1]; i--) { // loop in reverse order until time changes or zero index is reached
			n += delta[i] * weights[i]; // weight the living
			d += status[i] * weights[i]; // weight the dead
		}
		unique[j] = time2[i+1];
		if (n > 0) {
			surv[j--] = 1-d/n;
		} else {
			surv[j--] = 1;
		}
	}
	*istart = j+2; // first index of output vector in R's notation
	for (i = *istart; i < *len; i++) {
		surv[i] *= surv[i-1]; // compute survival probabilities
	}
	return;
} // WeightedKaplanMeierVector

/*
Author:
	Artur Araújo

Description:
	Sorts time2, status, weights and delta and then calls 'WeightedKaplanMeierVector'.

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	weights[inout]	pointer to weights first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, weights, and delta
	t[in]			pointer to time to compute the probability at (defines last index)
	unique[out]		pointer to unique first element (vector of unique times)
	surv[out]		pointer to survival probabilities vector
	istart[out]		pointer to start index of unique and surv vectors

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, weights and delta must have the same length.
	istart is the start index in R's notation (one more than C's notation)
*/

void WeightedKaplanMeierVectorSort(
	double *const time2,
	int *const status,
	double *const weights,
	int *const delta,
	const int *const len,
	const double *const t,
	double *const unique,
	double *const surv,
	int *const istart)
{
	int end = *len/2;
	sort_biv(time2, status, weights, delta, *len); // sort data
	if (time2[end] > *t) end = 0;
	for (; end < *len; end++) {
		if (time2[end] > *t) break; // determine last index
	}
	WeightedKaplanMeierVector(time2, status, weights, delta, len, &end, unique, surv, istart); // compute survival probabilities vector
	return;
} // WeightedKaplanMeierVectorSort

/*
Author:
	Artur Araújo

Description:
	Computes the Kaplan-Meier weights of a given delta vector.

Parameters:
	delta[in]		pointer to delta first element
	len[in]			pointer to length of delta
	end[in]			pointer to index limit
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.

Remarks:
	Vector weights has the same length as delta.
*/

void WeightsKaplanMeier(
	const int *const delta,
	const int *const len,
	const int *const end,
	double *const weights)
{
	register int i;
	double aux[2]; // declare auxiliary vector needed for the computation
	for (i = 0, aux[0] = 1; i < *end; i++) { // loop through the sample until last index is reached
		weights[i] = (double)delta[i]/(*len-i); // compute needed factor
		aux[1] = 1-weights[i]; // factor needed for the computation
		weights[i] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	return;
} // WeightsKaplanMeier

/*
Author:
	Artur Araújo

Description:
	Calls 'sort_univ_surv_index' to sort time1 and delta
		and then calls 'WeightsKaplanMeier'.

Parameters:
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time1 and delta
	t1[in]			pointer to time1 value limit
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time1 and delta must have the same length.
	Vector weights has the same length as time1 and delta.
*/

void WeightsKaplanMeierSort(
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t1,
	double *const weights)
{
	int i, index[*len];
	for (i = 0; i < *len; i++) index[i] = i;
	i = sort_univ_surv_index(time1, delta, index, *t1, *len); // sort time1 by increasing order with events first and censures last
	WeightsKaplanMeier(delta, len, &i, weights); // compute Kaplan-Meier weights
	sort_back_univ_weights(index, time1, delta, weights, *len); // put time1, delta and weights back to their previous order
	return;
} // WeightsKaplanMeierSort

/*
Author:
	Artur Araújo

Description:
	Computes the Kaplan-Meier weights of a given delta vector.

Parameters:
	delta[in]		pointer to delta first element
	len[in]			pointer to length of delta
	end[in]			pointer to index limit
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.

Remarks:
	Vector weights has the same length as delta.
*/

void WeightsKaplanMeierEx(
	const double *const delta,
	const int *const len,
	const int *const end,
	double *const weights)
{
	register int i;
	double aux[2]; // declare auxiliary vector needed for the computation
	for (i = 0, aux[0] = 1; i < *end; i++) { // loop through the sample until last index is reached
		weights[i] = delta[i]/(*len-i); // compute needed factor
		aux[1] = 1-weights[i]; // factor needed for the computation
		weights[i] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	return;
} // WeightsKaplanMeierEx

/*
Author:
	Artur Araújo

Description:
	Calls 'sort_univ_surv_index_double' to sort time1 and delta
		and then calls 'WeightsKaplanMeierEx'.

Parameters:
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time1 and delta
	t1[in]			pointer to time1 value limit
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time1 and delta must have the same length.
	Vector weights has the same length as time1 and delta.
*/

void WeightsKaplanMeierSortEx(
	double *const time1,
	double *const delta,
	const int *const len,
	const double *const t1,
	double *const weights)
{
	int i, index[*len];
	for (i = 0; i < *len; i++) index[i] = i;
	i = sort_univ_surv_index_double(time1, delta, index, *t1, *len); // sort time1 by increasing order with events first and censures last
	WeightsKaplanMeierEx(delta, len, &i, weights); // compute Kaplan-Meier weights
	sort_back_univ_weights_double(index, time1, delta, weights, *len); // put time1, delta and weights back to their previous order
	return;
} // WeightsKaplanMeierSortEx

/*
Author:
	Artur Araújo

Description:
	Computes nearest neighbour weights.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	span[in]		pointer to span
	window[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void WeightsNNE(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const span,
	const char *const *const window,
	double *const weights)
{
	char *window1 = "asymmetric";
	char *window2 = "symmetric";
	register int i;
	int index[*len], tr;
	double lambda1;
	for (i = 0; i < *len; i++) {
		weights[i] = time1[i]-*t1; // initialize weights vector
		index[i] = i; // initialize index vector
	}
	rsort_with_index(weights, index, *len); // use internal R sorting to sort weights
	if (strcmp(*window, window1) == 0) { // if window is asymmetric
		tr = *len * *span;
		for (i = *len-1; i > 0; i--) {
			if (weights[i] < 0) break; // compute index
		}
		if (i+tr+1 > *len-1) lambda1 = weights[*len-1];
		else lambda1 = weights[i+tr+1];
		for (i = 0; i < *len; i++) {
			weights[i] = (weights[i] >= 0 && weights[i] <= lambda1); // compute weights
		}
	}
	else if (strcmp(*window, window2) == 0) { // if window is symmetric
		double lambda0;
		tr = *len * *span / 2;
		for (i = *len-1; i > 0; i--) {
			if (weights[i] <= 0) break; // compute lower index
		}
		if (i-tr < 0) lambda0 = -fabs(weights[0]);
		else lambda0 = -fabs(weights[i-tr]);
		for (; i > 0; i--) {
			if (weights[i] < 0) break; // compute upper index
		}
		if (i+tr+1 > *len-1) lambda1 = weights[*len-1];
		else lambda1 = weights[i+tr+1];
		for (i = 0; i < *len; i++) {
			weights[i] = ( (weights[i] >= lambda0 && weights[i] <= 0) || (weights[i] >= 0 && weights[i] <= lambda1) ); // compute weights
		}
	}
	sort_back_double(index, weights, *len); // put weights back to their previous order
	return;
} // WeightsNNE

/*
Author:
	Artur Araújo

Description:
	Computes nearest neighbour Nadaraya-Watson weights.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	span[in]		pointer to span
	window[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void NWWeightsNNE(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const span,
	const char *const *const window,
	double *const weights)
{
	register int i;
	double sum;
	WeightsNNE(time1, len, t1, span, window, weights); // compute nearest neighbour density
	for (i = 0, sum = 0; i < *len; i++) {
		sum += weights[i]; // compute sum
	}
	for (i = 0; i < *len; i++) {
		weights[i] /= sum; // compute Nadaraya-Watson weight
	}
	return;
} // NWWeightsNNE

/*
Author:
	Artur Araújo

Description:
	Computes local linear nearest neighbour weights.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	span[in]		pointer to span
	window[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void LLWeightsNNE(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const span,
	const char *const *const window,
	double *const weights)
{
	register int i;
	double aux[2], sum[2] = {0, 0};
	NWWeightsNNE(time1, len, t1, span, window, weights); // compute nearest neighbour density
	for (i = 0; i < *len; i++) {
		aux[0] = time1[i]-*t1;
		aux[1] = aux[0]*weights[i];
		sum[0] += aux[1];
		sum[1] += aux[0]*aux[1];
	}
	for (i = 0; i < *len; i++) {
		weights[i] *= sum[1]-(time1[i]-*t1)*sum[0];
	}
	for (i = 0, sum[0] = 0; i < *len; i++) {
		sum[0] += weights[i]; // compute sum
	}
	for (i = 0; i < *len; i++) {
		weights[i] /= sum[0]; // compute local linear weight
	}
	return;
} // LLWeightsNNE

/*
Author:
	Artur Araújo

Description:
	Computes weights based on a kernel.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void WeightsKernel(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const weights)
{
	char *kernel1 = "gaussian";
	char *kernel2 = "epanechnikov";
	char *kernel3 = "tricube"; // also known as triweight
	char *kernel4 = "boxcar"; // also known as uniform
	char *kernel5 = "triangular";
	char *kernel6 = "quartic"; // also known as biweight
	char *kernel7 = "cosine";
	register int i;
	for (i = 0; i < *len; i++) {
		weights[i] = (time1[i]-*t1) / *lambda;
	}
	if (strcmp(*kernel, kernel1) == 0) { // if kernel is gaussian
		for (i = 0; i < *len; i++) {
			weights[i] = exp(-pow(weights[i], 2)/2);
		}
	}
	else if (strcmp(*kernel, kernel2) == 0) { // if kernel is epanechnikov
		for (i = 0; i < *len; i++) {
			weights[i] = ( 1-pow(weights[i], 2) )*(fabs(weights[i]) <= 1);
		}
	}
	else if (strcmp(*kernel, kernel3) == 0) { // if kernel is tricube
		for (i = 0; i < *len; i++) {
			weights[i] = fabs(weights[i]);
			weights[i] = pow(1-pow(weights[i], 3), 3)*(weights[i] <= 1);
		}
	}
	else if (strcmp(*kernel, kernel4) == 0) { // if kernel is boxcar
		for (i = 0; i < *len; i++) {
			weights[i] = (fabs(weights[i]) <= 1);
		}
	}
	else if (strcmp(*kernel, kernel5) == 0) { // if kernel is triangular
		for (i = 0; i < *len; i++) {
			weights[i] = fabs(weights[i]);
			weights[i] = (1-weights[i])*(weights[i] <= 1);
		}
	}
	else if (strcmp(*kernel, kernel6) == 0) { // if kernel is quartic
		for (i = 0; i < *len; i++) {
			weights[i] = pow(1-pow(weights[i], 2), 2)*(fabs(weights[i]) <= 1);
		}
	}
	else if (strcmp(*kernel, kernel7) == 0) { // if kernel is cosine
		for (i = 0; i < *len; i++) {
			weights[i] = cos(M_PI*weights[i]/2)*(fabs(weights[i]) <= 1);
		}
	}
	return;
} // WeightsKernel

/*
Author:
	Artur Araújo

Description:
	Computes the Nadaraya-Watson weights.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void NWWeightsKernel(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const weights)
{
	register int i;
	double sum;
	WeightsKernel(time1, len, t1, lambda, kernel, weights); // compute kernel density
	for (i = 0, sum = 0; i < *len; i++) {
		sum += weights[i]; // compute sum
	}
	for (i = 0; i < *len; i++) {
		weights[i] /= sum; // compute Nadaraya-Watson weight
	}
	return;
} // NWWeightsKernel

/*
Author:
	Artur Araújo

Description:
	Computes local linear weights based on a kernel.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void LLWeightsKernel(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const weights)
{
	register int i;
	double aux[2], sum[2] = {0, 0};
	WeightsKernel(time1, len, t1, lambda, kernel, weights); // compute kernel density
	for (i = 0; i < *len; i++) {
		aux[0] = time1[i]-*t1;
		aux[1] = aux[0]*weights[i];
		sum[0] += aux[1];
		sum[1] += aux[0]*aux[1];
	}
	for (i = 0; i < *len; i++) {
		weights[i] *= sum[1]-(time1[i]-*t1)*sum[0];
	}
	for (i = 0, sum[0] = 0; i < *len; i++) {
		sum[0] += weights[i]; // compute sum
	}
	for (i = 0; i < *len; i++) {
		weights[i] /= sum[0]; // compute local linear weight
	}
	return;
} // LLWeightsKernel

/*
Author:
	Artur Araújo

Description:
	Computes the conditional survival probability P(T2>t2|T1=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, time1 and delta
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	span[in]		pointer to span
	window[in]		pointer to pointer that points to a char array
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and delta must have the same length.
*/

void SurvBeranNNE(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t2,
	const double *const t1,
	const double *const span,
	const char *const *const window,
	double *const surv)
{
	double weights[*len];
	WeightsNNE(time1, len, t1, span, window, weights); // compute NNE weights
	WeightedKaplanMeierValueSort(time2, status, weights, delta, len, t2, surv); // compute conditional survival probability
	return;
} // SurvBeranNNE

/*
Author:
	Artur Araújo

Description:
	Computes the conditional survival probability P(T2>t2|T1=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, time1 and delta
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and delta must have the same length.
*/

void SurvBeranKernel(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t2,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const surv)
{
	double weights[*len];
	WeightsKernel(time1, len, t1, lambda, kernel, weights); // compute kernel weights
	WeightedKaplanMeierValueSort(time2, status, weights, delta, len, t2, surv); // compute conditional survival probability
	return;
} // SurvBeranKernel

/*
Author:
	Artur Araújo

Description:
	Computes the bivariate probability P(T2<=t2,T1<=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, time1 and delta
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	span[in]		pointer to span
	window[in]		pointer to pointer that points to a char array
	p[out]			pointer to probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and delta must have the same length.
*/

void BivDistBeranNNE(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t2,
	const double *const t1,
	const double *const span,
	const char *const *const window,
	double *const p)
{
	int end = *len/2;
	register int i;
	double W[*len], weights[*len], surv;
	sort_biv(time2, status, time1, delta, *len); // sort data
	if (time2[end] > *t2) end = 0;
	for (; end < *len; end++) {
		if (time2[end] > *t2) break; // determine last index
	}
	WeightsKaplanMeierSort(time1, delta, len, t1, W); // compute Kaplan-Meier weights
	for (i = 0, *p = 0; i < *len; i++) {
		if (time1[i] <= *t1) {
			WeightsNNE(time1, len, &time1[i], span, window, weights); // compute NNE weights
			WeightedKaplanMeierValue(time2, status, weights, delta, len, &end, &surv); // compute conditional survival probability
			*p += (1-surv)*W[i]; // compute bivariate probability
		}
	}
	return;
} // BivDistBeranNNE

/*
Author:
	Artur Araújo

Description:
	Computes the bivariate probability P(T2<=t2,T1<=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, time1 and delta
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	p[out]			pointer to probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and delta must have the same length.
*/

void BivDistBeranKernel(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t2,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const p)
{
	int end = *len/2;
	register int i;
	double W[*len], weights[*len], surv;
	sort_biv(time2, status, time1, delta, *len); // sort data
	if (time2[end] > *t2) end = 0;
	for (; end < *len; end++) {
		if (time2[end] > *t2) break; // determine last index
	}
	WeightsKaplanMeierSort(time1, delta, len, t1, W); // compute Kaplan-Meier weights
	for (i = 0, *p = 0; i < *len; i++) {
		if (time1[i] <= *t1) {
			WeightsKernel(time1, len, &time1[i], lambda, kernel, weights); // compute kernel weights
			WeightedKaplanMeierValue(time2, status, weights, delta, len, &end, &surv); // compute conditional survival probability
			*p += (1-surv)*W[i]; // compute bivariate probability
		}
	}
	return;
} // BivDistBeranKernel

/*
Author:
	Artur Araújo

Description:
	Computes the bivariate probability P(T2<=t2,T1<=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, time1 and delta
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	p[out]			pointer to probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and delta must have the same length.
*/

void BivDistCKM(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t2,
	const double *const t1,
	double *const p)
{
	register int i;
	int end = *len/2, j;
	double surv = 1;
	sort_biv(time1, delta, time2, status, *len); // sort data
	if (time1[end] > *t1) end = 0;
	for (; end < *len; end++) {
		if (time1[end] > *t1) break; // determine last index
	}
	KaplanMeierValue(time1, delta, len, &end, &surv); // compute survival probability
	for (i = 0, j = 0; i < end; i++) { // subset data
		if (delta[i]) {
			time2[j] = time2[i];
			status[j++] = status[i];
		}
	}
	KaplanMeierValueSort(time2, status, &j, t2, p); // compute conditional survival probability
	*p = (1-surv)*(1-*p); // compute bivariate probability
	return;
} // BivDistCKM

/*
Author:
	Artur Araújo

Description:
	Computes the bivariate probability P(T2<=t2,T1<=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	Stime[inout]	pointer to Stime first element
	len[in]			pointer to length of time2, status, time1 and Stime
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	p[out]			pointer to probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and Stime must have the same length.
*/

void BivDistKMW(
	double *const time2,
	int *const status,
	double *const time1,
	double *const Stime,
	const int *const len,
	const double *const t2,
	const double *const t1,
	double *const p)
{
	register int i;
	double aux[3];
	sort_biv_surv_stime(Stime, status, time1, time2, *len); // sort data
	for (i = 0, aux[0] = 1, *p = 0; i < *len; i++) { // loop through the sample until last index is reached
		aux[2] = (double)status[i]/(*len-i); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		*p += aux[2]*(time1[i] <= *t1 && time2[i] <= *t2); // compute bivariate probability
	}
	return;
} // BivDistKMW

/*
Author:
	Artur Araújo

Description:
	Computes the bivariate probability P(T2<=t2,T1<=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	m[inout]		pointer to m first element
	time1[inout]	pointer to time1 first element
	Stime[inout]	pointer to Stime first element
	len[in]			pointer to length of time2, m, time1 and Stime
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	p[out]			pointer to probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, m, time1 and Stime must have the same length.
*/

void BivDistKMPW(
	double *const time2,
	double *const m,
	double *const time1,
	double *const Stime,
	const int *const len,
	const double *const t2,
	const double *const t1,
	double *const p)
{
	register int i;
	double aux[3];
	sort_biv_surv(Stime, m, time1, time2, *len); // sort data
	for (i = 0, aux[0] = 1, *p = 0; i < *len; i++) { // loop through the sample until last index is reached
		aux[2] = m[i]/(*len-i); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		*p += aux[2]*(time1[i] <= *t1 && time2[i] <= *t2); // compute bivariate probability
	}
	return;
} // BivDistKMPW

/*
Author:
	Artur Araújo

Description:
	Computes the bivariate probability P(T2<=t2,T1<=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	Stime[inout]	pointer to Stime first element
	len[in]			pointer to length of time2, status, time1, delta and Stime
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	p[out]			pointer to probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1, delta and Stime must have the same length.
*/

void BivDistIPCW(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	double *const Stime,
	const int *const len,
	const double *const t2,
	const double *const t1,
	double *const p)
{
	register int i, j, k;
	int n, d;
	double surv, survS;
	sort_biv_time(time1, delta, time2, *len); // sort time1, delta and time2
	rsort_with_index(Stime, status, *len); // sort Stime and status
	for (i = 0, j = 0, k = 0, surv = 1, survS = 1, *p = 0; i < *len && time1[i] <= *t1; i++) { // loop through the sample
		if (j < *len && time1[j] == time1[i]) {
			n = *len-j; // count the living
			d = 1-delta[j]; // initialize dead count
			for (j++; j < *len && time1[j] == time1[j-1]; j++) { // loop until time changes or last index is reached
				d += 1-delta[j]; // count the dead
			}
			surv *= 1-(double)d/n; // compute survival probability
		}
		if (surv > 0) *p += (time2[i] > 0)/surv; // add term
		while (k < *len && Stime[k] <= time1[i]+*t2) {
			n = *len-k; // count the living
			d = 1-status[k]; // initialize dead count
			for (k++; k < *len && Stime[k] == Stime[k-1]; k++) { // loop until time changes or last index is reached
				d += 1-status[k]; // count the dead
			}
			survS *= 1-(double)d/n; // compute survival probability
		}
		if (survS > 0) *p -= (time2[i] > *t2)/survS; // subtract term
	}
	*p /= *len; // compute bivariate probability
	return;
} // BivDistIPCW
