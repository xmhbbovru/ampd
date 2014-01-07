// ampd-arm.cc rev. 06 Jan 2014 by Stuart Ambler.  Implements automatic
// multiscale-based peak detection (AMPD) algorithm as in An Efficient
// Algorithm for Automatic Peak Detection in Noisy Periodic and Quasi-Periodic
// Signals, by Felix Scholkmann, Jens Boss and Martin Wolf, Algorithms 2012, 5,
// 588-603.  Self-contained implementation and main function for test, except
// for dependency on argtable2.h for parsing command line arguments.
// Revised from ampd.cc for armadillo 4.000.0.
// Copyright (c) 2014 Stuart Ambler.
// Distributed under the Boost License in the accompanying file LICENSE.

// The algorithm, restated using 0-based array and matrix subscripts, and
// changed to add computation of minima as well as maxima.  In the process of
// doing this I came across something unexplained in the paper: it implicitly
// wants the number of windows evaluated for highest k to be >= 2.  It gives no
// reason why k >= 1 isn't good enough; perhaps it's so as to have one formula
// regardless of whether N is even or odd.  To see this, the number of windows
// at highest k, in the paper's notation, is N-khi+1 -(khi+2) + 1 = N-2*khi.
// Since the paper gives khi=L=ceil(N/2)-1, the number of windows is
// N-2*(ceil(N/2)-1) = N-2*ceil(N/2)+2 = 2 for even N, 3 for odd N.  If khi
// were increased by 1 regardless of the parity of N, the number of windows
// would be 0 for even N, which is no good, or 1 for odd N.

// Let x = {x0, x1, ... x(n-1)} be the sample of length n.  First calculate
// the least-squares straight line fit to x and subtract it.  Then,
// let el = ceil(n/2) - 1,, and for k=1,...,el, kix=k-1, do the following:
//   Not for calculation, let wk = 2*(k+1) be the window width, or more
//   precisely, comparison length to check if a given value of x, xi can be a
//   maximum or minimum vs. its neighbors, x(i+-k).  Consider the elements of x
//   too near its start or end to be able to check both such neighbors within x,
//   are not maxima or minima.  Construct two el x n matrices mpk, mtr, of
//   doubles with  alpha + random uniform number in [0, 1] (evaluated for each
//   element of mpk, mtr that's given this value) being assigned to the elements
//   of mpk that aren't (strict) maxima and of mtr that aren't (strict) minima,
//   and 0 assigned to those of mpk that are maxima, and of mtr that are minima,
//   as follows:
//     for i = 0,...,k-1=kix define not max/min since can't do both comparisons
//         i = n-k=n-kix-1,...,n-1 the same
//         i = k=kix+1,...,n-k-1=n-kix-2
//           xi > x(i-k)=x(i-kix-1) and xi > x(i+k)=x(i+kix+1) => max
//                                                             => mpk(kix,i) = 0
//           xi < x(i-k)=x(i-kix-1) and xi < x(i+k)=x(i+kix+1) => min
//                                                             => mtr(kix,i) = 0

// Then the remainder of the algorithm is applied to mpk and mtr separately for
// max and min, still adjusting for 0-base.  We refer to mpk/mtr as m here.
// Calculate gamma = {gamma0, gamma1, ..., gamma(el-1)} as
// gamma(kix) = sum of row k of m.  Define lamb (lambda) = the first kix for
// which gamma(kix) is a minimum among all the values of gamma, and remove rows
// kix+1,... from m, resulting in a (lamb + 1) x n matrix mr (or, simply don't
// use those rows).  Note that the method assumes that lamb > 0 in the form
// in the paper that computes something like the sample standard deviation (but
// no sqrt for 1/(lamb-1).  Since the more maxima/minima there are, the more
// zero entries there are in a row, the other entries being positive with
// expected value 1.5, the expected value of the sum of the row is
// 1.5 * (number of non-maxima/minima in the row).  Hence one would expect lamb
// to be the index of the row with the most maxima/minima in the row.
// Then for each column i of mr, calculate sumsqdev = sum over that column of
// the square of the value minus the mean of that column.  Those indices i for
// which sumsqdev is zero are returned in a vector of peaks/troughs, according
// to the algorithm in the paper (see below for optional modification).  (The
// algorithm actually takes the square root of sumsqdev and divides by lamb - 1,
// but I think this is a mistake though it may have been useful to the authors
// for their plots of sigma_i.  For the algorithm, it's unnecessary to take the
// square root, and dividing by lamb - 1 will blow up if lamb == 1.)

// I haven't figured out why the value of lamb is "correct", but it makes some
// sense intuitively.  If the signal is exactly periodic and monotone between
// peak and trough, peaks will be detected for scale values with windows up to
// the period.  Even if not "periodically monotone" there will be some region
// near the peak which falls away on either side, unless the sampling only gives
// one sample there.  The method fails to allow for lamb == 0, but this can be
// taken care of as a special case, just using the values in the first row of m.

// lamb is not always correct for all data.  For example it skips 5 of the 29
// peaks in data generated in the following (output in subdirectory ampd2)
// ../ampd -a 1. -b 1. -c .5 -d .1 -f 10. -g 70. -h 5. -i 5. -q 12. -s 0. -l 5. -n 100 -t 0. -u 0.1 -v 0.5 -w 0.083333 -z
// which calculates lamb == 1.  lamb == 0 works correctly for this data.
// I think the algorithm might work better for data with finer sampling of
// something smooth.

// I wasn't sure how close to zero to check sumsqdev; tried 8 * epsilon() and
// then 1.0e-16, both of which worked.  The paper says "equal", and I think
// that's what it should be.  In fact, I don't see why not just to check that
// all elements of the column (up through lamb) are zero.  If no randomness were
// used in m, then a column with *no* maxima/ minima at any scale would have
// sumsqdev zero and be chosen as having a peak/trough; not good.  Anyway it's
// simpler just to check all elements of the column zero, and so I added the
// -z command line option; I think it should always be set, and except that
// I was trying to implement the algorithm as it was, I would have made it
// the default and eliminated the option of setting it differently.

// Just to see if they mattered (so far I haven't found that they do), I added
// options to use normally distributed rather than uniformly distributed
// random numbers in the algorithm, and to change the mean and standard
// deviation of these distributions.

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <armadillo>
#include <argtable2.h>

using namespace std;
using namespace arma;

typedef pair<bool, double> UseIndexAndDevReturn;
typedef tuple<int, int, int, Col<int>, Col<int>> AmpdReturn;

// Print program usage.

void usage() {
  const char *s = "Test AMPD.\n"
"\n"
"Usage:\n"
" ./ampd [-a <double>] [-b <double>] [-c <double>] [-d <double>]\n"
"        [-f <double>] [-g <double>] [-h <double>] [-i <double>]\n"
"        [-q <double>] [-s <double>] [-l <double>] [-n <int>]\n"
"        [-t <double>] [-u <double>] [-v <double>] [-w <double>] [-o] [-z]\n"
"Defaults:\n"
" ./ampd -a 1. -b 1. -c .5 -d .1 -f 10. -g 70. -h 5. -i 5. -q 12. -s 0. -l 5.\n"
"        -n 1000 -t 0. -u 1.0 -v .5 -w 0.08333333\n"
"\n"
"Arguments (all optional):\n"
" -a --a         coefficient of freq. f1/fs term     [default  1.0]\n"
" -b --b         coefficient of freq. f2/fs term     [default  1.0]\n"
" -c --c         coefficient of freq. f3/fs term     [default  0.5]\n"
" -d --d         coefficient of random error term    [default  0.1]\n"
" -f --f1        frequency 1 (will be divided by fs) [default 10.0]\n"
" -g --f2        frequency 2 (will be divided by fs) [default 70.0]\n"
" -h --f3_start  frequency 3 starting value (/fs)    [default  5.0]\n"
" -i --f3_end    frequency 3 ending value   (/fs)    [default  5.0]\n"
" -q --fs        frequency divisor                   [default 12.0]\n"
" -s --start_t   starting time in 'seconds'          [default  0.0]\n"
" -l --len_t     time length   in 'seconds'          [default  5.0]\n"
" -n --n         number of samples of time series - twice its square times\n"
"                sizeof(double) must fit in memory   [default 1000]\n"
" -t --err_mean  mean of normal random error         [default  0.0]\n"
" -u --err_stdev standard dev of normal random error [default  1.0]\n"
" -v --alg_mean  mean of random nrs used in alg      [default  0.5]\n"
" -w --alg_stdev standard deviation rand nrs in alg  [default  1/12]\n"
"\n"
"Options:\n"
" -h --help   Show this help message and exit.\n"
" -o --normal Use normal rather than uniform dist for rand nrs in alg.\n"
" -z --zero   Test for column zero rather than zero variance.\n"
"\n";
  std::cout << s;
} 

// Get and process command-line arguments.  See usage().

void get_cfg(int argc, char *argv[],
             double *p_a,         double *p_b,
             double *p_c,         double *p_d,
             double *p_f1,        double *p_f2,
             double *p_f3_start,  double *p_f3_end, double *p_fs,
             double *p_start_t,   double *p_len_t,  int *p_n,
             double *p_err_mean,  double *p_err_stdev, 
             double *p_alg_mean,  double *p_alg_stdev, 
             bool   *p_help_only, bool   *p_normal,
             bool   *p_col_zero,  bool   *p_error) {
  struct arg_dbl *a =         arg_dbl0("a", "a", NULL, // "<a>",
                                       "coefficient of freq. f1/fs term");
  struct arg_dbl *b =         arg_dbl0("b", "b", NULL, // "<b>",
                                       "coefficient of freq. f2/fs term");
  struct arg_dbl *c =         arg_dbl0("c", "c", NULL, // "<c>",
                                       "coefficient of freq. f3/fs term");
  struct arg_dbl *d =         arg_dbl0("d", "d", NULL, // "<d>",
                                       "coefficient of random error term");
  struct arg_dbl *f1 =        arg_dbl0("f", "f1", NULL, // "<f1>",
                                       "frequency 1 (will be divided by fs");
  struct arg_dbl *f2 =        arg_dbl0("g", "f2", NULL, // "<f2>",
                                       "frequency 2 (will be divided by fs");
  struct arg_dbl *f3_start =  arg_dbl0("h", "f3_start", NULL, // "<f3_start>",
                                       "frequency 3 starting value (/fs");
  struct arg_dbl *f3_end =    arg_dbl0("i", "f3_end", NULL, // "<f3_end>",
                                       "frequency 3 ending value (/fs");
  struct arg_dbl *fs =        arg_dbl0("q", "fs", NULL, // "<fs>",
                                       "frequency divisor");
  struct arg_dbl *start_t =   arg_dbl0("s", "start_t", NULL, // "<start_t>",
                                       "starting time in 'seconds");
  struct arg_dbl *len_t =     arg_dbl0("l", "len_t", NULL, // "<len_t>",
                                       "time length   in 'seconds");
  struct arg_int *n =         arg_int0("n", "n", NULL, // "<n>",
                                       "number of samples of time series");
  struct arg_dbl *err_mean =  arg_dbl0("t", "err_mean", NULL, // "<err_mean>",
                                       "mean of normal random error");
  struct arg_dbl *err_stdev = arg_dbl0("u", "err_stdev", NULL, // "<err_stdev>",
                                       "standard dev of normal random error");
  struct arg_dbl *alg_mean =  arg_dbl0("v", "alg_mean", NULL, // "<alg_mean>",
                                       "mean of random nrs used in alg");
  struct arg_dbl *alg_stdev = arg_dbl0("w", "alg_stdev", NULL, // "<alg_stdev>",
                                       "standard deviation rand nrs in alg");
  struct arg_lit *help =      arg_lit0("h", "help",
                                      "Show help message and exit.");
  struct arg_lit *normal =    arg_lit0("o", "normal",
                                       "Use normal rather than uniform dist "
                                       "for rand nrs in alg");
  struct arg_lit *col_zero =  arg_lit0("z", "zero",
                                       "Test for column zero rather than zero "
                                       "variance.");
  struct arg_end *end  = arg_end(20);
  void *argtable[]     = { a, b, c, d, f1, f2, f3_start, f3_end, fs,
                           start_t, len_t, n, err_mean, err_stdev,
                           alg_mean, alg_stdev, help, normal, col_zero, end };
  if (arg_nullcheck(argtable) != 0) {
    std::cout << "Insufficient memory to parse command-line arguments.\n"
              << std::endl;
    *p_error = true;
    return;
  }
  int nr_errors = arg_parse(argc, argv, argtable);
  if (nr_errors == 0) {
    if (help->count > 0) {
      usage();
      *p_help_only = true;
      return;
    }
    if (a->count > 0)
      *p_a = a->dval[0];
    if (b->count > 0)
      *p_b = b->dval[0];
    if (c->count > 0)
      *p_c = c->dval[0];
    if (d->count > 0)
      *p_d = d->dval[0];
    if (f1->count > 0)
      *p_f1 = f1->dval[0];
    if (f2->count > 0)
      *p_f2 = f2->dval[0];
    if (f3_start->count > 0)
      *p_f3_start = f3_start->dval[0];
    if (f3_end->count > 0)
      *p_f3_end = f3_end->dval[0];
    if (fs->count > 0)
      *p_fs = fs->dval[0];
    if (start_t->count > 0)
      *p_start_t = start_t->dval[0];
    if (len_t->count > 0)
      *p_len_t = len_t->dval[0];
    if (n->count > 0)
      *p_n = n->ival[0];
    if (err_mean->count > 0)
      *p_err_mean =err_mean->dval[0];
    if (err_stdev->count > 0)
      *p_err_stdev = err_stdev->dval[0];
    if (alg_mean->count > 0)
      *p_alg_mean = alg_mean->dval[0];
    if (alg_stdev->count > 0)
      *p_alg_stdev = alg_stdev->dval[0];
    if (col_zero->count > 0)
      *p_col_zero = true;
    if (normal ->count > 0)
      *p_normal = true;

    if (*p_n <= 0) {
      std::cout << "n must be strictly positive" << std::endl;
      usage();
      *p_error = true;
      return;
    }
  } else {
    std::cout << "Incorrect usage, " << nr_errors << " errors." << std::endl;
    usage();
    *p_error = true;
    return;
  }
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
}

// Must be called with setup=true before other use!  The default initialization
// for stdev is no good.

double alg_rand(bool setup=false, bool arg_normal=false,
                double arg_mean=0.5, double arg_stdev=1.0/12.0) {
  static bool normal;
  static double mean;
  static double stdev;
  static default_random_engine rng;
  static uniform_real_distribution<double> uniform_dist(0.0, 1.0);
  static normal_distribution<double> normal_dist(0.0, 1.0);
  if (setup) {
    normal   = arg_normal;
    if (normal) {
      mean   = arg_mean;
      stdev  = arg_stdev;
      normal_dist = normal_distribution<double>(mean, stdev);
    } else {
      mean   = arg_mean;
      stdev  = arg_stdev;
      // using mean=(a+b)/2, stdev=(b-a)^2/12
      double a = mean - stdev * sqrt(3.0);
      double b = mean + stdev * sqrt(3.0);
      uniform_dist = uniform_real_distribution<double>(a, b);
    }
  }
  if (normal)
     return normal_dist(rng);
  else
    return uniform_dist(rng);
}

double calc_sum_sq_dev(vec v) {
  double mean = sum(v) / v.n_rows;
  v -= mean;
  v = v % v;
  return sum(v);
}

// Faster than calc_nr_nonzero, if don't need the graphs.

bool is_zero(vec v) {
  for (int i = 0; i < v.n_rows; i++) {
    if (v[i] != 0.0)
      return false;
  }
  return true;
}

int calc_nr_nonzero(vec v) {
  int nr_nonzero = 0;
  for (int i = 0; i < v.n_rows; i++) {
    if (v[i] != 0.0)
      nr_nonzero++;
  }
  return nr_nonzero;
}

// Faster than use_index_and_dev if don't need graphs.

bool use_index(vec v, bool col_zero) {
  static constexpr double sum_sq_dev_min = 1.0e-16;
  if (col_zero)
    return is_zero(v);
  else
    return (calc_sum_sq_dev(v) < sum_sq_dev_min);
}

UseIndexAndDevReturn use_index_and_dev(vec v, bool col_zero) {
  static constexpr double sum_sq_dev_min = 1.0e-16;
  if (col_zero) {
    int nr_nonzero = calc_nr_nonzero(v);
    return make_pair(nr_nonzero == 0,
                     double(nr_nonzero) / double(v.n_rows));
  } else {
    double sum_sq_dev = calc_sum_sq_dev(v);
    return make_pair(sum_sq_dev < sum_sq_dev_min,
                     sqrt(sum_sq_dev / double(v.n_rows)));
  }
}

AmpdReturn ampd(vec x, double alg_mean, double alg_stdev,
                    bool normal=false, bool col_zero=false,
                    bool write_files=false,
                    const char *lms_str=NULL, const char *gamma_str=NULL,
                    const char *peaks_str=NULL, const char *troughs_str=NULL)
{
  static const double  alpha = 1.0;
  int n  = x.n_rows;
  int el = ceil(n / 2) - 1;
  mat mpk(el, n, fill::zeros);
  mat mtr(el, n, fill::zeros);
  alg_rand(true, normal, alg_mean, alg_stdev); // setup=true before other use.
  for (int kix = 0; kix < el; kix++) {
    for (int i = 0; i <= kix; i++) {
      mpk(kix, i) = alpha + alg_rand();
      mtr(kix, i) = alpha + alg_rand();
    }
    for (int i = kix + 1; i < n - kix -1; i++) {
      if (x[i] <= x[i - kix - 1] || x[i] <= x[i + kix + 1])
        mpk(kix, i) = alpha + alg_rand();
      if (x[i] >= x[i - kix - 1] || x[i] >= x[i + kix + 1])
        mtr(kix, i) = alpha + alg_rand();
    }
    for (int i = n - kix - 1; i < n; i++) {
      mpk(kix, i) = alpha + alg_rand();
      mtr(kix, i) = alpha + alg_rand();
    }
  }

  vec pk_gamma(el);
  vec tr_gamma(el);
  for (int kix = 0; kix < el; kix++) {
    pk_gamma[kix] = sum(mpk.row(kix));
    tr_gamma[kix] = sum(mtr.row(kix));
  }
  double min_pk_gamma = pk_gamma[0];
  double min_tr_gamma = tr_gamma[0];
  int    min_pk_gamma_kix = 0;
  int    min_tr_gamma_kix = 0;
  for (int kix = 1; kix < el; kix++) {
    if (pk_gamma[kix] < min_pk_gamma) {
      min_pk_gamma = pk_gamma[kix];
      min_pk_gamma_kix = kix;
    }
    if (tr_gamma[kix] < min_tr_gamma) {
      min_tr_gamma = tr_gamma[kix];
      min_tr_gamma_kix = kix;
    }
  }
  int pk_lamb = min_pk_gamma_kix;
  int tr_lamb = min_tr_gamma_kix;

  Col<int> pk_zero_dev_ixs(n);
  Col<int> tr_zero_dev_ixs(n);
  int pk_zero_ixs_ix = 0;
  int tr_zero_ixs_ix = 0;
  vec pk_dev(n); // for writing files for graphs
  for (int i = 0; i < n; i++) {
    // use_index is faster if don't need graphs
    UseIndexAndDevReturn pk_ret = use_index_and_dev(mpk(span(0, pk_lamb), i),
                                                    col_zero);
    pk_dev[i] = pk_ret.second; // for writing files for graphs
    // could use use_index but want to use same code for pk and tr
    UseIndexAndDevReturn tr_ret = use_index_and_dev(mtr(span(0, tr_lamb), i),
                                                    col_zero);
    if (pk_lamb == 0 && mpk(0, i) == 0.0 || pk_lamb != 0 && pk_ret.first) {
      pk_zero_dev_ixs[pk_zero_ixs_ix] = i;
      pk_zero_ixs_ix++;
    }
    if (tr_lamb == 0 && mtr(0, i) == 0.0 || tr_lamb != 0 && tr_ret.first) {
      tr_zero_dev_ixs[tr_zero_ixs_ix] = i;
      tr_zero_ixs_ix++;
    }
  }
  
  if (write_files) {
    ofstream lms(lms_str);
    for (int k = 0; k < el; k++)
      for (int i = 0; i < n; i++)
        lms << i << " " << k << " " << mpk(k, i) << endl;
    lms.close();

    ofstream gamma(gamma_str);
    for (int k = 0; k < el; k++)
      gamma << k << " " << pk_gamma[k] / double(n) << endl;
    gamma.close();

    ofstream peaks(peaks_str);
    for (int j = 0; j < pk_zero_ixs_ix; j++) {
      int i = pk_zero_dev_ixs[j];
      peaks << i << " " << x[i] << endl;
    }
    peaks.close();

    ofstream troughs(troughs_str);
    for (int j = 0; j < tr_zero_ixs_ix; j++) {
      int i = tr_zero_dev_ixs[j];
      troughs << i << " " << x[i] << endl;
    }
    troughs.close();
}
  
  return make_tuple(el, pk_lamb, tr_lamb,
                    Col<int>(pk_zero_dev_ixs(span(0, pk_zero_ixs_ix))),
                    Col<int>(tr_zero_dev_ixs(span(0, tr_zero_ixs_ix))));
}

vec test_data(double a,  double b,  double c,  double d,
              double f1, double f2, double f3_start, double f3_end, double fs,
              double start_t, double len_t, int n,
              double err_mean, double err_stdev,
              bool write_files=false, const char *data_str=NULL) {
  static constexpr double twopi = 8.0 * atan(1.0);
  static default_random_engine rng;
  static normal_distribution<double> distribution(err_mean, err_stdev);
  vec v(n);
  double f1s       = f1 / fs;
  double f2s       = f2 / fs;
  double f3s_start = f3_start / fs;
  double f3s_end   = f3_end / fs;
  double t         = start_t;
  double t_inc     = len_t / double(n);
  double len_t_inv = 1.0 / len_t;

  for (int i = 0; i < n; i++, t += t_inc) {
    v[i] =   a * sin(twopi * f1s * t)
           + b * sin(twopi * f2s * t)
           + c * sin(twopi * (  f3s_start * (1.0 - t * len_t_inv)
                              + f3s_end   *        t * len_t_inv ) * t)
           + d * distribution(rng);
  }
  
  if (write_files) {
    ofstream data(data_str);
    for (int i = 0; i < n; i++)
      data << i << " " << v[i] << endl;
    data.close();
  }

  return v;
}

int main(int argc, char *argv[]) {
  double a         =  1.0;
  double b         =  1.0;
  double c         =  0.5;
  double d         =  0.1;
  double f1        = 10.0;
  double f2        = 70.0;
  double f3_start  =  5.0;
  double f3_end    =  5.0;
  double fs        = 12.0;
  double start_t   =  0.0;
  double len_t     =  5.0;
  int    n         = 1000;
  double err_mean  =  0.0;
  double err_stdev =  1.0;
  double alg_mean  =  0.5;
  double alg_stdev =  1.0 / 12.0;
  bool   help_only = false;
  bool   normal    = false;
  bool   col_zero  = false;
  bool   error     = false;
  get_cfg(argc, argv, &a, &b, &c, &d, &f1, &f2, &f3_start, &f3_end, &fs,
          &start_t, &len_t, &n, &err_mean, &err_stdev, &alg_mean, &alg_stdev,
          &help_only, &normal, &col_zero, &error);
  if (help_only)
    return 0;
  else if (error)
    return -1;

  ostringstream args;
  args << "ampd" << "-a" << a << "-b" << b <<"-c" << c << "-d" << d
       << "-f" << f1 << "-g" << f2 << "-h" << f3_start <<"-i" << f3_end
       << "-q" << fs << "-s" << start_t << "-l" << len_t <<"-n" << n
       << "-t" << err_mean << "-u" << err_stdev
       << "-v" << alg_mean <<"-w" << alg_stdev
       << "-o" << normal << "-z" << col_zero;
  string args_str(args.str());
  string cmds_str(args_str + ".cmd.txt");
  string data_str(args_str + ".data.txt");
  string lms_str(args_str + ".peaks.lms.txt");
  string gamma_str(args_str + ".peaks.gamma.txt");
  string peaks_str(args_str + ".peaks.txt");
  string troughs_str(args_str + ".troughs.txt");
  string png_str(args_str + ".png");
  
  vec v = test_data(a, b, c, d, f1, f2, f3_start, f3_end, fs, start_t, len_t,
                    n, err_mean, err_stdev, true, data_str.c_str());
  AmpdReturn calc =
      ampd(v, alg_mean, alg_stdev, normal, col_zero, true, lms_str.c_str(),
           gamma_str.c_str(), peaks_str.c_str(), troughs_str.c_str());
  int el = get<0>(calc);
  int pk_lamb = get<1>(calc);
  int tr_lamb = get<2>(calc); // unused below
  Col<int> peaks(get<3>(calc));   // unused below
  Col<int> troughs(get<4>(calc)); // unused below

  const int time_str_len = 20;
  char time_str[time_str_len];
  time_t curr_time = time(NULL);
  strftime(time_str, time_str_len, "%Y-%m-%d-%H:%M:%S", localtime(&curr_time));
  ofstream cmds(cmds_str);
  cmds << "# " << args_str << " " << time_str << endl;
  cmds << "reset" << endl;
  cmds << "set term png size 1300, 800" << endl;
  cmds << "set origin 0.0, 0.0" << endl;
  cmds << "set size 1.0, 1.0" << endl;
  cmds << "set output '" << png_str << "'" << endl;
  cmds << "set multiplot" << endl;
  cmds << "unset key" << endl;
  
  // Layout small plots on top, one large one underneath them.  The calculations
  // depend on this, are not more general.

  const int nr_large_plots_x = 1;
  const int nr_small_plots_x = 2;
  const int nr_plots_y       = 2;
  double margin        = 0.01;  // on all sides
  double height        =   (1.0 - (nr_plots_y + 1) * margin)
                         / double(nr_plots_y);
  double large_bottom  = margin;
  double small_bottom  = margin + height + margin;
  double large_width   =    (1.0 - (nr_large_plots_x + 1) * margin)
                          / double(nr_large_plots_x);
  double large_left    = margin;
  double small_width_fraction[nr_small_plots_x] = { 0.7, 0.3 };
  double small_width[nr_small_plots_x];
  double small_left[nr_small_plots_x];
  double running_small_left = margin;
  for (int i = 0; i < nr_small_plots_x; i++) {
    small_width[i] =   small_width_fraction[i]
                     * (1.0 - (nr_small_plots_x + 1) * margin);
    small_left[i] = running_small_left;
    running_small_left += (small_width[i] + margin);
  }
  
  cmds << "set origin " << small_left[0] << ", " << small_bottom << endl;
  cmds << "set size " << small_width[0] << ", " << height << endl;
  cmds << "set xrange[-1:" << n  << "]" << endl;
  cmds << "set yrange[-1:" << el << "]" << endl;
  cmds << "plot '" << lms_str << "' with image" << endl;

  cmds << "set title 'gamma; pk_lamb " << pk_lamb << "'" << endl;
  cmds << "set origin " << small_left[1] << ", " << small_bottom << endl;
  cmds << "set size " << small_width[1] << ", " << height << endl;
  cmds << "set xrange[-1:" << el << "]" << endl;
  cmds << "set yrange[0:3]" << endl;
  cmds << "plot '" << gamma_str << "' with lines" << endl;
  
  cmds << "set origin " << large_left << ", " << large_bottom << endl;
  cmds << "set size " << large_width << ", " << height << endl;
  cmds << "set xrange[-1:" << n << "]" << endl;
  cmds << "set yrange[-3:3]" << endl;
  cmds << "plot '"
       << troughs_str << "' with points pt 11, '"
       << peaks_str << "' with points pt 9, '"
       << data_str << "' with lines lt -1 lw 1" << endl;
  cmds << "unset multiplot" << endl;
  cmds.close();
  
  return 0;
}
