#include <Rcpp.h>
#include <complex>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include "mm_malloc.h"
using namespace Rcpp;

#define ALIGNMENT (32)

const int CONV_ROUNDING = 2048;
const int MINIMUM_SIZE_FOR_FFTW_CONVOLUTION = 80;

inline double* allocate_aligned_doubles(int n)
{
  return static_cast<double*>(_mm_malloc(n*sizeof(double), ALIGNMENT));
}

inline std::complex<double>* allocate_aligned_complexes(int n)
{
  return static_cast<std::complex<double>*>(_mm_malloc(n*sizeof(std::complex<double>), ALIGNMENT));
}


inline void free_aligned_mem(void* p)
{
  _mm_free(p);
}

template<class T>
class DoubleBuffer {
public:
  DoubleBuffer(int n, T value);
  std::vector<T>& get_src()  { return buf0_is_src ? buf0 : buf1; };
  std::vector<T>& get_dest() { return buf0_is_src ? buf1 : buf0; };
  void flip() { buf0_is_src = !buf0_is_src; };

private:
  std::vector<T> buf0;
  std::vector<T> buf1;
  bool buf0_is_src;
};

template<class T>
DoubleBuffer<T>::DoubleBuffer(int n, T value) :
  buf0(n, value), buf1(n, value), buf0_is_src(true)
{
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << std::string("[");
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << std::string("\b\b]");
  }
  return out;
}

// Computes the probability of a Poisson random variable with intensity lambda:
// Pr[Pois(lambda)=k] = e^-lambda * lambda^k / k!
inline double poisson_pmf(double lambda, int k)
{
  assert(k >= 0);
  assert(lambda >= 0.0);

  if (lambda == 0.0) {
    return k == 0 ? 1.0 : 0.0;
  }
  return std::exp(-lambda + k*std::log(lambda) - std::lgamma(k+1));
}


class PoissonPMFGenerator {
public:
  PoissonPMFGenerator(int max_k);
  ~PoissonPMFGenerator();
  double evaluate_pmf(double lambda, int k) const
  {
    assert(lambda >= 0);
    assert(k >= 0);

    if (lambda == 0) {
      return (k == 0) ? 1.0 : 0.0;
    }

    return exp(-lambda + k*log(lambda) - log_gamma_LUT[k+1]);
  };
  // Fills pmf_array_ptr with the PMF of a Poisson random variable:
  //     Pr[Pois(lambda) = 0], ..., Pr[Pois(lambda) = k]
  // Returns the smallest integer N such that for all indices i >= N we have that due to double-precision rounding
  //     Pr[Pois(lambda) = i] = 0
  void compute_array(int k, double lambda);
  const double* get_array() const {return pmf_array_ptr;}
private:
  int max_k;
  double* log_gamma_LUT;
  double* pmf_array_ptr;
};

class FFTWConvolver {
public:
  FFTWConvolver(int maximum_input_size);
  ~FFTWConvolver();
  void convolve_same_size(int size, const double* input_a, const double* input_b, double* output);
private:
  int maximum_input_size;

  std::complex<double>* tmp_complex;


  // The r2c plans perform, for various sizes, a real to complex FFT with input at r2c_in and output at r2c_out
  double* r2c_in;
  std::complex<double>* r2c_out;
  std::vector<fftw_plan> r2c_plans;
  fftw_plan memoized_r2c_plan(int rounded_size);

  // The c2r plans perform, for various sizes, a complex to real FFT with input at c2r_in and output at c2r_out
  std::complex<double>* c2r_in;
  double* c2r_out;
  std::vector<fftw_plan> c2r_plans;
  fftw_plan memoized_c2r_plan(int rounded_size);

};

void convolve_same_size(int size, const double* src0, const double* src1, double* dest)
{
  for (int j = 0; j < size; ++j) {
    double convolution_at_j = 0.0;
    for (int k = 0; k <= j; ++k) {
      convolution_at_j += src0[k] * src1[j-k];
    }
    dest[j] = convolution_at_j;
  }
}

PoissonPMFGenerator::PoissonPMFGenerator(int max_k)
{
  assert(max_k > 0);

  this->max_k = max_k;
  log_gamma_LUT = allocate_aligned_doubles(max_k+2);
  for (int i = 0; i < max_k+2; ++i) {
    log_gamma_LUT[i] = lgamma(i);
  }
  pmf_array_ptr = allocate_aligned_doubles(max_k+1);
  for (int i = 0; i < max_k+1; ++i) {
    pmf_array_ptr[i] = 0;
  }
}

PoissonPMFGenerator::~PoissonPMFGenerator()
{
  free_aligned_mem(pmf_array_ptr);
  free_aligned_mem(log_gamma_LUT);
}


void PoissonPMFGenerator::compute_array(int k, double lambda)
{
  assert(k >= 0);
  assert(k <= max_k);

  if (lambda < 0) {
    throw std::runtime_error("Expecting lambda>0 in PoissonPMFGenerator::compute_array()");
  }
  if (lambda == 0) {
    pmf_array_ptr[0] = 1;
    for (int i = 1; i < k+1; ++i) {
      pmf_array_ptr[i] = 0;
    }
    return;
  }

  double log_lambda = log(lambda);
  for (int i = 0; i < k+1; ++i) {
    pmf_array_ptr[i] = exp(-lambda + i*log_lambda - log_gamma_LUT[i+1]);
  }
}

int round_up(int n, int rounding)
{
  assert(rounding >= 0);
  return ((n+rounding-1)/rounding)*rounding;
}

FFTWConvolver::FFTWConvolver(int maximum_input_size) :
  maximum_input_size(maximum_input_size+CONV_ROUNDING-1),
  r2c_plans(round_up(2*maximum_input_size, CONV_ROUNDING)/CONV_ROUNDING, NULL),
  c2r_plans(round_up(2*maximum_input_size, CONV_ROUNDING)/CONV_ROUNDING, NULL)
{
  int maximum_padded_input_size = round_up(2*maximum_input_size, CONV_ROUNDING);

  r2c_in = allocate_aligned_doubles(maximum_padded_input_size);
  r2c_out = allocate_aligned_complexes(maximum_padded_input_size);

  c2r_in = allocate_aligned_complexes(maximum_padded_input_size);
  c2r_out = allocate_aligned_doubles(maximum_padded_input_size);

  tmp_complex = allocate_aligned_complexes(maximum_padded_input_size);
}

void convolve_same_size_naive(int size, const double* __restrict__ src0, const double* __restrict__ src1, double* __restrict__ dest)
{
  for (int j = 0; j < size; ++j) {
    double convolution_at_j = 0.0;
    for (int k = 0; k <= j; ++k) {
      convolution_at_j += src0[k] * src1[j-k];
    }
    dest[j] = convolution_at_j;

  }
}

void elementwise_complex_product(
    int size,
    const std::complex<double>* __restrict__ src0,
    const std::complex<double>* __restrict__  src1,
    std::complex<double>* __restrict__ dest,
    double multiplicative_constant)
{
  for (int i = 0; i < size; ++i) {
    dest[i] = multiplicative_constant*src0[i]*src1[i];
  }
}

fftw_plan FFTWConvolver::memoized_r2c_plan(int rounded_size)
{
  assert(rounded_size > 0);
  assert((rounded_size % CONV_ROUNDING) == 0);

  int index = rounded_size/CONV_ROUNDING - 1;
  assert(index < r2c_plans.size());

  if (r2c_plans[index] == NULL) {
    r2c_plans[index] = fftw_plan_dft_r2c_1d(rounded_size, r2c_in, reinterpret_cast<fftw_complex*>(r2c_out), FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
  }

  return r2c_plans[index];
}

fftw_plan FFTWConvolver::memoized_c2r_plan(int rounded_size)
{
  assert(rounded_size > 0);
  assert((rounded_size % CONV_ROUNDING) == 0);

  int index = rounded_size/CONV_ROUNDING - 1;
  assert(index < c2r_plans.size());

  if (c2r_plans[index] == NULL) {
    c2r_plans[index] = fftw_plan_dft_c2r_1d(rounded_size, reinterpret_cast<fftw_complex*>(c2r_in), c2r_out, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
  }

  return c2r_plans[index];
}

template<class T>
void copy_zero_padded(const T* src, T* dest, int src_size, int dest_size)
{
  memcpy(dest, src, sizeof(T)*src_size);
  memset(&dest[src_size], 0, sizeof(T)*(dest_size-src_size));
}

void FFTWConvolver::convolve_same_size(
    int size,
    const double* __restrict__ input_a,
    const double* __restrict__ input_b,
    double* __restrict__ output)
{
  if (size > maximum_input_size) {
    std::stringstream ss;
    ss << "FFTWConvolver::convolve_same_size received input of size " << size << ". This is bigger than maximum_input_size==" << maximum_input_size;
    throw std::runtime_error(ss.str());
  }
  if (size <= 0) {
    return; // Nothing to do
  }

  if (size < MINIMUM_SIZE_FOR_FFTW_CONVOLUTION) {
    convolve_same_size_naive(size, input_a, input_b, output);
    return;
  }

  int padded_size = round_up(2*size, CONV_ROUNDING);

  // tmp_complex <- FFT(zeropad(input_a));
  copy_zero_padded(input_a, r2c_in, size, padded_size);
  //fftw_execute(memoized_r2c_plan(padded_size));
  //memcpy(tmp_complex, r2c_out, sizeof(complex<double>)*(padded_size/2 + 1));
  fftw_execute_dft_r2c(memoized_r2c_plan(padded_size), r2c_in, reinterpret_cast<fftw_complex*>(tmp_complex));
  // Try this instead of last two lines:
  // fftw_execute_dft_r2c(memoized_r2c_plan(padded_length), r2c_in, tmp_complex);

  // r2c_out <- FFT(zeropad(input_b));
  copy_zero_padded(input_b, r2c_in, size, padded_size);
  fftw_execute(memoized_r2c_plan(padded_size));

  // Perform element-wise product of FFT(a) and FFT(b) and then compute inverse fourier transform.
  // FFTW returns unnormalized output. To normalize it one must divide each element of the result by the number of elements.
  elementwise_complex_product(padded_size/2 + 1, tmp_complex, r2c_out, c2r_in, 1.0/double(padded_size));
  fftw_execute(memoized_c2r_plan(padded_size));
  std::memcpy(output, c2r_out, size * sizeof(double));
}

FFTWConvolver::~FFTWConvolver()
{
  for (size_t i = 0; i < r2c_plans.size(); ++i) {
    if (r2c_plans[i] != NULL) {
      fftw_destroy_plan(r2c_plans[i]);
    }
  }

  for (size_t i = 0; i < c2r_plans.size(); ++i) {
    if (c2r_plans[i] != NULL) {
      fftw_destroy_plan(c2r_plans[i]);
    }
  }

  free_aligned_mem(r2c_in);
  free_aligned_mem(r2c_out);

  free_aligned_mem(c2r_in);
  free_aligned_mem(c2r_out);

  free_aligned_mem(tmp_complex);
}

enum BoundType {bSTEP, BSTEP, END};

struct Bound {
  double location;
  BoundType tag;
};

// Needed for using sort()
static bool operator<(Bound b0, Bound b1)
{
  return (b0.location < b1.location);
}

static std::vector<Bound> join_all_bounds(const std::vector<double>& b, const std::vector<double>& B)
{
  std::vector<Bound> bounds;
  bounds.reserve(b.size()+B.size()+1);

  Bound bound;

  for (int i = 0; i < (int)b.size(); ++i) {
    bound.location = b[i];
    bound.tag = bSTEP;
    bounds.push_back(bound);
  }

  for (int i = 0; i < (int)B.size(); ++i) {
    bound.location = B[i];
    bound.tag = BSTEP;
    bounds.push_back(bound);
  }

  sort(bounds.begin(), bounds.end());

  bound.location = 1.0;
  bound.tag = END;
  bounds.push_back(bound);

  return bounds;
}

void update_dest_buffer_and_step_counts(BoundType bound_tag, std::vector<double>& dest_buffer, int& b_step_count, int& B_step_count)
{
  if (bound_tag == bSTEP) {
    ++b_step_count;
    dest_buffer[b_step_count] = 0.0;
  } else if (bound_tag == BSTEP) {
    dest_buffer[B_step_count] = 0.0;
    ++B_step_count;
  } else {
    if (bound_tag != END) {
      throw std::runtime_error("Expecting END tag");
    }
  }
}

// TODO: Split function into 2 cases: with_fft and no_fft
std::vector<double> poisson_process_noncrossing_probability(int n, double intensity, const std::vector<double>& b, const std::vector<double>& B, bool use_fft)
{
  std::vector<Bound> bounds = join_all_bounds(b, B);

  DoubleBuffer<double> buffers(n+1, 0.0);
  buffers.get_src()[0] = 1.0;

  FFTWConvolver fftconvolver(n+1);
  PoissonPMFGenerator pmfgen(n+1);

  int b_step_count = 0;
  int B_step_count = 0;

  double prev_location = 0.0;

  for (unsigned int i = 0; i < bounds.size(); ++i) {
    int cur_size = b_step_count - B_step_count + 1;

    double lambda = intensity*(bounds[i].location-prev_location);
    if (lambda > 0) {
      pmfgen.compute_array(cur_size, lambda);
      if (use_fft) {
        fftconvolver.convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[B_step_count], &buffers.get_dest()[B_step_count]);
      } else {
        convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[B_step_count], &buffers.get_dest()[B_step_count]);
      }
      update_dest_buffer_and_step_counts(bounds[i].tag, buffers.get_dest(), b_step_count, B_step_count);
      buffers.flip();
    } else if (lambda==0) {
      // No need to convolve or copy anything -- just modify src buffer in place.
      update_dest_buffer_and_step_counts(bounds[i].tag, buffers.get_src(), b_step_count, B_step_count);
    } else {
      throw std::runtime_error("lambda<0 in poisson_process_noncrossing_probability(). This should never happen.");
    }
    prev_location = bounds[i].location;
  }
  return buffers.get_src();
}

// [[Rcpp::export]]
double fft_get_level_from_bounds_two_sided(const std::vector<double>& b, const std::vector<double>& B)
{
  int n = b.size();

  std::vector<double> poisson_nocross_probs = poisson_process_noncrossing_probability(n, n, b, B, true);

  return poisson_nocross_probs[n] / poisson_pmf(n, n);
}
