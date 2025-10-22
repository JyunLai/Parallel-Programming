#include "PPintrin.h"
#include <cstring>
#include <iostream>
#include <ostream>

// implementation of absSerial(), but it is vectorized using PP intrinsics
void absVector(float *values, float *output, int N)
{
  __pp_vec_float x;
  __pp_vec_float result;
  __pp_vec_float zero = _pp_vset_float(0.f);
  __pp_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {

    // All ones
    maskAll = _pp_init_ones();

    // All zeros
    maskIsNegative = _pp_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _pp_vload_float(x, values + i, maskAll); // x = values[i];

    // Set mask according to predicate
    _pp_vlt_float(maskIsNegative, x, zero, maskAll); // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _pp_vsub_float(result, zero, x, maskIsNegative); //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _pp_mask_not(maskIsNegative); // } else {

    // Execute instruction ("else" clause)
    _pp_vload_float(result, values + i, maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _pp_vstore_float(output + i, result, maskAll);
  }
}

void clampedExpVector(float *values, int *exponents, float *output, int N)
{
  //
  // PP STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //
  __pp_vec_float vecValues, vecOutput;
  __pp_vec_int vecExponents;
  __pp_vec_float vecMax = _pp_vset_float(9.999999f);
  __pp_vec_float vecMin = _pp_vset_float(-9.999999f);
  __pp_mask maskAll, maskGreaterThanMax, maskLessThanMin;

  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {
    int batch_size = (i + VECTOR_WIDTH <= N) ? VECTOR_WIDTH : N - i;
    maskAll = _pp_init_ones(batch_size);

    _pp_vload_float(vecValues, values+i, maskAll);
    _pp_vload_int(vecExponents, exponents+i, maskAll);

    _pp_vset_float(vecOutput, 1.0f, maskAll);

    for (int j = 0; j < VECTOR_WIDTH; ++j) {
      int exponent = vecExponents.value[j];
      float result = 1.0f;
      for (int k = 0; k < exponent; ++k) {
        result *= vecValues.value[j];
      }
      vecOutput.value[j] = result;
    }

    _pp_vgt_float(maskGreaterThanMax, vecOutput, vecMax, maskAll);
    _pp_vmove_float(vecOutput, vecMax, maskGreaterThanMax);

    _pp_vlt_float(maskLessThanMin, vecOutput, vecMin, maskAll);
    _pp_vmove_float(vecOutput, vecMin, maskLessThanMin);

    _pp_vstore_float(output+i, vecOutput, maskAll);
  }
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N)
{
  //
  // PP STUDENTS TODO: Implement your vectorized version of arraySumSerial here
  //
  __pp_vec_float sumVec, loaded;
  __pp_mask mask;

  mask = _pp_init_ones();
  _pp_vset_float(sumVec, 0.f, mask);

  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    int width = (i + VECTOR_WIDTH <= N) ? VECTOR_WIDTH : (N - i);
    mask = _pp_init_ones();
    
    if (width < VECTOR_WIDTH) {
      for (int j = width; j < VECTOR_WIDTH; ++j) {
        mask.value[j] = 0;
      }
    }

    _pp_vload_float(loaded, values+i, mask);
    _pp_vadd_float(sumVec, sumVec, loaded, mask);
  }

  int steps = (int)log2(VECTOR_WIDTH);
  for (int i = 1; i < steps; i++) {
    _pp_hadd_float(sumVec, sumVec);
    _pp_interleave_float(sumVec, sumVec);
  }
  _pp_hadd_float(sumVec, sumVec);

  float result[VECTOR_WIDTH];
  mask = _pp_init_ones();
  _pp_vstore_float(result, sumVec, mask);
  return result[0];
}