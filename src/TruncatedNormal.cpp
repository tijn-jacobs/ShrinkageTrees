/* MIT License
 * Copyright (c) 2024 Tijn Jacobs
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "TruncatedNormal.h"

double rtnorm(double mean, double truncation, double sd, Random& gen) {

  // Normalize the truncation point
  truncation = (truncation - mean) / sd;

  double z, lambda;

  if (truncation <= 0.0) {
    // Case 1: Truncation point is less than or equal to 0
    // Perform rejection sampling from a standard normal distribution
    do {
      z = gen.normal();
    } while (z < truncation);
  } else {
    // Case 2: Truncation point is greater than 0
    // Use exponential proposal distribution for efficient sampling

    // Calculate the optimal exponential rate parameter
    lambda = 0.5 * (truncation + sqrt(truncation * truncation + 4.0));

    // Perform rejection sampling
    do {
      z = gen.exp() / lambda + truncation;  // Sample from exponential distribution
    } while (gen.uniform() > exp(-0.5 * pow(z - lambda, 2.0)));  // Accept/reject step
  }

  // Scale back to the original distribution
  double sample = z * sd + mean;

  return sample;
}
