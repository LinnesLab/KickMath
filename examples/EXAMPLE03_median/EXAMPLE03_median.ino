/*
 * FILENAME: EXAMPLE03_median.ino
 * AUTHOR:   Orlando S. Hoilett
 * CONTACT:  orlandohoilett@gmail.com
 * VERSION:  1.0.0
 * 
 * 
 * AFFILIATIONS
 * Linnes Lab, Weldon School of Biomedical Engineering,
 * Purdue University, West Lafayette, IN 47907
 * 
 * 
 * DESCRIPTION
 * Basic test of the KickMath class to calculate median.
 * This is a templated, static class, so function calls must be preceded with
 * KickMath<variable_type>:: where variable_type should be replaced with
 * int16_t, int, float, etc.
 * 
 * 
 * UPDATES
 * Version 1.0.0
 * 2020/08/22:1651>
 *           - Initial creation.
 * 
 * 
 * DISCLAIMER
 * Linnes Lab code, firmware, and software is released under the
 * MIT License (http://opensource.org/licenses/MIT).
 * 
 * The MIT License (MIT)
 * 
 * Copyright (c) 2020 Linnes Lab, Purdue University
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 */


#include "KickMath.h"


const uint16_t samples = 100;
uint16_t input[samples] = { 3, 53, 70, 56, 18, 85, 27, 14, 37, 94, 9, 55, 40, 60, 52, 61, 15, 65, 13, 8, 57, 97, 69, 4, 35, 82, 22, 73, 59, 68, 78, 24, 21, 36, 71, 80, 74, 39, 17, 12, 29, 76, 49, 51, 30, 90, 88, 2, 84, 50, 62, 28, 77, 43, 5, 16, 58, 26, 32, 34, 1, 75, 66, 95, 38, 89, 67, 87, 100, 54, 92, 81, 25, 83, 46, 33, 23, 45, 96, 99, 79, 48, 11, 31, 7, 6, 19, 91, 93, 44, 47, 98, 86, 41, 63, 20, 72, 10, 42, 64 };

const uint16_t samples2 = 7;
uint16_t input2[samples] = { 3, 53, 70, 56, 18, 85, 27 };


void setup()
{
  Serial.begin(9600);
  while(!Serial); //will not run until Serial Monitor or Plotter is open


  uint16_t tmp[samples] = {0};
  uint16_t med = KickMath<uint16_t>::calcMedian(samples, input, tmp);
  Serial.print("Median: "); Serial.print(med); Serial.println();


  uint16_t tmp2[samples2] = {0};
  uint16_t med2 = KickMath<uint16_t>::calcMedian(samples2, input2, tmp2);
  Serial.print("Median: "); Serial.print(med2); Serial.println();
}


void loop()
{
}
