/*
 FILENAME:		KickMath.h
 AUTHOR:        Orlando S. Hoilett, Alyson S. Pickering, and Akio K. Fujita
 EMAIL:     	orlandohoilett@gmail.com
 
 
 Please see .cpp file for extended descriptions, instructions, and version updates
 
 
 DISCLAIMER
 Linnes Lab code, firmware, and software is released under the
 MIT License (http://opensource.org/licenses/MIT).
 
 The MIT License (MIT)
 
 Copyright (c) 2019 Linnes Lab, Purdue University
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 
 */



#ifndef KickMath_h
#define KickMath_h


//Standard Arduino libraries
#include <Arduino.h>

//Custom Kick Libraries
#include "ttestTable.h"


class KickMath
{

public:
	
	static int32_t calcSqrt(int32_t num);
	static uint16_t calcMagnitude(int16_t x, int16_t y);
	static uint16_t calcMagnitude(int16_t x, int16_t y, int16_t z);
	
	
	static int16_t getMax(uint16_t samples, const int16_t data[]);
	static int16_t getMin(uint16_t samples, const int16_t data[]);
	static int32_t getSum(uint16_t samples, const int16_t data[]);
	
	
	static uint16_t getMaxIndex(uint16_t samples, const int32_t data[]);
	static void getMaxIndex(uint16_t samples, const int32_t data[], uint16_t &max1, uint16_t &max2);
	static void getMaxIndex(uint16_t samples, const int32_t data[], uint16_t maxes[], const uint16_t num_maxes);
	
	static uint16_t getMinIndex(uint16_t samples, const int32_t data[]);

	
	static int16_t calcAverage(uint16_t samples, const int16_t data[]);
	static float calcStDev(uint16_t samples, const int32_t data[]);
	
	static int16_t calcPeaktoPeak(uint16_t samples, const int16_t data[]);
	//static int16_t calcPeaktoPeak(uint16_t samples, const int32_t data[]);
	static float calcRMS(uint16_t samples, const int16_t data[]);
	
	
	static void calcDerivative(uint16_t samples, const int16_t data[], int16_t ddx[]);
	static void calcDerivative(uint16_t samples, uint16_t dt, const int16_t data[], float ddx[]);
	static void calcDerivative(uint16_t samples, const int16_t data[], float ddx[], const uint32_t t[]);

	
	static float calcCentroid(float fs, uint16_t samples, const int32_t mag[], uint16_t fcenter, uint16_t width);
	
	
	static bool tTest(const int16_t data1[], const int16_t data2[], const uint8_t samples, const float alpha);
	
};


#endif /* KickMath_h */


