/*
 FILENAME:      KickMath.h
 AUTHOR:        Orlando S. Hoilett, Alyson S. Pickering, and Akio K. Fujita
 EMAIL:     	orlandohoilett@gmail.com
 VERSION:		3.0.0
 
 
 DESCRIPTION
 A library for performing a few simple mathematical calculations and for use
 with arrays. Functions include max and min detection, square root, centroid,
 derivatives, etc.
 
 This is a static class. Function calls must be preceded with the class name
 and scope resolution operator as follows KickMath::
 
 
 UPDATES
 2019/12/03:2157>
 			Initiated.
 2020/03/06:2300>
 			Added t-test function for matached paired, two-tailed test at
 			alpha = 0.05.
 VERSION 2.0.0
 2020/07/04:1018> (UTC-5)
 			- Added a lower processing power square root function.
 			- Added getMaxIndex and a getMinIndex functions.
 			- Moved Akio's peak detection functions to the KickPeaks class.
 2020/07/12:0658> (UTC-5)
 			- Updated comments.
 Version 3.0.0
 2020/08/18:1143> (UTC-5)
 			- moving to a templated class
 
 
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


template<typename Type>


class KickMath
{

public:
	
	static int32_t calcSqrt(Type num);
	static int32_t calcMagnitude(Type x, Type y);
	static int32_t calcMagnitude(Type x, Type y, Type z);
	
	
	static Type getMax(uint16_t samples, const Type data[]);
	static Type getMin(uint16_t samples, const Type data[]);
	static float getSum(uint16_t samples, const Type data[]);
	
	
	static uint16_t getMaxIndex(uint16_t samples, const Type data[]);
	static void getMaxIndex(uint16_t samples, const Type data[], uint16_t &max1, uint16_t &max2);
	static void getMaxIndex(uint16_t samples, const Type data[], uint16_t maxes[], const uint16_t num_maxes);
	
	static uint16_t getMinIndex(uint16_t samples, const Type data[]);

	
	static Type calcAverage(uint16_t samples, const Type data[]);
	static float calcStDev(uint16_t samples, const Type data[]);
	
	static Type calcPeaktoPeak(uint16_t samples, const Type data[]);
	static Type calcRMS(uint16_t samples, const Type data[]);
	
	
	static void calcDerivative(uint16_t samples, const Type data[], Type ddx[]);
	static void calcDerivative(uint16_t samples, uint16_t dt, const Type data[], float ddx[]);
	static void calcDerivative(uint16_t samples, const Type data[], float ddx[], const uint32_t t[]);

	
	static float calcCentroid(float fs, uint16_t samples, const Type mag[], uint16_t fcenter, uint16_t width);
	
	
	static bool tTest(const Type data1[], const Type data2[], const uint16_t samples, const float alpha);
	
};



//int32_t KickMath::isqrt(int32_t num)
//num to calculate square root of
//
//A PDF of the sources are also included in the "extras/references/" folder
//
//Source: https://en.wikipedia.org/wiki/Methods_of_computing_square_roots
//https://web.archive.org/web/20120306040058/http://medialab.freaknet.org/martin/src/sqrt/sqrt.c
template<typename Type>
int32_t KickMath<Type>::calcSqrt(Type val)
{
	int32_t num = (int32_t)val;
	
	if(num < 0) return 0;
	
	int32_t res = 0;
	int32_t bit = 1 << 30; // The second-to-top bit is set.
	// Same as ((unsigned) INT32_MAX + 1) / 2.
	
	// "bit" starts at the highest power of four <= the argument.
	while (bit > num)
		bit >>= 2;
	
	while (bit != 0) {
		if (num >= res + bit) {
			num -= res + bit;
			res = (res >> 1) + bit;
		} else
			res >>= 1;
		bit >>= 2;
	}
	return res;
}


//uint16_t KickMath::calcMagnitude(int16_t x, int16_t y)
//x				x component of the 2D vector
//y				y component of the 2D vector
//
//Calculates the magnitude of a 2D vector
template<typename Type>
int32_t KickMath<Type>::calcMagnitude(Type x, Type y)
{
	return calcSqrt(x*x + y*y);
}


//uint16_t KickMath::calcMagnitude(int16_t x, int16_t y, int16_t z)
//x				x component of the 3D vector
//y				y component of the 3D vector
//z				z component of the 3D vector
//
//Calculates the magnitude of a 3D vector
template<typename Type>
int32_t KickMath<Type>::calcMagnitude(Type x, Type y, Type z)
{
	return calcSqrt(x*x + y*y + z*z);
}


//int16_t KickMath::getMax(uint16_t samples, const int16_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//Finds the max value within an input array.
template<typename Type>
Type KickMath<Type>::getMax(uint16_t samples, const Type data[])
{
	//Fence post solution: assume the first value
	//is the max then compare and update from there
	Type max = data[0];
	
	for(uint16_t i = 1; i < samples; i++)
	{
		if (data[i] > max) max = data[i];
	}
	
	return max;
}


//uint16_t KickMath::getMaxIndex(uint16_t samples, const int32_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//Finds the index of the max value within an input array.
template<typename Type>
uint16_t KickMath<Type>::getMaxIndex(uint16_t samples, const Type data[])
{
	//Fence post solution: assume the first value
	//is the max then compare and update from there
	uint16_t imax = 0;
	
	for(uint16_t i = 1; i < samples; i++)
	{
		if (data[i] > data[imax]) imax = i;
	}
	
	return imax;
}


//void KickMath::getMaxIndex(uint16_t samples, const int32_t data[], uint16_t &max1, uint16_t &max2)
//samples		number of samples within the array
//data			input array containing signal
//max1			variable to store the index of the first highest value
//max2			variable to store the index of the second highest value
//
//Finds the index of the two highest values within an input array and stores
//them in the max1 and max2 variables.
//
//Could be used as a general function and return an array of up to 10 or so maxes (do hardcode a limit)
template<typename Type>
void KickMath<Type>::getMaxIndex(uint16_t samples, const Type data[], uint16_t &max1, uint16_t &max2)
{
	//Fence post solution: assume the first value
	//is the max then compare and update from there
	max1 = 0;
	max2 = 0;
	
	for(uint16_t i = 1; i < samples; i++)
	{
		if (data[i] > data[max1])
		{
			max2 = max1;
			max1 = i;
		}
	}
}


//void KickMath::getMaxIndex(uint16_t samples, const int32_t data[], uint16_t maxes[], const uint16_t num_maxes)
//samples		number of samples within the array
//data			input array containing signal
//maxes			array to store max values. C++ allows us to modify arrays used
//					as function parameters
//num_maxes		desired number of the top values requested by function caller
//
//Finds the indices of the max values in a function and returns the indices in an
//array. Indices are ordered from highest max at the 0th index to lowest max in
//the last index.
//
//C++ allows us to modify arrays used as function parameters so the maxes
//are stored in the "maxes" array.
template<typename Type>
void KickMath<Type>::getMaxIndex(uint16_t samples, const Type data[], uint16_t maxes[], const uint16_t num_maxes)
{
	//Fence post solution: assume the first value
	//is the max then compare and update from there
	for(uint16_t i = 0; i < num_maxes; i++)
	{
		maxes[i] = 0;
	}
	
	
	for(uint16_t i = 1; i < samples; i++)
	{
		if (data[i] > data[maxes[0]])
		{
			//shifting the maxes in the array if a new max is found
			for(uint16_t j = 1; j < num_maxes; j++)
			{
				maxes[num_maxes-j] = maxes[num_maxes-j-1];
			}
			
			//update index 0 which contains the index of the highest max.
			maxes[0] = i;
		}
	}
}


//int16_t KickMath::getMin(uint16_t samples, const int16_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//Finds the minimum value within an input array.
template<typename Type>
Type KickMath<Type>::getMin(uint16_t samples, const Type data[])
{
	//Fence post solution: assume the first value
	//is the min then compare and update from there
	Type min = data[0];
	
	for(uint16_t i = 1; i < samples; i++)
	{
		if (data[i] < min) min = data[i];
	}
	
	return min;
}


//uint16_t KickMath::getMinIndex(uint16_t samples, const int32_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//Finds the index of the min value within an input array.
template<typename Type>
uint16_t KickMath<Type>::getMinIndex(uint16_t samples, const Type data[])
{
	//Fence post solution: assume the first value
	//is the min then compare and update from there
	uint16_t imin = 0;
	
	for(uint16_t i = 1; i < samples; i++)
	{
		if (data[i] < data[imin]) imin = i;
	}
	
	return imin;
}


//int32_t KickMath::getSum(uint16_t samples, const int16_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//Returns the sum of values in a given input array.
template<typename Type>
float KickMath<Type>::getSum(uint16_t samples, const Type data[])
{
	float sum = 0;
	
	for(uint16_t i = 0; i < samples; i++) sum += data[i];
	
	return sum;
}


//int16_t KickMath::calcAverage(uint16_t samples, const int16_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//Calculates the average of the values in the given input array.
template<typename Type>
Type KickMath<Type>::calcAverage(uint16_t samples, const Type data[])
{
	return getSum(samples,data)/(float)samples;
}


//float KickMath::calcStDev(uint16_t samples, const int32_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//return		standard deviation of signal stored in data array
//
//Calculates the standard deviation of the values in the given input array.
template<typename Type>
float KickMath<Type>::calcStDev(uint16_t samples, const Type data[])
{
	//float avg = calcAverage(samples, data);
	float avg = 0;
	for(uint16_t i = 0; i < samples; i++)
	{
		avg += data[i];
	}
	avg = avg/samples;
	
	
	float num = 0;
	for(uint16_t i = 0; i < samples; i++)
	{
		num += pow(data[i]-avg, 2); //or just num*num? which is faster?
	}
	
	return sqrt(num/(samples-1)); 
}


//int16_t KickMath::calcPeaktoPeak(uint16_t samples, const int16_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//return		peak-to-peak value of signal stored in data array
//
//Calculates the peak-to-peak amplitude of the values in the given input array.
//It would be faster to avoid calling the getMax and getMin functions separately
//to avoid scanning through the array twice. Will update in a future version.
template<typename Type>
Type KickMath<Type>::calcPeaktoPeak(uint16_t samples, const Type data[])
{
	return getMax(samples, data) - getMin(samples, data);
}


//float KickMath::calcRMS(uint16_t samples, const int16_t data[])
//samples		number of samples within the array
//data			input array containing signal
//
//return		root-mean-squared of signal stored in data array
//
//Calculates the root mean square of the signal across a given number of
//samples in accordance to the equation as defined here:
//<https://en.wikipedia.org/wiki/Root_mean_square>
//This method first calculates the average and removes the mean from each data
//point, then calculates RMS.
template<typename Type>
Type KickMath<Type>::calcRMS(uint16_t samples, const Type data[])
{
	Type avg = calcAverage(samples, data);
	
	
	float sum_of_squares = 0;
	for(uint16_t i = 0; i < samples; i++)
	{
		sum_of_squares += pow(data[i]-avg, 2);
	}
	
	//maybe switch to faster square root method
	return sqrt(sum_of_squares/samples);
}


//void KickMath::calcDerivative(uint16_t samples, const int16_t data[], int16_t ddx[])
//samples		number of samples within the array
//data			input array containing signal
//ddx			the array for storing the derivatives. C++ allows us to modify
//					arrays used as function parameters
//
//Calculates the derivative of the given input signal assuming the change in time
//is 1. This simply becomes the difference between subsequent values in the
//given input array.
//
//C++ allows us to modify arrays used as function parameters so the derivatives
//are stored in the "ddx" array.
template<typename Type>
void KickMath<Type>::calcDerivative(uint16_t samples, const Type data[], Type ddx[])
{
	//The first value in the ddx array is zero since the derivative needs two
	//values before a valid derivaitve is can be determined. As a result, the
	//derivative array has samples-1 valid answers. The first value is not useful.
	ddx[0] = 0;
	for (uint16_t i = 1; i < samples; i++)
	{
		ddx[i] = (data[i]-data[i-1]);
	}
}


//void KickMath::calcDerivative(uint16_t samples, uint16_t dt, const int16_t data[], float ddx[])
//samples		number of samples within the array
//dt			the difference in time between samples. Assumes constant
//					sampling rate for the entire array.
//data			input array containing signal
//ddx			the array for storing the derivatives. C++ allows us to modify
//					arrays used as function parameters
//
//Calculates the derivative of the given input signal assuming the change in time
//is dt.
//
//C++ allows us to modify arrays used as function parameters so the derivatives
//are stored in the "ddx" array.
template<typename Type>
void KickMath<Type>::calcDerivative(uint16_t samples, uint16_t dt, const Type data[], float ddx[])
{
	//The first value in the ddx array is zero since the derivative needs two
	//values before a valid derivaitve is can be determined. As a result, the
	//derivative array has samples-1 valid answers. The first value is not useful.
	ddx[0] = 0;
	for (uint16_t i = 1; i < samples; i++)
	{
		ddx[i] = (data[i]-data[i-1]) / (float)dt;
	}
}


//void KickMath::calcDerivative(uint16_t samples, const int16_t data[], float ddx[], const uint32_t t[])
//samples		number of samples within the array
//data			input array containing signal
//ddx			the array for storing the derivatives. C++ allows us to modify
//					arrays used as function parameters
//t				time array. holds the time each sample in data was taken. t[1]
//					corresponds to the time at which data[1] was samples. t[10]
//					corresponds to data[10] and so on.
//
//Calculates the derivative of the given input signal assuming the change in time
//is dt.
//
//C++ allows us to modify arrays used as function parameters so the derivatives
//are stored in the "ddx" array.
template<typename Type>
void KickMath<Type>::calcDerivative(uint16_t samples, const Type data[], float ddx[], const uint32_t t[])
{
	//The first value in the ddx array is zero since the derivative needs two
	//values before a valid derivaitve is can be determined. As a result, the
	//derivative array has samples-1 valid answers. The first value is not useful.
	ddx[0] = 0;
	for (uint16_t i = 1; i < samples; i++)
	{
		ddx[i] = (data[i]-data[i-1]) / (float)(t[i]-t[i-1]);
	}
}


//float KickMath::calcCentroid(float fs, uint16_t samples, const int32_t mag[], uint16_t fcenter, uint16_t width)
//fs			sampling frequency for data array
//samples		number of samples within the array
//mag			frequency weights
//fcenter		the index of the mag array to calculate the centroid around.
//width			the number of indices around the index for determinig centroid
//
//This metods calculates the centroid of a certain section of the frequency spectrum.
//
//The centroid is calculated by perfoming a weigthed average.
//https://en.wikipedia.org/wiki/Spectral_centroid
template<typename Type>
float KickMath<Type>::calcCentroid(float fs, uint16_t samples, const Type mag[], uint16_t fcenter, uint16_t width)
{
	float denominator = 0;
	float numerator = 0;
	float fs_scale = fs/samples; //frequency increments
	
	
	//if width exceeds the bounds of the array, limit width
	if(fcenter+width > samples-1) width = samples-fcenter;
	if(fcenter-width < 0) width = fcenter-0;
	
	
	for(uint16_t i = fcenter-width; i <= fcenter+width; i++)
	{
		denominator += mag[i];
		numerator += mag[i]*i*fs_scale;
	}
	
	
	return numerator/denominator;
}


//bool KickMath::tTest(const int16_t data1[], const int16_t data2[], const uint8_t samples, const float alpha)
//data1			pre-data
//data2			post-data
//samples		number of observations
//alpha			alpha (significance) level...currently only for alpha = 0.05
//
//returns true if reject null hypothesis
//returns false if fail to reject null hypothesis
//
//Runs a matched paired t-test and returns true if the the null hypothesis is rejected
//Reference saved here: /extras/references/paired-t-test.pdf
template<typename Type>
bool KickMath<Type>::tTest(const Type data1[], const Type data2[], const uint16_t samples, const float alpha)
{
	//if (samples > 30) return -1;
	float diff[30] = {};
	const uint16_t df = samples-1; //degress of freedom
	
	
	float sum = 0;
	for (uint16_t i = 0; i < samples; i++)
	{
		diff[i]= data2[i] - data1[i];
		sum += diff[i];
	}
	
	
	float meanDiff = sum/samples;
	float SdNum = 0;
	for (int i = 0; i < samples; i++)
	{
		SdNum += sq(diff[i] - meanDiff);
	}
	float std = sqrt(SdNum/df); //should this be samples or degress of freedom??
	float SE = std/sqrt(samples); //standard error
	
	
	//t-statistic = meanDiff/SEfloat
	//a95 is look up table for up to 35 degrees of freedom at
	//alpha = 0.05 for two-tailed test
	return abs(meanDiff/SE) > a95[df];
}



#endif /* KickMath_h */


