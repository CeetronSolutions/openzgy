// Copyright 2017-2021, Schlumberger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#if defined(_MSC_VER) && !defined(_MANAGED)
#if defined( _M_X64 ) || (_M_IX86_FP >= 2) // we can use SSE2 if we are on a x64 platform or we are on x86 and the /arch flag is set to at least SSE2 
#define USE_SSE2_INSTRUCTIONS
#endif
#endif

#ifdef __GNUC__
#if defined(__x86_64) || (defined(__i386__) && defined(__SSE2__)) // we can use SSE2 if we are on an x64 platform or we are on x86 and the -m flag is set to at least -msse2
#define USE_SSE2_INSTRUCTIONS
#endif
#endif

#ifdef USE_SSE2_INSTRUCTIONS
#include <emmintrin.h> // intel SSE2 intrinsics used for min / max and summing calculations
#endif

#include "minmaxscan.h"
#include <limits>
#include <cstring>

namespace {

#ifdef USE_SSE2_INSTRUCTIONS

    inline void reduceMinMax(float& minValue, float& maxValue, __m128 minVec, __m128 maxVec)
    {
        // find the single minimum/maximum value
        // we do this by performing 2 shuffles and two min/max operations
        // this will set all the elements in minValue to the minimum actual value and all the elements in maxValue to all the maximum value

        minVec = _mm_min_ps(minVec, _mm_shuffle_ps(minVec, minVec, _MM_SHUFFLE(2, 1, 0, 3)));
        minVec = _mm_min_ps(minVec, _mm_shuffle_ps(minVec, minVec, _MM_SHUFFLE(1, 0, 3, 2)));

        maxVec = _mm_max_ps(maxVec, _mm_shuffle_ps(maxVec, maxVec, _MM_SHUFFLE(2, 1, 0, 3)));
        maxVec = _mm_max_ps(maxVec, _mm_shuffle_ps(maxVec, maxVec, _MM_SHUFFLE(1, 0, 3, 2)));

        _mm_store_ss(&minValue, minVec);
        _mm_store_ss(&maxValue, maxVec);
    }

#endif

    inline void stridedMinMax(const float* const values, const size_t size, const size_t stride, float& minValue, float& maxValue)
    {
        size_t index = 0;
#ifdef USE_SSE2_INSTRUCTIONS
        __m128 minVal = _mm_set_ss(std::numeric_limits<float>::infinity());
        __m128 maxVal = _mm_set_ss(-std::numeric_limits<float>::infinity());

        // loop through the array, one element at a time.
        for (size_t i = 0; i < size; i++)
        {
            __m128 temp = _mm_set_ss(values[index]);
            minVal = _mm_min_ss(minVal, temp);
            maxVal = _mm_max_ss(maxVal, temp);
            index += stride;
        }
        // store the values out of the register, into the minValue / maxValue array.
        _mm_store_ss(&minValue, minVal);
        _mm_store_ss(&maxValue, maxVal);
#else
        minValue = std::numeric_limits<float>::infinity();
        maxValue = -std::numeric_limits<float>::infinity();

        for (size_t i = 0; i < size; i++)
        {
            if (values[index] < minValue) minValue = values[index];
            if (values[index] > maxValue) maxValue = values[index];
            index += stride;
        }
#endif
    }

#ifdef USE_SSE2_INSTRUCTIONS
    // fills minValue and maxValue with the smallest / largest element from values. Generally faster for large arrays
    // size must be at least 4.
    inline void vectorMinMax(const float* const values, const size_t size, float& minValue, float& maxValue)
    {
        size_t unrollSize = (size / 4) * 4;

        __m128 minVals = _mm_loadu_ps(values);
        __m128 maxVals = minVals;
        for (size_t i = 4; i < unrollSize; i += 4)
        {
            __m128 temp(_mm_loadu_ps(&values[i]));

            minVals = _mm_min_ps(temp, minVals);
            maxVals = _mm_max_ps(temp, maxVals);
        }

        for (size_t i = unrollSize; i < size; i++)
        {
            __m128 temp = _mm_set_ss(values[i]);
            minVals = _mm_min_ss(minVals, temp);
            maxVals = _mm_max_ss(maxVals, temp);
        }

        reduceMinMax(minValue, maxValue, minVals, maxVals);
    }
#endif

    inline bool isNotNaN(const float val)
    {
        // Assume that sizeof(int) == sizeof(float) and that we are IEEE arithmetic
        //
        // Disable these asserts for now because they fail under the old versions of GCC.
        //
        //  static_assert(sizeof(float) == sizeof(unsigned int), "sizeof(float) is expected to equal sizeof(unsigned int)");
        //  static_assert(std::numeric_limits<float>::is_iec559, "Not IEEE float");

        unsigned int punnedVal;

        // type-punning via memcpy is the *only* safe way of doing direct conversion (without a cast)
        // if we actually had to do a memcpy, this would be quite slow, thankfully the compiler optimizes the call away
        memcpy(&punnedVal, &val, sizeof(float));

        // since we know that IEEE expects that inf / NaN has an all one exponent, we can simply mask out everything except the exponent
        // and then check that the integer reinterpretation of the float is not equal to an all one exponent (conveniently this is the same value). If that is the case, we are not inf / NaN
        return ((punnedVal & 0x7FFFFFFF) <= 0x7F800000);
    }

    inline void stridedMinMaxSafe(const float* const values, const size_t size, const size_t stride, float& minValue, float& maxValue)
    {
        size_t index = 0;
#ifdef USE_SSE2_INSTRUCTIONS
        __m128  pInfinity = _mm_set1_ps(std::numeric_limits<float>::infinity());
        __m128  nInfinity = _mm_set1_ps(-std::numeric_limits<float>::infinity());
        __m128  minVal = _mm_set_ss(std::numeric_limits<float>::infinity());
        __m128  maxVal = _mm_set_ss(-std::numeric_limits<float>::infinity());

        for (size_t i = 0; i < size; i++)
        {
            __m128 temp = _mm_set_ss(values[index]);
            __m128 isValid = _mm_and_ps(_mm_and_ps(_mm_cmpord_ps(temp, temp), _mm_cmpneq_ps(temp, pInfinity)), _mm_cmpneq_ps(temp, nInfinity));

            __m128 maskedMinTemp = _mm_or_ps(_mm_and_ps(isValid, temp), _mm_andnot_ps(isValid, minVal)); // any NaN values should be "masked" out
            __m128 maskedMaxTemp = _mm_or_ps(_mm_and_ps(isValid, temp), _mm_andnot_ps(isValid, maxVal)); // any NaN values should be "masked" out

            minVal = _mm_min_ss(minVal, maskedMinTemp);
            maxVal = _mm_max_ss(maxVal, maskedMaxTemp);

            index += stride;
        }

        _mm_store_ss(&minValue, minVal);
        _mm_store_ss(&maxValue, maxVal);
#else

        minValue = std::numeric_limits<float>::infinity();
        maxValue = -std::numeric_limits<float>::infinity();
        for (size_t i = 0; i < size; i++)
        {
            if (isNotNaN(values[index]))
            {
                if (values[index] < minValue) minValue = values[index];
                if (values[index] > maxValue) maxValue = values[index];
            }

            index += stride;
        }

#endif
    }

#ifdef USE_SSE2_INSTRUCTIONS

    // fills minValue and maxValue with the smallest / largest element from values, ignoring NaN values. Generally faster for large arrays
    // size must be at least 4.
    inline void vectorMinMaxSafe(const float* const values, const size_t size, float& minValue, float& maxValue)
    {
        size_t unrollSize = (size / 4) * 4;

        __m128  pInfinity = _mm_set1_ps(std::numeric_limits<float>::infinity());
        __m128  nInfinity = _mm_set1_ps(-std::numeric_limits<float>::infinity());
        __m128  minVals = _mm_set1_ps(std::numeric_limits<float>::infinity());
        __m128  maxVals = _mm_set1_ps(-std::numeric_limits<float>::infinity());

        for (size_t i = 0; i < unrollSize; i += 4)
        {
            __m128 temp(_mm_loadu_ps(&values[i]));
            
            __m128 isValid = _mm_and_ps(_mm_and_ps(_mm_cmpord_ps(temp, temp), _mm_cmpneq_ps(temp, pInfinity)), _mm_cmpneq_ps(temp, nInfinity));

            __m128 maskedMinTemp = _mm_or_ps(_mm_and_ps(isValid, temp), _mm_andnot_ps(isValid, minVals)); // any NaN values should be "masked" out
            __m128 maskedMaxTemp = _mm_or_ps(_mm_and_ps(isValid, temp), _mm_andnot_ps(isValid, maxVals)); // any NaN values should be "masked" out

            minVals = _mm_min_ps(maskedMinTemp, minVals);
            maxVals = _mm_max_ps(maskedMaxTemp, maxVals);
        }

        for (size_t i = unrollSize; i < size; i++)
        {
            __m128 temp = _mm_set_ss(values[i]);

            __m128 isValid = _mm_and_ps(_mm_and_ps(_mm_cmpord_ps(temp, temp), _mm_cmpneq_ps(temp, pInfinity)), _mm_cmpneq_ps(temp, nInfinity));

            __m128 maskedMinTemp = _mm_or_ps(_mm_and_ps(isValid, temp), _mm_andnot_ps(isValid, minVals)); // any NaN values should be "masked" out
            __m128 maskedMaxTemp = _mm_or_ps(_mm_and_ps(isValid, temp), _mm_andnot_ps(isValid, maxVals)); // any NaN values should be "masked" out

            minVals = _mm_min_ss(minVals, maskedMinTemp);
            maxVals = _mm_max_ss(maxVals, maskedMaxTemp);
        }

        reduceMinMax(minValue, maxValue, minVals, maxVals);
    }

#endif
}

namespace InternalZGY {
    void MinMaxScan::unsafeScanArray(const float* const values, const size_t size, const size_t stride, float& min, float& max)
    {
#ifdef USE_SSE2_INSTRUCTIONS
        if (stride == 1 && size > 8)
        {
            vectorMinMax(values, size, min, max);
        }
        else
        {
            stridedMinMax(values, size, stride, min, max);
        }
#else
        stridedMinMax(values, size, stride, min, max);
#endif
    }

    void MinMaxScan::scanArray(const float* const values, const size_t size, const size_t stride, float& min, float& max)
    {
#ifdef USE_SSE2_INSTRUCTIONS
        if (stride == 1 && size > 8)
        {
            vectorMinMaxSafe(values, size, min, max);
        }
        else
        {
            stridedMinMaxSafe(values, size, stride, min, max);
        }
#else
        stridedMinMaxSafe(values, size, stride, min, max);
#endif
    }

}
