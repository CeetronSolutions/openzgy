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

#pragma once

#include "declspec.h"
#include <cstddef>

namespace InternalZGY {
	/**
	 * Definition of Min-Max scanning method.
	 *
	 * Array scanning method which uses machine-specific instructions where possible.
	 * This potentially gives large speed-ups, particularly in the case where the size is large
	 * and the stride is 1.
	 */
	class MinMaxScan
	{
	public:
		/// @brief Find the minimum and maximum in the values array. Safe for all floating point values.
		/// @details NaN elements in values are ignored, +/- infinity elements are correctly honoured.
		/// @note Calling with size=0 will result in a range of [std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity()]
		static void scanArray(const float* const values, const size_t size, const size_t stride, float& min, float& max);

		/// @brief Find the minimum and maximum in the values array. Not NaN safe.
		/// @details While this version is faster than scanArray above, it is not safe if your values array potentially contains NaN.
		/// @warning If the values array contains NaN, the output of min and max is undefined.
		/// @note Calling with size=0 will result in a range of [std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity()]
		static void unsafeScanArray(const float* const values, const size_t size, const size_t stride, float& min, float& max);
	};
}
