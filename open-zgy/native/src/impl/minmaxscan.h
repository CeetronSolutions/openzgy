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

#include "../declspec.h"
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
		/// @details NaN and +/- Inf elements in values are ignored.
		/// @note Calling with size=0 or no valid samples will return with min > max.
		static void scanArray(const float* const values, const size_t size, const size_t stride, float& min, float& max);

		/// @brief Find the minimum and maximum in the values array. Not NaN safe.
		/// @details While this version is faster than scanArray above, it is not safe if your values array potentially contains NaN.
		/// @warning If the values array contains NaN, the output of min and max is undefined.
		/// @note Calling with size=0 will return with min > max.
		static void unsafeScanArray(const float* const values, const size_t size, const size_t stride, float& min, float& max);

		/// @brief Return true if sse2 intrinsics are available and enabled.
		/// @details Use this for logging performance data only.
                static bool has_sse2();

		/// @brief Return true unless callers are requested to not call us at all.
		/// @details Use this for testing performance only.
		/// @details Callers don't need to honor this setting.
                static bool use_sse2();
	};
}
