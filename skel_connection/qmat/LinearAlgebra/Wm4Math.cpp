#include "Wm4Math.h"

namespace Wm4
{
template<> const float Math<float>::EPSILON = FLT_EPSILON;
template<> const float Math<float>::ZERO_TOLERANCE = 1.0e-12;
template<> const float Math<float>::INVERSE_TOLERANCE = 1.0e-12;
template<> const float Math<float>::MAX_REAL = FLT_MAX;
template<> const float Math<float>::PI = (float)(4.0*atan(1.0));
template<> const float Math<float>::TWO_PI = 2.0f*Math<float>::PI;
template<> const float Math<float>::HALF_PI = 0.5f*Math<float>::PI;
template<> const float Math<float>::INV_PI = 1.0f/Math<float>::PI;
template<> const float Math<float>::INV_TWO_PI = 1.0f/Math<float>::TWO_PI;
template<> const float Math<float>::DEG_TO_RAD = Math<float>::PI/180.0f;
template<> const float Math<float>::RAD_TO_DEG = 180.0f/Math<float>::PI;
template<> const float Math<float>::LN_2 = Math<float>::Log(2.0f);
template<> const float Math<float>::LN_10 = Math<float>::Log(10.0f);
template<> const float Math<float>::INV_LN_2 = 1.0f/Math<float>::LN_2;
template<> const float Math<float>::INV_LN_10 = 1.0f/Math<float>::LN_10;

template<> const double Math<double>::EPSILON = DBL_EPSILON;
template<> const double Math<double>::ZERO_TOLERANCE = 1.0e-12;
template<> const double Math<double>::INVERSE_TOLERANCE = 1.0e-12;
template<> const double Math<double>::MAX_REAL = DBL_MAX;
template<> const double Math<double>::PI = 4.0*atan(1.0);
template<> const double Math<double>::TWO_PI = 2.0*Math<double>::PI;
template<> const double Math<double>::HALF_PI = 0.5*Math<double>::PI;
template<> const double Math<double>::INV_PI = 1.0/Math<double>::PI;
template<> const double Math<double>::INV_TWO_PI = 1.0/Math<double>::TWO_PI;
template<> const double Math<double>::DEG_TO_RAD = Math<double>::PI/180.0;
template<> const double Math<double>::RAD_TO_DEG = 180.0/Math<double>::PI;
template<> const double Math<double>::LN_2 = Math<double>::Log(2.0);
template<> const double Math<double>::LN_10 = Math<double>::Log(10.0);
template<> const double Math<double>::INV_LN_2 = 1.0/Math<double>::LN_2;
template<> const double Math<double>::INV_LN_10 = 1.0/Math<double>::LN_10;

}
