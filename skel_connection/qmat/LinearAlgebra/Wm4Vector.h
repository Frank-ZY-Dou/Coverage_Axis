// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2007
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
// The license applies to versions 0 through 4 of Wild Magic.
//
// Version: 4.0.3 (2007/03/07)

#ifndef WM4VECTOR_H
#define WM4VECTOR_H

#include "Wm4Math.h"




// vector2


namespace Wm4
{

template <class Real>
class Vector2
{
public:
    // construction
    Vector2 ();  // uninitialized
    Vector2 (Real fX, Real fY);
    Vector2 (const Real* afTuple);
    Vector2 (const Vector2& rkV);

    // coordinate access
    inline operator const Real* () const;
    inline operator Real* ();
    inline Real operator[] (int i) const;
    inline Real& operator[] (int i);
    inline Real X () const;
    inline Real& X ();
    inline Real Y () const;
    inline Real& Y ();

    // assignment
    inline Vector2& operator= (const Vector2& rkV);

    // comparison
    bool operator== (const Vector2& rkV) const;
    bool operator!= (const Vector2& rkV) const;
    bool operator<  (const Vector2& rkV) const;
    bool operator<= (const Vector2& rkV) const;
    bool operator>  (const Vector2& rkV) const;
    bool operator>= (const Vector2& rkV) const;

    // arithmetic operations
    inline Vector2 operator+ (const Vector2& rkV) const;
    inline Vector2 operator- (const Vector2& rkV) const;
    inline Vector2 operator* (Real fScalar) const;
    inline Vector2 operator/ (Real fScalar) const;
    inline Vector2 operator- () const;

    // arithmetic updates
    inline Vector2& operator+= (const Vector2& rkV);
    inline Vector2& operator-= (const Vector2& rkV);
    inline Vector2& operator*= (Real fScalar);
    inline Vector2& operator/= (Real fScalar);

    // vector operations
    inline Real Length () const;
    inline Real SquaredLength () const;
    inline Real Dot (const Vector2& rkV) const;
    inline Real Normalize ();
	inline Vector2 GetNormalized ();

    // returns (y,-x)
    inline Vector2 Perp () const;

    // returns (y,-x)/sqrt(x*x+y*y)
    inline Vector2 UnitPerp () const;

    // returns DotPerp((x,y),(V.x,V.y)) = x*V.y - y*V.x
    inline Real DotPerp (const Vector2& rkV) const;

    // Compute the barycentric coordinates of the point with respect to the
    // triangle <V0,V1,V2>, P = b0*V0 + b1*V1 + b2*V2, where b0 + b1 + b2 = 1.
    void GetBarycentrics (const Vector2& rkV0, const Vector2& rkV1,
        const Vector2& rkV2, Real afBary[3]) const;

    // Gram-Schmidt orthonormalization.  Take linearly independent vectors U
    // and V and compute an orthonormal set (unit length, mutually
    // perpendicular).
    static void Orthonormalize (Vector2& rkU, Vector2& rkV);

    // Input V must be a nonzero vector.  The output is an orthonormal basis
    // {U,V}.  The input V is normalized by this function.  If you know V is
    // already unit length, use U = V.Perp().
    static void GenerateOrthonormalBasis (Vector2& rkU, Vector2& rkV);

    // Compute the extreme values.
    static void ComputeExtremes (int iVQuantity, const Vector2* akPoint,
        Vector2& rkMin, Vector2& rkMax);

    // special vectors
    static const Vector2 ZERO;    // (0,0)
    static const Vector2 UNIT_X;  // (1,0)
    static const Vector2 UNIT_Y;  // (0,1)
    static const Vector2 ONE;     // (1,1)

private:
    // support for comparisons
    int CompareArrays (const Vector2& rkV) const;

    Real m_afTuple[2];
};

// arithmetic operations
template <class Real>
Vector2<Real> operator* (Real fScalar, const Vector2<Real>& rkV);

// debugging output
template <class Real>
std::ostream& operator<< (std::ostream& rkOStr, const Vector2<Real>& rkV);

// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2007
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
// The license applies to versions 0 through 4 of Wild Magic.
//
// Version: 4.0.3 (2007/03/07)

//----------------------------------------------------------------------------
template <class Real>
Vector2<Real>::Vector2 ()
{
    // uninitialized for performance in array construction
}
//----------------------------------------------------------------------------
template <class Real>
Vector2<Real>::Vector2 (Real fX, Real fY)
{
    m_afTuple[0] = fX;
    m_afTuple[1] = fY;
}
//----------------------------------------------------------------------------
template <class Real>
Vector2<Real>::Vector2 (const Real* afTuple)
{
    m_afTuple[0] = afTuple[0];
    m_afTuple[1] = afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
Vector2<Real>::Vector2 (const Vector2& rkV)
{
    m_afTuple[0] = rkV.m_afTuple[0];
    m_afTuple[1] = rkV.m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>::operator const Real* () const
{
    return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>::operator Real* ()
{
    return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::operator[] (int i) const
{
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector2<Real>::operator[] (int i)
{
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::X () const
{
    return m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector2<Real>::X ()
{
    return m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::Y () const
{
    return m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector2<Real>::Y ()
{
    return m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>& Vector2<Real>::operator= (const Vector2& rkV)
{
    m_afTuple[0] = rkV.m_afTuple[0];
    m_afTuple[1] = rkV.m_afTuple[1];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
int Vector2<Real>::CompareArrays (const Vector2& rkV) const
{
    return memcmp(m_afTuple,rkV.m_afTuple,2*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector2<Real>::operator== (const Vector2& rkV) const
{
    return CompareArrays(rkV) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector2<Real>::operator!= (const Vector2& rkV) const
{
    return CompareArrays(rkV) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector2<Real>::operator< (const Vector2& rkV) const
{
    return CompareArrays(rkV) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector2<Real>::operator<= (const Vector2& rkV) const
{
    return CompareArrays(rkV) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector2<Real>::operator> (const Vector2& rkV) const
{
    return CompareArrays(rkV) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector2<Real>::operator>= (const Vector2& rkV) const
{
    return CompareArrays(rkV) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::operator+ (const Vector2& rkV) const
{
    return Vector2(
        m_afTuple[0]+rkV.m_afTuple[0],
        m_afTuple[1]+rkV.m_afTuple[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::operator- (const Vector2& rkV) const
{
    return Vector2(
        m_afTuple[0]-rkV.m_afTuple[0],
        m_afTuple[1]-rkV.m_afTuple[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::operator* (Real fScalar) const
{
    return Vector2(
        fScalar*m_afTuple[0],
        fScalar*m_afTuple[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::operator/ (Real fScalar) const
{
    Vector2 kQuot;

    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        kQuot.m_afTuple[0] = fInvScalar*m_afTuple[0];
        kQuot.m_afTuple[1] = fInvScalar*m_afTuple[1];
    }
    else
    {
        kQuot.m_afTuple[0] = Math<Real>::MAX_REAL;
        kQuot.m_afTuple[1] = Math<Real>::MAX_REAL;
    }

    return kQuot;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::operator- () const
{
    return Vector2(
        -m_afTuple[0],
        -m_afTuple[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> operator* (Real fScalar, const Vector2<Real>& rkV)
{
    return Vector2<Real>(
        fScalar*rkV[0],
        fScalar*rkV[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>& Vector2<Real>::operator+= (const Vector2& rkV)
{
    m_afTuple[0] += rkV.m_afTuple[0];
    m_afTuple[1] += rkV.m_afTuple[1];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>& Vector2<Real>::operator-= (const Vector2& rkV)
{
    m_afTuple[0] -= rkV.m_afTuple[0];
    m_afTuple[1] -= rkV.m_afTuple[1];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>& Vector2<Real>::operator*= (Real fScalar)
{
    m_afTuple[0] *= fScalar;
    m_afTuple[1] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real>& Vector2<Real>::operator/= (Real fScalar)
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        m_afTuple[0] *= fInvScalar;
        m_afTuple[1] *= fInvScalar;
    }
    else
    {
        m_afTuple[0] = Math<Real>::MAX_REAL;
        m_afTuple[1] = Math<Real>::MAX_REAL;
    }

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::Length () const
{
    return Math<Real>::Sqrt(
        m_afTuple[0]*m_afTuple[0] +
        m_afTuple[1]*m_afTuple[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::SquaredLength () const
{
    return
        m_afTuple[0]*m_afTuple[0] +
        m_afTuple[1]*m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::Dot (const Vector2& rkV) const
{
    return
        m_afTuple[0]*rkV.m_afTuple[0] +
        m_afTuple[1]*rkV.m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::Normalize ()
{
    Real fLength = Length();

    if (fLength > Math<Real>::ZERO_TOLERANCE)
    {
        Real fInvLength = ((Real)1.0)/fLength;
        m_afTuple[0] *= fInvLength;
        m_afTuple[1] *= fInvLength;
    }
    else
    {
        fLength = (Real)0.0;
        m_afTuple[0] = (Real)0.0;
        m_afTuple[1] = (Real)0.0;
    }

    return fLength;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::GetNormalized ()
{
    Real fLength = Length();
	Real c[2];

    if (fLength > Math<Real>::ZERO_TOLERANCE)
    {
        Real fInvLength = ((Real)1.0)/fLength;
        c[0] = m_afTuple[0]*fInvLength;
        c[1] = m_afTuple[1]*fInvLength;
    }
    else
    {
        c[0] = (Real)0.0;
        c[1] = (Real)0.0;
    }

    return Vector2(c[0],c[1]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::Perp () const
{
    return Vector2(m_afTuple[1],-m_afTuple[0]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Vector2<Real>::UnitPerp () const
{
    Vector2 kPerp(m_afTuple[1],-m_afTuple[0]);
    kPerp.Normalize();
    return kPerp;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector2<Real>::DotPerp (const Vector2& rkV) const
{
    return m_afTuple[0]*rkV.m_afTuple[1] - m_afTuple[1]*rkV.m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
void Vector2<Real>::GetBarycentrics (const Vector2& rkV0, const Vector2& rkV1,
    const Vector2& rkV2, Real afBary[3]) const
{
    // compute the vectors relative to V2 of the triangle
    Vector2 akDiff[3] =
    {
        rkV0 - rkV2,
        rkV1 - rkV2,
        *this - rkV2
    };

    // If the vertices have large magnitude, the linear system of equations
    // for computing barycentric coordinates can be ill-conditioned.  To avoid
    // this, uniformly scale the triangle edges to be of order 1.  The scaling
    // of all differences does not change the barycentric coordinates.
    Real fMax = (Real)0.0;
    int i;
    for (i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            Real fValue = Math<Real>::FAbs(akDiff[i][j]);
            if (fValue > fMax)
            {
                fMax = fValue;
            }
        }
    }

    // scale down only large data
    if (fMax > (Real)1.0)
    {
        Real fInvMax = ((Real)1.0)/fMax;
        for (i = 0; i < 3; i++)
        {
            akDiff[i] *= fInvMax;
        }
    }

    Real fDet = akDiff[0].DotPerp(akDiff[1]);
    if (Math<Real>::FAbs(fDet) > Math<Real>::ZERO_TOLERANCE)
    {
        Real fInvDet = ((Real)1.0)/fDet;
        afBary[0] = akDiff[2].DotPerp(akDiff[1])*fInvDet;
        afBary[1] = akDiff[0].DotPerp(akDiff[2])*fInvDet;
        afBary[2] = (Real)1.0 - afBary[0] - afBary[1];
    }
    else
    {
        // The triangle is a sliver.  Determine the longest edge and
        // compute barycentric coordinates with respect to that edge.
        Vector2 kE2 = rkV0 - rkV1;
        Real fMaxSqrLength = kE2.SquaredLength();
        int iMaxIndex = 2;
        Real fSqrLength = akDiff[1].SquaredLength();
        if (fSqrLength > fMaxSqrLength)
        {
            iMaxIndex = 1;
            fMaxSqrLength = fSqrLength;
        }
        fSqrLength = akDiff[0].SquaredLength();
        if (fSqrLength > fMaxSqrLength)
        {
            iMaxIndex = 0;
            fMaxSqrLength = fSqrLength;
        }

        if (fMaxSqrLength > Math<Real>::ZERO_TOLERANCE)
        {
            Real fInvSqrLength = ((Real)1.0)/fMaxSqrLength;
            if (iMaxIndex == 0)
            {
                // P-V2 = t(V0-V2)
                afBary[0] = akDiff[2].Dot(akDiff[0])*fInvSqrLength;
                afBary[1] = (Real)0.0;
                afBary[2] = (Real)1.0 - afBary[0];
            }
            else if (iMaxIndex == 1)
            {
                // P-V2 = t(V1-V2)
                afBary[0] = (Real)0.0;
                afBary[1] = akDiff[2].Dot(akDiff[1])*fInvSqrLength;
                afBary[2] = (Real)1.0 - afBary[1];
            }
            else
            {
                // P-V1 = t(V0-V1)
                akDiff[2] = *this - rkV1;
                afBary[0] = akDiff[2].Dot(kE2)*fInvSqrLength;
                afBary[1] = (Real)1.0 - afBary[0];
                afBary[2] = (Real)0.0;
            }
        }
        else
        {
            // triangle is a nearly a point, just return equal weights
            afBary[0] = ((Real)1.0)/(Real)3.0;
            afBary[1] = afBary[0];
            afBary[2] = afBary[0];
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Vector2<Real>::Orthonormalize (Vector2& rkU, Vector2& rkV)
{
    // If the input vectors are v0 and v1, then the Gram-Schmidt
    // orthonormalization produces vectors u0 and u1 as follows,
    //
    //   u0 = v0/|v0|
    //   u1 = (v1-(u0*v1)u0)/|v1-(u0*v1)u0|
    //
    // where |A| indicates length of vector A and A*B indicates dot
    // product of vectors A and B.

    // compute u0
    rkU.Normalize();

    // compute u1
    Real fDot0 = rkU.Dot(rkV); 
    rkV -= rkU*fDot0;
    rkV.Normalize();
}
//----------------------------------------------------------------------------
template <class Real>
void Vector2<Real>::GenerateOrthonormalBasis (Vector2& rkU, Vector2& rkV)
{
    rkV.Normalize();
    rkU = rkV.Perp();
}
//----------------------------------------------------------------------------
template <class Real>
void Vector2<Real>::ComputeExtremes (int iVQuantity, const Vector2* akPoint,
    Vector2& rkMin, Vector2& rkMax)
{
    assert(iVQuantity > 0 && akPoint);

    rkMin = akPoint[0];
    rkMax = rkMin;
    for (int i = 1; i < iVQuantity; i++)
    {
        const Vector2<Real>& rkPoint = akPoint[i];
        for (int j = 0; j < 2; j++)
        {
            if (rkPoint[j] < rkMin[j])
            {
                rkMin[j] = rkPoint[j];
            }
            else if (rkPoint[j] > rkMax[j])
            {
                rkMax[j] = rkPoint[j];
            }
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
std::ostream& operator<< (std::ostream& rkOStr, const Vector2<Real>& rkV)
{
     return rkOStr << rkV.X() << ' ' << rkV.Y();
}
//----------------------------------------------------------------------------


typedef Vector2<float> Vector2f;
typedef Vector2<double> Vector2d;

}







namespace Wm4
{

template <class Real>
class Vector3
{
public:
    // construction
    Vector3 ();  // uninitialized
    Vector3 (Real fX, Real fY, Real fZ);
    Vector3 (const Real* afTuple);
    Vector3 (const Vector3& rkV);

    // coordinate access
    inline operator const Real* () const;
    inline operator Real* ();
    inline Real operator[] (int i) const;
    inline Real& operator[] (int i);
    inline Real X () const;
    inline Real& X ();
    inline Real Y () const;
    inline Real& Y ();
    inline Real Z () const;
    inline Real& Z ();

    // assignment
    inline Vector3& operator= (const Vector3& rkV);

    // comparison
    bool operator== (const Vector3& rkV) const;
    bool operator!= (const Vector3& rkV) const;
    bool operator<  (const Vector3& rkV) const;
    bool operator<= (const Vector3& rkV) const;
    bool operator>  (const Vector3& rkV) const;
    bool operator>= (const Vector3& rkV) const;

    // arithmetic operations
    inline Vector3 operator+ (const Vector3& rkV) const;
    inline Vector3 operator- (const Vector3& rkV) const;
    inline Vector3 operator* (Real fScalar) const;
    inline Vector3 operator/ (Real fScalar) const;
    inline Vector3 operator- () const;

    // arithmetic updates
    inline Vector3& operator+= (const Vector3& rkV);
    inline Vector3& operator-= (const Vector3& rkV);
    inline Vector3& operator*= (Real fScalar);
    inline Vector3& operator/= (Real fScalar);

    // vector operations
    inline Real Length () const;
    inline Real SquaredLength () const;
    inline Real Dot (const Vector3& rkV) const;
    inline Real Normalize ();

    // The cross products are computed using the right-handed rule.  Be aware
    // that some graphics APIs use a left-handed rule.  If you have to compute
    // a cross product with these functions and send the result to the API
    // that expects left-handed, you will need to change sign on the vector
    // (replace each component value c by -c).
    inline Vector3 Cross (const Vector3& rkV) const;
    inline Vector3 UnitCross (const Vector3& rkV) const;

    // Compute the barycentric coordinates of the point with respect to the
    // tetrahedron <V0,V1,V2,V3>, P = b0*V0 + b1*V1 + b2*V2 + b3*V3, where
    // b0 + b1 + b2 + b3 = 1.
    void GetBarycentrics (const Vector3& rkV0, const Vector3& rkV1,
        const Vector3& rkV2, const Vector3& rkV3, Real afBary[4]) const;

    // Gram-Schmidt orthonormalization.  Take linearly independent vectors
    // U, V, and W and compute an orthonormal set (unit length, mutually
    // perpendicular).
    static void Orthonormalize (Vector3& rkU, Vector3& rkV, Vector3& rkW);
    static void Orthonormalize (Vector3* akV);

    // Input W must be a nonzero vector. The output is an orthonormal basis
    // {U,V,W}.  The input W is normalized by this function.  If you know
    // W is already unit length, use GenerateComplementBasis to compute U
    // and V.
    static void GenerateOrthonormalBasis (Vector3& rkU, Vector3& rkV,
        Vector3& rkW);

    // Input W must be a unit-length vector.  The output vectors {U,V} are
    // unit length and mutually perpendicular, and {U,V,W} is an orthonormal
    // basis.
    static void GenerateComplementBasis (Vector3& rkU, Vector3& rkV,
        const Vector3& rkW);

    // Compute the extreme values.
    static void ComputeExtremes (int iVQuantity, const Vector3* akPoint,
        Vector3& rkMin, Vector3& rkMax);

    // special vectors
    static const Vector3 ZERO;    // (0,0,0)
    static const Vector3 UNIT_X;  // (1,0,0)
    static const Vector3 UNIT_Y;  // (0,1,0)
    static const Vector3 UNIT_Z;  // (0,0,1)
    static const Vector3 ONE;     // (1,1,1)

private:
    // support for comparisons
    int CompareArrays (const Vector3& rkV) const;

    Real m_afTuple[3];
};

// arithmetic operations
template <class Real>
Vector3<Real> operator* (Real fScalar, const Vector3<Real>& rkV);

// debugging output
template <class Real>
std::ostream& operator<< (std::ostream& rkOStr, const Vector3<Real>& rkV);

// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2007
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
//
// Version: 4.0.3 (2007/03/07)

//----------------------------------------------------------------------------
template <class Real>
Vector3<Real>::Vector3 ()
{
    // uninitialized for performance in array construction
}
//----------------------------------------------------------------------------
template <class Real>
Vector3<Real>::Vector3 (Real fX, Real fY, Real fZ)
{
    m_afTuple[0] = fX;
    m_afTuple[1] = fY;
    m_afTuple[2] = fZ;
}
//----------------------------------------------------------------------------
template <class Real>
Vector3<Real>::Vector3 (const Real* afTuple)
{
    m_afTuple[0] = afTuple[0];
    m_afTuple[1] = afTuple[1];
    m_afTuple[2] = afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
Vector3<Real>::Vector3 (const Vector3& rkV)
{
    m_afTuple[0] = rkV.m_afTuple[0];
    m_afTuple[1] = rkV.m_afTuple[1];
    m_afTuple[2] = rkV.m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>::operator const Real* () const
{
    return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>::operator Real* ()
{
    return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::operator[] (int i) const
{
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector3<Real>::operator[] (int i)
{
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::X () const
{
    return m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector3<Real>::X ()
{
    return m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::Y () const
{
    return m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector3<Real>::Y ()
{
    return m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::Z () const
{
    return m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector3<Real>::Z ()
{
    return m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>& Vector3<Real>::operator= (const Vector3& rkV)
{
    m_afTuple[0] = rkV.m_afTuple[0];
    m_afTuple[1] = rkV.m_afTuple[1];
    m_afTuple[2] = rkV.m_afTuple[2];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
int Vector3<Real>::CompareArrays (const Vector3& rkV) const
{
    return memcmp(m_afTuple,rkV.m_afTuple,3*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector3<Real>::operator== (const Vector3& rkV) const
{
    return CompareArrays(rkV) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector3<Real>::operator!= (const Vector3& rkV) const
{
    return CompareArrays(rkV) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector3<Real>::operator< (const Vector3& rkV) const
{
    return CompareArrays(rkV) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector3<Real>::operator<= (const Vector3& rkV) const
{
    return CompareArrays(rkV) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector3<Real>::operator> (const Vector3& rkV) const
{
    return CompareArrays(rkV) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector3<Real>::operator>= (const Vector3& rkV) const
{
    return CompareArrays(rkV) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::operator+ (const Vector3& rkV) const
{
    return Vector3(
        m_afTuple[0]+rkV.m_afTuple[0],
        m_afTuple[1]+rkV.m_afTuple[1],
        m_afTuple[2]+rkV.m_afTuple[2]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::operator- (const Vector3& rkV) const
{
    return Vector3(
        m_afTuple[0]-rkV.m_afTuple[0],
        m_afTuple[1]-rkV.m_afTuple[1],
        m_afTuple[2]-rkV.m_afTuple[2]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::operator* (Real fScalar) const
{
    return Vector3(
        fScalar*m_afTuple[0],
        fScalar*m_afTuple[1],
        fScalar*m_afTuple[2]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::operator/ (Real fScalar) const
{
    Vector3 kQuot;

    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        kQuot.m_afTuple[0] = fInvScalar*m_afTuple[0];
        kQuot.m_afTuple[1] = fInvScalar*m_afTuple[1];
        kQuot.m_afTuple[2] = fInvScalar*m_afTuple[2];
    }
    else
    {
        kQuot.m_afTuple[0] = Math<Real>::MAX_REAL;
        kQuot.m_afTuple[1] = Math<Real>::MAX_REAL;
        kQuot.m_afTuple[2] = Math<Real>::MAX_REAL;
    }

    return kQuot;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::operator- () const
{
    return Vector3(
        -m_afTuple[0],
        -m_afTuple[1],
        -m_afTuple[2]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> operator* (Real fScalar, const Vector3<Real>& rkV)
{
    return Vector3<Real>(
        fScalar*rkV[0],
        fScalar*rkV[1],
        fScalar*rkV[2]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>& Vector3<Real>::operator+= (const Vector3& rkV)
{
    m_afTuple[0] += rkV.m_afTuple[0];
    m_afTuple[1] += rkV.m_afTuple[1];
    m_afTuple[2] += rkV.m_afTuple[2];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>& Vector3<Real>::operator-= (const Vector3& rkV)
{
    m_afTuple[0] -= rkV.m_afTuple[0];
    m_afTuple[1] -= rkV.m_afTuple[1];
    m_afTuple[2] -= rkV.m_afTuple[2];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>& Vector3<Real>::operator*= (Real fScalar)
{
    m_afTuple[0] *= fScalar;
    m_afTuple[1] *= fScalar;
    m_afTuple[2] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real>& Vector3<Real>::operator/= (Real fScalar)
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        m_afTuple[0] *= fInvScalar;
        m_afTuple[1] *= fInvScalar;
        m_afTuple[2] *= fInvScalar;
    }
    else
    {
        m_afTuple[0] = Math<Real>::MAX_REAL;
        m_afTuple[1] = Math<Real>::MAX_REAL;
        m_afTuple[2] = Math<Real>::MAX_REAL;
    }

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::Length () const
{
    return Math<Real>::Sqrt(
        m_afTuple[0]*m_afTuple[0] +
        m_afTuple[1]*m_afTuple[1] +
        m_afTuple[2]*m_afTuple[2]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::SquaredLength () const
{
    return
        m_afTuple[0]*m_afTuple[0] +
        m_afTuple[1]*m_afTuple[1] +
        m_afTuple[2]*m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::Dot (const Vector3& rkV) const
{
    return
        m_afTuple[0]*rkV.m_afTuple[0] +
        m_afTuple[1]*rkV.m_afTuple[1] +
        m_afTuple[2]*rkV.m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector3<Real>::Normalize ()
{
    Real fLength = Length();

    if (fLength > Math<Real>::ZERO_TOLERANCE)
    {
        Real fInvLength = ((Real)1.0)/fLength;
        m_afTuple[0] *= fInvLength;
        m_afTuple[1] *= fInvLength;
        m_afTuple[2] *= fInvLength;
    }
    else
    {
        fLength = (Real)0.0;
        m_afTuple[0] = (Real)0.0;
        m_afTuple[1] = (Real)0.0;
        m_afTuple[2] = (Real)0.0;
    }

    return fLength;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::Cross (const Vector3& rkV) const
{
    return Vector3(
        m_afTuple[1]*rkV.m_afTuple[2] - m_afTuple[2]*rkV.m_afTuple[1],
        m_afTuple[2]*rkV.m_afTuple[0] - m_afTuple[0]*rkV.m_afTuple[2],
        m_afTuple[0]*rkV.m_afTuple[1] - m_afTuple[1]*rkV.m_afTuple[0]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Vector3<Real>::UnitCross (const Vector3& rkV) const
{
    Vector3 kCross(
        m_afTuple[1]*rkV.m_afTuple[2] - m_afTuple[2]*rkV.m_afTuple[1],
        m_afTuple[2]*rkV.m_afTuple[0] - m_afTuple[0]*rkV.m_afTuple[2],
        m_afTuple[0]*rkV.m_afTuple[1] - m_afTuple[1]*rkV.m_afTuple[0]);
    kCross.Normalize();
    return kCross;
}
//----------------------------------------------------------------------------
template <class Real>
void Vector3<Real>::GetBarycentrics (const Vector3<Real>& rkV0,
    const Vector3<Real>& rkV1, const Vector3<Real>& rkV2,
    const Vector3<Real>& rkV3, Real afBary[4]) const
{
    // compute the vectors relative to V3 of the tetrahedron
    Vector3<Real> akDiff[4] =
    {
        rkV0 - rkV3,
        rkV1 - rkV3,
        rkV2 - rkV3,
        *this - rkV3
    };

    // If the vertices have large magnitude, the linear system of
    // equations for computing barycentric coordinates can be
    // ill-conditioned.  To avoid this, uniformly scale the tetrahedron
    // edges to be of order 1.  The scaling of all differences does not
    // change the barycentric coordinates.
    Real fMax = (Real)0.0;
    int i;
    for (i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Real fValue = Math<Real>::FAbs(akDiff[i][j]);
            if (fValue > fMax)
            {
                fMax = fValue;
            }
        }
    }

    // scale down only large data
    if (fMax > (Real)1.0)
    {
        Real fInvMax = ((Real)1.0)/fMax;
        for (i = 0; i < 4; i++)
        {
            akDiff[i] *= fInvMax;
        }
    }

    Real fDet = akDiff[0].Dot(akDiff[1].Cross(akDiff[2]));
    Vector3<Real> kE1cE2 = akDiff[1].Cross(akDiff[2]);
    Vector3<Real> kE2cE0 = akDiff[2].Cross(akDiff[0]);
    Vector3<Real> kE0cE1 = akDiff[0].Cross(akDiff[1]);
    if (Math<Real>::FAbs(fDet) > Math<Real>::ZERO_TOLERANCE)
    {
        Real fInvDet = ((Real)1.0)/fDet;
        afBary[0] = akDiff[3].Dot(kE1cE2)*fInvDet;
        afBary[1] = akDiff[3].Dot(kE2cE0)*fInvDet;
        afBary[2] = akDiff[3].Dot(kE0cE1)*fInvDet;
        afBary[3] = (Real)1.0 - afBary[0] - afBary[1] - afBary[2];
    }
    else
    {
        // The tetrahedron is potentially flat.  Determine the face of
        // maximum area and compute barycentric coordinates with respect
        // to that face.
        Vector3<Real> kE02 = rkV0 - rkV2;
        Vector3<Real> kE12 = rkV1 - rkV2;
        Vector3<Real> kE02cE12 = kE02.Cross(kE12);
        Real fMaxSqrArea = kE02cE12.SquaredLength();
        int iMaxIndex = 3;
        Real fSqrArea = kE0cE1.SquaredLength();
        if (fSqrArea > fMaxSqrArea)
        {
            iMaxIndex = 0;
            fMaxSqrArea = fSqrArea;
        }
        fSqrArea = kE1cE2.SquaredLength();
        if (fSqrArea > fMaxSqrArea)
        {
            iMaxIndex = 1;
            fMaxSqrArea = fSqrArea;
        }
        fSqrArea = kE2cE0.SquaredLength();
        if (fSqrArea > fMaxSqrArea)
        {
            iMaxIndex = 2;
            fMaxSqrArea = fSqrArea;
        }

        if (fMaxSqrArea > Math<Real>::ZERO_TOLERANCE)
        {
            Real fInvSqrArea = ((Real)1.0)/fMaxSqrArea;
            Vector3<Real> kTmp;
            if (iMaxIndex == 0)
            {
                kTmp = akDiff[3].Cross(akDiff[1]);
                afBary[0] = kE0cE1.Dot(kTmp)*fInvSqrArea;
                kTmp = akDiff[0].Cross(akDiff[3]);
                afBary[1] = kE0cE1.Dot(kTmp)*fInvSqrArea;
                afBary[2] = (Real)0.0;
                afBary[3] = (Real)1.0 - afBary[0] - afBary[1];
            }
            else if (iMaxIndex == 1)
            {
                afBary[0] = (Real)0.0;
                kTmp = akDiff[3].Cross(akDiff[2]);
                afBary[1] = kE1cE2.Dot(kTmp)*fInvSqrArea;
                kTmp = akDiff[1].Cross(akDiff[3]);
                afBary[2] = kE1cE2.Dot(kTmp)*fInvSqrArea;
                afBary[3] = (Real)1.0 - afBary[1] - afBary[2];
            }
            else if (iMaxIndex == 2)
            {
                kTmp = akDiff[2].Cross(akDiff[3]);
                afBary[0] = kE2cE0.Dot(kTmp)*fInvSqrArea;
                afBary[1] = (Real)0.0;
                kTmp = akDiff[3].Cross(akDiff[0]);
                afBary[2] = kE2cE0.Dot(kTmp)*fInvSqrArea;
                afBary[3] = (Real)1.0 - afBary[0] - afBary[2];
            }
            else
            {
                akDiff[3] = *this - rkV2;
                kTmp = akDiff[3].Cross(kE12);
                afBary[0] = kE02cE12.Dot(kTmp)*fInvSqrArea;
                kTmp = kE02.Cross(akDiff[3]);
                afBary[1] = kE02cE12.Dot(kTmp)*fInvSqrArea;
                afBary[2] = (Real)1.0 - afBary[0] - afBary[1];
                afBary[3] = (Real)0.0;
            }
        }
        else
        {
            // The tetrahedron is potentially a sliver.  Determine the edge of
            // maximum length and compute barycentric coordinates with respect
            // to that edge.
            Real fMaxSqrLength = akDiff[0].SquaredLength();
            iMaxIndex = 0;  // <V0,V3>
            Real fSqrLength = akDiff[1].SquaredLength();
            if (fSqrLength > fMaxSqrLength)
            {
                iMaxIndex = 1;  // <V1,V3>
                fMaxSqrLength = fSqrLength;
            }
            fSqrLength = akDiff[2].SquaredLength();
            if (fSqrLength > fMaxSqrLength)
            {
                iMaxIndex = 2;  // <V2,V3>
                fMaxSqrLength = fSqrLength;
            }
            fSqrLength = kE02.SquaredLength();
            if (fSqrLength > fMaxSqrLength)
            {
                iMaxIndex = 3;  // <V0,V2>
                fMaxSqrLength = fSqrLength;
            }
            fSqrLength = kE12.SquaredLength();
            if (fSqrLength > fMaxSqrLength)
            {
                iMaxIndex = 4;  // <V1,V2>
                fMaxSqrLength = fSqrLength;
            }
            Vector3<Real> kE01 = rkV0 - rkV1;
            fSqrLength = kE01.SquaredLength();
            if (fSqrLength > fMaxSqrLength)
            {
                iMaxIndex = 5;  // <V0,V1>
                fMaxSqrLength = fSqrLength;
            }

            if (fMaxSqrLength > Math<Real>::ZERO_TOLERANCE)
            {
                Real fInvSqrLength = ((Real)1.0)/fMaxSqrLength;
                if (iMaxIndex == 0)
                {
                    // P-V3 = t*(V0-V3)
                    afBary[0] = akDiff[3].Dot(akDiff[0])*fInvSqrLength;
                    afBary[1] = (Real)0.0;
                    afBary[2] = (Real)0.0;
                    afBary[3] = (Real)1.0 - afBary[0];
                }
                else if (iMaxIndex == 1)
                {
                    // P-V3 = t*(V1-V3)
                    afBary[0] = (Real)0.0;
                    afBary[1] = akDiff[3].Dot(akDiff[1])*fInvSqrLength;
                    afBary[2] = (Real)0.0;
                    afBary[3] = (Real)1.0 - afBary[1];
                }
                else if (iMaxIndex == 2)
                {
                    // P-V3 = t*(V2-V3)
                    afBary[0] = (Real)0.0;
                    afBary[1] = (Real)0.0;
                    afBary[2] = akDiff[3].Dot(akDiff[2])*fInvSqrLength;
                    afBary[3] = (Real)1.0 - afBary[2];
                }
                else if (iMaxIndex == 3)
                {
                    // P-V2 = t*(V0-V2)
                    akDiff[3] = *this - rkV2;
                    afBary[0] = akDiff[3].Dot(kE02)*fInvSqrLength;
                    afBary[1] = (Real)0.0;
                    afBary[2] = (Real)1.0 - afBary[0];
                    afBary[3] = (Real)0.0;
                }
                else if (iMaxIndex == 4)
                {
                    // P-V2 = t*(V1-V2)
                    akDiff[3] = *this - rkV2;
                    afBary[0] = (Real)0.0;
                    afBary[1] = akDiff[3].Dot(kE12)*fInvSqrLength;
                    afBary[2] = (Real)1.0 - afBary[1];
                    afBary[3] = (Real)0.0;
                }
                else
                {
                    // P-V1 = t*(V0-V1)
                    akDiff[3] = *this - rkV1;
                    afBary[0] = akDiff[3].Dot(kE01)*fInvSqrLength;
                    afBary[1] = (Real)1.0 - afBary[0];
                    afBary[2] = (Real)0.0;
                    afBary[3] = (Real)0.0;
                }
            }
            else
            {
                // tetrahedron is a nearly a point, just return equal weights
                afBary[0] = (Real)0.25;
                afBary[1] = afBary[0];
                afBary[2] = afBary[0];
                afBary[3] = afBary[0];
            }
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Vector3<Real>::Orthonormalize (Vector3& rkU, Vector3& rkV, Vector3& rkW)
{
    // If the input vectors are v0, v1, and v2, then the Gram-Schmidt
    // orthonormalization produces vectors u0, u1, and u2 as follows,
    //
    //   u0 = v0/|v0|
    //   u1 = (v1-(u0*v1)u0)/|v1-(u0*v1)u0|
    //   u2 = (v2-(u0*v2)u0-(u1*v2)u1)/|v2-(u0*v2)u0-(u1*v2)u1|
    //
    // where |A| indicates length of vector A and A*B indicates dot
    // product of vectors A and B.

    // compute u0
    rkU.Normalize();

    // compute u1
    Real fDot0 = rkU.Dot(rkV); 
    rkV -= fDot0*rkU;
    rkV.Normalize();

    // compute u2
    Real fDot1 = rkV.Dot(rkW);
    fDot0 = rkU.Dot(rkW);
    rkW -= fDot0*rkU + fDot1*rkV;
    rkW.Normalize();
}
//----------------------------------------------------------------------------
template <class Real>
void Vector3<Real>::Orthonormalize (Vector3* akV)
{
    Orthonormalize(akV[0],akV[1],akV[2]);
}
//----------------------------------------------------------------------------
template <class Real>
void Vector3<Real>::GenerateOrthonormalBasis (Vector3& rkU, Vector3& rkV,
    Vector3& rkW)
{
    rkW.Normalize();
    GenerateComplementBasis(rkU,rkV,rkW);
}
//----------------------------------------------------------------------------
template <class Real>
void Vector3<Real>::GenerateComplementBasis (Vector3& rkU, Vector3& rkV,
    const Vector3& rkW)
{
    Real fInvLength;

    if (Math<Real>::FAbs(rkW.m_afTuple[0]) >=
        Math<Real>::FAbs(rkW.m_afTuple[1]) )
    {
        // W.x or W.z is the largest magnitude component, swap them
        fInvLength = Math<Real>::InvSqrt(rkW.m_afTuple[0]*rkW.m_afTuple[0] +
            rkW.m_afTuple[2]*rkW.m_afTuple[2]);
        rkU.m_afTuple[0] = -rkW.m_afTuple[2]*fInvLength;
        rkU.m_afTuple[1] = (Real)0.0;
        rkU.m_afTuple[2] = +rkW.m_afTuple[0]*fInvLength;
        rkV.m_afTuple[0] = rkW.m_afTuple[1]*rkU.m_afTuple[2];
        rkV.m_afTuple[1] = rkW.m_afTuple[2]*rkU.m_afTuple[0] -
            rkW.m_afTuple[0]*rkU.m_afTuple[2];
        rkV.m_afTuple[2] = -rkW.m_afTuple[1]*rkU.m_afTuple[0];
    }
    else
    {
        // W.y or W.z is the largest magnitude component, swap them
        fInvLength = Math<Real>::InvSqrt(rkW.m_afTuple[1]*rkW.m_afTuple[1] +
            rkW.m_afTuple[2]*rkW.m_afTuple[2]);
        rkU.m_afTuple[0] = (Real)0.0;
        rkU.m_afTuple[1] = +rkW.m_afTuple[2]*fInvLength;
        rkU.m_afTuple[2] = -rkW.m_afTuple[1]*fInvLength;
        rkV.m_afTuple[0] = rkW.m_afTuple[1]*rkU.m_afTuple[2] -
            rkW.m_afTuple[2]*rkU.m_afTuple[1];
        rkV.m_afTuple[1] = -rkW.m_afTuple[0]*rkU.m_afTuple[2];
        rkV.m_afTuple[2] = rkW.m_afTuple[0]*rkU.m_afTuple[1];
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Vector3<Real>::ComputeExtremes (int iVQuantity, const Vector3* akPoint,
    Vector3& rkMin, Vector3& rkMax)
{
    assert(iVQuantity > 0 && akPoint);

    rkMin = akPoint[0];
    rkMax = rkMin;
    for (int i = 1; i < iVQuantity; i++)
    {
        const Vector3<Real>& rkPoint = akPoint[i];
        for (int j = 0; j < 3; j++)
        {
            if (rkPoint[j] < rkMin[j])
            {
                rkMin[j] = rkPoint[j];
            }
            else if (rkPoint[j] > rkMax[j])
            {
                rkMax[j] = rkPoint[j];
            }
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
std::ostream& operator<< (std::ostream& rkOStr, const Vector3<Real>& rkV)
{
     return rkOStr << rkV.X() << ' ' << rkV.Y() << ' ' << rkV.Z();
}
//----------------------------------------------------------------------------


typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;
typedef Vector3<int> Vector3i;

}







namespace Wm4
{

template <class Real>
class Vector4
{
public:
    // construction
    Vector4 ();  // uninitialized
    Vector4 (Real fX, Real fY, Real fZ, Real fW);
    Vector4 (const Real* afTuple);
    Vector4 (const Vector4& rkV);

    // coordinate access
    inline operator const Real* () const;
    inline operator Real* ();
    inline Real operator[] (int i) const;
    inline Real& operator[] (int i);
    inline Real X () const;
    inline Real& X ();
    inline Real Y () const;
    inline Real& Y ();
    inline Real Z () const;
    inline Real& Z ();
    inline Real W () const;
    inline Real& W ();

    // assignment
    inline Vector4& operator= (const Vector4& rkV);

    // comparison
    bool operator== (const Vector4& rkV) const;
    bool operator!= (const Vector4& rkV) const;
    bool operator<  (const Vector4& rkV) const;
    bool operator<= (const Vector4& rkV) const;
    bool operator>  (const Vector4& rkV) const;
    bool operator>= (const Vector4& rkV) const;

    // arithmetic operations
    inline Vector4 operator+ (const Vector4& rkV) const;
    inline Vector4 operator- (const Vector4& rkV) const;
    inline Vector4 operator* (Real fScalar) const;
    inline Vector4 operator/ (Real fScalar) const;
    inline Vector4 operator- () const;

    // arithmetic updates
    inline Vector4& operator+= (const Vector4& rkV);
    inline Vector4& operator-= (const Vector4& rkV);
    inline Vector4& operator*= (Real fScalar);
    inline Vector4& operator/= (Real fScalar);

    // vector operations
    inline Real Length () const;
    inline Real SquaredLength () const;
    inline Real Dot (const Vector4& rkV) const;
    inline Real Normalize ();

    // special vectors
    static const Vector4 ZERO;
    static const Vector4 UNIT_X;  // (1,0,0,0)
    static const Vector4 UNIT_Y;  // (0,1,0,0)
    static const Vector4 UNIT_Z;  // (0,0,1,0)
    static const Vector4 UNIT_W;  // (0,0,0,1)
    static const Vector4 ONE;     // (1,1,1,1)

private:
    // support for comparisons
    int CompareArrays (const Vector4& rkV) const;

    Real m_afTuple[4];
};

// arithmetic operations
template <class Real>
Vector4<Real> operator* (Real fScalar, const Vector4<Real>& rkV);

// debugging output
template <class Real>
std::ostream& operator<< (std::ostream& rkOStr, const Vector4<Real>& rkV);

// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2007
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
// The license applies to versions 0 through 4 of Wild Magic.
//
// Version: 4.0.2 (2007/03/07)

//----------------------------------------------------------------------------
template <class Real>
Vector4<Real>::Vector4 ()
{
	m_afTuple[0] = (Real)0.0;
	m_afTuple[1] = (Real)0.0;
	m_afTuple[2] = (Real)0.0;
	m_afTuple[3] = (Real)0.0;
    // uninitialized for performance in array construction
}
//----------------------------------------------------------------------------
template <class Real>
Vector4<Real>::Vector4 (Real fX, Real fY, Real fZ, Real fW)
{
    m_afTuple[0] = fX;
    m_afTuple[1] = fY;
    m_afTuple[2] = fZ;
    m_afTuple[3] = fW;
}
//----------------------------------------------------------------------------
template <class Real>
Vector4<Real>::Vector4 (const Real* afTuple)
{
    m_afTuple[0] = afTuple[0];
    m_afTuple[1] = afTuple[1];
    m_afTuple[2] = afTuple[2];
    m_afTuple[3] = afTuple[3];
}
//----------------------------------------------------------------------------
template <class Real>
Vector4<Real>::Vector4 (const Vector4& rkV)
{
    m_afTuple[0] = rkV.m_afTuple[0];
    m_afTuple[1] = rkV.m_afTuple[1];
    m_afTuple[2] = rkV.m_afTuple[2];
    m_afTuple[3] = rkV.m_afTuple[3];
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>::operator const Real* () const
{
    return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>::operator Real* ()
{
    return m_afTuple;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::operator[] (int i) const
{
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector4<Real>::operator[] (int i)
{
    return m_afTuple[i];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::X () const
{
    return m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector4<Real>::X ()
{
    return m_afTuple[0];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::Y () const
{
    return m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector4<Real>::Y ()
{
    return m_afTuple[1];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::Z () const
{
    return m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector4<Real>::Z ()
{
    return m_afTuple[2];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::W () const
{
    return m_afTuple[3];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Vector4<Real>::W ()
{
    return m_afTuple[3];
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>& Vector4<Real>::operator= (const Vector4& rkV)
{
    m_afTuple[0] = rkV.m_afTuple[0];
    m_afTuple[1] = rkV.m_afTuple[1];
    m_afTuple[2] = rkV.m_afTuple[2];
    m_afTuple[3] = rkV.m_afTuple[3];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
int Vector4<Real>::CompareArrays (const Vector4& rkV) const
{
    return memcmp(m_afTuple,rkV.m_afTuple,4*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector4<Real>::operator== (const Vector4& rkV) const
{
    return CompareArrays(rkV) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector4<Real>::operator!= (const Vector4& rkV) const
{
    return CompareArrays(rkV) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector4<Real>::operator< (const Vector4& rkV) const
{
    return CompareArrays(rkV) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector4<Real>::operator<= (const Vector4& rkV) const
{
    return CompareArrays(rkV) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector4<Real>::operator> (const Vector4& rkV) const
{
    return CompareArrays(rkV) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Vector4<Real>::operator>= (const Vector4& rkV) const
{
    return CompareArrays(rkV) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> Vector4<Real>::operator+ (const Vector4& rkV) const
{
    return Vector4(
        m_afTuple[0]+rkV.m_afTuple[0],
        m_afTuple[1]+rkV.m_afTuple[1],
        m_afTuple[2]+rkV.m_afTuple[2],
        m_afTuple[3]+rkV.m_afTuple[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> Vector4<Real>::operator- (const Vector4& rkV) const
{
    return Vector4(
        m_afTuple[0]-rkV.m_afTuple[0],
        m_afTuple[1]-rkV.m_afTuple[1],
        m_afTuple[2]-rkV.m_afTuple[2],
        m_afTuple[3]-rkV.m_afTuple[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> Vector4<Real>::operator* (Real fScalar) const
{
    return Vector4(
        fScalar*m_afTuple[0],
        fScalar*m_afTuple[1],
        fScalar*m_afTuple[2],
        fScalar*m_afTuple[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> Vector4<Real>::operator/ (Real fScalar) const
{
    Vector4 kQuot;

    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        kQuot.m_afTuple[0] = fInvScalar*m_afTuple[0];
        kQuot.m_afTuple[1] = fInvScalar*m_afTuple[1];
        kQuot.m_afTuple[2] = fInvScalar*m_afTuple[2];
        kQuot.m_afTuple[3] = fInvScalar*m_afTuple[3];
    }
    else
    {
        kQuot.m_afTuple[0] = Math<Real>::MAX_REAL;
        kQuot.m_afTuple[1] = Math<Real>::MAX_REAL;
        kQuot.m_afTuple[2] = Math<Real>::MAX_REAL;
        kQuot.m_afTuple[3] = Math<Real>::MAX_REAL;
    }

    return kQuot;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> Vector4<Real>::operator- () const
{
    return Vector4(
        -m_afTuple[0],
        -m_afTuple[1],
        -m_afTuple[2],
        -m_afTuple[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> operator* (Real fScalar, const Vector4<Real>& rkV)
{
    return Vector4<Real>(
        fScalar*rkV[0],
        fScalar*rkV[1],
        fScalar*rkV[2],
        fScalar*rkV[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>& Vector4<Real>::operator+= (const Vector4& rkV)
{
    m_afTuple[0] += rkV.m_afTuple[0];
    m_afTuple[1] += rkV.m_afTuple[1];
    m_afTuple[2] += rkV.m_afTuple[2];
    m_afTuple[3] += rkV.m_afTuple[3];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>& Vector4<Real>::operator-= (const Vector4& rkV)
{
    m_afTuple[0] -= rkV.m_afTuple[0];
    m_afTuple[1] -= rkV.m_afTuple[1];
    m_afTuple[2] -= rkV.m_afTuple[2];
    m_afTuple[3] -= rkV.m_afTuple[3];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>& Vector4<Real>::operator*= (Real fScalar)
{
    m_afTuple[0] *= fScalar;
    m_afTuple[1] *= fScalar;
    m_afTuple[2] *= fScalar;
    m_afTuple[3] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real>& Vector4<Real>::operator/= (Real fScalar)
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        m_afTuple[0] *= fInvScalar;
        m_afTuple[1] *= fInvScalar;
        m_afTuple[2] *= fInvScalar;
        m_afTuple[3] *= fInvScalar;
    }
    else
    {
        m_afTuple[0] = Math<Real>::MAX_REAL;
        m_afTuple[1] = Math<Real>::MAX_REAL;
        m_afTuple[2] = Math<Real>::MAX_REAL;
        m_afTuple[3] = Math<Real>::MAX_REAL;
    }

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::Length () const
{
    return Math<Real>::Sqrt(
        m_afTuple[0]*m_afTuple[0] +
        m_afTuple[1]*m_afTuple[1] +
        m_afTuple[2]*m_afTuple[2] +
        m_afTuple[3]*m_afTuple[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::SquaredLength () const
{
    return
        m_afTuple[0]*m_afTuple[0] +
        m_afTuple[1]*m_afTuple[1] +
        m_afTuple[2]*m_afTuple[2] +
        m_afTuple[3]*m_afTuple[3];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::Dot (const Vector4& rkV) const
{
    return
        m_afTuple[0]*rkV.m_afTuple[0] +
        m_afTuple[1]*rkV.m_afTuple[1] +
        m_afTuple[2]*rkV.m_afTuple[2] +
        m_afTuple[3]*rkV.m_afTuple[3];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Vector4<Real>::Normalize ()
{
    Real fLength = Length();

    if (fLength > Math<Real>::ZERO_TOLERANCE)
    {
        Real fInvLength = ((Real)1.0)/fLength;
        m_afTuple[0] *= fInvLength;
        m_afTuple[1] *= fInvLength;
        m_afTuple[2] *= fInvLength;
        m_afTuple[3] *= fInvLength;
    }
    else
    {
        fLength = (Real)0.0;
        m_afTuple[0] = (Real)0.0;
        m_afTuple[1] = (Real)0.0;
        m_afTuple[2] = (Real)0.0;
        m_afTuple[3] = (Real)0.0;
    }

    return fLength;
}
//----------------------------------------------------------------------------
template <class Real>
std::ostream& operator<< (std::ostream& rkOStr, const Vector4<Real>& rkV)
{
     return rkOStr << rkV.X() << ' ' << rkV.Y() << ' ' << rkV.Z()
         << ' ' << rkV.W();
}
//----------------------------------------------------------------------------


typedef Vector4<float> Vector4f;
typedef Vector4<double> Vector4d;
typedef Vector4<int> Vector4i;

}






namespace Wm4
{

	template <class Real>
	class Vector5
	{
	public:
		// construction
		Vector5 ();  // uninitialized
		Vector5 (Real fX, Real fY, Real fZ, Real fW, Real fU);
		Vector5 (const Real* afTuple);
		Vector5 (const Vector5& rkV);

		// coordinate access
		inline operator const Real* () const;
		inline operator Real* ();
		inline Real operator[] (int i) const;
		inline Real& operator[] (int i);
		inline Real X () const;
		inline Real& X ();
		inline Real Y () const;
		inline Real& Y ();
		inline Real Z () const;
		inline Real& Z ();
		inline Real W () const;
		inline Real& W ();
		inline Real U () const;
		inline Real& U ();

		// assignment
		inline Vector5& operator= (const Vector5& rkV);

		// comparison
		bool operator== (const Vector5& rkV) const;
		bool operator!= (const Vector5& rkV) const;
		bool operator<  (const Vector5& rkV) const;
		bool operator<= (const Vector5& rkV) const;
		bool operator>  (const Vector5& rkV) const;
		bool operator>= (const Vector5& rkV) const;

		// arithmetic operations
		inline Vector5 operator+ (const Vector5& rkV) const;
		inline Vector5 operator- (const Vector5& rkV) const;
		inline Vector5 operator* (Real fScalar) const;
		inline Vector5 operator/ (Real fScalar) const;
		inline Vector5 operator- () const;

		// arithmetic updates
		inline Vector5& operator+= (const Vector5& rkV);
		inline Vector5& operator-= (const Vector5& rkV);
		inline Vector5& operator*= (Real fScalar);
		inline Vector5& operator/= (Real fScalar);

		// vector operations
		inline Real Length () const;
		inline Real SquaredLength () const;
		inline Real Dot (const Vector5& rkV) const;
		inline Real Normalize ();

		// special vectors
		static const Vector5 ZERO;
		static const Vector5 UNIT_X;  // (1,0,0,0,0)
		static const Vector5 UNIT_Y;  // (0,1,0,0,0)
		static const Vector5 UNIT_Z;  // (0,0,1,0,0)
		static const Vector5 UNIT_W;  // (0,0,0,1,0)
		static const Vector5 UNIT_U;  // (0,0,0,0,1)
		static const Vector5 ONE;     // (1,1,1,1,1)

	private:
		// support for comparisons
		int CompareArrays (const Vector5& rkV) const;

		Real m_afTuple[5];
	};

	// arithmetic operations
	template <class Real>
	Vector5<Real> operator* (Real fScalar, const Vector5<Real>& rkV);

	// debugging output
	template <class Real>
	std::ostream& operator<< (std::ostream& rkOStr, const Vector5<Real>& rkV);

	// Wild Magic Source Code
	// David Eberly
	// http://www.geometrictools.com
	// Copyright (c) 1998-2007
	//
	// This library is free software; you can redistribute it and/or modify it
	// under the terms of the GNU Lesser General Public License as published by
	// the Free Software Foundation; either version 2.1 of the License, or (at
	// your option) any later version.  The license is available for reading at
	// either of the locations:
	//     http://www.gnu.org/copyleft/lgpl.html
	//     http://www.geometrictools.com/License/WildMagicLicense.pdf
	// The license applies to versions 0 through 4 of Wild Magic.
	//
	// Version: 4.0.2 (2007/03/07)

	//----------------------------------------------------------------------------
	template <class Real>
	Vector5<Real>::Vector5 ()
	{
		m_afTuple[0] = (Real)0.0;
		m_afTuple[1] = (Real)0.0;
		m_afTuple[2] = (Real)0.0;
		m_afTuple[3] = (Real)0.0;
		m_afTuple[4] = (Real)0.0;
		// uninitialized for performance in array construction
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Vector5<Real>::Vector5 (Real fX, Real fY, Real fZ, Real fW, Real fU)
	{
		m_afTuple[0] = fX;
		m_afTuple[1] = fY;
		m_afTuple[2] = fZ;
		m_afTuple[3] = fW;
		m_afTuple[4] = fU;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Vector5<Real>::Vector5 (const Real* afTuple)
	{
		m_afTuple[0] = afTuple[0];
		m_afTuple[1] = afTuple[1];
		m_afTuple[2] = afTuple[2];
		m_afTuple[3] = afTuple[3];
		m_afTuple[4] = afTuple[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Vector5<Real>::Vector5 (const Vector5& rkV)
	{
		m_afTuple[0] = rkV.m_afTuple[0];
		m_afTuple[1] = rkV.m_afTuple[1];
		m_afTuple[2] = rkV.m_afTuple[2];
		m_afTuple[3] = rkV.m_afTuple[3];
		m_afTuple[4] = rkV.m_afTuple[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>::operator const Real* () const
	{
		return m_afTuple;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>::operator Real* ()
	{
		return m_afTuple;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::operator[] (int i) const
	{
		return m_afTuple[i];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Vector5<Real>::operator[] (int i)
	{
		return m_afTuple[i];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::X () const
	{
		return m_afTuple[0];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Vector5<Real>::X ()
	{
		return m_afTuple[0];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::Y () const
	{
		return m_afTuple[1];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Vector5<Real>::Y ()
	{
		return m_afTuple[1];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::Z () const
	{
		return m_afTuple[2];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Vector5<Real>::Z ()
	{
		return m_afTuple[2];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::W () const
	{
		return m_afTuple[3];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Vector5<Real>::W ()
	{
		return m_afTuple[3];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::U () const
	{
		return m_afTuple[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Vector5<Real>::U ()
	{
		return m_afTuple[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>& Vector5<Real>::operator= (const Vector5& rkV)
	{
		m_afTuple[0] = rkV.m_afTuple[0];
		m_afTuple[1] = rkV.m_afTuple[1];
		m_afTuple[2] = rkV.m_afTuple[2];
		m_afTuple[3] = rkV.m_afTuple[3];
		m_afTuple[4] = rkV.m_afTuple[4];
		return *this;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	int Vector5<Real>::CompareArrays (const Vector5& rkV) const
	{
		return memcmp(m_afTuple,rkV.m_afTuple,5*sizeof(Real));
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Vector5<Real>::operator== (const Vector5& rkV) const
	{
		return CompareArrays(rkV) == 0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Vector5<Real>::operator!= (const Vector5& rkV) const
	{
		return CompareArrays(rkV) != 0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Vector5<Real>::operator< (const Vector5& rkV) const
	{
		return CompareArrays(rkV) < 0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Vector5<Real>::operator<= (const Vector5& rkV) const
	{
		return CompareArrays(rkV) <= 0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Vector5<Real>::operator> (const Vector5& rkV) const
	{
		return CompareArrays(rkV) > 0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Vector5<Real>::operator>= (const Vector5& rkV) const
	{
		return CompareArrays(rkV) >= 0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> Vector5<Real>::operator+ (const Vector5& rkV) const
	{
		return Vector5(
			m_afTuple[0]+rkV.m_afTuple[0],
			m_afTuple[1]+rkV.m_afTuple[1],
			m_afTuple[2]+rkV.m_afTuple[2],
			m_afTuple[3]+rkV.m_afTuple[3],
			m_afTuple[4]+rkV.m_afTuple[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> Vector5<Real>::operator- (const Vector5& rkV) const
	{
		return Vector5(
			m_afTuple[0]-rkV.m_afTuple[0],
			m_afTuple[1]-rkV.m_afTuple[1],
			m_afTuple[2]-rkV.m_afTuple[2],
			m_afTuple[3]-rkV.m_afTuple[3],
			m_afTuple[4]-rkV.m_afTuple[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> Vector5<Real>::operator* (Real fScalar) const
	{
		return Vector5(
			fScalar*m_afTuple[0],
			fScalar*m_afTuple[1],
			fScalar*m_afTuple[2],
			fScalar*m_afTuple[3],
			fScalar*m_afTuple[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> Vector5<Real>::operator/ (Real fScalar) const
	{
		Vector5 kQuot;

		if (fScalar != (Real)0.0)
		{
			Real fInvScalar = ((Real)1.0)/fScalar;
			kQuot.m_afTuple[0] = fInvScalar*m_afTuple[0];
			kQuot.m_afTuple[1] = fInvScalar*m_afTuple[1];
			kQuot.m_afTuple[2] = fInvScalar*m_afTuple[2];
			kQuot.m_afTuple[3] = fInvScalar*m_afTuple[3];
			kQuot.m_afTuple[4] = fInvScalar*m_afTuple[4];
		}
		else
		{
			kQuot.m_afTuple[0] = Math<Real>::MAX_REAL;
			kQuot.m_afTuple[1] = Math<Real>::MAX_REAL;
			kQuot.m_afTuple[2] = Math<Real>::MAX_REAL;
			kQuot.m_afTuple[3] = Math<Real>::MAX_REAL;
			kQuot.m_afTuple[4] = Math<Real>::MAX_REAL;
		}

		return kQuot;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> Vector5<Real>::operator- () const
	{
		return Vector5(
			-m_afTuple[0],
			-m_afTuple[1],
			-m_afTuple[2],
			-m_afTuple[3],
			-m_afTuple[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> operator* (Real fScalar, const Vector5<Real>& rkV)
	{
		return Vector5<Real>(
			fScalar*rkV[0],
			fScalar*rkV[1],
			fScalar*rkV[2],
			fScalar*rkV[3],
			fScalar*rkV[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>& Vector5<Real>::operator+= (const Vector5& rkV)
	{
		m_afTuple[0] += rkV.m_afTuple[0];
		m_afTuple[1] += rkV.m_afTuple[1];
		m_afTuple[2] += rkV.m_afTuple[2];
		m_afTuple[3] += rkV.m_afTuple[3];
		m_afTuple[4] += rkV.m_afTuple[4];
		return *this;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>& Vector5<Real>::operator-= (const Vector5& rkV)
	{
		m_afTuple[0] -= rkV.m_afTuple[0];
		m_afTuple[1] -= rkV.m_afTuple[1];
		m_afTuple[2] -= rkV.m_afTuple[2];
		m_afTuple[3] -= rkV.m_afTuple[3];
		m_afTuple[4] -= rkV.m_afTuple[4];
		return *this;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>& Vector5<Real>::operator*= (Real fScalar)
	{
		m_afTuple[0] *= fScalar;
		m_afTuple[1] *= fScalar;
		m_afTuple[2] *= fScalar;
		m_afTuple[3] *= fScalar;
		m_afTuple[4] *= fScalar;
		return *this;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real>& Vector5<Real>::operator/= (Real fScalar)
	{
		if (fScalar != (Real)0.0)
		{
			Real fInvScalar = ((Real)1.0)/fScalar;
			m_afTuple[0] *= fInvScalar;
			m_afTuple[1] *= fInvScalar;
			m_afTuple[2] *= fInvScalar;
			m_afTuple[3] *= fInvScalar;
			m_afTuple[4] *= fInvScalar;
		}
		else
		{
			m_afTuple[0] = Math<Real>::MAX_REAL;
			m_afTuple[1] = Math<Real>::MAX_REAL;
			m_afTuple[2] = Math<Real>::MAX_REAL;
			m_afTuple[3] = Math<Real>::MAX_REAL;
			m_afTuple[4] = Math<Real>::MAX_REAL;
		}

		return *this;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::Length () const
	{
		return Math<Real>::Sqrt(
			m_afTuple[0]*m_afTuple[0] +
			m_afTuple[1]*m_afTuple[1] +
			m_afTuple[2]*m_afTuple[2] +
			m_afTuple[3]*m_afTuple[3] +
			m_afTuple[4]*m_afTuple[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::SquaredLength () const
	{
		return
			m_afTuple[0]*m_afTuple[0] +
			m_afTuple[1]*m_afTuple[1] +
			m_afTuple[2]*m_afTuple[2] +
			m_afTuple[3]*m_afTuple[3] +
			m_afTuple[4]*m_afTuple[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::Dot (const Vector5& rkV) const
	{
		return
			m_afTuple[0]*rkV.m_afTuple[0] +
			m_afTuple[1]*rkV.m_afTuple[1] +
			m_afTuple[2]*rkV.m_afTuple[2] +
			m_afTuple[3]*rkV.m_afTuple[3] +
			m_afTuple[4]*rkV.m_afTuple[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Vector5<Real>::Normalize ()
	{
		Real fLength = Length();

		if (fLength > Math<Real>::ZERO_TOLERANCE)
		{
			Real fInvLength = ((Real)1.0)/fLength;
			m_afTuple[0] *= fInvLength;
			m_afTuple[1] *= fInvLength;
			m_afTuple[2] *= fInvLength;
			m_afTuple[3] *= fInvLength;
			m_afTuple[4] *= fInvLength;
		}
		else
		{
			fLength = (Real)0.0;
			m_afTuple[0] = (Real)0.0;
			m_afTuple[1] = (Real)0.0;
			m_afTuple[2] = (Real)0.0;
			m_afTuple[3] = (Real)0.0;
			m_afTuple[4] = (Real)0.0;
		}

		return fLength;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	std::ostream& operator<< (std::ostream& rkOStr, const Vector5<Real>& rkV)
	{
		return rkOStr << rkV.X() << ' ' << rkV.Y() << ' ' << rkV.Z()
			<< ' ' << rkV.W() << ' ' << rkV.U();
	}
	//----------------------------------------------------------------------------


	typedef Vector5<float> Vector5f;
	typedef Vector5<double> Vector5d;
	typedef Vector5<int> Vector5i;

}




#endif
