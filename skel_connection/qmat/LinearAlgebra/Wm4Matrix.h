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

#ifndef WM4MATRIX_H
#define WM4MATRIX_H

// Matrix operations are applied on the left.  For example, given a matrix M
// and a vector V, matrix-times-vector is M*V.  That is, V is treated as a
// column vector.  Some graphics APIs use V*M where V is treated as a row
// vector.  In this context the "M" matrix is really a transpose of the M as
// represented in Wild Magic.  Similarly, to apply two matrix operations M0
// and M1, in that order, you compute M1*M0 so that the transform of a vector
// is (M1*M0)*V = M1*(M0*V).  Some graphics APIs use M0*M1, but again these
// matrices are the transpose of those as represented in Wild Magic.  You
// must therefore be careful about how you interface the transformation code
// with graphics APIS.
//
// For memory organization it might seem natural to chose Real[N][N] for the
// matrix storage, but this can be a problem on a platform/console that
// chooses to store the data in column-major rather than row-major format.
// To avoid potential portability problems, the matrix is stored as Real[N*N]
// and organized in row-major order.  That is, the entry of the matrix in row
// r (0 <= r < N) and column c (0 <= c < N) is stored at index i = c+N*r
// (0 <= i < N*N).

// Rotation matrices are of the form
//   R = cos(t) -sin(t)
//       sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the xy-plane.

#include "Wm4Vector.h"

namespace Wm4
{

template <class Real>
class Matrix2
{
public:
    // If bZero is true, create the zero matrix.  Otherwise, create the
    // identity matrix.
    Matrix2 (bool bZero = true);

    // copy constructor
    Matrix2 (const Matrix2& rkM);

    // input Mrc is in row r, column c.
    Matrix2 (Real fM00, Real fM01, Real fM10, Real fM11);

    // Create a matrix from an array of numbers.  The input array is
    // interpreted based on the Boolean input as
    //   true:  entry[0..3] = {m00,m01,m10,m11}  [row major]
    //   false: entry[0..3] = {m00,m10,m01,m11}  [column major]
    Matrix2 (const Real afEntry[4], bool bRowMajor);

    // Create matrices based on vector input.  The Boolean is interpreted as
    //   true: vectors are columns of the matrix
    //   false: vectors are rows of the matrix
    Matrix2 (const Vector2<Real>& rkU, const Vector2<Real>& rkV,
        bool bColumns);
    Matrix2 (const Vector2<Real>* akV, bool bColumns);

    // create a diagonal matrix
    Matrix2 (Real fM00, Real fM11);

    // create a rotation matrix (positive angle - counterclockwise)
    Matrix2 (Real fAngle);

    // create a tensor product U*V^T
    Matrix2 (const Vector2<Real>& rkU, const Vector2<Real>& rkV);

    // create various matrices
    void MakeZero ();
    void MakeIdentity ();
    void MakeDiagonal (Real fM00, Real fM11);
    void FromAngle (Real fAngle);
    void MakeTensorProduct (const Vector2<Real>& rkU,
        const Vector2<Real>& rkV);

    // member access
    inline operator const Real* () const;
    inline operator Real* ();
    inline const Real* operator[] (int iRow) const;
    inline Real* operator[] (int iRow);
    inline Real operator() (int iRow, int iCol) const;
    inline Real& operator() (int iRow, int iCol);
    void SetRow (int iRow, const Vector2<Real>& rkV);
    Vector2<Real> GetRow (int iRow) const;
    void SetColumn (int iCol, const Vector2<Real>& rkV);
    Vector2<Real> GetColumn (int iCol) const;
    void GetColumnMajor (Real* afCMajor) const;

    // assignment
    inline Matrix2& operator= (const Matrix2& rkM);

    // comparison
    bool operator== (const Matrix2& rkM) const;
    bool operator!= (const Matrix2& rkM) const;
    bool operator<  (const Matrix2& rkM) const;
    bool operator<= (const Matrix2& rkM) const;
    bool operator>  (const Matrix2& rkM) const;
    bool operator>= (const Matrix2& rkM) const;

    // arithmetic operations
    inline Matrix2 operator+ (const Matrix2& rkM) const;
    inline Matrix2 operator- (const Matrix2& rkM) const;
    inline Matrix2 operator* (const Matrix2& rkM) const;
    inline Matrix2 operator* (Real fScalar) const;
    inline Matrix2 operator/ (Real fScalar) const;
    inline Matrix2 operator- () const;

    // arithmetic updates
    inline Matrix2& operator+= (const Matrix2& rkM);
    inline Matrix2& operator-= (const Matrix2& rkM);
    inline Matrix2& operator*= (Real fScalar);
    inline Matrix2& operator/= (Real fScalar);

    // matrix times vector
    inline Vector2<Real> operator* (const Vector2<Real>& rkV) const;  // M * v

    // other operations
    Matrix2 Transpose () const;  // M^T
    Matrix2 TransposeTimes (const Matrix2& rkM) const;  // this^T * M
    Matrix2 TimesTranspose (const Matrix2& rkM) const;  // this * M^T
    Matrix2 Inverse () const;
    Matrix2 Adjoint () const;
    Real Determinant () const;
    Real QForm (const Vector2<Real>& rkU,
        const Vector2<Real>& rkV) const;  // u^T*M*v

    // The matrix must be a rotation for these functions to be valid.  The
    // last function uses Gram-Schmidt orthonormalization applied to the
    // columns of the rotation matrix.  The angle must be in radians, not
    // degrees.
    void ToAngle (Real& rfAngle) const;
    void Orthonormalize ();

    // The matrix must be symmetric.  Factor M = R * D * R^T where
    // R = [u0|u1] is a rotation matrix with columns u0 and u1 and
    // D = diag(d0,d1) is a diagonal matrix whose diagonal entries are d0 and
    // d1.  The eigenvector u[i] corresponds to eigenvector d[i].  The
    // eigenvalues are ordered as d0 <= d1.
    void EigenDecomposition (Matrix2& rkRot, Matrix2& rkDiag) const;
	
	// special vectors
	//static const Matrix2 J; // (0,1,-1,0)

private:
    // support for comparisons
    int CompareArrays (const Matrix2& rkM) const;

    // matrix stored in row-major order
    Real m_afEntry[4];
};

// c * M
template <class Real>
inline Matrix2<Real> operator* (Real fScalar, const Matrix2<Real>& rkM);

// v^T * M
template <class Real>
inline Vector2<Real> operator* (const Vector2<Real>& rkV,
    const Matrix2<Real>& rkM);

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
Matrix2<Real>::Matrix2 (bool bZero)
{
    if (bZero)
    {
        MakeZero();
    }
    else
    {
        MakeIdentity();
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (const Matrix2& rkM)
{
    m_afEntry[0] = rkM.m_afEntry[0];
    m_afEntry[1] = rkM.m_afEntry[1];
    m_afEntry[2] = rkM.m_afEntry[2];
    m_afEntry[3] = rkM.m_afEntry[3];
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (Real fM00, Real fM01, Real fM10, Real fM11)
{
    m_afEntry[0] = fM00;
    m_afEntry[1] = fM01;
    m_afEntry[2] = fM10;
    m_afEntry[3] = fM11;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (const Real afEntry[4], bool bRowMajor)
{
    if (bRowMajor)
    {
        m_afEntry[0] = afEntry[0];
        m_afEntry[1] = afEntry[1];
        m_afEntry[2] = afEntry[2];
        m_afEntry[3] = afEntry[3];
    }
    else
    {
        m_afEntry[0] = afEntry[0];
        m_afEntry[1] = afEntry[2];
        m_afEntry[2] = afEntry[1];
        m_afEntry[3] = afEntry[3];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (const Vector2<Real>& rkU, const Vector2<Real>& rkV,
    bool bColumns)
{
    if (bColumns)
    {
        m_afEntry[0] = rkU[0];
        m_afEntry[1] = rkV[0];
        m_afEntry[2] = rkU[1];
        m_afEntry[3] = rkV[1];
    }
    else
    {
        m_afEntry[0] = rkU[0];
        m_afEntry[1] = rkU[1];
        m_afEntry[2] = rkV[0];
        m_afEntry[3] = rkV[1];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (const Vector2<Real>* akV, bool bColumns)
{
    if (bColumns)
    {
        m_afEntry[0] = akV[0][0];
        m_afEntry[1] = akV[1][0];
        m_afEntry[2] = akV[0][1];
        m_afEntry[3] = akV[1][1];
    }
    else
    {
        m_afEntry[0] = akV[0][0];
        m_afEntry[1] = akV[0][1];
        m_afEntry[2] = akV[1][0];
        m_afEntry[3] = akV[1][1];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (Real fM00, Real fM11)
{
    MakeDiagonal(fM00,fM11);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (Real fAngle)
{
    FromAngle(fAngle);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real>::Matrix2 (const Vector2<Real>& rkU, const Vector2<Real>& rkV)
{
    MakeTensorProduct(rkU,rkV);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>::operator const Real* () const
{
    return m_afEntry;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>::operator Real* ()
{
    return m_afEntry;
}
//----------------------------------------------------------------------------
template <class Real>
inline const Real* Matrix2<Real>::operator[] (int iRow) const
{
    return &m_afEntry[2*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real* Matrix2<Real>::operator[] (int iRow)
{
    return &m_afEntry[2*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Matrix2<Real>::operator() (int iRow, int iCol) const
{
    return m_afEntry[iCol + 2*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Matrix2<Real>::operator() (int iRow, int iCol)
{
    return m_afEntry[iCol + 2*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::MakeZero ()
{
    m_afEntry[0] = (Real)0.0;
    m_afEntry[1] = (Real)0.0;
    m_afEntry[2] = (Real)0.0;
    m_afEntry[3] = (Real)0.0;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::MakeIdentity ()
{
    m_afEntry[0] = (Real)1.0;
    m_afEntry[1] = (Real)0.0;
    m_afEntry[2] = (Real)0.0;
    m_afEntry[3] = (Real)1.0;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::MakeDiagonal (Real fM00, Real fM11)
{
    m_afEntry[0] = fM00;
    m_afEntry[1] = (Real)0.0;
    m_afEntry[2] = (Real)0.0;
    m_afEntry[3] = fM11;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::FromAngle (Real fAngle)
{
    m_afEntry[0] = Math<Real>::Cos(fAngle);
    m_afEntry[2] = Math<Real>::Sin(fAngle);
    m_afEntry[1] = -m_afEntry[2];
    m_afEntry[3] =  m_afEntry[0];
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::MakeTensorProduct (const Vector2<Real>& rkU,
    const Vector2<Real>& rkV)
{
    m_afEntry[0] = rkU[0]*rkV[0];
    m_afEntry[1] = rkU[0]*rkV[1];
    m_afEntry[2] = rkU[1]*rkV[0];
    m_afEntry[3] = rkU[1]*rkV[1];
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::SetRow (int iRow, const Vector2<Real>& rkV)
{
    int i0 = 2*iRow ,i1 = i0+1;
    m_afEntry[i0] = rkV[0];
    m_afEntry[i1] = rkV[1];
}
//----------------------------------------------------------------------------
template <class Real>
Vector2<Real> Matrix2<Real>::GetRow (int iRow) const
{
    int i0 = 2*iRow ,i1 = i0+1;
    return Vector2<Real>(m_afEntry[i0],m_afEntry[i1]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::SetColumn (int iCol, const Vector2<Real>& rkV)
{
    m_afEntry[iCol] = rkV[0];
    m_afEntry[iCol+2] = rkV[1];
}
//----------------------------------------------------------------------------
template <class Real>
Vector2<Real> Matrix2<Real>::GetColumn (int iCol) const
{
    return Vector2<Real>(m_afEntry[iCol],m_afEntry[iCol+2]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::GetColumnMajor (Real* afCMajor) const
{
    afCMajor[0] = m_afEntry[0];
    afCMajor[1] = m_afEntry[2];
    afCMajor[2] = m_afEntry[1];
    afCMajor[3] = m_afEntry[3];
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>& Matrix2<Real>::operator= (const Matrix2& rkM)
{
    m_afEntry[0] = rkM.m_afEntry[0];
    m_afEntry[1] = rkM.m_afEntry[1];
    m_afEntry[2] = rkM.m_afEntry[2];
    m_afEntry[3] = rkM.m_afEntry[3];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
int Matrix2<Real>::CompareArrays (const Matrix2& rkM) const
{
    return memcmp(m_afEntry,rkM.m_afEntry,4*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix2<Real>::operator== (const Matrix2& rkM) const
{
    return CompareArrays(rkM) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix2<Real>::operator!= (const Matrix2& rkM) const
{
    return CompareArrays(rkM) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix2<Real>::operator<  (const Matrix2& rkM) const
{
    return CompareArrays(rkM) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix2<Real>::operator<= (const Matrix2& rkM) const
{
    return CompareArrays(rkM) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix2<Real>::operator>  (const Matrix2& rkM) const
{
    return CompareArrays(rkM) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix2<Real>::operator>= (const Matrix2& rkM) const
{
    return CompareArrays(rkM) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> Matrix2<Real>::operator+ (const Matrix2& rkM) const

{
    return Matrix2<Real>(
        m_afEntry[0] + rkM.m_afEntry[0],
        m_afEntry[1] + rkM.m_afEntry[1],
        m_afEntry[2] + rkM.m_afEntry[2],
        m_afEntry[3] + rkM.m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> Matrix2<Real>::operator- (const Matrix2& rkM) const
{
    return Matrix2<Real>(
        m_afEntry[0] - rkM.m_afEntry[0],
        m_afEntry[1] - rkM.m_afEntry[1],
        m_afEntry[2] - rkM.m_afEntry[2],
        m_afEntry[3] - rkM.m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> Matrix2<Real>::operator* (const Matrix2& rkM) const
{
    return Matrix2<Real>(
        m_afEntry[0]*rkM.m_afEntry[0] + m_afEntry[1]*rkM.m_afEntry[2],
        m_afEntry[0]*rkM.m_afEntry[1] + m_afEntry[1]*rkM.m_afEntry[3],
        m_afEntry[2]*rkM.m_afEntry[0] + m_afEntry[3]*rkM.m_afEntry[2],
        m_afEntry[2]*rkM.m_afEntry[1] + m_afEntry[3]*rkM.m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> Matrix2<Real>::operator* (Real fScalar) const
{
    return Matrix2<Real>(
        fScalar*m_afEntry[0],
        fScalar*m_afEntry[1],
        fScalar*m_afEntry[2],
        fScalar*m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> Matrix2<Real>::operator/ (Real fScalar) const
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        return Matrix2<Real>(
            fInvScalar*m_afEntry[0],
            fInvScalar*m_afEntry[1],
            fInvScalar*m_afEntry[2],
            fInvScalar*m_afEntry[3]);
    }

    return Matrix2<Real>(
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> Matrix2<Real>::operator- () const
{
    return Matrix2<Real>(
        -m_afEntry[0],
        -m_afEntry[1],
        -m_afEntry[2],
        -m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>& Matrix2<Real>::operator+= (const Matrix2& rkM)
{
    m_afEntry[0] += rkM.m_afEntry[0];
    m_afEntry[1] += rkM.m_afEntry[1];
    m_afEntry[2] += rkM.m_afEntry[2];
    m_afEntry[3] += rkM.m_afEntry[3];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>& Matrix2<Real>::operator-= (const Matrix2& rkM)
{
    m_afEntry[0] -= rkM.m_afEntry[0];
    m_afEntry[1] -= rkM.m_afEntry[1];
    m_afEntry[2] -= rkM.m_afEntry[2];
    m_afEntry[3] -= rkM.m_afEntry[3];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>& Matrix2<Real>::operator*= (Real fScalar)
{
    m_afEntry[0] *= fScalar;
    m_afEntry[1] *= fScalar;
    m_afEntry[2] *= fScalar;
    m_afEntry[3] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real>& Matrix2<Real>::operator/= (Real fScalar)
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        m_afEntry[0] *= fInvScalar;
        m_afEntry[1] *= fInvScalar;
        m_afEntry[2] *= fInvScalar;
        m_afEntry[3] *= fInvScalar;
    }
    else
    {
        m_afEntry[0] = Math<Real>::MAX_REAL;
        m_afEntry[1] = Math<Real>::MAX_REAL;
        m_afEntry[2] = Math<Real>::MAX_REAL;
        m_afEntry[3] = Math<Real>::MAX_REAL;
    }

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> Matrix2<Real>::operator* (const Vector2<Real>& rkV) const
{
    return Vector2<Real>(
        m_afEntry[0]*rkV[0] + m_afEntry[1]*rkV[1],
        m_afEntry[2]*rkV[0] + m_afEntry[3]*rkV[1]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real> Matrix2<Real>::Transpose () const
{
    return Matrix2<Real>(
        m_afEntry[0],
        m_afEntry[2],
        m_afEntry[1],
        m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real> Matrix2<Real>::TransposeTimes (const Matrix2& rkM) const
{
    // P = A^T*B
    return Matrix2<Real>(
        m_afEntry[0]*rkM.m_afEntry[0] + m_afEntry[2]*rkM.m_afEntry[2],
        m_afEntry[0]*rkM.m_afEntry[1] + m_afEntry[2]*rkM.m_afEntry[3],
        m_afEntry[1]*rkM.m_afEntry[0] + m_afEntry[3]*rkM.m_afEntry[2],
        m_afEntry[1]*rkM.m_afEntry[1] + m_afEntry[3]*rkM.m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real> Matrix2<Real>::TimesTranspose (const Matrix2& rkM) const
{
    // P = A*B^T
    return Matrix2<Real>(
        m_afEntry[0]*rkM.m_afEntry[0] + m_afEntry[1]*rkM.m_afEntry[1],
        m_afEntry[0]*rkM.m_afEntry[2] + m_afEntry[1]*rkM.m_afEntry[3],
        m_afEntry[2]*rkM.m_afEntry[0] + m_afEntry[3]*rkM.m_afEntry[1],
        m_afEntry[2]*rkM.m_afEntry[2] + m_afEntry[3]*rkM.m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real> Matrix2<Real>::Inverse () const
{
    Matrix2<Real> kInverse;

    Real fDet = m_afEntry[0]*m_afEntry[3] - m_afEntry[1]*m_afEntry[2];
    if (Math<Real>::FAbs(fDet) > Math<Real>::INVERSE_TOLERANCE)
    {
        Real fInvDet = ((Real)1.0)/fDet;
        kInverse.m_afEntry[0] =  m_afEntry[3]*fInvDet;
        kInverse.m_afEntry[1] = -m_afEntry[1]*fInvDet;
        kInverse.m_afEntry[2] = -m_afEntry[2]*fInvDet;
        kInverse.m_afEntry[3] =  m_afEntry[0]*fInvDet;
    }
    else
    {
        kInverse.m_afEntry[0] = (Real)0.0;
        kInverse.m_afEntry[1] = (Real)0.0;
        kInverse.m_afEntry[2] = (Real)0.0;
        kInverse.m_afEntry[3] = (Real)0.0;
    }

    return kInverse;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix2<Real> Matrix2<Real>::Adjoint () const
{
    return Matrix2<Real>(
        m_afEntry[3],
        -m_afEntry[1],
        -m_afEntry[2],
        m_afEntry[0]);
}
//----------------------------------------------------------------------------
template <class Real>
Real Matrix2<Real>::Determinant () const
{
    return m_afEntry[0]*m_afEntry[3] - m_afEntry[1]*m_afEntry[2];
}
//----------------------------------------------------------------------------
template <class Real>
Real Matrix2<Real>::QForm (const Vector2<Real>& rkU,
    const Vector2<Real>& rkV) const
{
    return rkU.Dot((*this)*rkV);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::ToAngle (Real& rfAngle) const
{
    // assert:  matrix is a rotation
    rfAngle = Math<Real>::ATan2(m_afEntry[2],m_afEntry[0]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::Orthonormalize ()
{
    // Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is
    // M = [m0|m1], then orthonormal output matrix is Q = [q0|q1],
    //
    //   q0 = m0/|m0|
    //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
    //
    // where |V| indicates length of vector V and A*B indicates dot
    // product of vectors A and B.

    // compute q0
    Real fInvLength = Math<Real>::InvSqrt(m_afEntry[0]*m_afEntry[0] +
        m_afEntry[2]*m_afEntry[2]);

    m_afEntry[0] *= fInvLength;
    m_afEntry[2] *= fInvLength;

    // compute q1
    Real fDot0 = m_afEntry[0]*m_afEntry[1] + m_afEntry[2]*m_afEntry[3];
    m_afEntry[1] -= fDot0*m_afEntry[0];
    m_afEntry[3] -= fDot0*m_afEntry[2];

    fInvLength = Math<Real>::InvSqrt(m_afEntry[1]*m_afEntry[1] +
        m_afEntry[3]*m_afEntry[3]);

    m_afEntry[1] *= fInvLength;
    m_afEntry[3] *= fInvLength;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix2<Real>::EigenDecomposition (Matrix2& rkRot, Matrix2& rkDiag) const
{
    Real fTrace = m_afEntry[0] + m_afEntry[3];
    Real fDiff = m_afEntry[0] - m_afEntry[3];
    Real fDiscr = Math<Real>::Sqrt(fDiff*fDiff +
        ((Real)4.0)*m_afEntry[1]*m_afEntry[1]);
    Real fEVal0 = ((Real)0.5)*(fTrace-fDiscr);
    Real fEVal1 = ((Real)0.5)*(fTrace+fDiscr);
    rkDiag.MakeDiagonal(fEVal0,fEVal1);

    Real fCos, fSin;
    if (fDiff >= (Real)0.0)
    {
        fCos = m_afEntry[1];
        fSin = fEVal0 - m_afEntry[0];
    }
    else
    {
        fCos = fEVal0 - m_afEntry[3];
        fSin = m_afEntry[1];
    }
    Real fTmp = Math<Real>::InvSqrt(fCos*fCos + fSin*fSin);
    fCos *= fTmp;
    fSin *= fTmp;

    rkRot.m_afEntry[0] = fCos;
    rkRot.m_afEntry[1] = -fSin;
    rkRot.m_afEntry[2] = fSin;
    rkRot.m_afEntry[3] = fCos;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix2<Real> operator* (Real fScalar, const Matrix2<Real>& rkM)
{
    return rkM*fScalar;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector2<Real> operator* (const Vector2<Real>& rkV,
    const Matrix2<Real>& rkM)
{
    return Vector2<Real>(
        rkV[0]*rkM[0][0] + rkV[1]*rkM[1][0],
        rkV[0]*rkM[0][1] + rkV[1]*rkM[1][1]);
}
//----------------------------------------------------------------------------
//template<> const Matrix2<float> Matrix2<float>::J(0,1,-1,0);
//template<> const Matrix2<double> Matrix2<double>::J(0,1,-1,0);


typedef Matrix2<float> Matrix2f;
typedef Matrix2<double> Matrix2d;

}





namespace Wm4
{

template <class Real>
class Matrix3
{
public:
    // If bZero is true, create the zero matrix.  Otherwise, create the
    // identity matrix.
    Matrix3 (bool bZero = true);

    // copy constructor
    Matrix3 (const Matrix3& rkM);

    // input Mrc is in row r, column c.
    Matrix3 (Real fM00, Real fM01, Real fM02,
             Real fM10, Real fM11, Real fM12,
             Real fM20, Real fM21, Real fM22);

    // Create a matrix from an array of numbers.  The input array is
    // interpreted based on the Boolean input as
    //   true:  entry[0..8]={m00,m01,m02,m10,m11,m12,m20,m21,m22} [row major]
    //   false: entry[0..8]={m00,m10,m20,m01,m11,m21,m02,m12,m22} [col major]
    Matrix3 (const Real afEntry[9], bool bRowMajor);

    // Create matrices based on vector input.  The Boolean is interpreted as
    //   true: vectors are columns of the matrix
    //   false: vectors are rows of the matrix
    Matrix3 (const Vector3<Real>& rkU, const Vector3<Real>& rkV,
        const Vector3<Real>& rkW, bool bColumns);
    Matrix3 (const Vector3<Real>* akV, bool bColumns);

    // create a diagonal matrix
    Matrix3 (Real fM00, Real fM11, Real fM22);

    // Create rotation matrices (positive angle - counterclockwise).  The
    // angle must be in radians, not degrees.
    Matrix3 (const Vector3<Real>& rkAxis, Real fAngle);

    // create a tensor product U*V^T
    Matrix3 (const Vector3<Real>& rkU, const Vector3<Real>& rkV);

    // create various matrices
    Matrix3& MakeZero ();
    Matrix3& MakeIdentity ();
    Matrix3& MakeDiagonal (Real fM00, Real fM11, Real fM22);
    Matrix3& FromAxisAngle (const Vector3<Real>& rkAxis, Real fAngle);
    Matrix3& MakeTensorProduct (const Vector3<Real>& rkU,
        const Vector3<Real>& rkV);

    // member access
    inline operator const Real* () const;
    inline operator Real* ();
    inline const Real* operator[] (int iRow) const;
    inline Real* operator[] (int iRow);
    inline Real operator() (int iRow, int iCol) const;
    inline Real& operator() (int iRow, int iCol);
    void SetRow (int iRow, const Vector3<Real>& rkV);
    Vector3<Real> GetRow (int iRow) const;
    void SetColumn (int iCol, const Vector3<Real>& rkV);
    Vector3<Real> GetColumn (int iCol) const;
    void GetColumnMajor (Real* afCMajor) const;

    // assignment
    inline Matrix3& operator= (const Matrix3& rkM);

    // comparison
    bool operator== (const Matrix3& rkM) const;
    bool operator!= (const Matrix3& rkM) const;
    bool operator<  (const Matrix3& rkM) const;
    bool operator<= (const Matrix3& rkM) const;
    bool operator>  (const Matrix3& rkM) const;
    bool operator>= (const Matrix3& rkM) const;

    // arithmetic operations
    inline Matrix3 operator+ (const Matrix3& rkM) const;
    inline Matrix3 operator- (const Matrix3& rkM) const;
    inline Matrix3 operator* (const Matrix3& rkM) const;
    inline Matrix3 operator* (Real fScalar) const;
    inline Matrix3 operator/ (Real fScalar) const;
    inline Matrix3 operator- () const;

    // arithmetic updates
    inline Matrix3& operator+= (const Matrix3& rkM);
    inline Matrix3& operator-= (const Matrix3& rkM);
    inline Matrix3& operator*= (Real fScalar);
    inline Matrix3& operator/= (Real fScalar);

    // matrix times vector
    inline Vector3<Real> operator* (const Vector3<Real>& rkV) const;  // M * v

    // other operations
    Matrix3 Transpose () const;  // M^T
    Matrix3 TransposeTimes (const Matrix3& rkM) const;  // this^T * M
    Matrix3 TimesTranspose (const Matrix3& rkM) const;  // this * M^T
    Matrix3 Inverse () const;
    Matrix3 Adjoint () const;
    Real Determinant () const;
    Real QForm (const Vector3<Real>& rkU,
        const Vector3<Real>& rkV) const;  // u^T*M*v
    Matrix3 TimesDiagonal (const Vector3<Real>& rkDiag) const;  // M*D
    Matrix3 DiagonalTimes (const Vector3<Real>& rkDiag) const;  // D*M

    // The matrix must be a rotation for these functions to be valid.  The
    // last function uses Gram-Schmidt orthonormalization applied to the
    // columns of the rotation matrix.  The angle must be in radians, not
    // degrees.
    void ToAxisAngle (Vector3<Real>& rkAxis, Real& rfAngle) const;
    void Orthonormalize ();

    // The matrix must be symmetric.  Factor M = R * D * R^T where
    // R = [u0|u1|u2] is a rotation matrix with columns u0, u1, and u2 and
    // D = diag(d0,d1,d2) is a diagonal matrix whose diagonal entries are d0,
    // d1, and d2.  The eigenvector u[i] corresponds to eigenvector d[i].
    // The eigenvalues are ordered as d0 <= d1 <= d2.
    void EigenDecomposition (Matrix3& rkRot, Matrix3& rkDiag) const;

    // Create rotation matrices from Euler angles.
    Matrix3& FromEulerAnglesXYZ (Real fXAngle, Real fYAngle, Real fZAngle);
    Matrix3& FromEulerAnglesXZY (Real fXAngle, Real fZAngle, Real fYAngle);
    Matrix3& FromEulerAnglesYXZ (Real fYAngle, Real fXAngle, Real fZAngle);
    Matrix3& FromEulerAnglesYZX (Real fYAngle, Real fZAngle, Real fXAngle);
    Matrix3& FromEulerAnglesZXY (Real fZAngle, Real fXAngle, Real fYAngle);
    Matrix3& FromEulerAnglesZYX (Real fZAngle, Real fYAngle, Real fXAngle);

    // Extract Euler angles from rotation matrices.  The return value is
    // 'true' iff the factorization is unique relative to certain angle
    // ranges.  That is, if (U,V,W) is some permutation of (X,Y,Z), the angle
    // ranges for the outputs from ToEulerAnglesUVW(uAngle,vAngle,wAngle) are
    // uAngle in [-pi,pi], vAngle in [-pi/2,pi/2], and wAngle in [-pi,pi].  If
    // the function returns 'false', wAngle is 0 and vAngle is either pi/2 or
    // -pi/2.
    bool ToEulerAnglesXYZ (Real& rfXAngle, Real& rfYAngle, Real& rfZAngle)
        const;
    bool ToEulerAnglesXZY (Real& rfXAngle, Real& rfZAngle, Real& rfYAngle)
        const;
    bool ToEulerAnglesYXZ (Real& rfYAngle, Real& rfXAngle, Real& rfZAngle)
        const;
    bool ToEulerAnglesYZX (Real& rfYAngle, Real& rfZAngle, Real& rfXAngle)
        const;
    bool ToEulerAnglesZXY (Real& rfZAngle, Real& rfXAngle, Real& rfYAngle)
        const;
    bool ToEulerAnglesZYX (Real& rfZAngle, Real& rfYAngle, Real& rfXAngle)
        const;

    // SLERP (spherical linear interpolation) without quaternions.  Computes
    // R(t) = R0*(Transpose(R0)*R1)^t.  If Q is a rotation matrix with
    // unit-length axis U and angle A, then Q^t is a rotation matrix with
    // unit-length axis U and rotation angle t*A.
    Matrix3& Slerp (Real fT, const Matrix3& rkR0, const Matrix3& rkR1);

    // Singular value decomposition, M = L*S*R, where L and R are orthogonal
    // and S is a diagonal matrix whose diagonal entries are nonnegative.
    void SingularValueDecomposition (Matrix3& rkL, Matrix3& rkS,
        Matrix3& rkR) const;
    void SingularValueComposition (const Matrix3& rkL, const Matrix3& rkS,
        const Matrix3& rkR);
	
	// Polar decomposition, M = Q*S, where Q is orthogonal and S is symmetric.
	// This uses the singular value decomposition:
	//   M = L*D*Transpose(R) = (L*Transpose(R))*(R*D*Transpose(R)) = Q*S
	// where Q = L*Transpose(R) and S = R*D*Transpose(R).
	void PolarDecomposition (Matrix3& rkQ, Matrix3& rkS);

    // factor M = Q*D*U with orthogonal Q, diagonal D, upper triangular U
    void QDUDecomposition (Matrix3& rkQ, Matrix3& rkD, Matrix3& rkU) const;

private:
    // Support for eigendecomposition.  The Tridiagonalize function applies
    // a Householder transformation to the matrix.  If that transformation
    // is the identity (the matrix is already tridiagonal), then the return
    // value is 'false'.  Otherwise, the transformation is a reflection and
    // the return value is 'true'.  The QLAlgorithm returns 'true' iff the
    // QL iteration scheme converged.
    bool Tridiagonalize (Real afDiag[3], Real afSubd[2]);
    bool QLAlgorithm (Real afDiag[3], Real afSubd[2]);

    // support for singular value decomposition
    static void Bidiagonalize (Matrix3& rkA, Matrix3& rkL, Matrix3& rkR);
    static void GolubKahanStep (Matrix3& rkA, Matrix3& rkL, Matrix3& rkR);

    // support for comparisons
    int CompareArrays (const Matrix3& rkM) const;

    Real m_afEntry[9];
};

// c * M
template <class Real>
inline Matrix3<Real> operator* (Real fScalar, const Matrix3<Real>& rkM);

// v^T * M
template <class Real>
inline Vector3<Real> operator* (const Vector3<Real>& rkV,
    const Matrix3<Real>& rkM);

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
// Version: 4.0.3 (2006/03/07)

//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (bool bZero)
{
    if (bZero)
    {
        MakeZero();
    }
    else
    {
        MakeIdentity();
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (const Matrix3& rkM)
{
    m_afEntry[0] = rkM.m_afEntry[0];
    m_afEntry[1] = rkM.m_afEntry[1];
    m_afEntry[2] = rkM.m_afEntry[2];
    m_afEntry[3] = rkM.m_afEntry[3];
    m_afEntry[4] = rkM.m_afEntry[4];
    m_afEntry[5] = rkM.m_afEntry[5];
    m_afEntry[6] = rkM.m_afEntry[6];
    m_afEntry[7] = rkM.m_afEntry[7];
    m_afEntry[8] = rkM.m_afEntry[8];
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (Real fM00, Real fM01, Real fM02, Real fM10, Real fM11,
    Real fM12, Real fM20, Real fM21, Real fM22)
{
    m_afEntry[0] = fM00;
    m_afEntry[1] = fM01;
    m_afEntry[2] = fM02;
    m_afEntry[3] = fM10;
    m_afEntry[4] = fM11;
    m_afEntry[5] = fM12;
    m_afEntry[6] = fM20;
    m_afEntry[7] = fM21;
    m_afEntry[8] = fM22;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (const Real afEntry[9], bool bRowMajor)
{
    if (bRowMajor)
    {
        m_afEntry[0] = afEntry[0];
        m_afEntry[1] = afEntry[1];
        m_afEntry[2] = afEntry[2];
        m_afEntry[3] = afEntry[3];
        m_afEntry[4] = afEntry[4];
        m_afEntry[5] = afEntry[5];
        m_afEntry[6] = afEntry[6];
        m_afEntry[7] = afEntry[7];
        m_afEntry[8] = afEntry[8];
    }
    else
    {
        m_afEntry[0] = afEntry[0];
        m_afEntry[1] = afEntry[3];
        m_afEntry[2] = afEntry[6];
        m_afEntry[3] = afEntry[1];
        m_afEntry[4] = afEntry[4];
        m_afEntry[5] = afEntry[7];
        m_afEntry[6] = afEntry[2];
        m_afEntry[7] = afEntry[5];
        m_afEntry[8] = afEntry[8];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (const Vector3<Real>& rkU, const Vector3<Real>& rkV,
    const Vector3<Real>& rkW, bool bColumns)
{
    if (bColumns)
    {
        m_afEntry[0] = rkU[0];
        m_afEntry[1] = rkV[0];
        m_afEntry[2] = rkW[0];
        m_afEntry[3] = rkU[1];
        m_afEntry[4] = rkV[1];
        m_afEntry[5] = rkW[1];
        m_afEntry[6] = rkU[2];
        m_afEntry[7] = rkV[2];
        m_afEntry[8] = rkW[2];
    }
    else
    {
        m_afEntry[0] = rkU[0];
        m_afEntry[1] = rkU[1];
        m_afEntry[2] = rkU[2];
        m_afEntry[3] = rkV[0];
        m_afEntry[4] = rkV[1];
        m_afEntry[5] = rkV[2];
        m_afEntry[6] = rkW[0];
        m_afEntry[7] = rkW[1];
        m_afEntry[8] = rkW[2];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (const Vector3<Real>* akV, bool bColumns)
{
    if (bColumns)
    {
        m_afEntry[0] = akV[0][0];
        m_afEntry[1] = akV[1][0];
        m_afEntry[2] = akV[2][0];
        m_afEntry[3] = akV[0][1];
        m_afEntry[4] = akV[1][1];
        m_afEntry[5] = akV[2][1];
        m_afEntry[6] = akV[0][2];
        m_afEntry[7] = akV[1][2];
        m_afEntry[8] = akV[2][2];
    }
    else
    {
        m_afEntry[0] = akV[0][0];
        m_afEntry[1] = akV[0][1];
        m_afEntry[2] = akV[0][2];
        m_afEntry[3] = akV[1][0];
        m_afEntry[4] = akV[1][1];
        m_afEntry[5] = akV[1][2];
        m_afEntry[6] = akV[2][0];
        m_afEntry[7] = akV[2][1];
        m_afEntry[8] = akV[2][2];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (Real fM00, Real fM11, Real fM22)
{
    MakeDiagonal(fM00,fM11,fM22);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (const Vector3<Real>& rkAxis, Real fAngle)
{
    FromAxisAngle(rkAxis,fAngle);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>::Matrix3 (const Vector3<Real>& rkU, const Vector3<Real>& rkV)
{
    MakeTensorProduct(rkU,rkV);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>::operator const Real* () const
{
    return m_afEntry;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>::operator Real* ()
{
    return m_afEntry;
}
//----------------------------------------------------------------------------
template <class Real>
inline const Real* Matrix3<Real>::operator[] (int iRow) const
{
    return &m_afEntry[3*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real* Matrix3<Real>::operator[] (int iRow)
{
    return &m_afEntry[3*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Matrix3<Real>::operator() (int iRow, int iCol) const
{
    return m_afEntry[iCol+3*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Matrix3<Real>::operator() (int iRow, int iCol)
{
    return m_afEntry[iCol+3*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::MakeZero ()
{
    m_afEntry[0] = (Real)0.0;
    m_afEntry[1] = (Real)0.0;
    m_afEntry[2] = (Real)0.0;
    m_afEntry[3] = (Real)0.0;
    m_afEntry[4] = (Real)0.0;
    m_afEntry[5] = (Real)0.0;
    m_afEntry[6] = (Real)0.0;
    m_afEntry[7] = (Real)0.0;
    m_afEntry[8] = (Real)0.0;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::MakeIdentity ()
{
    m_afEntry[0] = (Real)1.0;
    m_afEntry[1] = (Real)0.0;
    m_afEntry[2] = (Real)0.0;
    m_afEntry[3] = (Real)0.0;
    m_afEntry[4] = (Real)1.0;
    m_afEntry[5] = (Real)0.0;
    m_afEntry[6] = (Real)0.0;
    m_afEntry[7] = (Real)0.0;
    m_afEntry[8] = (Real)1.0;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::MakeDiagonal (Real fM00, Real fM11, Real fM22)
{
    m_afEntry[0] = fM00;
    m_afEntry[1] = (Real)0.0;
    m_afEntry[2] = (Real)0.0;
    m_afEntry[3] = (Real)0.0;
    m_afEntry[4] = fM11;
    m_afEntry[5] = (Real)0.0;
    m_afEntry[6] = (Real)0.0;
    m_afEntry[7] = (Real)0.0;
    m_afEntry[8] = fM22;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromAxisAngle (const Vector3<Real>& rkAxis,
    Real fAngle)
{
    Real fCos = Math<Real>::Cos(fAngle);
    Real fSin = Math<Real>::Sin(fAngle);
    Real fOneMinusCos = ((Real)1.0)-fCos;
    Real fX2 = rkAxis[0]*rkAxis[0];
    Real fY2 = rkAxis[1]*rkAxis[1];
    Real fZ2 = rkAxis[2]*rkAxis[2];
    Real fXYM = rkAxis[0]*rkAxis[1]*fOneMinusCos;
    Real fXZM = rkAxis[0]*rkAxis[2]*fOneMinusCos;
    Real fYZM = rkAxis[1]*rkAxis[2]*fOneMinusCos;
    Real fXSin = rkAxis[0]*fSin;
    Real fYSin = rkAxis[1]*fSin;
    Real fZSin = rkAxis[2]*fSin;
    
    m_afEntry[0] = fX2*fOneMinusCos+fCos;
    m_afEntry[1] = fXYM-fZSin;
    m_afEntry[2] = fXZM+fYSin;
    m_afEntry[3] = fXYM+fZSin;
    m_afEntry[4] = fY2*fOneMinusCos+fCos;
    m_afEntry[5] = fYZM-fXSin;
    m_afEntry[6] = fXZM-fYSin;
    m_afEntry[7] = fYZM+fXSin;
    m_afEntry[8] = fZ2*fOneMinusCos+fCos;

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::MakeTensorProduct (const Vector3<Real>& rkU,
    const Vector3<Real>& rkV)
{
    m_afEntry[0] = rkU[0]*rkV[0];
    m_afEntry[1] = rkU[0]*rkV[1];
    m_afEntry[2] = rkU[0]*rkV[2];
    m_afEntry[3] = rkU[1]*rkV[0];
    m_afEntry[4] = rkU[1]*rkV[1];
    m_afEntry[5] = rkU[1]*rkV[2];
    m_afEntry[6] = rkU[2]*rkV[0];
    m_afEntry[7] = rkU[2]*rkV[1];
    m_afEntry[8] = rkU[2]*rkV[2];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::SetRow (int iRow, const Vector3<Real>& rkV)
{
    int i0 = 3*iRow, i1 = i0+1, i2 = i1+1;
    m_afEntry[i0] = rkV[0];
    m_afEntry[i1] = rkV[1];
    m_afEntry[i2] = rkV[2];
}
//----------------------------------------------------------------------------
template <class Real>
Vector3<Real> Matrix3<Real>::GetRow (int iRow) const
{
    int i0 = 3*iRow, i1 = i0+1, i2 = i1+1;
    return Vector3<Real>(m_afEntry[i0],m_afEntry[i1],m_afEntry[i2]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::SetColumn (int iCol, const Vector3<Real>& rkV)
{
    m_afEntry[iCol] = rkV[0];
    m_afEntry[iCol+3] = rkV[1];
    m_afEntry[iCol+6] = rkV[2];
}
//----------------------------------------------------------------------------
template <class Real>
Vector3<Real> Matrix3<Real>::GetColumn (int iCol) const
{
    return Vector3<Real>(m_afEntry[iCol],m_afEntry[iCol+3],m_afEntry[iCol+6]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::GetColumnMajor (Real* afCMajor) const
{
    afCMajor[0] = m_afEntry[0];
    afCMajor[1] = m_afEntry[3];
    afCMajor[2] = m_afEntry[6];
    afCMajor[3] = m_afEntry[1];
    afCMajor[4] = m_afEntry[4];
    afCMajor[5] = m_afEntry[7];
    afCMajor[6] = m_afEntry[2];
    afCMajor[7] = m_afEntry[5];
    afCMajor[8] = m_afEntry[8];
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>& Matrix3<Real>::operator= (const Matrix3& rkM)
{
    m_afEntry[0] = rkM.m_afEntry[0];
    m_afEntry[1] = rkM.m_afEntry[1];
    m_afEntry[2] = rkM.m_afEntry[2];
    m_afEntry[3] = rkM.m_afEntry[3];
    m_afEntry[4] = rkM.m_afEntry[4];
    m_afEntry[5] = rkM.m_afEntry[5];
    m_afEntry[6] = rkM.m_afEntry[6];
    m_afEntry[7] = rkM.m_afEntry[7];
    m_afEntry[8] = rkM.m_afEntry[8];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
int Matrix3<Real>::CompareArrays (const Matrix3& rkM) const
{
    return memcmp(m_afEntry,rkM.m_afEntry,9*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::operator== (const Matrix3& rkM) const
{
    return CompareArrays(rkM) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::operator!= (const Matrix3& rkM) const
{
    return CompareArrays(rkM) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::operator<  (const Matrix3& rkM) const
{
    return CompareArrays(rkM) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::operator<= (const Matrix3& rkM) const
{
    return CompareArrays(rkM) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::operator>  (const Matrix3& rkM) const
{
    return CompareArrays(rkM) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::operator>= (const Matrix3& rkM) const
{
    return CompareArrays(rkM) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> Matrix3<Real>::operator+ (const Matrix3& rkM) const
{
    return Matrix3<Real>(
        m_afEntry[0] + rkM.m_afEntry[0],
        m_afEntry[1] + rkM.m_afEntry[1],
        m_afEntry[2] + rkM.m_afEntry[2],
        m_afEntry[3] + rkM.m_afEntry[3],
        m_afEntry[4] + rkM.m_afEntry[4],
        m_afEntry[5] + rkM.m_afEntry[5],
        m_afEntry[6] + rkM.m_afEntry[6],
        m_afEntry[7] + rkM.m_afEntry[7],
        m_afEntry[8] + rkM.m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> Matrix3<Real>::operator- (const Matrix3& rkM) const
{
    return Matrix3<Real>(
        m_afEntry[0] - rkM.m_afEntry[0],
        m_afEntry[1] - rkM.m_afEntry[1],
        m_afEntry[2] - rkM.m_afEntry[2],
        m_afEntry[3] - rkM.m_afEntry[3],
        m_afEntry[4] - rkM.m_afEntry[4],
        m_afEntry[5] - rkM.m_afEntry[5],
        m_afEntry[6] - rkM.m_afEntry[6],
        m_afEntry[7] - rkM.m_afEntry[7],
        m_afEntry[8] - rkM.m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> Matrix3<Real>::operator* (const Matrix3& rkM) const
{
    return Matrix3<Real>(
        m_afEntry[0]*rkM.m_afEntry[0] +
        m_afEntry[1]*rkM.m_afEntry[3] +
        m_afEntry[2]*rkM.m_afEntry[6],

        m_afEntry[0]*rkM.m_afEntry[1] +
        m_afEntry[1]*rkM.m_afEntry[4] +
        m_afEntry[2]*rkM.m_afEntry[7],

        m_afEntry[0]*rkM.m_afEntry[2] +
        m_afEntry[1]*rkM.m_afEntry[5] +
        m_afEntry[2]*rkM.m_afEntry[8],

        m_afEntry[3]*rkM.m_afEntry[0] +
        m_afEntry[4]*rkM.m_afEntry[3] +
        m_afEntry[5]*rkM.m_afEntry[6],

        m_afEntry[3]*rkM.m_afEntry[1] +
        m_afEntry[4]*rkM.m_afEntry[4] +
        m_afEntry[5]*rkM.m_afEntry[7],

        m_afEntry[3]*rkM.m_afEntry[2] +
        m_afEntry[4]*rkM.m_afEntry[5] +
        m_afEntry[5]*rkM.m_afEntry[8],

        m_afEntry[6]*rkM.m_afEntry[0] +
        m_afEntry[7]*rkM.m_afEntry[3] +
        m_afEntry[8]*rkM.m_afEntry[6],

        m_afEntry[6]*rkM.m_afEntry[1] +
        m_afEntry[7]*rkM.m_afEntry[4] +
        m_afEntry[8]*rkM.m_afEntry[7],

        m_afEntry[6]*rkM.m_afEntry[2] +
        m_afEntry[7]*rkM.m_afEntry[5] +
        m_afEntry[8]*rkM.m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> Matrix3<Real>::operator* (Real fScalar) const
{
    return Matrix3<Real>(
        fScalar*m_afEntry[0],
        fScalar*m_afEntry[1],
        fScalar*m_afEntry[2],
        fScalar*m_afEntry[3],
        fScalar*m_afEntry[4],
        fScalar*m_afEntry[5],
        fScalar*m_afEntry[6],
        fScalar*m_afEntry[7],
        fScalar*m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> Matrix3<Real>::operator/ (Real fScalar) const
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        return Matrix3<Real>(
            fInvScalar*m_afEntry[0],
            fInvScalar*m_afEntry[1],
            fInvScalar*m_afEntry[2],
            fInvScalar*m_afEntry[3],
            fInvScalar*m_afEntry[4],
            fInvScalar*m_afEntry[5],
            fInvScalar*m_afEntry[6],
            fInvScalar*m_afEntry[7],
            fInvScalar*m_afEntry[8]);
    }

    return Matrix3<Real>(
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> Matrix3<Real>::operator- () const
{
    return Matrix3<Real>(
        -m_afEntry[0],
        -m_afEntry[1],
        -m_afEntry[2],
        -m_afEntry[3],
        -m_afEntry[4],
        -m_afEntry[5],
        -m_afEntry[6],
        -m_afEntry[7],
        -m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>& Matrix3<Real>::operator+= (const Matrix3& rkM)
{
    m_afEntry[0] += rkM.m_afEntry[0];
    m_afEntry[1] += rkM.m_afEntry[1];
    m_afEntry[2] += rkM.m_afEntry[2];
    m_afEntry[3] += rkM.m_afEntry[3];
    m_afEntry[4] += rkM.m_afEntry[4];
    m_afEntry[5] += rkM.m_afEntry[5];
    m_afEntry[6] += rkM.m_afEntry[6];
    m_afEntry[7] += rkM.m_afEntry[7];
    m_afEntry[8] += rkM.m_afEntry[8];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>& Matrix3<Real>::operator-= (const Matrix3& rkM)
{
    m_afEntry[0] -= rkM.m_afEntry[0];
    m_afEntry[1] -= rkM.m_afEntry[1];
    m_afEntry[2] -= rkM.m_afEntry[2];
    m_afEntry[3] -= rkM.m_afEntry[3];
    m_afEntry[4] -= rkM.m_afEntry[4];
    m_afEntry[5] -= rkM.m_afEntry[5];
    m_afEntry[6] -= rkM.m_afEntry[6];
    m_afEntry[7] -= rkM.m_afEntry[7];
    m_afEntry[8] -= rkM.m_afEntry[8];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>& Matrix3<Real>::operator*= (Real fScalar)
{
    m_afEntry[0] *= fScalar;
    m_afEntry[1] *= fScalar;
    m_afEntry[2] *= fScalar;
    m_afEntry[3] *= fScalar;
    m_afEntry[4] *= fScalar;
    m_afEntry[5] *= fScalar;
    m_afEntry[6] *= fScalar;
    m_afEntry[7] *= fScalar;
    m_afEntry[8] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real>& Matrix3<Real>::operator/= (Real fScalar)
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        m_afEntry[0] *= fInvScalar;
        m_afEntry[1] *= fInvScalar;
        m_afEntry[2] *= fInvScalar;
        m_afEntry[3] *= fInvScalar;
        m_afEntry[4] *= fInvScalar;
        m_afEntry[5] *= fInvScalar;
        m_afEntry[6] *= fInvScalar;
        m_afEntry[7] *= fInvScalar;
        m_afEntry[8] *= fInvScalar;
    }
    else
    {
        m_afEntry[0] = Math<Real>::MAX_REAL;
        m_afEntry[1] = Math<Real>::MAX_REAL;
        m_afEntry[2] = Math<Real>::MAX_REAL;
        m_afEntry[3] = Math<Real>::MAX_REAL;
        m_afEntry[4] = Math<Real>::MAX_REAL;
        m_afEntry[5] = Math<Real>::MAX_REAL;
        m_afEntry[6] = Math<Real>::MAX_REAL;
        m_afEntry[7] = Math<Real>::MAX_REAL;
        m_afEntry[8] = Math<Real>::MAX_REAL;
    }

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> Matrix3<Real>::operator* (const Vector3<Real>& rkV) const
{
    return Vector3<Real>(
        m_afEntry[0]*rkV[0] + m_afEntry[1]*rkV[1] + m_afEntry[2]*rkV[2],
        m_afEntry[3]*rkV[0] + m_afEntry[4]*rkV[1] + m_afEntry[5]*rkV[2],
        m_afEntry[6]*rkV[0] + m_afEntry[7]*rkV[1] + m_afEntry[8]*rkV[2]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::Transpose () const
{
    return Matrix3<Real>(
        m_afEntry[0],
        m_afEntry[3],
        m_afEntry[6],
        m_afEntry[1],
        m_afEntry[4],
        m_afEntry[7],
        m_afEntry[2],
        m_afEntry[5],
        m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::TransposeTimes (const Matrix3& rkM) const
{
    // P = A^T*B
    return Matrix3<Real>(
        m_afEntry[0]*rkM.m_afEntry[0] +
        m_afEntry[3]*rkM.m_afEntry[3] +
        m_afEntry[6]*rkM.m_afEntry[6],

        m_afEntry[0]*rkM.m_afEntry[1] +
        m_afEntry[3]*rkM.m_afEntry[4] +
        m_afEntry[6]*rkM.m_afEntry[7],

        m_afEntry[0]*rkM.m_afEntry[2] +
        m_afEntry[3]*rkM.m_afEntry[5] +
        m_afEntry[6]*rkM.m_afEntry[8],

        m_afEntry[1]*rkM.m_afEntry[0] +
        m_afEntry[4]*rkM.m_afEntry[3] +
        m_afEntry[7]*rkM.m_afEntry[6],

        m_afEntry[1]*rkM.m_afEntry[1] +
        m_afEntry[4]*rkM.m_afEntry[4] +
        m_afEntry[7]*rkM.m_afEntry[7],

        m_afEntry[1]*rkM.m_afEntry[2] +
        m_afEntry[4]*rkM.m_afEntry[5] +
        m_afEntry[7]*rkM.m_afEntry[8],

        m_afEntry[2]*rkM.m_afEntry[0] +
        m_afEntry[5]*rkM.m_afEntry[3] +
        m_afEntry[8]*rkM.m_afEntry[6],

        m_afEntry[2]*rkM.m_afEntry[1] +
        m_afEntry[5]*rkM.m_afEntry[4] +
        m_afEntry[8]*rkM.m_afEntry[7],

        m_afEntry[2]*rkM.m_afEntry[2] +
        m_afEntry[5]*rkM.m_afEntry[5] +
        m_afEntry[8]*rkM.m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::TimesTranspose (const Matrix3& rkM) const
{
    // P = A*B^T
    return Matrix3<Real>(
        m_afEntry[0]*rkM.m_afEntry[0] +
        m_afEntry[1]*rkM.m_afEntry[1] +
        m_afEntry[2]*rkM.m_afEntry[2],

        m_afEntry[0]*rkM.m_afEntry[3] +
        m_afEntry[1]*rkM.m_afEntry[4] +
        m_afEntry[2]*rkM.m_afEntry[5],

        m_afEntry[0]*rkM.m_afEntry[6] +
        m_afEntry[1]*rkM.m_afEntry[7] +
        m_afEntry[2]*rkM.m_afEntry[8],

        m_afEntry[3]*rkM.m_afEntry[0] +
        m_afEntry[4]*rkM.m_afEntry[1] +
        m_afEntry[5]*rkM.m_afEntry[2],

        m_afEntry[3]*rkM.m_afEntry[3] +
        m_afEntry[4]*rkM.m_afEntry[4] +
        m_afEntry[5]*rkM.m_afEntry[5],

        m_afEntry[3]*rkM.m_afEntry[6] +
        m_afEntry[4]*rkM.m_afEntry[7] +
        m_afEntry[5]*rkM.m_afEntry[8],

        m_afEntry[6]*rkM.m_afEntry[0] +
        m_afEntry[7]*rkM.m_afEntry[1] +
        m_afEntry[8]*rkM.m_afEntry[2],

        m_afEntry[6]*rkM.m_afEntry[3] +
        m_afEntry[7]*rkM.m_afEntry[4] +
        m_afEntry[8]*rkM.m_afEntry[5],

        m_afEntry[6]*rkM.m_afEntry[6] +
        m_afEntry[7]*rkM.m_afEntry[7] +
        m_afEntry[8]*rkM.m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::Inverse () const
{
    // Invert a 3x3 using cofactors.  This is faster than using a generic
    // Gaussian elimination because of the loop overhead of such a method.

    Matrix3 kInverse;

    kInverse.m_afEntry[0] =
        m_afEntry[4]*m_afEntry[8] - m_afEntry[5]*m_afEntry[7];
    kInverse.m_afEntry[1] =
        m_afEntry[2]*m_afEntry[7] - m_afEntry[1]*m_afEntry[8];
    kInverse.m_afEntry[2] =
        m_afEntry[1]*m_afEntry[5] - m_afEntry[2]*m_afEntry[4];
    kInverse.m_afEntry[3] =
        m_afEntry[5]*m_afEntry[6] - m_afEntry[3]*m_afEntry[8];
    kInverse.m_afEntry[4] =
        m_afEntry[0]*m_afEntry[8] - m_afEntry[2]*m_afEntry[6];
    kInverse.m_afEntry[5] =
        m_afEntry[2]*m_afEntry[3] - m_afEntry[0]*m_afEntry[5];
    kInverse.m_afEntry[6] =
        m_afEntry[3]*m_afEntry[7] - m_afEntry[4]*m_afEntry[6];
    kInverse.m_afEntry[7] =
        m_afEntry[1]*m_afEntry[6] - m_afEntry[0]*m_afEntry[7];
    kInverse.m_afEntry[8] =
        m_afEntry[0]*m_afEntry[4] - m_afEntry[1]*m_afEntry[3];

    Real fDet =
        m_afEntry[0]*kInverse.m_afEntry[0] +
        m_afEntry[1]*kInverse.m_afEntry[3] +
        m_afEntry[2]*kInverse.m_afEntry[6];

    if (Math<Real>::FAbs(fDet) <= Math<Real>::INVERSE_TOLERANCE)
    {
        return Matrix3<Real>(true);
    }

    Real fInvDet = ((Real)1.0)/fDet;
    kInverse.m_afEntry[0] *= fInvDet;
    kInverse.m_afEntry[1] *= fInvDet;
    kInverse.m_afEntry[2] *= fInvDet;
    kInverse.m_afEntry[3] *= fInvDet;
    kInverse.m_afEntry[4] *= fInvDet;
    kInverse.m_afEntry[5] *= fInvDet;
    kInverse.m_afEntry[6] *= fInvDet;
    kInverse.m_afEntry[7] *= fInvDet;
    kInverse.m_afEntry[8] *= fInvDet;
    return kInverse;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::Adjoint () const
{
    return Matrix3<Real>(
        m_afEntry[4]*m_afEntry[8] - m_afEntry[5]*m_afEntry[7],
        m_afEntry[2]*m_afEntry[7] - m_afEntry[1]*m_afEntry[8],
        m_afEntry[1]*m_afEntry[5] - m_afEntry[2]*m_afEntry[4],
        m_afEntry[5]*m_afEntry[6] - m_afEntry[3]*m_afEntry[8],
        m_afEntry[0]*m_afEntry[8] - m_afEntry[2]*m_afEntry[6],
        m_afEntry[2]*m_afEntry[3] - m_afEntry[0]*m_afEntry[5],
        m_afEntry[3]*m_afEntry[7] - m_afEntry[4]*m_afEntry[6],
        m_afEntry[1]*m_afEntry[6] - m_afEntry[0]*m_afEntry[7],
        m_afEntry[0]*m_afEntry[4] - m_afEntry[1]*m_afEntry[3]);
}
//----------------------------------------------------------------------------
template <class Real>
Real Matrix3<Real>::Determinant () const
{
    Real fCo00 = m_afEntry[4]*m_afEntry[8] - m_afEntry[5]*m_afEntry[7];
    Real fCo10 = m_afEntry[5]*m_afEntry[6] - m_afEntry[3]*m_afEntry[8];
    Real fCo20 = m_afEntry[3]*m_afEntry[7] - m_afEntry[4]*m_afEntry[6];
    Real fDet = m_afEntry[0]*fCo00 + m_afEntry[1]*fCo10 + m_afEntry[2]*fCo20;
    return fDet;
}
//----------------------------------------------------------------------------
template <class Real>
Real Matrix3<Real>::QForm (const Vector3<Real>& rkU, const Vector3<Real>& rkV)
    const
{
    return rkU.Dot((*this)*rkV);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::TimesDiagonal (const Vector3<Real>& rkDiag) const
{
    return Matrix3<Real>(
        m_afEntry[0]*rkDiag[0],m_afEntry[1]*rkDiag[1],m_afEntry[2]*rkDiag[2],
        m_afEntry[3]*rkDiag[0],m_afEntry[4]*rkDiag[1],m_afEntry[5]*rkDiag[2],
        m_afEntry[6]*rkDiag[0],m_afEntry[7]*rkDiag[1],m_afEntry[8]*rkDiag[2]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real> Matrix3<Real>::DiagonalTimes (const Vector3<Real>& rkDiag) const
{
    return Matrix3<Real>(
        rkDiag[0]*m_afEntry[0],rkDiag[0]*m_afEntry[1],rkDiag[0]*m_afEntry[2],
        rkDiag[1]*m_afEntry[3],rkDiag[1]*m_afEntry[4],rkDiag[1]*m_afEntry[5],
        rkDiag[2]*m_afEntry[6],rkDiag[2]*m_afEntry[7],rkDiag[2]*m_afEntry[8]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::ToAxisAngle (Vector3<Real>& rkAxis, Real& rfAngle) const
{
    // Let (x,y,z) be the unit-length axis and let A be an angle of rotation.
    // The rotation matrix is R = I + sin(A)*P + (1-cos(A))*P^2 where
    // I is the identity and
    //
    //       +-        -+
    //   P = |  0 -z +y |
    //       | +z  0 -x |
    //       | -y +x  0 |
    //       +-        -+
    //
    // If A > 0, R represents a counterclockwise rotation about the axis in
    // the sense of looking from the tip of the axis vector towards the
    // origin.  Some algebra will show that
    //
    //   cos(A) = (trace(R)-1)/2  and  R - R^t = 2*sin(A)*P
    //
    // In the event that A = pi, R-R^t = 0 which prevents us from extracting
    // the axis through P.  Instead note that R = I+2*P^2 when A = pi, so
    // P^2 = (R-I)/2.  The diagonal entries of P^2 are x^2-1, y^2-1, and
    // z^2-1.  We can solve these for axis (x,y,z).  Because the angle is pi,
    // it does not matter which sign you choose on the square roots.

    Real fTrace = m_afEntry[0] + m_afEntry[4] + m_afEntry[8];
    Real fCos = ((Real)0.5)*(fTrace-(Real)1.0);
    rfAngle = Math<Real>::ACos(fCos);  // in [0,PI]

    if (rfAngle > (Real)0.0)
    {
        if (rfAngle < Math<Real>::PI)
        {
            rkAxis[0] = m_afEntry[7]-m_afEntry[5];
            rkAxis[1] = m_afEntry[2]-m_afEntry[6];
            rkAxis[2] = m_afEntry[3]-m_afEntry[1];
            rkAxis.Normalize();
        }
        else
        {
            // angle is PI
            Real fHalfInverse;
            if (m_afEntry[0] >= m_afEntry[4])
            {
                // r00 >= r11
                if (m_afEntry[0] >= m_afEntry[8])
                {
                    // r00 is maximum diagonal term
                    rkAxis[0] = ((Real)0.5)*Math<Real>::Sqrt(m_afEntry[0] -
                        m_afEntry[4] - m_afEntry[8] + (Real)1.0);
                    fHalfInverse = ((Real)0.5)/rkAxis[0];
                    rkAxis[1] = fHalfInverse*m_afEntry[1];
                    rkAxis[2] = fHalfInverse*m_afEntry[2];
                }
                else
                {
                    // r22 is maximum diagonal term
                    rkAxis[2] = ((Real)0.5)*Math<Real>::Sqrt(m_afEntry[8] -
                        m_afEntry[0] - m_afEntry[4] + (Real)1.0);
                    fHalfInverse = ((Real)0.5)/rkAxis[2];
                    rkAxis[0] = fHalfInverse*m_afEntry[2];
                    rkAxis[1] = fHalfInverse*m_afEntry[5];
                }
            }
            else
            {
                // r11 > r00
                if (m_afEntry[4] >= m_afEntry[8])
                {
                    // r11 is maximum diagonal term
                    rkAxis[1] = ((Real)0.5)*Math<Real>::Sqrt(m_afEntry[4] -
                        m_afEntry[0] - m_afEntry[8] + (Real)1.0);
                    fHalfInverse  = ((Real)0.5)/rkAxis[1];
                    rkAxis[0] = fHalfInverse*m_afEntry[1];
                    rkAxis[2] = fHalfInverse*m_afEntry[5];
                }
                else
                {
                    // r22 is maximum diagonal term
                    rkAxis[2] = ((Real)0.5)*Math<Real>::Sqrt(m_afEntry[8] -
                        m_afEntry[0] - m_afEntry[4] + (Real)1.0);
                    fHalfInverse = ((Real)0.5)/rkAxis[2];
                    rkAxis[0] = fHalfInverse*m_afEntry[2];
                    rkAxis[1] = fHalfInverse*m_afEntry[5];
                }
            }
        }
    }
    else
    {
        // The angle is 0 and the matrix is the identity.  Any axis will
        // work, so just use the x-axis.
        rkAxis[0] = (Real)1.0;
        rkAxis[1] = (Real)0.0;
        rkAxis[2] = (Real)0.0;
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::Orthonormalize ()
{
    // Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is
    // M = [m0|m1|m2], then orthonormal output matrix is Q = [q0|q1|q2],
    //
    //   q0 = m0/|m0|
    //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
    //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
    //
    // where |V| indicates length of vector V and A*B indicates dot
    // product of vectors A and B.

    // compute q0
    Real fInvLength = Math<Real>::InvSqrt(m_afEntry[0]*m_afEntry[0] +
        m_afEntry[3]*m_afEntry[3] + m_afEntry[6]*m_afEntry[6]);

    m_afEntry[0] *= fInvLength;
    m_afEntry[3] *= fInvLength;
    m_afEntry[6] *= fInvLength;

    // compute q1
    Real fDot0 = m_afEntry[0]*m_afEntry[1] + m_afEntry[3]*m_afEntry[4] +
        m_afEntry[6]*m_afEntry[7];

    m_afEntry[1] -= fDot0*m_afEntry[0];
    m_afEntry[4] -= fDot0*m_afEntry[3];
    m_afEntry[7] -= fDot0*m_afEntry[6];

    fInvLength = Math<Real>::InvSqrt(m_afEntry[1]*m_afEntry[1] +
        m_afEntry[4]*m_afEntry[4] + m_afEntry[7]*m_afEntry[7]);

    m_afEntry[1] *= fInvLength;
    m_afEntry[4] *= fInvLength;
    m_afEntry[7] *= fInvLength;

    // compute q2
    Real fDot1 = m_afEntry[1]*m_afEntry[2] + m_afEntry[4]*m_afEntry[5] +
        m_afEntry[7]*m_afEntry[8];

    fDot0 = m_afEntry[0]*m_afEntry[2] + m_afEntry[3]*m_afEntry[5] +
        m_afEntry[6]*m_afEntry[8];

    m_afEntry[2] -= fDot0*m_afEntry[0] + fDot1*m_afEntry[1];
    m_afEntry[5] -= fDot0*m_afEntry[3] + fDot1*m_afEntry[4];
    m_afEntry[8] -= fDot0*m_afEntry[6] + fDot1*m_afEntry[7];

    fInvLength = Math<Real>::InvSqrt(m_afEntry[2]*m_afEntry[2] +
        m_afEntry[5]*m_afEntry[5] + m_afEntry[8]*m_afEntry[8]);

    m_afEntry[2] *= fInvLength;
    m_afEntry[5] *= fInvLength;
    m_afEntry[8] *= fInvLength;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::EigenDecomposition (Matrix3& rkRot, Matrix3& rkDiag) const
{
    // Factor M = R*D*R^T.  The columns of R are the eigenvectors.  The
    // diagonal entries of D are the corresponding eigenvalues.
    Real afDiag[3], afSubd[2];
    rkRot = *this;
    bool bReflection = rkRot.Tridiagonalize(afDiag,afSubd);
    bool bConverged = rkRot.QLAlgorithm(afDiag,afSubd);
    assert(bConverged);

    // (insertion) sort eigenvalues in increasing order, d0 <= d1 <= d2
    int i;
    Real fSave;

    if (afDiag[1] < afDiag[0])
    {
        // swap d0 and d1
        fSave = afDiag[0];
        afDiag[0] = afDiag[1];
        afDiag[1] = fSave;

        // swap V0 and V1
        for (i = 0; i < 3; i++)
        {
            fSave = rkRot[i][0];
            rkRot[i][0] = rkRot[i][1];
            rkRot[i][1] = fSave;
        }
        bReflection = !bReflection;
    }

    if (afDiag[2] < afDiag[1])
    {
        // swap d1 and d2
        fSave = afDiag[1];
        afDiag[1] = afDiag[2];
        afDiag[2] = fSave;

        // swap V1 and V2
        for (i = 0; i < 3; i++)
        {
            fSave = rkRot[i][1];
            rkRot[i][1] = rkRot[i][2];
            rkRot[i][2] = fSave;
        }
        bReflection = !bReflection;
    }

    if (afDiag[1] < afDiag[0])
    {
        // swap d0 and d1
        fSave = afDiag[0];
        afDiag[0] = afDiag[1];
        afDiag[1] = fSave;

        // swap V0 and V1
        for (i = 0; i < 3; i++)
        {
            fSave = rkRot[i][0];
            rkRot[i][0] = rkRot[i][1];
            rkRot[i][1] = fSave;
        }
        bReflection = !bReflection;
    }

    rkDiag.MakeDiagonal(afDiag[0],afDiag[1],afDiag[2]);

    if (bReflection)
    {
        // The orthogonal transformation that diagonalizes M is a reflection.
        // Make the eigenvectors a right--handed system by changing sign on
        // the last column.
        rkRot[0][2] = -rkRot[0][2];
        rkRot[1][2] = -rkRot[1][2];
        rkRot[2][2] = -rkRot[2][2];
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromEulerAnglesXYZ (Real fXAngle, Real fYAngle,
    Real fZAngle)
{
    Real fCos, fSin;

    fCos = Math<Real>::Cos(fXAngle);
    fSin = Math<Real>::Sin(fXAngle);
    Matrix3 kXMat(
        (Real)1.0,(Real)0.0,(Real)0.0,
        (Real)0.0,fCos,-fSin,
        (Real)0.0,fSin,fCos);

    fCos = Math<Real>::Cos(fYAngle);
    fSin = Math<Real>::Sin(fYAngle);
    Matrix3 kYMat(
        fCos,(Real)0.0,fSin,
        (Real)0.0,(Real)1.0,(Real)0.0,
        -fSin,(Real)0.0,fCos);

    fCos = Math<Real>::Cos(fZAngle);
    fSin = Math<Real>::Sin(fZAngle);
    Matrix3 kZMat(
        fCos,-fSin,(Real)0.0,
        fSin,fCos,(Real)0.0,
        (Real)0.0,(Real)0.0,(Real)1.0);

    *this = kXMat*(kYMat*kZMat);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromEulerAnglesXZY (Real fXAngle, Real fZAngle,
    Real fYAngle)
{
    Real fCos, fSin;

    fCos = Math<Real>::Cos(fXAngle);
    fSin = Math<Real>::Sin(fXAngle);
    Matrix3 kXMat(
        (Real)1.0,(Real)0.0,(Real)0.0,
        (Real)0.0,fCos,-fSin,
        (Real)0.0,fSin,fCos);

    fCos = Math<Real>::Cos(fZAngle);
    fSin = Math<Real>::Sin(fZAngle);
    Matrix3 kZMat(
        fCos,-fSin,(Real)0.0,
        fSin,fCos,(Real)0.0,
        (Real)0.0,(Real)0.0,(Real)1.0);

    fCos = Math<Real>::Cos(fYAngle);
    fSin = Math<Real>::Sin(fYAngle);
    Matrix3 kYMat(
        fCos,(Real)0.0,fSin,
        (Real)0.0,(Real)1.0,(Real)0.0,
        -fSin,(Real)0.0,fCos);

    *this = kXMat*(kZMat*kYMat);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromEulerAnglesYXZ (Real fYAngle, Real fXAngle,
    Real fZAngle)
{
    Real fCos, fSin;

    fCos = Math<Real>::Cos(fYAngle);
    fSin = Math<Real>::Sin(fYAngle);
    Matrix3 kYMat(
        fCos,(Real)0.0,fSin,
        (Real)0.0,(Real)1.0,(Real)0.0,
        -fSin,(Real)0.0,fCos);

    fCos = Math<Real>::Cos(fXAngle);
    fSin = Math<Real>::Sin(fXAngle);
    Matrix3 kXMat(
        (Real)1.0,(Real)0.0,(Real)0.0,
        (Real)0.0,fCos,-fSin,
        (Real)0.0,fSin,fCos);

    fCos = Math<Real>::Cos(fZAngle);
    fSin = Math<Real>::Sin(fZAngle);
    Matrix3 kZMat(
        fCos,-fSin,(Real)0.0,
        fSin,fCos,(Real)0.0,
        (Real)0.0,(Real)0.0,(Real)1.0);

    *this = kYMat*(kXMat*kZMat);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromEulerAnglesYZX (Real fYAngle, Real fZAngle,
    Real fXAngle)
{
    Real fCos, fSin;


    fCos = Math<Real>::Cos(fYAngle);
    fSin = Math<Real>::Sin(fYAngle);
    Matrix3 kYMat(
        fCos,(Real)0.0,fSin,
        (Real)0.0,(Real)1.0,(Real)0.0,
        -fSin,(Real)0.0,fCos);

    fCos = Math<Real>::Cos(fZAngle);
    fSin = Math<Real>::Sin(fZAngle);
    Matrix3 kZMat(
        fCos,-fSin,(Real)0.0,
        fSin,fCos,(Real)0.0,
        (Real)0.0,(Real)0.0,(Real)1.0);

    fCos = Math<Real>::Cos(fXAngle);
    fSin = Math<Real>::Sin(fXAngle);
    Matrix3 kXMat(
        (Real)1.0,(Real)0.0,(Real)0.0,
        (Real)0.0,fCos,-fSin,
        (Real)0.0,fSin,fCos);

    *this = kYMat*(kZMat*kXMat);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromEulerAnglesZXY (Real fZAngle, Real fXAngle,
    Real fYAngle)
{
    Real fCos, fSin;

    fCos = Math<Real>::Cos(fZAngle);
    fSin = Math<Real>::Sin(fZAngle);
    Matrix3 kZMat(
        fCos,-fSin,(Real)0.0,
        fSin,fCos,(Real)0.0,
        (Real)0.0,(Real)0.0,(Real)1.0);

    fCos = Math<Real>::Cos(fXAngle);
    fSin = Math<Real>::Sin(fXAngle);
    Matrix3 kXMat(
        (Real)1.0,(Real)0.0,(Real)0.0,
        (Real)0.0,fCos,-fSin,
        (Real)0.0,fSin,fCos);

    fCos = Math<Real>::Cos(fYAngle);
    fSin = Math<Real>::Sin(fYAngle);
    Matrix3 kYMat(
        fCos,(Real)0.0,fSin,
        (Real)0.0,(Real)1.0,(Real)0.0,
        -fSin,(Real)0.0,fCos);

    *this = kZMat*(kXMat*kYMat);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::FromEulerAnglesZYX (Real fZAngle, Real fYAngle,
    Real fXAngle)
{
    Real fCos, fSin;

    fCos = Math<Real>::Cos(fZAngle);
    fSin = Math<Real>::Sin(fZAngle);
    Matrix3 kZMat(
        fCos,-fSin,(Real)0.0,
        fSin,fCos,(Real)0.0,
        (Real)0.0,(Real)0.0,(Real)1.0);

    fCos = Math<Real>::Cos(fYAngle);
    fSin = Math<Real>::Sin(fYAngle);
    Matrix3 kYMat(
        fCos,(Real)0.0,fSin,
        (Real)0.0,(Real)1.0,(Real)0.0,
        -fSin,(Real)0.0,fCos);

    fCos = Math<Real>::Cos(fXAngle);
    fSin = Math<Real>::Sin(fXAngle);
    Matrix3 kXMat(
        (Real)1.0,(Real)0.0,(Real)0.0,
        (Real)0.0,fCos,-fSin,
        (Real)0.0,fSin,fCos);

    *this = kZMat*(kYMat*kXMat);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::ToEulerAnglesXYZ (Real& rfXAngle, Real& rfYAngle,
    Real& rfZAngle) const
{
    // rot =  cy*cz          -cy*sz           sy
    //        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
    //       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

    if (m_afEntry[2] < (Real)1.0)
    {
        if (m_afEntry[2] > -(Real)1.0)
        {
            rfXAngle = Math<Real>::ATan2(-m_afEntry[5],m_afEntry[8]);
            rfYAngle = (Real)asin((double)m_afEntry[2]);
            rfZAngle = Math<Real>::ATan2(-m_afEntry[1],m_afEntry[0]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  XA - ZA = -atan2(r10,r11)
            rfXAngle = -Math<Real>::ATan2(m_afEntry[3],m_afEntry[4]);
            rfYAngle = -Math<Real>::HALF_PI;
            rfZAngle = (Real)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  XAngle + ZAngle = atan2(r10,r11)
        rfXAngle = Math<Real>::ATan2(m_afEntry[3],m_afEntry[4]);
        rfYAngle = Math<Real>::HALF_PI;
        rfZAngle = (Real)0.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::ToEulerAnglesXZY (Real& rfXAngle, Real& rfZAngle,
    Real& rfYAngle) const
{
    // rot =  cy*cz          -sz              cz*sy
    //        sx*sy+cx*cy*sz  cx*cz          -cy*sx+cx*sy*sz
    //       -cx*sy+cy*sx*sz  cz*sx           cx*cy+sx*sy*sz

    if (m_afEntry[1] < (Real)1.0)
    {
        if (m_afEntry[1] > -(Real)1.0)
        {
            rfXAngle = Math<Real>::ATan2(m_afEntry[7],m_afEntry[4]);
            rfZAngle = (Real)asin(-(double)m_afEntry[1]);
            rfYAngle = Math<Real>::ATan2(m_afEntry[2],m_afEntry[0]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  XA - YA = atan2(r20,r22)
            rfXAngle = Math<Real>::ATan2(m_afEntry[6],m_afEntry[8]);
            rfZAngle = Math<Real>::HALF_PI;
            rfYAngle = (Real)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  XA + YA = atan2(-r20,r22)
        rfXAngle = Math<Real>::ATan2(-m_afEntry[6],m_afEntry[8]);
        rfZAngle = -Math<Real>::HALF_PI;
        rfYAngle = (Real)0.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::ToEulerAnglesYXZ (Real& rfYAngle, Real& rfXAngle,
    Real& rfZAngle) const
{
    // rot =  cy*cz+sx*sy*sz  cz*sx*sy-cy*sz  cx*sy
    //        cx*sz           cx*cz          -sx
    //       -cz*sy+cy*sx*sz  cy*cz*sx+sy*sz  cx*cy

    if (m_afEntry[5] < (Real)1.0)
    {
        if (m_afEntry[5] > -(Real)1.0)
        {
            rfYAngle = Math<Real>::ATan2(m_afEntry[2],m_afEntry[8]);
            rfXAngle = (Real)asin(-(double)m_afEntry[5]);
            rfZAngle = Math<Real>::ATan2(m_afEntry[3],m_afEntry[4]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  YA - ZA = atan2(r01,r00)
            rfYAngle = Math<Real>::ATan2(m_afEntry[1],m_afEntry[0]);
            rfXAngle = Math<Real>::HALF_PI;
            rfZAngle = (Real)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  YA + ZA = atan2(-r01,r00)
        rfYAngle = Math<Real>::ATan2(-m_afEntry[1],m_afEntry[0]);
        rfXAngle = -Math<Real>::HALF_PI;
        rfZAngle = (Real)0.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::ToEulerAnglesYZX (Real& rfYAngle, Real& rfZAngle,
    Real& rfXAngle) const
{
    // rot =  cy*cz           sx*sy-cx*cy*sz  cx*sy+cy*sx*sz
    //        sz              cx*cz          -cz*sx
    //       -cz*sy           cy*sx+cx*sy*sz  cx*cy-sx*sy*sz

    if (m_afEntry[3] < (Real)1.0)
    {
        if (m_afEntry[3] > -(Real)1.0)
        {
            rfYAngle = Math<Real>::ATan2(-m_afEntry[6],m_afEntry[0]);
            rfZAngle = (Real)asin((double)m_afEntry[3]);
            rfXAngle = Math<Real>::ATan2(-m_afEntry[5],m_afEntry[4]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  YA - XA = -atan2(r21,r22);
            rfYAngle = -Math<Real>::ATan2(m_afEntry[7],m_afEntry[8]);
            rfZAngle = -Math<Real>::HALF_PI;
            rfXAngle = (Real)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  YA + XA = atan2(r21,r22)
        rfYAngle = Math<Real>::ATan2(m_afEntry[7],m_afEntry[8]);
        rfZAngle = Math<Real>::HALF_PI;
        rfXAngle = (Real)0.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::ToEulerAnglesZXY (Real& rfZAngle, Real& rfXAngle,
    Real& rfYAngle) const
{
    // rot =  cy*cz-sx*sy*sz -cx*sz           cz*sy+cy*sx*sz
    //        cz*sx*sy+cy*sz  cx*cz          -cy*cz*sx+sy*sz
    //       -cx*sy           sx              cx*cy

    if (m_afEntry[7] < (Real)1.0)
    {
        if (m_afEntry[7] > -(Real)1.0)
        {
            rfZAngle = Math<Real>::ATan2(-m_afEntry[1],m_afEntry[4]);
            rfXAngle = (Real)asin((double)m_afEntry[7]);
            rfYAngle = Math<Real>::ATan2(-m_afEntry[6],m_afEntry[8]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  ZA - YA = -atan(r02,r00)
            rfZAngle = -Math<Real>::ATan2(m_afEntry[2],m_afEntry[0]);
            rfXAngle = -Math<Real>::HALF_PI;
            rfYAngle = (Real)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  ZA + YA = atan2(r02,r00)
        rfZAngle = Math<Real>::ATan2(m_afEntry[2],m_afEntry[0]);
        rfXAngle = Math<Real>::HALF_PI;
        rfYAngle = (Real)0.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::ToEulerAnglesZYX (Real& rfZAngle, Real& rfYAngle,
    Real& rfXAngle) const
{
    // rot =  cy*cz           cz*sx*sy-cx*sz  cx*cz*sy+sx*sz
    //        cy*sz           cx*cz+sx*sy*sz -cz*sx+cx*sy*sz
    //       -sy              cy*sx           cx*cy

    if (m_afEntry[6] < (Real)1.0)
    {
        if (m_afEntry[6] > -(Real)1.0)
        {
            rfZAngle = Math<Real>::ATan2(m_afEntry[3],m_afEntry[0]);
            rfYAngle = (Real)asin(-(double)m_afEntry[6]);
            rfXAngle = Math<Real>::ATan2(m_afEntry[7],m_afEntry[8]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  ZA - XA = -atan2(r01,r02)
            rfZAngle = -Math<Real>::ATan2(m_afEntry[1],m_afEntry[2]);
            rfYAngle = Math<Real>::HALF_PI;
            rfXAngle = (Real)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  ZA + XA = atan2(-r01,-r02)
        rfZAngle = Math<Real>::ATan2(-m_afEntry[1],-m_afEntry[2]);
        rfYAngle = -Math<Real>::HALF_PI;
        rfXAngle = (Real)0.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix3<Real>& Matrix3<Real>::Slerp (Real fT, const Matrix3& rkR0,
    const Matrix3& rkR1)
{
    Vector3<Real> kAxis;
    Real fAngle;
    Matrix3 kProd = rkR0.TransposeTimes(rkR1);
    kProd.ToAxisAngle(kAxis,fAngle);
    FromAxisAngle(kAxis,fT*fAngle);
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::Tridiagonalize (Real afDiag[3], Real afSubd[2])
{
    // Householder reduction T = Q^t M Q
    //   Input:   
    //     mat, symmetric 3x3 matrix M
    //   Output:  
    //     mat, orthogonal matrix Q (a reflection)
    //     diag, diagonal entries of T
    //     subd, subdiagonal entries of T (T is symmetric)

    Real fM00 = m_afEntry[0];
    Real fM01 = m_afEntry[1];
    Real fM02 = m_afEntry[2];
    Real fM11 = m_afEntry[4];
    Real fM12 = m_afEntry[5];
    Real fM22 = m_afEntry[8];

    afDiag[0] = fM00;
    if (Math<Real>::FAbs(fM02) >= Math<Real>::ZERO_TOLERANCE)
    {
        afSubd[0] = Math<Real>::Sqrt(fM01*fM01+fM02*fM02);
        Real fInvLength = ((Real)1.0)/afSubd[0];
        fM01 *= fInvLength;
        fM02 *= fInvLength;
        Real fTmp = ((Real)2.0)*fM01*fM12+fM02*(fM22-fM11);
        afDiag[1] = fM11+fM02*fTmp;
        afDiag[2] = fM22-fM02*fTmp;
        afSubd[1] = fM12-fM01*fTmp;

        m_afEntry[0] = (Real)1.0;
        m_afEntry[1] = (Real)0.0;
        m_afEntry[2] = (Real)0.0;
        m_afEntry[3] = (Real)0.0;
        m_afEntry[4] = fM01;
        m_afEntry[5] = fM02;
        m_afEntry[6] = (Real)0.0;
        m_afEntry[7] = fM02;
        m_afEntry[8] = -fM01;
        return true;
    }
    else
    {
        afDiag[1] = fM11;
        afDiag[2] = fM22;
        afSubd[0] = fM01;
        afSubd[1] = fM12;

        m_afEntry[0] = (Real)1.0;
        m_afEntry[1] = (Real)0.0;
        m_afEntry[2] = (Real)0.0;
        m_afEntry[3] = (Real)0.0;
        m_afEntry[4] = (Real)1.0;
        m_afEntry[5] = (Real)0.0;
        m_afEntry[6] = (Real)0.0;
        m_afEntry[7] = (Real)0.0;
        m_afEntry[8] = (Real)1.0;
        return false;
    }
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix3<Real>::QLAlgorithm (Real afDiag[3], Real afSubd[2])
{
    // This is an implementation of the symmetric QR algorithm from the book
    // "Matrix Computations" by Gene H. Golub and Charles F. Van Loan, second
    // edition.  The algorithm is 8.2.3.  The implementation has a slight
    // variation to actually make it a QL algorithm, and it traps the case
    // when either of the subdiagonal terms s0 or s1 is zero and reduces the
    // 2-by-2 subblock directly.

    const int iMax = 32;
    for (int i = 0; i < iMax; i++)
    {
        Real fSum, fDiff, fDiscr, fEValue0, fEValue1, fCos, fSin, fTmp;
        int iRow;

        fSum = Math<Real>::FAbs(afDiag[0]) + Math<Real>::FAbs(afDiag[1]);
        if (Math<Real>::FAbs(afSubd[0]) + fSum == fSum)
        {
            // The matrix is effectively
            //       +-        -+
            //   M = | d0  0  0 |
            //       | 0  d1 s1 |
            //       | 0  s1 d2 |
            //       +-        -+

            // Compute the eigenvalues as roots of a quadratic equation.
            fSum = afDiag[1] + afDiag[2];
            fDiff = afDiag[1] - afDiag[2];
            fDiscr = Math<Real>::Sqrt(fDiff*fDiff +
                ((Real)4.0)*afSubd[1]*afSubd[1]);
            fEValue0 = ((Real)0.5)*(fSum - fDiscr);
            fEValue1 = ((Real)0.5)*(fSum + fDiscr);

            // Compute the Givens rotation.
            if (fDiff >= (Real)0.0)
            {
                fCos = afSubd[1];
                fSin = afDiag[1] - fEValue0;
            }
            else
            {
                fCos = afDiag[2] - fEValue0;
                fSin = afSubd[1];
            }
            fTmp = Math<Real>::InvSqrt(fCos*fCos + fSin*fSin);
            fCos *= fTmp;
            fSin *= fTmp;

            // Postmultiply the current orthogonal matrix with the Givens
            // rotation.
            for (iRow = 0; iRow < 3; iRow++)
            {
                fTmp = m_afEntry[2+3*iRow];
                m_afEntry[2+3*iRow] = fSin*m_afEntry[1+3*iRow] + fCos*fTmp;
                m_afEntry[1+3*iRow] = fCos*m_afEntry[1+3*iRow] - fSin*fTmp;
            }

            // Update the tridiagonal matrix.
            afDiag[1] = fEValue0;
            afDiag[2] = fEValue1;
            afSubd[0] = (Real)0.0;
            afSubd[1] = (Real)0.0;
            return true;
        }

        fSum = Math<Real>::FAbs(afDiag[1]) + Math<Real>::FAbs(afDiag[2]);
        if (Math<Real>::FAbs(afSubd[1]) + fSum == fSum)
        {
            // The matrix is effectively
            //       +-         -+
            //   M = | d0  s0  0 |
            //       | s0  d1  0 |
            //       | 0   0  d2 |
            //       +-         -+

            // Compute the eigenvalues as roots of a quadratic equation.
            fSum = afDiag[0] + afDiag[1];
            fDiff = afDiag[0] - afDiag[1];
            fDiscr = Math<Real>::Sqrt(fDiff*fDiff +
                ((Real)4.0)*afSubd[0]*afSubd[0]);
            fEValue0 = ((Real)0.5)*(fSum - fDiscr);
            fEValue1 = ((Real)0.5)*(fSum + fDiscr);

            // Compute the Givens rotation.
            if (fDiff >= (Real)0.0)
            {
                fCos = afSubd[0];
                fSin = afDiag[0] - fEValue0;
            }
            else
            {
                fCos = afDiag[1] - fEValue0;
                fSin = afSubd[0];
            }
            fTmp = Math<Real>::InvSqrt(fCos*fCos + fSin*fSin);
            fCos *= fTmp;
            fSin *= fTmp;

            // Postmultiply the current orthogonal matrix with the Givens
            // rotation.
            for (iRow = 0; iRow < 3; iRow++)
            {
                fTmp = m_afEntry[1+3*iRow];
                m_afEntry[1+3*iRow] = fSin*m_afEntry[0+3*iRow] + fCos*fTmp;
                m_afEntry[0+3*iRow] = fCos*m_afEntry[0+3*iRow] - fSin*fTmp;
            }

            // Update the tridiagonal matrix.
            afDiag[0] = fEValue0;
            afDiag[1] = fEValue1;
            afSubd[0] = (Real)0.0;
            afSubd[1] = (Real)0.0;
            return true;
        }

        // The matrix is
        //       +-        -+
        //   M = | d0 s0  0 |
        //       | s0 d1 s1 |
        //       | 0  s1 d2 |
        //       +-        -+

        // Set up the parameters for the first pass of the QL step.  The
        // value for A is the difference between diagonal term D[2] and the
        // implicit shift suggested by Wilkinson.
        Real fRatio = (afDiag[1]-afDiag[0])/(((Real)2.0f)*afSubd[0]);
        Real fRoot = Math<Real>::Sqrt((Real)1.0 + fRatio*fRatio);
        Real fB = afSubd[1];
        Real fA = afDiag[2] - afDiag[0];
        if (fRatio >= (Real)0.0)
        {
            fA += afSubd[0]/(fRatio + fRoot);
        }
        else
        {
            fA += afSubd[0]/(fRatio - fRoot);
        }

        // Compute the Givens rotation for the first pass.
        if (Math<Real>::FAbs(fB) >= Math<Real>::FAbs(fA))
        {
            fRatio = fA/fB;
            fSin = Math<Real>::InvSqrt((Real)1.0 + fRatio*fRatio);
            fCos = fRatio*fSin;
        }
        else
        {
            fRatio = fB/fA;
            fCos = Math<Real>::InvSqrt((Real)1.0 + fRatio*fRatio);
            fSin = fRatio*fCos;
        }

        // Postmultiply the current orthogonal matrix with the Givens
        // rotation.
        for (iRow = 0; iRow < 3; iRow++)
        {
            fTmp = m_afEntry[2+3*iRow];
            m_afEntry[2+3*iRow] = fSin*m_afEntry[1+3*iRow]+fCos*fTmp;
            m_afEntry[1+3*iRow] = fCos*m_afEntry[1+3*iRow]-fSin*fTmp;
        }

        // Set up the parameters for the second pass of the QL step.  The
        // values tmp0 and tmp1 are required to fully update the tridiagonal
        // matrix at the end of the second pass.
        Real fTmp0 = (afDiag[1] - afDiag[2])*fSin +
            ((Real)2.0)*afSubd[1]*fCos;
        Real fTmp1 = fCos*afSubd[0];
        fB = fSin*afSubd[0];
        fA = fCos*fTmp0 - afSubd[1];
        fTmp0 *= fSin;

        // Compute the Givens rotation for the second pass.  The subdiagonal
        // term S[1] in the tridiagonal matrix is updated at this time.
        if (Math<Real>::FAbs(fB) >= Math<Real>::FAbs(fA))
        {
            fRatio = fA/fB;
            fRoot = Math<Real>::Sqrt((Real)1.0 + fRatio*fRatio);
            afSubd[1] = fB*fRoot;
            fSin = ((Real)1.0)/fRoot;
            fCos = fRatio*fSin;
        }
        else
        {
            fRatio = fB/fA;
            fRoot = Math<Real>::Sqrt((Real)1.0 + fRatio*fRatio);
            afSubd[1] = fA*fRoot;
            fCos = ((Real)1.0)/fRoot;
            fSin = fRatio*fCos;
        }

        // Postmultiply the current orthogonal matrix with the Givens
        // rotation.
        for (iRow = 0; iRow < 3; iRow++)
        {
            fTmp = m_afEntry[1+3*iRow];
            m_afEntry[1+3*iRow] = fSin*m_afEntry[0+3*iRow]+fCos*fTmp;
            m_afEntry[0+3*iRow] = fCos*m_afEntry[0+3*iRow]-fSin*fTmp;
        }

        // Complete the update of the tridiagonal matrix.
        Real fTmp2 = afDiag[1] - fTmp0;
        afDiag[2] += fTmp0;
        fTmp0 = (afDiag[0] - fTmp2)*fSin + ((Real)2.0)*fTmp1*fCos;
        afSubd[0] = fCos*fTmp0 - fTmp1;
        fTmp0 *= fSin;
        afDiag[1] = fTmp2 + fTmp0;
        afDiag[0] -= fTmp0;
    }
    return false;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::Bidiagonalize (Matrix3& rkA, Matrix3& rkL, Matrix3& rkR)
{
    Real afV[3], afW[3];
    Real fLength, fSign, fT1, fInvT1, fT2;
    bool bIdentity;

    // map first column to (*,0,0)
    fLength = Math<Real>::Sqrt(rkA[0][0]*rkA[0][0] + rkA[1][0]*rkA[1][0] +
        rkA[2][0]*rkA[2][0]);
    if (fLength > (Real)0.0)
    {
        fSign = (rkA[0][0] > (Real)0.0 ? (Real)1.0 : -(Real)1.0);
        fT1 = rkA[0][0] + fSign*fLength;
        fInvT1 = ((Real)1.0)/fT1;
        afV[1] = rkA[1][0]*fInvT1;
        afV[2] = rkA[2][0]*fInvT1;

        fT2 = -((Real)2.0)/(((Real)1.0)+afV[1]*afV[1]+afV[2]*afV[2]);
        afW[0] = fT2*(rkA[0][0]+rkA[1][0]*afV[1]+rkA[2][0]*afV[2]);
        afW[1] = fT2*(rkA[0][1]+rkA[1][1]*afV[1]+rkA[2][1]*afV[2]);
        afW[2] = fT2*(rkA[0][2]+rkA[1][2]*afV[1]+rkA[2][2]*afV[2]);
        rkA[0][0] += afW[0];
        rkA[0][1] += afW[1];
        rkA[0][2] += afW[2];
        rkA[1][1] += afV[1]*afW[1];
        rkA[1][2] += afV[1]*afW[2];
        rkA[2][1] += afV[2]*afW[1];
        rkA[2][2] += afV[2]*afW[2];

        rkL[0][0] = ((Real)1.0)+fT2;
        rkL[0][1] = fT2*afV[1];
        rkL[1][0] = rkL[0][1];
        rkL[0][2] = fT2*afV[2];
        rkL[2][0] = rkL[0][2];
        rkL[1][1] = ((Real)1.0)+fT2*afV[1]*afV[1];
        rkL[1][2] = fT2*afV[1]*afV[2];
        rkL[2][1] = rkL[1][2];
        rkL[2][2] = ((Real)1.0)+fT2*afV[2]*afV[2];
        bIdentity = false;
    }
    else
    {
        rkL = Matrix3<Real>(false);
        bIdentity = true;
    }

    // map first row to (*,*,0)
    fLength = Math<Real>::Sqrt(rkA[0][1]*rkA[0][1]+rkA[0][2]*rkA[0][2]);
    if (fLength > (Real)0.0)
    {
        fSign = (rkA[0][1] > (Real)0.0 ? (Real)1.0 : -(Real)1.0);
        fT1 = rkA[0][1] + fSign*fLength;
        afV[2] = rkA[0][2]/fT1;

        fT2 = -((Real)2.0)/(((Real)1.0)+afV[2]*afV[2]);
        afW[0] = fT2*(rkA[0][1]+rkA[0][2]*afV[2]);
        afW[1] = fT2*(rkA[1][1]+rkA[1][2]*afV[2]);
        afW[2] = fT2*(rkA[2][1]+rkA[2][2]*afV[2]);
        rkA[0][1] += afW[0];
        rkA[1][1] += afW[1];
        rkA[1][2] += afW[1]*afV[2];
        rkA[2][1] += afW[2];
        rkA[2][2] += afW[2]*afV[2];
        
        rkR[0][0] = (Real)1.0;
        rkR[0][1] = (Real)0.0;
        rkR[1][0] = (Real)0.0;
        rkR[0][2] = (Real)0.0;
        rkR[2][0] = (Real)0.0;
        rkR[1][1] = ((Real)1.0)+fT2;
        rkR[1][2] = fT2*afV[2];
        rkR[2][1] = rkR[1][2];
        rkR[2][2] = ((Real)1.0)+fT2*afV[2]*afV[2];
    }
    else
    {
        rkR = Matrix3<Real>(false);
    }

    // map second column to (*,*,0)
    fLength = Math<Real>::Sqrt(rkA[1][1]*rkA[1][1]+rkA[2][1]*rkA[2][1]);
    if (fLength > (Real)0.0)
    {
        fSign = (rkA[1][1] > (Real)0.0 ? (Real)1.0 : -(Real)1.0);
        fT1 = rkA[1][1] + fSign*fLength;
        afV[2] = rkA[2][1]/fT1;

        fT2 = -((Real)2.0)/(((Real)1.0)+afV[2]*afV[2]);
        afW[1] = fT2*(rkA[1][1]+rkA[2][1]*afV[2]);
        afW[2] = fT2*(rkA[1][2]+rkA[2][2]*afV[2]);
        rkA[1][1] += afW[1];
        rkA[1][2] += afW[2];
        rkA[2][2] += afV[2]*afW[2];

        Real fA = ((Real)1.0)+fT2;
        Real fB = fT2*afV[2];
        Real fC = ((Real)1.0)+fB*afV[2];

        if (bIdentity)
        {
            rkL[0][0] = (Real)1.0;
            rkL[0][1] = (Real)0.0;
            rkL[1][0] = (Real)0.0;
            rkL[0][2] = (Real)0.0;
            rkL[2][0] = (Real)0.0;
            rkL[1][1] = fA;
            rkL[1][2] = fB;
            rkL[2][1] = fB;
            rkL[2][2] = fC;
        }
        else
        {
            for (int iRow = 0; iRow < 3; iRow++)
            {
                Real fTmp0 = rkL[iRow][1];
                Real fTmp1 = rkL[iRow][2];
                rkL[iRow][1] = fA*fTmp0+fB*fTmp1;
                rkL[iRow][2] = fB*fTmp0+fC*fTmp1;
            }
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::GolubKahanStep (Matrix3& rkA, Matrix3& rkL, Matrix3& rkR)
{
    Real fT11 = rkA[0][1]*rkA[0][1]+rkA[1][1]*rkA[1][1];
    Real fT22 = rkA[1][2]*rkA[1][2]+rkA[2][2]*rkA[2][2];
    Real fT12 = rkA[1][1]*rkA[1][2];
    Real fTrace = fT11+fT22;
    Real fDiff = fT11-fT22;
    Real fDiscr = Math<Real>::Sqrt(fDiff*fDiff+((Real)4.0)*fT12*fT12);
    Real fRoot1 = ((Real)0.5)*(fTrace+fDiscr);
    Real fRoot2 = ((Real)0.5)*(fTrace-fDiscr);

    // adjust right
    Real fY = rkA[0][0] - (Math<Real>::FAbs(fRoot1-fT22) <=
        Math<Real>::FAbs(fRoot2-fT22) ? fRoot1 : fRoot2);
    Real fZ = rkA[0][1];
    Real fInvLength = Math<Real>::InvSqrt(fY*fY+fZ*fZ);
    Real fSin = fZ*fInvLength;
    Real fCos = -fY*fInvLength;

    Real fTmp0 = rkA[0][0];
    Real fTmp1 = rkA[0][1];
    rkA[0][0] = fCos*fTmp0-fSin*fTmp1;
    rkA[0][1] = fSin*fTmp0+fCos*fTmp1;
    rkA[1][0] = -fSin*rkA[1][1];
    rkA[1][1] *= fCos;

    int iRow;
    for (iRow = 0; iRow < 3; iRow++)
    {
        fTmp0 = rkR[0][iRow];
        fTmp1 = rkR[1][iRow];
        rkR[0][iRow] = fCos*fTmp0-fSin*fTmp1;
        rkR[1][iRow] = fSin*fTmp0+fCos*fTmp1;
    }

    // adjust left
    fY = rkA[0][0];
    fZ = rkA[1][0];
    fInvLength = Math<Real>::InvSqrt(fY*fY+fZ*fZ);
    fSin = fZ*fInvLength;
    fCos = -fY*fInvLength;

    rkA[0][0] = fCos*rkA[0][0]-fSin*rkA[1][0];
    fTmp0 = rkA[0][1];
    fTmp1 = rkA[1][1];
    rkA[0][1] = fCos*fTmp0-fSin*fTmp1;
    rkA[1][1] = fSin*fTmp0+fCos*fTmp1;
    rkA[0][2] = -fSin*rkA[1][2];
    rkA[1][2] *= fCos;

    int iCol;
    for (iCol = 0; iCol < 3; iCol++)
    {
        fTmp0 = rkL[iCol][0];
        fTmp1 = rkL[iCol][1];
        rkL[iCol][0] = fCos*fTmp0-fSin*fTmp1;
        rkL[iCol][1] = fSin*fTmp0+fCos*fTmp1;
    }

    // adjust right
    fY = rkA[0][1];
    fZ = rkA[0][2];
    fInvLength = Math<Real>::InvSqrt(fY*fY+fZ*fZ);
    fSin = fZ*fInvLength;
    fCos = -fY*fInvLength;

    rkA[0][1] = fCos*rkA[0][1]-fSin*rkA[0][2];
    fTmp0 = rkA[1][1];
    fTmp1 = rkA[1][2];
    rkA[1][1] = fCos*fTmp0-fSin*fTmp1;
    rkA[1][2] = fSin*fTmp0+fCos*fTmp1;
    rkA[2][1] = -fSin*rkA[2][2];
    rkA[2][2] *= fCos;

    for (iRow = 0; iRow < 3; iRow++)
    {
        fTmp0 = rkR[1][iRow];
        fTmp1 = rkR[2][iRow];
        rkR[1][iRow] = fCos*fTmp0-fSin*fTmp1;
        rkR[2][iRow] = fSin*fTmp0+fCos*fTmp1;
    }

    // adjust left
    fY = rkA[1][1];
    fZ = rkA[2][1];
    fInvLength = Math<Real>::InvSqrt(fY*fY+fZ*fZ);
    fSin = fZ*fInvLength;
    fCos = -fY*fInvLength;

    rkA[1][1] = fCos*rkA[1][1]-fSin*rkA[2][1];
    fTmp0 = rkA[1][2];
    fTmp1 = rkA[2][2];
    rkA[1][2] = fCos*fTmp0-fSin*fTmp1;
    rkA[2][2] = fSin*fTmp0+fCos*fTmp1;

    for (iCol = 0; iCol < 3; iCol++)
    {
        fTmp0 = rkL[iCol][1];
        fTmp1 = rkL[iCol][2];
        rkL[iCol][1] = fCos*fTmp0-fSin*fTmp1;
        rkL[iCol][2] = fSin*fTmp0+fCos*fTmp1;
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::SingularValueDecomposition (Matrix3& rkL, Matrix3& rkS,
    Matrix3& rkR) const
{
    int iRow, iCol;

    Matrix3 kA = *this;
    Bidiagonalize(kA,rkL,rkR);
    rkS.MakeZero();

    const int iMax = 32;
    const Real fEpsilon = (Real)1e-04;
    for (int i = 0; i < iMax; i++)
    {
        Real fTmp, fTmp0, fTmp1;
        Real fSin0, fCos0, fTan0;
        Real fSin1, fCos1, fTan1;

        bool bTest1 = (Math<Real>::FAbs(kA[0][1]) <=
            fEpsilon*(Math<Real>::FAbs(kA[0][0]) +
            Math<Real>::FAbs(kA[1][1])));
        bool bTest2 = (Math<Real>::FAbs(kA[1][2]) <=
            fEpsilon*(Math<Real>::FAbs(kA[1][1]) +
            Math<Real>::FAbs(kA[2][2])));
        if (bTest1)
        {
            if (bTest2)
            {
                rkS[0][0] = kA[0][0];
                rkS[1][1] = kA[1][1];
                rkS[2][2] = kA[2][2];
                break;
            }
            else
            {
                // 2x2 closed form factorization
                fTmp = (kA[1][1]*kA[1][1] - kA[2][2]*kA[2][2] +
                    kA[1][2]*kA[1][2])/(kA[1][2]*kA[2][2]);
                fTan0 = ((Real)0.5)*(fTmp + Math<Real>::Sqrt(fTmp*fTmp +
                    ((Real)4.0)));
                fCos0 = Math<Real>::InvSqrt(((Real)1.0)+fTan0*fTan0);
                fSin0 = fTan0*fCos0;

                for (iCol = 0; iCol < 3; iCol++)
                {
                    fTmp0 = rkL[iCol][1];
                    fTmp1 = rkL[iCol][2];
                    rkL[iCol][1] = fCos0*fTmp0-fSin0*fTmp1;
                    rkL[iCol][2] = fSin0*fTmp0+fCos0*fTmp1;
                }
                
                fTan1 = (kA[1][2]-kA[2][2]*fTan0)/kA[1][1];
                fCos1 = Math<Real>::InvSqrt(((Real)1.0)+fTan1*fTan1);
                fSin1 = -fTan1*fCos1;

                for (iRow = 0; iRow < 3; iRow++)
                {
                    fTmp0 = rkR[1][iRow];
                    fTmp1 = rkR[2][iRow];
                    rkR[1][iRow] = fCos1*fTmp0-fSin1*fTmp1;
                    rkR[2][iRow] = fSin1*fTmp0+fCos1*fTmp1;
                }

                rkS[0][0] = kA[0][0];
                rkS[1][1] = fCos0*fCos1*kA[1][1] -
                    fSin1*(fCos0*kA[1][2]-fSin0*kA[2][2]);
                rkS[2][2] = fSin0*fSin1*kA[1][1] +
                    fCos1*(fSin0*kA[1][2]+fCos0*kA[2][2]);
                break;
            }
        }
        else 
        {
            if (bTest2)
            {
                // 2x2 closed form factorization 
                fTmp = (kA[0][0]*kA[0][0] + kA[1][1]*kA[1][1] -
                    kA[0][1]*kA[0][1])/(kA[0][1]*kA[1][1]);
                fTan0 = ((Real)0.5)*(-fTmp + Math<Real>::Sqrt(fTmp*fTmp +
                    ((Real)4.0)));
                fCos0 = Math<Real>::InvSqrt(((Real)1.0)+fTan0*fTan0);
                fSin0 = fTan0*fCos0;

                for (iCol = 0; iCol < 3; iCol++)
                {
                    fTmp0 = rkL[iCol][0];
                    fTmp1 = rkL[iCol][1];
                    rkL[iCol][0] = fCos0*fTmp0-fSin0*fTmp1;
                    rkL[iCol][1] = fSin0*fTmp0+fCos0*fTmp1;
                }
                
                fTan1 = (kA[0][1]-kA[1][1]*fTan0)/kA[0][0];
                fCos1 = Math<Real>::InvSqrt(((Real)1.0)+fTan1*fTan1);
                fSin1 = -fTan1*fCos1;

                for (iRow = 0; iRow < 3; iRow++)
                {
                    fTmp0 = rkR[0][iRow];
                    fTmp1 = rkR[1][iRow];
                    rkR[0][iRow] = fCos1*fTmp0-fSin1*fTmp1;
                    rkR[1][iRow] = fSin1*fTmp0+fCos1*fTmp1;
                }

                rkS[0][0] = fCos0*fCos1*kA[0][0] -
                    fSin1*(fCos0*kA[0][1]-fSin0*kA[1][1]);
                rkS[1][1] = fSin0*fSin1*kA[0][0] +
                    fCos1*(fSin0*kA[0][1]+fCos0*kA[1][1]);
                rkS[2][2] = kA[2][2];
                break;
            }
            else
            {
                GolubKahanStep(kA,rkL,rkR);
            }
        }
    }

    // positize diagonal
    for (iRow = 0; iRow < 3; iRow++)
    {
        if (rkS[iRow][iRow] < (Real)0.0)
        {
            rkS[iRow][iRow] = -rkS[iRow][iRow];
            for (iCol = 0; iCol < 3; iCol++)
                rkR[iRow][iCol] = -rkR[iRow][iCol];
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::SingularValueComposition (const Matrix3& rkL,
    const Matrix3& rkS, const Matrix3& rkR)
{
    *this = rkL*(rkS*rkR);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::PolarDecomposition (Matrix3& rkQ, Matrix3& rkS)
{
	// Decompose M = L*D*R^T.
	Matrix3 kL, kD, kRTranspose;
	SingularValueDecomposition(kL,kD,kRTranspose);

	// Compute Q = L*R^T.
	rkQ = kL*kRTranspose;

	// Compute S = R*D*R^T.
	rkS = kRTranspose.TransposeTimes(kD*kRTranspose);

	// Numerical round-off errors can cause S not to be symmetric in the
	// sense that S[i][j] and S[j][i] are slightly different for i != j.
	// Correct this by averaging S and Transpose(S).
	rkS[0][1] = ((Real)0.5)*(rkS[0][1] + rkS[1][0]);
	rkS[1][0] = rkS[0][1];
	rkS[0][2] = ((Real)0.5)*(rkS[0][2] + rkS[2][0]);
	rkS[2][0] = rkS[0][2];
	rkS[1][2] = ((Real)0.5)*(rkS[1][2] + rkS[2][1]);
	rkS[2][1] = rkS[1][2];
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix3<Real>::QDUDecomposition (Matrix3& rkQ, Matrix3& rkD,

    Matrix3& rkU) const
{
    // Factor M = QR = QDU where Q is orthogonal (rotation), D is diagonal
    // (scaling),  and U is upper triangular with ones on its diagonal
    // (shear).  Algorithm uses Gram-Schmidt orthogonalization (the QR
    // algorithm).
    //
    // If M = [ m0 | m1 | m2 ] and Q = [ q0 | q1 | q2 ], then
    //
    //   q0 = m0/|m0|
    //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
    //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
    //
    // where |V| indicates length of vector V and A*B indicates dot
    // product of vectors A and B.  The matrix R has entries
    //
    //   r00 = q0*m0  r01 = q0*m1  r02 = q0*m2
    //   r10 = 0      r11 = q1*m1  r12 = q1*m2
    //   r20 = 0      r21 = 0      r22 = q2*m2
    //
    // so D = diag(r00,r11,r22) and U has entries u01 = r01/r00,
    // u02 = r02/r00, and u12 = r12/r11.

    // build orthogonal matrix Q
    Real fInvLength = Math<Real>::InvSqrt(m_afEntry[0]*m_afEntry[0] +
        m_afEntry[3]*m_afEntry[3] + m_afEntry[6]*m_afEntry[6]);
    rkQ[0][0] = m_afEntry[0]*fInvLength;
    rkQ[1][0] = m_afEntry[3]*fInvLength;
    rkQ[2][0] = m_afEntry[6]*fInvLength;

    Real fDot = rkQ[0][0]*m_afEntry[1] + rkQ[1][0]*m_afEntry[4] +
        rkQ[2][0]*m_afEntry[7];
    rkQ[0][1] = m_afEntry[1]-fDot*rkQ[0][0];
    rkQ[1][1] = m_afEntry[4]-fDot*rkQ[1][0];
    rkQ[2][1] = m_afEntry[7]-fDot*rkQ[2][0];
    fInvLength = Math<Real>::InvSqrt(rkQ[0][1]*rkQ[0][1] +
        rkQ[1][1]*rkQ[1][1] + rkQ[2][1]*rkQ[2][1]);
    rkQ[0][1] *= fInvLength;
    rkQ[1][1] *= fInvLength;
    rkQ[2][1] *= fInvLength;

    fDot = rkQ[0][0]*m_afEntry[2] + rkQ[1][0]*m_afEntry[5] +
        rkQ[2][0]*m_afEntry[8];
    rkQ[0][2] = m_afEntry[2]-fDot*rkQ[0][0];
    rkQ[1][2] = m_afEntry[5]-fDot*rkQ[1][0];
    rkQ[2][2] = m_afEntry[8]-fDot*rkQ[2][0];
    fDot = rkQ[0][1]*m_afEntry[2] + rkQ[1][1]*m_afEntry[5] +
        rkQ[2][1]*m_afEntry[8];
    rkQ[0][2] -= fDot*rkQ[0][1];
    rkQ[1][2] -= fDot*rkQ[1][1];
    rkQ[2][2] -= fDot*rkQ[2][1];
    fInvLength = Math<Real>::InvSqrt(rkQ[0][2]*rkQ[0][2] +
        rkQ[1][2]*rkQ[1][2] + rkQ[2][2]*rkQ[2][2]);
    rkQ[0][2] *= fInvLength;
    rkQ[1][2] *= fInvLength;
    rkQ[2][2] *= fInvLength;

    // guarantee that orthogonal matrix has determinant 1 (no reflections)
    Real fDet = rkQ[0][0]*rkQ[1][1]*rkQ[2][2] + rkQ[0][1]*rkQ[1][2]*rkQ[2][0]
        +  rkQ[0][2]*rkQ[1][0]*rkQ[2][1] - rkQ[0][2]*rkQ[1][1]*rkQ[2][0]
        -  rkQ[0][1]*rkQ[1][0]*rkQ[2][2] - rkQ[0][0]*rkQ[1][2]*rkQ[2][1];

    if (fDet < (Real)0.0)
    {
        for (int iRow = 0; iRow < 3; iRow++)
        {
            for (int iCol = 0; iCol < 3; iCol++)
            {
                rkQ[iRow][iCol] = -rkQ[iRow][iCol];
            }
        }
    }

    // build "right" matrix R
    Matrix3 kR;
    kR[0][0] = rkQ[0][0]*m_afEntry[0] + rkQ[1][0]*m_afEntry[3] +
        rkQ[2][0]*m_afEntry[6];
    kR[0][1] = rkQ[0][0]*m_afEntry[1] + rkQ[1][0]*m_afEntry[4] +
        rkQ[2][0]*m_afEntry[7];
    kR[1][1] = rkQ[0][1]*m_afEntry[1] + rkQ[1][1]*m_afEntry[4] +
        rkQ[2][1]*m_afEntry[7];
    kR[0][2] = rkQ[0][0]*m_afEntry[2] + rkQ[1][0]*m_afEntry[5] +
        rkQ[2][0]*m_afEntry[8];
    kR[1][2] = rkQ[0][1]*m_afEntry[2] + rkQ[1][1]*m_afEntry[5] +
        rkQ[2][1]*m_afEntry[8];
    kR[2][2] = rkQ[0][2]*m_afEntry[2] + rkQ[1][2]*m_afEntry[5] +
        rkQ[2][2]*m_afEntry[8];

    // the scaling component
    rkD.MakeDiagonal(kR[0][0],kR[1][1],kR[2][2]);

    // the shear component
    Real fInvD0 = ((Real)1.0)/rkD[0][0];
    rkU[0][0] = (Real)1.0;
    rkU[0][1] = kR[0][1]*fInvD0;
    rkU[0][2] = kR[0][2]*fInvD0;
    rkU[1][0] = (Real)0.0;
    rkU[1][1] = (Real)1.0;
    rkU[1][2] = kR[1][2]/rkD[1][1];
    rkU[2][0] = (Real)0.0;
    rkU[2][1] = (Real)0.0;
    rkU[2][2] = (Real)1.0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix3<Real> operator* (Real fScalar, const Matrix3<Real>& rkM)
{
    return rkM*fScalar;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector3<Real> operator* (const Vector3<Real>& rkV,
    const Matrix3<Real>& rkM)
{
    return Vector3<Real>(
        rkV[0]*rkM[0][0] + rkV[1]*rkM[1][0] + rkV[2]*rkM[2][0],
        rkV[0]*rkM[0][1] + rkV[1]*rkM[1][1] + rkV[2]*rkM[2][1],
        rkV[0]*rkM[0][2] + rkV[1]*rkM[1][2] + rkV[2]*rkM[2][2]);
}
//----------------------------------------------------------------------------


typedef Matrix3<float> Matrix3f;
typedef Matrix3<double> Matrix3d;

}








namespace Wm4
{

template <class Real>
class Matrix4
{
public:
    // If bZero is true, create the zero matrix.  Otherwise, create the
    // identity matrix.
    Matrix4 (bool bZero = true);

    // copy constructor
    Matrix4 (const Matrix4& rkM);

    // input Mrc is in row r, column c.
    Matrix4 (Real fM00, Real fM01, Real fM02, Real fM03,
             Real fM10, Real fM11, Real fM12, Real fM13,
             Real fM20, Real fM21, Real fM22, Real fM23,
             Real fM30, Real fM31, Real fM32, Real fM33);

    // Create a matrix from an array of numbers.  The input array is
    // interpreted based on the Boolean input as
    //   true:  entry[0..15]={m00,m01,m02,m03,m10,m11,m12,m13,m20,m21,m22,
    //                        m23,m30,m31,m32,m33} [row major]
    //   false: entry[0..15]={m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,
    //                        m32,m03,m13,m23,m33} [col major]
    Matrix4 (const Real afEntry[16], bool bRowMajor);

    void MakeZero ();
    void MakeIdentity ();

    // member access
    inline operator const Real* () const;
    inline operator Real* ();
    inline const Real* operator[] (int iRow) const;
    inline Real* operator[] (int iRow);
    inline Real operator() (int iRow, int iCol) const;
    inline Real& operator() (int iRow, int iCol);
    void SetRow (int iRow, const Vector4<Real>& rkV);
    Vector4<Real> GetRow (int iRow) const;
    void SetColumn (int iCol, const Vector4<Real>& rkV);
    Vector4<Real> GetColumn (int iCol) const;
    void GetColumnMajor (Real* afCMajor) const;

    // assignment
    inline Matrix4& operator= (const Matrix4& rkM);

    // comparison
    bool operator== (const Matrix4& rkM) const;
    bool operator!= (const Matrix4& rkM) const;
    bool operator<  (const Matrix4& rkM) const;
    bool operator<= (const Matrix4& rkM) const;
    bool operator>  (const Matrix4& rkM) const;
    bool operator>= (const Matrix4& rkM) const;

    // arithmetic operations
    inline Matrix4 operator+ (const Matrix4& rkM) const;
    inline Matrix4 operator- (const Matrix4& rkM) const;
    inline Matrix4 operator* (const Matrix4& rkM) const;
    inline Matrix4 operator* (Real fScalar) const;
    inline Matrix4 operator/ (Real fScalar) const;
    inline Matrix4 operator- () const;

    // arithmetic updates
    inline Matrix4& operator+= (const Matrix4& rkM);
    inline Matrix4& operator-= (const Matrix4& rkM);
    inline Matrix4& operator*= (Real fScalar);
    inline Matrix4& operator/= (Real fScalar);

    // matrix times vector
    inline Vector4<Real> operator* (const Vector4<Real>& rkV) const;  // M * v

    // other operations
    Matrix4 Transpose () const;  // M^T
    Matrix4 TransposeTimes (const Matrix4& rkM) const;  // this^T * M
    Matrix4 TimesTranspose (const Matrix4& rkM) const;  // this * M^T
    Matrix4 Inverse () const;
    Matrix4 Adjoint () const;
    Real Determinant () const;
    Real QForm (const Vector4<Real>& rkU,
		const Vector4<Real>& rkV) const;  // u^T*M*v    
	
	Matrix4& MakeTensorProduct (const Vector4<Real>& rkU,
	const Vector4<Real>& rkV);

    // projection matrices onto a specified plane
    void MakeObliqueProjection (const Vector3<Real>& rkNormal,
        const Vector3<Real>& rkPoint, const Vector3<Real>& rkDirection);
    void MakePerspectiveProjection (const Vector3<Real>& rkNormal,
        const Vector3<Real>& rkPoint, const Vector3<Real>& rkEye);

    // reflection matrix through a specified plane
    void MakeReflection (const Vector3<Real>& rkNormal,
        const Vector3<Real>& rkPoint);


private:
    // support for comparisons
    int CompareArrays (const Matrix4& rkM) const;

    Real m_afEntry[16];
};

// c * M
template <class Real>
inline Matrix4<Real> operator* (Real fScalar, const Matrix4<Real>& rkM);

// v^T * M
template <class Real>
inline Vector4<Real> operator* (const Vector4<Real>& rkV,
    const Matrix4<Real>& rkM);

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
Matrix4<Real>::Matrix4 (bool bZero)
{
    if (bZero)
    {
        MakeZero();
    }
    else
    {
        MakeIdentity();
    }
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real>::Matrix4 (const Matrix4& rkM)
{
    m_afEntry[ 0] = rkM.m_afEntry[ 0];
    m_afEntry[ 1] = rkM.m_afEntry[ 1];
    m_afEntry[ 2] = rkM.m_afEntry[ 2];
    m_afEntry[ 3] = rkM.m_afEntry[ 3];
    m_afEntry[ 4] = rkM.m_afEntry[ 4];
    m_afEntry[ 5] = rkM.m_afEntry[ 5];
    m_afEntry[ 6] = rkM.m_afEntry[ 6];
    m_afEntry[ 7] = rkM.m_afEntry[ 7];
    m_afEntry[ 8] = rkM.m_afEntry[ 8];
    m_afEntry[ 9] = rkM.m_afEntry[ 9];
    m_afEntry[10] = rkM.m_afEntry[10];
    m_afEntry[11] = rkM.m_afEntry[11];
    m_afEntry[12] = rkM.m_afEntry[12];
    m_afEntry[13] = rkM.m_afEntry[13];
    m_afEntry[14] = rkM.m_afEntry[14];
    m_afEntry[15] = rkM.m_afEntry[15];
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real>::Matrix4 (Real fM00, Real fM01, Real fM02, Real fM03,
    Real fM10, Real fM11, Real fM12, Real fM13, Real fM20, Real fM21,
    Real fM22, Real fM23, Real fM30, Real fM31, Real fM32, Real fM33)
{
    m_afEntry[ 0] = fM00;
    m_afEntry[ 1] = fM01;
    m_afEntry[ 2] = fM02;
    m_afEntry[ 3] = fM03;
    m_afEntry[ 4] = fM10;
    m_afEntry[ 5] = fM11;
    m_afEntry[ 6] = fM12;
    m_afEntry[ 7] = fM13;
    m_afEntry[ 8] = fM20;
    m_afEntry[ 9] = fM21;
    m_afEntry[10] = fM22;
    m_afEntry[11] = fM23;
    m_afEntry[12] = fM30;
    m_afEntry[13] = fM31;
    m_afEntry[14] = fM32;
    m_afEntry[15] = fM33;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real>::Matrix4 (const Real afEntry[16], bool bRowMajor)
{
    if (bRowMajor)
    {
        m_afEntry[ 0] = afEntry[ 0];
        m_afEntry[ 1] = afEntry[ 1];
        m_afEntry[ 2] = afEntry[ 2];
        m_afEntry[ 3] = afEntry[ 3];
        m_afEntry[ 4] = afEntry[ 4];
        m_afEntry[ 5] = afEntry[ 5];
        m_afEntry[ 6] = afEntry[ 6];
        m_afEntry[ 7] = afEntry[ 7];
        m_afEntry[ 8] = afEntry[ 8];
        m_afEntry[ 9] = afEntry[ 9];
        m_afEntry[10] = afEntry[10];
        m_afEntry[11] = afEntry[11];
        m_afEntry[12] = afEntry[12];
        m_afEntry[13] = afEntry[13];
        m_afEntry[14] = afEntry[14];
        m_afEntry[15] = afEntry[15];
    }
    else
    {
        m_afEntry[ 0] = afEntry[ 0];
        m_afEntry[ 1] = afEntry[ 4];
        m_afEntry[ 2] = afEntry[ 8];
        m_afEntry[ 3] = afEntry[12];
        m_afEntry[ 4] = afEntry[ 1];
        m_afEntry[ 5] = afEntry[ 5];
        m_afEntry[ 6] = afEntry[ 9];
        m_afEntry[ 7] = afEntry[13];
        m_afEntry[ 8] = afEntry[ 2];
        m_afEntry[ 9] = afEntry[ 6];
        m_afEntry[10] = afEntry[10];
        m_afEntry[11] = afEntry[14];
        m_afEntry[12] = afEntry[ 3];
        m_afEntry[13] = afEntry[ 7];
        m_afEntry[14] = afEntry[11];
        m_afEntry[15] = afEntry[15];
    }
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>::operator const Real* () const
{
    return m_afEntry;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>::operator Real* ()
{
    return m_afEntry;
}
//----------------------------------------------------------------------------
template <class Real>
inline const Real* Matrix4<Real>::operator[] (int iRow) const
{
    return &m_afEntry[4*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real* Matrix4<Real>::operator[] (int iRow)
{
    return &m_afEntry[4*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real Matrix4<Real>::operator() (int iRow, int iCol) const
{
    return m_afEntry[iCol+4*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
inline Real& Matrix4<Real>::operator() (int iRow, int iCol)
{
    return m_afEntry[iCol+4*iRow];
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::MakeZero ()
{
    m_afEntry[ 0] = (Real)0.0;
    m_afEntry[ 1] = (Real)0.0;
    m_afEntry[ 2] = (Real)0.0;
    m_afEntry[ 3] = (Real)0.0;
    m_afEntry[ 4] = (Real)0.0;
    m_afEntry[ 5] = (Real)0.0;
    m_afEntry[ 6] = (Real)0.0;
    m_afEntry[ 7] = (Real)0.0;
    m_afEntry[ 8] = (Real)0.0;
    m_afEntry[ 9] = (Real)0.0;
    m_afEntry[10] = (Real)0.0;
    m_afEntry[11] = (Real)0.0;
    m_afEntry[12] = (Real)0.0;
    m_afEntry[13] = (Real)0.0;
    m_afEntry[14] = (Real)0.0;
    m_afEntry[15] = (Real)0.0;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::MakeIdentity ()
{
    m_afEntry[ 0] = (Real)1.0;
    m_afEntry[ 1] = (Real)0.0;
    m_afEntry[ 2] = (Real)0.0;
    m_afEntry[ 3] = (Real)0.0;
    m_afEntry[ 4] = (Real)0.0;
    m_afEntry[ 5] = (Real)1.0;
    m_afEntry[ 6] = (Real)0.0;
    m_afEntry[ 7] = (Real)0.0;
    m_afEntry[ 8] = (Real)0.0;
    m_afEntry[ 9] = (Real)0.0;
    m_afEntry[10] = (Real)1.0;
    m_afEntry[11] = (Real)0.0;
    m_afEntry[12] = (Real)0.0;
    m_afEntry[13] = (Real)0.0;
    m_afEntry[14] = (Real)0.0;
    m_afEntry[15] = (Real)1.0;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::SetRow (int iRow, const Vector4<Real>& rkV)
{
    int i0 = 4*iRow, i1 = i0+1, i2 = i1+1, i3 = i2+1;
    m_afEntry[i0] = rkV[0];
    m_afEntry[i1] = rkV[1];
    m_afEntry[i2] = rkV[2];
    m_afEntry[i3] = rkV[3];
}
//----------------------------------------------------------------------------
template <class Real>
Vector4<Real> Matrix4<Real>::GetRow (int iRow) const
{
    int i0 = 4*iRow, i1 = i0+1, i2 = i1+1, i3 = i2+1;
    return Vector4<Real>(m_afEntry[i0],m_afEntry[i1],m_afEntry[i2],
        m_afEntry[i3]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::SetColumn (int iCol, const Vector4<Real>& rkV)
{
    m_afEntry[iCol] = rkV[0];
    m_afEntry[iCol+4] = rkV[1];
    m_afEntry[iCol+8] = rkV[2];
    m_afEntry[iCol+12] = rkV[3];
}
//----------------------------------------------------------------------------
template <class Real>
Vector4<Real> Matrix4<Real>::GetColumn (int iCol) const
{
    return Vector4<Real>(m_afEntry[iCol],m_afEntry[iCol+4],m_afEntry[iCol+8],
        m_afEntry[iCol+12]);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::GetColumnMajor (Real* afCMajor) const
{
    afCMajor[ 0] = m_afEntry[ 0];
    afCMajor[ 1] = m_afEntry[ 4];
    afCMajor[ 2] = m_afEntry[ 8];
    afCMajor[ 3] = m_afEntry[12];
    afCMajor[ 4] = m_afEntry[ 1];
    afCMajor[ 5] = m_afEntry[ 5];
    afCMajor[ 6] = m_afEntry[ 9];
    afCMajor[ 7] = m_afEntry[13];
    afCMajor[ 8] = m_afEntry[ 2];
    afCMajor[ 9] = m_afEntry[ 6];
    afCMajor[10] = m_afEntry[10];
    afCMajor[11] = m_afEntry[14];
    afCMajor[12] = m_afEntry[ 3];
    afCMajor[13] = m_afEntry[ 7];
    afCMajor[14] = m_afEntry[11];
    afCMajor[15] = m_afEntry[15];
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>& Matrix4<Real>::operator= (const Matrix4& rkM)
{
    m_afEntry[ 0] = rkM.m_afEntry[ 0];
    m_afEntry[ 1] = rkM.m_afEntry[ 1];
    m_afEntry[ 2] = rkM.m_afEntry[ 2];
    m_afEntry[ 3] = rkM.m_afEntry[ 3];
    m_afEntry[ 4] = rkM.m_afEntry[ 4];
    m_afEntry[ 5] = rkM.m_afEntry[ 5];
    m_afEntry[ 6] = rkM.m_afEntry[ 6];
    m_afEntry[ 7] = rkM.m_afEntry[ 7];
    m_afEntry[ 8] = rkM.m_afEntry[ 8];
    m_afEntry[ 9] = rkM.m_afEntry[ 9];
    m_afEntry[10] = rkM.m_afEntry[10];
    m_afEntry[11] = rkM.m_afEntry[11];
    m_afEntry[12] = rkM.m_afEntry[12];
    m_afEntry[13] = rkM.m_afEntry[13];
    m_afEntry[14] = rkM.m_afEntry[14];
    m_afEntry[15] = rkM.m_afEntry[15];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
int Matrix4<Real>::CompareArrays (const Matrix4& rkM) const
{
    return memcmp(m_afEntry,rkM.m_afEntry,16*sizeof(Real));
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix4<Real>::operator== (const Matrix4& rkM) const
{
    return CompareArrays(rkM) == 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix4<Real>::operator!= (const Matrix4& rkM) const
{
    return CompareArrays(rkM) != 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix4<Real>::operator<  (const Matrix4& rkM) const
{
    return CompareArrays(rkM) < 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix4<Real>::operator<= (const Matrix4& rkM) const
{
    return CompareArrays(rkM) <= 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix4<Real>::operator>  (const Matrix4& rkM) const
{
    return CompareArrays(rkM) > 0;
}
//----------------------------------------------------------------------------
template <class Real>
bool Matrix4<Real>::operator>= (const Matrix4& rkM) const
{
    return CompareArrays(rkM) >= 0;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real> Matrix4<Real>::operator+ (const Matrix4& rkM) const
{
    return Matrix4<Real>(
        m_afEntry[ 0] + rkM.m_afEntry[ 0],
        m_afEntry[ 1] + rkM.m_afEntry[ 1],
        m_afEntry[ 2] + rkM.m_afEntry[ 2],
        m_afEntry[ 3] + rkM.m_afEntry[ 3],
        m_afEntry[ 4] + rkM.m_afEntry[ 4],
        m_afEntry[ 5] + rkM.m_afEntry[ 5],
        m_afEntry[ 6] + rkM.m_afEntry[ 6],
        m_afEntry[ 7] + rkM.m_afEntry[ 7],
        m_afEntry[ 8] + rkM.m_afEntry[ 8],
        m_afEntry[ 9] + rkM.m_afEntry[ 9],
        m_afEntry[10] + rkM.m_afEntry[10],
        m_afEntry[11] + rkM.m_afEntry[11],
        m_afEntry[12] + rkM.m_afEntry[12],
        m_afEntry[13] + rkM.m_afEntry[13],
        m_afEntry[14] + rkM.m_afEntry[14],
        m_afEntry[15] + rkM.m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real> Matrix4<Real>::operator- (const Matrix4& rkM) const
{
    return Matrix4<Real>(
        m_afEntry[ 0] - rkM.m_afEntry[ 0],
        m_afEntry[ 1] - rkM.m_afEntry[ 1],
        m_afEntry[ 2] - rkM.m_afEntry[ 2],
        m_afEntry[ 3] - rkM.m_afEntry[ 3],
        m_afEntry[ 4] - rkM.m_afEntry[ 4],
        m_afEntry[ 5] - rkM.m_afEntry[ 5],
        m_afEntry[ 6] - rkM.m_afEntry[ 6],
        m_afEntry[ 7] - rkM.m_afEntry[ 7],
        m_afEntry[ 8] - rkM.m_afEntry[ 8],
        m_afEntry[ 9] - rkM.m_afEntry[ 9],
        m_afEntry[10] - rkM.m_afEntry[10],
        m_afEntry[11] - rkM.m_afEntry[11],
        m_afEntry[12] - rkM.m_afEntry[12],
        m_afEntry[13] - rkM.m_afEntry[13],
        m_afEntry[14] - rkM.m_afEntry[14],
        m_afEntry[15] - rkM.m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real> Matrix4<Real>::operator* (const Matrix4& rkM) const
{
    return Matrix4<Real>(
        m_afEntry[ 0]*rkM.m_afEntry[ 0] +
        m_afEntry[ 1]*rkM.m_afEntry[ 4] +
        m_afEntry[ 2]*rkM.m_afEntry[ 8] +
        m_afEntry[ 3]*rkM.m_afEntry[12],

        m_afEntry[ 0]*rkM.m_afEntry[ 1] +
        m_afEntry[ 1]*rkM.m_afEntry[ 5] +
        m_afEntry[ 2]*rkM.m_afEntry[ 9] +
        m_afEntry[ 3]*rkM.m_afEntry[13],

        m_afEntry[ 0]*rkM.m_afEntry[ 2] +
        m_afEntry[ 1]*rkM.m_afEntry[ 6] +
        m_afEntry[ 2]*rkM.m_afEntry[10] +
        m_afEntry[ 3]*rkM.m_afEntry[14],

        m_afEntry[ 0]*rkM.m_afEntry[ 3] +
        m_afEntry[ 1]*rkM.m_afEntry[ 7] +
        m_afEntry[ 2]*rkM.m_afEntry[11] +
        m_afEntry[ 3]*rkM.m_afEntry[15],

        m_afEntry[ 4]*rkM.m_afEntry[ 0] +
        m_afEntry[ 5]*rkM.m_afEntry[ 4] +
        m_afEntry[ 6]*rkM.m_afEntry[ 8] +
        m_afEntry[ 7]*rkM.m_afEntry[12],

        m_afEntry[ 4]*rkM.m_afEntry[ 1] +
        m_afEntry[ 5]*rkM.m_afEntry[ 5] +
        m_afEntry[ 6]*rkM.m_afEntry[ 9] +
        m_afEntry[ 7]*rkM.m_afEntry[13],

        m_afEntry[ 4]*rkM.m_afEntry[ 2] +
        m_afEntry[ 5]*rkM.m_afEntry[ 6] +
        m_afEntry[ 6]*rkM.m_afEntry[10] +
        m_afEntry[ 7]*rkM.m_afEntry[14],

        m_afEntry[ 4]*rkM.m_afEntry[ 3] +
        m_afEntry[ 5]*rkM.m_afEntry[ 7] +
        m_afEntry[ 6]*rkM.m_afEntry[11] +
        m_afEntry[ 7]*rkM.m_afEntry[15],

        m_afEntry[ 8]*rkM.m_afEntry[ 0] +
        m_afEntry[ 9]*rkM.m_afEntry[ 4] +
        m_afEntry[10]*rkM.m_afEntry[ 8] +
        m_afEntry[11]*rkM.m_afEntry[12],

        m_afEntry[ 8]*rkM.m_afEntry[ 1] +
        m_afEntry[ 9]*rkM.m_afEntry[ 5] +
        m_afEntry[10]*rkM.m_afEntry[ 9] +
        m_afEntry[11]*rkM.m_afEntry[13],

        m_afEntry[ 8]*rkM.m_afEntry[ 2] +
        m_afEntry[ 9]*rkM.m_afEntry[ 6] +
        m_afEntry[10]*rkM.m_afEntry[10] +
        m_afEntry[11]*rkM.m_afEntry[14],

        m_afEntry[ 8]*rkM.m_afEntry[ 3] +
        m_afEntry[ 9]*rkM.m_afEntry[ 7] +
        m_afEntry[10]*rkM.m_afEntry[11] +
        m_afEntry[11]*rkM.m_afEntry[15],

        m_afEntry[12]*rkM.m_afEntry[ 0] +
        m_afEntry[13]*rkM.m_afEntry[ 4] +
        m_afEntry[14]*rkM.m_afEntry[ 8] +
        m_afEntry[15]*rkM.m_afEntry[12],

        m_afEntry[12]*rkM.m_afEntry[ 1] +
        m_afEntry[13]*rkM.m_afEntry[ 5] +
        m_afEntry[14]*rkM.m_afEntry[ 9] +
        m_afEntry[15]*rkM.m_afEntry[13],

        m_afEntry[12]*rkM.m_afEntry[ 2] +
        m_afEntry[13]*rkM.m_afEntry[ 6] +
        m_afEntry[14]*rkM.m_afEntry[10] +
        m_afEntry[15]*rkM.m_afEntry[14],

        m_afEntry[12]*rkM.m_afEntry[ 3] +
        m_afEntry[13]*rkM.m_afEntry[ 7] +
        m_afEntry[14]*rkM.m_afEntry[11] +
        m_afEntry[15]*rkM.m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real> Matrix4<Real>::operator* (Real fScalar) const
{
    return Matrix4<Real>(
        fScalar*m_afEntry[ 0],
        fScalar*m_afEntry[ 1],
        fScalar*m_afEntry[ 2],
        fScalar*m_afEntry[ 3],
        fScalar*m_afEntry[ 4],
        fScalar*m_afEntry[ 5],
        fScalar*m_afEntry[ 6],
        fScalar*m_afEntry[ 7],
        fScalar*m_afEntry[ 8],
        fScalar*m_afEntry[ 9],
        fScalar*m_afEntry[10],
        fScalar*m_afEntry[11],
        fScalar*m_afEntry[12],
        fScalar*m_afEntry[13],
        fScalar*m_afEntry[14],
        fScalar*m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real> Matrix4<Real>::operator/ (Real fScalar) const
{
    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        return Matrix4<Real>(
            fInvScalar*m_afEntry[ 0],
            fInvScalar*m_afEntry[ 1],
            fInvScalar*m_afEntry[ 2],
            fInvScalar*m_afEntry[ 3],
            fInvScalar*m_afEntry[ 4],
            fInvScalar*m_afEntry[ 5],
            fInvScalar*m_afEntry[ 6],
            fInvScalar*m_afEntry[ 7],
            fInvScalar*m_afEntry[ 8],
            fInvScalar*m_afEntry[ 9],
            fInvScalar*m_afEntry[10],
            fInvScalar*m_afEntry[11],
            fInvScalar*m_afEntry[12],
            fInvScalar*m_afEntry[13],
            fInvScalar*m_afEntry[14],
            fInvScalar*m_afEntry[15]);
    }

    return Matrix4<Real>(
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL,
        Math<Real>::MAX_REAL);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real> Matrix4<Real>::operator- () const
{
    return Matrix4<Real>(
        -m_afEntry[ 0],
        -m_afEntry[ 1],
        -m_afEntry[ 2],
        -m_afEntry[ 3],
        -m_afEntry[ 4],
        -m_afEntry[ 5],
        -m_afEntry[ 6],
        -m_afEntry[ 7],
        -m_afEntry[ 8],
        -m_afEntry[ 9],
        -m_afEntry[10],
        -m_afEntry[11],
        -m_afEntry[12],
        -m_afEntry[13],
        -m_afEntry[14],
        -m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>& Matrix4<Real>::operator+= (const Matrix4& rkM)
{
    m_afEntry[ 0] += rkM.m_afEntry[ 0];
    m_afEntry[ 1] += rkM.m_afEntry[ 1];
    m_afEntry[ 2] += rkM.m_afEntry[ 2];
    m_afEntry[ 3] += rkM.m_afEntry[ 3];
    m_afEntry[ 4] += rkM.m_afEntry[ 4];
    m_afEntry[ 5] += rkM.m_afEntry[ 5];
    m_afEntry[ 6] += rkM.m_afEntry[ 6];
    m_afEntry[ 7] += rkM.m_afEntry[ 7];
    m_afEntry[ 8] += rkM.m_afEntry[ 8];
    m_afEntry[ 9] += rkM.m_afEntry[ 9];
    m_afEntry[10] += rkM.m_afEntry[10];
    m_afEntry[11] += rkM.m_afEntry[11];
    m_afEntry[12] += rkM.m_afEntry[12];
    m_afEntry[13] += rkM.m_afEntry[13];
    m_afEntry[14] += rkM.m_afEntry[14];
    m_afEntry[15] += rkM.m_afEntry[15];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>& Matrix4<Real>::operator-= (const Matrix4& rkM)
{
    m_afEntry[ 0] -= rkM.m_afEntry[ 0];
    m_afEntry[ 1] -= rkM.m_afEntry[ 1];
    m_afEntry[ 2] -= rkM.m_afEntry[ 2];
    m_afEntry[ 3] -= rkM.m_afEntry[ 3];
    m_afEntry[ 4] -= rkM.m_afEntry[ 4];
    m_afEntry[ 5] -= rkM.m_afEntry[ 5];
    m_afEntry[ 6] -= rkM.m_afEntry[ 6];
    m_afEntry[ 7] -= rkM.m_afEntry[ 7];
    m_afEntry[ 8] -= rkM.m_afEntry[ 8];
    m_afEntry[ 9] -= rkM.m_afEntry[ 9];
    m_afEntry[10] -= rkM.m_afEntry[10];
    m_afEntry[11] -= rkM.m_afEntry[11];
    m_afEntry[12] -= rkM.m_afEntry[12];
    m_afEntry[13] -= rkM.m_afEntry[13];
    m_afEntry[14] -= rkM.m_afEntry[14];
    m_afEntry[15] -= rkM.m_afEntry[15];
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>& Matrix4<Real>::operator*= (Real fScalar)
{
    m_afEntry[ 0] *= fScalar;
    m_afEntry[ 1] *= fScalar;
    m_afEntry[ 2] *= fScalar;
    m_afEntry[ 3] *= fScalar;
    m_afEntry[ 4] *= fScalar;
    m_afEntry[ 5] *= fScalar;
    m_afEntry[ 6] *= fScalar;
    m_afEntry[ 7] *= fScalar;
    m_afEntry[ 8] *= fScalar;
    m_afEntry[ 9] *= fScalar;
    m_afEntry[10] *= fScalar;
    m_afEntry[11] *= fScalar;
    m_afEntry[12] *= fScalar;
    m_afEntry[13] *= fScalar;
    m_afEntry[14] *= fScalar;
    m_afEntry[15] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Matrix4<Real>& Matrix4<Real>::operator/= (Real fScalar)
{
    int i;

    if (fScalar != (Real)0.0)
    {
        Real fInvScalar = ((Real)1.0)/fScalar;
        m_afEntry[ 0] *= fInvScalar;
        m_afEntry[ 1] *= fInvScalar;
        m_afEntry[ 2] *= fInvScalar;
        m_afEntry[ 3] *= fInvScalar;
        m_afEntry[ 4] *= fInvScalar;
        m_afEntry[ 5] *= fInvScalar;
        m_afEntry[ 6] *= fInvScalar;
        m_afEntry[ 7] *= fInvScalar;
        m_afEntry[ 8] *= fInvScalar;
        m_afEntry[ 9] *= fInvScalar;
        m_afEntry[10] *= fInvScalar;
        m_afEntry[11] *= fInvScalar;
        m_afEntry[12] *= fInvScalar;
        m_afEntry[13] *= fInvScalar;
        m_afEntry[14] *= fInvScalar;
        m_afEntry[15] *= fInvScalar;
    }
    else
    {
        m_afEntry[ 0] = Math<Real>::MAX_REAL;
        m_afEntry[ 1] = Math<Real>::MAX_REAL;
        m_afEntry[ 2] = Math<Real>::MAX_REAL;
        m_afEntry[ 3] = Math<Real>::MAX_REAL;
        m_afEntry[ 4] = Math<Real>::MAX_REAL;
        m_afEntry[ 5] = Math<Real>::MAX_REAL;
        m_afEntry[ 6] = Math<Real>::MAX_REAL;
        m_afEntry[ 7] = Math<Real>::MAX_REAL;
        m_afEntry[ 8] = Math<Real>::MAX_REAL;
        m_afEntry[ 9] = Math<Real>::MAX_REAL;
        m_afEntry[10] = Math<Real>::MAX_REAL;
        m_afEntry[11] = Math<Real>::MAX_REAL;
        m_afEntry[12] = Math<Real>::MAX_REAL;
        m_afEntry[13] = Math<Real>::MAX_REAL;
        m_afEntry[14] = Math<Real>::MAX_REAL;
        m_afEntry[15] = Math<Real>::MAX_REAL;
    }

    return *this;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> Matrix4<Real>::operator* (const Vector4<Real>& rkV) const
{
    return Vector4<Real>(
        m_afEntry[ 0]*rkV[0] +
        m_afEntry[ 1]*rkV[1] +
        m_afEntry[ 2]*rkV[2] +
        m_afEntry[ 3]*rkV[3],

        m_afEntry[ 4]*rkV[0] +
        m_afEntry[ 5]*rkV[1] +
        m_afEntry[ 6]*rkV[2] +
        m_afEntry[ 7]*rkV[3],

        m_afEntry[ 8]*rkV[0] +
        m_afEntry[ 9]*rkV[1] +
        m_afEntry[10]*rkV[2] +
        m_afEntry[11]*rkV[3],

        m_afEntry[12]*rkV[0] +
        m_afEntry[13]*rkV[1] +
        m_afEntry[14]*rkV[2] +
        m_afEntry[15]*rkV[3]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real> Matrix4<Real>::Transpose () const
{
    return Matrix4<Real>(
        m_afEntry[ 0],
        m_afEntry[ 4],
        m_afEntry[ 8],
        m_afEntry[12],
        m_afEntry[ 1],
        m_afEntry[ 5],
        m_afEntry[ 9],
        m_afEntry[13],
        m_afEntry[ 2],
        m_afEntry[ 6],
        m_afEntry[10],
        m_afEntry[14],
        m_afEntry[ 3],
        m_afEntry[ 7],
        m_afEntry[11],
        m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real> Matrix4<Real>::TransposeTimes (const Matrix4& rkM) const
{
    // P = A^T*B
    return Matrix4<Real>(
        m_afEntry[ 0]*rkM.m_afEntry[ 0] +
        m_afEntry[ 4]*rkM.m_afEntry[ 4] +
        m_afEntry[ 8]*rkM.m_afEntry[ 8] +
        m_afEntry[12]*rkM.m_afEntry[12],

        m_afEntry[ 0]*rkM.m_afEntry[ 1] +
        m_afEntry[ 4]*rkM.m_afEntry[ 5] +
        m_afEntry[ 8]*rkM.m_afEntry[ 9] +
        m_afEntry[12]*rkM.m_afEntry[13],

        m_afEntry[ 0]*rkM.m_afEntry[ 2] +
        m_afEntry[ 4]*rkM.m_afEntry[ 6] +
        m_afEntry[ 8]*rkM.m_afEntry[10] +
        m_afEntry[12]*rkM.m_afEntry[14],

        m_afEntry[ 0]*rkM.m_afEntry[ 3] +
        m_afEntry[ 4]*rkM.m_afEntry[ 7] +
        m_afEntry[ 8]*rkM.m_afEntry[11] +
        m_afEntry[12]*rkM.m_afEntry[15],

        m_afEntry[ 1]*rkM.m_afEntry[ 0] +
        m_afEntry[ 5]*rkM.m_afEntry[ 4] +
        m_afEntry[ 9]*rkM.m_afEntry[ 8] +
        m_afEntry[13]*rkM.m_afEntry[12],

        m_afEntry[ 1]*rkM.m_afEntry[ 1] +
        m_afEntry[ 5]*rkM.m_afEntry[ 5] +
        m_afEntry[ 9]*rkM.m_afEntry[ 9] +
        m_afEntry[13]*rkM.m_afEntry[13],

        m_afEntry[ 1]*rkM.m_afEntry[ 2] +
        m_afEntry[ 5]*rkM.m_afEntry[ 6] +
        m_afEntry[ 9]*rkM.m_afEntry[10] +
        m_afEntry[13]*rkM.m_afEntry[14],

        m_afEntry[ 1]*rkM.m_afEntry[ 3] +
        m_afEntry[ 5]*rkM.m_afEntry[ 7] +
        m_afEntry[ 9]*rkM.m_afEntry[11] +
        m_afEntry[13]*rkM.m_afEntry[15],

        m_afEntry[ 2]*rkM.m_afEntry[ 0] +
        m_afEntry[ 6]*rkM.m_afEntry[ 4] +
        m_afEntry[10]*rkM.m_afEntry[ 8] +
        m_afEntry[14]*rkM.m_afEntry[12],

        m_afEntry[ 2]*rkM.m_afEntry[ 1] +
        m_afEntry[ 6]*rkM.m_afEntry[ 5] +
        m_afEntry[10]*rkM.m_afEntry[ 9] +
        m_afEntry[14]*rkM.m_afEntry[13],

        m_afEntry[ 2]*rkM.m_afEntry[ 2] +
        m_afEntry[ 6]*rkM.m_afEntry[ 6] +
        m_afEntry[10]*rkM.m_afEntry[10] +
        m_afEntry[14]*rkM.m_afEntry[14],

        m_afEntry[ 2]*rkM.m_afEntry[ 3] +
        m_afEntry[ 6]*rkM.m_afEntry[ 7] +
        m_afEntry[10]*rkM.m_afEntry[11] +
        m_afEntry[14]*rkM.m_afEntry[15],

        m_afEntry[ 3]*rkM.m_afEntry[ 0] +
        m_afEntry[ 7]*rkM.m_afEntry[ 4] +
        m_afEntry[11]*rkM.m_afEntry[ 8] +
        m_afEntry[15]*rkM.m_afEntry[12],

        m_afEntry[ 3]*rkM.m_afEntry[ 1] +
        m_afEntry[ 7]*rkM.m_afEntry[ 5] +
        m_afEntry[11]*rkM.m_afEntry[ 9] +
        m_afEntry[15]*rkM.m_afEntry[13],

        m_afEntry[ 3]*rkM.m_afEntry[ 2] +
        m_afEntry[ 7]*rkM.m_afEntry[ 6] +
        m_afEntry[11]*rkM.m_afEntry[10] +
        m_afEntry[15]*rkM.m_afEntry[14],

        m_afEntry[ 3]*rkM.m_afEntry[ 3] +
        m_afEntry[ 7]*rkM.m_afEntry[ 7] +
        m_afEntry[11]*rkM.m_afEntry[11] +
        m_afEntry[15]*rkM.m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real> Matrix4<Real>::TimesTranspose (const Matrix4& rkM) const
{
    // P = A*B^T
    return Matrix4<Real>(
        m_afEntry[ 0]*rkM.m_afEntry[ 0] +
        m_afEntry[ 1]*rkM.m_afEntry[ 1] +
        m_afEntry[ 2]*rkM.m_afEntry[ 2] +
        m_afEntry[ 3]*rkM.m_afEntry[ 3],

        m_afEntry[ 0]*rkM.m_afEntry[ 4] +
        m_afEntry[ 1]*rkM.m_afEntry[ 5] +
        m_afEntry[ 2]*rkM.m_afEntry[ 6] +
        m_afEntry[ 3]*rkM.m_afEntry[ 7],

        m_afEntry[ 0]*rkM.m_afEntry[ 8] +
        m_afEntry[ 1]*rkM.m_afEntry[ 9] +
        m_afEntry[ 2]*rkM.m_afEntry[10] +
        m_afEntry[ 3]*rkM.m_afEntry[11],

        m_afEntry[ 0]*rkM.m_afEntry[12] +
        m_afEntry[ 1]*rkM.m_afEntry[13] +
        m_afEntry[ 2]*rkM.m_afEntry[14] +
        m_afEntry[ 3]*rkM.m_afEntry[15],

        m_afEntry[ 4]*rkM.m_afEntry[ 0] +
        m_afEntry[ 5]*rkM.m_afEntry[ 1] +
        m_afEntry[ 6]*rkM.m_afEntry[ 2] +
        m_afEntry[ 7]*rkM.m_afEntry[ 3],

        m_afEntry[ 4]*rkM.m_afEntry[ 4] +
        m_afEntry[ 5]*rkM.m_afEntry[ 5] +
        m_afEntry[ 6]*rkM.m_afEntry[ 6] +
        m_afEntry[ 7]*rkM.m_afEntry[ 7],

        m_afEntry[ 4]*rkM.m_afEntry[ 8] +
        m_afEntry[ 5]*rkM.m_afEntry[ 9] +
        m_afEntry[ 6]*rkM.m_afEntry[10] +
        m_afEntry[ 7]*rkM.m_afEntry[11],

        m_afEntry[ 4]*rkM.m_afEntry[12] +
        m_afEntry[ 5]*rkM.m_afEntry[13] +
        m_afEntry[ 6]*rkM.m_afEntry[14] +
        m_afEntry[ 7]*rkM.m_afEntry[15],

        m_afEntry[ 8]*rkM.m_afEntry[ 0] +
        m_afEntry[ 9]*rkM.m_afEntry[ 1] +
        m_afEntry[10]*rkM.m_afEntry[ 2] +
        m_afEntry[11]*rkM.m_afEntry[ 3],

        m_afEntry[ 8]*rkM.m_afEntry[ 4] +
        m_afEntry[ 9]*rkM.m_afEntry[ 5] +
        m_afEntry[10]*rkM.m_afEntry[ 6] +
        m_afEntry[11]*rkM.m_afEntry[ 7],

        m_afEntry[ 8]*rkM.m_afEntry[ 8] +
        m_afEntry[ 9]*rkM.m_afEntry[ 9] +
        m_afEntry[10]*rkM.m_afEntry[10] +
        m_afEntry[11]*rkM.m_afEntry[11],

        m_afEntry[ 8]*rkM.m_afEntry[12] +
        m_afEntry[ 9]*rkM.m_afEntry[13] +
        m_afEntry[10]*rkM.m_afEntry[14] +
        m_afEntry[11]*rkM.m_afEntry[15],

        m_afEntry[12]*rkM.m_afEntry[ 0] +
        m_afEntry[13]*rkM.m_afEntry[ 1] +
        m_afEntry[14]*rkM.m_afEntry[ 2] +
        m_afEntry[15]*rkM.m_afEntry[ 3],

        m_afEntry[12]*rkM.m_afEntry[ 4] +
        m_afEntry[13]*rkM.m_afEntry[ 5] +
        m_afEntry[14]*rkM.m_afEntry[ 6] +
        m_afEntry[15]*rkM.m_afEntry[ 7],

        m_afEntry[12]*rkM.m_afEntry[ 8] +
        m_afEntry[13]*rkM.m_afEntry[ 9] +
        m_afEntry[14]*rkM.m_afEntry[10] +
        m_afEntry[15]*rkM.m_afEntry[11],

        m_afEntry[12]*rkM.m_afEntry[12] +
        m_afEntry[13]*rkM.m_afEntry[13] +
        m_afEntry[14]*rkM.m_afEntry[14] +
        m_afEntry[15]*rkM.m_afEntry[15]);
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real> Matrix4<Real>::Inverse () const
{
    Real fA0 = m_afEntry[ 0]*m_afEntry[ 5] - m_afEntry[ 1]*m_afEntry[ 4];
    Real fA1 = m_afEntry[ 0]*m_afEntry[ 6] - m_afEntry[ 2]*m_afEntry[ 4];
    Real fA2 = m_afEntry[ 0]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 4];
    Real fA3 = m_afEntry[ 1]*m_afEntry[ 6] - m_afEntry[ 2]*m_afEntry[ 5];
    Real fA4 = m_afEntry[ 1]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 5];
    Real fA5 = m_afEntry[ 2]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 6];
    Real fB0 = m_afEntry[ 8]*m_afEntry[13] - m_afEntry[ 9]*m_afEntry[12];
    Real fB1 = m_afEntry[ 8]*m_afEntry[14] - m_afEntry[10]*m_afEntry[12];
    Real fB2 = m_afEntry[ 8]*m_afEntry[15] - m_afEntry[11]*m_afEntry[12];
    Real fB3 = m_afEntry[ 9]*m_afEntry[14] - m_afEntry[10]*m_afEntry[13];
    Real fB4 = m_afEntry[ 9]*m_afEntry[15] - m_afEntry[11]*m_afEntry[13];
    Real fB5 = m_afEntry[10]*m_afEntry[15] - m_afEntry[11]*m_afEntry[14];

	Real fDet = fA0*fB5-fA1*fB4+fA2*fB3+fA3*fB2-fA4*fB1+fA5*fB0;
	Matrix4 kInv;
    if (Math<Real>::FAbs(fDet) <= Math<Real>::INVERSE_TOLERANCE)
    {
        return kInv;
    }

    kInv.m_afEntry[ 0] =
        + m_afEntry[ 5]*fB5 - m_afEntry[ 6]*fB4 + m_afEntry[ 7]*fB3;
    kInv.m_afEntry[ 4] =
        - m_afEntry[ 4]*fB5 + m_afEntry[ 6]*fB2 - m_afEntry[ 7]*fB1;
    kInv.m_afEntry[ 8] =
        + m_afEntry[ 4]*fB4 - m_afEntry[ 5]*fB2 + m_afEntry[ 7]*fB0;
    kInv.m_afEntry[12] =
        - m_afEntry[ 4]*fB3 + m_afEntry[ 5]*fB1 - m_afEntry[ 6]*fB0;
    kInv.m_afEntry[ 1] =
        - m_afEntry[ 1]*fB5 + m_afEntry[ 2]*fB4 - m_afEntry[ 3]*fB3;
    kInv.m_afEntry[ 5] =
        + m_afEntry[ 0]*fB5 - m_afEntry[ 2]*fB2 + m_afEntry[ 3]*fB1;
    kInv.m_afEntry[ 9] =
        - m_afEntry[ 0]*fB4 + m_afEntry[ 1]*fB2 - m_afEntry[ 3]*fB0;
    kInv.m_afEntry[13] =
        + m_afEntry[ 0]*fB3 - m_afEntry[ 1]*fB1 + m_afEntry[ 2]*fB0;
    kInv.m_afEntry[ 2] =
        + m_afEntry[13]*fA5 - m_afEntry[14]*fA4 + m_afEntry[15]*fA3;
    kInv.m_afEntry[ 6] =
        - m_afEntry[12]*fA5 + m_afEntry[14]*fA2 - m_afEntry[15]*fA1;
    kInv.m_afEntry[10] =
        + m_afEntry[12]*fA4 - m_afEntry[13]*fA2 + m_afEntry[15]*fA0;
    kInv.m_afEntry[14] =
        - m_afEntry[12]*fA3 + m_afEntry[13]*fA1 - m_afEntry[14]*fA0;
    kInv.m_afEntry[ 3] =
        - m_afEntry[ 9]*fA5 + m_afEntry[10]*fA4 - m_afEntry[11]*fA3;
    kInv.m_afEntry[ 7] =
        + m_afEntry[ 8]*fA5 - m_afEntry[10]*fA2 + m_afEntry[11]*fA1;
    kInv.m_afEntry[11] =
        - m_afEntry[ 8]*fA4 + m_afEntry[ 9]*fA2 - m_afEntry[11]*fA0;
    kInv.m_afEntry[15] =
        + m_afEntry[ 8]*fA3 - m_afEntry[ 9]*fA1 + m_afEntry[10]*fA0;

    Real fInvDet = ((Real)1.0)/fDet;
    kInv.m_afEntry[ 0] *= fInvDet;
    kInv.m_afEntry[ 1] *= fInvDet;
    kInv.m_afEntry[ 2] *= fInvDet;
    kInv.m_afEntry[ 3] *= fInvDet;
    kInv.m_afEntry[ 4] *= fInvDet;
    kInv.m_afEntry[ 5] *= fInvDet;
    kInv.m_afEntry[ 6] *= fInvDet;
    kInv.m_afEntry[ 7] *= fInvDet;
    kInv.m_afEntry[ 8] *= fInvDet;
    kInv.m_afEntry[ 9] *= fInvDet;
    kInv.m_afEntry[10] *= fInvDet;
    kInv.m_afEntry[11] *= fInvDet;
    kInv.m_afEntry[12] *= fInvDet;
    kInv.m_afEntry[13] *= fInvDet;
    kInv.m_afEntry[14] *= fInvDet;
    kInv.m_afEntry[15] *= fInvDet;

    return kInv;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real> Matrix4<Real>::Adjoint () const
{
    Real fA0 = m_afEntry[ 0]*m_afEntry[ 5] - m_afEntry[ 1]*m_afEntry[ 4];
    Real fA1 = m_afEntry[ 0]*m_afEntry[ 6] - m_afEntry[ 2]*m_afEntry[ 4];
    Real fA2 = m_afEntry[ 0]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 4];
    Real fA3 = m_afEntry[ 1]*m_afEntry[ 6] - m_afEntry[ 2]*m_afEntry[ 5];
    Real fA4 = m_afEntry[ 1]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 5];
    Real fA5 = m_afEntry[ 2]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 6];
    Real fB0 = m_afEntry[ 8]*m_afEntry[13] - m_afEntry[ 9]*m_afEntry[12];
    Real fB1 = m_afEntry[ 8]*m_afEntry[14] - m_afEntry[10]*m_afEntry[12];
    Real fB2 = m_afEntry[ 8]*m_afEntry[15] - m_afEntry[11]*m_afEntry[12];
    Real fB3 = m_afEntry[ 9]*m_afEntry[14] - m_afEntry[10]*m_afEntry[13];
    Real fB4 = m_afEntry[ 9]*m_afEntry[15] - m_afEntry[11]*m_afEntry[13];
    Real fB5 = m_afEntry[10]*m_afEntry[15] - m_afEntry[11]*m_afEntry[14];

    return Matrix4<Real>(
        + m_afEntry[ 5]*fB5 - m_afEntry[ 6]*fB4 + m_afEntry[ 7]*fB3,
        - m_afEntry[ 1]*fB5 + m_afEntry[ 2]*fB4 - m_afEntry[ 3]*fB3,
        + m_afEntry[13]*fA5 - m_afEntry[14]*fA4 + m_afEntry[15]*fA3,
        - m_afEntry[ 9]*fA5 + m_afEntry[10]*fA4 - m_afEntry[11]*fA3,
        - m_afEntry[ 4]*fB5 + m_afEntry[ 6]*fB2 - m_afEntry[ 7]*fB1,
        + m_afEntry[ 0]*fB5 - m_afEntry[ 2]*fB2 + m_afEntry[ 3]*fB1,
        - m_afEntry[12]*fA5 + m_afEntry[14]*fA2 - m_afEntry[15]*fA1,
        + m_afEntry[ 8]*fA5 - m_afEntry[10]*fA2 + m_afEntry[11]*fA1,
        + m_afEntry[ 4]*fB4 - m_afEntry[ 5]*fB2 + m_afEntry[ 7]*fB0,
        - m_afEntry[ 0]*fB4 + m_afEntry[ 1]*fB2 - m_afEntry[ 3]*fB0,
        + m_afEntry[12]*fA4 - m_afEntry[13]*fA2 + m_afEntry[15]*fA0,
        - m_afEntry[ 8]*fA4 + m_afEntry[ 9]*fA2 - m_afEntry[11]*fA0,
        - m_afEntry[ 4]*fB3 + m_afEntry[ 5]*fB1 - m_afEntry[ 6]*fB0,
        + m_afEntry[ 0]*fB3 - m_afEntry[ 1]*fB1 + m_afEntry[ 2]*fB0,
        - m_afEntry[12]*fA3 + m_afEntry[13]*fA1 - m_afEntry[14]*fA0,
        + m_afEntry[ 8]*fA3 - m_afEntry[ 9]*fA1 + m_afEntry[10]*fA0);
}
//----------------------------------------------------------------------------
template <class Real>
Real Matrix4<Real>::Determinant () const
{
    Real fA0 = m_afEntry[ 0]*m_afEntry[ 5] - m_afEntry[ 1]*m_afEntry[ 4];
    Real fA1 = m_afEntry[ 0]*m_afEntry[ 6] - m_afEntry[ 2]*m_afEntry[ 4];
    Real fA2 = m_afEntry[ 0]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 4];
    Real fA3 = m_afEntry[ 1]*m_afEntry[ 6] - m_afEntry[ 2]*m_afEntry[ 5];
    Real fA4 = m_afEntry[ 1]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 5];
    Real fA5 = m_afEntry[ 2]*m_afEntry[ 7] - m_afEntry[ 3]*m_afEntry[ 6];
    Real fB0 = m_afEntry[ 8]*m_afEntry[13] - m_afEntry[ 9]*m_afEntry[12];
    Real fB1 = m_afEntry[ 8]*m_afEntry[14] - m_afEntry[10]*m_afEntry[12];
    Real fB2 = m_afEntry[ 8]*m_afEntry[15] - m_afEntry[11]*m_afEntry[12];
    Real fB3 = m_afEntry[ 9]*m_afEntry[14] - m_afEntry[10]*m_afEntry[13];
    Real fB4 = m_afEntry[ 9]*m_afEntry[15] - m_afEntry[11]*m_afEntry[13];
    Real fB5 = m_afEntry[10]*m_afEntry[15] - m_afEntry[11]*m_afEntry[14];
    Real fDet = fA0*fB5-fA1*fB4+fA2*fB3+fA3*fB2-fA4*fB1+fA5*fB0;
    return fDet;
}
//----------------------------------------------------------------------------
template <class Real>
Real Matrix4<Real>::QForm (const Vector4<Real>& rkU,
    const Vector4<Real>& rkV) const
{
    return rkU.Dot((*this)*rkV);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::MakeObliqueProjection (const Vector3<Real>& rkNormal,
    const Vector3<Real>& rkPoint, const Vector3<Real>& rkDirection)
{
    // The projection plane is Dot(N,X-P) = 0 where N is a 3-by-1 unit-length
    // normal vector and P is a 3-by-1 point on the plane.  The projection
    // is oblique to the plane, in the direction of the 3-by-1 vector D.
    // Necessarily Dot(N,D) is not zero for this projection to make sense.
    // Given a 3-by-1 point U, compute the intersection of the line U+t*D
    // with the plane to obtain t = -Dot(N,U-P)/Dot(N,D).  Then
    //
    //   projection(U) = P + [I - D*N^T/Dot(N,D)]*(U-P)
    //
    // A 4-by-4 homogeneous transformation representing the projection is
    //
    //       +-                               -+
    //   M = | D*N^T - Dot(N,D)*I   -Dot(N,P)D |
    //       |          0^T          -Dot(N,D) |
    //       +-                               -+
    //
    // where M applies to [U^T 1]^T by M*[U^T 1]^T.  The matrix is chosen so
    // that M[3][3] > 0 whenever Dot(N,D) < 0 (projection is onto the
    // "positive side" of the plane).

    Real fNdD = rkNormal.Dot(rkDirection);
    Real fNdP = rkNormal.Dot(rkPoint);
    m_afEntry[ 0] = rkDirection[0]*rkNormal[0] - fNdD;
    m_afEntry[ 1] = rkDirection[0]*rkNormal[1];
    m_afEntry[ 2] = rkDirection[0]*rkNormal[2];
    m_afEntry[ 3] = -fNdP*rkDirection[0];
    m_afEntry[ 4] = rkDirection[1]*rkNormal[0];
    m_afEntry[ 5] = rkDirection[1]*rkNormal[1] - fNdD;
    m_afEntry[ 6] = rkDirection[1]*rkNormal[2];
    m_afEntry[ 7] = -fNdP*rkDirection[1];
    m_afEntry[ 8] = rkDirection[2]*rkNormal[0];
    m_afEntry[ 9] = rkDirection[2]*rkNormal[1];
    m_afEntry[10] = rkDirection[2]*rkNormal[2] - fNdD;
    m_afEntry[11] = -fNdP*rkDirection[2];
    m_afEntry[12] = 0.0f;
    m_afEntry[13] = 0.0f;
    m_afEntry[14] = 0.0f;
    m_afEntry[15] = -fNdD;
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::MakePerspectiveProjection (const Vector3<Real>& rkNormal,
    const Vector3<Real>& rkPoint, const Vector3<Real>& rkEye)
{
    //     +-                                                 -+
    // M = | Dot(N,E-P)*I - E*N^T    -(Dot(N,E-P)*I - E*N^T)*E |
    //     |        -N^t                      Dot(N,E)         |
    //     +-                                                 -+
    //
    // where E is the eye point, P is a point on the plane, and N is a
    // unit-length plane normal.

    Real fNdEmP = rkNormal.Dot(rkEye-rkPoint);

    m_afEntry[ 0] = fNdEmP - rkEye[0]*rkNormal[0];
    m_afEntry[ 1] = -rkEye[0]*rkNormal[1];
    m_afEntry[ 2] = -rkEye[0]*rkNormal[2];
    m_afEntry[ 3] = -(m_afEntry[0]*rkEye[0] + m_afEntry[1]*rkEye[1] +
        m_afEntry[2]*rkEye[2]);
    m_afEntry[ 4] = -rkEye[1]*rkNormal[0];
    m_afEntry[ 5] = fNdEmP - rkEye[1]*rkNormal[1];
    m_afEntry[ 6] = -rkEye[1]*rkNormal[2];
    m_afEntry[ 7] = -(m_afEntry[4]*rkEye[0] + m_afEntry[5]*rkEye[1] +
        m_afEntry[6]*rkEye[2]);
    m_afEntry[ 8] = -rkEye[2]*rkNormal[0];
    m_afEntry[ 9] = -rkEye[2]*rkNormal[1];
    m_afEntry[10] = fNdEmP- rkEye[2]*rkNormal[2];
    m_afEntry[11] = -(m_afEntry[8]*rkEye[0] + m_afEntry[9]*rkEye[1] +
        m_afEntry[10]*rkEye[2]);
    m_afEntry[12] = -rkNormal[0];
    m_afEntry[13] = -rkNormal[1];
    m_afEntry[14] = -rkNormal[2];
    m_afEntry[15] = rkNormal.Dot(rkEye);
}
//----------------------------------------------------------------------------
template <class Real>
void Matrix4<Real>::MakeReflection (const Vector3<Real>& rkNormal,
    const Vector3<Real>& rkPoint)
{
    //     +-                         -+
    // M = | I-2*N*N^T    2*Dot(N,P)*N |
    //     |     0^T            1      |
    //     +-                         -+
    //
    // where P is a point on the plane and N is a unit-length plane normal.

    Real fTwoNdP = ((Real)2.0)*(rkNormal.Dot(rkPoint));

    m_afEntry[ 0] = (Real)1.0 - ((Real)2.0)*rkNormal[0]*rkNormal[0];
    m_afEntry[ 1] = -((Real)2.0)*rkNormal[0]*rkNormal[1];
    m_afEntry[ 2] = -((Real)2.0)*rkNormal[0]*rkNormal[2];
    m_afEntry[ 3] = fTwoNdP*rkNormal[0];
    m_afEntry[ 4] = -((Real)2.0)*rkNormal[1]*rkNormal[0];
    m_afEntry[ 5] = (Real)1.0 - ((Real)2.0)*rkNormal[1]*rkNormal[1];
    m_afEntry[ 6] = -((Real)2.0)*rkNormal[1]*rkNormal[2];
    m_afEntry[ 7] = fTwoNdP*rkNormal[1];
    m_afEntry[ 8] = -((Real)2.0)*rkNormal[2]*rkNormal[0];
    m_afEntry[ 9] = -((Real)2.0)*rkNormal[2]*rkNormal[1];
    m_afEntry[10] = (Real)1.0 - ((Real)2.0)*rkNormal[2]*rkNormal[2];
    m_afEntry[11] = fTwoNdP*rkNormal[2];
    m_afEntry[12] = (Real)0.0;
    m_afEntry[13] = (Real)0.0;
    m_afEntry[14] = (Real)0.0;
    m_afEntry[15] = (Real)1.0;
}
//----------------------------------------------------------------------------
template <class Real>
Matrix4<Real>& Matrix4<Real>::MakeTensorProduct (const Vector4<Real>& rkU,
	const Vector4<Real>& rkV)
{
	m_afEntry[0] = rkU[0]*rkV[0];
	m_afEntry[1] = rkU[0]*rkV[1];
	m_afEntry[2] = rkU[0]*rkV[2];
	m_afEntry[3] = rkU[0]*rkV[3];
	m_afEntry[4] = rkU[1]*rkV[0];
	m_afEntry[5] = rkU[1]*rkV[1];
	m_afEntry[6] = rkU[1]*rkV[2];
	m_afEntry[7] = rkU[1]*rkV[3];
	m_afEntry[8] = rkU[2]*rkV[0];
	m_afEntry[9] = rkU[2]*rkV[1];
	m_afEntry[10] = rkU[2]*rkV[2];
	m_afEntry[11] = rkU[2]*rkV[3];
	m_afEntry[12] = rkU[3]*rkV[0];
	m_afEntry[13] = rkU[3]*rkV[1];
	m_afEntry[14] = rkU[3]*rkV[2];
	m_afEntry[15] = rkU[3]*rkV[3];
	return *this;
}
//----------------------------------------------------------------------------

template <class Real>
inline Matrix4<Real> operator* (Real fScalar, const Matrix4<Real>& rkM)
{
    return rkM*fScalar;
}
//----------------------------------------------------------------------------
template <class Real>
inline Vector4<Real> operator* (const Vector4<Real>& rkV,
    const Matrix4<Real>& rkM)
{
    return Vector4<Real>(
        rkV[0]*rkM[0][0]+rkV[1]*rkM[1][0]+rkV[2]*rkM[2][0]+rkV[3]*rkM[3][0],
        rkV[0]*rkM[0][1]+rkV[1]*rkM[1][1]+rkV[2]*rkM[2][1]+rkV[3]*rkM[3][1],
        rkV[0]*rkM[0][2]+rkV[1]*rkM[1][2]+rkV[2]*rkM[2][2]+rkV[3]*rkM[3][2],
        rkV[0]*rkM[0][3]+rkV[1]*rkM[1][3]+rkV[2]*rkM[2][3]+rkV[3]*rkM[3][3]);
}
//----------------------------------------------------------------------------


typedef Matrix4<float> Matrix4f;
typedef Matrix4<double> Matrix4d;

}






namespace Wm4
{

	template <class Real>
	class Matrix5
	{
	public:
		// If bZero is true, create the zero matrix.  Otherwise, create the
		// identity matrix.
		Matrix5 (bool bZero = true);

		// copy constructor
		Matrix5 (const Matrix5& rkM);

		// input Mrc is in row r, column c.
		Matrix5 (Real fM00, Real fM01, Real fM02, Real fM03, Real fM04,
			Real fM10, Real fM11, Real fM12, Real fM13, Real fM14,
			Real fM20, Real fM21, Real fM22, Real fM23, Real fM24,
			Real fM30, Real fM31, Real fM32, Real fM33, Real fM34,
			Real fM40, Real fM41, Real fM42, Real fM43, Real fM44);

		void MakeZero ();
		void MakeIdentity ();
		
		// member access
		inline operator const Real* () const;
		inline operator Real* ();
		inline const Real* operator[] (int iRow) const;
		inline Real* operator[] (int iRow);
		inline Real operator() (int iRow, int iCol) const;
		inline Real& operator() (int iRow, int iCol);

		void SetRow (int iRow, const Vector5<Real>& rkV);
		Vector5<Real> GetRow (int iRow) const;
		void SetColumn (int iCol, const Vector5<Real>& rkV);
		Vector5<Real> GetColumn (int iCol) const;

		// assignment
		inline Matrix5& operator= (const Matrix5& rkM);

		// arithmetic operations
		inline Matrix5 operator+ (const Matrix5& rkM) const;
		inline Matrix5 operator* (const Matrix5& rkM) const;
		inline Matrix5 operator* (Real fScalar) const;


		bool operator== (const Matrix5& rkM) const;

		// matrix times vector
		inline Vector5<Real> operator* (const Vector5<Real>& rkV) const;  // M * v

		// other operations
		Matrix5 Inverse () const;
		Matrix5 Inverse2() const;
		Real QForm (const Vector5<Real>& rkU,
			const Vector5<Real>& rkV) const;  // u^T*M*v

	private:
		// support for comparisons
		int CompareArrays (const Matrix5& rkM) const;
		
	Real m_afEntry[25];
	};

	// v^T * M
	template <class Real>
	inline Vector5<Real> operator* (const Vector5<Real>& rkV,
		const Matrix5<Real>& rkM);

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
	Matrix5<Real>::Matrix5 (bool bZero)
	{
		if (bZero)
		{
			MakeZero();
		}
		else
		{
			MakeIdentity();
		}
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Matrix5<Real>::Matrix5 (const Matrix5& rkM)
	{
		m_afEntry[ 0] = rkM.m_afEntry[ 0];
		m_afEntry[ 1] = rkM.m_afEntry[ 1];
		m_afEntry[ 2] = rkM.m_afEntry[ 2];
		m_afEntry[ 3] = rkM.m_afEntry[ 3];
		m_afEntry[ 4] = rkM.m_afEntry[ 4];
		m_afEntry[ 5] = rkM.m_afEntry[ 5];
		m_afEntry[ 6] = rkM.m_afEntry[ 6];
		m_afEntry[ 7] = rkM.m_afEntry[ 7];
		m_afEntry[ 8] = rkM.m_afEntry[ 8];
		m_afEntry[ 9] = rkM.m_afEntry[ 9];
		m_afEntry[10] = rkM.m_afEntry[10];
		m_afEntry[11] = rkM.m_afEntry[11];
		m_afEntry[12] = rkM.m_afEntry[12];
		m_afEntry[13] = rkM.m_afEntry[13];
		m_afEntry[14] = rkM.m_afEntry[14];
		m_afEntry[15] = rkM.m_afEntry[15];		
		m_afEntry[16] = rkM.m_afEntry[16];
		m_afEntry[17] = rkM.m_afEntry[17];
		m_afEntry[18] = rkM.m_afEntry[18];
		m_afEntry[19] = rkM.m_afEntry[19];
		m_afEntry[20] = rkM.m_afEntry[20];
		m_afEntry[21] = rkM.m_afEntry[21];
		m_afEntry[22] = rkM.m_afEntry[22];
		m_afEntry[23] = rkM.m_afEntry[23];
		m_afEntry[24] = rkM.m_afEntry[24];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Matrix5<Real>::Matrix5 (Real fM00, Real fM01, Real fM02, Real fM03, Real fM04,
		Real fM10, Real fM11, Real fM12, Real fM13, Real fM14,
		Real fM20, Real fM21, Real fM22, Real fM23, Real fM24,
		Real fM30, Real fM31, Real fM32, Real fM33, Real fM34,
		Real fM40, Real fM41, Real fM42, Real fM43, Real fM44)
	{
		m_afEntry[ 0] = fM00;
		m_afEntry[ 1] = fM01;
		m_afEntry[ 2] = fM02;
		m_afEntry[ 3] = fM03;
		m_afEntry[ 4] = fM04;
		m_afEntry[ 5] = fM10;
		m_afEntry[ 6] = fM11;
		m_afEntry[ 7] = fM12;
		m_afEntry[ 8] = fM13;
		m_afEntry[ 9] = fM14;
		m_afEntry[10] = fM20;
		m_afEntry[11] = fM21;
		m_afEntry[12] = fM22;
		m_afEntry[13] = fM23;
		m_afEntry[14] = fM24;
		m_afEntry[15] = fM30;
		m_afEntry[16] = fM31;
		m_afEntry[17] = fM32;
		m_afEntry[18] = fM33;
		m_afEntry[19] = fM34;
		m_afEntry[20] = fM40;
		m_afEntry[21] = fM41;
		m_afEntry[22] = fM42;
		m_afEntry[23] = fM43;
		m_afEntry[24] = fM44;
	}
	
	//----------------------------------------------------------------------------
	template <class Real>
	void Matrix5<Real>::MakeZero ()
	{
		m_afEntry[ 0] = (Real)0.0;
		m_afEntry[ 1] = (Real)0.0;
		m_afEntry[ 2] = (Real)0.0;
		m_afEntry[ 3] = (Real)0.0;
		m_afEntry[ 4] = (Real)0.0;
		m_afEntry[ 5] = (Real)0.0;
		m_afEntry[ 6] = (Real)0.0;
		m_afEntry[ 7] = (Real)0.0;
		m_afEntry[ 8] = (Real)0.0;
		m_afEntry[ 9] = (Real)0.0;
		m_afEntry[10] = (Real)0.0;
		m_afEntry[11] = (Real)0.0;
		m_afEntry[12] = (Real)0.0;
		m_afEntry[13] = (Real)0.0;
		m_afEntry[14] = (Real)0.0;
		m_afEntry[15] = (Real)0.0;
		m_afEntry[16] = (Real)0.0;
		m_afEntry[17] = (Real)0.0;
		m_afEntry[18] = (Real)0.0;
		m_afEntry[19] = (Real)0.0;
		m_afEntry[20] = (Real)0.0;
		m_afEntry[21] = (Real)0.0;
		m_afEntry[22] = (Real)0.0;
		m_afEntry[23] = (Real)0.0;
		m_afEntry[24] = (Real)0.0;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	void Matrix5<Real>::MakeIdentity ()
	{
		m_afEntry[ 0] = (Real)1.0;
		m_afEntry[ 1] = (Real)0.0;
		m_afEntry[ 2] = (Real)0.0;
		m_afEntry[ 3] = (Real)0.0;
		m_afEntry[ 4] = (Real)0.0;
		m_afEntry[ 5] = (Real)0.0;
		m_afEntry[ 6] = (Real)1.0;
		m_afEntry[ 7] = (Real)0.0;
		m_afEntry[ 8] = (Real)0.0;
		m_afEntry[ 9] = (Real)0.0;
		m_afEntry[10] = (Real)0.0;
		m_afEntry[11] = (Real)0.0;
		m_afEntry[12] = (Real)1.0;
		m_afEntry[13] = (Real)0.0;
		m_afEntry[14] = (Real)0.0;
		m_afEntry[15] = (Real)0.0;
		m_afEntry[16] = (Real)0.0;
		m_afEntry[17] = (Real)0.0;
		m_afEntry[18] = (Real)1.0;
		m_afEntry[19] = (Real)0.0;
		m_afEntry[20] = (Real)0.0;
		m_afEntry[21] = (Real)0.0;
		m_afEntry[22] = (Real)0.0;
		m_afEntry[23] = (Real)0.0;
		m_afEntry[24] = (Real)1.0;
	}
	
	//----------------------------------------------------------------------------
	template <class Real>
	inline Matrix5<Real>::operator const Real* () const
	{
		return m_afEntry;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Matrix5<Real>::operator Real* ()
	{
		return m_afEntry;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline const Real* Matrix5<Real>::operator[] (int iRow) const
	{
		return &m_afEntry[5*iRow];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real* Matrix5<Real>::operator[] (int iRow)
	{
		return &m_afEntry[5*iRow];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real Matrix5<Real>::operator() (int iRow, int iCol) const
	{
		return m_afEntry[iCol+5*iRow];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Real& Matrix5<Real>::operator() (int iRow, int iCol)
	{
		return m_afEntry[iCol+5*iRow];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	void Matrix5<Real>::SetRow (int iRow, const Vector5<Real>& rkV)
	{
		int i0 = 5*iRow, i1 = i0+1, i2 = i1+1, i3 = i2+1, i4 = i3+1;
		m_afEntry[i0] = rkV[0];
		m_afEntry[i1] = rkV[1];
		m_afEntry[i2] = rkV[2];
		m_afEntry[i3] = rkV[3];
		m_afEntry[i4] = rkV[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Vector5<Real> Matrix5<Real>::GetRow (int iRow) const
	{
		int i0 = 5*iRow, i1 = i0+1, i2 = i1+1, i3 = i2+1, i4 = i3+1;
		return Vector5<Real>(m_afEntry[i0],m_afEntry[i1],m_afEntry[i2],
			m_afEntry[i3],m_afEntry[i4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	void Matrix5<Real>::SetColumn (int iCol, const Vector5<Real>& rkV)
	{
		m_afEntry[iCol] = rkV[0];
		m_afEntry[iCol+5] = rkV[1];
		m_afEntry[iCol+10] = rkV[2];
		m_afEntry[iCol+15] = rkV[3];
		m_afEntry[iCol+20] = rkV[4];
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Vector5<Real> Matrix5<Real>::GetColumn (int iCol) const
	{
		return Vector5<Real>(m_afEntry[iCol],m_afEntry[iCol+5],m_afEntry[iCol+10],
			m_afEntry[iCol+15], m_afEntry[iCol+20]);
	}
	//----------------------------------------------------------------------------
	
	template <class Real>
	inline Matrix5<Real>& Matrix5<Real>::operator= (const Matrix5& rkM)
	{
		m_afEntry[ 0] = rkM.m_afEntry[ 0];
		m_afEntry[ 1] = rkM.m_afEntry[ 1];
		m_afEntry[ 2] = rkM.m_afEntry[ 2];
		m_afEntry[ 3] = rkM.m_afEntry[ 3];
		m_afEntry[ 4] = rkM.m_afEntry[ 4];
		m_afEntry[ 5] = rkM.m_afEntry[ 5];
		m_afEntry[ 6] = rkM.m_afEntry[ 6];
		m_afEntry[ 7] = rkM.m_afEntry[ 7];
		m_afEntry[ 8] = rkM.m_afEntry[ 8];
		m_afEntry[ 9] = rkM.m_afEntry[ 9];
		m_afEntry[10] = rkM.m_afEntry[10];
		m_afEntry[11] = rkM.m_afEntry[11];
		m_afEntry[12] = rkM.m_afEntry[12];
		m_afEntry[13] = rkM.m_afEntry[13];
		m_afEntry[14] = rkM.m_afEntry[14];
		m_afEntry[15] = rkM.m_afEntry[15];	
		m_afEntry[16] = rkM.m_afEntry[16];
		m_afEntry[17] = rkM.m_afEntry[17];
		m_afEntry[18] = rkM.m_afEntry[18];
		m_afEntry[19] = rkM.m_afEntry[19];
		m_afEntry[20] = rkM.m_afEntry[20];
		m_afEntry[21] = rkM.m_afEntry[21];
		m_afEntry[22] = rkM.m_afEntry[22];
		m_afEntry[23] = rkM.m_afEntry[23];
		m_afEntry[24] = rkM.m_afEntry[24];
		return *this;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Matrix5<Real> Matrix5<Real>::operator+ (const Matrix5& rkM) const
	{
		return Matrix5<Real>(
			m_afEntry[ 0] + rkM.m_afEntry[ 0],
			m_afEntry[ 1] + rkM.m_afEntry[ 1],
			m_afEntry[ 2] + rkM.m_afEntry[ 2],
			m_afEntry[ 3] + rkM.m_afEntry[ 3],
			m_afEntry[ 4] + rkM.m_afEntry[ 4],
			m_afEntry[ 5] + rkM.m_afEntry[ 5],
			m_afEntry[ 6] + rkM.m_afEntry[ 6],
			m_afEntry[ 7] + rkM.m_afEntry[ 7],
			m_afEntry[ 8] + rkM.m_afEntry[ 8],
			m_afEntry[ 9] + rkM.m_afEntry[ 9],
			m_afEntry[10] + rkM.m_afEntry[10],
			m_afEntry[11] + rkM.m_afEntry[11],
			m_afEntry[12] + rkM.m_afEntry[12],
			m_afEntry[13] + rkM.m_afEntry[13],
			m_afEntry[14] + rkM.m_afEntry[14],
			m_afEntry[15] + rkM.m_afEntry[15],
			m_afEntry[16] + rkM.m_afEntry[16],
			m_afEntry[17] + rkM.m_afEntry[17],
			m_afEntry[18] + rkM.m_afEntry[18],
			m_afEntry[19] + rkM.m_afEntry[19],
			m_afEntry[20] + rkM.m_afEntry[20],
			m_afEntry[21] + rkM.m_afEntry[21],
			m_afEntry[22] + rkM.m_afEntry[22],
			m_afEntry[23] + rkM.m_afEntry[23],
			m_afEntry[24] + rkM.m_afEntry[24]
		);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Matrix5<Real> Matrix5<Real>::operator* (const Matrix5& rkM) const
	{
		return Matrix5<Real>(
			m_afEntry[ 0]*rkM.m_afEntry[ 0] +
			m_afEntry[ 1]*rkM.m_afEntry[ 5] +
			m_afEntry[ 2]*rkM.m_afEntry[10] +
			m_afEntry[ 3]*rkM.m_afEntry[15] +
			m_afEntry[ 4]*rkM.m_afEntry[20],

			m_afEntry[ 0]*rkM.m_afEntry[ 1] +
			m_afEntry[ 1]*rkM.m_afEntry[ 6] +
			m_afEntry[ 2]*rkM.m_afEntry[11] +
			m_afEntry[ 3]*rkM.m_afEntry[16] +
			m_afEntry[ 4]*rkM.m_afEntry[21],

			m_afEntry[ 0]*rkM.m_afEntry[ 2] +
			m_afEntry[ 1]*rkM.m_afEntry[ 7] +
			m_afEntry[ 2]*rkM.m_afEntry[12] +
			m_afEntry[ 3]*rkM.m_afEntry[17] +
			m_afEntry[ 4]*rkM.m_afEntry[22],

			m_afEntry[ 0]*rkM.m_afEntry[ 3] +
			m_afEntry[ 1]*rkM.m_afEntry[ 8] +
			m_afEntry[ 2]*rkM.m_afEntry[13] +
			m_afEntry[ 3]*rkM.m_afEntry[18] +
			m_afEntry[ 4]*rkM.m_afEntry[23],

			m_afEntry[ 0]*rkM.m_afEntry[ 4] +
			m_afEntry[ 1]*rkM.m_afEntry[ 9] +
			m_afEntry[ 2]*rkM.m_afEntry[14] +
			m_afEntry[ 3]*rkM.m_afEntry[19] +
			m_afEntry[ 4]*rkM.m_afEntry[24],

			m_afEntry[ 5]*rkM.m_afEntry[ 0] +
			m_afEntry[ 6]*rkM.m_afEntry[ 5] +
			m_afEntry[ 7]*rkM.m_afEntry[10] +
			m_afEntry[ 8]*rkM.m_afEntry[15] +
			m_afEntry[ 9]*rkM.m_afEntry[20],

			m_afEntry[ 5]*rkM.m_afEntry[ 1] +
			m_afEntry[ 6]*rkM.m_afEntry[ 6] +
			m_afEntry[ 7]*rkM.m_afEntry[11] +
			m_afEntry[ 8]*rkM.m_afEntry[16] +
			m_afEntry[ 9]*rkM.m_afEntry[21],

			m_afEntry[ 5]*rkM.m_afEntry[ 2] +
			m_afEntry[ 6]*rkM.m_afEntry[ 7] +
			m_afEntry[ 7]*rkM.m_afEntry[12] +
			m_afEntry[ 8]*rkM.m_afEntry[17] +
			m_afEntry[ 9]*rkM.m_afEntry[22],

			m_afEntry[ 5]*rkM.m_afEntry[ 3] +
			m_afEntry[ 6]*rkM.m_afEntry[ 8] +
			m_afEntry[ 7]*rkM.m_afEntry[13] +
			m_afEntry[ 8]*rkM.m_afEntry[18] +
			m_afEntry[ 9]*rkM.m_afEntry[23],

			m_afEntry[ 5]*rkM.m_afEntry[ 4] +
			m_afEntry[ 6]*rkM.m_afEntry[ 9] +
			m_afEntry[ 7]*rkM.m_afEntry[14] +
			m_afEntry[ 8]*rkM.m_afEntry[19] +
			m_afEntry[ 9]*rkM.m_afEntry[24],

			m_afEntry[ 10]*rkM.m_afEntry[ 0] +
			m_afEntry[ 11]*rkM.m_afEntry[ 5] +
			m_afEntry[ 12]*rkM.m_afEntry[10] +
			m_afEntry[ 13]*rkM.m_afEntry[15] +
			m_afEntry[ 14]*rkM.m_afEntry[20],

			m_afEntry[ 10]*rkM.m_afEntry[ 1] +
			m_afEntry[ 11]*rkM.m_afEntry[ 6] +
			m_afEntry[ 12]*rkM.m_afEntry[11] +
			m_afEntry[ 13]*rkM.m_afEntry[16] +
			m_afEntry[ 14]*rkM.m_afEntry[21],
			
			m_afEntry[ 10]*rkM.m_afEntry[ 2] +
			m_afEntry[ 11]*rkM.m_afEntry[ 7] +
			m_afEntry[ 12]*rkM.m_afEntry[12] +
			m_afEntry[ 13]*rkM.m_afEntry[17] +
			m_afEntry[ 14]*rkM.m_afEntry[22],

			m_afEntry[ 10]*rkM.m_afEntry[ 3] +
			m_afEntry[ 11]*rkM.m_afEntry[ 8] +
			m_afEntry[ 12]*rkM.m_afEntry[13] +
			m_afEntry[ 13]*rkM.m_afEntry[18] +
			m_afEntry[ 14]*rkM.m_afEntry[23],

			m_afEntry[ 10]*rkM.m_afEntry[ 4] +
			m_afEntry[ 11]*rkM.m_afEntry[ 9] +
			m_afEntry[ 12]*rkM.m_afEntry[14] +
			m_afEntry[ 13]*rkM.m_afEntry[19] +
			m_afEntry[ 14]*rkM.m_afEntry[24],

			m_afEntry[ 15]*rkM.m_afEntry[ 0] +
			m_afEntry[ 16]*rkM.m_afEntry[ 5] +
			m_afEntry[ 17]*rkM.m_afEntry[10] +
			m_afEntry[ 18]*rkM.m_afEntry[15] +
			m_afEntry[ 19]*rkM.m_afEntry[20],

			m_afEntry[ 15]*rkM.m_afEntry[ 1] +
			m_afEntry[ 16]*rkM.m_afEntry[ 6] +
			m_afEntry[ 17]*rkM.m_afEntry[11] +
			m_afEntry[ 18]*rkM.m_afEntry[16] +
			m_afEntry[ 19]*rkM.m_afEntry[21],

			m_afEntry[ 15]*rkM.m_afEntry[ 2] +
			m_afEntry[ 16]*rkM.m_afEntry[ 7] +
			m_afEntry[ 17]*rkM.m_afEntry[12] +
			m_afEntry[ 18]*rkM.m_afEntry[17] +
			m_afEntry[ 19]*rkM.m_afEntry[22],

			m_afEntry[ 15]*rkM.m_afEntry[ 3] +
			m_afEntry[ 16]*rkM.m_afEntry[ 8] +
			m_afEntry[ 17]*rkM.m_afEntry[13] +
			m_afEntry[ 18]*rkM.m_afEntry[18] +
			m_afEntry[ 19]*rkM.m_afEntry[23],

			m_afEntry[ 15]*rkM.m_afEntry[ 4] +
			m_afEntry[ 16]*rkM.m_afEntry[ 9] +
			m_afEntry[ 17]*rkM.m_afEntry[14] +
			m_afEntry[ 18]*rkM.m_afEntry[19] +
			m_afEntry[ 19]*rkM.m_afEntry[24],

			m_afEntry[ 20]*rkM.m_afEntry[ 0] +
			m_afEntry[ 21]*rkM.m_afEntry[ 5] +
			m_afEntry[ 22]*rkM.m_afEntry[10] +
			m_afEntry[ 23]*rkM.m_afEntry[15] +
			m_afEntry[ 24]*rkM.m_afEntry[20],

			m_afEntry[ 20]*rkM.m_afEntry[ 1] +
			m_afEntry[ 21]*rkM.m_afEntry[ 6] +
			m_afEntry[ 22]*rkM.m_afEntry[11] +
			m_afEntry[ 23]*rkM.m_afEntry[16] +
			m_afEntry[ 24]*rkM.m_afEntry[21],

			m_afEntry[ 20]*rkM.m_afEntry[ 2] +
			m_afEntry[ 21]*rkM.m_afEntry[ 7] +
			m_afEntry[ 22]*rkM.m_afEntry[12] +
			m_afEntry[ 23]*rkM.m_afEntry[17] +
			m_afEntry[ 24]*rkM.m_afEntry[22],

			m_afEntry[ 20]*rkM.m_afEntry[ 3] +
			m_afEntry[ 21]*rkM.m_afEntry[ 8] +
			m_afEntry[ 22]*rkM.m_afEntry[13] +
			m_afEntry[ 23]*rkM.m_afEntry[18] +
			m_afEntry[ 24]*rkM.m_afEntry[23],

			m_afEntry[ 20]*rkM.m_afEntry[ 4] +
			m_afEntry[ 21]*rkM.m_afEntry[ 9] +
			m_afEntry[ 22]*rkM.m_afEntry[14] +
			m_afEntry[ 23]*rkM.m_afEntry[19] +
			m_afEntry[ 24]*rkM.m_afEntry[24]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Matrix5<Real> Matrix5<Real>::operator* (Real fScalar) const
	{
		return Matrix5<Real>(
			fScalar*m_afEntry[ 0],
			fScalar*m_afEntry[ 1],
			fScalar*m_afEntry[ 2],
			fScalar*m_afEntry[ 3],
			fScalar*m_afEntry[ 4],
			fScalar*m_afEntry[ 5],
			fScalar*m_afEntry[ 6],
			fScalar*m_afEntry[ 7],
			fScalar*m_afEntry[ 8],
			fScalar*m_afEntry[ 9],
			fScalar*m_afEntry[10],
			fScalar*m_afEntry[11],
			fScalar*m_afEntry[12],
			fScalar*m_afEntry[13],
			fScalar*m_afEntry[14],
			fScalar*m_afEntry[15],
			fScalar*m_afEntry[16],
			fScalar*m_afEntry[17],
			fScalar*m_afEntry[18],
			fScalar*m_afEntry[19],
			fScalar*m_afEntry[20],
			fScalar*m_afEntry[21],
			fScalar*m_afEntry[22],
			fScalar*m_afEntry[23],
			fScalar*m_afEntry[24]);
	}
	//----------------------------------------------------------------------------

	template <class Real>
	inline Vector5<Real> Matrix5<Real>::operator* (const Vector5<Real>& rkV) const
	{
		return Vector5<Real>(
			m_afEntry[ 0]*rkV[0] +
			m_afEntry[ 1]*rkV[1] +
			m_afEntry[ 2]*rkV[2] +
			m_afEntry[ 3]*rkV[3] +
			m_afEntry[ 4]*rkV[4],

			m_afEntry[ 5]*rkV[0] +
			m_afEntry[ 6]*rkV[1] +
			m_afEntry[ 7]*rkV[2] +
			m_afEntry[ 8]*rkV[3] +
			m_afEntry[ 9]*rkV[4],

			m_afEntry[10]*rkV[0] +
			m_afEntry[11]*rkV[1] +
			m_afEntry[12]*rkV[2] +
			m_afEntry[13]*rkV[3] +
			m_afEntry[14]*rkV[4],

			m_afEntry[15]*rkV[0] +
			m_afEntry[16]*rkV[1] +
			m_afEntry[17]*rkV[2] +
			m_afEntry[18]*rkV[3] +
			m_afEntry[19]*rkV[4],

			m_afEntry[20]*rkV[0] +
			m_afEntry[21]*rkV[1] +
			m_afEntry[22]*rkV[2] +
			m_afEntry[23]*rkV[3] +
			m_afEntry[24]*rkV[4]);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	int Matrix5<Real>::CompareArrays (const Matrix5& rkM) const
	{
		return memcmp(m_afEntry,rkM.m_afEntry,25*sizeof(Real));
	}
	//----------------------------------------------------------------------------
	template <class Real>
	bool Matrix5<Real>::operator== (const Matrix5& rkM) const
	{
		return CompareArrays(rkM) == 0;
	}

	//----------------------------------------------------------------------------
	template <class Real>
	Matrix5<Real> Matrix5<Real>::Inverse () const
	{
		Matrix5 kInM(*this);

		int actualsize = 5, maxsize = 5;
		for (int i = 1; i < actualsize; i++) 
			kInM.m_afEntry[i] /= kInM.m_afEntry[0]; // normalize row 0

		for (int i = 1; i < actualsize; i++)  { 
			for (int j = i; j < actualsize; j++)  { // do a column of L
				Real sum = 0.0;
				for (int k = 0; k < i; k++)  
					sum += kInM.m_afEntry[j * maxsize + k] * kInM.m_afEntry[k * maxsize + i];
				kInM.m_afEntry[j * maxsize + i] -= sum;
			}
			if (i == actualsize - 1) continue;
			for (int j = i + 1; j < actualsize; j++)  {  // do a row of U
				Real sum = 0.0;
				for (int k = 0; k < i; k++)
					sum += kInM.m_afEntry[i * maxsize + k]*kInM.m_afEntry[k * maxsize + j];
				kInM.m_afEntry[i*maxsize+j] = 
					(kInM.m_afEntry[i*maxsize+j]-sum) / kInM.m_afEntry[i*maxsize+i];
			}
		}

		for ( int i = 0; i < actualsize; i++ )  // invert L
			for ( int j = i; j < actualsize; j++ )  {
				Real x = 1.0;
				if ( i != j ) {
					x = 0.0;
					for ( int k = i; k < j; k++ ) 
						x -= kInM.m_afEntry[j * maxsize + k] * kInM.m_afEntry[k * maxsize + i];
				}
				kInM.m_afEntry[j * maxsize + i] = x / kInM.m_afEntry[j * maxsize + j];
			}

		for ( int i = 0; i < actualsize; i++ )   // invert U
			for ( int j = i; j < actualsize; j++ )  {
				if ( i == j ) continue;
				Real sum = 0.0;
				for ( int k = i; k < j; k++ )
					sum += kInM.m_afEntry[k * maxsize + j] * ( (i == k) ? 1.0 : kInM.m_afEntry[i*maxsize+k] );
				kInM.m_afEntry[i * maxsize + j] = -sum;
			}
			for ( int i = 0; i < actualsize; i++ )   // final inversion
				for ( int j = 0; j < actualsize; j++ )  {
					Real sum = 0.0;
					for ( int k = ((i > j) ? i : j); k < actualsize; k++ )  
						sum += ((j == k) ? 1.0 : kInM.m_afEntry[j * maxsize + k]) * kInM.m_afEntry[k * maxsize + i];
					kInM.m_afEntry[j * maxsize + i] = sum;
				}
		return kInM;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Matrix5<Real> Matrix5<Real>::Inverse2 () const
	{
		Matrix5 kInv;
		Matrix4<Real> A11(m_afEntry[6], m_afEntry[7], m_afEntry[8], m_afEntry[9], m_afEntry[11], m_afEntry[12], m_afEntry[13], m_afEntry[14], 
			m_afEntry[16], m_afEntry[17], m_afEntry[18], m_afEntry[19], m_afEntry[21], m_afEntry[22], m_afEntry[23], m_afEntry[24]);
		Matrix4<Real> A12(m_afEntry[5], m_afEntry[7], m_afEntry[8], m_afEntry[9], m_afEntry[10], m_afEntry[12], m_afEntry[13], m_afEntry[14], 
			m_afEntry[15], m_afEntry[17], m_afEntry[18], m_afEntry[19], m_afEntry[20], m_afEntry[22], m_afEntry[23], m_afEntry[24]);
		Matrix4<Real> A13(m_afEntry[5], m_afEntry[6], m_afEntry[8], m_afEntry[9], m_afEntry[10], m_afEntry[11], m_afEntry[13], m_afEntry[14], 
			m_afEntry[15], m_afEntry[16], m_afEntry[18], m_afEntry[19], m_afEntry[20], m_afEntry[21], m_afEntry[23], m_afEntry[24]);
		Matrix4<Real> A14(m_afEntry[5], m_afEntry[6], m_afEntry[7], m_afEntry[9], m_afEntry[10], m_afEntry[11], m_afEntry[12], m_afEntry[14], 
			m_afEntry[15], m_afEntry[16], m_afEntry[17], m_afEntry[19], m_afEntry[20], m_afEntry[21], m_afEntry[22], m_afEntry[24]);
		Matrix4<Real> A15(m_afEntry[5], m_afEntry[6], m_afEntry[7], m_afEntry[8], m_afEntry[10], m_afEntry[11], m_afEntry[12], m_afEntry[13], 
			m_afEntry[15], m_afEntry[16], m_afEntry[17], m_afEntry[18], m_afEntry[20], m_afEntry[21], m_afEntry[22], m_afEntry[23]);

		Real fDet = m_afEntry[0]*A11.Determinant() - m_afEntry[1]*A12.Determinant() + m_afEntry[2]*A13.Determinant()
			- m_afEntry[3]*A14.Determinant() + m_afEntry[4]*A15.Determinant();

		if (Math<Real>::FAbs(fDet) <= Math<Real>::INVERSE_TOLERANCE)
		{
			return kInv;
		}

		std::vector<Matrix4<Real>> adjoint_matrix;
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 5; col++)
			{
				Real entries[16];
				int index = 0;
				for (int i = 0; i < 5; i++)
					for (int j = 0; j < 5; j++)
						if (i != row && j != col)
							entries[index++] = m_afEntry[i + 5 * j];

				adjoint_matrix.push_back(
					Matrix4<Real>(entries[0], entries[1], entries[2], entries[3], 
					entries[4],	entries[5], entries[6], entries[7],
					entries[8], entries[9], entries[10], entries[11],
					entries[12], entries[13], entries[14], entries[15]));
			}
		
		for (int i = 0; i < 25; i++)
		{
			kInv.m_afEntry[i] = pow(-1.0, i) * adjoint_matrix[i].Determinant() / fDet;
		}

		return kInv;
	}
	//----------------------------------------------------------------------------
	template <class Real>
	Real Matrix5<Real>::QForm (const Vector5<Real>& rkU,
		const Vector5<Real>& rkV) const
	{
		return rkU.Dot((*this)*rkV);
	}
	//----------------------------------------------------------------------------
	template <class Real>
	inline Vector5<Real> operator* (const Vector5<Real>& rkV,
		const Matrix5<Real>& rkM)
	{
		return Vector5<Real>(
			rkV[0]*rkM[0][0]+rkV[1]*rkM[1][0]+rkV[2]*rkM[2][0]+rkV[3]*rkM[3][0]+rkV[4]*rkM[4][0],
			rkV[0]*rkM[0][1]+rkV[1]*rkM[1][1]+rkV[2]*rkM[2][1]+rkV[3]*rkM[3][1]+rkV[4]*rkM[4][1],
			rkV[0]*rkM[0][2]+rkV[1]*rkM[1][2]+rkV[2]*rkM[2][2]+rkV[3]*rkM[3][2]+rkV[4]*rkM[4][2],
			rkV[0]*rkM[0][3]+rkV[1]*rkM[1][3]+rkV[2]*rkM[2][3]+rkV[3]*rkM[3][3]+rkV[4]*rkM[4][3],
			rkV[0]*rkM[0][4]+rkV[1]*rkM[1][4]+rkV[2]*rkM[2][4]+rkV[3]*rkM[3][4]+rkV[4]*rkM[4][4]);
	}
	//----------------------------------------------------------------------------


	typedef Matrix5<float> Matrix5f;
	typedef Matrix5<double> Matrix5d;

}



// v0^T m0 v1
// when m0 is a symmetric matrix, the order of v0 and v1 do not matter.
inline double VMV(Wm4::Vector3d & v0, Wm4::Matrix3d & m0, Wm4::Vector3d & v1)
{
	return v0.Dot(m0*v1);
}

// v0^T m0 v1
// when m0 is a symmetric matrix, the order of v0 and v1 do not matter.
inline double VMV(Wm4::Vector3d & v0, Wm4::Vector3d & v1, Wm4::Matrix3d & m0)
{
	return v0.Dot(m0*v1);
}

// v0^T m0 v1
// when m0 is a symmetric matrix, the order of v0 and v1 do not matter.
inline double VMV(Wm4::Matrix3d & m0, Wm4::Vector3d & v0, Wm4::Vector3d & v1)
{
	return v0.Dot(m0*v1);
}





#endif
