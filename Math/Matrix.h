#ifndef Matrix_h
#define Matrix_h

#include "Vectors.h"
#include <array>
#include <cstddef>
#include <ostream>

namespace Math {

class Matrix33 {
  std::array<Vec3, 3> x{};

public:
  // Constructors
  Matrix33(Vec3 a, Vec3 b, Vec3 c);
  Matrix33() = default;
  Matrix33(Matrix33 &&m) = default;
  Matrix33(const Matrix33 &m) = default;
  Matrix33 &operator=(Matrix33 &&m) = default;
  Matrix33 &operator=(const Matrix33 &m) = default;
  ~Matrix33() = default;

  // Subscripting
  const Vec3 operator[](size_t i) const;
  Vec3 &operator[](size_t i);
  auto begin() { return x.begin(); }
  auto end()   { return x.end(); }
  auto begin() const { return x.cbegin(); }
  auto end()   const { return x.cend(); }

  // Basic operations
  Matrix33 &operator+=(const Matrix33 &m);
  Matrix33 &operator-=(const Matrix33 &m);
  Matrix33 &operator*=(const float &s);

  Matrix33 operator+(const Matrix33 &m) const;
  Matrix33 operator-(const Matrix33 &m) const;
  Matrix33 operator*(const float &s) const;

  // Matrix operations
  void Off_Diagonal();
  Vec3 VecSum() const;

  // Print
  friend std::ostream &operator<<(std::ostream &os, const Matrix33 &m);

  static Matrix33 Identity();
};
} // namespace Math

#endif /* Matrix_h */
