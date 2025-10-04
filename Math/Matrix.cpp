#include "Matrix.h"
#include "Vectors.h"
#include <algorithm>
#include <iomanip>

namespace Math {
Matrix33::Matrix33(Vec3 a, Vec3 b, Vec3 c)
    : x{a, b, c} {};
const Vec3 Matrix33::operator[](size_t i) const { return x[i]; }
Vec3 &Matrix33::operator[](size_t i) { return x[i]; }
Matrix33 &Matrix33::operator+=(const Matrix33 &m) {
  x[0] += m.x[0];
  x[1] += m.x[1];
  x[2] += m.x[2];
  return *this;
}
Matrix33 &Matrix33::operator-=(const Matrix33 &m) {
  x[0] -= m.x[0];
  x[1] -= m.x[1];
  x[2] -= m.x[2];
  return *this;
}
Matrix33 &Matrix33::operator*=(const float &s) {
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return *this;
}
Matrix33 Matrix33::operator+(const Matrix33 &m) const {
  Matrix33 res;
  res.x[0] = x[0] + m.x[0];
  res.x[1] = x[1] + m.x[1];
  res.x[2] = x[2] + m.x[2];
  return res;
}

Matrix33 Matrix33::operator-(const Matrix33 &m) const {
  Matrix33 res;
  res.x[0] = x[0] - m.x[0];
  res.x[1] = x[1] - m.x[1];
  res.x[2] = x[2] - m.x[2];
  return res;
}

Matrix33 Matrix33::operator*(const float &s) const {
  Matrix33 result(*this);
  std::transform(result.x.begin(), result.x.end(), result.x.begin(),
            [s](Vec3 m) { return m * s; });
  return result;
}

Matrix33 Vec3::Diagonal() const {
  return {Vec3{r[0], 0, 0}, Vec3{0, r[1], 0}, Vec3{0, 0, r[2]}};
}
void Matrix33::Off_Diagonal() {
  x[0][1] = 0.0f;
  x[0][2] = 0.0f;
  x[1][0] = 0.0f;
  x[1][2] = 0.0f;
  x[2][0] = 0.0f;
  x[2][1] = 0.0f;
}
Vec3 Matrix33::VecSum() const {
  return Vec3{x[0].sum(), x[1].sum(), x[2].sum()};
}

std::ostream &operator<<(std::ostream &os, const Matrix33 &m) {
  return os << std::setprecision(5) << std::fixed << '[' << m[0][0] << ' ' << m[1][0] << ' '
            << m[2][0] << '\n'
            << ' ' << m[0][1] << ' ' << m[1][1] << ' ' << m[2][1] << '\n'
            << ' ' << m[0][2] << ' ' << m[1][2] << ' ' << m[2][2] << "]";
}

Matrix33 Matrix33::Identity() {
  return Matrix33(Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f),
                  Vec3(0.0f, 0.0f, 1.0f));
};
} // namespace Math
