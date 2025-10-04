#include "Vectors.h"

using namespace std;

namespace Math {

float &Vec3::operator[](size_t i) { return r[i]; }
const float Vec3::operator[](size_t i) const { return r[i]; }
Vec3 &Vec3::operator+=(const Vec3 &v) {
  r[0] += v.r[0];
  r[1] += v.r[1];
  r[2] += v.r[2];
  return *this;
}
Vec3 &Vec3::operator-=(const Vec3 &v) {
  r[0] -= v.r[0];
  r[1] -= v.r[1];
  r[2] -= v.r[2];
  return *this;
}
Vec3 &Vec3::operator*=(const float s) {
  r[0] *= s;
  r[1] *= s;
  r[2] *= s;
  return *this;
};
Vec3 &Vec3::operator/=(const float s) {
  r[0] /= s;
  r[1] /= s;
  r[2] /= s;
  return *this;
};
Vec3 Vec3::operator+(const Vec3 &v) const {
  return Vec3(r[0] + v.r[0], r[1] + v.r[1], r[2] + v.r[2]);
}
Vec3 Vec3::operator-(const Vec3 &v) const {
  return Vec3(r[0] - v.r[0], r[1] - v.r[1], r[2] - v.r[2]);
}
Vec3 Vec3::operator*(const float s) const {
  return Vec3(s * r[0], s * r[1], s * r[2]);
}
Vec3 Vec3::operator/(const float s) const {
  return Vec3(r[0] / s, r[1] / s, r[2] / s);
}

float Vec3::dot(const Vec3 &v) const {
  return r[0] * v.r[0] + r[1] * v.r[1] + r[2] * v.r[2];
};
float Vec3::operator*(const Vec3 &v) const {
  return r[0] * v.r[0] + r[1] * v.r[1] + r[2] * v.r[2];
};
Vec3 Vec3::cross(const Vec3 &v) const {
  return Vec3(r[1] * v.r[2] - r[2] * v.r[1], r[2] * v.r[0] - r[0] * v.r[2],
              r[0] * v.r[1] - r[1] * v.r[0]);
};
float Vec3::norm() const {
  return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
};
float Vec3::norm2() const { return r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; };
float Vec3::sum() const { return r[0] + r[1] + r[2]; };
float Vec3::product() const { return r[0] * r[1] * r[2]; };
Vec3 Vec3::normalize() const {
  Vec3 u;
  float n = norm();
  if (n != 0) {
    float factor = 1.0 / n;
    u.r[0] = r[0] * factor;
    u.r[1] = r[1] * factor;
    u.r[2] = r[2] * factor;
  }

  return u;
}

istream &operator>>(istream &is, Vec3 &v) {
  return is >> v.r[0] >> v.r[1] >> v.r[2];
}
ostream &operator<<(ostream &os, const Vec3 &v) {
  return os << '(' << v.r[0] << ' ' << v.r[1] << ' ' << v.r[2] << ')';
}
Vec3 Vec3::Identity() { return Vec3(1, 1, 1); }
} // namespace Math
