#ifndef Vector_h
#define Vector_h

#include <array>
#include <cmath>
#include <iostream>

namespace Math {
class Matrix33;

class Vec3 {
private:
  std::array<float, 3> r{};

public:
  /* Constructors */
  Vec3(float rx, float ry, float rz) : r{rx, ry, rz} {};
  Vec3() = default;
  Vec3(Vec3 &&v) = default;
  Vec3(const Vec3 &v) = default;
  Vec3 &operator=(Vec3 &&v) = default;
  Vec3 &operator=(const Vec3 &v) = default;
  ~Vec3() = default;

  float &operator[](size_t i);
  float operator[](size_t i) const;
  auto begin() { return r.begin(); }
  auto end()   { return r.end(); }
  auto begin() const { return r.cbegin(); }
  auto end()   const { return r.cend(); }

  // Basic operations
  Vec3 &operator+=(const Vec3 &v);
  Vec3 &operator-=(const Vec3 &v);
  Vec3 &operator*=(const float s);
  Vec3 &operator/=(const float s);

  Vec3 operator+(const Vec3 &v) const;
  Vec3 operator-(const Vec3 &v) const;
  Vec3 operator*(const float s) const;
  Vec3 operator/(const float s) const;

  // Vector operations
  float dot(const Vec3 &v) const;
  float operator*(const Vec3 &v) const;
  Vec3 cross(const Vec3 &v) const;
  float norm() const;
  float norm2() const;
  float sum() const;
  float product() const;

  Vec3 normalize() const;
  Vec3 rotate(float uAngle, Vec3 uAxis);

  // Print/Write
  friend std::istream &operator>>(std::istream &is, Vec3 &v);
  friend std::ostream &operator<<(std::ostream &os, const Vec3 &v);

  static Vec3 Identity();
  Matrix33 Diagonal() const;
};
} // namespace Math

#endif /* Vector_h */
