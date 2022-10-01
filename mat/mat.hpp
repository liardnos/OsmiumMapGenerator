#pragma once

#include "../utils.hpp"

class Mat3 {
public:
    Mat3();

    float operator[](const int i) const;

    float &operator[](const int i);

    Mat3 operator*(const Mat3 &other) const;
    Vector2f operator*(const Vector2f &p1) const;

    void rx(float a);
    void ry(float a);
    void rz(float a);
    void rrx(float a);
    void rry(float a);
    void rrz(float a);

    void tx(float t);
    void ty(float t);
    void tz(float t);
    void ttx(float t);
    void tty(float t);
    void ttz(float t);

    void scale(float s);

    Mat3 inv();

    float mat[16] = {0};
};

std::ostream &operator<<(std::ostream &stream, const Mat3 &that);