#pragma once

#include "../utils.hpp"


template <typename T>
class Mat4 {
public:
    Mat4() {
        mat[0] = 1;
        mat[5] = 1;
        mat[10] = 1;
        mat[15] = 1;
    }

    T operator[](const int i) const {
        return mat[i];
    }

    T &operator[](const int i) {
        return mat[i];
    }

    Mat4 operator*(const Mat4 &other) const {
        Mat4 res;
        for (int c = 0; c < 4; c++)
            for (int r = 0; r < 4; r++)
                res[r*4 + c] = mat[r*4 + 0] * other[0*4 + c] + mat[r*4 + 1] * other[1*4 + c] + mat[r*4 + 2] * other[2*4 + c] + mat[r*4 + 3] * other[3*4 + c];
        return res;
    }

    NDVector<T, 4> operator*(NDVector<T, 4> const &p1) const {
        NDVector<T, 4> res;
        res[0] = p1[0]*mat[0*4+0] + p1[1]*mat[1*4+0] + p1[2]*mat[2*4+0] + p1[3]*mat[3*4+0];
        res[1] = p1[0]*mat[0*4+1] + p1[1]*mat[1*4+1] + p1[2]*mat[2*4+1] + p1[3]*mat[3*4+1];
        res[2] = p1[0]*mat[0*4+2] + p1[1]*mat[1*4+2] + p1[2]*mat[2*4+2] + p1[3]*mat[3*4+2];
        res[3] = p1[0]*mat[0*4+3] + p1[1]*mat[1*4+3] + p1[2]*mat[2*4+3] + p1[3]*mat[3*4+3];
        return res;
    }
    NDVector<T, 3> operator*(NDVector<T, 3> const &p1) const {
        NDVector<T, 3> res;
        res[0] = p1[0]*mat[0*4+0] + p1[1]*mat[1*4+0] + p1[2]*mat[2*4+0] + 1*mat[3*4+0];
        res[1] = p1[0]*mat[0*4+1] + p1[1]*mat[1*4+1] + p1[2]*mat[2*4+1] + 1*mat[3*4+1];
        res[2] = p1[0]*mat[0*4+2] + p1[1]*mat[1*4+2] + p1[2]*mat[2*4+2] + 1*mat[3*4+2];
        return res;
    }

    void pointAt(NDVector<T, 3> const &pos, NDVector<T, 3> const &forward, NDVector<T, 3> const &up) {
        NDVector<double, 3> newForward2 = forward.normalize();
        NDVector<T, 3> _forward = newForward2.cast<T>();

        NDVector<T, 3> a = _forward * up.dotProduct(_forward);
        NDVector<double, 3> newUp2 = (up - a).normalize();
        NDVector<T, 3> newUp = newUp2.cast<T>();

        NDVector<T, 3> newRight = newUp.crossProduct(_forward);


        mat[0*4+0] = newRight[0];
        mat[0*4+1] = newRight[1];
        mat[0*4+2] = newRight[2];
        mat[0*4+3] = 0.0f;

        mat[1*4+0] = newUp[0];
        mat[1*4+1] = newUp[1];
        mat[1*4+2] = newUp[2];
        mat[1*4+3] = 0.0f;

        mat[2*4+0] = newForward2[0];
        mat[2*4+1] = newForward2[1];
        mat[2*4+2] = newForward2[2];
        mat[2*4+3] = 0.0f;

        mat[3*4+0] = pos[0];
        mat[3*4+1] = pos[1];
        mat[3*4+2] = pos[2];
        mat[3*4+3] = 1.0f;
    }

    void rx(T a) {
        Mat4 tmpMat;
        tmpMat[1*4+1] = cosf(a);
        tmpMat[1*4+2] = sinf(a);
        tmpMat[2*4+1] = -sinf(a);
        tmpMat[2*4+2] = cosf(a);
        *this = tmpMat * *this;
    }

    void ry(T a) {
        Mat4 tmpMat;
        tmpMat[0*4+0] = cosf(a);
        tmpMat[0*4+2] = sinf(a);
        tmpMat[2*4+0] = -sinf(a);
        tmpMat[2*4+2] = cosf(a);
        *this = tmpMat * *this;
    }

    void rz(T a) {
        Mat4 tmpMat;
        tmpMat[0*4+0] = cosf(a);
        tmpMat[0*4+1] = sinf(a);
        tmpMat[1*4+0] = -sinf(a);
        tmpMat[1*4+1] = cosf(a);
        *this = tmpMat * *this;
    }

    void rrx(T a) {
        Mat4 tmpMat;
        tmpMat[1*4+1] = cosf(a);
        tmpMat[1*4+2] = sinf(a);
        tmpMat[2*4+1] = -sinf(a);
        tmpMat[2*4+2] = cosf(a);
        *this = *this * tmpMat;
    }

    void rry(T a) {
        Mat4 tmpMat;
        tmpMat[0*4+0] = cosf(a);
        tmpMat[0*4+2] = sinf(a);
        tmpMat[2*4+0] = -sinf(a);
        tmpMat[2*4+2] = cosf(a);
        *this = *this * tmpMat;
    }

    void rrz(T a) {
        Mat4 tmpMat;
        tmpMat[0*4+0] = cosf(a);
        tmpMat[0*4+1] = sinf(a);
        tmpMat[1*4+0] = -sinf(a);
        tmpMat[1*4+1] = cosf(a);
        *this = *this * tmpMat;
    }

    void t(T x, T y, T z) {
        Mat4 tmpMat;
        tmpMat[3*4+0] = x;
        tmpMat[3*4+1] = y;
        tmpMat[3*4+2] = z;
        *this = tmpMat * *this;
    }

    void tt(T x, T y, T z) {
        Mat4 tmpMat;
        tmpMat[3*4+0] = x;
        tmpMat[3*4+1] = y;
        tmpMat[3*4+2] = z;
        *this = *this * tmpMat;
    }

    /*void tx(T t) {
        Mat4 tmpMat;
        tmpMat[0*4+3] = t;
        *this = tmpMat * *this;
    }

    void ty(T t) {
        Mat4 tmpMat;
        tmpMat[3+1*4] = t;
        *this = tmpMat * *this;
    }

    void tz(T t) {
        Mat4 tmpMat;
        tmpMat[3+2*4] = t;
        *this = tmpMat * *this;
    }

    void ttx(T t) {
        Mat4 tmpMat;
        tmpMat[3+0*4] = t;
        *this = *this * tmpMat;
    }

    void tty(T t) {
        Mat4 tmpMat;
        tmpMat[3+1*4] = t;
        *this = *this * tmpMat;
    }

    void ttz(T t) {
        Mat4 tmpMat;
        tmpMat[3+2*4] = t;
        *this = *this * tmpMat;
    }*/

    void scale(T s) {
        Mat4 tmpMat;
        tmpMat[0] = s;
        tmpMat[5] = s;
        tmpMat[10] = s;
        tmpMat[15] = 1;
        *this = *this * tmpMat;
    }

    void scale(T x, T y, T z) {
        Mat4 tmpMat;
        tmpMat[0] = x;
        tmpMat[5] = y;
        tmpMat[10] = z;
        tmpMat[15] = 1;
        *this = *this * tmpMat;
    }
    
    void scalescale(T s) {
        Mat4 tmpMat;
        tmpMat[0] = s;
        tmpMat[5] = s;
        tmpMat[10] = s;
        tmpMat[15] = 1;
        *this = tmpMat * *this;
    }

    void scalescale(T x, T y, T z) {
        Mat4 tmpMat;
        tmpMat[0] = x;
        tmpMat[5] = y;
        tmpMat[10] = z;
        tmpMat[15] = 1;
        *this = tmpMat * *this;
    }

    Mat4 inv() {
        Mat4 matrix;
        T inv[16], det;

        inv[0]  =  mat[5] * mat[10] * mat[15] - mat[5] * mat[11] * mat[14] - mat[9] * mat[6] * mat[15] + mat[9] * mat[7] * mat[14] + mat[13] * mat[6] * mat[11] - mat[13] * mat[7] * mat[10];
        inv[4]  = -mat[4] * mat[10] * mat[15] + mat[4] * mat[11] * mat[14] + mat[8] * mat[6] * mat[15] - mat[8] * mat[7] * mat[14] - mat[12] * mat[6] * mat[11] + mat[12] * mat[7] * mat[10];
        inv[8]  =  mat[4] * mat[9]  * mat[15] - mat[4] * mat[11] * mat[13] - mat[8] * mat[5] * mat[15] + mat[8] * mat[7] * mat[13] + mat[12] * mat[5] * mat[11] - mat[12] * mat[7] * mat[9];
        inv[12] = -mat[4] * mat[9]  * mat[14] + mat[4] * mat[10] * mat[13] + mat[8] * mat[5] * mat[14] - mat[8] * mat[6] * mat[13] - mat[12] * mat[5] * mat[10] + mat[12] * mat[6] * mat[9];
        inv[1]  = -mat[1] * mat[10] * mat[15] + mat[1] * mat[11] * mat[14] + mat[9] * mat[2] * mat[15] - mat[9] * mat[3] * mat[14] - mat[13] * mat[2] * mat[11] + mat[13] * mat[3] * mat[10];
        inv[5]  =  mat[0] * mat[10] * mat[15] - mat[0] * mat[11] * mat[14] - mat[8] * mat[2] * mat[15] + mat[8] * mat[3] * mat[14] + mat[12] * mat[2] * mat[11] - mat[12] * mat[3] * mat[10];
        inv[9]  = -mat[0] * mat[9]  * mat[15] + mat[0] * mat[11] * mat[13] + mat[8] * mat[1] * mat[15] - mat[8] * mat[3] * mat[13] - mat[12] * mat[1] * mat[11] + mat[12] * mat[3] * mat[9];
        inv[13] =  mat[0] * mat[9]  * mat[14] - mat[0] * mat[10] * mat[13] - mat[8] * mat[1] * mat[14] + mat[8] * mat[2] * mat[13] + mat[12] * mat[1] * mat[10] - mat[12] * mat[2] * mat[9];
        inv[2]  =  mat[1] * mat[6]  * mat[15] - mat[1] * mat[7]  * mat[14] - mat[5] * mat[2] * mat[15] + mat[5] * mat[3] * mat[14] + mat[13] * mat[2] * mat[7]  - mat[13] * mat[3] * mat[6];
        inv[6]  = -mat[0] * mat[6]  * mat[15] + mat[0] * mat[7]  * mat[14] + mat[4] * mat[2] * mat[15] - mat[4] * mat[3] * mat[14] - mat[12] * mat[2] * mat[7]  + mat[12] * mat[3] * mat[6];
        inv[10] =  mat[0] * mat[5]  * mat[15] - mat[0] * mat[7]  * mat[13] - mat[4] * mat[1] * mat[15] + mat[4] * mat[3] * mat[13] + mat[12] * mat[1] * mat[7]  - mat[12] * mat[3] * mat[5];
        inv[14] = -mat[0] * mat[5]  * mat[14] + mat[0] * mat[6]  * mat[13] + mat[4] * mat[1] * mat[14] - mat[4] * mat[2] * mat[13] - mat[12] * mat[1] * mat[6]  + mat[12] * mat[2] * mat[5];
        inv[3]  = -mat[1] * mat[6]  * mat[11] + mat[1] * mat[7]  * mat[10] + mat[5] * mat[2] * mat[11] - mat[5] * mat[3] * mat[10] - mat[9]  * mat[2] * mat[7]  + mat[9]  * mat[3] * mat[6];
        inv[7]  =  mat[0] * mat[6]  * mat[11] - mat[0] * mat[7]  * mat[10] - mat[4] * mat[2] * mat[11] + mat[4] * mat[3] * mat[10] + mat[8]  * mat[2] * mat[7]  - mat[8]  * mat[3] * mat[6];
        inv[11] = -mat[0] * mat[5]  * mat[11] + mat[0] * mat[7]  * mat[9]  + mat[4] * mat[1] * mat[11] - mat[4] * mat[3] * mat[9]  - mat[8]  * mat[1] * mat[7]  + mat[8]  * mat[3] * mat[5];
        inv[15] =  mat[0] * mat[5]  * mat[10] - mat[0] * mat[6]  * mat[9]  - mat[4] * mat[1] * mat[10] + mat[4] * mat[2] * mat[9]  + mat[8]  * mat[1] * mat[6]  - mat[8]  * mat[2] * mat[5];

        det = mat[0] * inv[0] + mat[1] * inv[4] + mat[2] * inv[8] + mat[3] * inv[12];
        if (det == 0) 
            return matrix;
        det = 1.0 / det;

        int i;
        for (i = 0; i < 16; i++)
            matrix[i] = inv[i] * det;

        return matrix;
    }


    T mat[16] = {0};
};

template <class T>
std::ostream &operator<<(std::ostream &stream, const Mat4<T> &that) {
    stream << typeid(T).name() << std::endl;
    for (int i = 0; i < 16; i++) {
        stream << that[i] << " ";
        if (!((i+1) % 4)) 
            stream << std::endl;
    }
    return stream;
}