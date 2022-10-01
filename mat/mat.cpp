#include "mat.hpp"

Mat3::Mat3() {
    mat[0] = 1;
    mat[5] = 1;
    mat[10] = 1;
    mat[15] = 1;
}

float Mat3::operator[](const int i) const {
    return mat[i];
}

float &Mat3::operator[](const int i) {
    return mat[i];
}

Mat3 Mat3::operator*(const Mat3 &other) const {
    Mat3 res;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            res[r*4 + c] = mat[r*4 + 0] * other[0*4 + c] + mat[r*4 + 1] * other[1*4 + c] + mat[r*4 + 2] * other[2*4 + c] + mat[r*4 + 3] * other[3*4 + c];
    return res;
}

Vector2f Mat3::operator*(const Vector2f &p1) const {
    Vector2f res;
    res[0] = p1[0] * mat[0] + p1[1] * mat[4] + 1 * mat[12];
    res[1] = p1[0] * mat[1] + p1[1] * mat[5] + 1 * mat[13];
    return res;
}

void Mat3::rx(float a) {
    Mat3 tmpMat;
    tmpMat[5] = cosf(a);
    tmpMat[6] = sinf(a);
    tmpMat[9] = -sinf(a);
    tmpMat[10] = cosf(a);
    *this = tmpMat * *this;
}

void Mat3::ry(float a) {
    Mat3 tmpMat;
    tmpMat[0] = cosf(a);
    tmpMat[2] = sinf(a);
    tmpMat[8] = -sinf(a);
    tmpMat[10] = cosf(a);
    *this = tmpMat * *this;
}

void Mat3::rz(float a) {
    Mat3 tmpMat;
    tmpMat[0] = cosf(a);
    tmpMat[1] = sinf(a);
    tmpMat[4] = -sinf(a);
    tmpMat[5] = cosf(a);
    *this = tmpMat * *this;
}

void Mat3::rrx(float a) {
    Mat3 tmpMat;
    tmpMat[5] = cosf(a);
    tmpMat[6] = sinf(a);
    tmpMat[9] = -sinf(a);
    tmpMat[10] = cosf(a);
    *this = *this * tmpMat;
}

void Mat3::rry(float a) {
    Mat3 tmpMat;
    tmpMat[0] = cosf(a);
    tmpMat[2] = sinf(a);
    tmpMat[8] = -sinf(a);
    tmpMat[10] = cosf(a);
    *this = *this * tmpMat;
}

void Mat3::rrz(float a) {
    Mat3 tmpMat;
    tmpMat[0] = cosf(a);
    tmpMat[1] = sinf(a);
    tmpMat[4] = -sinf(a);
    tmpMat[5] = cosf(a);
    *this = *this * tmpMat;
}

void Mat3::tx(float t) {
    Mat3 tmpMat;
    tmpMat[12] = t;
    *this = tmpMat * *this;
}

void Mat3::ty(float t) {
    Mat3 tmpMat;
    tmpMat[13] = t;
    *this = tmpMat * *this;
}

void Mat3::tz(float t) {
    Mat3 tmpMat;
    tmpMat[14] = t;
    *this = tmpMat * *this;
}

void Mat3::ttx(float t) {
    Mat3 tmpMat;
    tmpMat[12] = t;
    *this = *this * tmpMat;
}

void Mat3::tty(float t) {
    Mat3 tmpMat;
    tmpMat[13] = t;
    *this = *this * tmpMat;
}

void Mat3::ttz(float t) {
    Mat3 tmpMat;
    tmpMat[14] = t;
    *this = *this * tmpMat;
}

void Mat3::scale(float s) {
    Mat3 tmpMat;
    tmpMat[0] = s;
    tmpMat[5] = s;
    tmpMat[10] = s;
    tmpMat[15] = 1;
    *this = *this * tmpMat;
}

Mat3 Mat3::inv() {
    Mat3 matrix;
    float inv[16], det;

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


std::ostream &operator<<(std::ostream &stream, const Mat3 &that) {
    stream << std::endl;
    stream << std::endl;
    for (int i = 0; i < 16; i++) {
        stream << that[i] << " ";
        if (!((i+1) % 4)) 
            stream << std::endl;
    }
    return stream;
}