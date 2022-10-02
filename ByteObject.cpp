/*
** EPITECH PROJECT, 2020
** YEP Indie Studio
** File description:
** Indie Object Class
*/

#include <ostream>
#include <fstream>
#include "ByteObject.hpp"

std::basic_istream<char> &operator>>(std::basic_istream<char> &stream, ByteObject &obj) {
    stream.read((char *)&obj._size, sizeof(obj._size));
    obj._data.resize((int)obj._size);
    stream.read((char *)&(obj._data[0]), obj._size);
    obj.setCursor(0);
    return stream;
}

std::basic_ostream<char> &operator<<(std::basic_ostream<char> &stream, const ByteObject &obj)
{
    stream.write((char *)&obj._size, sizeof(obj._size));
    stream.write((char *)&(obj._data[0]), obj._size);
    return stream;
}   