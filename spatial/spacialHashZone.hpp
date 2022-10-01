#pragma once

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <functional>
#include <vector>
#include <bitset>
#include <memory>
#include <unordered_map>
#include "../utils.hpp"

template <class T, class Zone, uint D>
class SpatialHashZone {
public:
    SpatialHashZone(Zone const &zone) {}

    ~SpatialHashZone() {}

    void addData(T *data) {}

    void removeData(T *data) {}

    std::shared_ptr<std::vector<T*>> getColides(Zone const &zone, std::shared_ptr<std::vector<T*>> res = std::make_shared<std::vector<T*>>()) const {}

    void getColidesR(Zone const &zone, std::shared_ptr<std::vector<T*>> &res) const {}

    void forEach(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, uint depth)> func, uint d = 0) const {
}

    void forEachReverse(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, uint depth)> func, uint d = 0) const {
        for (uint i = 0; i < (uint)std::pow(2, D); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
        func(*this, d);
    }

    std::vector<T *> _datas;
    UniTreeZone<T, Zone, D> *_nexts[(int)std::pow(2, D)] = {0};
};


