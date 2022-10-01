#pragma once

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <functional>
#include <vector>
#include <bitset>
#include <memory>
#include "../utils.hpp"

template <class T, class Zone, uint D>
class UniTreeZone {
public:
    UniTreeZone(Zone const &zone) :
        _zone(zone)
    {}

    ~UniTreeZone() {
        for (int i = 0; i < (int)std::pow(2, D); ++i)
            if (_nexts[i])
                delete _nexts[i];
    }

    void addData(T *data) {
        Zone z(_zone.pos, {0});
        bool intersect = data->intersectOnAnyD(z);
        if (intersect) { // add data to array
            _datas.push_back(data);
        } else { // pass data to next node 
            relativePos cmp = _zone.cmp(*data);
            if (!_nexts[cmp.off]) {
                Zone zone(_zone);
                zone.size /= 2;
             
                for (uint i = 0; i < D; ++i) {
                    zone.pos[i] += ((int)(((cmp.off >> i) & 1)*2)-1)*zone.size[i];
                }

                _nexts[cmp.off] = new UniTreeZone<T, Zone, D>(zone);
            }
            _nexts[cmp.off]->addData(data);
        }
    }

    void removeData(T *data) {
        Zone z(_zone.pos, {0});
        bool intersect = data->intersectOnAnyD(z);
        if (intersect) { // rm data from array
            _datas.erase(std::find(_datas.begin(), _datas.end(), data));
        } else { // pass to next node 
            relativePos cmp = _zone.cmp(*data);
            if (!_nexts[cmp.off]) {
                // this data is not in the tree
                int *p = 0;
                *p = 0;
            }
            _nexts[cmp.off]->removeData(data);
        }

    }

    std::shared_ptr<std::vector<T*>> getColides(Zone const &zone, std::shared_ptr<std::vector<T*>> res = std::make_shared<std::vector<T*>>()) const {
        getColidesR(zone, res);
        return res;
    }

    void getColidesR(Zone const &zone, std::shared_ptr<std::vector<T*>> &res) const {
        for (T *data : _datas) 
            if (data->intersect(zone))
                res->push_back(data);
        Zone z(_zone.pos, {0});
        relativePos cmp = z.cmp(zone);
        
        if (cmp.colide) {
            for (uint i = 0; i < (uint)std::pow(2, D); ++i)
                if (_nexts[i])
                    _nexts[i]->getColidesR(zone, res);
        } else
            if (_nexts[cmp.off])
                _nexts[cmp.off]->getColidesR(zone, res);
    }



    void draw(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, uint depth)> func, uint d = 0) const {
        for (uint i = 0; i < (int)std::pow(2, D); ++i)
            if (_nexts[i])
                _nexts[i]->draw(func, d+1);
        func(*this, d);
    }

    void forEach(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, uint depth)> func, uint d = 0) const {
        func(*this, d);
        for (uint i = 0; i < (uint)std::pow(2, D); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
    }

    void forEachReverse(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, uint depth)> func, uint d = 0) const {
        for (uint i = 0; i < (uint)std::pow(2, D); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
        func(*this, d);
    }

    void forEach(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, T const &data, uint depth)> func, uint d = 0) const {
        for (auto &_data : _datas)
            func(*this, *_data, d);
        for (uint i = 0; i < (uint)std::pow(2, D); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
    }

    void forEachReverse(std::function<void(UniTreeZone<T, Zone, D> const &treeNode, T const &data, uint depth)> func, uint d = 0) const {
        for (uint i = 0; i < (uint)std::pow(2, D); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
        for (auto &_data : _datas)
            func(*this, *_data, d);
    }

    Zone const _zone;

    std::vector<T *> _datas;
    UniTreeZone<T, Zone, D> *_nexts[(int)std::pow(2, D)] = {0};
};


