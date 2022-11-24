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


template <class TKey, class TData, uint TD>
class UniTreeZone {
public:
    class Storage {
    public:
        Storage(Zone<TKey, TD> key, TData *data) : _key(key), _data(data) {}
        
        Zone<TKey, TD> const _key;
        TData *_data;

        mutable bool isAlreadyCheck = false;
    };


    UniTreeZone(Zone<TKey, TD> const &zone) :
        _zone(zone)
    {}

    ~UniTreeZone() {
        for (int i = 0; i < (int)std::pow(2, TD); ++i)
            if (_nexts[i])
                delete _nexts[i];
    }

    void addData(Zone<TKey, TD> key, TData *data) {
        auto max = std::abs(key.size[0]);
        for (uint i = 1; i < TD; ++i)
            if (std::abs(key.size[i]) > max)
                max = std::abs(key.size[i]);

        Storage *storage = new Storage(key, data);
        _addData(max, storage);
    }


    void _addData(TKey max, Storage *storage) {
        if (!storage->_key.intersect(_zone))
            return;

        // Zone z(_zone.pos, {0});
        // relativePos cmp = z.cmp(*data);
        // std::cout << DEBUGVAR(max) << std::endl;

        // std::cout << DEBUGVAR(std::abs(_zone.size[0])) << std::endl;

        if (max*2 > std::abs(_zone.size[0])) { // add data to array
            _datas.push_back(storage);
        } else { // pass data to next node 
            { // find the 2**'n' node that need to be call
                // bool buf[std::pow(2, D)] = 0;
                // for (uint i = 0; i < D; ++i) {
                //     if ((cmp.colide >> i) & 1) {
                //         // split on two more
                //     } else {
                //         // go only on this side
                //         (cmp.off >> i) & 1
                //     }
                // }
                
                for (uint cmpoff = 0; cmpoff < std::pow(2, TD); ++cmpoff) { // offset the new zone}
                    if (!_nexts[cmpoff]) { // if next zone doesn't exist
                        Zone<TKey, TD> zone(_zone);
                        zone.size /= 2;
                    
                        for (uint i = 0; i < TD; ++i) { // offset the new zone
                            zone.pos[i] += ((int)(((cmpoff >> i) & 1)*2)-1)*zone.size[i];
                        }

                        _nexts[cmpoff] = new UniTreeZone<TKey, TData, TD>(zone);
                    }
                    _nexts[cmpoff]->_addData(max, storage);
                }
            }
        }
    }

    void removeData(Zone<TKey, 2> key, TData *data) {
        if (!key.intersect(_zone))
            return;


        auto it = _datas.begin();
        for (; it != _datas.end(); ++it)
            if ((*it)->_data == data) {
                _datas.erase(it);
                break;
            }

        // pass to next nodes
        for (uint cmpoff = 0; cmpoff < std::pow(2, TD); ++cmpoff)
            if (_nexts[cmpoff])
                _nexts[cmpoff]->removeData(key, data);
    }

    std::shared_ptr<std::vector<Storage *>> getColides(Zone<TKey, 2> const &zone, TKey minSize = 0, std::shared_ptr<std::vector<Storage *>> res = std::make_shared<std::vector<Storage *>>()) const {
        _getColides(zone, minSize, res);
        for (Storage const *data : *res)
            data->isAlreadyCheck = false;
        return res;
    }

    void _getColides(Zone<TKey, 2> const &zone, TKey minSize, std::shared_ptr<std::vector<Storage*>> &res) const {
        if (!_zone.intersect(zone) || _zone.size[0] < minSize)
            return;
        for (Storage *data : _datas)
            if (data->_key.intersect(zone) && data->_key.size[0] >= minSize && !data->isAlreadyCheck) {
                data->isAlreadyCheck = true;
                res->push_back(data);
            }
        for (uint cmpoff = 0; cmpoff < std::pow(2, TD); ++cmpoff)
            if (_nexts[cmpoff])
                _nexts[cmpoff]->_getColides(zone, minSize, res);
    }

    void forEach(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, uint depth)> func, uint d = 0) const {
        func(*this, d);
        for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
    }

    void forEachReverse(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, uint depth)> func, uint d = 0) const {
        for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
            if (_nexts[i])
                _nexts[i]->forEach(func, d+1);
        func(*this, d);
    }

    void forEachStorage(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, Storage const &data, uint depth)> func, uint d = 0) const {
        for (Storage const *data : _datas)
            func(*this, *data, d);
        for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
            if (_nexts[i])
                _nexts[i]->forEachStorage(func, d+1);
    }

    void forEach(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, Zone<TKey, 2> const &zone, TData const &data, uint depth)> func, uint d = 0) const {
        _forEach(func, d);
        forEachStorage([](UniTreeZone<TKey, TData, TD> const &treeNode, Storage const &storage, uint depth) {storage.isAlreadyCheck = false;});
    }

    void _forEach(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, Zone<TKey, 2> const &zone, TData const &data, uint depth)> func, uint d) const {
        for (Storage const *data : _datas)
            if (!data->isAlreadyCheck) {
                func(*this, data->_key, *data->_data, d);
                data->isAlreadyCheck = true;
            }
        for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
            if (_nexts[i])
                _nexts[i]->_forEach(func, d+1);
    }

    void forEachReverse(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, Zone<TKey, 2> const &zone, TData const &data, uint depth)> func, uint d = 0) const {
        _forEachReverse(func, d);
        forEachStorage([](UniTreeZone<TKey, TData, TD> const &treeNode, Storage const &storage, uint depth) {(void)treeNode;(void)depth;storage.isAlreadyCheck = false;});
    }

    void _forEachReverse(std::function<void(UniTreeZone<TKey, TData, TD> const &treeNode, Zone<TKey, 2> const &zone, TData const &data, uint depth)> func, uint d) const {
        for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
            if (_nexts[i])
                _nexts[i]->_forEachReverse(func, d+1);
        for (Storage const *data : _datas) {
            if (!data->isAlreadyCheck) {
               func(*this, data->_key, *data->_data, d);
                data->isAlreadyCheck = true;
            }
        }
    }

    Zone<TKey, TD> const _zone;

    std::vector<Storage *> _datas;
    UniTreeZone<TKey, TData, TD> *_nexts[(int)std::pow(2, TD)] = {0};
};