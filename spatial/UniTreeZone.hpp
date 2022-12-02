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
#include <shared_mutex>
#include <atomic>

template <class TKey, class TData, uint TD>
class UniTreeZone {
public:
    class Storage {
    public:
        Storage(Zone<TKey, TD> key, TData *data) : _key(key), _data(data) {}
        
        inline uint64_t isLockMask(uint64_t const &i) const {
            return i & isAlreadyCheck;
        }

        inline void lockMask(uint64_t const &i) const {
            isAlreadyCheck |= i;
        }

        inline void unlockMask(uint64_t const &i) const {
            isAlreadyCheck &= ~i;   
        }

        Zone<TKey, TD> const _key;
        TData * const _data;
    private:
        mutable std::atomic<uint64_t> isAlreadyCheck = 0;
    };
    class _UniTreeZone {
    public:
        

        _UniTreeZone(Zone<TKey, TD> const &zone) :
            _zone(zone)
        {}

        ~_UniTreeZone() {
            for (int i = 0; i < (int)std::pow(2, TD); ++i)
                if (_nexts[i])
                    delete _nexts[i];
        }

        void _addData(uint64_t const &threadMask, TKey const &max, std::shared_ptr<Storage> &storage) {
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

                            _nexts[cmpoff] = new _UniTreeZone(zone);
                        }
                        _nexts[cmpoff]->_addData(threadMask, max, storage);
                    }
                }
            }
        }

        void _removeData(uint64_t const &threadMask, Zone<TKey, TD> const &key, TData *data) {
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
                    _nexts[cmpoff]->_removeData(threadMask, key, data);
        }

        void _getColides(uint64_t const &threadMask, Zone<TKey, TD> const &zone, TKey const &minSize, std::shared_ptr<std::vector<std::shared_ptr<Storage>>> &res) const {
            if (!_zone.intersect(zone) || _zone.size[0] < minSize)
                return;
            for (std::shared_ptr<Storage> data : _datas)
                if (data->_key.intersect(zone) && data->_key.size[0] >= minSize && !data->isLockMask(threadMask)) {
                    data->lockMask(threadMask);
                    res->push_back(data);
                }
            for (uint cmpoff = 0; cmpoff < std::pow(2, TD); ++cmpoff)
                if (_nexts[cmpoff])
                    _nexts[cmpoff]->_getColides(threadMask, zone, minSize, res);
        }

        void _forEach(uint64_t const &threadMask, std::function<void(_UniTreeZone const &treeNode, uint depth)> const &func, uint d) const {
            func(*this, d);
            for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
                if (_nexts[i])
                    _nexts[i]->_forEach(threadMask, func, d+1);
        }

        void _forEachReverse(uint64_t const &threadMask, std::function<void(_UniTreeZone const &treeNode, uint depth)> const &func, uint d) const {
            for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
                if (_nexts[i])
                    _nexts[i]->_forEachReverse(threadMask, func, d+1);
            func(*this, d);
        }

        void _forEachStorage(uint64_t const &threadMask, std::function<void(_UniTreeZone const &treeNode, Storage const &data, uint depth)> const &func, uint d) const {
            for (std::shared_ptr<Storage> const data : _datas)
                func(*this, *data, d);
            for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
                if (_nexts[i])
                    _nexts[i]->_forEachStorage(threadMask, func, d+1);
        }


        void _forEach(uint64_t const &threadMask, std::function<void(_UniTreeZone const &treeNode, Zone<TKey, TD> const &zone, TData const &data, uint depth)> const &func, uint d) const {
            for (std::shared_ptr<Storage> const data : _datas)
                if (!data->isLockMask(threadMask)) {
                    func(*this, data->_key, *data->_data, d);
                    data->lockMask(threadMask);
                }
            for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
                if (_nexts[i])
                    _nexts[i]->_forEach(threadMask, func, d+1);
        }

        void _forEachReverse(uint64_t const &threadMask, std::function<void(_UniTreeZone const &treeNode, Zone<TKey, TD> const &zone, TData const &data, uint depth)> const &func, uint d) const {
            for (uint i = 0; i < (uint)std::pow(2, TD); ++i)
                if (_nexts[i])
                    _nexts[i]->_forEachReverse(threadMask, func, d+1);
            for (std::shared_ptr<Storage> const data : _datas) {
                if (!data->isLockMask(threadMask)) {
                    func(*this, data->_key, *data->_data, d);
                    data->lockMask(threadMask);
                }
            }
        }

        Zone<TKey, TD> const _zone;

        std::vector<std::shared_ptr<Storage>> _datas;
        _UniTreeZone *_nexts[(int)std::pow(2, TD)] = {0};
    };


    UniTreeZone(Zone<TKey, TD> zone)
    {
        TKey max = zone.size[0];
        for (uint i = 1; i < TD; ++i)
            if (zone.size[i] > max)
                max = zone.size[i];
        for (uint i = 1; i < TD; ++i)
                zone.size[i] = max;

        _tree = new _UniTreeZone(zone);
    }

    ~UniTreeZone() {
        delete _tree;
    }

    uint64_t lockThreadId() const {
        std::lock_guard<std::mutex> lock(_threadIdsMut);
        for (uint i = 0; 1 ; i++, i %= sizeof(_threadsIds)*8) {
            if (!(_threadsIds >> i & (((uint64_t)1) << i))) {
                _threadsIds |= (((uint64_t)1) << i);
                return (((uint64_t)1) << i);
            }
        }
    }

    void unlockThreadId(uint64_t const &mask) const {
        std::lock_guard<std::mutex> lock(_threadIdsMut);
        if (!(_threadsIds & mask))
            throw std::logic_error("thread mask already deleted");
        _threadsIds ^= mask;
    }


    void addData(Zone<TKey, TD> key, TData *data) {
        auto max = std::abs(key.size[0]);
        for (uint i = 1; i < TD; ++i)
            if (std::abs(key.size[i]) > max)
                max = std::abs(key.size[i]);

        std::shared_ptr<Storage> storage = std::make_shared<Storage>(key, data);
        std::unique_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_addData(threadMask, max, storage);
        unlockThreadId(threadMask);
    }

    void removeData(Zone<TKey, TD> key, TData *data) {
        std::unique_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_removeData(threadMask, key, data);
        unlockThreadId(threadMask);
    }


    std::shared_ptr<std::vector<std::shared_ptr<Storage>>> getColides(Zone<TKey, TD> const &zone, TKey minSize = 0, std::shared_ptr<std::vector<std::shared_ptr<Storage>>> res = std::make_shared<std::vector<std::shared_ptr<Storage>>>()) const {
        std::shared_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_getColides(threadMask, zone, minSize, res);
        for (std::shared_ptr<Storage> const data : *res)
            data->unlockMask(threadMask);
        unlockThreadId(threadMask);
        return res;
    }

    void forEachStorage(std::function<void(_UniTreeZone const &treeNode, Storage const &data, uint depth)> const &func, uint d = 0) const {
        std::shared_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_forEachStorage(threadMask, func, d);
        unlockThreadId(threadMask);
    }

    void forEach(std::function<void(_UniTreeZone const &treeNode, uint depth)> const &func, uint d = 0) const {
        std::shared_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_forEach(threadMask, func, d);
        unlockThreadId(threadMask);
    }

    void forEachReverse(std::function<void(_UniTreeZone const &treeNode, uint depth)> const &func, uint d = 0) const {
        std::shared_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_forEachReverse(threadMask, func, d);
        unlockThreadId(threadMask);
    }

    void forEach(std::function<void(_UniTreeZone const &treeNode, Zone<TKey, TD> const &zone, TData const &data, uint depth)> const &func, uint d = 0) const {
        std::shared_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_forEach(threadMask, func, d);
        _tree->_forEachStorage(threadMask, [&threadMask](_UniTreeZone const &, Storage const &storage, uint) {storage.unlockMask(threadMask);}, 0);
        unlockThreadId(threadMask);
    }

    void forEachReverse(std::function<void(_UniTreeZone const &treeNode, Zone<TKey, TD> const &zone, TData const &data, uint depth)> const &func, uint d = 0) const {
        std::shared_lock lock(_smut);
        uint64_t const threadMask = lockThreadId();
        _tree->_forEachReverse(threadMask, func, d);
        _tree->_forEachStorage(threadMask, [&threadMask](_UniTreeZone const &, Storage const &storage, uint) {;storage.unlockMask(threadMask);}, 0);
        unlockThreadId(threadMask);
    }

    mutable std::mutex _threadIdsMut;
    mutable uint64_t _threadsIds = 0;
    mutable std::shared_mutex _smut;
    _UniTreeZone *_tree;
};