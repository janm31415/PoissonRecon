#pragma once

#include <algorithm>
#include <list>
#include <unordered_map>
#include <vector>
#include <cassert>

namespace jtk
  {

  /////////////////////////////////////////
  // hashed_heap
  /////////////////////////////////////////

  namespace hashed_heap_details
    {

    template <typename _RandomAccessIterator, typename _Compare>
    inline void _up_heap(_RandomAccessIterator _first, _RandomAccessIterator _pos, _Compare _comp)
      {
      auto _index = _pos - _first;
      auto _parent = (_index - 1) / 2;
      auto _val = *_pos;

      while (_index > 0 && _comp(*(_first + _parent), _val)) {
        *(_first + _index) = *(_first + _parent);
        _index = _parent;
        _parent = (_parent - 1) / 2;
        }

      if (_pos != (_first + _index))
        *(_first + _index) = _val;
      }

    template <typename _RandomAccessIterator, typename _Compare>
    inline void _down_heap(_RandomAccessIterator _first,
      _RandomAccessIterator _last,
      _RandomAccessIterator _pos,
      _Compare _comp)
      {
      auto _len = _last - _first;
      auto _index = _pos - _first;
      auto _left = _index * 2 + 1;
      auto _right = _index * 2 + 2;
      auto _largest = _right;
      auto _val = *_pos;

      while (_index < _len) {
        if (_right < _len) {
          _largest = _comp(*(_first + _right), *(_first + _left)) ? _left : _right;
          }
        else if (_left < _len) {
          _largest = _left;
          }
        else {
          // Force termination
          _largest = _len;
          }

        if (_largest < _len && _comp(_val, *(_first + _largest))) {
          *(_first + _index) = *(_first + _largest);
          _index = _largest;
          _left = _index * 2 + 1;
          _right = _index * 2 + 2;
          }
        else
          break;
        }

      if (_pos != (_first + _index))
        *(_first + _index) = _val;
      }

    template <typename _RandomAccessIterator, typename _Compare>
    inline void _update_heap(_RandomAccessIterator _first,
      _RandomAccessIterator _last,
      _RandomAccessIterator _pos,
      _Compare _comp)
      {
      auto _index = (_pos - _first);
      auto _parent = (_index - 1) / 2;

      if (_index > 0 && _comp(*(_first + _parent), *(_pos)))
        _up_heap(_first, _pos, _comp);
      else
        _down_heap(_first, _last, _pos, _comp);
      }

    }

  template <class Key,
    class Data,
    class HashKey = std::hash<Key>,
    class KeyEquality = std::equal_to<Key>,
    class CompareData = std::less<Data>>
    class hashed_heap
    {
    private:
      using hash_type = std::unordered_map<Key, Data, HashKey, KeyEquality>;
      using heap_type = std::vector<Key>;

    public:
      using my_type = hashed_heap<Key, Data, HashKey, KeyEquality, CompareData>;
      using value_type = std::pair<Key, Data>;
      using const_reference = const std::pair<Key, Data>&;

    public:
      hashed_heap()
        : _hash()
        , _heap()
        {}

      hashed_heap(const hashed_heap& other)
        : _hash(other._hash)
        , _heap(other._heap)
        {}

      my_type& operator=(const my_type& right)
        {
        my_type temp(right);
        swap(temp);
        return (*this);
        }

      void swap(my_type& right)
        {
        _hash.swap(right._hash);
        _heap.swap(right._heap);
        }

      void reserve(std::size_t capacity) { _heap.reserve(capacity); }

      void push(const_reference val)
        {
        auto it = _hash.find(val.first);
        if (it != _hash.end() &&
          CompareData {}(val.second, it->second)) // element was already present, but with higher priority
          {
          it->second = val.second;
          auto iter =
            std::find_if(_heap.begin(), _heap.end(), [=](const Key& key) { return KeyEquality{}(key, val.first); });
          hashed_heap_details::_update_heap(_heap.begin(), _heap.end(), iter, Compare(_hash, _compare_data));
          }
        else if (it == _hash.end()) {
          _hash[val.first] = val.second;
          _heap.push_back(val.first);
          hashed_heap_details::_update_heap(_heap.begin(),
            _heap.end(),
            _heap.end() - 1,
            Compare(_hash, _compare_data)); // do not use std::push_heap, the STL does not guarantee that heaps
                                              // are implemented corresponding to _update_heap
          }
        }

      const hash_type& get_hash() const { return _hash; }

      const heap_type& get_heap() const { return _heap; }

      value_type top() const
        {
        std::pair<Key, Data> top_element(_heap.front(), _hash.at(_heap.front()));
        return top_element;
        }

      void pop() // do not use std::pop_heap, the STL does not guarantee that heaps are implemented corresponding to _update_heap
        {
        Key k = _heap.front();
        std::swap(_heap.front(), _heap.back());
        _heap.pop_back();
        if (!_heap.empty())
          hashed_heap_details::_update_heap(_heap.begin(), _heap.end(), _heap.begin(), Compare(_hash, _compare_data));
        _hash.erase(k);
        }

      bool empty() const { return _heap.empty(); }

    private:
      class Compare
        {
        private:
          hash_type& _hash;
          CompareData& _compare_data;

        public:

          Compare(const Compare& other)
            : _hash(other._hash)
            , _compare_data(other._compare_data)
            {}

          void swap(Compare& other)
            {
            std::swap(_hash, other._hash);
            std::swap(_compare_data, other._compare_data);
            }

          Compare& operator=(const Compare& other)
            {
            Compare temp(other);
            swap(temp);
            return *this;
            }

          Compare(hash_type& i_hash, CompareData& compare_data)
            : _hash(i_hash)
            , _compare_data(compare_data)
            {}

          bool operator()(const Key& left, const Key& right) const { return _compare_data(_hash[right], _hash[left]); }
        };

    private:
      hash_type _hash;
      heap_type _heap;
      CompareData _compare_data;
    };

  } // namespace jtk
