// sync with
//  https://github.com/cubao/headers/blob/main/include/packedrtree.hpp
//  https://github.com/cubao/nano-fmm/blob/master/3rdparty/packedrtree.hpp

/******************************************************************************
 *
 * Project:  FlatGeobuf
 * Purpose:  Packed RTree management
 * Author:   Björn Harrtell <bjorn at wololo dot org>
 *
 ******************************************************************************
 * Copyright (c) 2018-2020, Björn Harrtell <bjorn at wololo dot org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

// NOTE: The upstream of this file is in
// https://github.com/bjornharrtell/flatgeobuf/tree/master/src/cpp

#ifndef FLATGEOBUF_PACKEDRTREE_H_
#define FLATGEOBUF_PACKEDRTREE_H_

#include <cmath>
#include <deque>
#include <numeric>
#include <cstdint>
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <map>
#include <unordered_map>

namespace FlatGeobuf
{

struct NodeItem
{
    double minX;
    double minY;
    double maxX;
    double maxY;
    uint64_t offset;
    double width() const { return maxX - minX; }
    double height() const { return maxY - minY; }
    static NodeItem sum(NodeItem a, const NodeItem &b)
    {
        a.expand(b);
        return a;
    }
    static NodeItem create(uint64_t offset = 0)
    {
        return {std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(),
                -1 * std::numeric_limits<double>::infinity(),
                -1 * std::numeric_limits<double>::infinity(), offset};
    }
    const NodeItem &expand(const NodeItem &r)
    {
        if (r.minX < minX)
            minX = r.minX;
        if (r.minY < minY)
            minY = r.minY;
        if (r.maxX > maxX)
            maxX = r.maxX;
        if (r.maxY > maxY)
            maxY = r.maxY;
        return *this;
    }
    bool intersects(const NodeItem &r) const
    {
        if (maxX < r.minX)
            return false;
        if (maxY < r.minY)
            return false;
        if (minX > r.maxX)
            return false;
        if (minY > r.maxY)
            return false;
        return true;
    }
    std::vector<double> toVector() const
    {
        return std::vector<double>{minX, minY, maxX, maxY};
    }
};

inline bool operator==(const NodeItem &lhs, const NodeItem &rhs)
{
    return lhs.minX == rhs.minX && lhs.minY == rhs.minY &&
           lhs.maxX == rhs.maxX && lhs.maxY == rhs.maxY &&
           lhs.offset == rhs.offset;
}

struct Item
{
    NodeItem nodeItem;
};

struct SearchResultItem
{
    uint64_t offset;
    uint64_t index;
};

inline bool operator==(const SearchResultItem &lhs, const SearchResultItem &rhs)
{
    return lhs.index == rhs.index && lhs.offset == rhs.offset;
}

// Based on public domain code at
// https://github.com/rawrunprotected/hilbert_curves
inline uint32_t hilbert(uint32_t x, uint32_t y)
{
    uint32_t a = x ^ y;
    uint32_t b = 0xFFFF ^ a;
    uint32_t c = 0xFFFF ^ (x | y);
    uint32_t d = x & (y ^ 0xFFFF);

    uint32_t A = a | (b >> 1);
    uint32_t B = (a >> 1) ^ a;
    uint32_t C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
    uint32_t D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;

    a = A;
    b = B;
    c = C;
    d = D;
    A = ((a & (a >> 2)) ^ (b & (b >> 2)));
    B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));
    C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
    D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));

    a = A;
    b = B;
    c = C;
    d = D;
    A = ((a & (a >> 4)) ^ (b & (b >> 4)));
    B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));
    C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
    D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));

    a = A;
    b = B;
    c = C;
    d = D;
    C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
    D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));

    a = C ^ (C >> 1);
    b = D ^ (D >> 1);

    uint32_t i0 = x ^ y;
    uint32_t i1 = b | (0xFFFF ^ (i0 | a));

    i0 = (i0 | (i0 << 8)) & 0x00FF00FF;
    i0 = (i0 | (i0 << 4)) & 0x0F0F0F0F;
    i0 = (i0 | (i0 << 2)) & 0x33333333;
    i0 = (i0 | (i0 << 1)) & 0x55555555;

    i1 = (i1 | (i1 << 8)) & 0x00FF00FF;
    i1 = (i1 | (i1 << 4)) & 0x0F0F0F0F;
    i1 = (i1 | (i1 << 2)) & 0x33333333;
    i1 = (i1 | (i1 << 1)) & 0x55555555;

    uint32_t value = ((i1 << 1) | i0);

    return value;
}

inline uint32_t hilbert(const NodeItem &r, uint32_t hilbertMax,
                        const double minX, const double minY,
                        const double width, const double height)
{
    uint32_t x = 0;
    uint32_t y = 0;
    uint32_t v;
    if (width != 0.0)
        x = static_cast<uint32_t>(
            floor(hilbertMax * ((r.minX + r.maxX) / 2 - minX) / width));
    if (height != 0.0)
        y = static_cast<uint32_t>(
            floor(hilbertMax * ((r.minY + r.maxY) / 2 - minY) / height));
    v = hilbert(x, y);
    return v;
}

constexpr uint32_t HILBERT_MAX = (1 << 16) - 1;

inline NodeItem calcExtent(const std::vector<NodeItem> &nodes)
{
    return std::accumulate(
        nodes.begin(), nodes.end(), NodeItem::create(0),
        [](NodeItem a, const NodeItem &b) { return a.expand(b); });
}

inline NodeItem calcExtent(const std::vector<std::shared_ptr<Item>> &items)
{
    return std::accumulate(items.begin(), items.end(), NodeItem::create(0),
                           [](NodeItem a, const std::shared_ptr<Item> &b) {
                               return a.expand(b->nodeItem);
                           });
}

template <class ITEM_TYPE> void hilbertSort(std::deque<ITEM_TYPE> &items)
{
    NodeItem extent = calcExtent(items);
    const double minX = extent.minX;
    const double minY = extent.minY;
    const double width = extent.width();
    const double height = extent.height();
    std::sort(
        items.begin(), items.end(),
        [minX, minY, width, height](const ITEM_TYPE &a, const ITEM_TYPE &b) {
            uint32_t ha =
                hilbert(a.nodeItem, HILBERT_MAX, minX, minY, width, height);
            uint32_t hb =
                hilbert(b.nodeItem, HILBERT_MAX, minX, minY, width, height);
            return ha > hb;
        });
}

inline void hilbertSort(std::vector<NodeItem> &items, const NodeItem &extent)
{
    const double minX = extent.minX;
    const double minY = extent.minY;
    const double width = extent.width();
    const double height = extent.height();
    std::sort(
        items.begin(), items.end(),
        [minX, minY, width, height](const NodeItem &a, const NodeItem &b) {
            uint32_t ha = hilbert(a, HILBERT_MAX, minX, minY, width, height);
            uint32_t hb = hilbert(b, HILBERT_MAX, minX, minY, width, height);
            return ha > hb;
        });
}
inline void hilbertSort(std::vector<NodeItem> &items)
{
    hilbertSort(items, calcExtent(items));
}

inline void hilbertSort(std::vector<std::shared_ptr<Item>> &items, const NodeItem &extent)
{
    const double minX = extent.minX;
    const double minY = extent.minY;
    const double width = extent.width();
    const double height = extent.height();
    std::sort(items.begin(), items.end(),
              [minX, minY, width, height](std::shared_ptr<Item> a,
                                          std::shared_ptr<Item> b) {
                  uint32_t ha = hilbert(a->nodeItem, HILBERT_MAX, minX, minY,
                                        width, height);
                  uint32_t hb = hilbert(b->nodeItem, HILBERT_MAX, minX, minY,
                                        width, height);
                  return ha > hb;
              });
}

inline void hilbertSort(std::vector<std::shared_ptr<Item>> &items) {
    return hilbertSort(items, calcExtent(items));
}


/**
 * Packed R-Tree
 * Based on https://github.com/mourner/flatbush
 */
class PackedRTree
{
    NodeItem _extent;
    std::vector<NodeItem> _nodeItems;
    uint64_t _numItems;
    uint64_t _numNodes;
    uint16_t _nodeSize;
    std::vector<std::pair<uint64_t, uint64_t>> _levelBounds;
    void init(const uint16_t nodeSize)
    {
        if (nodeSize < 2)
            throw std::invalid_argument("Node size must be at least 2");
        if (_numItems == 0)
            throw std::invalid_argument("Cannot create empty tree");
        _nodeSize = std::min(std::max(nodeSize, static_cast<uint16_t>(2)),
                             static_cast<uint16_t>(65535));
        _levelBounds = generateLevelBounds(_numItems, _nodeSize);
        _numNodes = _levelBounds.front().second;
        _nodeItems.resize(_numNodes);
    }
    void generateNodes()
    {
        for (uint32_t i = 0; i < _levelBounds.size() - 1; i++) {
            auto pos = _levelBounds[i].first;
            auto end = _levelBounds[i].second;
            auto newpos = _levelBounds[i + 1].first;
            while (pos < end) {
                NodeItem node = NodeItem::create(pos);
                for (uint32_t j = 0; j < _nodeSize && pos < end; j++)
                    node.expand(_nodeItems[pos++]);
                _nodeItems[newpos++] = node;
            }
        }
    }

    void fromData(const void *data)
    {
        auto buf = reinterpret_cast<const uint8_t *>(data);
        const NodeItem *pn = reinterpret_cast<const NodeItem *>(buf);
        for (uint64_t i = 0; i < _numNodes; i++) {
            NodeItem n = *pn++;
            _nodeItems[i] = n;
            _extent.expand(n);
        }
    }

  public:
    PackedRTree(const std::vector<NodeItem> &nodes, const NodeItem &extent,
                const uint16_t nodeSize = 16)
        : _extent(extent), _numItems(nodes.size())
    {
        init(nodeSize);
        for (size_t i = 0; i < _numItems; i++)
            _nodeItems[_numNodes - _numItems + i] = nodes[i];
        generateNodes();
    }

    PackedRTree(const std::vector<std::shared_ptr<Item>> &items,
                const NodeItem &extent, const uint16_t nodeSize = 16)
        : _extent(extent), _numItems(items.size())
    {
        init(nodeSize);
        for (size_t i = 0; i < _numItems; i++)
            _nodeItems[_numNodes - _numItems + i] = items[i]->nodeItem;
        generateNodes();
    }

    PackedRTree(const void *data, const uint64_t numItems,
                const uint16_t nodeSize)
        : _extent(NodeItem::create(0)), _numItems(numItems)
    {
        init(nodeSize);
        fromData(data);
    }

    std::vector<SearchResultItem> search(double minX, double minY, double maxX,
                                         double maxY) const
    {
        uint64_t leafNodesOffset = _levelBounds.front().first;
        NodeItem n{minX, minY, maxX, maxY, 0};
        std::vector<SearchResultItem> results;
        std::unordered_map<uint64_t, uint64_t> queue;
        queue.insert(std::pair<uint64_t, uint64_t>(0, _levelBounds.size() - 1));
        while (queue.size() != 0) {
            auto next = queue.begin();
            uint64_t nodeIndex = next->first;
            uint64_t level = next->second;
            queue.erase(next);
            bool isLeafNode = nodeIndex >= _numNodes - _numItems;
            // find the end index of the node
            uint64_t end =
                std::min(static_cast<uint64_t>(nodeIndex + _nodeSize),
                         _levelBounds[static_cast<size_t>(level)].second);
            // search through child nodes
            for (uint64_t pos = nodeIndex; pos < end; pos++) {
                auto nodeItem = _nodeItems[static_cast<size_t>(pos)];
                if (!n.intersects(nodeItem))
                    continue;
                if (isLeafNode)
                    results.push_back({nodeItem.offset, pos - leafNodesOffset});
                else
                    queue.insert(std::pair<uint64_t, uint64_t>(nodeItem.offset,
                                                               level - 1));
            }
        }
        return results;
    }
    static std::vector<std::pair<uint64_t, uint64_t>>
    generateLevelBounds(const uint64_t numItems, const uint16_t nodeSize)
    {
        if (nodeSize < 2)
            throw std::invalid_argument("Node size must be at least 2");
        if (numItems == 0)
            throw std::invalid_argument(
                "Number of items must be greater than 0");
        if (numItems >
            std::numeric_limits<uint64_t>::max() - ((numItems / nodeSize) * 2))
            throw std::overflow_error("Number of items too large");

        // number of nodes per level in bottom-up order
        std::vector<uint64_t> levelNumNodes;
        uint64_t n = numItems;
        uint64_t numNodes = n;
        levelNumNodes.push_back(n);
        do {
            n = (n + nodeSize - 1) / nodeSize;
            numNodes += n;
            levelNumNodes.push_back(n);
        } while (n != 1);

        // bounds per level in reversed storage order (top-down)
        std::vector<uint64_t> levelOffsets;
        n = numNodes;
        for (auto size : levelNumNodes)
            levelOffsets.push_back(n -= size);
        std::vector<std::pair<uint64_t, uint64_t>> levelBounds;
        for (size_t i = 0; i < levelNumNodes.size(); i++)
            levelBounds.push_back(std::pair<uint64_t, uint64_t>(
                levelOffsets[i], levelOffsets[i] + levelNumNodes[i]));
        return levelBounds;
    }
    uint64_t size() const { return _numNodes * sizeof(NodeItem); }
    static uint64_t size(const uint64_t numItems, const uint16_t nodeSize = 16)
    {
        if (nodeSize < 2)
            throw std::invalid_argument("Node size must be at least 2");
        if (numItems == 0)
            throw std::invalid_argument(
                "Number of items must be greater than 0");
        const uint16_t nodeSizeMin =
            std::min(std::max(nodeSize, static_cast<uint16_t>(2)),
                     static_cast<uint16_t>(65535));
        // limit so that resulting size in bytes can be represented by uint64_t
        if (numItems > static_cast<uint64_t>(1) << 56)
            throw std::overflow_error("Number of items must be less than 2^56");
        uint64_t n = numItems;
        uint64_t numNodes = n;
        do {
            n = (n + nodeSizeMin - 1) / nodeSizeMin;
            numNodes += n;
        } while (n != 1);
        return numNodes * sizeof(NodeItem);
    }

    void streamWrite(const std::function<void(uint8_t *, size_t)> &writeData)
    {
        writeData(reinterpret_cast<uint8_t *>(_nodeItems.data()),
                  static_cast<size_t>(_numNodes * sizeof(NodeItem)));
    }

    NodeItem getExtent() const { return _extent; }
};

} // namespace FlatGeobuf

#endif
