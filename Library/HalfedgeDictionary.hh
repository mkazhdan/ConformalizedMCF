////////////////////////////////////////////////////////////////////////////////
// HalfedgeDictionary.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Implements a custom hashmap-like data structure for speeding up
//      halfedge queries from (start vertex, end vertex) pairs.
//
//      Assumes that vertices don't have high incident edge count:
//      Lookups are O(max_{v \in vertices} incidentCount(v))
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//
//  Created:  07/19/2011 16:35:46
//  Revision History:
//      07/19/2011  Julian Panetta    Initial Revision
////////////////////////////////////////////////////////////////////////////////
#ifndef HALFEDGE_DICTIONARY_HH
#define HALFEDGE_DICTIONARY_HH
#include <list>
#include <vector>
#include <cassert>

class HalfedgeDictionary
{
public:
    /** Returned for failed queries */
    enum {INVALID_HALFEDGE = -1, INVALID_RECORD = -1, INVALID_VERTEX = -1};
private:
    // NOTE: IF THIS IS CHANGED, THE getBucket FUNCTION MUST BE UPDATED
    static const unsigned int BUCKETS_PER_VERTEX = 8;

    struct EdgeRecord {
        unsigned int endVertex;
        unsigned int halfedge;
        unsigned int nextRecord;

        EdgeRecord()
            : endVertex(INVALID_VERTEX), halfedge(INVALID_HALFEDGE),
              nextRecord(INVALID_RECORD)
        { }

        EdgeRecord(unsigned int end, unsigned int halfedge)
            : endVertex(end), halfedge(halfedge),
              nextRecord(INVALID_RECORD)
        { }
    };

    typedef unsigned int    EdgeBucket;
    typedef unsigned int    EdgeRecordIdx;

    // Indices of the first EdgeRecord
    std::vector<unsigned int> m_halfedgeBuckets;
    std::vector<EdgeRecord>   m_halfedges;

    EdgeBucket &getBucket(unsigned int start, unsigned int end)
    {
        return m_halfedgeBuckets[(start << 3) + (end & 7)];
    }

    EdgeBucket getBucket(unsigned int start, unsigned int end) const
    {
        return m_halfedgeBuckets[(start << 3) + (end & 7)];
    }

public:
    ////////////////////////////////////////////////////////////////////////////
    /*! Constructor allocates all space ever needed by HE structure.
    //  @param[in]  numVertices     the number of vertices
    //  @param[in]  numHalfedges    the number of halfedges
    *///////////////////////////////////////////////////////////////////////////
    HalfedgeDictionary(unsigned int numVertices, unsigned int numHalfedges)
    {
        m_halfedgeBuckets.resize(BUCKETS_PER_VERTEX * numVertices,
                                 INVALID_RECORD);
        m_halfedges.reserve(numHalfedges);
    }

    ////////////////////////////////////////////////////////////////////////////
    /*! Inserts a halfedge (start --> end)
    //  @param[in]  start       halfedge start vertex index
    //  @param[in]  end         halfedge end vertex index
    //  @param[in]  halfedge    halfedge index
    *///////////////////////////////////////////////////////////////////////////
    void insert(unsigned int start, unsigned int end, unsigned int halfedge)
    {
        EdgeBucket &bucket = getBucket(start, end);
        EdgeRecordIdx curr = bucket;
        EdgeRecordIdx prev = INVALID_RECORD;
        while (curr != INVALID_RECORD && (m_halfedges[curr].endVertex != end)) {
            prev = curr;
            curr = m_halfedges[curr].nextRecord;
        };

        if (curr == INVALID_RECORD) {
            // Insert a new halfedge, not exceeding the known HE count
            curr = HalfedgeDictionary::EdgeRecordIdx( m_halfedges.size() );
            assert(curr < m_halfedges.capacity());
            m_halfedges.push_back(EdgeRecord(end, halfedge));

            if (prev != INVALID_RECORD) {
                // Insert at end of list
                assert(m_halfedges[prev].nextRecord == INVALID_RECORD);
                m_halfedges[prev].nextRecord = curr;
            }
            else    {
                // Insert into empty list
                bucket = curr;
            }
        }
        else    {
            // Update an existing record
            EdgeRecord &rec = m_halfedges[curr];
            if (rec.halfedge != halfedge)    {
                printf("WARNING: changing index for halfedge %i --> %i to %i\n",
                        start, end, halfedge);
                rec.halfedge = halfedge;
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    /*! Retrieves the index for the halfedge (start --> end)
    //  @param[in]  start       halfedge start vertex index
    //  @param[in]  end         halfedge end vertex index
    //  @return     index of (start --> end), or
    //              INVALID_HALFEDGE if none is found
    *///////////////////////////////////////////////////////////////////////////
    unsigned int halfedge(unsigned int start, unsigned int end) const
    {
        EdgeRecordIdx curr = getBucket(start, end);
        while (curr != INVALID_RECORD) {
            if (m_halfedges[curr].endVertex == end)
                break;
            curr = m_halfedges[curr].nextRecord;
        };

        return (curr != INVALID_RECORD) ? m_halfedges[curr].halfedge
                                        : INVALID_HALFEDGE;
    }
};

#endif
