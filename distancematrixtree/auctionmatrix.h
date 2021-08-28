//
//  auctionmatrix.h
//  iqtree
//
//  Created by James Barbetti on 14/11/20.
//

#ifndef auctionmatrix_h
#define auctionmatrix_h

#include "rapidnj.h"

namespace StartTree {
    template <class T=NJFloat, class M=NJMatrix<T>>
    class AuctionMatrix: public BoundingMatrix<T, M>  {
        //
        //A reverse auction algorithm: each row "bids" to be clustered
        //(lowest bid wins) based on the difference between the best join
        //involving it (and another lower-numbered cluster) and the
        //second-best.  if there aren't two lower-numbered clusters,
        //no bid is placed (infiniteDistance is returned).
        //
        //This is a little faster than actually looking at all of the
        //distances!  And it practice, yields answers that aren't that
        //much worse.
        //
    public:
        typedef StartTree::BoundingMatrix<T, M> super;
        using super::clusterToRow;
        using super::entryToCluster;
        using super::entriesSorted;
        using super::purgeRow;
        virtual std::string getAlgorithmName() const {
            return "Auction" + M::getAlgorithmName();
        }
        StartTree::Position<T> getRowMinimum(size_t row, T maxTot, T qBest) const {
            StartTree::Position<T> pos(row, 0, infiniteDistance, 0);
            const T*   rowData   = entriesSorted.rows[row];
            const int* toCluster = entryToCluster.rows[row];
            size_t i = 0;
            int   bestOtherRow;
            for (i=0; rowData[i]<infiniteDistance; ++i) {
                bestOtherRow = clusterToRow[toCluster[i]];
                if (bestOtherRow != notMappedToRow) {
                    break;
                }
            }
            if (rowData[i] < infiniteDistance) {
                T     bestV = rowData[i];
                int   secondBestOtherRow;
                for (; rowData[i]<infiniteDistance; ++i) {
                    secondBestOtherRow = clusterToRow[toCluster[i]];
                    if (secondBestOtherRow != notMappedToRow) {
                        break;
                    }
                }
                if (rowData[i] < infiniteDistance) {
                    pos.row   = bestOtherRow;
                    pos.value = bestV - rowData[secondBestOtherRow];
                }
                if (1<i) {
                    purgeRow(row);
                }
            }
            return pos;
        }
    };

    typedef AuctionMatrix<NJFloat, NJMatrix<NJFloat>> DistanceAuctionMatrix;
}

#endif /* auctionmatrix_h */
