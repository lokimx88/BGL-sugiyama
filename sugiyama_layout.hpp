// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// this project is translated from this project : digraph [https://github.com/damian0815/digraph]

//  Author:  MingXing Liu
//  Date:    2017-12-25
//  version:
//  v0.1
//  1.init verion

#ifndef BOOST_GRAPH_SUGIYAMA_LAYOUT_HPP
#define BOOST_GRAPH_SUGIYAMA_LAYOUT_HPP

#include <math.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/simple_point.hpp>
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/static_assert.hpp>
#include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>

namespace boost {

static const int    MAX_SWEEPS = 100;

template<typename GraphType>
static void splitIntoLayers(const GraphType& g,
                            std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                            std::map<typename graph_traits<GraphType>::vertex_iterator, int> &nodemap) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;
    typedef typename graph_traits<GraphType>::vertex_descriptor VertexDescriptor;
    typedef typename graph_traits<GraphType>::out_edge_iterator OutEdgeIter;

    typedef std::deque< VertexDescriptor > V_VEC;
    V_VEC sorted;
    topological_sort(g, std::front_inserter(sorted));

    std::map<VertexDescriptor, int> lmap;
    for (typename V_VEC::value_type &vt:sorted){
        lmap[vt] = 0;
    }

    int h = 1;
    for (typename V_VEC::iterator it = sorted.begin(); it != sorted.end(); ++it )
    {
        VertexDescriptor n1 = *it;

        OutEdgeIter out, out_end;
        tie(out, out_end) = out_edges(n1, g);
        for(; out != out_end; ++out){
            VertexDescriptor n2=target(*out,g);

            int inc = 1;
            lmap[n2] = std::max(lmap[n1] + inc, lmap[n2]);
            h = std::max(h, lmap[n2] + 1);
        }
    }

    //stack_init(h, (int)sorted.size());
    layers.resize(h);

    Viter n, nend;
    tie(n,nend) = vertices(g);
    for (typename V_VEC::value_type &vt:sorted){
        //stack_add(n+vt, lmap[vt]);
        int layerIndex=lmap[vt];

        layers[layerIndex].push_back(n+vt);
        nodemap[n+vt] = layerIndex;
    }
}

template<typename GraphType,typename IndexMap>
static void setOrderedIndexes(IndexMap &index,
                       std::vector<typename graph_traits<GraphType>::vertex_iterator>& ln) {
    for (int i = 0; i < ln.size(); i++) {
        index[*(ln[i])]=i;
    }
}

template<typename GraphType,typename RectMap>
static double maxHeight(RectMap &rect,
                 const std::vector<typename graph_traits<GraphType>::vertex_iterator> &ln) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    double mh = 0;
    for (typename std::vector<Viter>::const_iterator it = ln.begin(); it != ln.end(); ++it ) {
        Viter n = (*it);
        mh = std::max(mh, (double)(rect[*n][1]));
    }
    return mh;
}

template<typename GraphType,typename PositionMap,typename RectMap>
static double avgX(PositionMap &position,
            RectMap &rect,
            const std::vector<typename graph_traits<GraphType>::vertex_iterator>& ln) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    double m = 0;
    for (typename std::vector<Viter>::const_iterator it = ln.begin(); it != ln.end(); ++it ) {
        Viter n = *it;
        m += (position[*n][0]+rect[*n][0]/2);
    }
    return m / ln.size();
}

template<typename GraphType,typename IndexMap>
static int barycenter( IndexMap &index,
                const std::vector<typename graph_traits<GraphType>::vertex_iterator>& ln) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    if (ln.size() == 0) {
        return 0;
    }

    double bc = 0;
    for (typename std::vector<Viter>::const_iterator it = ln.begin(); it != ln.end(); ++it ) {
        Viter n = *it;

        bc += index[*n];
    }
    return (int) floor((bc / ln.size())+0.5f);
}

template<typename GraphType>
std::vector<typename graph_traits<GraphType>::vertex_iterator>
static stack_getConnectedTo(const GraphType& g,
                     std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                     typename graph_traits<GraphType>::vertex_iterator n1,
                     int layerIndex) {
    typedef typename graph_traits<GraphType>::vertex_iterator        Viter;
    typedef typename boost::graph_traits<GraphType>::edge_descriptor EdgeDescriptor;

    std::vector<Viter> ln;
    if (layerIndex < layers.size() && layerIndex >= 0) {
        for (typename std::vector<Viter>::iterator it = layers.at(layerIndex).begin(); it != layers.at(layerIndex).end(); ++it )
        {
            Viter n2 = *it;
            std::pair<EdgeDescriptor, bool> test1 = edge(*n1, *n2, g);
            if (test1.second == true) { // 判断vd1和vd2所指的两节点是否相邻
                ln.push_back(n2);
            }

            /*std::pair<EdgeDescriptor, bool> test2 = edge( *n2,*n1, g);
              if (test2.second == true) { // 判断vd1和vd2所指的两节点是否相邻
                 ln.push_back(n2);
             }*/
        }
    }
    return ln;
}

template<typename GraphType,typename IndexMap>
static void stack_reduceCrossings2L(const GraphType& g,
                             std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                             IndexMap &index,
                             int staticIndex, int flexIndex) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    assert( flexIndex < layers.size() );
    std::vector<Viter>& flex = layers[flexIndex];
    for (typename std::vector<Viter>::const_iterator it = flex.begin(); it != flex.end(); ++it )
    {
        Viter n = (*it);
        std::vector<Viter> neighbors = stack_getConnectedTo(g,layers,n, staticIndex);

        index[*n]=barycenter<GraphType>(index,neighbors);
    }

    sort( flex.begin(), flex.end(),  [&index]( Viter& n1,  Viter& n2 )->bool{
        return index[*n1]< index[*n2];
    });

    setOrderedIndexes<GraphType>(index,flex);
}

template<typename GraphType,typename PositionMap,typename RectMap>
static void stack_layerHeights(std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                        PositionMap &position,
                        RectMap &rect,
                        int &yspacing) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    double offset = 0;
    for (int l = 0; l < layers.size(); l++) {
        std::vector<Viter>& ln = layers[l];
        double maxh = maxHeight<GraphType>(rect,ln);
        for (typename std::vector<Viter>::iterator it = ln.begin(); it != ln.end(); ++it ) {
            Viter n = (*it);
            //if (n->isVirtual()) {
            if (false) {
                position[*n][1]=offset+maxHeight<GraphType>(rect,ln)/2.0;
            } else {
                position[*n][1]=offset;
            }
        }
        offset += maxh + yspacing;
    }
}

template<typename GraphType, typename PositionMap,typename RectMap>
static void stack_xPosPack(std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                    PositionMap &position,
                    RectMap &rect,
                    int flexIndex,
                    int &xspacing) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    assert( flexIndex < layers.size() );
    std::vector<Viter>& flex = layers[flexIndex];
    double offset = 0;
    for (typename std::vector<Viter>::const_iterator it = flex.begin(); it != flex.end(); ++it ) {
        Viter n = (*it);
        position[*n][0]=offset;
        offset = position[*n][0] + rect[*n][0] + xspacing;
    }
}

template<typename GraphType, typename PositionMap,typename RectMap>
static  void stack_xPosDown(const GraphType& g,
                    PositionMap &position,
                    RectMap &rect,
                    std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                    int staticIndex, int flexIndex,
                    int &xspacing) {
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    assert( flexIndex < layers.size() );
    std::vector<Viter> flex = layers[flexIndex];
    for (int i = 0; i < flex.size(); i++) {
        Viter n = flex[i];
        std::vector<Viter> neighbors = stack_getConnectedTo(g,layers,n, staticIndex);
        double avg = avgX<GraphType>(position,rect,neighbors);
        double min = (i > 0) ? (position[*(flex[i-1])][0] + rect[*(flex[i-1])][0] + xspacing) : (-std::numeric_limits<double>::max());
        if (!isnan(avg)) {
            position[*n][0]= std::max(min, avg -  rect[*n][0]/2.0);
        }
    }
}

template<typename GraphType, typename PositionMap,typename RectMap>
static  void stack_xPosUp(const GraphType& g,
                  PositionMap &position,
                  RectMap &rect,
                  std::vector<std::vector<typename graph_traits<GraphType>::vertex_iterator>> &layers,
                  int staticIndex, int flexIndex,
                  int &xspacing) {//staticIndex=flexIndex+1
    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    assert( flexIndex < layers.size() );
    std::vector<Viter> flex = layers[flexIndex];
    for (int i = flex.size() - 1; i >= 0; i--) {
        Viter n = flex[i];
        std::vector<Viter> neighbors = stack_getConnectedTo(g,layers,n, staticIndex);

        //for ( int i=0; i<neighbors.size() ;i++ )
        double avg = avgX<GraphType>(position,rect,neighbors);
        // calculate min, max
        double min = (i>0)?(position[*(flex[i-1])][0] + rect[*(flex[i-1])][0] + xspacing) : -std::numeric_limits<double>::max();
        double max =std::numeric_limits<double>::max();
        if(i<flex.size()-1){
            max=position[*(flex[i+1])][0] - rect[*n][0] - xspacing;
        }

        if (!isnan(avg)) {
            position[*n][0]=std::max(min, std::min(max, avg - rect[*n][0] / 2.0));
        }
    }
}


template<typename GraphType, typename PositionMap,typename RectMap,typename IndexMap>
void sugiyama_graph_layout(const GraphType& g,
                           PositionMap &position,
                           RectMap &rect,
                           IndexMap &index,
                           int xspacing=20,
                           int yspacing=30)
{
    BOOST_STATIC_ASSERT (property_traits<PositionMap>::value_type::dimensions >= 2);
    BOOST_STATIC_ASSERT (property_traits<RectMap>::value_type::dimensions >= 2);

    typedef typename graph_traits<GraphType>::vertex_iterator   Viter;

    typedef  std::vector<Viter>     Viter_VEC;
    typedef  std::vector<Viter_VEC> Viter_VEC_VEC;
    Viter_VEC_VEC        layers;
    std::map<Viter, int> nodemap;

    //removeCycles
    //removeCycles();

    //splitIntoLayers
    splitIntoLayers(g,layers,nodemap);

    //insertDummies
    //insertDummies();

    //initIndexes
    for (typename Viter_VEC_VEC::iterator it = layers.begin(); it != layers.end(); ++it )
    {
        Viter_VEC& l = *it;

        /*sort( l.begin(), l.end(), [&getXYWH]( Viter n1, Viter n2 )->bool{
              return getXYWH[*n1]->id < getXYWH[*n2]->id;
          });*/
        /*
          @todo figure out how to do this
          Collections.sort(l, new Comparator<Viter>() {
              public int compare(Viter n1, Viter n2) {
                  return n1.getIndex() - n2.getIndex();
              }
          });*/

        setOrderedIndexes<GraphType>(index,l);
    }

    //reduceCrossings
    for (int round = 0; round < MAX_SWEEPS; round++) {
        if (round % 2 == 0) {
            for (int l = 0; l < layers.size() - 1; l++) {
                stack_reduceCrossings2L(g,layers,index,l, l + 1);
            }
        } else {
            for (int l = layers.size() - 1; l > 0; l--) {
                stack_reduceCrossings2L(g,layers,index,l, l - 1);
            }
        }
    }

    //undoRemoveCycles
    //undoRemoveCycles();

    //layerHeights
    stack_layerHeights<GraphType>(layers,position,rect,yspacing);

    //xPos
    for (int l = 0; l < layers.size(); l++) {
        stack_xPosPack<GraphType>(layers,position,rect,l,xspacing);
    }

    printf("2\n");
    for (int l = 0; l < layers.size() - 1; l++) {
        stack_xPosDown(g,position,rect,layers,l, l + 1,xspacing);
    }
    for (int l = layers.size() - 1; l >0; l--) {
        stack_xPosUp(g,position,rect,layers,l, l - 1,xspacing);
    }
}


}

#endif
