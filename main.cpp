#include "sugiyama_layout.hpp"

enum vertex_xy_t { vertex_xy };
enum vertex_wh_t { vertex_wh };

namespace boost {
    BOOST_INSTALL_PROPERTY(vertex, xy);
    BOOST_INSTALL_PROPERTY(vertex, wh);
}

int main(){
    typedef boost::rectangle_topology<>::point_type            Point;
    typedef boost::rectangle_topology<>::point_difference_type Rect;

    //x0 y1 w0 h1
    typedef boost::property<boost::vertex_index1_t, int,
                boost::property<boost::vertex_index2_t, int,
                    boost::property<vertex_xy_t, Point,
                        boost::property<vertex_wh_t, Rect> > > > V_PROPERTY;

    typedef boost::adjacency_list<boost::vecS,
                                    boost::vecS,
                                    boost::bidirectionalS,
                                    V_PROPERTY> Graph;

    typedef  boost::property_map<Graph, boost::vertex_index1_t>::type Index1Map;
    typedef  boost::property_map<Graph, boost::vertex_index2_t>::type Index2Map;
    typedef  boost::property_map<Graph, vertex_xy_t>::type            PointMap;
    typedef  boost::property_map<Graph, vertex_wh_t>::type            RectMap;

    typedef  boost::graph_traits<Graph>::vertex_descriptor  VertexDescriptor;
    typedef  boost::graph_traits<Graph>::vertex_iterator    Viter;

    Graph g;

    Index1Map vidxmap =boost::get(boost::vertex_index1, g);
    Index2Map vidx2map=boost::get(boost::vertex_index2, g);
    PointMap  ptmap   =boost::get(vertex_xy, g);
    RectMap   rtmap   =boost::get(vertex_wh, g);


    VertexDescriptor n6=boost::add_vertex(g);
    vidxmap[n6]=6;
    VertexDescriptor n5=boost::add_vertex(g);
    vidxmap[n5]=5;
    VertexDescriptor n4=boost::add_vertex(g);
    vidxmap[n4]=4;
    VertexDescriptor n1=boost::add_vertex(g);
    vidxmap[n1]=1;
    vidx2map[n1]=0;
    ptmap[n1][0]=0;
    ptmap[n1][1]=0;
    rtmap[n1][0]=50.0;
    rtmap[n1][1]=20.0;
    VertexDescriptor n2=boost::add_vertex(g);
    vidxmap[n2]=2;
    VertexDescriptor n3=boost::add_vertex(g);
    vidxmap[n3]=3;
    VertexDescriptor n7=boost::add_vertex(g);
    vidxmap[n7]=7;
    VertexDescriptor n8=boost::add_vertex(g);
    vidxmap[n8]=8;


    Viter ui, uiend;
    boost::tie(ui, uiend) = boost::vertices(g);
    for (; ui != uiend; ++ui) {
        vidx2map[*ui]=0;
        ptmap[*ui][0]=0.0;
        ptmap[*ui][1]=0.0;
        rtmap[*ui][0]=50.0;
        rtmap[*ui][1]=20.0;
    }

    boost::add_edge(n3, n6, g);
    boost::add_edge(n2, n5, g);
    boost::add_edge(n2, n4, g);
    boost::add_edge(n1, n2, g);
    boost::add_edge(n1, n3, g);

    boost::sugiyama_graph_layout(g,ptmap,rtmap,vidx2map);

    boost::tie(ui, uiend) = boost::vertices(g);
    for (; ui != uiend; ++ui) {
        printf("%d %d %f %f \n",vidxmap[*ui],vidx2map[*ui],
                ptmap[*ui][0],ptmap[*ui][1]) ;
    }
}
