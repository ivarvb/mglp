#ifndef GRAPH_H
#define	GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <map>
#include <typeinfo>
#include <cassert>

#include "../com/cp/image/Imagecpp.h"

#include "VectorMatrix.h"
#include "DinamicFeatures.h"


typedef std::map<int, double> MapGraph;
typedef std::queue<int> ClusterNode;

class Graph{

public:
    int *LABELS;
    int *U;
    int *V;
    double *H;
    int *ORDER;
    int *ORDER2;
    MapGraph *GRAPH;
    MapGraph *GAUX;
    ClusterNode *CLUS;
    int SIZE;
    int N;
    int ITER;
    double threshold;

    Graph(int n, Image *img){
        LABELS = new int[n]();
        U = new int[n]();
        V = new int[n]();
        H = new double[n]();
        ORDER = new int[n]();
        ORDER2 = new int[n]();
        GRAPH = new MapGraph[n];
        GAUX = new MapGraph[n];
        CLUS = new ClusterNode[n];
        SIZE = n;
        N = n;
        ITER = 1;
        constructionfromimage(img);
    }  
    /**
     * @brief Destructor.
     */    
    ~Graph(){
        delete []LABELS;
        delete []U;
        delete []V;
        delete []H;
        delete []ORDER;
        delete []ORDER2;
        delete []GRAPH;
        delete []CLUS;
    }

    void constructionfromimage(Image *img){
        ITER++;
        int cc=0;
        int xp[]={1,1,0,-1,-1,-1, 0, 1};
        int yp[]={0,1,1, 1, 0,-1,-1,-1};

        int height = img->height();
        int width = img->width();
        int i,j,yj,xj,lbi,lbj,x,y,z;

        for (y=0; y<height-1; ++y) {
            for (x=1; x<width-1; ++x) {
                i = y*width+x;
                lbi = img->label(i);
                if (V[lbi]!=ITER){
                    V[lbi] = ITER;
                    ORDER[cc]=lbi;cc++;
                }
                for (z=0; z<5; ++z) {
                    yj = y+yp[z];
                    xj = x+xp[z];
                    j = yj*width+xj;
                    lbj = img->label(j);
                    if (lbi!=lbj){
                        GAUX[lbi][lbj] = 1.0;
                        GAUX[lbj][lbi] = 1.0;
                    }
                    if (V[lbj]!=ITER){
                        V[lbj] = ITER;
                        ORDER[cc]=lbj;cc++;
                    }
                }
            }
        }
        N = img->size_labels();
    }
    
    double proximity(VectorMatrix *vm, int di, int dj){
        double d = 0.0;
        double a = 0.0;
        int i;
        for(i=0; i<vm->cols(); ++i){            
            a = vm->at(di, i) - vm->at(di, i);
            d += a*a;
        }
        d = std::sqrt(d);
        return exp( -1.0*( (d*d)/(2.0*(0.5*0.5)) ) );
    }

    //make edges
    void makegraph(VectorMatrix *vertex){
        double df = -1.0;
        int i,j,li,lj;

        for (i=0;i<N;++i){            
            li = ORDER[i];
            GRAPH[li].clear();           
        }

        for (i=0;i<N;++i){
            li = ORDER[i];
            for (auto &e : GAUX[li]){
                lj = e.first;
                if(li<lj){
                    df = proximity(vertex, li, lj);
                    //df = 0.0;
                    // if(df <= threshold){
                    //     df = 0.0;
                    // }                            

                    GRAPH[li][lj] = df;
                    GRAPH[lj][li] = df;
                }
            }
        }
    }

    void remakegraph(DinamicFeatures* dfe, VectorMatrix* vertex){        
        ITER++;
        int i, j, id, li, lj, cc = 0;

        // nuevo orden de labels
        i = ORDER[0];
        
        for (i=0;i<N;++i){
            id = ORDER[i];
            ORDER2[i] = id;
            li = LABELS[id];
            GAUX[li].clear();
        }

        for (id=0;id<N;++id){
            i = ORDER2[id];
            li = LABELS[i];
            for (auto &ne : GRAPH[i]){
                j = ne.first;
                lj = LABELS[j];
                //join vertex percorrer o grafo atual
                if (U[j] != ITER){
                    U[j]  = ITER;
//                    Q.push_back(j);
                    dfe->add(lj, vertex, j);
                    if(lj!=j){
                        mergecluster(lj, j);
                    }
                }
                //to new order of labels
                if (V[lj] != ITER){
                    V[lj]  = ITER;
                    ORDER[cc]=lj;cc++;
                }                
                ////////build new edges
                if (li!=lj){
                    GAUX[li][lj] = 1;
                    GAUX[lj][li] = 1;
                }
            }
        }      
        N = cc;

        dfe->transfer(vertex, N, ORDER);
        makegraph(vertex);
    }

    void mergecluster(int labeli, int labelj){
        int e;
        ClusterNode &aux1 = CLUS[labeli];
        ClusterNode &aux2 = CLUS[labelj];
        while (!aux2.empty()){
            e = aux2.front();
            aux2.pop();            
            aux1.push(e);
            LABELS[e] = labeli;
        }
    }

    static VectorMatrix* makevertex(Image* image){
        int i, sizeS = image->size_labels();
        //std::cout<<"sizeS"<<sizeS<<std::endl;
        VectorMatrix* features = new VectorMatrix(sizeS, 9);
        vx_image_norm_zscore(image->imageio_t());
        vx_image_to_colormoments_features(image->imageio_t(), features->data());
        return features;
    }
};
#endif