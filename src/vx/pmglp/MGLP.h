#ifndef MGLP_H
#define	MGLP_H

extern "C"{
    #include "../com/cc/superpixel/SNICT.h"
}

#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <map>
#include <typeinfo>
#include <cassert>


#include "../com/cp/image/Imagecpp.h"


#include "VectorMatrix.h"
#include "Graph.h"
#include "DinamicFeatures.h"
#include "Model.h"

class MGLP: public Model{
public:
    MGLP(){

    }

    ~MGLP(){

    }

    void setTheshold(double t){
        threshold = t;
    }

    void setTraversal(int t){
        traversal = t;
    }

    void update_labels(Image * img, Graph *gml ){
        gml->ITER++;
        int i,c=0;
        for(i=0; i<gml->SIZE; ++i){
            if (gml->V[gml->LABELS[i]]!=gml->ITER){
                gml->V[gml->LABELS[i]]=gml->ITER;
                gml->U[gml->LABELS[i]]=c;
                c++;
            }    
            gml->ORDER[i] = gml->U[gml->LABELS[i]];
        }
        //img->updateFromSP(gml->ORDER, c);
    }

    int sp_center(Image *img){
        int i, mx, my;
        mx=img->width()/2;
        my=img->height()/2;
        return img->label( (my)*img->width()+(mx) );        
    }
    
    int sp_init(Image *img){
        return img->label(0);
    }

    int highFreqLabel(Graph* gml, int olabel, int i){
        gml->ITER++;

        int j, nlabel=olabel;
        double maxcont = 0;
        double w, we;
        for (auto &e : gml->GRAPH[i]){
            j = gml->LABELS[e.first];
            w = e.second;
            if(w>=threshold){
                if(gml->V[j]!= gml->ITER){
                    gml->V[j] = gml->ITER;
                    gml->H[j] = 0.0;
                }
                gml->H[j] += 1.0;
                
                we = gml->H[j];
                if(we > maxcont){
                    maxcont = we;
                    nlabel = j;
                }
            }
        }
        return nlabel;
    }

    void ILP(Graph *glm, int iter){
        int olabel, nlabel, cc, deg, id, it, i;
        for(it=0; it<iter; ++it){
            cc = 0;
            for(i=0; i<glm->N; ++i){
                id = glm->ORDER[i];
                deg = glm->GRAPH[id].size();
                if(deg > 0){
                    olabel = glm->LABELS[id];
                    nlabel = highFreqLabel(glm, olabel, id);
                    glm->LABELS[id] = nlabel;
                    if(olabel != nlabel){
                        cc++;
                    }                        
                }
            }
            if(cc == 0){
                break;
            }
        }
    }

    void makeorder(int idsp, Graph *gml){   
        gml->ITER++;
        
        int i,j,co=0;
        i = gml->LABELS[idsp];

        std::deque<int> Q;
        Q.push_back(i);

        while (!Q.empty()){
            i = Q.front();
            Q.pop_front();
            if (gml->V[i] != gml->ITER){
                gml->V[i] = gml->ITER;
                gml->ORDER[co] = i;co++;

                for (auto &en : gml->GRAPH[i]){
                    j = en.first;
                    Q.push_back(j);
                }
            }
        }
        if(gml->N != co)
            std::cout<<">>>SPIRAL:"<<(gml->N)<<":"<<co<<std::endl;
        assert(gml->N==co);
    }

    void spiralx(
        int idsp,
        Graph *gml)
    {   
        gml->ITER++;
        int i,j,co=0;
        i = gml->LABELS[idsp];
        std::deque<int> Q;
        Q.push_back(i);
        
        gml->V[i] = gml->ITER;
        gml->ORDER[co] = i;co++;
        
        while (!Q.empty()){
            i = Q.front();
            Q.pop_front();
            for (auto &en : gml->GRAPH[i]){
                j = en.first;
                if (gml->V[j] != gml->ITER){
                    gml->V[j] = gml->ITER;
                    gml->ORDER[co] = j;co++;
//                    Q.push_front(jh);
                    Q.push_back(j);
                }
            }
        }
        assert(gml->N==co);
    }


    void execute(){
        Image* img = new Image(in_filename);

        //  super-pixels_t extraction;
        vx_superpixels_snic(img->imageio_t(), sp_side);
        vx_image_draw_regions(img->imageio_t());

        //img->write(ou_filename);
        VectorMatrix* vertex = Graph::makevertex(img);

        DinamicFeatures *dfe = new DinamicFeatures(vertex->rows(), vertex->cols());

        Graph *gml = new Graph(vertex->rows(), img);

        int i, idsp, aux;
        for (i=0; i<vertex->rows(); ++i){
            gml->LABELS[i] = i;
            gml->ORDER[i] = i;
            gml->CLUS[i].push(i);
        }

        // contruction graph with edges
        gml->makegraph(vertex);
        
        //convolution
        if (traversal==1){
            idsp = sp_init(img);
        }
        //espiral
        else if(traversal==2){
            idsp = sp_center(img);
        }
        
        int iteration = 0;
        bool nextlevel = true;
        while(nextlevel){
            makeorder(idsp, gml);
            aux = gml->N;
            ILP(gml, 1);
            gml->remakegraph(dfe,vertex);
            iteration++;
            nextlevel = (aux==gml->N) ? false : true; 
        }

        //update_labels(img, gml);
        
        delete gml;        
        delete dfe;
        delete vertex;
        delete img;
    }
};

#endif

