#ifndef DINAMICFEATURES_H
#define	DINAMICFEATURES_H

#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <map>
#include <typeinfo>
#include <cassert>

typedef std::map<int, double*> MapFeatures;

class DinamicFeatures{
public:
    MapFeatures *mdata;
    MapFeatures::iterator it;
    int *size;
    int cols;
    int rows;
    DinamicFeatures(int rows, int cols){
        mdata = new MapFeatures();
        size = new int[rows]();
        this->cols = cols;
        this->rows = rows;
    }
    ~DinamicFeatures(){
        for(auto & e : *mdata){
            delete []e.second;
        }
        delete mdata;
        delete[]size;
    }

    void insert(int k){
        it = mdata->find(k);
        if (it == mdata->end()){
            double *d = new double[cols]();
            int i = 0;
            //initialization with zero
            for(i=0;i<cols;++i){
                d[i] = 0.0;
            }
            (*mdata)[k] = d;
            size[k] = 0;
        }
    }

    void add(int k, VectorMatrix *vm, int j){
        insert(k);
        int i;
        double *aux = (*mdata)[k];
        for(i=0; i<cols; ++i){
            aux[i] += vm->at(j, i);
        }
        size[k] += 1;
    }

    void transfer(VectorMatrix *vm, int N, int* ORDER){
        int i,c,k;
        double *aux;
        for(i=0; i<N; ++i){
            k = ORDER[i];
            aux = (*mdata)[k];
            for(c=0; c<cols; ++c){
                aux[c] /= size[k];
                vm->at(k, c) = aux[c];
                aux[c] = 0.0;
            }
            size[k] = 0;
        }
    }
};

#endif
