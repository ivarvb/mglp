#ifndef VECTORMATRIX_H
#define	VECTORMATRIX_H

#include <iostream>

class VectorMatrix{
private:                
    double* _data;
    int _rows;
    int _cols;
public:

    /**
     * @brief Constructor.
     * */
    VectorMatrix(int rows, int cols){
        _rows = rows;
        _cols = cols;

        _data = new double[_rows*_cols];
        int i;
        for(i=0; i<_rows*_cols; ++i){
            _data[i] = 0.0;
        }
    }

    /**
     * @brief Destructor.
     */    
    ~VectorMatrix(){
        delete []_data;
    }
    
    int rows(){
        return _rows;
    }
    
    int cols(){
        return _cols;
    }

    double* data(){
        return _data;
    }

    inline double& at(int row, int col){
        return _data[row*(_cols)+col];
    }
};
#endif