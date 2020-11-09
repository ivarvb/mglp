#ifndef ASNICI_H
#define	ASNICI_H

//#include <mex.h>
#include <cmath>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <queue>

#include "../image/Imagecpp.h"

using namespace std;

class ElementOrderi{
public:
    int i;
    double d;
    int k;
    ElementOrderi()
    {
        i = 0;
        d = 0.0;
        k = 0.0;
    }

    ElementOrderi(int _i, double _d, int _k)
    {
        i = _i;
        d = _d;
        k = _k;
    }

    bool operator<(const ElementOrderi& obj) const
    {
//        //ordena de mayor hasta menor
//        return d < obj.d;
        //ordena de menor a maior 
        return d > obj.d;
    }
};





class SNICI{
public:
   
    double proximity(double cl, double ca, double cb, double l, double a, double b){
        double d;
        d = (   ((cl-l)*(cl-l))+
                ((ca-a)*(ca-a))+
                ((cb-b)*(cb-b)) 
            );
//        ai = cl*ca*cb;
//        bi = a*b*c;
//        d = std::sqrt(d/20.98);
//        d = std::sqrt(d);
//        r = (((cl-l)*(cl-l))+((ca-a)*(ca-a))+((cb-b)*(cb-b)));
        
//        d = exp((-1.0*(d*d))/(2.0*(10.75*10.75)));
//        r = d;
//        return 1.0-d;
        return d;
    }
    void execute(ImageIO * image, int slide){
        double _time;
        RunTime *timeru = new RunTime;   
        timeru->start();
        ImageIO_t* imgt = image->imageio_t();
        int size = imgt->size_pixels;
        int i, j, s, k, x, y, f, xi, yi, xf, yf;
        int n = image->size();
        int w = image->width();
        int h = image->height();

        const double M = 5.0;//10.0;
//        const double M = 10.0;//10.0;
//        const double M = 2.0;//10.0;
        double superpixelSize = slide*slide;
        double innumk = double(image->size())/double(superpixelSize);
        const double invwt = (M*M*innumk)/double(image->size());
        
        double d;
        double lf;
        double af;
        double bf;
        double colordist;
        double xydist;
                            
//        int *TARGETS = new int[n]{};
//        std::cout<<"11"<<std::endl;
        for(i = 0; i < image->size(); ++i){
            imgt->targets[i]=-1;
            imgt->pixel[i].size_region = 0;
        }
        
        std::priority_queue<ElementOrderi> Q;
        std::vector<double>KL;
        std::vector<double>KA;
        std::vector<double>KB;
        std::vector<double>KX;
        std::vector<double>KY;
        std::vector<double>KZ;

//        std::cout<<"22"<<std::endl;
        int MXK=0;
        for(x=0; x<w; x=x+slide){
            for(y=0; y<h; y=y+slide){
                i = y*w+x;
                Q.push({i, 0, MXK});
                KL.push_back(0);
                KA.push_back(0);
                KB.push_back(0);
                KZ.push_back(0);
                KX.push_back(0);
                KY.push_back(0);
                MXK++;
            }
        }
        int *verif = new int[MXK];
        ElementOrderi node;
        while (!Q.empty()) { 
            node = Q.top();
            Q.pop();
            i = node.i;
            k = node.k;
            if (imgt->targets[i]==-1){
                imgt->targets[i] = k;
                KL[k] += imgt->pixel_color[0*size+i];
                KA[k] += imgt->pixel_color[1*size+i];
                KB[k] += imgt->pixel_color[2*size+i];
                KX[k] += imgt->pixel[i].x;
                KY[k] += imgt->pixel[i].y;
                KZ[k] += 1.0;
                verif[k] = 1;

                for(f = 0; f < 8; ++f){
                    x = image->x(i)+x8[f];
                    y = image->y(i)+y8[f];
                    if(x>=0 && x<image->width() && y>=0 && y<image->height()){
                        j = y*w+x;
                        if(imgt->targets[j]==-1){
                            lf = KL[k] - imgt->pixel_color[0*size+j]*KZ[k];
                            af = KA[k] - imgt->pixel_color[1*size+j]*KZ[k];
                            bf = KB[k] - imgt->pixel_color[2*size+j]*KZ[k];
                            xf = KX[k] - imgt->pixel[j].x*KZ[k];
                            yf = KY[k] - imgt->pixel[j].y*KZ[k];

                            colordist   = lf*lf + af*af + bf*bf;
                            xydist      = xf*xf + yf*yf;
                            d    = (colordist + xydist*invwt)/(KZ[k]*KZ[k]);//late normalization by ksize[k], to have only one division operation
//                            d    = (colordist)/(KZ[k]*KZ[k]);//late normalization by ksize[k], to have only one division operation
                            
//                            d = proximity(KL[k], KA[k], KB[k], lf, af, bf);
//                            d /= 0.005;
                            Q.push({j,d,k});                            
                        }
                    }
                }
            }
        }
//        if(labels[0] < 0) labels[0] = 0;
//        for(int y = 1; y < height; y++)
//        {
//            for(int x = 1; x < width; x++)
//            {
//                int i = y*width+x;
//                if(labels[i] < 0)//find an adjacent label
//                {
//                    if(labels[i-1] >= 0) labels[i] = labels[i-1];
//                    else if(labels[i-width] >= 0) labels[i] = labels[i-width];
//                }//if labels[i] < 0 ends
//            }
//        }

        timeru->end();
        _time = timeru->getRunTime();
        timeru->printRunTime("Alg. SNICI");
        delete timeru;

        KL.clear();
        KB.clear();
        KA.clear();
        KX.clear();
        KY.clear();
        KZ.clear();
        
//        Color3RGB *colors = makecolor(MXK);
//        int *verx = new int[MXK];
        
        int t, maxlabel = 0;
        for(i = 0; i < MXK; ++i){
            if (verif[i] == 1){
                verif[i] = maxlabel;
                maxlabel++;
            }
        }
        for(i = 0; i < image->size(); ++i){
            t = verif[imgt->targets[i]];
            imgt->targets[i] = t;
            
            imgt->pixel[i].next = &imgt->pixel[imgt->pixel[t].region];
            imgt->pixel[t].region = imgt->pixel[i].i;
            imgt->pixel[t].size_region += 1;

        }
        imgt->size_targets = maxlabel;
        
        std::cout<<"sp size:: "<<image->sizeSegments()<<" - "<<MXK<<std::endl;
        delete []verif;

    }

};
#endif