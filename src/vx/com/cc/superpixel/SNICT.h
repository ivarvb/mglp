#ifndef SNICT_H
#define	SNICT_H

#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 

#include "../image/Imagec.h"
#include "../algorithm/Priority.h"


void vx_superpixels_snic(Image_t * image, int side);


void vx_superpixels_snic(Image_t * imgt, int side){
        clock_t start = clock();
        clock_t end;
        //ImageIO_t* imgt = image->imageio_t();
        int N = imgt->size_pixels;
        int i, j, s, k, x, y, f, xi, yi, xf, yf;
        int w = imgt->width;
        int h = imgt->height;
        const int dn8[8] = {-1, -w, 1, w, -1-w,1-w,1+w,-1+w};

        //const double M = 5.0;//10.0;
        const double M = 10.0;//10.0;
//        const double M = 10.0;//10.0;
//        const double M = 2.0;//10.0;
        double superpixelSize = side*side;
        double innumk = (double)(N/(double)(superpixelSize));
        const double invwt = (M*M*innumk)/((double)(N));
        
        double d;
        double lf;
        double af;
        double bf;
        double colordist;
        double xydist;
                            
        for(i = 0; i < N; ++i){
            imgt->pixel[i].label=-1;
        }

        Heap *Q = new_heap(1,10);
        int MXK=0;
        for(x=0; x<w; x=x+side){
            for(y=0; y<h; y=y+side){
                i = y*w+x;
                HeapPair e = {0.0, i, MXK};
                heap_push(Q, e);
                MXK++;
            }
        }
        double *KL = (double*)calloc(MXK,sizeof(double));
        double *KA = (double*)calloc(MXK,sizeof(double));
        double *KB = (double*)calloc(MXK,sizeof(double));

        double *KX = (double*)calloc(MXK,sizeof(double));
        double *KY = (double*)calloc(MXK,sizeof(double));
        double *KZ = (double*)calloc(MXK,sizeof(double));

        //int *verif = (int*)calloc(MXK,sizeof(int));

        HeapPair node;
        int chan = imgt->channels;
        int z = imgt->size_pixels;
        //while (!Q.empty()) { 
        while (Q->count != 0) {
            node = heap_pop(Q);
            //Q.pop();
            i = node.i;
            k = node.k;
            
            if (imgt->pixel[i].label==-1){
                imgt->pixel[i].label = k;
                KL[k] += imgt->pixel_lab[0*z+i];
                KA[k] += imgt->pixel_lab[1*z+i];
                KB[k] += imgt->pixel_lab[2*z+i];
                KX[k] += imgt->pixel[i].x;
                KY[k] += imgt->pixel[i].y;
                KZ[k] += 1.0;
                //verif[k] = 1;

                for(f = 0; f < 8; ++f){
                    x = imgt->pixel[i].x+x8[f];
                    y = imgt->pixel[i].y+y8[f];
                    if(x>=0 && x<w && y>=0 && y<h){
                        j = y*w+x; 
                        //j = i+dn8[f];
                        if(imgt->pixel[j].label ==-1){
                            lf = KL[k] - imgt->pixel_lab[0*z+i]*KZ[k];
                            af = KA[k] - imgt->pixel_lab[1*z+i]*KZ[k];
                            bf = KB[k] - imgt->pixel_lab[2*z+i]*KZ[k];
                            xf = KX[k] - x*KZ[k];
                            yf = KY[k] - y*KZ[k];

                            colordist   = lf*lf + af*af + bf*bf;
                            xydist      = xf*xf + yf*yf;
                            d    = (colordist + xydist*invwt)/(KZ[k]*KZ[k]);
                            HeapPair ei = {d, (y*w+x), k};
                            heap_push(Q, ei);
                            //Q.push({j,d,k});                            
                            //printf("(%f, %d, %d)", d, j, k);
                            //printf("(%d)", Q->count);
                        }
                    }
                }
            }
        }
        end = clock();
        printf("time superpixels SNIC: %f\n", (double)(end - start)/CLOCKS_PER_SEC);

        free(KL); free(KB); free(KA);free(KX);free(KY);free(KZ);
        vx_image_update_regions(imgt);
        printf("SNIC LABELS  %d\n",imgt->size_labels);
        //imgt->size_labels = maxlabel;
        // std::cout<<"sp size:: "<<image->sizeSegments()<<" - "<<MXK<<std::endl;
        //free(verif);
        free_Heap(Q);
}
#endif