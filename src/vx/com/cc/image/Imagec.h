/** 
/* Copyright (C) 2012-2020 Ivar Vargas Belizario
/*
/*
/*
 */

#ifndef IMAGEC_H
#define IMAGEC_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "Color3RGB.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb_image/stb_image_write.h"

const int x8[8] = {-1, 0, 1, 0, -1, 1, 1, -1}; //for 4 or 8 connectivity
const int y8[8] = {0, -1, 0, 1, -1, -1, 1, 1}; //for 4 or 8 connectivity

// declare types and methods
typedef struct Region_t Region_t;
typedef struct Pixel_t Pixel_t;
typedef struct Image_t Image_t;

Image_t* vx_image_read(const char filename[]);
Image_t *vx_image_create(int width, int height, int init_intensity, unsigned char *fdata, int fdatai);
void vx_image_setcolor(Image_t *img, int index, Color3RGB color);
Color3RGB vx_image_getcolor(Image_t *img, int index);
void vx_image_copy(unsigned char *datafrom, unsigned char *datato, int size);
void vx_image_write(Image_t* img, const char filename[]);
void vx_image_draw_regions(Image_t* img);
void vx_image_free(Image_t *img);
//void Image_t_create_lab_color(Image_t *image);
void vx_image_update_regions(Image_t *image);

void rgb2lab(
    unsigned char *sR,
    unsigned char *sG,
    unsigned char *sB,
    double *lval,
    double *aval,
    double *bval);
void rgb2xyz(
    unsigned char *sR,
    unsigned char *sG,
    unsigned char *sB,
    double *X,
    double *Y,
    double *Z);

void RGB2LABT(  unsigned char *R_v,
                unsigned char *G_v,
                unsigned char *B_v,
                double *Lv, double *Bv, double *Av);
// implementing
typedef struct Region_t
{
    Pixel_t *next;
    int size;
} Region_t;

typedef struct Pixel_t
{
    int i;
    int x;
    int y;
    int label;
    Pixel_t *next;
} Pixel_t;

typedef struct Image_t
{
    int width, height, channels, size_pixels, size_labels;
    Pixel_t *pixel;
    unsigned char *pixel_rgb;
    double *pixel_lab;

    Region_t *region;

} Image_t;

Image_t *vx_image_create(   int width, int height, int init_intensity, 
                            unsigned char *fdata, int fdatai){

    printf("fdatai:: %d",fdatai);

    clock_t start = clock();
    clock_t end;

    Image_t *img = (Image_t *)malloc(sizeof(Image_t));
    img->channels = 3;
    img->width = width;
    img->height = height;
    img->size_pixels = width * height;
    img->pixel_rgb = (unsigned char *)malloc(img->size_pixels * 3 * sizeof(unsigned char));
    img->pixel_lab = (double *)malloc(img->size_pixels * 3 * sizeof(double));
    img->pixel = (Pixel_t *)malloc(img->size_pixels * sizeof(Pixel_t));
    img->region = (Region_t *)malloc(img->size_pixels * sizeof(Region_t));

    //img->pixel_rgb = fdata;
    // initializing attributes
    int x, y, i;
    int z = img->size_pixels;
    unsigned char R, G, B;
    double L_c, A_c, B_c;


    if (fdata==NULL){
        for (y = 0; y < img->height; ++y){
            for (x = 0; x < img->width; ++x){
                i = y * img->width + x;

                img->pixel[i].i = i;
                img->pixel[i].x = x;
                img->pixel[i].y = y;
                img->pixel[i].label = i;

                img->region[i].size = 1;
                // img->region[i].next = &img->pixel[i];

                img->pixel_rgb[0*z+i] = init_intensity;
                img->pixel_rgb[1*z+i] = init_intensity;
                img->pixel_rgb[2*z+i] = init_intensity;

                // img->pixel_rgb[i*img->channels+0] = init_intensity;
                // img->pixel_rgb[i*img->channels+1] = init_intensity;
                // img->pixel_rgb[i*img->channels+2] = init_intensity;

                // img->pixel_lab[i*chan+0] = L_c;
                // img->pixel_lab[i*chan+1] = A_c;
                // img->pixel_lab[i*chan+2] = B_c;

            }
        }
    }
    else{
        if (fdatai==1){
            for (y = 0; y < img->height; ++y){
                for (x = 0; x < img->width; ++x){
                    i = y * img->width + x;
                    //i = x * img->height + y;
                    img->pixel[i].i = i;
                    img->pixel[i].x = x;
                    img->pixel[i].y = y;
                    img->pixel[i].label = i;

                    img->region[i].size = 1;
                    // img->region[i].next = &img->pixel[i];

                    img->pixel_rgb[0*z+i] = fdata[0*z+i];
                    img->pixel_rgb[1*z+i] = fdata[0*z+i];
                    img->pixel_rgb[2*z+i] = fdata[0*z+i];

                    // img->pixel_rgb[0*z+i] = fdata[0*z+i];
                    // img->pixel_rgb[1*z+i] = fdata[0*z+i];
                    // img->pixel_rgb[2*z+i] = fdata[0*z+i];


                    R = img->pixel_rgb[0*z+i];
                    G = img->pixel_rgb[1*z+i];
                    B = img->pixel_rgb[2*z+i];
                    L_c = 0.0;A_c = 0.0;B_c = 0.0;
                    rgb2lab(&R, &G, &B, &L_c, &A_c, &B_c);
                    //printf("%f, %f, %f",L_c, A_c, B_c);
                    img->pixel_lab[0*z+i] = L_c;
                    img->pixel_lab[1*z+i] = A_c;
                    img->pixel_lab[2*z+i] = B_c;
                }
            }
        }
        else if (fdatai==3){
            for (y = 0; y < img->height; ++y){
                for (x = 0; x < img->width; ++x){
                    i = y * img->width + x;
                    //i = x * img->height + y;
                    img->pixel[i].i = i;
                    img->pixel[i].x = x;
                    img->pixel[i].y = y;
                    img->pixel[i].label = i;

                    img->region[i].size = 1;
                    // img->region[i].next = &img->pixel[i];

                    img->pixel_rgb[0*z+i] = fdata[0*z+i];
                    img->pixel_rgb[1*z+i] = fdata[1*z+i];
                    img->pixel_rgb[2*z+i] = fdata[2*z+i];

                    L_c = 0.0;A_c = 0.0;B_c = 0.0;
                    RGB2LABT(    &img->pixel_rgb[0*z+i],
                                &img->pixel_rgb[1*z+i], 
                                &img->pixel_rgb[2*z+i],
                                &img->pixel_lab[0*z+i],
                                &img->pixel_lab[1*z+i],
                                &img->pixel_lab[2*z+i]
                                );
                   //printf("%f, %f, %f",L_c, A_c, B_c);
                }
            }
        }
    }






    //create_lab_color(img);
    img->size_labels = img->size_pixels;
    end = clock();
    printf("time creando image: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    return img;
}

Image_t* vx_image_read(const char filename[]){
    clock_t start = clock();
    clock_t end;

    Image_t *img = (Image_t*)malloc(sizeof(Image_t));

    unsigned char * _rgb = stbi_load(filename, &img->width, &img->height, &img->channels, 0);
    if(_rgb == NULL){printf("Error in loading the image\n");exit(1);}
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", img->width, img->height, img->channels);

    img->size_pixels = img->width*img->height;
    img->pixel_rgb = (unsigned char *)malloc(img->size_pixels * 3 * sizeof(unsigned char));
    img->pixel_lab = (double*)malloc(img->size_pixels*img->channels*sizeof(double));
    img->pixel = (Pixel_t*)malloc(img->size_pixels*sizeof(Pixel_t));
    img->region = (Region_t*)malloc(img->size_pixels*sizeof(Region_t));

    // initializing attributes
    int x, y, i;
    int z = img->size_pixels;
    int chan = img->channels;
    //printf("img->channels %d", img->channels);
    unsigned char R, G, B;
    double L_c, A_c, B_c;
    for (y=0;y<img->height;++y){
        for (x=0;x<img->width;++x){
            i =y*img->width+x;
            img->pixel[i].i = i;
            img->pixel[i].x = x;
            img->pixel[i].y = y;
            img->pixel[i].label = i;

            img->region[i].size = 1;
            img->region[i].next = &img->pixel[i];
            
            
            R = _rgb[i*chan+0];
            G = _rgb[i*chan+1];
            B = _rgb[i*chan+2];
            img->pixel_rgb[0*z+i] = R;
            img->pixel_rgb[1*z+i] = G;
            img->pixel_rgb[2*z+i] = B;

            L_c = 0.0;
            A_c = 0.0;
            B_c = 0.0;
            rgb2lab(&R, &G, &B, &L_c, &A_c, &B_c);
            //printf("%f, %f, %f",L_c, A_c, B_c);
            img->pixel_lab[0*z+i] = L_c;
            img->pixel_lab[1*z+i] = A_c;
            img->pixel_lab[2*z+i] = B_c;
        }
    }
    img->size_labels = img->size_pixels;
    //create_lab_color(img);
    end = clock();
    printf("time read image: %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    free(_rgb);
    return img;
}

void vx_image_free(Image_t *img){
    free(img->pixel_rgb);
    free(img->pixel_lab);
    free(img->pixel);
    free(img->region);
    free(img);
}

// void create_lab_color(Image_t *image){
//     int i;
//     int size = image->size_pixels;
//     int chan = image->channels;

//     unsigned char R, G, B;
//     double L_c, A_c, B_c;
//     for(i=0; i<size; ++i){
//         R = image->pixel_rgb[i*chan+0];
//         G = image->pixel_rgb[i*chan+1];
//         B = image->pixel_rgb[i*chan+2];

//         L_c = 0.0;
//         A_c = 0.0;
//         B_c = 0.0;
//         rgb2lab(&R, &G, &B, &L_c, &A_c, &B_c);
//         //printf("%f, %f, %f",L_c, A_c, B_c);
//         image->pixel_lab[i*chan+0] = L_c;
//         image->pixel_lab[i*chan+1] = A_c;
//         image->pixel_lab[i*chan+2] = B_c;
//     }
// }

void vx_image_update_regions(Image_t *img){
    clock_t start = clock();
    clock_t end;

    int *aux = (int *)malloc(img->size_pixels * sizeof(int));
    int i, j, label;
    int count = 0;

    for (i = 0; i < img->size_pixels; ++i){
        aux[img->pixel[i].label] = -1;
        // img->region[i].next = NULL;
        img->region[i].size = 0;
    }

    for (i = 0; i < img->size_pixels; ++i){
        /* code */
        label = img->pixel[i].label;

        if (aux[label] == -1){
            aux[label] = count;
            count++;
        }
        j = aux[label];
        // img->pixel[i].next = &*(img->region[j].next);

        // img->region[j].next = &img->pixel[i];
        img->region[j].size++;

        img->pixel[i].label = j;
    }
    img->size_labels = count;

    free(aux);

    end = clock();
    printf("time update regions: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
}

 void vx_image_copy(unsigned char *datafrom, unsigned char *datato, int size){
    clock_t start = clock();
    clock_t end;

    int i;
    for (i=0; i<size; ++i){
        datato[i] = datafrom[i];
    }

    end = clock();
    printf("time copy image: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
 }

void vx_image_write(Image_t* img, const char filename[]){
    clock_t start = clock();
    clock_t end;

    //stbi_write_jpg(filename, img->width, img->height, img->channels, img->pixel_rgb, 0);

    // c.R = img->pixel_rgb[index * img->channels + 0];
    // c.R = img->pixel_rgb[index * img->channels + 1];
    // c.R = img->pixel_rgb[index * img->channels + 2];
    unsigned char *_rgb = (unsigned char *)malloc(img->size_pixels * 3 * sizeof(unsigned char));
    int i;
    for (i=0; i<img->size_pixels; ++i){
        _rgb[i*img->channels+0] = img->pixel_rgb[0*img->size_pixels+i];
        _rgb[i*img->channels+1] = img->pixel_rgb[1*img->size_pixels+i];
        _rgb[i*img->channels+2] = img->pixel_rgb[2*img->size_pixels+i];
    }
    stbi_write_png(filename, img->width, img->height, img->channels, _rgb, 0);
    free(_rgb);

    end = clock();
    printf("time write image: %f\n", (double)(end - start)/CLOCKS_PER_SEC);
}


void vx_image_setcolor(Image_t *img, int index, const Color3RGB color)
{
    // img->pixel_rgb[index * img->channels + 0] = color.R;
    // img->pixel_rgb[index * img->channels + 1] = color.G;
    // img->pixel_rgb[index * img->channels + 2] = color.B;

    img->pixel_rgb[0 * img->size_pixels + index] = color.R;
    img->pixel_rgb[1 * img->size_pixels + index] = color.G;
    img->pixel_rgb[2 * img->size_pixels + index] = color.B;
}

Color3RGB vx_image_getcolor(Image_t *img, int index)
{
    Color3RGB c;
    // c.R = img->pixel_rgb[index * img->channels + 0];
    // c.R = img->pixel_rgb[index * img->channels + 1];
    // c.R = img->pixel_rgb[index * img->channels + 2];
    c.R = img->pixel_rgb[0 * img->size_pixels + index];
    c.G = img->pixel_rgb[1 * img->size_pixels + index];
    c.B = img->pixel_rgb[2 * img->size_pixels + index];
    return c;
}


void vx_image_draw_regions(Image_t* img){
    clock_t start = clock();
    clock_t end;

    //Image_t* imgp = vx_image_create(img->width, img->height, 0, NULL, 0);
    Color3RGB* colors = makecolor(img->size_labels);
    Color3RGB c;
    int i;
    for (i = 0; i < img->size_pixels; ++i){
        c = colors[img->pixel[i].label];

        //c = vx_image_getcolor(img, i);
        //vx_image_setcolor(img, i, c);
        vx_image_setcolor(img, i, c);
        //printf("(%d, %d, %d, %d)",c.R,c.G,c.B, img->pixel[i].label);
        //printf("(%d)\n",img->pixel[i].label);
    }
    //vx_image_write(imgp, filename);
    free(colors);
    //vx_image_free(imgp);

    end = clock();
    printf("time write regions: %f\n", (double)(end - start)/CLOCKS_PER_SEC);

}


void rgb2lab(
    unsigned char *sR,
    unsigned char *sG,
    unsigned char *sB,
    double *lval,
    double *aval,
    double *bval)
{
    //------------------------
    // sRGB to XYZ conversion
    //------------------------
    double X, Y, Z;
    rgb2xyz(sR, sG, sB, &X, &Y, &Z);

    //------------------------
    // XYZ to LAB conversion
    //------------------------
    double epsilon = 0.008856; //actual CIE standard
    double kappa = 903.3;      //actual CIE standard

    double Xr = 0.950456; //reference white
    double Yr = 1.0;      //reference white
    double Zr = 1.088754; //reference white

    double xr = X / Xr;
    double yr = Y / Yr;
    double zr = Z / Zr;

    double fx, fy, fz;
    if (xr > epsilon)
        fx = pow(xr, 1.0 / 3.0);
    else
        fx = (kappa * xr + 16.0) / 116.0;
    if (yr > epsilon)
        fy = pow(yr, 1.0 / 3.0);
    else
        fy = (kappa * yr + 16.0) / 116.0;
    if (zr > epsilon)
        fz = pow(zr, 1.0 / 3.0);
    else
        fz = (kappa * zr + 16.0) / 116.0;

    *lval = 116.0 * fy - 16.0;
    *aval = 500.0 * (fx - fy);
    *bval = 200.0 * (fy - fz);
}

void rgb2xyz(
    unsigned char *sR,
    unsigned char *sG,
    unsigned char *sB,
    double *X,
    double *Y,
    double *Z)
{

    double R = *sR / 255.0;
    double G = *sG / 255.0;
    double B = *sB / 255.0;

    double r, g, b;

    if (R <= 0.04045)
        r = R / 12.92;
    else
        r = pow((R + 0.055) / 1.055, 2.4);
    if (G <= 0.04045)
        g = G / 12.92;
    else
        g = pow((G + 0.055) / 1.055, 2.4);
    if (B <= 0.04045)
        b = B / 12.92;
    else
        b = pow((B + 0.055) / 1.055, 2.4);

    *X = r * 0.4124564 + g * 0.3575761 + b * 0.1804375;
    *Y = r * 0.2126729 + g * 0.7151522 + b * 0.0721750;
    *Z = r * 0.0193339 + g * 0.1191920 + b * 0.9503041;
}









inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

inline double fastPow2(double a, double b) {
    return exp(b*log(a));
}

inline double fastPow3(double a, double b) {
  // calculate approximation with fraction of the exponent

  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;

  // exponentiation by squaring with the exponent's integer part
  // double r = u.d makes everything much slower, not sure why
  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }

  return r * u.d;
}

double H(double q)
{
	double value;
	if ( q > 0.008856 ) {
		// value = pow ( q, 0.333333 );
		// return value;

        //value = exp(log(q) * 0.333333);
        //value = fastPow2(q, 0.333333);
        value = fastPow3(q, 0.333333);
        //value = fastPow(q, 0.333333);
        return value;

		// return 0.5;
	}
	else {
		value = 7.787*q + 0.137931;
		return value;
		
        // return 0.5;
	}
    //return 0.5;
}


void RGB2LABT(  unsigned char *R_v,
                unsigned char *G_v,
                unsigned char *B_v,
                double *Lv, double *Bv, double *Av){
	double RGB[3];
	double XYZ[3];
	double Lab[3];
	double adapt[3];

	adapt[0] = 0.950467;
	adapt[1] = 1.000000;
	adapt[2] = 1.088969;

	RGB[0] = *R_v * 0.003922;
	RGB[1] = *G_v * 0.003922;
	RGB[2] = *B_v * 0.003922;

	XYZ[0] = 0.412424 * RGB[0] + 0.357579 * RGB[1] + 0.180464 * RGB[2];
	XYZ[1] = 0.212656 * RGB[0] + 0.715158 * RGB[1] + 0.0721856 * RGB[2];
	XYZ[2] = 0.0193324 * RGB[0] + 0.119193 * RGB[1] + 0.950444 * RGB[2];

	*Lv = 116 * H( XYZ[1] / adapt[1] ) - 16;
	*Bv = 500 * ( H( XYZ[0] / adapt[0] ) - H ( XYZ[1] / adapt[1] ) );
	*Av = 200 * ( H( XYZ[1] / adapt[1] ) - H ( XYZ[2] / adapt[2] ) );
}




void vx_image_norm_zscore(Image_t* img){
    double a;
    double mean = 0.0;
    double sigma = 0.0;
    double eps = 1e-6;
    double sx = 0.0;


    int i, j;
    int n = img->size_pixels;
    for(j=0; j<3; ++j){
        sx = 0.0;
        for(i=0; i<n; ++i){
            sx += img->pixel_lab[j*n+i];
        }
        
        mean = sx/n;
        
        for(i=0; i<n; ++i){
            a = img->pixel_lab[j*n+i] - mean;
            sigma += a*a;
        }
        sigma = sqrt( sigma/(n-1) );

        if (sigma < eps){
            sigma = 1.0;
        }
        for(i=0; i<n; ++i){
            img->pixel_lab[j*n+i] = (img->pixel_lab[j*n+i] - mean)/(sigma);
        }
    }
}
void vx_image_to_colormoments_features(Image_t* image, double *features){
    int sr = image->size_labels;
    int a, i, j, l;
    int size = image->size_pixels;
    int chanel = 3;
    //explore each region
    for (i=0; i<size; ++i){
        l = image->pixel[i].label;
        features[(l*9)+0] += image->pixel_lab[0*size+i];
        features[(l*9)+1] += image->pixel_lab[1*size+i];
        features[(l*9)+2] += image->pixel_lab[2*size+i];
    }
    
    //medias
    for (i=0; i<sr; ++i){
        features[(i*9)+0] /= image->region[i].size;
        features[(i*9)+1] /= image->region[i].size;
        features[(i*9)+2] /= image->region[i].size;
    }

    //explore each region
    for (i=0; i<size; ++i){
        l = image->pixel[i].label;
        features[(l*9)+3] += pow(abs(image->pixel_lab[0*size+i] - features[(l*9)+0]), 2.0);
        features[(l*9)+4] += pow(abs(image->pixel_lab[1*size+i] - features[(l*9)+1]), 2.0);
        features[(l*9)+5] += pow(abs(image->pixel_lab[2*size+i] - features[(l*9)+2]), 2.0);

        features[(l*9)+6] += pow(abs(image->pixel_lab[0*size+i] - features[(l*9)+0]), 3.0);
        features[(l*9)+7] += pow(abs(image->pixel_lab[1*size+i] - features[(l*9)+1]), 3.0);
        features[(l*9)+8] += pow(abs(image->pixel_lab[2*size+i] - features[(l*9)+2]), 3.0);
    }

    //medias
    double e2 = 1.0/2.0;
    double e4 = 1.0/3.0;
    for (i=0; i<sr; ++i){
        
        features[(i*9)+3] /= image->region[i].size;
        features[(i*9)+4] /= image->region[i].size;
        features[(i*9)+5] /= image->region[i].size;

        features[(i*9)+6] /= image->region[i].size;
        features[(i*9)+7] /= image->region[i].size;
        features[(i*9)+8] /= image->region[i].size;


        features[(i*9)+3] = pow(abs(features[(i*9)+3]), e2);
        features[(i*9)+4] = pow(abs(features[(i*9)+4]), e2);
        features[(i*9)+5] = pow(abs(features[(i*9)+5]), e2);

        features[(i*9)+6] = pow(abs(features[(i*9)+6]), e4);
        features[(i*9)+7] = pow(abs(features[(i*9)+7]), e4);
        features[(i*9)+8] = pow(abs(features[(i*9)+8]), e4);
    }
}

#endif