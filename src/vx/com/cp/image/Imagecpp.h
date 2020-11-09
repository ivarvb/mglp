
extern "C"{
    #include "../../cc/image/Imagec.h"
}

#ifndef IMAGECPP_H
#define IMAGECPP_H

#include <iostream>
#include <fstream>
#include <sstream>

#include "CImg.h"
using namespace cimg_library;
typedef CImg<PixelData> ImageCImg;

typedef unsigned char PixelData;

class Image{
private:
    // ImageCImg *_buff;
    Image_t *_img;
    std::string _img_name;

public:
    //from file
    Image(std::string filename){
        std::string ext = filename.substr(filename.find_last_of(".") + 1);
        if (ext == "seg"){
            //read seg file
        }
        else if (ext == "bou"){
            //read boud file
        }
        else{

            clock_t start = clock(); clock_t end;
            ImageCImg *datar = new ImageCImg();
            datar->load(filename.c_str());
            _img = vx_image_create(datar->width(), datar->height(),
                                   0, datar->_data, datar->spectrum());
            delete datar;
            end = clock();
            printf("time read image: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
            
            // _img = vx_image_read(filename.c_str());

        }
    }

    //from especific size
    Image(int width, int height, int intensity = 0){
        _img = vx_image_create(width, height, intensity, NULL, 3);
        // _buff = new ImageCImg(width, height, intensity, 3, 1);
    }

    ~Image(){
        // delete _buff;
        vx_image_free(_img);
    }

    Image_t *imageio_t(){
        return _img;
    };

    std::string &name(){
        return this->_img_name;
    }

    inline int i(int index){
        return _img->pixel[index].i;
    }

    inline int x(int index){
        return _img->pixel[index].x;
    }

    inline int y(int index){
        return _img->pixel[index].y;
    }

    void save(std::string path){
        //_img->save_png(path.c_str());
    }

    inline int width(){
        return _img->width;
    }

    inline int height(){
        return _img->height;
    }

    inline int size(){
        return _img->size_pixels;
    }

    inline int& label(int i){
        return _img->pixel[i].label;
    }

    inline int& size_labels(){
        return _img->size_labels;
    }


    // inline Color3RGB getPixel(int x, int y){
    //     Color3RGB c;
    //     c.R = _img->pixel_rgb[0*size()+(y*width()+x)];
    //     c.G = _img->pixel_rgb[1*size()+(y*width()+x)];
    //     c.B = _img->pixel_rgb[2*size()+(y*width()+x)];
    //     return c;
    // }

    // inline void sePixel(int x, int y, Color3RGB c){
    //     _img->pixel_rgb[0*size()+(y*width()+x)] = c.R;
    //     _img->pixel_rgb[1*size()+(y*width()+x)] = c.G;
    //     _img->pixel_rgb[2*size()+(y*width()+x)] = c.B;
    // }

    // inline void setPixel(int i, Color3RGB c){
    //     _img->pixel_rgb[0*size()+i] = c.R;
    //     _img->pixel_rgb[1*size()+i] = c.G;
    //     _img->pixel_rgb[2*size()+i] = c.B;
    // }

    // //////////////////////////////////////////////////
    // //////////////////////////////////////////////////
    // inline Color3RGB getBuffPixel(int x, int y){
    //     Color3RGB c;
    //     c.R = _buff->_data[0*size()+(y*width()+x)];
    //     c.G = _buff->_data[1*size()+(y*width()+x)];
    //     c.B = _buff->_data[2*size()+(y*width()+x)];
    //     return c;
    // }

    // inline void setBuffPixel(int x, int y, Color3RGB c){
    //     _buff->_data[0*size()+(y*width()+x)] = c.R;
    //     _buff->_data[1*size()+(y*width()+x)] = c.G;
    //     _buff->_data[2*size()+(y*width()+x)] = c.B;
    // }

    // inline void setBuffPixel(int i, Color3RGB c){
    //     _buff->_data[0*size()+i] = c.R;
    //     _buff->_data[1*size()+i] = c.G;
    //     _buff->_data[2*size()+i] = c.B;
    // }


    void display(){
        // ImageCImg *_buff = new ImageCImg(width(), height(), 1, 3, 1);
        // vx_image_copy(_img->pixel_rgb, _buff->_data, this->size() * 3);
        // //std::memcpy(_buff->_data, _img->pixel_rgb, _img->size_pixels);

        // std::stringstream stm;
        // stm << this->width() << " x " << this->height();
        // _buff->display(stm.str().c_str(), false);





        // std::string path = "outcimg.jpg";
        // _buff->save_jpeg(path.c_str());

        // // _buff->display(stm.str().c_str(), false);
        // delete _buff;

        // std::string pathf = "outxxxx.jpg";
        // vx_image_write(_img, pathf.c_str());
    }


    void write(std::string out_file){
        ImageCImg *_buff = new ImageCImg(width(), height(), 1, 3, 1);
        vx_image_copy(_img->pixel_rgb, _buff->_data, this->size() * 3);
        //std::memcpy(_buff->_data, _img->pixel_rgb, _img->size_pixels*3);

        // std::stringstream stm;
        // stm << this->width() << " x " << this->height();
        // _buff->display(stm.str().c_str(), false);


        _buff->save_png(out_file.c_str());

        // _buff->display(stm.str().c_str(), false);
        delete _buff;

        // std::string pathf = "outxxxx.jpg";
        // vx_image_write(_img, pathf.c_str());
    }

    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
};
#endif
