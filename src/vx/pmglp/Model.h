

#ifndef MODEL_H
#define	MODEL_h

//class Configure;

class Model{
public:
    std::string in_filename;
    std::string ou_filename;
    int ag_type;// algorithm
    int sg_type;// type segmentation: automatic(0) or interactive (1)
    int sp_side;// superpixels side
    int fd_type;// feature descriptor
    int px_type;// proximity
    double threshold;// threshold
    int traversal;// proximity

};
#endif

