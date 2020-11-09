//#include "vx/com/cp/image/Imagecpp.h"

#include "vx/pmglp/MGLP.h"


int main(int argc, char* *argv) {
    MGLP *mglp = new MGLP();
    mglp->in_filename = std::string(argv[1]);
    mglp->ou_filename = std::string(argv[2]);
    mglp->sp_side = atoi(argv[3]);
    mglp->threshold = atof(argv[4]);
    mglp->traversal = atof(argv[5]);
    
    mglp->execute();
    delete mglp;


    // Image* img = new Image(pathimginput);
    // img->display();
    // delete img;

}