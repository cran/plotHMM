#include <math.h>//for INFINITY.

// http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf

double elnsum(double elnx, double elny){
  if(elnx == -INFINITY){
    return elny;
  }else if(elny == -INFINITY){
    return elnx;
  }else if(elny < elnx){
    return elnx + log(1+exp(elny-elnx));
  }else{
    return elny + log(1+exp(elnx-elny));
  }
}

double elnproduct(double elnx, double elny){
  if(elnx == -INFINITY || elny == -INFINITY){
    return -INFINITY;
  }else{
    return elnx+elny;
  }
}
