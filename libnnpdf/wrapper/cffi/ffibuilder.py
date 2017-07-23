""" ffibuilder.py
Gernerate sources for CFFI style bindings
"""
import argparse
import cffi

I = cffi.FFI()

#TODO: Maybe move to independent files?
I.set_source("NNPDF.ffi",
"""

#include <cmath>
#include "NNPDF/lhapdfset.h"

extern "C"{

static const int all_flavours[] = {-6,-5,-4,-3,-2,-1,21,1,2,3,4,5,6,22};

double xfxQ(void *pdf , double x, double Q, int member, int fl){
    try{
        return ((NNPDF::LHAPDFSet*)pdf)->xfxQ(x,Q,member,fl);
    }catch(...){
        return NAN;
    }

}

struct along_x_user_data{
    void *pdf;
    double Q;
    int member;
    int fl;
};

double xfxQ_along_x(double x, void *user_data){
    auto dt = (struct along_x_user_data*)user_data;
    return xfxQ(dt->pdf, x, dt->Q, dt->member, dt->fl);
}

double xfxQ_along_x_sum_all(double x, void* user_data){
    double res = 0;
    //TODO: This should probably reset user_data->fl
    auto dt = (struct along_x_user_data*)user_data;
    for (auto fl : all_flavours){
        dt->fl = fl;
        res += xfxQ_along_x(x, dt);
    }
    return res;
}

double xfxQ_valence_along_x(double x, void *user_data){
    auto dt = (struct along_x_user_data*)user_data;
    return xfxQ(dt->pdf, x, dt->Q, dt->member, dt->fl) -
           xfxQ(dt->pdf, x, dt->Q, dt->member, -dt->fl);
}

double fxQ_along_x(double x, void *user_data){
    return xfxQ_along_x(x, user_data)/x;
}

double fxQ_valence_along_x(double x, void *user_data){
    return xfxQ_valence_along_x(x, user_data)/x;
}

}// extern C
"""
,source_extension='.cpp')

I.cdef("""
double xfxQ(void *pdf , double x, double Q, int member, int fl);
double xfxQ_along_x(double x, void *user_data);
double xfxQ_along_x_sum_all(double x, void* user_data);
double xfxQ_valence_along_x(double x, void *user_data);
double fxQ_along_x(double x, void *user_data);
double fxQ_valence_along_x(double x, void *user_data);
struct along_x_user_data{
    void *pdf;
    double Q;
    int member;
    int fl;
};
""")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("path")
    args = p.parse_args()
    I.emit_c_code(args.path)
