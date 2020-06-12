#ifndef GIRIH_STENCILS_
#define GIRIH_STENCILS_

#include <cstdint>

///
/// \brief Abstract class using as an interface for all stencils
///
template<typename real_t>
class StencilsStrategy {
public:
    virtual ~StencilsStrategy() {};
    virtual void ComputeStencils(const int xmin, const int xmax,
                                 const int ymin, const int ymax,
                                 const int tmin, const int tmax,const int dt,
                                 const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                                 real_t* __restrict__ u,
                                 const real_t* __restrict__ v,
                                 const real_t* __restrict__ coef,
                                 const real_t* __restrict__ roc) const = 0;
};

///
/// \brief Context of strategy design pattern for selecting stencil method
///
template<typename real_t>
class StencilContext {
private:
    StencilsStrategy<real_t> *strategy_;

public:
    StencilContext(StencilsStrategy<real_t> *strategy = nullptr) : strategy_(strategy) {};
    ~StencilContext() {
        delete this->strategy_;
    }

    void set_strategy(StencilsStrategy<real_t> *strategy) {
        delete this->strategy_;
        this->strategy_ = strategy;
    }

    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax, const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const {
        this->strategy_->ComputeStencils(xmin,xmax,ymin,ymax,tmin,tmax,dt,nnx,nnxy,nnxyz,u,v,coef,roc);
    } 
};

///
/// \brief ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
///
template<typename real_t>
class space_8_time_2_constant_coef : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax, const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override{

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u[1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v[1ULL*k*nnxy + j*nnx]);
        
                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = (real_t (2.0)) * vx[i] - ux[i] 
                            + roc[i] * ( coef[0] * vx[i] 
                                        + coef[1] * (vx[i+1]       + vx[i-1]     )          
                                        + coef[1] * (vx[i+nnx]     + vx[i-nnx]   )       
                                        + coef[1] * (vx[i+nnxy]    + vx[i-nnxy]  )     
                                        + coef[2] * (vx[i+2]       + vx[i-2]     )         
                                        + coef[2] * (vx[i+2*nnx]   + vx[i-2*nnx] )     
                                        + coef[2] * (vx[i+2*nnxy]  + vx[i-2*nnxy])   
                                        + coef[3] * (vx[i+3]       + vx[i-3]     )         
                                        + coef[3] * (vx[i+3*nnx]   + vx[i-3*nnx] )     
                                        + coef[3] * (vx[i+3*nnxy]  + vx[i-3*nnxy])   
                                        + coef[4] * (vx[i+4]       + vx[i-4]     )          
                                        + coef[4] * (vx[i+4*nnx]   + vx[i-4*nnx] )      
                                        + coef[4] * (vx[i+4*nnxy]  + vx[i-4*nnxy]) 
                                        );
                }
            }
        }
    }
};

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
///
template<typename real_t>
class space_2_time_1_constant_coef : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax, const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override {

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u[1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v[1ULL*k*nnxy + j*nnx]);
        
                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = coef[0] * vx[i] + coef[1] * (vx[i+1]    + vx[i-1]   ) 
                                            + coef[1] * (vx[i+nnx]  + vx[i-nnx] ) 
                                            + coef[1] * (vx[i+nnxy] + vx[i-nnxy]);
                }
            }
        }
    }
};

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
///
template<typename real_t>
class space_2_time_1_variable_coef : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax,const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override {

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ cx = &(coef[1ULL*k*nnxy + j*nnx]);

                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = cx[i] * vx[i] + cx[i+nnxyz] * (vx[i+1]    + vx[i-1]   ) 
                                          + cx[i+nnxyz] * (vx[i+nnx]  + vx[i-nnx] ) 
                                          + cx[i+nnxyz] * (vx[i+nnxy] + vx[i-nnxy]);
                }
            }
        }
    }
};

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with variable axis symmetric coefficient
///
template<typename real_t>
class space_2_time_1_variable_coef_axis_symm : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax,const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override {

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ cx = &(coef[1ULL*k*nnxy + j*nnx]);

                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = cx[i] * vx[i] + cx[i+nnxyz  ] * (vx[i+1]    + vx[i-1]   ) 
                                          + cx[i+nnxyz*2] * (vx[i+nnx]  + vx[i-nnx] ) 
                                          + cx[i+nnxyz*3] * (vx[i+nnxy] + vx[i-nnxy]);
                }
            }
        }
    }
};

///
/// \brief ISO stencil 8th-order-in-space-1st-order-in-time with variable axis symmetric coefficient
///
template<typename real_t>
class space_8_time_1_variable_coef_axis_symm : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax,const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override {

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ cx = &(coef[1ULL*k*nnxy + j*nnx]);

                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = cx[i] * vx[i] + cx[i+nnxyz   ] * (vx[i+1]      + vx[i-1]     )
                                          + cx[i+nnxyz*2 ] * (vx[i+  nnx]  + vx[i-  nnx] ) 
                                          + cx[i+nnxyz*3 ] * (vx[i+  nnxy] + vx[i-  nnxy])
                                          + cx[i+nnxyz*4 ] * (vx[i+2]      + vx[i-2]     )
                                          + cx[i+nnxyz*5 ] * (vx[i+2*nnx]  + vx[i-2*nnx] ) 
                                          + cx[i+nnxyz*6 ] * (vx[i+2*nnxy] + vx[i-2*nnxy])
                                          + cx[i+nnxyz*7 ] * (vx[i+3]      + vx[i-3]     )
                                          + cx[i+nnxyz*8 ] * (vx[i+3*nnx]  + vx[i-3*nnx] ) 
                                          + cx[i+nnxyz*9 ] * (vx[i+3*nnxy] + vx[i-3*nnxy])
                                          + cx[i+nnxyz*10] * (vx[i+4]      + vx[i-4]     )
                                          + cx[i+nnxyz*11] * (vx[i+4*nnx]  + vx[i-4*nnx] ) 
                                          + cx[i+nnxyz*12] * (vx[i+4*nnxy] + vx[i-4*nnxy]);
                }
            }
        }
    }
};

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with variable no symmetry coefficient 
///
template<typename real_t>
class space_2_time_1_variable_coef_no_symm : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax,const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override {

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v   [1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ cx = &(coef[1ULL*k*nnxy + j*nnx]);

                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = cx[i        ] * vx[i] 
                          + cx[i+nnxyz*1] * vx[i-1]
                          + cx[i+nnxyz*2] * vx[i+1]
                          + cx[i+nnxyz*3] * vx[i-nnx]
                          + cx[i+nnxyz*4] * vx[i+nnx]
                          + cx[i+nnxyz*5] * vx[i-nnxy] 
                          + cx[i+nnxyz*6] * vx[i+nnxy];
                }
            }
        }
    }
};

///
/// \brief ISO Box stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
///
template<typename real_t>
class box_space_8_time_1_variable_coef_axis_symm : public StencilsStrategy<real_t> {
public :
     
    void ComputeStencils(const int xmin, const int xmax,
                         const int ymin, const int ymax,
                         const int tmin, const int tmax, const int dt,
                         const int nnx, const uint64_t nnxy, const uint64_t nnxyz,
                         real_t* __restrict__ u,
                         const real_t* __restrict__ v,
                         const real_t* __restrict__ coef,
                         const real_t* __restrict__ roc) const override{

        for(int k = tmin; k < tmax; k += dt) {
            for(int j = ymin; j < ymax; j++) {
                real_t * __restrict__ ux = &(u[1ULL*k*nnxy + j*nnx]);
                real_t * __restrict__ vx = &(v[1ULL*k*nnxy + j*nnx]);
        
                #pragma simd
                for(int i = xmin; i < xmax; i++) {
                    ux[i] = coef[0] * vx[i] + coef[1] * (vx[i+1         ] + vx[i-1         ])
                                            + coef[1] * (vx[i  +nnx     ] + vx[i  -nnx     ]) 
                                            + coef[1] * (vx[i      +nnxy] + vx[i      -nnxy])
                                            + coef[2] * (vx[i+1    -nnxy] + vx[i-1    -nnxy])
                                            + coef[2] * (vx[i  +nnx-nnxy] + vx[i  -nnx-nnxy]) 
                                            + coef[2] * (vx[i+1+nnx     ] + vx[i-1-nnx     ]) 
                                            + coef[2] * (vx[i+1-nnx     ] + vx[i-1+nnx     ])
                                            + coef[2] * (vx[i+1    +nnxy] + vx[i-1    +nnxy])
                                            + coef[2] * (vx[i  +nnx+nnxy] + vx[i  -nnx+nnxy]) 
                                            + coef[3] * (vx[i+1+nnx+nnxy] + vx[i-1-nnx-nnxy])
                                            + coef[3] * (vx[i+1-nnx+nnxy] + vx[i-1+nnx-nnxy])
                                            + coef[3] * (vx[i-1-nnx+nnxy] + vx[i+1+nnx-nnxy])
                                            + coef[3] * (vx[i-1+nnx+nnxy] + vx[i+1-nnx-nnxy]);
                }
            }
        }
    }
};

#endif