#ifndef GIRIH_STENCILS_
#define GIRIH_STENCILS_

///
/// \brief ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
///
template<typename Real>
inline void space_8_time_2_constant_coef(const int i,
                                         const int nnx,
                                         const int nnxy,
                                         const int nnxyz,
                                         Real* __restrict__ ux,
                                         const Real* __restrict__ vx,
                                         const Real* __restrict__ coef,
                                         const Real* __restrict__ roc) {
    ux[i] = (Real (2.0)) * vx[i] - ux[i] 
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

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
///
template<typename Real>
inline void space_2_time_1_constant_coef(const int i,
                                         const int nnx,
                                         const int nnxy,
                                         const int nnxyz,
                                         Real* __restrict__ ux,
                                         const Real* __restrict__ vx,
                                         const Real* __restrict__ coef,
                                         const Real* __restrict__ roc) {
    ux[i] = coef[0] * vx[i] + coef[1] * (vx[i+1]    + vx[i-1]   ) 
                            + coef[1] * (vx[i+nnx]  + vx[i-nnx] ) 
                            + coef[1] * (vx[i+nnxy] + vx[i-nnxy]);
}

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
///
template<typename Real>
inline void space_2_time_1_variable_coef(const int i,
                                         const int nnx,
                                         const int nnxy,
                                         const int nnxyz,
                                         Real* __restrict__ ux,
                                         const Real* __restrict__ vx,
                                         const Real* __restrict__ coef,
                                         const Real* __restrict__ roc) {
    ux[i] = coef[i] * vx[i] + coef[i+nnxyz] * (vx[i+1]    + vx[i-1]   ) 
                            + coef[i+nnxyz] * (vx[i+nnx]  + vx[i-nnx] ) 
                            + coef[i+nnxyz] * (vx[i+nnxy] + vx[i-nnxy]);
}

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with variable axis symmetric coefficient
///
template<typename Real>
inline void space_2_time_1_variable_coef_axis_symm(const int i,
                                                   const int nnx,
                                                   const int nnxy,
                                                   const int nnxyz,
                                                   Real* __restrict__ ux,
                                                   const Real* __restrict__ vx,
                                                   const Real* __restrict__ coef,
                                                   const Real* __restrict__ roc) {
    ux[i] = coef[i] * vx[i] + coef[i+nnxyz  ] * (vx[i+1]    + vx[i-1]   ) 
                            + coef[i+nnxyz*2] * (vx[i+nnx]  + vx[i-nnx] ) 
                            + coef[i+nnxyz*3] * (vx[i+nnxy] + vx[i-nnxy]);
}

///
/// \brief ISO stencil 8th-order-in-space-1st-order-in-time with variable axis symmetric coefficient
///
template<typename Real>
inline void space_8_time_1_variable_coef_axis_symm(const int i,
                                                   const int nnx,
                                                   const int nnxy,
                                                   const int nnxyz,
                                                   Real* __restrict__ ux,
                                                   const Real* __restrict__ vx,
                                                   const Real* __restrict__ coef,
                                                   const Real* __restrict__ roc) {
    ux[i] = coef[i] * vx[i] + coef[i+nnxyz   ] * (vx[i+1]      + vx[i-1]     )
                            + coef[i+nnxyz*2 ] * (vx[i+  nnx]  + vx[i-  nnx] ) 
                            + coef[i+nnxyz*3 ] * (vx[i+  nnxy] + vx[i-  nnxy])
                            + coef[i+nnxyz*4 ] * (vx[i+2]      + vx[i-2]     )
                            + coef[i+nnxyz*5 ] * (vx[i+2*nnx]  + vx[i-2*nnx] ) 
                            + coef[i+nnxyz*6 ] * (vx[i+2*nnxy] + vx[i-2*nnxy])
                            + coef[i+nnxyz*7 ] * (vx[i+3]      + vx[i-3]     )
                            + coef[i+nnxyz*8 ] * (vx[i+3*nnx]  + vx[i-3*nnx] ) 
                            + coef[i+nnxyz*9 ] * (vx[i+3*nnxy] + vx[i-3*nnxy])
                            + coef[i+nnxyz*10] * (vx[i+4]      + vx[i-4]     )
                            + coef[i+nnxyz*11] * (vx[i+4*nnx]  + vx[i-4*nnx] ) 
                            + coef[i+nnxyz*12] * (vx[i+4*nnxy] + vx[i-4*nnxy]);
}

///
/// \brief ISO stencil 2nd-order-in-space-1st-order-in-time with variable no symmetry coefficient 
///
template<typename Real>
inline void space_2_time_1_variable_coef_no_symm(const int i,
                                                 const int nnx,
                                                 const int nnxy,
                                                 const int nnxyz,
                                                 Real* __restrict__ ux,
                                                 const Real* __restrict__ vx,
                                                 const Real* __restrict__ coef,
                                                 const Real* __restrict__ roc) {
    ux[i] = coef[i        ] * vx[i] 
          + coef[i+nnxyz*1] * vx[i-1]
          + coef[i+nnxyz*2] * vx[i+1]
          + coef[i+nnxyz*3] * vx[i-nnx]
          + coef[i+nnxyz*4] * vx[i+nnx]
          + coef[i+nnxyz*5] * vx[i-nnxy] 
          + coef[i+nnxyz*6] * vx[i+nnxy];
}

///
/// \brief ISO Box stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
///
template<typename Real>
inline void box_space_8_time_1_variable_coef_axis_symm(const int i,
                                                       const int nnx,
                                                       const int nnxy,
                                                       const int nnxyz,
                                                       Real* __restrict__ ux,
                                                       const Real* __restrict__ vx,
                                                       const Real* __restrict__ coef,
                                                       const Real* __restrict__ roc) {
    ux[i] = coef[i] * vx[i] + coef[i+nnxyz  ] * (vx[i+1         ] + vx[i-1         ])
                            + coef[i+nnxyz  ] * (vx[i  +nnx     ] + vx[i  -nnx     ]) 
                            + coef[i+nnxyz  ] * (vx[i      +nnxy] + vx[i      -nnxy])
                            + coef[i+nnxyz*2] * (vx[i+1    -nnxy] + vx[i-1    -nnxy])
                            + coef[i+nnxyz*2] * (vx[i  +nnx-nnxy] + vx[i  -nnx-nnxy]) 
                            + coef[i+nnxyz*2] * (vx[i+1+nnx     ] + vx[i-1-nnx     ]) 
                            + coef[i+nnxyz*2] * (vx[i+1-nnx     ] + vx[i-1+nnx     ])
                            + coef[i+nnxyz*2] * (vx[i+1    +nnxy] + vx[i-1    +nnxy])
                            + coef[i+nnxyz*2] * (vx[i  +nnx+nnxy] + vx[i  -nnx+nnxy]) 
                            + coef[i+nnxyz*3] * (vx[i+1+nnx+nnxy] + vx[i-1-nnx-nnxy])
                            + coef[i+nnxyz*3] * (vx[i+1-nnx+nnxy] + vx[i-1+nnx-nnxy])
                            + coef[i+nnxyz*3] * (vx[i-1-nnx+nnxy] + vx[i+1+nnx-nnxy])
                            + coef[i+nnxyz*3] * (vx[i-1+nnx+nnxy] + vx[i+1-nnx-nnxy]);
}

#endif