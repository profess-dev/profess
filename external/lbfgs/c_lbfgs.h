#ifndef __C_LBFGS_H__
#define __C_LBFGS_H__

extern "C"
{
void c_lbfgs_step(int*, int*, double*, double*, double*, int*, double*,
                  int*, double*, double*, double*, int*);
}

#endif // __C_LBFGS_H

