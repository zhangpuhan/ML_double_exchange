//
// Created by Puhan Zhang on 12/16/19.
//

#ifndef DE_C_TORCH_HELP_FUNCTION_H
#define DE_C_TORCH_HELP_FUNCTION_H

#include <iostream>
#include "lattice_base.h"

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

#endif //DE_C_TORCH_HELP_FUNCTION_H
