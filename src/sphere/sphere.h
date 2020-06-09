#ifndef SPHERE_H
#define SPHERE_H

#include "data_utility.h"
#include "search.h"
#include "operation.h"

#include "lp.h"

// The complete Sphere algorithm
point_set_t* sphereWSImpLP(point_set_t* point_set, int k);

#endif