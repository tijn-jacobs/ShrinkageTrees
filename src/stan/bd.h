
#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "Info.h"
#include "StanTree.h"
#include "StanTreeFunctions.h"
#include "StanForestFunctions.h"

bool bd(StanTree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
	std::vector<size_t>& nv, std::vector<double>& pv, bool aug, Random& random);

#endif
