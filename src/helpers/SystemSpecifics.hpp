/*
 * SystemSpecifics.hpp
 *
 *  Created on: Dec 10, 2019
 *      Author: jsteinbu
 */

#ifndef SYSTEMSPECIFICS_HPP_
#define SYSTEMSPECIFICS_HPP_

#include <sys/resource.h>

/*int setDataLimitInMB(rlim_t value){
   rlimit rlim;
   rlim.rlim_cur = value*(1<<20);
   rlim.rlim_max = value*(1<<20);
   return setrlimit(RLIMIT_DATA, &rlim);
}
*/

#endif /* SYSTEMSPECIFICS_HPP_ */
