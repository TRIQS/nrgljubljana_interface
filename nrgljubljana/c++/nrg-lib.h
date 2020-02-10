#ifndef _nrg_lib_h_
#define _nrg_lib_h_

// This header is included for both executable and library.

#define NRG_COMMON

#include <string>

void run_nrg_master();
void run_nrg_slave();
void set_workdir(const std::string &workdir);
namespace time_mem {
  void timing_report();
  void memory_report();
}

#endif

