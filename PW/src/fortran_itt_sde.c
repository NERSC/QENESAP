#ifdef __USE_VTUNE
#include "ittnotify.h"
#endif

void fortran_sde_start(int mark)
{
  // __SSC_MARK(mark);
}  

void fortran_sde_stop(int mark)
{
  // __SSC_MARK(mark);
}  

void fortran_itt_resume()
{
#ifdef __USE_VTUNE
  __itt_resume();
#endif
}

void fortran_itt_pause()
{
#ifdef __USE_VTUNE
  __itt_pause();
#endif
}
