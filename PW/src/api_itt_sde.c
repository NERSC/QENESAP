#ifdef __USE_VTUNE
#include "ittnotify.h"
#endif

void fortran_sde_start()
{
#ifdef __USE_SDE
   __SSC_MARK(0x111);
#else
#warning Not using SDE markers
#endif
}  

void fortran_sde_stop()
{
#ifdef __USE_SDE
  __SSC_MARK(0x222);
#endif
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
