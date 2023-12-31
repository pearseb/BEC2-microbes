      function get_frac_year(days_elapsed)

      implicit none
      real days_elapsed, get_frac_year
#include "cppdefs.h"

! This function computes the number of years elapsed after days_elapsed
! days. The result is a real number.

#ifdef USE_REAL_YEAR
      get_frac_year = days_elapsed/365.25
#else
      get_frac_year = days_elapsed/360
#endif
      end
