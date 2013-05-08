      REAL sflux(NPTS,NSFLXS,5,0:NJDT),VAF(NPTS,NSFLXSP2)
      common/ flx_sfc   / sflux,VAF

