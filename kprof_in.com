c       REAL ustar(npts),B0(npts),buoy(npts,NZP1tmax),B0sol(npts)
c       common/ kprof_in  / ustar,B0,buoy,B0sol
	REAL ustar,B0,B0sol
