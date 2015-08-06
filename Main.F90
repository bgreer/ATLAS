
! ATLAS CODE

PROGRAM FRACK ! ignore name..

	USE ParseInput
	USE Communication_Library
	USE Track
	USE Grid
	USE File_Management
	USE Timing
	USE Fit
	USE omp_lib
	IMPLICIT NONE

	INTEGER :: ierr, ii, ij, it, tilecount, currsize, ind
	INTEGER :: currtile, ix, iy, completed, tilesperproc, starttile
	INTEGER :: endtile, num, nproc, myid
	DOUBLE PRECISION :: starttime, endtime
	INTEGER :: ctest

	! Initialize MPI and timer
		! Communication_Library.F90
	CALL Initialize_Communication(myid, nproc)
	starttime = MPI_Wtime()
		! Timing.F90
	CALL InitTimer(myid,10)
	CALL AddTime(1) ! adds tagged timing to add-on library

	! Debug: print number of active openmp threads
	!$OMP PARALLEL
	!PRINT*, 'OMP_NUM_THREADS = ', OMP_GET_NUM_THREADS()
	!$OMP END PARALLEL

	! Set defaults
		! ParseInput.F90
	CALL SetDefaults()
	! Sort out command line arguments
	CALL ReadCommandLine()
	CALL AddTime(2)

	! Output some stats to stdout
	CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
	IF (verbose .AND. myid .EQ. 0) &
		WRITE(*,'(I3,A)') nproc," processors checking in."

	IF (verbose .AND. myid .EQ. 0) CALL PrintDetails()

	! Allocate space for the backframe
	ALLOCATE(backframe(doppix,doppix))
	backframe = 0D0

	! every processor read background file into memory
	IF (backgroundset .EQ. 1) THEN
		CALL Read_Background(myid, background, backframe, doppix)
	ENDIF

	! Set up tiling scheme, either from file or pure command-line
	IF (use_gridfile) THEN
		! when using grid file, can only do one tile size
		ij = 0
		DO ii=6,1,-1
			IF (dotilesize(ii)) THEN
				ij = ij + 1
				ind = ii
			ENDIF
		ENDDO
		IF (ij .NE. 1) THEN
			PRINT*, "ERROR: When using grid file, specify a single tile size using -ts#"
			STOP
		ENDIF

		! load grid file
		CALL ReadGrid_head(myid, gridfile, tilecount, numtiles, ind, verbose)
		ALLOCATE(tilesize(tilecount))
		ALLOCATE(lon(tilecount))
		ALLOCATE(lat(tilecount))
		ALLOCATE(delta_rot(tilecount))
		tilesize(:) = 2**(ind-1)
		CALL ReadGrid(myid, gridfile, tilecount, lon, lat, verbose)
		CALL SetTrackingRate(delta_rot,a0,a2,a4,lat,tilecount)
	ELSE
		! Count number of tiles needed for each tilesize
		CALL CountTiles(myid, tilecount, dotilesize, numtiles, lonrn, latrn, apode, verbose)
		! Allocate info foreach tile
		ALLOCATE(tilesize(tilecount))
		ALLOCATE(lon(tilecount))
		ALLOCATE(lat(tilecount))
		ALLOCATE(delta_rot(tilecount))
	
		! Go back through and load details of each tile
		currtile = 1
		DO ii=6,1,-1
			IF (dotilesize(ii)) THEN
				currsize = 2**(ii-1)
				tilesize(currtile:currtile+numtiles(ii)-1) = currsize
					! Grid.F90
				CALL GenerateGrid(currsize,lonrn,latrn,clon,clat,apode,&
					lon(currtile:currtile+numtiles(ii)-1),lat(currtile:currtile+numtiles(ii)-1))
					! Grid.F90
				CALL SetTrackingRate(delta_rot,a0,a2,a4,lat,numtiles(ii))
				currtile = currtile + numtiles(ii)
			ENDIF
		ENDDO
	ENDIF
	CALL AddTime(3)

	! Let proc0 read the master dop list
	CALL AddTime(4)
	IF (myid .EQ. 0) THEN
		IF (verbose) WRITE(*,'(A)') "Reading Master List.."
			! File_Management.F90
		CALL Load_Dopplergram_List(myid, nsteps, masterlist, dopfname, doptime, dopinterp)
	ENDIF
	CALL AddTime(5)

	CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

	! Communicate master dop list to all procs
	IF (myid .EQ. 0 .AND. verbose) &
		WRITE(*,'(A)') "Broadcasting dopplergram list to all procs.."
	CALL MPI_BCAST(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	IF (myid .NE. 0) THEN
		ALLOCATE(doptime(nsteps))
		ALLOCATE(dopfname(nsteps))
		ALLOCATE(dopinterp(nsteps))
	ENDIF
	ALLOCATE(dop_p_angle(nsteps))
	dop_p_angle = 0.0
	ALLOCATE(dop_cen_lon(nsteps))
	dop_cen_lon = 0.0
	ALLOCATE(dop_cen_lat(nsteps))
	dop_cen_lat = 0.0
	ALLOCATE(dop_cen_xpix(nsteps))
	dop_cen_xpix = 0.0
	ALLOCATE(dop_cen_ypix(nsteps))
	dop_cen_ypix = 0.0
	ALLOCATE(dop_r_sun_pix(nsteps))
	dop_r_sun_pix = 0.0
	CALL MPI_BCAST(doptime, nsteps, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!	CALL MPI_BCAST(doptime_ends, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	CALL MPI_BCAST(dopfname, nsteps*200, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
!	CALL MPI_BCAST(dopfname_ends, 400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
	CALL MPI_BCAST(dopinterp, nsteps, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
	CALL AddTime(6)

	! Can't interpolate first or last dopplergram
	IF (dopinterp(1) .OR. dopinterp(nsteps)) THEN
		WRITE(*,'(A)') "FATAL ERROR: Can't interpolate first or last dopplergram :("
		STOP
	ENDIF

	! Main loop through tile size
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (verbose .AND. myid .EQ. 0) PRINT*, "Entering Main Loop"
	currtile = 1
	DO ii=6,1,-1
		CALL AddTime(1000+ii)
		IF (dotilesize(ii)) THEN
			completed = 0 ! completed of this tilesize
			currsize = 2**(ii-1)

			IF (verbose .AND. myid .EQ. 0) WRITE(*,'(A,I0,A,I0,A)') &
					"Tracking tilesize ",currsize," (",numtiles(ii)," tiles)"

			! fill up each processor and run tracking
			! until no tiles are left
!			maxtiles = FLOOR((memlimittotal/nproc) &
!			/memusage(tilesize(currtile), nsteps, doppix))
			
			DO WHILE (completed .LT. numtiles(ii))
				
				! fill up each processor with tiles
				num = MIN(maxtiles*nproc, numtiles(ii)-completed)
				
				! call main tracking
					! Track.F90
				CALL TrackTiles(myid,nproc,num,currsize,nsteps,completed,&
					dopfname,doptime,dopfname_ends,doptime_ends,dopinterp,&
					loaddops,outdir,backframe,doppix,verbose,makepspec,&
					fname_guess,extended)

				completed = completed + num
				currtile = currtile + num

				CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
			
			ENDDO
		ENDIF
	ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL AddTime(7)

	! Free ALL the things
	DEALLOCATE(tilesize)
	DEALLOCATE(lon)
	DEALLOCATE(lat)
	DEALLOCATE(delta_rot)
	DEALLOCATE(doptime)
	DEALLOCATE(dopfname)
	DEALLOCATE(dopinterp)
	DEALLOCATE(dop_p_angle)
	DEALLOCATE(dop_cen_lon)
	DEALLOCATE(dop_cen_lat)
	DEALLOCATE(dop_cen_xpix)
	DEALLOCATE(dop_cen_ypix)
	DEALLOCATE(dop_r_sun_pix)
	DEALLOCATE(backframe)

	! End everything
	IF (verbose .AND. myid .EQ. 0) PRINT*, "Main loop done, Finalizing MPI.."
	CALL AddTime(8)
	IF (dotiming) THEN
		CALL SaveTimer(myid, nproc, "timing")
	ENDIF
	CALL KillTimer()
	endtime = MPI_Wtime()
	IF (verbose .AND. myid .EQ. 0) THEN
		WRITE(*,'(A,I0,A)') "Total Time = ",INT(endtime-starttime)," seconds"
	ENDIF
	CALL Finalize_Communication()

END PROGRAM
