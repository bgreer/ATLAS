
! This file handles things having to do with the gridding scheme
! it will decide how many tiles are needed, where they should be placed, etc

MODULE Grid
	USE ParseInput
	USE Communication_Library
	IMPLICIT NONE

CONTAINS

	! Parallel
	SUBROUTINE ReadGrid_head (myid, gridfile, tilecount, numtiles, ind, verbose)
		IMPLICIT NONE

		INTEGER :: myid, tilecount, ind, stat, nlon, nlat, ierr
		CHARACTER(LEN=*) :: gridfile
		INTEGER :: numtiles(6)
		LOGICAL :: verbose

		IF (myid .EQ. 0) THEN
			! read file
			OPEN(2,FILE=gridfile)
			READ(2,*,IOSTAT=stat) tilecount, nlon, nlat
			CLOSE(2)
			numtiles(ind) = tilecount
			IF (verbose) PRINT*, "    Tiles in Grid File: ",tilecount
		ENDIF
		CALL MPI_BCAST(tilecount, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		CALL MPI_BCAST(numtiles, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	END SUBROUTINE ReadGrid_head

	SUBROUTINE ReadGrid (myid, gridfile, tilecount, lon, lat, verbose)
		IMPLICIT NONE

		INTEGER :: myid, tilecount, ii, ierr, stat
		INTEGER :: temp1, temp2, temp3
		REAL :: lon(tilecount), lat(tilecount)
		CHARACTER(LEN=*) :: gridfile
		LOGICAL :: verbose

		IF (myid .EQ. 0) THEN
			OPEN(2,FILE=gridfile)
			! don't need the header right now
			READ(2,*,IOSTAT=stat) temp1, temp2, temp3
			DO ii=1,tilecount
				READ(2,*,IOSTAT=stat) lon(ii), lat(ii)
			ENDDO
			CLOSE(2)
		ENDIF

		CALL MPI_BCAST(lon, tilecount, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
		CALL MPI_BCAST(lat, tilecount, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

	END SUBROUTINE ReadGrid

	! Generate central tile coords for all tiles in a region
	! INPUT: ts, lonrn, latrn, clon, clat, apode
	! OUTPUT: lon(), lat()  assumed to be allocated with correct size
	! Only one kind of grid supported, dense-pack
	SUBROUTINE GenerateGrid (ts, lonrn, latrn, clon, clat, apode, lon, lat)
		IMPLICIT NONE

		INTEGER :: ts, numx, numy, ii, ij, currtile
		REAL :: lonrn, latrn, clon, clat, apode, space
		REAL :: lonmin, latmin
		REAL :: lon(:), lat(:)

		space = ts*densepack*apode

		numx = INT(lonrn/space+1)
		numy = INT(latrn/space+1)

		lonmin = clon - 0.5*lonrn
		latmin = clat - 0.5*latrn

!		IF (verbose .AND. myid .EQ. 0) WRITE(*,'(A,F10.5)') "Generating Grid with spacing: ", space
		currtile = 1
		DO ii=1,numx
			DO ij=1,numy
				lon(currtile) = lonmin+(ii-1)*space
				lat(currtile) = latmin+(ij-1)*space
				
				! Wrap longitude around
				IF (lon(currtile) .GT. 360D0) &
					lon(currtile) = lon(currtile) - 360D0
				IF (lon(currtile) .LT. 0D0) &
					lon(currtile) = lon(currtile) + 360D0
				currtile = currtile + 1
			ENDDO
		ENDDO

	END SUBROUTINE GenerateGrid

	! Set delta_rot for each tile
	! INPUT: a0,a2,a4, lat(num), num
	! OUTPUT: delta_rot(num)
	SUBROUTINE SetTrackingRate (delta_rot, a0, a2, a4, lat, num)
		IMPLICIT NONE
		INTEGER :: num, ii
		REAL :: delta_rot(num), lat(num)
		REAL :: sinlat, a0, a2, a4

		DO ii=1,num
			sinlat = sin(lat(ii)*0.01745329251D0) ! forgot to convert to radians earlier... oops..
			sinlat = sinlat * sinlat
			! units of degrees per second
			delta_rot(ii) = 57.2957795131D0*(a0 + sinlat * (a2 + a4*sinlat))*1E-6
		ENDDO
	END SUBROUTINE SetTrackingRate

	! Given a lat/lon range and tilesize, give number of tiles
	INTEGER FUNCTION GetTileCount (ts, lonrn, latrn, apode)
		IMPLICIT NONE
		INTEGER :: ts
		REAL :: lonrn, latrn, space, apode
		space = ts*densepack*apode
		GetTileCount = INT(lonrn/space+1)*INT(latrn/space+1)
	END FUNCTION

	
	SUBROUTINE CountTiles(myid, tilecount, dotilesize, numtiles, lonrn, latrn, apode, verbose)
		IMPLICIT NONE
		INTEGER :: tilecount, myid, currsize, ii
		INTEGER :: numtiles(6)
		LOGICAL :: verbose, dotilesize(6)
		REAL :: lonrn, latrn, apode

		tilecount = 0
		DO ii=6,1,-1
			IF (dotilesize(ii)) THEN
				currsize = 2**(ii-1)
				numtiles(ii) = GetTileCount(currsize,lonrn,latrn,apode)
				tilecount = tilecount + numtiles(ii)
				IF (verbose .AND. myid .EQ. 0) &
					WRITE(*,'(A,I6,A,I2,A)') "Adding ",numtiles(ii)," ", &
						currsize," degree tiles"
			ENDIF
		ENDDO
	END SUBROUTINE CountTiles

END MODULE
