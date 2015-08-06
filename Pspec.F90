
! Taken from PSPEC code
! input a tracked data cube, make a power spectrum data cube, send it off for
! fitting

MODULE Pspec
	USE Interpolation
	USE File_Management
	USE Fit

	IMPLICIT NONE
	
	INCLUDE "fftw3.f"


CONTAINS

	SUBROUTINE Make_Pspec(myid, arr, dim1, dim2, dim3, verbose, outfile, fname_guess, outdir, lon, lat)
		IMPLICIT NONE
		
		LOGICAL :: verbose
		INTEGER :: dim1, dim2, dim3, myid
		REAL, DIMENSION(:,:,:) :: arr
		CHARACTER(LEN=*) :: outfile, outdir
		CHARACTER(LEN=*) :: fname_guess
		REAL :: lon, lat

		REAL, DIMENSION(:,:), ALLOCATABLE :: mask2d
		REAL, DIMENSION(:), ALLOCATABLE :: mask1d
		REAL :: apod_min, apod_min_t, r, cx, cy, x, y, x0, y0, dtheta, tht, pi
		INTEGER :: ii, ij, ik
		INTEGER*8 :: plan
		DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: inarr2
		DOUBLE COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: outarr, outarr2
		REAL, DIMENSION(:,:,:), ALLOCATABLE :: power, power_shift
		INTEGER :: ntheta, nk, ntheta_sub, reducefactor
		DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: power_cart, power_sub
		DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: linein, lineout, subin, subout
		REAL :: val

		pi = ACOS(-1D0)

		! Background subtraction
!$OMP PARALLEL DO PRIVATE(ii)
		DO ii=1,dim3
			arr(:,:,ii) = arr(:,:,ii) - SUM(arr(:,:,ii)) / (float(dim1)*float(dim2))
		ENDDO
!$OMP END PARALLEL DO

		! Create Apodization Masks
		ALLOCATE(mask2d(dim1,dim2))
		ALLOCATE(mask1d(dim3))
		mask1d = 1.0
		mask2d = 1.0
		apod_min = 0.9375
		apod_min_t = 0.96875
		cx = dim1/2.+0.5
		cy = dim2/2.+0.5
		DO ii=1,dim3 ! Time
			IF (FLOAT(ii-1)/FLOAT(dim3) .LT. 1.0-apod_min_t) THEN
				mask1d(ii) = (1.-(1.-(FLOAT(ii-1)/FLOAT(dim3))/(1.0-apod_min_t))**2.)**2.
			ENDIF
			IF (FLOAT(ii-1)/FLOAT(dim3) .GT. apod_min_t) THEN
				mask1d(ii) = (1.-(1.-(FLOAT(dim3-ii)/FLOAT(dim3))/(1.-apod_min_t))**2.)**2.
			ENDIF
		ENDDO
		DO ii=1,dim2 ! Space
			DO ij=1,dim1
				r = MIN(SQRT((ij-cx)**2.+(ii-cy)**2.) / cx, 1.0)
				IF (r .GT. apod_min) THEN
					mask2d(ii,ij) = (1.-((r-apod_min)/(1.-apod_min))**2.)**2.
				ENDIF
			ENDDO
		ENDDO

		! Apply Apodization to Data
!$OMP PARALLEL DO PRIVATE(ii)
		DO ii=1,dim3
			arr(:,:,ii) = arr(:,:,ii)*mask1d(ii)*mask2d
		ENDDO
!$OMP END PARALLEL DO

		! Free masks
		DEALLOCATE(mask2d)
		DEALLOCATE(mask1d)

		! transform because i'm lazy
		ALLOCATE(inarr2(dim3, dim2, dim1))
		! copy inarr to inarr2
!$OMP PARALLEL DO PRIVATE(ii,ij)
		DO ii=1,dim3
			DO ij=1,dim1
				inarr2(ii,:,ij) = arr(ij,:,ii)
			ENDDO
		ENDDO
!$OMP END PARALLEL DO

		ALLOCATE(outarr2(dim3/2+1, dim2, dim1))

		plan = 0
		! Set up FFT of inarr->outarr
		CALL DFFTW_PLAN_DFT_R2C_3D(plan,dim3,dim2,dim1,inarr2,outarr2,FFTW_ESTIMATE)
		CALL DFFTW_EXECUTE(plan)
		CALL DFFTW_DESTROY_PLAN(plan)
		CALL DFFTW_CLEANUP()

		DEALLOCATE(inarr2)

		! copy outarr2 to outarr
		ALLOCATE(outarr(dim1,dim2,dim3/2+1))
!$OMP PARALLEL DO PRIVATE(ii,ij)
		DO ii=1,dim3/2+1
			DO ij=1,dim1
				outarr(ij,:,ii) = outarr2(ii,:,ij)
			ENDDO
		ENDDO
!$OMP END PARALLEL DO
		DEALLOCATE(outarr2)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! POWER SPECTRUM
		ALLOCATE(power(dim1,dim2,dim3/2+1))
		! Process result
!$OMP PARALLEL DO PRIVATE(ii)
		DO ii=1,dim3/2+1
			power(:,:,ii) = 0.5*LOG(REAL(outarr(:,:,ii))**2. + &
				AIMAG(outarr(:,:,ii))**2.)
		ENDDO
!$OMP END PARALLEL DO
		DEALLOCATE(outarr)


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CART TO POLAR

		! Shift pixels
		ALLOCATE(power_shift(dim1,dim2,dim3/2+1))
		power_shift(dim1/2+1:dim1,:,:) = power(1:dim1/2,:,:)
		power_shift(1:dim1/2,:,:) = power(dim1/2+1:dim1,:,:)
		power(:,dim2/2+1:dim2,:) = power_shift(:,1:dim2/2,:)
		power(:,1:dim2/2,:) = power_shift(:,dim2/2+1:dim2,:)
		DEALLOCATE(power_shift)

		nk = dim1/2
		ntheta = 200
		SELECT CASE (dim1)
			CASE (768) ! 32 degrees
				ntheta = 512
			CASE (384) ! 16 degree
				ntheta = 256
			CASE (192) ! 8 degree
				ntheta = 128
			CASE (96) ! 4 degree
				ntheta = 64
			CASE (48) ! 2 degree
				ntheta = 32
		END SELECT
		ALLOCATE(power_cart(ntheta,nk,dim3/2))
		x0 = FLOAT(dim1/2)+1.
		y0 = FLOAT(dim2/2)+1.
		dtheta = 2.0*pi/FLOAT(ntheta)
!$OMP PARALLEL DO PRIVATE(ii,ij,ik,tht,x,y)
		DO ii=1,dim3/2
			DO ij=1,nk
				DO ik=1,ntheta
					tht = (ik-1)*dtheta
					x = FLOAT(ij)*COS(tht)+x0
					y = FLOAT(ij)*SIN(tht)+y0
					CALL INTERP(x,y,power(:,:,ii),dim1,val)
					power_cart(ik,ij,ii) = EXP(val)
				ENDDO
			ENDDO
		ENDDO
!$OMP END PARALLEL DO
		DEALLOCATE(power)

		ntheta_sub = INT(MAX(ntheta/4,16))
		reducefactor = INT(ntheta/ntheta_sub)
		ALLOCATE(power_sub(ntheta_sub,nk,dim3/2))
		ALLOCATE(linein(ntheta))
		ALLOCATE(lineout(ntheta))
		ALLOCATE(subin(ntheta))
		ALLOCATE(subout(ntheta))

		DO ii=1,dim3/2
			DO ij=1,nk
				! copy spectrum -> linein
				linein = CMPLX(power_cart(:,ij,ii),0.0)
				! ft linein -> lineout
				CALL DFFTW_PLAN_DFT_1D(plan,ntheta,linein,lineout,FFTW_FORWARD,FFTW_ESTIMATE)
				CALL DFFTW_EXECUTE_DFT(plan, linein, lineout)
				CALL DFFTW_DESTROY_PLAN(plan)
				! copy lineout -> subin
				subin = lineout
				subin(ntheta_sub+1:ntheta-ntheta_sub+1) = CMPLX(0.0,0.0)
				! ft subin -> subout
				CALL DFFTW_PLAN_DFT_1D(plan,ntheta,subin,subout,FFTW_BACKWARD,FFTW_ESTIMATE)
				CALL DFFTW_EXECUTE_DFT(plan, subin, subout)
				CALL DFFTW_DESTROY_PLAN(plan)
				! copy subout > spectrum
				DO ik=1,ntheta_sub
					power_sub(ik,ij,ii) = REAL(subout((ik-1)*reducefactor+1))
				ENDDO
			ENDDO
		ENDDO
		DEALLOCATE(linein)
		DEALLOCATE(lineout)
		DEALLOCATE(subin)
		DEALLOCATE(subout)
		DEALLOCATE(power_cart)

		CALL FitPspec(myid, lon, lat, power_sub, ntheta_sub, nk, dim3/2, 360D0/(384*696D0/24D0), 1D6/(45D0*dim3), fname_guess, verbose, outdir)
		!CALL Save_Tile(power_sub(:,:,:), ntheta_sub, nk, dim3/2, dim3, TRIM(outfile), verbose)
		DEALLOCATE(power_sub)
	END SUBROUTINE Make_Pspec


END MODULE
