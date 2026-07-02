MODULE AirfoilInfo_RotCor

   USE                                             AirfoilInfo_Types
   USE                                          :: NWTC_LAPACK
   USE                                          :: NWTC_Library_Types

   IMPLICIT NONE

   PRIVATE

   PUBLIC                                       :: AFI_PreCalcRotCorTables
   PUBLIC                                       :: AFI_ComputeAirfoilCoefsRotCor1D
   PUBLIC                                       :: AFI_ComputeAirfoilCoefsRotCor2D
   PUBLIC                                       :: AFI_ComputeUACoefsRotCor1D
   PUBLIC                                       :: AFI_ComputeUACoefsRotCor2D
   PUBLIC                                       :: AFI_CalcSnel

   integer, parameter                           :: MaxNumAFCoeffs = 7 !cl,cd,cm,cpMin, UA:f_st, FullySeparate, FullyAttached

   ABSTRACT INTERFACE
      SUBROUTINE CalculateUACoeffs_Proc(CalcDefaults,p,ColCl,ColCd,ColCm,ColUAf,UAMod)
         USE AirfoilInfo_Types
         TYPE (AFI_UA_BL_Default_Type),INTENT(IN):: CalcDefaults
         TYPE (AFI_Table_Type),    INTENT(INOUT) :: p
         INTEGER(IntKi),           INTENT(IN   ) :: ColCl
         INTEGER(IntKi),           INTENT(IN   ) :: ColCd
         INTEGER(IntKi),           INTENT(IN   ) :: ColCm
         INTEGER(IntKi),           INTENT(IN   ) :: ColUAf
         INTEGER(IntKi),           INTENT(IN   ) :: UAMod
      END SUBROUTINE CalculateUACoeffs_Proc
   END INTERFACE

CONTAINS

SUBROUTINE AFI_PreCalcRotCorTables(p, iTable, InitInp, ErrStat, ErrMsg, RoutineName, CalculateUACoeffsCB)
    ! Description:
    ! This subroutine pre-calculates the Snel-corrected airfoil tables for all
    ! Snel factors and stores them in the rotCorTables component of the
    ! derived type. This helps avoid runtime calculations when the Snel effect is active.
    
    ! Arguments:
    IMPLICIT NONE
    
	TYPE (AFI_ParameterType), INTENT(INOUT)   :: p            ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
    INTEGER(IntKi),          INTENT(IN)       :: iTable       ! Index of the table to process
    TYPE (AFI_InitInputType), INTENT(IN)      :: InitInp                       ! This structure stores values that are set by the calling routine during the initialization phase.
    INTEGER(IntKi),          INTENT(OUT)      :: ErrStat      ! Error status flag
    CHARACTER(*),            INTENT(OUT)      :: ErrMsg       ! Error message string
    CHARACTER(*),            INTENT(IN)       :: RoutineName  ! Name of the calling routine for error handling
   PROCEDURE(CalculateUACoeffs_Proc)         :: CalculateUACoeffsCB
   TYPE (AFI_UA_BL_Default_Type) :: CalcDefaults            ! Whether to calculate default UA params

    ! Local Variables
    INTEGER(IntKi)            :: rotCorTabIdx      ! Loop index for rotational correction tables
    REAL(ReKi)                :: snel_factor, schepers_factor  ! The correction factors
    REAL(ReKi),    ALLOCATABLE:: Cl_vec(:),Cd_vec(:)    ! Local copy of Cl for interpolation
    REAL(ReKi)                :: Cl_0         ! Lift coefficient at zero angle of attack
	REAL(ReKi)                :: max_snel_factor         ! Maximum snel_factor to tabulate
	INTEGER(IntKi)            :: num_rotCor_tables         ! Number of tables to create. 
    INTEGER(IntKi)            :: iLo          ! Lower index for interpolation
    INTEGER(IntKi)            :: ErrStat2     ! Local error status for splines
    CHARACTER(300)            :: ErrMsg2      ! Local error message for splines
    
    TYPE(AFI_Table_Type)      :: tempRotCorTable ! Temporary table for deep copy operations

    ! Initialize error handling
    ErrStat = ErrID_None
    ErrMsg  = ""	
	
    ! --- Main Logic ---
	
	max_snel_factor = 1.0_ReKi
	num_rotCor_tables = 200

   if ( p%Table(iTable)%ConstData ) then
      ! Constant tables are identical for all snel factors, so store only one table.
      if (allocated(p%Table(iTable)%rotCorTables)) deallocate(p%Table(iTable)%rotCorTables)
      allocate(p%Table(iTable)%rotCorTables(1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, 'Error allocating RotCor tables.', ErrStat, ErrMsg, RoutineName)
         return
      end if

      tempRotCorTable%UserProp   = p%Table(iTable)%UserProp
      tempRotCorTable%Re         = p%Table(iTable)%Re
      tempRotCorTable%NumAlf     = p%Table(iTable)%NumAlf
      tempRotCorTable%ConstData  = p%Table(iTable)%ConstData
      tempRotCorTable%InclUAdata = p%Table(iTable)%InclUAdata
      tempRotCorTable%UA_BL      = p%Table(iTable)%UA_BL

      if ( allocated(p%Table(iTable)%Alpha) ) then
         if (allocated(tempRotCorTable%Alpha)) deallocate(tempRotCorTable%Alpha)
         allocate(tempRotCorTable%Alpha(size(p%Table(iTable)%Alpha)))
         tempRotCorTable%Alpha = p%Table(iTable)%Alpha
      end if
      if ( allocated(p%Table(iTable)%Coefs) ) then
         if (allocated(tempRotCorTable%Coefs)) deallocate(tempRotCorTable%Coefs)
         allocate(tempRotCorTable%Coefs(size(p%Table(iTable)%Coefs,1), size(p%Table(iTable)%Coefs,2)))
         tempRotCorTable%Coefs = p%Table(iTable)%Coefs
      end if

      p%Table(iTable)%rotCorTables(1)%snel_factor = 0.0_ReKi
      p%Table(iTable)%rotCorTables(1)%UserProp    = tempRotCorTable%UserProp
      p%Table(iTable)%rotCorTables(1)%Re          = tempRotCorTable%Re
      p%Table(iTable)%rotCorTables(1)%NumAlf      = tempRotCorTable%NumAlf
      p%Table(iTable)%rotCorTables(1)%ConstData   = tempRotCorTable%ConstData
      p%Table(iTable)%rotCorTables(1)%InclUAdata  = tempRotCorTable%InclUAdata
      p%Table(iTable)%rotCorTables(1)%UA_BL       = tempRotCorTable%UA_BL

      if (allocated(p%Table(iTable)%rotCorTables(1)%Alpha)) deallocate(p%Table(iTable)%rotCorTables(1)%Alpha)
      allocate(p%Table(iTable)%rotCorTables(1)%Alpha(size(tempRotCorTable%Alpha)))
      p%Table(iTable)%rotCorTables(1)%Alpha = tempRotCorTable%Alpha

      if (allocated(p%Table(iTable)%rotCorTables(1)%Coefs)) deallocate(p%Table(iTable)%rotCorTables(1)%Coefs)
      allocate(p%Table(iTable)%rotCorTables(1)%Coefs(size(tempRotCorTable%Coefs,1), size(tempRotCorTable%Coefs,2)))
      p%Table(iTable)%rotCorTables(1)%Coefs = tempRotCorTable%Coefs

      if (allocated(p%Table(iTable)%rotCorTables(1)%SplineCoefs)) deallocate(p%Table(iTable)%rotCorTables(1)%SplineCoefs)

      if (allocated(tempRotCorTable%Alpha)) deallocate(tempRotCorTable%Alpha)
      if (allocated(tempRotCorTable%Coefs)) deallocate(tempRotCorTable%Coefs)
      if (allocated(tempRotCorTable%SplineCoefs)) deallocate(tempRotCorTable%SplineCoefs)
      return
   end if

    ! The rotCorTables array must be allocated before it can be used.
   if (allocated(p%Table(iTable)%rotCorTables)) deallocate(p%Table(iTable)%rotCorTables)
    allocate(p%Table(iTable)%rotCorTables(num_rotCor_tables+1), STAT=ErrStat2)
    if (ErrStat2 /= 0) then
       call SetErrStat(ErrID_Fatal, 'Error allocating RotCor tables.', ErrStat, ErrMsg, RoutineName)
       return
    end if
    
    DO rotCorTabIdx = 1, num_rotCor_tables+1 ! Loop over snel factor values
        ! This formula creates evenly spaced values from 0.0 to max_snel_factor.
        snel_factor = ( (rotCorTabIdx-1.0_ReKi)/(num_rotCor_tables) ) * max_snel_factor
        
      ! Copy only the source table fields needed for correction.
      tempRotCorTable%UserProp   = p%Table(iTable)%UserProp
      tempRotCorTable%Re         = p%Table(iTable)%Re
      tempRotCorTable%NumAlf     = p%Table(iTable)%NumAlf
      tempRotCorTable%ConstData  = p%Table(iTable)%ConstData
      tempRotCorTable%InclUAdata = p%Table(iTable)%InclUAdata
      tempRotCorTable%UA_BL      = p%Table(iTable)%UA_BL
        
        ! Check if allocatable arrays exist before attempting to allocate and copy.
        if ( allocated(p%Table(iTable)%Alpha) ) then
            if (allocated(tempRotCorTable%Alpha)) deallocate(tempRotCorTable%Alpha)
            allocate(tempRotCorTable%Alpha(size(p%Table(iTable)%Alpha)))
            tempRotCorTable%Alpha = p%Table(iTable)%Alpha
        end if
        if ( allocated(p%Table(iTable)%Coefs) ) then
            if (allocated(tempRotCorTable%Coefs)) deallocate(tempRotCorTable%Coefs)
            allocate(tempRotCorTable%Coefs(size(p%Table(iTable)%Coefs,1), size(p%Table(iTable)%Coefs,2)))
            tempRotCorTable%Coefs = p%Table(iTable)%Coefs
        end if
        
        if ( p%Table(iTable)%ConstData ) then
            ! For constant tables, copy the entire table payload directly.
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%snel_factor = snel_factor
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%UserProp    = tempRotCorTable%UserProp
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%Re          = tempRotCorTable%Re
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%NumAlf      = tempRotCorTable%NumAlf
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%ConstData   = tempRotCorTable%ConstData
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%InclUAdata  = tempRotCorTable%InclUAdata
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%UA_BL       = tempRotCorTable%UA_BL

            if (allocated(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha)) deallocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha)
            allocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha(size(tempRotCorTable%Alpha)))
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha = tempRotCorTable%Alpha

            if (allocated(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs)) deallocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs)
            allocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs(size(tempRotCorTable%Coefs,1), size(tempRotCorTable%Coefs,2)))
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs = tempRotCorTable%Coefs

            if (allocated(p%Table(iTable)%rotCorTables(rotCorTabIdx)%SplineCoefs)) deallocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%SplineCoefs)
        else
            ! The Cl_vec needs to be a local copy to avoid a cumulative error.
            Cl_vec = tempRotCorTable%Coefs(:, p%ColCl)
			Cd_vec = tempRotCorTable%Coefs(:, p%ColCd)
            
            ! Find Cl at alpha=0 using InterpBinReal
            iLo = 0
            Cl_0 = InterpBinReal( 0.0_ReKi, tempRotCorTable%Alpha, Cl_vec, iLo, size(tempRotCorTable%Alpha) )
            
			if (p%RotCorParams%RotCor==1) then ! apply Snel Cl correction
               call AFI_ApplySnel(tempRotCorTable%Alpha, Cl_vec, snel_factor, Cl_0, tempRotCorTable%UA_BL%C_lalpha)
			else if (p%RotCorParams%RotCor==2) then ! apply Schepers corrections
			   schepers_factor =  1.72*snel_factor**0.7 ! 3.8_ReKi * (chord_over_r**1.4) * (cos(twist+pitch))**6			   
			   call AFI_ApplySchepers(tempRotCorTable%Alpha, Cl_vec, Cd_vec, schepers_factor, tempRotCorTable%UA_BL%C_lalpha, Cl_0)
			else if (p%RotCorParams%RotCor==3) then ! apply Snel Cl and Gertz Cd correction
			   call AFI_ApplySnel(tempRotCorTable%Alpha, Cl_vec, snel_factor, Cl_0, tempRotCorTable%UA_BL%C_lalpha)			   
			   call AFI_ApplyGertzCd(tempRotCorTable%Alpha, Cl_vec, Cd_vec, snel_factor)			   
			end if

			! Copy the corrected values back into the tempRotCorTable%Coefs array
            tempRotCorTable%Coefs(:, p%ColCl) = Cl_vec			
			tempRotCorTable%Coefs(:, p%ColCd) = Cd_vec
			
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%snel_factor = snel_factor
            
            ! Explicitly allocate the allocatable components before assigning to them.
            if (allocated(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha)) deallocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha)
            allocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha(size(tempRotCorTable%Alpha)))
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha = tempRotCorTable%Alpha
            
            if (allocated(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs)) deallocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs)
            allocate(p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs(size(tempRotCorTable%Coefs,1), size(tempRotCorTable%Coefs,2)))

   	        ! Allocate spline coefficients on the temporary corrected table.
            if (allocated(tempRotCorTable%SplineCoefs)) deallocate(tempRotCorTable%SplineCoefs)
            allocate(tempRotCorTable%SplineCoefs( p%Table(iTable)%NumAlf-1, size(tempRotCorTable%Coefs,2), 0:3 ), STAT=ErrStat2)
            if (ErrStat2 /= 0) then
               call SetErrStat(ErrID_Fatal, 'Error allocating temporary RotCor spline coefficients.', ErrStat, ErrMsg, RoutineName)
               return
            end if
			
            if (InitInp%UAMod>0) then			
                ! For UA_Flag, we pass the local temporary table to the function.
               call CalculateUACoeffsCB(CalcDefaults, tempRotCorTable, p%ColCl, p%ColCd, p%ColCm, p%ColUAf, InitInp%UAMod)
            end if
            
            if ( p%InterpOrd == 3_IntKi ) then
                call CubicSplineInitM(tempRotCorTable%Alpha, tempRotCorTable%Coefs, tempRotCorTable%SplineCoefs, ErrStat2, ErrMsg2)
                call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            else if ( p%InterpOrd == 1_IntKi ) then
                call CubicLinSplineInitM(tempRotCorTable%Alpha, tempRotCorTable%Coefs, tempRotCorTable%SplineCoefs, ErrStat2, ErrMsg2)
                call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            end if
            
            ! Now copy everything to the final table
			p%Table(iTable)%rotCorTables(rotCorTabIdx)%Alpha = tempRotCorTable%Alpha
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%Coefs = tempRotCorTable%Coefs			
			p%Table(iTable)%rotCorTables(rotCorTabIdx)%SplineCoefs = tempRotCorTable%SplineCoefs
            
            ! Copy all other relevant properties
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%UserProp = tempRotCorTable%UserProp
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%Re = tempRotCorTable%Re
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%NumAlf = tempRotCorTable%NumAlf
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%ConstData = tempRotCorTable%ConstData
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%InclUAdata = tempRotCorTable%InclUAdata
            p%Table(iTable)%rotCorTables(rotCorTabIdx)%UA_BL = tempRotCorTable%UA_BL			
			
        end if
    ENDDO ! rotCorTabIdx loop
    ! Deallocate temporary arrays
    if (allocated(Cl_vec)) deallocate(Cl_vec)
   if (allocated(Cd_vec)) deallocate(Cd_vec)
    if (allocated(tempRotCorTable%Alpha)) deallocate(tempRotCorTable%Alpha)
    if (allocated(tempRotCorTable%Coefs)) deallocate(tempRotCorTable%Coefs)
    if (allocated(tempRotCorTable%SplineCoefs)) deallocate(tempRotCorTable%SplineCoefs)

END SUBROUTINE AFI_PreCalcRotCorTables    
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine computes airfoil coefficients from a single rotationally corrected table
subroutine AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, RotCorTable, p, AFI_interp, errStat, errMsg )

   real(ReKi),                             intent(in   ) :: AOA
   type(RotCor_RotCorTableType),           intent(in   ) :: RotCorTable
   TYPE (AFI_ParameterType),               intent(in   ) :: p
   type(AFI_OutputType),                   intent(  out) :: AFI_interp
   integer(IntKi),                         intent(  out) :: errStat
   character(*),                           intent(  out) :: errMsg
   
   ! Local variables
   real                                    :: IntAFCoefs(MaxNumAFCoeffs)
   real(ReKi)                              :: Alpha
   integer                                 :: s1
   
   errStat = ErrID_None
   errMsg  = ""

   if (.not. allocated(RotCorTable%Coefs)) then
      errStat = ErrID_Fatal
      errMsg  = 'RotCor table coefficients are not allocated.'
      return
   end if
   
   IntAFCoefs = 0.0_ReKi ! initialize
   s1 = size(RotCorTable%Coefs,2)
   
   if (RotCorTable%ConstData) then
      IntAFCoefs(1:s1) = RotCorTable%Coefs(1,:)   ! all the rows are constant
   else
      if (.not. allocated(RotCorTable%Alpha)) then
         errStat = ErrID_Fatal
         errMsg  = 'RotCor table alpha values are not allocated.'
         return
      end if
      if (.not. allocated(RotCorTable%SplineCoefs)) then
         errStat = ErrID_Fatal
         errMsg  = 'RotCor table spline coefficients are not allocated.'
         return
      end if
      Alpha = AOA
      call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
      
      ! Spline interpolation based on requested AOA
      CALL CubicSplineInterpM( Alpha, RotCorTable%Alpha, RotCorTable%Coefs, &
                               RotCorTable%SplineCoefs, IntAFCoefs(1:s1) )
						 
   end if
  
   AFI_interp%Cl    = IntAFCoefs(p%ColCl)
   AFI_interp%Cd    = IntAFCoefs(p%ColCd)

   if ( p%ColCm > 0 ) then
      AFI_interp%Cm = IntAFCoefs(p%ColCm)
   else
      AFI_interp%Cm    = 0.0_ReKi
   end if
   
   if ( p%ColCpmin > 0 ) then
      AFI_interp%Cpmin = IntAFCoefs(p%ColCpmin)
   else
      AFI_interp%Cpmin = 0.0_ReKi
   end if

   if ( p%ColUAf > 0 ) then
      AFI_interp%f_st          = IntAFCoefs(p%ColUAf)
      AFI_interp%fullySeparate = IntAFCoefs(p%ColUAf+1)
      AFI_interp%fullyAttached = IntAFCoefs(p%ColUAf+2)
   else
      AFI_interp%f_st          = 0.0_ReKi
      AFI_interp%fullySeparate = 0.0_ReKi
      AFI_interp%fullyAttached = 0.0_ReKi
   end if
   
   ! needed if using UnsteadyAero:
   if (RotCorTable%InclUAdata) then
      AFI_interp%Cd0 = RotCorTable%UA_BL%Cd0
      AFI_interp%Cm0 = RotCorTable%UA_BL%Cm0
   else
      AFI_interp%Cd0 = 0.0_ReKi
      AFI_interp%Cm0 = 0.0_ReKi
   end if

end subroutine AFI_ComputeAirfoilCoefsFromRotCorTable
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine handles 1D airfoil coefficient calculation with rotation correction
subroutine AFI_ComputeAirfoilCoefsRotCor1D( AOA, p, AFI_interp, errStat, errMsg, iTable )

   real(ReKi),                             intent(in   ) :: AOA
   TYPE (AFI_ParameterType),               intent(in   ) :: p
   type(AFI_OutputType),                   intent(  out) :: AFI_interp
   integer(IntKi),                         intent(  out) :: errStat
   character(*),                           intent(  out) :: errMsg
   integer(IntKi),                         intent(in   ) :: iTable
   
   ! Local variables
   integer(IntKi)                          :: rotCorTabIdxLo, rotCorTabIdxHi
   real(ReKi)                              :: rotCorInterpFrac
   type(AFI_OutputType)                    :: AFI_interpLo, AFI_interpHi
   integer(IntKi)                          :: ErrStat2
   character(300)                          :: ErrMsg2
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Handle constant data case - similar to AFI_ComputeAirfoilCoefs1D
   if (p%Table(iTable)%ConstData) then
      ! For constant data, use the precomputed RotCor table directly.
      if (allocated(p%Table(iTable)%rotCorTables) .and. size(p%Table(iTable)%rotCorTables) > 0) then
         call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(iTable)%rotCorTables(1), p, AFI_interp, errStat, errMsg )
      else
         ! Defensive fallback in case precomputed tables are unavailable.
         AFI_interp%Cl = p%Table(iTable)%Coefs(1,p%ColCl)
         AFI_interp%Cd = p%Table(iTable)%Coefs(1,p%ColCd)

         if ( p%ColCm > 0 ) then
            AFI_interp%Cm = p%Table(iTable)%Coefs(1,p%ColCm)
         else
            AFI_interp%Cm = 0.0_ReKi
         end if

         if ( p%ColCpmin > 0 ) then
            AFI_interp%Cpmin = p%Table(iTable)%Coefs(1,p%ColCpmin)
         else
            AFI_interp%Cpmin = 0.0_ReKi
         end if

         if ( p%ColUAf > 0 ) then
            AFI_interp%f_st          = p%Table(iTable)%Coefs(1,p%ColUAf)
            AFI_interp%fullySeparate = p%Table(iTable)%Coefs(1,p%ColUAf+1)
            AFI_interp%fullyAttached = p%Table(iTable)%Coefs(1,p%ColUAf+2)
         else
            AFI_interp%f_st          = 0.0_ReKi
            AFI_interp%fullySeparate = 0.0_ReKi
            AFI_interp%fullyAttached = 0.0_ReKi
         end if

         if (p%Table(iTable)%InclUAdata) then
            AFI_interp%Cd0 = p%Table(iTable)%UA_BL%Cd0
            AFI_interp%Cm0 = p%Table(iTable)%UA_BL%Cm0
         else
            AFI_interp%Cd0 = 0.0_ReKi
            AFI_interp%Cm0 = 0.0_ReKi
         end if
      end if
      return
   end if
   
   ! Find the table indices for interpolation
   call FindRotCorTableIndices(p%RotCorParams%current_snel_factor, p%Table(iTable), &
                             rotCorTabIdxLo, rotCorTabIdxHi, rotCorInterpFrac, errStat2, errMsg2)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeAirfoilCoefsRotCor1D')
   if (errStat >= AbortErrLev) return
   
   ! Get coefficients from lower table
   call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(iTable)%rotCorTables(rotCorTabIdxLo), p, AFI_interpLo, errStat2, errMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeAirfoilCoefsRotCor1D')
   if (errStat >= AbortErrLev) return
   
   if (rotCorTabIdxLo == rotCorTabIdxHi) then
      ! No interpolation needed
      AFI_interp = AFI_interpLo
   else
      ! Get coefficients from higher table
      call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(iTable)%rotCorTables(rotCorTabIdxHi), p, AFI_interpHi, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeAirfoilCoefsRotCor1D')
      if (errStat >= AbortErrLev) return
      
      ! Interpolate between the two results
      AFI_interp%Cl = AFI_interpLo%Cl + rotCorInterpFrac * (AFI_interpHi%Cl - AFI_interpLo%Cl)
      AFI_interp%Cd = AFI_interpLo%Cd + rotCorInterpFrac * (AFI_interpHi%Cd - AFI_interpLo%Cd)
      AFI_interp%Cm = AFI_interpLo%Cm + rotCorInterpFrac * (AFI_interpHi%Cm - AFI_interpLo%Cm)
      AFI_interp%Cd0 = AFI_interpLo%Cd0 + rotCorInterpFrac * (AFI_interpHi%Cd0 - AFI_interpLo%Cd0)
      AFI_interp%Cm0 = AFI_interpLo%Cm0 + rotCorInterpFrac * (AFI_interpHi%Cm0 - AFI_interpLo%Cm0)
      AFI_interp%Cpmin = AFI_interpLo%Cpmin + rotCorInterpFrac * (AFI_interpHi%Cpmin - AFI_interpLo%Cpmin)
      AFI_interp%f_st = AFI_interpLo%f_st + rotCorInterpFrac * (AFI_interpHi%f_st - AFI_interpLo%f_st)
      AFI_interp%FullySeparate = AFI_interpLo%FullySeparate + rotCorInterpFrac * (AFI_interpHi%FullySeparate - AFI_interpLo%FullySeparate)
      AFI_interp%FullyAttached = AFI_interpLo%FullyAttached + rotCorInterpFrac * (AFI_interpHi%FullyAttached - AFI_interpLo%FullyAttached)
   end if

end subroutine AFI_ComputeAirfoilCoefsRotCor1D
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine handles 2D airfoil coefficient calculation with rotation correction
subroutine AFI_ComputeAirfoilCoefsRotCor2D( AOA, SecondProp, p, AFI_interp, errStat, errMsg )

   real(ReKi),                             intent(in   ) :: AOA
   real(ReKi),                             intent(in   ) :: SecondProp              ! Re or UserProp
   TYPE (AFI_ParameterType),               intent(in   ) :: p
   type(AFI_OutputType),                   intent(  out) :: AFI_interp
   integer(IntKi),                         intent(  out) :: errStat
   character(*),                           intent(  out) :: errMsg
   
   ! Local variables
   integer(IntKi)                          :: rotCorTabIdxLo, rotCorTabIdxHi
   real(ReKi)                              :: rotCorInterpFrac
   type(AFI_OutputType)                    :: AFI_interpLo, AFI_interpHi
   integer(IntKi)                          :: ErrStat2
   character(300)                          :: ErrMsg2
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Handle constant data case - similar to AFI_ComputeAirfoilCoefs1D
   if (p%Table(1)%ConstData) then
      ! For constant data, use RotCor interpolation with the sole precomputed table.
      call AFI_Compute2DInterpolationRotCor( AOA, SecondProp, p, 1_IntKi, AFI_interp, errStat, errMsg )
      return
   end if
   
   ! For 2D case, we need to do 2D interpolation on Re/UserProp first, then rot cor interpolation
   ! This requires more complex logic to handle multiple table interpolation
   
   ! Find the table indices for interpolation based on current_snel_factor
   call FindRotCorTableIndices(p%RotCorParams%current_snel_factor, p%Table(1), &
                             rotCorTabIdxLo, rotCorTabIdxHi, rotCorInterpFrac, errStat2, errMsg2)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeAirfoilCoefsRotCor2D')
   if (errStat >= AbortErrLev) return
   
   ! Get coefficients using 2D interpolation on the lower snel factor tables
   call AFI_Compute2DInterpolationRotCor( AOA, SecondProp, p, rotCorTabIdxLo, &
                                          AFI_interpLo, errStat2, errMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeAirfoilCoefsRotCor2D')
   if (errStat >= AbortErrLev) return
   
   if (rotCorTabIdxLo == rotCorTabIdxHi) then
      ! No snel interpolation needed
      AFI_interp = AFI_interpLo
   else
      ! Get coefficients using 2D interpolation on the higher snel factor tables
      call AFI_Compute2DInterpolationRotCor( AOA, SecondProp, p, rotCorTabIdxHi, &
                                             AFI_interpHi, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeAirfoilCoefsRotCor2D')
      if (errStat >= AbortErrLev) return
      
      ! Interpolate between the two snel factor results
      AFI_interp%Cl = AFI_interpLo%Cl + rotCorInterpFrac * (AFI_interpHi%Cl - AFI_interpLo%Cl)
      AFI_interp%Cd = AFI_interpLo%Cd + rotCorInterpFrac * (AFI_interpHi%Cd - AFI_interpLo%Cd)
      AFI_interp%Cm = AFI_interpLo%Cm + rotCorInterpFrac * (AFI_interpHi%Cm - AFI_interpLo%Cm)
      AFI_interp%Cd0 = AFI_interpLo%Cd0 + rotCorInterpFrac * (AFI_interpHi%Cd0 - AFI_interpLo%Cd0)
      AFI_interp%Cm0 = AFI_interpLo%Cm0 + rotCorInterpFrac * (AFI_interpHi%Cm0 - AFI_interpLo%Cm0)
      AFI_interp%Cpmin = AFI_interpLo%Cpmin + rotCorInterpFrac * (AFI_interpHi%Cpmin - AFI_interpLo%Cpmin)
      AFI_interp%f_st = AFI_interpLo%f_st + rotCorInterpFrac * (AFI_interpHi%f_st - AFI_interpLo%f_st)
      AFI_interp%FullySeparate = AFI_interpLo%FullySeparate + rotCorInterpFrac * (AFI_interpHi%FullySeparate - AFI_interpLo%FullySeparate)
      AFI_interp%FullyAttached = AFI_interpLo%FullyAttached + rotCorInterpFrac * (AFI_interpHi%FullyAttached - AFI_interpLo%FullyAttached)
   end if

end subroutine AFI_ComputeAirfoilCoefsRotCor2D
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine handles 1D UA coefficient calculation with rotation correction
subroutine AFI_ComputeUACoefsRotCor1D( p, UA_BL, errStat, errMsg )

   type(AFI_ParameterType), intent(in   ) :: p
   type(AFI_UA_BL_Type),    intent(  out) :: UA_BL
   integer(IntKi),          intent(  out) :: errStat
   character(*),            intent(  out) :: errMsg
   
   ! Local variables
   integer(IntKi)                          :: rotCorTabIdxLo, rotCorTabIdxHi
   real(ReKi)                              :: rotCorInterpFrac
   type(AFI_UA_BL_Type)                    :: UA_BL_Lo, UA_BL_Hi
   integer(IntKi)                          :: ErrStat2
   character(300)                          :: ErrMsg2
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Handle constant data case - similar to AFI_ComputeAirfoilCoefs1D
   if (p%Table(1)%ConstData) then
      ! For constant data, just copy the UA_BL from the first table
      call AFI_CopyUA_BL_Type( p%Table(1)%UA_BL, UA_BL, MESH_NEWCOPY, errStat, errMsg )
      return
   end if
   
   ! Find the table indices for interpolation
   call FindRotCorTableIndices(p%RotCorParams%current_snel_factor, p%Table(1), &
                             rotCorTabIdxLo, rotCorTabIdxHi, rotCorInterpFrac, errStat2, errMsg2)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUACoefsRotCor1D')
   if (errStat >= AbortErrLev) return
   
   ! Copy UA_BL from lower table
   call AFI_CopyUA_BL_Type( p%Table(1)%rotCorTables(rotCorTabIdxLo)%UA_BL, UA_BL_Lo, MESH_NEWCOPY, errStat2, errMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUACoefsRotCor1D')
   if (errStat >= AbortErrLev) return
   
   if (rotCorTabIdxLo == rotCorTabIdxHi) then
      ! No interpolation needed
      UA_BL = UA_BL_Lo
   else
      ! Copy UA_BL from higher table
      call AFI_CopyUA_BL_Type( p%Table(1)%rotCorTables(rotCorTabIdxHi)%UA_BL, UA_BL_Hi, MESH_NEWCOPY, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUACoefsRotCor1D')
      if (errStat >= AbortErrLev) return
      
      ! Interpolate UA parameters (copy structure from lower, then interpolate values)
      UA_BL = UA_BL_Lo
      call InterpolateUABLType(UA_BL_Lo, UA_BL_Hi, rotCorInterpFrac, UA_BL)
      
      ! Clean up
      call AFI_DestroyUA_BL_Type(UA_BL_Hi, errStat2, errMsg2)
   end if
   
   ! Clean up
   call AFI_DestroyUA_BL_Type(UA_BL_Lo, errStat2, errMsg2)

end subroutine AFI_ComputeUACoefsRotCor1D
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine handles 2D UA coefficient calculation with rotation correction
subroutine AFI_ComputeUACoefsRotCor2D( SecondProp, p, UA_BL, errStat, errMsg )

   real(ReKi),              intent(in   ) :: SecondProp                     ! Re or UserProp
   type(AFI_ParameterType), intent(in   ) :: p
   type(AFI_UA_BL_Type),    intent(  out) :: UA_BL
   integer(IntKi),          intent(  out) :: errStat
   character(*),            intent(  out) :: errMsg
   
   ! Local variables
   integer(IntKi)                          :: rotCorTabIdxLo, rotCorTabIdxHi
   real(ReKi)                              :: rotCorInterpFrac
   type(AFI_UA_BL_Type)                    :: UA_BL_Lo, UA_BL_Hi
   integer(IntKi)                          :: ErrStat2
   character(300)                          :: ErrMsg2
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Handle constant data case - similar to AFI_ComputeAirfoilCoefs1D
   if (p%Table(1)%ConstData) then
      ! For constant data, just copy the UA_BL from the first table
      call AFI_CopyUA_BL_Type( p%Table(1)%UA_BL, UA_BL, MESH_NEWCOPY, errStat, errMsg )
      return
   end if
   
   ! Find the table indices for interpolation
   call FindRotCorTableIndices(p%RotCorParams%current_snel_factor, p%Table(1), &
                             rotCorTabIdxLo, rotCorTabIdxHi, rotCorInterpFrac, errStat2, errMsg2)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUACoefsRotCor2D')
   if (errStat >= AbortErrLev) return
   
   ! Get UA coefficients using 2D interpolation on the lower factor tables
   call AFI_ComputeUA2DInterpolationRotCor( SecondProp, p, rotCorTabIdxLo, &
                                            UA_BL_Lo, errStat2, errMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUACoefsRotCor2D')
   if (errStat >= AbortErrLev) return
   
   if (rotCorTabIdxLo == rotCorTabIdxHi) then
      ! No interpolation needed
      UA_BL = UA_BL_Lo
   else
      ! Get UA coefficients using 2D interpolation on the higher factor tables
      call AFI_ComputeUA2DInterpolationRotCor( SecondProp, p, rotCorTabIdxHi, &
                                               UA_BL_Hi, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUACoefsRotCor2D')
      if (errStat >= AbortErrLev) return
      
      ! Interpolate UA parameters between the two snel factor results
      UA_BL = UA_BL_Lo
      call InterpolateUABLType(UA_BL_Lo, UA_BL_Hi, rotCorInterpFrac, UA_BL)
      
      ! Clean up
      call AFI_DestroyUA_BL_Type(UA_BL_Hi, errStat2, errMsg2)
   end if
   
   ! Clean up
   call AFI_DestroyUA_BL_Type(UA_BL_Lo, errStat2, errMsg2)

end subroutine AFI_ComputeUACoefsRotCor2D
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine finds the table indices for interpolation
subroutine FindRotCorTableIndices(current_snel_factor, Table, rotCorTabIdxLo, rotCorTabIdxHi, rotCorInterpFrac, errStat, errMsg)

   real(ReKi),              intent(in   ) :: current_snel_factor
   type(AFI_Table_Type),    intent(in   ) :: Table
   integer(IntKi),          intent(  out) :: rotCorTabIdxLo, rotCorTabIdxHi
   real(ReKi),              intent(  out) :: rotCorInterpFrac
   integer(IntKi),          intent(  out) :: errStat
   character(*),            intent(  out) :: errMsg
   
   ! Local variables
   real(ReKi)                              :: max_snel_factor = 1.0_ReKi
   integer(IntKi)                          :: num_rotCor_tables
   real(ReKi)                              :: snel_factor_step
   
   errStat = ErrID_None
   errMsg  = ""

   if (.not. allocated(Table%rotCorTables)) then
      errStat = ErrID_Fatal
      errMsg  = 'RotCor tables are not allocated.'
      return
   end if

   if (size(Table%rotCorTables) == 1) then
      rotCorTabIdxLo = 1
      rotCorTabIdxHi = 1
      rotCorInterpFrac = 0.0_ReKi
      return
   end if

   num_rotCor_tables = size(Table%rotCorTables) - 1
   if (num_rotCor_tables < 1) then
      errStat = ErrID_Fatal
      errMsg  = 'RotCor table count is invalid.'
      return
   end if
   
   ! Calculate step size (matching AFI_PrecalculateRotCorTables)
   snel_factor_step = max_snel_factor / real(num_rotCor_tables, ReKi)
   
   ! Clamp current_snel_factor to valid range
   if (current_snel_factor <= 0.0_ReKi) then
      rotCorTabIdxLo = 1
      rotCorTabIdxHi = 1
      rotCorInterpFrac = 0.0_ReKi
   elseif (current_snel_factor >= max_snel_factor) then
      rotCorTabIdxLo = num_rotCor_tables + 1
      rotCorTabIdxHi = num_rotCor_tables + 1
      rotCorInterpFrac = 0.0_ReKi
   else
      ! Find the indices for interpolation
      rotCorTabIdxLo = int(current_snel_factor / snel_factor_step) + 1
      rotCorTabIdxHi = rotCorTabIdxLo + 1
      
      ! Ensure indices are within bounds
      rotCorTabIdxLo = max(1, min(rotCorTabIdxLo, num_rotCor_tables + 1))
      rotCorTabIdxHi = max(1, min(rotCorTabIdxHi, num_rotCor_tables + 1))
      
      ! Calculate interpolation fraction
      if (rotCorTabIdxHi > rotCorTabIdxLo) then
         rotCorInterpFrac = (current_snel_factor - (rotCorTabIdxLo - 1) * snel_factor_step) / snel_factor_step
      else
         rotCorInterpFrac = 0.0_ReKi
      end if
   end if

end subroutine FindRotCorTableIndices
!----------------------------------------------------------------------------------------------------------------------------------
subroutine FindBoundingTables(p, secondProp, lowerTable, upperTable, xVals)

   TYPE (AFI_ParameterType),               intent(in   ) :: p
   real(ReKi),                             intent(in   ) :: secondProp
   integer(IntKi),                         intent(  out) :: lowerTable, upperTable
   real(ReKi),                             intent(  out) :: xVals(2)

   integer(IntKi)                                        :: i

   lowerTable = 1
   upperTable = p%NumTabs

   do i = 1, p%NumTabs-1
      if ( secondProp >= p%secondVals(i) .and. secondProp <= p%secondVals(i+1) ) then
         lowerTable = i
         upperTable = i+1
         exit
      end if
   end do

   xVals(1) = p%secondVals(lowerTable)
   xVals(2) = p%secondVals(upperTable)

end subroutine FindBoundingTables
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine handles 2D interpolation across Re/UserProp for a specific table index
subroutine AFI_Compute2DInterpolationRotCor( AOA, SecondProp, p, rotCorTabIdx, AFI_interp, errStat, errMsg )

   real(ReKi),                             intent(in   ) :: AOA
   real(ReKi),                             intent(in   ) :: SecondProp
   TYPE (AFI_ParameterType),               intent(in   ) :: p
   integer(IntKi),                         intent(in   ) :: rotCorTabIdx
   type(AFI_OutputType),                   intent(  out) :: AFI_interp
   integer(IntKi),                         intent(  out) :: errStat
   character(*),                           intent(  out) :: errMsg
   
   ! Local variables
   integer                                 :: lowerTable, upperTable
   real(ReKi)                              :: xVals(2)
   type(AFI_OutputType)                    :: AFI_lower, AFI_upper
   integer(IntKi)                          :: ErrStat2
   character(300)                          :: ErrMsg2
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Handle boundary conditions first
   IF ( SecondProp <= p%secondVals( 1 ) )  THEN
      ! Use the first table's table
      call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(1)%rotCorTables(rotCorTabIdx), &
                                                 p, AFI_interp, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_Compute2DInterpolationRotCor')
      return
   ELSE IF ( SecondProp >= p%secondVals( p%NumTabs ) ) THEN
      ! Use the last table's table
      call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(p%NumTabs)%rotCorTables(rotCorTabIdx), &
                                                 p, AFI_interp, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_Compute2DInterpolationRotCor')
      return
   END IF
   
   ! Find bounding tables for interpolation
   call FindBoundingTables(p, SecondProp, lowerTable, upperTable, xVals)
   
   ! Get coefficients from lower table's table
   call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(lowerTable)%rotCorTables(rotCorTabIdx), &
                                              p, AFI_lower, errStat2, errMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_Compute2DInterpolationRotCor')
   if (errStat >= AbortErrLev) return
   
   ! Get coefficients from upper table's table
   call AFI_ComputeAirfoilCoefsFromRotCorTable( AOA, p%Table(upperTable)%rotCorTables(rotCorTabIdx), &
                                              p, AFI_upper, errStat2, errMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_Compute2DInterpolationRotCor')
   if (errStat >= AbortErrLev) return

   ! Linearly interpolate between the two results
   call AFI_Output_ExtrapInterp1(AFI_lower, AFI_upper, xVals, AFI_interp, SecondProp, ErrStat2, ErrMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_Compute2DInterpolationRotCor')

end subroutine AFI_Compute2DInterpolationRotCor
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine handles 2D UA interpolation for a specific table index
subroutine AFI_ComputeUA2DInterpolationRotCor( SecondProp, p, rotCorTabIdx, UA_BL, errStat, errMsg )

   real(ReKi),                             intent(in   ) :: SecondProp
   TYPE (AFI_ParameterType),               intent(in   ) :: p
   integer(IntKi),                         intent(in   ) :: rotCorTabIdx
   type(AFI_UA_BL_Type),                   intent(  out) :: UA_BL
   integer(IntKi),                         intent(  out) :: errStat
   character(*),                           intent(  out) :: errMsg
   
   ! Local variables
   real(ReKi)                              :: xVals(2)
   integer                                 :: lowerTable, upperTable
   integer(IntKi)                          :: ErrStat2
   character(300)                          :: ErrMsg2
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Handle boundary conditions first
   IF ( SecondProp <= p%secondVals( 1 ) )  THEN
      call AFI_CopyUA_BL_Type( p%Table(1)%rotCorTables(rotCorTabIdx)%UA_BL, UA_BL, MESH_NEWCOPY, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUA2DInterpolationRotCor')
      return
   ELSE IF ( SecondProp >= p%secondVals( p%NumTabs ) ) THEN
      call AFI_CopyUA_BL_Type( p%Table(p%NumTabs)%rotCorTables(rotCorTabIdx)%UA_BL, UA_BL, MESH_NEWCOPY, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUA2DInterpolationRotCor')
      return
   END IF

   call FindBoundingTables(p, SecondProp, lowerTable, upperTable, xVals)

   ! Linearly interpolate UA parameters between tables
   call AFI_UA_BL_Type_ExtrapInterp1(p%Table(lowerTable)%rotCorTables(rotCorTabIdx)%UA_BL, &
                                     p%Table(upperTable)%rotCorTables(rotCorTabIdx)%UA_BL, &
                                     xVals, UA_BL, SecondProp, ErrStat2, ErrMsg2 )
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AFI_ComputeUA2DInterpolationRotCor')

end subroutine AFI_ComputeUA2DInterpolationRotCor
!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine interpolates between two UA_BL_Type structures
subroutine InterpolateUABLType(UA_BL_Lo, UA_BL_Hi, InterpFrac, UA_BL_Out)

   type(AFI_UA_BL_Type),    intent(in   ) :: UA_BL_Lo, UA_BL_Hi
   real(ReKi),              intent(in   ) :: InterpFrac
   type(AFI_UA_BL_Type),    intent(inout) :: UA_BL_Out
   
   ! Linear interpolation of UA parameters
   ! Note: UA_BL_Out should already be initialized (copied from UA_BL_Lo)
   
   UA_BL_Out%alpha0            = UA_BL_Lo%alpha0     + InterpFrac * (UA_BL_Hi%alpha0     - UA_BL_Lo%alpha0)				
   UA_BL_Out%alpha1            = UA_BL_Lo%alpha1     + InterpFrac * (UA_BL_Hi%alpha1     - UA_BL_Lo%alpha1)                
   UA_BL_Out%alpha2            = UA_BL_Lo%alpha2     + InterpFrac * (UA_BL_Hi%alpha2     - UA_BL_Lo%alpha2)                
   UA_BL_Out%eta_e             = UA_BL_Lo%eta_e      + InterpFrac * (UA_BL_Hi%eta_e      - UA_BL_Lo%eta_e)       
   UA_BL_Out%C_nalpha          = UA_BL_Lo%C_nalpha   + InterpFrac * (UA_BL_Hi%C_nalpha   - UA_BL_Lo%C_nalpha)       
   UA_BL_Out%C_lalpha          = UA_BL_Lo%C_lalpha   + InterpFrac * (UA_BL_Hi%C_lalpha   - UA_BL_Lo%C_lalpha)                             
   UA_BL_Out%T_f0              = UA_BL_Lo%T_f0       + InterpFrac * (UA_BL_Hi%T_f0       - UA_BL_Lo%T_f0)    																											  
   UA_BL_Out%T_V0              = UA_BL_Lo%T_V0       + InterpFrac * (UA_BL_Hi%T_V0       - UA_BL_Lo%T_V0)                                                                                                                
   UA_BL_Out%T_p               = UA_BL_Lo%T_p        + InterpFrac * (UA_BL_Hi%T_p        - UA_BL_Lo%T_p)                                                                                                               
   UA_BL_Out%T_VL              = UA_BL_Lo%T_VL       + InterpFrac * (UA_BL_Hi%T_VL       - UA_BL_Lo%T_VL)                                                                                                              
   UA_BL_Out%b1                = UA_BL_Lo%b1         + InterpFrac * (UA_BL_Hi%b1         - UA_BL_Lo%b1)                                                                                                               
   UA_BL_Out%b2                = UA_BL_Lo%b2         + InterpFrac * (UA_BL_Hi%b2         - UA_BL_Lo%b2)                                                                                                               
   UA_BL_Out%b5                = UA_BL_Lo%b5         + InterpFrac * (UA_BL_Hi%b5         - UA_BL_Lo%b5)                                                                                                               
   UA_BL_Out%A1                = UA_BL_Lo%A1         + InterpFrac * (UA_BL_Hi%A1         - UA_BL_Lo%A1)                                                                                                               
   UA_BL_Out%A2                = UA_BL_Lo%A2         + InterpFrac * (UA_BL_Hi%A2         - UA_BL_Lo%A2)                                                                                                          
   UA_BL_Out%A5                = UA_BL_Lo%A5         + InterpFrac * (UA_BL_Hi%A5         - UA_BL_Lo%A5)
   UA_BL_Out%S1                = UA_BL_Lo%S1         + InterpFrac * (UA_BL_Hi%S1         - UA_BL_Lo%S1)
   UA_BL_Out%S2                = UA_BL_Lo%S2         + InterpFrac * (UA_BL_Hi%S2         - UA_BL_Lo%S2)
   UA_BL_Out%S3                = UA_BL_Lo%S3         + InterpFrac * (UA_BL_Hi%S3         - UA_BL_Lo%S3)
   UA_BL_Out%S4                = UA_BL_Lo%S4         + InterpFrac * (UA_BL_Hi%S4         - UA_BL_Lo%S4)
   UA_BL_Out%Cn1               = UA_BL_Lo%Cn1        + InterpFrac * (UA_BL_Hi%Cn1        - UA_BL_Lo%Cn1)
   UA_BL_Out%Cn2               = UA_BL_Lo%Cn2        + InterpFrac * (UA_BL_Hi%Cn2        - UA_BL_Lo%Cn2)
   UA_BL_Out%St_sh             = UA_BL_Lo%St_sh      + InterpFrac * (UA_BL_Hi%St_sh      - UA_BL_Lo%St_sh)
   UA_BL_Out%Cd0               = UA_BL_Lo%Cd0        + InterpFrac * (UA_BL_Hi%Cd0        - UA_BL_Lo%Cd0)
   UA_BL_Out%Cm0               = UA_BL_Lo%Cm0        + InterpFrac * (UA_BL_Hi%Cm0        - UA_BL_Lo%Cm0)   
   UA_BL_Out%k0                = UA_BL_Lo%k0         + InterpFrac * (UA_BL_Hi%k0         - UA_BL_Lo%k0)
   UA_BL_Out%k1                = UA_BL_Lo%k1         + InterpFrac * (UA_BL_Hi%k1         - UA_BL_Lo%k1)
   UA_BL_Out%k2                = UA_BL_Lo%k2         + InterpFrac * (UA_BL_Hi%k2         - UA_BL_Lo%k2)
   UA_BL_Out%k3                = UA_BL_Lo%k3         + InterpFrac * (UA_BL_Hi%k3         - UA_BL_Lo%k3)
   UA_BL_Out%k1_hat            = UA_BL_Lo%k1_hat     + InterpFrac * (UA_BL_Hi%k1_hat     - UA_BL_Lo%k1_hat)
   UA_BL_Out%x_cp_bar          = UA_BL_Lo%x_cp_bar   + InterpFrac * (UA_BL_Hi%x_cp_bar   - UA_BL_Lo%x_cp_bar)
   UA_BL_Out%UACutout          = UA_BL_Lo%UACutout   + InterpFrac * (UA_BL_Hi%UACutout   - UA_BL_Lo%UACutout)
   UA_BL_Out%UACutout_delta    = UA_BL_Lo%UACutout_delta + InterpFrac * (UA_BL_Hi%UACutout_delta - UA_BL_Lo%UACutout_delta) 
   UA_BL_Out%UACutout_blend    = UA_BL_Lo%UACutout_blend + InterpFrac * (UA_BL_Hi%UACutout_blend - UA_BL_Lo%UACutout_blend)  
   UA_BL_Out%filtCutOff        = UA_BL_Lo%filtCutOff + InterpFrac * (UA_BL_Hi%filtCutOff - UA_BL_Lo%filtCutOff)
   UA_BL_Out%alphaUpper        = UA_BL_Lo%alphaUpper + InterpFrac * (UA_BL_Hi%alphaUpper - UA_BL_Lo%alphaUpper) 
   UA_BL_Out%alphaLower        = UA_BL_Lo%alphaLower + InterpFrac * (UA_BL_Hi%alphaLower - UA_BL_Lo%alphaLower)  
   UA_BL_Out%c_alphaLower      = UA_BL_Lo%c_alphaLower + InterpFrac * (UA_BL_Hi%c_alphaLower - UA_BL_Lo%c_alphaLower)  
   UA_BL_Out%c_alphaUpper      = UA_BL_Lo%c_alphaUpper + InterpFrac * (UA_BL_Hi%c_alphaUpper - UA_BL_Lo%c_alphaUpper)  
   UA_BL_Out%alpha0ReverseFlow = UA_BL_Lo%alpha0ReverseFlow + InterpFrac * (UA_BL_Hi%alpha0ReverseFlow - UA_BL_Lo%alpha0ReverseFlow)  
   UA_BL_Out%alphaBreakUpper   = UA_BL_Lo%alphaBreakUpper + InterpFrac * (UA_BL_Hi%alphaBreakUpper - UA_BL_Lo%alphaBreakUpper)  
   UA_BL_Out%CnBreakUpper      = UA_BL_Lo%CnBreakUpper + InterpFrac * (UA_BL_Hi%CnBreakUpper - UA_BL_Lo%CnBreakUpper)  
   UA_BL_Out%alphaBreakLower   = UA_BL_Lo%alphaBreakLower + InterpFrac * (UA_BL_Hi%alphaBreakLower - UA_BL_Lo%alphaBreakLower)  
   UA_BL_Out%CnBreakLower      = UA_BL_Lo%CnBreakLower + InterpFrac * (UA_BL_Hi%CnBreakLower - UA_BL_Lo%CnBreakLower)  

end subroutine InterpolateUABLType
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AFI_CalcSnel(AFInfo, snel_factor)
   implicit none
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
   real(ReKi),                  intent(  out)   :: snel_factor
   
   ! Local variables
   real(ReKi) :: r_over_R, chord_over_R, tsr, tsr_local
   real(ReKi) :: chord, rLocal, rMax
   
   ! Get parameters from AFInfo structure
   tsr = AFInfo%RotCorParams%tsr
   rLocal = AFInfo%RotCorParams%rLocal
   rMax = AFInfo%RotCorParams%rMax
   chord = AFInfo%RotCorParams%chord
   
   ! Calculate normalized ratios with safe divide
   r_over_R = merge(rLocal/rMax, 0.0_ReKi, rMax > tiny(1.0_ReKi))
   chord_over_r = merge(chord/rLocal, 0.0_ReKi, rLocal > tiny(1.0_ReKi))
   tsr_local = tsr * r_over_R
   
   ! Calculate Snel factor 
   ! H. Snel, R. Houwink, and W. J. Piers. 1993. Linenberg, 2003.
   snel_factor = 0.0_ReKi
   if (r_over_R < 0.80_ReKi ) then
      snel_factor = (3.1_ReKi * tsr_local**2 / (1.0_ReKi + tsr_local**2)) * (chord_over_r**2) 
   end if    
   
end subroutine AFI_CalcSnel
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AFI_ApplySnel(AOA_vec, Cl_vec, snel_factor, Cl_0, slope)
    implicit none
    real(ReKi),      intent(in)     :: AOA_vec(:)
    real(ReKi),      intent(inout)  :: Cl_vec(:)
    real(ReKi),      intent(in)     :: snel_factor
    real(ReKi),      intent(in)     :: Cl_0, slope
    
    ! Local variables
    integer(IntKi)  :: num_aoa
    real(ReKi), allocatable :: alpha_deg(:), g(:)
    real(ReKi), allocatable :: cl_lin(:), delta_cl(:), correction(:)
    logical, allocatable :: mask1(:), mask2(:)
    
    ! Get the size of the input arrays for allocation
    num_aoa = size(AOA_vec)
    
    ! Allocate all local arrays based on the number of AOA points
    allocate(alpha_deg(num_aoa), g(num_aoa))
    allocate(cl_lin(num_aoa), delta_cl(num_aoa), correction(num_aoa))
    allocate(mask1(num_aoa), mask2(num_aoa))

    ! Vectorized calculations
    alpha_deg = AOA_vec * 180.0_ReKi / Pi_D
    cl_lin = slope * AOA_vec + Cl_0
    delta_cl = cl_lin - Cl_vec
    
    ! Vectorized blending factor calculation
	! From QBlade documentation. https://docs.qblade.org/src/theory/aerodynamics/secondary_effects/himmelskamp.html#himmelskamp-effect
    mask1 = (alpha_deg > 0.0_ReKi) .AND. (alpha_deg < 30.0_ReKi)
    mask2 = (alpha_deg >= 30.0_ReKi) .AND. (alpha_deg < 60.0_ReKi)
    
    g = 0.0_ReKi  ! default for alpha_deg >= 60
    where (mask1)
        g = 1.0_ReKi
    end where
    where (mask2)
        g = 0.5_ReKi * (1.0_ReKi + cos(D2R*(6.0_ReKi*alpha_deg - 180.0_ReKi)))
    end where
    
    ! Vectorized final correction
    correction = snel_factor * g * delta_cl
    Cl_vec = Cl_vec + correction
    
    deallocate(alpha_deg, g, cl_lin, delta_cl, correction, mask1, mask2)

end subroutine AFI_ApplySnel
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AFI_ApplySchepers(AOA_vec, Cl_vec, Cd_vec, schepers_factor, cl_slope, cl_0)
    implicit none
    ! Inputs
    real(ReKi), intent(in)    :: AOA_vec(:)          ! AoA [rad]
    real(ReKi), intent(in)    :: schepers_factor     ! f_cl at 20° (amplitude incl. TSR/phi/c/r)
    real(ReKi), intent(in)    :: cl_slope            ! dCl/dalpha [per rad]
    real(ReKi), intent(in)    :: cl_0                ! Cl at alpha=0
    ! In/out
    real(ReKi), intent(inout) :: Cl_vec(:), Cd_vec(:)

    integer(IntKi)               :: n
    real(ReKi)                   :: D2R, eps
    real(ReKi)                   :: a20_rad, a25_deg
    real(ReKi), allocatable      :: alpha_deg(:), abs_a(:)

    ! ===== LIFT =====
    real(ReKi), allocatable      :: cl_lin(:), dcl(:), fcl(:), w(:)
    logical,   allocatable       :: m0_20(:), m20_25(:), m25_60(:)

    real(ReKi)                   :: cl20, cllin20, denom_cl20

    ! ===== DRAG =====
    real(ReKi), allocatable      :: fcd(:), t(:)
    logical,   allocatable       :: m0_15(:), m15_20(:), m20_30(:), m30_45(:)
    real(ReKi)                   :: cdmin_used, cd15, denom_cd15, fcd15, fcd20

    n   = size(AOA_vec)
    D2R = 180.0_ReKi / Pi_D
    eps = 1.0e-12_ReKi

    a20_rad = 20.0_ReKi / D2R
    a25_deg = 25.0_ReKi

    allocate(alpha_deg(n), abs_a(n))
    alpha_deg = AOA_vec * D2R
    abs_a     = abs(alpha_deg)

    ! =========================
    ! 1) LIFT correction (paper §1)
    ! =========================
    allocate(cl_lin(n), dcl(n), fcl(n), w(n))
    allocate(m0_20(n), m20_25(n), m25_60(n))

    ! potential/linear lift and deficit
    cl_lin = cl_slope*AOA_vec + cl_0
    dcl    = cl_lin - Cl_vec

    ! normalization at |alpha|=20°
    cl20      = interp_on_abs(alpha_deg, Cl_vec, 20.0_ReKi)
    cllin20   = cl_slope*a20_rad + cl_0
    denom_cl20= MAX(eps, ABS(cllin20 - cl20))

    ! masks
    m0_20  = (abs_a <= 20.0_ReKi)
    m20_25 = (abs_a >  20.0_ReKi) .AND. (abs_a <= 25.0_ReKi)
    m25_60 = (abs_a >  25.0_ReKi) .AND. (abs_a <= 60.0_ReKi)

    fcl = 0.0_ReKi
    w   = 0.0_ReKi

    ! 0..20° : proportional to dCl relative to value at 20°
    where (m0_20)
        fcl = schepers_factor * MAX(0.0_ReKi, MIN(1.0_ReKi, ABS(dcl)/denom_cl20))
    end where

    ! 20..25° : constant at schepers_factor
    where (m20_25)
        fcl = schepers_factor
    end where

    ! 25..60° : linear drop to zero at 60°
    where (m25_60)
        w   = (60.0_ReKi - abs_a) / (60.0_ReKi - a25_deg)   ! 1 at 25°, 0 at 60°
        w   = MAX(0.0_ReKi, MIN(1.0_ReKi, w))
        fcl = schepers_factor * w
    end where

    ! Apply: Cl_3D = Cl_2D + fcl*(Cl_lin - Cl_2D)
    Cl_vec = Cl_vec + fcl * dcl

    ! =========================
    ! 2) DRAG correction (paper §2)
    ! =========================
    allocate(fcd(n), t(n))
    allocate(m0_15(n), m15_20(n), m20_30(n), m30_45(n))

    ! Cd_min from low-alpha if available (|α|<=15°), else global min
    if (ANY(abs_a <= 15.0_ReKi)) then
        cdmin_used = MINVAL(Cd_vec, MASK=(abs_a <= 15.0_ReKi))
    else
        cdmin_used = MINVAL(Cd_vec)
    end if

    cd15        = interp_on_abs(alpha_deg, Cd_vec, 15.0_ReKi)
    denom_cd15  = MAX(eps, cd15 - cdmin_used)

    ! anchors
    fcd15 = 6.0_ReKi * schepers_factor
    fcd20 = 2.5_ReKi * schepers_factor

    ! masks
    m0_15  = (abs_a <= 15.0_ReKi)
    m15_20 = (abs_a >  15.0_ReKi) .AND. (abs_a <= 20.0_ReKi)
    m20_30 = (abs_a >  20.0_ReKi) .AND. (abs_a <= 30.0_ReKi)
    m30_45 = (abs_a >  30.0_ReKi) .AND. (abs_a <= 45.0_ReKi)

    fcd = 0.0_ReKi
    t   = 0.0_ReKi

    ! 0..15° : proportional to local excess drag
    where (m0_15)
        fcd = fcd15 * MAX(0.0_ReKi, MIN(1.0_ReKi, (Cd_vec - cdmin_used)/denom_cd15))
    end where

    ! 15..20° : linear fcd15 -> fcd20
    where (m15_20)
        t   = (abs_a - 15.0_ReKi)/5.0_ReKi
        fcd = (1.0_ReKi - t)*fcd15 + t*fcd20
    end where

    ! 20..30° : plateau
    where (m20_30)
        fcd = fcd20
    end where

    ! 30..45° : linear drop to zero
    where (m30_45)
        t   = (45.0_ReKi - abs_a)/15.0_ReKi
        t   = MAX(0.0_ReKi, MIN(1.0_ReKi, t))
        fcd = fcd20 * t
    end where

    ! Apply: Cd_3D = Cd_2D + fcd*(Cd_2D - Cd_min)
    Cd_vec = Cd_vec + fcd * (Cd_vec - cdmin_used)

    ! cleanup
    deallocate(alpha_deg, abs_a, cl_lin, dcl, fcl, w, fcd, t, m0_20, m20_25, m25_60, m0_15, m15_20, m20_30, m30_45)

contains
    ! Interpolate y at target |alpha| (deg) using absolute(alpha_deg)
    pure function interp_on_abs(alpha_deg_all, y_all, target_abs_deg) result(yq)
        real(ReKi), intent(in) :: alpha_deg_all(:), y_all(:), target_abs_deg
        real(ReKi)             :: yq
        integer(IntKi)         :: i, n, ilo, ihi
        real(ReKi)             :: xa, xlo, xhi, ylo, yhi, tt
        n   = size(alpha_deg_all)
        ilo = -1; ihi = -1
        xlo = -HUGE(1.0_ReKi); xhi = HUGE(1.0_ReKi)
        do i = 1, n
            xa = abs(alpha_deg_all(i))
            if (xa <= target_abs_deg .AND. xa > xlo) then
                xlo = xa; ilo = i
            end if
            if (xa >= target_abs_deg .AND. xa < xhi) then
                xhi = xa; ihi = i
            end if
        end do
        if (ilo < 0 .AND. ihi < 0) then
            yq = y_all(1)
        else if (ilo < 0) then
            yq = y_all(ihi)
        else if (ihi < 0) then
            yq = y_all(ilo)
        else if (ABS(xhi - xlo) <= 1.0e-12_ReKi) then
            yq = 0.5_ReKi*(y_all(ilo) + y_all(ihi))
        else
            ylo = y_all(ilo); yhi = y_all(ihi)
            tt  = (target_abs_deg - xlo)/(xhi - xlo)
            yq  = ylo + tt*(yhi - ylo)
        end if
    end function interp_on_abs
end subroutine AFI_ApplySchepers
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AFI_ApplyGertzCd(AOA_vec, Cl_vec, Cd_vec, snel_factor)
   implicit none
   real(ReKi),      intent(in)     :: AOA_vec(:), Cl_vec(:)
   real(ReKi),      intent(inout)  :: Cd_vec(:)
   real(ReKi),      intent(in)     :: snel_factor

   integer(IntKi)               :: i, n
   real(ReKi), allocatable      :: alpha_deg(:), alpha_rad(:)
   real(ReKi), allocatable      :: Cd_orig(:), Cd_fp(:)
   real(ReKi)                   :: D2R, tiny
   real(ReKi)                   :: xk(6), yk(6), y2(6), Cd_fp_cont_pt_1,Cd_fp_cont_pt_2

   n    = size(AOA_vec)
   D2R  = Pi_D/180.0_ReKi
   tiny = 1.0e-12_ReKi

   allocate(alpha_deg(n), alpha_rad(n))
   allocate(Cd_orig(n), Cd_fp(n))

   alpha_rad = AOA_vec
   alpha_deg = AOA_vec * 180.0_ReKi / Pi_D

   ! --- preserve original Cd for control points ---
   Cd_orig = Cd_vec

   ! ================== CD: spline(15..60) with requested control points ==================
   if (snel_factor > 0.0_ReKi) then ! no correction unless Snel factor is active
   
      ! Flat-plate target using CURRENT Cl (after Snel): L/D_fp = cot(a) => Cd_fp = |Cl|*|tan(a)|
      do i = 1, n
         Cd_fp(i) = abs(Cl_vec(i)) * abs(tan(alpha_rad(i)))
      end do
   
      Cd_fp_cont_pt_1 = 25.0_ReKi
      Cd_fp_cont_pt_2 = 30.0_ReKi
      ! Control-point abscissae [deg]
      xk = (/ 10.0_ReKi, 15.0_ReKi, Cd_fp_cont_pt_1, Cd_fp_cont_pt_2, 60.0_ReKi, 70.0_ReKi /)
   
      ! Control-point ordinates:
      ! y(10)=Cd_orig(10), y(15)=Cd_orig(15), y(20)=Cd_fp(20), y(30)=Cd_fp(30), y(60)=Cd_orig(60), y(70)=Cd_orig(70)
      yk(1) = linear_interp(alpha_deg, Cd_orig, 10.0_ReKi)
      yk(2) = linear_interp(alpha_deg, Cd_orig, 15.0_ReKi)
      yk(3) = linear_interp(alpha_deg, Cd_fp,   Cd_fp_cont_pt_1)*0.9
      yk(4) = linear_interp(alpha_deg, Cd_fp,   Cd_fp_cont_pt_2)*0.9
      yk(5) = linear_interp(alpha_deg, Cd_orig, 60.0_ReKi)
      yk(6) = linear_interp(alpha_deg, Cd_orig, 70.0_ReKi)
   
      ! Build natural cubic spline on those 6 control points
      call cubic_spline_build(xk, yk, y2)
   
      ! Apply ONLY between 15° and 60°; keep original elsewhere
      do i = 1, n
         if (alpha_deg(i) > 15.0_ReKi .and. alpha_deg(i) < 60.0_ReKi) then
            Cd_vec(i) = cubic_spline_eval(xk, yk, y2, alpha_deg(i))
         else
            Cd_vec(i) = Cd_orig(i)
         end if
      end do

   end if ! snel_factor > 0
   
   deallocate(alpha_deg, alpha_rad, Cd_orig, Cd_fp)

contains
   !---------------- linear interpolation on a 1D grid (clamps to ends) ----------------
   pure function linear_interp(x, y, xq) result(yq)
      real(ReKi), intent(in) :: x(:), y(:), xq
      real(ReKi)             :: yq
      integer(IntKi)         :: j, nloc
      real(ReKi)             :: x1, x2, y1, y2, t

      nloc = size(x)
      if (nloc == 1) then
         yq = y(1); return
      end if

      ! clamp
      if (xq <= x(1)) then
         yq = y(1); return
      elseif (xq >= x(nloc)) then
         yq = y(nloc); return
      end if

      ! find bracket (assumes mostly increasing; still safe if not perfectly uniform)
      do j = 1, nloc-1
         if ( (xq >= x(j) .and. xq <= x(j+1)) .or. (xq >= x(j+1) .and. xq <= x(j)) ) then
            x1 = x(j);   x2 = x(j+1)
            y1 = y(j);   y2 = y(j+1)
            if (abs(x2 - x1) < 1.0e-12_ReKi) then
               yq = 0.5_ReKi*(y1 + y2)
            else
               t  = (xq - x1) / (x2 - x1)
               yq = y1 + t*(y2 - y1)
            end if
            return
         end if
      end do

      ! fallback (shouldn’t hit)
      yq = y(nloc)
   end function linear_interp

   !-------------------- Natural cubic spline (Numerical-Recipes style) -------------------
   subroutine cubic_spline_build(x, y, y2)
      real(ReKi), intent(in)  :: x(:), y(:)
      real(ReKi), intent(out) :: y2(:)
      integer(IntKi)          :: n, i, k
      real(ReKi), allocatable :: u(:)
      real(ReKi)              :: sig, p

      n = size(x)
      if (size(y) /= n .or. size(y2) /= n) stop 'spline_build: size mismatch'
      allocate(u(n-1))
      y2(1) = 0.0_ReKi
      u(1)  = 0.0_ReKi

      do i = 2, n-1
         sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
         p   = sig*y2(i-1) + 2.0_ReKi
         y2(i) = (sig - 1.0_ReKi) / p
         u(i)  = ((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1)))
         u(i)  = (6.0_ReKi*u(i)/(x(i+1)-x(i-1)) - sig*u(i-1)) / p
      end do

      y2(n) = 0.0_ReKi
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1) + u(k)
      end do
      deallocate(u)
   end subroutine cubic_spline_build

   function cubic_spline_eval(x, y, y2, xx) result(yy)
      real(ReKi), intent(in) :: x(:), y(:), y2(:), xx
      real(ReKi)             :: yy
      integer(IntKi)         :: klo, khi, k, n
      real(ReKi)             :: h, a, b

      n = size(x)
      ! clamp
      if (xx <= x(1)) then
         yy = y(1); return
      elseif (xx >= x(n)) then
         yy = y(n); return
      end if

      klo = 1; khi = n
      do while (khi - klo > 1)
         k = (khi + klo)/2
         if (x(k) > xx) then
            khi = k
         else
            klo = k
         end if
      end do

      h = x(khi) - x(klo)
      if (h <= 0.0_ReKi) then
         yy = 0.5_ReKi*(y(klo)+y(khi)); return
      end if
      a = (x(khi) - xx)/h
      b = (xx - x(klo))/h
      yy = a*y(klo) + b*y(khi) + ((a**3 - a)*y2(klo) + (b**3 - b)*y2(khi))*(h*h)/6.0_ReKi
   end function cubic_spline_eval

end subroutine AFI_ApplyGertzCd
!=============================================================================
   

END MODULE AirfoilInfo_RotCor

