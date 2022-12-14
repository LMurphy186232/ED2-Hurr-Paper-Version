!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 09/30/08 Including the cumulus/radiation feedback.                                 !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
!   This is the main driver for the radiation effect on the thermodynamic properties.      !
!------------------------------------------------------------------------------------------!
subroutine radiate(mzp,mxp,myp,ia,iz,ja,jz,mynum)
   use mem_tend   ,  only: tend_g         ! ! intent(inout)
   use mem_grid   ,  only: ngrid          & ! intent(in)
                          ,time           & ! intent(in)
                          ,dtlt           & ! intent(in)
                          ,itimea         & ! intent(in)
                          ,nzg            & ! intent(in)
                          ,nzs            & ! intent(in)
                          ,npatch         & ! intent(in)
                          ,grid_g         & ! intent(in)
                          ,nnzp           & ! intent(in)
                          ,if_adap        & ! intent(in)
                          ,zm             & ! intent(in)
                          ,zt             ! ! intent(in)
   use mem_leaf   ,  only: leaf_g         ! ! intent(inout)
   use mem_radiate,  only: ilwrtyp        & ! intent(in)
                          ,iswrtyp        & ! intent(in)
                          ,icumfdbk       & ! intent(in)
                          ,radiate_g      & ! intent(inout)
                          ,radfrq         & ! intent(in)
                          ,ncall_i        & ! intent(inout)
                          ,prsnz          & ! intent(in)
                          ,prsnzp         & ! intent(in)
                          ,jday           & ! intent(inout)
                          ,solfac         & ! intent(inout)
                          ,nadd_rad       & ! intent(in)
                          ,ncrad          ! ! intent(in)
   use mem_basic  ,  only: basic_g        ! ! intent(in)
   use mem_scratch,  only: vctr1          & ! intent(inout)
                          ,vctr2          & ! intent(inout)
                          ,vctr3          & ! intent(inout)
                          ,vctr4          & ! intent(inout)
                          ,vctr5          & ! intent(inout)
                          ,vctr6          & ! intent(inout)
                          ,vctr7          & ! intent(inout)
                          ,vctr8          & ! intent(inout)
                          ,vctr9          & ! intent(inout)
                          ,vctr10         & ! intent(inout)
                          ,vctr11         & ! intent(inout)
                          ,vctr12         & ! intent(inout)
                          ,scratch        ! ! intent(inout)
   use mem_micro   , only: micro_g        ! ! intent(in)
   use therm_lib   , only: uint2tl        & ! subroutine
                          ,level          & ! intent(in)
                          ,cloud_on       & ! intent(in)
                          ,bulk_on        ! ! intent(in)
   use micphys     , only: availcat       ! ! intent(in)
   use mem_cuparm  , only: cuparm_g       & ! intent(in)
                          ,nnqparm        & ! intent(in)
                          ,nclouds        ! ! intent(in)
   use rad_carma   , only: radcomp_carma  ! ! subroutine
   use catt_start  , only: catt           ! ! intent(in)
   use mem_scalar  , only: scalar_g       ! ! intent(in)
   use rconstants  , only: hr_sec         ! ! intent(in)

   implicit none
  
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)  :: mzp
   integer                     , intent(in)  :: mxp
   integer                     , intent(in)  :: myp
   integer                     , intent(in)  :: ia
   integer                     , intent(in)  :: iz
   integer                     , intent(in)  :: ja
   integer                     , intent(in)  :: jz
   integer                     , intent(in)  :: mynum
   !----- Local variables -----------------------------------------------------------------!
   integer                                   :: koff
   integer                                   :: nrad
   integer                                   :: i
   integer                                   :: j
   integer                                   :: k
   integer                                   :: ka
   integer                                   :: kz
   integer                                   :: icld
   real                                      :: solc
   real, dimension(:,:)        , allocatable :: st_rain
   real, dimension(:,:,:)      , allocatable :: cb_rain
   real, dimension(:,:,:)      , allocatable :: lwl
   real, dimension(:,:,:)      , allocatable :: iwl
   real                                      :: tcoal
   real                                      :: fracliq
   real                                      :: max_albedt
   real                                      :: max_rlongup
   real                                      :: time_rfrq
   !---------------------------------------------------------------------------------------!

   !----- Leave this driver if the user turned off the radiation --------------------------!
   if (ilwrtyp + iswrtyp == 0) return



   kz = mzp-1

   !---------------------------------------------------------------------------------------!
   !     Transfering some variables to scratch arrays, or filling the scratch arrays with  !
   ! zeroes in case the variables are unavailable.                                         !
   !---------------------------------------------------------------------------------------!
   call rad_copy2scratch(mzp,mxp,myp,npatch)

   !----- Updating the temperature tendency -----------------------------------------------!
   call tend_accum(mzp,mxp,myp,ia,iz,ja,jz,tend_g(ngrid)%tht,radiate_g(ngrid)%fthrd)

   !----- If this is the first time Harrington is called, run initialization. -------------!  
   if ((iswrtyp == 3 .or. ilwrtyp == 3) .and. ncall_i == 0) then
      !------------------------------------------------------------------------------------!
      !    If first call for this node, initialize several quantities & mclatchy sounding  !
      ! data.                                                                              !
      !------------------------------------------------------------------------------------!
      call harr_radinit(nnzp(1))
      ncall_i = ncall_i + 1
   end if

   !----- Checking whether this is time to update the radiative forcing. ------------------!
   time_rfrq = real(dmod(time + 0.001,dble(radfrq)))

   if ( time_rfrq  < dtlt .or. time < 0.001) then

      !----- Reset the radiative forcing --------------------------------------------------!
      call azero(mzp*mxp*myp,radiate_g(ngrid)%fthrd     )
      call azero(mzp*mxp*myp,radiate_g(ngrid)%fthrd_lw  )
      call azero(mxp*myp,radiate_g(ngrid)%par_beam      )
      call azero(mxp*myp,radiate_g(ngrid)%par_diffuse   )
      call azero(mxp*myp,radiate_g(ngrid)%nir_beam      )
      call azero(mxp*myp,radiate_g(ngrid)%nir_diffuse   )
      call azero(mxp*myp,radiate_g(ngrid)%rshort        )
      call azero(mxp*myp,radiate_g(ngrid)%rshort_diffuse)
      call azero(mxp*myp,radiate_g(ngrid)%rshort_top    )
      call azero(mxp*myp,radiate_g(ngrid)%rshortup_top  )
      call azero(mxp*myp,radiate_g(ngrid)%rlongup_top   )
      call azero(mxp*myp,radiate_g(ngrid)%rlong         )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute solar zenith angle.                                                     !
      !------------------------------------------------------------------------------------!
      call zen(mxp,myp,ia,iz,ja,jz,iswrtyp,ilwrtyp,grid_g(ngrid)%glon,grid_g(ngrid)%glat   &
              ,radiate_g(ngrid)%cosz)
      !------------------------------------------------------------------------------------!



      !----- CARMA radiation --------------------------------------------------------------!
      if (ilwrtyp == 4 .or. iswrtyp == 4) then
         allocate(st_rain(    mxp,myp)        )
         allocate(cb_rain(    mxp,myp,nclouds))
         allocate(lwl    (mzp,mxp,myp        ))
         allocate(iwl    (mzp,mxp,myp        ))

         call azero(mzp*mxp*myp        ,lwl    )
         call azero(mzp*mxp*myp        ,iwl    )
         call azero(    mxp*myp        ,st_rain)
         call azero(    mxp*myp*nclouds,cb_rain)


         !------ Copy precipitation -------------------------------------------------------!
         if (nnqparm(ngrid) > 0 ) call atob(mxp*myp*nclouds,cuparm_g(ngrid)%conprr,cb_rain)
         if (bulk_on)             call atob(mxp*myp        ,micro_g(ngrid)%pcpg   ,st_rain)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Add hydrometeors to lwl and iwl. This is safe even when some or all pure     !
         ! liquid or pure ice hydrometeors are unavailable, because in this case the       !
         ! arrays have been set to zero.                                                   !
         !---------------------------------------------------------------------------------!
         call ae1p1(mzp*mxp*myp,lwl,lwl,scratch%vt3da) ! Cloud
         call ae1p1(mzp*mxp*myp,lwl,lwl,scratch%vt3db) ! Rain
         call ae1p1(mzp*mxp*myp,iwl,iwl,scratch%vt3dc) ! Pristine ice
         call ae1p1(mzp*mxp*myp,iwl,iwl,scratch%vt3dd) ! Snow
         call ae1p1(mzp*mxp*myp,iwl,iwl,scratch%vt3de) ! Aggregates
         !---------------------------------------------------------------------------------!


         !----- Graupel: if available, we need to split between ice and liquid. -----------!
         if (availcat(6)) then
            do j=ja,jz
               do i=ia,iz
                  ka = nint(grid_g(ngrid)%flpw(i,j))
                  do k=ka,kz
                     call uint2tl(micro_g(ngrid)%q6(k,i,j),tcoal,fracliq)
                     lwl(k,i,j) = lwl(k,i,j) + fracliq*micro_g(ngrid)%rgp(k,i,j)
                     iwl(k,i,j) = iwl(k,i,j) + (1.-fracliq)*micro_g(ngrid)%rgp(k,i,j)
                  end do
               end do
            end do
         end if
         !---------------------------------------------------------------------------------!


         !----- Hail: if available, we need to split between ice and liquid. --------------!
         if (availcat(7)) then
            do j=ja,jz
               do i=ia,iz
                  ka = nint(grid_g(ngrid)%flpw(i,j))
                  do k=ka,kz
                     call uint2tl(micro_g(ngrid)%q7(k,i,j),tcoal,fracliq)
                     lwl(k,i,j) = lwl(k,i,j) + fracliq*micro_g(ngrid)%rhp(k,i,j)
                     iwl(k,i,j) = iwl(k,i,j) + (1.-fracliq)*micro_g(ngrid)%rhp(k,i,j)
                  end do
               end do
            end do
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!




         !----- Call CARMA radiation driver -----------------------------------------------!
         call radcomp_carma(mzp,mxp,myp,npatch,nclouds,ncrad,ia,iz,ja,jz,mynum,iswrtyp     &
                           ,ilwrtyp,icumfdbk,solfac                                        &
                           ,basic_g(ngrid)%theta           ,basic_g(ngrid)%pi0             &
                           ,basic_g(ngrid)%pp              ,basic_g(ngrid)%rv              &
                           ,st_rain,cb_rain,lwl,iwl        ,scratch%vt4da                  &
                           ,scratch%vt4dc                  ,scratch%vt3dr                  &
                           ,scratch%vt3ds                  ,basic_g(ngrid)%dn0             &
                           ,basic_g(ngrid)%rtp             ,radiate_g(ngrid)%fthrd         &
                           ,grid_g(ngrid)%rtgt             ,grid_g(ngrid)%f13t             &
                           ,grid_g(ngrid)%f23t             ,grid_g(ngrid)%glat             &
                           ,grid_g(ngrid)%glon             ,radiate_g(ngrid)%par_beam      &
                           ,radiate_g(ngrid)%par_diffuse   ,radiate_g(ngrid)%nir_beam      &
                           ,radiate_g(ngrid)%nir_diffuse   ,radiate_g(ngrid)%rshort        &
                           ,radiate_g(ngrid)%rshort_diffuse,radiate_g(ngrid)%rshort_top    &
                           ,radiate_g(ngrid)%rshortup_top  ,radiate_g(ngrid)%rlong         &
                           ,radiate_g(ngrid)%rlongup_top   ,radiate_g(ngrid)%albedt        &
                           ,radiate_g(ngrid)%cosz          ,radiate_g(ngrid)%rlongup       &
                           ,grid_g(ngrid)%fmapt            ,scalar_g(3,ngrid)%sclp         &
                           ,leaf_g(ngrid)%patch_area       )
         !---------------------------------------------------------------------------------!



         !----- Free memory. --------------------------------------------------------------!
         deallocate(st_rain)
         deallocate(cb_rain)
         deallocate(lwl    )
         deallocate(iwl    )
         !---------------------------------------------------------------------------------!
      end if

      if (ilwrtyp <= 2 .or. iswrtyp <= 2) then
         !----- If using Mahrer-Pielke and/or Chen-Cotton radiation, call radcomp. --------!
         call radcomp(mzp,mxp,myp,ngrid,ia,iz,ja,jz                                        &
                           ,basic_g(ngrid)%theta          ,basic_g(ngrid)%pi0              &
                           ,basic_g(ngrid)%pp             ,basic_g(ngrid)%rv               &
                           ,basic_g(ngrid)%dn0            ,basic_g(ngrid)%rtp              &
                           ,scratch%vt3do                 ,radiate_g(ngrid)%fthrd          &
                           ,grid_g(ngrid)%rtgt            ,grid_g(ngrid)%f13t              &
                           ,grid_g(ngrid)%f23t            ,grid_g(ngrid)%glon              &
                           ,radiate_g(ngrid)%rshort       ,radiate_g(ngrid)%rlong          &
                           ,radiate_g(ngrid)%albedt       ,radiate_g(ngrid)%cosz           &
                           ,radiate_g(ngrid)%rlongup      ,radiate_g(ngrid)%fthrd_lw       &
                           ,mynum                         )
      end if

      !----- Using Harrington radiation ---------------------------------------------------!
      if (iswrtyp == 3 .or. ilwrtyp == 3) then
         call harr_raddriv(mzp,mxp,myp,nclouds,ncrad,ngrid,if_adap,time,dtlt,ia,iz,ja,jz   &
                          ,nadd_rad,iswrtyp,ilwrtyp,icumfdbk                               &
                          ,grid_g(ngrid)%flpw             ,grid_g(ngrid)%topt              &
                          ,grid_g(ngrid)%glon             ,grid_g(ngrid)%glat              &
                          ,grid_g(ngrid)%rtgt             ,basic_g(ngrid)%pi0              &
                          ,basic_g(ngrid)%pp              ,basic_g(ngrid)%dn0              &
                          ,basic_g(ngrid)%theta           ,basic_g(ngrid)%rv               &
                          ,scratch%vt3do                  ,radiate_g(ngrid)%rshort         &
                          ,radiate_g(ngrid)%rshort_diffuse,radiate_g(ngrid)%rlong          &
                          ,radiate_g(ngrid)%fthrd         ,radiate_g(ngrid)%rlongup        &
                          ,radiate_g(ngrid)%cosz          ,radiate_g(ngrid)%albedt         &
                          ,radiate_g(ngrid)%rshort_top    ,radiate_g(ngrid)%rshortup_top   &
                          ,radiate_g(ngrid)%rlongup_top   ,radiate_g(ngrid)%fthrd_lw       &
                          ,scratch%vt3da                  ,scratch%vt3db                   &
                          ,scratch%vt3dc                  ,scratch%vt3dd                   &
                          ,scratch%vt3de                  ,scratch%vt3df                   &
                          ,scratch%vt3dg                  ,scratch%vt3dh                   &
                          ,scratch%vt3di                  ,scratch%vt3dj                   &
                          ,scratch%vt3dk                  ,scratch%vt3dl                   &
                          ,scratch%vt3dm                  ,scratch%vt3dn                   &
                          ,scratch%vt4da                  ,scratch%vt4dc                   &
                          ,scratch%vt3dr                  ,scratch%vt3ds                   &
                          ,mynum                          )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If shortwave radiation isn't CARMA, split the components using the Weiss and   !
      ! Norman (1985) method.                                                              !
      !------------------------------------------------------------------------------------!
      select case (iswrtyp)
      case (1,2,3)
         do j=ja,jz
            do i=ia,iz
               call wn85_rshort_split( radiate_g(ngrid)%rshort         (i,j)               &
                                     , radiate_g(ngrid)%cosz           (i,j)               &
                                     , radiate_g(ngrid)%par_beam       (i,j)               &
                                     , radiate_g(ngrid)%par_diffuse    (i,j)               &
                                     , radiate_g(ngrid)%nir_beam       (i,j)               &
                                     , radiate_g(ngrid)%nir_diffuse    (i,j)               &
                                     , radiate_g(ngrid)%rshort_diffuse (i,j)               )
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      case default
         continue
      end select
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!
   return
end subroutine radiate
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies some variables that may or may not be assigned depending on    !
! the user choices. This is essentially the microphysics, cumulus scheme, and TEB vari-    !
! ables. In case they are not assigned, the corresponding scratch arrays are filled with   !
! zeroes.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine rad_copy2scratch(mzp,mxp,myp,mpp)
   use mem_basic     , only : basic_g   & ! intent(in)
                            , co2_on    & ! intent(in)
                            , co2con    ! ! intent(in)
   use mem_grid      , only : ngrid     ! ! intent(in)
   use micphys       , only : availcat  & ! intent(in)
                             ,progncat  ! ! intent(in)
   use mem_micro     , only : micro_g   ! ! intent(in)
   use mem_cuparm    , only : nnqparm   & ! intent(in)
                             ,nclouds   & ! intent(in)
                             ,cuparm_g  ! ! intent(in)
   use mem_radiate   , only : icumfdbk  ! ! intent(in)
   use mem_scratch   , only : scratch   ! ! intent(inout)
   use mem_leaf      , only : leaf_g    ! ! intent(in)
   use mem_teb_common, only : tebc_g    ! ! intent(in)
   use teb_spm_start , only : teb_spm   ! ! intent(in)
   use rconstants    , only : mmcod1em6 ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: mzp
   integer, intent(in) :: mxp
   integer, intent(in) :: myp
   integer, intent(in) :: mpp
   !----- Local variables. ----------------------------------------------------------------!
   real                :: rco2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Resetting scratch arrays for Harrington/Bulk microphysics/Cumulus interface. The  !
   ! actual values will be copied depending on the microphysics level and whether the      !
   ! cumulus is on.                                                                        !
   !---------------------------------------------------------------------------------------!
   !----- Cloud droplets ------------------------------------------------------------------!
   if (availcat(1)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rcp,scratch%vt3da)
   else
      call azero(mzp*mxp*myp,scratch%vt3da)
   end if
   !----- Rain drops ----------------------------------------------------------------------!
   if (availcat(2)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rrp,scratch%vt3db)
   else
      call azero(mzp*mxp*myp,scratch%vt3db)
   end if
   !----- Pristine ice --------------------------------------------------------------------!
   if (availcat(3)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rpp,scratch%vt3dc)
   else
      call azero(mzp*mxp*myp,scratch%vt3dc)
   end if
   !----- Snow  ---------------------------------------------------------------------------!
   if (availcat(4)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rsp,scratch%vt3dd)
   else
      call azero(mzp*mxp*myp,scratch%vt3dd)
   end if
   !----- Aggregates ----------------------------------------------------------------------!
   if (availcat(5)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rap,scratch%vt3de)
   else
      call azero(mzp*mxp*myp,scratch%vt3de)
   end if
   !----- Aggregates ----------------------------------------------------------------------!
   if (availcat(6)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rgp,scratch%vt3df)
   else
      call azero(mzp*mxp*myp,scratch%vt3df)
   end if
   !----- Aggregates ----------------------------------------------------------------------!
   if (availcat(7)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%rhp,scratch%vt3dg)
   else
      call azero(mzp*mxp*myp,scratch%vt3dg)
   end if
   !----- Cloud droplets ------------------------------------------------------------------!
   if (progncat(1)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%ccp,scratch%vt3dh)
   else
      call azero(mzp*mxp*myp,scratch%vt3dh)
   end if
   !----- Rain drops ----------------------------------------------------------------------!
   if (progncat(2)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%crp,scratch%vt3di)
   else
      call azero(mzp*mxp*myp,scratch%vt3di)
   end if
   !----- Pristine ice --------------------------------------------------------------------!
   if (progncat(3)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%cpp,scratch%vt3dj)
   else
      call azero(mzp*mxp*myp,scratch%vt3dj)
   end if
   !----- Snow  ---------------------------------------------------------------------------!
   if (progncat(4)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%csp,scratch%vt3dk)
   else
      call azero(mzp*mxp*myp,scratch%vt3dk)
   end if
   !----- Aggregates ----------------------------------------------------------------------!
   if (progncat(5)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%cap,scratch%vt3dl)
   else
      call azero(mzp*mxp*myp,scratch%vt3dl)
   end if
   !----- Aggregates ----------------------------------------------------------------------!
   if (progncat(6)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%cgp,scratch%vt3dm)
   else
      call azero(mzp*mxp*myp,scratch%vt3dm)
   end if
   !----- Aggregates ----------------------------------------------------------------------!
   if (progncat(7)) then
      call atob(mzp*mxp*myp,micro_g(ngrid)%chp,scratch%vt3dn)
   else
      call azero(mzp*mxp*myp,scratch%vt3dn)
   end if

   !---------------------------------------------------------------------------------------!
   !     Copy CO2 array to scratch in case CO2 is prognosed, otherwise, dump CO2CON.  The  !
   ! radiation schemes expect CO2 mixing ratio in kg_CO2/kg_air, so we must convert the    !
   ! units here (BRAMS default CO2 unit is ppm, or ?mol_CO2/mol_air).                      !
   !---------------------------------------------------------------------------------------!
   if (co2_on) then
      call ae1t0(mzp*mxp*myp,scratch%vt3do,basic_g(ngrid)%co2p,mmcod1em6)
   else
      rco2 = co2con(1) * mmcod1em6
      call ae0(mzp*mxp*myp,scratch%vt3do,rco2)
   end if

   !---------------------------------------------------------------------------------------!
   !     Checking whether the user wants the cumulus clouds to interact with the radiation !
   ! scheme. This will happen only if the user is running Harrington scheme, though.       !
   !---------------------------------------------------------------------------------------!
   call azero(mxp*myp*nclouds,scratch%vt3dr)
   if (icumfdbk == 1 .and. nnqparm(ngrid) > 0) then
      call atob (mzp*mxp*myp*nclouds,cuparm_g(ngrid)%cuprliq,scratch%vt4da)
      call atob (mzp*mxp*myp*nclouds,cuparm_g(ngrid)%cuprice,scratch%vt4dc)
      call accum(    mxp*myp*nclouds,scratch%vt3dr,cuparm_g(ngrid)%areadn)
      call accum(    mxp*myp*nclouds,scratch%vt3dr,cuparm_g(ngrid)%areaup)
      call atob (    mxp*myp*nclouds,cuparm_g(ngrid)%xierr,scratch%vt3ds)
   else
      call azero(mzp*mxp*myp*nclouds,scratch%vt4da)
      call azero(mzp*mxp*myp*nclouds,scratch%vt4dc)
      call aone (    mxp*myp*nclouds,scratch%vt3ds)
   end if

   !---------------------------------------------------------------------------------------!
   !     Using scratch arrays to deal with TEB. The old-style interface was causing compi- !
   ! lation problems with ifort 10 and full interface so it was removed.                   !
   !---------------------------------------------------------------------------------------!
   if (teb_spm == 1) then
      call atob(mxp*myp    ,tebc_g(ngrid)%emis_town,scratch%vt2da)
      call atob(mxp*myp    ,tebc_g(ngrid)%alb_town ,scratch%vt2db)
      call atob(mxp*myp    ,tebc_g(ngrid)%ts_town  ,scratch%vt2dc)
      call atob(mxp*myp*mpp,leaf_g(ngrid)%g_urban  ,scratch%vt3dp)
   else
      call azero(mxp*myp    ,scratch%vt2da)
      call azero(mxp*myp    ,scratch%vt2db)
      call azero(mxp*myp    ,scratch%vt2dc)
      call azero(mxp*myp*mpp,scratch%vt3dp)
   end if

   return
end subroutine rad_copy2scratch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine accumulates the tendency due to one process to the overall tendency.  !
!------------------------------------------------------------------------------------------!
subroutine tend_accum(m1,m2,m3,ia,iz,ja,jz,totaltend,onetend)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1,m2,m3,ia,iz,ja,jz
   real, dimension(m1,m2,m3), intent(in)    :: onetend
   real, dimension(m1,m2,m3), intent(inout) :: totaltend
   !----- Local variables -----------------------------------------------------------------!
   integer                                  :: i,j,k
   !---------------------------------------------------------------------------------------!

   do j = ja,jz
      do i = ia,iz
         do k = 1,m1
            totaltend(k,i,j) = totaltend(k,i,j) + onetend(k,i,j)
         end do
      end do
   end do

   return
end subroutine tend_accum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This is the wrapper for Mahrer-Pielke and Chen-Cotton radiation methods. This sub-   !
! routine is different from the original because now it uses the new variables computed at !
! the "zen" subroutine to find the hour angle. Also, most sin/cos calculations are done in !
! double precision.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine radcomp(m1,m2,m3,ifm,ia,iz,ja,jz,theta,pi0,pp,rv,dn0,rtp,co2p,fthrd,rtgt,f13t   &
                  ,f23t,glon,rshort,rlong,albedt,cosz,rlongup,fthrd_lw,mynum)
   use mem_grid   , only : dzm             & ! intent(in)
                         , dzt             & ! intent(in)
                         , itopo           & ! intent(in)
                         , plonn           & ! intent(in)
                         , ngrid           & ! intent(in)
                         , time            & ! intent(in)
                         , itimea          & ! intent(in)
                         , centlon         & ! intent(in)
                         , dtlt            ! ! intent(in)
   use mem_scratch, only : scratch         ! ! intent(in)
   use mem_radiate, only : ilwrtyp         & ! intent(in)
                         , iswrtyp         & ! intent(in)
                         , lonrad          & ! intent(in)
                         , cdec            & ! intent(in)
                         , jday            & ! intent(in)
                         , solfac          & ! intent(in)
                         , sun_longitude   & ! intent(in)
                         , rad_cosz_min    ! ! intent(in)
   use rconstants , only : cpdryi          & ! intent(in)
                         , cpor            & ! intent(in)
                         , p00             & ! intent(in)
                         , stefan          & ! intent(in)
                         , solar           & ! intent(in)
                         , pio1808         & ! intent(in)
                         , pi1             & ! intent(in)
                         , halfpi          ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                          , intent(in)   :: m1,m2,m3,ifm,ia,iz,ja,jz,mynum
   real        , dimension(m1,m2,m3), intent(in)   :: theta,pi0,pp,rv,dn0,rtp,co2p
   real        , dimension(   m2,m3), intent(in)   :: rtgt,f13t,f23t,glon
   real        , dimension(   m2,m3), intent(in)   :: cosz,albedt,rlongup
   real        , dimension(m1,m2,m3), intent(inout):: fthrd,fthrd_lw
   real        , dimension(   m2,m3), intent(inout):: rshort,rlong
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k,kk
   real(kind=8)                                :: dzsdx,dzsdy,dlon,a1,a2,hrangl,sinz
   real(kind=8)                                :: sazmut,slazim,slangl,cosi,gglon 
   real        , dimension(m1)                 :: rvr,rtr,co2r,dn0r,pird,prd
   real        , dimension(m1)                 :: dzmr,dztr,temprd,fthrl,fthrs
   !----- Constants -----------------------------------------------------------------------!
   real(kind=8), parameter :: offset=1.d-20
   !---------------------------------------------------------------------------------------!


   !----- Constant longitude if lonrad is set to 0 ----------------------------------------!
   gglon=dble(centlon(1))

   do j = ja,jz
      do i = ia,iz
         !------ Initialising the fluxes, this way they will be zero if not called --------!
         do k = 1,m1
            fthrl(k) = 0.
            fthrs(k) = 0.
         end do
         do k = 1,m1
            !---- Compute some basic thermodynamic variables (pressure, temperature). -----!
            pird(k) = (pp(k,i,j) + pi0(k,i,j)) * cpdryi
            temprd(k) = theta(k,i,j) * pird(k)
            rvr(k)  = max(0.,rv(k,i,j))
            rtr(k)  = max(rvr(k),rtp(k,i,j))
            co2r(k) = co2p(k,i,j)
            !----- Convert the next 4 variables to cgs for now. ---------------------------!
            prd(k)  = pird(k) ** cpor * p00 * 10.
            dn0r(k) = dn0(k,i,j) * 1.e-3
            dzmr(k) = dzm(k) / rtgt(i,j) * 1.e-2
            dztr(k) = dzt(k) / rtgt(i,j) * 1.e-2

            fthrl(k) = 0.
            fthrs(k) = 0.

         end do
         temprd(1) = sqrt(sqrt(rlongup(i,j) / stefan))
         !---------------------------------------------------------------------------------!

         !----- Sanity check --------------------------------------------------------------!
         do k=1,m1
            if (prd(k)  <   0. .or. dn0r(k)   <   0. .or.  rtr(k)    <   0. .or.           &
                co2r(k) <   0. .or. temprd(k) < 160.                             ) then
               !---------------------------------------------------------------------------!
               ! TL(k) < 160.: This is -113 C, which is much colder than the Vostok,       !
               !               Antarctica world record and should also be colder than any  !
               !               atmospheric temperature, so to avoid frostbites we will     !
               !               halt the run. The other conditions are non-physical values  !
               !               so the run is already in trouble.                           !
               !---------------------------------------------------------------------------!
               write (unit=*,fmt='(a)') '================================================='
               write (unit=*,fmt='(a)') ' ERROR - rad_comp!!!'
               write (unit=*,fmt='(a)') '         The model is about to stop!'
               write (unit=*,fmt='(2(a,1x,i5,a))') ' - Node:',mynum,' Grid: ',ifm
               write (unit=*,fmt='(3(a,1x,i5,a))') ' - k = ',k,' i = ',i,' j = ',j
               write (unit=*,fmt='(a)') ' - Either the temperature is too low, or some'
               write (unit=*,fmt='(a)') '   negative density, mixing ratio or pressure!'
               write (unit=*,fmt='(a)') ' - Sanity check at Chen-Cotton/Mahrer-Pielke:'
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               write (unit=*,fmt='(a3,1x,5(a12,1x))')                                      &
                                        'LEV','  MIX. RATIO','        CO_2','     DENSITY' &
                                             ,'    PRESSURE',' TEMPERATURE'
               do kk=1,m1
                  write (unit=*,fmt='(i3,1x,5(es12.3,1x))')                                &
                                       kk, rtr(kk), co2r(kk), dn0r(kk), prd(kk), temprd(kk)
               enddo
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               write (unit=*,fmt='(a)') ' '
               write (unit=*,fmt='(a)') '================================================='
               call abort_run ('Weird thermodynamic values, caught at radiation'           &
                              ,'radcomp','rad_driv.f90')
            end if
         end do

         !----- Call the longwave parameterizations. --------------------------------------!
         select case (ilwrtyp)
         case (1) !----- Chen-Cotton (1983) -----------------------------------------------!
            call lwradc(m1,rvr,rtr,co2r,dn0r,temprd,prd,dztr,fthrl,rlong(i,j))
         case (2) !----- Mahrer-Pielke (1977) ---------------------------------------------!
            call lwradp(m1,temprd,rvr,co2r,dn0r,dztr,pird,scratch%vt3dq,fthrl,rlong(i,j))
         end select

         !---------------------------------------------------------------------------------!
         !     The shortwave parameterizations are only valid if the cosine of the zenith  !
         ! angle is greater than cosz_min.  Otherwise this is nighttime or polar night,    !
         ! and the shortwave is close to zero anyway.                                      !
         !---------------------------------------------------------------------------------!
         if (cosz(i,j) > rad_cosz_min) then
            select case (iswrtyp)
            case (1) !----- Chen-Cotton (1983) --------------------------------------------!
               call shradc(m1,rvr,rtr,dn0r,dztr,prd,albedt(i,j),solar*1.e3*solfac          &
                          ,cosz(i,j),fthrs,rshort(i,j))
            case (2) !----- Mahrer-Pielke (1977) ------------------------------------------!
               call shradp(m1,rvr,dn0r,dzmr,scratch%vt3dq,pird,cosz(i,j),albedt(i,j)       &
                          ,solar*1e3*solfac,fthrs,rshort(i,j))
            end select

            !------------------------------------------------------------------------------!
            !     Modify the downward surface shortwave flux by considering the slope of   !
            ! the topography.
            !------------------------------------------------------------------------------!
            if (itopo == 1) then                   
               dzsdx = dble(f13t(i,j) * rtgt(i,j)) 
               dzsdy = dble(f23t(i,j) * rtgt(i,j)) 
               !---------------------------------------------------------------------------!
               !     The y- and x-directions must be true north and east for this correc-  !
               ! tion. the following rotates the model y/x to the true north/east.         !
               !---------------------------------------------------------------------------!
               dlon = (dble(plonn(ngrid)) - dble(glon(i,j))) * pio1808
               a1    = dzsdx*dcos(dlon) + dzsdy * dsin(dlon)
               a2    = -dzsdx*dsin(dlon) + dzsdy * dcos(dlon)
               dzsdx = a1
               dzsdy = a2

               if (lonrad == 1) gglon = dble(glon(i,j))
               hrangl = (sun_longitude-gglon)*pio1808
               !---------------------------------------------------------------------------!
               !    sinz tends to zero when cosz(i,j) approaches 1. To avoid this singu-   !
               ! larity, include an offset.                                                !
               !---------------------------------------------------------------------------!
               sinz   = max(offset,sqrt(dble(1.) - dble(cosz(i,j)) ** 2))
               sazmut = dasin(max(dble(-1.),min(dble(1.),cdec*dsin(hrangl)/dble(sinz))))
              
               !----- Offset but preserving the sign... -----------------------------------!
               if (abs(dzsdx) < offset) dzsdx = sign(offset,dzsdx)
               if (abs(dzsdy) < offset) dzsdy = sign(offset,dzsdy)

               slazim      = halfpi - datan2(dzsdy,dzsdx)
               slangl      = datan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
               cosi        = dcos(slangl)*dble(cosz(i,j))                                       &
                           + dsin(slangl)*sinz*dcos(sazmut-slazim)
               rshort(i,j) = rshort(i,j) * sngl(cosi) / cosz(i,j)
            end if

         else
            !----- Make all shortwave to be zero (night time) -----------------------------!
            do k = 1,m1
               fthrs(k) = 0.
            end do
            rshort(i,j) = 0.
         end if

         fthrd(1,i,j)    = fthrd(2,i,j)
         fthrd_lw(1,i,j) = fthrd_lw(2,i,j)

         !----- Convert the downward flux at the ground to SI. ----------------------------!
         rshort(i,j)     = rshort(i,j) * 1.e-3 / (1. - albedt(i,j))
         rlong(i,j)      = rlong(i,j) * 1.e-3
      end do
   end do
   return
end subroutine radcomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zen(m2,m3,ia,iz,ja,jz,iswrtyp,ilwrtyp,glon,glat,cosz)

   use mem_grid   , only : nzpmax          & ! intent(in)
                         , imontha         & ! intent(in)
                         , idatea          & ! intent(in)
                         , iyeara          & ! intent(in)
                         , time            & ! intent(in)
                         , itimea          & ! intent(in)
                         , centlat         & ! intent(in)
                         , centlon         ! ! intent(in)
   use mem_radiate, only : lonrad          & ! intent(in)
                         , jday            & ! intent(out)
                         , solfac          & ! intent(out)
                         , sdec            & ! intent(out)
                         , cdec            & ! intent(out)
                         , declin          & ! intent(out) 
                         , sun_longitude   ! ! intent(out)
   use rconstants , only : pio1808         & ! intent(in)
                         , twopi8          & ! intent(in)
                         , day_sec8        & ! intent(in)
                         , day_hr8         & ! intent(in)
                         , hr_sec8         & ! intent(in)
                         , min_sec8        ! ! intent(in)
   use mem_mclat,   only : mclat_spline    ! ! subroutine
   use mem_harr,    only : nsolb           & ! intent(in)
                         , solar0          & ! intent(in)
                         , solar1          ! ! intent(out)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer               , intent(in)    :: m2,m3           ! Grid dimensions
   integer               , intent(in)    :: ia,iz,ja,jz     ! Node dimensions
   integer               , intent(in)    :: iswrtyp,ilwrtyp ! Radiation scheme flag
   real, dimension(m2,m3), intent(in)    :: glon,glat       ! Grid coordinates
   real, dimension(m2,m3), intent(inout) :: cosz
   !----- Local variables -----------------------------------------------------------------!
   integer                 :: i,j                ! Grid counters
   integer                 :: is                 ! Solar band counter
   real(kind=8)            :: dayhr              ! Hour of day without longitude correction
   real(kind=8)            :: dayhrr             ! Hour of day with longitude correction
   real(kind=8)            :: d0                 ! coefficient for solfac computation
   real(kind=8)            :: d02                ! coefficient for solfac computation
   real(kind=8)            :: radlat             ! latitude in radians
   real(kind=8)            :: hrangl             ! Hour angle in radians
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter :: capri    = -2.35d1 ! Tropic of Capricornium latitude
   real(kind=8), parameter :: ndaysyr  =  3.65d2 ! Number of days of year (no leap years)
   integer     , parameter :: shsummer = -9      ! Julian day of Southern Hemisphere summer
   !----- External functions. -------------------------------------------------------------!
   integer, external   :: julday                 ! Function to compute day of year
   !---------------------------------------------------------------------------------------!


   !----- Find current Julian day ---------------------------------------------------------!
   jday = julday(imontha,idatea,iyeara)
   jday = jday + floor(time/day_sec8)

   !---------------------------------------------------------------------------------------!
   !     Solfac is a multiplier of the solar constant to correct for Earth's varying       !
   ! distance to the sun.                                                                  !
   !---------------------------------------------------------------------------------------!
   d0 =  twopi8 * dble(jday-1) / ndaysyr
   d02 = d0 * 2.d0
   solfac = sngl(1.000110d0 + 3.4221d-2 * dcos (d0) + 1.280d-3 * dsin(d0)                  &
                            + 7.1900d-4 * dcos(d02) + 7.700d-5 * dsin(d02))

   !----- Check whether Harrington shortwave or longwave radiation is used ----------------!
   if (iswrtyp == 3 .or. ilwrtyp == 3) then
      !------------------------------------------------------------------------------------!
      !    Adjust solar fluxes at top of atmosphere for current Earth-Sun distance for     !
      ! Harrington shortwave radiation.                                                    !
      !------------------------------------------------------------------------------------!
      do is = 1,nsolb
         solar1(is) = solar0(is) * solfac
      end do
      !------------------------------------------------------------------------------------!
      !      Interpolate Mclatchy soundings between summer and winter values, and prepare  !
      ! spline coefficients for interpolation by latitude.                                 !
      !------------------------------------------------------------------------------------!
      call mclat_spline(jday)
   end if


   !----- Declin is the solar latitude in degrees (aka declination) -----------------------!
   declin = capri * cos(twopi8 * dble(jday - shsummer) / ndaysyr) * pio1808
   sdec = dsin(declin) !-----  sdec - sine   of declination -------------------------------!
   cdec = dcos(declin) !-----  cdec - cosine of declination -------------------------------!

   !----- Convert the time into hour of the day. ------------------------------------------!
   dayhr = time / hr_sec8 + dble(itimea/100) + dble(mod(itimea,100)) / min_sec8

  
   !----- Compute the cosine of zenith angle ----------------------------------------------!
   if (lonrad == 0) then
      !----- No longitude variation, use centlon to define the hour angle -----------------!
      dayhrr = mod(dayhr + centlon(1)/1.5d1 + day_hr8,day_hr8)
      hrangl = 1.5d1 * (dayhrr - 1.2d1) *pio1808
      do j=ja,jz
         do i=ia,iz
            radlat=dble(glat(i,j))*pio1808
            radlat=dble(glat(i,j))*pio1808
            cosz(i,j) = sngl(dsin(radlat)*sdec+dcos(radlat)*cdec*dcos(hrangl))
            !----- Making sure that it is bounded -----------------------------------------!
            cosz(i,j) = max(-1.0,min(1.0,cosz(i,j)))
         end do
      end do
   else
      !----- Use actual position to find the hour angle -----------------------------------!
      do j=ja,jz
         do i=ia,iz
            dayhrr = mod(dayhr + dble(glon(i,j))/1.5d1 + day_hr8,day_hr8)
            hrangl = 1.5d1 * (dayhrr - 1.2d1) *pio1808
            radlat=dble(glat(i,j))*pio1808
            if (radlat == declin) radlat = radlat + 1.d-5
            cosz(i,j) = sngl(dsin(radlat)*sdec+dcos(radlat)*cdec*dcos(hrangl))
            !----- Making sure that it is bounded -----------------------------------------!
            cosz(i,j) = max(-1.,min(1.,cosz(i,j)))
         end do
      end do
   end if

   return
end subroutine zen
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine computes the split between direct and diffuse radiation, and        !
! between visible and near-infrared radiation using the method suggested by:               !
!                                                                                          !
! Weiss, A., J. M. Norman, 1985: Partitioning solar radiation into direct and diffuse,     !
!     visible and near-infrared components.  Agric. For. Meteorol., 34, 205-213. (WN85)    !
!------------------------------------------------------------------------------------------!
subroutine wn85_rshort_split(rshort_full,cosz,par_beam,par_diff,nir_beam,nir_diff          &
                            ,rshort_diff)
   use rconstants , only : solar         & ! intent(in)
                         , prefsea       & ! intent(in)
                         , twothirds     ! ! intent(in)
   use mem_radiate, only : rad_cosz_min  & ! intent(in)
                         , ref_fvis_beam & ! intent(in)
                         , ref_fvis_diff & ! intent(in)
                         , ref_fnir_beam & ! intent(in)
                         , ref_fnir_diff ! ! intent(in)
   use leaf_coms  , only : atm_prss      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in)    :: rshort_full          ! Incident SW radiation   (total)   [  W/m?]
   real, intent(in)    :: cosz                 ! cos(zenith distance)              [   ---]
   real, intent(out)   :: par_beam             ! Incident PAR            (direct ) [  W/m?]
   real, intent(out)   :: par_diff             ! Incident PAR            (diffuse) [  W/m?]
   real, intent(out)   :: nir_beam             ! Incident near-infrared  (direct ) [  W/m?]
   real, intent(out)   :: nir_diff             ! Incident near-infrared  (diffuse) [  W/m?]
   real, intent(out)   :: rshort_diff          ! Incident SW radiation   (diffuse) [  W/m?]
   !----- Local variables. ----------------------------------------------------------------!
   real                :: par_beam_top         ! PAR at the top of atmosphere      [  W/m?]
   real                :: nir_beam_top         ! NIR at the top of atmosphere      [  W/m?]
   real                :: par_beam_pot         ! Potential incident PAR  (direct)  [  W/m?]
   real                :: par_diff_pot         ! Potential incident PAR  (diffuse) [  W/m?]
   real                :: par_full_pot         ! Potential incident PAR  (total)   [  W/m?]
   real                :: nir_beam_pot         ! Potential  PAR          (direct)  [  W/m?]
   real                :: nir_diff_pot         ! Potential  PAR          (diffuse) [  W/m?]
   real                :: nir_full_pot         ! Potential  PAR          (total)   [  W/m?]
   real                :: par_full             ! Actual incident PAR     (total)   [  W/m?]
   real                :: nir_full             ! Actual incident NIR     (total)   [  W/m?]
   real                :: fvis_beam_act        ! Actual beam fraction - PAR        [   ---]
   real                :: fnir_beam_act        ! Actual beam fraction - NIR        [   ---]
   real                :: fvis_diff_act        ! Actual diffuse fraction - PAR     [   ---]
   real                :: fnir_diff_act        ! Actual diffuse fraction - NIR     [   ---]
   real                :: ratio                ! Ratio between obs. and expected   [   ---]
   real                :: aux_par              ! Auxiliary variable                [   ---]
   real                :: aux_nir              ! Auxiliary variable                [   ---]
   real                :: secz                 ! sec(zenith distance)              [   ---]
   real                :: log10secz            ! log10[sec(zenith distance)]       [   ---]
   real                :: w10                  ! Minimum Water absorption          [  W/m?]
   !---------------------------------------------------------------------------------------!
   !    Local constants.                                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Extinction coefficient. (equations 1 and 4 of WN85) -----------------------------!
   real              , parameter :: par_beam_expext  = -0.185
   real              , parameter :: nir_beam_expext  = -0.060
   !----- This is the typical conversion of diffuse radiation in sunny days. --------------!
   real              , parameter :: par2diff_sun = 0.400
   real              , parameter :: nir2diff_sun = 0.600
   !----- Coefficients for various equations in WN85. -------------------------------------!
   real, dimension(3), parameter :: wn85_06 = (/ -1.1950, 0.4459, -0.0345 /)
   real, dimension(2), parameter :: wn85_11 = (/  0.90, 0.70 /)
   real, dimension(2), parameter :: wn85_12 = (/  0.88, 0.68 /)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     First thing to check is whether this is daytime or "night-time".  If the zenith   !
   ! angle is too close to horizon, we assume it's dawn/dusk and all radiation goes to     !
   ! diffuse.                                                                              !
   !---------------------------------------------------------------------------------------!
   if (cosz <= rad_cosz_min) then
      rshort_diff  = rshort_full
      par_beam     = 0.0
      nir_beam     = 0.0
      par_diff     = ref_fvis_diff * rshort_diff
      nir_diff     = ref_fnir_diff * rshort_diff
   else
      !----- Save 1/cos(zen), which is the secant.  We will use this several times. -------!
      secz      = 1. / cosz
      log10secz = log10(secz)
      !------------------------------------------------------------------------------------!


      !----- Total radiation at the top [  W/m?], using ED defaults. ----------------------!
      par_beam_top = ref_fvis_beam * solar
      nir_beam_top = ref_fnir_beam * solar
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the potential PAR components (beam, diffuse, total), using equations 1, 3, !
      ! and 9 of WN85.                                                                     !
      !------------------------------------------------------------------------------------!
      par_beam_pot = par_beam_top                                                          &
                   * exp ( par_beam_expext * (atm_prss / prefsea) * secz) * cosz
      par_diff_pot = par2diff_sun * (par_beam_top - par_beam_pot) * cosz
      par_full_pot = par_beam_pot + par_diff_pot
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6. !
      !------------------------------------------------------------------------------------!
      w10 = solar * 10 ** (wn85_06(1) + log10secz * (wn85_06(2) + wn85_06(3) * log10secz))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the potential direct and diffuse near-infrared radiation, using equations !
      ! 4, 5, and 10 of WN85.                                                              !
      !------------------------------------------------------------------------------------!
      nir_beam_pot = ( nir_beam_top                                                        &
                     * exp ( nir_beam_expext * (atm_prss / prefsea) * secz) - w10 ) * cosz
      nir_diff_pot = nir2diff_sun * ( nir_beam_top - nir_beam_pot - w10 ) * cosz
      nir_full_pot = nir_beam_pot + nir_diff_pot
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the actual total for PAR and NIR, using equations 7 and 8.                !
      !------------------------------------------------------------------------------------!
      ratio    = rshort_full / (par_full_pot + nir_full_pot)
      par_full = ratio * par_full_pot
      nir_full = ratio * nir_full_pot
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the fraction of PAR and NIR that stays as beam, using equations 11 and 12 !
      ! of WN85.                                                                           !
      !------------------------------------------------------------------------------------!
      !----- Make sure that the ratio is bounded. -----------------------------------------!
      aux_par       = min(wn85_11(1), max(0., ratio))
      aux_nir       = min(wn85_12(1), max(0., ratio))
      fvis_beam_act = min(1., max(0., par_beam_pot                                         &
                                 * (1. - ((wn85_11(1) - aux_par)/wn85_11(2)) ** twothirds) &
                                 / par_full_pot ) )
      fnir_beam_act = min(1., max(0., nir_beam_pot                                         &
                                 * (1. - ((wn85_12(1) - aux_nir)/wn85_12(2)) ** twothirds) &
                                 / nir_full_pot ) )
      fvis_diff_act = 1. - fvis_beam_act
      fnir_diff_act = 1. - fnir_beam_act
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the radiation components.                                                 !
      !------------------------------------------------------------------------------------!
      par_beam    = fvis_beam_act * par_full
      par_diff    = fvis_diff_act * par_full
      nir_beam    = fnir_beam_act * nir_full
      nir_diff    = fnir_diff_act * nir_full
      rshort_diff = par_diff + nir_diff
      !------------------------------------------------------------------------------------!
   end if

   return
end subroutine wn85_rshort_split
!==========================================================================================!
!==========================================================================================!
