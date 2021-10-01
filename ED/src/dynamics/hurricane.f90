!==========================================================================================!
!==========================================================================================!
!> \brief Implement scheduled hurricanes with the potential to damage and kill cohorts.
!> \details Hurricanes occur at pre-determined times, input by the user. Depending on their
!> severity, they will remove some amount of biomass, and possibly kill some of the cohort
!> as well.
!
!> A mild storm might remove leaves only. A stronger storm would remove some branches,
!> up to crown snap. To replicate this, the cohort will lose some degree of height.  In
!> conjunction with the tolerance for cohorts to be off allometry, this will simulate
!> damage that a cohort must repair.
!>
!> This module checks monthly to see if a hurricane is scheduled. The hurricane schedule
!> is in a separate file that is read in at start up, and whose name is passed in ED2IN.
!> While it is certainly possible in the real world for two hurricanes to occur in the same
!> month, in this module it is not allowed.
!>
!> Hurricanes cannot be used if ED2 is in bigleaf mode.
!>
!> If IDETAILED = 64, indicating disturbance details are to be written, this will write a
!> detailed file of storm activity called "hurricane_report.txt"
!>
!> Possible future work: randomly generate a storm regime; further decouple biomass pools
!> to allow preferential removal of, say, leaves (right now doing so causes the cohorts
!> to fail budget checks in structural growth, since part of the check is that the
!> pools are on allometry).
!> \author Lora Murphy, September 2021
!> \todo Parameter copying and passing in MPI situations
!> \todo Figure out how to deal with lianas
!------------------------------------------------------------------------------------------!
module hurricane

   !=======================================================================================!
   !=======================================================================================!

   contains

   !=======================================================================================!
   !=======================================================================================!
   !> \brief Main hurricane driver.
   !> \param cgrid Main grid.
   !---------------------------------------------------------------------------------------!
   subroutine apply_hurricane(cgrid)
     use ed_max_dims          , only : n_pft
      use disturb_coms        , only : include_hurricanes          &
                                     , hurricane_db_list           &
                                     , hurricane_db_list_len       &
                                     , hurricane_report
      use ed_misc_coms        , only : current_time
      use ed_state_vars       , only : edtype                      & ! structure
                                     , polygontype                 & ! structure
                                     , sitetype                    & ! structure
                                     , patchtype                   ! ! structure
      use pft_coms            , only : c2n_leaf                    & ! intent(in)
                                     , c2n_storage                 & ! intent(in)
                                     , c2n_stem                    & ! intent(in)
                                     , l2n_stem                    & ! intent(in)
                                    !, negligible_nplant           & ! intent(in)
                                    !, is_grass                    & ! intent(in)
                                     , agf_bs                      & ! intent(in)
                                    !, q                           & ! intent(in)
                                    !, storage_reflush_times       & ! intent(in)
                                    !, is_liana                    & ! intent(in)
                                    !, cbr_severe_stress           & ! intent(in)
                                    !, h_edge                      & ! intent(in)
                                     , f_labile_leaf               & ! intent(in)
                                     , f_labile_stem               ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap               & ! function
                                     , update_veg_energy_cweh      ! ! function
      use plant_hydro         , only : rwc2tw                      & ! subroutine
                                     , twi2twe                     ! ! subroutine
      use detailed_coms       , only : idetailed
      use allometry           , only : size2bd                     &
                                     , size2bl

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                    , target      :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      real, dimension(n_pft)        :: agb_removed      ! Amount of AGB removed per PFT
      real, dimension(n_pft)        :: nplant_removed   ! Amount of NPLANT removed per PFT

      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ihu
      integer                       :: ipft
      !integer                       :: prev_month
      !integer                       :: prev_year
      !integer                       :: prev_ndays
      !integer                       :: imonth
      !integer                       :: phenstatus_in
      !integer                       :: krdepth_in
      !real                          :: tor_fact
      !real                          :: fgc_in_in
      !real                          :: fsc_in_in
      !real                          :: stgc_in_in
      !real                          :: stsc_in_in
      !real                          :: nplant_loss
      !real                          :: pat_balive_in
      !real                          :: pat_bdead_in
      !real                          :: pat_bstorage_in
      !real                          :: pat_mortality
      !real                          :: pat_carbon_miss
      !real                          :: carbon_miss
      real                          :: bleaf_in
      real                          :: broot_in
      real                          :: bsapwooda_in
      real                          :: bsapwoodb_in
      real                          :: bbarka_in
      real                          :: bbarkb_in
      real                          :: balive_in
      real                          :: bdeada_in
      real                          :: bdeadb_in
      !real                          :: bevery_in
      !real                          :: hite_in
      !real                          :: dbh_in
      !real                          :: nplant_in
      !real                          :: bstorage_in
      !real                          :: bstorage_reserve
      !real                          :: agb_in
      !real                          :: lai_in
      !real                          :: wai_in
      !real                          :: cai_in
      !real                          :: ba_in
      !real                          :: bag_in
      !real                          :: bam_in
      !real                          :: vm_bar_in
      !real                          :: sla_in
      !real                          :: psi_open_in
      !real                          :: psi_closed_in
      !real                          :: cb_act
      !real                          :: cb_lightmax
      !real                          :: cb_moistmax
      !real                          :: cb_mlmax
      !real                          :: cbr_light
      !real                          :: cbr_moist
      !real                          :: cbr_ml
      !real                          :: f_bseeds
      !real                          :: f_bdeada
      !real                          :: f_bdeadb
      !real                          :: f_growth
      !real                          :: f_bstorage
      real                          :: a_bfast_litter
      real                          :: b_bfast_litter
      real                          :: a_bstruct_litter
      real                          :: b_bstruct_litter
      real                          :: a_bstorage_litter
      real                          :: b_bstorage_litter
      real                          :: a_bfast
      real                          :: b_bfast
      real                          :: a_bstruct
      real                          :: b_bstruct
      real                          :: a_bstorage
      real                          :: b_bstorage
      !real                          :: maxh !< maximum patch height
      !real                          :: mort_litter
      !real                          :: bseeds_mort_litter
      !real                          :: net_seed_N_uptake
      !real                          :: net_stem_N_uptake
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: old_leaf_water
      real                          :: old_wood_water
      real                          :: old_leaf_water_im2
      real                          :: old_wood_water_im2
      real                          :: nplant_in
      real                          :: nplant_loss
      real                          :: bleaf_loss
      real                          :: severity
      logical                       :: hurricane_time
      logical                       :: print_detailed
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !       Proceed only if hurricanes are enabled and a hurricane is due this month.    !
      !------------------------------------------------------------------------------------!
      if (include_hurricanes .eq. 0) return
      hurricane_time = .false.

      do ihu=1,hurricane_db_list_len
         if (hurricane_db_list(ihu)%year  .eq. current_time%year  .and. &
             hurricane_db_list(ihu)%month .eq. current_time%month) then

            severity = hurricane_db_list(ihu)%severity
            hurricane_time = .true.
         end if
      end do

      if (.not.hurricane_time) return
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Get ready: announce the hurricane, empty the reporting variables                !
      !------------------------------------------------------------------------------------!
      write (unit=*,fmt='(a,1x,es12.5)')  'Hurricane occurring of severity ', severity
      agb_removed    = 0
      nplant_removed = 0

      !------------------------------------------------------------------------------------!
      !    Loop over everything down to the cohort level.                                  !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               cohortloop: do ico = 1,cpatch%ncohorts

                  !----- Assigning an alias for PFT type. ---------------------------------!
                  ipft    = cpatch%pft(ico)

                  !------------------------------------------------------------------------!
                  !      Save original cohort biomass fractions. We will use these         !
                  ! later to determine how much biomass was removed from where, so we      !
                  ! can correctly assign litter                                            !
                  !------------------------------------------------------------------------!
                  !balive_in       = cpatch%balive          (ico)
                  bdeada_in       = cpatch%bdeada          (ico)
                  bdeadb_in       = cpatch%bdeadb          (ico)
                  bleaf_in        = cpatch%bleaf           (ico)
                  broot_in        = cpatch%broot           (ico)
                  bsapwooda_in    = cpatch%bsapwooda       (ico)
                  bsapwoodb_in    = cpatch%bsapwoodb       (ico)
                  bbarka_in       = cpatch%bbarka          (ico)
                  bbarkb_in       = cpatch%bbarkb          (ico)

                  !------------------------------------------------------------------------!
                  !      Save original heat capacity and water content for both leaves     !
                  ! and wood.  These are used to track changes in energy and water         !
                  ! storage due to vegetation dynamics.                                    !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap      = cpatch%leaf_hcap     (ico)
                  old_wood_hcap      = cpatch%wood_hcap     (ico)
                  old_leaf_water     = cpatch%leaf_water    (ico)
                  old_wood_water     = cpatch%wood_water    (ico)
                  old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
                  old_wood_water_im2 = cpatch%wood_water_im2(ico)



                  !------------------------------------------------------------------------!
                  !       Storm mortality                                                  !
                  !------------------------------------------------------------------------!
                  nplant_in = cpatch%nplant(ico)
                  cpatch%nplant(ico) = cpatch%nplant(ico) * 0.9
                  nplant_loss = nplant_in - cpatch%nplant(ico)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !       Storm damage                                                     !
                  ! We reduce the height of the cohort in an amount corresponding to the   !
                  ! level of damage this storm inflicts.
                  !------------------------------------------------------------------------!
                  hite_in = cpatch%hite(ico)
                  cpatch%hite(ico) = cpatch%hite(ico) * 0.9

                  !----- Have the cohort update itself ------------------------------------!
                  call update_cohort_derived_props(cpatch,ico,cpoly%lsl(isi),.false.       &
                                                  ,cpoly%llspan_toc(ipft,isi)              &
                                                  ,cpoly%vm_bar_toc(ipft,isi)              &
                                                  ,cpoly%rd_bar_toc(ipft,isi)              &
                                                  ,cpoly%sla_toc   (ipft,isi) )
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      Split biomass components that are labile or structural. ----------!
                  a_bfast    = f_labile_leaf(ipft) * cpatch%bleaf(ico)                     &
                             + f_labile_stem(ipft)                                         &
                             * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)                &
                             + cpatch%bdeada   (ico) )
                  b_bfast    = f_labile_leaf(ipft) * cpatch%broot(ico)                     &
                             + f_labile_stem(ipft)                                         &
                             * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)                &
                             + cpatch%bdeadb   (ico) )
                  a_bstruct  = (1.0 - f_labile_leaf(ipft)) * cpatch%bleaf(ico)             &
                             + (1.0 - f_labile_stem(ipft))                                 &
                             * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)                &
                             + cpatch%bdeada   (ico) )
                  b_bstruct  = (1.0 - f_labile_leaf(ipft)) * cpatch%broot(ico)             &
                             + (1.0 - f_labile_stem(ipft))                                 &
                             * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)                &
                             + cpatch%bdeadb   (ico) )
                  a_bstorage =        agf_bs(ipft)  * cpatch%bstorage(ico)
                  b_bstorage = (1.0 - agf_bs(ipft)) * cpatch%bstorage(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Get the amount of each biomass fraction lost, to send to litter   !
                  !------------------------------------------------------------------------!

                  ! LEM: This seems like it wouldn't add up...

                  !----- Amount lost due to storm damage - original amount minus current --!
                  bleaf_loss     = bleaf_in     - cpatch%bleaf    (ico)
                  bdeada_loss    = bdeada_in    - cpatch%bdeada   (ico)
                  bdeadb_loss    = bdeadb_in    - cpatch%bdeadb   (ico)
                  broot_loss     = broot_in     - cpatch%broot    (ico)
                  bbarka_loss    = bbarka_in    - cpatch%bbarka   (ico)
                  bbarkb_loss    = bbarkb_in    - cpatch%bbarkb   (ico)
                  bsapwooda_loss = bsapwooda_in - cpatch%bsapwooda(ico)
                  bsapwoodb_loss = bsapwoodb_in - cpatch%bsapwoodb(ico)

                  !----- Amount lost to storm mortality -----------------------------------!
                  bleaf_loss     = bleaf_loss     + (bleaf_in     * nplant_loss)
                  bdeada_loss    = bdeada_loss    + (bdeada_in    * nplant_loss)
                  bdeadb_loss    = bdeadb_loss    + (bdeadb_in    * nplant_loss)
                  broot_loss     = broot_loss     + (broot_in     * nplant_loss)
                  bbarka_loss    = bbarka_loss    + (bbarka_in    * nplant_loss)
                  bbarkb_loss    = bbarkb_loss    + (bbarkb_in    * nplant_loss)
                  bsapwooda_loss = bsapwooda_loss + (bsapwooda_in * nplant_loss)
                  bsapwoodb_loss = bsapwoodb_loss + (bsapwoodb_in * nplant_loss)

                  !----- Determine the amount of above and belowground structural lost ----!
                  !new_bdead = size2bd(cpatch%dbh(ico), cpatch%hite(ico), ipft)
                  !bdeada_loss = bdeada_in - (agf_bs(ipft) * new_bdead)
                  !bdeadb_loss = bdeadb_in - ((1.0 - agf_bs(ipft)) * new_bdead)

                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !       Update our reporting variables                                   !
                  ! LEM: I feel like this isn't right. Adjust loss for nplant?
                  !------------------------------------------------------------------------!
                  nplant_removed(ipft) = nplant_removed(ipft) +                            &
                                         (nplant_loss * csite%area(ipa))
                  agb_removed(ipft)    = agb_removed(ipft)                                 &
                                       + ((cpatch%agb(ico) * nplant_loss)                  &
                                       + (bleaf_loss + bdeada_loss + bbarka_loss
                                       +  bsapwooda_loss)) * csite%area(ipa)

                  !------------------------------------------------------------------------!
                  !     Add storm-killed trees to litter                                   !
                  ! Use approach from structural growth                                    !
                  !------------------------------------------------------------------------!


                  !----- Multiply by mortality to get litter inputs -----------------------!
                  a_bfast_litter    = a_bfast    * nplant_loss
                  b_bfast_litter    = b_bfast    * nplant_loss
                  a_bstruct_litter  = a_bstruct  * nplant_loss
                  b_bstruct_litter  = b_bstruct  * nplant_loss
                  a_bstorage_litter = a_bstorage * nplant_loss
                  b_bstorage_litter = b_bstorage * nplant_loss

                  !----- Damage removed biomass -------------------------------------------!
                  a_bfast_litter   = a_bfast_litter +                                      &
                                     bleaf_loss * f_labile_leaf(ipft)
                  a_bstruct_litter = a_bstruct_litter +                                    &
                                     bleaf_loss * (1.0 - f_labile_leaf(ipft))

                  !------------------------------------------------------------------------!
                  !     Finalize litter inputs.                                            !
                  !------------------------------------------------------------------------!
                  csite%fgc_in (ipa) = csite%fgc_in(ipa) + a_bfast_litter                  &
                                     + a_bstorage_litter
                  csite%fsc_in (ipa) = csite%fsc_in(ipa) + b_bfast_litter                  &
                                     + b_bstorage_litter
                  csite%fgn_in (ipa) = csite%fgn_in(ipa)                                   &
                                     + a_bfast_litter    / c2n_leaf   (ipft)               &
                                     + a_bstorage_litter / c2n_storage
                  csite%fsn_in (ipa) = csite%fsn_in(ipa)                                   &
                                     + b_bfast_litter    / c2n_leaf   (ipft)               &
                                     + b_bstorage_litter / c2n_storage
                  csite%stgc_in(ipa) = csite%stgc_in(ipa) + a_bstruct_litter
                  csite%stsc_in(ipa) = csite%stsc_in(ipa) + b_bstruct_litter
                  csite%stgl_in(ipa) = csite%stgl_in(ipa)                                  &
                                     + a_bstruct_litter * l2n_stem / c2n_stem(ipft)
                  csite%stsl_in(ipa) = csite%stsl_in(ipa)                                  &
                                     + b_bstruct_litter * l2n_stem / c2n_stem(ipft)
                  csite%stgn_in(ipa) = csite%stgn_in(ipa)                                  &
                                     + a_bstruct_litter  / c2n_stem   (ipft)
                  csite%stsn_in(ipa) = csite%stsn_in(ipa)                                  &
                                     + b_bstruct_litter  / c2n_stem   (ipft)
                  !------------------------------------------------------------------------!

                  ! Terminate cohorts

                  !------------------------------------------------------------------------!
                  !  Update the heat capacity and the vegetation internal energy, again,   !
                  !  following structural growth                                           !
                  !------------------------------------------------------------------------!
                  call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)                  &
                                    ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)              &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
                  call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                    &
                             ,cpatch%bleaf(ico),cpatch%bsapwooda(ico)                      &
                             ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico),cpatch%bdeadb(ico)  &
                             ,cpatch%broot(ico),cpatch%dbh(ico),cpatch%pft(ico)            &
                             ,cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico))
                  call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)       &
                              ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)               &
                              ,cpatch%wood_water_im2(ico))
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap    &
                                             ,old_leaf_water,old_wood_water                &
                                             ,old_leaf_water_im2,old_wood_water_im2        &
                                             ,.true.,.false.)
                  !------------------------------------------------------------------------!
               end do cohortloop
            end do patchloop
         end do siteloop
      end do polyloop
      !---------------------------------------------------------------------------!

      ! Prune lianas

      !------------------------------------------------------------------------------------!
      !       Find out whether to print detailed information on screen.                    !
      !------------------------------------------------------------------------------------!
      print_detailed = btest(idetailed,6)
      if (print_detailed) then
        open (unit=79,file=trim(hurricane_report),form='formatted',access='append'          &
             ,status='old')

        write (unit=79,fmt='(a,i4,a,i4,a,f12.5)')  'Hurricane occurring in '          &
             ,current_time%month, '/', current_time%year, ' of severity ', severity
        write (unit=79,fmt='(a)') 'NPLANT removed:'
        do ipft=1,n_pft
           if (nplant_removed(ipft) > 0) then
              write (unit=79,fmt='(a,i4,a,es12.5)') 'PFT ', ipft, ': ', nplant_removed(ipft)
           end if
        end do
        write (unit=79,fmt='(a,1x)') 'AGB removed:'
        do ipft=1,n_pft
           if (agb_removed(ipft) > 0) then
              write (unit=79,fmt='(a,i4,a,es12.5)') 'PFT ', ipft, ': ', agb_removed(ipft)
           end if
        end do
        close(unit=79,status='keep')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine apply_hurricane
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !> \brief Reads in the file of scheduled hurricanes.
   !> \details The hurricane filename is specified in the ED2IN file as "HURRICANE_DB".
   !> The format is a space-delimited text file with a header and three columns:
   !> year, month, and severity. This throws a fatal error if the file is not found, or
   !> multiple hurricanes are supposed to take place at the same time. A warning is written
   !> if the file is longer than the maximum number of allowed records.
   !> \param hurricane_db Name of hurricane schedule.
   !---------------------------------------------------------------------------------------!
   subroutine read_hurricane_db()
      use ed_max_dims         , only : str_len                 &
                                     , max_hurricanes
      use disturb_coms        , only : hurricane_db            &
                                     , hurricane_db_list       &
                                     , hurricane_db_list_len   &
                                     , hurricane_report
      use detailed_coms       , only : idetailed
      implicit none

      !----- Local variables. -------------------------------------------------------------!
      integer  :: ferr  ! error flag
      integer  :: record_counter
      integer  :: ihu
      logical  :: l1
      logical  :: remove_entry

      !------------------------------------------------------------------------------------!
      !      Make sure the hurricane schedule file exists.                                 !
      !------------------------------------------------------------------------------------!

      inquire(file=trim(hurricane_db),exist=l1)
      if (.not. l1) then
         write (unit=*,fmt='(a)') 'File '//trim(hurricane_db)//' not found!'
         write (unit=*,fmt='(a)') 'Specify HURRICANE_DB properly in ED namelist.'
         call fatal_error('HURRICANE_DB not found!','read_hurricane_times' &
                         ,'hurricane.f90')
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loading the hurricane times                                                   !
      !  File is in the format year month severity                                         !
      !------------------------------------------------------------------------------------!
      open(unit=12,file=trim(hurricane_db),form='formatted',status='old')
      read(unit=12,fmt=*)  ! skip header
      hurricane_db_list_len = 0 ! initialize the length of hurricane schedule
      ferr = 0

      do while ((ferr .eq. 0) .and. (hurricane_db_list_len .lt. max_hurricanes))
         record_counter = hurricane_db_list_len + 1

         !----- Read this hurricane's information -----------------------------------------!
         read(unit=12,fmt=*,iostat=ferr) hurricane_db_list(record_counter)%year     &
                                        ,hurricane_db_list(record_counter)%month    &
                                        ,hurricane_db_list(record_counter)%severity


         !----- Error trapping: can't be the same as an existing hurricane ----------------!
         hurricane_loop: do ihu=1,hurricane_db_list_len
            if (hurricane_db_list(record_counter)%year  .eq. hurricane_db_list(ihu)%year  .and. &
                hurricane_db_list(record_counter)%month .eq. hurricane_db_list(ihu)%month) then

               close(unit=12)
               write (unit=*,fmt='(a)') 'Two hurricanes cannot occur at the same time.'
               call fatal_error('Two hurricanes cannot occur at the same time.','read_hurricane_times' &
                          ,'hurricane.f90')

            end if
         end do hurricane_loop

         if (hurricane_db_list(record_counter)%year     .gt. 0 .and. &
             hurricane_db_list(record_counter)%month    .gt. 0 .and. &
             hurricane_db_list(record_counter)%severity .gt. 0) then
             hurricane_db_list_len = record_counter
         end if
      end do

      ! Close the file
      close(unit=12)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Error trapping                                                                !
      !  Doing this as a separate step so the input file is nicely closed.                 !
      !                                                                                    !
      !  Month must be between 1 and 12; storm severity must be between 0 and 1. It is     !
      !  possible to specify a bad year and hurricanes won't occur. I'm deliberately not   !
      !  trapping for that because someone may have generated a schedule for a long run    !
      !  and this way they won't have to edit it down for shorter test runs.               !
      !------------------------------------------------------------------------------------!
      hurricane_error_trap: do ihu=1,hurricane_db_list_len

         !----- Error trapping: month must be between 1 and 12 ----------------------------!
         if (hurricane_db_list(ihu)%month < 1   .or. &
             hurricane_db_list(ihu)%month > 12) then

            write (unit=*,fmt='(a,1x,i12,1x)')  'Cannot understand hurricane month value: ', &
                            hurricane_db_list(ihu)%month, &
                            'Hurricane months must be between 1 and 12.'
            call fatal_error('Cannot understand hurricane month value.','read_hurricane_times' &
                             ,'hurricane.f90')
         end if

         !----- Error trapping: storm severity must be between 0 and 1. --------------------!
         if (hurricane_db_list(ihu)%severity < 0  .or. &
             hurricane_db_list(ihu)%severity > 1) then

            write (unit=*,fmt='(a,1x,i12,1x)')  'Cannot understand hurricane severity value: ', &
                 hurricane_db_list(ihu)%severity, '. Hurricane severity must be between 0 and 1.'
            call fatal_error('Cannot understand hurricane severity value.','read_hurricane_times' &
                             ,'hurricane.f90')
         end if
      end do hurricane_error_trap
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      If we had too many observations, warn the user                                !
      !------------------------------------------------------------------------------------!
      if (hurricane_db_list_len .eq. max_hurricanes) then
         write (unit=*,fmt='(a,i4,a)') &
               'Too many entries in input HURRICANE_DB. Using only the first ', max_hurricanes, '...'
      end if
      !------------------------------------------------------------------------------------!
      !      If detailed output is requested, set up the hurricane report file             !
      !------------------------------------------------------------------------------------!
      hurricane_report = 'hurricane_report.txt'
      if (btest(idetailed,6)) then
        open  (unit=79,file=trim(hurricane_report),form='formatted',status='replace')
        write (unit=79,fmt='(a,1x)') 'Hurricane report'
        close (unit=79,status='keep')
      end if

end subroutine read_hurricane_db
!=======================================================================================!
!=======================================================================================!


end module hurricane
!==========================================================================================!
!==========================================================================================!
