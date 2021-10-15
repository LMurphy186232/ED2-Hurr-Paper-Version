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
!> The basic model for the cumulative probability of damage level j for stem i of species s
!> is logit(p<sub>isj</sub>) = a<sub>j</sub> + c<sub>s</sub> * S * dbh ^<sup>b<sub>s</sub></sup>
!> where S is the storm intensity.
!>
!> The cohort is treated as an ensemble of individuals; some of whom are lightly damaged,
!> some medium damaged, and some heavily damaged. Of the heavily damaged individuals, some
!> will die.
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
      use hurricane_coms      , only : include_hurricanes          &
                                     , hurricane_db_list           &
                                     , hurricane_db_list_len       &
                                     , agb_to_height_lut           &
                                     , height_for_agb_lut          &
                                     , ge_mask_hurr_lut            &
                                     , le_mask_hurr_lut            &
                                     , hurricane_report
      use update_derived_utils, only : update_cohort_derived_props
      use ed_misc_coms        , only : current_time
      use ed_state_vars       , only : edtype                      & ! structure
                                     , polygontype                 & ! structure
                                     , sitetype                    & ! structure
                                     , patchtype                   ! ! structure
      use pft_coms            , only : c2n_leaf                    & ! intent(in)
                                     , c2n_storage                 & ! intent(in)
                                     , c2n_stem                    & ! intent(in)
                                     , l2n_stem                    & ! intent(in)
                                     , qsw                         &
                                     , qbark                       &
                                     , hurr_a1                     &
                                     , hurr_a2                     &
                                     , hurr_b                      &
                                     , hurr_c                      &
                                     , hurr_g                      &
                                     , hurr_h                      &
                                     , min_hurr_dbh                &
                                     , is_grass                    & ! intent(in)
                                     , is_liana                    &
                                     , agf_bs                      & ! intent(in)
                                     , f_labile_leaf               & ! intent(in)
                                     , f_labile_stem               ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap               & ! function
                                     , update_veg_energy_cweh      ! ! function
      use plant_hydro         , only : rwc2tw                      & ! subroutine
                                     , twi2twe                     ! ! subroutine
      use detailed_coms       , only : idetailed
      use allometry           , only : size2bd                     &
                                     , size2bl                     &
                                     , h2dbh

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
      real                          :: bleaf_in
      real                          :: broot_in
      real                          :: bsapwooda_in
      real                          :: bsapwoodb_in
      real                          :: bbarka_in
      real                          :: bbarkb_in
      real                          :: balive_in
      real                          :: bdeada_in
      real                          :: bdeadb_in
      real                          :: bstorage_in
      real                          :: bleaf_loss
      real                          :: broot_loss
      real                          :: bsapwooda_loss
      real                          :: bsapwoodb_loss
      real                          :: bbarka_loss
      real                          :: bbarkb_loss
      real                          :: balive_loss
      real                          :: bdeada_loss
      real                          :: bdeadb_loss
      real                          :: bstorage_loss
      real                          :: a_bfast
      real                          :: b_bfast
      real                          :: a_bstruct
      real                          :: b_bstruct
      real                          :: a_bstorage
      real                          :: b_bstorage
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: old_leaf_water
      real                          :: old_wood_water
      real                          :: old_leaf_water_im2
      real                          :: old_wood_water_im2
      real                          :: nplant_in
      real                          :: nplant_loss
      real                          :: dbh_aim
      real                          :: severity
      real                          :: prob_light
      real                          :: prob_med
      real                          :: prob_heavy
      real                          :: prob_mort
      real                          :: agb_loss
      real                          :: new_agb
      real                          :: amt
      real                          :: new_hite
      real                          :: hite_in
      real                          :: h1
      real                          :: h2
      real                          :: agb1
      real                          :: agb2
      real                          :: agb_in
      integer                       :: ilwr
      integer                       :: iupr
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

      print_detailed = btest(idetailed,6)
      if (print_detailed) then
        open (unit=79,file=trim(hurricane_report),form='formatted',access='append'         &
             ,status='old')

        write (unit=79,fmt='(a,i4,a,i4,a,f12.5)')  'Hurricane occurring in '               &
             ,current_time%month, '/', current_time%year, ' of severity ', severity

        write (unit=79,fmt='(a)') 'PFT, DBH, Beginning height, Ending height, Beginning AGB, Ending AGB, Beginning NPLANT, Ending NPLANT'

      end if

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

                  !----- Hurricanes don't affect grass and small cohorts ------------------!
                  if (.not. is_grass(ipft) .and.                                           &
                      .not. is_liana(ipft) .and.                                           &
                      cpatch%dbh(ico)      .ge. min_hurr_dbh) then


                     !---------------------------------------------------------------------!
                     !      Save original cohort biomass fractions. We will use these      !
                     ! later to determine how much biomass was removed from where, so we   !
                     ! can correctly assign litter                                         !
                     !---------------------------------------------------------------------!
                     bdeada_in       = cpatch%bdeada          (ico)
                     bdeadb_in       = cpatch%bdeadb          (ico)
                     bleaf_in        = cpatch%bleaf           (ico)
                     broot_in        = cpatch%broot           (ico)
                     bsapwooda_in    = cpatch%bsapwooda       (ico)
                     bsapwoodb_in    = cpatch%bsapwoodb       (ico)
                     bbarka_in       = cpatch%bbarka          (ico)
                     bbarkb_in       = cpatch%bbarkb          (ico)
                     bstorage_in     = cpatch%bstorage        (ico)

                     !---------------------------------------------------------------------!
                     !      Save original heat capacity and water content for both leaves  !
                     ! and wood.  These are used to track changes in energy and water      !
                     ! storage due to vegetation dynamics.                                 !
                     !---------------------------------------------------------------------!
                     old_leaf_hcap      = cpatch%leaf_hcap     (ico)
                     old_wood_hcap      = cpatch%wood_hcap     (ico)
                     old_leaf_water     = cpatch%leaf_water    (ico)
                     old_wood_water     = cpatch%wood_water    (ico)
                     old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
                     old_wood_water_im2 = cpatch%wood_water_im2(ico)


                     !---------------------------------------------------------------------!
                     !       Calculate storm damage probabilities                          !
                     !---------------------------------------------------------------------!
                     !----- Probability of light damage -----------------------------------!
                     prob_light = exp(hurr_a1(ipft) + hurr_c(ipft) * severity              &
                                                    * cpatch%DBH(ico) ** hurr_b(ipft))     &
                          /  (1 + exp(hurr_a1(ipft) + hurr_c(ipft) * severity              &
                                                    * cpatch%DBH(ico) ** hurr_b(ipft)))

                     !----- Probability of light plus medium damage -----------------------!
                     prob_med = exp(hurr_a2(ipft) + hurr_c(ipft) * severity                &
                                                    * cpatch%DBH(ico) ** hurr_b(ipft))     &
                        /  (1 + exp(hurr_a2(ipft) + hurr_c(ipft) * severity                &
                                                    * cpatch%DBH(ico) ** hurr_b(ipft)))

                     !----- Subtract out light damage to make it medium alone -------------!
                     prob_med = prob_med - prob_light

                     !----- Make sure a very small medium probability doesn't end negative !
                     if (prob_med .lt. 0) then
                       prob_med = 0
                     end if

                     !----- Heavy damage probability is what's left -----------------------!
                     prob_heavy = 1 - (prob_light + prob_med)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !       Calculate amount of AGB to remove                             !
                     ! Use random numbers to choose final amount to remove                 !
                     !---------------------------------------------------------------------!
                     agb_loss = 0

                     !----- Light damage: biomass loss between 0 and 3.8 percent ----------!
                     agb_loss = agb_loss + (prob_light * (rand() * 0.038))

                     !----- Medium damage: loss between 20 and 53.8% ----------------------!
                     amt = (rand() * (0.538 - 0.2)) + 0.2
                     agb_loss = agb_loss + (prob_med   * amt)

                     !----- Heavy damage: loss between 51.3 and 83% -----------------------!
                     amt = (rand() * (0.83 - 0.513)) + 0.513
                     agb_loss = agb_loss + (prob_heavy * amt)
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !       Calculate storm mortality                                     !
                     ! Our function is the probability that a heavily damaged tree         !
                     ! survives; flip to mortality                                         !
                     !---------------------------------------------------------------------!
                     prob_mort = exp(hurr_g(ipft) + hurr_h(ipft) * cpatch%dbh(ico)) /      &
                            (1 + exp(hurr_g(ipft) + hurr_h(ipft) * cpatch%dbh(ico)))
                     prob_mort = 1.0 - prob_mort
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !       Apply storm mortality                                         !
                     ! The amount of NPLANT lost is the fraction that is heavily damaged   !
                     ! times the mortality probability for heavy damage                    !
                     !---------------------------------------------------------------------!
                     nplant_in = cpatch%nplant(ico)
                     cpatch%nplant(ico) = cpatch%nplant(ico)                               &
                                        * (1 - (prob_heavy * prob_mort))
                     nplant_loss = nplant_in - cpatch%nplant(ico)
                     !---------------------------------------------------------------------!




                     !---------------------------------------------------------------------!
                     !       Apply storm damage                                            !
                     ! Use our lookup table to determine how much height is lost based     !
                     ! on amount of AGB loss                                               !
                     !---------------------------------------------------------------------!

                     !----- Current AGB ---------------------------------------------------!
                     agb_in = cpatch%bleaf (ico) + cpatch%bbarka   (ico)                   &
                            + cpatch%bdeada(ico) + cpatch%bsapwooda(ico)

                     !----- New AGB -------------------------------------------------------!
                     new_agb = agb_in * (1 - agb_loss)

                     !----- Get the new height corresponding to this AGB ------------------!
                     hite_in = cpatch%hite(ico)
                     le_mask_hurr_lut(:) = new_agb <= agb_to_height_lut(:,ipft)
                     ge_mask_hurr_lut(:) = new_agb >= agb_to_height_lut(:,ipft)
                     ilwr  = maxloc (agb_to_height_lut(:,ipft),dim=1,mask=le_mask_hurr_lut)
                     iupr  = minloc (agb_to_height_lut(:,ipft),dim=1,mask=ge_mask_hurr_lut)

                     !----- Linear interpolate if necessary -------------------------------!
                     if (ilwr == iupr) then
                        new_hite = height_for_agb_lut(ilwr,ipft)
                     else
                        h2   = agb_to_height_lut (iupr, ipft)
                        h1   = agb_to_height_lut (ilwr, ipft)
                        agb2 = agb_to_height_lut (iupr, ipft)
                        agb1 = agb_to_height_lut (ilwr, ipft)

                        new_hite = agb1 + (new_agb - agb1) * ((h2 - h1) / (agb2 - agb1))
                     end if

                     cpatch%bleaf(ico) = cpatch%bleaf(ico) * 0.9

                     !----- Don't let height grow -----------------------------------------!
                     if (new_hite < cpatch%hite(ico)) then

                        cpatch%hite(ico) = new_hite

                        !----- What is the dbh that matches this new height? --------------!
                        dbh_aim = h2dbh(cpatch%hite(ico), ipft)

                        !----- Use these size measures to get other ABG values ------------!
                        cpatch%bdeada   (ico) = size2bd(dbh_aim, cpatch%hite(ico), ipft) * agf_bs(ipft)
                        cpatch%bleaf    (ico) = size2bl(dbh_aim, cpatch%hite(ico), cpatch%sla(ico), ipft)
                        cpatch%bbarka   (ico) = cpatch%bleaf(ico) * qbark(ipft) * cpatch%hite(ico) * agf_bs(ipft)
                        cpatch%bsapwooda(ico) = cpatch%bleaf(ico) * qsw  (ipft) * cpatch%hite(ico) * agf_bs(ipft)

                        !----- Have the cohort update itself ------------------------------!
                        call update_cohort_derived_props(cpatch,ico,cpoly%lsl(isi),.false. &
                                                        ,cpoly%llspan_toc(ipft,isi)        &
                                                        ,cpoly%vm_bar_toc(ipft,isi)        &
                                                        ,cpoly%rd_bar_toc(ipft,isi)        &
                                                        ,cpoly%sla_toc   (ipft,isi) )
                        !------------------------------------------------------------------!
                     end if



                     !---------------------------------------------------------------------!
                     !      Get the amount of each biomass fraction lost, to send to       !
                     ! litter. Biomass fractions are reported in kgC / plant. nplant is in !
                     ! plants / m2. Multiplying by nplant gives us kgC / m2, which is the  !
                     ! units of the litter pools                                           !
                     !---------------------------------------------------------------------!
                     !----- Amount lost due to storm damage - original amount minus current!
                     if (new_hite < cpatch%hite(ico)) then
                        bleaf_loss     = (bleaf_in     - cpatch%bleaf    (ico)) * cpatch%nplant(ico)
                        bdeada_loss    = (bdeada_in    - cpatch%bdeada   (ico)) * cpatch%nplant(ico)
                        bdeadb_loss    = (bdeadb_in    - cpatch%bdeadb   (ico)) * cpatch%nplant(ico)
                        broot_loss     = (broot_in     - cpatch%broot    (ico)) * cpatch%nplant(ico)
                        bbarka_loss    = (bbarka_in    - cpatch%bbarka   (ico)) * cpatch%nplant(ico)
                        bbarkb_loss    = (bbarkb_in    - cpatch%bbarkb   (ico)) * cpatch%nplant(ico)
                        bsapwooda_loss = (bsapwooda_in - cpatch%bsapwooda(ico)) * cpatch%nplant(ico)
                        bsapwoodb_loss = (bsapwoodb_in - cpatch%bsapwoodb(ico)) * cpatch%nplant(ico)
                     else
                        bleaf_loss     = 0.0
                        bdeada_loss    = 0.0
                        bdeadb_loss    = 0.0
                        broot_loss     = 0.0
                        bbarka_loss    = 0.0
                        bbarkb_loss    = 0.0
                        bsapwooda_loss = 0.0
                        bsapwoodb_loss = 0.0
                     end if


                     !----- Add amount lost to storm mortality -----------------------------!
                     bleaf_loss     = bleaf_loss     + (bleaf_in     * nplant_loss)
                     bdeada_loss    = bdeada_loss    + (bdeada_in    * nplant_loss)
                     bdeadb_loss    = bdeadb_loss    + (bdeadb_in    * nplant_loss)
                     broot_loss     = broot_loss     + (broot_in     * nplant_loss)
                     bbarka_loss    = bbarka_loss    + (bbarka_in    * nplant_loss)
                     bbarkb_loss    = bbarkb_loss    + (bbarkb_in    * nplant_loss)
                     bsapwooda_loss = bsapwooda_loss + (bsapwooda_in * nplant_loss)
                     bsapwoodb_loss = bsapwoodb_loss + (bsapwoodb_in * nplant_loss)
                     ! There shouldn't be a change of storage due to damage - only mortality
                     bstorage_loss  =                   bstorage_in  * nplant_loss
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !      Split biomass lost into fast and slow to add to appropriate    !
                     ! litter pools.                                                       !
                     !---------------------------------------------------------------------!
                     !---------------------------------------------------------------------!
                     !      Split biomass components that are labile or structural. -------!
                     a_bfast    = f_labile_leaf(ipft) * bleaf_loss                         &
                                + f_labile_stem(ipft)                                      &
                                * (bsapwooda_loss + bbarka_loss + bdeada_loss)
                     b_bfast    = f_labile_leaf(ipft) * broot_loss                         &
                                + f_labile_stem(ipft)                                      &
                                * (bsapwoodb_loss + bbarkb_loss + bdeadb_loss)
                     a_bstruct  = (1.0 - f_labile_leaf(ipft)) * bleaf_loss                 &
                                + (1.0 - f_labile_stem(ipft))                              &
                                * (bsapwooda_loss + bbarka_loss + bdeada_loss)
                     b_bstruct  = (1.0 - f_labile_leaf(ipft)) * broot_loss                 &
                                + (1.0 - f_labile_stem(ipft))                              &
                                * (bsapwoodb_loss + bbarkb_loss + bdeadb_loss)
                     a_bstorage =        agf_bs(ipft)  * bstorage_loss
                     b_bstorage = (1.0 - agf_bs(ipft)) * bstorage_loss
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !       Update our reporting variables                                !
                     !---------------------------------------------------------------------!
                     nplant_removed(ipft) = nplant_removed(ipft) + nplant_loss
                     agb_removed(ipft)    = agb_removed(ipft)                              &
                                          + (bleaf_loss + bdeada_loss + bbarka_loss        &
                                          +  bsapwooda_loss)

                     amt = cpatch%bbarka(ico) + cpatch%bleaf(ico) + cpatch%bdeada(ico)     &
                         + cpatch%bsapwooda(ico)
                     !write (unit=79,fmt='(i4,a,f12.5,a,f12.5,a,f12.5)')                    &
                    !   ipft, ',', cpatch%dbh(ico), ',', hite_in, ','                 &
                    !            ,cpatch%hite(ico)

                    write (unit=79,fmt='(i4,a,f12.5,a,f12.5,a,f12.5,a,f12.5,a,f12.5,a,f12.5,a,f12.5)')&
                             ipft, ',', cpatch%dbh(ico), ',', hite_in, ','                 &
                             , cpatch%hite(ico), ',', agb_in, ','                          &
                             , amt, ',', nplant_in, ',', cpatch%nplant(ico)
                    !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Finalize litter inputs.                                         !
                     !---------------------------------------------------------------------!
                     csite%fgc_in (ipa) = csite%fgc_in(ipa) + a_bfast + a_bstorage
                     csite%fsc_in (ipa) = csite%fsc_in(ipa) + b_bfast + b_bstorage
                     csite%fgn_in (ipa) = csite%fgn_in(ipa)                                &
                                        + a_bfast    / c2n_leaf   (ipft)                   &
                                        + a_bstorage / c2n_storage
                     csite%fsn_in (ipa) = csite%fsn_in(ipa)                                &
                                        + b_bfast    / c2n_leaf   (ipft)                   &
                                        + b_bstorage / c2n_storage
                     csite%stgc_in(ipa) = csite%stgc_in(ipa) + a_bstruct
                     csite%stsc_in(ipa) = csite%stsc_in(ipa) + b_bstruct
                     csite%stgl_in(ipa) = csite%stgl_in(ipa)                               &
                                        + a_bstruct * l2n_stem / c2n_stem(ipft)
                     csite%stsl_in(ipa) = csite%stsl_in(ipa)                               &
                                        + b_bstruct * l2n_stem / c2n_stem(ipft)
                     csite%stgn_in(ipa) = csite%stgn_in(ipa)                               &
                                        + a_bstruct  / c2n_stem   (ipft)
                     csite%stsn_in(ipa) = csite%stsn_in(ipa)                               &
                                        + b_bstruct  / c2n_stem   (ipft)
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !  Update the heat capacity and the vegetation internal energy,       !
                     ! again, following structural growth                                  !
                     !---------------------------------------------------------------------!
                     call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)               &
                                       ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)           &
                                       ,cpatch%nplant(ico),cpatch%pft(ico)                 &
                                       ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
                     call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                 &
                                ,cpatch%bleaf(ico),cpatch%bsapwooda(ico)                   &
                                ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico)                  &
                                ,cpatch%bdeadb(ico),cpatch%broot(ico),cpatch%dbh(ico)      &
                                ,cpatch%pft(ico),cpatch%leaf_water_int(ico),               &
                                cpatch%wood_water_int(ico))
                     call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)    &
                                 ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)            &
                                 ,cpatch%wood_water_im2(ico))
                     call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap &
                                                ,old_leaf_water,old_wood_water             &
                                                ,old_leaf_water_im2,old_wood_water_im2     &
                                                ,.true.,.false.)
                     !---------------------------------------------------------------------!
                   end if  ! Eliminating grasses
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
        !open (unit=79,file=trim(hurricane_report),form='formatted',access='append'          &
        !     ,status='old')

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
      use hurricane_coms      , only : hurricane_db            &
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
!=========================================================================================!
!=========================================================================================!


!=========================================================================================!
!=========================================================================================!
!> \brief Calculates lookup table relating height to AGB.
!> \details Since allometry is strongly dependent on user settings, including those in
!> an XML config file, this will calculate the height-AGB relationship based on how
!> allometry is set up and put the results in a lookup table. This will allow the
!> hurricane module to quickly determine how much height needs to be removed from a
!> cohort to remove a desired amount of AGB.
!-----------------------------------------------------------------------------------------!
subroutine init_hurricane_derived_params_after_xml()
  use allometry     , only : size2bd                     &
                           , size2bl                     &
                           , h2dbh
  use ed_max_dims   , only : n_pft
  use pft_coms      , only : hgt_max                     &
                           , agf_bs                      &
                           , SLA                         &
                           , qbark                       &
                           , qsw
  use hurricane_coms, only : agb_to_height_lut           &
                           , height_for_agb_lut          &
                           , ge_mask_hurr_lut            &
                           , le_mask_hurr_lut            &
                           , nbt_agb_to_height_lut

  !----- Local variables ------------------------------------------------------------------!
  real                    :: height_chunk
  real                    :: height_now
  real                    :: dbh_now
  real                    :: bleaf_now
  real                    :: bbarka_now
  real                    :: bdeada_now
  real                    :: bsapwooda_now
  integer                 :: ipft
  integer                 :: ilut

  !----- Set the number of height divisions, and allocate the table array -----------------!
  nbt_agb_to_height_lut = 100.0
  allocate(agb_to_height_lut (nbt_agb_to_height_lut,n_pft))
  allocate(height_for_agb_lut(nbt_agb_to_height_lut,n_pft))
  allocate(le_mask_hurr_lut  (nbt_agb_to_height_lut))
  allocate(ge_mask_hurr_lut  (nbt_agb_to_height_lut))

  !------ Divide height up to max height, and get AGB for each ----------------------------!
  do ipft=1, n_pft

    !----- How big will the divisions be for this PFT? ------------------------------------!
    height_chunk = hgt_max(ipft) / nbt_agb_to_height_lut
    height_now   = height_chunk

    !----- Calculate AGB for each height increment ----------------------------------------!
    do ilut=1, nbt_agb_to_height_lut
      height_for_agb_lut(ilut, ipft) = height_now

      !----- AGB is leaf, aboveground sapwood, aboveground bdead, and aboveground bark ----!

      !----- What is the dbh that matches this current height? ----------------------------!
      dbh_now = h2dbh(height_now, ipft)

      !----- Calculate the fractions for this size ----------------------------------------!
      bdeada_now    = size2bd(dbh_now, height_now, ipft) * agf_bs(ipft)
      bleaf_now     = size2bl(dbh_now, height_now, SLA(ipft), ipft)
      bbarka_now    = bleaf_now * qbark(ipft) * height_now * agf_bs(ipft)
      bsapwooda_now = bleaf_now * qsw(ipft) * height_now * agf_bs(ipft)

      agb_to_height_lut(ilut, ipft) = bdead_now + bleaf_now + bbarka_now + bsapwooda_now

      height_now = height_now + height_chunk
    end do

  end do

end subroutine init_hurricane_derived_params_after_xml
!==========================================================================================!
!==========================================================================================!


end module hurricane
!==========================================================================================!
!==========================================================================================!
