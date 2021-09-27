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
!> Possible future work: randomly generate a storm regime.
!> \author Lora Murphy, September 2021
!> \todo Check include hurricane flag when executing
!> \todo Parameter copying and passing in MPI situations
!> \todo Flag error if in big leaf situations
!> \todo Figure out how to deal with lianas
!> \todo ed_opspec_misc for validating the hurricane choice (0 or 1)
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
      use ed_state_vars       , only : edtype                     & ! structure
                                     , polygontype                & ! structure
                                     , sitetype                   & ! structure
                                     , patchtype                  ! ! structure
      use detailed_coms       , only : idetailed

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                    , target      :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      logical                                       :: hurricane_time
      logical                                       :: print_detailed
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !       Proceed only if a hurricane is due this month.                               !
      !------------------------------------------------------------------------------------!
      hurricane_time = .false.
      if (.not.hurricane_time) return

      ! No big leaf
      ! Prune lianas

      !------------------------------------------------------------------------------------!
      !       Find out whether to print detailed information on screen.                    !
      !------------------------------------------------------------------------------------!
      print_detailed = btest(idetailed,6)
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
      use ed_max_dims           , only : str_len                 &
                                       , max_hurricanes
      use disturb_coms          , only : hurricane_db            &
                                       , hurricane_db_list       &
                                       , hurricane_db_list_len
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

end subroutine read_hurricane_db
!=======================================================================================!
!=======================================================================================!


end module hurricane
!==========================================================================================!
!==========================================================================================!
