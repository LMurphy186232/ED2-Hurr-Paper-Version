!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This module contains some variables used in LEAF-3.                                   !
!------------------------------------------------------------------------------------------!
module leaf_coms

   use grid_dims
   use rconstants, only: grav        & ! intent(in)
                       , vonk        & ! intent(in)
                       , alvl3       & ! intent(in)
                       , onethird    & ! intent(in)
                       , twothirds   & ! intent(in)
                       , umol_2_mol  & ! intent(in)
                       , t00         & ! intent(in)
                       , cliq        & ! intent(in)
                       , mmcod       & ! intent(in)
                       , prefsea     & ! intent(in)
                       , mmo2        & ! intent(in)
                       , mmdryi      & ! intent(in)
                       , Watts_2_Ein ! ! intent(in)

   !----- Parameters that are initialised from RAMSIN. ------------------------------------! 
   real    :: ugbmin          ! Minimum leaf-level velocity                     [      m/s]
   real    :: ubmin           ! Minimum velocity                                [      m/s]
   real    :: ustmin          ! Minimum ustar                                   [      m/s]
   real    :: gamm            ! Gamma used by Businger et al. (1971) - momentum.
   real    :: gamh            ! Gamma used by Businger et al. (1971) - heat.
   real    :: tprandtl        ! Turbulent Prandtl number.
   real    :: ribmax          ! Maximum bulk Richardson number
   real    :: leaf_maxwhc     ! Leaf maximum water holding capacity             [kg/m2leaf]
   real    :: min_patch_area  !  Minimum patch area to consider
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Commons used by LEAF-3.                                                           !
   !---------------------------------------------------------------------------------------!
   integer :: niter_leaf3      ! ! number of LEAF-3 timesteps
   real    :: ndtvegi          ! ! inverse of the vegetation time steps
   integer :: flag_sfcwater    ! ! flag to determine the pounding water stability.
   logical :: resolvable       ! ! Flag to determine whether to resolve vegetation or not.

   real    :: dtl3             & ! LEAF-3 time step for most pools (ground, CAS, snow)
            , dtl3_factor      & ! Timestep factor for most pools (dtgc/dtlong)
            , dtvg             & ! LEAF-3 time step for vegetation
            , dtl3owcc         & ! leaf timestep   / (can_depth * can_rhos)
            , dtl3ohcc         & ! leaf timestep   / (can_depth * can_rhos)
            , dtl3occc         & ! mmdry * leaf timestep   / (can_depth * can_rhos)

            , atm_up           & ! U velocity at top of surface layer            [     m/s]
            , atm_vp           & ! V velocity at top of surface layer            [     m/s]
            , atm_thil         & ! ice-liquid pot. temp. at top of surface layer [       K]
            , atm_theta        & ! potential temperature at top of surface layer [       K]
            , atm_temp         & ! temperature at top of surface layer           [       K]
            , atm_temp_zcan    & ! air temperature just above the canopy air     [       K]
            , atm_enthalpy     & ! specific enthalpy above the canopy air space  [    J/kg]
            , atm_rvap         & ! vapor mixing ratio at top of surface layer    [   kg/kg]
            , atm_rtot         & ! total mixing ratio at top of surface layer    [   kg/kg]
            , atm_shv          & ! specific humidity at top of surface layer     [   kg/kg]
            , atm_co2          & ! CO2 mixing ratio at top of surface layer      [?mol/mol]
            , atm_theiv        & ! atmospheric ice-vapour equiv. potential temp. [       K]
            , atm_vpdef        & ! atmospheric vapour pressure deficit           [      Pa]
            , atm_rhos         & ! air density                                   [   kg/m?]
            , geoht            & ! height at top of surface layer                [       m]
            , atm_exner        & ! "Exner" function at surface (Exner/cp)        [     ---]
            , atm_prss         & ! pressure at surface                           [      Pa]
            , atm_vels         & ! wind speed at top of surface layer            [       m]
            , pcpgl            & ! precip mass in leaf timestep                  [   kg/m?]
            , qpcpgl           & ! precip energy in leaf timestep                [    J/m?]
            , dpcpgl           & ! precip depth in leaf timestep                 [       m]

            , snowfac          & ! frac. of veg. height covered by sfcwater      [     ---]
            , can_exner        & ! canopy air Exner function                     [  J/kg/K]
            , can_temp         & ! canopy air temperature                        [       K]
            , can_enthalpy     & ! canopy air specific enthalpy                  [    J/kg]
            , can_cp           & ! canopy air specific heat at constant pressure [  J/kg/K]
            , can_rhos         & ! canopy air density                            [   kg/m?]
            , can_rsat         & ! canopy air saturation mixing ratio            [   kg/kg]
            , can_shv          & ! canopy air specific humidity                  [   kg/kg]
            , can_rhv          & ! canopy air relative humidity                  [     ---]
            , can_depth        & ! canopy depth                                  [       m]
            , veg_temp         & ! vegetation temperature                        [       K]
            , veg_fliq         & ! liquid fraction of vegetation surface water   [      --]
            , veg_wind         & ! wind speed at the vegetation height           [     m/s]
            , estar            & ! theta_Eiv characteristic friction scale       [       K]
            , qstar            & ! spec. hum. characteristic friction scale      [   kg/kg]
            , timefac_sst      & ! time interpolation factor for SST             [     ---]
            , timefac_ndvi     & ! time interpolation factor for NDVI or LAI     [     ---]
            
            , gbh              & ! leaf boundary layer heat conductance          [J/K/m2/s]
            , gbw              & ! leaf boundary layer water conductance         [ kg/m2/s]
            , gsw              & ! stomatal conductance                          [ kg/m2/s]
            , ggnet            & ! net ground heat/water conductance             [     m/s]
            , ggbare           & ! heat/water conductance, bare ground           [     m/s]
            , ggveg            & ! heat/water conductance,  vegetated ground     [     m/s]
            , ggsoil           & ! water conductance only, soil moisture effect  [     m/s]
            , rho_ustar        & ! canopy density time friction velocity
            , rshort_g         & ! net SW radiation absorbed by grnd
            , rshort_v         & ! net SW radiation absorbed by veg
            , rshort_a         & ! net SW rad. reflected to atm by veg +grnd
            , rlong_v          & ! net LW rad. absorbed by veg
            , rlong_g          & ! net LW rad. absorbed by ground
            , rlong_s          & ! net LW rad. absorbed by snow

            , hflxsc           & ! sensible heat from snowpack to canopy        [   J/m?/s]
            , wflxsc           & ! water vapor from snowpack to canopy          [  kg/m?/s]
            , qwflxsc          & ! latent heat from snowpack to canopy          [   J/m?/s]
            , hflxgc           & ! sensible heat from ground to canopy          [   J/m?/s]
            , wflxgc           & ! water vapor from ground to canopy            [  kg/m?/s]
            , qwflxgc          & ! latent heat from ground to canopy            [   J/m?/s]
            , cflxgc           & ! carbon from ground to canopy                 [?mol/m?/s]
            , eflxac           & ! enthalpy flux from atmosphere to canopy      [   J/m?/s]
            , hflxac           & ! sens. heat flux from atmosphere to canopy    [   J/m?/s]
            , wflxac           & ! water flux from atmosphere to canopy         [  kg/m?/s]
            , qwflxac          & ! water flux from atmosphere to canopy         [   J/m?/s]
            , cflxac           & ! carbon flux from atmosphere to canopy        [?mol/m?/s]
            , hflxvc           & ! sensible heat from vegetation to canopy      [   J/m?/s]
            , wflxvc           & ! water vapor from veg. to canopy (evap.)      [  kg/m?/s]
            , qwflxvc          & ! latent heat from veg. to canopy (evap.)      [   J/m?/s]
            , cflxvc           & ! carbon from vegetation to canopy             [?mol/m?/s]
            , gpp              & ! water flux due to transpiration              [  kg/m?/s]
            , transp           & ! water flux due to transpiration              [  kg/m?/s]
            , qtransp          & ! latent heat flux due to transpiration        [  kg/m?/s]
            , wshed            & ! mass flux due to canopy dripping             [  kg/m?/s]
            , qwshed           & ! energy flux due to canopy dripping           [   J/m?/s]
            , dwshed           & ! depth flux due to canopy dripping            [      m/s]
            , hflxvc_tot       & ! sensible heat from vegetation to canopy      [   J/m?/s]
            , wflxvc_tot       & ! water vapor from veg. to canopy (evap.)      [  kg/m?/s]
            , qwflxvc_tot      & ! latent heat from veg. to canopy (evap.)      [   J/m?/s]
            , cflxvc_tot       & ! carbon from vegetation to canopy             [?mol/m?/s]
            , cflxgc_tot       & ! carbon from vegetation to canopy             [?mol/m?/s]
            , transp_tot       & ! water flux due to transpiration              [  kg/m?/s]
            , qtransp_tot      & ! latent heat flux due to transpiration        [   J/m?/s]
            , wshed_tot        & ! water shed from vegetation to ground         [  kg/m?/s]
            , qwshed_tot       & ! energy from shed water                       [   J/m?/s]
            , dwshed_tot       & ! total depth of water shed                    [      m/s]
            , throughfall_tot  & ! Total throughfall water                      [  kg/m?/s]
            , qthroughfall_tot & ! Total throughfall internal energy            [  kg/m?/s]
            , dthroughfall_tot & ! Total throughfall depth                      [  kg/m?/s]
            , intercepted_tot  & ! Total intercepted water                      [  kg/m?/s]
            , qintercepted_tot & ! Total intercepted internal energy            [  kg/m?/s]
            , dintercepted_tot & ! Total intercepted depth                      [  kg/m?/s]
            , dewgnd_tot       & ! dew formation on ground                      [  kg/m?/s]
            , qdewgnd_tot      & ! energy from dew formation on ground          [   J/m?/s]
            , ddewgnd_tot      & ! depth gain from dew formation on ground      [      m/s]
            , virtual_energy   & ! "virtual layer" internal energy              [     J/m?]
            , virtual_water    & ! "virtual layer" water mass                   [    kg/m?]
            , virtual_depth    & ! "virtual layer" depth                        [        m]

            , emis_town        & ! Urban emissivity (for TEB)                   [     ----]
            , alb_town         & ! Urban albedo (for TEB)                       [     ----]
            , ts_town          & ! Urban temperature                            [        K]
            , g_urban          & ! Urban something...                           [     ????]

            , transp_o         & ! Transpiration rate                           [  kg/m?/s]
            , gpp_o            & ! Gross primary productivity                   [?mol/m?/s]
            , leaf_resp_o      & ! Leaf respiration rate                        [?mol/m?/s]
            , root_resp_o      & ! Root respiration rate                        [?mol/m?/s]
            , het_resp_o       ! ! Heterotrophic respiration rate               [?mol/m?/s]


   !---------------------------------------------------------------------------------------!
   !    For LEAF-4, we use sun and shade variables for vegetation, then we aggregate the   !
   ! variables for output.                                                                 !
   !---------------------------------------------------------------------------------------!
   !----- Variables that must be split into sun and shade. --------------------------------!
   real(kind=4), dimension(2) :: sla_ss
   real(kind=4), dimension(2) :: lai_ss
   real(kind=4), dimension(2) :: par_l_ss
   real(kind=4), dimension(2) :: vm0_ss
   real(kind=4), dimension(2) :: rd0_ss
   real(kind=4), dimension(2) :: gpp_ss
   real(kind=4), dimension(2) :: leaf_resp_ss
   real(kind=4), dimension(2) :: transp_ss
   !---------------------------------------------------------------------------------------!


   !----- These are used for the soil bottom boundary condition. --------------------------!
   real    :: soil_water_0
   real    :: soil_tempk_0
   real    :: soil_fracliq_0

   real, allocatable, dimension(:) ::  &
              dslz                 & ! soil layer thickness at T point
            , dslzi                & ! (1. / soil layer thickness at T point)
            , dslzidt              & ! (dtll / soil layer thickness at T point)
            , slzt                 & ! soil depth at T point
            , dslzt                & ! soil layer thickness at M point
            , dslzti               & ! (1. / soil layer thickness at M point)
            , dslztidt             & ! (dtll / soil layer thickness at M point)

            , rshort_s             & ! net SW radiation absorbed by snow layers (mzs)
            , sfcwater_energy_ext  & ! Extensive sfc. water internal energy (mzs)
            , sfcwater_tempk       & ! diagnosed temp (K) of surface water
            , sfcwater_fracliq     & ! diagnosed liquid fraction surface water
            , soil_tempk           & ! diagnosed temp (K) of soil
            , soil_fracliq         & ! diagnosed liquid fraction of soil water

            , psiplusz             & ! soil water potential plus geopotential [      m]
            , hydcond              & ! hydraulic conductivity                 [    m/s]
            , th_cond_s            & ! Thermal conductivity                   [  W/m/K]
            , th_cond_p            & ! Thermal conductivity                   [  W/m/K]
            , h_flux_g             & ! Sensible heat flux at staggered layer  [   W/m2]
            , h_flux_s             & ! Sensible heat flux at staggered layer  [   W/m2]
            , w_flux_g             & ! Water flux at staggered layers         [kg/m2/s]
            , qw_flux_g            ! ! Energy flux at staggered layers        [kg/m2/s]

   logical, dimension(:), allocatable :: drysoil ! Soil is too dry
   logical, dimension(:), allocatable :: satsoil ! Soil is too wet

   !----- Variables to define the snow layers. --------------------------------------------!
   real, dimension(:,:) , allocatable   :: thick
   real, dimension(:)   , allocatable   :: thicknet


   !---------------------------------------------------------------------------------------!



   !----- Some dimensions -----------------------------------------------------------------!
   integer, parameter :: nstyp     = 17 ! # of soil types
   integer, parameter :: nvtyp     = 20 ! # of land use types
   integer, parameter :: nvtyp_teb =  1 ! # of TEB extra land use types (21 - Very urban).
   integer, parameter :: nscol     = 21 ! # of soil colour types
   !---------------------------------------------------------------------------------------!



   !----- Soil properties -----------------------------------------------------------------!
   real, dimension(nstyp)           :: slden,slcpd,slbs,slcond,sfldcap,slcons,slmsts,slpots
   real, dimension(nstyp)           :: ssand,sclay,sorgan,sporo,soilwp,soilcp,slfc,emisg
   real, dimension(nstyp)           :: slcons00,slcons0,fhydraul,xsilt,xsand,xclay
   real, dimension(nstyp)           :: psild,psiwp,psifc
   real, dimension(nstyp)           :: thcond0,thcond1,thcond2,thcond3
   real, dimension(0:nzgmax,nstyp)  :: slcons1
   real, dimension(nscol)           :: alb_vis_dry,alb_vis_wet
   real, dimension(nscol)           :: alb_nir_dry,alb_nir_wet
   !---------------------------------------------------------------------------------------!



   !----- Plant properties ----------------------------------------------------------------!
   real   , dimension(nvtyp+nvtyp_teb)        :: albv_green,albv_brown,emisv,sr_max,tai_max
   real   , dimension(nvtyp+nvtyp_teb)        :: sai,veg_clump,veg_frac,veg_ht,glai_max
   real   , dimension(nvtyp+nvtyp_teb)        :: gsw_max,dead_frac,leaf_width,stom_side
   integer, dimension(nvtyp+nvtyp_teb)        :: kroot,phenology
   real   , dimension(nzgmax,nvtyp+nvtyp_teb) :: root
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Optical properties of leaves and wood.                                            !
   !                                                                                       !
   ! CLUMPING FACTOR - factor indicating the degree of clumpiness of leaves.               !
   ! ORIENT_FACTOR   - mean leaf orientation.                                              !
   !                     0 -- leaves are randomly oriented                                 !
   !                     1 -- all leaves are perfectly horizontal                          !
   !                    -1 -- all leaves are perfectly vertical.                           !
   ! PHI1            - The phi1 term from the CLM technical manual                         !
   ! PHI2            - The phi2 term from the CLM technical manual                         !
   ! MU_BAR          - average cosine of incidence angle for hemispheric (diffuse)         !
   !                   radiation (for both short wave and long wave)                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: clumping_factor
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: orient_factor
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: phi1
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: phi2
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: mu_bar
   !---------------------------------------------------------------------------------------!
   !     Reflectance coefficients.                                                         !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: leaf_reflect_vis
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: wood_reflect_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: leaf_reflect_nir
   real(kind=8), dimension(nvtyp+nvtyp_teb)   :: wood_reflect_nir
   !---------------------------------------------------------------------------------------!
   !     Transmittance coefficients.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_trans_vis
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_trans_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_trans_nir
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_trans_nir
   !----- Emissivity of the vegetation (TIR). ---------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_emiss_tir
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_emiss_tir
   !---------------------------------------------------------------------------------------!
   !     Scattering coefficients.                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_scatter_vis
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_scatter_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_scatter_nir
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_scatter_nir
   !----- Thermal infrared. ---------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_scatter_tir
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_scatter_tir
   !---------------------------------------------------------------------------------------!
   !     Fraction of diffuse radiation that is upscattered.                                !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_backscatter_vis
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_backscatter_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_backscatter_nir
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_backscatter_nir
   !----- Backscattering of thermal infrared. ---------------------------------------------!
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: leaf_backscatter_tir
   real(kind=8), dimension(nvtyp+nvtyp_teb) :: wood_backscatter_tir
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Photosynthesis parameters.                                                       !
   !---------------------------------------------------------------------------------------!
   integer     , dimension(nvtyp+nvtyp_teb) :: pathway       ! photosyn_pathway
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: leaf_supply    ! d0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: dr_gamma       ! dark_respiration_factor
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: gsw_0          ! cuticular_cond
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: gsw_m          ! stomatal_slope
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: quantum_yield  ! quantum_efficiency
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: cn_l           ! CLM-4 CN_L
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: f_lnr          ! CLM-4 F_LNR
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: fun_nitro      ! CLM-4 f(N)
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: sla_0          ! SLA_0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: sla_m          ! SLA_m
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: vm0_qten       ! Q10 for Vm0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: vm0_dec        ! Decay factor for Vm0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: rd0_qten       ! Q10 for Rd0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: rd0_dec        ! Decay factor for Rd0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: rr0_0          ! Reference value for Rr0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: rr0_qten       ! Q10 for Rr0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: rr0_dec        ! Decay factor for Rr0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: gr_factor      ! Growth respiration factor
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: phys_low_temp  ! Decay factor for Vm0
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: phys_high_temp ! Decay factor for Vm0
   !----- Derived quantities. -------------------------------------------------------------!
   real(kind=4), dimension(nvtyp+nvtyp_teb) :: vm0_0         ! Vm0 coefficient
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!
   !    Variables needed for photosynthesis.  Suggested references:                        !
   !                                                                                       !
   ! - M09 - Medvigy, D.M., S. C. Wofsy, J. W. Munger, D. Y. Hollinger, P. R. Moorcroft,   !
   !         2009: Mechanistic scaling of ecosystem function and dynamics in space and     !
   !         time: Ecosystem Demography model version 2.  J. Geophys. Res., 114, G01002,   !
   !         doi:10.1029/2008JG000812.                                                     !
   ! - M06 - Medvigy, D.M., 2006: The State of the Regional Carbon Cycle: results from a   !
   !         constrained coupled ecosystem-atmosphere model, 2006.  Ph.D. dissertation,    !
   !         Harvard University, Cambridge, MA, 322pp.                                     !
   ! - M01 - Moorcroft, P. R., G. C. Hurtt, S. W. Pacala, 2001: A method for scaling       !
   !         vegetation dynamics: the ecosystem demography model, Ecological Monographs,   !
   !         71, 557-586.                                                                  !
   ! - F96 - Foley, J. A., I. Colin Prentice, N. Ramankutty, S. Levis, D. Pollard,         !
   !         S. Sitch, A. Haxeltime, 1996: An integrated biosphere model of land surface   !
   !         processes, terrestrial carbon balance, and vegetation dynamics. Glob.         !
   !         Biogeochem. Cycles, 10, 603-602.                                              !
   ! - L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf      !
   !         nitrogen, photosynthesis, conductance, and transpiration: scaling from leaves !
   !         to canopies. Plant, Cell, and Environ., 18, 1183-1200.                        !
   ! - F80 - Farquhar, G. D., S. von Caemmerer, J. A. Berry, 1980: A biochemical model of  !
   !         photosynthetic  CO2 assimilation in leaves of C3 species. Planta, 149, 78-90. !
   ! - C91 - Collatz, G. J., J. T. Ball, C. Grivet, J. A. Berry, 1991: Physiology and      !
   !         environmental regulation of stomatal conductance, photosynthesis and          !
   !         transpiration: a model that includes a laminar boundary layer, Agric. and     !
   !         Forest Meteorol., 54, 107-136.                                                !
   ! - C92 - Collatz, G. J., M. Ribas-Carbo, J. A. Berry, 1992: Coupled photosynthesis-    !
   !         stomatal conductance model for leaves of C4 plants.  Aust. J. Plant Physiol., !
   !         19, 519-538.                                                                  !
   ! - E78 - Ehleringer, J. R., 1978: Implications of quantum yield differences on the     !
   !         distributions of C3 and C4 grasses.  Oecologia, 31, 255-267.                  !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Bounds for internal carbon and water stomatal conductance. ----------------------!
   real(kind=4), parameter :: c34smin_lint_co2 = 0.500 * umol_2_mol
   real(kind=4), parameter :: c34smax_lint_co2 = 1200. * umol_2_mol
   real(kind=4), parameter :: c34smax_gsw      = 100.
   !----- Reference temperature for the Q10 functions. ------------------------------------!
   real(kind=4), parameter :: temp0_q10        = 25.0 + t00
   real(kind=4), parameter :: slope_q10        = 0.10
   !----- Oxygen concentration. -----------------------------------------------------------!
   real(kind=4), parameter :: o2_ref           = 0.209
   !----- Michaelis-Mentel constants for CO2, O2, and compensation point. -----------------!
   real(kind=4), parameter :: kco2_q10         = 2.1
   real(kind=4), parameter :: ko2_q10          = 1.2
   real(kind=4), parameter :: kco2_refval      = 30. * mmcod / prefsea
   real(kind=4), parameter :: ko2_refval       = 30000. * mmo2 * mmdryi / prefsea
   real(kind=4), parameter :: tau_refval_c91   =  2600.
   real(kind=4), parameter :: tau_q10          =  0.57 
   real(kind=4), parameter :: compp_refval     = o2_ref * tau_q10 / (2. * tau_refval_c91)
   real(kind=4), parameter :: compp_q10        = 1. / tau_q10
   !----- Coefficient for near saturated conditions for C4 photosynthesis. ----------------!
   real(kind=4), parameter :: klowco2          = 4000. * mmcod
   !----- Ratios between different conductances, based on L95 and C91. --------------------!
   real(kind=4), parameter :: gbh_2_gbw        = 1.075
   real(kind=4), parameter :: gbw_2_gbc        = 1.0 / 1.4
   real(kind=4), parameter :: gsw_2_gsc        = 1.0 / 1.6
   real(kind=4), parameter :: gsc_2_gsw        = 1./gsw_2_gsc
   !----- Minimum PAR to even think about opening stomata. --------------------------------!
   real(kind=4), parameter :: par_twilight_min = 0.5 * Watts_2_Ein
   !===== Double precision version of the parameters. =====================================!
   real(kind=8), parameter :: c34smin_lint_co28 = dble(c34smin_lint_co2)
   real(kind=8), parameter :: c34smax_lint_co28 = dble(c34smax_lint_co2)
   real(kind=8), parameter :: c34smax_gsw8      = dble(c34smax_gsw     )
   real(kind=8), parameter :: temp0_q108        = dble(temp0_q10       )
   real(kind=8), parameter :: slope_q108        = dble(slope_q10       )
   real(kind=8), parameter :: o2_ref8           = dble(o2_ref          )
   real(kind=8), parameter :: kco2_q108         = dble(kco2_q10        )
   real(kind=8), parameter :: ko2_q108          = dble(ko2_q10         )
   real(kind=8), parameter :: kco2_refval8      = dble(kco2_refval     )
   real(kind=8), parameter :: ko2_refval8       = dble(ko2_refval      )
   real(kind=8), parameter :: tau_refval_c918   = dble(tau_refval_c91  )
   real(kind=8), parameter :: tau_q108          = dble(tau_q10         )
   real(kind=8), parameter :: compp_refval8     = dble(compp_refval    )
   real(kind=8), parameter :: compp_q108        = dble(compp_q10       )
   real(kind=8), parameter :: klowco28          = dble(klowco2         )
   real(kind=8), parameter :: gbh_2_gbw8        = dble(gbh_2_gbw       )
   real(kind=8), parameter :: gbw_2_gbc8        = dble(gbw_2_gbc       )
   real(kind=8), parameter :: gsw_2_gsc8        = dble(gsw_2_gsc       )
   real(kind=8), parameter :: gsc_2_gsw8        = dble(gsc_2_gsw       )
   real(kind=8), parameter :: par_twilight_min8 = dble(par_twilight_min)
   !---------------------------------------------------------------------------------------!

   !----- Other variables -----------------------------------------------------------------!
   real                             :: cmin,corg,cwat,cair,cka,ckw
   !---------------------------------------------------------------------------------------!

   !----- Roughness -----------------------------------------------------------------------!
   real(kind=4), parameter :: z0fac_water     = .016 / grav ! Coefficient before ustar?
   real(kind=4), parameter :: min_waterrough  = .0001       ! Min. water roughness height
   real(kind=4), parameter :: waterrough      = .0001       ! Water roughness height
   real(kind=4), parameter :: snowrough       = .0024       ! Snow roughness height
   !----- Double precision version of some variables above. -------------------------------!
   real(kind=8), parameter :: z0fac_water8    = dble(z0fac_water   )
   real(kind=8), parameter :: min_waterrough8 = dble(min_waterrough)
   !---------------------------------------------------------------------------------------!

   !----- Maximum transpiration allowed. --------------------------------------------------!
   real, parameter :: transp_max     = 400. / alvl3
   !---------------------------------------------------------------------------------------!

   !----- Is super-saturation fine? -------------------------------------------------------!
   logical, parameter :: supersat_ok    = .false.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Original LEAF-3 parameters for heat capacity of vegetation.                       !
   !---------------------------------------------------------------------------------------!
   real, parameter :: hcapveg_ref  = 3.e3      ! [J/m2/K]
   real, parameter :: hcapveg_hmin = 1.5
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     New parameters for heat capacity, based on ED-2.  AGB is a simple estimate, based !
   ! on lidar allometric equations for tropical environments.                              !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Gu, L., T. Meyers, S. G. Pallardy, 2007: Influences of biomass heat and biochemical   !
   !      energy storages on the land surface fluxes and radiative temperature.            !
   !      J. Geophys. Res., v. 112, doi: 10.1029/2006JD007425.                             !
   !                                                                                       !
   ! Longo, M., 2013: Amazon forest response to changes in rainfall regime: results from   !
   !    an individual-based dynamic vegetation model.   Ph.D. dissertation, Harvard        !
   !    University, Cambridge, MA, Dec 2013.                                               !
   !                                                                                       !
   ! Asner, G. P. J. Mascaro, 2014: Mapping tropical forest carbon: calibrating plot       !
   !    estimates to a simple LiDar metric.  Remote Sens. Environ., v. 140, 614-624,       !
   !    doi: 10.1016/j.rse.2013.09.023.                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: agb_am14_a      = 0.685
   real, parameter :: agb_am14_b      = 0.952
   real, parameter :: gu_tref         = t00 + 15.0
   real, parameter :: gu_c_leaf_dry   = 3218.
   real, parameter :: gu_wat2dry_leaf = 1.85
   real, parameter :: gu_wat2dry_wood = 0.7
   real, parameter :: gu_c_wood_dry   = 103.1 + 3.867 * gu_tref
   real, parameter :: gu_delta_c      = 1.e5 * gu_wat2dry_wood                             &
                                      * ( - 0.06191                                        &
                                          + 2.36e-4 * gu_tref                              &
                                          - 0.0133  * gu_wat2dry_wood )
   real, parameter :: gu_spheat_leaf  = ( gu_c_leaf_dry + gu_wat2dry_leaf * cliq )         &
                                      / ( 1. + gu_wat2dry_leaf )
   real, parameter :: gu_spheat_wood  = ( gu_c_wood_dry + gu_wat2dry_wood * cliq )         &
                                      / (1. + gu_wat2dry_wood ) + gu_delta_c
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Minimum canopy air space depth.                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: can_depth_min = 5.0 ! [J/m2/K]
   !---------------------------------------------------------------------------------------!

   !----- Some constants to ensure the model good behaviour -------------------------------!
   real, parameter :: min_sfcwater_mass     = 1.e-6
   real, parameter :: min_sfcwater_depth    = 1.e-9
   real, parameter :: water_stab_thresh     = 3.0  
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Constants for surface layer models.                                              !
   !---------------------------------------------------------------------------------------!
   real, parameter :: vh2vr    = 0.13          ! Vegetation roughness:vegetation hgt ratio
   real, parameter :: vh2dh    = 0.63          ! Displacement height:vegetation hgt ratio
   !----- Louis (1979) model. -------------------------------------------------------------!
   real, parameter :: bl79     = 5.0 ! b prime parameter                                       
   real, parameter :: csm      = 7.5 ! C* for momentum (eqn. 20, not co2 char. scale)          
   real, parameter :: csh      = 5.0 ! C* for heat (eqn.20, not co2 char. scale)               
   real, parameter :: dl79     = 5.0 ! ???                                                     
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   real, parameter :: beta_s   = 5.0    ! Beta used by Businger et al. (1971)
   !----- Beljaars and Holtslag (1991) model. ---------------------------------------------!
   real, parameter :: abh91       = -1.00         ! -a from equation  (28) and (32)
   real, parameter :: bbh91       = -twothirds    ! -b from equation  (28) and (32)
   real, parameter :: cbh91       =  5.0          !  c from equations (28) and (32)
   real, parameter :: dbh91       =  0.35         !  d from equations (28) and (32)
   real, parameter :: ebh91       = -twothirds    ! the 2/3 factor in equation (32)
   real, parameter :: fbh91       =  1.50         ! exponent in equation (32)
   real, parameter :: cod         = cbh91/dbh91   ! c/d
   real, parameter :: bcod        = bbh91 * cod   ! b*c/d
   real, parameter :: fm1         = fbh91 - 1.0   ! f-1
   real, parameter :: ate         = abh91 * ebh91 ! a * e
   real, parameter :: atetf       = ate   * fbh91 ! a * e * f
   real, parameter :: z0moz0h     = 1.0           ! z0(M)/z0(h)
   real, parameter :: z0hoz0m     = 1. / z0moz0h  ! z0(M)/z0(h)
   !    Modified CLM (2004) model.   These will be initialised later. -------------------!
   real   :: beta_vs     ! Beta for the very stable case (CLM eq. 5.30)
   real   :: chim        ! CLM coefficient for very unstable case (momentum)
   real   :: chih        ! CLM coefficient for very unstable case (heat)
   real   :: zetac_um    ! critical zeta below which it becomes very unstable (momentum)
   real   :: zetac_uh    ! critical zeta below which it becomes very unstable (heat)
   real   :: zetac_sm    ! critical zeta above which it becomes very stable   (momentum)
   real   :: zetac_sh    ! critical zeta above which it becomes very stable   (heat)
   real   :: zetac_umi   ! 1. / zetac_umi
   real   :: zetac_uhi   ! 1. / zetac_uhi
   real   :: zetac_smi   ! 1. / zetac_smi
   real   :: zetac_shi   ! 1. / zetac_shi
   real   :: zetac_umi16 ! 1/(-zetac_umi)^(1/6)
   real   :: zetac_uhi13 ! 1/(-zetac_umi)^(1/6)
   real   :: psimc_um    ! psim evaluation at zetac_um
   real   :: psihc_uh    ! psih evaluation at zetac_uh
   !---------------------------------------------------------------------------------------!



   !----- Parameters that used to be in LEAF-3 (leaftw) -----------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Parameters used to compute the aerodynamic resistance.                            !
   !---------------------------------------------------------------------------------------!
   real, parameter                   :: exar       =  2.50
   real, parameter                   :: covr       =  2.16
   real, parameter                   :: rasveg_min =  1.e-6
   real, parameter                   :: taumin     =  1.e-6
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar))                                 !
   !     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is                        !
   !     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.                   !
   !     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the                    !
   !     total expression inside the radical in Eq. 3.37.                                  !
   !---------------------------------------------------------------------------------------!
   real, parameter                   :: c1=116.6
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   real, parameter                   :: brad =   196.0,  srad =   0.047
   real, parameter                   :: btlo =   281.5,  stlo =   0.26
   real, parameter                   :: bthi =   310.1,  sthi =  -0.124
   real, parameter                   :: bvpd =  4850.0,  svpd =  -0.0051
   real, parameter                   :: bsmp = -1.07e6,  ssmp =   7.42e-6
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Parameters for the aerodynamic resistance between the leaf and the canopy air    !
   ! space.  These are the A, B, n, and m parameters that define the Nusselt number for    !
   ! forced and free convection, at equations 10.7 and 10.9.  The parameters are found at  !
   ! the appendix A.5(a) and A.5(b).                                                       !
   !                                                                                       !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !---------------------------------------------------------------------------------------!
   real(kind=4), parameter :: aflat_lami = 0.600    ! A (forced convection), laminar   flow
   real(kind=4), parameter :: nflat_lami = 0.500    ! n (forced convection), laminar   flow
   real(kind=4), parameter :: aflat_turb = 0.032    ! A (forced convection), turbulent flow
   real(kind=4), parameter :: nflat_turb = 0.800    ! n (forced convection), turbulent flow
   real(kind=4), parameter :: bflat_lami = 0.500    ! B (free   convection), laminar   flow
   real(kind=4), parameter :: mflat_lami = 0.250    ! m (free   convection), laminar   flow
   real(kind=4), parameter :: bflat_turb = 0.190    ! B (free   convection), turbulent flow
   real(kind=4), parameter :: mflat_turb = onethird ! m (free   convection), turbulent flow
   !---------------------------------------------------------------------------------------!



   !----- Constants from  equation E27 (Medvigy 2007) -------------------------------------!
   real(kind=4), dimension(6), parameter :: ss = (/ 1.093e-3, 2.800e-2, 3.000e-2           &
                                                  , 3.030e-4,-1.770e-7, 2.250e-9 /)  
   !---------------------------------------------------------------------------------------!



   !-----  Exponent in the frozen soil hydraulic conductivity correction. -----------------!
   real(kind=4)              , parameter :: freezeCoef = 7.0    
   !---------------------------------------------------------------------------------------!



   !-----  Bounds for solving the vegetation. ---------------------------------------------!
   real(kind=4)              , parameter :: tai_min     = 0.1
   real(kind=8)              , parameter :: tai_min8    = dble(tai_min)
   real(kind=4)              , parameter :: snowfac_max = 0.9
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   Soil conductance terms, from:                                                       !
   !                                                                                       !
   ! Passerat de Silans, A., 1986: Transferts de masse et de chaleur dans un sol stratifi? !
   !     soumis ? une excitation amtosph?rique naturelle. Comparaison: Mod?les-exp?rience. !
   !     Thesis, Institut National Polytechnique de Grenoble. (P86)                        !
   !                                                                                       !
   ! retrieved from:                                                                       !
   ! Mahfouf, J. F., J. Noilhan, 1991: Comparative study of various formulations of        !
   !     evaporation from bare soil using in situ data. J. Appl. Meteorol., 30, 1354-1365. !
   !     (MN91)                                                                            !
   !                                                                                       !
   !     Please notice that the values are inverted because we compute conductance, not    !
   ! resistance.                                                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4)              , parameter :: ggsoil0 = 1. / 38113.
   real(kind=4)              , parameter :: kksoil  = 13.515
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for fraction covered with snow, which is based on:                     !
   !                                                                                       !
   ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover     !
   !    fraction and its evaluation over large North American river basins,                !
   !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                            !
   !                                                                                       !
   !    These are the parameters in equation 4.  Fresh snow density is defined at          !
   ! consts_coms.f90                                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4)              , parameter :: ny07_eq04_a  = 2.5
   real(kind=4)              , parameter :: ny07_eq04_m  = 1.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Ocean optical parameters.                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4)              , parameter :: alb_oc_inter = -0.0139
   real(kind=4)              , parameter :: alb_oc_slope =  0.0467
   real(kind=4)              , parameter :: alb_oc_min   =  0.0300
   real(kind=4)              , parameter :: alb_oc_max   =  0.9999
   real(kind=4)              , parameter :: emiss_oc     =  0.97
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Optical properties for snow.  Values are a first guess, and a more thorough snow  !
   ! model that takes snow age and snow melt into account (like in CLM-4 or ECHAM-5) are   !
   ! very welcome.                                                                         !
   !                                                                                       !
   !  References for current snow values:                                                  !
   !  Roesch, A., et al., 2002: Comparison of spectral surface albedos and their           !
   !      impact on the general circulation model simulated surface climate.  J.           !
   !      Geophys. Res.-Atmosph., 107(D14), 4221, 10.1029/2001JD000809.                    !
   !      Average between minimum and maximum snow albedo on land, af = 0. and af=1.       !
   !                                                                                       !
   !  Oleson, K.W., et al., 2010: Technical description of version 4.0 of the              !
   !      Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR.                 !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4)              , parameter :: alb_snow_par = 0.518
   real(kind=4)              , parameter :: alb_snow_nir = 0.430
   real(kind=4)              , parameter :: alb_snow     = 0.500
   real(kind=4)              , parameter :: emiss_snow   = 0.970
   !----- Water pounding parameter (LEAF-3 only). -----------------------------------------!
   real(kind=4)              , parameter :: alb_damp     = 0.14
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Soil heterotrophic respiration parameters.                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=4)              , parameter :: rhmax_0         = 1.20
   real(kind=4)              , parameter :: rhmax_m         = 2.58
   real(kind=4)              , parameter :: depth_hetresp   = -0.20
   integer                               :: k_hetresp
   real(kind=4)              , parameter :: decay_low_rh    = 0.24
   real(kind=4)              , parameter :: decay_high_rh   = 0.60
   real(kind=4)              , parameter :: low_temp_rh     = 18.0 + t00
   real(kind=4)              , parameter :: high_temp_rh    = 45.0 + t00
   real(kind=4)              , parameter :: decay_dry_rh    = 12.0 ! 18.0
   real(kind=4)              , parameter :: decay_wet_rh    = 36.0 ! 36.0
   real(kind=4)              , parameter :: dry_smoist_rh   = 0.48 ! 0.36
   real(kind=4)              , parameter :: wet_smoist_rh   = 0.98 ! 0.96
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine allocates some specific scratch arrays used by LEAF-3.              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_leafcol(nzg,nzs)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: nzg
      integer, intent(in) :: nzs
      !----- Local variables. -------------------------------------------------------------!
      integer :: k
      integer :: kzs
      real    :: thick_1st
      real    :: stretch
      !------------------------------------------------------------------------------------!

      !----- Allocate leaf column arrays. -------------------------------------------------!
      allocate (dslz                (0:nzg)    )
      allocate (dslzi               (0:nzg)    )
      allocate (dslzidt             (0:nzg)    )
      allocate (slzt                (0:nzg)    )
      allocate (dslzt               (nzg)      )
      allocate (dslzti              (nzg)      )
      allocate (dslztidt            (nzg)      )
      allocate (soil_tempk          (nzg)      )
      allocate (soil_fracliq        (nzg)      )

      allocate (rshort_s            (nzs)      )
      allocate (sfcwater_energy_ext (nzs)      )
      allocate (sfcwater_tempk      (nzs)      )
      allocate (sfcwater_fracliq    (nzs)      )

      allocate (psiplusz            (0:nzg)    )
      allocate (hydcond             (0:nzg)    )
      allocate (th_cond_s           (0:nzg)    )
      allocate (th_cond_p           (0:nzs)    )
      allocate (drysoil             (0:nzg)    )
      allocate (satsoil             (0:nzg)    )
      allocate (h_flux_g            (nzg+1)    )
      allocate (h_flux_s            (nzs+1)    )
      allocate (w_flux_g            (nzg+1)    )
      allocate (qw_flux_g           (nzg+1)    )

      allocate (thick             (nzs+1,nzs+1))
      allocate (thicknet                (nzs+1))

      !----- Initialise snow thickness scaling array. -------------------------------------!
      stretch = 2.0
      do kzs = 1,nzs
         thick_1st = 1.0
         thicknet(kzs) = 0.0
         do k = 1,(kzs+1)/2
            thick(k,kzs) = thick_1st
            thick(kzs+1-k,kzs) = thick_1st
            thicknet(kzs) = thicknet(kzs) + 2. * thick_1st
            thick_1st = thick_1st * stretch
         end do
         if ((kzs+1)/2 /= kzs/2) thicknet(kzs) = thicknet(kzs) - thick_1st/stretch
         do k = 1,kzs
            thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
         end do
      end do

               
      return
   end subroutine alloc_leafcol   
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine flushes to zero all scratch variables that should be recycled by !
   ! LEAF within a call.  The flag is used to determine whether which variables should be  !
   ! flushed:                                                                              !
   ! 0 - Local, only variables that should be recalculated every iteration are flushed.    !
   ! 1 - Partial, variables that are updated once for each patch are flushed, plus all     !
   !     variables that were flushed in option zero.                                       !
   ! 2 - Partial, variables that are updated once for each grid point are flushed, plus    !
   !     all variables that were flushed in options zero and one.                          !
   ! 3 - Total, all variables are flushed.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine flush_leaf_coms(idel)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*), intent(in) :: idel
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following variables are reset only once per leaf3 call.  Internal          !
      ! constants or variables that depend only on the grid called should be flushed here. !
      !------------------------------------------------------------------------------------!
      select case (trim(idel))
      case ('INITIAL')
         niter_leaf3   = 0
         ndtvegi       = 0.
         dtl3          = 0.
         dtl3_factor   = 0.
         dtvg          = 0.
         timefac_sst   = 0.
         timefac_ndvi  = 0.

         emis_town     = 0.
         alb_town      = 0.
         ts_town       = 0.
         g_urban       = 0.

         dslz    (:)   = 0
         dslzi   (:)   = 0.
         dslzidt (:)   = 0.
         slzt    (:)   = 0.
         dslzt   (:)   = 0.
         dslzti  (:)   = 0.
         dslztidt(:)   = 0.

      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following variables are reset only once per grid point call.  Variables    !
      ! like the atmospheric forcing should be placed here.                                !
      !------------------------------------------------------------------------------------!
      select case (trim(idel))
      case ('INITIAL','GRID_POINT')
         atm_up        = 0
         atm_vp        = 0.
         atm_thil      = 0.
         atm_theta     = 0
         atm_temp      = 0.
         atm_temp_zcan = 0.
         atm_enthalpy  = 0.
         atm_rvap      = 0.
         atm_rtot      = 0
         atm_shv       = 0.
         atm_co2       = 0.
         atm_theiv     = 0
         atm_vpdef     = 0
         atm_rhos      = 0.
         geoht         = 0.
         atm_exner     = 0
         atm_prss      = 0.
         atm_vels      = 0.
         pcpgl         = 0
         qpcpgl        = 0.
         dpcpgl        = 0.
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following variables are reset every time a new patch is called.            !
      ! Internal prognostic variables and variables that should remain constant during the !
      ! time steps should be placed here.                                                  !
      !------------------------------------------------------------------------------------!
      select case (trim(idel))
      case ('INITIAL','GRID_POINT','PATCH')
         resolvable              = .false.

         can_enthalpy            = 0.
         can_exner               = 0
         can_rhos                = 0.
         can_depth               = 0.
         flag_sfcwater           = 0

         sfcwater_energy_ext (:) = 0.

      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following variables are reset every time.  These are mostly diagnostic     !
      ! variables that should be computed from scratch every time the time step is called. !
      !------------------------------------------------------------------------------------!
      dtl3owcc               = 0
      dtl3ohcc               = 0.
      dtl3occc               = 0.
      snowfac                = 0.

      can_cp                 = 0.
      can_temp               = 0
      can_rsat               = 0.
      can_shv                = 0.
      can_rhv                = 0.
      veg_temp               = 0
      veg_fliq               = 0.
      veg_wind               = 0.
      estar                  = 0.
      qstar                  = 0.

      soil_water_0           = 0.
      soil_tempk_0           = 0.
      soil_fracliq_0         = 0.

      gbh                    = 0.
      gbw                    = 0
      ggnet                  = 0.
      ggbare                 = 0.
      ggsoil                 = 0.
      ggveg                  = 0.
      rho_ustar              = 0.

      hflxsc                 = 0
      wflxsc                 = 0.
      qwflxsc                = 0.
      hflxgc                 = 0
      wflxgc                 = 0.
      qwflxgc                = 0.
      cflxgc                 = 0.
      eflxac                 = 0
      hflxac                 = 0.
      wflxac                 = 0.
      qwflxac                = 0.
      cflxac                 = 0.
      hflxvc                 = 0.
      wflxvc                 = 0
      qwflxvc                = 0.
      cflxvc                 = 0.
      transp                 = 0
      qtransp                = 0.
      gpp                    = 0.
      wshed                  = 0
      qwshed                 = 0.
      dwshed                 = 0.

      hflxvc_tot             = 0.
      wflxvc_tot             = 0
      qwflxvc_tot            = 0.
      cflxvc_tot             = 0.
      cflxgc_tot             = 0.
      transp_tot             = 0.
      qtransp_tot            = 0.

      wshed_tot              = 0.
      qwshed_tot             = 0.
      dwshed_tot             = 0
      throughfall_tot        = 0.
      qthroughfall_tot       = 0.
      dthroughfall_tot       = 0.
      intercepted_tot        = 0
      qintercepted_tot       = 0.
      dintercepted_tot       = 0.
      dewgnd_tot             = 0.
      qdewgnd_tot            = 0.
      ddewgnd_tot            = 0.

      virtual_energy         = 0.
      virtual_water          = 0.
      virtual_depth          = 0.

      rshort_s         (:)   = 0
      sfcwater_tempk   (:)   = 0.
      sfcwater_fracliq (:)   = 0.
      soil_tempk       (:)   = 0
      soil_fracliq     (:)   = 0.

      psiplusz         (:)   = 0.
      hydcond          (:)   = 0.
      th_cond_s        (:)   = 0.
      th_cond_p        (:)   = 0.
      drysoil          (:)   = .false.
      satsoil          (:)   = .false.
      h_flux_g         (:)   = 0.
      h_flux_s         (:)   = 0.
      w_flux_g         (:)   = 0.
      qw_flux_g        (:)   = 0.

      rshort_g               = 0.
      rshort_v               = 0.
      rshort_a               = 0.
      rlong_v                = 0.
      rlong_g                = 0.
      rlong_s                = 0.

      sla_ss           (:)   = 0.
      lai_ss           (:)   = 0.
      par_l_ss         (:)   = 0.
      vm0_ss           (:)   = 0.
      rd0_ss           (:)   = 0.
      gpp_ss           (:)   = 0.
      leaf_resp_ss     (:)   = 0.
      transp_ss        (:)   = 0.

      gsw                    = 0.
      leaf_resp_o            = 0.
      root_resp_o            = 0.
      het_resp_o             = 0.
      gpp_o                  = 0.
      transp_o               = 0.
      !------------------------------------------------------------------------------------!

      select case (trim(idel))
      case ('PATCH')
         continue
      end select

      return
   end subroutine flush_leaf_coms
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine initialises several parameters for the surface layer model.      !
   !---------------------------------------------------------------------------------------!
   subroutine sfclyr_init_params()
      use rconstants, only : onesixth
      implicit none
      !----- External functions. ----------------------------------------------------------!
      real, external :: cbrt ! Cubic root
      !------------------------------------------------------------------------------------!
     
      !----- Similar to CLM (2004), but with different phi_m for very unstable case. ------!
      zetac_um    = -1.5
      zetac_uh    = -0.5
      zetac_sm    =  1.0
      zetac_sh    =  zetac_sm
      !----- Define chim and chih so the functions are continuous. ------------------------!
      chim        = (-zetac_um) ** onesixth / sqrt(sqrt(1.0 - gamm * zetac_um))
      chih        = cbrt(-zetac_uh) / sqrt(1.0 - gamh * zetac_uh)
      beta_vs     = 1.0 - (1.0 - beta_s) * zetac_sm
      !----- Define derived values to speed up the code a little. -------------------------!
      zetac_umi   = 1.0 / zetac_um
      zetac_uhi   = 1.0 / zetac_uh
      zetac_smi   = 1.0 / zetac_sm
      zetac_shi   = 1.0 / zetac_sh
      zetac_umi16 = 1.0 / (- zetac_um) ** onesixth
      zetac_uhi13 = 1.0 / cbrt(-zetac_uh)

      !------------------------------------------------------------------------------------!
      !     Initialise these values with dummies, it will be updated after we define the   !
      ! functions.                                                                         !
      !------------------------------------------------------------------------------------!
      psimc_um  = 0.
      psimc_um  = psim(zetac_um,.false.)
      psihc_uh  = 0.
      psihc_uh  = psih(zetac_uh,.false.)
      !------------------------------------------------------------------------------------!

      return
   end subroutine sfclyr_init_params
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for momentum.            !
   !---------------------------------------------------------------------------------------!
   real function psim(zeta,stable)
      use rconstants, only : halfpi   & ! intent(in)
                           , onesixth ! ! intent(in)
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: xx
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psim = - beta_s * zeta 
         case (3) !----- Beljaars and Holtslag (1991). ------------------------------------!
            psim = abh91 * zeta                                                            &
                 + bbh91 * (zeta - cod) * exp(max(-38.,-dbh91 * zeta))                     &
                 + bcod
         case (4) !----- CLM (2004) (including neglected terms). --------------------------!
            if (zeta > zetac_sm) then
               !----- Very stable case. ---------------------------------------------------!
               psim = (1.0 - beta_vs) * log(zeta * zetac_smi)                              &
                    + (1.0 - beta_s ) * zetac_sm - zeta
            else
               !----- Normal stable case. -------------------------------------------------!
               psim = - beta_s * zeta
            end if
         end select
      else
         select case (istar)
         case (2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). -----!
            xx   = sqrt(sqrt(1.0 - gamm * zeta))
            psim = log(0.125 * (1.0+xx) * (1.0+xx) * (1.0 + xx*xx)) - 2.0*atan(xx) + halfpi
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               psim = log(zeta * zetac_umi)                                                &
                    + 6.0 * chim * ((- zeta) ** (-onesixth) - zetac_umi16)                 &
                    + psimc_um
            else
               !----- Normal unstable case. -----------------------------------------------!
               xx   = sqrt(sqrt(1.0 - gamm * zeta))
               psim = log(0.125 * (1.0+xx) * (1.0+xx) * (1.0 + xx*xx))                     &
                    - 2.0*atan(xx) + halfpi
            end if
         end select
      end if
      return
   end function psim
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for heat (and vapour,    !
   ! and carbon dioxide too.)                                                              !
   !---------------------------------------------------------------------------------------!
   real function psih(zeta,stable)
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: yy
      !----- External functions. ----------------------------------------------------------!
      real   , external   :: cbrt
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psih = - beta_s * zeta 
         case (3) !----- Beljaars and Holtslag (1991). ------------------------------------!
            psih = 1.0 - (1.0 + ate * zeta)**fbh91                                         &
                 + bbh91 * (zeta - cod) * exp(max(-38.,-dbh91 * zeta)) + bcod
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sh) then
               !----- Very stable case. ---------------------------------------------------!
               psih = (1.0 - beta_vs) * log(zeta * zetac_shi)                              &
                    + (1.0 - beta_s ) * zetac_sh - zeta
            else
               !----- Normal stable case. -------------------------------------------------!
               psih = - beta_s * zeta 
            end if
         end select
      else
         select case (istar)
         case (2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). -----!
            yy   = sqrt(1.0 - gamh * zeta)
            psih = log(0.25 * (1.0+yy) * (1.0+yy))
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               psih = log(zeta * zetac_uhi)                                                &
                    + 3.0 * chih * (1./cbrt(-zeta) - zetac_uhi13)                          &
                    + psihc_uh
            else
               !----- Normal unstable case. -----------------------------------------------!
               yy   = sqrt(1.0 - gamh * zeta)
               psih = log(0.25 * (1.0+yy) * (1.0+yy))
            end if
         end select
      end if
      return
   end function psih
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! momentum with respect to zeta.                                                        !
   !---------------------------------------------------------------------------------------!
   real function dpsimdzeta(zeta,stable)
      use mem_leaf  , only : istar
      use rconstants, only : onesixth
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: xx
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsimdzeta = - beta_s 
         case (3) !----- Beljaars and Holtslag (1991). ------------------------------------!
            dpsimdzeta = abh91 + bbh91 * (1.0 - dbh91 * zeta + cbh91)                      &
                               * exp(max(-38.,-dbh91 * zeta))
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sm) then
               !----- Very stable case. ---------------------------------------------------!
               dpsimdzeta = (1.0 - beta_vs) / zeta - 1.0
            else
               !----- Normal stable case. -------------------------------------------------!
               dpsimdzeta = - beta_s 
            end if
         end select
      else
         select case (istar)
         case (2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). -----!
            xx         = sqrt(sqrt(1.0 - gamm * zeta))
            dpsimdzeta = - gamm / (xx * (1.0+xx) * (1.0 + xx*xx)) 
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               dpsimdzeta = (1.0 - chim * (-zeta)**onesixth) / zeta
            else
               !----- Normal unstable case. -----------------------------------------------!
               xx         = sqrt(sqrt(1.0 - gamm * zeta))
               dpsimdzeta = - gamm / (xx * (1.0+xx) * (1.0 + xx*xx))
            end if
         end select
      end if

      return
   end function dpsimdzeta
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! heat/moisture/CO2 with respect to zeta.                                               !
   !---------------------------------------------------------------------------------------!
   real function dpsihdzeta(zeta,stable)
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: yy
      !----- External functions. ----------------------------------------------------------!
      real   , external   :: cbrt
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsihdzeta = - beta_s
         case (3) !----- Beljaars and Holtslag (1991). ------------------------------------!
            dpsihdzeta = - atetf * (1.0 + ate * zeta)**fm1                                 &
                         + bbh91 * (1.0 - dbh91 * zeta + cbh91)                            &
                         * exp(max(-38.,-dbh91 * zeta))
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sh) then
               !----- Very stable case. ---------------------------------------------------!
               dpsihdzeta = (1.0 - beta_vs) / zeta - 1.0
            else
               !----- Normal stable case. -------------------------------------------------!
               dpsihdzeta = - beta_s
            end if
         end select
      else
         select case (istar)
         case (2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). -----!
            yy   = sqrt(1.0 - gamh * zeta)
            dpsihdzeta = -gamh / (yy * (1.0 + yy))
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               dpsihdzeta = (1.0 + chih / cbrt(zeta)) / zeta
            else
               !----- Normal unstable case. -----------------------------------------------!
               yy   = sqrt(1.0 - gamh * zeta)
               dpsihdzeta = -gamh / (yy * (1.0 + yy))
            end if
         end select
      end if

      return
   end function dpsihdzeta
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the value of zeta for a given Richardson number, reference    !
   ! height and the roughness scale.  This is solved by using the definition of Obukhov    !
   ! length scale as stated in Louis (1979) equation (10), modified to define z/L rather   !
   ! than L.  The solution is found  iteratively since it's not a simple function to       !
   ! invert.  It tries to use Newton's method, which should take care of most cases.  In   !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real function zoobukhov(rib,zstar,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
      use therm_lib, only : toler  & ! intent(in)
                          , maxfpo & ! intent(in)
                          , maxit  ! ! intent(in)
      use mem_leaf , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: rib       ! Bulk Richardson number                    [   ---]
      real   , intent(in) :: zstar     ! Reference height                          [     m]
      real   , intent(in) :: rough     ! Roughness length scale                    [     m]
      real   , intent(in) :: zoz0m     ! zstar/roughness(momentum)                 [   ---]
      real   , intent(in) :: lnzoz0m   ! ln[zstar/roughness(momentum)]             [   ---]
      real   , intent(in) :: zoz0h     ! zstar/roughness(heat)                     [   ---]
      real   , intent(in) :: lnzoz0h   ! ln[zstar/roughness(heat)]                 [   ---]
      logical, intent(in) :: stable    ! Flag... This surface layer is stable      [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: ribuse    ! Richardson number to use                  [   ---]
      real                :: fm        ! lnzoz0m - psim(zeta) + psim(zeta0m)       [   ---]
      real                :: fh        ! lnzoz0h - psih(zeta) + psih(zeta0h)       [   ---]
      real                :: dfmdzeta  ! d(fm)/d(zeta)                             [   ---]
      real                :: dfhdzeta  ! d(fh)/d(zeta)                             [   ---]
      real                :: z0moz     ! Roughness(momentum) / Reference height    [   ---]
      real                :: zeta0m    ! Roughness(momentum) / Obukhov length      [   ---]
      real                :: z0hoz     ! Roughness(heat) / Reference height        [   ---]
      real                :: zeta0h    ! Roughness(heat) / Obukhov length          [   ---]
      real                :: zetaa     ! Smallest guess (or previous guess)        [   ---]
      real                :: zetaz     ! Largest guess (or new guess in Newton's)  [   ---]
      real                :: deriv     ! Function Derivative                       [   ---]
      real                :: fun       ! Function for which we seek a root.        [   ---]
      real                :: funa      ! Smallest guess function.                  [   ---]
      real                :: funz      ! Largest guess function.                   [   ---]
      real                :: delta     ! Aux. var --- 2nd guess for bisection      [   ---]
      real                :: coeff     ! RiB * zstar / (Pr * (zstar - z0))         [   ---]
      real                :: zetamin   ! Minimum zeta for stable case.             [   ---]
      real                :: zetamax   ! Maximum zeta for unstable case.           [   ---]
      real                :: zetasmall ! Zeta dangerously close to zero            [   ---]
      integer             :: itb       ! Iteration counters                        [   ---]
      integer             :: itn       ! Iteration counters                        [   ---]
      integer             :: itp       ! Iteration counters                        [   ---]
      logical             :: converged ! Flag... The method converged!             [   T|F]
      logical             :: zside     ! Flag... I'm on the z-side.                [   T|F]
      !------------------------------------------------------------------------------------!



      !----- Define some values that won't change during the iterative method. ------------!
      z0moz = 1. / zoz0m
      z0hoz = 1. / zoz0h
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     First thing, check whether this is a stable case and we are running methods 2  !
      ! or 4.  In these methods, there is a singularity that must be avoided.              !
      !------------------------------------------------------------------------------------!
      select case (istar)
      case (2,4)
         ribuse = min(rib, (1.0 - toler) * tprandtl / beta_s)

         !---------------------------------------------------------------------------------!
         !    Stable case, use Oncley and Dudhia, we can solve it analytically.            !
         !---------------------------------------------------------------------------------!
         if (stable .and. istar == 2) then
            zoobukhov = ribuse * zstar * min(lnzoz0m,lnzoz0h)                              &
                      / ( (zstar-rough) * (tprandtl - beta_s * ribuse) )
            return
         end if
         !---------------------------------------------------------------------------------!
      case default
         ribuse = rib
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Define the coefficient Ri * zstar / [Pr * (zstar-z0)]                          !
      !------------------------------------------------------------------------------------!
      coeff = ribuse * zstar / (tprandtl * (zstar - rough))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the bulk Richardson number is zero or almost zero, then we rather just      !
      ! assign z/L to be the one similar to Oncley and Dudhia (1995).  This saves time and !
      ! also avoids the risk of having zeta with the opposite sign.                        !
      !------------------------------------------------------------------------------------!
      zetasmall = coeff * min(lnzoz0m,lnzoz0h)
      if (ribuse <= 0. .and. zetasmall > - z0moz0h * toler) then
         zoobukhov = zetasmall
         return
      elseif (ribuse > 0. .and. zetasmall < z0moz0h * toler) then
         zoobukhov = zetasmall / (1.0 - beta_s * ribuse / tprandtl)
         return
      else
         zetamin    =  toler
         zetamax    = -toler
      end if
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(5(a,1x,f11.4,1x),a,l1)')                                         &
      !   'Input values: Rib =',rib,'zstar=',zstar,'rough=',rough,'zoz0=',zoz0              &
      !           ,'lnzoz0=',lnzoz0,'stable=',stable
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!



      !------------------------------------------------------------------------------------!
      !     First guess, using Oncley and Dudhia (1995) approximation for unstable case.   !
      ! We won't use the stable case to avoid FPE or zeta with opposite sign when          !
      ! Ri is too positive.                                                                !
      !------------------------------------------------------------------------------------!
      zetaa = zetasmall
      !------------------------------------------------------------------------------------!



      !----- Find the function and its derivative. ----------------------------------------!
      zeta0m   = zetaa * z0moz
      zeta0h   = zetaa * z0hoz
      fm       = lnzoz0m - psim(zetaa,stable) + psim(zeta0m,stable)
      fh       = lnzoz0h - psih(zetaa,stable) + psih(zeta0h,stable)
      dfmdzeta = z0moz * dpsimdzeta(zeta0m,stable) - dpsimdzeta(zetaa,stable)
      dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaa,stable)
      funa     = coeff * fm * fm / fh - zetaa
      deriv    = coeff * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / (fh * fh) - 1.
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration. ---------------------------!
      zetaz = zetaa
      fun   = funa
      !------------------------------------------------------------------------------------!


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
      !   '1STGSS: itn=',0,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm          &
      !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Enter Newton's method loop. --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1, maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Newton's method converges fast when it's on the right track, but there are  !
         ! cases in which it becomes ill-behaved.  Two situations are known to cause       !
         ! trouble:                                                                        !
         ! 1.  If the derivative is tiny, the next guess can be too far from the actual    !
         !     answer;                                                                     !
         ! 2.  For this specific problem, when zeta is too close to zero.  In this case    !
         !     the derivative will tend to infinity at this point and Newton's method is   !
         !     not going to perform well and can potentially enter in a weird behaviour or !
         !     lead to the wrong answer.  In any case, so we rather go with bisection.     !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler) then
            exit newloop
         elseif(stable .and. (zetaz - fun/deriv < zetamin)) then
            exit newloop
         elseif((.not. stable) .and. (zetaz - fun/deriv > zetamax)) then
            exit newloop
         end if

         !----- Copy the previous guess ---------------------------------------------------!
         zetaa = zetaz
         funa  = fun
         !---------------------------------------------------------------------------------!


         !----- New guess, its function and derivative evaluation -------------------------!
         zetaz    = zetaa - fun/deriv
         zeta0m   = zetaz * z0moz
         zeta0h   = zetaz * z0hoz
         fm       = lnzoz0m - psim(zetaz,stable) + psim(zeta0m,stable)
         fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
         dfmdzeta = z0moz * dpsimdzeta(zeta0m,stable) - dpsimdzeta(zetaz,stable)
         dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaz,stable)
         fun      = coeff * fm * fm / fh - zetaz
         deriv    = coeff * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / (fh * fh) - 1.
         !---------------------------------------------------------------------------------!



         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'NEWTON: itn=',itn,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm        &
         !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         converged = abs(zetaz-zetaa) < toler * abs(zetaz)

         if (converged) then
            zoobukhov = 0.5 * (zetaa+zetaz)
            return
         elseif (fun == 0.0) then !---- Converged by luck. --------------------------------!
            zoobukhov = zetaz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     If we reached this point then it's because Newton's method failed or it has    !
      ! become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a  !
      ! fancier bisection.  For this we need two guesses, and the guesses must have        !
      ! opposite signs.                                                                    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.0) then
         funz  = fun
         zside = .true. 
      else
         if (abs(fun-funa) < 100. * toler * abs(zetaa)) then
            if (stable) then
               delta = max(0.5 * abs(zetaa-zetamin),100. * toler * abs(zetaa))
            else
               delta = max(0.5 * abs(zetaa-zetamax),100. * toler * abs(zetaa))
            end if
         else
            if (stable) then
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,100. * toler * abs(zetaa)                                       &
                          ,0.5 * abs(zetaa-zetamin))
            else
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,100. * toler * abs(zetaa)                                       &
                          ,0.5 * abs(zetaa-zetamax))
            end if
         end if
         if (stable) then
            zetaz = max(zetamin,zetaa + delta)
         else
            zetaz = min(zetamax,zetaa + delta)
         end if
         zside = .false.
         zgssloop: do itp=1,maxfpo
            if (stable) then
               zetaz    = max(zetamin,zetaa + real((-1)**itp * (itp+3)/2) * delta)
            else
               zetaz    = min(zetamax,zetaa + real((-1)**itp * (itp+3)/2) * delta)
            end if
            zeta0m   = zetaz * z0moz
            zeta0h   = zetaz * z0hoz
            fm       = lnzoz0m - psim(zetaz,stable) + psim(zeta0m,stable)
            fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
            funz     = coeff * fm * fm / fh - zetaz
            zside    = funa * funz < 0.0
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                 &
            !   '2NDGSS: itp=',itp,'zside=',zside,'zetaa=',zetaa,'zetaz=',zetaz             &
            !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(a)') '    No second guess for you...'
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zstar  =',zstar  ,'rough  =',rough
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'lnzoz0m=',lnzoz0m,'lnzoz0h=',lnzoz0h
            write (unit=*,fmt='(1(a,1x,es14.7,1x))') 'rib    =',rib    ,'ribuse =',ribuse
            write (unit=*,fmt='(1(a,1x,l1,1x))')     'stable =',stable
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'fun    =',fun    ,'delta  =',delta
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaa  =',zetaa  ,'funa   =',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaz  =',zetaz  ,'funz   =',funz
            call abort_run('Failed finding the second guess for regula falsi'              &
                            ,'zoobukhov','leaf_coms.f90')
         end if
      end if

      !----- Now we are ready to start the regula falsi method. ---------------------------!
      bisloop: do itb=itn,maxfpo
         zoobukhov = (funz*zetaa-funa*zetaz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(zoobukhov-zetaa) < toler * abs(zoobukhov)
         if (converged) exit bisloop

         !------ Update function evaluation. ----------------------------------------------!
         zeta0m   = zoobukhov * z0moz
         zeta0h   = zoobukhov * z0hoz
         fm       = lnzoz0m - psim(zoobukhov,stable) + psim(zeta0m,stable)
         fh       = lnzoz0h - psih(zoobukhov,stable) + psih(zeta0h,stable)
         fun      = coeff * fm * fm / fh - zoobukhov
         !---------------------------------------------------------------------------------!

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'REGULA: itn=',itb,'bisection=',.true.,'zetaa=',zetaa,'zetaz=',zetaz,'fun=',fun   &
         !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

         !------ Define new interval based on the intermediate value theorem. -------------!
         if (fun*funa < 0. ) then
            zetaz = zoobukhov
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 0.5
            !----- We just updated zside, set zside to true. ------------------------------!
            zside = .true.
         else
            zetaa = zoobukhov
            funa  = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 0.5
            !----- We just updated aside, set aside to true. ------------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Zeta finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rib             [   ---] =',rib
         write (unit=*,fmt='(a,1x,f12.4)' ) 'ribuse          [   ---] =',ribuse
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zstar           [     m] =',zstar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rough           [     m] =',rough
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0m           [   ---] =',zoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0m         [   ---] =',lnzoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0h           [   ---] =',zoz0h
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0h         [   ---] =',lnzoz0h
         write (unit=*,fmt='(a,1x,l1)'    ) 'stable          [   T|F] =',stable
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaa           [   ---] =',zetaa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaz           [   ---] =',zetaz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [   ---] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fm              [   ---] =',fm
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fh              [   ---] =',fh
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [   ---] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [   ---] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [   ---] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [   ---] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [   ---] ='                   &
                                                            ,abs(zetaz-zetaa)/abs(zetaz)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoobukhov       [   ---] =',zoobukhov
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

         call abort_run('Zeta didn''t converge, giving up!!!'                              &
                         ,'zoobukhov','leaf_coms.f90')
      end if

      return
   end function zoobukhov
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil moisture to soil matric potential.                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function leaf3_matric_potential(nsoil,soil_water)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: nsoil      ! Soil texture                         [  idx]
      real(kind=4), intent(in) :: soil_water ! Soil moisture                        [m3/m3]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmoist   ! Relative soil moisture               [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmoist      = min(soil_water/slmsts(nsoil),1.0)
      !------------------------------------------------------------------------------------!



      !----- Find the matric potential. ---------------------------------------------------!
      leaf3_matric_potential = slpots(nsoil) / relmoist ** slbs(nsoil)
      !------------------------------------------------------------------------------------!

      return
   end function leaf3_matric_potential
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil matric potential to soil moisture.                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function leaf3_matric_potential_inv(nsoil,smpot)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: nsoil      ! Soil texture                         [  idx]
      real(kind=4), intent(in) :: smpot      ! Soil moisture                        [    m]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmatpot  ! Relative matric potential            [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmatpot = max(smpot/slpots(nsoil),1.0)
      !------------------------------------------------------------------------------------!



      !----- Find the matric potential. ---------------------------------------------------!
      leaf3_matric_potential_inv = slmsts(nsoil) / relmatpot ** (1. / slbs(nsoil))
      !------------------------------------------------------------------------------------!

      return
   end function leaf3_matric_potential_inv
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil moisture to hydraulic conductivity.                  !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function leaf3_hydr_conduct(k,nsoil,soil_water,soil_fracliq)
      use rconstants, only : lnexp_min ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: k            ! Layer index                        [  idx]
      integer     , intent(in) :: nsoil        ! Soil texture                       [  idx]
      real(kind=4), intent(in) :: soil_water   ! Soil moisture                      [m3/m3]
      real(kind=4), intent(in) :: soil_fracliq ! Liquid fraction                    [  ---]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmoist     ! Relative soil moisture             [  ---]
      real(kind=4)             :: fzcorr       ! Freezing correction                [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find correction for frozen soils. -------------------------------------------!
      fzcorr = exp( max( lnexp_min, - freezecoef * (1.0 - soil_fracliq) ) )
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmoist = min(soil_water/slmsts(nsoil),1.0)
      !------------------------------------------------------------------------------------!



      !----- Find the hydraulic conductivity. ---------------------------------------------!
      leaf3_hydr_conduct = fzcorr * slcons1(k,nsoil) * relmoist ** (2. * slbs(nsoil) + 3.)
      !------------------------------------------------------------------------------------!


      return
   end function leaf3_hydr_conduct
   !=======================================================================================!
   !=======================================================================================!
end module leaf_coms
!==========================================================================================!
!==========================================================================================!



