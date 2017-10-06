!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Fri Oct  6 16:24:59 2017
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E O H N C 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 42, species_electrons = 42, species_length = 9, reactions_max = 67, reactions_length = 21
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-03
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 40, bolsig_species_length = 9, bolsig_rates_max = 67, &
                                               bolsig_addsect_max = 33 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0, &
                                               bolsig_addsect(2,bolsig_addsect_max) 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1/
  data species_name(1:species_max) &
  /"O2       ","O2(V1)   ","O2(V2)   ","O2(V3)   ","O2(V4)   ","O2(A1DG) ","N2       ","N2(V1)   ","N2(V2)   ","N2(V3)   ",&
   "N2(V4)   ","N2(V5)   ","N2(V6)   ","N2(V7)   ","N2(V8)   ","H2       ","H2(V1)   ","H2(V2)   ","H2(V3)   ","CO2      ",&
   "CO2(V010)","CO2(V020)","CO2(V100)","CO       ","CO(V1)   ","CO(V2)   ","CO(V3)   ","CO(V4)   ","CO(V5)   ","CO(V6)   ",&
   "CO(V7)   ","CO(V8)   ","CO(V9)   ","CO(V10)  ","H2O      ","H2O(V010)","H2O(V101)","H2O(ROT) ","CH4      ","CH4(V13) ",&
   "CH4(V24) ","E        "/
  data reaction_sign(1:reactions_max) &
  /"bolsig:O2->O2(A1DG)  ","bolsig:N2->N2(V1)    ","bolsig:N2->N2(V2)    ","bolsig:N2->N2(V3)    ","bolsig:N2->N2(V4)    ",&
   "bolsig:N2->N2(V5)    ","bolsig:N2->N2(V6)    ","bolsig:N2->N2(V7)    ","bolsig:N2->N2(V8)    ","bolsig:N2(V1)->N2    ",&
   "bolsig:N2(V2)->N2    ","bolsig:N2(V3)->N2    ","bolsig:N2(V4)->N2    ","bolsig:N2(V5)->N2    ","bolsig:N2(V6)->N2    ",&
   "bolsig:N2(V7)->N2    ","bolsig:N2(V8)->N2    ","bolsig:O2->O2(V1)    ","bolsig:O2->O2(V2)    ","bolsig:O2->O2(V3)    ",&
   "bolsig:O2->O2(V4)    ","bolsig:O2(V1)->O2    ","bolsig:O2(V2)->O2    ","bolsig:O2(V3)->O2    ","bolsig:O2(V4)->O2    ",&
   "bolsig:CH4->CH4(V13) ","bolsig:CH4->CH4(V24) ","bolsig:CH4(V13)->CH4 ","bolsig:CH4(V24)->CH4 ","bolsig:CO2->CO2(V010)",&
   "bolsig:CO2->CO2(V020)","bolsig:CO2->CO2(V100)","bolsig:CO2(V010)->CO2","bolsig:CO2(V020)->CO2","bolsig:CO2(V100)->CO2",&
   "bolsig:CO->CO(V1)    ","bolsig:CO->CO(V2)    ","bolsig:CO->CO(V3)    ","bolsig:CO->CO(V4)    ","bolsig:CO->CO(V5)    ",&
   "bolsig:CO->CO(V6)    ","bolsig:CO->CO(V7)    ","bolsig:CO->CO(V8)    ","bolsig:CO->CO(V9)    ","bolsig:CO->CO(V10)   ",&
   "bolsig:CO(V1)->CO    ","bolsig:CO(V2)->CO    ","bolsig:CO(V3)->CO    ","bolsig:CO(V4)->CO    ","bolsig:CO(V5)->CO    ",&
   "bolsig:CO(V6)->CO    ","bolsig:CO(V7)->CO    ","bolsig:CO(V8)->CO    ","bolsig:CO(V9)->CO    ","bolsig:CO(V10)->CO   ",&
   "bolsig:H2O->H2O(V010)","bolsig:H2O->H2O(V101)","bolsig:H2O->H2O(ROT) ","bolsig:H2O(V010)->H2O","bolsig:H2O(V101)->H2O",&
   "bolsig:H2O(ROT)->H2O ","bolsig:H2->H2(V1)    ","bolsig:H2->H2(V2)    ","bolsig:H2->H2(V3)    ","bolsig:H2(V1)->H2    ",&
   "bolsig:H2(V2)->H2    ","bolsig:H2(V3)->H2    "/
  data bolsig_species(1:bolsig_species_max) &
  /"N2       ","N2(V1)   ","N2(V2)   ","N2(V3)   ","N2(V4)   ","N2(V5)   ","N2(V6)   ","N2(V7)   ","N2(V8)   ","O2       ",&
   "O2(V1)   ","O2(V2)   ","O2(V3)   ","O2(V4)   ","CH4      ","CH4(V13) ","CH4(V24) ","CO2      ","CO2(V010)","CO2(V020)",&
   "CO2(V100)","CO       ","CO(V1)   ","CO(V2)   ","CO(V3)   ","CO(V4)   ","CO(V5)   ","CO(V6)   ","CO(V7)   ","CO(V8)   ",&
   "CO(V9)   ","CO(V10)  ","H2O      ","H2O(V010)","H2O(V101)","H2O(ROT) ","H2       ","H2(V1)   ","H2(V2)   ","H2(V3)   "/
  data bolsig_addsect(1:2,1:bolsig_addsect_max) &
  / 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9,10,11,10,12,10,13,10,14,15,16,15,17,18,19,18,20,18,21,22,23,22,24,22,25,22,26,&
   22,27,22,28,22,29,22,30,22,31,22,32,33,36,33,34,33,35,37,38,37,39,37,40/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A9)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i2,1x,A21)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(251x,42(1x,i9))
312 format(A2,1x,A21,1x,42(1x,A9))
313 format(i2,1x,A21,1x,42(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(67(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,42(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,67(111x,i2.2))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,42(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(68(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(01,01) = - reac_rate_local(01) 
  reac_source_local(06,01) = + reac_rate_local(01) 
  reac_source_local(07,02) = - reac_rate_local(02) 
  reac_source_local(08,02) = + reac_rate_local(02) 
  reac_source_local(07,03) = - reac_rate_local(03) 
  reac_source_local(09,03) = + reac_rate_local(03) 
  reac_source_local(07,04) = - reac_rate_local(04) 
  reac_source_local(10,04) = + reac_rate_local(04) 
  reac_source_local(07,05) = - reac_rate_local(05) 
  reac_source_local(11,05) = + reac_rate_local(05) 
  reac_source_local(07,06) = - reac_rate_local(06) 
  reac_source_local(12,06) = + reac_rate_local(06) 
  reac_source_local(07,07) = - reac_rate_local(07) 
  reac_source_local(13,07) = + reac_rate_local(07) 
  reac_source_local(07,08) = - reac_rate_local(08) 
  reac_source_local(14,08) = + reac_rate_local(08) 
  reac_source_local(07,09) = - reac_rate_local(09) 
  reac_source_local(15,09) = + reac_rate_local(09) 
  reac_source_local(07,10) = + reac_rate_local(10) 
  reac_source_local(08,10) = - reac_rate_local(10) 
  reac_source_local(07,11) = + reac_rate_local(11) 
  reac_source_local(09,11) = - reac_rate_local(11) 
  reac_source_local(07,12) = + reac_rate_local(12) 
  reac_source_local(10,12) = - reac_rate_local(12) 
  reac_source_local(07,13) = + reac_rate_local(13) 
  reac_source_local(11,13) = - reac_rate_local(13) 
  reac_source_local(07,14) = + reac_rate_local(14) 
  reac_source_local(12,14) = - reac_rate_local(14) 
  reac_source_local(07,15) = + reac_rate_local(15) 
  reac_source_local(13,15) = - reac_rate_local(15) 
  reac_source_local(07,16) = + reac_rate_local(16) 
  reac_source_local(14,16) = - reac_rate_local(16) 
  reac_source_local(07,17) = + reac_rate_local(17) 
  reac_source_local(15,17) = - reac_rate_local(17) 
  reac_source_local(01,18) = - reac_rate_local(18) 
  reac_source_local(02,18) = + reac_rate_local(18) 
  reac_source_local(01,19) = - reac_rate_local(19) 
  reac_source_local(03,19) = + reac_rate_local(19) 
  reac_source_local(01,20) = - reac_rate_local(20) 
  reac_source_local(04,20) = + reac_rate_local(20) 
  reac_source_local(01,21) = - reac_rate_local(21) 
  reac_source_local(05,21) = + reac_rate_local(21) 
  reac_source_local(01,22) = + reac_rate_local(22) 
  reac_source_local(02,22) = - reac_rate_local(22) 
  reac_source_local(01,23) = + reac_rate_local(23) 
  reac_source_local(03,23) = - reac_rate_local(23) 
  reac_source_local(01,24) = + reac_rate_local(24) 
  reac_source_local(04,24) = - reac_rate_local(24) 
  reac_source_local(01,25) = + reac_rate_local(25) 
  reac_source_local(05,25) = - reac_rate_local(25) 
  reac_source_local(39,26) = - reac_rate_local(26) 
  reac_source_local(40,26) = + reac_rate_local(26) 
  reac_source_local(39,27) = - reac_rate_local(27) 
  reac_source_local(41,27) = + reac_rate_local(27) 
  reac_source_local(39,28) = + reac_rate_local(28) 
  reac_source_local(40,28) = - reac_rate_local(28) 
  reac_source_local(39,29) = + reac_rate_local(29) 
  reac_source_local(41,29) = - reac_rate_local(29) 
  reac_source_local(20,30) = - reac_rate_local(30) 
  reac_source_local(21,30) = + reac_rate_local(30) 
  reac_source_local(20,31) = - reac_rate_local(31) 
  reac_source_local(22,31) = + reac_rate_local(31) 
  reac_source_local(20,32) = - reac_rate_local(32) 
  reac_source_local(23,32) = + reac_rate_local(32) 
  reac_source_local(20,33) = + reac_rate_local(33) 
  reac_source_local(21,33) = - reac_rate_local(33) 
  reac_source_local(20,34) = + reac_rate_local(34) 
  reac_source_local(22,34) = - reac_rate_local(34) 
  reac_source_local(20,35) = + reac_rate_local(35) 
  reac_source_local(23,35) = - reac_rate_local(35) 
  reac_source_local(24,36) = - reac_rate_local(36) 
  reac_source_local(25,36) = + reac_rate_local(36) 
  reac_source_local(24,37) = - reac_rate_local(37) 
  reac_source_local(26,37) = + reac_rate_local(37) 
  reac_source_local(24,38) = - reac_rate_local(38) 
  reac_source_local(27,38) = + reac_rate_local(38) 
  reac_source_local(24,39) = - reac_rate_local(39) 
  reac_source_local(28,39) = + reac_rate_local(39) 
  reac_source_local(24,40) = - reac_rate_local(40) 
  reac_source_local(29,40) = + reac_rate_local(40) 
  reac_source_local(24,41) = - reac_rate_local(41) 
  reac_source_local(30,41) = + reac_rate_local(41) 
  reac_source_local(24,42) = - reac_rate_local(42) 
  reac_source_local(31,42) = + reac_rate_local(42) 
  reac_source_local(24,43) = - reac_rate_local(43) 
  reac_source_local(32,43) = + reac_rate_local(43) 
  reac_source_local(24,44) = - reac_rate_local(44) 
  reac_source_local(33,44) = + reac_rate_local(44) 
  reac_source_local(24,45) = - reac_rate_local(45) 
  reac_source_local(34,45) = + reac_rate_local(45) 
  reac_source_local(24,46) = + reac_rate_local(46) 
  reac_source_local(25,46) = - reac_rate_local(46) 
  reac_source_local(24,47) = + reac_rate_local(47) 
  reac_source_local(26,47) = - reac_rate_local(47) 
  reac_source_local(24,48) = + reac_rate_local(48) 
  reac_source_local(27,48) = - reac_rate_local(48) 
  reac_source_local(24,49) = + reac_rate_local(49) 
  reac_source_local(28,49) = - reac_rate_local(49) 
  reac_source_local(24,50) = + reac_rate_local(50) 
  reac_source_local(29,50) = - reac_rate_local(50) 
  reac_source_local(24,51) = + reac_rate_local(51) 
  reac_source_local(30,51) = - reac_rate_local(51) 
  reac_source_local(24,52) = + reac_rate_local(52) 
  reac_source_local(31,52) = - reac_rate_local(52) 
  reac_source_local(24,53) = + reac_rate_local(53) 
  reac_source_local(32,53) = - reac_rate_local(53) 
  reac_source_local(24,54) = + reac_rate_local(54) 
  reac_source_local(33,54) = - reac_rate_local(54) 
  reac_source_local(24,55) = + reac_rate_local(55) 
  reac_source_local(34,55) = - reac_rate_local(55) 
  reac_source_local(35,56) = - reac_rate_local(56) 
  reac_source_local(36,56) = + reac_rate_local(56) 
  reac_source_local(35,57) = - reac_rate_local(57) 
  reac_source_local(37,57) = + reac_rate_local(57) 
  reac_source_local(35,58) = - reac_rate_local(58) 
  reac_source_local(38,58) = + reac_rate_local(58) 
  reac_source_local(35,59) = + reac_rate_local(59) 
  reac_source_local(36,59) = - reac_rate_local(59) 
  reac_source_local(35,60) = + reac_rate_local(60) 
  reac_source_local(37,60) = - reac_rate_local(60) 
  reac_source_local(35,61) = + reac_rate_local(61) 
  reac_source_local(38,61) = - reac_rate_local(61) 
  reac_source_local(16,62) = - reac_rate_local(62) 
  reac_source_local(17,62) = + reac_rate_local(62) 
  reac_source_local(16,63) = - reac_rate_local(63) 
  reac_source_local(18,63) = + reac_rate_local(63) 
  reac_source_local(16,64) = - reac_rate_local(64) 
  reac_source_local(19,64) = + reac_rate_local(64) 
  reac_source_local(16,65) = + reac_rate_local(65) 
  reac_source_local(17,65) = - reac_rate_local(65) 
  reac_source_local(16,66) = + reac_rate_local(66) 
  reac_source_local(18,66) = - reac_rate_local(66) 
  reac_source_local(16,67) = + reac_rate_local(67) 
  reac_source_local(19,67) = - reac_rate_local(67) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(43)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(01) = rrt(01) * density(01) * density(42) 
  rrt(02) = rrt(02) * density(07) * density(42) 
  rrt(03) = rrt(03) * density(07) * density(42) 
  rrt(04) = rrt(04) * density(07) * density(42) 
  rrt(05) = rrt(05) * density(07) * density(42) 
  rrt(06) = rrt(06) * density(07) * density(42) 
  rrt(07) = rrt(07) * density(07) * density(42) 
  rrt(08) = rrt(08) * density(07) * density(42) 
  rrt(09) = rrt(09) * density(07) * density(42) 
  rrt(10) = rrt(10) * density(08) * density(42) 
  rrt(11) = rrt(11) * density(09) * density(42) 
  rrt(12) = rrt(12) * density(10) * density(42) 
  rrt(13) = rrt(13) * density(11) * density(42) 
  rrt(14) = rrt(14) * density(12) * density(42) 
  rrt(15) = rrt(15) * density(13) * density(42) 
  rrt(16) = rrt(16) * density(14) * density(42) 
  rrt(17) = rrt(17) * density(15) * density(42) 
  rrt(18) = rrt(18) * density(01) * density(42) 
  rrt(19) = rrt(19) * density(01) * density(42) 
  rrt(20) = rrt(20) * density(01) * density(42) 
  rrt(21) = rrt(21) * density(01) * density(42) 
  rrt(22) = rrt(22) * density(02) * density(42) 
  rrt(23) = rrt(23) * density(03) * density(42) 
  rrt(24) = rrt(24) * density(04) * density(42) 
  rrt(25) = rrt(25) * density(05) * density(42) 
  rrt(26) = rrt(26) * density(39) * density(42) 
  rrt(27) = rrt(27) * density(39) * density(42) 
  rrt(28) = rrt(28) * density(40) * density(42) 
  rrt(29) = rrt(29) * density(41) * density(42) 
  rrt(30) = rrt(30) * density(20) * density(42) 
  rrt(31) = rrt(31) * density(20) * density(42) 
  rrt(32) = rrt(32) * density(20) * density(42) 
  rrt(33) = rrt(33) * density(21) * density(42) 
  rrt(34) = rrt(34) * density(22) * density(42) 
  rrt(35) = rrt(35) * density(23) * density(42) 
  rrt(36) = rrt(36) * density(24) * density(42) 
  rrt(37) = rrt(37) * density(24) * density(42) 
  rrt(38) = rrt(38) * density(24) * density(42) 
  rrt(39) = rrt(39) * density(24) * density(42) 
  rrt(40) = rrt(40) * density(24) * density(42) 
  rrt(41) = rrt(41) * density(24) * density(42) 
  rrt(42) = rrt(42) * density(24) * density(42) 
  rrt(43) = rrt(43) * density(24) * density(42) 
  rrt(44) = rrt(44) * density(24) * density(42) 
  rrt(45) = rrt(45) * density(24) * density(42) 
  rrt(46) = rrt(46) * density(25) * density(42) 
  rrt(47) = rrt(47) * density(26) * density(42) 
  rrt(48) = rrt(48) * density(27) * density(42) 
  rrt(49) = rrt(49) * density(28) * density(42) 
  rrt(50) = rrt(50) * density(29) * density(42) 
  rrt(51) = rrt(51) * density(30) * density(42) 
  rrt(52) = rrt(52) * density(31) * density(42) 
  rrt(53) = rrt(53) * density(32) * density(42) 
  rrt(54) = rrt(54) * density(33) * density(42) 
  rrt(55) = rrt(55) * density(34) * density(42) 
  rrt(56) = rrt(56) * density(35) * density(42) 
  rrt(57) = rrt(57) * density(35) * density(42) 
  rrt(58) = rrt(58) * density(35) * density(42) 
  rrt(59) = rrt(59) * density(36) * density(42) 
  rrt(60) = rrt(60) * density(37) * density(42) 
  rrt(61) = rrt(61) * density(38) * density(42) 
  rrt(62) = rrt(62) * density(16) * density(42) 
  rrt(63) = rrt(63) * density(16) * density(42) 
  rrt(64) = rrt(64) * density(16) * density(42) 
  rrt(65) = rrt(65) * density(17) * density(42) 
  rrt(66) = rrt(66) * density(18) * density(42) 
  rrt(67) = rrt(67) * density(19) * density(42) 
  ydot(01) = -rrt(01)-rrt(18)-rrt(19)-rrt(20)-rrt(21)+rrt(22)+rrt(23)+rrt(24)+rrt(25) 
  ydot(02) = +rrt(18)-rrt(22) 
  ydot(03) = +rrt(19)-rrt(23) 
  ydot(04) = +rrt(20)-rrt(24) 
  ydot(05) = +rrt(21)-rrt(25) 
  ydot(06) = +rrt(01) 
  ydot(07) = -rrt(02)-rrt(03)-rrt(04)-rrt(05)-rrt(06)-rrt(07)-rrt(08)-rrt(09)+rrt(10)+rrt(11)+rrt(12)+rrt(13)+rrt(14)+rrt(15)&
             +rrt(16)+rrt(17) 
  ydot(08) = +rrt(02)-rrt(10) 
  ydot(09) = +rrt(03)-rrt(11) 
  ydot(10) = +rrt(04)-rrt(12) 
  ydot(11) = +rrt(05)-rrt(13) 
  ydot(12) = +rrt(06)-rrt(14) 
  ydot(13) = +rrt(07)-rrt(15) 
  ydot(14) = +rrt(08)-rrt(16) 
  ydot(15) = +rrt(09)-rrt(17) 
  ydot(16) = -rrt(62)-rrt(63)-rrt(64)+rrt(65)+rrt(66)+rrt(67) 
  ydot(17) = +rrt(62)-rrt(65) 
  ydot(18) = +rrt(63)-rrt(66) 
  ydot(19) = +rrt(64)-rrt(67) 
  ydot(20) = -rrt(30)-rrt(31)-rrt(32)+rrt(33)+rrt(34)+rrt(35) 
  ydot(21) = +rrt(30)-rrt(33) 
  ydot(22) = +rrt(31)-rrt(34) 
  ydot(23) = +rrt(32)-rrt(35) 
  ydot(24) = -rrt(36)-rrt(37)-rrt(38)-rrt(39)-rrt(40)-rrt(41)-rrt(42)-rrt(43)-rrt(44)-rrt(45)+rrt(46)+rrt(47)+rrt(48)+rrt(49)&
             +rrt(50)+rrt(51)+rrt(52)+rrt(53)+rrt(54)+rrt(55) 
  ydot(25) = +rrt(36)-rrt(46) 
  ydot(26) = +rrt(37)-rrt(47) 
  ydot(27) = +rrt(38)-rrt(48) 
  ydot(28) = +rrt(39)-rrt(49) 
  ydot(29) = +rrt(40)-rrt(50) 
  ydot(30) = +rrt(41)-rrt(51) 
  ydot(31) = +rrt(42)-rrt(52) 
  ydot(32) = +rrt(43)-rrt(53) 
  ydot(33) = +rrt(44)-rrt(54) 
  ydot(34) = +rrt(45)-rrt(55) 
  ydot(35) = -rrt(56)-rrt(57)-rrt(58)+rrt(59)+rrt(60)+rrt(61) 
  ydot(36) = +rrt(56)-rrt(59) 
  ydot(37) = +rrt(57)-rrt(60) 
  ydot(38) = +rrt(58)-rrt(61) 
  ydot(39) = -rrt(26)-rrt(27)+rrt(28)+rrt(29) 
  ydot(40) = +rrt(26)-rrt(28) 
  ydot(41) = +rrt(27)-rrt(29) 
  ydot(42) = 0.0d0
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(43) = 0.0d0
  if( lgas_heating ) then
    ydot(43) = ( ZDPlasKin_cfg(14)/k_B + ydot(43) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(43) = ydot(43) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(43)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) - rrt(01) * density(42) 
  pd(01,42) = pd(01,42) - rrt(01) * density(01) 
  pd(06,01) = pd(06,01) + rrt(01) * density(42) 
  pd(06,42) = pd(06,42) + rrt(01) * density(01) 
  pd(07,07) = pd(07,07) - rrt(02) * density(42) 
  pd(07,42) = pd(07,42) - rrt(02) * density(07) 
  pd(08,07) = pd(08,07) + rrt(02) * density(42) 
  pd(08,42) = pd(08,42) + rrt(02) * density(07) 
  pd(07,07) = pd(07,07) - rrt(03) * density(42) 
  pd(07,42) = pd(07,42) - rrt(03) * density(07) 
  pd(09,07) = pd(09,07) + rrt(03) * density(42) 
  pd(09,42) = pd(09,42) + rrt(03) * density(07) 
  pd(07,07) = pd(07,07) - rrt(04) * density(42) 
  pd(07,42) = pd(07,42) - rrt(04) * density(07) 
  pd(10,07) = pd(10,07) + rrt(04) * density(42) 
  pd(10,42) = pd(10,42) + rrt(04) * density(07) 
  pd(07,07) = pd(07,07) - rrt(05) * density(42) 
  pd(07,42) = pd(07,42) - rrt(05) * density(07) 
  pd(11,07) = pd(11,07) + rrt(05) * density(42) 
  pd(11,42) = pd(11,42) + rrt(05) * density(07) 
  pd(07,07) = pd(07,07) - rrt(06) * density(42) 
  pd(07,42) = pd(07,42) - rrt(06) * density(07) 
  pd(12,07) = pd(12,07) + rrt(06) * density(42) 
  pd(12,42) = pd(12,42) + rrt(06) * density(07) 
  pd(07,07) = pd(07,07) - rrt(07) * density(42) 
  pd(07,42) = pd(07,42) - rrt(07) * density(07) 
  pd(13,07) = pd(13,07) + rrt(07) * density(42) 
  pd(13,42) = pd(13,42) + rrt(07) * density(07) 
  pd(07,07) = pd(07,07) - rrt(08) * density(42) 
  pd(07,42) = pd(07,42) - rrt(08) * density(07) 
  pd(14,07) = pd(14,07) + rrt(08) * density(42) 
  pd(14,42) = pd(14,42) + rrt(08) * density(07) 
  pd(07,07) = pd(07,07) - rrt(09) * density(42) 
  pd(07,42) = pd(07,42) - rrt(09) * density(07) 
  pd(15,07) = pd(15,07) + rrt(09) * density(42) 
  pd(15,42) = pd(15,42) + rrt(09) * density(07) 
  pd(07,08) = pd(07,08) + rrt(10) * density(42) 
  pd(07,42) = pd(07,42) + rrt(10) * density(08) 
  pd(08,08) = pd(08,08) - rrt(10) * density(42) 
  pd(08,42) = pd(08,42) - rrt(10) * density(08) 
  pd(07,09) = pd(07,09) + rrt(11) * density(42) 
  pd(07,42) = pd(07,42) + rrt(11) * density(09) 
  pd(09,09) = pd(09,09) - rrt(11) * density(42) 
  pd(09,42) = pd(09,42) - rrt(11) * density(09) 
  pd(07,10) = pd(07,10) + rrt(12) * density(42) 
  pd(07,42) = pd(07,42) + rrt(12) * density(10) 
  pd(10,10) = pd(10,10) - rrt(12) * density(42) 
  pd(10,42) = pd(10,42) - rrt(12) * density(10) 
  pd(07,11) = pd(07,11) + rrt(13) * density(42) 
  pd(07,42) = pd(07,42) + rrt(13) * density(11) 
  pd(11,11) = pd(11,11) - rrt(13) * density(42) 
  pd(11,42) = pd(11,42) - rrt(13) * density(11) 
  pd(07,12) = pd(07,12) + rrt(14) * density(42) 
  pd(07,42) = pd(07,42) + rrt(14) * density(12) 
  pd(12,12) = pd(12,12) - rrt(14) * density(42) 
  pd(12,42) = pd(12,42) - rrt(14) * density(12) 
  pd(07,13) = pd(07,13) + rrt(15) * density(42) 
  pd(07,42) = pd(07,42) + rrt(15) * density(13) 
  pd(13,13) = pd(13,13) - rrt(15) * density(42) 
  pd(13,42) = pd(13,42) - rrt(15) * density(13) 
  pd(07,14) = pd(07,14) + rrt(16) * density(42) 
  pd(07,42) = pd(07,42) + rrt(16) * density(14) 
  pd(14,14) = pd(14,14) - rrt(16) * density(42) 
  pd(14,42) = pd(14,42) - rrt(16) * density(14) 
  pd(07,15) = pd(07,15) + rrt(17) * density(42) 
  pd(07,42) = pd(07,42) + rrt(17) * density(15) 
  pd(15,15) = pd(15,15) - rrt(17) * density(42) 
  pd(15,42) = pd(15,42) - rrt(17) * density(15) 
  pd(01,01) = pd(01,01) - rrt(18) * density(42) 
  pd(01,42) = pd(01,42) - rrt(18) * density(01) 
  pd(02,01) = pd(02,01) + rrt(18) * density(42) 
  pd(02,42) = pd(02,42) + rrt(18) * density(01) 
  pd(01,01) = pd(01,01) - rrt(19) * density(42) 
  pd(01,42) = pd(01,42) - rrt(19) * density(01) 
  pd(03,01) = pd(03,01) + rrt(19) * density(42) 
  pd(03,42) = pd(03,42) + rrt(19) * density(01) 
  pd(01,01) = pd(01,01) - rrt(20) * density(42) 
  pd(01,42) = pd(01,42) - rrt(20) * density(01) 
  pd(04,01) = pd(04,01) + rrt(20) * density(42) 
  pd(04,42) = pd(04,42) + rrt(20) * density(01) 
  pd(01,01) = pd(01,01) - rrt(21) * density(42) 
  pd(01,42) = pd(01,42) - rrt(21) * density(01) 
  pd(05,01) = pd(05,01) + rrt(21) * density(42) 
  pd(05,42) = pd(05,42) + rrt(21) * density(01) 
  pd(01,02) = pd(01,02) + rrt(22) * density(42) 
  pd(01,42) = pd(01,42) + rrt(22) * density(02) 
  pd(02,02) = pd(02,02) - rrt(22) * density(42) 
  pd(02,42) = pd(02,42) - rrt(22) * density(02) 
  pd(01,03) = pd(01,03) + rrt(23) * density(42) 
  pd(01,42) = pd(01,42) + rrt(23) * density(03) 
  pd(03,03) = pd(03,03) - rrt(23) * density(42) 
  pd(03,42) = pd(03,42) - rrt(23) * density(03) 
  pd(01,04) = pd(01,04) + rrt(24) * density(42) 
  pd(01,42) = pd(01,42) + rrt(24) * density(04) 
  pd(04,04) = pd(04,04) - rrt(24) * density(42) 
  pd(04,42) = pd(04,42) - rrt(24) * density(04) 
  pd(01,05) = pd(01,05) + rrt(25) * density(42) 
  pd(01,42) = pd(01,42) + rrt(25) * density(05) 
  pd(05,05) = pd(05,05) - rrt(25) * density(42) 
  pd(05,42) = pd(05,42) - rrt(25) * density(05) 
  pd(39,39) = pd(39,39) - rrt(26) * density(42) 
  pd(39,42) = pd(39,42) - rrt(26) * density(39) 
  pd(40,39) = pd(40,39) + rrt(26) * density(42) 
  pd(40,42) = pd(40,42) + rrt(26) * density(39) 
  pd(39,39) = pd(39,39) - rrt(27) * density(42) 
  pd(39,42) = pd(39,42) - rrt(27) * density(39) 
  pd(41,39) = pd(41,39) + rrt(27) * density(42) 
  pd(41,42) = pd(41,42) + rrt(27) * density(39) 
  pd(39,40) = pd(39,40) + rrt(28) * density(42) 
  pd(39,42) = pd(39,42) + rrt(28) * density(40) 
  pd(40,40) = pd(40,40) - rrt(28) * density(42) 
  pd(40,42) = pd(40,42) - rrt(28) * density(40) 
  pd(39,41) = pd(39,41) + rrt(29) * density(42) 
  pd(39,42) = pd(39,42) + rrt(29) * density(41) 
  pd(41,41) = pd(41,41) - rrt(29) * density(42) 
  pd(41,42) = pd(41,42) - rrt(29) * density(41) 
  pd(20,20) = pd(20,20) - rrt(30) * density(42) 
  pd(20,42) = pd(20,42) - rrt(30) * density(20) 
  pd(21,20) = pd(21,20) + rrt(30) * density(42) 
  pd(21,42) = pd(21,42) + rrt(30) * density(20) 
  pd(20,20) = pd(20,20) - rrt(31) * density(42) 
  pd(20,42) = pd(20,42) - rrt(31) * density(20) 
  pd(22,20) = pd(22,20) + rrt(31) * density(42) 
  pd(22,42) = pd(22,42) + rrt(31) * density(20) 
  pd(20,20) = pd(20,20) - rrt(32) * density(42) 
  pd(20,42) = pd(20,42) - rrt(32) * density(20) 
  pd(23,20) = pd(23,20) + rrt(32) * density(42) 
  pd(23,42) = pd(23,42) + rrt(32) * density(20) 
  pd(20,21) = pd(20,21) + rrt(33) * density(42) 
  pd(20,42) = pd(20,42) + rrt(33) * density(21) 
  pd(21,21) = pd(21,21) - rrt(33) * density(42) 
  pd(21,42) = pd(21,42) - rrt(33) * density(21) 
  pd(20,22) = pd(20,22) + rrt(34) * density(42) 
  pd(20,42) = pd(20,42) + rrt(34) * density(22) 
  pd(22,22) = pd(22,22) - rrt(34) * density(42) 
  pd(22,42) = pd(22,42) - rrt(34) * density(22) 
  pd(20,23) = pd(20,23) + rrt(35) * density(42) 
  pd(20,42) = pd(20,42) + rrt(35) * density(23) 
  pd(23,23) = pd(23,23) - rrt(35) * density(42) 
  pd(23,42) = pd(23,42) - rrt(35) * density(23) 
  pd(24,24) = pd(24,24) - rrt(36) * density(42) 
  pd(24,42) = pd(24,42) - rrt(36) * density(24) 
  pd(25,24) = pd(25,24) + rrt(36) * density(42) 
  pd(25,42) = pd(25,42) + rrt(36) * density(24) 
  pd(24,24) = pd(24,24) - rrt(37) * density(42) 
  pd(24,42) = pd(24,42) - rrt(37) * density(24) 
  pd(26,24) = pd(26,24) + rrt(37) * density(42) 
  pd(26,42) = pd(26,42) + rrt(37) * density(24) 
  pd(24,24) = pd(24,24) - rrt(38) * density(42) 
  pd(24,42) = pd(24,42) - rrt(38) * density(24) 
  pd(27,24) = pd(27,24) + rrt(38) * density(42) 
  pd(27,42) = pd(27,42) + rrt(38) * density(24) 
  pd(24,24) = pd(24,24) - rrt(39) * density(42) 
  pd(24,42) = pd(24,42) - rrt(39) * density(24) 
  pd(28,24) = pd(28,24) + rrt(39) * density(42) 
  pd(28,42) = pd(28,42) + rrt(39) * density(24) 
  pd(24,24) = pd(24,24) - rrt(40) * density(42) 
  pd(24,42) = pd(24,42) - rrt(40) * density(24) 
  pd(29,24) = pd(29,24) + rrt(40) * density(42) 
  pd(29,42) = pd(29,42) + rrt(40) * density(24) 
  pd(24,24) = pd(24,24) - rrt(41) * density(42) 
  pd(24,42) = pd(24,42) - rrt(41) * density(24) 
  pd(30,24) = pd(30,24) + rrt(41) * density(42) 
  pd(30,42) = pd(30,42) + rrt(41) * density(24) 
  pd(24,24) = pd(24,24) - rrt(42) * density(42) 
  pd(24,42) = pd(24,42) - rrt(42) * density(24) 
  pd(31,24) = pd(31,24) + rrt(42) * density(42) 
  pd(31,42) = pd(31,42) + rrt(42) * density(24) 
  pd(24,24) = pd(24,24) - rrt(43) * density(42) 
  pd(24,42) = pd(24,42) - rrt(43) * density(24) 
  pd(32,24) = pd(32,24) + rrt(43) * density(42) 
  pd(32,42) = pd(32,42) + rrt(43) * density(24) 
  pd(24,24) = pd(24,24) - rrt(44) * density(42) 
  pd(24,42) = pd(24,42) - rrt(44) * density(24) 
  pd(33,24) = pd(33,24) + rrt(44) * density(42) 
  pd(33,42) = pd(33,42) + rrt(44) * density(24) 
  pd(24,24) = pd(24,24) - rrt(45) * density(42) 
  pd(24,42) = pd(24,42) - rrt(45) * density(24) 
  pd(34,24) = pd(34,24) + rrt(45) * density(42) 
  pd(34,42) = pd(34,42) + rrt(45) * density(24) 
  pd(24,25) = pd(24,25) + rrt(46) * density(42) 
  pd(24,42) = pd(24,42) + rrt(46) * density(25) 
  pd(25,25) = pd(25,25) - rrt(46) * density(42) 
  pd(25,42) = pd(25,42) - rrt(46) * density(25) 
  pd(24,26) = pd(24,26) + rrt(47) * density(42) 
  pd(24,42) = pd(24,42) + rrt(47) * density(26) 
  pd(26,26) = pd(26,26) - rrt(47) * density(42) 
  pd(26,42) = pd(26,42) - rrt(47) * density(26) 
  pd(24,27) = pd(24,27) + rrt(48) * density(42) 
  pd(24,42) = pd(24,42) + rrt(48) * density(27) 
  pd(27,27) = pd(27,27) - rrt(48) * density(42) 
  pd(27,42) = pd(27,42) - rrt(48) * density(27) 
  pd(24,28) = pd(24,28) + rrt(49) * density(42) 
  pd(24,42) = pd(24,42) + rrt(49) * density(28) 
  pd(28,28) = pd(28,28) - rrt(49) * density(42) 
  pd(28,42) = pd(28,42) - rrt(49) * density(28) 
  pd(24,29) = pd(24,29) + rrt(50) * density(42) 
  pd(24,42) = pd(24,42) + rrt(50) * density(29) 
  pd(29,29) = pd(29,29) - rrt(50) * density(42) 
  pd(29,42) = pd(29,42) - rrt(50) * density(29) 
  pd(24,30) = pd(24,30) + rrt(51) * density(42) 
  pd(24,42) = pd(24,42) + rrt(51) * density(30) 
  pd(30,30) = pd(30,30) - rrt(51) * density(42) 
  pd(30,42) = pd(30,42) - rrt(51) * density(30) 
  pd(24,31) = pd(24,31) + rrt(52) * density(42) 
  pd(24,42) = pd(24,42) + rrt(52) * density(31) 
  pd(31,31) = pd(31,31) - rrt(52) * density(42) 
  pd(31,42) = pd(31,42) - rrt(52) * density(31) 
  pd(24,32) = pd(24,32) + rrt(53) * density(42) 
  pd(24,42) = pd(24,42) + rrt(53) * density(32) 
  pd(32,32) = pd(32,32) - rrt(53) * density(42) 
  pd(32,42) = pd(32,42) - rrt(53) * density(32) 
  pd(24,33) = pd(24,33) + rrt(54) * density(42) 
  pd(24,42) = pd(24,42) + rrt(54) * density(33) 
  pd(33,33) = pd(33,33) - rrt(54) * density(42) 
  pd(33,42) = pd(33,42) - rrt(54) * density(33) 
  pd(24,34) = pd(24,34) + rrt(55) * density(42) 
  pd(24,42) = pd(24,42) + rrt(55) * density(34) 
  pd(34,34) = pd(34,34) - rrt(55) * density(42) 
  pd(34,42) = pd(34,42) - rrt(55) * density(34) 
  pd(35,35) = pd(35,35) - rrt(56) * density(42) 
  pd(35,42) = pd(35,42) - rrt(56) * density(35) 
  pd(36,35) = pd(36,35) + rrt(56) * density(42) 
  pd(36,42) = pd(36,42) + rrt(56) * density(35) 
  pd(35,35) = pd(35,35) - rrt(57) * density(42) 
  pd(35,42) = pd(35,42) - rrt(57) * density(35) 
  pd(37,35) = pd(37,35) + rrt(57) * density(42) 
  pd(37,42) = pd(37,42) + rrt(57) * density(35) 
  pd(35,35) = pd(35,35) - rrt(58) * density(42) 
  pd(35,42) = pd(35,42) - rrt(58) * density(35) 
  pd(38,35) = pd(38,35) + rrt(58) * density(42) 
  pd(38,42) = pd(38,42) + rrt(58) * density(35) 
  pd(35,36) = pd(35,36) + rrt(59) * density(42) 
  pd(35,42) = pd(35,42) + rrt(59) * density(36) 
  pd(36,36) = pd(36,36) - rrt(59) * density(42) 
  pd(36,42) = pd(36,42) - rrt(59) * density(36) 
  pd(35,37) = pd(35,37) + rrt(60) * density(42) 
  pd(35,42) = pd(35,42) + rrt(60) * density(37) 
  pd(37,37) = pd(37,37) - rrt(60) * density(42) 
  pd(37,42) = pd(37,42) - rrt(60) * density(37) 
  pd(35,38) = pd(35,38) + rrt(61) * density(42) 
  pd(35,42) = pd(35,42) + rrt(61) * density(38) 
  pd(38,38) = pd(38,38) - rrt(61) * density(42) 
  pd(38,42) = pd(38,42) - rrt(61) * density(38) 
  pd(16,16) = pd(16,16) - rrt(62) * density(42) 
  pd(16,42) = pd(16,42) - rrt(62) * density(16) 
  pd(17,16) = pd(17,16) + rrt(62) * density(42) 
  pd(17,42) = pd(17,42) + rrt(62) * density(16) 
  pd(16,16) = pd(16,16) - rrt(63) * density(42) 
  pd(16,42) = pd(16,42) - rrt(63) * density(16) 
  pd(18,16) = pd(18,16) + rrt(63) * density(42) 
  pd(18,42) = pd(18,42) + rrt(63) * density(16) 
  pd(16,16) = pd(16,16) - rrt(64) * density(42) 
  pd(16,42) = pd(16,42) - rrt(64) * density(16) 
  pd(19,16) = pd(19,16) + rrt(64) * density(42) 
  pd(19,42) = pd(19,42) + rrt(64) * density(16) 
  pd(16,17) = pd(16,17) + rrt(65) * density(42) 
  pd(16,42) = pd(16,42) + rrt(65) * density(17) 
  pd(17,17) = pd(17,17) - rrt(65) * density(42) 
  pd(17,42) = pd(17,42) - rrt(65) * density(17) 
  pd(16,18) = pd(16,18) + rrt(66) * density(42) 
  pd(16,42) = pd(16,42) + rrt(66) * density(18) 
  pd(18,18) = pd(18,18) - rrt(66) * density(42) 
  pd(18,42) = pd(18,42) - rrt(66) * density(18) 
  pd(16,19) = pd(16,19) + rrt(67) * density(42) 
  pd(16,42) = pd(16,42) + rrt(67) * density(19) 
  pd(19,19) = pd(19,19) - rrt(67) * density(42) 
  pd(19,42) = pd(19,42) - rrt(67) * density(19) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(43,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(43,:) = pd(43,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  call ZDPlasKin_bolsig_rates()
  rrt(01) = bolsig_rates(bolsig_pointer(1))
  rrt(02) = bolsig_rates(bolsig_pointer(2))
  rrt(03) = bolsig_rates(bolsig_pointer(3))
  rrt(04) = bolsig_rates(bolsig_pointer(4))
  rrt(05) = bolsig_rates(bolsig_pointer(5))
  rrt(06) = bolsig_rates(bolsig_pointer(6))
  rrt(07) = bolsig_rates(bolsig_pointer(7))
  rrt(08) = bolsig_rates(bolsig_pointer(8))
  rrt(09) = bolsig_rates(bolsig_pointer(9))
  rrt(10) = bolsig_rates(bolsig_pointer(10))
  rrt(11) = bolsig_rates(bolsig_pointer(11))
  rrt(12) = bolsig_rates(bolsig_pointer(12))
  rrt(13) = bolsig_rates(bolsig_pointer(13))
  rrt(14) = bolsig_rates(bolsig_pointer(14))
  rrt(15) = bolsig_rates(bolsig_pointer(15))
  rrt(16) = bolsig_rates(bolsig_pointer(16))
  rrt(17) = bolsig_rates(bolsig_pointer(17))
  rrt(18) = bolsig_rates(bolsig_pointer(18))
  rrt(19) = bolsig_rates(bolsig_pointer(19))
  rrt(20) = bolsig_rates(bolsig_pointer(20))
  rrt(21) = bolsig_rates(bolsig_pointer(21))
  rrt(22) = bolsig_rates(bolsig_pointer(22))
  rrt(23) = bolsig_rates(bolsig_pointer(23))
  rrt(24) = bolsig_rates(bolsig_pointer(24))
  rrt(25) = bolsig_rates(bolsig_pointer(25))
  rrt(26) = bolsig_rates(bolsig_pointer(26))
  rrt(27) = bolsig_rates(bolsig_pointer(27))
  rrt(28) = bolsig_rates(bolsig_pointer(28))
  rrt(29) = bolsig_rates(bolsig_pointer(29))
  rrt(30) = bolsig_rates(bolsig_pointer(30))
  rrt(31) = bolsig_rates(bolsig_pointer(31))
  rrt(32) = bolsig_rates(bolsig_pointer(32))
  rrt(33) = bolsig_rates(bolsig_pointer(33))
  rrt(34) = bolsig_rates(bolsig_pointer(34))
  rrt(35) = bolsig_rates(bolsig_pointer(35))
  rrt(36) = bolsig_rates(bolsig_pointer(36))
  rrt(37) = bolsig_rates(bolsig_pointer(37))
  rrt(38) = bolsig_rates(bolsig_pointer(38))
  rrt(39) = bolsig_rates(bolsig_pointer(39))
  rrt(40) = bolsig_rates(bolsig_pointer(40))
  rrt(41) = bolsig_rates(bolsig_pointer(41))
  rrt(42) = bolsig_rates(bolsig_pointer(42))
  rrt(43) = bolsig_rates(bolsig_pointer(43))
  rrt(44) = bolsig_rates(bolsig_pointer(44))
  rrt(45) = bolsig_rates(bolsig_pointer(45))
  rrt(46) = bolsig_rates(bolsig_pointer(46))
  rrt(47) = bolsig_rates(bolsig_pointer(47))
  rrt(48) = bolsig_rates(bolsig_pointer(48))
  rrt(49) = bolsig_rates(bolsig_pointer(49))
  rrt(50) = bolsig_rates(bolsig_pointer(50))
  rrt(51) = bolsig_rates(bolsig_pointer(51))
  rrt(52) = bolsig_rates(bolsig_pointer(52))
  rrt(53) = bolsig_rates(bolsig_pointer(53))
  rrt(54) = bolsig_rates(bolsig_pointer(54))
  rrt(55) = bolsig_rates(bolsig_pointer(55))
  rrt(56) = bolsig_rates(bolsig_pointer(56))
  rrt(57) = bolsig_rates(bolsig_pointer(57))
  rrt(58) = bolsig_rates(bolsig_pointer(58))
  rrt(59) = bolsig_rates(bolsig_pointer(59))
  rrt(60) = bolsig_rates(bolsig_pointer(60))
  rrt(61) = bolsig_rates(bolsig_pointer(61))
  rrt(62) = bolsig_rates(bolsig_pointer(62))
  rrt(63) = bolsig_rates(bolsig_pointer(63))
  rrt(64) = bolsig_rates(bolsig_pointer(64))
  rrt(65) = bolsig_rates(bolsig_pointer(65))
  rrt(66) = bolsig_rates(bolsig_pointer(66))
  rrt(67) = bolsig_rates(bolsig_pointer(67))
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
