program artes

  use omp_lib

  implicit none

  ! Constants
  integer,  parameter       :: dp       = kind(0.d0)                  ! Double precision
  real(dp), parameter       :: pi       = 4._dp*datan(1._dp)          ! Surface area of a disk with radius unity
  real(dp), parameter       :: k_b      = 1.3806488e-23_dp            ! Boltzmann constant [m2 kg s-2 K-1]
  real(dp), parameter       :: sb       = 5.670373e-8_dp              ! Stefan-Boltzmann constant [J s-1 m-2 K-4]
  real(dp), parameter       :: hh       = 6.62606957e-34_dp           ! Planck constant [m2 kg s-1]
  real(dp), parameter       :: cc       = 2.99792458e8_dp             ! Speed of light [m s-1]
  real(dp), parameter       :: r_sun    = 6.95500e8_dp                ! Solar radius [m]
  real(dp), parameter       :: pc       = 3.08572e16_dp               ! Parsec [m]
  real(dp), parameter       :: au       = 1.49598e11_dp               ! Astronomical Unit [m]

  ! Input files
  logical                   :: log_file                               ! Log file on or output to the screen
  integer                   :: photon_source                          ! Photon emission point (1=star, 2=planet)
  logical                   :: phase_curve
  logical                   :: spectrum
  logical                   :: imaging_mono
  logical                   :: imaging_broad
  real(dp)                  :: wavelength                             ! Photon wavelength
  integer(16)               :: packages                               ! Number of photon packages
  real(dp)                  :: orbit                                  ! Star-planet distance [m]
  real(dp)                  :: distance_planet                        ! Distance from the observer to the planet [m]
  real(dp)                  :: fstop                                  ! Photon stop parameter [0:1]
  real(dp)                  :: photon_minimum                         ! Minimum photon energy possible (0<=photon_minimum<=1), otherwise removed
  logical                   :: thermal_weight                         ! Weight cell luminosity
  logical                   :: photon_scattering                      ! Scattering on or off
  integer                   :: photon_emission                        ! Isotropic (1) or biased (2) emission
  real(dp)                  :: photon_bias                            ! Bias parameter (0 <= bias < 1)
  real(dp)                  :: t_star                                 ! Stellar effective temperature [K]
  real(dp)                  :: r_star                                 ! Stellar radius [m]
  logical                   :: stellar_direction                      ! Stellar direction set manual (on) or automatic (off)
  real(dp)                  :: theta_star                             ! Stellar location, theta direction
  real(dp)                  :: phi_star                               ! Stellar location, phi direction
  real(dp)                  :: surface_albedo                         ! Surface albedo [0,1] (0=black, 1=Lambertian)
  real(dp)                  :: oblateness                             ! Planet oblateness (0-1)
  real(dp)                  :: oblate_x, oblate_y, oblate_z           ! Planet aspect_ratio (=R_equator/R_poles)
  integer                   :: nr                                     ! Number of grid points in r-direction
  integer                   :: ntheta                                 ! Number of grid points in theta-direction
  integer                   :: nphi                                   ! Number of grid points in phi-direction
  real(dp)                  :: det_theta, det_phi                     ! Detector direction
  real(dp)                  :: theta_phase, phi_phase                 ! Theta and phi coordinate of phase curve detector
  integer                   :: nx                                     ! Number of pixels in x-direction of high-res detector
  integer                   :: ny                                     ! Number of pixels in y-direction of high-res detector
  real(dp)                  :: x_max                                  ! Image plane size ranges from [-x_max:x_max]
  real(dp)                  :: y_max                                  ! Image plane size ranges from [-y_max:y_max]
  character(100)            :: email                                  ! Send email when job is finished
  logical                   :: flow_global                            ! Global energy transport output (on/off)
  logical                   :: flow_theta                             ! Latitudinal energy transport output (on/off)
  logical                   :: ring                                   ! Ring system (on/off)

  ! Arrays
  real(dp), allocatable     :: rfront(:)                              ! Radial front point of the grid cell
  real(dp), allocatable     :: thetafront(:)                          ! Theta front point of the grid cell
  integer,  allocatable     :: thetaplane(:)                          ! Type of theta surface (1=normal/cone, 2=flat/xy-plane)
  real(dp), allocatable     :: phifront(:)                            ! Phi front point of the grid cell
  real(dp), allocatable     :: cell_density(:,:,:)                    ! Cell density [kg/m3]
  real(dp), allocatable     :: cell_temperature(:,:,:)                ! Cell temperature [K]
  real(dp), allocatable     :: cell_scattering_opacity(:,:,:,:)       ! Cell scattering opacity [m-1]
  real(dp), allocatable     :: cell_absorption_opacity(:,:,:,:)       ! Cell absorption opacity [m-1]
  real(dp), allocatable     :: cell_opacity(:,:,:,:)                  ! Cell total opacity [m-1]
  real(dp), allocatable     :: cell_luminosity(:,:,:)                 ! Cell thermal luminosity [W]
  real(dp), allocatable     :: cell_weight(:,:,:)                     ! Cell luminosity weight
  real(dp), allocatable     :: cell_scatter_matrix(:,:,:,:,:,:)       ! Cell scatter matrix number
  real(dp), allocatable     :: cell_albedo(:,:,:,:)                   ! Cell albedo
  real(dp), allocatable     :: cell_volume(:,:,:)                     ! Cell volume [m3]
  real(dp), allocatable     :: cell_p11_int(:,:,:,:)                  ! Integral of the P11 scattering matrix element
  real(dp), allocatable     :: cell_p12_int(:,:,:,:)                  ! Integral of the P12 scattering matrix element
  real(dp), allocatable     :: cell_p13_int(:,:,:,:)                  ! Integral of the P13 scattering matrix element
  real(dp), allocatable     :: cell_p14_int(:,:,:,:)                  ! Integral of the P14 scattering matrix element
  real(dp), allocatable     :: theta_grid_cos(:)                      ! Cosine of theta boundaries
  real(dp), allocatable     :: theta_grid_tan(:)                      ! Tangent of theta boundaries
  real(dp), allocatable     :: phi_grid_sin(:)                        ! Sine of phi boundaries
  real(dp), allocatable     :: phi_grid_cos(:)                        ! Cosine of phi boundaries
  real(dp), allocatable     :: wavelengths(:)                         ! Wavelength points [micron]
  real(dp), allocatable     :: cell_flow_global(:,:,:,:,:)            ! Global energy transport through each cell
  real(dp), allocatable     :: cell_flow(:,:,:,:,:)                   ! Energy transport through cell boundary (1=upward, 2=downward, 3=south, 4=north)
  real(dp), allocatable     :: emissivity_cumulative(:,:,:)           ! Cumulative emissivity in each grid cell CDF
  real(dp), allocatable     :: detector(:,:,:,:)                      ! Array(flux/flux^2/photons, Stokes I/Q/U/V, nx, ny)
  real(dp), allocatable     :: detector_thread(:,:,:,:,:)             ! Array(thread_id, flux/flux^2/photons, Stokes I/Q/U/V, nx, ny)
  real(dp), allocatable     :: flux_emitted(:)                        ! Thermal flux emitted in the atmosphere
  real(dp), allocatable     :: flux_exit(:)                           ! Flux exiting from the atmosphere
  real(dp)                  :: sinbeta(360), cosbeta(360)             ! sin(beta), cos(beta)
  real(dp)                  :: sin2beta(360), cos2beta(360)           ! sin(2*beta), cos(2*beta)
  real(dp)                  :: photometry(11)                         ! Detector integrated values
  real(dp)                  :: det_dir(5)

  ! Global variables
  real(dp)                  :: phase_observer                         ! Phase angle of observation [degrees]
  real(dp)                  :: sin_det_theta, cos_det_theta
  real(dp)                  :: sin_det_phi, cos_det_phi
  integer                   :: out_unit                               ! Output unit, 10=logfile, 6=screen
  integer                   :: cores                                  ! Number of computer cores
  integer                   :: threads                                ! Number of threads in use
  integer                   :: cells                                  ! Number of grid cells
  integer                   :: cell_depth                             ! Maximum cell depth to emit thermal photons
  real(dp)                  :: t_start, t_end
  character(100)            :: output_name
  character(100)            :: atmosphere                             ! Atmosphere name
  character(100)            :: atmosphere_directory                   ! Directory with atmosphere input files
  character(100)            :: input_file                             ! Input parameter file
  real(dp)                  :: package_energy                         ! Photon package energy [J]
  real(dp)                  :: x_fov, y_fov                           ! Field of view [mas]
  real(dp)                  :: pixel_scale                            ! Pixel size [mas pixel-1]
  integer                   :: n_wavelength                           ! Number of wavelengths
  integer                   :: wl_count                               ! Wavelength point
  character(300)            :: error_log                              ! Filename of error log
  integer, parameter        :: ns = 4                                 ! Random number generator
  integer, parameter        :: default_seed(ns) = (/ 521288629, 362436069, 16163801, 1131199299 /)
  integer, allocatable      :: state(:,:)                             ! Random number generator

  call run

contains

  subroutine run

    integer  :: i
    
    call initialize
    call python
    call grid_initialize(1)
    call output(1)

    wl_count = 1

    if (spectrum) then

       do

          wavelength = wavelengths(wl_count)

          call array_start
          call grid_initialize(2)

          if (wl_count.eq.1) call output(2)

          write (6,fmt="(a1,a,t13,f7.3,a)", advance="no") achar(13), &
               "Wavelength: ", wavelength*1.e6_dp, " micron"

          call radiative_transfer
          call write_output

          call grid_finished(1)
          call grid_finished(2)

          if (wl_count.eq.size(wavelengths)) then

             write (6,fmt="(a1,a,t13,f6.3,a)", advance="yes") achar(13), ""
             if (.not.log_file) write (6,'(a)') ""
             call output(3)
             exit

          end if

          wl_count = wl_count + 1

       end do

       call grid_finished(3)
       call output(4)

    else if (imaging_broad) then

       do

          wavelength = wavelengths(wl_count)
          call grid_initialize(2)

          if (wl_count.eq.1) then

             call array_start
             call output(2)

          end if

          write (6,fmt="(a1,a,t13,f6.3,a)", advance="no") achar(13), &
               "Wavelength: ", wavelength*1.e6_dp, " micron"

          call radiative_transfer
          call grid_finished(1)

          if (wl_count.eq.size(wavelengths)) then

             write (6,fmt="(a1,a,t13,f6.3,a)", advance="yes") achar(13), ""
             if (.not.log_file) write (6,'(a)') ""
             exit

          end if

          wl_count = wl_count + 1

       end do

       call write_output
       call output(3)
       call grid_finished(2)
       call grid_finished(3)
       call output(4)

    else if (phase_curve.or.imaging_mono) then

       wavelength = wavelengths(wl_count)
       
       call grid_initialize(2)
       call output(2)

       if (phase_curve) then

          do i = 1,74

             if (i.eq.1) then
                det_phi = 1.e-5_dp*pi/180._dp
             else if (i.eq.2) then
                det_phi = 2.5_dp*pi/180._dp
             else if (i.eq.73) then
                det_phi = (180._dp-1e-5_dp)*pi/180._dp
             else if (i.eq.74) then
                write (6,fmt="(a1,a,t14,f6.1,a)", advance="yes") achar(13), ""
                if (.not.log_file) write (6,'(a)') ""
                exit
             else
                det_phi = det_phi + 2.5_dp*pi/180._dp
             end if

             det_dir(5) = det_phi
             sin_det_phi = sin(det_dir(5))
             cos_det_phi = cos(det_dir(5))

             call spherical_cartesian(1._dp, det_theta, det_phi, det_dir(1), det_dir(2), det_dir(3))

             write (6,fmt="(a1,a,t14,f6.1,a)", advance="no") achar(13), &
                  "Phase angle: ", det_phi*180._dp/pi, " degrees"

             call array_start
             call radiative_transfer
             call write_output
             call grid_finished(2)

          end do

          call output(3)
          call grid_finished(1)
          call grid_finished(3)
          call output(4)

       else if (imaging_mono) then

          call array_start
          call radiative_transfer
          call write_output
          call output(3)
          call grid_finished(1)
          call grid_finished(2)
          call grid_finished(3)
          call output(4)

       end if

    end if

  end subroutine run

  subroutine initialize

    integer                     :: i, j, nlines
    logical                     :: file_exists
    character(1)                :: first_char
    character(100)              :: key_word, key_value
    character(500), allocatable :: data(:)
    real(dp)                    :: xi

    ! ARTES input file
    
    log_file = .false.
    email = ""

    photon_source = 1
    packages = 100000
    fstop = 1.e-5_dp
    photon_minimum = 1.e-20_dp
    thermal_weight = .true.
    photon_scattering = .true.
    photon_emission = 1
    photon_bias = 0.8

    t_star = 5800._dp
    r_star = r_sun
    stellar_direction = .false.
    theta_star = pi/2._dp
    phi_star = 0._dp

    surface_albedo = 0._dp
    oblateness = 0._dp
    orbit = 5._dp*au
    ring = .false.
    
    phase_curve = .false.
    spectrum = .false.
    imaging_mono = .false.
    imaging_broad = .false.
    det_theta = 90._dp
    det_phi = 90._dp
    nx = 25
    ny = 25
    distance_planet = 10._dp*pc

    flow_global = .false.
    flow_theta = .false.

    ! Other parameters
    
    det_dir = 0._dp
    phase_observer = 0._dp
    atmosphere = ""
    atmosphere_directory = ""
    input_file = ""
    package_energy = 0._dp
    i = 0
    j = 0
    nlines = 0
    first_char = ""
    file_exists = .false.
    wavelength = 0._dp
    x_max = 0._dp
    y_max = 0._dp
    oblate_x = 1._dp
    oblate_y = 1._dp
    oblate_z = 1._dp
    cell_depth = 0
    photometry = 0._dp

    inquire(file="/proc/cpuinfo", exist=file_exists)

    if (file_exists) then

       ! Check number of cores on a Linux computer

       call system("grep -c ^processor /proc/cpuinfo > cores.txt")
       open (100, file='cores.txt')
       read (100,*) cores
       close (100, status='delete')

    else if (.not.file_exists) then

       ! Check number of cores on a Apple computer

       call system("sysctl hw.ncpu | awk '{print $2}' > cores.txt")
       open (100, file='cores.txt')
       read (100,*) cores
       close (100, status='delete')

    end if
    
    ! Number of threads in use
    !$omp parallel
    threads = omp_get_num_threads()
    !$omp end parallel

    ! CPU start time
    t_start = omp_get_wtime()

    ! Get input directory and number of photons
    call argument_input(atmosphere)
    atmosphere_directory = 'input/'//trim(atmosphere)
    input_file = trim(atmosphere_directory)//'/artes.in'

    ! Check if input file exists
    inquire (file=input_file, exist=file_exists)
    if (.not.file_exists) then
       write (6,'(a)') "Input file does not exist!"
       call exit(0)
    end if

    ! Read the input file
    call readfile(input_file, data, nlines)

    ! Get the keywords and values, ignore comment lines
    do i=1,nlines

       ! Get first character of string
       first_char = data(i)(1:1)

       ! Check if the line is not a comment line
       if (first_char.ne."*".and.first_char.ne."-".and.first_char.ne."=".and.len_trim(data(i)).ne.0) then

          call get_key_value(data(i), key_word, key_value)
          call input_parameters(key_word, key_value)

       end if

    enddo

    ! Get command line keywords
    call argument_keywords
    
    ! Sine and cosine arrays

    cosbeta  = 0._dp
    sinbeta  = 0._dp
    cos2beta = 0._dp
    sin2beta = 0._dp

    do i=1,180

       cosbeta(i)      = ( cos(dble(i)*pi/180._dp) + cos(dble(i-1)*pi/180._dp) ) / 2._dp
       cosbeta(i+180)  = -cosbeta(i)
       sinbeta(i)      = ( sin(dble(i)*pi/180._dp) + sin(dble(i-1)*pi/180._dp) ) / 2._dp
       sinbeta(i+180)  = -sinbeta(i)
       cos2beta(i)     = ( cos(2._dp*dble(i)*pi/180._dp) + cos(2._dp*dble(i-1)*pi/180._dp) ) / 2._dp
       cos2beta(i+180) = cos2beta(i)
       sin2beta(i)     = ( sin(2._dp*dble(i)*pi/180._dp) + sin(2._dp*dble(i-1)*pi/180._dp) ) / 2._dp
       sin2beta(i+180) = sin2beta(i)

    end do

    ! Get atmospheric structure
    call get_atmosphere
    
    ! Open error log
    error_log = 'output/'//trim(output_name)//'/error.log'
    open (11, file=trim(error_log), status="new")
    close (11)

    ! Open output log
    if (log_file) open (10, file='output/'//trim(output_name)//'/output.log', status='new')

    ! Random number generator

    call init_random_seed

    allocate (state(threads,ns))
    state = 0.
    
    do j=1,ns
       do i=1,threads
          if (j.eq.1) then
             call random_number(xi)
             state(i,j) = int(xi*1.e6_dp)
          else
             state(i,j) = default_seed(j)
          end if
       end do
    end do

    ! Detector

    if (spectrum) then

       nx = 1
       ny = 1

    else if (phase_curve) then

       nx = 1
       ny = 1
       det_theta = pi/2._dp
       det_phi = 1.e-5_dp

    end if
    
    ! Oblateness

    oblate_x = 1._dp / (1._dp-oblateness)
    oblate_y = oblate_x
    oblate_z = 1._dp

    ! Detector field of view

    x_max = 1.3_dp*rfront(nr)
    y_max = 1.3_dp*rfront(nr)

    x_max = ( oblateness + 1._dp ) * x_max
    y_max = ( oblateness + 1._dp ) * y_max

    ! Field of view [mas]
 
    x_fov = 2._dp*atan(x_max/distance_planet)*3600._dp*180._dp/pi*1000._dp
    y_fov = 2._dp*atan(y_max/distance_planet)*3600._dp*180._dp/pi*1000._dp
    
    ! Pixel size [mas pixel-1]
    
    pixel_scale = x_fov/nx

    ! Detectors

    if (abs(det_phi).lt.1.e-3_dp.or.det_phi.gt.2._dp*pi-1.e-3_dp) det_phi = 1.e-3_dp
    if (det_phi.gt.pi-1.e-3_dp.and.det_phi.lt.pi+1.e-3_dp) det_phi = pi-1.e-3_dp

    call spherical_cartesian(1._dp, det_theta, det_phi, det_dir(1), det_dir(2), det_dir(3))
    det_dir(4) = det_theta
    det_dir(5) = det_phi

    sin_det_theta = sin(det_dir(4))
    cos_det_theta = cos(det_dir(4))
    sin_det_phi = sin(det_dir(5))
    cos_det_phi = cos(det_dir(5))
    
    ! Observer phase angle

    if (.not.phase_curve) then

       phase_observer = sin(theta_star)*cos(phi_star)*sin(det_theta)*cos(det_phi) + &
            sin(theta_star)*sin(phi_star)*sin(det_theta)*sin(det_phi) + &
            cos(theta_star)*cos(det_theta)

       phase_observer = acos(phase_observer) * 180._dp/pi

    end if
    
  end subroutine initialize

  subroutine radiative_transfer

    integer(16) :: i
    integer     :: j, k, l, thread_id, cell(3), cell_out(3), current_face(2), next_face(2)
    integer     :: current_face_check(2), cell_check(3), buffer(13), status
    logical     :: grid_exit, surface, check1, check2, check3, check4, cell_error
    real(dp)    :: face_distance, tau, tau_cell, tau_run, xi, x_out, y_out, z_out, tau_first
    real(dp)    :: stokes(4), stokes_new(4), alpha, beta, scatter(4,4), gamma
    real(dp)    :: direction(3), direction_new(3), x, y, z, s, x_check, y_check, z_check
    real(dp)    :: r_disk, phi_disk, dpi, dummy, bias_weight
    
    check1 = .false.
    check2 = .false.
    check3 = .false.
    check4 = .false.

    !$omp parallel do default(none) private(surface) &

    !$omp& private(x, y, z, s, x_out, y_out, z_out, stokes, stokes_new, direction, direction_new) &
    !$omp& private(current_face, next_face, face_distance, alpha, beta, cell, cell_out, tau, tau_run, tau_first) &
    !$omp& private(scatter, xi, theta_phase, phi_phase, grid_exit, tau_cell, cell_error, thread_id, bias_weight) &
    !$omp& private(x_check, y_check, z_check, current_face_check, cell_check, buffer, status, r_disk, phi_disk) &

    !$omp& shared(cell_opacity, cell_absorption_opacity, cell_albedo, nr, gamma, photon_source, log_file, photon_minimum) &
    !$omp& shared(surface_albedo, packages, threads, cell_depth, fstop, spectrum, imaging_mono, imaging_broad) &
    !$omp& shared(check1, check2, check3, check4, det_phi, det_theta, output_name, wl_count, rfront, cell_scattering_opacity) &
    !$omp& shared(flow_global, flow_theta, package_energy, cell_weight, flux_emitted, flux_exit, error_log, photon_scattering)
    
    do i=1,packages

       ! Check if error file is not becoming larger than 100 MB
       call stat('output/'//trim(output_name)//'/error.log', buffer, status)
       if (buffer(8).gt.1.e8_dp) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 001"
          close (11)
          call exit(0)
       end if

       surface = .false.
       grid_exit = .false.
       cell_error = .false.
       cell = 0
       cell_out = 0
       current_face = 0
       next_face = 0
       x = 0._dp
       y = 0._dp
       z = 0._dp
       bias_weight = 0._dp

       thread_id = OMP_GET_THREAD_NUM()

       if (thread_id.eq.0.and.imaging_mono) then

          if (dble(i)/(dble(packages)/dble(threads)).gt.0.2_dp.and..not.check1) then
             write (6,fmt="(a3)", advance="no") "20%"
             check1 = .true.
          else if (dble(i)/(dble(packages)/dble(threads)).gt.0.4_dp.and..not.check2) then
             write (6,fmt="(a9)", advance="no") "  --  40%"
             check2 = .true.
          else if (dble(i)/(dble(packages)/dble(threads)).gt.0.6_dp.and..not.check3) then
             write (6,fmt="(a9)", advance="no") "  --  60%"
             check3 = .true.
          else if (dble(i)/(dble(packages)/dble(threads)).gt.0.8_dp.and..not.check4) then
             write (6,fmt="(a9)", advance="no") "  --  80%"
             check4 = .true.
          else if (i.eq.packages/threads) then
             write (6,fmt='(a10)', advance="yes") "  --  100%"
             if (.not.log_file) write (6,'(a)') ""
          end if

       end if

       call emit_photon(thread_id, x, y, z, direction, current_face, cell, r_disk, phi_disk, bias_weight)

       stokes(1) = 1._dp
       stokes(2) = 0._dp
       stokes(3) = 0._dp
       stokes(4) = 0._dp

       if (photon_source.eq.2) then

          ! Correction of the photon energy for the cell emission probability
          ! stokes(1) = stokes(1) * exp( cell_absorption_opacity(cell(1),cell(2),cell(3),wl_count) * &
          !     ( rfront(cell(1)+1) - rfront(cell(1)) ) )

          stokes(1) = stokes(1) * bias_weight / cell_weight(cell(1),cell(2),cell(3))

          flux_emitted(thread_id+1) = flux_emitted(thread_id+1) + stokes(1)

          call peel_thermal(thread_id, x, y, z, stokes, cell, current_face, cell_error)

          if (cell_error) then

             open (11, file=trim(error_log), position="append")
             write (11,*) "error 047"
             close (11)

             cycle
             
          end if
          
       end if

       ! Get optical depth to grid boundary or surface

       tau_first = 0._dp

       x_check = x
       y_check = y
       z_check = z
       current_face_check = current_face
       cell_check = cell

       do

          call cell_face(x_check, y_check, z_check, direction, current_face_check, &
               next_face, face_distance, grid_exit, cell_check, cell_out, cell_error)

          if (cell_error) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 002"
             close (11)
          end if
          
          tau_cell = face_distance*cell_opacity(cell_check(1), cell_check(2), cell_check(3), wl_count)
          tau_first = tau_first + tau_cell

          x_check = x_check + face_distance*direction(1)
          y_check = y_check + face_distance*direction(2)
          z_check = z_check + face_distance*direction(3)

          if (cell_error.or.grid_exit.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

          current_face_check = next_face
          cell_check = cell_out

       end do
       
       ! First optical depth

       if (tau_first.lt.1.e-6_dp.and..not.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) then

          ! New photon in case optical depth is zero and not crossing the surface
          
          cycle

       else if (tau_first.lt.1.e-6_dp.and.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) then

          ! Zero optical depth, photon hits planet surface

          call random(thread_id,xi)
          tau = -log(1._dp-xi)

       else

          ! Sample optical depth, force first interaction, weight photon

          call random(thread_id,xi)
          if (tau_first.lt.50._dp) then
             tau = -log(1._dp-xi*(1._dp-exp(-tau_first)))
             stokes = stokes * (1._dp-exp(-tau_first))
          else
             tau = -log(1._dp-xi)
          end if

       end if

       ! Do loop to start crossing cells till next interaction point

       tau_run = 0._dp

       do

          call cell_face(x, y, z, direction, current_face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

          if (cell_error) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 003"
             close (11)
          end if

          tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)

          ! Check if next interaction happens in the upcoming cell

          if (tau_run+tau_cell.gt.tau) then
             
             grid_exit = .false.

             s = (tau-tau_run) / cell_opacity(cell(1), cell(2), cell(3), wl_count)

             x = x + s*direction(1)
             y = y + s*direction(2)
             z = z + s*direction(3)

             if (flow_global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), s, cell)
             
             current_face = 0
             next_face = 0

             exit

          else

             x = x + face_distance*direction(1)
             y = y + face_distance*direction(2)
             z = z + face_distance*direction(3)

             if (flow_global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), face_distance, cell)

             if (flow_theta) then
                if (next_face(1).eq.1) then
                   if (cell_out(1).gt.cell(1)) then
                      call add_flow(thread_id, 1, stokes(1), cell)
                   else if (cell_out(1).lt.cell(1)) then
                      call add_flow(thread_id, 2, stokes(1), cell)
                   end if
                else if (next_face(1).eq.2) then
                   if (cell_out(2).gt.cell(2)) then
                      call add_flow(thread_id, 3, stokes(1), cell)
                   else if (cell_out(2).lt.cell(2)) then
                      call add_flow(thread_id, 4, stokes(1), cell)
                   end if
                end if
             end if
             
             current_face = next_face
             cell = cell_out

             if (grid_exit) exit

          end if
          
          ! Photon is crossing the planet surface

          if (next_face(1).eq.1.and.next_face(2).eq.cell_depth) then

             call random(thread_id,xi)

             if (xi.gt.surface_albedo) then

                surface = .true.
                exit

             else if (xi.le.surface_albedo) then

                call lambertian(thread_id, x, y, z, stokes, stokes_new, direction)
                call peel_surface(thread_id, x, y, z, stokes, cell, current_face)

                stokes = stokes_new
                cell(1) = cell(1)+1

             end if

          end if

          tau_run = tau_run + tau_cell

       end do

       if (grid_exit.and.photon_source.eq.2) flux_exit(thread_id+1) = flux_exit(thread_id+1) + stokes(1)
       
       ! Emit new photon when current photon leaves the atmosphere or is absorbed by the planet surface
       
       if (surface.or.grid_exit.or.cell_error) cycle

       ! Let the photon scatter in the atmosphere until it either exits or is absorbed
       
       do

          ! No scattering
          if (.not.photon_scattering) exit
          
          call random(thread_id,xi)

          if (xi.lt.fstop) then

             exit

          else

             if (cell_albedo(cell(1),cell(2),cell(3),wl_count).lt.1._dp.and. &
                  cell_albedo(cell(1),cell(2),cell(3),wl_count).gt.0._dp) then

                gamma = cell_albedo(cell(1),cell(2),cell(3),wl_count)/(1._dp-fstop)
                stokes = gamma*stokes

             end if

             ! Remove photon when Stokes I becomes too small
             if (stokes(1).le.photon_minimum) then
                cell_error = .true.
                exit
             end if

             call peel_photon(thread_id, x, y, z, stokes, direction, cell, current_face, cell_error)

             if (cell_error) exit
             
             call scatter_photon(thread_id, direction, stokes, cell, direction_new, alpha, beta, scatter)

             if (abs(alpha).lt.1._dp) then

                call polarization_rotation(alpha, beta, stokes, scatter, direction, direction_new, stokes_new, .false.)
                
                stokes = stokes_new
                direction = direction_new

             else
                
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 050"
                write (11,*) alpha, beta
                write (11,*) stokes
                write (11,*) direction
                close (11)

                cell_error = .true.
                exit

             end if

          end if

          ! Sample optical depth
          call random(thread_id,xi)
          tau = -log(1._dp-xi)

          tau_run = 0._dp

          do

             call cell_face(x, y, z, direction, current_face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

             if (cell_error) then
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 004"
                close (11)
             end if

             tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)

             if (tau_run+tau_cell.gt.tau) then

                ! Next interaction point is in the same grid cell

                grid_exit = .false.

                s = (tau-tau_run) / cell_opacity(cell(1), cell(2), cell(3), wl_count)

                x = x + s*direction(1)
                y = y + s*direction(2)
                z = z + s*direction(3)

                if (flow_global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), s, cell)
                
                current_face = 0
                next_face = 0

                exit

             else

                ! Photon crosses a cell face before next interaction point

                x = x + face_distance*direction(1)
                y = y + face_distance*direction(2)
                z = z + face_distance*direction(3)

                if (flow_global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), face_distance, cell)

                if (flow_theta) then
                   if (next_face(1).eq.1) then
                      if (cell_out(1).gt.cell(1)) then
                         call add_flow(thread_id, 1, stokes(1), cell)
                      else if (cell_out(1).lt.cell(1)) then
                         call add_flow(thread_id, 2, stokes(1), cell)
                      end if
                   else if (next_face(1).eq.2) then
                      if (cell_out(2).gt.cell(2)) then
                         call add_flow(thread_id, 3, stokes(1), cell)
                      else if (cell_out(2).lt.cell(2)) then
                         call add_flow(thread_id, 4, stokes(1), cell)
                      end if
                   end if
                end if
                
                current_face = next_face
                cell = cell_out

                if (grid_exit) exit

             end if

             ! Photon crosses the surface

             if (next_face(1).eq.1.and.next_face(2).eq.cell_depth) then

                call random(thread_id,xi)

                if (xi.gt.surface_albedo) then

                   surface = .true.
                   exit

                else if (xi.le.surface_albedo) then

                   call lambertian(thread_id, x, y, z, stokes, stokes_new, direction)

                   call peel_surface(thread_id, x, y, z, stokes, cell, current_face)

                   stokes = stokes_new
                   cell(1) = cell(1)+1

                end if

             end if

             tau_run = tau_run + tau_cell
             current_face = next_face

          end do

          if (cell_error) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 005"
             close (11)
          end if

          if (grid_exit.or.surface.or.cell_error) exit

       end do

       if (grid_exit.and.photon_source.eq.2) flux_exit(thread_id+1) = flux_exit(thread_id+1) + stokes(1)

    end do
    
    call photon_package
    
    do l = 1,3
       do k = 1,4
          do j = 1,ny
             do i = 1,nx

                if (l.eq.1) then
                   detector(i,j,k,l) = sum(detector_thread(i,j,k,l,:)) * package_energy
                else if (l.eq.2) then
                   detector(i,j,k,l) = sum(detector_thread(i,j,k,l,:)) * package_energy * package_energy
                else if (l.eq.3) then
                   detector(i,j,k,l) = sum(detector_thread(i,j,k,l,:))
                end if

             end do
          end do
       end do
    end do

    photometry(1) = sum(detector(:,:,1,1)) ! Stokes I
    photometry(3) = sum(detector(:,:,2,1)) ! Stokes Q
    photometry(5) = sum(detector(:,:,3,1)) ! Stokes U
    photometry(7) = sum(detector(:,:,4,1)) ! Stokes V
    photometry(9) = sqrt( sum(detector(:,:,2,1))**2 + sum(detector(:,:,3,1))**2 ) ! PI
    photometry(10) = photometry(9) / photometry(1) ! PI/I

    do i=1,4

       if (sum(detector(:,:,i,3)).gt.0._dp) then
       
          dummy = ( sum(detector(:,:,i,2)) / sum(detector(:,:,i,3)) ) - ( sum(detector(:,:,i,1)) / sum(detector(:,:,i,3)) )**2
          if (dummy.gt.0._dp) photometry(i*2) = sqrt(dummy) * sqrt(sum(detector(:,:,i,3)))

       end if
       
    end do

    ! Degree of polarization error

    if (photometry(3)**2+photometry(5)**2.gt.0._dp) then
    
       dpi = sqrt( ( (photometry(3)*photometry(4))**2 + (photometry(5)*photometry(6))**2 ) / &
            ( 2._dp*(photometry(3)**2+photometry(5)**2 ) ) )

       photometry(11) = photometry(10) * sqrt( (dpi/photometry(9))**2 + (photometry(2)/photometry(1))**2  )

    end if

  end subroutine radiative_transfer

  subroutine emit_photon(thread_id, x_emission, y_emission, z_emission, emission_direction, &
       face_emission, cell_emission, r_disk, phi_disk, bias_weight)

    ! Photon emission

    integer,  intent(in)  :: thread_id
    integer,  intent(out) :: cell_emission(3), face_emission(2)
    integer               :: i, j, k
    real(dp), intent(out) :: x_emission, y_emission, z_emission, emission_direction(3), r_disk, phi_disk, bias_weight
    real(dp)              :: alpha, beta, rot_matrix(3,3), disk(2), emissivity_cumulative_sampled, cos_phi_sampled
    real(dp)              :: r_sampled, phi_sampled, phase_obs, sin_theta_sampled, cos_theta_sampled, emis_prev, xi
    real(dp)              :: theta_direction, phi_direction, x_temp, y_temp, z_temp, cos_beta, sin_beta, sin_phi_sampled
    real(dp)              :: theta_sampled, radial_unit(3), y_bias, norm
    logical               :: solution, grid_exit, cell_ok

    solution = .false.
    grid_exit = .false.
    bias_weight = 1._dp

    if (photon_source.eq.1) then

       ! Stellar photons

       face_emission(1) = 1
       face_emission(2) = nr

       solution = .false.

       r_disk = 0._dp
       phi_disk = 0._dp

       do while (.not.solution)

          if (phase_curve.and.det_phi*180._dp/pi.ge.170._dp) then

             do

                call random(thread_id,xi)
                r_disk = sqrt(xi)
                if (r_disk.gt.0.9_dp) exit

             end do

             call random(thread_id,xi)
             phi_disk = 2._dp * pi * xi

             disk(1) = rfront(nr) * r_disk * sin(phi_disk)
             disk(2) = rfront(nr) * r_disk * cos(phi_disk)

          else

             call random(thread_id,xi)
             r_disk = sqrt(xi)

             call random(thread_id,xi)
             phi_disk = 2._dp * pi * xi

             disk(1) = rfront(nr) * r_disk * sin(phi_disk)
             disk(2) = rfront(nr) * r_disk * cos(phi_disk)

          end if
          
          solution = .true.

          emission_direction(1) = -1._dp
          emission_direction(2) = 0._dp
          emission_direction(3) = 0._dp

          x_emission = sqrt( rfront(nr)*rfront(nr) - disk(1)*disk(1) - disk(2)*disk(2) )
          y_emission = disk(1)
          z_emission = disk(2)

          if (stellar_direction) then

             ! Rotate coordinates
             
             ! Rotate coordinates around y-axis
             call rotation_matrix(2, -(pi/2._dp-theta_star), rot_matrix)

             x_temp = x_emission*rot_matrix(1,1) + y_emission*rot_matrix(1,2) + z_emission*rot_matrix(1,3)
             y_temp = x_emission*rot_matrix(2,1) + y_emission*rot_matrix(2,2) + z_emission*rot_matrix(2,3)
             z_temp = x_emission*rot_matrix(3,1) + y_emission*rot_matrix(3,2) + z_emission*rot_matrix(3,3)
             
             ! Rotate coordinates around z-axis
             call rotation_matrix(3, phi_star, rot_matrix)

             x_emission = x_temp*rot_matrix(1,1) + y_temp*rot_matrix(1,2) + z_temp*rot_matrix(1,3)
             y_emission = x_temp*rot_matrix(2,1) + y_temp*rot_matrix(2,2) + z_temp*rot_matrix(2,3)
             z_emission = x_temp*rot_matrix(3,1) + y_temp*rot_matrix(3,2) + z_temp*rot_matrix(3,3)

             ! Rotate photon direction

             theta_direction = pi - theta_star
             phi_direction = pi + phi_star

             if (theta_direction.lt.0._dp) theta_direction = theta_direction + 2._dp*pi
             if (theta_direction.gt.2._dp*pi) theta_direction = theta_direction - 2._dp*pi
             if (phi_direction.lt.0._dp) phi_direction = phi_direction + 2._dp*pi
             if (phi_direction.gt.2._dp*pi) phi_direction = phi_direction - 2._dp*pi

             call spherical_cartesian(1._dp, theta_direction, phi_direction, emission_direction(1), &
                  emission_direction(2), emission_direction(3))
             
          end if

       end do

       call initial_cell(x_emission, y_emission, z_emission, cell_emission)

    else if (photon_source.eq.2) then

       ! Thermal photon emission from the planet

       face_emission(1) = 0
       face_emission(2) = 0

       ! Sample grid cell

       call random(thread_id,xi)
       emissivity_cumulative_sampled = xi*emissivity_cumulative(nr-1,ntheta-1,nphi-1)

       emis_prev = 0._dp
       cell_ok = .false.
       
       do i=cell_depth,nr-1
          do j=0,ntheta-1
             do k=0,nphi-1

                if ( emissivity_cumulative_sampled.ge.emis_prev.and. &
                     emissivity_cumulative_sampled.le.emissivity_cumulative(i,j,k) ) then

                   cell_emission(1) = i
                   cell_emission(2) = j
                   cell_emission(3) = k

                   cell_ok = .true.
                   
                   exit

                end if

                emis_prev = emissivity_cumulative(i,j,k)

             end do
             if (cell_ok) exit
          end do
          if (cell_ok) exit
       end do

       ! Sample location in grid cell

       call random(thread_id,xi)
       r_sampled = xi * ( rfront(cell_emission(1)+1) - rfront(cell_emission(1)) )
       r_sampled = rfront(cell_emission(1)) + r_sampled

       call random(thread_id,xi)
       cos_theta_sampled = xi * ( theta_grid_cos(cell_emission(2)+1) - theta_grid_cos(cell_emission(2)) )
       cos_theta_sampled = theta_grid_cos(cell_emission(2)) + cos_theta_sampled
       sin_theta_sampled = sqrt(1._dp-cos_theta_sampled*cos_theta_sampled)

       if (nphi.eq.1) then

          call random(thread_id,xi)
          phi_sampled = 2._dp*pi*xi
          
       else if (nphi.gt.1) then

          if (cell_emission(3).lt.nphi-1) then

             call random(thread_id,xi)
             phi_sampled = xi * ( phifront(cell_emission(3)+1) - phifront(cell_emission(3)) )
             phi_sampled = phifront(cell_emission(3)) + phi_sampled

          else if (cell_emission(3).eq.nphi-1) then

             call random(thread_id,xi)
             phi_sampled = xi * ( 2._dp*pi - phifront(cell_emission(3)) )
             phi_sampled = phifront(cell_emission(3)) + phi_sampled

          end if

       end if

       cos_phi_sampled = cos(phi_sampled)
       sin_phi_sampled = sqrt(1._dp-cos_phi_sampled*cos_phi_sampled)
       if (phi_sampled.gt.pi) sin_phi_sampled = -sin_phi_sampled
       
       ! Emission coordinates
       
       x_emission = r_sampled*sin_theta_sampled*cos_phi_sampled
       y_emission = r_sampled*sin_theta_sampled*sin_phi_sampled
       z_emission = r_sampled*cos_theta_sampled

       ! Scale spherical coordinates to oblate spheroid coordinates

       x_emission = oblate_x * x_emission
       y_emission = oblate_y * y_emission
       z_emission = oblate_z * z_emission

       phase_obs = sin_theta_sampled*cos_phi_sampled*sin_det_theta*cos_det_phi + &
            sin_theta_sampled*sin_phi_sampled*sin_det_theta*sin_det_phi + cos_theta_sampled*cos_det_theta

       ! Sample emission direction

       if (photon_emission.eq.1) then

          ! Isotropic
          
          call random(thread_id,xi)
          alpha  = 2._dp*xi - 1._dp
          call random(thread_id,xi)
          beta = 2._dp*pi*xi

          cos_beta = cos(beta)
          sin_beta = sqrt(1._dp-cos_beta*cos_beta)
          if (beta.gt.pi) sin_beta = -sin_beta

          emission_direction(1) = sqrt(1._dp-alpha*alpha)*cos_beta
          emission_direction(2) = sqrt(1._dp-alpha*alpha)*sin_beta
          emission_direction(3) = alpha

       else if (photon_emission.eq.2) then

          ! Biased upward (Gordon 1987)
          ! http://www.oceanopticsbook.info/view/monte_carlo_simulation/importance_sampling

          call random(thread_id,xi)
          y_bias = (1._dp+photon_bias) * tan(pi*xi/2._dp) / sqrt(1._dp-photon_bias*photon_bias)
          theta_sampled = acos( (1._dp-y_bias*y_bias) / (1._dp+y_bias*y_bias) )

          call random(thread_id,xi)
          beta = 2._dp*pi*xi

          ! Local radial vector, normal vector to sphere
          radial_unit(1) = x_emission / ( oblate_x * oblate_x )
          radial_unit(2) = y_emission / ( oblate_y * oblate_y )
          radial_unit(3) = z_emission / ( oblate_z * oblate_z )
          norm = sqrt(radial_unit(1)*radial_unit(1)+radial_unit(2)*radial_unit(2)+radial_unit(3)*radial_unit(3))
          radial_unit = radial_unit / norm

          ! Photon emission direction
          call direction_cosine(cos(pi-theta_sampled), beta, radial_unit, emission_direction)

          ! Weight factor
          bias_weight = ( pi * sin(theta_sampled) * (1._dp+photon_bias*cos(theta_sampled)) ) / &
               ( 2._dp * sqrt(1._dp-photon_bias*photon_bias) )

       end if

       if (abs(emission_direction(3)).ge.1._dp) then

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 054"
          write (11,*) alpha, beta, xi
          close (11)

       end if

    end if

  end subroutine emit_photon
  
  subroutine rotation_matrix(rotation_axis, rotation_angle, matrix)

    ! Create rotation matrix for rotation aroung x-axis (1), y-axis (2) or z-axis (3)

    integer,  intent(in)  :: rotation_axis
    real(dp), intent(in)  :: rotation_angle
    real(dp), intent(out) :: matrix(3,3)
    real(dp)              :: cos_rot_angle, sin_rot_angle

    cos_rot_angle = cos(rotation_angle)
    sin_rot_angle = sin(rotation_angle)
    
    if (rotation_axis.eq.1) then

       matrix(1,1) = 1._dp
       matrix(1,2) = 0._dp
       matrix(1,3) = 0._dp

       matrix(2,1) = 0._dp
       matrix(2,2) = cos_rot_angle
       matrix(2,3) = -sin_rot_angle

       matrix(3,1) = 0._dp
       matrix(3,2) = sin_rot_angle
       matrix(3,3) = cos_rot_angle

    else if (rotation_axis.eq.2) then

       matrix(1,1) = cos_rot_angle
       matrix(1,2) = 0._dp
       matrix(1,3) = sin_rot_angle

       matrix(2,1) = 0._dp
       matrix(2,2) = 1._dp
       matrix(2,3) = 0._dp

       matrix(3,1) = -sin_rot_angle
       matrix(3,2) = 0._dp
       matrix(3,3) = cos_rot_angle

    else if (rotation_axis.eq.3) then

       matrix(1,1) = cos_rot_angle
       matrix(1,2) = -sin_rot_angle
       matrix(1,3) = 0._dp

       matrix(2,1) = sin_rot_angle
       matrix(2,2) = cos_rot_angle
       matrix(2,3) = 0._dp

       matrix(3,1) = 0._dp
       matrix(3,2) = 0._dp
       matrix(3,3) = 1._dp

    end if

  end subroutine rotation_matrix
  
  subroutine python

    ! Write input parameters for Python scripts

    character(100) :: dummy

    open (100, file='output/'//trim(output_name)//'/plot.dat', status='replace')
    write (100,'(a)') "[plot]"
    write (dummy,'(i10)') photon_source
    write (100,'(2a)') "photon_source=", adjustl(trim(dummy))
    write (dummy,'(es14.7)') distance_planet
    write (100,'(2a)') "distance=", adjustl(trim(dummy))
    write (dummy,'(es14.7)') rfront(0)
    write (100,'(2a)') "planet_radius=", adjustl(trim(dummy))
    write (dummy,'(i10)') ntheta
    write (100,'(2a)') "ntheta=", adjustl(trim(dummy))
    write (dummy,'(es14.7)') x_fov
    write (100,'(2a)') "fov=", adjustl(trim(dummy))
    close (100)

  end subroutine python

  subroutine planck_function(temperature, flux)

    real(dp), intent(in)  :: temperature
    real(dp), intent(out) :: flux

    if (photon_source.eq.1) then

       ! Planck function [W m-2 m-1]
       flux = (2._dp*pi*hh*cc*cc/(wavelength**5._dp)) / ( exp(hh*cc/(wavelength*k_b*temperature)) - 1._dp )
       
    else if (photon_source.eq.2) then

       ! Planck function [W m-2 m-1 sr-1]
       flux = (2._dp*hh*cc*cc/(wavelength**5._dp)) / ( exp(hh*cc/(wavelength*k_b*temperature)) - 1._dp )

    end if

  end subroutine planck_function

  subroutine lambertian(thread_id, x, y, z, stokes, stokes_new, direction)

    ! Isotropic Lambertian reflection

    integer,  intent(in)  :: thread_id
    real(dp), intent(in)  :: x, y, z, stokes(4)
    real(dp), intent(out) :: stokes_new(4), direction(3)
    real(dp)              :: xi, alpha, beta, surface_normal(3), norm

    ! Surface normal on spherical surface at point of reflection and scale for oblate surface
    surface_normal(1) = x / ( oblate_x * oblate_x )
    surface_normal(2) = y / ( oblate_y * oblate_y )
    surface_normal(3) = z / ( oblate_z * oblate_z )

    ! Make the vector a unit vector
    norm = sqrt(surface_normal(1)*surface_normal(1)+surface_normal(2)*surface_normal(2)+surface_normal(3)*surface_normal(3))
    surface_normal = surface_normal / norm

    ! Sample emission direction
    call random(thread_id,xi)
    alpha = sqrt(xi)
    call random(thread_id,xi)
    beta = 2._dp*pi*xi

    ! Calculate photon emission direction on a spherical planet
    call direction_cosine(alpha, beta, surface_normal, direction)

    ! Lambertian surface depolarizes the light 100%
    stokes_new(1) = stokes(1)
    stokes_new(2) = 0._dp
    stokes_new(3) = 0._dp
    stokes_new(4) = 0._dp

  end subroutine lambertian

  subroutine cartesian_spherical(x, y, z, r, theta, phi)

    ! Transform from Cartesian to spherical coordinates
    ! Note that atan2 returns a value in the range [-pi:pi]

    real(dp), intent(in)  :: x, y, z
    real(dp), intent(out) :: r, theta, phi

    r = sqrt(x*x+y*y+z*z)

    theta = acos(z/r)
    phi = atan2(y,x)

    if (phi.lt.0._dp) phi = phi+2._dp*pi

  end subroutine cartesian_spherical

  subroutine spherical_cartesian(r,theta,phi,x,y,z)

    ! Transform from spherical to Cartesian coordinates

    real(dp), intent(in)  :: r, theta, phi
    real(dp), intent(out) :: x, y, z

    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)

  end subroutine spherical_cartesian

  subroutine scatter_photon(thread_id, direction, stokes, cell, direction_new, alpha, beta, scatter)

    ! Sample scattering angle, calculate new direction and get the scatterings matrix

    integer               :: i, angle_upper, angle_lower
    integer,  intent(in)  :: cell(3), thread_id
    real(dp), intent(in)  :: direction(3), stokes(4)
    real(dp), intent(out) :: direction_new(3), alpha, beta, scatter(4,4)
    real(dp)              :: scatter_dummy(16), x0(16), x1(16), y0, y1, y_inter, acos_alpha

    call scattering_angle_sampling(thread_id, stokes, alpha, beta, cell)

    call direction_cosine(alpha, beta, direction, direction_new)

    acos_alpha = acos(alpha)
    
    if (mod(acos_alpha*180._dp/pi,1._dp).gt.0.5_dp) then

       angle_upper = int(acos_alpha*180._dp/pi) + 2
       angle_lower = int(acos_alpha*180._dp/pi) + 1

    else

       angle_upper = int(acos_alpha*180._dp/pi) + 1
       angle_lower = int(acos_alpha*180._dp/pi)

    end if

    if (angle_upper.eq.1) then

       scatter(1,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,1,1)
       scatter(1,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,2,1)
       scatter(1,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,3,1)
       scatter(1,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,4,1)
       scatter(2,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,5,1)
       scatter(2,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,6,1)
       scatter(2,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,7,1)
       scatter(2,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,8,1)
       scatter(3,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,9,1)
       scatter(3,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,10,1)
       scatter(3,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,11,1)
       scatter(3,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,12,1)
       scatter(4,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,13,1)
       scatter(4,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,14,1)
       scatter(4,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,15,1)
       scatter(4,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,16,1)

    else if (angle_lower.eq.180) then

       scatter(1,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,1,180)
       scatter(1,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,2,180)
       scatter(1,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,3,180)
       scatter(1,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,4,180)
       scatter(2,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,5,180)
       scatter(2,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,6,180)
       scatter(2,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,7,180)
       scatter(2,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,8,180)
       scatter(3,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,9,180)
       scatter(3,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,10,180)
       scatter(3,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,11,180)
       scatter(3,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,12,180)
       scatter(4,1) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,13,180)
       scatter(4,2) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,14,180)
       scatter(4,3) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,15,180)
       scatter(4,4) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,16,180)

    else

       do i=1,16

          x0(i) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,i,angle_lower)
          x1(i) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,i,angle_upper)
          y0 = dble(angle_lower) - 0.5_dp
          y1 = dble(angle_upper) - 0.5_dp
          y_inter = acos_alpha*180._dp/pi
          scatter_dummy(i) = (x1(i)-x0(i)) * (y_inter-y0) / (y1-y0) + x0(i)

       end do

       scatter(1,1) = scatter_dummy(1)
       scatter(1,2) = scatter_dummy(2)
       scatter(1,3) = scatter_dummy(3)
       scatter(1,4) = scatter_dummy(4)
       scatter(2,1) = scatter_dummy(5)
       scatter(2,2) = scatter_dummy(6)
       scatter(2,3) = scatter_dummy(7)
       scatter(2,4) = scatter_dummy(8)
       scatter(3,1) = scatter_dummy(9)
       scatter(3,2) = scatter_dummy(10)
       scatter(3,3) = scatter_dummy(11)
       scatter(3,4) = scatter_dummy(12)
       scatter(4,1) = scatter_dummy(13)
       scatter(4,2) = scatter_dummy(14)
       scatter(4,3) = scatter_dummy(15)
       scatter(4,4) = scatter_dummy(16)

    end if

  end subroutine scatter_photon

  subroutine scattering_angle_sampling(thread_id, stokes, alpha, beta, cell)

    integer               :: i
    integer,  intent(in)  :: cell(3), thread_id
    real(dp), intent(in)  :: stokes(4)
    real(dp), intent(out) :: alpha, beta
    real(dp)              :: intensity(180), intensity_cumulative(0:180), c2b, s2b
    real(dp)              :: intensity_cumulative_sampled, x0, x1, y0, y1, xi

    ! Azimuthal angle sampling

    intensity_cumulative(0) = 0._dp

    do i=1,180

       ! Calculate the intensity for different azimuthal angles
       
       intensity(i) = cell_p11_int(cell(1),cell(2),cell(3),wl_count)*stokes(1) + &
            cell_p12_int(cell(1),cell(2),cell(3),wl_count)*stokes(2)*cos2beta(i) + &
            cell_p12_int(cell(1),cell(2),cell(3),wl_count)*stokes(3)*sin2beta(i) - &
            cell_p13_int(cell(1),cell(2),cell(3),wl_count)*stokes(2)*sin2beta(i) + &
            cell_p13_int(cell(1),cell(2),cell(3),wl_count)*stokes(3)*cos2beta(i) + &
            cell_p14_int(cell(1),cell(2),cell(3),wl_count)*stokes(4)

       intensity_cumulative(i) = intensity_cumulative(i-1) + intensity(i)

    end do

    call random(thread_id,xi)
    intensity_cumulative_sampled = xi*intensity_cumulative(180)

    do i=1,180

       if ( intensity_cumulative_sampled.ge.intensity_cumulative(i-1) .and. &
            intensity_cumulative_sampled.le.intensity_cumulative(i) ) then

          x0 = dble(i-1)
          x1 = dble(i)
          y0 = intensity_cumulative(i-1)
          y1 = intensity_cumulative(i)
          beta = (x1-x0) * (intensity_cumulative_sampled-y0) / (y1-y0) + x0
          beta = beta*pi/180._dp

          exit

       end if

       if (i.eq.180) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 006"
          close (11)
       end if

    end do

    call random(thread_id,xi)
    if (xi.gt.0.5_dp) beta = beta + pi

    if (beta.ge.2._dp*pi) beta = 2._dp*pi - 1.e-10_dp
    if (beta.le.0._dp) beta = -2._dp*pi + 1.e-10_dp
    
    ! Scattering angle sampling

    c2b = cos(2._dp*beta)
    s2b = sqrt(1._dp-c2b*c2b)

    if (beta.gt.pi/2._dp.and.beta.lt.pi) then
       s2b = -s2b
    else if (beta.gt.3._dp*pi/2._dp.and.beta.lt.2._dp*pi) then
       s2b = -s2b
    else if (beta.gt.-pi/2._dp.and.beta.lt.0._dp) then
       s2b = -s2b
    else if (beta.gt.-2._dp*pi.and.beta.lt.-3._dp*pi/2._dp) then
       s2b = -s2b
    end if
    
    do i=1,180

       intensity(i) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,1,i)*stokes(1) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,2,i)*c2b*stokes(2) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,2,i)*s2b*stokes(3) - &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,3,i)*s2b*stokes(2) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,3,i)*c2b*stokes(3) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,4,i)*stokes(4)

       intensity(i) = intensity(i) * sinbeta(i) * pi/180._dp

       intensity_cumulative(i) = intensity_cumulative(i-1) + intensity(i)

    end do

    call random(thread_id,xi)
    intensity_cumulative_sampled = xi*intensity_cumulative(180)

    do i=1,180

       if ( intensity_cumulative_sampled.ge.intensity_cumulative(i-1) .and. &
            intensity_cumulative_sampled.le.intensity_cumulative(i) ) then

          x0 = dble(i-1)
          x1 = dble(i)
          y0 = intensity_cumulative(i-1)
          y1 = intensity_cumulative(i)
          alpha = (x1-x0) * (intensity_cumulative_sampled-y0) / (y1-y0) + x0
          alpha = cos(alpha*pi/180._dp)

          if (abs(alpha).ge.1._dp) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 056"
             close (11)
          end if

          exit

          if (i.eq.180) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 007"
             close (11)
          end if

       end if

    end do

    if (alpha.ge.1._dp) alpha = 1._dp - 1.e-10_dp
    if (alpha.le.-1._dp) alpha = -1._dp + 1.e-10_dp

  end subroutine scattering_angle_sampling

  subroutine polarization_rotation(alpha, beta, stokes_meridian_in, scatter, direction, direction_new, &
       stokes_meridian_out, peeling)

    ! Calculate the new Stokes vector

    logical,  intent(in)  :: peeling
    real(dp), intent(in)  :: alpha, beta, stokes_meridian_in(4), scatter(4,4)
    real(dp), intent(in)  :: direction_new(3), direction(3)
    real(dp), intent(out) :: stokes_meridian_out(4)
    real(dp)              :: beta2, mueller_matrix(2,2), stokes_rotation(4), stokes_scattered(4), norm, num_check_2
    ! real(dp)              :: num_check_1, beta_check

    norm = 0._dp

    ! Check if sampled beta and beta_check are the same

    ! if (abs(alpha).lt.1._dp.and.abs(direction(3)).lt.1._dp) then

    !    num_check_1 = (direction_new(3) - direction(3)*alpha) / &
    !         ( sqrt(1._dp-alpha*alpha)*sqrt(1._dp-direction(3)*direction(3)) )

    !    if (abs(num_check_1).le.1._dp) then

    !       beta_check = acos( num_check_1 )

    !    else if (num_check_1.gt.1._dp.and.num_check_1.lt.1.0001_dp) then

    !       beta_check = 0._dp

    !    else if (num_check_1.lt.-1._dp.and.num_check_1.gt.-1.0001_dp) then

    !       beta_check = pi

    !    else

    !       open (11, file=trim(error_log), position="append")
    !       write (11,*) "error 008"
    !       write (11,*) num_check_1
    !       close (11)

    !    end if

    !    if (beta.ge.pi.and.beta.le.2._dp*pi) beta_check = 2._dp*pi - beta_check

    !    if (abs(beta-beta_check)*180._dp/pi.gt.1.e-2_dp) then

    !       open (11, file=trim(error_log), position="append")
    !       write (11,*) "error 009"
    !       write (11,*) abs(beta-beta_check)*180._dp/pi
    !       close (11)

    !    end if

    ! else

    !    open (11, file=trim(error_log), position="append")
    !    write (11,*) "error 010"
    !    write (11,*) alpha, beta
    !    write (11,*) direction
    !    close (11)
       
    ! end if

    ! Rotation from scattering plane to meridiane plane with the spherical cosine rule

    if (abs(alpha).lt.1._dp.and.abs(direction_new(3)).lt.1._dp) then

       num_check_2 = (direction(3) - direction_new(3)*alpha) / &
            ( sqrt(1._dp-alpha*alpha) * sqrt(1._dp-direction_new(3)*direction_new(3)) )

       if (abs(num_check_2).le.1._dp) then

          beta2 = acos( num_check_2 )

       else if (num_check_2.gt.1._dp.and.num_check_2.lt.1.00001_dp) then

          beta2 = 0._dp

       else if (num_check_2.lt.-1._dp.and.num_check_2.gt.-1.00001_dp) then

          beta2 = pi

       else

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 011"
          close (11)

       end if

       call mueller_matrix_filler(beta, mueller_matrix)

       ! Rotate the Stokes vector from the meridian plane to the scattering plane

       stokes_rotation(1) = stokes_meridian_in(1)
       stokes_rotation(2) = mueller_matrix(1,1)*stokes_meridian_in(2) + mueller_matrix(1,2)*stokes_meridian_in(3)
       stokes_rotation(3) = mueller_matrix(2,1)*stokes_meridian_in(2) + mueller_matrix(2,2)*stokes_meridian_in(3)
       stokes_rotation(4) = stokes_meridian_in(4)

       ! Re-normalized the polarized intensity to keep it constant before and after the Stokes rotation

       if (sqrt(stokes_rotation(2)**2+stokes_rotation(3)**2+stokes_rotation(4)**2).gt.0._dp) then

          norm = sqrt(stokes_meridian_in(2)**2+stokes_meridian_in(3)**2+stokes_meridian_in(4)**2) / &
               sqrt(stokes_rotation(2)**2+stokes_rotation(3)**2+stokes_rotation(4)**2)

       else

          norm = 1._dp

       end if

       if (norm.lt.1._dp.or.norm.gt.1._dp) then

          stokes_rotation(2) = stokes_rotation(2) * norm
          stokes_rotation(3) = stokes_rotation(3) * norm
          stokes_rotation(4) = stokes_rotation(4) * norm

       end if

       ! Apply scattering matrix

       stokes_scattered(1) = scatter(1,1)*stokes_rotation(1) + scatter(1,2)*stokes_rotation(2) + &
            scatter(1,3)*stokes_rotation(3) + scatter(1,4)*stokes_rotation(4)

       stokes_scattered(2) = scatter(2,1)*stokes_rotation(1) + scatter(2,2)*stokes_rotation(2) + &
            scatter(2,3)*stokes_rotation(3) + scatter(2,4)*stokes_rotation(4)

       stokes_scattered(3) = scatter(3,1)*stokes_rotation(1) + scatter(3,2)*stokes_rotation(2) + &
            scatter(3,3)*stokes_rotation(3) + scatter(3,4)*stokes_rotation(4)

       stokes_scattered(4) = scatter(4,1)*stokes_rotation(1) + scatter(4,2)*stokes_rotation(2) + &
            scatter(4,3)*stokes_rotation(3) + scatter(4,4)*stokes_rotation(4)

       ! Normalize the stokes vector such that intensity is conserved

       if (.not.peeling) then

          if (stokes_scattered(1).gt.0._dp) then

             norm = stokes_rotation(1) / stokes_scattered(1)
             stokes_scattered = norm * stokes_scattered

          else

             open (11, file=trim(error_log), position="append")
             write (11,*) "error 012"
             close (11)

          end if

       end if

       ! Obtain Mueller matrix

       if (beta.ge.0._dp.and.beta.lt.pi) then

          call mueller_matrix_filler(beta2, mueller_matrix)

       else if (beta.ge.pi.and.beta.lt.2._dp*pi) then

          call mueller_matrix_filler(-beta2, mueller_matrix)

       end if

       ! Rotate the new Stokes vector into the meridian plane

       stokes_meridian_out(1) = stokes_scattered(1)
       stokes_meridian_out(2) = mueller_matrix(1,1)*stokes_scattered(2) + mueller_matrix(1,2)*stokes_scattered(3)
       stokes_meridian_out(3) = mueller_matrix(2,1)*stokes_scattered(2) + mueller_matrix(2,2)*stokes_scattered(3)
       stokes_meridian_out(4) = stokes_scattered(4)

       ! Re-normalized the polarized intensity to keep it constant before and after the Stokes rotation

       if (sqrt(stokes_meridian_out(2)**2+stokes_meridian_out(3)**2+stokes_meridian_out(4)**2).gt.0._dp) then

          norm = sqrt(stokes_scattered(2)**2+stokes_scattered(3)**2+stokes_scattered(4)**2) / &
               sqrt(stokes_meridian_out(2)**2+stokes_meridian_out(3)**2+stokes_meridian_out(4)**2)

       else

          norm = 1._dp

       end if

       if (norm.lt.1._dp.or.norm.gt.1._dp) then

          stokes_meridian_out(2) = stokes_meridian_out(2) * norm
          stokes_meridian_out(3) = stokes_meridian_out(3) * norm
          stokes_meridian_out(4) = stokes_meridian_out(4) * norm

       end if

    else if (alpha.ge.1._dp.and.alpha.lt.1.0001_dp) then

       ! Forward scattering

       stokes_meridian_out = stokes_meridian_in

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 013"
       write (11,*) alpha, beta
       close (11)

    else if (alpha.le.-1._dp.and.alpha.gt.-1.0001_dp) then

       ! Backward scattering

       stokes_scattered(1) = scatter(1,1)*stokes_meridian_in(1) + scatter(1,2)*stokes_meridian_in(2) + &
            scatter(1,3)*stokes_meridian_in(3) + scatter(1,4)*stokes_meridian_in(4)

       stokes_scattered(2) = scatter(2,1)*stokes_meridian_in(1) + scatter(2,2)*stokes_meridian_in(2) + &
            scatter(2,3)*stokes_meridian_in(3) + scatter(2,4)*stokes_meridian_in(4)

       stokes_scattered(3) = scatter(3,1)*stokes_meridian_in(1) + scatter(3,2)*stokes_meridian_in(2) + &
            scatter(3,3)*stokes_meridian_in(3) + scatter(3,4)*stokes_meridian_in(4)

       stokes_scattered(4) = scatter(4,1)*stokes_meridian_in(1) + scatter(4,2)*stokes_meridian_in(2) + &
            scatter(4,3)*stokes_meridian_in(3) + scatter(4,4)*stokes_meridian_in(4)

       if (peeling) then

          stokes_meridian_out = stokes_scattered

       else if (.not.peeling) then

          ! Normalize the stokes vector such that intensity is conserved

          if (stokes_scattered(1).gt.0._dp) then

             norm = stokes_meridian_in(1) / stokes_scattered(1)
             stokes_meridian_out = norm * stokes_scattered

          else

             stokes_meridian_out = 0._dp

             open (11, file=trim(error_log), position="append")
             write (11,*) "error 014"
             close (11)

          end if

       end if

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 015"
       write (11,*) alpha, beta
       write (11,*) direction
       close (11)

    else

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 016"
       close (11)

    end if

    ! if (sqrt(stokes_meridian_out(2)*stokes_meridian_out(2)+ &
    !      stokes_meridian_out(3)*stokes_meridian_out(3)).gt.stokes_meridian_out(1)) then

    !    open (11, file=trim(error_log), position="append")
    !    write (11,*) "error 017"
    !    write (11,*) peeling, alpha, beta
    !    close (11)

    ! end if

  end subroutine polarization_rotation

  subroutine mueller_matrix_filler(psi, mueller)

    ! Mueller matrix (row, column)

    real(dp), intent(in)  :: psi
    real(dp), intent(out) :: mueller(2,2)
    real(dp)              :: c2p, s2p

    c2p = cos(2._dp*psi)
    s2p = sqrt(1._dp-c2p*c2p)

    if (psi.gt.pi/2._dp.and.psi.lt.pi) then
       s2p = -s2p
    else if (psi.gt.3._dp*pi/2._dp.and.psi.lt.2._dp*pi) then
       s2p = -s2p
    else if (psi.gt.-pi/2._dp.and.psi.lt.0._dp) then
       s2p = -s2p
    else if (psi.gt.-2._dp*pi.and.psi.lt.-3._dp*pi/2._dp) then
       s2p = -s2p
    end if

    mueller(1,1) = c2p
    mueller(2,1) = -s2p
    mueller(1,2) = s2p
    mueller(2,2) = c2p

  end subroutine mueller_matrix_filler

  subroutine direction_cosine(alpha, beta, direction, direction_new)

    real(dp), intent(in)  :: alpha, beta, direction(3)
    real(dp), intent(out) :: direction_new(3)
    real(dp)              :: cos_theta_old, sin_theta_old, cos_theta_new, sin_theta_new
    real(dp)              :: cos_phi_new, sin_phi_new, phi_old, phi_new, num_check_2

    ! Initital photon direction
    
    cos_theta_old = direction(3)/sqrt(direction(1)*direction(1)+direction(2)*direction(2)+direction(3)*direction(3))
    sin_theta_old = sqrt(1._dp-cos_theta_old*cos_theta_old)

    phi_old = atan2(direction(2),direction(1))
    if (phi_old.lt.0._dp) phi_old = phi_old+2._dp*pi
    
    ! New theta direction from spherical cosine rule

    if (beta.ge.pi.and.beta.lt.2._dp*pi) then

       cos_theta_new = cos_theta_old*alpha + sin_theta_old*sqrt(1._dp-alpha*alpha)*cos(2._dp*pi-beta)

    else if (beta.ge.0._dp.and.beta.lt.pi) then

       cos_theta_new = cos_theta_old*alpha + sin_theta_old*sqrt(1._dp-alpha*alpha)*cos(beta)

    else

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 018"
       close (11)

    end if

    sin_theta_new = sqrt(1._dp-cos_theta_new*cos_theta_new)

    ! New phi direction from spherical cosine rule

    num_check_2 = (alpha-cos_theta_new*cos_theta_old) / (sin_theta_new*sin_theta_old)

    if (num_check_2.ge.1._dp) then
       num_check_2 = 1._dp - 1.e-10_dp
    else if (num_check_2.le.-1._dp) then
       num_check_2 = -1._dp + 1.e-10_dp
    end if
    
    if (abs(num_check_2).le.1._dp) then

       if (beta.ge.pi.and.beta.lt.2._dp*pi) then
          phi_new = phi_old - acos( num_check_2  )
       else if (beta.ge.0._dp.and.beta.lt.pi) then
          phi_new = phi_old + acos( num_check_2 )
       else
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 019"
          close (11)
       end if

    else

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 020"
       write (11,*) alpha, beta, num_check_2
       write (11,*) direction
       close (11)

    end if

    ! Shift new phi direction by 2*pi if needed

    if (phi_new.lt.0._dp) phi_new = phi_new + 2._dp*pi
    if (phi_new.gt.2._dp*pi) phi_new = phi_new - 2._dp*pi

    ! New Cartesian photon direction vector

    cos_phi_new = cos(phi_new)
    
    if (phi_new.ge.0._dp.and.phi_new.lt.pi) then
       sin_phi_new = sqrt(1._dp-cos_phi_new*cos_phi_new)
    else if (phi_new.ge.pi.and.phi_new.le.2._dp*pi) then
       sin_phi_new = -sqrt(1._dp-cos_phi_new*cos_phi_new)
    else
       open (11, file=trim(error_log), position="append")
       write (11,*) "error 021"
       close (11)
    end if

    direction_new(1) = sin_theta_new*cos_phi_new
    direction_new(2) = sin_theta_new*sin_phi_new
    direction_new(3) = cos_theta_new
    
  end subroutine direction_cosine

  subroutine get_atmosphere

    integer               :: unit, status, readwrite, hdutype, nfound, blocksize, i, j, k, m, n
    integer, allocatable  :: naxes(:)
    real(dp), allocatable :: temp(:)
    logical               :: anynul
    character(100)        :: fitsfile

    status = 0
    readwrite = 0

    fitsfile  = trim(atmosphere_directory)//"/atmosphere.fits"

    call ftgiou(unit, status)
    call ftopen(unit, fitsfile, readwrite, blocksize, status)
    call ftmahd(unit, 1, hdutype, status)

    ! Radial grid [m]

    call ftgknj(unit,'NAXIS',1,1,nr,nfound,status)
    allocate (temp(nr))
    temp = 0._dp
    call ftgpvd(unit,1,1,nr,-999._dp,temp,anynul,status)
    nr = nr-1
    allocate (rfront(0:nr))
    rfront = 0._dp
    do i = 0,nr
       rfront(i) = temp(i+1)
    end do
    deallocate (temp)

    ! Theta grid [rad]

    call ftmrhd(unit,1,hdutype,status)    
    call ftgknj(unit,'NAXIS',1,1,ntheta,nfound,status)
    allocate (temp(ntheta))
    temp = 0._dp
    call ftgpvd(unit,1,1,ntheta,-999._dp,temp,anynul,status)
    ntheta = ntheta-1
    allocate (thetafront(0:ntheta))
    thetafront = 0._dp
    allocate (thetaplane(0:ntheta))
    thetaplane = 0._dp
    do i = 0,ntheta
       thetafront(i) = temp(i+1)
       if (thetafront(i).lt.90._dp-1.e-6_dp.or.thetafront(i).gt.90._dp+1.e-6_dp) then
          thetaplane(i) = 1
       else
          thetaplane(i) = 2
       end if
    end do
    deallocate (temp)
    thetafront = thetafront*pi/180.

    ! Phi grid [rad]

    call ftmrhd(unit,1,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,nphi,nfound,status)
    allocate (temp(nphi))
    temp = 0._dp
    call ftgpvd(unit,1,1,nphi,-999._dp,temp,anynul,status)
    nphi = nphi
    allocate (phifront(0:nphi-1))
    phifront = 0._dp
    do i = 0,nphi-1
       phifront(i) = temp(i+1)
    end do
    deallocate (temp)
    phifront = phifront*pi/180.

    ! Wavelengths [micron]

    call ftmrhd(unit,1,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,n_wavelength,nfound,status)
    allocate (wavelengths(n_wavelength))
    wavelengths = 0._dp
    call ftgpvd(unit,1,1,n_wavelength,-999._dp,wavelengths,anynul,status)
    wavelengths = wavelengths * 1.e-6_dp ! [micron] to [m]

    ! Density [kg m-3]

    allocate (cell_density(0:nr-1,0:ntheta-1,0:nphi-1))
    cell_density = 0._dp
    allocate (naxes(3))
    naxes = 0
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3), -999._dp, cell_density, anynul, status)
    deallocate (cell_density)
    
    ! Temperature [K]

    allocate (cell_temperature(0:nr-1,0:ntheta-1,0:nphi-1))
    cell_temperature = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3), -999._dp, cell_temperature, anynul, status)
    deallocate (naxes)
    
    ! Scattering opacity [m-1]

    allocate (naxes(4))
    naxes = 0
    allocate (cell_scattering_opacity(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_scattering_opacity = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4), -999._dp, cell_scattering_opacity, anynul, status)

    ! Absorption opacity [m-1]

    allocate (cell_absorption_opacity(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_absorption_opacity = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4), -999._dp, cell_absorption_opacity, anynul, status)
    deallocate (naxes)

    ! Opacity [m-1]

    allocate (cell_opacity(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_albedo(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_opacity = 0._dp
    cell_albedo = 0._dp
    do i=0,nr-1
       do j=0,ntheta-1
          do k=0,nphi-1
             do m=1,n_wavelength
                cell_opacity(i,j,k,m) = cell_scattering_opacity(i,j,k,m) + cell_absorption_opacity(i,j,k,m)
                if (cell_opacity(i,j,k,m).gt.0._dp) cell_albedo(i,j,k,m) = cell_scattering_opacity(i,j,k,m) / cell_opacity(i,j,k,m)
                if (cell_albedo(i,j,k,m).lt.1.e-20_dp) cell_albedo(i,j,k,m) = 1.e-20_dp
             end do
          end do
       end do
    end do
    
    ! Scattering matrix
    
    allocate (naxes(6))
    naxes = 0
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,6,naxes,nfound,status)
    allocate (cell_scatter_matrix(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength,16,180))
    cell_scatter_matrix = 0._dp
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6), -999._dp, cell_scatter_matrix, anynul, status)

    call ftclos(unit, status)
    call ftfiou(unit, status)

    ! Integral of first row elements
    
    allocate (cell_p11_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_p12_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_p13_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_p14_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))

    cell_p11_int = 0._dp
    cell_p12_int = 0._dp
    cell_p13_int = 0._dp
    cell_p14_int = 0._dp

    do i=0,nr-1
       do j=0,ntheta-1
          do k=0,nphi-1
             do m=1,n_wavelength
                do n=1,180

                   cell_p11_int(i,j,k,m) = cell_p11_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,1,n)*sinbeta(n)*pi/180._dp
                   cell_p12_int(i,j,k,m) = cell_p12_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,2,n)*sinbeta(n)*pi/180._dp
                   cell_p13_int(i,j,k,m) = cell_p13_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,3,n)*sinbeta(n)*pi/180._dp
                   cell_p14_int(i,j,k,m) = cell_p14_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,4,n)*sinbeta(n)*pi/180._dp

                end do
             end do
          end do
       end do
    end do

    ! Total number of cells
    cells = nr*ntheta*nphi

  end subroutine get_atmosphere

  subroutine grid_initialize(mode)

    ! Atmospheric grid initialization

    integer, intent(in) :: mode
    integer             :: i, j, k, cell_max, grid_out
    real(dp)            :: planck_flux, optical_depth_total1, optical_depth_total2, optical_depth_total
    real(dp)            :: weight_norm, emissivity_total
    logical             :: exist

    if (mode.eq.1) then

       ! Wavelength independent arrays

       allocate (theta_grid_cos(0:ntheta))
       theta_grid_cos = 0._dp
       allocate (theta_grid_tan(0:ntheta))
       theta_grid_tan = 0._dp
       allocate (phi_grid_cos(0:nphi-1))
       phi_grid_cos = 0._dp
       allocate (phi_grid_sin(0:nphi-1))
       phi_grid_sin = 0._dp

       ! Set sine and cosine values of theta faces
       do i=0,ntheta
          theta_grid_cos(i) = cos(thetafront(i))
          theta_grid_tan(i) = tan(thetafront(i))
       end do

       ! Set sine and cosine values of phi faces
       do i=0,nphi-1
          phi_grid_cos(i) = cos(phifront(i))
          phi_grid_sin(i) = sin(phifront(i))
       end do

       ! Cell volume

       allocate (cell_volume(0:nr-1,0:ntheta-1,0:nphi-1))
       cell_volume = 0._dp

       do k=0,nphi-1
          do j=0,ntheta-1
             do i=0,nr-1

                ! Cell volume
                ! Volume = a * b * c * (1/3) * (r_out^3-r_in^3) * (cos(theta_in)-cos(theta_out)) * (phi_out-phi_in)

                if (nphi.eq.1) then
                
                   cell_volume(i,j,k) = oblate_x * oblate_y * oblate_z * (1._dp/3._dp) * &
                        (rfront(i+1)**3-rfront(i)**3) * (theta_grid_cos(j)-theta_grid_cos(j+1)) * 2._dp*pi

                else if (nphi.gt.1) then

                   if (k.lt.nphi-1) then

                      cell_volume(i,j,k) = oblate_x * oblate_y * oblate_z * (1._dp/3._dp) * &
                           (rfront(i+1)**3-rfront(i)**3) * (theta_grid_cos(j)-theta_grid_cos(j+1)) * (phifront(k+1)-phifront(k))

                   else if (k.eq.nphi-1) then

                      cell_volume(i,j,k) = oblate_x * oblate_y * oblate_z * (1._dp/3._dp) * &
                           (rfront(i+1)**3-rfront(i)**3) * (theta_grid_cos(j)-theta_grid_cos(j+1)) * (2._dp*pi-phifront(k))

                   end if
                      
                end if

             end do
          end do
       end do

       ! Cell energy transport
       
       if (flow_global) then
       
          allocate (cell_flow_global(threads,3,0:nr-1,0:ntheta-1,0:nphi-1))
          cell_flow_global = 0._dp

       end if

       if (flow_theta) then

          allocate (cell_flow(threads,4,0:nr-1,0:ntheta-1,0:nphi-1))
          cell_flow = 0._dp

       end if

    else if (mode.eq.2) then

       ! Wavelength dependent arrays

       if (photon_source.eq.1) then

          ! Check at which radial depth the optical thickness is 30 or larger
          ! Determine at which theta/phi location this is deepest in the atmosphere

          cell_max = 1000000

          do j=0,ntheta-1
             do k=0,nphi-1

                optical_depth_total = 0._dp

                do i=0,nr-1

                   optical_depth_total = optical_depth_total + cell_opacity(nr-i-1,j,k,wl_count) * &
                        ( rfront(nr-i)-rfront(nr-i-1) )

                   cell_depth = nr-i-1

                   if (optical_depth_total.gt.30._dp) exit

                end do

                if (cell_depth.lt.cell_max) cell_max = cell_depth

             end do
          end do

          cell_depth = cell_max

       else if (photon_source.eq.2) then

          ! Check at which radial depth the absorption optical thickness is 5 or larger
          ! Determine at which theta/phi location this is deepest in the atmosphere

          cell_max = 1000000

          if (ring) then
             grid_out = 2
          else if (.not.ring) then
             grid_out = 0
          end if
          
          do j=0,ntheta-1
             do k=0,nphi-1

                optical_depth_total = 0._dp

                do i=grid_out,nr-1

                   optical_depth_total = optical_depth_total + cell_absorption_opacity(nr-i-1,j,k,wl_count) * &
                        ( rfront(nr-i)-rfront(nr-i-1) )

                   cell_depth = nr-i-1

                   if (optical_depth_total.gt.5._dp) exit

                end do

                if (cell_depth.lt.cell_max) cell_max = cell_depth

             end do
          end do

          cell_depth = cell_max

          ! Thermal cell luminosity
          
          weight_norm = 0._dp

          do i=cell_depth,nr-1
             do j=0,ntheta-1
                do k=0,nphi-1

                   if (cell_temperature(i,j,k).gt.0._dp) then
                      
                      call planck_function(cell_temperature(i,j,k), planck_flux)
                      weight_norm = weight_norm + cell_absorption_opacity(i,j,k,wl_count)*planck_flux*cell_volume(i,j,k)

                   end if

                end do
             end do
          end do

          allocate (cell_weight(0:nr-1,0:ntheta-1,0:nphi-1))
          cell_weight = 0._dp

          allocate (cell_luminosity(0:nr-1,0:ntheta-1,0:nphi-1))
          cell_luminosity = 0._dp

          allocate (emissivity_cumulative(0:nr-1,0:ntheta-1,0:nphi-1))
          emissivity_cumulative = 0._dp

          emissivity_total = 0._dp
          
          do i=cell_depth,nr-1
             do j=0,ntheta-1
                do k=0,nphi-1

                   if (cell_temperature(i,j,k).gt.0._dp.and.cell_absorption_opacity(i,j,k,wl_count).gt.0._dp) then
                      
                      call planck_function(cell_temperature(i,j,k), planck_flux)

                      if (thermal_weight) then
                         cell_weight(i,j,k) = weight_norm / (cell_volume(i,j,k)*cell_absorption_opacity(i,j,k,wl_count)*planck_flux)
                      else if (.not.thermal_weight) then
                         cell_weight(i,j,k) = 1._dp
                      end if
                      
                      cell_luminosity(i,j,k) = 4._dp*pi*cell_volume(i,j,k)*cell_absorption_opacity(i,j,k,wl_count)*planck_flux ! [W m-1]

                      emissivity_cumulative(i,j,k) = emissivity_total + cell_luminosity(i,j,k) * cell_weight(i,j,k)

                      emissivity_total = emissivity_cumulative(i,j,k)

                   else

                      emissivity_cumulative(i,j,k) = emissivity_total

                   end if

                end do
             end do
          end do

       end if

       ! Optical depth output

       if (imaging_broad.or.spectrum) then

          inquire(file='output/'//adjustl(trim(output_name))//'/output/optical_depth.dat', exist=exist)

          if (exist) then

             open (100, file='output/'//adjustl(trim(output_name))//'/output/optical_depth.dat', &
                  status='old', position='append', action='write')

          else if (.not.exist) then

             open (100, file='output/'//adjustl(trim(output_name))//'/output/optical_depth.dat', &
                  status='new', action='write')

             write (100,*) "# Wavelength [micron] - Total optical depth - Absorption optical depth - Scattering optical depth"
             write (100,*)

          end if

          optical_depth_total = 0._dp
          optical_depth_total1 = 0._dp
          optical_depth_total2 = 0._dp

          do i=0,nr-1
             optical_depth_total = optical_depth_total + ( (rfront(i+1)-rfront(i)) * cell_opacity(i,0,0,wl_count) )
             optical_depth_total1 = optical_depth_total1 + ( (rfront(i+1)-rfront(i)) * cell_scattering_opacity(i,0,0,wl_count) )
             optical_depth_total2 = optical_depth_total2 + ( (rfront(i+1)-rfront(i)) * cell_absorption_opacity(i,0,0,wl_count) )
          end do

          ! Wavelength [micron] - Total radial optical depth total - absorption - scattering
          write (100,*) wavelength*1.e6_dp, optical_depth_total, optical_depth_total2, optical_depth_total1

          close (100)

       end if

       if (photon_source.eq.2) then
       
          allocate (flux_emitted(threads))
          flux_emitted = 0._dp
          
          allocate (flux_exit(threads))
          flux_exit = 0._dp

       end if
          
    end if

  end subroutine grid_initialize

  subroutine photon_package

    ! Photon package energy

    real(dp) :: planck_flux

    if (photon_source.eq.1) then
    
       ! L_star = 4 * pi * R_star^2 * F_lambda = 4 * pi * D^2 * F_lambda,planet
       ! F_lambda,planet = F_lambda * R_star^2/D^2
       ! F_lambda,planet * pi * R_planet^2 = d^2 * F_obs
       ! F_obs = pi * R_planet^2 * R_star^2 * F_lambda / ( D^2 * d^2 )
    
       call planck_function(t_star, planck_flux)

       package_energy = pi * planck_flux * rfront(nr) * rfront(nr) * r_star * r_star / &
            ( orbit * orbit * distance_planet * distance_planet * dble(packages) )

       if (phase_curve.and.det_phi*180._dp/pi.ge.170._dp) then
          
          package_energy = package_energy * (pi*r_star*r_star-0.9_dp*0.9_dp*pi*r_star*r_star)/(pi*r_star*r_star)

       end if

    else if (photon_source.eq.2) then
       
       package_energy = emissivity_cumulative(nr-1,ntheta-1,nphi-1) / ( distance_planet * distance_planet * dble(packages) )

    end if

  end subroutine photon_package

  subroutine array_start

    allocate (detector(nx,ny,4,3))
    detector = 0._dp
    allocate (detector_thread(nx,ny,4,3,threads))
    detector_thread = 0._dp

  end subroutine array_start

  subroutine grid_finished(i)

    integer, intent(in) :: i

    if (i.eq.1) then

       ! Wavelength dependent arrays

       if (photon_source.eq.2) then
       
          deallocate (cell_luminosity)
          deallocate (cell_weight)
          deallocate (emissivity_cumulative)
          deallocate (flux_emitted)
          deallocate (flux_exit)

       end if
       
    else if (i.eq.2) then

       ! Detector arrays

       deallocate (detector)
       deallocate (detector_thread)

    else if (i.eq.3) then

       ! Global grid arrays

       deallocate (cell_absorption_opacity)
       deallocate (cell_scattering_opacity)
       deallocate (cell_opacity)
       deallocate (cell_albedo)
       deallocate (cell_scatter_matrix)
       deallocate (cell_p11_int)
       deallocate (cell_p12_int)
       deallocate (cell_p13_int)
       deallocate (cell_p14_int)
       deallocate (rfront)
       deallocate (thetafront)
       deallocate (thetaplane)
       deallocate (phifront)
       deallocate (theta_grid_cos)
       deallocate (theta_grid_tan)
       deallocate (phi_grid_cos)
       deallocate (phi_grid_sin)
       deallocate (cell_volume)
       deallocate (wavelengths)
       if (flow_global) deallocate (cell_flow_global)
       if (flow_theta) deallocate (cell_flow)
       
    end if

  end subroutine grid_finished

  subroutine initial_cell(x, y, z, cell)

    ! Initial grid cell that is crossed

    integer               :: j
    real(dp), intent(in)  :: x, y, z
    real(dp)              :: r, theta, phi
    integer,  intent(out) :: cell(3)

    cell = 0

    call cartesian_spherical(x, y, z, r, theta, phi)

    ! Cell(1)

    if (photon_source.eq.1) then

       cell(1) = nr-1

    else if (photon_source.eq.2) then

       cell(1) = 0

    end if

    ! Cell(2)

    do j=0,ntheta-1

       if (theta.gt.thetafront(j).and.theta.lt.thetafront(j+1)) then

          cell(2) = j
          exit

       end if

    end do

    ! Cell(3)

    do j=0,nphi-1

       if (j.lt.nphi-1) then

          if (phi.gt.phifront(j).and.phi.lt.phifront(j+1)) then

             cell(3) = j
             exit

          end if

       else if (j.eq.nphi-1) then

          if (phi.gt.phifront(j).and.phi.lt.2._dp*pi) then

             cell(3) = j
             exit

          end if

       end if

    end do

  end subroutine initial_cell

  subroutine next_cell(current_face, next_face, cell_in, cell_out)

    integer, intent(in)  :: current_face(2), next_face(2), cell_in(3)
    integer, intent(out) :: cell_out(3)

    cell_out = 0

    ! Next face is radial
    if (next_face(1).eq.1) then

       if (current_face(1).eq.1.and.next_face(2).eq.current_face(2)) then

          ! Next cell is radially outward, photon crosses the same radial surface
          cell_out(1) = cell_in(1)+1
          cell_out(2) = cell_in(2)
          cell_out(3) = cell_in(3)

       else if (next_face(2).eq.cell_in(1)) then

          ! Next cell is radially inward
          cell_out(1) = cell_in(1)-1
          cell_out(2) = cell_in(2)
          cell_out(3) = cell_in(3)

       else if (next_face(2).eq.cell_in(1)+1) then

          ! Next cell is radially outward
          cell_out(1) = cell_in(1)+1
          cell_out(2) = cell_in(2)
          cell_out(3) = cell_in(3)

       else

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 022"
          close (11)

       end if

    end if

    ! Next face is polar
    if (next_face(1).eq.2) then

       if (current_face(1).eq.2.and.next_face(2).eq.current_face(2).and.thetafront(next_face(2)).lt.pi/2._dp) then

          ! Next cell is theta outward, photon crosses the same theta surface
          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)+1
          cell_out(3) = cell_in(3)

       else if (current_face(1).eq.2.and.next_face(2).eq.current_face(2).and.thetafront(next_face(2)).gt.pi/2._dp) then

          ! Next cell is theta inward, photon crosses the same theta surface
          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)-1
          cell_out(3) = cell_in(3)

       else if (next_face(2).eq.cell_in(2)) then

          ! Next cell is theta inward
          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)-1
          cell_out(3) = cell_in(3)

       else if (next_face(2).eq.cell_in(2)+1) then

          ! Next cell is theta outward
          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)+1
          cell_out(3) = cell_in(3)

       else

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 023"
          close (11)

       end if

    end if

    ! Change in phi grid cell
    if (next_face(1).eq.3) then

       if (cell_in(3).eq.nphi-1.and.next_face(2).eq.0) then

          ! Next cell is azimuthally outward

          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)
          cell_out(3) = 0

       else if (cell_in(3).eq.0.and.next_face(2).eq.0) then

          ! Next cell is azimuthally inward

          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)
          cell_out(3) = nphi-1

       else if (next_face(2).eq.cell_in(3)+1) then

          ! Next cell is azimuthally outward

          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)
          cell_out(3) = cell_in(3)+1

       else if (next_face(2).eq.cell_in(3)) then

          ! Next cell is azimuthally inward

          cell_out(1) = cell_in(1)
          cell_out(2) = cell_in(2)
          cell_out(3) = cell_in(3)-1

       else

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 024"
          close (11)

       end if

    end if

  end subroutine next_cell

  subroutine cell_face(x, y, z, n, current_face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

    integer,  intent(in)  :: current_face(2), cell(3)
    integer,  intent(out) :: next_face(2), cell_out(3)
    integer               :: face(3,3), face_loc(2), i, j, buffer(13), status
    real(dp), intent(in)  :: x, y, z, n(3)
    real(dp), intent(out) :: face_distance
    real(dp)              :: a, b, c, qa, qb, qc, distance(3,3), z_test, tan_theta
    real(dp)              :: solutions_r(2), solutions_t(2), solutions_p(2), r, theta, phi
    logical,  intent(out) :: grid_exit, cell_error

    if (cell(1).lt.cell_depth) then
       open (11, file=trim(error_log), position="append")
       write (11,*) "error 025"
       write (11,*) cell_depth
       write (11,*) current_face
       close (11)
    end if

    ! Check if error file is not becoming larger than 100 MB
    call stat('output/'//adjustl(trim(output_name))//'/error.log', buffer, status)
    if (buffer(8).gt.1.e8_dp) then
       open (11, file=trim(error_log), position="append")
       write (11,*) "error 026"
       close (11)
       call exit(0)
    end if

    solutions_r = 0._dp
    solutions_t = 0._dp
    solutions_p = 0._dp
    
    grid_exit = .false.
    cell_error = .false.

    distance = 0._dp
    face_loc = 0

    a = 1._dp/oblate_x
    b = 1._dp/oblate_y
    c = 1._dp/oblate_z

    face(1,1) = cell(1)
    face(1,2) = cell(1) + 1
    face(1,3) = -999

    face(2,1) = cell(2)
    face(2,2) = cell(2) + 1
    face(2,3) = -999

    face(3,1) = cell(3)
    face(3,2) = cell(3) + 1
    face(3,3) = -999
    if (face(3,2).eq.nphi) face(3,2) = 0

    if (current_face(1).eq.1) then

       face(1,1) = current_face(2)-1
       face(1,2) = current_face(2)+1
       face(1,3) = current_face(2)

    else if (current_face(1).eq.2) then

       face(2,1) = current_face(2)-1
       face(2,2) = current_face(2)+1
       face(2,3) = current_face(2)

    else if (current_face(1).eq.3) then

       if (current_face(2).eq.0) then
          face(3,1) = nphi-1
       else
          face(3,1) = current_face(2)-1
       end if
       
       if (current_face(2).eq.nphi-1) then
          face(3,2) = 0
       else
          face(3,2) = current_face(2)+1
       end if

    end if

    ! Current cell face is a radial face
    
    if (current_face(1).eq.1) then

       if (cell(1).eq.current_face(2)-1) then

          ! Radial inward face, photon comes from radial outward cell

          qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) + c*c*n(3)*n(3)
          qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) + c*c*z*n(3))
          qc = a*a*x*x + b*b*y*y + c*c*z*z - rfront(face(1,1))*rfront(face(1,1))

          call quadratic_equation(qa, qb, qc, solutions_r)

          if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).le.1.e-15_dp.and.solutions_r(1).lt.1.e100_dp) then
             distance(1,1) = solutions_r(1)
          else if (solutions_r(2).gt.1.e-15_dp.and.solutions_r(1).le.1.e-15_dp.and.solutions_r(2).lt.1.e100_dp) then
             distance(1,1) = solutions_r(2)
          else if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).gt.1.e-15_dp) then
             if (solutions_r(1).lt.1.e100_dp.and.solutions_r(1).lt.solutions_r(2)) then
                distance(1,1) = solutions_r(1)
             else if (solutions_r(2).lt.1.e100_dp.and.solutions_r(2).lt.solutions_r(1)) then
                distance(1,1) = solutions_r(2)
             end if
          end if

       else if (cell(1).eq.current_face(2)) then

          ! Radial outward face, photon comes from radial inward cell

          qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) + c*c*n(3)*n(3)
          qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) + c*c*z*n(3))
          qc = a*a*x*x + b*b*y*y + c*c*z*z - rfront(face(1,2))*rfront(face(1,2))

          call quadratic_equation(qa, qb, qc, solutions_r)

          if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).le.1.e-15_dp.and.solutions_r(1).lt.1.e100_dp) then
             distance(1,2) = solutions_r(1)
          else if (solutions_r(2).gt.1.e-15_dp.and.solutions_r(1).le.1.e-15_dp.and.solutions_r(2).lt.1.e100_dp) then
             distance(1,2) = solutions_r(2)
          else if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).gt.1.e-15_dp) then
             if (solutions_r(1).lt.1.e100_dp.and.solutions_r(1).lt.solutions_r(2)) then
                distance(1,2) = solutions_r(1)
             else if (solutions_r(2).lt.1.e100_dp.and.solutions_r(2).lt.solutions_r(1)) then
                distance(1,2) = solutions_r(2)
             end if
          end if

       end if

       if (cell(1).eq.current_face(2)-1) then

          ! Radial same face
          ! Only possible for radially inward moving photon

          qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) + c*c*n(3)*n(3)
          qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) + c*c*z*n(3))
          qc = a*a*x*x + b*b*y*y + c*c*z*z - rfront(face(1,3))*rfront(face(1,3))

          call quadratic_equation(qa, qb, qc, solutions_r)

          if (solutions_r(1).gt.1.e-3_dp.and.solutions_r(2).le.1.e-3_dp.and.solutions_r(1).lt.1.e100_dp) then
             distance(1,3) = solutions_r(1)
          else if (solutions_r(2).gt.1.e-3_dp.and.solutions_r(1).le.1.e-3_dp.and.solutions_r(2).lt.1.e100_dp) then
             distance(1,3) = solutions_r(2)
          else if (solutions_r(1).gt.1.e-3_dp.and.solutions_r(2).gt.1.e-3_dp) then
             if (solutions_r(1).lt.1.e100_dp.and.solutions_r(1).lt.solutions_r(2)) then
                distance(1,3) = solutions_r(1)
             else if (solutions_r(2).lt.1.e100_dp.and.solutions_r(2).lt.solutions_r(1)) then
                distance(1,3) = solutions_r(2)
             end if
          end if

       else if (cell(1).eq.current_face(2)) then

       else

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 027"
          close (11)

       end if

    else if (current_face(1).ne.1) then

       ! Cell face is not a radial cell face

       ! Radial inward face

       qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) + c*c*n(3)*n(3)
       qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) + c*c*z*n(3))
       qc = a*a*x*x + b*b*y*y + c*c*z*z - rfront(face(1,1))*rfront(face(1,1))

       call quadratic_equation(qa, qb, qc, solutions_r)

       if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).le.1.e-15_dp.and.solutions_r(1).lt.1.e100_dp) then
          distance(1,1) = solutions_r(1)
       else if (solutions_r(2).gt.1.e-15_dp.and.solutions_r(1).le.1.e-15_dp.and.solutions_r(2).lt.1.e100_dp) then
          distance(1,1) = solutions_r(2)
       else if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).gt.1.e-15_dp) then
          if (solutions_r(1).lt.1.e100_dp.and.solutions_r(1).lt.solutions_r(2)) then
             distance(1,1) = solutions_r(1)
          else if (solutions_r(2).lt.1.e100_dp.and.solutions_r(2).lt.solutions_r(1)) then
             distance(1,1) = solutions_r(2)
          end if
       end if

       ! Radial outward face

       qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) + c*c*n(3)*n(3)
       qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) + c*c*z*n(3))
       qc = a*a*x*x + b*b*y*y + c*c*z*z - rfront(face(1,2))*rfront(face(1,2))

       call quadratic_equation(qa, qb, qc, solutions_r)

       if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).le.1.e-15_dp.and.solutions_r(1).lt.1.e100_dp) then
          distance(1,2) = solutions_r(1)
       else if (solutions_r(2).gt.1.e-15_dp.and.solutions_r(1).le.1.e-15_dp.and.solutions_r(2).lt.1.e100_dp) then
          distance(1,2) = solutions_r(2)
       else if (solutions_r(1).gt.1.e-15_dp.and.solutions_r(2).gt.1.e-15_dp) then
          if (solutions_r(1).lt.1.e100_dp.and.solutions_r(1).lt.solutions_r(2)) then
             distance(1,2) = solutions_r(1)
          else if (solutions_r(2).lt.1.e100_dp.and.solutions_r(2).lt.solutions_r(1)) then
             distance(1,2) = solutions_r(2)
          end if
       end if

    end if

    ! Current cell face is a theta face

    if (current_face(1).eq.2.and.face(2,3).eq.-999) then
       open (11, file=trim(error_log), position="append")
       write (11,*) "error 028"
       close (11)
    end if

    if (current_face(1).eq.2) then

       if (cell(2).eq.current_face(2)-1.and.face(2,1).ne.0) then

          ! Theta inward face, photon coming from theta outward face

          if (thetaplane(face(2,1)).eq.1) then
          
             tan_theta = theta_grid_tan(face(2,1))

             qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) - c*c*n(3)*n(3)*tan_theta*tan_theta
             qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) - c*c*z*n(3)*tan_theta*tan_theta)
             qc = a*a*x*x + b*b*y*y - c*c*z*z*tan_theta*tan_theta

             call quadratic_equation(qa, qb, qc, solutions_t)

             if (solutions_t(1).gt.1.e-15_dp) then

                z_test = z+solutions_t(1)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,1)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,1)).lt.pi/2._dp) ) solutions_t(1) = 0._dp

             end if

             if (solutions_t(2).gt.1.e-15_dp) then

                z_test = z+solutions_t(2)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,1)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,1)).lt.pi/2._dp) ) solutions_t(2) = 0._dp

             end if

             if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).le.1.e-15_dp.and.solutions_t(1).lt.1.e100_dp) then
                distance(2,1) = solutions_t(1)
             else if (solutions_t(2).gt.1.e-15_dp.and.solutions_t(1).le.1.e-15_dp.and.solutions_t(2).lt.1.e100_dp) then
                distance(2,1) = solutions_t(2)
             else if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).gt.1.e-15_dp) then
                if (solutions_t(1).lt.1.e100_dp.and.solutions_t(1).lt.solutions_t(2)) then
                   distance(2,1) = solutions_t(1)
                else if (solutions_t(2).lt.1.e100_dp.and.solutions_t(2).lt.solutions_t(1)) then
                   distance(2,1) = solutions_t(2)
                end if
             end if

          else if (thetaplane(face(2,1)).eq.2) then

             if (-z/n(3).gt.0._dp.and.n(3).gt.1.e-15_dp) distance(2,1) = -z/n(3)
             
          end if
          
       else if (cell(2).eq.current_face(2).and.face(2,2).ne.ntheta) then

          ! Theta outward face, photon coming from inward face

          if (thetaplane(face(2,2)).eq.1) then

             tan_theta = theta_grid_tan(face(2,2))

             qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) - c*c*n(3)*n(3)*tan_theta*tan_theta
             qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) - c*c*z*n(3)*tan_theta*tan_theta)
             qc = a*a*x*x + b*b*y*y - c*c*z*z*tan_theta*tan_theta

             call quadratic_equation(qa, qb, qc, solutions_t)

             if (solutions_t(1).gt.1.e-15_dp) then

                z_test = z+solutions_t(1)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,2)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,2)).lt.pi/2._dp) ) solutions_t(1) = 0._dp

             end if

             if (solutions_t(2).gt.1.e-15_dp) then

                z_test = z+solutions_t(2)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,2)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,2)).lt.pi/2._dp) ) solutions_t(2) = 0._dp

             end if

             if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).le.1.e-15_dp.and.solutions_t(1).lt.1.e100_dp) then
                distance(2,2) = solutions_t(1)
             else if (solutions_t(2).gt.1.e-15_dp.and.solutions_t(1).le.1.e-15_dp.and.solutions_t(2).lt.1.e100_dp) then
                distance(2,2) = solutions_t(2)
             else if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).gt.1.e-15_dp) then
                if (solutions_t(1).lt.1.e100_dp.and.solutions_t(1).lt.solutions_t(2)) then
                   distance(2,2) = solutions_t(1)
                else if (solutions_t(2).lt.1.e100_dp.and.solutions_t(2).lt.solutions_t(1)) then
                   distance(2,2) = solutions_t(2)
                end if
             end if

          else if (thetaplane(face(2,2)).eq.2) then

             if (-z/n(3).gt.0._dp.and.n(3).lt.-1.e-15_dp) distance(2,2) = -z/n(3)
             
          end if

       end if

       if ( (thetafront(face(2,3)).lt.pi/2._dp.and.cell(2).eq.current_face(2)-1) .or. &
            (thetafront(face(2,3)).gt.pi/2._dp.and.cell(2).eq.current_face(2)) ) then

          ! Theta same face

          if (thetaplane(face(2,3)).eq.1) then

             tan_theta = theta_grid_tan(face(2,3))
          
             qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) - c*c*n(3)*n(3)*tan_theta*tan_theta
             qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) - c*c*z*n(3)*tan_theta*tan_theta)
             qc = a*a*x*x + b*b*y*y - c*c*z*z*tan_theta*tan_theta

             call quadratic_equation(qa, qb, qc, solutions_t)

             if (solutions_t(1).gt.1.e-15_dp) then

                z_test = z+solutions_t(1)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,3)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,3)).lt.pi/2._dp) ) solutions_t(1) = 0._dp

             end if

             if (solutions_t(2).gt.1.e-15_dp) then

                z_test = z+solutions_t(2)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,3)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,3)).lt.pi/2._dp) ) solutions_t(2) = 0._dp

             end if

             if (solutions_t(1).gt.1.e-3_dp.and.solutions_t(2).le.1.e-3_dp.and.solutions_t(1).lt.1.e100_dp) then
                distance(2,3) = solutions_t(1)
             else if (solutions_t(2).gt.1.e-3_dp.and.solutions_t(1).le.1.e-3_dp.and.solutions_t(2).lt.1.e100_dp) then
                distance(2,3) = solutions_t(2)
             else if (solutions_t(1).gt.1.e-3_dp.and.solutions_t(2).gt.1.e-3_dp) then
                if (solutions_t(1).lt.1.e100_dp.and.solutions_t(1).lt.solutions_t(2)) then
                   distance(2,3) = solutions_t(1)
                else if (solutions_t(2).lt.1.e100_dp.and.solutions_t(2).lt.solutions_t(1)) then
                   distance(2,3) = solutions_t(2)
                end if
             end if

          end if

       end if

    else if (current_face(1).ne.2) then

       ! Current face is not a theta face

       ! Check if the inward theta face exists

       if (face(2,1).lt.0.or.face(2,1).gt.ntheta) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 029"
          close (11)
          face(2,1) = 0
       end if

       if (face(2,1).ne.0) then

          ! Theta inward face

          if (thetaplane(face(2,1)).eq.1) then

             tan_theta = theta_grid_tan(face(2,1))

             qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) - c*c*n(3)*n(3)*tan_theta*tan_theta
             qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) - c*c*z*n(3)*tan_theta*tan_theta)
             qc = a*a*x*x + b*b*y*y - c*c*z*z*tan_theta*tan_theta

             call quadratic_equation(qa, qb, qc, solutions_t)

             if (solutions_t(1).gt.1.e-15_dp) then

                z_test = z+solutions_t(1)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,1)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,1)).lt.pi/2._dp) ) solutions_t(1) = 0._dp

             end if

             if (solutions_t(2).gt.1.e-15_dp) then

                z_test = z+solutions_t(2)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,1)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,1)).lt.pi/2._dp) ) solutions_t(2) = 0._dp

             end if

             if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).le.1.e-15_dp.and.solutions_t(1).lt.1.e100_dp) then
                distance(2,1) = solutions_t(1)
             else if (solutions_t(2).gt.1.e-15_dp.and.solutions_t(1).le.1.e-15_dp.and.solutions_t(2).lt.1.e100_dp) then
                distance(2,1) = solutions_t(2)
             else if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).gt.1.e-15_dp) then
                if (solutions_t(1).lt.1.e100_dp.and.solutions_t(1).lt.solutions_t(2)) then
                   distance(2,1) = solutions_t(1)
                else if (solutions_t(2).lt.1.e100_dp.and.solutions_t(2).lt.solutions_t(1)) then
                   distance(2,1) = solutions_t(2)
                end if
             end if

          else if (thetaplane(face(2,1)).eq.2) then
             
             if (-z/n(3).gt.0._dp.and.n(3).gt.1.e-15_dp) distance(2,1) = -z/n(3)

          end if

       end if

       if (face(2,2).ne.ntheta) then

          ! Theta outward face

          if (thetaplane(face(2,2)).eq.1) then

             tan_theta = theta_grid_tan(face(2,2))

             qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) - c*c*n(3)*n(3)*tan_theta*tan_theta
             qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) - c*c*z*n(3)*tan_theta*tan_theta)
             qc = a*a*x*x + b*b*y*y - c*c*z*z*tan_theta*tan_theta

             call quadratic_equation(qa, qb, qc, solutions_t)

             if (solutions_t(1).gt.1.e-15_dp) then

                z_test = z+solutions_t(1)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,2)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,2)).lt.pi/2._dp) ) solutions_t(1) = 0._dp

             end if

             if (solutions_t(2).gt.1.e-15_dp) then

                z_test = z+solutions_t(2)*n(3)

                if ( (z_test.gt.0._dp.and.thetafront(face(2,2)).gt.pi/2._dp) .or. &
                     (z_test.lt.0._dp.and.thetafront(face(2,2)).lt.pi/2._dp) ) solutions_t(2) = 0._dp

             end if

             if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).le.1.e-15_dp.and.solutions_t(1).lt.1.e100_dp) then
                distance(2,2) = solutions_t(1)
             else if (solutions_t(2).gt.1.e-15_dp.and.solutions_t(1).le.1.e-15_dp.and.solutions_t(2).lt.1.e100_dp) then
                distance(2,2) = solutions_t(2)
             else if (solutions_t(1).gt.1.e-15_dp.and.solutions_t(2).gt.1.e-15_dp) then
                if (solutions_t(1).lt.1.e100_dp.and.solutions_t(1).lt.solutions_t(2)) then
                   distance(2,2) = solutions_t(1)
                else if (solutions_t(2).lt.1.e100_dp.and.solutions_t(2).lt.solutions_t(1)) then
                   distance(2,2) = solutions_t(2)
                end if
             end if

          else if (thetaplane(face(2,2)).eq.2) then

             if (-z/n(3).gt.0._dp.and.n(3).lt.-1.e-15_dp) distance(2,2) = -z/n(3)

          end if

       end if

    end if

    if (current_face(1).eq.3) then

       ! Current face is a phi face

       if (cell(3).eq.current_face(2)-1.or.(cell(3).eq.nphi-1.and.current_face(2).eq.0)) then

          ! Photon is travelling inward in phi direction

          if (abs(b*n(2)*phi_grid_cos(face(3,1))-a*n(1)*phi_grid_sin(face(3,1))).gt.0._dp) then
          
             solutions_p(1) = (a*x*phi_grid_sin(face(3,1))-b*y*phi_grid_cos(face(3,1))) / &
                  (b*n(2)*phi_grid_cos(face(3,1))-a*n(1)*phi_grid_sin(face(3,1)))

             if (solutions_p(1).gt.1.e-15_dp.and.solutions_p(1).lt.1.e100_dp) distance(3,1) = solutions_p(1)
             
          end if

       else if (cell(3).eq.current_face(2)) then

          ! Photon is travelling outward in phi direction

          if (abs((b*n(2)*phi_grid_cos(face(3,2))-a*n(1)*phi_grid_sin(face(3,2)))).gt.0._dp) then

             solutions_p(2) = (a*x*phi_grid_sin(face(3,2))-b*y*phi_grid_cos(face(3,2))) / &
                  (b*n(2)*phi_grid_cos(face(3,2))-a*n(1)*phi_grid_sin(face(3,2)))

             if (solutions_p(2).gt.1.e-15_dp.and.solutions_p(1).lt.1.e100_dp) distance(3,2) = solutions_p(2)
             
          end if

       end if

    else if (current_face(1).ne.3.and.nphi.gt.1) then

       ! Current face in not a phi face
       
       ! Inward phi direction

       if (abs(b*n(2)*phi_grid_cos(face(3,1))-a*n(1)*phi_grid_sin(face(3,1))).gt.0._dp) then

          solutions_p(1) = (a*x*phi_grid_sin(face(3,1))-b*y*phi_grid_cos(face(3,1))) / &
               (b*n(2)*phi_grid_cos(face(3,1))-a*n(1)*phi_grid_sin(face(3,1)))

          if (solutions_p(1).gt.1.e-15_dp.and.solutions_p(1).lt.1.e100_dp) distance(3,1) = solutions_p(1)
          
       end if

       ! Outward phi direction

       if (abs(n(2)*phi_grid_cos(face(3,2))-n(1)*phi_grid_sin(face(3,2))).gt.0._dp) then

          solutions_p(2) = (a*x*phi_grid_sin(face(3,2))-b*y*phi_grid_cos(face(3,2))) / &
               (b*n(2)*phi_grid_cos(face(3,2))-a*n(1)*phi_grid_sin(face(3,2)))

          if (solutions_p(2).gt.1.e-15_dp.and.solutions_p(1).lt.1.e100_dp) distance(3,2) = solutions_p(2)
          
       end if

    end if

    ! Get the next cell face by checking the shortest distance
    ! next_face(1) --> 1 = radial, 2 = theta, 3 = phi
    ! next_face(2) --> cell face number

    ! First try to find a 'large' distance solution

    face_distance = 1.e100_dp
    
    do j=1,3
       do i=1,3
          
          if (distance(i,j).gt.1.e-9_dp.and.distance(i,j).lt.face_distance) then
             
             face_distance = distance(i,j)
             face_loc(1) = i
             face_loc(2) = j

          end if
          
       end do
    end do

    if (face_loc(1).eq.0) then

       ! If no 'large' solution than take 'small solution

       face_distance = 1.e100_dp

       do j=1,3
          do i=1,3

             if (distance(i,j).gt.1.e-12_dp.and.distance(i,j).lt.face_distance) then

                face_distance = distance(i,j)
                face_loc(1) = i
                face_loc(2) = j

             end if

          end do
       end do
       
       ! face_distance = minval(distance, mask=(distance > 1.e-12_dp))
       ! face_loc = minloc(distance, mask=(distance > 1.e-12_dp))

       if (face_loc(1).eq.0) then

          call cartesian_spherical(x,y,z,r,theta,phi)
          
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 031"
          write (11,*) solutions_r
          write (11,*) solutions_t
          write (11,*) solutions_p
          write (11,*) distance
          write (11,*) current_face
          write (11,*) cell
          write (11,*) face
          write (11,*) r,theta,phi
          write (11,*) n
          close (11)

          cell_error = .true.

       end if

    end if

    next_face(1) = face_loc(1)
    next_face(2) = face(face_loc(1),face_loc(2))

    if (next_face(2).eq.-999) then

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 033"
       close (11)
       cell_error = .true.

    else

       call next_cell(current_face, next_face, cell, cell_out)

    end if

    if (next_face(1).eq.1.and.next_face(2).eq.nr) grid_exit = .true.

    if (current_face(1).eq.1.and.current_face(2).eq.cell_depth.and.&
         next_face(1).eq.1.and.next_face(2).eq.cell_depth) then

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 034"
       close (11)
       cell_error = .true.

    else if (cell_out(1).eq.nr.and..not.grid_exit) then

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 035"
       close (11)
       cell_error = .true.

    else if (cell_out(2).eq.ntheta) then

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 036"
       close (11)
       cell_error = .true.

    else if (cell(1).eq.cell_out(1).and.cell(2).eq.cell_out(2).and.&
         cell(3).eq.cell_out(3).and..not.grid_exit) then

       open (11, file=trim(error_log), position="append")
       write (11,*) "error 037"
       close (11)
       cell_error = .true.

    end if

  end subroutine cell_face

  subroutine write_output

    integer               :: i, j, k, m, n
    real(dp)              :: planck_flux, error(nx,ny,5), jansky, norm, e_pack, dummy, pol, dpol
    real(dp), allocatable :: transport(:,:,:,:)
    logical               :: exist
    
    error = 0._dp

    ! Stokes error
    
    do k=1,4
       do j=1,ny
          do i=1,nx

             ! Error: stdev = sqrt( sum_x2/n - mean^2 )
             ! mean = sum_x / n

             if (detector(i,j,k,3).gt.0._dp) then

                dummy = (detector(i,j,k,2)/detector(i,j,k,3)) - (detector(i,j,k,1)/detector(i,j,k,3))**2
                if (dummy.gt.0._dp) error(i,j,k) = sqrt(dummy) * sqrt(detector(i,j,k,3))

             end if

          end do
       end do
    end do

    ! Degree of polarization error
    
    do j=1,ny
       do i=1,nx

          if (detector(i,j,2,1)**2+detector(i,j,3,1)**2.gt.0._dp) then
          
             pol = sqrt(detector(i,j,2,1)**2+detector(i,j,3,1)**2)

             dpol = sqrt( ( (detector(i,j,2,1)*error(i,j,2))**2 + (detector(i,j,3,1)*error(i,j,3))**2 ) / &
                  ( 2._dp*(detector(i,j,2,1)**2+detector(i,j,3,1)**2) ) )

          end if

          if (detector(i,j,1,1).gt.0._dp) error(i,j,5) = &
               (pol/detector(i,j,1,1)) * sqrt( (dpol/pol)**2 + (error(i,j,1)/detector(i,j,1,1))**2  )
          
       end do
    end do

    if (phase_curve) then

       ! Phase curve

       inquire(file='output/'//adjustl(trim(output_name))//'/output/phase.dat', exist=exist)

       if (exist) then

          open (100, file='output/'//adjustl(trim(output_name))//'/output/phase.dat', &
               status='old', position='append', action='write')

       else if (.not.exist) then

          open (100, file='output/'//adjustl(trim(output_name))//'/output/phase.dat', &
               status='new', action='write')
          write (100,*) "# Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]"
          write (100,*)

       end if

       ! Phase [deg] - Stokes I, I error, Q, Q error, U, U error, V, V error [W m-2 micron-1]
       
       if (det_phi*180._dp/pi.lt.1._dp) then

          write (100,*) 0._dp, detector(1,1,1,1)*1.e-6_dp, error(1,1,1)*1.e-6_dp, &
               detector(1,1,2,1)*1.e-6_dp, error(1,1,2)*1.e-6_dp, detector(1,1,3,1)*1.e-6_dp, &
               error(1,1,3)*1.e-6_dp, detector(1,1,4,1)*1.e-6_dp, error(1,1,4)*1.e-6_dp

       else if (det_phi*180._dp/pi.gt.179_dp) then

          write (100,*) 180._dp, detector(1,1,1,1)*1.e-6_dp, error(1,1,1)*1.e-6_dp, &
               detector(1,1,2,1)*1.e-6_dp, error(1,1,2)*1.e-6_dp, detector(1,1,3,1)*1.e-6_dp, &
               error(1,1,3)*1.e-6_dp, detector(1,1,4,1)*1.e-6_dp, error(1,1,4)*1.e-6_dp
          
       else
          
          write (100,*) det_phi*180._dp/pi, detector(1,1,1,1)*1.e-6_dp, error(1,1,1)*1.e-6_dp, &
               detector(1,1,2,1)*1.e-6_dp, error(1,1,2)*1.e-6_dp, detector(1,1,3,1)*1.e-6_dp, &
               error(1,1,3)*1.e-6_dp, detector(1,1,4,1)*1.e-6_dp, error(1,1,4)*1.e-6_dp

       end if

       close (100)

    else if (imaging_mono.or.imaging_broad) then
          
       ! Write FITS images

       call write_fits_3D('stokes.fits', detector(:,:,:,1)*1.e-6_dp/(pixel_scale*pixel_scale), nx, ny, 4)
       call write_fits_3D('error.fits', error, nx, ny, 5)
       
       ! Write photometry

       if (imaging_mono) then

          open (100, file='output/'//adjustl(trim(output_name))//'/output/photometry.dat', &
               status='new', action='write')

          write (100,*) "# Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]"
          write (100,*)

          ! Wavelength [micron] - Stokes I, I error, Q, Q error, U, U error, V, V error [W m-2 micron-1]
          write (100,*) wavelength*1.e6_dp, 1.e-6_dp*photometry(1), 1.e-6_dp*photometry(2), 1.e-6_dp*photometry(3), &
               1.e-6_dp*photometry(4), 1.e-6_dp*photometry(5), 1.e-6_dp*photometry(6), &
               1.e-6_dp*photometry(7), 1.e-6_dp*photometry(8)

          close (100)

       end if
       
    else if (spectrum) then

       ! Spectrum
       
       ! [W m-2 m-1] -> [Jy]
       jansky = (wavelength * wavelength * 1.e26) / cc

       inquire(file='output/'//adjustl(trim(output_name))//'/output/spectrum.dat', exist=exist)

       if (exist) then

          open (100, file='output/'//adjustl(trim(output_name))//'/output/spectrum.dat', &
               status='old', position='append', action='write')

       else if (.not.exist) then

          open (100, file='output/'//adjustl(trim(output_name))//'/output/spectrum.dat', &
               status='new', action='write')

          write (100,*) "# Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]"
          write (100,*)

       end if

       ! Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]
       write (100,*) wavelength*1.e6_dp, 1.e-6_dp*detector(1,1,1,1), 1.e-6_dp*detector(1,1,2,1), &
            1.e-6_dp*detector(1,1,3,1), 1.e-6_dp*detector(1,1,4,1)

       close (100)

    end if

    if (photon_source.eq.1) then

       ! Flux normalization constants

       if ((phase_curve.and.det_phi.lt.pi/180._dp).or..not.phase_curve) then

          inquire(file='output/'//adjustl(trim(output_name))//'/output/normalization.dat', exist=exist)

          if (exist) then

             open (100, file='output/'//adjustl(trim(output_name))//'/output/normalization.dat', &
                  status='old', position='append', action='write')

          else if (.not.exist) then

             open (100, file='output/'//adjustl(trim(output_name))//'/output/normalization.dat', &
                  status='new', action='write')

          end if

          call planck_function(t_star, planck_flux)

          ! Wavelength [micron] - Norm1 [W m-2 micron-1] - Norm2 [W m-2 micron-1]
          write (100,*) wavelength*1.e6_dp, 1.e-6_dp*planck_flux*r_star*r_star / &
               (distance_planet*distance_planet), 1.e-6_dp*planck_flux*rfront(nr)*rfront(nr)*r_star*r_star / &
               (orbit*orbit*distance_planet*distance_planet)

          close (100)

       end if

    else if (photon_source.eq.2) then

       ! Write cell luminosity
       
       if (imaging_mono) call write_fits_3D("cell_luminosity.fits", cell_luminosity, nr, ntheta, nphi)

       ! Write emitted and emergent luminosity
       
       e_pack = emissivity_cumulative(nr-1,ntheta-1,nphi-1)/dble(packages)

       inquire(file='output/'//trim(output_name)//'/output/luminosity.dat', exist=exist)

       if (exist) then
          
          open (100, file='output/'//trim(output_name)//'/output/luminosity.dat', &
               status='old', position='append', action='write')
          
       else if (.not.exist) then

          open (100, file='output/'//trim(output_name)//'/output/luminosity.dat', status='new', action='write')

          write (100,*) "# Wavelength [deg] - Emitted luminosity [W micron-1] - &
               & Emergent luminosity [W micron-1] - Emergent luminosity [a.u.]"
          write (100,*)
          
       end if

       write (100,*) wavelength, sum(flux_emitted)*e_pack*1.e-6_dp, sum(flux_exit)*e_pack*1.e-6_dp, sum(flux_exit)

       close (100)
       
    end if

    ! Write cell depth

    if (imaging_mono.or.spectrum) then
    
       inquire(file='output/'//adjustl(trim(output_name))//'/output/cell_depth.dat', exist=exist)

       if (exist) then
          
          open (100, file='output/'//adjustl(trim(output_name))//'/output/cell_depth.dat', &
               status='old', position='append', action='write')
          
       else if (.not.exist) then

          open (100, file='output/'//adjustl(trim(output_name))//'/output/cell_depth.dat', status='new', action='write')

          write (100,*) "# Wavelength [micron] - Cell depth"
          write (100,*)
          
       end if

       write (100,*) wavelength*1.e6_dp, cell_depth

       close (100)

    end if
    
    ! Write flow output

    if (flow_global) then
    
       allocate (transport(3,0:nr-1,0:ntheta-1,0:nphi-1))
       transport = 0._dp

       do k=0,nphi-1
          do j=0,ntheta-1
             do i=cell_depth,nr-1
                do m=1,3
                   do n=1,threads
                      transport(m,i,j,k) = transport(m,i,j,k) + cell_flow_global(n,m,i,j,k)
                   end do
                end do
                norm = sqrt(transport(1,i,j,k)**2+transport(2,i,j,k)**2+transport(3,i,j,k)**2)
                if (norm.gt.0._dp) then
                   transport(1,i,j,k) = transport(1,i,j,k)/norm
                   transport(2,i,j,k) = transport(2,i,j,k)/norm
                   transport(3,i,j,k) = transport(3,i,j,k)/norm
                end if
             end do
          end do
       end do

       call write_fits_4D("flow_global.fits", transport, 3, nr, ntheta, nphi)

       deallocate (transport)

    end if

    if (flow_theta) then

       ! Write latitudinal energy transport, normalized to the total emergent flux

       allocate (transport(4,0:nr-1,0:ntheta-1,0:nphi-1))
       transport = 0._dp

       do k=0,nphi-1
          do j=0,ntheta-1
             do i=cell_depth,nr-1

                transport(1,i,j,k) = sum(cell_flow(:,1,i,j,k))
                transport(2,i,j,k) = sum(cell_flow(:,2,i,j,k))
                transport(3,i,j,k) = sum(cell_flow(:,3,i,j,k))
                transport(4,i,j,k) = sum(cell_flow(:,4,i,j,k))
                
             end do
          end do
       end do

       transport = transport / sum(flux_exit)

       call write_fits_4D("flow_latitudinal.fits", transport, 4, nr, ntheta, nphi)

       deallocate (transport)

    end if

  end subroutine write_output

  subroutine write_fits_3D(filename, array, n1, n2, n3)

    ! Write 2D array to FITS file

    integer,      intent(in) :: n1, n2, n3
    integer                  :: status, unit, blocksize, bitpix, naxis, naxes(3), group, fpixel
    real(dp),     intent(in) :: array(n1,n2,n3)
    character(*), intent(in) :: filename
    logical                  :: simple, extend
    character(300)           :: fits_path

    status = 0
    blocksize = 1
    simple = .true.
    naxis = 3
    naxes(1) = n1
    naxes(2) = n2
    naxes(3) = n3
    extend = .true.
    group = 1
    fpixel = 1
    bitpix = -64

    write (fits_path, "('!output/',(a),'/output/',(a))") trim(output_name), trim(filename)

    call ftgiou(unit, status)
    call ftinit(unit, trim(fits_path), blocksize, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, n1*n2*n3, array, status)
    call ftclos(unit, status)
    call ftfiou(unit, status)

  end subroutine write_fits_3D

  subroutine write_fits_4D(filename, array, n1, n2, n3, n4)

    ! Write 4D array to FITS file

    integer,      intent(in) :: n1, n2, n3, n4
    integer                  :: status, unit, blocksize, bitpix, naxis, naxes(4), group, fpixel
    real(dp),     intent(in) :: array(n1,n2,n3,n4)
    character(*), intent(in) :: filename
    logical                  :: simple, extend
    character(300)           :: fits_path

    status = 0
    blocksize = 1
    simple = .true.
    naxis = 4
    naxes(1) = n1
    naxes(2) = n2
    naxes(3) = n3
    naxes(4) = n4
    extend = .true.
    group = 1
    fpixel = 1
    bitpix = -64

    write (fits_path, "('!output/',(a),'/output/',(a))") trim(output_name), trim(filename)

    call ftgiou(unit, status)
    call ftinit(unit, trim(fits_path), blocksize, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, n1*n2*n3*n4, array, status)
    call ftclos(unit, status)
    call ftfiou(unit, status)

  end subroutine write_fits_4D
  
  subroutine output(log_number)

    integer, intent(in) :: log_number
    integer             :: i, j, k, hours, minutes, seconds, log_size
    real(dp)            :: angular_size, total_optical_depth, cpu_total, e_pack, norm, norm2, planck_flux
    logical             :: file_exists
    character(200)      :: hostname

    if (log_number.eq.1) then

       angular_size = 2._dp*atan(rfront(nr)/distance_planet)*3600._dp*180._dp/pi*1000._dp

       write (out_unit,'(a)') "########################################################"
       write (out_unit,'(a)') "              _     ___   _____   ___   ___             "
       write (out_unit,'(a)') "             /_\   | _ \ |_   _| | __| / __|            "
       write (out_unit,'(a)') "            / _ \  |   /   | |   | _|  \__ \            "
       write (out_unit,'(a)') "           /_/ \_\ |_|_\   |_|   |___| |___/            "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "  Atmospheric Radiative Transfer for Exoplanet Science  "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "                      Developed by:                     "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "                      Tomas Stolker                     "
       write (out_unit,'(a)') "                  T.Stolker [at] uva.nl                 "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "        Anton Pannekoek Institute for Astronomy         "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)')
       write (out_unit,'(a)') "       Please cite Stolker et al. 2017 whenever         "
       write (out_unit,'(a)') "       ARTES results are used in a publication.         "
       write (out_unit,'(a)')
       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)')
       write (out_unit,'(a)') "--> Initialization"
       write (out_unit,'(a)')
       write (out_unit,'(a,a)') "Atmosphere: ", atmosphere
       write (out_unit,'(a,i0)') "Computer cores: ", cores
       write (out_unit,'(a,i0)') "Threads in use: ", threads
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--> Build planet atmosphere"
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Planet radius [km]: ", rfront(nr)/1000._dp
       write (out_unit,'(a,es8.2)') "Atmosphere height [km]: ", (rfront(nr)-rfront(0))/1000._dp
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Oblateness: ", oblateness
       if (oblateness.gt.0._dp) write (out_unit,'(a,es8.2)') "Eccentricity: ", &
            sqrt( 1._dp - (oblate_z*oblate_z) / (oblate_x*oblate_x) )
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Surface albedo: ", surface_albedo
       write (out_unit,'(a)') ""
       write (out_unit,'(a,i0)') "Radial grid cells: ", nr
       write (out_unit,'(a,i0)') "Latitudinal grid cells: ", ntheta
       write (out_unit,'(a,i0)') "Longitudial grid cells: ", nphi
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2,a,es8.2)') "Field of view [mas x mas]: ", x_fov, " x ", y_fov
       write (out_unit,'(a,es8.2,a,es8.2)') "Planet angular diameter [mas]: ", angular_size
       write (out_unit,'(a,es8.2,a,es8.2)') "Pixel scale [mas pixel-1]: ", pixel_scale
       write (out_unit,'(a)') ""

    else if (log_number.eq.2) then

       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--> Photon transfer"
       write (out_unit,'(a)') ""
       if (photon_source.eq.1) then
          write (out_unit,'(a,a)') "Photon source: star"
       else if (photon_source.eq.2) then
          write (out_unit,'(a,a)') "Photon source: planet"
       end if
       write (out_unit,'(a,es8.2)') "Emitted photons: ", dble(packages)
       write (out_unit,'(a)') ""
       if (photon_source.eq.1.and..not.phase_curve) then
          write (out_unit,'(a,es8.2)') "Phase angle [deg]: ", phase_observer
          write (out_unit,'(a)') ""
       end if
       if (phase_curve.or.imaging_mono) then
          write (out_unit,'(a,es8.2)') "Monochromatic wavelength [micron]: ", wavelength*1.e6_dp
       else if (imaging_broad) then
          write (out_unit,'(a,es8.2,a,es8.2)') "Broadband wavelength range [micron]: ", wavelengths(1)*1.e6, &
               " - ", wavelengths(size(wavelengths))*1.e6
       else if (spectrum) then
          write (out_unit,'(a,es8.2,a,es8.2)') "Spectrum wavelength range [micron]: ", wavelengths(1)*1.e6, &
               " - ", wavelengths(size(wavelengths))*1.e6
       end if
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Stellar luminosity [W]: ", 4._dp*pi*r_star*r_star*sb*t_star**4
       write (out_unit,'(a)') ""
       if (.not.spectrum) then
          write (out_unit,'(a)') "Total optical depth:"
          do i=0,ntheta-1
             do j=0,nphi-1
                total_optical_depth = 0._dp
                do k=cell_depth,nr-1
                   total_optical_depth = total_optical_depth + ( (rfront(k+1)-rfront(k)) * &
                        cell_opacity(k,i,j,wl_count) )
                end do
                write (out_unit,'(a,i0,a,i0,a,es10.4)') "[Theta, phi] = [", i, ", ", j, "] --> ", total_optical_depth
             end do
          end do
          write (out_unit,'(a)') ""
          write (out_unit,'(a)') "Scattering optical depth:"
          do i=0,ntheta-1
             do j=0,nphi-1
                total_optical_depth = 0._dp
                do k=cell_depth,nr-1
                   total_optical_depth = total_optical_depth + ( (rfront(k+1)-rfront(k)) * &
                        cell_scattering_opacity(k,i,j,wl_count) )
                end do
                write (out_unit,'(a,i0,a,i0,a,es10.4)') "[Theta, phi] = [", i, ", ", j, "] --> ", total_optical_depth
             end do
          end do
          write (out_unit,'(a)') ""
          write (out_unit,'(a)') "Absorption optical depth:"
          do i=0,ntheta-1
             do j=0,nphi-1
                total_optical_depth = 0._dp
                do k=cell_depth,nr-1
                   total_optical_depth = total_optical_depth + ( (rfront(k+1)-rfront(k)) * &
                        cell_absorption_opacity(k,i,j,wl_count) )
                end do
                write (out_unit,'(a,i0,a,i0,a,es10.4)') "[Theta, phi] = [", i, ", ", j, "] --> ", total_optical_depth
             end do
          end do
          write (out_unit,'(a)') ""
       end if

    else if (log_number.eq.3) then

       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""

       if (imaging_mono.or.imaging_broad) then

          if (photon_source.eq.1) then

             call planck_function(t_star, planck_flux)

             norm = planck_flux*rfront(nr)*rfront(nr)*r_star*r_star / (orbit*orbit*distance_planet*distance_planet)
             norm2 = planck_flux*r_star*r_star/(distance_planet*distance_planet)

             if (photometry(1).gt.0._dp) then

                write (out_unit,'(a)') "Planet integrated flux"
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Stokes I [W m-2 micron-1]: ", photometry(1)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes Q [W m-2 micron-1]: ", photometry(3)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes U [W m-2 micron-1]: ", photometry(5)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes V [W m-2 micron-1]: ", photometry(7)*1.e-6_dp
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Normalized Stokes I: ", photometry(1)/norm
                write (out_unit,'(a,es9.2)') "Normalized Stokes Q: ", photometry(3)/norm
                write (out_unit,'(a,es9.2)') "Normalized Stokes U: ", photometry(5)/norm
                write (out_unit,'(a,es9.2)') "Normalized Stokes V: ", photometry(7)/norm
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes I: ", photometry(1)/norm2
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes Q: ", photometry(3)/norm2
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes U: ", photometry(5)/norm2
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes V: ", photometry(7)/norm2
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "-Q/I: ", -photometry(3)/photometry(1)
                write (out_unit,'(a,es9.2)') " U/I: ", photometry(5)/photometry(1)
                write (out_unit,'(a,es9.2)') " V/I: ", photometry(7)/photometry(1)
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es8.2,a,es8.2)') "Degree of polarization [%]: ", &
                     100._dp * photometry(10), " +/- ", 100._dp * photometry(11)
                write (out_unit,'(a,es9.2)') "Direction of polarization [deg]: ", &
                     0.5*atan2(photometry(5), photometry(3)) * 180._dp/pi
                write (out_unit,'(a)') ""

             else

                write (out_unit,'(a)') "Error in Stokes I"
                write (out_unit,'(a)') ""

             end if

             write (out_unit,'(a)') "--------------------------------------------------------"
             write (out_unit,'(a)') ""

          else if (photon_source.eq.2) then

             if (photometry(1).gt.0._dp) then

                write (out_unit,'(a)') "Planet integrated flux"
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Stokes I [W m-2 micron-1]: ", photometry(1)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes Q [W m-2 micron-1]: ", photometry(3)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes U [W m-2 micron-1]: ", photometry(5)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes V [W m-2 micron-1]: ", photometry(7)*1.e-6_dp
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "-Q/I: ", -photometry(3)/photometry(1)
                write (out_unit,'(a,es9.2)') " U/I: ", photometry(5)/photometry(1)
                write (out_unit,'(a,es9.2)') " V/I: ", photometry(7)/photometry(1)
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es8.2,a,es8.2)') "Degree of polarization [%]: ", &
                     100._dp * photometry(10), " +/- ", 100._dp * photometry(11)
                write (out_unit,'(a,es9.2)') "Direction of polarization [deg]: ", &
                     0.5*atan2(photometry(5), photometry(3)) * 180._dp/pi
                write (out_unit,'(a)') ""

             else

                write (out_unit,'(a)') "Error: Stokes I is zero"
                write (out_unit,'(a)') ""

             end if

             write (out_unit,'(a)') "--------------------------------------------------------"
             write (out_unit,'(a)') ""

             if (imaging_mono.and.photon_source.eq.2) then

                e_pack = emissivity_cumulative(nr-1,ntheta-1,nphi-1)/dble(packages)
                
                write (out_unit,'(a,es9.2)') "Total emitted flux [W micron-1] =", sum(flux_emitted)*e_pack*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Total emergent flux [W micron-1] =", sum(flux_exit)*e_pack*1.e-6_dp
                write (out_unit,'(a)') ""
                write (out_unit,'(a)') "--------------------------------------------------------"
                write (out_unit,'(a)') ""

             end if

          end if

       end if

    else if (log_number.eq.4) then

       t_end = omp_get_wtime()
       cpu_total = t_end - t_start

       hours   = int( cpu_total / 3600._dp )
       minutes = int( cpu_total / 60._dp - dble(hours)*60._dp )
       seconds = int( cpu_total - dble(hours)*3600._dp - dble(minutes)*60._dp )

       inquire (file=trim(error_log),size=log_size)

       write (out_unit,'(a,i0.2,a,i0.2,a,i0.2)') "CPU time [hour:min:sec]: ", hours, ":", minutes, ":", seconds
       write (out_unit,'(a)') ""
       if (log_size.ne.0) then
          write (out_unit,'(a)') "WARNING: check error log!"
          write (out_unit,'(a)') ""
       end if
       write (out_unit,'(a)') "########################################################"

       if (log_file) close (10)

       if (len(trim(email)).ne.0) then

          call system('hostname > hostname.txt')
          open (100, file='hostname.txt')
          read (100,*) hostname
          close (100, status='delete')

          inquire(file="/usr/bin/mail", exist=file_exists)

          if (file_exists) then

             open (100, file='mail.txt')
             write (100,'(a,a,a,a,a)') "Job with input ", adjustl(trim(atmosphere)), " on ", &
                  adjustl(trim(hostname)), " is finished."
             write (100,'(a)')
             write (100,'(a)') "Have a nice day!"

             call system('mail -s "ARTES is finished" ' // adjustl(trim(email)) // ' < mail.txt')

             close (100, status='delete')

          else if (.not.file_exists) then

             inquire(file="/usr/sbin/ssmtp", exist=file_exists)

             if (file_exists) then

                open (100, file='mail.txt')

                write (100,'(a,a)') "To:" // adjustl(trim(email))
                write (100,'(a)') "From:ARTES"
                write (100,'(a)') "Subject: ARTES is finished"
                write (100,'(a)') ""
                write (100,'(a,a,a,a,a)') "Job with input ", adjustl(trim(input_file)), " on ", &
                     adjustl(trim(hostname)), " is finished."
                write (100,'(a)')
                write (100,'(a)') "Have a nice day!"

                call system('ssmtp ' // adjustl(trim(email)) // ' < mail.txt')

                close (100, status='delete')

             else

                open (11, file=trim(error_log), position="append")
                write (11,*) "error 038"
                close (11)

             end if

          end if

       end if

       close (11)

    end if

  end subroutine output

  subroutine quadratic_equation(a, b, c, solutions)

    ! Solve a quadratic equation

    real(dp), intent(in)  :: a, b, c
    real(dp), intent(out) :: solutions(2)
    real(dp)              :: discriminant, q

    solutions = 0._dp
    discriminant = b*b - 4._dp*a*c

    if (discriminant.ge.0._dp) then

       q = -0.5_dp * ( b + sign(1._dp,b) * sqrt(discriminant)  )
       if (abs(a).gt.1.e-100_dp) solutions(1) = q/a
       if (abs(q).gt.1.e-100_dp) solutions(2) = c/q

    end if

  end subroutine quadratic_equation

  subroutine init_random_seed

    ! Initialize a random seed for the random number generator

    integer               :: clock, i, n
    integer, allocatable  :: seed(:)

    call random_seed(size = n)

    allocate (seed(n))
    seed = 0

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)

    call random_seed(put = seed)

    deallocate (seed)

  end subroutine init_random_seed

  subroutine random(thread_id, xi)

    ! Marsaglia & Zaman (1994)

    integer,  intent(in)  :: thread_id
    integer               :: imz
    real(dp), intent(out) :: xi

    imz = state(thread_id+1,1) - state(thread_id+1,3)

    if (imz.lt.0) imz = imz + 2147483579

    state(thread_id+1,1) = state(thread_id+1,2)
    state(thread_id+1,2) = state(thread_id+1,3)
    state(thread_id+1,3) = imz
    state(thread_id+1,4) = 69069 * state(thread_id+1,4) + 1013904243

    imz = imz + state(thread_id+1,4)

    xi = 0.5_dp + 0.23283064e-9_dp * imz

    if (xi.gt.0._dp.and.xi.lt.1._dp) then
    else
       open (11, file=trim(error_log), position="append")
       write (11,*) "error 055"
       write (11,*) thread_id, imz, xi
       write (11,*) state(thread_id+1,1)
       write (11,*) state(thread_id+1,2)
       write (11,*) state(thread_id+1,3)
       write (11,*) state(thread_id+1,4)
       close (11)
    end if
       
  end subroutine random

  subroutine argument_input(atmosphere)

    ! Check the command line arguments

    character(100), intent(out) :: atmosphere
    character(100)              :: dummy_char
    real(dp                   ) :: dummy_dble

    atmosphere = ""

    if (iargc().le.1) then
       
       write (6,'(a)') "How to run ARTES:"
       write (6,'(a)') "./bin/ARTES [inputDirectory] [photons] -o [outputDirectory] -k [keyWord]=[value]"
       call exit(0)
    
    else

       call getarg(1, atmosphere)
       call getarg(2, dummy_char)

       read (dummy_char,*) dummy_dble
       packages = int(dummy_dble,kind=16)

    end if
    
  end subroutine argument_input

  subroutine argument_keywords

    integer        :: i
    character(100) :: arg, keyword_dummy, key_word, key_value
    logical        :: exists

    do i = 1, command_argument_count()

       call get_command_argument(i, arg)

       select case (arg)
       case ('-o')
          
          call get_command_argument(i+1, output_name)

          ! Make output directories
          call system("rm -rf output/" // trim(output_name))
          call system('mkdir -p output/' // trim(output_name) )
          call system('mkdir -p output/' // trim(output_name) // '/input' )
          call system('mkdir -p output/' // trim(output_name) // '/output' )
          call system('mkdir -p output/' // trim(output_name) // '/plot' )

          ! Copy input files to output directory
          call system('cp '//adjustl(trim(input_file))//' output/'//adjustl(trim(output_name))//'/input/')
          call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'atmosphere.in')// &
               ' output/'//adjustl(trim(output_name))//'/input/')
          call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'atmosphere.fits')// &
               ' output/'//adjustl(trim(output_name))//'/input/')
          inquire(file=trim(input_file(1:len_trim(input_file)-8))//'atmosphere.dat', exist=exists)
          if (exists) call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'atmosphere.dat')// &
               ' output/'//adjustl(trim(output_name))//'/input/')
          inquire(file=trim(input_file(1:len_trim(input_file)-8))//'pressureTemperature.dat', exist=exists)
          if (exists) call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'pressureTemperature.dat')// &
               ' output/'//adjustl(trim(output_name))//'/input/')
          
       case ('-k')
          
          call get_command_argument(i+1, keyword_dummy)
          call get_key_value(keyword_dummy, key_word, key_value)
          call input_parameters(key_word, key_value)

          open (100, file='output/'//trim(output_name)//'/input/artes.in', status='old', position='append', action='write')
          write (100,*) new_line('a')//trim(keyword_dummy)
          close (100)
          
       end select
       
    end do

  end subroutine argument_keywords
  
  subroutine readfile(input_file, data, nlines)

    ! Check number of lines in a data file and read into array

    integer,                     intent(out) :: nlines
    integer                                  :: j, line_counter, ios
    integer, parameter                       :: max_lines = 10000
    character(100),              intent(in)  :: input_file 
    character(500), allocatable, intent(out) :: data(:)
    character(100)                           :: junk

    line_counter = 0

    open (100,file=input_file)

    do j=1,max_lines

       read (100,'(a500)',iostat=ios) junk

       if (ios.ne.0) exit

       if (j.eq.max_lines) then

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 039"
          close (11)
          stop

       end if

       line_counter = line_counter + 1

    end do

    rewind(100)

    allocate (data(line_counter))

    do j=1,line_counter

       read (100,'(a500)') data(j)

    enddo

    nlines = line_counter

    close (100)

  end subroutine readfile

  subroutine input_parameters(key_word, key_value)

    integer                     :: key_value_int
    real(dp)                    :: key_value_double
    character(100), intent(in)  :: key_word, key_value
    character(100)              :: key_value_char

    key_value_int = 0
    key_value_double = 0._dp
    key_value_char = ""

    select case(key_word)

    case("general:log")
       if (key_value.eq."on") then
          log_file = .true.
          out_unit = 10
       else if (key_value.eq."off") then
          log_file = .false.
          out_unit = 6
       end if
    case("general:email")
       read (key_value, '(a100)') key_value_char
       email = key_value_char
    case("photon:source")
       if (key_value.eq."star") then
          photon_source = 1
       else if (key_value.eq."planet") then
          photon_source = 2
       end if
    case("photon:fstop")
       read (key_value, '(e50.0)') key_value_double
       fstop = key_value_double
    case("photon:minimum")
       read (key_value, '(e50.0)') key_value_double
       photon_minimum = key_value_double
    case("photon:weight")
       if (key_value.eq."on") then
          thermal_weight = .true.
       else if (key_value.eq."off") then
          thermal_weight = .false.
       end if
    case("photon:scattering")
       if (key_value.eq."on") then
          photon_scattering = .true.
       else if (key_value.eq."off") then
          photon_scattering = .false.
       end if
    case("photon:emission")
       if (key_value.eq."isotropic") then
          photon_emission = 1
       else if (key_value.eq."biased") then
          photon_emission = 2
       end if
    case("photon:bias")
       read (key_value, '(e50.0)') key_value_double
       photon_bias = key_value_double
    case("star:temperature")
       read (key_value, '(e50.0)') key_value_double
       t_star = key_value_double
    case("star:radius")
       read (key_value, '(e50.0)') key_value_double
       r_star = key_value_double*r_sun
    case("star:direction")
       if (key_value.eq."on") then
          stellar_direction = .true.
       else if (key_value.eq."off") then
          stellar_direction = .false.
       end if
    case("star:theta")
       if (stellar_direction) then
          read (key_value, '(e50.0)') key_value_double
          theta_star = key_value_double*pi/180._dp
          if (theta_star.lt.1.e-3) theta_star = 1.e-3_dp
          if (theta_star.gt.pi-1.e-3) theta_star = pi-1.e-3
       end if
    case("star:phi")
       if (stellar_direction) then
          read (key_value, '(e50.0)') key_value_double
          phi_star = key_value_double*pi/180._dp
       end if
    case("planet:surface_albedo")
       read (key_value, '(e50.0)') key_value_double
       surface_albedo = key_value_double
    case("planet:oblateness")
       read (key_value, '(e50.0)') key_value_double
       oblateness = key_value_double
    case("planet:orbit")
       read (key_value, '(e50.0)') key_value_double
       orbit = key_value_double*au
    case("planet:ring")
       if (key_value.eq."on") then
          ring = .true.
       else if (key_value.eq."off") then
          ring = .false.
       end if
    case("detector:type")
       if (key_value.eq."phase") then
          phase_curve = .true.
       else if (key_value.eq."spectrum") then
          spectrum = .true.
       else if (key_value.eq."imaging_mono") then
          imaging_mono = .true.
       else if (key_value.eq."imaging_broad") then
          imaging_broad = .true.
       end if
    case("detector:theta")
       read (key_value, '(e50.0)') key_value_double
       det_theta = key_value_double*pi/180._dp
       if (det_theta.lt.1.e-3) det_theta = 1.e-3_dp
       if (det_theta.gt.pi-1.e-3) det_theta = pi-1.e-3
    case("detector:phi")
       read (key_value, '(e50.0)') key_value_double
       det_phi = key_value_double*pi/180._dp
    case("detector:pixel")
       read (key_value, '(i50)') key_value_int
       nx = key_value_int
       ny = key_value_int
    case("detector:distance")
       read (key_value, '(e50.0)') key_value_double
       distance_planet = key_value_double*pc
    case("output:flow_global")
       if (key_value.eq."on") then
          flow_global = .true.
       else if (key_value.eq."off") then
          flow_global = .false.
       end if
    case("output:flow_latitudinal")
       if (key_value.eq."on") then
          flow_theta = .true.
       else if (key_value.eq."off") then
          flow_theta = .false.
       end if
    case default
       write (6,*) "Wrong keyword found in input file: ", key_word
       call exit(0)

    end select

  end subroutine input_parameters

  subroutine get_key_value(line,key_word,key_value)

    ! Seperate a key_word and key_value component of a string given key_word=key_value syntax

    character(100), intent(in)  :: line
    character(100), intent(out) :: key_word, key_value

    key_word = ""
    key_value = ""

    key_word=line(1:index(line,'=')-1)
    key_value=line(index(line,'=')+1:len_trim(line))

    if(key_value(1:1).eq.'"'.or.key_value(1:1).eq."'") key_value=key_value(2:len_trim(key_value)-1)

  end subroutine get_key_value

  subroutine peel_thermal(thread_id, x_photon, y_photon, z_photon, stokes_peel_in, cell_in, face_in, cell_error)

    integer,  intent(in) :: cell_in(3), face_in(2), thread_id
    integer              :: cell(3), cell_out(3), current_face(2), next_face(2), ix, iy
    real(dp), intent(in) :: x_photon, y_photon, z_photon, stokes_peel_in(4)
    real(dp)             :: x, y, z, face_distance, photon_weight, tau_cell, tau_total, x_im, y_im, det(3)
    logical              :: grid_exit
    logical, intent(out) :: cell_error

    x = x_photon
    y = y_photon
    z = z_photon

    det(1) = det_dir(1)
    det(2) = det_dir(2)
    det(3) = det_dir(3)

    current_face = face_in
    cell = cell_in

    tau_total = 0._dp
    photon_weight = 0._dp

    do

       call cell_face(x, y, z, det, current_face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)
       
       if (cell_error) then
          
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 046"
          write (11,*) grid_exit
          close (11)
          
          exit
          
       end if

       tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)
       tau_total = tau_total + tau_cell

       x = x + face_distance*det(1)
       y = y + face_distance*det(2)
       z = z + face_distance*det(3)

       if (grid_exit.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

       current_face = next_face
       cell = cell_out

    end do

    if (grid_exit.and.tau_total.lt.50._dp) then

       photon_weight = exp(-tau_total) / (4._dp*pi)

       x_im = y_photon*cos_det_phi - x_photon*sin_det_phi
       y_im = z_photon*sin_det_theta - y_photon*cos_det_theta*sin_det_phi - x_photon*cos_det_theta*cos_det_phi

       ix = int(nx*(x_im+x_max)/(2._dp*x_max))+1
       iy = int(ny*(y_im+y_max)/(2._dp*y_max))+1

       if (photon_weight*stokes_peel_in(1).gt.0._dp.and.photon_weight*stokes_peel_in(1).lt.1.e100_dp) then
       
          detector_thread(ix,iy,1,1,thread_id+1) = detector_thread(ix,iy,1,1,thread_id+1) + photon_weight*stokes_peel_in(1)
          detector_thread(ix,iy,1,2,thread_id+1) = detector_thread(ix,iy,1,2,thread_id+1) + (photon_weight*stokes_peel_in(1))**2
          detector_thread(ix,iy,1,3,thread_id+1) = detector_thread(ix,iy,1,3,thread_id+1) + 1._dp

       else

          open (11, file=trim(error_log), position="append")
          write (11,*) "error 051"
          write (11,*) photon_weight, stokes_peel_in(1)
          close (11)

       end if
       
    end if

  end subroutine peel_thermal

  subroutine peel_surface(thread_id, x_photon, y_photon, z_photon, stokes_peel_in, cell_in, face_in)

    integer,  intent(in) :: cell_in(3), face_in(2), thread_id
    integer              :: cell(3), cell_out(3), current_face(2), next_face(2), ix, iy
    real(dp), intent(in) :: x_photon, y_photon, z_photon, stokes_peel_in(4)
    real(dp)             :: x, y, z, face_distance, photon_weight, surface_normal(3), det_sphere(3), normal_sphere(3)
    real(dp)             :: tau_cell, tau_total, x_im, y_im, det(3), cos_angle, norm
    logical              :: grid_exit, cell_error

    ! Surface normal on spherical surface at point of reflection and scale for oblate surface
    surface_normal(1) = x_photon / ( oblate_x * oblate_x )
    surface_normal(2) = y_photon / ( oblate_y * oblate_y )
    surface_normal(3) = z_photon / ( oblate_z * oblate_z )

    ! Make the vector a unit vector
    norm = sqrt(surface_normal(1)*surface_normal(1)+surface_normal(2)*surface_normal(2)+surface_normal(3)*surface_normal(3))
    surface_normal = surface_normal / norm

    x = x_photon
    y = y_photon
    z = z_photon

    det(1) = det_dir(1)
    det(2) = det_dir(2)
    det(3) = det_dir(3)

    ! Determine angle between surface normal and detector direction
    
    call cartesian_spherical(det(1), det(2), det(3), det_sphere(1), det_sphere(2), det_sphere(3))
    call cartesian_spherical(surface_normal(1), surface_normal(2), surface_normal(3), &
         normal_sphere(1), normal_sphere(2), normal_sphere(3))

    cos_angle = sin(det_sphere(2))*cos(det_sphere(3))*sin(normal_sphere(2))*cos(normal_sphere(3)) + &
         sin(det_sphere(2))*sin(det_sphere(3))*sin(normal_sphere(2))*sin(normal_sphere(3)) + &
         cos(det_sphere(2))*cos(normal_sphere(2))

    ! Angle should be smaller than 90 degrees otherwise photon is directed into the surface

    if (cos_angle.gt.0._dp) then

       current_face = face_in

       cell(1) = cell_in(1) + 1
       cell(2) = cell_in(2)
       cell(3) = cell_in(3)

       ! Calculate optical depth

       tau_total = 0._dp
       photon_weight = 0._dp

       do

          call cell_face(x, y, z, det, current_face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

          if (cell_error) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 042"
             close (11)
          end if

          tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)
          tau_total = tau_total + tau_cell

          x = x + face_distance*det(1)
          y = y + face_distance*det(2)
          z = z + face_distance*det(3)

          if (grid_exit.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

          current_face = next_face
          cell = cell_out

       end do

       if (grid_exit.and.tau_total.lt.50._dp) then

          ! Lambertian scattering probability: P(mu,phi) dmu dphi = 2mu dmu dphi/(2pi) -> P(mu,phi) = mu/pi

          ! Weight photon for optical depth and scattering angle

          photon_weight = exp(-tau_total) * cos_angle / pi

          x_im = y_photon*cos_det_phi - x_photon*sin_det_phi
          y_im = z_photon*sin_det_theta - y_photon*cos_det_theta*sin_det_phi - x_photon*cos_det_theta*cos_det_phi

          ix = int(nx*(x_im+x_max)/(2._dp*x_max))+1
          iy = int(ny*(y_im+y_max)/(2._dp*y_max))+1

          if (photon_weight*stokes_peel_in(1).gt.0._dp.and.photon_weight*stokes_peel_in(1).lt.1.e100_dp) then

             detector_thread(ix,iy,1,1,thread_id+1) = detector_thread(ix,iy,1,1,thread_id+1) + photon_weight*stokes_peel_in(1)
             detector_thread(ix,iy,1,2,thread_id+1) = detector_thread(ix,iy,1,2,thread_id+1) + (photon_weight*stokes_peel_in(1))**2
             detector_thread(ix,iy,1,3,thread_id+1) = detector_thread(ix,iy,1,3,thread_id+1) + 1._dp
             
          else

             open (11, file=trim(error_log), position="append")
             write (11,*) "error 052"
             write (11,*) photon_weight, stokes_peel_in(1)
             close (11)

          end if

       end if

    end if

  end subroutine peel_surface

  subroutine peel_photon(thread_id, x_photon, y_photon, z_photon, stokes_peel_in, direction_photon, cell_in, face_in, cell_error)

    integer,  intent(in) :: cell_in(3), face_in(2), thread_id
    integer              :: cell(3), cell_out(3), current_face(2), next_face(2)
    integer              :: i, angle_upper, angle_lower, ix, iy
    real(dp), intent(in) :: x_photon, y_photon, z_photon, stokes_peel_in(4), direction_photon(3)
    real(dp)             :: x, y, z, face_distance, stokes_out(4), x_im, y_im, phi_old
    real(dp)             :: mu_scatter, phi_scatter, scatter_dummy(16), scatter(4,4), num_check, phi_new
    real(dp)             :: x0(16), x1(16), y0, y1, y_inter, tau_cell, tau_total, photon_weight, det(3), acos_mu_scatter
    logical, intent(out) :: cell_error
    logical              :: grid_exit

    x = x_photon
    y = y_photon
    z = z_photon

    det(1) = det_dir(1)
    det(2) = det_dir(2)
    det(3) = det_dir(3)

    current_face = face_in
    cell = cell_in

    ! Calculate optical depth

    tau_total = 0._dp
    photon_weight = 0._dp
    cell_error = .false.

    do

       call cell_face(x, y, z, det, current_face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

       if (cell_error) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 043"
          close (11)
       end if

       tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)
       tau_total = tau_total + tau_cell

       x = x + face_distance*det(1)
       y = y + face_distance*det(2)
       z = z + face_distance*det(3)

       if (grid_exit.or.cell_error.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

       current_face = next_face
       cell = cell_out

    end do

    if (.not.cell_error) then

       if (grid_exit.and.tau_total.lt.50._dp) then

          ! Weight photon for optical depth
          photon_weight = exp(-tau_total)

          mu_scatter = direction_photon(1)*det(1) + direction_photon(2)*det(2) + direction_photon(3)*det(3)

          if (mu_scatter.ge.1._dp) then
             mu_scatter = 1._dp - 1.e-10_dp
          else if (mu_scatter.le.-1._dp) then
             mu_scatter = -1._dp + 1.e-10_dp
          end if

          ! Scattering matrix

          acos_mu_scatter = acos(mu_scatter)

          if (mod(acos_mu_scatter*180._dp/pi,1._dp).gt.0.5_dp) then

             angle_upper = int(acos_mu_scatter*180._dp/pi) + 2
             angle_lower = int(acos_mu_scatter*180._dp/pi) + 1

          else

             angle_upper = int(acos_mu_scatter*180._dp/pi) + 1
             angle_lower = int(acos_mu_scatter*180._dp/pi)

          end if

          if (angle_upper.eq.1) then

             scatter(1,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,1,1)
             scatter(1,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,2,1)
             scatter(1,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,3,1)
             scatter(1,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,4,1)
             scatter(2,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,5,1)
             scatter(2,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,6,1)
             scatter(2,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,7,1)
             scatter(2,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,8,1)
             scatter(3,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,9,1)
             scatter(3,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,10,1)
             scatter(3,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,11,1)
             scatter(3,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,12,1)
             scatter(4,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,13,1)
             scatter(4,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,14,1)
             scatter(4,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,15,1)
             scatter(4,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,16,1)

          else if (angle_lower.eq.180) then

             scatter(1,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,1,180)
             scatter(1,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,2,180)
             scatter(1,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,3,180)
             scatter(1,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,4,180)
             scatter(2,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,5,180)
             scatter(2,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,6,180)
             scatter(2,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,7,180)
             scatter(2,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,8,180)
             scatter(3,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,9,180)
             scatter(3,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,10,180)
             scatter(3,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,11,180)
             scatter(3,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,12,180)
             scatter(4,1) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,13,180)
             scatter(4,2) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,14,180)
             scatter(4,3) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,15,180)
             scatter(4,4) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,16,180)

          else

             do i=1,16

                x0(i) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,i,angle_lower)
                x1(i) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,i,angle_upper)
                y0 = dble(angle_lower) - 0.5_dp
                y1 = dble(angle_upper) - 0.5_dp
                y_inter = acos_mu_scatter*180._dp/pi
                scatter_dummy(i) = (x1(i)-x0(i)) * (y_inter-y0) / (y1-y0) + x0(i)

                scatter(1,1) = scatter_dummy(1)
                scatter(1,2) = scatter_dummy(2)
                scatter(1,3) = scatter_dummy(3)
                scatter(1,4) = scatter_dummy(4)
                scatter(2,1) = scatter_dummy(5)
                scatter(2,2) = scatter_dummy(6)
                scatter(2,3) = scatter_dummy(7)
                scatter(2,4) = scatter_dummy(8)
                scatter(3,1) = scatter_dummy(9)
                scatter(3,2) = scatter_dummy(10)
                scatter(3,3) = scatter_dummy(11)
                scatter(3,4) = scatter_dummy(12)
                scatter(4,1) = scatter_dummy(13)
                scatter(4,2) = scatter_dummy(14)
                scatter(4,3) = scatter_dummy(15)
                scatter(4,4) = scatter_dummy(16)

             end do

          end if

          ! Check phi difference between old and new direction

          phi_old = atan2(direction_photon(2),direction_photon(1))

          if (phi_old.lt.0._dp) phi_old = phi_old + 2._dp*pi
          if (phi_old.gt.2._dp*pi) phi_old = phi_old - 2._dp*pi

          phi_new = atan2(det(2),det(1))

          if (phi_new.lt.0._dp) phi_new = phi_new + 2._dp*pi
          if (phi_new.gt.2._dp*pi) phi_new = phi_new - 2._dp*pi

          ! Azimuthal scattering angle

          if (abs(direction_photon(3)).lt.1._dp) then

             num_check = (det(3) - direction_photon(3)*mu_scatter) / &
                  ( sqrt(1._dp-mu_scatter*mu_scatter)*sqrt(1._dp-direction_photon(3)*direction_photon(3)) )

             if (abs(num_check).lt.1._dp) then

                phi_scatter = acos(num_check)

             else if (num_check.ge.1._dp) then

                phi_scatter = 0._dp + 1.e-10_dp

             else if (num_check.le.-1._dp) then

                phi_scatter = pi - 1.e-10_dp

             else

                open (11, file=trim(error_log), position="append")
                write (11,*) "error 044"
                write (11,*) num_check
                close (11)

             end if

             if (phi_old-phi_new.ge.0._dp.and.phi_old-phi_new.lt.pi) then

                phi_scatter = 2._dp*pi - phi_scatter

             end if

             if (2._dp*pi+phi_old-phi_new.ge.0._dp.and.2._dp*pi+phi_old-phi_new.lt.pi) then

                phi_scatter = 2._dp*pi - phi_scatter

             end if

             if (phi_scatter.lt.0._dp) phi_scatter = phi_scatter + 2._dp*pi

             if (abs(mu_scatter).lt.1._dp) then

                call polarization_rotation(mu_scatter, phi_scatter, stokes_peel_in, scatter, direction_photon, &
                     det, stokes_out, .true.)

             else
                
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 049"
                write (11,*) mu_scatter, phi_scatter
                write (11,*) stokes_peel_in
                write (11,*) direction_photon
                close (11)

                cell_error = .true.

             end if

          else

             open (11, file=trim(error_log), position="append")
             write (11,*) "error 045"
             write (11,*) direction_photon
             close (11)

          end if

          if (.not.cell_error) then

             x_im = y_photon*cos_det_phi - x_photon*sin_det_phi
             y_im = z_photon*sin_det_theta - y_photon*cos_det_theta*sin_det_phi - x_photon*cos_det_theta*cos_det_phi

             ix = int(nx*(x_im+x_max)/(2._dp*x_max))+1
             iy = int(ny*(y_im+y_max)/(2._dp*y_max))+1

             if (photon_weight*stokes_out(1).gt.0._dp.and.photon_weight*stokes_out(1).lt.1.e100_dp) then

                detector_thread(ix,iy,1,1,thread_id+1) = detector_thread(ix,iy,1,1,thread_id+1) + photon_weight*stokes_out(1)
                detector_thread(ix,iy,2,1,thread_id+1) = detector_thread(ix,iy,2,1,thread_id+1) - photon_weight*stokes_out(2)
                detector_thread(ix,iy,3,1,thread_id+1) = detector_thread(ix,iy,3,1,thread_id+1) + photon_weight*stokes_out(3)
                detector_thread(ix,iy,4,1,thread_id+1) = detector_thread(ix,iy,4,1,thread_id+1) + photon_weight*stokes_out(4)

                detector_thread(ix,iy,1,2,thread_id+1) = &
                     detector_thread(ix,iy,1,2,thread_id+1) + (photon_weight*stokes_out(1))**2
                detector_thread(ix,iy,2,2,thread_id+1) = &
                     detector_thread(ix,iy,2,2,thread_id+1) + (photon_weight*stokes_out(2))**2
                detector_thread(ix,iy,3,2,thread_id+1) = &
                     detector_thread(ix,iy,3,2,thread_id+1) + (photon_weight*stokes_out(3))**2
                detector_thread(ix,iy,4,2,thread_id+1) = &
                     detector_thread(ix,iy,4,2,thread_id+1) + (photon_weight*stokes_out(4))**2

                detector_thread(ix,iy,1,3,thread_id+1) = detector_thread(ix,iy,1,3,thread_id+1) + 1._dp
                detector_thread(ix,iy,2,3,thread_id+1) = detector_thread(ix,iy,2,3,thread_id+1) + 1._dp
                detector_thread(ix,iy,3,3,thread_id+1) = detector_thread(ix,iy,3,3,thread_id+1) + 1._dp
                detector_thread(ix,iy,4,3,thread_id+1) = detector_thread(ix,iy,4,3,thread_id+1) + 1._dp
             
             else

                open (11, file=trim(error_log), position="append")
                write (11,*) "error 053"
                write (11,*) photon_weight
                write (11,*) stokes_out
                close (11)

             end if

          end if

       end if

    end if

  end subroutine peel_photon

  subroutine add_flow_global(thread_id, x, y, z, direction, energy, face_distance, cell)

    integer,  intent(in) :: thread_id, cell(3)
    real(dp), intent(in) :: x, y, z, face_distance, direction(3), energy
    real(dp)             :: theta, phi, r_dir, theta_dir, phi_dir

    theta = acos(z/sqrt(x*x+y*y+z*z))
    phi = atan2(y,x)

    r_dir = sin(theta)*cos(phi)*direction(1) + sin(theta)*sin(phi)*direction(2) + cos(theta)*direction(3)
    theta_dir = cos(theta)*cos(phi)*direction(1) + cos(theta)*sin(phi)*direction(2) - sin(theta)*direction(3)
    phi_dir = -sin(phi)*direction(1) + cos(phi)*direction(2)

    cell_flow_global(thread_id+1,1,cell(1),cell(2),cell(3)) = &
         cell_flow_global(thread_id+1,1,cell(1),cell(2),cell(3)) + r_dir*face_distance*energy

    cell_flow_global(thread_id+1,2,cell(1),cell(2),cell(3)) = &
         cell_flow_global(thread_id+1,2,cell(1),cell(2),cell(3)) + theta_dir*face_distance*energy

    cell_flow_global(thread_id+1,3,cell(1),cell(2),cell(3)) = &
         cell_flow_global(thread_id+1,3,cell(1),cell(2),cell(3)) + phi_dir*face_distance*energy

  end subroutine add_flow_global

  subroutine add_flow(thread_id, direction, energy, cell)

    integer,  intent(in) :: thread_id, direction, cell(3)
    real(dp), intent(in) :: energy

    if (direction.eq.1) then

       ! Energy flowing upward

       cell_flow(thread_id+1,1,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,1,cell(1),cell(2),cell(3)) + energy

    else if (direction.eq.2) then

       ! Energy flowing downward

       cell_flow(thread_id+1,2,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,2,cell(1),cell(2),cell(3)) + energy

    else if (direction.eq.3) then

       ! Energy flowing towards the south
       
       cell_flow(thread_id+1,3,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,3,cell(1),cell(2),cell(3)) + energy

    else if (direction.eq.4) then

       ! Energy flowing towards the north

       cell_flow(thread_id+1,4,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,4,cell(1),cell(2),cell(3)) + energy

    end if

  end subroutine add_flow
  
end program artes
