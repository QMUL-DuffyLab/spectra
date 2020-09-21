program coupling_calc
  use iso_fortran_env
  implicit none
  integer, parameter :: sp = REAL64
  logical :: verbose
  character(100) :: coord_fmt, ei_file, dipoles_file,&
  lambda_file, gnt_file, lifetimes_file, g_i_count
  character(100) :: input_dir, temp_string, tau_string
  character(200) :: input_file, jij_file, pop_file,&
    eigvecs_file, eigvals_file, mu_i_file, mu_n_file, lambda_i_file,&
    gamma_i_file, spectra_input_file, aw_output_file, fw_output_file
  character(100), dimension(:), allocatable :: coord_files,&
  gnt_files
  integer :: i, j, k, coord_stat, control_len, tau
  integer, dimension(:), allocatable :: coord_lengths
  real :: start_time, end_time
  real(sp) :: r, e2kb, kc, temperature
  real(sp), dimension(:), allocatable :: work,&
  eigvals, ei, lambda, lifetimes, dipoles
  real(sp), dimension(:,:), allocatable :: coords_i, coords_j,&
  Jij, Jeig, mu, mu_ex
  complex(sp), dimension(:,:), allocatable :: gnt

  verbose = .true.
  call cpu_time(start_time)
  coord_fmt = '(E016.8 1X E016.8 1X E016.8 1X E016.8)'

  if (command_argument_count().ne.4) then
    write (*,*) "Wrong number of arguments. Try again."
  else
    call get_command_argument(1, input_file)
    call get_command_argument(2, input_dir)
    call get_command_argument(3, temp_string)
    call get_command_argument(4, tau_string)
    read(temp_string, '(F16.1)') temperature
    read(tau_string, '(I5)') tau
  end if

  if (verbose) then
    write(*,*) "Input file = ", input_file
    write(*,*) "Directory = ", input_dir
    write(*,*) "Temperature = ", temperature
    write(*,*) "Tau = ", tau
  end if

  jij_file        = trim(adjustl(input_dir)) // "/J_ij.out"
  eigvecs_file    = trim(adjustl(input_dir)) // "/eigvecs.out"
  eigvals_file    = trim(adjustl(input_dir)) // "/eigvals.out"
  mu_n_file       = trim(adjustl(input_dir)) // "/mu_site.out"
  mu_i_file       = trim(adjustl(input_dir)) // "/mu_exciton.out"
  lambda_i_file   = trim(adjustl(input_dir)) // "/lambda_exciton.out"
  gamma_i_file    = trim(adjustl(input_dir)) // "/lifetimes_exciton.out"
  aw_output_file  = trim(adjustl(input_dir)) // "/aw.dat"
  fw_output_file  = trim(adjustl(input_dir)) // "/fw.dat"
  pop_file        = trim(adjustl(input_dir)) // "/populations.dat"
  spectra_input_file = "in/input_spectra.dat"

  ! first number is e_c^2 / 1.98E-23 * 1E-10, for conversion
  ! the 1.98E-23 isn't kB, it's some conversion factor;
  ! came from Chris. should probably sit down and figure it out lol
  ! second is Coulomb constant 1 / 4_pi_e0
  ! the 0.5 is because e_r = 2 for proteins
  e2kb = 1.2955E-5
  kc = 8.988E9 * 0.5

  ! these are the ones we read in from
  ei_file         = trim(adjustl(input_dir)) // "/ei.txt"
  lambda_file     = trim(adjustl(input_dir)) // "/lambda.txt"
  gnt_file        = trim(adjustl(input_dir)) // "/gnt.txt"
  lifetimes_file  = trim(adjustl(input_dir)) // "/lifetimes.txt"
  dipoles_file    = trim(adjustl(input_dir)) // "/dipoles.txt"

  ! this way we automatically deal with varying numbers of pigments
  control_len = get_file_length(input_file)

  if (verbose) then
    write(*,*) "Number of pigments = ", control_len
  end if

  ! coord/tresp_length arrays allow for varying numbers of
  ! atoms and tresp charges (not necessarily the same thing?)
  ! per pigment; we need to know these numbers to allocate
  ! the coordinate and tresp arrays later on in the main loop
  allocate(coord_files(control_len))
  allocate(gnt_files(control_len))
  allocate(coord_lengths(control_len))
  allocate(Jij(control_len, control_len))
  allocate(Jeig(control_len, control_len))
  allocate(mu(control_len, 3))
  allocate(mu_ex(control_len, 3)) ! fortran is column major
  allocate(eigvals(control_len))
  allocate(lambda(control_len))
  allocate(lifetimes(control_len))
  allocate(dipoles(control_len))
  allocate(ei(control_len))
  allocate(gnt(control_len, tau))

  coord_lengths = 0
  Jij = 0.0
  Jeig = 0.0
  mu = 0.0
  mu_ex = 0.0
  eigvals = 0.0
  lambda = 0.0
  lifetimes = 0.0
  dipoles = 0.0
  ei = 0.0
  gnt = (0.0, 0.0)

  ! in principle we could probably do the reading in as
  ! part of the main loop below, but the time taken to
  ! read in the control file and parse the filenames is
  ! negligible and it reads much nicer this way, I think.
  open(unit=10, file=input_file)
  open(unit=11, file=ei_file)
  open(unit=12, file=lambda_file)
  open(unit=13, file=gnt_file)
  open(unit=14, file=lifetimes_file)
  open(unit=15, file=dipoles_file)
  do i = 1, control_len

    read(11, *) ei(i)
    read(12, *) lambda(i)
    read(13, '(a)') gnt_files(i)
    read(14, *) lifetimes(i)
    read(15, *) dipoles(i)
    read(10, '(a)') coord_files(i)

    coord_lengths(i) = get_file_length(coord_files(i))

    if (verbose) then
      write(*, *) trim(adjustl(coord_files(i))), coord_lengths(i)
    end if

  end do
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)

  ! we multiply the whole Jij matrix later on to convert
  ! to wavenumbers; the read in excitation energies are already
  ! in wavenumbers so do that conversion backwards here
  ei = ei / (e2kb * kc)

  i_loop: do i = 1, control_len

    ! number of atoms per pigment isn't known at
    ! compile time, so allocate at run time
    allocate(coords_i(4, coord_lengths(i)))
    open(unit=10, file=trim(adjustl(coord_files(i))))
    open(unit=12, file=trim(adjustl(gnt_files(i))))
    do k = 1, tau
      read (12, '(F18.10, 1X, F18.10, 1X, F18.10)') r, gnt(i, k)
    end do
    close(12)

    do k = 1, coord_lengths(i)
      read(10, fmt=coord_fmt) coords_i(1, k), coords_i(2, k),&
                              coords_i(3, k), coords_i(4, k)
    end do
    close(10)

    j_loop: do j = i, control_len

      allocate(coords_j(4, coord_lengths(j)))
      open(unit=10, file=trim(adjustl(coord_files(j))))
      do k = 1, coord_lengths(j)
        read(10, fmt=coord_fmt) coords_j(1, k), coords_j(2, k),&
                                coords_j(3, k), coords_j(4, k)
      end do
      close(10)

      ! diagonal elements of Jij are the excitation energies
      ! we also want to calculate transition dipole moments
      if (i.eq.j) then
        Jij(i, j) = ei(i)
        call mu_calc(coords_i, coord_lengths(i), \
                 dipoles(i), mu(i, :), d_raw)
        write(*,*) "i = ", i, "mu(i) = ", mu(i, :),&
        "mu^2 = ", sum(mu(i, :)**2), "e_i = ", Jij(i,j)
      else
        Jij(i, j) = J_calc(coords_i, coords_j,&
                    coord_lengths(i), coord_lengths(j))
      end if

      deallocate(coords_j)
      
    end do j_loop
    deallocate(coords_i)
  end do i_loop

  Jij = Jij * e2kb * kc
  Jeig = Jij

  if (verbose) then
    write(*,*) Jij
  end if

  ! DSYEV does eigendecomposition of a real symmetric matrix
  ! this first one is the query: find optimal size of work array
  ! using lwork = -1 sets r to the optimal work size
  call dsyev('V', 'U', control_len, Jeig, control_len,&
              eigvals, r, -1, coord_stat)
  if (verbose) then
    write(*,*) "Optimal size of work array = ", int(r)
  end if
  allocate(work(int(r)))
  call dsyev('V', 'U', control_len, Jeig, control_len,&
              eigvals, work, int(r), coord_stat)
  if (verbose) then
    write(*,*) "LAPACK INFO (should be 0) = ", coord_stat
  end if

  if (verbose) then
    ! write(*,*) mu
    write(*,*)
    write(*,*) Jeig
  end if

  lifetimes = 1.0 / lifetimes ! mix rates not lifetimes!!

  ! transpose Jeig because we want to multiply the columns of
  ! Jeig with our quantitites (columns are eigenvectors)
  Jeig = transpose(Jeig)
  mu_ex     = matmul(Jeig, mu) ! mix transition dipole moments
  gnt       = matmul(Jeig**4, gnt) ! mix lineshape functions
  lambda    = matmul(Jeig**4, lambda) ! mix reorganisation energies
  lifetimes = matmul(Jeig**2, lifetimes) ! mix relaxation times
  ! transpose back to ensure the C code calculates Redfield rates
  ! correctly (it goes down the columns as you'd expect)
  Jeig = transpose(Jeig)

  lifetimes = 1.0 / lifetimes ! get back lifetimes

  open(unit=20, file=spectra_input_file)
  ! stuff to read into spectra.c
  write(20, *) control_len
  write(20, *) tau
  write(20, *) temperature
  write(20, '(a)') adjustl(trim(adjustl(eigvecs_file)))
  write(20, '(a)') adjustl(trim(adjustl(eigvals_file)))
  write(20, '(a)') adjustl(trim(adjustl(mu_i_file)))
  write(20, '(a)') adjustl(trim(adjustl(lambda_i_file)))
  write(20, '(a)') adjustl(trim(adjustl(gamma_i_file)))
  write(20, '(a)') adjustl(trim(adjustl(aw_output_file)))
  write(20, '(a)') adjustl(trim(adjustl(fw_output_file)))
  write(20, '(a)') adjustl(trim(adjustl(pop_file)))

  open(unit=10, file=jij_file)
  open(unit=11, file=eigvecs_file)
  open(unit=12, file=eigvals_file)
  open(unit=13, file=mu_n_file)
  open(unit=14, file=mu_i_file)
  open(unit=15, file=lambda_i_file)
  open(unit=16, file=gamma_i_file)

  if (verbose) then
    write(*,'(a)') "site"
    do i = 1, control_len
      do j = 1, 3
        write(*, '(2I2, F18.10, 1X)') i, j, mu(i, j)
      end do
    end do
    write(*,*)
    write(*,'(a)') "exciton"
    do i = 1, control_len
      do j = 1, 3
        write(*, '(2I2, F18.10, 1X)') i, j, mu_ex(i, j)
      end do
    end do
  end if

  do i = 1, control_len
    do j = 1, control_len
      ! can write these with implied do loops
      ! (Jij(i, j), j = 1, control_len)
      ! but it made the code slower?
      if (i.gt.j) then 
        write(10, '(17F16.5)', advance='no') Jij(j, i)
      else
        write(10, '(17F16.5)', advance='no') Jij(i, j)
      end if
      write(11, '(17F16.5)', advance='no') Jeig(i, j)
    end do
    write(10,*) ! blank line between rows!
    write(11,*)
    write(13, '(F18.10, 1X, F18.10, 1X, F18.10)') mu(i, 1),&
      mu(i, 2), mu(i, 3)
    write(14, '(F18.10, 1X, F18.10, 1X, F18.10)') mu_ex(i, 1),&
      mu_ex(i, 2), mu_ex(i, 3)
    write(12, '(F18.10)') eigvals(i)
    write(15, '(F18.10)') lambda(i)
    write(16, '(F18.10)') lifetimes(i)

    ! now write out all the g_i(tau)s
    write(unit=g_i_count,fmt='(I0.2)') i
    open(unit=17, file=trim(adjustl(input_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat")
    write(20, '(a)') trim(adjustl(input_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat"
    do j = 1, tau
      write(17, '(F18.10, 1X, F18.10)') gnt(i, j)
    end do
    close(17)

  end do
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)

  deallocate(coord_files)
  deallocate(gnt_files)
  deallocate(coord_lengths)
  deallocate(Jij)
  deallocate(Jeig)
  deallocate(work)
  deallocate(eigvals)
  deallocate(ei)
  deallocate(mu)
  deallocate(mu_ex)
  deallocate(lambda)
  deallocate(lifetimes)
  deallocate(gnt)

  call cpu_time(end_time)
  if (verbose) then
    write (*,*) "Time taken: ", end_time - start_time, " seconds."
  end if

  stop

  contains

  function get_file_length(buffer) result(res)
    implicit none
    logical :: file_check
    character(len=*), intent(in) :: buffer
    integer :: n, res, stat

    inquire(file=trim(adjustl(buffer)),exist=file_check)
    if (file_check) then
      open(unit=99, file=trim(adjustl(buffer)))
      ! write(*,*) "Opened file ", trim(adjustl(buffer))
    else
      write (*,*) "File ",trim(adjustl(buffer)), &
        " doesn't exist. Check and try again."
      stop
    end if

    n = 0
    do
      read(99, *, iostat=stat)
      if (stat == iostat_end) then
        exit
      else
        n = n + 1
      end if
    end do
    close(99)
    res = n

  end function get_file_length

  function J_calc(p1, p2, len1, len2) result(res)
    implicit none
    integer, parameter :: sp = real64
    integer, intent(in) :: len1, len2
    integer :: i, j
    real(sp), dimension(4, len1) :: p1
    real(sp), dimension(4, len2) :: p2
    real(sp) :: s, rx, ry, rz, r, res

    s = 0.0
    do i = 1, len1
      do j = 1, len2
        rx = p1(1, i) - p2(1, j)
        ry = p1(2, i) - p2(2, j)
        rz = p1(3, i) - p2(3, j)
        r = sqrt(rx**2 + ry**2 + rz**2)
        ! hack - so i can set homodimer parameters easier
        if (r.eq.0.0) then
          r = 1.0
        end if
        s = s + (p1(4, i) * p2(4, j)) / r
      end do
    end do
    res = s
    
  end function J_calc

  subroutine mu_calc(p, len, osc, mu, d_raw)
    implicit none
    integer, parameter :: sp = real64
    integer, intent(in) :: len
    real(sp), intent(in) :: osc
    real(sp), dimension(4, len), intent(in) :: p
    real(sp), dimension(3), intent(out) :: mu
    real(sp), intent(out) :: d_raw
    integer :: i, j
    real(sp), dimension(3) :: mu

    mu = 0.0
    do i = 1, len
      do j = 1, 3
        mu(j) = mu(j) + (p(j, i) * p(4, i))
      end do
    end do

    d_raw = 0.0
    do i = 1, 3
      d_raw = d_raw + mu(i)**2
    end do
    mu = mu * osc / sqrt(d_raw)
    
  end subroutine mu_calc

end program coupling_calc
