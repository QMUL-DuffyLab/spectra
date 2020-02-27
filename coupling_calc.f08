program coupling_calc
  use iso_fortran_env
  use c_interface
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: sp = REAL64
  logical :: verbose
  character(50) :: coord_fmt, ei_file,&
  lambda_file, gnt_file, lifetimes_file, g_i_count
  character(100) :: output_dir
  character(200) :: input_file, jij_file,&
    eigvecs_file, eigvals_file, mu_i_file, mu_n_file, lambda_i_file,&
    gamma_i_file, spectra_input_file
  character(100), dimension(:), allocatable :: coord_files,&
  gnt_files
  integer :: i, j, k, coord_stat, control_len, tau
  integer, dimension(:), allocatable :: coord_lengths
  real :: start_time, end_time
  real(sp) :: r, e2kb, kc
  real(sp), dimension(:), allocatable :: work,&
  eigvals, ei, lambda, lifetimes
  real(sp), dimension(:,:), allocatable :: coords_i, coords_j,&
  Jij, Jeig, mu, mu_ex
  complex(sp), dimension(:,:), allocatable :: gnt

  verbose = .true.
  call cpu_time(start_time)
  coord_fmt = '(E016.8 1X E016.8 1X E016.8 1X E016.8)'

  if (command_argument_count().ne.2) then
    write (*,*) "Wrong number of arguments. Try again."
  else
    call get_command_argument(1, input_file)
    call get_command_argument(2, output_dir)
  end if

  if (verbose) then
    write(*,*) "Input file = ", input_file
    write(*,*) "Output dir = ", output_dir
  end if

  jij_file      = trim(adjustl(output_dir)) // "/J_ij.out"
  eigvecs_file  = trim(adjustl(output_dir)) // "/eigvecs.out"
  eigvals_file  = trim(adjustl(output_dir)) // "/eigvals.out"
  mu_n_file     = trim(adjustl(output_dir)) // "/mu_site.out"
  mu_i_file     = trim(adjustl(output_dir)) // "/mu_exciton.out"
  lambda_i_file = trim(adjustl(output_dir)) // "/lambda_exciton.out"
  gamma_i_file  = trim(adjustl(output_dir)) // "/lifetimes_exciton.out"
  spectra_input_file = "spectra_input.dat"

  ! first number is e_c^2 / 1.98E-23 * 1E-10, for conversion
  ! the 1.98E-23 isn't kB, it's some conversion factor;
  ! came from Chris. should probably sit down and figure it out lol
  ! second is Coulomb constant 1 / 4_pi_e0
  ! the 0.5 is because e_r = 2 for proteins
  e2kb = 1.2955E-5
  kc = 8.988E9 * 0.5

 ! tau is the number of femtoseconds i calculated lineshapes for.
 ! don't like hardcoding but it shouldn't need to be dynamic, really
  tau = 2000

  ! these are the ones we read in from
  ei_file = "ei.txt"
  lambda_file = "lambda.txt"
  gnt_file = "gnt.txt"
  lifetimes_file = "lifetimes.txt"

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
  allocate(mu(3, control_len)) ! fortran is column major
  allocate(mu_ex(3, control_len)) ! fortran is column major
  allocate(eigvals(control_len))
  allocate(lambda(control_len))
  allocate(lifetimes(control_len))
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
  do i = 1, control_len

    read(11, *) ei(i)
    read(12, *) lambda(i)
    read(13, '(a)') gnt_files(i)
    read(14, *) lifetimes(i)
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
        mu(:,i) = mu_calc(coords_i, coord_lengths(i))
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

  write(*,*) Jij

  ! DSYEV does eigendecomposition of a real symmetric matrix
  ! this first one is the query: find optimal size of work array
  ! using lwork = -1 sets r to the optimal work size
  call dsyev('V', 'U', control_len, Jeig, control_len,&
              eigvals, r, -1, coord_stat)
  write(*,*) "Optimal size of work array = ", int(r)
  allocate(work(int(r)))
  call dsyev('V', 'U', control_len, Jeig, control_len,&
              eigvals, work, int(r), coord_stat)
  write(*,*) "LAPACK INFO (should be 0) = ", coord_stat

  write(*,*) Jeig

  mu_ex = matmul(mu, Jeig) ! mix transition dipole moments
  gnt = matmul(Jeig**4, gnt) ! mix lineshape functions
  lambda = matmul(Jeig**4, lambda) ! mix reorganisation energies
  lifetimes = matmul(Jeig, lifetimes) ! mix relaxation times

  ! do i = 1, control_len
  !   do j = 1, control_len
  !     wij = (eigvals(i) - lambda(i)) - (eigvals(j) - lambda(j))
  !     do n = 1, control_len
  !       ! this should be C_n(wij) but C_n?????
  !       ! should be possible to call the C functions from
  !       ! fortran but would require building the structs
  !       k(i, j) = Jeig(i, n)**2 * Jeig(j, n)**2 * wij
  !   end do
  ! end do

  open(unit=20, file=spectra_input_file)
  ! stuff to read into spectra.c
  write(20, *) control_len
  write(20, *) trim(adjustl(eigvecs_file))
  write(20, *) trim(adjustl(eigvals_file))
  write(20, *) trim(adjustl(mu_i_file))
  write(20, *) trim(adjustl(lambda_i_file))
  write(20, *) trim(adjustl(gamma_i_file))

  open(unit=10, file=jij_file)
  open(unit=11, file=eigvecs_file)
  open(unit=12, file=eigvals_file)
  open(unit=13, file=mu_n_file)
  open(unit=14, file=mu_i_file)
  open(unit=15, file=lambda_i_file)
  open(unit=16, file=gamma_i_file)

  do i = 1, control_len
    do j = 1, control_len
      ! can write these with implied do loops
      ! (Jij(i, j), j = 1, control_len)
      ! but it made the code slower?
      write(10, '(17F16.5)', advance='no') Jij(i, j)
      write(11, '(17F16.5)', advance='no') Jeig(i, j)
    end do
    write(10,*) ! blank line between rows!
    write(11,*)
    write(12, *) eigvals(i)
    write(13, *) mu(1, i), mu(2, i), mu(3, i)
    write(14, *) mu_ex(1, i), mu_ex(2, i), mu_ex(3, i)
    write(15, *) lambda(i)
    write(16, *) lifetimes(i)

    ! now write out all the g_i(tau)s
    write(unit=g_i_count,fmt='(I0.2)') i
    write (*,*) trim(adjustl(output_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat"
    open(unit=17, file=trim(adjustl(output_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat")
    write(20, *) trim(adjustl(output_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat"
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


  ! here would be the place to use iso_c_binding and call the GSL
  ! integration routines; might as well just do it now instead of
  ! writing something separate to read everything back in again

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
  write (*,*) "Time taken: ", end_time - start_time, " seconds."

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

  function parse_tresp_line(buffer) result(res)
    implicit none
    character(len=*), intent(in) :: buffer
    character(len=len(buffer)) :: line
    integer :: pos
    real(kind=8) :: res
    line = buffer
    ! for some reason read didn't like hard tab delimiters and
    ! kept breaking, so this is a bit of a hack: the decimal point
    ! is the only full stop on each line, so find that and work back
    pos = scan(line, '.')
    line = line(pos - 2:)
    read(line, *) res

  end function parse_tresp_line

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
        s = s + (p1(4, i) * p2(4, j)) / r
      end do
    end do
    res = s
    
  end function J_calc

  function mu_calc(p, len) result(res)
    implicit none
    integer, parameter :: sp = real64
    integer, intent(in) :: len
    integer :: i, j
    real(sp), dimension(4, len) :: p
    real(sp), dimension(3) :: mu, res

    mu = 0.0
    do i = 1, len
      do j = 1, 3
        mu(j) = p(j, i) * p(4, i)
      end do
    end do
    res = mu
    
  end function mu_calc

end program coupling_calc
