program coupling_calc
  use iso_fortran_env
  implicit none
  logical :: verbose
  character(31) :: pdb_temp
  character(50) :: coord_fmt, control_file, ei_file,&
  lambda_file, gnt_file, g_i_count
  character(200) :: line
  character(100), dimension(:), allocatable :: coord_files, tresp_files,&
  gnt_files
  integer :: i, j, k, l, posit, coord_stat, control_len, tau
  integer, dimension(:), allocatable :: coord_lengths, tresp_lengths
  real :: start_time, end_time
  real(kind=8) :: rx, ry, rz, r, e2kb, kc
  real(kind=8), dimension(:), allocatable :: tresp_i, tresp_j, work,&
  eigvals, ei, lambda
  real(kind=8), dimension(:,:), allocatable :: coords_i, coords_j,&
  Jij, Jeig, mu, mu_ex
  complex(kind=16), dimension(:,:), allocatable :: gnt

  ! thinking about this for future use
  ! state not pigment bc e.g. qy has different params than qx
  type State
    character(3) :: lig
    real(kind=8) :: eigvals_n
    real(kind=8) :: gamma_n
    real(kind=8) :: mu_n
    complex(kind=8), dimension(:), allocatable :: g_n
  end type State

  verbose = .false.
  call cpu_time(start_time)
  coord_fmt = '(A31 F8.3 F8.3 F8.3 A)'

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

  ! this way we automatically deal with varying numbers of pigments
  control_file = "J_control.txt"
  ei_file = "ei.txt"
  lambda_file = "lambda.txt"
  gnt_file = "gnt.txt"
  control_len = get_file_length(control_file)

  if (verbose) then
    write(*,*) "Control length = ", control_len
  end if

  ! coord/tresp_length arrays allow for varying numbers of
  ! atoms and tresp charges (not necessarily the same thing?)
  ! per pigment; we need to know these numbers to allocate
  ! the coordinate and tresp arrays later on in the main loop
  allocate(coord_files(control_len))
  allocate(tresp_files(control_len))
  allocate(gnt_files(control_len))
  allocate(coord_lengths(control_len))
  allocate(tresp_lengths(control_len))
  allocate(Jij(control_len, control_len))
  allocate(Jeig(control_len, control_len))
  allocate(mu(3, control_len)) ! fortran is column major
  allocate(mu_ex(3, control_len))
  allocate(eigvals(control_len))
  allocate(lambda(control_len))
  allocate(ei(control_len))
  allocate(gnt(control_len, tau))

  coord_lengths = 0
  tresp_lengths = 0
  Jij = 0.0
  Jeig = 0.0
  mu = 0.0
  mu_ex = 0.0
  eigvals = 0.0
  lambda = 0.0
  ei = 0.0
  gnt = (0.0, 0.0)

  ! in principle we could probably do the reading in as
  ! part of the main loop below, but the time taken to
  ! read in the control file and parse the filenames is
  ! negligible and it reads much nicer this way, I think.
  open(unit=10, file=control_file)
  open(unit=11, file=ei_file)
  open(unit=12, file=lambda_file)
  open(unit=13, file=gnt_file)
  do i = 1, control_len

    read(11, *) ei(i)
    read(12, *) lambda(i)
    read(13, '(a)') gnt_files(i)
    read(10, '(a)') line
    posit = scan(line, ' ')
    coord_files(i) = line(1:posit - 1)
    tresp_files(i) = line(posit:)

    coord_lengths(i) = get_file_length(coord_files(i))
    tresp_lengths(i) = get_file_length(tresp_files(i))

    if (verbose) then
      write(*, *) trim(adjustl(coord_files(i))), coord_lengths(i)
      write(*, *) trim(adjustl(tresp_files(i))), tresp_lengths(i)
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

  ! we can also speed this up by only taking the j loop from
  ! i, control_len instead of 1, control_len, since the
  ! matrix J_ij is symmetric by construction. It'll require an
  ! extra few lines when writing output which I haven't bothered
  ! to write yet, hence why I'm still doing it the slow way here.
  i_loop: do i = 1, control_len

    ! number of atoms/tresp charges per pigment isn't known at
    ! compile time, so allocate at run time once we've got them
    allocate(coords_i(3, coord_lengths(i)))
    allocate(tresp_i(tresp_lengths(i)))
    open(unit=10, file=trim(adjustl(coord_files(i))))
    open(unit=11, file=trim(adjustl(tresp_files(i))))
    open(unit=12, file=trim(adjustl(gnt_files(i))))
    do k = 1, tau
      read (12, '(F18.10, 1X, F18.10, 1X, F18.10)') r, gnt(i, k)
    end do
    close(12)

    read(11, *) ! first line in tresp files is a comment
    do k = 1, coord_lengths(i)
      read(10, fmt=coord_fmt) pdb_temp,&
      coords_i(1, k), coords_i(2, k), coords_i(3, k), pdb_temp
    end do
    do k = 1, tresp_lengths(i) - 1
      read(11, '(a)') line
      ! had to write a separate function to parse the tresp input
      ! because fortran was being funny about the hard tabs
      tresp_i(k) = parse_tresp_line(line)
    end do
    close(10)
    close(11)

    j_loop: do j = 1, control_len

      allocate(coords_j(3, coord_lengths(j)))
      allocate(tresp_j(tresp_lengths(j)))
      open(unit=10, file=trim(adjustl(coord_files(j))))
      open(unit=11, file=trim(adjustl(tresp_files(j))))
      read(11, *)
      do k = 1, coord_lengths(j)
        read(10, fmt=coord_fmt) pdb_temp, &
        coords_j(1, k), coords_j(2, k), coords_j(3, k), pdb_temp
      end do
      do k = 1, tresp_lengths(j) - 1
        read(11, '(a)') line
        tresp_j(k) = parse_tresp_line(line)
      end do
      close(10)
      close(11)

      ! only need to do this once; doesn't have to be in kl loops
      ! diagonal elements are the excitation energies
      if (i.eq.j) then
        Jij(i, j) = ei(i)
      end if

      k_loop: do k = 1, coord_lengths(i)
        l_loop: do l = 1, coord_lengths(j)

          if (i.ne.j) then
            rx = coords_i(1, k) - coords_j(1, l)
            ry = coords_i(2, k) - coords_j(2, l)
            rz = coords_i(3, k) - coords_j(3, l)
            r = sqrt(rx**2 + ry**2 + rz**2)
            Jij(i, j) = Jij(i, j) + ((tresp_i(k) * tresp_j(l)) / (r))
          else
            ! calculate transition dipole moments for each pigment
            mu(1, i) = mu(1, i) + tresp_i(k) * coords_i(1, k)
            mu(2, i) = mu(2, i) + tresp_i(k) * coords_i(2, k)
            mu(3, i) = mu(3, i) + tresp_i(k) * coords_i(3, k)
          end if

        end do l_loop
      end do k_loop

      deallocate(coords_j)
      deallocate(tresp_j)
      
    end do j_loop

    deallocate(coords_i)
    deallocate(tresp_i)

  end do i_loop

  ! not sure if this is correct or not
  Jij = Jij * e2kb * kc
  Jeig = Jij
  mu_ex = matmul(mu, Jeig) ! mix transition dipole moments
  gnt = matmul(Jeig, gnt) ! mix lineshape functions
  lambda = matmul(Jeig, lambda) ! mix lineshape functions

  ! DSYEV does eigendecomposition of a real symmetric matrix
  ! this first one is the query: find optimal size of work array
  ! using lwork = -1 sets r to the optimal work size
  call dsyev('V', 'U', control_len, Jeig, control_len,&
              eigvals, r, -1, coord_stat)
  write(*,*) "Optimal size of work array = ", int(r)
  allocate(work(int(r)))
  call dsyev('V', 'U', control_len, Jeig, control_len,&
              eigvals, work, int(r), coord_stat)
  write(*,*) "LAPACK INFO (should be 0) = ",coord_stat
  write(*, *) Jeig(1,1), Jeig(1,2)

  open(unit=10, file="out/J_ij.out")
  open(unit=11, file="out/J_eig.out")
  open(unit=12, file="out/eigvals.out")
  open(unit=13, file="out/mu_site.out")
  open(unit=14, file="out/mu_exciton.out")
  open(unit=15, file="out/lambda_exciton.out")
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

    ! now write out all the g_i(tau)s
    write(unit=g_i_count,fmt='(I0.2)') i
    open(unit=16, file="out/g_i_" // trim(adjustl(g_i_count)) // ".dat")
    do j = 1, tau
      write(16, '(F18.10, 1X, F18.10)') gnt(i, j)
    end do
    close(16)

  end do
  close(10)
  close(11)
  close(12)

  deallocate(coord_files)
  deallocate(tresp_files)
  deallocate(gnt_files)
  deallocate(coord_lengths)
  deallocate(tresp_lengths)
  deallocate(Jij)
  deallocate(Jeig)
  deallocate(work)
  deallocate(eigvals)
  deallocate(ei)
  deallocate(mu)
  deallocate(mu_ex)
  deallocate(lambda)
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

end program coupling_calc
