program coupling_calc
  use iso_fortran_env
  implicit none
  logical :: file_check, verbose
  character(31) :: pdb_temp
  character(50) :: coord_fmt, control_file
  character(200) :: line
  character(100), dimension(:), allocatable :: coord_files, tresp_files
  integer :: i, j, k, l, posit, coord_stat, tresp_stat, control_len
  integer, dimension(:), allocatable :: coord_lengths, tresp_lengths
  real :: start_time, end_time
  real(kind=8) :: rx, ry, rz, r, conv
  real(kind=8), dimension(:), allocatable :: tresp_i, tresp_j, work, lambda
  real(kind=8), dimension(:,:), allocatable :: coords_i, coords_j, Jij, Jeig

  verbose = .false.
  call cpu_time(start_time)
  coord_fmt = '(A31 F8.3 F8.3 F8.3 A)'

  ! NB: some of the fucking around here could probably be avoided
  ! if all the tresp files are the same length as the corresponding PDB
  ! files, but this doesn't seem to be the case? i have no idea why

  ! control_file is the name of the control file
  ! this way we automatically deal with varying numbers of pigments
  control_file = "J_control.txt"
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
  allocate(coord_lengths(control_len))
  allocate(tresp_lengths(control_len))
  allocate(Jij(control_len, control_len))
  allocate(Jeig(control_len, control_len))
  allocate(lambda(control_len))

  coord_lengths = 0
  tresp_lengths = 0
  Jij = 0.0

  ! in principle we could probably do the reading in as
  ! part of the main loop below, but the time taken to
  ! read in the control file and parse the filenames is
  ! negligible and it reads much nicer this way, I think.
  open(unit=10, file=control_file)
  do i = 1, control_len

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
      read(11, *) ! first line in tresp files is a comment
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

      k_loop: do k = 1, coord_lengths(i)
        l_loop: do l = 1, coord_lengths(j)

          if (i.ne.j) then
            rx = coords_i(1, k) - coords_j(1, l)
            ry = coords_i(2, k) - coords_j(2, l)
            rz = coords_i(3, k) - coords_j(3, l)
            r = sqrt(rx**2 + ry**2 + rz**2)
            ! again, tresp length must equal pdb length, surely?????
            Jij(i, j) = Jij(i, j) + ((tresp_i(k) * tresp_j(l)) / (r))
          else
            ! not sure this is correct yet but it should be
            ! possible to get oscillator strengths this way
            rx = coords_i(1, k)
            ry = coords_i(2, k)
            rz = coords_i(3, k)
            r = sqrt(rx**2 + ry**2 + rz**2)
            Jij(i, j) = Jij(i, j) + tresp_i(k) * tresp_i(k) * r
          end if

        end do l_loop
      end do k_loop

      ! first number is e_c^2 / k_B * 1E-10, for conversion
      ! although NB: this seems to use KB = 1.98E-23, it
      ! should be 1.38E23. check this with Chris!
      ! second is Coulomb constant 1 / 4_pi_e0
      ! third is e_r = 2 for proteins
      conv = 1.295E-5 * 8.988E9 * 0.5
      ! write(*, *) trim(adjustl(coord_files(i)))," ",&
      !             trim(adjustl(tresp_files(i)))," ",&
      !             trim(adjustl(coord_files(j)))," ",&
      !             trim(adjustl(tresp_files(j)))," ",&
      !             Jij(i, j) * conv

      deallocate(coords_j)
      deallocate(tresp_j)
      
    end do j_loop

    deallocate(coords_i)
    deallocate(tresp_i)

  end do i_loop

  call cpu_time(end_time)
  write (*,*) "Time taken: ", end_time - start_time, " seconds."

  Jij = Jij * conv
  Jeig = Jij
  ! SSYEV does eigendecomposition of a real symmetric matrix
  ! this first one is the query: find optimal size of work array
  ! using lwork = -1 sets r to the optimal work size
  call ssyev('V', 'U', control_len, Jeig, control_len, lambda,&
              r, -1, coord_stat)
  write(*,*) "Optimal size of work array = ", int(r)
  allocate(work(int(r)))
  call ssyev('V', 'U', control_len, Jeig, control_len, lambda,&
              work, 4 * control_len, coord_stat)
  ! now Jeig and lambda contain eigenvectors and eigenvalues of J

  open(unit=10, file="out/J_ij.out")
  open(unit=11, file="out/J_eig.out")
  open(unit=12, file="out/lambda.out")
  do i = 1, control_len
    do j = 1, control_len
      ! can write these with implied do loops
      ! (Jij(i, j), j = 1, control_len)
      ! but it made the code slower?
      write(10, '(17F18.10)', advance='no') Jij(i, j)
      write(11, '(17F18.10)', advance='no') Jeig(i, j)
    end do
    write(10,*) ! blank line between rows!
    write(11,*)
    write(12, '(F18.10)') lambda(i)
  end do
  close(10)
  close(11)
  close(12)

  deallocate(coord_files)
  deallocate(tresp_files)
  deallocate(coord_lengths)
  deallocate(tresp_lengths)
  deallocate(Jij)
  deallocate(Jeig)
  deallocate(lambda)

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
