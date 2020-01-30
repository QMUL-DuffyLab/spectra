program coupling_calc
  use iso_fortran_env
  implicit none
  logical :: file_check
  character(3) :: tresp_temp
  character(31) :: pdb_temp
  character(50) :: coord_fmt, tresp_fmt, control_file
  character(200) :: line
  character(100), dimension(:), allocatable :: coord_files, tresp_files
  integer :: i, j, k, l, posit, coord_stat, tresp_stat, control_len
  integer, dimension(:), allocatable :: coord_lengths, tresp_lengths
  real(kind=8) :: rx, ry, rz, r, conv
  real(kind=8), dimension(:), allocatable :: tresp_i, tresp_j
  real(kind=8), dimension(:,:), allocatable :: coords_i, coords_j, Jij

  coord_fmt = '(A31 F8.3 F8.3 F8.3 A)'
  tresp_fmt = '(A 1X A 1X F10.6)'

  ! NB: a lot of fucking around here can be avoided if all the
  ! tresp files are the same length as the corresponding PDB
  ! files, but this doesn't seem to be the case? check with Chris

  ! control_file is the name of the control file
  control_file = "J_control.txt"
  open(unit=10, file=control_file)
  control_len = 0
  do
    read(10, *, iostat=coord_stat)
    if (coord_stat == iostat_end) then
      exit
    else
      control_len = control_len + 1
    end if
  end do
  close(10)

  write(*,*) "Control length = ", control_len

  allocate(coord_files(control_len))
  allocate(tresp_files(control_len))
  allocate(coord_lengths(control_len))
  allocate(tresp_lengths(control_len))
  allocate(Jij(control_len, control_len))

  coord_lengths = 0
  tresp_lengths = 0
  Jij = 0.0

  open(unit=10, file=control_file)
  do i = 1, control_len
    read(10, '(a)') line
    posit = scan(line, ' ')
    coord_files(i) = line(1:posit - 1)
    tresp_files(i) = line(posit:)
  end do
  close(10)

  do i = 1, control_len
    inquire(file=trim(adjustl(coord_files(i))),exist=file_check)
    if (file_check) then
      open(unit=10, file=trim(adjustl(coord_files(i))))
      write(*,*) "Opened file ", trim(adjustl(coord_files(i)))
    else
      write (*,*) "File ",trim(adjustl(coord_files(i))), &
        " doesn't exist. Check and try again."
      stop
    end if

    do
      read(10, '(a)', iostat=coord_stat) line
      if (coord_stat == iostat_end) then
        exit
      else
        coord_lengths(i) = coord_lengths(i) + 1
      end if
    end do
    close(10)

  end do

  do i = 1, control_len

    inquire(file=trim(adjustl(tresp_files(i))),exist=file_check)
    if (file_check) then
      open(unit=10, file=trim(adjustl(tresp_files(i))))
      write(*,*) "Opened file ", trim(adjustl(tresp_files(i)))
    else
      write (*,*) "File ",trim(adjustl(tresp_files(i))), &
        " doesn't exist. Check and try again."
      stop
    end if

    do
      read(10, '(a)', iostat=tresp_stat) line
      if (tresp_stat == iostat_end) then
        exit
      else
        tresp_lengths(i) = tresp_lengths(i) + 1
      end if
    end do
    close(10)

  end do

  do i = 1, control_len
    write(*, *) trim(coord_files(i)), coord_lengths(i)
    write(*, *) trim(tresp_files(i)), tresp_lengths(i)
  end do

  i_loop: do i = 1, control_len

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
        ! write(*, *) tresp_j(k)
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
      write(*, *) trim(adjustl(coord_files(i)))," ",&
                  trim(adjustl(tresp_files(i)))," ",&
                  trim(adjustl(coord_files(j)))," ",&
                  trim(adjustl(tresp_files(j)))," ",&
                  Jij(i, j) * conv

      deallocate(coords_j)
      deallocate(tresp_j)
      
    end do j_loop

    deallocate(coords_i)
    deallocate(tresp_i)

  end do i_loop

  contains

  function parse_tresp_line(buffer) result(res)
    implicit none
    character(len=*), intent(in) :: buffer
    character(len=len(buffer)) :: line
    integer :: pos
    real(kind=8) :: res
    line = buffer
    pos = scan(line, '.')
    line = line(pos - 2:)
    read(line, *) res

  end function parse_tresp_line

end program coupling_calc

