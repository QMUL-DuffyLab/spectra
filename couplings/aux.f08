module aux
  use iso_fortran_env
  implicit none
  integer, parameter :: dp = REAL64

  type pigment
    character(len = 3) :: code
    real(dp) :: D
  end type pigment

  contains
    
    function get_file_length(buffer) result(res)
      implicit none
      logical :: file_check
      character(len=*), intent(in) :: buffer
      integer :: n, res, stat

      inquire(file=trim(adjustl(buffer)),exist=file_check)
      if (file_check) then
        open(unit=99, file=trim(adjustl(buffer)))
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
      integer, parameter :: dp = real64
      integer, intent(in) :: len1, len2
      integer :: i, j
      real(dp), dimension(4, len1) :: p1
      real(dp), dimension(4, len2) :: p2
      real(dp) :: s, rx, ry, rz, r, res

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
      integer, intent(in) :: len
      real(dp), intent(in) :: osc
      real(dp), dimension(4, len), intent(in) :: p
      real(dp), dimension(3), intent(out) :: mu
      real(dp), intent(out) :: d_raw
      integer :: i, j

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

    function get_pigment_code(gnt_files, i) result(code)
      implicit none
      character(100), dimension(:), intent(in) :: gnt_files
      integer, intent(in) :: i
      integer :: j
      character(3) :: code

      ! assumes the pigment code is at the start of the filename,
      ! straight after the directory separator, otherwise we'd have to
      ! come up with a list of possible codes and compare to that
      j = index(gnt_files(i), "/")
      code = gnt_files(i)(j + 1 : j + 3)

    end function get_pigment_code

    function get_unique_pigments(gnt_files) result(codes_uniq)
      implicit none
      character(100), dimension(:), intent(in) :: gnt_files
      integer :: i, j, num_uniq, gnt_length
      logical :: uniq
      character(3), dimension(:), allocatable :: codes, codes_uniq

      gnt_length = size(gnt_files)
      write(*,*) gnt_length
      allocate(codes(gnt_length))
      num_uniq = 0

      do i = 1, gnt_length
        codes(i) = get_pigment_code(gnt_files, i)
        uniq = .true.

        do j = 1, i
          ! not the most efficient but should be fine.
          ! note that we can't assume the list is sorted
          if (codes(i).eq.codes(j)) then
            uniq = .false.
            exit
          end if
        end do

        if (uniq) then
          codes_uniq(num_uniq) = codes(i)
          num_uniq = num_uniq + 1
        end if

        write(*, '(a)') codes(i)
      end do

      write(*, *) num_uniq
      do i = 1, num_uniq
        write(*, '(a)') codes_uniq(i)
      end do

    end function get_unique_pigments

    function which_pigment(gnt_files, unique_pigments, i) result(n)
      ! returns the relevant index in unique_pigments
      implicit none
      character(100), dimension(:), intent(in) :: gnt_files
      character(3), dimension(:), intent(in) :: unique_pigments
      character(3) :: code
      integer, intent(in) :: i
      integer :: j, uniq, n
      uniq = size(unique_pigments)
      code = get_pigment_code(gnt_files, i)
      n = 0 ! fortran is 1-valued so we can use 0 as an error value

      do j = 1, size(unique_pigments)
        if (index(unique_pigments(j), code).ne.0) then
          n = j
        end if
      end do

      if (n.eq.0) then
        write(*,*) "which_pigment did not find match. Details:"
        write(*,*) "code = ", code
        stop
      end if

    end function which_pigment

end module aux
