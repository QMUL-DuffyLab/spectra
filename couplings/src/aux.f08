module aux
  use iso_fortran_env
  implicit none
  integer, parameter :: dp = REAL64
  character(3), dimension(6) :: pigment_list = (/"CLA", "CHL", "LUT",&
                                                 "NEX", "XAT", "A86"/)

  type pigment
    character(len = 3) :: code
    real(dp) :: D, energy, reorg, lifetime
    real(dp), dimension(3) :: mu
    real(dp), dimension(:,:), allocatable :: coords
  end type pigment

  contains

    function cross(a, b)
      real(dp), dimension(3) :: cross
      real(dp), dimension(3), intent(in) :: a, b

      cross(1) = (a(2) * b(3)) - (a(3) * b(2))
      cross(2) = (a(3) * b(1)) - (a(1) * b(3))
      cross(3) = (a(1) * b(2)) - (a(2) * b(1))

    end function cross

    function CD(mu, r, eig, wi, gn)
      implicit none
      real(dp), dimension(:),    intent(in) :: wi
      real(dp), dimension(:, :), intent(in) :: mu, r, eig, gn
      real(dp), dimension(:), allocatable :: CD, gco
      real(dp) :: dd, w, pi
      integer :: Nt, Ns, n, m, k, step
      pi = 4 * atan(1.0_dp)
      Nt = size(mu, 2) 
      Ns = size(gn, 2)
      allocate(gco(Ns))
      do n = 1, Nt
        do m = 1, Nt
          dd = dot_product(cross(mu(n, :), mu(m, :)), r(n, :) - r(m, :))
          do k = 1, Nt
            do step = 1, Ns
              w = step * 2 * PI / Ns
              ! not right yet - need to think about the gn part
              gco(step) = eig(k, n) * eig(k, m) *&
                  (1 / sqrt(2 * pi)) * exp(-(w - wi(k))**2)
            end do
          end do
        end do
      end do
      deallocate(gco)

    end function CD
    
    function get_file_length(buffer) result(res)
      implicit none
      logical :: file_check
      character(len=*), intent(in) :: buffer
      character(99) :: line
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
        read(99, *, iostat=stat) line
        if (stat == iostat_end) then
          exit
        ! now check for the TER line; after that is connection info
        else if (index(line,"TER").ne.0) then
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

      s = 0.0_dp
      do i = 1, len1
        do j = 1, len2
          rx = p1(1, i) - p2(1, j)
          ry = p1(2, i) - p2(2, j)
          rz = p1(3, i) - p2(3, j)
          r = sqrt(rx**2 + ry**2 + rz**2)
          ! hack - so i can set homodimer parameters easier
          if (r.eq.0.0) then
            r = 1.0_dp
          end if
          s = s + (p1(4, i) * p2(4, j)) / r
        end do
      end do
      res = s
      
    end function J_calc

    ! subroutine because we want the dipole and the raw oscillator
    ! strength; we need to average those over each unique pigment
    subroutine mu_calc(p, len, mu, d_raw)
      implicit none
      integer, intent(in) :: len
      real(dp), dimension(4, len), intent(in) :: p
      real(dp), dimension(3), intent(out) :: mu
      real(dp), intent(out) :: d_raw
      integer :: i, j

      mu = 0.0_dp
      do i = 1, len
        do j = 1, 3
          mu(j) = mu(j) + (p(j, i) * p(4, i))
        end do
      end do

      d_raw = 0.0_dp
      do i = 1, 3
        d_raw = d_raw + mu(i)**2
      end do
      
    end subroutine mu_calc

    function angle(mu1, mu2) result(theta)
      implicit none
      real(dp), dimension(3), intent(in) :: mu1, mu2
      real(dp) :: theta, abs1, abs2, dot
      integer :: i

      abs1 = 0.0_dp
      abs2 = 0.0_dp
      do i = 1, 3
        abs1 = abs1 + mu1(i)**2
        abs2 = abs2 + mu2(i)**2
      end do

      dot = 0.0_dp
      do i = 1, 3
        dot = dot + (mu1(i) / sqrt(abs1)) * (mu2(i) / sqrt(abs2))
      end do
      ! for i = j dot can come out as 1.000000000002 or something;
      ! fix that
      if (dot.gt.1.0_dp.and.abs(dot - 1.0_dp).lt.1e-6) then
        dot = 1.0_dp
      end if
      theta = acos(dot)
      
    end function angle

    function k_factor(mu1, mu2, r1, r2) result(k12)
      implicit none
      real(dp), dimension(3), intent(in) :: mu1, mu2, r1, r2
      real(dp), dimension(3) :: rij
      real(dp) :: k12, abs1, abs2, absrij
      integer :: i

      ! normalise the dipoles and get r_ij
      abs1 = 0.0_dp
      abs2 = 0.0_dp
      absrij = 0.0_dp
      do i = 1, 3
        abs1 = abs1 + mu1(i)**2
        abs2 = abs2 + mu2(i)**2
        rij(i) = r2(i) - r1(i)
        absrij = absrij + rij(i)**2
      end do
      ! if i = j r_ij will be = 0; in that case don't do
      ! rij / sqrt(absrij) because 0.0 / 0.0 returns NAN
      if (absrij.gt.1E-10) then
        rij = rij / sqrt(absrij)
      end if

      k12 = 0.0_dp
      do i = 1, 3
        k12 = k12 + (mu1(i) / sqrt(abs1)) * (mu2(i) / sqrt(abs2))&
        - 3.0_dp * (mu1(i) * rij(i)) * (mu2(i) * rij(i))
      end do
      ! k12 = k12 * (sqrt(absrij))**3
      
    end function k_factor

    function get_pigment_code(gnt_files, i) result(code)
      implicit none
      character(100), dimension(:), intent(in) :: gnt_files
      integer, intent(in) :: i
      integer :: j, k
      character(3) :: code
      code = achar(0)

      ! check the list of pigments we're interested in
      do j = 1, size(pigment_list)
        k = index(gnt_files(i), pigment_list(j))
        if (k.eq.0) then
          cycle
        else
          code = pigment_list(j)
        end if
      end do

      if (index(code, achar(0)).ne.0) then
        write(*,*) "pigment code not found. details:"
        write(*,*) "code = ", code
        write(*,*) "gnt_file = ", gnt_files(i)
        write(*,*) "pigment_list = ", pigment_list
        stop
      end if

    end function get_pigment_code

    function get_unique_pigments(gnt_files) result(codes_uniq)
      implicit none
      character(100), dimension(:), intent(in) :: gnt_files
      integer :: i, j, num_uniq, gnt_length
      logical :: uniq
      character(3) :: code
      character(3), dimension(:), allocatable :: codes, uniq_temp, codes_uniq

      gnt_length = size(gnt_files)
      allocate(codes(gnt_length))
      allocate(uniq_temp(gnt_length))
      num_uniq = 0

      do i = 1, gnt_length
        code = get_pigment_code(gnt_files, i)
        codes(i) = code
        uniq = .true.

        if (i.gt.1) then
          do j = 1, i
            ! not the most efficient but should be fine.
            ! note that we can't assume the list is sorted
            if (uniq_temp(j).eq.code) then
              uniq = .false.
              exit
            end if
          end do
        end if

        if (uniq) then
          num_uniq = num_uniq + 1
          uniq_temp(num_uniq) = code
        end if

      end do

      allocate(codes_uniq(num_uniq))
      do i = 1, num_uniq
        codes_uniq(i) = uniq_temp(i)
      end do

    end function get_unique_pigments

    function count_pigments(gnt_files, unique_pigments) result(counts)
      implicit none
      character(100), dimension(:), intent(in) :: gnt_files
      character(3), dimension(:), intent(in) :: unique_pigments
      integer, dimension(:), allocatable :: counts
      integer :: i, n

      allocate(counts(size(unique_pigments)))
      counts = 0

      do i = 1, size(gnt_files)
        n = which_pigment(gnt_files, unique_pigments, i)
        counts(n) = counts(n) + 1
      end do

    end function count_pigments

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
        else
          cycle
        end if
      end do

      if (n.eq.0) then
        write(*,*) "which_pigment did not find match. Details:"
        write(*,*) "code = ", code
        stop
      end if

    end function which_pigment

    ! the raw partial charges won't give correct oscillator strengths,
    ! so we keep a running total of the raw numbers for each pigment
    ! type. here we normalise that total so we can then fix the average
    ! oscillator strength for each pigment to match the expected values
    function normalise_dipoles(raw_strengths, counts) result(norms)
      implicit none
      integer, dimension(:), intent(in) :: counts
      real(dp), dimension(:), intent(in) :: raw_strengths
      real(dp), dimension(:), allocatable :: norms
      integer :: i

      if (size(raw_strengths).ne.size(counts)) then
        write(*,*) "raw strengths array not the same length as counts"&
          " array. No idea how this happened?"
        stop
      end if

      allocate(norms(size(counts)))
      do i = 1, size(counts)
        norms(i) = raw_strengths(i) / float(counts(i))
      end do

    end function normalise_dipoles

    function pick_chls(gnt_files, unique_pigments, n)&
        result(chl_indices)
      implicit none
      ! real(dp), dimension(:,:), intent(in) :: Jij
      character(100), dimension(:) :: gnt_files
      character(3), dimension(:) :: unique_pigments
      integer :: i, n, n_chl, unique_index
      integer, dimension(:), allocatable :: chl_indices
      n_chl = 0

      do i = 1, n
        unique_index = which_pigment(gnt_files, unique_pigments, i)
        if (index(unique_pigments(unique_index), "CLA").ne.0) then
          n_chl = n_chl + 1
        else if (index(unique_pigments(unique_index), "CHL").ne.0) then
          n_chl = n_chl + 1
        end if
      end do

      allocate(chl_indices(n_chl))

      if (allocated(chl_indices)) then
        chl_indices = 0
        do i = 1, n
          unique_index = which_pigment(gnt_files, unique_pigments, i)
          if (index(unique_pigments(unique_index), "CLA").ne.0) then
            chl_indices(i) = i
          else if (index(unique_pigments(unique_index), "CHL").ne.0) then
            chl_indices(i) = i
          end if
        end do
      else
        write(*,*) "chl_incides not allocated"
      end if

    end function pick_chls

    ! note - this currently modifies the Jij passed to it.
    ! alternatively could give it another argument Jij_new
    ! and keep the pre-Redfield couplings as well.
    subroutine block_diag(Jij, indices, bloc, eigvals)
      implicit none
      real(dp), dimension(:,:) :: Jij
      integer, dimension(:) :: indices
      real(dp), dimension(:,:), allocatable :: bloc, redfield_jij
      real(dp) :: r
      real(dp), dimension(:), allocatable :: work, eigvals
      integer :: i, j, n, block_size, info
      n = floor(sqrt(size(Jij)+0.001)) ! prob don't need the 0.001 but
      block_size = size(indices)
      allocate(redfield_jij(n, n))
      ! allocate(bloc(block_size, block_size))

      if (allocated(bloc)) then
        do i = lbound(bloc, 1), ubound(bloc, 1)
          do j = lbound(bloc, 2), ubound(bloc, 2)
            bloc(i, j) = Jij(indices(i), indices(j))
          end do
        end do
        
        call dsyev('V', 'U', block_size, bloc, block_size,&
                    eigvals, r, -1, info)
        allocate(work(int(r)))
        call dsyev('V', 'U', block_size, bloc, block_size,&
                    eigvals, work, int(r), info)
        if (info.ne.0) then
          write(*, *) "diagonalisation of chl block failed."
          write(*, *) "Info = ", info
          do i = 1, block_size
            do j = 1, block_size
              bloc(i, j) = 0.0_dp
            end do
          end do
        end if

        ! note - we could modify the couplings here. might be more
        ! consistent to do do!!
        ! do i = 1, block_size
        !   do j = 1, block_size
        !     elem = 0.0_dp
        !     do k = 1, block_size
        !       elem = elem + bloc(i, k)**2 * Jij(indices(k), indices(j))
        !     end do
        !     Jij(indices(i), indices(j)) = elem
        !   end do
        ! end do
      else
        write (*,*) "block array not allocated for block_diag"
      end if
      deallocate(work)

    end subroutine block_diag

end module aux
