program coupling_calc
  use iso_fortran_env
  use aux
  implicit none
  integer, parameter :: cdp = REAL128
  logical :: verbose, print_jij, print_jij_diag
  character(100) :: coord_fmt, ei_file, dipoles_file,&
  lambda_file, gnt_file, lifetimes_file, g_i_count
  character(100) :: input_dir, temp_string, tau_string
  character(200) :: input_file, jij_file, pop_file,&
    eigvecs_file, eigvals_file, mu_i_file, mu_n_file, lambda_i_file,&
    gamma_i_file, spectra_input_file, aw_output_file, fw_output_file,&
    kappa_file, theta_file, com_file
  character(3), dimension(:), allocatable :: unique_pigments
  character(5) :: atom_code, pigment_code
  character(100), dimension(:), allocatable :: coord_files,&
    gnt_files
  integer :: i, j, k, control_len, tau,&
    num_unique_pigments, unique_index, myunit
  integer, dimension(:), allocatable :: coord_lengths, pigment_counts,&
    block_indices
  real :: start_time, end_time
  real(dp) :: r, e2kb, kc, temperature, d_raw, ratio
  real(dp), dimension(:), allocatable :: eigvals, ei, lambda,&
    lifetimes, dipoles, raw_osc, norm_osc, osc_check, block_eigvals
  real(dp), dimension(:,:), allocatable :: coords_i, coords_j,&
    Jij, Jeig, mu, mu_ex, r_charge, kappa, theta, bloc
  complex(cdp), dimension(:,:), allocatable :: gnt

  verbose = .true.
  print_jij = .false.
  print_jij_diag = .false.
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
  kappa_file  = trim(adjustl(input_dir)) // "/kappa.dat"
  theta_file  = trim(adjustl(input_dir)) // "/theta.dat"
  com_file  = trim(adjustl(input_dir)) // "/c_o_m.dat"

  ! first number is e_c^2 / 1.98E-23 * 1E-10, for conversion
  ! the 1.98E-23 isn't kB, it's some conversion factor;
  ! came from Chris. should probably sit down and figure it out lol
  ! second is Coulomb constant 1 / 4_pi_e0
  ! the 0.5 is because e_r = 2 for proteins
  e2kb = 1.2955E-5_dp
  kc = 8.988E9_dp * 0.5_dp

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
  allocate(kappa(control_len, control_len))
  allocate(theta(control_len, control_len))
  allocate(mu(control_len, 3))
  allocate(mu_ex(control_len, 3)) ! fortran is column major
  allocate(r_charge(control_len, 3))
  allocate(eigvals(control_len))
  allocate(lambda(control_len))
  allocate(lifetimes(control_len))
  allocate(dipoles(control_len))
  allocate(ei(control_len))
  allocate(gnt(control_len, tau))

  coord_lengths = 0
  Jij       = 0.0_dp
  Jeig      = 0.0_dp
  mu        = 0.0_dp
  mu_ex     = 0.0_dp
  eigvals   = 0.0_dp
  lambda    = 0.0_dp
  lifetimes = 0.0_dp
  dipoles   = 0.0_dp
  ei        = 0.0_dp
  gnt       = (0.0_dp, 0.0_dp)

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

  ! we want the list of unique pigments and how many of each there are
  unique_pigments = get_unique_pigments(gnt_files)
  pigment_counts = count_pigments(gnt_files, unique_pigments)
  num_unique_pigments = size(pigment_counts)
  if (verbose) then
    write(*, '(a)') "Unique pigments: ", unique_pigments
    write(*, *) "Pigment counts: ", pigment_counts
    write(*, *) "Number of unique pigments: ", num_unique_pigments
  end if

  allocate(raw_osc(num_unique_pigments))
  raw_osc = 0.0_dp

  ! we multiply the whole Jij matrix later on to convert
  ! to wavenumbers; the read in excitation energies are already
  ! in wavenumbers so do that conversion backwards here
  ei = ei / (e2kb * kc)

  open(newunit=myunit, file=trim(adjustl(com_file)))

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
    pigment_code = get_pigment_code(gnt_files, i)

    do k = 1, coord_lengths(i)
      read(10, fmt='(A5, 1X)', advance='no') atom_code
      read(10, fmt=coord_fmt) coords_i(1, k), coords_i(2, k),&
                              coords_i(3, k), coords_i(4, k)
      if (trim(adjustl(pigment_code)).eq."CHL") then
        if (trim(adjustl(atom_code)).eq."MG") then
          r_charge(i, :) = coords_i(1:3, k)
        end if
      else if (trim(adjustl(pigment_code)).eq."CLA") then
        if (trim(adjustl(atom_code)).eq."MG") then
          r_charge(i, :) = coords_i(1:3, k)
        end if
      else if (trim(adjustl(pigment_code)).eq."LUT") then
        if (trim(adjustl(atom_code)).eq."C20") then
          r_charge(i, :) = coords_i(1:3, k)
        end if
      end if
    end do
    close(10)
    write(myunit, '(F18.10, 1X, F18.10, 1X, F18.10)') r_charge(i, :)

    j_loop: do j = i, control_len

      allocate(coords_j(4, coord_lengths(j)))
      open(unit=10, file=trim(adjustl(coord_files(j))))
      do k = 1, coord_lengths(j)
        read(10, fmt='(A5, 1X)', advance='no') atom_code
        read(10, fmt=coord_fmt) coords_j(1, k), coords_j(2, k),&
                                coords_j(3, k), coords_j(4, k)
      end do
      close(10)
      
      ! diagonal elements of Jij are the excitation energies
      ! we also want to calculate transition dipole moments
      if (i.eq.j) then
        Jij(i, j) = ei(i)
        unique_index = which_pigment(gnt_files, unique_pigments, i)
        call mu_calc(coords_i, coord_lengths(i),  mu(i, :), d_raw)
        ! keep running total of the raw oscillator strengths by pigment
        raw_osc(unique_index) = raw_osc(unique_index) + d_raw
      else
        Jij(i, j) = J_calc(coords_i, coords_j,&
                    coord_lengths(i), coord_lengths(j))
      end if
      ! write(*, *) i, j, Jij(i, j)

      ! NB: after discussion with Chris - couplings less than 1 cm^-1
      ! are too small to worry about
      if (abs(Jij(i, j) * e2kb * kc).lt.1.0_dp) then
        Jij(i, j) = 0.0_dp
      end if

      deallocate(coords_j)
      
    end do j_loop
    deallocate(coords_i)
  end do i_loop
  close(myunit)

  do i = 1, control_len
    do j = 1, control_len
      if (i.eq.j) then
        kappa(i, j) = 1.0_dp
        theta(i, j) = 0.0_dp
      else
        theta(i, j) = angle(mu(i, :), mu(j, :))
        kappa(i, j) = k_factor(mu(i, :), mu(j, :), &
                      r_charge(i, :), r_charge(j, :))
      end if
      if (verbose) then
        write(*,*) i, j, kappa(i, j), theta(i, j)
      end if
    end do
  end do

  Jij = Jij * e2kb * kc

  ! can do away with the final part once this works
  block_indices = pick_chls(gnt_files, unique_pigments, control_len)
  ! alternatively, could go back to all Redfield approach by doing
  ! allocate(block_indices(control_len))
  ! forall(j = 1:control_len) block_indices(j) = j

  allocate(bloc(size(block_indices), size(block_indices)))
  allocate(block_eigvals(size(block_indices)))
  call block_diag(Jij, block_indices, bloc, block_eigvals)
  if (verbose) then
    do i = 1, size(block_indices)
      write(*, '(14F10.6)') bloc(i, :)
    end do
    do i = 1, size(block_indices)
      write(*, '(I4, F10.3)') i, block_eigvals(i)
    end do
  end if

  ! set Jeig to the identity
  Jeig = 0.0_dp
  eigvals = ei
  forall(j = 1:control_len) Jeig(j, j) = 1.0_dp
  ! now overwite the Redfield block with the diagonalised bit
  ! i think the block_indices(i) is correct??
  do i = 1, size(block_indices)
    eigvals(block_indices(i)) = block_eigvals(i)
    do j = 1, size(block_indices)
      Jeig(block_indices(i), j) = bloc(i, j)
    end do
  end do
  deallocate(bloc)
  deallocate(block_eigvals)
  deallocate(block_indices)

  norm_osc = normalise_dipoles(raw_osc, pigment_counts)
  allocate(osc_check(num_unique_pigments))
  osc_check = 0.0_dp
  do i = 1, control_len
    unique_index = which_pigment(gnt_files, unique_pigments, i)
    ! this is maybe a bit ugly, but it works:
    ! for each pigment, check which type it is and what the raw
    ! average oscillator strength was for that pigment type.
    ! then divide by the expected average oscillator strength to
    ! get a ratio and divide by the sqrt of that ratio. this ensures
    ! that the average oscillator strength for each pigment type will
    ! match the expected value (i.e. for LHCII, the average |mu^2|
    ! over all eight CLAs will equal whatever number I put in for the
    ! oscillator strength of CLA - currently 4.0 Debye
    ratio = norm_osc(unique_index) / dipoles(i)**2
    mu(i, :) = mu(i, :) / sqrt(ratio)
    osc_check(unique_index) = osc_check(unique_index) + sum(mu(i, :)**2)
  end do

  if (verbose) then
    do i = 1, num_unique_pigments
      write(*, *) "pigment type: ", unique_pigments(i)
      write(*, *) "count: ", pigment_counts(i)
      write(*, *) "<D>, <|D^2|> ",&
        sqrt(osc_check(i) / float(pigment_counts(i))),&
        osc_check(i) / float(pigment_counts(i))
    end do
  end if

  if (print_jij) then
    write(*,*) Jij
  end if

  ! DSYEV does eigendecomposition of a real symmetric matrix
  ! this first one is the query: find optimal size of work array
  ! using lwork = -1 sets r to the optimal work size
  ! call dsyev('V', 'U', control_len, Jeig, control_len,&
  !             eigvals, r, -1, coord_stat)
  ! if (verbose) then
  !   write(*,*) "Optimal size of work array = ", int(r)
  ! end if
  ! allocate(work(int(r)))
  ! call dsyev('V', 'U', control_len, Jeig, control_len,&
  !             eigvals, work, int(r), coord_stat)
  ! if (verbose) then
  !   write(*,*) "LAPACK INFO (should be 0) = ", coord_stat
  ! end if

  if (print_jij_diag) then
    ! write(*,*) mu
    write(*,*)
    write(*,*) Jeig
  end if

  lifetimes = 1.0_dp / lifetimes ! mix rates not lifetimes!!

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

  lifetimes = 1.0_dp / lifetimes ! get back lifetimes

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
  write(20, '(a)') adjustl(trim(adjustl(jij_file)))

  open(unit=10, file=jij_file)
  open(unit=11, file=eigvecs_file)
  open(unit=12, file=eigvals_file)
  open(unit=13, file=mu_n_file)
  open(unit=14, file=mu_i_file)
  open(unit=15, file=lambda_i_file)
  open(unit=16, file=gamma_i_file)
  open(unit=17, file=kappa_file)
  open(unit=18, file=theta_file)

  do i = 1, control_len
    do j = 1, control_len
      ! can write these with implied do loops
      ! (Jij(i, j), j = 1, control_len)
      ! but it made the code slower?
      if (i.gt.j) then 
        write(10, '(20F16.5)', advance='no') Jij(j, i)
      else
        write(10, '(20F16.5)', advance='no') Jij(i, j)
      end if
      write(11, '(20F16.5)', advance='no') Jeig(i, j)
      write(17, '(20F16.5)', advance='no') kappa(i, j)
      write(18, '(20F16.5)', advance='no') theta(i, j)
    end do
    write(10,*) ! blank line between rows!
    write(11,*)
    write(17,*)
    write(18,*)
    write(13, '(F18.10, 1X, F18.10, 1X, F18.10)') mu(i, 1),&
      mu(i, 2), mu(i, 3)
    write(14, '(F18.10, 1X, F18.10, 1X, F18.10)') mu_ex(i, 1),&
      mu_ex(i, 2), mu_ex(i, 3)
    write(12, '(F18.10)') eigvals(i)
    write(15, '(F18.10)') lambda(i)
    write(16, '(F18.10)') lifetimes(i)

    ! now write out all the g_i(tau)s
    write(unit=g_i_count,fmt='(I0.2)') i
    open(unit=19, file=trim(adjustl(input_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat")
    write(20, '(a)') trim(adjustl(input_dir)) // "/g_i_" // trim(adjustl(g_i_count)) // ".dat"
    do j = 1, tau
      write(19, '(F18.10, 1X, F18.10)') gnt(i, j)
    end do
    close(19)

  end do
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)

  deallocate(coord_files)
  deallocate(gnt_files)
  deallocate(coord_lengths)
  deallocate(Jij)
  deallocate(Jeig)
  deallocate(kappa)
  deallocate(theta)
  deallocate(eigvals)
  deallocate(ei)
  deallocate(mu)
  deallocate(mu_ex)
  deallocate(r_charge)
  deallocate(lambda)
  deallocate(lifetimes)
  deallocate(gnt)

  call cpu_time(end_time)
  if (verbose) then
    write (*,*) "Time taken: ", end_time - start_time, " seconds."
  end if

  stop

end program coupling_calc
