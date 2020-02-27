module c_interface
  use, intrinsic :: iso_c_binding
  implicit none

  ! this has to correspond to the C struct params!
  ! otherwise I can't call the C spectral density
  ! functions to get rate constants for the
  ! fluorescence spectrum. Really I should've just
  ! written the whole thing in C, but I didn't, because I hate C.
  type, bind(C) :: Params_f
    real(kind=C_DOUBLE) :: s0, s1, s2, g0, g1, g2,&
      l0, l1, l2, w1, w2, ti, T
    type(C_FUNPTR) :: cw
    integer(kind=C_INT) :: ligand
    character(kind=C_CHAR), dimension(200) :: aw_file, gt_file,&
      fw_file, lambda_file
  end type Params_f

  interface

    real(kind=C_DOUBLE) function cw(w, par) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE) :: w
      TYPE(C_PTR), value :: par
    end function cw

    type(C_PTR) function fortran_wrapper(ligand) bind(C, name="fortran_wrapper")
      use, intrinsic :: iso_c_binding
      import :: Params_f
      type(Params_f) :: p
      integer(kind=C_INT) :: ligand
    end function fortran_wrapper
    
  end interface

  contains

  function get_ligand_code(filename) result(res)
    character(len=*), intent(in) :: filename
    character(3) :: code
    integer :: i, res
    i = index(filename, "_CSV")
    code = filename(i - 6 : i - 4)

    ! 
    ! lookup table:
    ! CLA (Chlorophyll a)  = 0
    ! CHL (Chlorophyll b)  = 1
    ! KC1 (Chlorophyll c1) = 2
    ! KC2 (Chlorophyll c2) = 3
    ! A86 (fucoxanthin)    = 4
    ! DD6 (diadinoxanthin) = 5
    ! LUT (lutein)         = 6
    ! 

    select case (code)
      case ('CLA')
        res = 0
      case('CHL')
        res = 1
      case ('KC1')
        res = 2
      case ('KC2')
        res = 3
      case('A86')
        res = 4
      case ('DD6')
        res = 5
      case ('LUT')
        res = 6
      case default
        res = -1
        write (*,*) "Invalid ligand code.", code, ". Try again"
        stop
      end select

  end function get_ligand_code

end module c_interface

program test_c_interface
  use, intrinsic :: iso_c_binding
  use c_interface
  implicit none
  character(200) :: test_csv
  integer(kind=C_INT) :: ligand
  type(C_PTR) :: cptr
  type(Params_f), pointer :: par

  test_csv = "/Users/cgray/code/couplings/FCP/CLA401_CSV/frame1.csv"
  ligand = get_ligand_code(trim(adjustl(test_csv)))
  write (*,*) "Ligand number is ", ligand
  cptr = fortran_wrapper(ligand)
  call C_F_POINTER(cptr, par)

  write (*,*) "Params l0 = ",par%l0
  write (*,*) "Params g0 = ",par%g0
  write (*,*) "Params s0 = ",par%s0

  stop

  ! idea: jut write out the list of files, read that in,
  ! read in the eigvecs array to C, then we can do all the spectra
  ! shit in C. Can just include lineshapes/src/functions.c etc.

end program test_c_interface
