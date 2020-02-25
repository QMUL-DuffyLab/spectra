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

    subroutine get_parameters(filename) bind(C, name="get_parameters")
      use, intrinsic :: iso_c_binding
      import :: Params_f
      character(len=1, kind=C_CHAR), intent(in) :: filename(*)
      ! type(C_PTR) :: filename
    end subroutine get_parameters

  end interface

  contains

  function get_ligand_code(filename) result(res)
    character(len=*), intent(in) :: filename
    character(3) :: res
    integer :: i
    i = index(filename, "_CSV")
    res = filename(i - 6 : i - 4)
  end function get_ligand_code

  ! taken from http://fortranwiki.org/fortran/show/Generating+C+Interfaces
  ! apparently necessary because doing str // C_NULL_CHAR is bugged
  PURE FUNCTION F_C_STRING_FUNC (F_STRING) RESULT (C_STRING)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: F_STRING
    CHARACTER(LEN=1,KIND=C_CHAR) :: C_STRING(LEN_TRIM(F_STRING)+1)
    INTEGER                      :: N, I

    N = LEN_TRIM(F_STRING)
    DO I = 1, N
      C_STRING(I) = F_STRING(I:I)
    END DO
    C_STRING(N + 1) = C_NULL_CHAR

  END FUNCTION F_C_STRING_FUNC

end module c_interface

program test_c_interface
  use, intrinsic :: iso_c_binding
  use c_interface
  implicit none
  character(200) :: test_csv
  character(3) :: code
  type(Params_f) :: params
  character(len=20) :: test_str
  character(kind=C_CHAR), target :: c_str(21)
  type(c_ptr) :: str_ptr

  test_csv = "/Users/cgray/code/couplings/FCP/CLA401_CSV/frame1.csv"
  code = get_ligand_code(trim(adjustl(test_csv)))
  test_str = "lineshape/in/"//code//".def"
  c_str = F_C_STRING_FUNC(test_str)
  write (*,*) "Ligand code is ", code
  write (*,*) "test_str = ", test_str
  write (*,*) "c_str = ", c_str
  call get_parameters(c_str)

  write (*,*) "Params l0 = ",params%l0
  write (*,*) "Params g0 = ",params%g0
  write (*,*) "Params s0 = ",params%s0

  stop

end program test_c_interface
