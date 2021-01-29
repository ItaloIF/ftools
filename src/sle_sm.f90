module sle_sm
    implicit none
    type sm
        integer(8) :: n
        integer(8) :: m
        integer(8), allocatable :: ia(:)   ! ia(n+1)
        integer(8), allocatable :: ja(:)   ! ja(m)
        real(8), allocatable ::  a(:)      ! a(m)
    end type

    contains

    subroutine pardiso_solver(m,nth)
        type(sm), intent(in out) :: m
        integer, intent(in) :: nth

        

    end subroutine

end module sle_sm