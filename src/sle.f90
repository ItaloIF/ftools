! Solvers System of Lineal Equations
! $II-022020-2
!		sle_u	: solver upper triangular matrix
!		sle_l	: solver lower triangular matrix
!		sle_lu	: solver  matrix LU method
!		des_lu	: matrix descomposition A = L·U
!		sle_gss	: solver  matrix gauss method
!		sle_ch	: solver  matrix cholesky method
!		des_ch	: matrix descomposition A = L·Lt
!		sle_gc	: solver conjugate gradient method
!		sle_gs	: solver cg method with gauss-seidel precondition
!		gs_lu	: gauss-seidel precondition

module sle
    implicit none
    contains

    subroutine sle_u(u,x)
        ! solver upper triangular matrix
        implicit none
        real(8), intent(in) :: u(:,:)
        real(8), intent(in out) :: x(:)
        
        integer(8) :: n
        integer(8) :: i
        n = ubound(u,1)
        x(n) = x(n)/u(n,n)
        do i = 1,n-1
            x(n-i) = (x(n-i) - dot_product(u(n-i,n-i+1:n),x(n-i+1:n)))/u(n-i,n-i)
        end do
    end subroutine sle_u

    subroutine sle_l(l,x)
        ! solver lower triangular matrix
        implicit none
        real(8), intent(in) :: l(:,:)
        real(8), intent(in out) :: x(:)

        integer(8) :: n
        integer(8) :: i
        n = ubound(l,1)
        x(1) = x(1)/l(1,1)
        do i = 2,n
            x(i) = (x(i) - dot_product(l(i,1:i-1),x(1:i-1)))/l(i,i)
        end do
    end subroutine sle_l

    subroutine sle_lu(a,x)
        ! solver  matrix A = L·U method
        ! this subroutine doesn't change matrix A
        implicit none
        real(8), intent(in) :: a(:,:)
        real(8), intent(in out) :: x(:)

        integer(8) :: n
        integer(8) :: i
        real(8), allocatable :: t(:)
        real(8), allocatable :: at(:,:)
        n = ubound(a,1)
        allocate(t(n))
        allocate(at(n,n))
        at = a
        call des_lu(at)
        forall(i = 1:n) t(i) = at(i,i)
        forall(i = 1:n) at(i,i) = 1.0_8
        call sle_l(at,x)
        forall(i = 1:n) at(i,i) = t(i)
        call sle_u(at,x)
    end subroutine sle_lu

    subroutine des_lu(a)
        ! matrix descomposition A = L·U
        ! this subroutine changes matrix A (L\U)
        implicit none
        real(8), intent(in out) :: a(:,:)

        integer(8) :: n
        integer(8) :: i
        n = ubound(a,1)
        a(2:n,1) = a(2:n,1)/a(1,1)
        do i = 2,n
            a(i,i) = a(i,i) - dot_product(a(i,1:i-1),a(1:i-1,i))
            a(i,i+1:n) = a(i,i+1:n) - matmul(a(i,1:i-1),a(1:i-1,i+1:n))
            a(i+1:n,i) = (a(i+1:n,i) - matmul(a(i+1:n,1:i-1),a(1:i-1,i)))/a(i,i)
        end do
    end subroutine des_lu

    subroutine sle_gss(a,x)
        ! solver  matrix gauss method
        implicit none
        real(8), intent(in) :: a(:,:)
        real(8), intent(in out) :: x(:)

        integer(8) :: n
        integer(8) :: i
        integer(8) :: j
        real(8), allocatable :: at(:,:)
        n = ubound(a,1)
        allocate(at(n,n))
        at = a
        do i = 1,n-1
            do j = 1+1,n
                x(j) = x(j) - at(j,i)/at(i,i)*x(i)
                at(j,i:n) = at(j,i:n)-at(j,i)/at(i,i)*at(i,i:n)
            end do
        end do
        call sle_u(at,x)
    end subroutine sle_gss

    subroutine sle_ch(a,x)
        ! solver  matrix cholesky method
        implicit none
        real(8), intent(in) :: a(:,:)
        real(8), intent(in out) :: x(:)

        integer(8) :: n
        real(8), allocatable :: at(:,:)
        n = ubound(a,1)
        allocate(at(n,n))
        at = a
        call des_ch(at)
        call sle_l(at,x)
        call sle_u(at,x)
    end subroutine sle_ch

    subroutine des_ch(a)
        ! matrix descomposition A = L·Lt
        ! only for a Symetric and Definite Positive Matrix
        implicit none
        real(8), intent(in out) :: a(:,:)
        
        integer(8) :: n
        integer(8) :: i
        n = ubound(a,1)
        a(1,1) = a(1,1)**0.5
        a(1,2:n) = a(1,2:n)/a(1,1)
        a(2:n,1) = a(1,2:n)
        do i = 2,n
            a(i,i) = (a(i,i) - dot_product(a(1:i-1,i),a(1:i-1,i)))**0.5
            a(i,i+1:n) = (a(i,i+1:n) - matmul(a(1:i-1,i),a(1:i-1,i+1:n)))/a(i,i)
            a(i,1:i-1) = a(1:i-1,i)
        end do
    end subroutine des_ch

    subroutine sle_gc(a,x,tol,imax)
        ! solver conjugate gradient method
        ! relajation method
        ! important tolerance and max number of iterations
        implicit none
        real(8), intent(in) :: a(:,:)
        real(8), intent(in out) :: x(:)
        real(8), intent(in) :: tol
        integer(8), intent(in) ::  imax

        integer(8) :: n
        integer(8) :: i
        real(8) :: er
        real(8) :: alpha
        real(8) :: betha
        real(8), allocatable :: r(:)
        real(8), allocatable :: s(:)
        real(8), allocatable :: q(:)

        n = ubound(a,1)
        allocate(r(n))
        allocate(s(n))
        allocate(q(n))
        r = x
        s = x
        x = 0
        i = 1
        ! iterations
        do
            q = matmul(a,s)
            alpha = dot_product(r,r)/dot_product(s,q)
            x = x + alpha*s
            er = sqrt(dot_product(alpha*s,alpha*s))
            if (er.lt.tol.or.i.gt.imax) exit
            i = i + 1
            r = r - alpha*q
            betha = dot_product(r,q)/dot_product(s,q)
            s = r - betha*s
        end do 
    end subroutine sle_gc

    subroutine sle_gs(a,x,tol,imax)
        ! solver conjugate gradient method with Gauss-Seidel precondition
        ! relajation method
        ! important tolerance and max number of iterations
        implicit none
        real(8), intent(in) :: a(:,:)
        real(8), intent(in out) :: x(:)
        real(8), intent(in) :: tol
        integer(8), intent(in) ::  imax

        integer(8) :: n
        integer(8) :: i
        real(8) :: er
        real(8) :: alpha
        real(8) :: betha
        real(8), allocatable :: r(:)
        real(8), allocatable :: s(:)
        real(8), allocatable :: q(:)
        real(8), allocatable :: p(:)

        n = ubound(a,1)
        allocate(r(n))
        allocate(s(n))
        allocate(q(n))
        allocate(p(n))
        r = x
        s = x
        call gs_lu(a,p)
        s = p
        x = 0
        do
            q = matmul(a,s)
            alpha = dot_product(p,r)/dot_product(s,q)
            x = x + alpha*s
            if (er.lt.tol.or.i.gt.imax) exit
            i = i + 1
            betha = 1/dot_product(p,r)   ! check order
            r = r - alpha*q
            p = r
            call  gs_lu(a,p)
            betha = betha*dot_product(p,r)
            s = p + betha*s
            q =  matmul(a,s)
            alpha = dot_product(p,r)/dot_product(s,q)
            x = x + alpha*s
        end do
    end subroutine sle_gs

    subroutine gs_lu(a,p)
        ! Gauss-Seidel precondition
        implicit none
        real(8), intent(in) :: a(:,:)
        real(8), intent(in out) :: p(:)
        
        integer(8) :: n
        integer(8) :: i
        real(8), allocatable :: m(:,:)

        n = ubound(a,1)
        allocate(m(n,n))
        do i = 1,n-1
            m(i+1:n,i) = a(i+1:n,i)/a(i,i)
        end do
        forall(i = 1:n) m(i,i) = m(i,i) + 1
        call sle_u(m,p)
        call sle_l(a,p)
    end subroutine gs_lu
end module sle