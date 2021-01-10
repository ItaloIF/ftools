! Fast Fourier Transform - module
! $II-012020-1

module fftm
    implicit none
    contains

    subroutine pw2(n,m)
        ! power of two greater than n
        integer(8) :: n
        integer(4) :: m
        m = 0
        do while (n.gt.2**m)
            m = m + 1
        end do
    end subroutine pw2

    subroutine fft(m,x)
        ! fast fourier transform - 1D
        integer(8), intent(in) :: m
        complex(8), intent(in out) :: x(2**m)
        call ctalg(m,x)
        call calfft(m,x)
    end subroutine fft

    subroutine ctalg(m,x)
        ! Cooley-Tukey FFT algorithm - 1D
        integer(8), intent(in) :: m
        complex(8), intent(in out) :: x(2**m)
        integer(8) :: i
        integer(8) :: j
        integer(8) :: n
        integer(8) :: r
        complex(8) :: tp
        n = 2**m
        j = 0
        do i = 1,n-1
            if ((i-1).le.j) then
                tp = x(j+1)
                x(j+1) = x(i)
                x(i) = tp
            end if
            r = n/2
            do while (j.ge.r)
                j = j - r
                r = r/2
            end do
            j = j + r
        end do
    end subroutine ctalg



end module fftm
