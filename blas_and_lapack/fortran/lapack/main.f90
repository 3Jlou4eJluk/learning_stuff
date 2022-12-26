program calling_lapack_svd
    implicit none

    integer, parameter :: m = 3, n = 3
    integer :: i, j
    integer :: info, lda, ldu, ldv, lwork, sigma_size

    double precision, allocatable       :: a_copy(:, :), b(:), s(:), u(:, :), vs(:, :), work(:)
    double precision, dimension(m, n)   :: mat
    integer, allocatable       :: iwork(:)
    double precision                    :: stupid_thing(1, 1)

    if (m > n) then
        sigma_size = n
    else
        sigma_size = m
    end if

    allocate (s(sigma_size), vs(n, n), u(m, m), iwork(8 * m))

    read(*,*) ((mat(i, j), j = 1, m), i = 1, m)

    lwork = -1
    lda = m
    ldu = m
    ldv = n




    call dgesdd('A', m, n, mat, lda, s, u, ldu, vs, ldv, stupid_thing, lwork, iwork, info)

    lwork = stupid_thing(1, 1)
    allocate (work(lwork))

    call dgesdd('A', m, n, mat, lda, s, u, ldu, vs, ldv, work, lwork, iwork, info)


    write(*,*) "U is "

    call print_matrix(u, m, m)

    write(*,*) "sigma is "

    call printArray(s, sigma_size)

    write(*,*) "V is "

    call print_matrix(vs, n, n)
    
    contains

    subroutine print_matrix(matrix, matrix_vsize, matrix_hsize)
        implicit none
        ! Data section Start
        integer                                                     :: i, j
        integer                                                     :: matrix_vsize, matrix_hsize
        double precision, dimension(:, :)     :: matrix
        ! Data section End

        do i = 1, matrix_vsize
            do j = 1, matrix_hsize
                if (j < matrix_hsize) then
                    write(*, 100) matrix(i, j)
                    100 format (D10.3, ' ', $)
                else
                    write(*, 101) matrix(i, j)
                    101 format (D10.3, ' ')
                end if

            end do
        end do
    end subroutine

    subroutine printArray(a, len)
        integer :: len
        double precision, dimension (:) :: a  
        integer::i
   
        do i = 1, len
            Print *, a(i)
        end do
   
    end subroutine printArray
end program calling_lapack_svd