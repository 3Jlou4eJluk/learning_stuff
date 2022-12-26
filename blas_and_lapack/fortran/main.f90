program calling_blas_multiplication
    implicit none

    integer, parameter :: m = 3
    integer, parameter :: n = 3
    integer, parameter :: k = 3

    integer :: i, j
    double precision, dimension(m, n) :: mat1
    double precision, dimension(n, k) :: mat2
    double precision, dimension(m, k) :: mat3

    read(*,*) ((mat1(i, j), j = 1, m), i = 1, m)
    read(*,*) ((mat2(i, j), j = 1, n), i = 1, k)

    call dgemm('N', 'N', m, n, k,  1D0, mat1, n, mat2, k, 0D0, mat3, k)

    call print_matrix(mat3, m, k)

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

end program calling_blas_multiplication 
