import numpy as np
import mpi4py


def init_matrix(matrix: np.ndarray, n: int, first_row: int,
                last_row: int, max_rows: int, my_rank: int,
                p: int, path: str, comm):
    i, j, k = 0, 0, 0
    a = matrix
    buf = np.zeros((max_rows, n))
    tmp = 0
    first_row_k, last_row_k, err = 0, 0, 0
    if my_rank == 0:
        try:
            fp = open(path, "r")
        except:
            err = 1

    # comm.Bcast(err, root=0)
    if err:
        return None
    for k in range(p):
        if my_rank == 0:
            if k == 0:
                for i in range(first_row, last_row + 1):
                    row = np.array(list(map(float, fp.readline().split(' '))))
                    a[i - first_row, :] = row
            else:
                buf.fill(0)
                first_row_k = (n * k) // p
                last_row_k = n * (k + 1) // p - 1
                for i in range(first_row_k, last_row_k + 1):
                    buf[i - first_row_k, :] = np.array(list(map(float, fp.readline().split(' '))))
                comm.Send(buf, dest=k, tag=0)
        else:
            if my_rank == k:
                comm.Recv(buf, source=0, tag=0, status=mpi4py.MPI.Status())
                for i in range(last_row - first_row + 1):
                    a[i, :] = buf[i, :]

    if my_rank == 0:
        fp.close()



def init_vector(vector: np.ndarray, n: int, first_row: int, 
                last_row: int, max_rows: int, my_rank: int, 
                p: int, path: str, comm):
    i, k = 0, 0
    a = vector
    buf = np.zeros(max_rows)
    tmp = 0
    first_row_k, last_row_k, err = 0, 0, 0
    fp = None
    if my_rank == 0:
        try:
            fp = open(path, 'r')
        except:
            err = 1
    # comm.Bcast(err, root=0)
    if err:
        return None
    for k in range(p):
        if my_rank == 0:
            if k == 0:
                for i in range(first_row, last_row + 1):
                    try:
                        num = list(map(float, fp.readline().split()))[0]
                    except:
                        return None
                    a[i - first_row] = num
            else:
                buf.fill(0)
                first_row_k = n * k // p
                last_row_k = n * (k + 1) // p - 1
                for i in range(first_row_k, last_row_k + 1):
                    try:
                        num = list(map(float, fp.readline().split()))[0]
                    except:
                        return None
                    buf[i - first_row_k] = num
                comm.Send(buf, dest=k, tag=0)
        else:
            if my_rank == k:
                comm.Recv(buf, source=0, tag=0, status=mpi4py.MPI.Status())
                a[0:last_row - first_row + 1] = buf[0:last_row - first_row + 1]
    
    if my_rank == 0:
        fp.close()


def print_matrix_full(matr, first_row, last_row, my_rank, size, n, max_rows, comm):
    #rows1 = np.zeros(last_rows - first_rows + 1, dtype=np.int64)
    rows1 = np.array([last_row - first_row + 1]).astype(np.int64)
    buf = np.zeros((max_rows, n), dtype=np.double)

    if my_rank == 0:
        print(matr[:last_row - first_row + 1])

    for i in range(1, size):
        if my_rank == 0:
            comm.Recv(rows1, source=i, tag=0)
            comm.Recv(buf, source=i, tag=0)
            print(buf[:int(rows1[0])])
        elif my_rank == i:
            comm.Send(rows1, dest=0, tag=0)
            comm.Send(matr, dest=0, tag=0)


def check_float(s: str):
    try:
        tmp = float(s)
    except:
        return False
    return True

def check_int(s: str):
    try:
        tmp = int(s)
    except:
        return False
    return True


def matrix_mult_vector(a: np.ndarray, b: np.ndarray, c: np.ndarray, n: int, first_row: int, last_row: int, my_rank: int, p: int):
    i, j, k, l = 0, 0, 0, 0
    rows = last_row - first_row + 1
    max_rows = 0
    first_row_k, last_row_k, rows_k = tuple(np.zeros(3))
    s = 0
    dest, source, tag = tuple(np.zeros(3))
    status = MPI.Status()
    
    max_rows = n // p + n % p
    c.fill(0)

    source = (my_rank + 1) % p
    if my_rank == 0:
        dest = p - 1
    else:
        dest = my_rank - 1
    
    for l in range(p):
        k = (my_rank + l) % p
        first_row_k = n * k // p
        last_row_k = n * (k + 1) // p - 1
        rows_k = last_row_k - first_row_k + 1

        for i in range(rows):
            s = 0
            s += np.sum(a[i, first_row_k:first_row_k + rows_k] * b[:rows_k])
            # print(s, "my_rank is ", my_rank, 'rows is ', rows, "| a is ", a[i, first_row_k:first_row_k + rows_k], " | b is ", b[:rows_k])
            c[i] += s
        comm.Sendrecv_replace(b, dest, tag, source, tag, status)


def matrix_mult_vector_t(a: np.ndarray, b: np.ndarray, c: np.ndarray, n: int, first_row: int, last_row: int, my_rank: int, p: int):
    i, j, k, l = 0, 0, 0, 0
    rows = last_row - first_row + 1
    max_rows = 0
    first_row_k, last_row_k, rows_k = tuple(np.zeros(3))
    s = 0
    dest, source, tag = tuple(np.zeros(3))
    status = MPI.Status()
    
    max_rows = n // p + n % p
    c.fill(0)

    source = (my_rank + 1) % p
    if my_rank == 0:
        dest = p - 1
    else:
        dest = my_rank - 1
    
    for l in range(p):
        k = (my_rank + l) % p
        first_row_k = n * k // p
        last_row_k = n * (k + 1) // p - 1
        rows_k = last_row_k - first_row_k + 1

        for i in range(rows):
            s = a[i, first_row_k:first_row_k + rows_k] * b[i]
            c[:rows_k] += s
            # print(s, "my_rank is ", my_rank, 'rows is ', rows, "| a is ", a[i, first_row_k:first_row_k + rows_k], " | b is ", b[:rows_k])
        comm.Sendrecv_replace(c, dest, tag, source, tag, status)


def bicg(matr, b, x, max_iter, first_row, last_row, max_rows, n, my_rank, size, eps):
    rows = first_row - last_row + 1
    rows = max_rows
    i, j = tuple(np.zeros(2, dtype='int64'))
    err = np.zeros(1, dtype=np.int64)
    alpha, beta = tuple(np.zeros(2, dtype='float64'))
    r = np.zeros(rows)
    r_current = np.zeros(rows)
    p = np.zeros(rows)
    p_current = np.zeros(rows)
    z = np.zeros(rows)
    s = np.zeros(rows)
    result = np.zeros(rows)
    result1 = np.zeros(rows)


    matrix_mult_vector(matr, x, result, n, first_row, last_row, my_rank, size)
    b[first_row:last_row + 1] -= result[first_row:last_row + 1]
    r[:rows] = b[:rows]
    p[:rows] = r[:rows]
    z[:rows] = r[:rows]
    s[:rows] = r[:rows]

    sum = np.zeros(1)
    res = np.zeros(1)
    res1 = np.zeros(1)
    norm = np.zeros(1)

    for j in range(max_iter):
        sum = np.array([np.sum(p[:last_row - first_row + 1] * r[:last_row - first_row + 1])])
        sum1 = np.array([np.sum(r[:last_row - first_row + 1] ** 2)])
        comm.Allreduce(sum, res, MPI.SUM)
        comm.Allreduce(sum1, norm, MPI.SUM)

        if my_rank == 0:
            norm = np.sqrt(norm)
            if norm[0] < eps:
                print('BICG Converged. Euclidean norm is ', norm)
                err = np.array([1])
        
        comm.Bcast(err, 0)
        # print('my_rank is ', my_rank, 'my err is ', err)
        if err[0] == 1:
            # print('NIger number', my_rank, 'is done')
            return None
    
        result.fill(0)
        matrix_mult_vector(matr, z, result, n, first_row, last_row, my_rank, size)
        sum = np.array([np.sum(s[:last_row - first_row + 1] * result[:last_row - first_row + 1])])
        comm.Allreduce(sum, res1, MPI.SUM)
        matrix_mult_vector_t(matr, s, result1, n, first_row, last_row, my_rank, size)
        alpha = res[0] / res1[0]
        x[:rows] += alpha * z[:rows]
        r_current[:rows] = r[:rows] - alpha * result[:rows]
        p_current[:rows] = p[:rows] - alpha * result1[:rows]
        comm.Barrier()
        res = np.zeros(1)
        res1 = np.zeros(1)
    
        sum = np.array([np.sum(p_current[:last_row - first_row + 1] * r_current[:last_row - first_row + 1])])
        comm.Allreduce(sum, res1, MPI.SUM)
        beta = res[0] / res1[0]
        z[:rows] = r_current[:rows] + beta * z[:rows]
        s[:rows] = p_current[:rows] + beta * s[:rows]
        r[:rows] = r_current[:rows]
        p[:rows] = p_current[:rows]
    comm.Barrier()

from mpi4py import MPI
import sys

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

# Проверка аргументов командной строки
if not (len(sys.argv) == 6 and
        check_int(sys.argv[1]) and
        check_int(sys.argv[2]) and
        check_float(str(sys.argv[3]))):

    if my_rank == 0:
        print("Usage:", sys.argv[0], "n max_iter eps [file] [file]")
    comm.Barrier()  # Синхронизация всех процессов
    MPI.Finalize()  # Завершение MPI
    sys.exit(0)

n = int(sys.argv[1])
max_iter = int(sys.argv[2])
eps = float(sys.argv[3])

path_to_matrix = str(sys.argv[4])
path_to_vec = str(sys.argv[5])

# Теперь у нас есть переменные n, max_iter и eps с правильными значениями, а также проверка на количество аргументов. # Продолжайте свой код отсюда и используйте переменные n, max_iter и eps в вашей программе. 
# Не забудьте синхронизироваться перед завершением программы. comm.Barrier()

first_row = n * my_rank // p
last_row = n * (my_rank + 1) // p - 1

max_rows = n // p + n % p

try:
    matrix = np.zeros((max_rows, n))
    vector = np.zeros(max_rows)
    result = np.zeros(max_rows)
except:
    print('Ошибка выделения памяти')

init_matrix(matrix, n, first_row, last_row, max_rows, my_rank, p, path_to_matrix, comm=comm)
#print_matrix_full(matrix, first_row, last_row, my_rank, p, n, max_rows, comm)
init_vector(vector, n, first_row, last_row, max_rows, my_rank, p, path_to_vec, comm=comm)
#print_matrix_full(vector, first_row, last_row, my_rank, p, 1, max_rows, comm)

# matrix_mult_vector_t(matrix, vector, result, n, first_row, last_row, my_rank, p)
bicg(matrix, vector, result, max_iter, first_row, last_row, max_rows, n, my_rank, p, eps)
comm.Barrier()
print_matrix_full(result, first_row, last_row, my_rank, p, 1, max_rows, comm)

MPI.Finalize()
