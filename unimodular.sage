import numpy as np

R.<x> = QQ[]
der = R.derivation()

A.<D> = R['D',der]

def null_space_vector(M):
    k = kernel(M)
    return (list(k.basis_matrix()[0]))

def find_t_index(v, delta, n):
    appr_delta = [-1 if v[i] == 0 else delta[i] for i in range(n)]

    return appr_delta.index(max(appr_delta))

def frMat_and_rowRanks(M):
    rrs = []
    fm = []
    for row in M:
        rr = max(map(lambda e: 0 if e == 0 else e.degree(), row))
        rrs.append(rr)
        fr = list(map(lambda e: e.leading_coefficient() if e.degree() == rr else 0, row))
        fm.append(fr)
    return (matrix(fm), rrs)

def rowRed(M):
    n = M.nrows()
    L = copy(M)
    U = identity_matrix(R, n)
    U = U.change_ring(A)
    U1 = copy(U)
    Us = []
    Us.append(U1)
#     print(Us)
    frM, rr = frMat_and_rowRanks(L)

    while frM.rank() < n:
#     for i in range(4):
        v = null_space_vector(frM)
        t = find_t_index(v, rr, n)
        L_t_th_row = sum ([ (0 if v[i] == 0 else v[i]*D^(rr[t]-rr[i])) * L[i] for i in range(n)])
        U_t_th_row = sum ([ (0 if v[i] == 0 else v[i]*D^(rr[t]-rr[i])) * U[i] for i in range(n)])
        L.set_row(t,L_t_th_row)
        U.set_row(t,U_t_th_row)
        U1 = copy(U)
        frM, rr = frMat_and_rowRanks(L)
        Us.append(U1)
#         print(Us)
#         print(U)
#         print(L)
#         print(frM)
#         print('v: ', v)
#         print('t: ', t)
#         print('delta: ',rr)
#         print('-------')
    return (L, Us, max(rr))

def isUnimodular(M):
    Mrr, Us, ordM = rowRed(M)
    return True if ordM == 0 else False

