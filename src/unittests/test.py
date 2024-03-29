import numpy as np
from scipy.linalg import qr
def test_svd():
    a = (8.79, 6.11, -9.15, 9.57, -3.49, 9.84, 9.93, 6.91, -7.93, 1.64, 4.02, 0.15,
         9.83, 5.04, 4.86, 8.83, 9.80, -8.99, 5.45, -0.27, 4.85, 0.74, 10.00,
         -6.02, 3.16, 7.98, 3.01, 5.80, 4.27, -5.31)
    a = np.reshape(a, (6,5))
    U, S, VT=np.linalg.svd(a)
    print("U:", U)
    print("S:", S)
    print("VT:", VT)

    #Q, R, piv=qr(U.T.conj(), mode='full', pivoting=True)
    Q, R, piv=qr(a, mode='full', pivoting=True)
    print(piv)

test_svd()
