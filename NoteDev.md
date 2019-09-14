# 

# SCDMk-phonon procedure:
  * dephase all the eigenvectors for each k-point
  * select anchor point
  * sum up the projections to the anchor points for each eigenvector
  * multiply by weight functions
  * select column
  * SVD -> Amn(k)
  * Amn(k) & psi -> wannier functions for each k
  * H0(k) and Amn(k) -> H(k)
  * 
