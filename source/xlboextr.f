      subroutine xlboextr
      use qmmm
      implicit none
      integer i, imat(8), k
      real*8 kappa, coeff(8), alpha, two
      save kappa, coeff, alpha, two
      data kappa/1.86d0/, alpha/0.0016d0/, (coeff(i), i = 1, 8)    
     $  /-36.0d0,99.0d0,-88.0d0,11.0d0,32.0d0,-25.0d0,8.0d0,-1.0d0/,
     $  two/2.0d0/
c
c     assemble a guess density following niklasson's extended lagrangian
c     scheme.
c
      if (nmat.le.7) return
c
c     get the indices of the density matrices so that the last is
c     the eighth:
c
      do i = 1, 8
        imat(i) = mod(nmat-i+1,8) 
        if (imat(i).eq.0) imat(i) = 8
      end do
c     write(6,*) 'nmat, imat(1):',   nmat,   imat(1)
c     write(6,*) 'nmat-1, imat(2):', nmat-1, imat(2)
c     write(6,*) 'nmat-2, imat(3):', nmat-2, imat(3)
c     write(6,*) 'nmat-3, imat(4):', nmat-3, imat(4)
c     write(6,*) 'nmat-4, imat(5):', nmat-4, imat(5)
c     write(6,*) 'nmat-5, imat(6):', nmat-5, imat(6)
c     write(6,*) 'nmat-6, imat(7):', nmat-6, imat(7)
c     write(6,*) 'nmat-7, imat(8):', nmat-7, imat(8)
c
c     guess = two*pguess(8) - pguess(7) + kappa*(pconv - pguess(8)) +
c           + alpha*sum(i=1,8) c(i)*pguess(i)
c
      guess = two*pguess(:,imat(1)) - pguess(:,imat(2))
      guess = guess + kappa*(pconv - pguess(:,imat(1)))
      do i = 1, 8
        guess = guess + alpha*coeff(i)*pguess(:,imat(i))
      end do
c
c     write(6,*) 'guess density:'
c     write(6,'(7f12.6)') (guess(k*(k+1)/2),k=1,7)
c
c     update the pguess array by overwriting the oldest density:
c
      pguess(:,imat(8)) = guess
      return
      end
