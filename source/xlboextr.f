      subroutine xlboextr
      use qmmm
      implicit none
      integer i, imat(8), k, iter, j, ij
      real*8 kappa, coeff(8), alpha, two
      real*8 :: err_rms, err_max, thres, thres1
      real*8, allocatable :: p2(:,:), p3(:,:), pscr(:,:), perr(:,:)
      logical conv
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
      allocate (p2(nbasis,nbasis),p3(nbasis,nbasis),pscr(nbasis,nbasis))
      allocate (perr(nbasis,nbasis))
c
      thres  = 1.0d-8
      thres1 = 1.0d-7
c     write(6,*) 'nbasis=', nbasis
      ij = 0
      do i = 1, nbasis
        do j = 1, i
          ij = ij + 1
          pscr(i,j) = guess(ij)
          pscr(j,i) = guess(ij)
        end do
      end do
      pscr = 0.50d0*pscr
c
      conv = .false.
      do iter = 1, 11
        p2 = matmul(pscr,pscr)
c       write(6,*) 'p^2'
c       write(6,'(12f8.4)') p2
        perr = abs(p2 - pscr)
        err_rms = 0.0d0
        err_max = 0.0d0
        do i = 1, nbasis
          do j = 1, nbasis
            err_max = max(err_max,perr(i,j))
            err_rms = err_rms + perr(i,j)**2
          end do
        end do
        err_rms = sqrt(err_rms/float(nbasis))
c       write(6,1000) iter-1, err_max, err_rms
        if (err_rms.lt.thres .and. err_max.lt.thres1) then
          conv = .true.
        end if
        p3 = matmul(p2,pscr)
        pscr = 3.0d0*p2 - 2.0d0*p3
        if(conv) goto 999
      end do
 999  continue
 1000 format(t3,'it = ', i3, ' error (rms,max): ', 2d9.2)
      if (conv) then
        pscr = 2.0d0*pscr
        ij = 0
        do i = 1, nbasis
          do j = 1, i
            ij = ij + 1
            guess(ij) = pscr(i,j)
          end do
        end do
      end if
      deallocate (p2,p3,pscr,perr)
c
c     write(6,*) 'guess density:'
c     write(6,'(7f12.6)') (guess(k*(k+1)/2),k=1,7)
c
c     update the pguess array by overwriting the oldest density:
c
      pguess(:,imat(8)) = guess
      return
      end
