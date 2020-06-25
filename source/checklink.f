c
c     "checkling" checks for bonds between the QM and MM parts
c
      subroutine checklink
      use usage
      use couple
      use iounit
      implicit none
      integer i, j, k, l
      integer maxla, nlinkatoms
      parameter (maxla=100)
      integer, allocatable :: linkatoms(:,:)
      
c
      write(6,*) 'Initializing link atoms'
      allocate(linkatoms(2,maxla))
c     allocate(qmneedff(4*maxla))
      nlinkatoms = 0
      linkatoms = 0
c
c     loop over qm atoms and check for mm neighbors.
c
      do i = 1, nqmatoms
        j = qmlist(i)
        do k = 1, n12(j)
          if (.not.qmatoms(i12(k,j))) then
            nlinkatoms = nlinkatoms + 1
            linkatoms(1,nlinkatoms) = j
            linkatoms(2,nlinkatoms) = i12(k,j)
          end if
        end do
      end do
c
      if (nlinkatoms.eq.0) then
        write(iout,100) 
        return
      end if
  100 format(/,' there is no covalent bond between the qm and mm part')
      write(iout,101) nlinkatoms
  101 format(/,' there are ',i3,' covalent bonds between the qm and',
     $  ' mm subsystems')
c
      do i = 1, nlinkatoms
        write(6,*) linkatoms(1,i), linkatoms(2,i)
      end do
c
c    check for multiple bonds between qm and mm interesting
c    the same atoms.
c
      do i = 1, nlinkatoms
        do j = i + 1, nlinkatoms
          if (linkatoms(1,i).eq.linkatoms(1,j)) then
            write (iout,110)
  110       format (/,' INITLA  --  The same QM atom cannot be linked',
     &                ' to more than one MM atom.')
            call fatal
          else if (linkatoms(2,i).eq.linkatoms(2,j)) then
  120       format (/,' INITLA  --  The same MM atom cannot be linked',
     &                ' to more than one QM atom.')
            write (iout,120)
            call fatal
          end if
        end do
      end do
c
c     set the list of qm atoms which need the force field
c     1-2 from the qm-mm bonds.
c
c     nqmff = 0
c     do i = 1, nlinkatoms
c       j = linkatoms(1,i)
c       nqmff = nqmff + 1
c       qmneedff(nqmff) = j
c       do k = 1, n12(j)
c         if (qmatoms(i12(k,j))) then
c           nqmff = nqmff + 1
c           qmneedff(nqmff) = i12(k,j)
c         end if
c       end do
c     end do
c     write(6,*) 'qm need ff'
c     do i = 1, nqmff
c       write(6,*) qmneedff(i)
c     end do
      return
      end
