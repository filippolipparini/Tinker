      subroutine qmwrite(iunit,istep,dt,wtstep)
      use qmmm
      implicit none
c
      integer iunit, istep, lenfmt
      real*8  dt, wtstep, psday
      character*(40) fmtu
c
      open (unit=iunit,file='qminfo.log',status='unknown',
     $  access='append')
 1000 format(i6,2x,f18.4,i5,f18.4,3f10.4,2d10.2,f9.2,f12.4)
      psday = 86400.0d0*dt/(wtstep)
      write(iunit,1000) istep,escf,itscf,etd,qmdip(1),qmdip(2),qmdip(3),
     $  maxfor(1),maxfor(2),wtstep,psday
c
      if (mod(istep,100).eq.0) write(iunit,*)
      close (100)
c
c     if there are properties available, dump them on file, too:
c
      write(fmtu,'(a,i4,a)') '(',lenprop,'f14.8)'
      lenfmt = len(trim(fmtu))
      if (nprops.ne.0) then
        open (unit=100,file='qmprops.log',status='unknown',
     $    access='append')
        write(100,fmtu(1:lenfmt)) proparray
        close (100)
      end if
      return
      end
