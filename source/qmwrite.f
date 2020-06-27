      subroutine qmwrite(iunit,istep,dt,wtstep)
      use qmmm
      implicit none
c
      integer iunit, istep
      real*8  dt, wtstep, psday
c
      open (unit=iunit,file='qminfo.log',status='unknown',
     $  access='append')
 1000 format(i6,2x,2f18.4,3f10.4,2d10.2,f9.2,f12.4)
      psday = 86400.0d0*dt/(wtstep)
      write(iunit,1000) istep,escf,etd,qmdip(1),qmdip(2),qmdip(3),
     $  maxfor(1),maxfor(2),wtstep,psday
c
      if (mod(istep,100).eq.0) write(iunit,*)
      close (100)
c
c     if there are properties available, dump them on file, too:
c
      if (nprops.ne.0) then
        open (unit=100,file='qmprops.log',status='unknown',
     $    access='append')
        write(100,*) proparray
        close (100)
      end if
      return
      end
