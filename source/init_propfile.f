      subroutine init_propfile
      use qmmm
      implicit none
      integer i
      character*10 dname(100), pname(100)
      save dname, pname
      data (pname(i),i=1,3)
     $  /'Dipole    ','Dipole Vel','r x Del   '/
      data (dname(i),i=1,3)
     $  /'       SCF','     Total','Transition'/
      if (nprops.gt.0) then
        open (unit=100,file='qmprops.log',status='unknown',
     $    access='sequential')
        write(100,*) ' dumping the following QM properties '
        write(100,*) ' all the properties are consecutive floats in',
     $               ' one line per md step.'
        do i = 1, nprops
          write(100,*) dname(idens(i)),' ',pname(iprop(i))
        end do
        close (100)
      end if
      return
      end
