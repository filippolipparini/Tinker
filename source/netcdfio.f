c
c     "netcdfio" writes out a set of Cartesian coordinates
c     to an external disk file in the amber netcdf format
c
      subroutine netcdfio(init,istep,dt)
#ifdef USE_NETCDF
      use netcdf
      use sizes
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use titles
      implicit none
      integer i,j,k,ixyz
      integer*4 size,crdsiz
      integer*4 ncid,status,nf90_create_
      integer*4 dimframe,dimspatial,dimatoms,dimcells,dimcella
      integer*4 dimidchar(1),dimidcellspatial(1),dimidcellangular(1)
      integer*4 dimidframe(1),dimidcoordinates(3)
      integer*4 dimidcelll(2),dimidcella(2)
      integer*4 varid,nloc
      integer istep
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 xyzfile
      real, allocatable :: temp(:,:)
      real*8 dt
      real tempcell(3)
      real pico
      integer*4 starttime(1),starttab(3),starttab1(2),counttab(3)
      integer*4 counttab1(2)
      logical init,exist
      character*240 netcdffile
c     integer nf90_def_dim, nf90_create, nf90_clobber, nf90_unlimited
c     integer nf90_def_var, nf90_char, nf90_float, nf90_put_att
c     integer nf90_global, nf90_enddef, nf90_open, nf90_write
c     integer nf90_inq_varid, nf90_put_var, nf90_close
c     external nf90_def_dim, nf90_create, nf90_clobber, nf90_unlimited
c     external nf90_def_var, nf90_char, nf90_float, nf90_put_att
c     external nf90_global, nf90_enddef, nf90_open, nf90_write
c     external nf90_inq_varid, nf90_put_var, nf90_close
c

      nloc = int(n,4)
c
      netcdffile = filename(1:leng)//'.nc'
      init = istep .eq. iwrite
      if (init) then
c  
c     create the netcdf dataset
c
        status = nf90_create(netcdffile,nf90_clobber,ncid)
c
c     create the netcdf dimensions
c
        status = nf90_def_dim(ncid,"frame",NF90_unlimited,dimframe)
        status = nf90_def_dim(ncid,"spatial",3,dimspatial)
        status = nf90_def_dim(ncid,"atom",nloc,dimatoms)
        status = nf90_def_dim(ncid,"cell_spatial",3,dimcells)
        status = nf90_def_dim(ncid,"cell_angular",3,dimcella)
        dimidchar(1) = int(dimspatial,4)
        dimidcellspatial(1) = int(dimcells,4)
        dimidcellangular(1) = int(dimcella,4)
        dimidframe(1) = int(dimframe,4)
        dimidcoordinates(1) = int(dimspatial,4)
        dimidcoordinates(2) = int(dimatoms,4)
        dimidcoordinates(3) = int(dimframe,4)
        dimidcelll(1) = int(dimcells,4)
        dimidcelll(2) = int(dimframe,4)
        dimidcella(1) = int(dimcella,4)
        dimidcella(2) = int(dimframe,4)
c
c     create the netcdf variables
c
        status = nf90_def_var(ncid,"spatial",NF90_CHAR,dimidchar,varid)

        status = nf90_def_var(ncid,"cell_spatial",NF90_CHAR,
     $     dimidcellspatial,varid)

        status = nf90_def_var(ncid,"cell_angular",NF90_CHAR,
     $     dimidcellangular,varid)

        status = nf90_def_var(ncid,"time",NF90_FLOAT,dimidframe,varid)
        status = nf90_put_att(ncid,varid,"units","picosecond")

        status = nf90_def_var(ncid,"coordinates",NF90_FLOAT,
     $     dimidcoordinates,varid)
        status = nf90_put_att(ncid,varid,"units","angstrom")

        status = nf90_def_var(ncid,"cell_lengths",NF90_FLOAT,
     $     dimidcelll,varid)
        status = nf90_put_att(ncid,varid,"units","angstrom")

        status = nf90_def_var(ncid,"cell_angles",NF90_FLOAT,
     $     dimidcella,varid)
        status = nf90_put_att(ncid,varid,"units","degree")

        status = nf90_put_att(ncid,NF90_GLOBAL,"title","amber output")
        status = nf90_put_att(ncid,NF90_GLOBAL,"application","tinker")
        status = nf90_put_att(ncid,NF90_GLOBAL,"program","tinker")
        status = nf90_put_att(ncid,NF90_GLOBAL,"programVersion","8.2.1")
        status = nf90_put_att(ncid,NF90_GLOBAL,"Conventions","AMBER")
        status = nf90_put_att(ncid,NF90_GLOBAL,"ConventionVersion",
     $   "1.0")
         
        status = nf90_enddef(ncid)
      else
c
c     open the netcdf dataset
c
        status = nf90_open(netcdffile,nf90_write,ncid)
      endif
c
c     append the netcdf file
c
      status = nf90_inq_varid(ncid,"time",varid)
      pico = real(istep) * real(dt)
      starttime(1) = int(istep,4)
      status = nf90_put_var(ncid,varid,pico,starttime)

      status = nf90_inq_varid(ncid,"cell_lengths",varid)
      tempcell(1) = real(xbox)
      tempcell(2) = real(ybox)
      tempcell(3) = real(zbox)
      starttab1(1) = int(1,4)
      starttab1(2) = int(istep ,4)
      counttab1(1) = int(3,4)
      counttab1(2) = int(1,4)
      status = nf90_put_var(ncid,varid,tempcell,starttab1,counttab1)

      status = nf90_inq_varid(ncid,"coordinates",varid)
      allocate (temp(3,n))
      do i = 1, n
        temp(1,i) = real(x(i))
        temp(2,i) = real(y(i))
        temp(3,i) = real(z(i))
      end do
      starttab(1) = int(1 ,4)
      starttab(2) = int(1 ,4)
      starttab(3) = int(istep ,4)
      counttab(1) = int(3,4)
      counttab(2) = int(n,4)
      counttab(3) = int(1,4)
      status = nf90_put_var(ncid,varid,temp,starttab,counttab)
      deallocate (temp)

      status = nf90_inq_varid(ncid,"cell_angles",varid)
      tempcell(1) = real(alpha)
      tempcell(2) = real(beta)
      tempcell(3) = real(gamma) 
      starttab1(1) = int(1,4)
      starttab1(2) = int(istep ,4)
      counttab1(1) = int(3,4)
      counttab1(2) = int(1,4)
      status = nf90_put_var(ncid,varid,tempcell,starttab1,counttab1)
c
c     close the netcdf dataset
c
      status = nf90_close(ncid)
      return
#endif
      end
