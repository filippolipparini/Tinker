c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kqmmm  --  qmmm term parameter assignment    ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kextra" assigns parameters to any additional user defined
c     potential energy contribution
c
c
      subroutine kqmmm
      use sizes
      use atoms
      use atomid
      use keys
      use qmmm
      use usage
      implicit none
      integer i,k,next,istat,nattmp,status
      real*8 rd
      logical header,havemat
      character*240 keyword
      character*240 value
      character*240 record
      character*240 string
      character*240 command
c
c     read the name of the gaussian input and matrix elements file.
c     note that the program expects a gau_name.com and mat_name.mat
c     file to be present. 
c     if no file is explicitly given in input, the program assumes
c     that the files are called gau.com and gau.mel
c
      gau_name  = 'gau.com'
      mat_name  = 'gau.mat'
      lgname    = 7
      lmname    = 7
      called    = .false.
      nmat      = 0
      nprops    = 0
      qmmm_post = .false.
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:13) .eq. 'QMMM-GAUFILE ') then
            call gettext(record,value,next)
            lgname = len(trim(value))
            gau_name(1:lgname) = trim(value)
         else if (keyword(1:13) .eq. 'QMMM-MATFILE ') then
            call gettext (record,value,next)
            lmname = len(trim(value))
            mat_name(1:lmname) = trim(value)
         else if (keyword(1:10) .eq. 'QMMM-XLBO ') then
            xlbo = .true.
         else if (keyword(1:8) .eq. 'QMPROPS ') then
            read(string,*,err=10,end=10) (idens(k), iprop(k), 
     $        k = nprops+1,100)
  10        continue
            do while (idens(nprops+1).ne.0 .and. iprop(nprops+1).ne.0)
              nprops = nprops + 1
            end do
         else if (keyword(1:12) .eq. 'POSTPROCESS ') then
            call getword (record,value,next)
            lscrname = len(trim(value))
            scr_name(1:lscrname) = trim(value)
            qmmm_post = .true.
         end if
      end do
c
c     If doing xlbo, create a second gaussian file with the appropriate
c     keywords.
c
      if (xlbo) then
        write(command,*) 'cat ', gau_name(1:lgname),
     $    ' | sed -e s/useoao/readoao/gI > ',gau_name(1:lgname-4),
     $    '_xlbo.com'
        status = system(command)
        if (status.ne.0) then
          write(6,*) 'Could not prepare XLBO .com file.'
          call fatal
        end if
      end if
c
c     allocate memory for qmmm:
c
      allocate (ian(n), iattyp(n), ibfatm(maxbas), ibftyp(maxbas),
     $  atmchg(n), cgau(3,n), atmwgt(n), stat=istat)
      if (istat.ne.0) then
         write(6,*) ' allocation error in kqmmm.'
         call fatal
      end if
c
c     If doing full qm use gaussian to create a matrix element file.
c
      if (nqmatoms.eq.n) then
        inquire(file=mat_name(1:lmname),exist=havemat)
        if (.not.havemat) then
          write(command,*) 'cat ',gau_name(1:lgname-4),
     $      '_init.com > bldmat.com'
          status = system(command)
          if (status.ne.0) then
            write(6,100)
  100       format(' A full QM calculation requires either a matrix',
     $        ' element file',/,
     $        ' or a .com (basename_init.com) containing',
     $        ' root, title, charge and spin multiplicity.')
            call fatal
          end if 
          open (unit=100,file='bldmat.com',status='Unknown',
     $      form='formatted',access='Append')
          do i = 1, n
            write(100,'(1i3,3x,3f16.8)') atomic(i),x(i),y(i),z(i)
          end do
          write(100,*)
          write(100,*) mat_name(1:lmname)
          write(100,*)
          close(100)
          write(command,*) 'gdvtest bldmat.com'
          status = system(command)
        end if
      end if
c
c     create an unformatted chk file for gaussian.
c     this will store all the information that stays unaltered through the
c     dynamics.
c
      write(command,*) 'unfchk -mat ',mat_name(1:lmname)
      status = system(command)
c
c     if required, allocate an array to store qm properties and compute its
c     length
c
      if (nprops.ne.0) then
c
        do_dip  = .false.
        do_del  = .false.
        do_rdel = .false.
c  
        do i = 1, nprops
          if (iprop(i).eq.1) do_dip  = .true.
          if (iprop(i).eq.2) do_del  = .true.
          if (iprop(i).eq.3) do_rdel = .true.
        end do
c
c       at the moment, only dipoles can be handled
c
        if (do_del .or. do_rdel) then
          write(6,*) ' only dipole-related properties are available'
          call fatal
        end if
        lenprop = 0
      end if
c
c     open the matrix element file and read some information that will 
c     be kept in memory for the whole process, and that is needed to
c     write the header of the file at each step:
c
      call open_read(mat_name(1:lmname),iumat,labfil,ivers,nlab,gvers,
     $  title,nattmp,nbasis,nbsuse,icharg,multip,ne,len12l,len4l,iopcl,
     $  icgu)
      if (iumat.le.0) then
         write(6,*) 'mat name:', mat_name(1:lmname)
         write(6,*) 'len12l, len4l:', len12l, len4l
         write(6,*) 'iunit=', iumat
         write(6,*) ' failed to open matrix element file in kqmmm.'
         call fatal
      end if
c
c     read the header and close the file.
c
      call rd_head(iumat,nlab,nattmp,nbasis,ian,iattyp,atmchg,cgau,
     &  ibfatm,ibftyp,atmwgt,nfc,nfv,itran,idum9,nshao,nprao,nshdb,
     &  nprdb,nbtot)
      call close_matf(.false.,iumat)
      return
      end
