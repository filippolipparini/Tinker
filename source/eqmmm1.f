c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine eqmmm1  --  qmmm with gaussian potentials  ##
c     ##                                                        ##
c     ############################################################
c
c     "eqmm1" calculates qm/mm energy and gradients with gaussian
c
      subroutine eqmmm1
      use atoms
      use atomid
      use deriv
      use energi
      use qmmm
      use units
      use usage
      implicit none
      integer i,iglob,j,system,nattmp,ni,nr,ntot,lenbuf
      integer n1,n2,n3,n4,n5,nri,lcbuf,lr,lenarr,lnz,ioff,k
      integer status, istat, typea
      logical yesno,eof
      real*8, allocatable :: temp(:,:)
      real*8  dipole(3), vmax, fac
      character*120 command
c
      lenbuf=4000
      allocate (temp(3,n))
      temp = 0.0d0
c
c     zero out the qmmm energy term and first derivatives
c
      eqmmm = 0.0d0
      do i = 1, n
         deqmmm(1,i) = 0.0d0
         deqmmm(2,i) = 0.0d0
         deqmmm(3,i) = 0.0d0
      end do
c
c     create the new (empty) matrix element file and write the coordinates
c
      call open_write(mat_name(1:lmname),iumat,labfil,gvers,title,n,
     $  nbasis,nbsuse,icharg,multip,ne,iopcl,icgu)
c
      if (iumat.lt.1) then
        write(6,*) ' failed to open matrix elements file for writing'
        call fatal
      end if
c
c     update the coordinates and write the header:
c
      do i = 1, n
        cgau(1,i) = x(i)/bohr
        cgau(2,i) = y(i)/bohr
        cgau(3,i) = z(i)/bohr
      end do
c
      call wr_head(iumat,n,3*n,nbasis,ian,iattyp,atmchg,cgau,ibfatm,
     $  ibftyp,atmwgt,nfc,nfv,itran,idum9,nshao,nprao,nshdb,nprdb,nbtot)
c
      if (xlbo) call xlboextr
c
      if (xlbo .and. nmat.ge.8) then
        call wr_lrbuf(iumat,'ALPHA OAO DENSITY MATRIX',1,lenbuf,-nbasis,
     $    nbasis,0,0,0,.false.,guess)
      end if
      call close_matf(.true.,iumat)
c
c     launch the qm/mm calculation by gaussian
c
      write(command,*) 'gdvtest ',gau_name(1:lgname)
      if (xlbo .and. nmat.ge.9) 
     $  write(command,*) 'gdvtest ',gau_name(1:lgname-4),'_xlbo.com'
      status = system(command)
c
c     if the execution of a post-processing script has been required,
c     run it now:
c
      if (qmmm_post) status = system(scr_name(1:lscrname))
c
c     open the matrix element file and read the energy and forces:
c
      call open_read(mat_name(1:lmname),iumat,labfil,ivers,nlab,gvers,
     $  title,nattmp,nbasis,nbsuse,icharg,multip,ne,len12l,len4l,iopcl,
     $  icgu)
c
      if (iumat.le.0) then
         write(6,*) ' failed to open matrix element file for reading'
         call fatal
      end if
c
c     if this is the first gaussian run, recover some important
c     information that is needed in order to allocate memory for
c     the xlbo procedure:
c
      if (.not. called) then
         called = .true.
         nbasxl = nbasis
         nttxl  = (nbasis*(nbasis+1))/2
c
c        allocate memory for the xlbo extrapolation:
c
         allocate (pguess(nttxl,8), pconv(nttxl), guess(nttxl),
     $     stat=istat)
         if (istat.ne.0) then
           write(6,*) ' failed to allocate memory for xlbo.'
           call fatal
         end if
      end if
c
      call rd_head(iumat,nlab,nattmp,nbasis,ian,iattyp,atmchg,cgau,
     $  ibfatm,ibftyp,atmwgt,nfc,nfv,itran,idum9,nshao,nprao,nshdb,
     $  nprdb,nbtot)
c
c     read the force/energy file for QM produced by gaussian
c
      eof = .false.
      do while(.not. eof)
        cbuf = ' '
        call rd_labl(iumat,ivers,cbuf,ni,nr,ntot,lenbuf,n1,n2,n3,n4,n5,
     $    typea,nri,eof)
        lcbuf = len_trim(cbuf)
        if (cbuf(1:lcbuf).eq.'GAUSSIAN SCALARS') then
          lr = lenarr(n1,n2,n3,n4,n5)
          call rd_rind(iumat,nr,lr,ntot,lenbuf,gen) 
          escf = gen(32)*hartree
          etd  = gen(25)*hartree
          if (etd.ne.0.0d0) then
             eqmmm = escf + etd
          else
             eqmmm = escf
          endif
        else if(cbuf(1:lcbuf).eq.'NUCLEAR GRADIENT') then
          call rd_rbuf(iumat,nr*ntot,nr*lenbuf,temp)
          maxfor(1) = 0.0d0
          maxfor(2) = 0.0d0
c
          do i = 1, n
            if (qmatoms(i)) then
              do j = 1, 3
                maxfor(1) = max(maxfor(1),abs(temp(j,i)))
              end do
            else
              do j = 1, 3
                maxfor(2) = max(maxfor(2),abs(temp(j,i)))
              end do
            end if
          end do
        else if(cbuf(1:lcbuf).eq.'ELECTRIC DIPOLE MOMENT') then
          call rd_rbuf(iumat,nr*ntot,nr*lenbuf,qmdip)
          do j = 1, 3
              qmdip(j) = qmdip(j)*debye*bohr
          end do
        else if(cbuf(1:lcbuf).eq.'ALPHA OAO DENSITY MATRIX') then
          if (xlbo) then
            nmat = nmat + 1
            call rd_rbuf(iumat,nr*ntot,nr*lenbuf,pconv)
            if (nmat.le.8) pguess(:,nmat) = pconv
          end if
        end if
      end do
c
c     close the matrix elements file:
c
      call close_matf(.false.,iumat)
c
c     read the force/energy file for MM produced by gaussian
c
      do i = 1, n
        do j = 1, 3
          deqmmm(j,i) = temp(j,i)*hartree/bohr
        end do
      end do
c
      deallocate (temp)
c
c     if required, compute some properties:
c
      if (nprops.ne.0) call qmproperties
      return
      end
