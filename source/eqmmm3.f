c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine eqmmm3  --  qmmm with gaussian potentials  ##
c     ##                                                        ##
c     ############################################################
c
c     "eqmm3" calculates qm/mm energy and gradients with gaussian
c
      subroutine eqmmm3
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
      integer status, istat
      logical yesno,eof,asym
      real*8  dipole(3), vmax
      character*120 command
c
      lenbuf=4000
c
c     zero out the qmmm energy term
c
      eqmmm = 0.0d0
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
      call close_matf(iumat)
c
c     launch the qm/mm calculation by gaussian
c
      write(command,*) 'gdvtest ',gau_name(1:lgname)
      status = system(command)
c
c     open the matrix element file and read the energy and forces:
c
      call open_read(mat_name(1:lmname),iumat,labfil,ivers,nlab,gvers,
     $  title,nattmp,nbasis,nbsuse,icharg,multip,ne,len12l,len4l,iopcl,
     $  icgu)
c
c     if this is the first gaussian run, recover some important
c     information that is needed in order to allocate memory for
c     the xlbo procedure:
c
      if (iumat.le.0) then
         write(6,*) ' failed to open matrix element file for reading'
         call fatal
      end if
c
      if (.not. called) then
         called = .true.
         nbasxl = nbasis
         nttxl  = (nbasis*(nbasis+1))/2
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
     $    asym,nri,eof)
        lcbuf = len_trim(cbuf)
        if (cbuf(1:lcbuf).eq.'GAUSSIAN SCALARS') then
          lr = lenarr(n1,n2,n3,n4,n5)
          call rd_rind(iumat,nr,lr,ntot,lenbuf,lnz,gen) 
          escf = gen(32)*hartree
          etd  = gen(25)*hartree
          if (etd.ne.0.0d0) then
             eqmmm = escf + etd
          else
             eqmmm = escf
          endif
        else if(cbuf(1:lcbuf).eq.'ELECTRIC DIPOLE MOMENT') then
          call rd_rbuf(iumat,nr*ntot,nr*lenbuf,qmdip)
          do j = 1, 3
              qmdip(j) = qmdip(j)*debye*bohr
          end do
        end if
      end do
c
c     close the matrix elements file:
c
      call close_matf(iumat)
      return
      end
