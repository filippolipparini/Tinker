      subroutine qmproperties
      use atoms
      use qmmm
      use units
      use usage
      implicit none
c
      integer istate,ndens,inext
      integer i,iglob,j,system,nattmp,ni,nr,ntot,lenbuf
      integer n1,n2,n3,n4,n5,nri,lcbuf,lr,lenarr,lnz,ioff,k
      integer status, istat, typea
      integer ntt, nbsq, nspin
      logical yesno,eof,haveb
      logical do_pa, do_pb, do_pta, do_ptb, do_tra, do_trb
      real*8  dip_nuc(3), dip_scf(3), dip_tot(3), dip_tr(3,100)
c
      real*8, allocatable :: p_scf(:,:), p_tr(:,:,:), p_rel(:,:)
      real*8, allocatable :: dipint(:,:), dipvel(:,:), rdel(:,:)
      real*8, allocatable :: sqint(:,:)
c
      real*8 zero, two, sq2
      save zero, two
      data zero/0.0d0/, two/2.0d0/
c
c
c     open the matrix element file and gather the required integrals
c     and density matrices to compute various properties.
c
      call open_read(mat_name(1:lmname),iumat,labfil,ivers,nlab,gvers,
     $  title,nattmp,nbasis,nbsuse,icharg,multip,ne,len12l,len4l,iopcl,
     $  icgu)
c
      ntt  = (nbasis*(nbasis+1))/2
      nbsq = nbasis*nbasis
c
      call rd_head(iumat,nlab,nattmp,nbasis,ian,iattyp,atmchg,cgau,
     $  ibfatm,ibftyp,atmwgt,nfc,nfv,itran,idum9,nshao,nprao,nshdb,
     $  nprdb,nbtot)
c
      nspin = 1
      if (iopcl.eq.1) nspin = 2
      haveb = nspin.eq.2
c
      if (iumat.le.0) then
         write(6,*) ' failed to open matrix element file for reading'
         call fatal
      end if
c
      do_pa  = .false.
      do_pta = .false.
      do_tra = .false.
      do_pb  = .false.
      do_ptb = .false.
      do_trb = .false.
      do i = 1, nprops
        if (idens(i).eq.1) do_pa  = .true.
        if (idens(i).eq.2) do_pta = .true.
        if (idens(i).eq.3) do_tra = .true.
      end do
c
      if (do_pa  .and. haveb) do_pb  = .true.
      if (do_pta .and. haveb) do_ptb = .true.
      if (do_tra .and. haveb) do_trb = .true.
c
      if (do_pa)  allocate (p_scf(ntt,nspin))
      if (do_pta) allocate (p_rel(ntt,nspin))
c
      if (do_dip)  allocate (dipint(ntt,3))
      if (do_del)  allocate (dipvel(ntt,3))
      if (do_rdel) allocate (rdel(ntt,3))   
c
      if (do_tra) allocate (sqint(nbsq,3))
c
      eof = .false.
      do while(.not. eof)
        cbuf = ' '
        call rd_labl(iumat,ivers,cbuf,ni,nr,ntot,lenbuf,n1,n2,n3,n4,n5,
     $    typea,nri,eof)
        lcbuf = len_trim(cbuf)
c
c     start scanning the records for interesting information:
c
        if (cbuf(1:lcbuf).eq.'ALPHA SCF DENSITY MATRIX') then
          if (do_pa) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,p_scf(:,1))
        else if (cbuf(1:lcbuf).eq.'BETA SCF DENSITY MATRIX') then
          if (do_pb) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,p_scf(:,2))
        else if (cbuf(1:lcbuf).eq.'ALPHA CI DENSITY MATRIX' .or. 
     $           cbuf(1:lcbuf).eq.'ALPHA MP2 DENSITY MATRIX') then
          if (do_pta) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,p_rel(:,1))
        else if (cbuf(1:lcbuf).eq.'BETA CI DENSITY MATRIX' .or. 
     $           cbuf(1:lcbuf).eq.'BETA MP2 DENSITY MATRIX') then
          if (do_pta) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,p_rel(:,2))
        else if (cbuf(1:lcbuf).eq.
     $    'GROUND TO EXCITED STATE DENSITIES') then
c
c     n5 contains the number of states, so we can finally allocate this.
c     we will read alpha and beta densities in any case.
c
          if (do_tra) then
            nstates = n5
            allocate (p_tr(n1*n2,2,nstates))
            call rd_rbuf(iumat,nr*ntot,nr*lenbuf,p_tr)
          end if
c
c     now read the required integrals:
c
        else if (cbuf(1:lcbuf).eq.'DIPOLE INTEGRALS') then
          if (do_dip) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,dipint)
        else if (cbuf(1:lcbuf).eq.'DIP VEL INTEGRALS') then
          if (do_del) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,dipvel)
        else if (cbuf(1:lcbuf).eq.'R X DEL INTEGRALS') then
          if (do_del) call rd_rbuf(iumat,nr*ntot,nr*lenbuf,rdel)
        end if
      end do
c
c     close the matrix elements file:
c
      call close_matf(.false.,iumat)
c
c     if this is the first time that we pass here, compute a total lenght of the
c     properties array and allocate it.
c
      if (lenprop.eq.0) then
c
c        count how many densities we are treating:
c
         ndens = 0
         if (do_pa)  ndens = ndens + 1
         if (do_pta) ndens = ndens + 1
         if (do_tra) ndens = ndens + nstates
c
c        compute the number of properties:
c
         if (do_dip)  lenprop = lenprop + 3*ndens
         if (do_del)  lenprop = lenprop + 3*ndens
         if (do_rdel) lenprop = lenprop + 3*ndens
c
         allocate (proparray(lenprop))
      end if
c
c     now, compute the required properties:
c
      if (do_dip) then
c
c       scf dipole:
c
        if (do_pa) then
          dip_scf = zero
          call trace_lt(.true.,nbasis,ntt,3,p_scf(:,1),dipint,dip_scf)
          if (do_pb) then
            call trace_lt(.false.,nbasis,ntt,3,p_scf(:,2),dipint,
     $        dip_scf)
          else
            dip_scf = two*dip_scf
          end if
c
          dip_scf = - dip_scf
c
c       add the nuclear contribution:
c
          dip_nuc = zero
c
          do i = 1, n
            if (qmatoms(i)) then
              dip_nuc(1) = dip_nuc(1) + atmchg(i)*x(i)/bohr
              dip_nuc(2) = dip_nuc(2) + atmchg(i)*y(i)/bohr
              dip_nuc(3) = dip_nuc(3) + atmchg(i)*z(i)/bohr
            end if
          end do
c
          write(6,100) dip_scf, dip_nuc, dip_scf+dip_nuc
 100  format(' electronic scf dipole: ',3f12.6,/,
     $       ' nuclear dipole:        ',3f12.6,/,
     $       ' total scf dipole:      ',3f12.6)
          dip_scf = dip_scf + dip_nuc
        end if
        if (do_pta) then
          call trace_lt(.true.,nbasis,ntt,3,p_rel(:,1),dipint,dip_tot)
          if (do_pb) then
            call trace_lt(.false.,nbasis,ntt,3,p_rel(:,2),dipint,
     $        dip_tot)
          else
            dip_tot = two*dip_tot
          end if
          dip_tot = - dip_tot
c
c       add the nuclear contribution:
c
          dip_nuc = zero
c
          do i = 1, n
            if (qmatoms(i)) then
              dip_nuc(1) = dip_nuc(1) + atmchg(i)*x(i)/bohr
              dip_nuc(2) = dip_nuc(2) + atmchg(i)*y(i)/bohr
              dip_nuc(3) = dip_nuc(3) + atmchg(i)*z(i)/bohr
            end if
          end do
c
          write(6,110) dip_tot, dip_nuc, dip_tot+dip_nuc
 110  format(/,' electronic dipole:     ',3f12.6,/,
     $       ' nuclear dipole:        ',3f12.6,/,
     $       ' total dipole:          ',3f12.6,/)
          dip_tot = dip_tot + dip_nuc
        end if
        if (do_tra) then
c
c         expand the dipole integrals to a square matrix:
c
          do i = 1, 3
            call lt2sq(nbasis,dipint(:,i),sqint(:,i))
          end do
          sq2 = sqrt(two)
          do istate = 1, nstates
            do i = 1, 3
              dip_tr(i,istate) = 
     $          - (dot_product(sqint(:,i),p_tr(:,1,istate)) + 
     $             dot_product(sqint(:,i),p_tr(:,2,istate)))/sq2
            end do
            write(6,120) istate, dip_tr(:,istate)
            tdip(:,istate) = dip_tr(:,istate)
          end do
 120  format(' transition dipole for state ',i3,': ',3f12.6)
        end if
c       ...
      end if
c
c     fill the array:
c
      inext = 0
      if (do_dip) then
        if (do_pa) then
          do i = 1, 3
            proparray(inext+i) = dip_scf(i)
          end do
          inext = inext + 3
        end if
        if (do_pta) then
          do i = 1, 3
            proparray(inext+i) = dip_tot(i)
          end do
          inext = inext + 3
        end if
        if (do_tra) then
          do istate = 1, nstates
            do i = 1, 3
              proparray(inext+i) = dip_tr(i,istate)
            end do
            inext = inext + 3
          end do
        end if
      end if
      return
      end
c
      subroutine trace_lt(init,nbasis,ntt,nints,p,x,prop)
      implicit none
      logical init
      integer nbasis, ntt, nints
      real*8  p(ntt), x(ntt,nints), prop(nints)
c
      integer iint, i, j, ij
      real*8  tr1, tr2
c
      if (init) prop = 0.0d0
c
      do iint = 1, nints
        tr1 = 0.0d0
        tr2 = 0.0d0
        do i = 1, ntt
          tr1 = tr1 + x(i,iint)*p(i)
        end do
        ij = 0
        do i = 1, nbasis
          ij = ij + i
          tr2 = tr2 + x(ij,iint)*p(ij)
        end do
        prop(iint) = tr1 + tr1 - tr2
      end do
      return
      end
c
      subroutine lt2sq(n,lt,sq)
      implicit none
      integer n
      real*8  lt(*), sq(n,*)
c
      integer i, j, ij
c
      ij = 0
      do i = 1, n
        do j = 1, i
          ij = ij + 1
          sq(i,j) = lt(ij)
          sq(j,i) = sq(i,j)
        end do
      end do
      return
      end
