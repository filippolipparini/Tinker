c
c
c     ###############################################################
c     ##  COPYRIGHT (C)  2000  by Filippo Lipparini, Daniele Loco, ##
c     ##  Louis Lagardere, Michele Nottoli, Jean-Philip Piquemal,  ##
c     ##  Benedetta Mennucci and Jay William Ponder                ##
c     ##              All Rights Reserved                          ##
c     ###############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module qmmmstuf  --  qmmm controls                         ##
c     ##                                                             ##
c     #################################################################
c
c
c     escf    scf energy at the current step
c     qmdip   qm dipole at the current step
c     maxfor  |deqmmm|_inf at the current step on the qm and mm atoms.
c
      module qmmm
      implicit none
c
c     quantities related to the gau_open i/o system:
c
      character*(64) gau_name*80, mat_name*80, scr_name*80
      character*(64) title,gvers,labfil,cbuf 
      integer ivers,nlab,ne,len12l,len4l
      logical qmmm_post
      integer lscrname
      integer lgname,lmname,iumat,nbasis,nbsuse,icharg,multip,iopcl
      integer icgu,nfc,nfv,itran,idum9,nshao,nprao,nshdb,nprdb,nbtot
      integer, parameter   :: maxbas=20000
      integer, allocatable :: ian(:),iattyp(:),ibfatm(:),ibftyp(:)
      real*8,  allocatable :: atmchg(:),cgau(:,:),atmwgt(:)
!
!     xlbo stuff
!
      logical  xlbo, called
      integer  nmat, nbasxl, nttxl
      real*8,  allocatable :: pguess(:,:), pconv(:), guess(:)
c
c     properties
c
      integer nprops, nstates, idens(100), iprop(100)
      real*8 escf, etd
      real*8 qmdip(3), tdip(3,100), gen(1000)
      real*8 maxfor(2)
      end
