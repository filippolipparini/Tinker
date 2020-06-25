      subroutine checkqmmm
      use fields
      use iounit
      use potent
      use usage
      implicit none
      logical ok_qmmm
c
c     only amber, amoeba, opls and charmm force field are compatible
c     with qmmm. 
c     issue an error and abort if a different force field is used.
c
      ok_qmmm = .false.
      if (use_qmmm) then
        if (forcefield(1:5).eq.'AMBER') then
          ok_qmmm = .true.
        else if (forcefield(1:6).eq.'AMOEBA') then
          ok_qmmm = forcefield(1:11).ne.'AMOEBA-PLUS'
        else if (forcefield(1:4).eq.'OPLS') then
          ok_qmmm = .true.
        else if (forcefield(1:6).eq.'CHARMM') then
          ok_qmmm = .true.
        end if
        if (.not. ok_qmmm) then
          write(6,*) ' qmmm and ', forcefield,' not supported.'
          call fatal
        end if
      end if
c
      if (nqmatoms.ne.0 .and. .not. use_qmmm) then
        write(6,*) ' qmatoms, but not qmmm'
        call fatal
      end if
      if (nqmatoms.eq.0 .and. use_qmmm) then
        write(6,*) ' no qmatoms, but qmmm'
        call fatal
      end if
c
c     check whether there are link atoms, and if these are well 
c     defined.
c
      call checklink
      return
      end
