#
#
#  ################################################################
#  ##                                                            ##
#  ##  debug.make  --  compile Tinker routines with debug flags  ##
#  ##              (GNU gfortran for Linux Version)              ##
#  ##                                                            ##
#  ################################################################
#
#
#  compile all the modules; "sizes" must be first since it is used
#  to set static array dimensions in many of the other modules
#
#
gfortran -c -Wall sizes.f
gfortran -c -Wall action.f
gfortran -c -Wall align.f
gfortran -c -Wall analyz.f
gfortran -c -Wall angang.f
gfortran -c -Wall angbnd.f
gfortran -c -Wall angpot.f
gfortran -c -Wall angtor.f
gfortran -c -Wall argue.f
gfortran -c -Wall ascii.f
gfortran -c -Wall atmlst.f
gfortran -c -Wall atomid.f
gfortran -c -Wall atoms.f
gfortran -c -Wall bath.f
gfortran -c -Wall bitor.f
gfortran -c -Wall bndpot.f
gfortran -c -Wall bndstr.f
gfortran -c -Wall bound.f
gfortran -c -Wall boxes.f
gfortran -c -Wall cell.f
gfortran -c -Wall cflux.f
gfortran -c -Wall charge.f
gfortran -c -Wall chgpen.f
gfortran -c -Wall chgpot.f
gfortran -c -Wall chgtrn.f
gfortran -c -Wall chrono.f
gfortran -c -Wall chunks.f
gfortran -c -Wall couple.f
gfortran -c -Wall ctrpot.f
gfortran -c -Wall deriv.f
gfortran -c -Wall dipole.f
gfortran -c -Wall disgeo.f
gfortran -c -Wall disp.f
gfortran -c -Wall dma.f
gfortran -c -Wall domega.f
gfortran -c -Wall dsppot.f
gfortran -c -Wall energi.f
gfortran -c -Wall ewald.f
gfortran -c -Wall faces.f
gfortran -c -Wall fft.f
gfortran -c -Wall fields.f
gfortran -c -Wall files.f
gfortran -c -Wall fracs.f
gfortran -c -Wall freeze.f
gfortran -c -Wall gkstuf.f
gfortran -c -Wall group.f
gfortran -c -Wall hescut.f
gfortran -c -Wall hessn.f
gfortran -c -Wall hpmf.f
gfortran -c -Wall ielscf.f
gfortran -c -Wall improp.f
gfortran -c -Wall imptor.f
gfortran -c -Wall inform.f
gfortran -c -Wall inter.f
gfortran -c -Wall iounit.f
gfortran -c -Wall kanang.f
gfortran -c -Wall kangs.f
gfortran -c -Wall kantor.f
gfortran -c -Wall katoms.f
gfortran -c -Wall kbonds.f
gfortran -c -Wall kcflux.f
gfortran -c -Wall kchrge.f
gfortran -c -Wall kcpen.f
gfortran -c -Wall kctrn.f
gfortran -c -Wall kdipol.f
gfortran -c -Wall kdsp.f
gfortran -c -Wall keys.f
gfortran -c -Wall khbond.f
gfortran -c -Wall kiprop.f
gfortran -c -Wall kitors.f
gfortran -c -Wall kmulti.f
gfortran -c -Wall kopbnd.f
gfortran -c -Wall kopdst.f
gfortran -c -Wall korbs.f
gfortran -c -Wall kpitor.f
gfortran -c -Wall kpolr.f
gfortran -c -Wall krepl.f
gfortran -c -Wall kstbnd.f
gfortran -c -Wall ksttor.f
gfortran -c -Wall ktorsn.f
gfortran -c -Wall ktrtor.f
gfortran -c -Wall kurybr.f
gfortran -c -Wall kvdwpr.f
gfortran -c -Wall kvdws.f
gfortran -c -Wall light.f
gfortran -c -Wall limits.f
gfortran -c -Wall linmin.f
gfortran -c -Wall math.f
gfortran -c -Wall mdstuf.f
gfortran -c -Wall merck.f
gfortran -c -Wall minima.f
gfortran -c -Wall molcul.f
gfortran -c -Wall moldyn.f
gfortran -c -Wall moment.f
gfortran -c -Wall mplpot.f
gfortran -c -Wall mpole.f
gfortran -c -Wall mrecip.f
gfortran -c -Wall mutant.f
gfortran -c -Wall neigh.f
gfortran -c -Wall nonpol.f
gfortran -c -Wall nucleo.f
gfortran -c -Wall omega.f
gfortran -c -Wall opbend.f
gfortran -c -Wall opdist.f
gfortran -c -Wall openmp.f
gfortran -c -Wall orbits.f
gfortran -c -Wall output.f
gfortran -c -Wall params.f
gfortran -c -Wall paths.f
gfortran -c -Wall pbstuf.f
gfortran -c -Wall pdb.f
gfortran -c -Wall phipsi.f
gfortran -c -Wall piorbs.f
gfortran -c -Wall pistuf.f
gfortran -c -Wall pitors.f
gfortran -c -Wall pme.f
gfortran -c -Wall polar.f
gfortran -c -Wall polgrp.f
gfortran -c -Wall polopt.f
gfortran -c -Wall polpcg.f
gfortran -c -Wall polpot.f
gfortran -c -Wall poltcg.f
gfortran -c -Wall potent.f
gfortran -c -Wall potfit.f
gfortran -c -Wall ptable.f
gfortran -c -Wall qmstuf.f
gfortran -c -Wall refer.f
gfortran -c -Wall repel.f
gfortran -c -Wall reppot.f
gfortran -c -Wall resdue.f
gfortran -c -Wall restrn.f
gfortran -c -Wall rgddyn.f
gfortran -c -Wall rigid.f
gfortran -c -Wall ring.f
gfortran -c -Wall rotbnd.f
gfortran -c -Wall rxnfld.f
gfortran -c -Wall rxnpot.f
gfortran -c -Wall scales.f
gfortran -c -Wall sequen.f
gfortran -c -Wall shunt.f
gfortran -c -Wall socket.f
gfortran -c -Wall solute.f
gfortran -c -Wall stodyn.f
gfortran -c -Wall strbnd.f
gfortran -c -Wall strtor.f
gfortran -c -Wall syntrn.f
gfortran -c -Wall tarray.f
gfortran -c -Wall titles.f
gfortran -c -Wall torpot.f
gfortran -c -Wall tors.f
gfortran -c -Wall tortor.f
gfortran -c -Wall tree.f
gfortran -c -Wall units.f
gfortran -c -Wall uprior.f
gfortran -c -Wall urey.f
gfortran -c -Wall urypot.f
gfortran -c -Wall usage.f
gfortran -c -Wall valfit.f
gfortran -c -Wall vdw.f
gfortran -c -Wall vdwpot.f
gfortran -c -Wall vibs.f
gfortran -c -Wall virial.f
gfortran -c -Wall warp.f
gfortran -c -Wall xtals.f
gfortran -c -Wall zclose.f
gfortran -c -Wall zcoord.f
#
#  now compile separately each of the Fortran source files
#
gfortran -c -Wall active.f
gfortran -c -Wall alchemy.f
gfortran -c -Wall alterchg.f
gfortran -c -Wall analysis.f
gfortran -c -Wall analyze.f
gfortran -c -Wall angles.f
gfortran -c -Wall anneal.f
gfortran -c -Wall archive.f
gfortran -c -Wall attach.f
gfortran -c -Wall baoab.f
gfortran -c -Wall bar.f
gfortran -c -Wall basefile.f
gfortran -c -Wall beeman.f
gfortran -c -Wall bicubic.f
gfortran -c -Wall bitors.f
gfortran -c -Wall bonds.f
gfortran -c -Wall born.f
gfortran -c -Wall bounds.f
gfortran -c -Wall bussi.f
gfortran -c -Wall calendar.f
gfortran -c -Wall center.f
gfortran -c -Wall chkpole.f
gfortran -c -Wall chkring.f
gfortran -c -Wall chkxyz.f
gfortran -c -Wall cholesky.f
gfortran -c -Wall clock.f
gfortran -c -Wall cluster.f
gfortran -c -Wall column.f
gfortran -c -Wall command.f
gfortran -c -Wall connect.f
gfortran -c -Wall connolly.f
gfortran -c -Wall control.f
gfortran -c -Wall correlate.f
gfortran -c -Wall crystal.f
gfortran -c -Wall cspline.f
gfortran -c -Wall cutoffs.f
gfortran -c -Wall damping.f
gfortran -c -Wall dcflux.f
gfortran -c -Wall deflate.f
gfortran -c -Wall delete.f
gfortran -c -Wall diagq.f
gfortran -c -Wall diffeq.f
gfortran -c -Wall diffuse.f
gfortran -c -Wall distgeom.f
gfortran -c -Wall document.f
gfortran -c -Wall dynamic.f
gfortran -c -Wall eangang.f
gfortran -c -Wall eangang1.f
gfortran -c -Wall eangang2.f
gfortran -c -Wall eangang3.f
gfortran -c -Wall eangle.f
gfortran -c -Wall eangle1.f
gfortran -c -Wall eangle2.f
gfortran -c -Wall eangle3.f
gfortran -c -Wall eangtor.f
gfortran -c -Wall eangtor1.f
gfortran -c -Wall eangtor2.f
gfortran -c -Wall eangtor3.f
gfortran -c -Wall ebond.f
gfortran -c -Wall ebond1.f
gfortran -c -Wall ebond2.f
gfortran -c -Wall ebond3.f
gfortran -c -Wall ebuck.f
gfortran -c -Wall ebuck1.f
gfortran -c -Wall ebuck2.f
gfortran -c -Wall ebuck3.f
gfortran -c -Wall echarge.f
gfortran -c -Wall echarge1.f
gfortran -c -Wall echarge2.f
gfortran -c -Wall echarge3.f
gfortran -c -Wall echgdpl.f
gfortran -c -Wall echgdpl1.f
gfortran -c -Wall echgdpl2.f
gfortran -c -Wall echgdpl3.f
gfortran -c -Wall echgtrn.f
gfortran -c -Wall echgtrn1.f
gfortran -c -Wall echgtrn2.f
gfortran -c -Wall echgtrn3.f
gfortran -c -Wall edipole.f
gfortran -c -Wall edipole1.f
gfortran -c -Wall edipole2.f
gfortran -c -Wall edipole3.f
gfortran -c -Wall edisp.f
gfortran -c -Wall edisp1.f
gfortran -c -Wall edisp2.f
gfortran -c -Wall edisp3.f
gfortran -c -Wall egauss.f
gfortran -c -Wall egauss1.f
gfortran -c -Wall egauss2.f
gfortran -c -Wall egauss3.f
gfortran -c -Wall egeom.f
gfortran -c -Wall egeom1.f
gfortran -c -Wall egeom2.f
gfortran -c -Wall egeom3.f
gfortran -c -Wall ehal.f
gfortran -c -Wall ehal1.f
gfortran -c -Wall ehal2.f
gfortran -c -Wall ehal3.f
gfortran -c -Wall eimprop.f
gfortran -c -Wall eimprop1.f
gfortran -c -Wall eimprop2.f
gfortran -c -Wall eimprop3.f
gfortran -c -Wall eimptor.f
gfortran -c -Wall eimptor1.f
gfortran -c -Wall eimptor2.f
gfortran -c -Wall eimptor3.f
gfortran -c -Wall elj.f
gfortran -c -Wall elj1.f
gfortran -c -Wall elj2.f
gfortran -c -Wall elj3.f
gfortran -c -Wall embed.f
gfortran -c -Wall emetal.f
gfortran -c -Wall emetal1.f
gfortran -c -Wall emetal2.f
gfortran -c -Wall emetal3.f
gfortran -c -Wall emm3hb.f
gfortran -c -Wall emm3hb1.f
gfortran -c -Wall emm3hb2.f
gfortran -c -Wall emm3hb3.f
gfortran -c -Wall empole.f
gfortran -c -Wall empole1.f
gfortran -c -Wall empole2.f
gfortran -c -Wall empole3.f
gfortran -c -Wall energy.f
gfortran -c -Wall eopbend.f
gfortran -c -Wall eopbend1.f
gfortran -c -Wall eopbend2.f
gfortran -c -Wall eopbend3.f
gfortran -c -Wall eopdist.f
gfortran -c -Wall eopdist1.f
gfortran -c -Wall eopdist2.f
gfortran -c -Wall eopdist3.f
gfortran -c -Wall epitors.f
gfortran -c -Wall epitors1.f
gfortran -c -Wall epitors2.f
gfortran -c -Wall epitors3.f
gfortran -c -Wall epolar.f
gfortran -c -Wall epolar1.f
gfortran -c -Wall epolar2.f
gfortran -c -Wall epolar3.f
gfortran -c -Wall erepel.f
gfortran -c -Wall erepel1.f
gfortran -c -Wall erepel2.f
gfortran -c -Wall erepel3.f
gfortran -c -Wall erf.f
gfortran -c -Wall erxnfld.f
gfortran -c -Wall erxnfld1.f
gfortran -c -Wall erxnfld2.f
gfortran -c -Wall erxnfld3.f
gfortran -c -Wall esolv.f
gfortran -c -Wall esolv1.f
gfortran -c -Wall esolv2.f
gfortran -c -Wall esolv3.f
gfortran -c -Wall estrbnd.f
gfortran -c -Wall estrbnd1.f
gfortran -c -Wall estrbnd2.f
gfortran -c -Wall estrbnd3.f
gfortran -c -Wall estrtor.f
gfortran -c -Wall estrtor1.f
gfortran -c -Wall estrtor2.f
gfortran -c -Wall estrtor3.f
gfortran -c -Wall etors.f
gfortran -c -Wall etors1.f
gfortran -c -Wall etors2.f
gfortran -c -Wall etors3.f
gfortran -c -Wall etortor.f
gfortran -c -Wall etortor1.f
gfortran -c -Wall etortor2.f
gfortran -c -Wall etortor3.f
gfortran -c -Wall eurey.f
gfortran -c -Wall eurey1.f
gfortran -c -Wall eurey2.f
gfortran -c -Wall eurey3.f
gfortran -c -Wall evcorr.f
gfortran -c -Wall extra.f
gfortran -c -Wall extra1.f
gfortran -c -Wall extra2.f
gfortran -c -Wall extra3.f
gfortran -c -Wall fatal.f
gfortran -c -Wall fft3d.f
gfortran -c -Wall fftpack.f
gfortran -c -Wall field.f
gfortran -c -Wall final.f
gfortran -c -Wall flatten.f
gfortran -c -Wall freeunit.f
gfortran -c -Wall gda.f
gfortran -c -Wall geometry.f
gfortran -c -Wall getarc.f
gfortran -c -Wall getint.f
gfortran -c -Wall getkey.f
gfortran -c -Wall getmol.f
gfortran -c -Wall getmol2.f
gfortran -c -Wall getnumb.f
gfortran -c -Wall getpdb.f
gfortran -c -Wall getprm.f
gfortran -c -Wall getref.f
gfortran -c -Wall getstring.f
gfortran -c -Wall gettext.f
gfortran -c -Wall getword.f
gfortran -c -Wall getxyz.f
gfortran -c -Wall ghmcstep.f
gfortran -c -Wall gradient.f
gfortran -c -Wall gradrgd.f
gfortran -c -Wall gradrot.f
gfortran -c -Wall groups.f
gfortran -c -Wall grpline.f
gfortran -c -Wall gyrate.f
gfortran -c -Wall hessian.f
gfortran -c -Wall hessrgd.f
gfortran -c -Wall hessrot.f
gfortran -c -Wall hybrid.f
gfortran -c -Wall image.f
gfortran -c -Wall impose.f
gfortran -c -Wall induce.f
gfortran -c -Wall inertia.f
gfortran -c -Wall initatom.f
gfortran -c -Wall initial.f
gfortran -c -Wall initprm.f
gfortran -c -Wall initres.f
gfortran -c -Wall initrot.f
gfortran -c -Wall insert.f
gfortran -c -Wall intedit.f
gfortran -c -Wall intxyz.f
gfortran -c -Wall invbeta.f
gfortran -c -Wall invert.f
gfortran -c -Wall jacobi.f
gfortran -c -Wall kangang.f
gfortran -c -Wall kangle.f
gfortran -c -Wall kangtor.f
gfortran -c -Wall katom.f
gfortran -c -Wall kbond.f
gfortran -c -Wall kcharge.f
gfortran -c -Wall kchgflx.f
gfortran -c -Wall kchgtrn.f
gfortran -c -Wall kdipole.f
gfortran -c -Wall kdisp.f
gfortran -c -Wall kewald.f
gfortran -c -Wall kextra.f
gfortran -c -Wall kgeom.f
gfortran -c -Wall kimprop.f
gfortran -c -Wall kimptor.f
gfortran -c -Wall kinetic.f
gfortran -c -Wall kmetal.f
gfortran -c -Wall kmpole.f
gfortran -c -Wall kopbend.f
gfortran -c -Wall kopdist.f
gfortran -c -Wall korbit.f
gfortran -c -Wall kpitors.f
gfortran -c -Wall kpolar.f
gfortran -c -Wall krepel.f
gfortran -c -Wall ksolv.f
gfortran -c -Wall kstrbnd.f
gfortran -c -Wall kstrtor.f
gfortran -c -Wall ktors.f
gfortran -c -Wall ktortor.f
gfortran -c -Wall kurey.f
gfortran -c -Wall kvdw.f
gfortran -c -Wall lattice.f
gfortran -c -Wall lbfgs.f
gfortran -c -Wall lights.f
gfortran -c -Wall makeint.f
gfortran -c -Wall makeref.f
gfortran -c -Wall makexyz.f
gfortran -c -Wall maxwell.f
gfortran -c -Wall mdinit.f
gfortran -c -Wall mdrest.f
gfortran -c -Wall mdsave.f
gfortran -c -Wall mdstat.f
gfortran -c -Wall mechanic.f
gfortran -c -Wall merge.f
gfortran -c -Wall minimize.f
gfortran -c -Wall minirot.f
gfortran -c -Wall minrigid.f
gfortran -c -Wall mol2xyz.f
gfortran -c -Wall molecule.f
gfortran -c -Wall molxyz.f
gfortran -c -Wall moments.f
gfortran -c -Wall monte.f
gfortran -c -Wall mutate.f
gfortran -c -Wall nblist.f
gfortran -c -Wall newton.f
gfortran -c -Wall newtrot.f
gfortran -c -Wall nextarg.f
gfortran -c -Wall nexttext.f
gfortran -c -Wall nose.f
gfortran -c -Wall nspline.f
gfortran -c -Wall nucleic.f
gfortran -c -Wall number.f
gfortran -c -Wall numeral.f
gfortran -c -Wall numgrad.f
gfortran -c -Wall ocvm.f
gfortran -c -Wall openend.f
gfortran -c -Wall optimize.f
gfortran -c -Wall optinit.f
gfortran -c -Wall optirot.f
gfortran -c -Wall optrigid.f
gfortran -c -Wall optsave.f
gfortran -c -Wall orbital.f
gfortran -c -Wall orient.f
gfortran -c -Wall orthog.f
gfortran -c -Wall overlap.f
gfortran -c -Wall path.f
gfortran -c -Wall pdbxyz.f
gfortran -c -Wall picalc.f
gfortran -c -Wall pmestuf.f
gfortran -c -Wall pmpb.f
gfortran -c -Wall polarize.f
gfortran -c -Wall poledit.f
gfortran -c -Wall polymer.f
gfortran -c -Wall potential.f
gfortran -c -Wall pressure.f
gfortran -c -Wall prmedit.f
gfortran -c -Wall prmkey.f
gfortran -c -Wall promo.f
gfortran -c -Wall protein.f
gfortran -c -Wall prtdyn.f
gfortran -c -Wall prterr.f
gfortran -c -Wall prtint.f
gfortran -c -Wall prtmol2.f
gfortran -c -Wall prtpdb.f
gfortran -c -Wall prtprm.f
gfortran -c -Wall prtseq.f
gfortran -c -Wall prtxyz.f
gfortran -c -Wall pss.f
gfortran -c -Wall pssrigid.f
gfortran -c -Wall pssrot.f
gfortran -c -Wall qrfact.f
gfortran -c -Wall quatfit.f
gfortran -c -Wall radial.f
gfortran -c -Wall random.f
gfortran -c -Wall rattle.f
gfortran -c -Wall readdyn.f
gfortran -c -Wall readgau.f
gfortran -c -Wall readgdma.f
gfortran -c -Wall readint.f
gfortran -c -Wall readmol.f
gfortran -c -Wall readmol2.f
gfortran -c -Wall readpdb.f
gfortran -c -Wall readprm.f
gfortran -c -Wall readseq.f
gfortran -c -Wall readxyz.f
gfortran -c -Wall replica.f
gfortran -c -Wall respa.f
gfortran -c -Wall rgdstep.f
gfortran -c -Wall rings.f
gfortran -c -Wall rmsfit.f
gfortran -c -Wall rotlist.f
gfortran -c -Wall rotpole.f
gfortran -c -Wall saddle.f
gfortran -c -Wall scan.f
gfortran -c -Wall sdstep.f
gfortran -c -Wall search.f
gfortran -c -Wall server.f
gfortran -c -Wall shakeup.f
gfortran -c -Wall sigmoid.f
gfortran -c -Wall simplex.f
gfortran -c -Wall sktstuf.f
gfortran -c -Wall sniffer.f
gfortran -c -Wall sort.f
gfortran -c -Wall spacefill.f
gfortran -c -Wall spectrum.f
gfortran -c -Wall square.f
gfortran -c -Wall suffix.f
gfortran -c -Wall superpose.f
gfortran -c -Wall surface.f
gfortran -c -Wall surfatom.f
gfortran -c -Wall switch.f
gfortran -c -Wall tcgstuf.f
gfortran -c -Wall temper.f
gfortran -c -Wall testgrad.f
gfortran -c -Wall testhess.f
gfortran -c -Wall testpair.f
gfortran -c -Wall testpol.f
gfortran -c -Wall testrot.f
gfortran -c -Wall testvir.f
gfortran -c -Wall timer.f
gfortran -c -Wall timerot.f
gfortran -c -Wall tncg.f
gfortran -c -Wall torphase.f
gfortran -c -Wall torque.f
gfortran -c -Wall torsfit.f
gfortran -c -Wall torsions.f
gfortran -c -Wall trimtext.f
gfortran -c -Wall unitcell.f
gfortran -c -Wall valence.f
gfortran -c -Wall verlet.f
gfortran -c -Wall version.f
gfortran -c -Wall vibbig.f
gfortran -c -Wall vibrate.f
gfortran -c -Wall vibrot.f
gfortran -c -Wall volume.f
gfortran -c -Wall xtalfit.f
gfortran -c -Wall xtalmin.f
gfortran -c -Wall xyzatm.f
gfortran -c -Wall xyzedit.f
gfortran -c -Wall xyzint.f
gfortran -c -Wall xyzmol2.f
gfortran -c -Wall xyzpdb.f
gfortran -c -Wall zatom.f
