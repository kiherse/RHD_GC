!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!     This subroutine calculates the fluxes in the direction of the
!     sweep.
!     If IGEOM = 3,4,5, the x coordinate will be a radial coordinate and
!     the index j will always refer to the x-direction.
!     Adapted from Prometheus.
!===============================================================================

subroutine hydrow3( tid, fznp, nznp, igodu, small, smlrho, epsiln, &
                    omg1, omg2, me, mp, x, dx, areal, p, rho,      &
                    rhoe, rhopr, vel, velt, veltt, f,              &
                    rflx, reflx, mflx, mtflx, mttflx, eflx, fflx, xyswp, j, k )

  USE workarrays
  implicit none
# include "types.h"

!------------------------------------------------------- Input Variables

  integer(is):: tid, fznp, nznp, igodu

  real(rs)   :: small, smlrho, epsiln, omg1, omg2, me, mp

  ! size -4:mn5
  real(rs)   :: x(-4:*), dx(-4:*), areal(-4:*)

!------------------------------------------------------- Output Variables

  ! size -4:mn5
  real(rs)   :: p(-4:*), rho(-4:*), rhoe(-4:*), rhopr(-4:*), vel(-4:*), &
                velt(-4:*), veltt(-4:*), f(-4:*)

  ! size -4:mn6
  real(rs)   :: rflx(-4:*), reflx(-4:*), mflx(-4:*), mtflx(-4:*), &
                mttflx(-4:*), eflx(-4:*), fflx(-4:*)

!------------------------------------------------------- Local Variables

  integer(is):: i, j, k, xyswp

  real(rs)   :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, norm

  real(rs), pointer :: scrch1(:), scrch2(:)

  real(rs), pointer :: wm(:), wp(:)

  real(rs), pointer :: rpflx(:)

  real(rs), pointer :: pl(:), pr(:), rhol(:), rhor(:), &
                       rhoel(:), rhoer(:),             &
                       rhoprl(:), rhoprr(:),           &
                       xel(:), xer(:), teml(:),        &
                       temr(:), vell(:), velr(:),      &
                       veltl(:), veltr(:),             &
                       velttl(:), velttr(:), wl(:),    &
                       wr(:), ul(:), ur(:), csl(:),    &
                       csr(:), hl(:), hr(:),           &
                       dpdrhl(:), dpdrhr(:),           &
                       dpdul(:), dpdur(:),             &
                       gaml(:), gamr(:),               &
                       rl(:), rr(:), ml(:),            &
                       mr(:), mtl(:), mtr(:),          &
                       mttl(:), mttr(:),               &
                       el(:), er(:), fl(:), fr(:)

  real(rs), pointer :: w(:), wvel(:), wvelt(:), wveltt(:)          

  real(rs), pointer :: coeff1(:), coeff2(:), coeff3(:), &
                       coeff4(:), coeff5(:)

  real(rs), pointer :: pm(:), pp(:), rhom(:), rhop(:), &
                       rhoem(:), rhoep(:),             &
                       rhopm(:), rhopp(:),             &
                       velm(:), velp(:), veltm(:),     &
                       veltp(:), velttm(:), velttp(:), &
                       fm(:), fp(:)

  real(rs), pointer :: delp(:), delrho(:), delrhoe(:), &
                       delrhop(:), delvel(:),          &
                       delvlt(:), delvltt(:),          &
                       delf(:)

  real(rs), pointer :: dp(:), p6(:), drho(:), rho6(:), &
                       drhoe(:), rhoe6(:),             &
                       drhop(:), rhop6(:),             &
                       dvel(:), vel6(:), dvelt(:),     &
                       velt6(:), dveltt(:), veltt6(:), &
                       df(:), f6(:)

  real(rs), pointer :: flatn(:), flatn1(:)

!--------------------------------------------------------------------------BEGIN

  ! -4:mn6 vectors
  rpflx  => wa(tid)%wa189
  scrch1 => wa(tid)%wa106
  scrch2 => wa(tid)%wa107

  wm => wa(tid)%wa108
  wp => wa(tid)%wa109

  pl     => wa(tid)%wa110
  pr     => wa(tid)%wa111
  rhol   => wa(tid)%wa112
  rhor   => wa(tid)%wa113
  rhoel  => wa(tid)%wa190
  rhoer  => wa(tid)%wa191
  rhoprl => wa(tid)%wa192
  rhoprr => wa(tid)%wa193
  xel    => wa(tid)%wa194
  xer    => wa(tid)%wa195
  teml   => wa(tid)%wa196
  temr   => wa(tid)%wa197
  vell   => wa(tid)%wa114
  velr   => wa(tid)%wa115
  veltl  => wa(tid)%wa116
  veltr  => wa(tid)%wa117
  velttl => wa(tid)%wa118
  velttr => wa(tid)%wa119
  wl     => wa(tid)%wa120
  wr     => wa(tid)%wa121
  ul     => wa(tid)%wa122
  ur     => wa(tid)%wa123
  csl    => wa(tid)%wa124
  csr    => wa(tid)%wa125
  hl     => wa(tid)%wa198
  hr     => wa(tid)%wa199
  dpdrhl => wa(tid)%wa1101
  dpdrhr => wa(tid)%wa1102
  dpdul  => wa(tid)%wa1103
  dpdur  => wa(tid)%wa1104
  gaml   => wa(tid)%wa1105
  gamr   => wa(tid)%wa1106
  rl     => wa(tid)%wa126
  rr     => wa(tid)%wa127
  ml     => wa(tid)%wa128
  mr     => wa(tid)%wa129
  mtl    => wa(tid)%wa130
  mtr    => wa(tid)%wa131
  mttl   => wa(tid)%wa132
  mttr   => wa(tid)%wa133
  el     => wa(tid)%wa134
  er     => wa(tid)%wa135
  fl     => wa(tid)%wa136
  fr     => wa(tid)%wa137

  ! -4:mn5 vectors
  w      => wa(tid)%wa224
  wvel   => wa(tid)%wa225
  wvelt  => wa(tid)%wa226
  wveltt => wa(tid)%wa227

  coeff1 => wa(tid)%wa228
  coeff2 => wa(tid)%wa229
  coeff3 => wa(tid)%wa230
  coeff4 => wa(tid)%wa231
  coeff5 => wa(tid)%wa232

  pm     => wa(tid)%wa233
  pp     => wa(tid)%wa234
  rhom   => wa(tid)%wa235
  rhop   => wa(tid)%wa236
  rhoem  => wa(tid)%wa237
  rhoep  => wa(tid)%wa238
  rhopm  => wa(tid)%wa239
  rhopp  => wa(tid)%wa240 
  velm   => wa(tid)%wa241
  velp   => wa(tid)%wa242
  veltm  => wa(tid)%wa243
  veltp  => wa(tid)%wa244
  velttm => wa(tid)%wa245
  velttp => wa(tid)%wa246
  fm     => wa(tid)%wa247
  fp     => wa(tid)%wa248

  delp    => wa(tid)%wa249
  delrho  => wa(tid)%wa250
  delrhoe => wa(tid)%wa251
  delrhop => wa(tid)%wa252
  delvel  => wa(tid)%wa253
  delvlt  => wa(tid)%wa254
  delvltt => wa(tid)%wa255
  delf    => wa(tid)%wa256  

  dp     => wa(tid)%wa257
  p6     => wa(tid)%wa258
  drho   => wa(tid)%wa259
  rho6   => wa(tid)%wa260
  drhoe  => wa(tid)%wa261
  rhoe6  => wa(tid)%wa262
  drhop  => wa(tid)%wa263
  rhop6  => wa(tid)%wa264
  dvel   => wa(tid)%wa265
  vel6   => wa(tid)%wa266
  dvelt  => wa(tid)%wa267
  velt6  => wa(tid)%wa268
  dveltt => wa(tid)%wa269
  veltt6 => wa(tid)%wa270
  df     => wa(tid)%wa271
  f6     => wa(tid)%wa272

  flatn  => wa(tid)%wa273
  flatn1 => wa(tid)%wa274

  !----------------------------------------------------- PPM reconstruction

  call coef( tid, fznp, nznp, dx, coeff1, coeff2, coeff3, coeff4, coeff5 )

  !Definimos las Variables wvel y wvelt para reconstruccion.
  !17/12/01. Perucho
  do i= fznp-5, nznp+5
    w(i)     = 1.0/SQRT(1.0-vel(i)*vel(i)-velt(i)*velt(i)-veltt(i)*veltt(i))
    wvel(i)  = w(i)*vel(i)
    wvelt(i) = w(i)*velt(i)
    wveltt(i)= w(i)*veltt(i)
  enddo

  ! Interpolation for finding the variables in the surfaces of the cells
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               p, pm, pp, delp  )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               rho, rhom, rhop, delrho  )
    ! Delect detect in stationary jets
    call detect( tid, fznp, nznp, smlrho, small,     &
                rho, delrho, x, dx, p, rho, rhom, rhop )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               rhoe, rhoem, rhoep, delrhoe  )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               rhopr, rhopm, rhopp, delrhop  )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               wvel, velm,  velp, delvel  )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               wvelt, veltm, veltp, delvlt  )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               wveltt, velttm, velttp, delvltt  )
  call interp( tid, fznp, nznp, small, coeff1, coeff2, coeff3, coeff4, coeff5, &
               f, fm, fp, delf  )

  ! Flattening
  call flaten( tid, fznp, nznp, igodu, epsiln, omg1, omg2, &
               p, vel, flatn, flatn1 )
  
  do i= fznp-1, nznp+1
    tt1 = flatn(i)*rho(i)
    tt7 = flatn(i)*rhoe(i)
    tt8 = flatn(i)*rhopr(i)
    tt2 = flatn(i)*wvel(i)
    tt3 = flatn(i)*wvelt(i)
    tt4 = flatn(i)*wveltt(i)
    tt5 = flatn(i)*p(i)
    tt6 = flatn(i)*f(i)

    rhom(i)  = tt1 + flatn1(i)*rhom(i)
    rhop(i)  = tt1 + flatn1(i)*rhop(i)
    rhoem(i) = tt7 + flatn1(i)*rhoem(i)
    rhoep(i) = tt7 + flatn1(i)*rhoep(i)
    rhopm(i) = tt8 + flatn1(i)*rhopm(i)
    rhopp(i) = tt8 + flatn1(i)*rhopp(i) 
    velm(i)  = tt2 + flatn1(i)*velm(i)
    velp(i)  = tt2 + flatn1(i)*velp(i)
    veltm(i) = tt3 + flatn1(i)*veltm(i)
    veltp(i) = tt3 + flatn1(i)*veltp(i)
    velttm(i)= tt4 + flatn1(i)*velttm(i)
    velttp(i)= tt4 + flatn1(i)*velttp(i)
    pm(i)    = tt5 + flatn1(i)*pm(i)
    pp(i)    = tt5 + flatn1(i)*pp(i)
    fm(i)    = tt6 + flatn1(i)*fm(i)
    fp(i)    = tt6 + flatn1(i)*fp(i)
  enddo

  call monot( tid, fznp, nznp, small, p,      pm,     pp,     dp,     p6     )
  call monot( tid, fznp, nznp, small, rho,    rhom,   rhop,   drho,   rho6   )
  call monot( tid, fznp, nznp, small, rhoe,   rhoem,  rhoep,  drhoe,  rhoe6   )
  call monot( tid, fznp, nznp, small, rhopr,  rhopm,  rhopp,  drhop,  rhop6   )  
  call monot( tid, fznp, nznp, small, wvel,   velm,   velp,   dvel,   vel6   )
  call monot( tid, fznp, nznp, small, wvelt,  veltm,  veltp,  dvelt,  velt6  )
  call monot( tid, fznp, nznp, small, wveltt, velttm, velttp, dveltt, veltt6 )
  call monot( tid, fznp, nznp, small, f,      fm,     fp,     df,     f6     )

  do i= fznp-1,nznp+1
    SCRCH1(i) = velm(i)*velm(i) + veltm(i)*veltm(i) + velttm(i)*velttm(i)
    wm(i)     = SQRT( 1.0+ SCRCH1(i) )
    tt1 = 1.0 / wm(i)

    SCRCH1(i) = velp(i)*velp(i) + veltp(i)*veltp(i) + velttp(i)*velttp(i)
    wp(i)     = SQRT( 1.0+ SCRCH1(i) )
    tt2 = 1.0 / wp(i)

    velm(i)   = velm(i)*tt1
    veltm(i)  = veltm(i)*tt1
    velp(i)   = velp(i)*tt2
    veltp(i)  = veltp(i)*tt2
    velttm(i) = velttm(i)*tt1
    velttp(i) = velttp(i)*tt2
  enddo

  do i= fznp, nznp+1
    rhol(i)  = rhop(i-1)
    rhoel(i) = rhoep(i-1)
    rhoprl(i)= rhopp(i-1)
    pl(i)    = pp(i-1)
    tt1      = pl(i) / rhol(i)
    vell(i)  = velp(i-1)
    veltl(i) = veltp(i-1)
    velttl(i)= velttp(i-1)
    wl(i)    = wp(i-1)
    xel(i)   = rhoel(i)/rhol(i)
    teml(i)  = pl(i) / (rhol(i)*( xel(i)/me + (1.0-xel(i))/mp ))   


    rhor(i)  = rhom(i)
    rhoer(i) = rhoem(i)
    rhoprr(i)= rhopm(i) 
    pr(i)    = pm(i)
    tt3      = pr(i) / rhor(i)
    velr(i)  = velm(i)
    veltr(i) = veltm(i)
    velttr(i)= velttm(i)
    wr(i)    = wm(i)
    xer(i)   = rhoer(i)/rhor(i)
    temr(i)  = pr(i) / (rhor(i)*( xer(i)/me + (1.0-xer(i))/mp ))

    if (xyswp==1 .and. mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. &
    i == 324 .and. j == 217 .and. k == 217) then
        print*,'------------------------------ Hydro x: before Synge -----------------------------'
        print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell ',i,j,k
        print*,'me, mp, rhol, xel, teml, ul: ',me, mp, rhol(i), xel(i), teml(i), ul(i)
        print*,'pl, gaml, hl, csl, dpdrhl, dpdu: ',pl(i), gaml(i), hl(i), csl(i), dpdrhl(i), dpdur(i)                                     
        print*,'-----------------------------------------------------------------------------------------'
    endif

  enddo

  call eos_synge( 1, tid, fznp, nznp+1, me, mp, rhol, xel, teml, ul,   &
                  pl, gaml, hl, csl, dpdrhl, dpdul)	
  call eos_synge( 1, tid, fznp, nznp+1, me, mp, rhor, xer, temr, ur,   &
                  pr, gamr, hr, csr, dpdrhr, dpdur)

  do i= fznp, nznp+1
    rl(i)    = rhol(i)*wl(i)
    tt2      = rl(i)*hl(i)*wl(i)
    ml(i)    = tt2*vell(i)
    mtl(i)   = tt2*veltl(i)
    mttl(i)  = tt2*velttl(i)
    el(i)    = tt2 - rl(i) - pl(i)
    fl(i)    = fp(i-1)

    rr(i)    = rhor(i)*wr(i)
    tt4      = rr(i)*hr(i)*wr(i)
    mr(i)    = tt4*velr(i)
    mtr(i)   = tt4*veltr(i)
    mttr(i)  = tt4*velttr(i)
    er(i)    = tt4 - rr(i) - pr(i)
    fr(i)    = fm(i)

    if (xyswp==1 .and. mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. &
    i == 324 .and. j == 217 .and. k == 217) then
        print*,'------------------------------ Hydro x sweep: Riemann parameters -----------------------------'
        print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell ',i,j,k
        print*,'vell, velr, veltl, veltr, velttl, velttr, wl, wr: ', vell(i), velr(i), veltl(i), veltr(i), velttl(i), velttr(i), wl(i), wr(i) 
        print*,'csl, csr, pl, pr, hl, hr, dpdrhl, dpdrhr, dpdul: ',csl(i), csr(i), pl(i), pr(i), hl(i), hr(i), dpdrhl(i), dpdrhr(i), dpdul(i)
        print*,'dpdur, rl, rr, ml, mr, mtl, mtr: ',dpdur(i), rl(i), rr(i), ml(i), mr(i), mtl(i), mtr(i)
        print*,'mttl, mttr, el, er: ',mttl(i), mttr(i), el(i), er(i)                                       
        print*,'-----------------------------------------------------------------------------------------'
    endif

  enddo

  ! Approximate Riemann solver
   call riem3dp( tid, fznp, nznp+1, rhol, rhor, ul, ur,            &
                 vell, velr, veltl, veltr, velttl, velttr, wl, wr, &
                 csl, csr, pl, pr, hl, hr, dpdrhl, dpdrhr, dpdul,  &
                 dpdur, rl, rr, ml, mr, mtl, mtr,                  &
                 mttl, mttr, el, er,                               &
                 rflx, mflx, mtflx, mttflx, eflx )

  do i=fznp,nznp+1
  
    if (xyswp==1 .and. mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. &
    i == 324 .and. j == 217 .and. k == 217) then
      print*,'---------- Hydro: after Riemann ------------'
      print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell ',i,j,k
      print*,'Fluxes i+1 and i: ',mflx(i+1),mflx(i) 
      print*,'dx(i): ',dx(i)                                            
      print*,'----------------------------------------------------------'
      stop
    endif

    if (rflx(i) > 0.0) then  
      reflx(i) = rhoel(i)*wl(i)*rflx(i)/rl(i)
      rpflx(i) = rhoprl(i)*wl(i)*rflx(i)/rl(i)
    else
      reflx(i) = rhoer(i)*wr(i)*rflx(i)/rr(i)
      rpflx(i) = rhoprr(i)*wr(i)*rflx(i)/rr(i)
    end if
  enddo

! renormalize fluxes:
  do i=fznp,nznp+1
    if (rflx(i) /= 0.0) then
      norm = (reflx(i)+rpflx(i))/rflx(i)
      reflx(i) = reflx(i)/norm
    endif
  enddo

  do i= fznp, nznp+1
    rflx(i)  = rflx(i)  *areal(i)
    reflx(i) = reflx(i) *areal(i)
    mflx(i)  = mflx(i)  *areal(i)
    mtflx(i) = mtflx(i) *areal(i)
    mttflx(i)= mttflx(i)*areal(i)
    eflx(i)  = eflx(i)  *areal(i)
    if (rflx(i) > 0.0) then
      fflx(i) = fl(i)*rflx(i)
    else
      fflx(i) = fr(i)*rflx(i)
    endif

  enddo

!----------------------------------------------------------------------------END

  END subroutine hydrow3

