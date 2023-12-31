c 23 Aug 2019, tcr not used in piriform.f, so this program acts as
c a placeholder
c 12 Sept 2006, start with /isoldeepVFOK/integrate_tcrx.f & add GABA-B
! 17 Feb. 2003: TCR cell looks much better after adjusting gCaT distrib. per
! Williams & Stuart (peak in dendritic stems); also adjust gAR, which has
! signif. effect on number of fast spikes/LTS.
! 13 Feb. 2003: Diego Contreras felt old TCR cell showed too much attenuation of spikes
! in a burst, and too many spikes/burst at 10 Hz.  Following changes undertaken:
! decrease gCaT density to 5 from 7.5 in layers 3, 4
! shift mnaf -3 mV and hnaf -7 mV
! gKA density x 0.2 everywhere
! gKDR density x 0.45 axon, soma, level 2.

! Integration program for tcr cells. 
! From dienf.                  

c      SUBROUTINE integrate_tcrxB (O, time, numcell, V, curr,
       SUBROUTINE integrate_placeholder5 (O, time, numcell, V, curr,
     &  initialize, firstcell, lastcell,
     &  gAMPA, gNMDA, gGABA_A, gGABA_B,
     &  Mg, gapcon, totaxgj, gjtable, dt,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar)

           SAVE

       integer, parameter:: numcomp = 137  ! should be compat. with calling prog

       integer numcell, totaxgj, gjtable(totaxgj,4)
       integer initialize, firstcell, lastcell
       INTEGER  I, J, L, L1, O
c numcell = total number of tcr cells in system;
c L = cell number relative to entire system
        real*8 Mg, gapcon, dt, time

c CINV is 1/C, i.e. inverse capacitance
      real*8 mnap(numcomp,numcell), persistentNa_shift, fastNa_shift
      real*8 v(numcomp,numcell), chi(numcomp,numcell), 
     x cinv(numcomp), c(numcomp),
     x mnaf(numcomp,numcell), hnaf(numcomp,numcell),
     x mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),
     x mk2(numcomp,numcell),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell),betchi(numcomp),
     x mar(numcomp,numcell),jacob(numcomp,numcomp),
     x gam(0:numcomp,0:numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x gampa(numcomp,numcell),gnmda(numcomp,numcell),
     x curr(numcomp,numcell), ggaba_b(numcomp,numcell),
     x ggaba_a(numcomp,numcell),cafor(numcomp)
       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)

c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),fmnaf(numcomp),
     x fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmk2(numcomp),
     x fhk2(numcomp),fmnap(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),
     x fmcal(numcomp),fmar(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp), 
     x  dfv_dmnaf(numcomp),
     x  dfv_dmnap(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),
     x dfhnaf_dhnaf(numcomp),
     x dfmnap_dmnap(numcomp), dfmnap_dv(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),
     x dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),
     x dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),
     x dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),
     x dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2, outrcd(20)

      REAL*8 vL(numcomp),vk(numcomp)
      REAL*8 vna,var,vca,vgaba_a,Z,Z1,Z2
         INTEGER  K0, K1, K2, NEIGH(numcomp,11), NNUM(numcomp)
       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp)
       REAL*8 gca_tot(numcomp),gar_tot(numcomp)
       REAL*8 gca_high(numcomp), A, BB1, BB2
c this will be gCa conductance corresponding to high-thresh channels

c        if (O.eq.1) then
         if (initialize.eq.0) then

c Program assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

       CALL   TCR_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL TCRMAJ (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)

          do i = 1, numcomp
             cinv(i) = 1.d0 / c(i)
          end do

        vL = -70.d0
        vK = -95.d0
        VNA = 50.d0
        VCA = 125.d0
        VAR = -43.d0
        VAR = -35.d0
c -43 mV from Huguenard & McCormick
        VGABA_A = -81.d0

c ? initialize membrane state variables?
         do i = 1, numcomp
         do j = 1, numcell
        v(i,j) = VL(i)
         end do
         end do

         chi = 0.d0

        k1 = idnint (4.d0 * (-85.d0 + 120.d0))
c       k1 = idnint (4.d0 * (v(1) + 120.d0))

       mnaf = 0.d0
       mnap = 0.d0
       mkdr = 0.d0
       mka = 0.d0
       mk2 = 0.d0
       mkm = 0.d0
       mkc = 0.d0
       mkahp = 0.d0
       mcat = 0.d0
       mcal = 0.d0
       open = 0.d0

      hnaf = alphah_naf(k1)/(alphah_naf(k1)+betah_naf(k1))
      hka = alphah_ka(k1)/(alphah_ka(k1)+betah_ka(k1))
      hk2 = alphah_k2(k1)/(alphah_k2(k1)+betah_k2(k1))
      hcat=alphah_cat(k1)/(alphah_cat(k1)+betah_cat(k1))
c     mar=alpham_ar(k1)/(alpham_ar(k1)+betam_ar(k1))
      mar= .25d0

c                 gnaf = 0.d0
c                 gnap = 0.d0
c                 gkdr = 0.d0
c                 gka = 0.d0
c                 gk2 = 0.d0
                  gkm = 0.d0
c                 gkc = 0.d0
                  gkahp = 0.d0
c                 gcat = 0.d0
c                 gcaL = 0.d0
c                 gar = 0.d0

! 13 Sept. 2006: for messing about with gCaT
            do i = 1, 137
             gCaT(i) = 1.d0 * gCaT(i)
            end do

! 13 Sept. 2006: for messing about with h current
            do i = 1, 137
             gAR(i) = 1.d0 * gAR(i)
            end do

            do i = 1, 137
             gKA(i) = 0.2d0 * gKA(i)
            end do

            do i = 132, 137 ! axon
             gKDR(i) = 0.45d0 * gKDR(i)
            end do
            gKDR(1) = 0.45d0 * gKDR(1)
            do i = 2, 119, 13  ! level 2
             gKDR(i) = 0.45d0 * gKDR(i)
            end do

c End initialization
            goto 6000

                endif


c          do 4000, L = 1, numcell    
           do 4000, L = firstcell, lastcell

       DO 301, I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL(i)) * cinv(i)
          DO 302, J = 1, NNUM(I)
             K = NEIGH(I,J)
302   FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
301    CONTINUE


        CALL FNMDA (V, OPEN, numcell, numcomp, MG, L, 
     &    A, BB1, BB2)

      DO 421, I = 1, numcomp
       FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) 
     X   - ggaba_b(I,L)*(V(I,L)-VK(i)  ) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV
421      continue

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m

       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
       end do

      DO 88, I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
c    x     gnap(i) * (mnaf(i,L)**3)
     x     gnap(i) * mnap(i,L)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


       FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK(i))
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
88           continue

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + ggaba_b(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

           do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK(i))
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L)       ) * (v(i,L) - vna)
        dfv_dmnap(i) = - cinv(i) *
     X    (               gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK(i))
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK(i))
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK(i))
       dfv_dmk2(i) = - cinv(i) * gk2(i) * hk2(i,L) * (v(i,L)-vK(i))
       dfv_dhk2(i) = - cinv(i) * gk2(i) * mk2(i,L) * (v(i,L)-vK(i))
        dfv_dmkm(i)  = - cinv(i) * gkm(i) * (v(i,L) - vK(i))
      dfv_dmkc(i) = - cinv(i) * gkc(i) * gamma(i) * (v(i,L)-vK(i))
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i,L) - vK(i))
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
            end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, 137
        alpham_ahp(i) = dmin1(0.2d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
          alpham_ahp_prime(i) = 0.2d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, 137
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
     x                  -.001d0 * mkahp(i,L)
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp

       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

       persistentNa_shift = 10.d0
       K2 = IDNINT ( 4.d0 * (V(I,L)+persistentNa_shift+ 120.d0))
       IF (K2.GT.640) K2 = 640
       IF (K2.LT.  0) K2 =   0

             fastNa_shift = -2.5d0
       K0 = IDNINT ( 4.d0 * (V(I,L)+    fastNa_shift+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0

        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fmnap(i) = alpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X              betam_naf(k2) * mnap(i,L)
        fhnaf(i) = alphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k1) * hnaf(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L))-
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfmnap_dv(i) = dalpham_naf(k2) * (1.d0 - mnap(i,L))-
     X                  dbetam_naf(k2) * mnap(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L))-
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L))-
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L))-
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L))-
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L))-
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfmnap_dmnap(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, 137
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, 137
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        mnap(i,L) = mnap(i,L) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
         end do

4000         CONTINUE
c all tcr cells on this node integrated

6000          END

C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE TCR_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr, shift_mnaf,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble (I)
          V = (V / 4.d0) - 120.d0

c gNa
           shift_mnaf = -3.d0
           V = V + shift_mnaf
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)
           V = V - shift_mnaf

            shift_hnaf =  -7.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do 2, i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
2      CONTINUE
       END

        SUBROUTINE TCRMAJ
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 10 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

        integer, parameter:: numcomp = 137

        REAL*8 C(numcomp),GL(numcomp),
     &   GAM(0:numcomp,0:numcomp),GNAF(numcomp),GCAT(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp)
        REAL*8 GKC(numcomp),GKAHP(numcomp),GCAL(numcomp)
        REAL*8 GK2(numcomp),GKM(numcomp)
        REAL*8 GNAP(numcomp),GAR(numcomp), CDENS, RM_AXON
        REAL*8 JACOB(numcomp,numcomp),RI_SD,RI_AXON,RM_SD
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:4), GCAT_DENS(0:4), GKDR_DENS(0:4)
        REAL*8 GKA_DENS(0:4), GKC_DENS(0:4), GKAHP_DENS(0:4)
        REAL*8 GCAL_DENS(0:4), GK2_DENS(0:4), GKM_DENS(0:4)
        REAL*8 GNAP_DENS(0:4), GAR_DENS(0:4)
        REAL*8 RES, RINPUT, ELEN(numcomp)
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp),LEN(numcomp),GAM1,GAM2
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI, Z
        INTEGER NEIGH(numcomp,11), NNUM(numcomp)
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS

        RI_SD = 175.d0
        RM_SD = 26400.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

        gnaf_dens(0) = 400.d0
        gnaf_dens(1) = 100.d0
        gnaf_dens(2) = 100.d0
        gnaf_dens(3) =   5.d0
        gnaf_dens(4) =   5.d0
c       gnaf_dens(3) =  40.d0                 
c       gnaf_dens(4) =  40.d0
c       gnaf_dens(3) =  20.d0                
c       gnaf_dens(4) =  20.d0

        gkdr_dens(0) = 400.d0
        gkdr_dens(1) =  75.d0
c       gkdr_dens(2) =  75.d0
c       gkdr_dens(3) =   2.d0
c       gkdr_dens(4) =   2.d0
        gkdr_dens(2) =  50.d0
        gkdr_dens(3) =   0.d0
        gkdr_dens(4) =   0.d0

        do i = 1, 4
          gnap_dens(i) = 0.002d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.01d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.10d0 * gnaf_dens(i)
        end do

c       gnap_dens(1) = 1.d0
c       gnap_dens(2) = 1.d0
c       gnap_dens(3) = 1.d0
c       gnap_dens(4) = 1.d0

c         gcat_dens(1) = 3.00d0
c         gcat_dens(2) = 3.00d0
c         gcat_dens(3) =  6.0d0
c         gcat_dens(4) =  6.0d0
!         gcat_dens(1) = 0.50d0
          gcat_dens(1) = 3.00d0
          gcat_dens(2) = 6.00d0
c         gcat_dens(3) =  3.0d0
c         gcat_dens(4) =  0.5d0
          gcat_dens(3) =  5.0d0
          gcat_dens(4) =  2.5d0

c         gcat_dens(3) =  7.5d0
c         gcat_dens(4) =  7.5d0

          gcal_dens(1) = 0.5d0
          gcal_dens(2) = 0.5d0
c         gcal_dens(3) = 1.0d0
c         gcal_dens(4) = 1.0d0
          gcal_dens(3) = 0.25d0
          gcal_dens(4) = 0.25d0

        gka_dens(0) = 1.d0
        gka_dens(1) = 30.d0
        gka_dens(2) = 30.d0
        gka_dens(3) =  1.d0
        gka_dens(4) =  1.d0

         gkc_dens(1) = 12.00d0
         gkc_dens(2) = 12.00d0
         gkc_dens(3) = 20.00d0
         gkc_dens(4) = 20.00d0

        do i = 1, 4
         gkm_dens(i) = 0.50d0
        end do

        gk2_dens(0) = .5d0
c       gk2_dens(1) = .5d0
c       gk2_dens(2) = .5d0
c       gk2_dens(3) = 5.d0
c       gk2_dens(4) = 5.d0
!       gk2_dens(1) =  2.d0
!       gk2_dens(2) =  2.d0
!       gk2_dens(3) =  2.d0
!       gk2_dens(4) =  2.d0
        gk2_dens(1) =  0.1d0
        gk2_dens(2) =  0.1d0
        gk2_dens(3) =  0.1d0
        gk2_dens(4) =  0.1d0

        do i = 1, 4
         gkahp_dens(i) = 0.05d0
        end do

c        gar_dens(1) = 2.5d0
c        gar_dens(2) = 2.5d0
c        gar_dens(3) =  5.0d0
c        gar_dens(4) =  5.0d0

c        gar_dens(0) = 1.50d0
c Run Q uses 1.5
c        gar_dens(1) = 1.5d0
c        gar_dens(2) = 1.5d0
c        gar_dens(3) = 1.50d0
c        gar_dens(4) = 1.50d0

c        gar_dens(1) = 0.5d0
c        gar_dens(2) = 0.5d0
c        gar_dens(3) = 0.50d0
c        gar_dens(4) = 0.50d0

         gar_dens(1) = 0.25d0
         gar_dens(2) = 0.50d0
         gar_dens(3) = 0.30d0
         gar_dens(4) = 0.30d0

c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(T)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 4
c         WRITE (6,9990) I, gnaf_dens(i), gcat_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1
        do i = 2, 119, 13
         level(i) = 2
        end do
        do i = 3, 120, 13
           level(i) = 3
           level(i+1) = 3
           level(i+2) = 3
        end do
        do i = 6, 123, 13
           level(i) = 4
           level(i+1) = 4
           level(i+2) = 4
           level(i+3) = 4
           level(i+4) = 4
           level(i+5) = 4
           level(i+6) = 4
           level(i+7) = 4
           level(i+8) = 4
        end do

        do i = 132, 137
         level(i) = 0
        end do

c connectivity of axon
        nnum(132) = 2
        nnum(133) = 3
        nnum(134) = 3
        nnum(136) = 3
        nnum(135) = 1
        nnum(137) = 1
         neigh(132,1) =  1
         neigh(132,2) = 133
         neigh(133,1) = 132
         neigh(133,2) = 134
         neigh(133,3) = 136
         neigh(134,1) = 133
         neigh(134,2) = 135
         neigh(134,3) = 136
         neigh(136,1) = 133
         neigh(136,2) = 134
         neigh(136,3) = 137
         neigh(135,1) = 134
         neigh(137,1) = 136

c connectivity of SD part
          nnum(1) = 11
          neigh(1,1) = 132
          neigh(1,2) =  2
          neigh(1,3) = 15
          neigh(1,4) = 28
          neigh(1,5) = 41
          neigh(1,6) = 54
          neigh(1,7) = 67
          neigh(1,8) = 80
          neigh(1,9) = 93
          neigh(1,10) = 106
          neigh(1,11) = 119

          do i = 2, 119, 13
           nnum(i) = 4
           neigh(i,1) = 1
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
           neigh(i,4) = i + 3
          end do

          do i = 3, 120, 13
           nnum(i) = 4
           neigh(i,1) = i - 1
           neigh(i,2) = i + 3
           neigh(i,3) = i + 4
           neigh(i,4) = i + 5
          end do

          do i = 4, 121, 13
           nnum(i) = 4
           neigh(i,1) = i - 2
           neigh(i,2) = i + 5
           neigh(i,3) = i + 6
           neigh(i,4) = i + 7
          end do

          do i = 5, 122, 13
           nnum(i) = 4
           neigh(i,1) = i - 3
           neigh(i,2) = i + 7
           neigh(i,3) = i + 8
           neigh(i,4) = i + 9
          end do

          do i = 6, 123, 13
           nnum(i) = 3
            neigh(i,1) = i - 3
            neigh(i,2) = i + 1
            neigh(i,3) = i + 2
          end do

          do i = 7, 124, 13
           nnum(i) = 3
           neigh(i,1) = i - 4
           neigh(i,2) = i - 1
           neigh(i,3) = i + 1
          end do

          do i = 8, 125, 13
           nnum(i) = 3
           neigh(i,1) = i - 5
           neigh(i,2) = i - 2
           neigh(i,3) = i - 1
          end do

          do i = 9, 126, 13
           nnum(i) = 3
           neigh(i,1) = i - 5
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
          end do

          do i = 10, 127, 13
           nnum(i) = 3
           neigh(i,1) = i - 6
           neigh(i,2) = i - 1
           neigh(i,3) = i + 1
          end do

          do i = 11, 128, 13
           nnum(i) = 3
           neigh(i,1) = i - 7
           neigh(i,2) = i - 2
           neigh(i,3) = i - 1
          end do

          do i = 12, 129, 13
           nnum(i) = 3
           neigh(i,1) = i - 7
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
          end do

          do i = 13, 130, 13
           nnum(i) = 3
           neigh(i,1) = i - 8
           neigh(i,2) = i - 1
           neigh(i,3) = i + 1
          end do

          do i = 14, 131, 13
           nnum(i) = 3
           neigh(i,1) = i - 9
           neigh(i,2) = i - 2
           neigh(i,3) = i - 1
          end do

         DO 332, I = 1, 137
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5),NEIGH(I,6),NEIGH(I,7),NEIGH(I,8),NEIGH(I,9),
c    X NEIGH(I,10),NEIGH(I,11)
3330     FORMAT(2X,12I5)
332      CONTINUE
          DO 858, I = 1, 137
           DO 858, J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
            DO 859, L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
              STOP
             ENDIF
858       CONTINUE

c length and radius of axonal compartments
          do i = 132, 137
            len(i) = 50.d0
          end do
          rad(132) = 0.80d0
          rad(133) = 0.7d0
          do i = 134, 137
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 42.d0
          rad(1) = 10.d0

          do i = 2, 119, 13
           len(i) = 20.d0
          end do

          do i = 3, 131
           if (level(i).ge.3) len(i) = 57.5d0
          end do

          do i = 2, 131
           if (level(i).eq.2) then
               rad(i) = 0.73d0
           else if (level(i).eq.3) then
               rad(i) = 0.8d0 * 0.73d0
           else
               rad(i) = 0.6d0 * 0.73d0
           endif
          end do


c       WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, 137
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 137
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
C NO CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
         DO 1019, I = 2, 131
           Z = Z + AREA(I)
1019     CONTINUE
c        WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

        DO 140, I = 1, numcomp
        DO 140, K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )

         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
140     CONTINUE
c gam computed in microS

        DO 299, I = 1, numcomp
299       BETCHI(I) = .05d0
        BETCHI( 1) =  .02d0

        DO 300, I = 1, numcomp
c300     D(I) = 2.D-4
300     D(I) = 5.D-4
        DO 301, I = 1, numcomp
         IF (LEVEL(I).EQ.1) D(I) = 1.D-3
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


       DO 160, I = 1, numcomp
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
C     NOTE CORRECTION

        do 200, i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
C     TO GO FROM MICROF TO NF.

      DO 909, I = 1, numcomp
       JACOB(I,I) = - GL(I)
      DO 909, J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

       DO 500, I = 1, numcomp
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F7.4)
500     CONTINUE
        END

