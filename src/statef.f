**==STATEF.FOR
      SUBROUTINE STATEF(FL, TL)
      IMPLICIT REAL*8(A-H, N-Z)
      PARAMETER (EPS=1.0D-30)
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /STAT1 / CSX(10), CS(121,191,10), CNIU(60,41,2), W(18400), 
     :                JCSX
      COMMON /STAT2 / PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA, CP, 
     :                CH, S, PR, PG, PF, PT, EN, WR(41)
      COMMON /STATFD/ RE, PE, QE, RET, PET, QET, REF, PEF, QEF, FDD(7)
      COMMON /ABUND / XA(10), NA(10), NEO, NIO, NZZ, AVM
      COMMON /ATDATA/ DH2, D1, D2, D3, CHI(26,9), OMG(27), AM(10), BN(10),
     &                IZ(10)
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CC, CD, CG, CR, CT, CEVB,
     &                CEN, CPL, CW(6)
      COMMON /IONISE/ NSTORE(8,5), NSTORE2(5)
      COMMON /OPDAT / cbase,obase,opT(141),opR(31),ZS
      COMMON /DEGEN/ PSI, VAL_PSI, V_X, V_Y
      COMMON/OTHR/ VAL_MAIN, VAL_X, VAL_Y 
      COMMON/OTHR2/ DPAA, DPBB, GAMMA, CV, TST
      common/massfracs/XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XHE3
      FXP(VX) = DEXP(MAX(-50.0D0,MIN(50.0D0,VX)))
      DIMENSION HA(26), VA(26), tkappa(2,2),Xcompos(5),COcompos(8)
      data Xcompos /0.0d0,0.03d0,0.1d0,0.35d0,0.7d0/
      data COcompos /0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.0d0/
      DATA JT, JX /2, 2/
* Evaluate Fermi-Dirac integrals according to Eggleton, Faulkner &
* Flannery (1973): RE, PE, SE and UE correspond to rho*, P*, S* and U*. 
* PSI is the usual degeneracy parameter
c     XA(10) are as follows: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XHE3
      XA(1)=XH
      XA(2)=XHE
      XA(3)=XC
      XA(4)=XN
      XA(5)=XO
      XA(6)=XNE
      XA(7)=XMG
      XA(8)=XSI
      XA(9)=XFE
      XA(10)=XHE3	
      F = EXP(FL)
c      T = EXP(TL) !I changed this to match with OPAL
      T = 1e6*TL	!TL = T6 
      UF = F/(1.0+F)
      WF = SQRT(1.0+F)
      PSI = FL + 2.0*(WF-LOG(1.0+WF))
      P2 = PSI    
      G = CT*T*WF
      CALL FDIRAC(F, G)
      PE = G*PE
      PET = PET + 1.0
      PEF = PEF + 0.5*UF
      QE = QE/(RE*WF)
      SE = QE + 2.0*WF - PSI
      SEF = QE*(QEF-REF-0.5*UF) - 1.0/WF
      SET = QE*(QET-RET)
      UE = SE + PSI - PE/(RE*CT*T)
* Some quantities that do not depend on the state of ionization:
* the NA are the element number densities (per baryon); AVM is the average
* mass per baryon (in amu); NEO and NIO are the numbers of electrons and
* ions assuming complete ionization.
      NEO = 0.0
      NIO = 0.0
      NZZ = 0.0
      AVM = 0.0
      DO I = 1, 10
         NA(I) = XA(I)/BN(I)
         AVM = AVM + AM(I)*NA(I)
         NIO = NIO + NA(I)
         NEO = NEO + IZ(I)*NA(I)
         NZZ = NZZ + IZ(I)*IZ(I)*NA(I)
      END DO
* PRESSI gives a crude model for pressure ionization and (if ICL=1) 
* a model for Coulomb interactions, returning corrections to the electron
* chemical potential, pressure, entropy and internal energy.
* TI is 1eV/kT, DE is 1 amu * number of electrons/cm3
      TI = CEVB/T
      DE = RE*CD
      CALL PRESSI(ICL, TI, DE, REF, RET, F, DC, DVT, DVF, DPA, DPAT,
     &     DPAF, DSA, DSAT, DSAF, DUA)
      DV = DC - PSI
      DVF = DVF - WF
* Contributions of the completely ionized species
      NE = 0.0
      NEF = 0.0
      NET = 0.0
      SI = 0.0
      SIF = 0.0
      SIT = 0.0
      UI = 0.0
      DO I = ION+1, 9
         NE = NE + IZ(I)*NA(I)
         VM = AM(I)*SQRT(AM(I))
         SI = SI - NA(I)*LOG(NA(I)/VM + EPS)
      END DO
* Calculate ionization of the first ION elements.
      DO I = ION, 1, -1
         SHA = 1.0
         SJHA = 0.0
         SCHA = 0.0
* compute potentials VA and number ratios HA of ionization state J
* relative to the ground state
         DO J = 1, IZ(I)
            VA(J) = -CHI(J,I)*TI + J*DV
            IF(J.EQ.1) THEN
               HA(J) = FXP(VA(J))*OMG(IZ(I))/OMG(IZ(I)+1)
            ELSE
               HA(J) = HA(J-1)*FXP(VA(J) - VA(J-1))*OMG(IZ(I)+1-J)
     &              /OMG(IZ(I)+2-J)
            END IF
            SHA = SHA + HA(J)
            SJHA = SJHA + J*HA(J)
            SCHA = SCHA + CHI(J,I)*TI*HA(J)
            NSTORE(J,I) = HA(J)
         END DO
         NSTORE2(I) = SHA
         VM = AM(I)*SQRT(AM(I))
         SI = SI + NA(I)*LOG(VM)
         IF(I.GT.1) THEN
* contributions to electron number density NE, entropy SI and
* internal energy UI for Helium and heavier
            VX = NA(I)/SHA
            SI = SI - VX*LOG(VX/OMG(IZ(I)+1) + EPS)
            DO J = 1, IZ(I)
               NX = HA(J)*VX
               NXF = NX*DVF*(J - SJHA/SHA)
               NXT = NXF*DVT/DVF + NX*(CHI(J,I)*TI - SCHA/SHA)
               NE = NE + J*NX
               NEF = NEF + J*NXF
               NET = NET + J*NXT
               SIF = SIF - VA(J)*NXF
               SIT = SIT - VA(J)*NXT
               SI = SI - NX*LOG(NX/OMG(IZ(I)+1-J) + EPS)
               UI = UI + CHI(J,I)*NX
            END DO
         END IF
      END DO
* Ionization and molecular dissciation of Hydrogen.
* partition function for H2 from Vardya (1960), Webbink (1975)
      DH2TI = DH2*TI	!D_H2/T
      D1TI = D1*TI	!D1/T
      D2TI = (D2*TI)**2 !(D2/T)^2
      D3TI = (D3*TI)**3 !(D3/T)^3
      ZET = 1.0 - (1.0 + DH2TI)*EXP(-DH2TI) !eq 3 PTE 
      DZH2T = -DH2TI**2*EXP(-DH2TI)*TI  !NOT /ZET, BUT /T
C      DZH2TT = (DH2TI - 2.0 - DZH2T)*DZH2T
      DZH2TT = EXP(-DH2TI)*(3*DH2TI**2 - DH2TI**3)*TI
      ZH2 = 6608.8*ZET*DH2TI**(-2.5)*EXP(-D1TI - D2TI - D3TI) !eq 2 PTE 
      ZH2T = 2.5 + D1TI + 2.0*D2TI + 3.0*D3TI + DZH2T
      ZH2TT = -D1TI - 4.0*D2TI - 9.0*D3TI + DZH2TT
      ZH2S = ZH2*SQRT(8D0)/VM
      H2A = CEN*(ZH2S/4.0)*(DE)/(T*SQRT(T))*EXP(DH2TI)
      H2BT = DH2TI + 1.5 - ZH2T
      H2AT = RET - H2BT
* solve for densities of H+, H, and H2
      AA = 2*H2A + HA(1)*(1.0 + HA(1))
      AB = NE + HA(1)*(NE - NA(1))
      AC = NA(1)*NE
      HG = 2*AC/(SQRT(AB*AB + 4*AA*AC) + AB)
      IF (HG.NE.HG) HG = 0.0	
      HI = HA(1)*HG
      NE = NE + HI
      EN = 1.0/NE
      H2 = H2A*HG*HG*EN	!eq 13 PTE
      IF (H2.NE.H2) H2=0.0
      NI = NIO - H2
* derivatives w.r.t. F and T
      AA = NE + 4*HG*H2A !
      IF (AA.NE.AA) AA = NE
      AB = HA(1)*(NE - 2*H2)
      AD = 1/(AA + AB)
C      IF (AD.NE.AD) AD=0.0
      AC = 2*H2*AD
C      IF (AC.NE.AC) AC=0.0
      AF = (NEF - NE*REF)*AC
      AT = (NET - NE*H2AT)*AC !
      AE = HG*AD
C      IF (AE.NE.AE) AE=0.0	
      BF = DVF*AE
      BT = (CHI(1,1)*TI + DVT)*AE
      HGF = AF - AB*BF
      HGT = AT - AB*BT
      HIF = HA(1)*(AF + AA*BF)
      HIT = HA(1)*(AT + AA*BT) !
      NEF = NEF + HIF
      NET = NET + HIT !
      H2F = H2*REF + EN*(2*H2A*HG*HGF - H2*NEF)
      IF (H2F.NE.H2F) H2F = 0.0
      H2T = H2*H2AT + EN*(2*H2A*HG*HGT - H2*NET)
      IF (H2T.NE.H2T) H2T = 0.0

* hydrogen contribution to entropy, internal energy
      SIF = SIF - VA(1)*HIF - H2BT*H2F
      SIT = SIT - VA(1)*HIT - H2BT*H2T + H2*(ZH2T + ZH2TT)
      SI = SI - HI*LOG(HI/OMG(1) + EPS) - HG*LOG(HG/OMG(2) + EPS)
     &        - H2*(LOG(H2/ZH2S + EPS) + ZH2T)
      UI = UI + CHI(1,1)*HI + 0.5*DH2*(HI + HG)
* DB is 1 amu * number of baryons/cm3; RHO is the mass density in g/cm3
      DB = EN*DE
      DL = LOG(DB)
      RHO = DB*AVM
      RL = LOG(RHO)
      RT = RET - NET/NEO !
      RF = REF - NEF/NEO ! 
* second call to PRESSI compensates for spurious pressure and entropy terms
      DE = DB*NEO
      CALL PRESSI(0, TI, DE, RF, RT, F, DC, DVT, DVF, DPB, DPBT, DPBF, 
     :            DSB, DSBT, DSBF, DUB)
* pressure terms
      PE = CB*PE
      TCR = T*CR
      P0 = TCR*DB
      PI = NI*P0
      T4 = T*T*T*T
      PR = CA*T4/3.0
      B = 4.0*PR/P0
      PG = PE + PI + TCR*(DPA-DPB)
      DPAA= DPA*TCR	
      DPBB= DPB*TCR 
      P = PG + PR
      PF = (PE*PEF      + PI*RF - H2F*P0 + TCR*(DPAF-DPBF))/P
      PT = (PE*PET + PI + PI*RT - H2T*P0 + TCR*(DPAT-DPBT) + PR*4.0)/P
      PL = LOG(P)
* entropy, in erg/g/K
      DSF = NEF*DSA + NE*DSAF - NEO*DSBF - RF*B
      DST = NET*DSA + NE*DSAT - NEO*DSBT - (RT - 3.0)*B
      SF = CR*(-NI*RF       + NEF*SE + NE*SEF + SIF + DSF)/AVM
      ST = CR*( NI*(1.5-RT) + NET*SE + NE*SET + SIT + DST)/AVM
      S = CR*(SE*NE + DSA*NE - DSB*NEO + B + (1.5*TL - DL + 2.5
     :        - LOG(CEN))*NI + SI)/AVM
* internal energy, in erg/g
      U = TCR*(UE*NE + DUA*NE - DUB*NEO + 1.5*NI + ZH2T*H2 + 0.75*B
     :         + TI*UI)/AVM
* other thermodynamic quantities
      Q = PT*SF - PF*ST
      CP = -Q/PF
      GRADA = SF/Q
      GAMMA = Q/(RT*SF-RF*ST)			!Gamma = CP/CV
      CV = -(RT*SF-RF*ST)/PF
c      CV = H2T
      ZT = SQRT(ABS((NE*REF/WF + NZZ)/NI))
* TST ought to be zero, if all the above programming is correct
      TST = SF/CR - P*(RT*PF-RF*PT)/(TCR*RHO)
      RETURN
      END
