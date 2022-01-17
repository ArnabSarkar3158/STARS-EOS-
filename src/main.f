	PROGRAM FD
	IMPLICIT REAL*8 (A-H,L,M,O-Z)
c	COMMON /STATFD/ D(3,3), DTT, DFT, DFF, DTTT, DFTT, DFFT, DFFF
       COMMON /STATFD/ RE, PE, QE, RET, PET, QET, REF, PEF, QEF, FDD(7)!FDD are the second order/third order derivatives of density, the rest (answer and their 1st derivative) are D(3,3)
       COMMON /STAT2 / PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA, CP, CH, S, PR, PG, PF, PT, EN, WR(41)
      COMMON /DEGEN/ PSI
      COMMON /OTHERS/ DMU, TGA, TGB, TX, TY	
      COMMON/OTHERS2/ H2O, DK, DK2
      COMMON /OT3/ EP, EU, SE
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON/OTHR/ VAL_MAIN
      COMMON/OTHR2/ DPAA, DPBB, GAMMA, CV, TST
      common/massfracs/XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XHE3

      CALL CONSTS

      ICL=1!TURNS ON/OFF (1/0) COULOMB CORRECTIONS
      ION=5!CALCULATES IONIZATION FOR ALL 5 ELEMENTS (3, 4, 5 => till C, N, O). DEFAULT IS 2 (=> H AND HE).

******Set X, Y, Z here:
      XH=0.0
      XHE=9.8478066E-01 
      XC=3.4796335179533172E-003
      XN= 0.0000000000000000
      XO=9.6975041099999996E-003
      XNE=0.0000000000000000
      XMG=3.8993723720466508E-003
      XSI=1.6000000000000001E-003
      XFE=1.4399999999999999E-003
      XHE3= 7.8040609560643153E-018	    
******
99002   FORMAT (E14.6,2E14.6, 3E14.6, E14.6, I4)
	OPEN(1, FILE="STARS_EOS")
	write(1,*) 'X = ', XH, 'Y = ', XHE 
	WRITE(1,*) 'Column 1: T [K] '
	write(1,*) 'Column 2: Density [g/cc]'
	write(1,*) 'Column 3: Gas Pressure [dyne/cm^2] with ICL:', ICL
	write(1,*) 'Column 4: Radiation Pressure [dyne/cm^2]'
	write(1,*) 'Column 5: Entropy [erg/g/K]'
	write(1,*) 'Column 6: Int. Energy [erg/g]'
	WRITE(1,*) 'Column 7: Cp/Cv'
	WRITE(1,*) '-1 in the last column indicates that the gas pressure numbers are unreliable'
	write(1,*) ' '
	DO I=-80,88,1 !(loge_f goes from -10.61 to 7.82 according to Eggleton's paper)
		DO J=20,-100,-1		!temperature goes from 10^8 to 10^-4 K (in 10^6 units). 
			FL = I/2.0	!-40 to 44 taken to iterate thorugh more density
			TL = 10**(J/10.0)
			CALL STATEF(FL,TL)
			if      (T.ge.1D3.and.T.le.1D8.and.RHO.ge.1D-12.and.RHO.LE.1e8) then
				if(PG.LE.0d0)WRITE(1,99002) T, RHO, PG, PR, S, U, GAMMA,-1!if coulomb interactions have made gas pressure negative
				if(PG.GE.0d0)WRITE(1,99002) T, RHO, PG, PR, S, U, GAMMA,+1
			end if
		END DO
	END DO
	STOP
	END

