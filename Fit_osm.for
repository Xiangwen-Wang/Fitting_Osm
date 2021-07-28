C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  Modules and routines
C  --------------------
C  
C  Fit_related   - this should not need altering.
C  
C
C  Model_related - solute information (**different for each solute**)
C
C  
C  main program  - does the fit and writes the results
C
C  
C  LSFUN1 - this subroutine is called by NAG routine E04FYF. It calls CalcOsm
C           and then assigns the 'residuals' (observed - fitted) for all the 
C           data.
C  
C
C  CalcOsm - this is the routine that solves the model (the bisection routine
C            is inside)
C
C  
C  Func - the function called by the bisection code to determine h.
C  
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
C
      MODULE Fit_related
C
      IMPLICIT NONE
C
C     ++++++
      PUBLIC
C     ++++++
C
!     ==================================================================
!   ..this module contains the data to be fitted, and most thermodynamic
!     and model constants.
!     ==================================================================
C
!   ..max numbers of data points and parameters to be fitted (actual values 
!     should always be less than or equal to the values of these constants).
!     ----------------------------------------------------------------------
      INTEGER, PARAMETER :: nValuesMax = 100, nParasMax = 10, Nmax = 25
C
C
!   ..data points to be fitted should be placed in these arrays.
!     -------------------------------------------------------------------   
      DOUBLE PRECISION :: mSolute(nValuesMax), Csolute(nValuesMax), 
     >                    Y(nValuesMax), T(nValuesMax),
     >                    VSol(nValuesMax), VH2O(nValuesMax), 
     >                    Fitted(nValuesMax), Wt(nValuesMax)
      INTEGER :: id(nValuesMax)
C
C
!   ..stored values of h from each result
!     -----------------------------------
      DOUBLE PRECISION :: h_all(nValuesmax)
C
C
      END MODULE Fit_related
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      MODULE Model_related
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C
!   ..the charges on the cation (zC) and anion (zA), 
!     and ion size (aIon in nm)
!     ----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: zC = 2.D0, zA = 1.D0, aIon = 0.4D0
C
!   ..the stoichiometric number of the solute (v) and
!     individual ions in the molecule (vC and vA)
!     ----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: vC = 1.D0, vA = 2.D0, v = vC+vA
C
C
!   ..Avogadro's number, and PI
!     -------------------------
      DOUBLE PRECISION, PARAMETER :: NAvo = 6.022D23, PI = 3.14159D0
C
C
!   ..the number of hydration stages (read from the data file)
!     --------------------------------------------------------
      INTEGER :: N
C
C
      END MODULE
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      PROGRAM Fit
C
C     +++
      USE Fit_related
C     +++
      USE Model_related
C     +++
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
!   ..this is used to set the size of an array for a NAG routine
!     ----------------------------------------------------------     
      INTEGER, PARAMETER :: LW = 10 + 7*nParasMax + nParasMax*nParasMax
     >                          + 2*nValuesMax*nParasMax + 3*nValuesMax 
     >                          + nParasMax*(nParasMax-1)/2
C 
!   ..user arrays (free to use for anything)
!     --------------------------------------
      INTEGER :: iUser(nValuesMax)
      DOUBLE PRECISION :: User(nValuesMax)
C
!   ..arrays of parameters, residuals, and three NAG arrays
!     --------------------------------------------
      DOUBLE PRECISION ::  W(LW), CJ(nParasMax), Work(nParasMax), 
     >                     Resid(nValuesMax), Params(nParasMax)
C
      DOUBLE PRECISION :: K_big, K_small, Ki, Cratio(Nmax), 
     >                    CiC0(nValuesMax,Nmax)
C
      EXTERNAL LSFUN1
C
C
C
!   ..open two files, one to read the data from, 
!     and the second to write the results to.
!     ------------------------------------------
      OPEN(1, FILE='Fit_osm.dat', STATUS='OLD')       ! input data file
      OPEN(2, FILE='Fit_osm.res', STATUS='UNKNOWN')   ! results file
C
C
!   ..1st line: read the total number of data points (inc. all
!     zero weighted points, and any dummy values for extrapolation),
!     and the number of fitted parameters.
!     --------------------------------------------------------------
      READ(1,*) nValues, nParas
C
C
C   ..2nd line: the number of hydration stages (N),
C     followed by first estimates of the parameters
C     Param(1) = K, Param(2) = k
!     ---------------------------------------------
      READ(1,*) N, (Params(I),I=1,nParas)
C
C
C   ..skip the line with the headers for the columns 
!     ----------------------------------------------
      READ(1,*)
C
C
C   ..read all the data, one line at a time. You could also 
!     do any necessary data processing in this loop. 'Y'
!     are the experimental osmotic coefficients that we 
!     are trying to fit.
!     -----------------------------------------------------
      DO I=1,nValues
C
        READ(1,*) mSolute(I), Csolute(I), VH2O(I), VSol(I), 
     >            Y(I), T(I), Wt(I), id(I)
C
!     ..change units to litres (decimetres^3) per mol
!       ---------------------------------------------
        VSol(I) = 0.001D0 * VSol(I)
        VH2O(I) = 0.001D0 * VH2O(I)
C
      ENDDO
C
C
C
!     ================================================
!    | Now fit the parameters                         |
!     ================================================
      iFail = 1
      Call E04FYF(nValues, nParas, LSFUN1, Params, FSumSq, W, LW,
     >            iUser, User, iFail)
C
C   ..Output the sum of squares for the completed calculation,  
C     and check the value of iFail, which usually will be 0 or 
C     5 for a successful calculation.
!     --------------------------------------------------------
      IF((iFail.NE.1) .AND. (iFail.NE.2)) THEN
        WRITE(2,100) iFail
        WRITE(2,200) FSUMSQ
      ELSE
        WRITE(2,300) iFail
C       ****
C     STOP
C       ****
      ENDIF
C
C
C    ===============================================================
C   | Compute estimates of the variances of the sample regression   |
C   | coefficients at the final point.                              |
C   |                                                               |
C   | E04YCF(JOB, M, N, FSumSq, S, V, LV, CJ, Work, iFail           |
C   |                                                               |
C   | M = number of observations (here nValues)                     |
C   |                                                               |
C   | N = number of fitted parameters (here nParas)                 |
C   |                                                               |
C   | S(N) is the array of singular values of the Jacobian returned |
C   |      by E04FYF (in W, starting at W(NS), where NS is defined  |
C   |      further below).                                          |
C   |                                                               |
C   | V(LV,N) is the N x N right-hand orthogonal matrix of          |        
C   |       J as returned by E04FDF. When V is passed in the        |
C   |       workspace array W (argument W(NV)) following E04FYF     |
C   |       LV must be the value N.                                 |
C   |                                                               |
C   | CJ(N) when Job = 0, CJ returns the N diagonal elements of C.  |
C   |       That is to say CJ(1->N) contains the variances of       |
C   |       fitted parameters 1 -> N.                               |
C   | WORK  is a work array (not used for anything special).        |
C   |                                                               |
C   | ------------------------------------------------------        |
C   | So, following E04FYF, the routine is called like this:        |
C   |                                                               |
C   | NS = 6*N + 2*M + M*N + 1 + MAX(1,(N*(N-1))/2)                 |
C   | NV = NS + N                                                   |
C   |                                                               |                                                                 
C   | iFail = 1                                                     |
C   | Call E04YCF(0, M, N, FSumSq, W(NS), W(NV), N, CJ, Work, iFail)|
C   |                                                               |
C    ===============================================================
C     
      NS = 6*nParas + 2*nValues + nValues*nParas + 1 
     >     + MAX(1,(nParas*(nParas-1))/2)
      NV = NS + nParas
C     
C     
C   ..Compute the uncertainties
C     -------------------------
      iFail = 1
      Call E04YCF(0, nValues, nParas, FSumSQ, W(NS), W(NV), nParas, 
     >            CJ, Work, iFail)
C
C
C   ..Output the parameter values, their standard errors, and 
C     the ratio (a value of 4 or higher is good).
!     -----------------------------------------------------------     
      IF((iFail.NE.1) .AND. (iFail.NE.2)) THEN
        WRITE(2,400) iFail
        WRITE(2,500) (I, Params(I), SQRT(CJ(I)), Params(I)/SQRT(CJ(I)),
     >                I=1,nParas)
      ELSE
        WRITE(2,600) iFail
C       ****
        STOP
C       ****
      ENDIF
C
C
C   ..Calculate the residuals, and output the results
C     -----------------------------------------------
      Resid(:nValues) = Fitted(:nValues) - Y(:nValues)
C
      WRITE(2,700) (I, mSolute(I), Y(I), Fitted(I), Resid(I), Wt(I),
     >              T(I), id(I), I=1,nValues)
C
C
C
C     ------------------------------------------------
C   ..Calculate Cratio (C(i)/C(i-1)), and C(i)/C(0), 
C     using the stored values of h (h_all) for the 
C     calculation.
C     ------------------------------------------------
      K_big   = Params(1)
      K_small = Params(2)
C
      WRITE(2,900) (I,I=1,N)  ! header
C
      DO iVal=1,nValues
C
        YConst = EXP(Csolute(iVal)*(VSol(iVal) 
     >               + h_all(iVal)*VH2O(iVal) 
     >               - v*VH2O(iVal)))
C
        aW = EXP(-0.0180152D0*v*mSolute(iVal)*Y(iVal))
C
!     ..now calculate the concentration (molarity) ratios
!       for each stepwise equilibrium (eq 19), where I=1,N
!       --------------------------------------------------
        DO I = 1,N
          Ki = K_big * k_small**(I-1)
          Cratio(I) = (Ki / YConst) * aW
        ENDDO
C
        DO I = 1,N
          IF(I .EQ. 1) THEN
            CiC0(iVal,I) = Cratio(I)
          ELSE
            CiC0(iVal,I) = Cratio(I) * CiC0(iVal,I-1)
          ENDIF
        ENDDO
C
        WRITE(2,910) mSolute(iVal), h_all(iVal), (CiC0(iVal,I),I=1,N)
C
      ENDDO
C
C
!     -----------------------------------------
!   ..now write the stepwise hydration profiles
!     that are shown in Figure 2.
!     -----------------------------------------
      WRITE(2,920) (I,I=0,N)  ! header
C
      DO iVal=1,nValues
C
!     ..calculate C0 for this solution
!       ------------------------------
        C0 = CSolute(iVal) / (1.D0 + SUM(CiC0(iVal,:N)))
C
!     ..values of Ci/C
!       --------------
        WRITE(2,930) mSolute(iVal), C0/CSolute(iVal), 
     >               (C0*CiC0(iVal,I)/CSolute(iVal),I=1,N)
      ENDDO
C
C
C     ****
      STOP
C     ****
C
100   FORMAT(//' Exit OK from E04FYF: iFail = ',I3)
200   FORMAT(1X,'On exit, the weighted sum of squares = ',F12.6)
300   FORMAT(/' !Error on exit, NAG failure: iFail= ',I5,/
     >        ' *** TERMINATED ***')
400   FORMAT(/' Exit OK from E04YCF: IFAIL = ',I3)
500   FORMAT(/' Parameters and their standard errors:',//
     >'        Parameter       Std. Error      Ratio'/
     >(1X,T2,I2,T6,E14.6,T23,E12.4,T39,F8.4))
600   FORMAT(//' !Error exit from E04YCF:  iFail = ',I3,/
     >         ' - see routine document')
700   FORMAT(//' ----------------- Results (1) -------------------',//
     >'   i        m            osm           fit         resid  ',
     >'     wt      T      id'/
     >(1X,I3,2X,3(E12.4,2X),F10.6,2X,F7.3,2X,F6.2,2X,I3))
900   FORMAT(//1X,'      m         h         ',<N>('C',I0,'C0',8X))
910   FORMAT(1X,F9.5,30(1X,ES11.4))
920   FORMAT(//1X,'      m        ',<N+1>('C',I0,'C',9X))
930   FORMAT(1X,F9.5,30(1X,ES11.4))
C
      END Program
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      Subroutine LSFUN1(nValues, nParas, Params, Resid, iUser, User)
C
C     +++
      USE Fit_related
C     +++
      USE Model_related
C     +++
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
!  -- arguments --
      DOUBLE PRECISION :: Resid(nValues), Params(nParas), User(*)
      INTEGER :: iUser(*)
C
!  -- local --
      DOUBLE PRECISION :: K_big, K_small
C
C
!   ..the two hydration parameters we are fitting
!     -------------------------------------------  
      K_big   = Params(1)
      K_small = Params(2)
C
C
      DO iVal = 1, nValues
C
        osmExpt = Y(iVal)
C
        Call CalcOsm(mSolute(iVal), Csolute(iVal), VSol(iVal), 
     >               VH2O(iVal), K_big, k_small, v, osmExpt, 
     >               osmCalc, h, N)

C
C     ..fitted quantity, and residual required by E04FYF
C       ------------------------------------------------
        Fitted(iVal) = osmCalc
        Resid(iVal)  = SQRT(Wt(iVal)) * (Fitted(iVal)-Y(iVal))
C
!     ..save the value of h for each solution
!       -------------------------------------
        h_all(iVal)  = h
C
      ENDDO
C
C
C     ****** 
      RETURN
C     ******
C
      END Subroutine
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CalcOsm(mSolute, Csolute, VSol, VH2O, K_big, k_small, 
     >                   v, osmExpt, osmCalc, h, N)
C
C     +++
      USE Model_related, ONLY : vC, vA, zC, zA, aIon, PI, NAvo
C     +++
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     --------------------------------------------------------------   
C     This routine receives the properties of a single solution,  
C     and calculates the osmotic coefficient according to the model. 
C
C     Inputs:
C     ------
C     mSolute = solute molality
C     Csolute = solute molarity 
C     VSol    = partial molar volume of the solute 
C     VH2O,   = partial molar volume of water
C     K_big   = K 
C     k_small = k 
C     N       = number of stepwise hydrations
C     osmExpt = experimental osmotic coefficient
C
C     Output:
C     ------
C     osmCalc = calculated osmotic coefficient
C     h       = calculated average number of waters of hydration
C     --------------------------------------------------------------
C
!   ..parameters for the bisection solution of the equations
!     ------------------------------------------------------
      INTEGER, PARAMETER :: Jmax = 250
      DOUBLE PRECISION, PARAMETER :: x_acc = 1.D-10, Hfrac = 0.999D0
C
!  -- arguments --
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN)  :: VSol, VH2O, K_big, k_small,  
     >                                 mSolute, Csolute, osmExpt, v
      DOUBLE PRECISION, INTENT(OUT) :: osmCalc, h
C
!  -- local --
      DOUBLE PRECISION :: ln_aW, Ionstr, Kappa
C
C
C
!     ===========================
!   ..Solve for h using bisection
!     ===========================
C
      aW = EXP(-0.0180152D0*v*mSolute*osmExpt)
C
C
!   ..assign the two guesses
!     ----------------------
      x1 = 0.0001D0  ! low value
C
!   ..the first part of the expression below is the 
!     number of water molecules per "molecule" of 
!     solute. We assume the upper bound of possible
!     h is equal to this number * some arbitrary (but
!     large) frction Hfrac
!     -----------------------------------------------
      x2 = (1.D0/0.0180152D0)/mSolute * Hfrac  ! high value
C
C
!   ..initial estimates
!     -----------------
      fmid = func(x2, Csolute, VSol, VH2O, aW, K_big, k_small, v, N)
      f    = func(x1, Csolute, VSol, VH2O, aW, K_big, k_small, v, N)
C
      IF(f*fmid .GE. 0.D0) THEN
        WRITE(*,*) 'root must be bracketed in rtbis', mSolute
C       ****
        STOP
C       ****
      ENDIF
C
      IF(f .LT. 0.D0) then
        rtbis = x1
        dx = x2 - x1
      ELSE
        rtbis = x2
        dx = x1 - x2
      ENDIF
C
      DO j=1,Jmax
        dx = dx * 0.5D0
        xmid = rtbis + dx
C
        fmid = func(xmid, Csolute, VSol, VH2O, aW, K_big, k_small, v, N)
        IF(fmid .LE. 0.D0) rtbis = xmid
C
        IF(ABS(dx).LT.x_acc .or. fmid .EQ. 0.D0) EXIT
C
        IF(j .EQ. Jmax) THEN
          WRITE(*,*) 'too many bisections. xmid = ', xmid
C         ****
          STOP
C         ****
        ENDIF
      ENDDO
C
!   ..rtbis is the value of the average hydration number h
!     -----------------------------------------------------
      h = rtbis
C
C
C
!     ===================================================
!     Now the water activity (aW) and osmotic coefficient
!     ===================================================
C
!   ..Debye Huckel term (expression for kappa below is for 
!     25 oC and comes from Wikipedia. The full expression
!     is more complex.)
!     ---------------------------------------------------
      Ionstr = 0.5D0*(vC*Csolute*zC**2 + vA*Csolute*zA**2)
      Kappa = SQRT(Ionstr) / 0.304D0
C
!   ..note: units of kappa are 1/nm and aIon are
!     nm so that the product is dimensionless.
!     ------------------------------------------
      dum = Kappa * aIon
      Sfunc = 3.D0/dum**3 
     >       * (1 + dum - 1.D0/(1.D0 + dum) - 2*LOG(1.D0 + dum))
C
!   ..note the 1.D8 factor converts the unit of Kappa
!     from 1/nm to 1/dm (decimetres) which is the same 
!     unit as VH2O (dm^3 per mol)
!     ------------------------------------------------
      Q = (1.D8*Kappa)**3/(24*PI*NAvo*CSolute**1.5D0)      

      DH_term = Q * VH2O * cSolute**1.5D0 * Sfunc
C
C
!   ..finally, the water activity and osmotic coefficient
!     ---------------------------------------------------
      V_h = VSol + h*VH2O
C
      ln_aW = LOG(1.D0 - CSolute*V_h) + Csolute*(V_h - v*VH2O) + DH_term
C
      OsmCalc = -ln_aw /(0.0180152D0 * v * mSolute)
C
C     ******
      RETURN
C     ******
C
      END Subroutine
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION func(h_guess, Csolute, VSol, VH2O, aW, K_big, k_small, 
     >              v, N)
C
C     ----------------------------------------------------------
C     This function calculates the difference between a value of 
C     h calculated from eq (23) of the Stokes paper, and an
C     initial guess (h_guess).
C     ----------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
!  -- arguments --
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN) :: h_guess, VSol, VH2O, K_big, 
     >                                k_small, aW, Csolute, v
C
!  -- local --
      DOUBLE PRECISION :: Cratio(N), CiC0(N), Ki, numerator 
C
C
C
!   ..value of Y from eq(18)
!     ----------------------
      Y = EXP(Csolute*(VSol + h_guess*VH2O - v*VH2O))
C
C
!   ..now calculate the concentration (molarity) ratios
!     for each stepwise equilibrium (eq 19), where I=1,N
!     --------------------------------------------------
      DO I = 1,N
C
        Ki = K_big * k_small**(I-1)
C
        Cratio(I) = (Ki / Y) * aW
C
      ENDDO
C
C
!   ..calculate the value of all ci/c0 in eq(20)
!     where i=1,N
!     ----------------------------------------------
      DO I = 1,N
C
        IF(I .EQ. 1) THEN
          CiC0(I) = Cratio(I)
        ELSE
          CiC0(I) = Cratio(I) * CiC0(I-1)
        ENDIF
C       	
      ENDDO
C
C
!   ..calculate the numerator and denominator in eq(23)
!     -------------------------------------------------
      numerator   = 0.D0
      denominator = 0.D0
C
      DO I = 0,N
C
        IF(I .EQ. 0) THEN
          numerator   = numerator   + 0.D0
          denominator = denominator + 1.D0
        ELSE
          numerator   = numerator   + I * CiC0(I)
          denominator = denominator + CiC0(I)
        ENDIF
C
      ENDDO
C
C
!   ..the value of h from eq (23)
!     ---------------------------
      h_calc = numerator / denominator
C
C
!   ..the function value
!     ------------------
      func = h_guess - h_calc
C
      END Function func
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
