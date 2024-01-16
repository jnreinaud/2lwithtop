      SUBROUTINE C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      
CVD$R NOVECTOR
C     .. Scalar Arguments ..
      INTEGER           IERROR, M, N, NQ
      CHARACTER*1       INIT
C     .. Array Arguments ..
      real*8  TRIG(2*N)
      INTEGER           Q(30)
C     .. Local Scalars ..
      INTEGER           NCHECK
C     .. External Subroutines ..
      EXTERNAL          C06FPY, C06FPZ
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     .. Save statement ..
      SAVE              NCHECK
C     .. Data statements ..
      DATA              NCHECK/-1/
C     .. Executable Statements ..
      IERROR = 0
C
      IF (M.LT.1) THEN
         IERROR = 1
         RETURN
      ELSE IF (N.LT.1) THEN
         IERROR = 2
         RETURN
      END IF
      IF (INIT.NE.'I' .AND. INIT.NE.'i' .AND. INIT.NE.'S' .AND. INIT.NE.
     *    's' .AND. INIT.NE.'R' .AND. INIT.NE.'r') THEN
         IERROR = 3
         RETURN
      END IF
      IF (INIT.EQ.'S' .OR. INIT.EQ.'s') THEN
         IF (NCHECK.EQ.-1) THEN
            IERROR = 4
            RETURN
         ELSE IF (NINT(TRIG(N)).NE.N .OR. NINT(TRIG(2*N)).NE.N) THEN
            IERROR = 5
            RETURN
         END IF
      END IF
      IF (INIT.EQ.'R' .OR. INIT.EQ.'r') THEN
         IF (NINT(TRIG(N)).NE.N .OR. NINT(TRIG(2*N)).NE.N) THEN
            IERROR = 5
            RETURN
         END IF
      END IF
      CALL C06FPZ(N,NQ,Q)
      IF (INIT.EQ.'I' .OR. INIT.EQ.'i') THEN
         CALL C06FPY(N,NQ,Q,TRIG(1),TRIG(N+1))
      END IF
      NCHECK = N
      RETURN
      END
	
      SUBROUTINE C06FPR(A,B,P,Q,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Real to Hermitian fast Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:Q-1,0:R-1), B(0:P-1,0:R-1,0:Q-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1)
C     .. Local Scalars ..
      real*8  AI, AR, TEMP, TEMP1, TEMP2, TEMP3, TEMPI, TEMPR,
     *                  TWOPI
      INTEGER           I, INDX, J, K, KP, L, Q2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      Q2 = (Q-1)/2
      TWOPI = 2.0D0 * acos(-1.d0)
C
C     Code for K=0 --
C
      DO 40 J = 1, Q2
         DO 20 I = 0, P - 1
            TEMP = A(I,J,0)
            A(I,J,0) = TEMP + A(I,Q-J,0)
            A(I,Q-J,0) = TEMP - A(I,Q-J,0)
   20    CONTINUE
   40 CONTINUE
      DO 120 L = 1, Q2
         DO 60 I = 0, P - 1
            B(I,0,L) = A(I,0,0)
            B(I,0,Q-L) = 0.0D0
   60    CONTINUE
         DO 100 J = 1, Q2
            INDX = MOD(J*L,Q)
            DO 80 I = 0, P - 1
               B(I,0,L) = B(I,0,L) + A(I,J,0)*COSINE(0,INDX)
               B(I,0,Q-L) = B(I,0,Q-L) + A(I,Q-J,0)*SINE(0,INDX)
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE
      DO 140 I = 0, P - 1
         B(I,0,0) = A(I,0,0)
  140 CONTINUE
      DO 180 J = 1, Q2
         DO 160 I = 0, P - 1
            B(I,0,0) = B(I,0,0) + A(I,J,0)
  160    CONTINUE
  180 CONTINUE
C
C     Code for general K --
C
      DO 460 K = 1, (R-1)/2
         KP = R - K
         DO 220 J = 1, Q - 1
            DO 200 I = 0, P - 1
               AR = A(I,J,K)
               AI = A(I,J,KP)
               A(I,J,K) = COSINE(K,J)*AR - SINE(K,J)*AI
               A(I,J,KP) = COSINE(K,J)*AI + SINE(K,J)*AR
  200       CONTINUE
  220    CONTINUE
         DO 260 J = 1, Q2
            DO 240 I = 0, P - 1
               TEMPR = A(I,J,K)
               TEMPI = A(I,J,KP)
               A(I,J,K) = TEMPR + A(I,Q-J,K)
               A(I,J,KP) = TEMPI + A(I,Q-J,KP)
               A(I,Q-J,K) = TEMPR - A(I,Q-J,K)
               A(I,Q-J,KP) = TEMPI - A(I,Q-J,KP)
  240       CONTINUE
  260    CONTINUE
         DO 340 L = 1, Q2
            DO 280 I = 0, P - 1
               B(I,K,L) = A(I,0,K)
               B(I,KP,Q-L-1) = A(I,0,KP)
               B(I,KP,L-1) = 0.0D0
               B(I,K,Q-L) = 0.0D0
  280       CONTINUE
            DO 320 J = 1, Q2
               INDX = MOD(J*L,Q)
               DO 300 I = 0, P - 1
                  B(I,K,L) = B(I,K,L) + A(I,J,K)*COSINE(0,INDX)
                  B(I,KP,Q-L-1) = B(I,KP,Q-L-1) + A(I,J,KP)*COSINE(0,
     *                            INDX)
                  B(I,KP,L-1) = B(I,KP,L-1) - A(I,Q-J,K)*SINE(0,INDX)
                  B(I,K,Q-L) = B(I,K,Q-L) + A(I,Q-J,KP)*SINE(0,INDX)
  300          CONTINUE
  320       CONTINUE
  340    CONTINUE
         DO 360 I = 0, P - 1
            B(I,K,0) = A(I,0,K)
            B(I,KP,Q-1) = A(I,0,KP)
  360    CONTINUE
         DO 400 J = 1, Q2
            DO 380 I = 0, P - 1
               B(I,K,0) = B(I,K,0) + A(I,J,K)
               B(I,KP,Q-1) = B(I,KP,Q-1) + A(I,J,KP)
  380       CONTINUE
  400    CONTINUE
         DO 440 L = 1, Q2
            DO 420 I = 0, P - 1
               TEMP1 = B(I,K,L)
               TEMP2 = B(I,KP,Q-L-1)
               TEMP3 = B(I,KP,L-1)
               B(I,K,L) = B(I,K,L) - B(I,K,Q-L)
               B(I,KP,Q-L-1) = B(I,KP,Q-L-1) - B(I,KP,L-1)
               B(I,KP,L-1) = TEMP1 + B(I,K,Q-L)
               B(I,K,Q-L) = -TEMP2 - TEMP3
  420       CONTINUE
  440    CONTINUE
  460 CONTINUE
C
C     Code for K=R/2 when R is even not needed - prime factors
C     always treated last, so R = 1 or R is odd!
C
      RETURN
      END
      SUBROUTINE C06FPS(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix six real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      real*8  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:5,0:R-1), B(0:P-1,0:R-1,0:5),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5)
C     .. Local Scalars ..
      real*8  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1, T1I, T1R, T2, T2I, T2R, T3, T3I, T3R, T4,
     *                  T5, T6, U0, U0I, U0R, U1I, U1R, U2I, U2R, UI,
     *                  UR, V0, V0I, V0R, V1I, V1R, V2I, V2R, VI, VR,
     *                  X1P, X2P, X3P, X4P, X5P, Y1P, Y2P, Y3P, Y4P, Y5P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,2,0) + A(I,4,0)
         UR = A(I,0,0) - 0.5D0*T1
         UI = -SIN60*(A(I,2,0)-A(I,4,0))
         U0 = A(I,0,0) + T1
         T1 = A(I,5,0) + A(I,1,0)
         VR = A(I,3,0) - 0.5D0*T1
         VI = -SIN60*(A(I,5,0)-A(I,1,0))
         V0 = A(I,3,0) + T1
         B(I,0,0) = U0 + V0
         B(I,0,1) = UR - VR
         B(I,0,2) = UR + VR
         B(I,0,3) = U0 - V0
         B(I,0,4) = -UI - VI
         B(I,0,5) = UI - VI
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               X4P = COSINE(K,4)*A(I,4,K) - SINE(K,4)*A(I,4,R-K)
               Y4P = COSINE(K,4)*A(I,4,R-K) + SINE(K,4)*A(I,4,K)
               X5P = COSINE(K,5)*A(I,5,K) - SINE(K,5)*A(I,5,R-K)
               Y5P = COSINE(K,5)*A(I,5,R-K) + SINE(K,5)*A(I,5,K)
               T1R = X2P + X4P
               T1I = Y2P + Y4P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,R-K) - 0.5D0*T1I
               T3R = SIN60*(X2P-X4P)
               T3I = SIN60*(Y2P-Y4P)
               U0R = A(I,0,K) + T1R
               U0I = A(I,0,R-K) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X5P + X1P
               T1I = Y5P + Y1P
               T2R = X3P - 0.5D0*T1R
               T2I = Y3P - 0.5D0*T1I
               T3R = SIN60*(X5P-X1P)
               T3I = SIN60*(Y5P-Y1P)
               V0R = X3P + T1R
               V0I = Y3P + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               B(I,K,0) = U0R + V0R
               B(I,R-K,0) = U2R - V2R
               B(I,K,1) = U1R - V1R
               B(I,R-K,1) = U1R + V1R
               B(I,K,2) = U2R + V2R
               B(I,R-K,2) = U0R - V0R
               B(I,K,3) = -U0I + V0I
               B(I,R-K,3) = U2I + V2I
               B(I,K,4) = -U1I - V1I
               B(I,R-K,4) = U1I - V1I
               B(I,K,5) = -U2I + V2I
               B(I,R-K,5) = U0I + V0I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            C5K = COSINE(K,5)
            S5K = SINE(K,5)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               X4P = C4K*A(I,4,K) - S4K*A(I,4,KP)
               Y4P = C4K*A(I,4,KP) + S4K*A(I,4,K)
               X5P = C5K*A(I,5,K) - S5K*A(I,5,KP)
               Y5P = C5K*A(I,5,KP) + S5K*A(I,5,K)
               T1R = X2P + X4P
               T1I = Y2P + Y4P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,KP) - 0.5D0*T1I
               T3R = SIN60*(X2P-X4P)
               T3I = SIN60*(Y2P-Y4P)
               U0R = A(I,0,K) + T1R
               U0I = A(I,0,KP) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X5P + X1P
               T1I = Y5P + Y1P
               T2R = X3P - 0.5D0*T1R
               T2I = Y3P - 0.5D0*T1I
               T3R = SIN60*(X5P-X1P)
               T3I = SIN60*(Y5P-Y1P)
               V0R = X3P + T1R
               V0I = Y3P + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               B(I,K,0) = U0R + V0R
               B(I,KP,0) = U2R - V2R
               B(I,K,1) = U1R - V1R
               B(I,KP,1) = U1R + V1R
               B(I,K,2) = U2R + V2R
               B(I,KP,2) = U0R - V0R
               B(I,K,3) = -U0I + V0I
               B(I,KP,3) = U2I + V2I
               B(I,K,4) = -U1I - V1I
               B(I,KP,4) = U1I - V1I
               B(I,K,5) = -U2I + V2I
               B(I,KP,5) = U0I + V0I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = A(I,2,R2) - A(I,4,R2)
            T2 = A(I,0,R2) + 0.5D0*T1
            T3 = SIN60*(A(I,2,R2)+A(I,4,R2))
            T4 = A(I,1,R2) + A(I,5,R2)
            T5 = -A(I,3,R2) - 0.5D0*T4
            T6 = SIN60*(A(I,1,R2)-A(I,5,R2))
            B(I,R2,0) = T2 + T6
            B(I,R2,1) = A(I,0,R2) - T1
            B(I,R2,2) = T2 - T6
            B(I,R2,3) = T5 + T3
            B(I,R2,4) = A(I,3,R2) - T4
            B(I,R2,5) = T5 - T3
  120    CONTINUE
      END IF
C
      RETURN
      END
	
      SUBROUTINE C06FPT(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix five real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      real*8  R54, SIN72, S36S72
      PARAMETER         (R54=0.559016994374947424102293417182819D0,
     *                  SIN72=0.951056516295153572116439333379382D0,
     *                  S36S72=0.618033988749894848204586834365638D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:4,0:R-1), B(0:P-1,0:R-1,0:4),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4)
C     .. Local Scalars ..
      real*8  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T1, T10I,
     *                  T10R, T11I, T11R, T1I, T1R, T2, T2I, T2R, T3,
     *                  T3I, T3R, T4, T4I, T4R, T5, T5I, T5R, T6, T6I,
     *                  T6R, T7, T7I, T7R, T8I, T8R, T9I, T9R, X1P, X2P,
     *                  X3P, X4P, Y1P, Y2P, Y3P, Y4P
      INTEGER           I, K, KP
c                              , R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,1,0) + A(I,4,0)
         T2 = A(I,2,0) + A(I,3,0)
         T3 = SIN72*(A(I,1,0)-A(I,4,0))
         T4 = SIN72*(A(I,2,0)-A(I,3,0))
         T5 = T1 + T2
         T6 = R54*(T1-T2)
         T7 = A(I,0,0) - 0.25D0*T5
         B(I,0,0) = A(I,0,0) + T5
         B(I,0,1) = T7 + T6
         B(I,0,2) = T7 - T6
         B(I,0,3) = -S36S72*T3 + T4
         B(I,0,4) = -T3 - S36S72*T4
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               X4P = COSINE(K,4)*A(I,4,K) - SINE(K,4)*A(I,4,R-K)
               Y4P = COSINE(K,4)*A(I,4,R-K) + SINE(K,4)*A(I,4,K)
               T1R = X1P + X4P
               T1I = Y1P + Y4P
               T2R = X2P + X3P
               T2I = Y2P + Y3P
               T3R = SIN72*(X1P-X4P)
               T3I = SIN72*(Y1P-Y4P)
               T4R = SIN72*(X2P-X3P)
               T4I = SIN72*(Y2P-Y3P)
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,0,K) - 0.25D0*T5R
               T7I = A(I,0,R-K) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               B(I,K,0) = A(I,0,K) + T5R
               B(I,R-K,0) = T8R - T10I
               B(I,K,1) = T8R + T10I
               B(I,R-K,1) = T9R - T11I
               B(I,K,2) = T9R + T11I
               B(I,R-K,2) = T9I - T11R
               B(I,K,3) = -T9I - T11R
               B(I,R-K,3) = T8I - T10R
               B(I,K,4) = -T8I - T10R
               B(I,R-K,4) = A(I,0,R-K) + T5I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               X4P = C4K*A(I,4,K) - S4K*A(I,4,KP)
               Y4P = C4K*A(I,4,KP) + S4K*A(I,4,K)
               T1R = X1P + X4P
               T1I = Y1P + Y4P
               T2R = X2P + X3P
               T2I = Y2P + Y3P
               T3R = SIN72*(X1P-X4P)
               T3I = SIN72*(Y1P-Y4P)
               T4R = SIN72*(X2P-X3P)
               T4I = SIN72*(Y2P-Y3P)
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,0,K) - 0.25D0*T5R
               T7I = A(I,0,KP) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               B(I,K,0) = A(I,0,K) + T5R
               B(I,KP,0) = T8R - T10I
               B(I,K,1) = T8R + T10I
               B(I,KP,1) = T9R - T11I
               B(I,K,2) = T9R + T11I
               B(I,KP,2) = T9I - T11R
               B(I,K,3) = -T9I - T11R
               B(I,KP,3) = T8I - T10R
               B(I,K,4) = -T8I - T10R
               B(I,KP,4) = A(I,0,KP) + T5I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = SIN72*(A(I,1,R2)+A(I,4,R2))
C           T2 = SIN72*(A(I,2,R2)+A(I,3,R2))
C           T3 = A(I,1,R2) - A(I,4,R2)
C           T4 = A(I,2,R2) - A(I,3,R2)
C           T5 = T4 - T3
C           T6 = R54*(T4+T3)
C           T7 = A(I,0,R2) - 0.25D0*T5
C           B(I,R2,0) = T7 + T6
C           B(I,R2,1) = T7 - T6
C           B(I,R2,2) = A(I,0,R2) + T5
C           B(I,R2,3) = -T1 + S36S72*T2
C           B(I,R2,4) = -S36S72*T1 - T2
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FPU(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix four real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      real*8  ROOT2I
      PARAMETER         (ROOT2I=0.707106781186547524400844362104849D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:3,0:R-1), B(0:P-1,0:R-1,0:3),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3)
C     .. Local Scalars ..
      real*8  C2K, C3K, CK, S2K, S3K, SK, T1, T1I, T1R, T2,
     *                  T2I, T2R, T3I, T3R, T4I, T4R, X1P, X2P, X3P,
     *                  Y1P, Y2P, Y3P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,0) + A(I,2,0)
         T2 = A(I,1,0) + A(I,3,0)
         B(I,0,0) = T1 + T2
         B(I,0,1) = A(I,0,0) - A(I,2,0)
         B(I,0,2) = T1 - T2
         B(I,0,3) = -A(I,1,0) + A(I,3,0)
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               T1R = A(I,0,K) + X2P
               T1I = A(I,0,R-K) + Y2P
               T2R = X1P + X3P
               T2I = Y1P + Y3P
               T3R = A(I,0,K) - X2P
               T3I = A(I,0,R-K) - Y2P
               T4R = X1P - X3P
               T4I = Y1P - Y3P
               B(I,K,0) = T1R + T2R
               B(I,R-K,0) = T3R - T4I
               B(I,K,1) = T3R + T4I
               B(I,R-K,1) = T1R - T2R
               B(I,K,2) = T2I - T1I
               B(I,R-K,2) = T3I - T4R
               B(I,K,3) = -T3I - T4R
               B(I,R-K,3) = T1I + T2I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               T1R = A(I,0,K) + X2P
               T1I = A(I,0,KP) + Y2P
               T2R = X1P + X3P
               T2I = Y1P + Y3P
               T3R = A(I,0,K) - X2P
               T3I = A(I,0,KP) - Y2P
               T4R = X1P - X3P
               T4I = Y1P - Y3P
               B(I,K,0) = T1R + T2R
               B(I,KP,0) = T3R - T4I
               B(I,K,1) = T3R + T4I
               B(I,KP,1) = T1R - T2R
               B(I,K,2) = T2I - T1I
               B(I,KP,2) = T3I - T4R
               B(I,K,3) = -T3I - T4R
               B(I,KP,3) = T1I + T2I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = ROOT2I*(A(I,1,R2)-A(I,3,R2))
            T2 = ROOT2I*(A(I,1,R2)+A(I,3,R2))
            B(I,R2,0) = A(I,0,R2) + T1
            B(I,R2,1) = A(I,0,R2) - T1
            B(I,R2,2) = A(I,2,R2) - T2
            B(I,R2,3) = -A(I,2,R2) - T2
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FPV(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix three Real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      real*8  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:2,0:R-1), B(0:P-1,0:R-1,0:2),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2)
C     .. Local Scalars ..
      real*8  C2K, CK, S2K, SK, T1, T1I, T1R, T2I, T2R, T3I,
     *                  T3R, X1P, X2P, Y1P, Y2P
      INTEGER           I, K, KP
c                               , R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,1,0) + A(I,2,0)
         B(I,0,0) = A(I,0,0) + T1
         B(I,0,1) = A(I,0,0) - 0.5D0*T1
         B(I,0,2) = -SIN60*(A(I,1,0)-A(I,2,0))
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               T1R = X1P + X2P
               T1I = Y1P + Y2P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,R-K) - 0.5D0*T1I
               T3R = SIN60*(X1P-X2P)
               T3I = SIN60*(Y1P-Y2P)
               B(I,K,0) = A(I,0,K) + T1R
               B(I,R-K,0) = T2R - T3I
               B(I,K,1) = T2R + T3I
               B(I,R-K,1) = T2I - T3R
               B(I,K,2) = -(T2I+T3R)
               B(I,R-K,2) = A(I,0,R-K) + T1I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               T1R = X1P + X2P
               T1I = Y1P + Y2P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,KP) - 0.5D0*T1I
               T3R = SIN60*(X1P-X2P)
               T3I = SIN60*(Y1P-Y2P)
               B(I,K,0) = A(I,0,K) + T1R
               B(I,KP,0) = T2R - T3I
               B(I,K,1) = T2R + T3I
               B(I,KP,1) = T2I - T3R
               B(I,K,2) = -(T2I+T3R)
               B(I,KP,2) = A(I,0,KP) + T1I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,1,R2) - A(I,2,R2)
C           B(I,R2,0) = A(I,0,R2) + 0.5D0*T1
C           B(I,R2,1) = A(I,0,R2) - T1
C           B(I,R2,2) = -SIN60*(A(I,1,R2)+A(I,2,R2))
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FPW(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix two Real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:1,0:R-1), B(0:P-1,0:R-1,0:1),
     *                  COSINE(0:R-1), SINE(0:R-1)
C     .. Local Scalars ..
      real*8  CK, SK, X1, Y1
      INTEGER           I, K, KP
c                               , R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         B(I,0,0) = A(I,0,0) + A(I,1,0)
         B(I,0,1) = A(I,0,0) - A(I,1,0)
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1 = COSINE(K)*A(I,1,K) - SINE(K)*A(I,1,R-K)
               Y1 = COSINE(K)*A(I,1,R-K) + SINE(K)*A(I,1,K)
               B(I,K,0) = A(I,0,K) + X1
               B(I,R-K,0) = A(I,0,K) - X1
               B(I,K,1) = -A(I,0,R-K) + Y1
               B(I,R-K,1) = A(I,0,R-K) + Y1
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K)
            SK = SINE(K)
            DO 80 I = 0, P - 1
               X1 = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1 = CK*A(I,1,KP) + SK*A(I,1,K)
               B(I,K,0) = A(I,0,K) + X1
               B(I,KP,0) = A(I,0,K) - X1
               B(I,K,1) = -A(I,0,KP) + Y1
               B(I,KP,1) = A(I,0,KP) + Y1
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           B(I,R2,0) = A(I,0,R2)
C           B(I,R2,1) = -A(I,1,R2)
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FPX(A,B,M,N,Q,NQ,TRIG)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Real to Hermitian Fast Fourier Transform Kernel Driver
C
C     Mixed-radix, self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NQ
C     .. Array Arguments ..
      real*8  A(0:M*N-1), B(0:M*N-1), TRIG(0:2*N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      real*8  FACTOR
      INTEGER           I, P, QI, R
      LOGICAL           INA
C     .. External Subroutines ..
      EXTERNAL          C06FPR, C06FPS, C06FPT, C06FPU, C06FPV, C06FPW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = N
      R = 1
      IF (N.EQ.1) RETURN
      DO 20 I = NQ, 1, -1
         QI = Q(I)
         P = P/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FPW(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FPV(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FPU(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FPT(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FPS(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FPR(A,B,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FPW(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FPV(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FPU(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FPT(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FPS(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FPR(B,A,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         END IF
         INA = .NOT. INA
         R = R*QI
   20 CONTINUE
C
      FACTOR = 1.0D0/SQRT(DBLE(N))
      IF (INA) THEN
         DO 40 I = 0, M*N - 1
            A(I) = A(I)*FACTOR
   40    CONTINUE
      ELSE
         DO 60 I = 0, M*N - 1
            A(I) = B(I)*FACTOR
   60    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FPY(N,NQ,Q,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Trig function initialisation subroutine
C
C     .. Scalar Arguments ..
      INTEGER           N, NQ
C     .. Array Arguments ..
      real*8  COSINE(0:N-1), SINE(0:N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      real*8  TWOPI, Z
      INTEGER           I, J, K, L, L1, QI, R
C     .. Intrinsic Functions ..
      INTRINSIC         COS, SIN, DBLE
C     .. Executable Statements ..
      TWOPI = 2.0D0 * acos(-1.d0)
      Z = TWOPI/DBLE(N)
C
      R = N
      L = 0
C
      DO 80 I = 1, NQ
         QI = Q(I)
         R = R/QI
         L1 = L
         DO 40 J = 1, QI - 1
            DO 20 K = 0, R - 1
               COSINE(L) = Z*DBLE(J*K)
               L = L + 1
   20       CONTINUE
   40    CONTINUE
         IF (QI.GE.7) THEN
            L = L1
            DO 60 J = 1, QI - 1
               COSINE(L) = Z*J*R
               L = L + R
   60       CONTINUE
         END IF
         Z = Z*QI
   80 CONTINUE
C
      DO 100 I = 0, N - 2
         SINE(I) = -SIN(COSINE(I))
         COSINE(I) = COS(COSINE(I))
  100 CONTINUE
C
C     Check on consistency of N and TRIG array --
C
      COSINE(N-1) = DBLE(N)
      SINE(N-1) = DBLE(N)
C
      RETURN
      END
      SUBROUTINE C06FPZ(N,NQ,Q)
      
CVD$R NOVECTOR
C     .. Scalar Arguments ..
      INTEGER           N, NQ
C     .. Array Arguments ..
      INTEGER           Q(30)
C     .. Local Scalars ..
      INTEGER           I, K, L, NN
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      NN = N
      K = 0
C
C     Trap the special case N = 1
C
      IF (N.EQ.1) THEN
         NQ = 1
         Q(1) = 1
         RETURN
      END IF
C
C     Factors of 6 --
C
   20 IF (MOD(NN,6).NE.0) GO TO 40
      K = K + 1
      Q(K) = 6
      NN = NN/6
      IF (NN.EQ.1) GO TO 160
      GO TO 20
C
C     Factors of 4 --
C
   40 IF (MOD(NN,4).NE.0) GO TO 60
      K = K + 1
      Q(K) = 4
      NN = NN/4
      IF (NN.EQ.1) GO TO 160
      GO TO 40
C
C     Factors of 2 --
C
   60 IF (MOD(NN,2).NE.0) GO TO 80
      K = K + 1
      Q(K) = 2
      NN = NN/2
      IF (NN.EQ.1) GO TO 160
      GO TO 60
C
C     Factors of 3 --
C
   80 IF (MOD(NN,3).NE.0) GO TO 100
      K = K + 1
      Q(K) = 3
      NN = NN/3
      IF (NN.EQ.1) GO TO 160
      GO TO 80
C
C     Remaining odd factors --
C
  100 L = 5
      I = 2
C
C     I is alternatively 2 or 4 --
C
  120 IF (MOD(NN,L).NE.0) GO TO 140
      K = K + 1
      Q(K) = L
      NN = NN/L
      IF (NN.EQ.1) GO TO 160
      GO TO 120
  140 L = L + I
      I = 6 - I
      GO TO 120
  160 NQ = K
C
      RETURN
      END

      SUBROUTINE C06HAF(M,N,X,INIT,TRIG,WORK,IFAIL)
      
C
C     C06HAF computes multiple Fourier sine transforms of sequences
C     of real data using the multiple real transform kernel C06FPX,
C     and pre- and post-processing steps described by Swarztrauber.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06HAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      real*8  TRIG(2*N), WORK(1:M,0:N-1), X(1:M,1:N)
C     .. Local Scalars ..
      real*8  PIBYN, ROOT2, ROOTN
c                                ,PI
      INTEGER           I, IERROR, J, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FPX
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE, SIN, SQRT
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0 .AND. N.GT.1) THEN
         PIBYN = acos(-1.d0) /DBLE(N)
         ROOTN = SQRT(DBLE(N))
         ROOT2 = SQRT(2.0D0)
C
C        Calculate the auxilliary array D and store it in the WORK array
C
C        Set the first column of the array to be transformed to zero
C
         DO 20 I = 1, M
            WORK(I,0) = 0.0D0
   20    CONTINUE
C
C        Set up the rest of the D array
C
         DO 60 J = 1, N - 1
            DO 40 I = 1, M
               WORK(I,J) = 0.5D0*(X(I,J)-X(I,N-J)) + SIN(DBLE(J)*PIBYN)
     *                     *(X(I,J)+X(I,N-J))
   40       CONTINUE
   60    CONTINUE
C
C        Transform the D array, stored in WORK, using X as the workspace
C        array
C
         CALL C06FPX(WORK,X,M,N,Q,NQ,TRIG)
C
C        Now extract the sine transform from the Fourier transform in
C        WORK, putting the required transform in X
C
         DO 80 I = 1, M
            X(I,1) = WORK(I,0)/ROOT2
   80    CONTINUE
         IF (MOD(N,2).EQ.0) THEN
            DO 140 J = 1, N/2 - 1
               DO 100 I = 1, M
                  X(I,2*J) = -ROOT2*WORK(I,N-J)
  100          CONTINUE
               DO 120 I = 1, M
                  X(I,2*J+1) = ROOT2*WORK(I,J) + X(I,2*J-1)
  120          CONTINUE
  140       CONTINUE
         ELSE IF (MOD(N,2).EQ.1) THEN
            DO 180 J = 1, (N-1)/2 - 1
               DO 160 I = 1, M
                  X(I,2*J) = -ROOT2*WORK(I,N-J)
                  X(I,2*J+1) = ROOT2*WORK(I,J) + X(I,2*J-1)
  160          CONTINUE
  180       CONTINUE
            DO 200 I = 1, M
               X(I,N-1) = -ROOT2*WORK(I,(N+1)/2)
  200       CONTINUE
         END IF
C
C        Set the arbitrary last element of the array to zero
C
         DO 220 I = 1, M
            X(I,N) = 0.0D0
  220    CONTINUE
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
***   IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,' but TRIG array never initialised')
99995 FORMAT (' ** INIT = ',A1,' but N and TRIG array incompatible')
      END



      SUBROUTINE C06HBF(M,N,X,INIT,TRIG,WORK,IFAIL)
      
C
C     C06HBF computes multiple Fourier cosine transforms of sequences
C     of real data using the multiple real transform kernel C06FPX,
C     and pre- and post-processing steps described by Swarztrauber.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06HBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      real*8  TRIG(2*N), WORK(1:M,0:N-1), X(1:M,0:N)
C     .. Local Scalars ..
      real*8  PIBYN, ROOT2, ROOTN, SUM
c                                     , PI
      INTEGER           I, IERROR, J, N2, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FPX
C     .. Intrinsic Functions ..
      INTRINSIC         COS, MOD, DBLE, SIN, SQRT
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0 .AND. N.GT.1) THEN
         PIBYN = acos(-1.d0) / DBLE(N)
         ROOTN = SQRT(DBLE(N))
         ROOT2 = SQRT(2.0D0)
C
C        Calculate the auxilliary array D and store it in the WORK array
C
C        Set the first column of the array to be transformed to the
C        required value computed from the input array
C
         DO 20 I = 1, M
            WORK(I,0) = 0.5D0*(X(I,0)+X(I,N))
   20    CONTINUE
C
C        Set up the rest of the D array
C
         DO 60 J = 1, N - 1
            DO 40 I = 1, M
               WORK(I,J) = 0.5D0*(X(I,J)+X(I,N-J)) - SIN(DBLE(J)*PIBYN)
     *                     *(X(I,J)-X(I,N-J))
   40       CONTINUE
   60    CONTINUE
C
C        The first element of the transform is generated here and
C        stored in X(1:M,N), which is not used as workspace by C06FPX
C
         DO 100 I = 1, M
            SUM = 0.0D0
            SUM = SUM + 0.5D0*X(I,0)
            DO 80 J = 1, N - 1
               SUM = SUM + X(I,J)*COS(DBLE(J)*PIBYN)
   80       CONTINUE
            SUM = SUM - 0.5D0*X(I,N)
            X(I,N) = ROOT2*SUM/ROOTN
  100    CONTINUE
C
C        Transform the D array, stored in WORK, using X as the workspace
C        array
C
         CALL C06FPX(WORK,X,M,N,Q,NQ,TRIG)
C
C        Now extract the sine transform from the Fourier transform in
C        WORK, putting the required transform in X
C
         DO 120 I = 1, M
            X(I,0) = ROOT2*WORK(I,0)
  120    CONTINUE
         DO 140 I = 1, M
            X(I,1) = X(I,N)
  140    CONTINUE
C
         IF (MOD(N,2).EQ.0) THEN
            DO 180 J = 1, N/2 - 1
               DO 160 I = 1, M
                  X(I,2*J) = ROOT2*WORK(I,J)
                  X(I,2*J+1) = X(I,2*J-1) - ROOT2*WORK(I,N-J)
  160          CONTINUE
  180       CONTINUE
            N2 = N/2
            DO 200 I = 1, M
               X(I,N) = ROOT2*WORK(I,N2)
  200       CONTINUE
         ELSE IF (MOD(N,2).EQ.1) THEN
            DO 240 J = 1, (N-1)/2
               DO 220 I = 1, M
                  X(I,2*J) = ROOT2*WORK(I,J)
                  X(I,2*J+1) = X(I,2*J-1) - ROOT2*WORK(I,N-J)
  220          CONTINUE
  240       CONTINUE
         END IF
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
***   IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,' but TRIG array never initialised')
99995 FORMAT (' ** INIT = ',A1,' but N and TRIG array incompatible')
      END

      SUBROUTINE C06FPF(M,N,X,INIT,TRIG,WORK,IFAIL)
      
CVD$R NOVECTOR
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FPF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      real*8  TRIG(2*N), WORK(M*N), X(M*N)
C     .. Local Scalars ..
      INTEGER           IERROR, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FPX
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0) THEN
         CALL C06FPX(X,WORK,M,N,Q,NQ,TRIG)
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
***   IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,', but TRIG array never initialized')
99995 FORMAT (' ** INIT = ',A1,', but N and TRIG array incompatible')
      END
      SUBROUTINE C06FQF(M,N,X,INIT,TRIG,WORK,IFAIL)
      
CVD$R NOVECTOR
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FQF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      real*8  TRIG(2*N), WORK(M*N), X(M*N)
C     .. Local Scalars ..
      INTEGER           IERROR, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FQX
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0) THEN
         CALL C06FQX(X,WORK,M,N,Q,NQ,TRIG)
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
***   IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,', but TRIG array never initialized')
99995 FORMAT (' ** INIT = ',A1,', but N and TRIG array incompatible')
      END
      SUBROUTINE C06GQF(M,N,X,IFAIL)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06GQF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      real*8  X(M*N)
C     .. Local Scalars ..
      INTEGER           I, IERROR, N2, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. Executable Statements ..
      IF (M.LE.0) GO TO 40
      IF (N.LE.0) GO TO 60
      IERROR = 0
      N2 = (N+4)/2
      DO 20 I = (N2-1)*M + 1, N*M
         X(I) = -X(I)
   20 CONTINUE
      GO TO 80
   40 IERROR = 1
      WRITE (REC(1),FMT=99999) M
      GO TO 80
   60 IERROR = 2
      WRITE (REC(1),FMT=99998) N
   80 NREC = 1
***   IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
      END


***********************

      SUBROUTINE C06FQQ(A,M,N)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      real*8  A(0:M-1,0:N-1)
C     .. Local Scalars ..
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 L = 0, M - 1
         A(L,0) = 0.5D0*A(L,0)
   20 CONTINUE
      IF (MOD(N,2).EQ.0) THEN
         DO 40 L = 0, M - 1
            A(L,N/2) = 0.5D0*A(L,N/2)
   40    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FQR(A,B,P,Q,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Hermitian to real fast Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:R-1,0:Q-1), B(0:P-1,0:Q-1,0:R-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1)
C     .. Local Scalars ..
      real*8  BI, BR, TEMP, TEMP1, TEMPI, TEMPR, TWOPI
      INTEGER           I, INDX, J, K, KP, L, Q2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      Q2 = (Q-1)/2
      TWOPI = 2.0D0 * acos(-1.d0)
C
C     Code for K=0 --
C
      DO 80 L = 1, Q2
         DO 20 I = 0, P - 1
            B(I,L,0) = A(I,0,0)
            B(I,Q-L,0) = 0.0D0
   20    CONTINUE
         DO 60 J = 1, Q2
            INDX = MOD(J*L,Q)
            DO 40 I = 0, P - 1
               B(I,L,0) = B(I,L,0) + A(I,0,J)*COSINE(0,INDX)
               B(I,Q-L,0) = B(I,Q-L,0) - A(I,0,Q-J)*SINE(0,INDX)
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      DO 100 I = 0, P - 1
         B(I,0,0) = A(I,0,0)
  100 CONTINUE
      DO 140 J = 1, Q2
         DO 120 I = 0, P - 1
            B(I,0,0) = B(I,0,0) + A(I,0,J)
  120    CONTINUE
  140 CONTINUE
      DO 180 J = 1, Q2
         DO 160 I = 0, P - 1
            TEMP = B(I,J,0)
            B(I,J,0) = B(I,J,0) + B(I,Q-J,0)
            B(I,Q-J,0) = TEMP - B(I,Q-J,0)
  160    CONTINUE
  180 CONTINUE
C
C     Code for general K --
C
      DO 460 K = 1, (R-1)/2
         KP = R - K
         DO 220 J = 1, Q2
            DO 200 I = 0, P - 1
               TEMPR = A(I,K,J)
               TEMPI = A(I,KP,Q-J-1)
               A(I,K,J) = TEMPR + A(I,KP,J-1)
               A(I,KP,Q-J-1) = TEMPI - A(I,K,Q-J)
               A(I,KP,J-1) = TEMPR - A(I,KP,J-1)
               A(I,K,Q-J) = -TEMPI - A(I,K,Q-J)
  200       CONTINUE
  220    CONTINUE
         DO 300 L = 1, Q2
            DO 240 I = 0, P - 1
               B(I,L,K) = A(I,K,0)
               B(I,L,KP) = A(I,KP,Q-1)
               B(I,Q-L,K) = 0.0D0
               B(I,Q-L,KP) = 0.0D0
  240       CONTINUE
            DO 280 J = 1, Q2
               INDX = MOD(J*L,Q)
               DO 260 I = 0, P - 1
                  B(I,L,K) = B(I,L,K) + A(I,K,J)*COSINE(0,INDX)
                  B(I,L,KP) = B(I,L,KP) + A(I,KP,Q-J-1)*COSINE(0,INDX)
                  B(I,Q-L,K) = B(I,Q-L,K) + A(I,KP,J-1)*SINE(0,INDX)
                  B(I,Q-L,KP) = B(I,Q-L,KP) - A(I,K,Q-J)*SINE(0,INDX)
  260          CONTINUE
  280       CONTINUE
  300    CONTINUE
         DO 320 I = 0, P - 1
            B(I,0,K) = A(I,K,0)
            B(I,0,KP) = A(I,KP,Q-1)
  320    CONTINUE
         DO 360 J = 1, Q2
            DO 340 I = 0, P - 1
               B(I,0,K) = B(I,0,K) + A(I,K,J)
               B(I,0,KP) = B(I,0,KP) + A(I,KP,Q-J-1)
  340       CONTINUE
  360    CONTINUE
         DO 400 J = 1, Q2
            DO 380 I = 0, P - 1
               TEMPR = B(I,J,K)
               TEMPI = B(I,J,KP)
               B(I,J,K) = TEMPR - B(I,Q-J,KP)
               B(I,J,KP) = TEMPI + B(I,Q-J,K)
               TEMP1 = B(I,Q-J,K)
               B(I,Q-J,K) = TEMPR + B(I,Q-J,KP)
               B(I,Q-J,KP) = TEMPI - TEMP1
  380       CONTINUE
  400    CONTINUE
         DO 440 J = 1, Q - 1
            DO 420 I = 0, P - 1
               BR = B(I,J,K)
               BI = B(I,J,KP)
               B(I,J,K) = COSINE(K,J)*BR - SINE(K,J)*BI
               B(I,J,KP) = COSINE(K,J)*BI + SINE(K,J)*BR
  420       CONTINUE
  440    CONTINUE
  460 CONTINUE
C
C     Code for K=R/2 when R is even not needed - prime factors
C     always treated last, so R = 1 or R is odd!
C
      RETURN
      END
      SUBROUTINE C06FQS(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     Radix six Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      real*8  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:R-1,0:5), B(0:P-1,0:5,0:R-1),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5)
C     .. Local Scalars ..
      real*8  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1, T1I, T1R, T2, T2I, T2R, T3, T3I, T3R, T4,
     *                  T5, T6, U0, U0I, U0R, U1, U1I, U1R, U2, U2I,
     *                  U2R, V0, V0I, V0R, V1, V1I, V1R, V2, V2I, V2R,
     *                  X0P, X1P, X2P, X3P, X4P, X5P, Y0P, Y1P, Y2P,
     *                  Y3P, Y4P, Y5P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,2)
         T2 = A(I,0,0) - 0.5D0*T1
         T3 = SIN60*A(I,0,4)
         U0 = A(I,0,0) + T1
         U1 = T2 + T3
         U2 = T2 - T3
         T1 = A(I,0,1)
         T2 = A(I,0,3) - 0.5D0*T1
         T3 = -SIN60*A(I,0,5)
         V0 = A(I,0,3) + T1
         V1 = T2 + T3
         V2 = T2 - T3
         B(I,0,0) = U0 + V0
         B(I,1,0) = U1 - V1
         B(I,2,0) = U2 + V2
         B(I,3,0) = U0 - V0
         B(I,4,0) = U1 + V1
         B(I,5,0) = U2 - V2
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,2) + A(I,R-K,1)
               T1I = A(I,R-K,3) - A(I,K,4)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,R-K,5) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,2)-A(I,R-K,1))
               T3I = SIN60*(A(I,R-K,3)+A(I,K,4))
               U0R = A(I,K,0) + T1R
               U0I = A(I,R-K,5) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = A(I,R-K,0) + A(I,K,1)
               T1I = -A(I,K,5) + A(I,R-K,4)
               T2R = A(I,R-K,2) - 0.5D0*T1R
               T2I = -A(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(A(I,R-K,0)-A(I,K,1))
               T3I = SIN60*(-A(I,K,5)-A(I,R-K,4))
               V0R = A(I,R-K,2) + T1R
               V0I = -A(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               X0P = U0R + V0R
               Y0P = U0I + V0I
               X1P = U1R - V1R
               Y1P = U1I - V1I
               X2P = U2R + V2R
               Y2P = U2I + V2I
               X3P = U0R - V0R
               Y3P = U0I - V0I
               X4P = U1R + V1R
               Y4P = U1I + V1I
               X5P = U2R - V2R
               Y5P = U2I - V2I
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
               B(I,4,K) = COSINE(K,4)*X4P - SINE(K,4)*Y4P
               B(I,4,R-K) = COSINE(K,4)*Y4P + SINE(K,4)*X4P
               B(I,5,K) = COSINE(K,5)*X5P - SINE(K,5)*Y5P
               B(I,5,R-K) = COSINE(K,5)*Y5P + SINE(K,5)*X5P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            C5K = COSINE(K,5)
            S5K = SINE(K,5)
            DO 80 I = 0, P - 1
               T1R = A(I,K,2) + A(I,KP,1)
               T1I = A(I,KP,3) - A(I,K,4)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,KP,5) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,2)-A(I,KP,1))
               T3I = SIN60*(A(I,KP,3)+A(I,K,4))
               U0R = A(I,K,0) + T1R
               U0I = A(I,KP,5) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = A(I,KP,0) + A(I,K,1)
               T1I = -A(I,K,5) + A(I,KP,4)
               T2R = A(I,KP,2) - 0.5D0*T1R
               T2I = -A(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(A(I,KP,0)-A(I,K,1))
               T3I = SIN60*(-A(I,K,5)-A(I,KP,4))
               V0R = A(I,KP,2) + T1R
               V0I = -A(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               X0P = U0R + V0R
               Y0P = U0I + V0I
               X1P = U1R - V1R
               Y1P = U1I - V1I
               X2P = U2R + V2R
               Y2P = U2I + V2I
               X3P = U0R - V0R
               Y3P = U0I - V0I
               X4P = U1R + V1R
               Y4P = U1I + V1I
               X5P = U2R - V2R
               Y5P = U2I - V2I
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
               B(I,4,K) = C4K*X4P - S4K*Y4P
               B(I,4,KP) = C4K*Y4P + S4K*X4P
               B(I,5,K) = C5K*X5P - S5K*Y5P
               B(I,5,KP) = C5K*Y5P + S5K*X5P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = A(I,R2,0) + A(I,R2,2)
            T2 = A(I,R2,5) + A(I,R2,3)
            T3 = A(I,R2,1) - 0.5D0*T1
            T4 = A(I,R2,4) + 0.5D0*T2
            T5 = SIN60*(A(I,R2,0)-A(I,R2,2))
            T6 = SIN60*(A(I,R2,5)-A(I,R2,3))
            B(I,0,R2) = A(I,R2,1) + T1
            B(I,1,R2) = T4 + T5
            B(I,2,R2) = T6 - T3
            B(I,3,R2) = T2 - A(I,R2,4)
            B(I,4,R2) = T3 + T6
            B(I,5,R2) = T4 - T5
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FQT(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix five Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      real*8  R54, SIN36, SIN72, S36S72
      PARAMETER         (R54=0.559016994374947424102293417182819D0,
     *                  SIN36=0.587785252292473129168705954639073D0,
     *                  SIN72=0.951056516295153572116439333379382D0,
     *                  S36S72=0.618033988749894848204586834365638D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:R-1,0:4), B(0:P-1,0:4,0:R-1),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4)
C     .. Local Scalars ..
      real*8  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T1, T10,
     *                  T10I, T10R, T11, T11I, T11R, T1I, T1R, T2, T2I,
     *                  T2R, T3, T3I, T3R, T4, T4I, T4R, T5, T5I, T5R,
     *                  T6, T6I, T6R, T7, T7I, T7R, T8, T8I, T8R, T9,
     *                  T9I, T9R, X0P, X1P, X2P, X3P, X4P, Y0P, Y1P,
     *                  Y2P, Y3P, Y4P
      INTEGER           I, K, KP
c                               , R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,1)
         T2 = A(I,0,2)
         T3 = SIN72*A(I,0,4)
         T4 = SIN72*A(I,0,3)
         T5 = T1 + T2
         T6 = R54*(T1-T2)
         T7 = A(I,0,0) - 0.25D0*T5
         T8 = T7 + T6
         T9 = T7 - T6
         T10 = T3 + S36S72*T4
         T11 = S36S72*T3 - T4
         B(I,0,0) = A(I,0,0) + T5
         B(I,1,0) = T8 + T10
         B(I,2,0) = T9 + T11
         B(I,3,0) = T9 - T11
         B(I,4,0) = T8 - T10
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,1) + A(I,R-K,0)
               T1I = A(I,R-K,3) - A(I,K,4)
               T2R = A(I,K,2) + A(I,R-K,1)
               T2I = A(I,R-K,2) - A(I,K,3)
               T3R = SIN72*(A(I,K,1)-A(I,R-K,0))
               T3I = SIN72*(A(I,R-K,3)+A(I,K,4))
               T4R = SIN72*(A(I,K,2)-A(I,R-K,1))
               T4I = SIN72*(A(I,R-K,2)+A(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,K,0) - 0.25D0*T5R
               T7I = A(I,R-K,4) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               X0P = A(I,K,0) + T5R
               Y0P = A(I,R-K,4) + T5I
               X1P = T8R + T10I
               Y1P = T8I - T10R
               X2P = T9R + T11I
               Y2P = T9I - T11R
               X3P = T9R - T11I
               Y3P = T9I + T11R
               X4P = T8R - T10I
               Y4P = T8I + T10R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
               B(I,4,K) = COSINE(K,4)*X4P - SINE(K,4)*Y4P
               B(I,4,R-K) = COSINE(K,4)*Y4P + SINE(K,4)*X4P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            DO 80 I = 0, P - 1
               T1R = A(I,K,1) + A(I,KP,0)
               T1I = A(I,KP,3) - A(I,K,4)
               T2R = A(I,K,2) + A(I,KP,1)
               T2I = A(I,KP,2) - A(I,K,3)
               T3R = SIN72*(A(I,K,1)-A(I,KP,0))
               T3I = SIN72*(A(I,KP,3)+A(I,K,4))
               T4R = SIN72*(A(I,K,2)-A(I,KP,1))
               T4I = SIN72*(A(I,KP,2)+A(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,K,0) - 0.25D0*T5R
               T7I = A(I,KP,4) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               X0P = A(I,K,0) + T5R
               Y0P = A(I,KP,4) + T5I
               X1P = T8R + T10I
               Y1P = T8I - T10R
               X2P = T9R + T11I
               Y2P = T9I - T11R
               X3P = T9R - T11I
               Y3P = T9I + T11R
               X4P = T8R - T10I
               Y4P = T8I + T10R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
               B(I,4,K) = C4K*X4P - S4K*Y4P
               B(I,4,KP) = C4K*Y4P + S4K*X4P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,R2,0) + A(I,R2,1)
C           T2 = 0.25D0*T1 - A(I,R2,2)
C           T3 = R54*(A(I,R2,0)-A(I,R2,1))
C           T4 = SIN36*A(I,R2,4) + SIN72*A(I,R2,3)
C           T5 = SIN72*A(I,R2,4) - SIN36*A(I,R2,3)
C           T6 = T2 + T3
C           T7 = T2 - T3
C           B(I,0,R2) = T1 + A(I,R2,2)
C           B(I,1,R2) = T4 + T6
C           B(I,2,R2) = T5 - T7
C           B(I,3,R2) = T5 + T7
C           B(I,4,R2) = T4 - T6
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FQU(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix four Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      real*8  ROOT2I
      PARAMETER         (ROOT2I=0.707106781186547524400844362104849D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:R-1,0:3), B(0:P-1,0:3,0:R-1),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3)
C     .. Local Scalars ..
      real*8  C2K, C3K, CK, S2K, S3K, SK, T1, T1I, T1R, T2,
     *                  T2I, T2R, T3, T3I, T3R, T4, T4I, T4R, X0P, X1P,
     *                  X2P, X3P, Y0P, Y1P, Y2P, Y3P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,0) + A(I,0,2)
         T2 = A(I,0,1)
         T3 = A(I,0,0) - A(I,0,2)
         T4 = A(I,0,3)
         B(I,0,0) = T1 + T2
         B(I,1,0) = T3 + T4
         B(I,2,0) = T1 - T2
         B(I,3,0) = T3 - T4
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,0) + A(I,R-K,1)
               T1I = A(I,R-K,3) - A(I,K,2)
               T2R = A(I,K,1) + A(I,R-K,0)
               T2I = A(I,R-K,2) - A(I,K,3)
               T3R = A(I,K,0) - A(I,R-K,1)
               T3I = A(I,R-K,3) + A(I,K,2)
               T4R = A(I,K,1) - A(I,R-K,0)
               T4I = A(I,R-K,2) + A(I,K,3)
               X0P = T1R + T2R
               Y0P = T1I + T2I
               X1P = T3R + T4I
               Y1P = T3I - T4R
               X2P = T1R - T2R
               Y2P = T1I - T2I
               X3P = T3R - T4I
               Y3P = T3I + T4R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            DO 80 I = 0, P - 1
               T1R = A(I,K,0) + A(I,KP,1)
               T1I = A(I,KP,3) - A(I,K,2)
               T2R = A(I,K,1) + A(I,KP,0)
               T2I = A(I,KP,2) - A(I,K,3)
               T3R = A(I,K,0) - A(I,KP,1)
               T3I = A(I,KP,3) + A(I,K,2)
               T4R = A(I,K,1) - A(I,KP,0)
               T4I = A(I,KP,2) + A(I,K,3)
               X0P = T1R + T2R
               Y0P = T1I + T2I
               X1P = T3R + T4I
               Y1P = T3I - T4R
               X2P = T1R - T2R
               Y2P = T1I - T2I
               X3P = T3R - T4I
               Y3P = T3I + T4R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            B(I,0,R2) = A(I,R2,0) + A(I,R2,1)
            B(I,2,R2) = A(I,R2,3) - A(I,R2,2)
            T3 = A(I,R2,0) - A(I,R2,1)
            T4 = A(I,R2,3) + A(I,R2,2)
            B(I,1,R2) = ROOT2I*(T3+T4)
            B(I,3,R2) = -ROOT2I*(T3-T4)
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FQV(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix three Hermitian to real Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      real*8  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:R-1,0:2), B(0:P-1,0:2,0:R-1),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2)
C     .. Local Scalars ..
      real*8  C2K, CK, S2K, SK, T1, T1I, T1R, T2, T2I, T2R,
     *                  T3, T3I, T3R, X0P, X1P, X2P, Y0P, Y1P, Y2P
      INTEGER           I, K, KP
c                               , R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,1)
         T2 = A(I,0,0) - 0.5D0*T1
         T3 = SIN60*A(I,0,2)
         B(I,0,0) = A(I,0,0) + T1
         B(I,1,0) = T2 + T3
         B(I,2,0) = T2 - T3
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,1) + A(I,R-K,0)
               T1I = A(I,R-K,1) - A(I,K,2)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,R-K,2) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,1)-A(I,R-K,0))
               T3I = SIN60*(A(I,R-K,1)+A(I,K,2))
               X0P = A(I,K,0) + T1R
               Y0P = A(I,R-K,2) + T1I
               X1P = T2R + T3I
               Y1P = T2I - T3R
               X2P = T2R - T3I
               Y2P = T2I + T3R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = SINE(K,1)*X1P + COSINE(K,1)*Y1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = SINE(K,2)*X2P + COSINE(K,2)*Y2P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            DO 80 I = 0, P - 1
               T1R = A(I,K,1) + A(I,KP,0)
               T1I = A(I,KP,1) - A(I,K,2)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,KP,2) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,1)-A(I,KP,0))
               T3I = SIN60*(A(I,KP,1)+A(I,K,2))
               X0P = A(I,K,0) + T1R
               Y0P = A(I,KP,2) + T1I
               X1P = T2R + T3I
               Y1P = T2I - T3R
               X2P = T2R - T3I
               Y2P = T2I + T3R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = SK*X1P + CK*Y1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = S2K*X2P + C2K*Y2P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,R2,0) + A(I,R2,1)
C           T2 = 0.5D0*A(I,R2,0) - A(I,R2,1)
C           T3 = SIN60*A(I,R2,2)
C           B(I,0,R2) = T1
C           B(I,1,R2) = T2 + T3
C           B(I,2,R2) = -T2 + T3
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FQW(A,B,P,R,COSINE,SINE)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Radix two Hermitian to real Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      real*8  A(0:P-1,0:R-1,0:1), B(0:P-1,0:1,0:R-1),
     *                  COSINE(0:R-1), SINE(0:R-1)
C     .. Local Scalars ..
      real*8  CK, SK, X1HAT, Y1HAT
      INTEGER           I, K, KP
c                               ,R2      
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         B(I,0,0) = A(I,0,0) + A(I,0,1)
         B(I,1,0) = A(I,0,0) - A(I,0,1)
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1HAT = A(I,K,0) - A(I,R-K,0)
               Y1HAT = A(I,R-K,1) + A(I,K,1)
               B(I,0,K) = A(I,K,0) + A(I,R-K,0)
               B(I,0,R-K) = A(I,R-K,1) - A(I,K,1)
               B(I,1,K) = COSINE(K)*X1HAT - SINE(K)*Y1HAT
               B(I,1,R-K) = COSINE(K)*Y1HAT + SINE(K)*X1HAT
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K)
            SK = SINE(K)
            DO 80 I = 0, P - 1
               X1HAT = A(I,K,0) - A(I,KP,0)
               Y1HAT = A(I,KP,1) + A(I,K,1)
               B(I,0,K) = A(I,K,0) + A(I,KP,0)
               B(I,0,KP) = A(I,KP,1) - A(I,K,1)
               B(I,1,K) = CK*X1HAT - SK*Y1HAT
               B(I,1,KP) = CK*Y1HAT + SK*X1HAT
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           B(I,0,R2) = A(I,R2,0)
C           B(I,1,R2) = A(I,R2,1)
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FQX(A,B,M,N,Q,NQ,TRIG)
      
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C
C     Hermitian to Real Fast Fourier Transform Kernel Driver
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NQ
C     .. Array Arguments ..
      real*8  A(0:M*N-1), B(0:M*N-1), TRIG(0:2*N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      real*8  FACTOR
      INTEGER           I, P, QI, R
      LOGICAL           INA
C     .. External Subroutines ..
      EXTERNAL          C06FQQ, C06FQR, C06FQS, C06FQT, C06FQU, C06FQV,
     *                  C06FQW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = 1
      R = N
      IF (N.EQ.1) RETURN
      CALL C06FQQ(A,M,N)
      DO 20 I = 1, NQ
         QI = Q(I)
         R = R/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FQW(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FQV(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FQU(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FQT(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FQS(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FQR(A,B,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FQW(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FQV(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FQU(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FQT(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FQS(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FQR(B,A,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         END IF
         INA = .NOT. INA
         P = P*QI
   20 CONTINUE
C
      FACTOR = 2.0D0/SQRT(DBLE(N))
      IF (INA) THEN
         DO 40 I = 0, M*N - 1
            A(I) = A(I)*FACTOR
   40    CONTINUE
      ELSE
         DO 60 I = 0, M*N - 1
            A(I) = B(I)*FACTOR
   60    CONTINUE
      END IF
C
      RETURN
      END

