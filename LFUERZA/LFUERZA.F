C   ********************************************************************
C   * LFUERZA.  LECTURA Y REDUCCION TIME-HISTORY: THICOMY              *
C   *           THICOMZ TRANVA BALDAG RE5FOR TRANVA1   REV. 2.2        *
C   ********************************************************************
      COMMON/CALCUL/XMAX,XMIN,YMAX,YMIN,NRA,NTRA,NMXM
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      CHARACTER*4 IXX(8),IXY(8),IYX(8),ICLAV(8),IXXX(8),IXXY(8)
      CHARACTER*4 W,ICLA
      CHARACTER*8  FIC
      character*20 filee,files
      character*1 sn
      COMPLEX*16 FDW(35000)
      DIMENSION WORK(50000),X(20,50000),Y(20,50000),
     *Y1(50000),X1(50000),F1(35000),FY1(35000),XM(35000),YM(35000)
      DIMENSION IRA(888),ITRA(888,999),NTRAT(888),IC3(999),IYY(15)
      DATA IXX/ 'TIEM', 'PO E', 'N SE', 'GUND', 'OS  ', '    ', '    ',
     *'    '/
      DATA IXY/ 'FUER', 'ZA D', 'E ON', 'DA  ', 'IC1=', '    ',
     * '    ', '    '/
      DATA IYX/ 'F. D', 'E DE', 'SCAR', 'GA  ', 'IC1=', '    ',
     * '    ', '    '/
      DATA IXXX/ 'FREC', 'UENC', 'IA E', 'N HZ', '.   ', '    ', '    ',
     *'    '/
      DATA IXXY/ 'DENS', 'IDAD', ' ESP', 'ECT.', 'IC1=', '    ',
     * '    ', '    '/
      DATA W/ 'PLOT'/
      DATA ICLAV/ 'LTHY','LTHZ','LTRA','LBAL','LR5F','LTR1','LTHX',
     +            'LR5C'/
c     CALL ERRSET(215,0,256,2,0,0)
c     CALL ERRSET(226,256,-1,1,1)
c     CALL ERRSET(213,256,-1,1,1)
c     CALL ERRSET(208,256,-1,1,1)
      write(*,2540)
2540  format (' ¨Fichero de entrada ? ',\)
      read (*,*) filee
      write(*,3540)
3540  format (' ¨Fichero de salida ? ',\)
      read (*,*) files
132   write(*,4540)
4540  format (' ¨Quiere ficheros individuales? s/n ',\)
      read (*,'(a)') sn
      if (sn.ne.'s'.and.sn.ne.'S'.and.sn.ne.'n'.and.sn.ne.'N') goto 132
      OPEN(8,FILE='file8',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(12,FILE='file12',FORM='FORMATTED')
      OPEN(15,FILE=filee ,FORM='FORMATTED')
      OPEN(16,FILE=files ,FORM='FORMATTED')
      READ(15,23) IYY,ICLA
  23  FORMAT(18A4)
C---------------------------------------------------------------------C
C     LEE LOS DOS PRIMEROS REGISTROS PARA INICIALIZAR VALORES DE      C
C     DIMENSIONADO                                                    C
C---------------------------------------------------------------------C
      ICAN=0
      IF(ICLA.EQ.ICLAV(1))CALL LTHY(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(2))CALL LTHZ(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(3))CALL LTRA(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(4))CALL LBAL(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(5))CALL LR5F(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(6))CALL LTR1(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(7))CALL LTHX(ICAN,ICAN1)
      IF(ICLA.EQ.ICLAV(8))CALL LR5C(ICAN,ICAN1)
      IF(ICAN.NE.0)  GO TO 1
      WRITE(16,2)
   2  FORMAT(///' ******** ERROR EN LA DEFINICION DE ICLA *******')
      STOP
   1  WRITE(16,1010)
 1010 FORMAT(1H+,55X,'LFUERZA REVISION 2.2'/1H0)
      WRITE(16,1000) HEAD,IYY,NTR,NTRAM
 1000 FORMAT(1X,15A4/1X,15A4/' NO. TOTAL DE RAMAS=',I4/1X,
     *'NO. TOTAL DE TRAMOS=',I4)
C---------------------------------------------------------------------C
C     LEE INDICES DE CONTROL                                          C
C---------------------------------------------------------------------C
      READ(15,20)W1,W2,I1,IR,I2,NRAT,IESC,FMAX
  20  FORMAT(2F10.3,5I5,F10.3)
  21  FORMAT(/////,10X,'   TIEMPOS PUERTA    I.GRAF.   I.FFT.     NPUNTO
     *S    NRAMAS    I.ESC.     F.MAX'
     *//8X,2F10.3,2X,5(I6,4X),F10.3)
      WRITE(16,21)W1,W2,I1,IR,I2,NRAT,IESC,FMAX
      IF(FMAX.EQ.0.)FMAX=1.E+20
      IESC=IESC-1
      IF(IESC.LT.0)IESC=0
      IF(I2.GT.200)THEN
       NGRUPO=I2/200
       NGR1=NGRUPO*200
       IF(NGR1.NE.I2)NGRUPO=NGRUPO+1
       I2=NGRUPO*200-3*NGRUPO+3
      END IF
      NPRA=I2
      IF(W3.EQ.0.) W3=0.001
      NMXM=0
C---------------------------------------------------------------------C
C     LEE LA NUMERACION DE RAMAS Y TRAMOS REQUERIDOS                  C
C---------------------------------------------------------------------C
      IF(NRAT.EQ.0) GO TO 11
      IF(NRAT.GT.NTR) GO TO 3
      I3=0
      DO 12 I=1,NRAT
  14  READ(15,13) IRA(I),NTRAT(I),(ITRA(I,I3+II),II=1,14)
  13  FORMAT(16I5)
      I3=I3+14
      IF(NTRAT(I).LE.I3) GO TO 12
      GO TO 14
  12  I3=0
      GO TO 17
  11  NRAT=NTR
      DO 18 I=1,NTR
      IRA(I)=I
      NTRAT(I)=M2(I)-M1(I)+1
      I3=NTRAT(I)
      DO 18 II=1,I3
  18  ITRA(I,II)=II
C---------------------------------------------------------------------C
C     ESCRIBE LOS DATOS LEIDOS                                        C
C---------------------------------------------------------------------C
  17  I3=0
      WRITE(16,16)
  16  FORMAT(///,2X,'RAMA TTRA  ------ TRAMOS NO.------',/)
      DO 22 I=1,NRAT
  24  WRITE(16,13) IRA(I),NTRAT(I),(ITRA(I,I3+II),II=1,NTRAT(I))
      I3=I3+14
      IF(NTRAT(I).LE.I3) GO TO 22
      GO TO 24
  22  I3=0
C---------------------------------------------------------------------C
C     HALLA LOS INDICES DE LOS TRAMOS SEGUN NUMERACION TOTAL          C
C---------------------------------------------------------------------C
      DO 4001 I=1,NRAT
      IF(NTRAT(I).NE.0) GO TO 9
      NTRAT(I)=M2(IRA(I))-M1(IRA(I))+1
      LTRAT=NTRAT(I)
      DO 19 I4=1,LTRAT
  19  ITRA(I,I4)=I4
   9  LTRAT=NTRAT(I)
      DO 4001 II=1,LTRAT
      I3=I3+1
 4001 IC3(I3)=M1(IRA(I))+ITRA(I,II)-1
      GO TO 5
C---------------------------------------------------------------------C
C     HALLA LOS INDICES DE LOS TRAMOS DE DESCARGA                     C
C---------------------------------------------------------------------C
   3  I3=NTS
      DO 4 I=1,NTS
   4  IC3(I)=M2(NTR)+I
   5  IC2=(I3-1)/20+1
      IC4=I3-20*(IC2-1)
      IC5=0
C---------------------------------------------------------------------C
C     BUCLE EXTENDIDO  AL TOTAL DE LOS TRAMOS                         C
C---------------------------------------------------------------------C
      INP=0
      DO 4000 IC=1,IC2
      IF(IC4.EQ.0) GO TO 4000
C---------------------------------------------------------------------C
C     BUCLE DE LECTURA DE LAS FUERZAS PARA CADA TIEMPO                C
C     (MAXIMO DE 20 TRAMOS)                                           C
C---------------------------------------------------------------------C
      JJ=0
 2000 IF(ICAN1.NE.1)THEN
       IF(ICLA.NE.'LR5C')THEN
        READ(8,END=3000)(WORK(J),J=1,LEN)
C       WRITE(*,*)' &&1 ',(WORK(J),J=1,LEN)
       ELSE
        IF(LEN.LE.8)THEN
         READ(8,999,END=3000)(WORK(J),J=1,LEN)
C        WRITE(*,*)' &&2 ',(WORK(J),J=1,LEN)
999      FORMAT(3X,8E15.5)
        ELSE
         READ(8,999,END=3000)(WORK(J),J=1,8)
C        WRITE(*,*)' &&3 ',(WORK(J),J=1,8)
         READ(8,999,END=3000)(WORK(J),J=9,LEN)
C         WRITE(*,*)' &&4 ',(WORK(J),J=9,LEN)
        END IF
       END IF
       GO TO 2010
      ELSE
       READ(8,END=3000) FIC,(WORK(J),J=3,LEN)
C       WRITE(*,*)' &&5 ',FIC,(WORK(J),J=3,LEN)
      END IF
      IF( FIC(1:4).NE.W) GO TO 2000
 2010 IND=ICAN1+2
      IF(WORK(IND).LT.W1) GO TO 2000
      IF(WORK(IND).GT.W2) GO TO 3000
      JJ=JJ+1
      DO 30 L=1,IC4
      N=IC3(L+IC5)+ICAN
      X(L,JJ)=WORK(IND)-W1
  30  Y(L,JJ)=WORK(N)
      IF(JJ.EQ.1) GO TO 2000
      IF(X(1,JJ).LE.X(1,JJ-1)) JJ=JJ-1
      GO TO 2000
 3000 CONTINUE
C---------------------------------------------------------------------C
C     ESCRIBE TRAMO SELECCIONADO                                      C
C---------------------------------------------------------------------C
      DO 4002 L=1,IC4
      DEL0 = 0
      DEL1=DEL0
      YMI0 = 0
      YMI1=YMI0
      IF(NRAT.LE.NTR) GO TO 6
      NRA=0
      NTRA=IC3(L+IC5)
      GO TO 7
  6   DO 4003 IC11=1,NTR
      IF((M2(IC11).GE.IC3(L+IC5)).AND.(M1(IC11).LE.IC3(L+IC5))) NRA=IC11
 4003 CONTINUE
      NTRA=IC3(L+IC5)-M1(NRA)+1
C 7   WRITE(12,1101)HEAD,NTRA,NRA,IC3(L+IC5)
C1101 FORMAT(5X,15A4/5X,'TRAMO NO.',I4,' DE LA RAMA NO.',I4/5X,
C    *'(IC1=',I4,')')
   7  DO 31 L1=1,JJ
      Y1(L1)=Y(L,L1)
  31  X1(L1)=X(L,L1)
C---------------------------------------------------------------------C
C     HALLA LOS VALORES MAXIMOS Y MINIMOS                             C
C---------------------------------------------------------------------C
      YMAX=-1.E+20
      YMIN=1.E+20
      DO 50 I=1,JJ
      YMAX=AMAX1(YMAX,Y(L,I))
   50 YMIN=AMIN1(YMIN,Y(L,I))
      YMAX=YMAX+0.01*ABS(YMAX)
      YMIN=YMIN-0.01*ABS(YMIN)
C       WRITE(*,*) JJ
      XMAX=X1(JJ)
      XMIN=X1(1)
C---------------------------------------------------------------------C
C     DIBUJA NUBE DE PUNTOS POR IMPRESORA                             C
C---------------------------------------------------------------------C
      WRITE(16,2022)
 2022 FORMAT(1H1,55X,'LFUERZA REVISION 2.2'/1H0)
      WRITE(16,2222) HEAD,NRA,NTRA,JJ,IC3(L+IC5)
 2222 FORMAT(1H ,1X,15A4/'  RAMA',I3,', TRAMO',I3,' CON  ',I4,
     *' PARES DE VALORES'
     */'  (IC1=',I4,')',
     */'  FUNCION TEMPORAL ORIGINAL')
      IF(I1.EQ.1)THEN
       IF(NRAT.LE.NTR)THEN
        CALL GRAF(X1,Y1,JJ,IXY,IXX)
       ELSE
        CALL GRAF(X1,Y1,JJ,IYX,IXX)
       END IF
      END IF
C---------------------------------------------------------------------C
C     HALLA DENSIDAD ESPECTRAL DE POTENCIA                            C
C---------------------------------------------------------------------C
      IF(IR.NE.0)THEN
       PUN=JJ
       L2N=INT(1.4427*ALOG(PUN))+1
       MPUN=2**L2N
       IF(JJ*2.LE.MPUN)THEN
        L2N=L2N-1
        MPUN=MPUN/2
       END IF
       INDICE=0
       CALL INTERP(X1,Y1,XM,YM,JJ,MPUN,INDICE)
       DW=1./XM(2)/MPUN
       F1(1)=0.
       KPUN=MPUN/2+1
       DO 574 IPUN=2,KPUN
        F1(IPUN)=F1(IPUN-1)+DW
        IF(F1(IPUN).GT.FMAX)THEN
         LPUN=IPUN-1
         GO TO 314
        END IF
574    CONTINUE
       LPUN=KPUN
314    DO 584 IPUN=1,MPUN
        FDW(IPUN)=CMPLX(YM(IPUN),0.)
584    CONTINUE
       CALL FFT(FDW,MPUN,L2N)
       DO 432 IPUN=1,LPUN
        FY1(IPUN)=REAL(FDW(IPUN))**2+DIMAG(FDW(IPUN))**2
432    CONTINUE
       XMAX1=XMAX
       YMAX1=YMAX
       XMIN1=XMIN
       YMIN1=YMIN
       XMIN=0.
       XMAX=F1(LPUN)
       YMAX=-1.E+20
       YMIN=1.E+20
       DO 675 IPUN=1,LPUN
        YMAX=AMAX1(YMAX,FY1(IPUN))
        YMIN=AMIN1(YMIN,FY1(IPUN))
675    CONTINUE
       WRITE(16,2032)
2032   FORMAT(1H1,55X,'LFUERZA REVISION 2.2'/1H0)
       WRITE(16,2332) HEAD,NRA,NTRA,JJ,IC3(L+IC5)
 2332  FORMAT(1H ,1X,15A4/'  RAMA',I3,', TRAMO',I3,' CON  ',I4,
     *' PARES DE VALORES'
     */'  (IC1=',I4,')',
     */'  DENSIDAD ESPECTRAL DE POTENCIA DE LA FUNCION TEMPORAL ORIGINAL
     *')
       CALL GRAF(F1,FY1,LPUN,IXXY,IXXX)
       XMAX=XMAX1
       YMAX=YMAX1
       XMIN=XMIN1
       YMIN=YMIN1
      END IF
C---------------------------------------------------------------------C
C     REDUCE NUMERO DE PUNTOS Y DIBUJA GRAFICA DE CURVA REDUCIDA      C
C---------------------------------------------------------------------C
      CALL REDUCE(X1,Y1,JJ,I2,ER1)
      IF(I1.NE.1) GO TO 33
      Y1(NPRA+1)=YMI1
      Y1(NPRA+2)=DEL1
      X1(NPRA+1)=X(L,1)
      DX = 0
      X1(NPRA+2)=DX
      IF(I1.EQ.1)THEN
       WRITE(16,2042)
2042   FORMAT(1H1,55X,'LFUERZA REVISION 2.2'/1H0)
       WRITE(16,2442) HEAD,NRA,NTRA,JJ,IC3(L+IC5)
 2442  FORMAT(1H0,1X,15A4/'  RAMA',I3,', TRAMO',I3,' CON  ',I4,
     *' PARES DE VALORES'
     */'  (IC1=',I4,')',
     */'  FUNCION TEMPORAL REDUCIDA')
       IF(NRAT.LE.NTR)THEN
        CALL GRAF(X1,Y1,I2,IXY,IXX)
       ELSE
        CALL GRAF(X1,Y1,I2,IYX,IXX)
       END IF
      END IF
C---------------------------------------------------------------------C
C     HALLA DENSIDAD ESPECTRAL DE POTENCIA DE LA CURVA REDUCIDA       C
C---------------------------------------------------------------------C
      IF(IR.NE.0)THEN
       INDICE=1
       CALL INTERP(X1,Y1,XM,YM,I2,MPUN,INDICE)
       DO 1584 IPUN=1,MPUN
        FDW(IPUN)=CMPLX(YM(IPUN),0.)
1584   CONTINUE
       CALL FFT(FDW,MPUN,L2N)
       DO 1432 IPUN=1,LPUN
        FY1(IPUN)=REAL(FDW(IPUN))**2+DIMAG(FDW(IPUN))**2
1432   CONTINUE
       XMAX1=XMAX
       YMAX1=YMAX
       XMIN1=XMIN
       YMIN1=YMIN
       XMIN=0.
       XMAX=F1(LPUN)
       YMAX=-1.E+20
       YMIN=1.E+20
       DO 673 IPUN=1,LPUN
        YMAX=AMAX1(YMAX,FY1(IPUN))
        YMIN=AMIN1(YMIN,FY1(IPUN))
673    CONTINUE
       WRITE(16,2052)
2052   FORMAT(1H1,55X,'LFUERZA REVISION 2.2'/1H0)
       WRITE(16,2552) HEAD,NRA,NTRA,JJ,IC3(L+IC5)
 2552  FORMAT(1H ,1X,15A4/'  RAMA',I3,', TRAMO',I3,' CON  ',I4,
     *' PARES DE VALORES'
     */'  (IC1=',I4,')',
     */'  DENSIDAD ESPECTRAL DE POTENCIA DE LA FUNCION TEMPORAL REDUCIDA
     *')
       CALL GRAF(F1,FY1,LPUN,IXXY,IXXX)
       XMAX=XMAX1
       YMAX=YMAX1
       XMIN=XMIN1
       YMIN=YMIN1
      END IF
C---------------------------------------------------------------------C
C     ESCRIBE VALORES TIEMPO-FUERZA                                   C
C---------------------------------------------------------------------C
      WRITE(16,2023)
2023  FORMAT(1H1,55X,'LFUERZA REVISION 2.2'/1H0)
      WRITE(16,2223) HEAD,NRA,NTRA,JJ
 2223 FORMAT(1H ,15A4/' DATOS DE ENTRADA PARA AGPIPE'//17X,'  RAMA',I5,
     *' ,TRAMO',I5,' CON',I6,' PUNTOS ORIGINALES'/)
  33  IF(I2.LE.200)THEN
       KK=NPRA/3
       K1=3*KK
       K2=NPRA-K1
       LLL=L+IC5
       NLIN=KK
      ELSE
       K2=2
       NLIN=66*NGRUPO
       TT=0.
       XX=0.
      END IF
      ILIN=0
      IESC=IESC+1
      LL=3
      LLL=0
      DO 510 K=1,NLIN
      IF(LL.GT.200)THEN
       IF(X1(LLL+1).EQ.0.)GO TO 4005
       WRITE(12,1102) IESC,(X1(LLL+I),Y1(LLL+I),I=1,2)
       WRITE(16,1103) IESC,(X1(LLL+I),Y1(LLL+I),I=1,2)
       IF(X1(LLL+2).EQ.0.)GO TO 4005
       IESC=IESC+1
       TT1=X1(LLL+2)+.00001
       WRITE(12,1102) IESC,TT,XX,X1(LLL+2),XX,TT1,Y1(LLL+2)
       WRITE(16,1103) IESC,TT,XX,X1(LLL+2),XX,TT1,Y1(LLL+2)
       LL=3
       LLL=LLL+2
      ELSE
       IF(X1(LLL+1).EQ.0..AND.X1(LLL+2).EQ.0.)GO TO 4005
       WRITE(12,1102) IESC,(X1(LLL+I),Y1(LLL+I),I=1,3)
       WRITE(16,1103) IESC,(X1(LLL+I),Y1(LLL+I),I=1,3)
       LL=LL+3
       LLL=LLL+3
      END IF
  510 CONTINUE
 1102 FORMAT('TIME FORCE',5X,I5,3(F10.5,F10.0))
 1103 FORMAT(1X,'TIME FORCE',5X,I5,3(F10.5,F10.0))
      IF(K2.EQ.0.OR.X1(LLL+1).EQ.0.) GO TO 4005
      WRITE(12,1102) IESC,(X1(LLL+I),Y1(LLL+I),I=1,K2)
      WRITE(16,1103) IESC,(X1(LLL+I),Y1(LLL+I),I=1,K2)
 4005 WRITE(16,6007) NPRA,ER1
      IF(NMXM.GT.0) WRITE(16,6009) NMXM
 6009 FORMAT(' NO. DE MAX. Y MIN. SIGNIFICATIVOS:',I5)
 6007 FORMAT(1H0//,' NO. DE PUNTOS:',I3/' EL ERROR EN LA REDUCCION ES ME
     *NOR DE:',E10.3)
 4002 CONTINUE
      IF(IC2.EQ.IC) GO TO 4004
      REWIND 8
      IF(ICAN.NE.1) THEN
        READ(8) FIC
C       WRITE(*,*)' &&6 ',FIC
      END IF
      IF(ICLA.EQ.'LTHZ')THEN
       READ(8)FIC
C       WRITE(*,*)' &&7 ',FIC
       READ(8)FIC
C       WRITE(*,*)' &&8 ',FIC
      END IF
      READ(8) FIC
C       WRITE(*,*)' &&9 ',FIC
      IC5=IC5+IC4
 4000 IC4=20
 4004 END FILE 8
      IF (sn.eq.'s'.or.sn.eq.'S')CALL TRANS
      STOP
      END
      SUBROUTINE GRAF(XMAT,YMAT,NPUN,IXY,IXX)
C   ******************************************************************
C   *   GRAF. REPRESENTACION GRAFICA POR IMPRESORA DE UNA NUBE
C   *   DE PUNTOS.     EPTISA-GHESA-TRSA (EMPRESARIOS AGRUPADOS)
C   ******************************************************************
      COMMON/CALCUL/XMAX,XMIN,YMAX,YMIN,NRA,NTRA,NMXM
      CHARACTER*2 A(101,51),X,PTO,C
      CHARACTER*4 IXY,IXX
      DIMENSION XMAT(1),YMAT(1),IXY(8),XS(11),YS(6),IXX(8)
      IF(NPUN.LT.2)RETURN
      DATA C/     '  '/,X/'X '/,PTO/'. '/
      DO 5 I=1,101
      DO 5 J=1,51
    5 A(I,J)=C
      DO 10 I=1,101
      A(I,1)=PTO
   10 A(I,51)=PTO
      DO 20 I=1,51
      A(1,I)=PTO
   20 A(101,I)=PTO
      DO 30 I=11,91,10
      A(I,2)=PTO
   30 A(I,50)=PTO
      DO 40 I=11,41,10
      A(2,I)=PTO
   40 A(100,I)=PTO
      DELTX=(XMAX-XMIN)/100.
      DELTY=(YMAX-YMIN)/50.
      IF(DELTX.LT.1.E-15.OR.DELTY.LT.1.E-15)RETURN
      DO 60 K=1,NPUN
      I=(XMAT(K)-XMIN)/DELTX+1.5
      J=52-(YMAT(K)-YMIN)/DELTY -0.5
   60 A(I,J)=X
      DO 70 I=1,11
   70 XS(I)=XMIN+(XMAX-XMIN)*(I-1)/10.
      DO 80 I=1,6
   80 YS(I)=YMAX-(YMAX-YMIN)*(I-1)/5.
      WRITE(16,1)
    1 FORMAT(1H0//)
      WRITE(16,6001)YS(1),(A(I,1),I=1,101)
 6001 FORMAT(1X,1PE20.3,101A1)
      DO 90 J=2,10
   90 WRITE(16,6002)(A(I,J),I=1,101)
 6002 FORMAT(21X,101A1)
      WRITE(16,6001)YS(2),(A(I,11),I=1,101)
      DO 100 J=12,20
  100 WRITE(16,6002)(A(I,J),I=1,101)
      WRITE(16,6001)YS(3),(A(I,21),I=1,101)
      DO 110 J=22,24
  110 WRITE(16,6002)(A(I,J),I=1,101)
      WRITE(16,6003)(IXY(I),I=1,4),(A(I,25),I=1,101)
      WRITE(16,6004)(A(I,26),I=1,101)
      WRITE(16,6004)(A(I,27),I=1,101)
 6003 FORMAT(1X,4A4,4X,101A1)
 6004 FORMAT(21X,101A1)
      DO 120 J=28,30
  120 WRITE(16,6002)(A(I,J),I=1,101)
      WRITE(16,6001)YS(4),(A(I,31),I=1,101)
      DO 130 J=32,40
  130 WRITE(16,6002)(A(I,J),I=1,101)
      WRITE(16,6001)YS(5),(A(I,41),I=1,101)
      DO 140 J=42,50
  140 WRITE(16,6002)(A(I,J),I=1,101)
      WRITE(16,6001)YS(6),(A(I,51),I=1,101)
      WRITE(16,6005)(XS(I),I=1,11)
      WRITE(16,6006)(IXX(I),I=1,8)
 6005 FORMAT(16X,11(1PE10.3))
 6006 FORMAT(1H0,50X,8A4)
      RETURN
      END
      SUBROUTINE LTHZ(ICAN,ICAN1)
      CHARACTER*8  ALFIN(50000)
      DIMENSION NUMBIN(50000)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      READ(8)(HEAD(I),I=1,12)
C     WRITE(*,*)' &&10 ',(HEAD(I),I=1,12)
      READ(8)IFIC,IFIC,IFIC,LEN
C     WRITE(*,*)' &&11 ',IFIC,IFIC,IFIC,LEN
      IF(LEN.GT.10001)THEN
       WRITE(16,'(1H ,''LONGITUD DE REGISTRO NO PERMITIDA'')')
       STOP
      END IF
      NVAR=LEN-1
      LEN=LEN+1
      READ(8)IFIC,IFIC,( ALFIN(I),I=1,NVAR)
C     WRITE(*,*)' &&12 ',IFIC,IFIC,( ALFIN(I),I=1,NVAR)
      READ(8)IFIC,IFIC,(IFIC,NUMBIN(I),I=1,NVAR)
C     WRITE(*,*)' &&13 ',IFIC,IFIC,(IFIC,NUMBIN(I),I=1,NVAR)
      READ(8)IFIC,IFIC,NTR
C     WRITE(*,*)' &&14 ',IFIC,IFIC,NTR
      DO 1 I=1,NVAR
       IF( ALFIN(I).EQ.'F. ONDA '.OR. ALFIN(I).EQ.'F. O    ')GO TO 2
1     CONTINUE
      WRITE(16,'(1H ,''ERROR'')')
      STOP
2     ICAN=I+1
      ICAN1=1
      NTRAM=0
      M1(1)=1
      M2(1)=0
      ITR=1
      II=ICAN-2
3     II=II+1
      IF(II.GT.NVAR.OR. ALFIN(II)(:4).NE.'F. O')RETURN
      NTRAM=NTRAM+1
      J=NUMBIN(II)/1000
      IF(J.EQ.ITR)THEN
       M2(ITR)=M2(ITR)+1
      ELSE
       ITR=ITR+1
       M1(ITR)=M2(ITR-1)+1
       M2(ITR)=M2(ITR-1)+1
      END IF
      GO TO 3
      END
      SUBROUTINE LTHY(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      NTRAM=0
      READ(8)HEAD
C     WRITE(*,*)' &&15 ',HEAD
      READ(8)NTR,NTC,NTB,NTST,NTV,TMAX,K1,(N1(I),I=1,NTR)
     *,(N2(I),I=1,NTR)
C     WRITE(*,*)' &&16 ',NTR,NTC,NTB,NTST,NTV,TMAX,K1,(N1(I),I=1,NTR)
C    *,(N2(I),I=1,NTR)
      IF(NTB.EQ.0) NTB=1
      IF(NTC.EQ.0) NTC=1
      IF(NTV.EQ.0) NTV=1
      IF(NTST.EQ.0) NTST=1
      DO 10 I=1,NTR
      NTRAM=NTRAM+N2(I)-1
      M2(I)=NTRAM
 10   M1(I)=M2(I)-N2(I)+2
      LEN=3+5*NTC+3*(NTB+NTST)+2*NTR+NTRAM+2*K1+NTV
      ICAN=LEN-K1-NTV-NTRAM
      ICAN1=1
      RETURN
      END
C***********************************************************************
      SUBROUTINE LTHX(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      READ(8)HEAD
C     WRITE(*,*)' &&17 ',HEAD
      READ(8)NTR,NTB,(N1(I),I=1,NTR),(N2(I),I=1,NTR)
C     WRITE(*,*)' &&18 ',NTR,NTB,(N1(I),I=1,NTR),(N2(I),I=1,NTR)
      IF(NTB.EQ.0) NTB=1
      NNJ=0
      NTRAM=0
      DO 1 I=1,NTR
       NNJ=NNJ+N1(I)
       NTRAM=NTRAM+N2(I)-1
       M2(I)=NTRAM
       M1(I)=NTRAM-N2(I)+2
1     CONTINUE
      ICAN=3+2*NTR+2*NTB+NNJ
      LEN=ICAN+NTRAM
      ICAN1=1
      RETURN
      END
C***********************************************************************
      SUBROUTINE LTRA(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      READ(8)HEAD
C     WRITE(*,*)' &&19 ',HEAD
      READ(8)NTR,NNOD,(N1(I),I=1,NTR),(N2(I),I=1,NTR),(M1(I),I=1,NTR),
     *(M2(I),I=1,NTR)
C     WRITE(*,*)' &&20 ',NTR,NNOD,(N1(I),I=1,NTR),(N2(I),I=1,NTR),
C    *(M1(I),I=1,NTR),
C    *(M2(I),I=1,NTR)
      NTRAM=M2(NTR)
      LEN=NTRAM+4*NNOD+2
      ICAN=2+4*NNOD
      ICAN1=0
      RETURN
      END
      SUBROUTINE LTR1(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      READ(8)HEAD
C     WRITE(*,*)' &&21 ',HEAD
      READ(8)NTR,NNOD,NTS,(N1(I),I=1,NTR),(N2(I),I=1,NTR),
     *(M1(I),I=1,NTR),(M2(I),I=1,NTR)
C     WRITE(*,*)' &&22 ',NTR,NNOD,NTS,(N1(I),I=1,NTR),(N2(I),I=1,NTR),
C    *(M1(I),I=1,NTR),(M2(I),I=1,NTR)
      NTRAM=M2(NTR)
      LEN=NTRAM+4*NNOD+2+NTS
      ICAN=2+4*NNOD
      ICAN1=0
      RETURN
      END
      SUBROUTINE LBAL(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      READ(8)HEAD
C     WRITE(*,*)' &&23 ',HEAD
      READ(8)NTR,NNOD,(N1(I),I=1,NTR),(N2(I),I=1,NTR),(M1(I),I=1,NTR),
     *(M2(I),I=1,NTR)
C     WRITE(*,*)' &&24 ',NTR,NNOD,(N1(I),I=1,NTR),(N2(I),I=1,NTR),
C    *(M1(I),I=1,NTR),
C    *(M2(I),I=1,NTR)
      NTRAM=M2(NTR)
      NTS=1
      LEN=NTRAM+4*NNOD+9
      ICAN=8+4*NNOD
      ICAN1=-1
      RETURN
      END
      SUBROUTINE LR5F(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      DIMENSION IIC(5)
      READ(8) (HEAD(I),I=1,15),IIC,NTRAM
C     WRITE(*,*)' &&25 ',(HEAD(I),I=1,15),IIC,NTRAM
      DO 10 I=1,NTRAM
      M1(I)=I
  10  M2(I)=I
      LEN=NTRAM+1
      NTR=NTRAM
      ICAN=1
      ICAN1=-1
      RETURN
      END
      SUBROUTINE LR5C(ICAN,ICAN1)
      COMMON/VARIAB/NTRAM,NTS,NTR,LEN,M1(999),M2(999),HEAD(15),N1(9000)
     *,N2(9000)
      READ(15,101)NTRAM
101   FORMAT(I5)
      DO 10 I=1,NTRAM
      M1(I)=I
  10  M2(I)=I
      LEN=NTRAM+1
      NTR=NTRAM
      ICAN=1
      ICAN1=-1
      RETURN
      END
*********************************************************************
*     SUBRUTINA DE INTERPOLACION
*********************************************************************
      SUBROUTINE INTERP(XA,YA,XM,YM,IA,IM,INDICE)
      DIMENSION XA(1),YA(1),XM(1),YM(1)
      IF(INDICE.EQ.0)THEN
       XM(1)=XA(1)
       XM(IM)=XA(IA)
       DX=(XM(IM)-XM(1))/(IM-1)
       DO 1 IPUN=2,IM-1
        XM(IPUN)=XM(IPUN-1)+DX
1      CONTINUE
      END IF
      YM(1)=YA(1)
      YM(IM)=YA(IA)
      IANT=2
      DO 4 IPUN=2,IM-1
       DO 2 JPUN=IANT,IA
        IF(XM(IPUN).LE.XA(JPUN))GO TO 3
2      CONTINUE
3      IANT=JPUN
       YM(IPUN)=YA(IANT-1)+(YA(IANT)-YA(IANT-1))/(XA(IANT)-XA(IANT-1))*
     +          (XM(IPUN)-XA(IANT-1))
4     CONTINUE
      RETURN
      END
****************************************************
*     CALCULO DE LA TRANSFORMADA DE FOURIER
****************************************************
      SUBROUTINE FFT(FDW,MPUN,L2N)
      COMPLEX*16 FDW,U,W,T
      REAL*8 PI,CO,SI
      DIMENSION FDW(35000)
      PI=3.141592653589793D0
      DO 1 J=1,MPUN
       FDW(J)=FDW(J)/MPUN
1     CONTINUE
      NBD2=MPUN/2
      NBD1=MPUN-1
      J=1
      DO 5 L=1,NBD1
       IF(L.GE.J)GO TO 2
       T=FDW(J)
       FDW(J)=FDW(L)
       FDW(L)=T
2      K=NBD2
3      IF(K.GE.J)GO TO 4
       J=J-K
       K=K/2
       GO TO 3
4      J=J+K
5     CONTINUE
      DO 8 M=1,L2N
       U=(1.0,0.0)
       ME=2**M
       K=ME/2
       CO=DCOS(PI/K)
       SI=-DSIN(PI/K)
       W=DCMPLX(CO,SI)
       DO 7 J=1,K
        DO 6 L=J,MPUN,ME
         LPK=L+K
         T=FDW(LPK)*U
         FDW(LPK)=FDW(L)-T
         FDW(L)=FDW(L)+T
6       CONTINUE
        U=U*W
7      CONTINUE
8     CONTINUE
      RETURN
      END
*******************************************************
*     SUBRUTINAS DE REDUCCION
*******************************************************
      SUBROUTINE REDUCE(TIMETH,ACCTH,NPUNTH,NNN,AMENOR)
      DIMENSION TIMETH(1),ACCTH(1)
      DIMENSION T(50000),AC(50000),TIME(50000),ACC(50000)
      AMENOR = 0.
      IF(NNN.GE.NPUNTH)GO TO 1000
      NPUN = NPUNTH - 1
      TIME(1)  = TIMETH(1)
      TIME(2)  = TIMETH(NPUNTH)
      ACC(1) = ACCTH(1)
      ACC  (2)= ACCTH(NPUNTH)
      N=2
      DO 600 I=3,NNN
      DO 20 J=1,N
      T(J)=TIME(J)
      AC(J)=ACC(J)
   20 CONTINUE
      II=N
      CALL INTER(TIMETH,NPUNTH,T,AC,II )
      AM=0.
      DO 100 K=2,NPUN
      A= ABS(ACCTH(K)-AC(K))
      IF(A.LE.AM)GO TO 100
      AM=A
      IM=K
  100 CONTINUE
	IM = MAX0(IM,1)
      DO 200 K=1,N
      IF(TIMETH(IM).LT.TIME(K))GO TO 225
  200 CONTINUE
  225 A1=TIME(K)
      A2=ACC(K)
      DO 300 J=K,N
      J1=J+1
      B1=TIME(J1)
      B2=ACC(J1)
      TIME(J1)=A1
      ACC(J1)=A2
      A1=B1
      A2=B2
  300 CONTINUE
      TIME(K)=TIMETH(IM)
      ACC(K)=ACCTH(IM)
      N=N+1
      AMENOR=AM
      IF(AM.LE.0.01)GO TO 625
  600 CONTINUE
  625 CONTINUE
      DO 855 K=1,N
      TIMETH(K) = TIME(K)
  855 ACCTH(K) = ACC(K)
      IF(N.GE.NNN)GO TO 999
      N1 = N+1
      DO 900 I = N1,NNN
      TIMETH(I) = 0.
  900 ACCTH(I) = 0.
  999 CONTINUE
 1000 RETURN
      END
      SUBROUTINE INTER(X1,NT1,X2,Y2,NT2)
      DIMENSION X1(1),X2(1),Y1(50000),Y2(1)
      IF(NT2.EQ.0)GO TO 250
      I1 =1
      DO 200 I=1,NT1
      DO 100 K=I1,NT2
      IF(X2(K).EQ.X1(I))GO TO 190
      IF(X1(I).GT.X2(K))GO TO 100
      GO TO 125
  100 CONTINUE
      K   = NT2
      GO TO 190
  125 K1 = K-1
      IF (K1.EQ.0)GO TO 190
      K2 = K1 +1
  180 Y1(I)= Y2(K1) + (X1(I)-X2(K1))*(Y2(K2)-Y2(K1))/(X2(K2)-X2(K1))
      GO TO 195
  190  Y1(I) = Y2(K)
  195 I1 = K
  200 CONTINUE
  250 DO 300 I=1,NT1
      Y2(I)=Y1(I)
      IF(NT2.EQ.0)Y2(I)= 0.
      X2(I)=X1(I)
  300 CONTINUE
      NT2=NT1
      RETURN
      END
C
      SUBROUTINE TRANS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
c
c  a partir de un file12 de salida de lfuerza crea ficheros de cada
c  rama del modelo (hasta 100 ramas) individualizando el fichero file12
c  en tantos ficheros como ramas del modelo
c
      CHARACTER*82 A
      CHARACTER*3 FIC(100)
      CHARACTER*10 FICT,TTT,UUU
      DIMENSION T(300),F(300)
      close(12)
C     WRITE(*,*) ' FICHERO TOTAL '
      FICT = 'file12'
c     OPEN(2,FILE='pp',FORM='FORMATTED')
      OPEN(12,FILE=FICT,FORM='FORMATTED')
      I0 = 0
1     READ(12,'(A82)',END = 900) A
      READ(A(19:20),'(I2.2)') I1
      IF(I1.NE.I0) THEN
         I0 = I1
        CLOSE(2)
        FIC(I1)(1:1) = 'T'
        WRITE(FIC(I1)(2:3),'(I2.2)') I1
        OPEN(2,FILE=FIC(I1),FORM='FORMATTED')
      ENDIF
      WRITE(2,'(A82)') A
      GOTO 1
      CLOSE(12)
 900  CONTINUE
      CLOSE(2)
      DO 50 I=1,I1
       OPEN(1,FILE=FIC(I),FORM='FORMATTED')
       II=1
   3   CONTINUE
       READ(1,'(20X,6(F10.0))',END=999)T(II),F(II),T(II+1),F(II+1),
     * T(II+2),F(II+2)
C      PRINT *,T(II),F(II),II
       II= II+3
       GOTO 3
 999   CONTINUE
       CLOSE(1)
       OPEN(1,FILE=FIC(I) ,FORM='FORMATTED')
        WRITE(1,'(A3)') FIC(I)
       DO 30 J=1,II
       IF (T(J).EQ.0..AND.F(J).EQ.0..AND.J.NE.1) GOTO 30
       WRITE(TTT,'(F10.5)') T(J)
       WRITE(UUU,'(F10.0)') F(J)
C      TE=' '
       LF = 11
 70    LF=LF-1
       IF(TTT(LF:LF).EQ.' ')GOTO 70
       LP=0
 75    LP = LP+1
       IF(TTT(LP:LP).EQ.' ')GOTO 75
       MF = 11
 60    MF=MF-1
       IF (UUU(MF:MF).EQ.' ')GOTO 60
       MP=0
 65    MP = MP+1
       IF(UUU(MP:MP).EQ.' ')GOTO 65
       WRITE(1,'(A,1X,A)') TTT(LP:LF),UUU(MP:MF)
  30   CONTINUE
       WRITE(1,'(A)') ' '
C      PRINT *, FIC(I)
       CLOSE(1)
  50   CONTINUE
       RETURN
       END
