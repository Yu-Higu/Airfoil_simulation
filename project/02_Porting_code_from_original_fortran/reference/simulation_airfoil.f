C=======================================================================
C/                                              2010/10/20  VER 0.20   /
C/                                                                     /
C/         SIMULATION OF TWO DIMENSIONAL FLOW PAST AN AIRFOIL          /
C/                    BY FINITE DIFFERENCE METHOD                      /
C/                                                                     /
C/    <CFD COURSE OF EXPERIMENTS IN MECHANICAL ENGINEERING FOR 3ME>    /
C/                                                                     /
C=======================================================================
C/                                                                     /
C/   PROGRAMER : TAKAHIRO TSUKAHARA (TOKYO UNVERSITY OF SCIENCE)       /
C/                                                                     /
C/   SMAC METHOD ON STAGGERED GRID (NON-UNIFORM MESHES IN X,Y DIR.)    /
C/      2ND-ORDER CENTRAL SCHEME                                       /
C/      2ND-ORDER ADAMS-BASHFORTH METHOD FOR TIME ADVANCEMENT          /
C/      SOR METHOD FOR POISSON SOLVER                                  /
C/                                                                     /
C=======================================================================
C                                                                       
      PROGRAM MAIN                                                      
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
      CHARACTER*14 NWING                                                
C                                                                       
C=================================================================      
C                                                                       
C------------------------------------ AIRFOIL DATA                      
C                                                                       
      OPEN(10, FILE='clark_y.dat', STATUS="OLD")                      
C                                                                       
C------------------------------------ SET PARAMETERS                    
C                                                                       
      WRITE(*,*)"REYNOLDS NUMBER, RE ="                                 
      READ(*,*) RE                                                      
      WRITE(*,*)"ANGLE OF ATTACK, Degree ="                             
      READ(*,*) DEG                                                     
C                                                                       
C----------------- TIME STEP                                            
C                                                                       
      DT=0.002500D0                                                      
C                                                                       
C----------------- NUMBER OF STEPS                                      
C                                                                       
      NC=10000                                                            
      NG=5                                                              
C                                                                       
C----------------- SIZE OF COMPUTATIONAL DOMAIN                         
C                                                                       
      XL=8.0D0                                                          
      YL=3.0D0                                                          
C                                                                       
C------------------------------------ INITIALIZATION                    
C                                                                       
      NN=0                                                              
C------------------------------------ VEROCITY SEED TIME STEP COUNTER   
      NUM=0                                                             
      TIME=0.0D0                                                        
C------------------------------------ AVERAGE DATA STEP COUNTER         
      NAVE=0                                                            
      AVETIM=0.0D0                                                      
C                                                                       
      IERR=0                                                            
C                                                                       
C------------------------------------ INPUT SEED                        
C      CALL SEEDIN(IERR)                                                 
      CALL CALBETA(IERR,NWING)                                          
      IF(IERR.NE.0) GO TO 1                                             
C                                                                       
C------------------------------------                                   
C                                                                       
      CALL UBC(DIV,CFL,USUM,CL,CD)                                      
C                                                                       
C-- FIRST STEP                                                          
C                                                                       
      NN=NN+1                                                           
      NUM=NUM+1                                                         
      TIME=TIME+DT                                                      
      AVETIM=AVETIM+DT                                                  
C                                                                       
      IF(NUM.EQ.1) CALL FS0                                             
      IF(NUM.NE.1) CALL FS1                                             
      CALL FS2(NCOUNT)                                                  
      CALL UBC(DIV,CFL,USUM,CL,CD)                                      
C                                                                       
      OPEN(1, FILE='TEXTDATA/time-series.dat', STATUS="UNKNOWN")        
      WRITE(1,2000)                                                     
      WRITE(6,2000)                                                     
      WRITE(1,1000) NN,NUM,TIME,DIV,CFL,USUM,CL,CD,NCOUNT               
      WRITE(6,1000) NN,NUM,TIME,DIV,CFL,USUM,CL,CD,NCOUNT               
C                                                                       
C-- SECOND STEP                                                         
C                                                                       
      DO NST=2,NC                                                       
       NN=NN+1                                                          
       NUM=NUM+1                                                        
       TIME=TIME+DT                                                     
       AVETIM=AVETIM+DT                                                 
C                                                                       
       CALL FS1                                                         
       CALL FS2(NCOUNT)                                                 
       CALL UBC(DIV,CFL,USUM,CL,CD)                                     
       IF(MOD(NN,NG).EQ.0) THEN                                         
        WRITE(1,1000) NN,NUM,TIME,DIV,CFL,USUM,CL,CD,NCOUNT             
        WRITE(6,1000) NN,NUM,TIME,DIV,CFL,USUM,CL,CD,NCOUNT             
       ENDIF                                                            
      ENDDO                                                             
C                                                                       
C------------------------------------ OUTPUT SEED                       
      CALL SEEDOUT                                                      
C                                                                       
      CALL OUTPUT(NWING)                                                
C                                                                       
 1000 FORMAT(I6,I8,F10.4,3E12.4,2F9.4,I8)                               
 2000 FORMAT('    NN     NUM    TIME    DIV         CFL         AVE.U   
     &     C_L      C_D     NCOUNT')                                    
C                                                                       
      CLOSE(1)                                                          
C                                                                       
    1 CONTINUE                                                          
      STOP                                                              
      END                                                               
C                                                                       
C====================================================================   
C     CALC. TEMPORARY VELOCIYT AT FIRST STEP (VELOCITY FIELD)           
C     EULER EXPLICIT METHOD                                             
      SUBROUTINE FS0                                                    
C     2010.10.13                                                        
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
C=================================================================      
C-- PRE-CALC.                                                           
      PRE=1.0D0/RE                                                      
      PDT=1.0D0/DT                                                      
C                                                                       
      DO J=0,JG+1                                                       
       DO I=0,IG+1                                                      
        U(I,J,1)=1.0D0*BETA(I,J,1)                                      
        U(I,J,2)=0.0D0*BETA(I,J,2)                                      
        P(I,J)=0.0D0                                                    
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
C-------------------------------------------------------------          
C--- STREAMWISE DIRECTION                                               
      DO J=1,JG                                                         
       DO I=1,IG                                                        
        IF(BETA(I,J,1).NE.1.0) GOTO 100                                 
        UDUDX=0.5D0                                                     
     &        *((U(I+1,J  ,1)+U(I  ,J,1))                               
     &         *(U(I+1,J  ,1)-U(I  ,J,1))*PDX(I+1)*DX(I)                
     &         +(U(I  ,J  ,1)+U(I-1,J,1))                               
     &         *(U(I  ,J  ,1)-U(I-1,J,1))*PDX(I)*DX(I+1)                
     &         )*PDX2(I)                                                
        VDUDY=0.5D0                                                     
     &        *((U(I+1,J  ,2)*DX(I)+U(I,J  ,2)*DX(I+1))                 
     &         *(U(I  ,J+1,1)-U(I,J  ,1))*PDYY(J+1)                     
     &         +(U(I+1,J-1,2)*DX(I)+U(I,J-1,2)*DX(I+1))                 
     &         *(U(I  ,J  ,1)-U(I,J-1,1))*PDYY(J  )                     
     &         )*PDX2(I)                                                
        U_CON=UDUDX+VDUDY                                               
C                                                                       
        UVISX=(DX(I+1)*U(I-1,J,1)-DX2(I)*U(I,J,1)+DX(I)*U(I+1,J,1)      
     &        )*2.0D0*PDX3(I)*PRE                                       
        UVISY=(DYY(J+1)*U(I,J-1,1)-DYY2(J)*U(I,J,1)+DYY(J)*U(I,J+1,1)   
     &        )*2.0D0*PDYY3(J)*PRE                                      
        U_VIS=UVISX+UVISY                                               
C                                                                       
        DPDX=(P(I+1,J)-P(I,J))*PDXX(I+1)                                
C                                                                       
        UT(I,J,1)=(U(I,J,1)-DT*(U_CON-U_VIS+DPDX))                      
C                                                                       
        CNX_OLD(I,J)=U_CON-U_VIS+DPDX                                   
  100  ENDDO                                                            
      ENDDO                                                             
C                                                                       
C--- VERTICAL DIRECTION                                                 
      DO J=1,JG-1                                                       
       DO I=1,IG                                                        
        IF(BETA(I,J,2).NE.1.0) GOTO 200                                 
        UDVDX=0.5D0                                                     
     &      *((U(I  ,J+1,1)*DY(J)+U(I  ,J,1)*DY(J+1))                   
     &       *(U(I+1,J  ,2)-U(I  ,J,2))*PDXX(I+1)                       
     &       +(U(I-1,J+1,1)*DY(J)+U(I-1,J,1)*DY(J+1))                   
     &       *(U(I  ,J  ,2)-U(I-1,J,2))*PDXX(I  )                       
     &       )*PDY2(J)                                                  
        VDVDY=0.5D0                                                     
     &      *((U(I  ,J+1,2)+U(I,J  ,2))                                 
     &       *(U(I  ,J+1,2)-U(I,J  ,2))*PDY(J+1)*DY(J)                  
     &       +(U(I  ,J  ,2)+U(I,J-1,2))                                 
     &       *(U(I  ,J  ,2)-U(I,J-1,2))*PDY(J)*DY(J+1)                  
     &       )*PDY2(J)                                                  
        V_CON=UDVDX+VDVDY                                               
C                                                                       
        VVISX=(DXX(I+1)*U(I-1,J,2)-DXX2(I)*U(I,J,2)+DXX(I)*U(I+1,J,2)   
     &        )*2.0D0*PDXX3(I)*PRE                                      
        VVISY=(DY(J+1)*U(I,J-1,2)-DY2(J)*U(I,J,2)+DY(J)*U(I,J+1,2)      
     &        )*2.0D0*PDY3(J)*PRE                                       
        V_VIS=VVISX+VVISY                                               
C                                                                       
        DPDY=(P(I,J+1)-P(I,J))*PDYY(J+1)                                
C                                                                       
        UT(I,J,2)=(U(I,J,2)-DT*(V_CON-V_VIS+DPDY))*BETA(I,J,2)          
C                                                                       
        CNY_OLD(I,J)=V_CON-V_VIS+DPDY                                   
  200  ENDDO                                                            
      ENDDO                                                             
C                                                                       
      DO I=1,IG                                                         
       DO J=1,JG                                                        
        IF(IDX(I,J).EQ.3) UT(I,J,1)=UT(I,J+1,1)*GAMMA(I,J,1)            
        IF(IDX(I,J).EQ.5) UT(I,J,1)=UT(I,J-1,1)*GAMMA(I,J,1)            
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     CALC. TEMPORARY VELOCIYT FROM SECOND STEPS (VELOCITY FIELD)       
C     ADAMS BASHFORTH METHOD                                            
      SUBROUTINE FS1                                                    
C     2010.10.13                                                        
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
C=================================================================      
C-- PRE-CALC.                                                           
      PRE=1.0D0/RE                                                      
      PDT=1.0D0/DT                                                      
C                                                                       
C-------------------------------------------------------------          
C--- STREAMWISE DIRECTION                                               
      DO J=1,JG                                                         
       DO I=1,IG                                                        
        IF(BETA(I,J,1).NE.1.0) GOTO 100                                 
        UDUDX=0.5D0                                                     
     &        *((U(I+1,J  ,1)+U(I  ,J,1))                               
     &         *(U(I+1,J  ,1)-U(I  ,J,1))*PDX(I+1)*DX(I)                
     &         +(U(I  ,J  ,1)+U(I-1,J,1))                               
     &         *(U(I  ,J  ,1)-U(I-1,J,1))*PDX(I)*DX(I+1)                
     &         )*PDX2(I)                                                
        VDUDY=0.5D0                                                     
     &        *((U(I+1,J  ,2)*DX(I)+U(I,J  ,2)*DX(I+1))                 
     &         *(U(I  ,J+1,1)-U(I,J  ,1))*PDYY(J+1)                     
     &         +(U(I+1,J-1,2)*DX(I)+U(I,J-1,2)*DX(I+1))                 
     &         *(U(I  ,J  ,1)-U(I,J-1,1))*PDYY(J  )                     
     &         )*PDX2(I)                                                
        U_CON=UDUDX+VDUDY                                               
C                                                                       
        UVISX=(DX(I+1)*U(I-1,J,1)-DX2(I)*U(I,J,1)+DX(I)*U(I+1,J,1)      
     &        )*2.0D0*PDX3(I)*PRE                                       
        UVISY=(DYY(J+1)*U(I,J-1,1)-DYY2(J)*U(I,J,1)+DYY(J)*U(I,J+1,1)   
     &        )*2.0D0*PDYY3(J)*PRE                                      
        U_VIS=UVISX+UVISY                                               
C                                                                       
        DPDX=(P(I+1,J)-P(I,J))*PDXX(I+1)                                
C                                                                       
        UT(I,J,1)=(U(I,J,1)                                             
     &             -DT*((U_CON-U_VIS+DPDX)*1.5D0-CNX_OLD(I,J)*0.5D0)    
     &            )*BETA(I,J,1)                                         
C                                                                       
        CNX_OLD(I,J)=U_CON-U_VIS+DPDX                                   
  100  ENDDO                                                            
      ENDDO                                                             
C                                                                       
C--- VERTICAL DIRECTION                                                 
      DO J=1,JG-1                                                       
       DO I=1,IG                                                        
        IF(BETA(I,J,2).NE.1.0) GOTO 200                                 
        UDVDX=0.5D0                                                     
     &      *((U(I  ,J+1,1)*DY(J)+U(I  ,J,1)*DY(J+1))                   
     &       *(U(I+1,J  ,2)-U(I  ,J,2))*PDXX(I+1)                       
     &       +(U(I-1,J+1,1)*DY(J)+U(I-1,J,1)*DY(J+1))                   
     &       *(U(I  ,J  ,2)-U(I-1,J,2))*PDXX(I  )                       
     &       )*PDY2(J)                                                  
        VDVDY=0.5D0                                                     
     &      *((U(I  ,J+1,2)+U(I,J  ,2))                                 
     &       *(U(I  ,J+1,2)-U(I,J  ,2))*PDY(J+1)*DY(J)                  
     &       +(U(I  ,J  ,2)+U(I,J-1,2))                                 
     &       *(U(I  ,J  ,2)-U(I,J-1,2))*PDY(J)*DY(J+1)                  
     &       )*PDY2(J)                                                  
        V_CON=UDVDX+VDVDY                                               
C                                                                       
        VVISX=(DXX(I+1)*U(I-1,J,2)-DXX2(I)*U(I,J,2)+DXX(I)*U(I+1,J,2)   
     &        )*2.0D0*PDXX3(I)*PRE                                      
        VVISY=(DY(J+1)*U(I,J-1,2)-DY2(J)*U(I,J,2)+DY(J)*U(I,J+1,2)      
     &        )*2.0D0*PDY3(J)*PRE                                       
        V_VIS=VVISX+VVISY                                               
C                                                                       
        DPDY=(P(I,J+1)-P(I,J))*PDYY(J+1)                                
C                                                                       
        UT(I,J,2)=(U(I,J,2)                                             
     &             -DT*((V_CON-V_VIS+DPDY)*1.5D0-CNY_OLD(I,J)*0.5D0)    
     &            )*BETA(I,J,2)                                         
C                                                                       
        CNY_OLD(I,J)=V_CON-V_VIS+DPDY                                   
  200  ENDDO                                                            
      ENDDO                                                             
C                                                                       
      DO I=1,IG                                                         
       DO J=1,JG                                                        
        IF(IDX(I,J).EQ.3) UT(I,J,1)=UT(I,J+1,1)*GAMMA(I,J,1)            
        IF(IDX(I,J).EQ.5) UT(I,J,1)=UT(I,J-1,1)*GAMMA(I,J,1)            
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C=================================================== POISSON EQ. FOR NWT
C     CALC. POISSON'S EQUATION (PRESSURE FIELD)                         
C     SOR (SUCCESSIVE OVER RELAXATION ITERATIVE) METHOD                 
      SUBROUTINE FS2(NCOUNT)                                            
C     2010.10.13                                                        
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
      DIMENSION PHI_NEW(0:IG+1,0:JG+1)                                  
C                                                                       
C=================================================================      
C
C
      IGF=6
C---- FRINGE CHODEWISE VEL.                                           
C
      DO J=1,JG                                                       
       DO I=IG*(IGF-1)/IGF+1,IG                                                        
       UT(I,J,1)=(UT(I,J,1)-1.0D0)
     &*(-DTANH(4.5D0*(X(I)-XL)))
     &          +1.0D0
      ENDDO                                                             
      ENDDO                                                             
C
C---- FRINGE SPANWISE VEL.                                           
C
      DO J=1,JG                                                       
       DO I=IG*(IGF-1)/IGF+1,IG                                                        
       UT(I,J,2)=UT(I,J,2)
     &*(-DTANH(4.5D0*(X(I)-XL)))
      ENDDO                                                             
      ENDDO                                                             
C
C---- VERTICAL DIRECTION B.C.                                           
      DO I=1,IG                                                         
       UT(I,   0,1)= UT(I,1,1)                                          
       UT(I,JG+1,1)= UT(I,JG,1)                                         
       UT(I,   0,2)= 0.0D0                                              
       UT(I,JG  ,2)= 0.0D0                                              
      ENDDO                                                             
C                                                                       
C---- STREAMWISE DIRECTION B.C.                                         
      DO J=1,JG                                                         
       UT(   0,J,1)= 1.0D0                                              
       UT(IG+1,J,1)= UT(IG,J,1)                                         
       UT(   0,J,2)=-UT( 1,J,2)                                         
       UT(IG+1,J,2)=-UT(IG,J,2)                                         
      ENDDO                                                             
C                                                                       
C------------------------------------ SOR                               
C                                                                       
      PDT=1.0D0/DT                                                      
      OMG=1.4D0                                                         
      EPS=1.0E-4                                                        
C                                                                       
      DO NCOUNT=1,50000                                                 
C                                                                       
      ERMAX=0.0D0                                                       
C                                                                       
      DO I=1,IG                                                         
       DO J=1,JG                                                        
        PHI_NEW(I,J)=(                                                  
     &                 D_W(I)*PHI(I-1,J)                                
     &                +D_E(I)*PHI(I+1,J)                                
     &                +D_S(J)*PHI(I,J-1)                                
     &                +D_N(J)*PHI(I,J+1)                                
     &                -                                                 
     &                PDT*((UT(I,J,1)-UT(I-1,J,1))*PDX(I)               
     &                    +(UT(I,J,2)-UT(I,J-1,2))*PDY(J))              
     &               )*PDPHI(I,J)                                       
        PHIDIF=PHI_NEW(I,J)-PHI(I,J)                                    
        PHI(I,J)=PHI(I,J)+OMG*PHIDIF                                    
        ERMAX=DMAX1(ERMAX,DABS(PHIDIF))                                 
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
C---- PHI B.C.                                                          
      DO I=1,IG                                                         
       PHI(I,   0)=PHI(I,   1)                                          
       PHI(I,JG+1)=PHI(I,JG  )                                          
      ENDDO                                                             
      DO J=1,JG                                                         
       PHI(   0,J)=PHI(1   ,J)                                          
       PHI(IG+1,J)=PHI(IG  ,J)                                          
      ENDDO                                                             
C                                                                       
      IF(ERMAX.LT.EPS) GO TO 100                                        
C                                                                       
      ENDDO                                                             
  100 CONTINUE                                                          
C                                                                       
C-------------------------------------- VELOCITY  NEW <- OLD            
C                                                                       
      DO J=1,JG                                                         
       DO I=1,IG                                                        
        U(I,J,1)=UT(I,J,1)-DT*(PHI(I+1,J)-PHI(I,J))*PDXX(I+1)           
        U(I,J,2)=UT(I,J,2)-DT*(PHI(I,J+1)-PHI(I,J))*PDYY(J+1)           
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
C-------------------------------------- PRESSURE NORMALIZATION          
C                                                                       
      P_STN=P(1,1)+PHI(1,1)                                             
      DO J=1,JG                                                         
       DO I=1,IG                                                        
        P(I,J)=P(I,J)+PHI(I,J)-P_STN                                    
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     BOUNDARY CONDITION FOR VELOCITY FIELD                             
      SUBROUTINE UBC(DIV,CFL,USUM,CL,CD)                                
C     2010.10.13                                                        
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
      CHARACTER*4 NSTEP                                                 
C                                                                       
C==============================================================         
C---- VERTICAL DIRECTION B.C.                                           
      DO I=1,IG                                                         
       U(I,   0,1)= U(I, 1,1)                                           
       U(I,JG+1,1)= U(I,JG,1)                                           
       U(I,   0,2)= 0.0D0                                               
       U(I,JG  ,2)= 0.0D0                                               
      ENDDO                                                             
C                                                                       
C---- STREAMWISE DIRECTION B.C.                                         
      DO J=1,JG                                                         
       U(   0,J,1)= 1.0D0                                               
       U(IG+1,J,1)= U(IG,J,1)                                           
       U(   0,J,2)=-U( 1,J,2)                                           
       U(IG+1,J,2)=-U(IG,J,2)                                           
      ENDDO                                                             
C                                                                       
C---- PRESSURE B.C.                                                     
C                                                                       
      DO 620 I=0,IG                                                     
        P(I,   0)=P(I,   1)                                             
        P(I,JG+1)=P(I,JG  )                                             
  620 CONTINUE                                                          
      DO 55 J=1,JG                                                      
        P(   0,J)=P(   1,J)                                             
        P(IG+1,J)=P(IG  ,J)                                             
   55 CONTINUE                                                          
C                                                                       
C------------------------------------ CONTINUITY CHECKER                
      DIV=0.0D0                                                         
C                                                                       
      DO 40 J=1,JG                                                      
      DO 40 I=1,IG                                                      
        QQQQ= (U(I,J,1)-U(I-1,J,1))*PDX(I)                              
     &       +(U(I,J,2)-U(I,J-1,2))*PDY(J)                              
        QQQ=DABS(QQQQ)                                                  
        DIV=DMAX1(DIV,QQQ)                                              
   40 CONTINUE                                                          
C                                                                       
C------------------------------------ CFL CHECKER                       
      CFL=0.0D0                                                         
C                                                                       
      DO 50 J=1,JG                                                      
      DO 50 I=1,IG                                                      
        CFLL=(DABS(U(I,J,1)+U(I-1,J,1))*PDX(I)                          
     &       +DABS(U(I,J,2)+U(I,J-1,2))*PDY(J)                          
     &       )*0.5D0*DT                                                 
        CFL=DMAX1(CFL,CFLL)                                             
   50 CONTINUE                                                          
C                                                                       
C------------------------------------ MEAN VELOCITY CHECKER             
      USUM=0.0D0                                                        
C                                                                       
      DO 60 J=1,JG                                                      
      DO 60 I=1,IG                                                      
        USUM=USUM+U(I,J,1)*DY(J)                                        
   60 CONTINUE                                                          
      USUM=USUM/DBLE(IG)/YL                                             
C                                                                       
C------------------------------------ OUTPUT INSTANTANEOUS VEL. FIELD   
C                                                                       
      IF(MOD(NUM,5).EQ.0) THEN                                        
       NONE=MOD(NUM,10)                                                 
       NTEN=MOD(NUM,100)-MOD(NUM,10)                                    
       NHUN=MOD(NUM,1000)-MOD(NUM,100)                                  
       NTHU=MOD(NUM,10000)-MOD(NUM,1000)                                
       NSTEP=CHAR(48+NTHU/1000)//CHAR(48+NHUN/100)//CHAR(48+NTEN/10)    
     & //CHAR(48+NONE)                                                  
C                                                                       
       OPEN(20, FILE='TEXTDATA/puv_field'//NSTEP//'.dat',               
     & STATUS="UNKNOWN")                                                
       DO I=1,IG                                                       
        DO J=1,JG                                                        
         WRITE(20,'(8E13.5)') XX(I),YY(J)                               
     &  ,P(I,J)                                                         
     &  ,(U(I,J,1)+U(I-1,J,1))*0.5D0                                    
     &  ,(U(I,J,2)+U(I,J-1,2))*0.5D0                                    
     &  ,(BETA(I,J,2)+BETA(I,J-1,2))*0.5D0                              
     &  ,(BETA(I,J,1)+BETA(I-1,J,1))*0.5D0                              
     &  ,-((U(I,J,1)+U(I-1,J,1))*0.5D0-(U(I,J-1,1)+U(I-1,J-1,1))*0.5D0)  
     &   /DY(J)
     &  +((U(I,J,2)+U(I,J-1,2))*0.5D0-(U(I-1,J,2)+U(I-1,J-1,2))*0.5D0)  
     &   /DX(I)                              
        ENDDO                                                           
       ENDDO                                                            
       CLOSE(20)                                                        
      ENDIF                                                             
C                                                                       
C------------------------------------ CALC. CL AND CD                   
C                                                                       
      P_SUM=0.0D0                                                       
      DO I=1,IG                                                         
       P_TOP=0.0D0                                                      
       P_BTM=0.0D0                                                      
       DO J=2,JG                                                        
C        IF(BETA(I,J,2)-BETA(I,J-4,2).GE. 0.50D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,2)-BETA(I,J-3,2).GE. 0.50D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,2)-BETA(I,J-2,2).GE. 0.50D0) P_TOP=P(I,J)            
        IF(BETA(I,J,2)-BETA(I,J-1,2).GE. 0.50D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,2)-BETA(I,J-4,2).LE.-0.50D0) P_BTM=P(I,J-4)            
C        IF(BETA(I,J,2)-BETA(I,J-3,2).LE.-0.50D0) P_BTM=P(I,J-3)            
C        IF(BETA(I,J,2)-BETA(I,J-2,2).LE.-0.50D0) P_BTM=P(I,J-2)            
        IF(BETA(I,J,2)-BETA(I,J-1,2).LE.-0.50D0) P_BTM=P(I,J-1)              
       ENDDO                                                            
       P_SUM=P_SUM+DX(I)*(P_BTM-P_TOP)                                  
      ENDDO                                                             
      CL=P_SUM*2.0D0                                                    
C                                                                       
      P_SUM =0.0D0                                                      
      DO J=1,JG                                                         
       P_TOP=0.0D0                                                      
       P_BTM=0.0D0                                                      
       DO I=2,IG                                                        
C        IF(BETA(I,J,1)-BETA(I-4,J,1).GE. 0.5D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,1)-BETA(I-3,J,1).GE. 0.5D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,1)-BETA(I-2,J,1).GE. 0.5D0) P_TOP=P(I,J)            
        IF(BETA(I,J,1)-BETA(I-1,J,1).GE. 0.5D0) P_TOP=P(I+1,J)            
C        IF(BETA(I,J,1)-BETA(I-4,J,1).LE.-0.5D0) P_BTM=P(I-4,J)            
C        IF(BETA(I,J,1)-BETA(I-3,J,1).LE.-0.5D0) P_BTM=P(I-3,J)            
C        IF(BETA(I,J,1)-BETA(I-2,J,1).LE.-0.5D0) P_BTM=P(I-2,J)            
        IF(BETA(I,J,1)-BETA(I-1,J,1).LE.-0.5D0) P_BTM=P(I-2,J)            
       ENDDO                                                            
       P_SUM=P_SUM+DY(J)*(P_BTM-P_TOP)                                  
      ENDDO                                                             
      CD=P_SUM*2.0D0                                                    
C                                                                       
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     MESH GENERATION AND BETA CALCULATION                              
      SUBROUTINE CALBETA(IERR,NWING)                                    
C     2010.10.13                                                        
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
      CHARACTER*14 NWING                                                
C                                                                       
      DIMENSION WNG_XT(IGS), WNG_XB(IGS)                                
      DIMENSION WNG_YT(IGS), WNG_YB(IGS)                                
      DIMENSION X_TEMP(-3:LX+5)                                         
      DIMENSION Y_TEMP(-3:LY+5)                                         
C                                                                       
C==============================================================         
C                                                                       
      THETA=PI*DEG/180.0D0                                              
 1000 FORMAT(A14)                                                       
 1010 FORMAT(3F9.5)                                                     
C                                                                       
c      OPEN(10, FILE='clark_y.dat', STATUS="OLD")                        
      READ(10,1000) NWING                                               
      READ(10,1010) (WNG_XT(II), WNG_YT(II), WNG_YB(II), II=1,IGS)      
      CLOSE(10)                                                         
C                                                                       
      DO II=1,IGS                                                       
        TWNG_XT=WNG_XT(II)*DCOS(-THETA)-WNG_YT(II)*DSIN(-THETA)         
        TWNG_YT=WNG_XT(II)*DSIN(-THETA)+WNG_YT(II)*DCOS(-THETA)         
        TWNG_XB=WNG_XT(II)*DCOS(-THETA)-WNG_YB(II)*DSIN(-THETA)         
        TWNG_YB=WNG_XT(II)*DSIN(-THETA)+WNG_YB(II)*DCOS(-THETA)         
C                                                                       
        WNG_XT(II)=TWNG_XT+XL*0.5D0-DCOS(THETA)*0.5D0                   
        WNG_XB(II)=TWNG_XB+XL*0.5D0-DCOS(THETA)*0.5D0                   
        WNG_YT(II)=TWNG_YT+YL*0.5D0+DSIN(THETA)*0.5D0-0.06D0            
        WNG_YB(II)=TWNG_YB+YL*0.5D0+DSIN(THETA)*0.5D0-0.06D0            
      ENDDO                                                             
C                                                                       
C------------------------------------------------------------- NEW MESH 
C                                                                       
      XALP1=1.0D0                                                       
      XALP2=3.0D0                                                       
      DO I=1,IG                                                         
        X1=2.0D0*DBLE(I)/DBLE(IG)-1.0D0                                 
        X_TEMP(I)=XL*(                                                  
     &                ((XALP1*(X1**XALP2)+X1)*0.5D0)                    
     &                /(XALP1+1.0D0)+0.5D0                              
     &               )                                                  
      ENDDO                                                             
      YALP1=1.0D0                                                       
      YALP2=5.0D0                                                       
      DO J=1,JG                                                         
        Y1=2.0D0*DBLE(J)/DBLE(JG)-1.0D0                                 
        Y_TEMP(J)=YL*(                                                  
     &                ((YALP1*(Y1**YALP2)+Y1)*0.5D0)                    
     &                /(YALP1+1.0D0)+0.5D0                              
     &               )                                                  
      ENDDO                                                             
C                                                                       
      X_TEMP(   0)=0.0D0                                                
      X_TEMP(  -1)=-X_TEMP(1)                                           
      X_TEMP(  -2)=-X_TEMP(2)                                           
      X_TEMP(  -3)=-X_TEMP(3)                                           
      X_TEMP(IG+1)=2.0D0*XL-X_TEMP(IG-1)                                
      X_TEMP(IG+2)=2.0D0*XL-X_TEMP(IG-2)                                
      X_TEMP(IG+3)=2.0D0*XL-X_TEMP(IG-3)                                
      X_TEMP(IG+4)=2.0D0*XL-X_TEMP(IG-4)                                
      X_TEMP(IG+5)=2.0D0*XL-X_TEMP(IG-5)                                
      Y_TEMP(   0)=0.0D0                                                
      Y_TEMP(  -1)=-Y_TEMP(1)                                           
      Y_TEMP(  -2)=-Y_TEMP(2)                                           
      Y_TEMP(  -3)=-Y_TEMP(3)                                           
      Y_TEMP(JG+1)=2.0D0*YL-Y_TEMP(JG-1)                                
      Y_TEMP(JG+2)=2.0D0*YL-Y_TEMP(JG-2)                                
      Y_TEMP(JG+3)=2.0D0*YL-Y_TEMP(JG-3)                                
      Y_TEMP(JG+4)=2.0D0*YL-Y_TEMP(JG-4)                                
      Y_TEMP(JG+5)=2.0D0*YL-Y_TEMP(JG-5)                                
C                                                                       
      DO I=-2,LX+5                                                      
        X(I)  =X_TEMP(I)                                                
        DX(I) =X(I)-X(I-1)                                              
        XX(I) =(X(I)+X(I-1))*0.5D0                                      
        DXX(I)=XX(I)-XX(I-1)                                            
      ENDDO                                                             
      DO J=-2,LY+5                                                      
        Y  (J)=Y_TEMP(J)                                                
        DY (J)=Y(J)-Y(J-1)                                              
        YY (J)=(Y(J)+Y(J-1))*0.5D0                                      
        DYY(J)=YY(J)-YY(J-1)                                            
      ENDDO                                                             
C                                                                       
C------------------------------------------------------------- CAL BETA 
C                                                                       
      DO I=0,IG                                                         
       DO J=0,JG                                                        
        BETA(I,J,1)=1.0D0                                               
        BETA(I,J,2)=1.0D0                                               
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
      DO I=0,IG                                                         
       IB=0                                                             
       IT=0                                                             
       DO II=1,IGS                                                      
        IF((X(I).LT.WNG_XB(II+1)).AND.                                  
     &     (X(I).GE.WNG_XB(II))) IB=II                                  
        IF((X(I).LT.WNG_XT(II+1)).AND.                                  
     &     (X(I).GE.WNG_XT(II))) IT=II                                  
       ENDDO                                                            
       IF(IB*IT.EQ.0) GO TO 110                                         
       XTA=X(I)-WNG_XT(IT)                                              
       XBA=X(I)-WNG_XB(IB)                                              
       XTB=WNG_XT(IT+1)-X(I)                                            
       XBB=WNG_XB(IB+1)-X(I)                                            
       YT=(WNG_YT(IT)*XTB+WNG_YT(IT+1)*XTA)/(XTA+XTB)                   
       YB=(WNG_YB(IB)*XBB+WNG_YB(IB+1)*XBA)/(XBA+XBB)                   
       DO J=0,JG                                                        
        IF((YY(J).GE.YB).AND.(YY(J).LT.YT)) BETA(I,J,1)=0.0D0           
        IF((YY(J).LE.YB).AND.(YY(J+1).GT.YB)) THEN                      
         GAMMA(I,J,1)=(YB-YY(J))/(YB-YY(J-1))                           
        ENDIF                                                           
        IF((YY(J).GE.YT).AND.(YY(J-1).LT.YT)) THEN                      
         GAMMA(I,J,1)=(YY(J)-YT)/(YY(J+1)-YT)                           
        ENDIF                                                           
       ENDDO                                                            
  110  CONTINUE                                                         
      ENDDO                                                             
C                                                                       
      DO I=0,IG                                                         
       IB=0                                                             
       IT=0                                                             
       DO II=1,IGS                                                      
        IF((XX(I).LT.WNG_XB(II+1)).AND.                                 
     &     (XX(I).GE.WNG_XB(II))) IB=II                                 
        IF((XX(I).LT.WNG_XT(II+1)).AND.                                 
     &     (XX(I).GE.WNG_XT(II))) IT=II                                 
       ENDDO                                                            
       IF(IB*IT.EQ.0) GO TO 120                                         
       XTA=XX(I)-WNG_XT(IT)                                             
       XBA=XX(I)-WNG_XB(IB)                                             
       XTB=WNG_XT(IT+1)-XX(I)                                           
       XBB=WNG_XB(IB+1)-XX(I)                                           
       YT=(WNG_YT(IT)*XTB+WNG_YT(IT+1)*XTA)/(XTA+XTB)                   
       YB=(WNG_YB(IB)*XBB+WNG_YB(IB+1)*XBA)/(XBA+XBB)                   
       DO J=0,JG                                                        
        IF((Y(J).GE.YB).AND.(Y(J).LT.YT)) BETA(I,J,2)=0.0D0             
       ENDDO                                                            
  120  CONTINUE                                                         
      ENDDO                                                             
C                                                                       
      WRITE(*,*)"GRID GENERATION OK!"                                   
C                                                                       
C-------------------------------------------------------------          
C                                                                       
      DO I=1,IG                                                         
       DO J=1,JG-1                                                      
        IDX(I,J)=4*IDINT(BETA(I,J-1,1))                                 
     &          +2*IDINT(BETA(I,J  ,1))                                 
     &          +1*IDINT(BETA(I,J+1,1))                                 
        IDY(I,J)=4*IDINT(BETA(I-1,J,2))                                 
     &          +2*IDINT(BETA(I  ,J,2))                                 
     &          +1*IDINT(BETA(I+1,J,2))                                 
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
      DO I=0,IG+1                                                       
       DX2(I)=DX(I)+DX(I+1)                                             
       DXX2(I)=DXX(I)+DXX(I+1)                                          
       PDX(I)=1.0D0/DX(I)                                               
       PDXX(I)=1.0D0/DXX(I)                                             
       PDX2(I)=1.0D0/(DX(I)+DX(I+1))                                    
       PDX3(I)=1.0D0/(DX(I)*DX(I+1)*(DX(I)+DX(I+1)))                    
       PDXX3(I)=1.0D0/(DXX(I)*DXX(I+1)*(DXX(I)+DXX(I+1)))               
      ENDDO                                                             
      DO J=0,JG+1                                                       
       DY2(J)=DY(J)+DY(J+1)                                             
       DYY2(J)=DYY(J)+DYY(J+1)                                          
       PDY(J)=1.0D0/DY(J)                                               
       PDYY(J)=1.0D0/DYY(J)                                             
       PDY2(J)=1.0D0/(DY(J)+DY(J+1))                                    
       PDY3(J)=1.0D0/(DY(J)*DY(J+1)*(DY(J)+DY(J+1)))                    
       PDYY3(J)=1.0D0/(DYY(J)*DYY(J+1)*(DYY(J)+DYY(J+1)))               
      ENDDO                                                             
C                                                                       
      DO I=1,IG                                                         
       DO J=1,JG                                                        
        D_W(I)=PDX(I)*PDXX(I)                                           
        D_E(I)=PDX(I)*PDXX(I+1)                                         
        D_S(J)=PDY(J)*PDYY(J)                                           
        D_N(J)=PDY(J)*PDYY(J+1)                                         
        PDPHI(I,J)=1.0D0/(D_W(I)+D_E(I)+D_S(J)+D_N(J))                  
       ENDDO                                                            
      ENDDO                                                             
C                                                                       
C------------------------------------ DIFFUSION NUMBER CHECKER          
      DFN=0.0D0                                                         
C                                                                       
      DO I=1,IG                                                         
       DO J=1,JG                                                        
        DFNN=DT/RE*(1.0D0/(DX(I)**2)+1.0D0/(DY(J)**2))                  
        DFN=DMAX1(DFN,DFNN)                                             
       ENDDO                                                            
      ENDDO                                                             
      WRITE(*,'(A20,E11.4)')"DIFFUSION NUMBER = ",DFN                    
C                                                                       
    1 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C=======================================================================
C     OUTPUT TIME-AVERAGE DATA (VELOCITY FIELD)                         
      SUBROUTINE OUTPUT(NWING)                                          
C     2010.10.13                                                        
C=======================================================================
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
      CHARACTER*4 NSTEP                                                 
      CHARACTER*14 NWING                                                
      CHARACTER*71 CPRM(99),CMTD(3)                                     
C                                                                       
C=======================================================================
C------------------------------------ OUTPUT INSTANTANEOUS VEL. FIELD   
C                                                                       
      IF(MOD(NUM,5).EQ.0) THEN                                        
       NONE=MOD(NUM,10)                                                 
       NTEN=MOD(NUM,100)-MOD(NUM,10)                                    
       NHUN=MOD(NUM,1000)-MOD(NUM,100)                                  
       NTHU=MOD(NUM,10000)-MOD(NUM,1000)                                
       NSTEP=CHAR(48+NTHU/1000)//CHAR(48+NHUN/100)//CHAR(48+NTEN/10)    
     & //CHAR(48+NONE)                                                  
C                                                                       
       OPEN(20, FILE='TEXTDATA/puv_field'//NSTEP//'.dat',               
     & STATUS="UNKNOWN")                                                
       DO I=1,IG                                                       
        DO J=1,JG                                                        
         WRITE(20,'(8E13.5)') XX(I),YY(J)                               
     &  ,P(I,J)                                                         
     &  ,(U(I,J,1)+U(I-1,J,1))*0.5D0                                    
     &  ,(U(I,J,2)+U(I,J-1,2))*0.5D0                                    
     &  ,(BETA(I,J,2)+BETA(I,J-1,2))*0.5D0                              
     &  ,(BETA(I,J,1)+BETA(I-1,J,1))*0.5D0                              
     &  ,-((U(I,J,1)+U(I-1,J,1))*0.5D0-(U(I,J-1,1)+U(I-1,J-1,1))*0.5D0)  
     &   /DY(J)
     &  +((U(I,J,2)+U(I,J-1,2))*0.5D0-(U(I-1,J,2)+U(I-1,J-1,2))*0.5D0)  
     &   /DX(I)                              
        ENDDO                                                           
       ENDDO                                                            
       CLOSE(20)                                                        
      ENDIF                                                             
C                                                                       
C------------------------------------ OUTPUT INSTANTANEOUS VEL. FIELD   
C                                                                       
      CMTD( 1)='#======================================================'
      CMTD( 2)='# SIMULATION OF FLOW PAST AN 2D AIRFOIL : '//NWING      
      CMTD( 3)='# STATISTICS IN VELOCITY FIELD'                         
      CPRM( 1)='# CALC.NUM.     CALC.TIME'                              
      CPRM( 2)='# NX  NY  TIMESTEP  REYNOLDS'                           
      CPRM( 3)='# X-LENGTH  Y-LENGTH AOA(DEG)'                          
      CPRM( 4)='# I X DX XX DXX'                                        
      CPRM( 5)='# J Y DY YY DYY'                                        
      CPRM( 6)='# X DX P(UPPER-SURFACE)  P(BOTTOM SURFACE)'             
      CPRM( 7)='# Y DY P(DOWNSTREAM-SURFACE)  P(UPSTREAM SURFACE)'      
C                                                                       
 1000 FORMAT(I10,E18.10)                                                
 1010 FORMAT(A71)                                                       
 1020 FORMAT(I4,I4,F11.6,F9.1)                                          
 1030 FORMAT(2F10.6,F7.2)                                               
 1040 FORMAT(5F10.7)                                                    
 1050 FORMAT(2F10.7)                                                    
 1060 FORMAT(I3,4E14.6)                                                 
 1070 FORMAT(4E14.6)                                                    
C                                                                       
      OPEN(14,FILE='TEXTDATA/output.dat',STATUS="UNKNOWN")              
C                                                                       
      WRITE(14,1010) CMTD(1)                                            
      WRITE(14,1010) CMTD(2)                                            
      WRITE(14,1010) CMTD(1)                                            
      WRITE(14,1010) CPRM(1)                                           
      WRITE(14,1000) NUM,TIME                               
      WRITE(14,1010) CPRM(2)                                            
      WRITE(14,1020) IG,JG,DT,RE                                  
      WRITE(14,1010) CPRM(3)                                            
      WRITE(14,1030) XL,YL,DEG                                   
      WRITE(14,1010) CMTD(3)                                            
C
      WRITE(14,1010) CPRM(4)                                            
      DO I=0,IG+1                                                       
        WRITE(14,1060) I,X(I),DX(I),XX(I),DXX(I)                        
      ENDDO                                                             
C                                                                       
      WRITE(14,1010) CPRM(5)                                            
      DO J=0,JG+1                                                       
        WRITE(14,1060) J,Y(J),DY(J),YY(J),DYY(J)                        
      ENDDO                                                             
C                                                                       
      WRITE(14,1010) CPRM(6)                                            
C------------------------------------ CALC. CL AND CD                   
C                                                                       
      P_SUM=0.0D0                                                       
      DO I=1,IG                                                         
       P_TOP=0.0D0                                                      
       P_BTM=0.0D0                                                      
       DO J=2,JG                                                        
C        IF(BETA(I,J,2)-BETA(I,J-4,2).GE. 0.50D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,2)-BETA(I,J-3,2).GE. 0.50D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,2)-BETA(I,J-2,2).GE. 0.50D0) P_TOP=P(I,J)            
        IF(BETA(I,J,2)-BETA(I,J-1,2).GE. 0.50D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,2)-BETA(I,J-4,2).LE.-0.50D0) P_BTM=P(I,J-4)            
C        IF(BETA(I,J,2)-BETA(I,J-3,2).LE.-0.50D0) P_BTM=P(I,J-3)            
C        IF(BETA(I,J,2)-BETA(I,J-2,2).LE.-0.50D0) P_BTM=P(I,J-2)            
        IF(BETA(I,J,2)-BETA(I,J-1,2).LE.-0.50D0) P_BTM=P(I,J-1)              
       ENDDO                                                            
       IF(P_TOP.NE.0.0D0) WRITE(14,1070) X(I),DX(I),P_TOP,P_BTM         
       P_SUM=P_SUM+DX(I)*(P_BTM-P_TOP)                                  
      ENDDO                                                             
      CL=P_SUM*2.0D0                                                    
C                                                                       
      WRITE(14,1010) CPRM(7)                                            
      P_SUM =0.0D0                                                      
      DO J=1,JG                                                         
       P_TOP=0.0D0                                                      
       P_BTM=0.0D0                                                      
       DO I=2,IG                                                        
C        IF(BETA(I,J,1)-BETA(I-4,J,1).GE. 0.5D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,1)-BETA(I-3,J,1).GE. 0.5D0) P_TOP=P(I,J)            
C        IF(BETA(I,J,1)-BETA(I-2,J,1).GE. 0.5D0) P_TOP=P(I,J)            
        IF(BETA(I,J,1)-BETA(I-1,J,1).GE. 0.5D0) P_TOP=P(I+1,J)            
C        IF(BETA(I,J,1)-BETA(I-4,J,1).LE.-0.5D0) P_BTM=P(I-4,J)            
C        IF(BETA(I,J,1)-BETA(I-3,J,1).LE.-0.5D0) P_BTM=P(I-3,J)            
C        IF(BETA(I,J,1)-BETA(I-2,J,1).LE.-0.5D0) P_BTM=P(I-2,J)            
        IF(BETA(I,J,1)-BETA(I-1,J,1).LE.-0.5D0) P_BTM=P(I-2,J)            
       ENDDO                                                            
       IF(P_TOP.NE.0.0D0) WRITE(14,1070) Y(J),DY(J),P_TOP,P_BTM         
       P_SUM=P_SUM+DY(J)*(P_BTM-P_TOP)                                  
      ENDDO                                                             
      CD=P_SUM*2.0D0                                                    
C                                                                       
      WRITE(*,'(A6,F6.3,A8,F6.3)') 'C_L =',CL,', C_D =',CD              
C                                                                       
      CLOSE(14)                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     SEED DATA (VELOCITY FIELD) : INPUT & OUTPUT                       
C====================================================================   
C                                                                       
      SUBROUTINE SEEDIN(IERR)                                           
C                                                                       
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
      CHARACTER*71 CMTD                                                 
C                                                                       
C=====================================================================  
C                                                                       
      OPEN(31, FILE="seed_u.dat", STATUS="OLD")                         
      OPEN(32, FILE="seed_p.dat", STATUS="OLD")                         
      OPEN(33, FILE="seed_old.dat", STATUS="OLD")                       
C                                                                       
C------------------------------- DATAIN                                 
C                                                                       
      READ(31,*) NUM,TIME                                               
      READ(31,*) IGD,JGD,XLD,YLD,DT,RE                                  
 1010 FORMAT(A71)                                                       
C                                                                       
      IERR=0                                                            
      IF ((IGD.NE.IG).OR.(JGD.NE.JG)) THEN                              
        IERR=1000                                                       
        CMTD=' MESH DATA ERROR'                                         
        WRITE(6,1010) CMTD                                              
        GOTO 1                                                          
      ENDIF                                                             
      IF ((XLD.NE.XL).OR.(XLD.NE.XL)) THEN                              
        IERR=1000                                                       
        CMTD=' DOMAIN DATA ERROR'                                         
        WRITE(6,1010) CMTD                                              
        GOTO 1                                                          
      ENDIF                                                             
C                                                                       
      READ(31,*) U                                                      
      READ(32,*) P                                                      
      READ(33,*) CNX_OLD,CNY_OLD                                        
C                                                                       
      CLOSE(31)                                                         
      CLOSE(32)                                                         
      CLOSE(33)                                                         
C                                                                       
      WRITE(*,*)"READ VELOCITY-SEED DATA OK!"                           
C                                                                       
    1 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C                                                                       
      SUBROUTINE SEEDOUT                                                
C                                                                       
C====================================================================   
C                                                                       
      INCLUDE "./INCFILE/incinit.f"                                     
      INCLUDE "./INCFILE/incdata.f"                                     
      INCLUDE "./INCFILE/incmesh.f"                                     
C                                                                       
C=====================================================================  
C                                                                       
      OPEN(31, FILE="seed_u.dat", STATUS="UNKNOWN")                     
      OPEN(32, FILE="seed_p.dat", STATUS="UNKNOWN")                     
      OPEN(33, FILE="seed_old.dat", STATUS="UNKNOWN")                   
C                                                                       
C------------------------------- DATAOUT                                
C                                                                       
      WRITE(31,*) NUM,TIME                                              
      WRITE(31,*) IG,JG,XL,YL,DT,RE                                           
      WRITE(31,*) U                                                     
      WRITE(32,*) P                                                     
      WRITE(33,*) CNX_OLD,CNY_OLD                                       
C                                                                       
      CLOSE(31)                                                         
      CLOSE(32)                                                         
      CLOSE(33)                                                         
C                                                                       
      RETURN                                                            
      END                                                               
