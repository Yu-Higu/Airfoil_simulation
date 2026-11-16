      COMMON /DAT1/ RE,DT,TIME                                          
      COMMON /DAT2/ XL,YL,DEG                                           
      COMMON /STEP/ NDRAW,NFIN,NUM                                      
C                                                                       
      COMMON /VELO/ U(0:IG+1,0:JG+1,2)                                  
      COMMON /TEMP/ UT(0:IG+1,0:JG+1,2)                                 
      COMMON /OLDC/ CNX_OLD(0:IG,0:JG),CNY_OLD(0:IG,0:JG)               
      COMMON /SCAL/ P(0:IG+1,0:JG+1)                                    
      COMMON /PHIP/ PHI(0:IG+1,0:JG+1)                                  
      COMMON /BIBM/ BETA(0:IG,0:JG,2)                                   
