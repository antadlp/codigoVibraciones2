      program rota123
      implicit double precision (a-h,o-z)
      parameter (m1=12)

      common/constants/zero,pi
      CHARACTER (LEN = 400) ich
      CHARACTER (LEN = 400) nmOut 
      CHARACTER (LEN = 400) nmIn
      character title*30
      dimension iswm(m1), angm(m1)
      INTEGER i, frames

      zero= 0.d0
      pi= 3.141592654
      frames = 5000 

c      strMDCYC = ADJUSTL(strMDCYC)
C      strSCFCYC = ADJUSTL(strSCFCYC)

      DO i=1,frames
         
         WRITE(ich,*) i 
        
         ich = ADJUSTL(ich)
!         nmIn = TRIM('/home/toshiba/out001/frame'//TRIM(ich)//'.xyzz')
!         OPEN(82,FILE=TRIM(nmIn),STATUS='NEW')
!         PRINT*, nmIn
!         PRINT*, ich

         call getData (i)
         call trackrota (hmx,hmy,hmz,iswm,angm,ich)
         call inveRota  (hmx,hmy,hmz,iswm,angm)
!

      END DO  

      STOP
      end
*
