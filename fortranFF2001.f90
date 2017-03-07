      PROGRAM FFT2

         common/constants/zero,pi
         COMPLEX, PARAMETER :: im = (0,1)
         
         REAL(8), DIMENSION(:),     ALLOCATABLE   :: X
         REAL(8), DIMENSION(:),     ALLOCATABLE   :: Y
         REAL(8), DIMENSION(:,:),   ALLOCATABLE   :: Z
         COMPLEX(8), DIMENSION(:,:),   ALLOCATABLE   :: FZ
         COMPLEX(8), DIMENSION(:,:),   ALLOCATABLE   :: FRZ

         REAL(8) fr1, fr2
         INTEGER i, j, k, l, N, u, v

         N = 5

         ALLOCATE(X(N))
         ALLOCATE(Y(N))
         ALLOCATE(Z(N,N))
         ALLOCATE(FZ(N,N))
         ALLOCATE(FRZ(N,N))

         zero= 0.d0
         pi= 3.141592654

         fr1 = 0.1
         fr2 = 0.5

         DO i=1,N
            
            X(i)=i
            Y(i)=i

         END DO

         
         DO i=1,N
            DO j=1,N
               
               Z(i,j) = 2*SIN(2*pi*fr1*X(i)) + 3*SIN(2*pi*fr2*Y(j))
     
            END DO
         END DO


         DO u=1,N
            DO v=1,N
               FZ(u,v) = 0.0
               DO i=1,N
                  DO j=1,N
                     
                     FZ(u,v) = FZ(u,v) + Z(i,j)                         &
                     *(COS(2*pi*(u-1)*(i-1)/N + 2*pi*(v-1)*(j-1)/N)     &
                     -im*SIN(2*pi*(u-1)*(i-1)/N + 2*pi*(v-1)*(j-1)/N))  
                  
                  END DO
               END DO
            END DO
         END DO


         DO i=1,N
            DO j=1,N
               FRZ(i,j) = 0.0
               DO u=1,N
                  DO v=1,N
                     
                     FRZ(i,j) = FRZ(i,j) + (1/(REAL(N*N,8)))*FZ(u,v)    &            
                     *(COS(2*pi*(u-1)*(i-1)/N + 2*pi*(v-1)*(j-1)/N)     &
                     +im*SIN(2*pi*(u-1)*(i-1)/N + 2*pi*(v-1)*(j-1)/N))  
                  
                  END DO
               END DO
            END DO
         END DO

         DO i=1,N
            DO j=1,N
               PRINT*, FRZ(i,j), Z(i,j)
            END DO
         END DO

         DEALLOCATE(X)
         DEALLOCATE(Y)
         DEALLOCATE(Z)
         DEALLOCATE(FZ)
         DEALLOCATE(FRZ)


      END PROGRAM FFT2
