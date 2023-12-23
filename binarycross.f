! creates binary vectors from raster data, constructs cross correlations 
C-----------------------------------------------------------------------
         PROGRAM binarycross   

      INTEGER*4 I J, K, L, ictr /0/
      real*8 time
      integer LOTvec (38:55, 500), fanvec (38:55, 500)
      integer cross (306,2) ! 306 = 17 x 18; entries = LOT cross, fan
c cross
C-----------------------------------------------------------------------
c initialize binary vectors
       do i = 38, 55
       do j = 1, 500 
        LOTvec (i,j) = 0
        fanvec (i,j) = 0
       end do
       end do

      OPEN(2,FILE='LEC38.LOTrast') 
       do i = 1, 3251  
         READ(2,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (38,K) = 1
       end do

      OPEN(3,FILE='LEC38.LECfanrast') 
       do i = 1, 2312  
         READ(3,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (38,K) = 1
       end do

      OPEN(4,FILE='LEC39.LOTrast') 
       do i = 1, 3226  
         READ(4,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (39,K) = 1
       end do

      OPEN(5,FILE='LEC39.LECfanrast') 
       do i = 1, 2220  
         READ(5,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (39,K) = 1
       end do

      OPEN(7,FILE='LEC40.LOTrast') 
       do i = 1, 3256  
         READ(7,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (40,K) = 1
       end do

      OPEN(8,FILE='LEC40.LECfanrast') 
       do i = 1, 2320  
         READ(8,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (40,K) = 1
       end do


      OPEN(9,FILE='LEC41.LOTrast') 
       do i = 1, 3116  
         READ(9,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (41,K) = 1
       end do

      OPEN(10,FILE='LEC41.LECfanrast') 
       do i = 1, 1827  
         READ(10,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (41,K) = 1
       end do


      OPEN(11,FILE='LEC42.LOTrast') 
       do i = 1, 3596  
         READ(11,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (42,K) = 1
       end do

      OPEN(12,FILE='LEC42.LECfanrast') 
       do i = 1, 4120  
         READ(12,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (42,K) = 1
       end do

      OPEN(13,FILE='LEC43.LOTrast') 
       do i = 1, 2786  
         READ(13,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (43,K) = 1
       end do

      OPEN(14,FILE='LEC43.LECfanrast') 
       do i = 1, 1432  
         READ(14,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (43,K) = 1
       end do

      OPEN(15,FILE='LEC44.LOTrast') 
       do i = 1, 2881  
         READ(15,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (44,K) = 1
       end do

      OPEN(16,FILE='LEC44.LECfanrast') 
       do i = 1, 1441  
         READ(16,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (44,K) = 1
       end do

      OPEN(17,FILE='LEC45.LOTrast') 
       do i = 1, 2976  
         READ(17,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (45,K) = 1
       end do

      OPEN(18,FILE='LEC45.LECfanrast') 
       do i = 1, 1511  
         READ(18,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (45,K) = 1
       end do

      OPEN(19,FILE='LEC46.LOTrast') 
       do i = 1, 3056  
         READ(19,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (46,K) = 1
       end do

      OPEN(20,FILE='LEC46.LECfanrast') 
       do i = 1, 1565  
         READ(20,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (46,K) = 1
       end do

      OPEN(21,FILE='LEC47.LOTrast') 
       do i = 1, 3196  
         READ(21,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (47,K) = 1
       end do

      OPEN(22,FILE='LEC47.LECfanrast') 
       do i = 1, 2077  
         READ(22,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (47,K) = 1
       end do

      OPEN(23,FILE='LEC48.LOTrast') 
       do i = 1, 3436  
         READ(23,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (48,K) = 1
       end do

      OPEN(24,FILE='LEC48.LECfanrast') 
       do i = 1, 3063  
         READ(24,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (48,K) = 1
       end do

      OPEN(25,FILE='LEC49.LOTrast') 
       do i = 1, 3391  
         READ(25,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (49,K) = 1
       end do

      OPEN(26,FILE='LEC49.LECfanrast') 
       do i = 1, 2805  
         READ(26,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (49,K) = 1
       end do

      OPEN(27,FILE='LEC50.LOTrast') 
       do i = 1, 3296  
         READ(27,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (50,K) = 1
       end do

      OPEN(28,FILE='LEC50.LECfanrast') 
       do i = 1, 2524  
         READ(28,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (50,K) = 1
       end do

      OPEN(29,FILE='LEC51.LOTrast') 
       do i = 1, 3511  
         READ(29,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (51,K) = 1
       end do

      OPEN(30,FILE='LEC51.LECfanrast') 
       do i = 1, 3521  
         READ(30,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (51,K) = 1
       end do

      OPEN(31,FILE='LEC52.LOTrast') 
       do i = 1, 3046  
         READ(31,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (52,K) = 1
       end do

      OPEN(32,FILE='LEC52.LECfanrast') 
       do i = 1, 1584  
         READ(32,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (52,K) = 1
       end do

      OPEN(33,FILE='LEC53.LOTrast') 
       do i = 1, 3166  
         READ(33,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (53,K) = 1
       end do

      OPEN(34,FILE='LEC53.LECfanrast') 
       do i = 1, 1890  
         READ(34,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (53,K) = 1
       end do

      OPEN(35,FILE='LEC54.LOTrast') 
       do i = 1, 3031  
         READ(35,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (54,K) = 1
       end do

      OPEN(36,FILE='LEC54.LECfanrast') 
       do i = 1, 1553  
         READ(36,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (54,K) = 1
       end do

      OPEN(37,FILE='LEC55.LOTrast') 
       do i = 1, 3091  
         READ(37,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) LOTvec (55,K) = 1
       end do

      OPEN(38,FILE='LEC55.LECfanrast') 
       do i = 1, 1590  
         READ(38,*,IOSTAT=IOTEST) time, K
        if (time.ge.1200.d0) fanvec (55,K) = 1
       end do

c compute the binary cross correlations
      do k = 38, 55
      do L = 38, 55
       if (K.ne.L) then
          ictr = ictr + 1
        do i = 1, 500
         if ((LOTvec(k,i).eq.1).and.(LOTvec(L,i).eq.1))  ! binary inner
c product
     x      cross(ictr,1) = cross(ictr,1) + 1 
         if ((fanvec(k,i).eq.1).and.(fanvec(L,i).eq.1))  ! binary inner
c product
     x      cross(ictr,2) = cross(ictr,2) + 1 
        end do
       endif
      end do
      end do

       do i = 1, ictr
      write(6,33) cross(i,1), cross(i,2)
33     format(2i7)
       end do

      END
