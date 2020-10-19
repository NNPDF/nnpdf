      program split
      implicit none

      integer idat,i,mbin,irun
      double precision yin,y(132),y2(132)
      double precision xsec_nor(132),xsec_abs(132),mll(6)

      open(unit=34,file="CMS-DY2D11.data",status="old")
      open(unit=35,file="CMS-DY2D11-ABS.data",status="unknown")
      open(unit=36,file="CMS-DY2D11-NOR.data",status="unknown")

      yin = -0.05
      do idat = 1,24
         yin = yin + 0.1
         y(idat) = yin 
      enddo
      yin = -0.1
      do idat = 1,12
         yin = yin + 0.2
         y2(idat) = yin 
      enddo

      mll(1) = 25d0
      mll(2) = 37.5d0
      mll(3) = 52.5d0
      mll(4) = 90d0
      mll(5) = 160d0
      mll(6) = 850d0

      do idat = 1,132

         read(34,*) i,xsec_nor(idat),xsec_abs(idat)
         
         if(idat.le.24) then
            mbin = 1
            irun = idat-24*(mbin-1)
            write(35,301)idat,y(irun),mll(mbin),xsec_abs(idat)
            write(36,302)idat,y(irun),mll(mbin),xsec_nor(idat)
         elseif(idat.ge.25.and.idat.le.48) then
            mbin = 2
            irun = idat-24*(mbin-1)
            write(35,301)idat,y(irun),mll(mbin),xsec_abs(idat)
            write(36,302)idat,y(irun),mll(mbin),xsec_nor(idat)
         elseif(idat.ge.49.and.idat.le.72) then
            mbin = 3
            irun = idat-24*(mbin-1)
            write(35,301)idat,y(irun),mll(mbin),xsec_abs(idat)
            write(36,302)idat,y(irun),mll(mbin),xsec_nor(idat)
         elseif(idat.ge.73.and.idat.le.96) then
            mbin = 4
            irun = idat-24*(mbin-1)
            write(35,301)idat,y(irun),mll(mbin),xsec_abs(idat)
            write(36,302)idat,y(irun),mll(mbin),xsec_nor(idat)
         elseif(idat.ge.97.and.idat.le.120) then
            mbin = 5
            irun = idat-24*(mbin-1)
            write(35,301)idat,y(irun),mll(mbin),xsec_abs(idat)
            write(36,302)idat,y(irun),mll(mbin),xsec_nor(idat)
         elseif(idat.ge.121.and.idat.le.132) then
            mbin = 6
            irun = idat-24*(mbin-1)
            write(35,301)idat,y2(irun),mll(mbin),xsec_abs(idat)
            write(36,302)idat,y2(irun),mll(mbin),xsec_nor(idat)
         endif

      enddo

      close(34)
      close(35)
      close(36)

 301  format(I4,2(3X,F7.3),3X,F13.3)
 302  format(I4,2(3X,F7.3),3X,ES15.8)

      stop
      end
