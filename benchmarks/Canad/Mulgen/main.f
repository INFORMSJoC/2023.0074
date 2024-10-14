      include 'dim.par'
      include 'mulgen.cmn'

      call infile
      call inipar
      call getpar
      call mulgen
      call outdat

      stop
      end 

      subroutine infile

      include 'dim.par'
      include 'mulgen.cmn'

      character*40 fparam

      if(iargc().lt.1)then
        print*, 'Usage: mulgen f (''f'' is a parameter file)'
        stop
      endif
      call getarg(1,fparam)
      open(uparam,file=fparam)

      return
      end

      subroutine inipar

      include 'dim.par'
      include 'mulgen.cmn'

        inseed=101
        dim1=2
        dim2=2
        comm=1
        addarcs=0
        minsources=1
        maxsources=1
        minsinks=1
        maxsinks=1
        minsup=0
        maxsup=10000
        mincost=0
        maxcost=10000
        mincap=1
        maxcap=10000
        minfcost=0
        maxfcost=10000
        mincca=0
        maxcca=10000
        outfile='outfile'

        dow=.false.
        noparallel=.false.
        iter=0
	ctight=.false.
	itight=2
	fixvar=.false.
	ifixvar=2
        mpsflg=.false.
        nulfcst=.false.
        topcap=.false.
        nulcom=.false.
        topcom=.false.
        single=.false.
        rnet=.false.
        capbin=.false.
        ccabin=.false.
        nocca=.false.
        fich4=.false.
        mrnet=.false.
        grid=.false.
        nosfil=.false.
        LP1=.false.
 
      return
      end

      subroutine getpar

      include 'dim.par'
      include 'mulgen.cmn'

 11    read(uparam,2001),var

       if(var.eq.'inseed')then
          backspace(1)
          read(uparam,2000) var,inseed
          goto 11
       endif

       if (var .eq. 'scalex') then
          backspace(1)
          read(uparam,2000) var,dim1
          if (dim1 .lt. 1) then
              print *,'scalex must be greater than 0'
              print *,'EXECUTION ABORTED'
              stop
          endif
          goto 11
       endif

       if (var .eq. 'scaley') then
          backspace(1)
          read(uparam,2000) var,dim2
          if (dim2 .lt. 1) then
              print *,'scaley must be greater than 0'
              print *,'EXECUTION ABORTED'
              stop
          endif
          goto 11
        endif
 
       if (var .eq.'commod') then
          backspace(1)
          read(uparam,2000) var,comm
          if (comm .gt. MAXCOMM) then
              print *,'maximum commodities = ',MAXCOMM
              print *,'EXECUTION ABORTED'
              stop
          endif
          goto 11
        endif
 
       if (var .eq.'addarc') then
          backspace(1)
          read(uparam,2000) var,addarcs
          if (addarcs .gt. MAXARCS) then
              print *,'maximum number of arcs = ',MAXARCS
              print *,'EXECUTION ABORTED'
              stop
          endif
           goto 11
        endif
 
       if (var .eq. 'minsrc') then
          backspace(1)
          read(uparam,2000) var,minsources
          if(minsources .le. 0) then
             print *,'minsrc must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif
 
       if (var .eq. 'maxsrc') then
          backspace(1)
          read(uparam,2000) var,maxsources
          if(maxsources .le. 0) then
             print *,'maxsrc must be greater than zero'
             print *,'EXECUTION ABORTED'
            stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'minsnk') then
          backspace(1)
          read(uparam,2000) var,minsinks
          if(minsinks .le. 0) then
             print *,'minsnk must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'maxsnk') then
          backspace(1)
          read(uparam,2000) var,maxsinks
          if(maxsinks .le. 0) then
             print *,'maxsnk must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'minsup') then
          backspace(1)
          read(uparam,2000) var, minsup
          if(minsup .le. 0) then
             print *,'minsup must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'maxsup') then
          backspace(1)
          read(uparam,2000) var, maxsup
          if(maxsup .le. 0) then
             print *,'maxsup must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'mincst') then
          backspace(1)
          read(uparam,2000) var,mincost
          if(mincost .lt. 0) then
             print *,'mincst must be greater or equal to zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'maxcst') then
          backspace(1)
          read(uparam,2000) var,maxcost
          if(maxcost .lt. 0) then
             print *,'maxcst must be greater or equal to zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'mincap') then
          backspace(1)
          read(uparam,2000) var,mincap
          if(mincap .le. 0) then
             print *,'mincap must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'maxcap') then
          backspace(1)
          read(uparam,2000) var,maxcap
          if(maxcap .le. 0) then
             print *,'maxcap must be greater than zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'minfct') then
          backspace(1)
          read(uparam,2000) var,minfcost
          if(minfcost .lt. 0) then
             print *,'minfct must be greater or equal to zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'maxfct') then
          backspace(1)
          read(uparam,2000) var,maxfcost
          if(maxfcost .lt. 0) then
             print *,'maxfct must be greater or equal to zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'mincca') then
          backspace(1)
          read(uparam,2000) var,mincca
          if(mincca .lt. 0) then
             print *,'mincca must be greater or equal to zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'maxcca') then
          backspace(1)
          read(uparam,2000) var,maxcca
          if(maxcca .lt. 0) then
             print *,'maxcca must be greater or equal to zero'
             print *,'EXECUTION ABORTED'
             stop
          else
             goto 11
          endif
       endif

       if (var .eq. 'outfil') then
          backspace(1)
          read(uparam,2002) var,outfile
          goto 11
       endif


       if (var .eq. 'end') then

          call traitenom(0)

          if((dim1*dim2) .gt. MAXNODES) then
              print *,'dim1*dim2 greater than ',MAXNODES
              print *,'EXECUTION ABORTED'
              stop
          endif

c  On doit verifier que les capacites partielles sont bien
c  inferieures a la capacite par arc ou au supply
c  ainsi que la sommes des capacites partielles soit
c  plus grande ou egale a la capacite de l'arc

        if((comm.gt.1).and.(.not.dow).and.(.not.nocca))then
c  -1-
          if(maxcca .gt. minsup) maxcca=minsup
          if(maxcca .gt. mincap) maxcca=mincap
          if(mincca .gt. maxcca) mincca=maxcca

c  -2-
          if(maxcap .gt. comm*mincca) maxcap=comm*mincca
          if(mincap .gt. maxcap) mincap=maxcap
c  -3-
          if(maxcca .gt. mincap) maxcca=mincap
        endif

          if(minsources .gt. maxsources) then
              print *,'minsrc greater than maxsrc!'
              print *,'EXECUTION ABORTED'
              stop
          endif

          if(minsinks .gt. maxsinks) then
              print *,'minsnk greater than maxsnk!'
              print *,'EXECUTION ABORTED'
              stop
          endif

          if(minsup .gt. maxsup) then
              print *,'minsup greater than maxsup!'
              print *,'EXECUTION ABORTED'
              stop
          endif

          if(mincost .gt. maxcost) then
              print *,'mincst greater than maxcst!'
              print *,'EXECUTION ABORTED'
              stop
          endif

          if(mincap .gt. maxcap) then
              print *,'mincap greater than maxcap!'
              print *,'EXECUTION ABORTED'
              stop
          endif

          if(minfct .gt. maxfct) then
              print *,'minfct greater than maxfct!'
              print *,'EXECUTION ABORTED'
              stop
          endif

          if(mincca .gt. maxcca) then
              print *,'mincca greater than maxcca!'
              print *,'EXECUTION ABORTED'
              stop
          endif
 
       if ((maxsinks+maxsources) .gt. (dim1*dim2)) then
          print *,'Sinks + sources > number of nodes'
          print *,'Check parameters maxsrc and maxsnk.'
          print *,'EXECUTION ABORTED'
          stop
       endif

        if ((dim2 .eq. 1) .and. (dim1 .eq. 1)) then
          print *,'Cannot create 1x1 network'
          print *,'EXECUTION ABORTED'
          stop
        endif

        if(rnet) comm=1

          return
       endif

C   Parameters...

       if (var .eq. 'dow') then
          dow=.true.
          goto 11
       endif

       if (var .eq. 'nulcst') then
          nulfcst=.true.
          backspace(1)
          read(uparam,2000) var, fcstprop
          if ((fcstprop .gt. 100) .or. (fcstprop .lt. 0)) then
                print *,'Fix cost proportion invalid [0,100]'
                print *,'EXECUTION ABORTED'
                stop
          endif
          write(*,3000)'with nulcst parameter at ',fcstprop,' %'
          goto 11
       endif

       if (var .eq. 'topcap') then
          topcap=.true.
          backspace(1)
          read(uparam,2000) var, capprop
          if ((capprop .gt. 100) .or. (capprop .lt. 0)) then
                print *,'Top capacity proportion invalid [0,100]'
                print *,'EXECUTION ABORTED'
                stop
          endif
          write(*,3000)'with topcap parameter at ',capprop,' %'
          goto 11
       endif

       if (var .eq. 'nulcom') then
          nulcom=.true.
          backspace(1)
          read(uparam,2000) var, nulcomprop
          if ((nulcomprop .gt. 100) .or. (nulcomprop .lt. 0)) then
                print *,'Nul commodity proportion invalid [0,100]'
                print *,'EXECUTION ABORTED'
                stop
          endif
          write(*,3000)'with nulcom parameter at ',nulcomprop,' %'
        goto 11
        endif

       if (var .eq. 'topcom') then
          topcom=.true.
          backspace(1)
          read(uparam,2000) var, topcomprop
          if ((topcomprop .gt. 100) .or. (topcomprop .lt. 0)) then
                print *,'Top commodity proportion invalid [0,100]'
                print *,'EXECUTION ABORTED'
                stop
          endif
          write(*,3000)'with topcom parameter at ',topcomprop,' %'
        goto 11
        endif

       if (var .eq. 'noprll') then
          noparallel=.true.
          backspace(1)
          read(uparam,2000) var, iter
         goto 11
       endif

       if (var .eq. 'ctight') then
          ctight=.true.
          backspace(1)
          read(uparam,2004) var,xtight 
         goto 11
       endif

       if (var .eq. 'fixvar') then
          fixvar=.true.
          backspace(1)
          read(uparam,2004) var,xfixvar 
         goto 11
       endif

       if (var .eq. 'mps') then
          mpsflg=.true.
c         backspace(1)
c         read(uparam,2003) var, mpsfile
c         print *,'with mps output in file :',mpsfile
          goto 11
       endif

       if (var .eq. 'capbin') then
          capbin=.true.
          print *,'with capbin'
         goto 11
       endif

       if (var .eq. 'ccabin') then
          ccabin=.true.
          print *,'with ccabin'
         goto 11
       endif

       if (var .eq. 'rnet') then
          rnet=.true.
         goto 11
       endif

       if (var .eq. 'single') then
          single=.true.
          print *,'with single od'
         goto 11
       endif

       if (var .eq. 'cplex') then
          cplex=.true.
         goto 11
       endif

       if (var .eq. 'nocca') then
          nocca=.true.
          print *,'no commodity capacity'
         goto 11
       endif

       if (var .eq. 'grdarc') then
          backspace(1)
          read(uparam,2000) var, itmp
          if(itmp .eq. 0)  then
                grid=.false.
                cgrid=.false.
                goto 11
          endif
          if(itmp .eq. 1)  then
                grid=.true.
                cgrid=.false.
                goto 11
          endif
          if(itmp .eq. 2)  then
                cgrid=.true.
                grid=.false.
                goto 11
          endif
c         print *,'Wrong grdarc value [0-1-2]'
c         stop
       endif

       if (var .eq. 'nosfil') then
          nosfil=.true.
         goto 11
       endif

       if (var .eq. '4fich') then
          fich4=.true.
         goto 11
       endif

 
       if(var .eq. 'mrnet') then
          mrnet=.true.
          goto 11
       endif

       if(var .eq. 'lp1') then
          LP1=.true.
          goto 11
       endif

c A. Frangioni 23/01/2013
c the following format statements used to have ' ' in place
c of the 1x, but this is no longer properly recognized by
c current fortran compilers, it seems

 2000  format(a6,1x,i8)
 2001  format(a6)
 2002  format(a6,1x,a40)
 2003  format(a6,1x,a8)
 2004  format(a6,1x,f6.2)
 3000  format(a25,i2,a2)

       print *,'ERROR: Bad command or no end statement'
       print *,'variable name: ',var
       stop
       end

      block data

      include 'dim.par'
      include 'mulgen.cmn'

      data uparam,udataf /1,2/
c      data ldataf /.false./

      end
