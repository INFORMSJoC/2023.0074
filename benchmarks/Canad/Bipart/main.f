      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'

      call infile
      call inipar
      call getpar
      write(udataf,1000)commod
1000  format('before calling bipart',i5)
      call bipart
      call outdat

      stop
      end

      subroutine infile
 
      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'

      character*40 fparam

      if(iargc().lt.1)then
        print*, 'Usage: bipart f (''f'' is a parameter file)'
        stop
      endif
      call getarg(1,fparam)
      open(uparam,file=fparam)

      return
      end

      subroutine inipar

      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'
 
      origin=10
      destin=10
      commod=1
      odarcs=10
      minsup=100
      maxsup=100
      mincap=1
      maxcap=100
      minfct=1
      maxfct=100
      mintct=1 
      maxtct=100
      iseed=101
      ctight=.false.
      xtight=2.
      fixvar=.false.
      ifixvar=2
      caploc=.false.

      return
      end

      subroutine getpar

      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'
 
      character*6 par
      character*80 name

  10  format(a6)
  20  format(7x,a80)
  30  format(7x,i5)
  40  format(7x,f6.2)

      read(uparam,10,end=77)par
      if(par.eq.'datafn')then
        backspace(1)
        read(uparam,20)name
        close(udataf)
        open(udataf,file=name)
      else
        goto 77
      endif

   1  read(uparam,10,end=99)par

      if(par.eq.'origin')then
        backspace(1)
        read(uparam,30)origin
        if(origin.lt.1)then
          print*,'ERROR: ORIGIN LESS THAN 1'
          stop
        endif
        goto 1
      endif

      if(par.eq.'destin')then
        backspace(1)
        read(uparam,30)destin
        if(destin.lt.1)then
          print*,'ERROR: DESTIN LESS THAN 1'
          stop
        endif
        goto 1
      endif

      if(par.eq.'odarcs')then
        backspace(1)
        read(uparam,30)odarcs
        if(odarcs.lt.1)then
          print*,'ERROR: ODARCS LESS THAN 1'
          stop
        endif
        goto 1
      endif

      if(par.eq.'commod')then
        backspace(1)
        read(uparam,30)commod
        if(commod.lt.1)then 
          print*,'ERROR: COMMOD LESS THAN 1'
          stop
        endif
        goto 1
      endif 

      if(par.eq.'minsup')then
        backspace(1)
        read(uparam,30)minsup
        if(minsup.lt.1)then
          print*,'ERROR: MINSUP LESS THAN 1'
          stop
        endif
        goto 1
      endif

      if(par.eq.'maxsup')then
        backspace(1)  
        read(uparam,30)maxsup
        if(maxsup.lt.1)then 
          print*,'ERROR: MAXSUP LESS THAN 1'
          stop
        endif
        goto 1
      endif

      if(par.eq.'mincap')then
        backspace(1)
        read(uparam,30)mincap
        if(mincap.lt.1)then
          print*,'ERROR: MINCAP LESS THAN 1'
          stop
        endif
        goto 1
      endif

      if(par.eq.'maxcap')then 
        backspace(1)  
        read(uparam,30)maxcap
        if(maxcap.lt.1)then 
          print*,'ERROR: MAXCAP LESS THAN 1'
          stop 
        endif 
        goto 1 
      endif  
 
      if(par.eq.'minfct')then
        backspace(1)
        read(uparam,30)minfct
        if(minfct.lt.0)then
          print*,'ERROR: MINFCT LESS THAN 0'
          stop
        endif
        goto 1
      endif

      if(par.eq.'maxfct')then
        backspace(1)  
        read(uparam,30)maxfct
        if(maxfct.lt.0)then 
          print*,'ERROR: MAXFCT LESS THAN 0' 
          stop 
        endif 
        goto 1
      endif

      if(par.eq.'mintct')then
        backspace(1)  
        read(uparam,30)mintct
        if(mintct.lt.0)then 
          print*,'ERROR: MINTCT LESS THAN 0' 
          stop 
        endif 
        goto 1
      endif 
 
      if(par.eq.'maxtct')then 
        backspace(1)   
        read(uparam,30)maxtct 
        if(maxtct.lt.0)then  
          print*,'ERROR: MAXTCT LESS THAN 0'  
          stop
        endif
        goto 1
      endif

      if(par.eq.'iseed ')then 
        backspace(1)   
        read(uparam,30)iseed 
        goto 1
      endif

      if(par.eq.'ctight')then
        ctight=.true.
        backspace(1)   
        read(uparam,40)xtight
        goto 1
      endif
 
      if(par.eq.'fixvar')then
        fixvar=.true.
        backspace(1)
        read(uparam,40)xfixvar
        goto 1
      endif

      if(par.eq.'caploc')then
        caploc=.true.
        goto 1
      endif

      print*,'UNKNOWN PARAMETER'
      stop
 
  77  print*,'DATA FILE NOT SPECIFIED'
      stop

  99  continue

      if(origin+destin.gt.maxnod)then
        print*,'NUMBER OF NODES LARGER THAN MAXIMUM: ',maxnod
        stop
      endif

      if(origin*odarcs.gt.maxarc)then
        print*,'NUMBER OF ARCS LARGER THAN MAXIMUM: ',maxarc
        stop
      endif

      if(minsup.gt.maxsup)then
        print*,'ERROR: MINSUP GREATER THAN MAXSUP'
        stop
      endif
 
      if(mincap.gt.maxcap)then
        print*,'ERROR: MINCAP GREATER THAN MAXCAP'
        stop
      endif

      if(minfct.gt.maxfct)then 
        print*,'ERROR: MINFCT GREATER THAN MAXFCT'
        stop 
      endif

      if(mintct.gt.maxtct)then 
        print*,'ERROR: MINTCT GREATER THAN MAXTCT'
        stop
      endif 

      return
      end

      block data

      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'

      data uparam,udataf /1,2/
      data ldataf /.false./
 
      end
