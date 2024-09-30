	subroutine outdat

	include 'dim.par'
	include 'mulgen.cmn'

        integer arc,com,compt
        integer start,end,qty
	character chf

C   ***** Writing output file(s) *******************
C   With dow parameter -----------------------------

       if(dow) then
          print *,'with dow parameter'
	  open(3,file=dowfile)
          write(3,*), 'MULTIGEN.DAT:'
          write(3,1040) n,na,comm
        
	do 67 noeud=1,n 
          do 43 i=1,na
	   if (startn(i) .eq. noeud) then
           write(3,1050) startn(i),endn(i),cost(1,i),
     &  u(i), c(i), 1, i
	   endif
 43       continue
 67     continue
	  do 45 compt=1,comm
            do 44 i=1,n
              if(b(compt,i) .gt. 0) then
                start=i
		qty=b(compt,i)
	      endif
              if(b(compt,i) .lt. 0) end=i
 44         continue
          write(3,1060) start, end, qty
 45	  continue

	close(3)
	endif

 1040  format(3i8)
 1050  format(7i8)
c 1051  format(2i8,f12.4,4i8)
 1060  format(3i8)
C --- end of dow parameter --------------------------



C ------- RNET parameter ----------------------------
	 if(rnet) then
	  open(2,file=rnetfile)
          print *,'with rnet parameter'
          write(2,1041) n,na,sources(1)+sinks(1)

	  do 133 noeud=1,n
	    do 131 i=1,na
	      if (startn(i) .eq. noeud)
     &	        write(2,1042) startn(i),endn(i),cost(1,i),u(i)
 131	    continue
 133	  continue


	    do 132 noeud=1,n
	      if (b(1,noeud) .ne. 0) write(2,1043) noeud,b(1,noeud)
 132	  continue
	 endif

 1041	format(3i8)
 1042	format(4i8)
 1043	format(2i8)
C ---- fin du parametre RNET ------------------------

	if(mrnet) then
	  print *,'Multi Rnet'
	  do 287 com=1,comm
	     rnetfile=outfile
	     call int2char(com,chf,.false.)
	     call append(rnetfile,chf//'.rnet')
	     open(2,file=rnetfile)
	     write(2,1041) n,na,sources(com)+sinks(com)
	     do 288 node=1,n
		do 289 arc=1,na
		  if(startn(arc) .eq. node) then
		    if(comm .eq. 1) then
                    write(2,1042) startn(arc),endn(arc),cost(com,arc),
     &                            u(arc)
		    else
                    write(2,1042) startn(arc),endn(arc),cost(com,arc),
     &                            cap(com,arc)
	            endif
                  endif
 289		continue
 288	     continue

	    do 290 node=1,n
              if (b(com,node).ne.0) write(2,1043) node,b(com,node)
 290	    continue
 287	  continue
	endif

C ---------------------------------------------------------------
           

       print *,'*************************************'
       print *,'Writing output to file ',outfile

	if(nocca) then
	   do 431 arc=1,na
	      do 432 com=1,comm
		cap(com,arc)=min(u(arc),totsup(com))
 432	      continue
 431	   continue
	endif

	if(dow) then
	   do 434 arc=1,na
	      do 435 com=1,comm
		 cap(com,arc)=min(u(arc),totsup(com))
		 cost(com,arc)=cost(1,arc)
 435	      continue
 434	   continue	
	endif


       if(fich4) then
c 1)
	  ccabin=.true.
	  capbin=.true.
	  call traitenom(1)
c	  call lp()
          call mps
c 2)
	  ccabin=.true.
	  capbin=.false.
	  call traitenom(2)
c	  call lp()
          call mps
c 3)
	  ccabin=.false.
	  capbin=.true.
	  call traitenom(3)
c	  call lp()
          call mps
c 4)
	  ccabin=.false.
	  capbin=.false.
	  call traitenom(4)
c	  call lp()
          call mps

       else

       if(mpsflg) call mps()
       if(cplex) call lp()

       endif

c   -- if LP1 option --
	if(LP1) then
	   print *,'Option LP1'
	   capbin=.false.
	   ccabin=.false.
c	   do 455 com=1,comm
c		do 455 arc=1,na
c		   cost(com,arc)=cost(com,arc)+
c     &			(float(c(arc))/float(u(arc)))
c 455	   continue 
	mpsfile=outfile
	call append(mpsfile,'LP1.mps')
	call mps()
	endif
c   -- end LP1 --

c   Writing the standard output file %%%%%%%%%%%%%%%%%%%%

       if(.not. nosfil) then

       open(13,file=outfile)
       rewind(13)

C   Write number of nodes, arcs and commodities

       write(13,1010) n,na,comm


C   Write start,end,cost,capacity of each arc
C   and cost capacity for each commodity
c   Let's count the number of commodities per arc

       do 670 noeud=1,n
       do 40 i=1,na
	 if(startn(i) .eq. noeud) then

         nbcomm=0
         do 122 j=1,comm
            if(cap(j,i) .ne. 0) nbcomm=nbcomm+1 
 122     continue

         write(13,1020) startn(i),endn(i),c(i),u(i),nbcomm
c AJOUTE PAR BG
         if(comm.eq.1)then
           write(13,1010) 1,cost(1,i),u(i)
         else

         do 41 j=1,comm
           if (cap(j,i) .ne. 0)
     &         write(13,1010) j,cost(j,i),cap(j,i)
 41      continue

         endif

	 endif
 40    continue
 670   continue

C   Write supply of nodes

       do 51 compt=1,comm
       do 50 i=1,n
       if(b(compt,i) .ne. 0) write(13,1010) compt,i,b(compt,i)
 50    continue
 51    continue


       endfile(13)
       rewind(13)
       print *,'End of writing.'

       endif

 1000  format(2i8)
 1010  format(3i8)
c 1011  format(i8,f8.2,i8)
 1020  format(5i8)
 1030  format(a23,a8)
 1039  format(i8)

       return
       end
