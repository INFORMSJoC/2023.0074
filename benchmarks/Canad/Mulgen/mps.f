	subroutine mps()

	include 'dim.par'
        include 'mulgen.cmn'

        integer arc,com,compt
	integer comfield, nodfield, arcfield
        real xxx

	comfield=int(log10(real(comm)))+1 
	nodfield=int(log10(real(n)))+1
	arcfield=int(log10(real(na)))+1

	if((comfield+nodfield) .gt. 7) then
	   print *,'COMMODITIES+NODES > 7'
	   print *,'FATAL ERROR: EXECUTION ABORTED'
	   stop
	endif

	if((comfield+arcfield) .gt. 7) then
	   print *,'COMMODITIES+ARCS > 7'
	   print *,'FATAL ERROR: EXECUTION ABORTED'
	   stop
	endif

	print *,'n na comm',n,na,comm
	print *,'nod arc com',nodfield,arcfield,comfield

c	write(*,*) 'mpsfile=',mpsfile
	open(9,file=mpsfile)
c	write(9,*)'n na comm',n,na,comm	
	write(9,*) '* Data generated by Gridgen (c)1993'
	write(9,3000) outfile

c   ***************** ECRITURE DES RANGEES ******************

	write(9,*) 'ROWS'

c   ** On alloue une variable a l'objectif ***

	write(9,*) ' N  cost'

c   ** Equations de continuite ( supply ) ***

	do 2200 noeud=1,n
	  do 2205 com=1,comm
             write(9,3002) com*10**nodfield+noeud
 2205	  continue
 2200	continue

 3000	format('NAME',t15,a8,'.mps')
 3002	format(' E  ','N',i7)

c   ** Equations de capacite de l'arc ***

	if((comm .gt. 1) .or. 
     &    ((comm .eq. 1) .and. (capbin .or. ccabin))) then
	do 2210 i=1,na
             write(9,3102) i
 2210	continue
	endif

 3102	format(' L  ','A',i7)

c   ** Equations de capacite par produit ***

	if((comm .gt. 1) .and. (ccabin)) then
	do 2220 arc=1,na
	  do 2225 com=1,comm
            write(9,3202) com*10**arcfield+arc
 2225     continue
 2220	continue
	endif

 3202	format(' L  ','P',i7)

c   *************** ECRITURE DES COLONNES ******************

	write(9,*) 'COLUMNS'

c   ** Ecriture pour les variables en x ********

c       * On pose l'equation avec l'objectif *

	do 2300 arc=1,na
	  do 2310 com=1,comm
           if(lp1)then
              if(dow)then
                xxx=cost(1,arc)+float(c(arc))/float(u(arc))
                write(9,4101)com*10**arcfield+arc,xxx
              else
                xxx=cost(com,arc)+float(c(arc))/float(u(arc))
                write(9,4101)com*10**arcfield+arc,xxx
              endif

           else

	   if(dow) then
               write(9,4100) com*10**arcfield+arc, cost(1,arc)
	   else
               write(9,4100) com*10**arcfield+arc, cost(com,arc)
	   endif

           endif

c       * On ecrit la contrainte de continuite *

	   do 2330 node=1,n
	      if(startn(arc) .eq. node) then
              write(9,4203) com*10**arcfield+arc,com*10**nodfield+node,1
	      else
	      if(endn(arc) .eq. node) then
             write(9,4203) com*10**arcfield+arc,com*10**nodfield+node,-1
	      endif
	      endif
 2330      continue

c        * On ecrit l'equation de capacite de l'arc *

	if((comm .gt. 1) .or. 
     &    ((comm .eq. 1) .and. (capbin .or. ccabin))) then
             write(9,4300) com*10**arcfield+arc,arc,1
	endif

c        * On ecrit l'equation de capacite par produit *

	if((comm .gt. 1) .and. ccabin) then
             write(9,4400) com*10**arcfield+arc,com*10**arcfield+arc,1
        endif

 2310	  continue
 2300	continue

 4100	format(t5,'x',i7,t15,'cost    ','  ',i12)
 4101	format(t5,'x',i7,t15,'cost    ','  ',f12.4)
 4203	format(t5,'x',i7,t15,'N',i7,t25,i12)
 4300	format(t5,'x',i7,t15,'A',i7,t25,i12)
 4400	format(t5,'x',i7,t15,'P',i7,t25,i12)

c   *** Ecriture pour les variables y **************

c	write(9,1234)

        do 3457 arc=1,na

c     * Objectif (cout fixe) *

	if(capbin .or. ccabin) then
	write(9,4103) arc,c(arc)
	endif

c     * Capacite par arc *
	if((comm .eq. 1) .and. (capbin .or. ccabin) .or.
     &     (comm .gt. 1) .and. capbin) then
	   if(u(arc) .ne. 0) then
	     write(9,4310) arc,arc,-u(arc)
	   endif
	endif

c     * Capacite par produit *
	if((comm .gt. 1) .and. ccabin) then
	    do 852 com=1,comm
	      write(9,4311) arc,com*10**arcfield+arc,-cap(com,arc)
 852	    continue
	endif

 3457   continue

c	write(9,1235)

c	endif

 4103	format(t5,'y',i7,t15,'cost    ','  ',i12)
 4310   format('    y',i7,t15,'A',i7,t25,i12)
 4311   format(t5,'y',i7,t15,'P',i7,t25,i12)
c 1234   format(t5,'MARK0000',2x,'''MARKER''',t40,'''INTORG''') 
c 1235   format(t5,'MARK0001',2x,'''MARKER''',t40,'''INTEND''') 


c   ***************** RIGHT-HAND SIDE *********************

	write(9,*) 'RHS'

c   ** Equations de continuite *** 

	compt=1
	do 2500 node=1,n
	  do 2510 com=1,comm
           if(b(com,node) .ne. 0) then
             write(9,4500) com*10**nodfield+node,b(com,node)
           endif
 2510     continue 
 2500	continue

 4500	format('    rhs',t15,'N',i7,t25,i12)

c   ** Equations de capacite par arc *** 

	if((comm .gt. 1) .and. (.not. capbin)) then
	do 2511 arc=1,na
	   write(9,4501) arc,u(arc)
 2511	continue
	endif

 4501	format('    rhs',t15,'A',i7,t25,i12)
	

c   ** Equations de capacite par produit ***
c	(pas rapport!)

c	if(strong) then
c	do 2600 arc=1,na
c	  do 2610 com=1,comm
c	    if(cap(com,arc) .ne. 0) then
c	      write(9,4600) com*10**arcfield+arc,cap(com,arc)
c	    endif
c 2610	  continue
c 2600	continue
c	endif

c 4600	format('    rhs',t15,'P',i7,t25,i12)

c	write(9,*) 'RANGES'

c   **************** ECRITURE DES BOUNDS ********************

	write(9,*) 'BOUNDS'
	
	if((comm .eq. 1) .and. .not.(capbin .or. ccabin)) then
	do 2692 arc=1,na
	      write(9,4655) 10**arcfield+arc,u(arc)
 2692	continue
	endif

	if((comm .gt. 1) .and. (.not. ccabin) .and.
c     &  (.not. dow) .and. (.not. nocca)) then
     &  (.not. dow)) then
	do 2690 arc=1,na
	   do 2695 com=1,comm
	      write(9,4655) com*10**arcfield+arc,cap(com,arc)
 2695	   continue
 2690	continue
	endif

	if(ccabin .or. capbin) then
	do 2700 arc=1,na
	   write(9,4650) arc
 2700	continue
	endif

 4650	format(' BV BOUND',t15,'y',i7)
 4655	format(' UP BOUND',t15,'x',i7,t25,i12)
	
	
	write(9,*) 'ENDATA'

        end
