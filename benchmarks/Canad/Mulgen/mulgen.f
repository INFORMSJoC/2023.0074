c   ****************************************************************
c   Network Problem Generator
c   by Dimitri P. Bertsekas from Linear Network Optimization
c   first version: september 1990
c   last revision: august 1993
c   by C. Beauchemin & B. Gendron
c   ****************************************************************

       subroutine mulgen

       include 'ran2.cmn'

       include 'dim.par'
       include 'mulgen.cmn'

       real ran2

       integer a,arc,gridarc
       logical mark(MAXCOMM, MAXNODES)
       logical omark(MAXNODES), dmark(MAXNODES)
       logical odmark(MAXNODES,MAXNODES)
       integer citer
       integer csum
       integer fsum
       real q


       n=dim1*dim2
       arc=0
       gridarc=0

       totsupply=0

c      a preliminary call to the random generator
       k= int(ran2(inseed))

       do 66 i=1,comm
          totsup(i)=minsup+int(ran2(inseed)*(maxsup-minsup))
          totsupply=totsupply+totsup(i)
  66   continue

c ****************** Grid arcs generation *********************

       if(grid) then

       do 10 j=1,dim2
         do 15 i=1,dim1
           node=i+(j-1)*dim1

c  Create horizontal arcs

	if (dim1 .gt. 1) then  
           if((i .gt. 1) .and. (i .lt. dim1)) then
              arc=arc+1
	      gridarc=gridarc+1
              startn(arc)=node
              endn(arc)=node+1
              c(arc)=maxfcost
              u(arc)=totsupply
	      do 311 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 311	      continue

              arc=arc+1
	      gridarc=gridarc+1
              startn(arc)=node
              endn(arc)=node-1
              c(arc)=maxfcost
              u(arc)=totsupply
	      do 312 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 312	      continue

           else
           if(i .eq. 1) then
              arc=arc+1
	      gridarc=gridarc+1
              startn(arc)=node
              endn(arc)=node+1
              c(arc)= maxfcost
              u(arc)= totsupply
	      do 313 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 313	      continue

 	   
	   if(dim1 .gt. 2) then
              arc=arc+1
	      gridarc=gridarc+1
              startn(arc)=node
              endn(arc)=node+(dim1-1)
              c(arc)=maxfcost
              u(arc)= totsupply
	      do 314 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 314	      continue

	   endif

           else

	   if (dim1 .gt. 2) then
              arc=arc+1 
	      gridarc=gridarc+1
              startn(arc)=node
              endn(arc)=node-(dim1-1)
              c(arc)=maxfcost
              u(arc)=totsupply
	      do 315 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 315	      continue

	   endif

              arc=arc+1
	      gridarc=gridarc+1
              startn(arc)=node
              endn(arc)=node-1
              c(arc)=maxfcost
              u(arc)=totsupply
	      do 316 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 316	      continue

           endif
           endif
      endif 

C   Create vertical arcs

	if(dim2 .gt. 1) then
           if((j .gt. 1) .and. (j .lt. dim2)) then
             arc=arc+1
	     gridarc=gridarc+1
             startn(arc)=node
             endn(arc)=node+dim1
             c(arc)=maxfcost
             u(arc)=totsupply
	      do 317 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 317	      continue

             arc=arc+1
	     gridarc=gridarc+1
             startn(arc)=node
             endn(arc)=node-dim1
             c(arc)=maxfcost
             u(arc)=totsupply
	      do 318 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 318	      continue

           else

           if(j .eq. 1) then
             arc=arc+1
	     gridarc=gridarc+1
             startn(arc)=node
             endn(arc)=node+dim1
             c(arc)=maxfcost
             u(arc)=totsupply
	      do 319 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 319	      continue


	   if(dim2 .gt. 2) then
             arc=arc+1
	     gridarc=gridarc+1
             startn(arc)=node
             endn(arc)=node+dim1*(dim2-1)
             c(arc)=maxfcost
             u(arc)=totsupply
	      do 321 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 321	      continue

	   endif

           else
             arc=arc+1
	     gridarc=gridarc+1
             startn(arc)=node
             endn(arc)=i
             c(arc)=maxfcost 
             u(arc)=totsupply
	      do 322 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 322	      continue



	   if(j .gt. 2) then
             arc=arc+1
	     gridarc=gridarc+1
             startn(arc)=node
             endn(arc)=node-dim1
             c(arc)=maxfcost
             u(arc)=totsupply
	      do 323 k=1,comm
		cost(k,arc)=maxcost
		cap(k,arc)=totsup(k)
 323	      continue

	   endif
           endif
           endif
        endif 

C   Create 2 arcs with random end node 

           do 8 k=1,2
	     citer=iter
             arc=arc+1
             startn(arc)=node
  7          endn(arc)=1+int(ran2(inseed)*n)
             if(endn(arc) .eq. node) goto 7

c   -- if noparallel parameter ---
             if (noparallel) then
                do 9 a=1,arc-1
                   if((startn(a) .eq. node) .and. 
     &                (endn(a) .eq. endn(arc))) then
			if(citer .gt. 0) then
			   citer=citer-1
			   goto 7
			else
 			   arc=arc-1
			   goto 8
			endif
		   endif
  9             continue
             endif
c   -- end of noparallel parameter --- 

             c(arc)=minfcost+ran2(inseed)*(maxfcost-minfcost)
             u(arc)=mincap+ran2(inseed)*(maxcap-mincap)
             do 33 kk=1,comm
	        cost(kk,arc)=mincost+ran2(inseed)*(maxcost-mincost)
		cap(kk,arc)=mincca+ran2(inseed)*(maxcca-mincca)
 33	     continue

  8        continue
 15      continue
 10    continue

	endif

c   ************** end; grid arc generation ***************
 
c   *********** Circle arc grid generation *****************
	if(cgrid) then

       do 117 node=1,n

              arc=arc+1
            gridarc=gridarc+1
              startn(arc)=node
              if(node .eq. n) then
              endn(arc)=1
            else
              endn(arc)=node+1
            endif
              c(arc)=maxfcost
              u(arc)=totsupply
            do 423 k=1,comm
              cost(k,arc)=maxcost
              cap(k,arc)=totsup(k)
 423        continue

 117	continue

	endif
c   ******************************************************** 

C   Create addarcs new arcs with random start and end node
        
       do 12 k=1,addarcs
         arc=arc+1
	 citer=iter
 130     node=1+int(ran2(inseed)*n)
         startn(arc)=node
 13      endn(arc)=1+int(ran2(inseed)*n)
         if (endn(arc) .eq. node) goto 13

c   -- if noparallel parameter ---
         if (noparallel) then
  
            do 14 a=1,arc-1
               if((startn(a) .eq. node) .and. 
     &            (endn(a) .eq. endn(arc))) then
		   if(citer .gt. 0) then
		      citer=citer-1
		      goto 130
		   else
                      arc=arc-1   
		      goto 12
		   endif
               endif  
 14         continue
         endif
c   -- end of noparallel parameter --- 

             c(arc)=minfcost+ran2(inseed)*(maxfcost-minfcost)
             u(arc)=mincap+ran2(inseed)*(maxcap-mincap)
             do 36 kk=1,comm
	        cost(kk,arc)=mincost+ran2(inseed)*(maxcost-mincost)
		cap(kk,arc)=mincca+ran2(inseed)*(maxcca-mincca)
 36	     continue

 12    continue

       na=arc

C   Adjust arcs for parameters nulfcst & topcap

c   -- if nulfcst parameter ----
	 if(nulfcst) then
	    do 11 km=na,na-(int(fcstprop*(na-gridarc)/100)-1),-1  
	       c(km)=0
 11	    continue
	 endif
c   -- end nulfcst parameter ---

c   -- if topcap parameter ----
	 if(topcap) then
	    do 34 lo=na,na-(int(capprop*(na-gridarc)/100)-1),-1  
	       u(lo)=totsupply
 34	    continue
	 endif
c   -- end nulcap parameter ---


c   -- if nulcom parameter ----
c      on choisit une proportion des arcs non-grid qu'on met a zero
	 if(nulcom) then
	   do 411 kk=1,int(nulcomprop*(na-gridarc)/100)
	   kkk=int(ran2(inseed)*comm)+1
	   arc=int(ran2(inseed)*(na-gridarc))+gridarc+1
	   cap(kkk,arc)=0
 411	   continue
	 endif
c   -- end nulcom parameter ---

c   -- if topcom parameter ----
c      on choisit une proportion des arcs non-grid qu'on met a maxcap
	 if(topcom) then
	   do 422 kk=1,int(topcomprop*(na-gridarc)/100)
	     kkk=int(ran2(inseed)*comm)+1
	     arc=int(ran2(inseed)*(na-gridarc))+gridarc+1
	     cap(kkk,arc)=maxcap
 422	   continue
	 endif
c   -- end topcom parameter ---


C   Generate the sinks and sources for each commodity

       do 167 i=1,comm 
       if(.not. dow) then
	if(.not. single) then
         sinks(i)=minsinks+ran2(inseed)*(maxsinks-minsinks)
         sources(i)=minsources+ran2(inseed)*(maxsources-minsources)
	else
	 ssinks=minsinks+ran2(inseed)*(maxsinks-minsinks)
	 ssources=minsources+ran2(inseed)*(maxsources-minsources)
	endif
       else
	sinks(i)=1
	sources(i)=1
       endif
 167   continue

c ----------- parametre single ----------------------------------
	if(single) then

	do 212 i=1,comm
	  do 213 nod=1,n
	  b(nod,i)=0
 213	  continue
	  omark(i)=.false.
	  dmark(i)=.false.
 212	continue

	do 214 i=1,ssources
 232	  node=1+int(ran2(inseed)*n)
	  if(omark(node)) then
	     goto 232
	  else
	     omark(node)=.true.
	  endif
 214	continue

	do 229 i=1,ssinks
 233	  node=1+int(ran2(inseed)*n)
	  if(dmark(node) .or. omark(node)) then
	     goto 233
	  else
	     dmark(node)=.true.
	  endif
 229	continue

	oindex=0
	dindex=0
	do 215 kk=1,comm
	   rest=0
	   supl=0
	   do 216 nod=1,n
	     if(omark(nod)) then
	       if(oindex .eq. 0) oindex=nod	
               temp=1+int(ran2(inseed)*(totsup(kk)/ssources))
	       b(kk,nod)=temp+rest
	       supl=supl+b(kk,nod)
	       rest=int(totsup(kk)/ssources)-temp
	     endif
 216	   continue
	   b(kk,oindex)=b(kk,oindex)+rest
	   supl=supl+rest

	   dem=0
	   rest=0

	   do 217 nod=1,n
	     if(dmark(nod)) then
	       if(dindex .eq. 0) dindex=nod	
               temp=1+int(ran2(inseed)*(totsup(kk)/ssinks))
	       b(kk,nod)=-(temp+rest)
	       dem=dem-b(kk,nod)
	       rest=int(totsup(kk)/ssources)-temp
	     endif
 217	   continue
	   b(kk,dindex)=b(kk,dindex)-rest
	   dem=dem+rest

C   Equalize supply and demand

       if (supl .gt. dem) then
          do 241 nod=1,n
          if(b(kk,nod) .lt. 0) then
            b(kk,nod)=b(kk,nod)-(supl-dem)
            goto 215
          endif
 241      continue
       else
       if(supl .lt. dem) then
          do 242 nod=1,n
          if(b(kk,nod) .gt. 0) then
            b(kk,nod)=b(kk,nod)+(dem-supl)
            goto 215
          endif
 242      continue
       endif
       endif
 
 215	continue
c ---------- fin du parametre single ---------------------------

	else

       if(dow)then
 
         do 333 i=1,n
         do 333 j=1,n
           odmark(i,j)=.false.
333      continue

         do 21 k=1,comm
         do 32 i=1,n
           b(k,i)=0
 32      continue
 230     onode=1+int(ran2(inseed)*n)
 23      dnode=1+int(ran2(inseed)*n)
         if(dnode.eq.onode)goto23
         if(odmark(onode,dnode))goto230
         odmark(onode,dnode)=.true.
         b(k,onode)=totsup(k)
         b(k,dnode)=-totsup(k)
 21      continue

       else

       do 200 kk=1,comm
         do 17 i=1,n
           b(kk,i)=0
           mark(kk,i)=.false.
 17      continue
 200   continue

       do 300 kk=1,comm
            
       supl=0
       rest=0

       do 20 i=1,sources(kk)
 18    node=1+int(ran2(inseed)*n)
       if(mark(kk,node)) then
         goto 18
       else
         mark(kk,node)=.true.
         if (i .eq. 1) index=node
         temp=1+int(ran2(inseed)*(totsup(kk)/sources(kk)))
         b(kk,node)=temp+rest
         supl=supl+b(kk,node)
         rest=int(totsup(kk)/sources(kk))-temp
       endif
 20    continue

C on donne a la premiere source le reste...
       b(kk,index)=b(kk,index)+rest
       supl=supl+rest

       dem=0
       rest=0

       do 25 i=1,sinks(kk)
 22    node=1+int(ran2(inseed)*n)
       if (mark(kk,node)) then
         goto 22
       else
         mark(kk,node)=.true.
         if(i .eq. 1) index=node
         temp=(1+int(ran2(inseed)*(totsup(kk)/sinks(kk))))
         b(kk,node)=-(temp+rest)
         dem=dem-b(kk,node)
         rest=int(totsup(kk)/sinks(kk))-temp
       endif
 25    continue

C on donne au premier puits le reste...
       b(kk,index)=b(kk,index)-rest
       dem=dem+rest

 
C   Equalize supply and demand

       if (supl .gt. dem) then
          do 30 node=1,n
          if(b(kk,node) .lt. 0) then
            b(kk,node)=b(kk,node)-(supl-dem)
            goto 300
          endif
 30       continue
       else
       if(supl .lt. dem) then
          do 35 node=1,n
          if(b(kk,node) .gt. 0) then
            b(kk,node)=b(kk,node)+(dem-supl)
            goto 300
          endif
 35       continue
       endif
       endif
 300   continue
       endif
       endif

       if(ctight)then
         csum=0
         do 111 i=1,na
           csum=csum+u(i)
111      continue
         q=float(na*totsupply)/float(csum) 

         do 112 i=1,na
           u(i)=min(int((q/xtight)*u(i))+1,totsupply)
112      continue

c        csum=0
c        do 113 i=1,na
c          csum=csum+u(i)
c13      continue
       endif

       if(fixvar)then
         fsum=0
         do 222 i=1,na
           fsum=fsum+c(i)
222      continue
         csum=0
         do 223 k=1,comm
         do 223 i=1,na
            csum=csum+cost(k,i)
223      continue
         csum=totsupply*(csum/comm)
         q=float(fsum)/float(csum)

         do 225 i=1,na
           c(i)=int((xfixvar/q)*c(i))+1
225      continue

c        fsum=0
c        do 226 i=1,na
c          fsum=fsum+c(i)
c26      continue
       endif

       return
       end
