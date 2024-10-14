	subroutine lp()

	include 'dim.par'
	include 'mulgen.cmn'

	character*80 ligne
	character lignevide
	character*8 chf
	character*20 block
	character*7  va
	logical sign,nosign

	integer comfield, nodfield, arcfield
        integer arc,com,compt

	comfield=int(log10(real(comm)))+1 
	nodfield=int(log10(real(n)))+1
	arcfield=int(log10(real(na)))+1

	sign=.true.
	nosign=.false.

	lignevide=' '
	block=' '
	open(10,file=lpfile)
	write(10,*) '\\Problem name: '//lpfile
	write(10,*)

c ********* Fonction objective ********************

	write(10,*) 'Minimize'

	ligne=' cost: '	
	compt=1
	do 10 arc=1,na
	  do 20 com=1,comm
	   if(dow) then
	     call int2char(cost(1,arc),chf,sign)
	   else
	     call int2char(cost(com,arc),chf,sign)
	   endif
	   call int2char(com*10**arcfield+arc,va,nosign) 
	   block=chf//' x'//va
	   call append(ligne,block)
	   if(mod(compt,4) .eq. 0) then
		write(10,*) ligne
		ligne=lignevide
	   endif
	   compt=compt+1
 20	  continue
 10	continue

	if(capbin .or. ccabin) then
	  do 25 arc=1,na
	    call int2char(c(arc),chf,sign)
	    call int2char(arc,va,nosign)
	    block=chf//' y'//va
	    call append(ligne,block)
	   if(mod(compt,4) .eq. 0) then
		write(10,*) ligne
		ligne=lignevide
	   endif
	   compt=compt+1
 25	  continue
	endif
	write(10,*) ligne

c ********* Contraintes ****************************

	write(10,*) 'Subject To'

c	--- contraintes de continuite ---
	compt=1
	do 30 noeud=1,n
	  do 40 com=1,comm
	     call int2char(com*10**nodfield+noeud,va,nosign)
	     ligne=' N'//va//':'
	     do 50 arc=1,na
		if(startn(arc) .eq. noeud) then
		  call int2char(com*10**arcfield+arc,va,nosign)
		  call append(ligne,' +x'//va)
		  compt=compt+1
		else
		if(endn(arc) .eq. noeud) then
		  call int2char(com*10**arcfield+arc,va,nosign)
		  call append(ligne,' -x'//va)
		  compt=compt+1
		endif
		endif
	   	if(mod(compt,8) .eq. 0) then
		  write(10,*) ligne
		  ligne=lignevide
		  compt=compt+1
		endif
 50	     continue
	     call int2char(b(com,noeud),va,sign)
	     call append(ligne,' = '//va)
	     write(10,*) ligne
	     compt=1
 40	  continue
 30	continue

c	--- contraintes de capacite de l'arc ---

	if((comm .gt. 1) .or. ((comm .eq. 1) .and.
     &  (capbin .or. ccabin))) then
	do 60 arc=1,na
	   call int2char(arc,va,nosign)
	   ligne=' A'//va//': '
	   do 70 com=1,comm
	     call int2char(com*10**arcfield+arc,va,nosign)
	     call append(ligne,' +x'//va)
	     if(mod(com,6) .eq. 0) then
		write(10,*) ligne
		ligne=lignevide
	     endif
 70	   continue
	   if((capbin .or. ((comm .eq. 1) .and. ccabin))
     &         .and. (u(arc) .ne. 0)) then
		call int2char(arc,va,nosign)
	        call int2char(-u(arc),chf,sign)
		block=chf//' y'//va
		call append(ligne,block)
		call append(ligne,' <= 0')
	   else
	   if((.not. capbin) .and. (comm .gt. 1)) then
		call int2char(u(arc),chf,sign)
		call append(ligne,' <= '//chf)
	   endif
	   endif
	   write(10,*) ligne
	   ligne=lignevide
 60	continue
	endif

c	--- capacites par produit ---

	if((comm .gt. 1) .and. ccabin) then
	do 81 arc=1,na
	   do 91 com=1,comm
	     call int2char(com*10**arcfield+arc,va,nosign)
	     ligne=' P'//va//': x'//va
	     call int2char(-cap(com,arc),chf,sign)
	     call append(ligne,chf//' y'//va//' <= 0')
	     write(10,*) ligne
 91	   continue
 81	continue
	endif
	ligne=lignevide

c  ************** BOUNDS ****************************
	write(10,*) 'Bounds'

	if((comm .eq. 1) .and. (.not. (capbin .or. 
     &  ccabin) .or. dow)) then
	  do 99 arc=1,na
	    call int2char(10**arcfield+arc,va,nosign)
	    call int2char(u(arc),chf,sign)
	    write(10,*) ' 0 <= x'//va//' <= '//chf
 99	  continue
	endif

	if((.not. ccabin) .and. (comm .gt. 1) .and.
     &  (.not. dow)) then
	do 100 arc=1,na
	  do 110 com=1,comm
	     call int2char(com*10**arcfield+arc,va,nosign)
	     call int2char(cap(com,arc),chf,sign)
	     write(10,*) ' 0 <= x'//va//' <= '//chf
 110	  continue
 100	continue
	endif

        if((capbin).or.(ccabin))then

	  do 120 arc=1,na
	    call int2char(arc,va,nosign)
	    write(10,*) ' 0 <= y'//va//' <= 1'
 120	  continue

	  write(10,*) 'Integers'
	  do 130 arc=1,na
	    call int2char(arc,va,nosign)
	    call append(ligne,' y'//va)
	    if(mod(compt,10) .eq. 0) then
		write(10,*) ligne
		ligne=lignevide
	    endif
	    compt=compt+1
 130	  continue
	write(10,*) ligne

	endif

	write(10,*) 'End'

	end


