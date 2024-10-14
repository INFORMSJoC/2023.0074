	subroutine int2char(nb,chr,flg)

	integer nb, temp, rest, nbr
	character*(*) chr
	character a
	logical flg

	if((nb .eq. 0) .and. .not. flg) then
	   chr='0'
	   return
	endif

	if((nb .eq. 0) .and. flg) then
	   chr='+0'
	   return
        endif

	if(flg) then
	if(nb .lt. 0) then
	   chr(1:1)='-'
	else
	   chr(1:1)='+'
	endif
	endif

	nbr=iabs(nb)
	long=log10(float(nbr))+1
	rest=nbr
	do 100 i=1,long
	   temp=rest/10**(long-i)
	   a=char(temp+48)
	   if(flg) then
		chr(i+1:i+1)=a
	   else
		chr(i:i)=a
	   endif  
	   rest=rest-temp*10**(long-i)
 100	continue

	if(flg) then
	do 150 i=long+2,len(chr)
	  chr(i:i)=' '
 150	continue

	else
	do 200 i=long+1,len(chr)
	   chr(i:i)=' '
 200	continue
	endif

	return

	end

