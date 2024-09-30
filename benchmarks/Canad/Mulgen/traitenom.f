c **************************************************
c Sous-routine pour trouver le nom du fichier de sortie
c ****************************************************

	subroutine traitenom(no)

	include 'dim.par'
	include "mulgen.cmn"

	integer no
	character chf
	
	mpsfile=outfile
	rnetfile=outfile
	dowfile=outfile
	lpfile=outfile

	if(no .eq. 0) then

	do 100 i=1,len(outfile)
	   if(outfile(i:i) .eq. char(32)) then
		mpsfile(i:i+4)='.mps'
		rnetfile(i:i+5)='.rnet'
		dowfile(i:i+4)='.dow'
		lpfile(i:i+3)='.lp'
		return
	   endif
 100	continue

	else

	call int2char(no,chf,.false.)
	do 200 i=1,len(outfile)
	   if(outfile(i:i) .eq. char(32)) then
		lpfile(i:i+5)='.'//chf//'.lp'
	        mpsfile(i:i+6)='.'//chf//'.mps'
		rnetfile(i:i+7)='.'//chf//'.rnet'
		dowfile(i:i+6)='.'//chf//'.dow'
		return
	   endif
 200	continue

	endif
	end

