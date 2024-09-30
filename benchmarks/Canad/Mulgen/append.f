	subroutine append(un,deux)

	character*(*) un,deux
	integer index

	do 100 i=len(un),1,-1
	   if((un(i:i) .gt. char(33))
     &  .and. (un(i:i) .lt. char(126))) then
		index=i+1
		un(index:index+len(deux))=deux
	   	return
	   endif
 100	continue

	un(2:len(deux)+2)=deux
	return
	end

