      subroutine bipart

      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'

      real tmpran
      real ran2

      real q
      integer fsum,csum
      integer msum,osum,dsum
      integer totsupply

      integer isupl, irest, itsup, itemp, idemd, imax
      integer jdest
      integer iarc, isply, idmnd

      if(caploc)commod=1

c a premiliminary call to the random generator
      k= int(ran2(iseed))

c generate supplies and demands

      totsupply=0
      do 2 k=1,commod
      tmpran= ran2(iseed)
      itsup=minsup+int(tmpran*(maxsup-minsup))
      totsupply=totsupply+itsup

      isupl=0 
      irest=0
      do 3 i=1,origin
         tmpran= ran2(iseed)
         itemp=1+int(tmpran*(itsup/origin))
         sup(i,k)=itemp+irest
         isupl=isupl+itemp+irest
         irest=int(itsup/origin)-itemp
   3  continue

      sup(1,k)=sup(1,k)+irest 
      isupl=isupl+irest

      idemd=0
      irest=0
      do 4 i=origin+1,origin+destin
         tmpran= ran2(iseed)
         itemp=1+int(tmpran*(itsup/destin)) 
         sup(i,k)=-itemp-irest 
         idemd=idemd+itemp+irest 
         irest=int(itsup/destin)-itemp 
   4  continue
      sup(origin+1,k)=sup(origin+1,k)-irest
      idemd=idemd+irest 

      if(isupl.gt.idemd)then
        sup(origin+1,k)=sup(origin+1,k)-(isupl-idemd)
      else
      if(isupl.lt.idemd)then 
        sup(1,k)=sup(1,k)+(idemd-isupl)
      endif
      endif

   2  continue


c generate arcs between origins and destinations

      iarc=0
      do 5 i=1,origin
         jdest=origin+1+int(ran2(iseed)*(destin-1))
         do 6 j=1,odarcs
            iarc=iarc+1
            from(iarc)=i
            to(iarc)=jdest
            fc(iarc)=minfct+int(ran2(iseed)*(maxfct-minfct))
            tcap(iarc)=mincap+int(ran2(iseed)*(maxcap-mincap))
            isply=0
            idmnd=0
            do 7 k=1,commod
               isply=isply+sup(i,k)
               idmnd=idmnd-sup(jdest,k)
               tc(iarc,k)=mintct+int(ran2(iseed)*(maxtct-mintct))
               pcap(iarc,k)=min0(tcap(iarc),sup(i,k),-sup(jdest,k))
   7        continue             
            tcap(iarc)=min0(tcap(iarc),isply,idmnd)
            jdest=mod(jdest,destin)+1+origin
   6     continue
   5  continue

       if(ctight)then
         csum=0
         msum=0
         do 111 i=1,iarc
           csum=csum+tcap(i)
           ifrom=from(i)
           jdest=to(i)
           osum=0
           dsum=0
           do 115 k=1,commod
             osum=osum+sup(ifrom,k)
             dsum=dsum-sup(jdest,k)
115        continue
           imax=min0(osum,dsum)
           msum=msum+imax           
111      continue
         q=dfloat(msum)/dfloat(csum)

         do 112 i=1,iarc
           tcap(i)=int((q/xtight)*tcap(i))+1
112      continue

         csum=0
         do 113 i=1,iarc
           csum=csum+tcap(i)
113      continue
       endif

       if(fixvar)then
         fsum=0
         do 222 i=1,iarc
           fsum=fsum+fc(i)
222      continue
         csum=0
         do 223 k=1,commod
         do 223 i=1,iarc
            csum=csum+tcap(i)*tc(i,k)/commod
223      continue
         q=dfloat(fsum)/dfloat(csum)
 
         do 225 i=1,iarc
           fc(i)=int((xfixvar/q)*fc(i))+1
225      continue
 
         fsum=0
         do 226 i=1,iarc
           fsum=fsum+fc(i)
226      continue
       endif
         
      if(caploc)then
        sup(origin+destin+1,1)=itsup
        do 10 i=1,origin
           sup(i,1)=0
           iarc=iarc+1
           from(iarc)=origin+destin+1
           to(iarc)=i
           fc(iarc)=minfct+int(ran2(iseed)*(maxfct-minfct))
           tcap(iarc)=mincap+int(ran2(iseed)*(maxcap-mincap))
           tcap(iarc)=sup(i,1)+tcap(iarc)
           tc(iarc,1)=0
           pcap(iarc,1)=tcap(iarc)
  10    continue
        do 20 ia=1,origin*odarcs
           fc(ia)=0
           tcap(iarc)=maxsup+maxcap
           pcap(iarc,1)=maxsup+maxcap
  20    continue
      endif     

      return
      end
