      subroutine outdat

      include 'dim.par'
      include 'bipart.cmn'
      include 'ran2.cmn'

      if(caploc)then
      iarcs=origin*(odarcs+1)
      nodes=origin+destin+1
      icoms=1
      else
      iarcs=origin*odarcs 
      nodes=origin+destin 
      icoms=commod
      endif

      write(udataf,1000) nodes,iarcs,icoms

      do 1 ia=1,iarcs
         write(udataf,1001) from(ia),to(ia),fc(ia),tcap(ia),icoms
         do 2 k=1,icoms
            write(udataf,1000) k,tc(ia,k),pcap(ia,k)
   2     continue
   1  continue

      do 3 k=1,icoms
      do 3 i=1,nodes
         write(udataf,1000) k,i,sup(i,k)
   3  continue

1000  format(3i8)
1001  format(5i8)

      return
      end
