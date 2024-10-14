#
# declaration des repertoires
#
SOURCE = .
EXEC = .

#
# declaration des fichiers-source
#
SRC =	$(SOURCE)/main.f\
	$(SOURCE)/bipart.f\
	$(SOURCE)/outdat.f\
	$(SOURCE)/ran2.f

#
# declaration des fichiers-objet
#
OBJ =	$(EXEC)/main.o\
	$(EXEC)/bipart.o\
	$(EXEC)/outdat.o\
	$(EXEC)/ran2.o

#
# declaration des fichiers de dimension
#
DIM = $(SOURCE)/dim.par
 
#
# declaration des fichiers de common
#
CMN = $(SOURCE)/bipart.cmn $(SOURCE)/ran2.cmn

F77= f77

LIB= +U77

default : $(EXEC)/bipart

clean : chdate default

chdate :
	touch -m $(SRC)

$(SRC) : $(DIM) $(CMN)
	touch -m $(SRC)

$(EXEC)/bipart : $(OBJ)
	$(F77) -o $@ -C -O $(OBJ) $(LIB)

$(EXEC)/main.o : $(SOURCE)/main.f
	$(F77) -c -O $? $(LIB)
$(EXEC)/bipart.o : $(SOURCE)/bipart.f
	$(F77) -c -O $? 
$(EXEC)/outdat.o : $(SOURCE)/outdat.f
	$(F77) -c -O $? 
$(EXEC)/ran2.o : $(SOURCE)/ran2.f
	$(F77) -c -O $? 
