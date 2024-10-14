#
# declaration des repertoires
#
SOURCE = .
EXEC = .

#
# declaration des fichiers-source
#
SRC =	$(SOURCE)/main.f\
	$(SOURCE)/mulgen.f\
	$(SOURCE)/ran2.f\
	$(SOURCE)/outdat.f\
	$(SOURCE)/mps.f\
	$(SOURCE)/lp.f\
	$(SOURCE)/append.f\
	$(SOURCE)/traitenom.f\
	$(SOURCE)/int2char.f 

#
# declaration des fichiers-objet
#
OBJ =	$(EXEC)/main.o\
	$(EXEC)/mulgen.o\
	$(EXEC)/ran2.o\
	$(EXEC)/outdat.o\
	$(EXEC)/mps.o\
	$(EXEC)/lp.o\
	$(EXEC)/append.o\
	$(EXEC)/traitenom.o\
	$(EXEC)/int2char.o 

#
# declaration des fichiers de dimension
#
DIM = $(SOURCE)/dim.par

#
# declaration des fichiers de common
#
CMN = $(SOURCE)/mulgen.cmn

F77= f77

LIB= +U77

#
# declaration des cibles principales
#
default : $(EXEC)/mulgen

clean : chdate default

chdate : 
	touch -m $(SRC)

$(SRC) : $(DIM) $(CMN)
	touch -m $(SRC)

$(EXEC)/mulgen : $(OBJ)
	$(F77) -o $@ -C -O $(OBJ) $(LIB)

$(EXEC)/main.o : $(SOURCE)/main.f
	$(F77) -c -O $? $(LIB)
$(EXEC)/mulgen.o : $(SOURCE)/mulgen.f
	$(F77) -c -O $?
$(EXEC)/ran2.o : $(SOURCE)/ran2.f
	$(F77) -c -O $?
$(EXEC)/outdat.o : $(SOURCE)/outdat.f
	$(F77) -c -O $?
$(EXEC)/mps.o : $(SOURCE)/mps.f
	$(F77) -c -O $?
$(EXEC)/lp.o : $(SOURCE)/lp.f
	$(F77) -c -O $?
$(EXEC)/append.o : $(SOURCE)/append.f
	$(F77) -c -O $?
$(EXEC)/traitenom.o : $(SOURCE)/traitenom.f
	$(F77) -c -O $?
$(EXEC)/int2char.o : $(SOURCE)/int2char.f
	$(F77) -c -O $?
