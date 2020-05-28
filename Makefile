#
# Les variables d'environnement libG2X, incG2X
# sont definies dans le fichier ~/.bashrc par le script ../install.sh
#
#compilateur
CC = gcc
#compil en mode 'debug' ou optmis�e (-O2)
DBG = yes

ifeq ($(DBG),yes) #en mode debug
  CFLAGS = -g -Wpointer-arith -Wall -ansi
else              #en mode normal
  CFLAGS = -O2 -ansi
endif

# assemblage des infos de lib. et inc.
lib =   $(libG3X) $(libG2X)
# fichiers *.c locaux
src = src/
# fichiers *.h locaux et lib.
inc = -I./include $(incG3X) $(incG2X)

# r�gle de compilation g�n�rique des objets
%.o : $(src)%.c
	@echo "module $@"
	@$(CC) $(CFLAGS) $(inc) -c $< -o $@
	
# r�gle de compilation g�n�rique
% : %.o 
	@echo "edition de lien $^ -> $@"
	@$(CC) $^ $(lib) -o $@
	
SimpleOsc : PMat.o Link.o SimpleOsc.o	
	@echo "edition de lien $^ -> $@"
	@$(CC) $^ $(lib) -o $@

Corde1 : PMat.o Link.o Corde1.o	
	@echo "edition de lien $^ -> $@"
	@$(CC) $^ $(lib) -o $@

Corde2 : PMat.o Link.o Corde2.o	
	@echo "edition de lien $^ -> $@"
	@$(CC) $^ $(lib) -o $@

.PHONY : clean cleanall

clean : 
	rm -f *.o .tmp*
cleanall : 
	rm -f *.o .tmp* $(EXEC)
