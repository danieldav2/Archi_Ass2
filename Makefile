all: root

root: root.o
	gcc -g -Wall -o root root.o
	
root.o: root.s
	nasm -g -f elf64 -w+all -o root.o root.s

.PHONY: clean


clean: 
	rm -f *.o root

        