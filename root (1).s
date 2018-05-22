section .data  

        result: db "root = %.17g %.17g" , 10, 0
        epsilonInput: db "epsilon = %lf ", 0
        orderInput: db "order = %d ", 0
        coeffInput: db "coeff %d ", 0
        coeffValueInput: db "= %lf %lf ", 0
        firstGuess: db "initial = %lf %lf", 0
        guessReal: DQ 0.0
        guessImg: DQ 0.0
        nirmul: DQ 0.0
        nirmulTemp: DQ 0.0
        DivTemp: DQ 0.0
        DerivTemp: DQ 0.0
        mallocation: DQ 0.0
        originalDegree: DQ 0.0
        mallocation2: DQ 0.0
        epsilonTwice: DQ 0.0
        
section .bss

        epsilon: resq 1  
        length: resq 1
        index: resq 1
        myPoly: resq 1
        
        
        
section .text 

        global main
        extern malloc
        extern printf
        extern __isoc99_scanf
        
        
%macro initialize 0
	push rbp
	mov rbp, rsp
	finit
%endmacro

%macro terminate 0
	mov rsp, rbp
	pop rbp
	ret
%endmacro

%macro myMalloc 1
        mov rdi, %1
        call malloc
%endmacro

%macro twoPop 2
        pop %1
        pop %2
%endmacro

%macro threePop 3
        pop %1
        pop %2
        pop %3
%endmacro

%macro fourPop 4
        pop %1
        pop %2
        pop %3
        pop %4
%endmacro

%macro popAll 7
    pop %1
    pop %2
    pop %3
    pop %4
    pop %5
    pop %6
    pop %7
%endmacro

%macro twoPush 2
    push %1
    push %2
%endmacro

%macro pushAll 7
    push %1
    push %2
    push %3
    push %4
    push %5
    push %6
    push %7
%endmacro

%macro checkOrderZero 0
    mov r10, qword [length]
    cmp r10, 0
    je orderZero
%endmacro           
        
    main:
        initialize
        ; receiving the epsilon from the input
        mov rdi, epsilonInput  
        mov rsi, epsilon    
        mov rax, 0
        call __isoc99_scanf
        ; receiving the epsilon from the input
        mov rdi, orderInput 
        mov rsi, length   
        mov rax, 0
        call __isoc99_scanf
        ; receiving the polynom from the input
        mov rax, qword [length]
        mov qword[originalDegree], rax
        mov rax, 16
        mul qword[originalDegree]        
        add rax, 16   
        mov rdi, rax
        call malloc
        mov qword[myPoly], rax
        mov r11, qword[originalDegree] 
        inc r11     
    
    coeffLoop: ;receiving the coefficient numbers in loop until we reaches the last degree
        mov rdi, coeffInput    
        mov rsi, index     
        mov rax, 0
        call __isoc99_scanf
        
        mov rsi, qword[myPoly]  
        mov rax, 16  
        mov r12, qword [length]
        sub r12, qword [index]
        mul r12
        add rsi, rax         
        mov rdi, coeffValueInput
        lea rsi, [rsi] 
        lea rdx, [rsi + 8]   
        mov rax, 0
        call __isoc99_scanf
        
        dec r11
        cmp r11, 0
        jg coeffLoop
        
        ; receiving the first guess from the input
        mov rdi, firstGuess
        lea rsi, [guessReal]     
        lea rdx, [guessImg]     
        mov rax, 0
        call __isoc99_scanf
        
        checkOrderZero ;checks if a given polynom is with order zero 
        ; putting input arguments in the order that we need 
        mov rdi, qword [myPoly]
        mov rsi, guessReal
        mov rdx, guessImg
        mov rcx, qword [originalDegree]
        mov r8, epsilon
        
        
        call rootEval ; the call for the core of our program
        terminate
        
                        
        rootEval: ; searches for the roots of the polynom
                initialize
                ; multiplication of epsilon by itself
                fld qword [r8]
                fst st1
                fmul st1
                fst qword [epsilonTwice] 
                ; evaluating the darivative of the polynom for later usage
                mov qword [originalDegree], rcx
                mov r14,rsi 
                mov r15,rdx 
                twoPush r14, r15
                mov rdx, rcx
                call calcDeriv
                ; evaluating the guess in the polynom
                twoPush rdi, rsi
                myMalloc 64
                mov qword [mallocation2], rax
                add rax, 48
                mov qword [mallocation], rax
                sub rax, 48
                fourPop rsi, rdi, r15, r14
                mov r13, rsi 
                push r13
                mov r9, rdi
                mov rdi, r14
                mov rsi, r15
                mov r8, qword [originalDegree]
                mov rdx, qword [mallocation2]
                mov rcx, qword [mallocation2]
                add rcx, 8
                push r15
                call calcEval
                ; evaluating the guess in the darivative of the polynom
                pop r15
                add rcx,16
                add rdx,16
                pop r13
                push r9
                mov r9, r13
                mov r8, qword [originalDegree]
                dec r8
                mov rdi, r14
                mov rsi, r15
                twoPush r13, r15
                call calcEval
                ; dividing the evaluation of the the guess in the the polynom by the guess in the darivative of the polynom
                threePop r15, r13, r9
                sub rdx, 16
                pushAll rdi, rsi, rdx, rcx, r8, r9, r13         
                mov rdi,rdx
                add rdx, 8
                mov rsi, rdx
                add rdx, 8
                add rcx, 8
                mov r8, rcx
                add rcx, 8
                mov r9, rcx
                sub rcx, 16
                call division  
                
                popAll r13, r9, r8, rcx, rdx, rsi, rdi
                sub rcx, 16
                
                finit
                fld qword [rdx]
                fst st1
                fmul
                fst qword [DivTemp]
                fld qword [rcx]
                fst st1
                fmul
                fld qword [DivTemp]
                fadd
                fst qword [nirmul]
                fld qword [nirmul]
                fcomp qword [epsilonTwice]
                fstsw ax
                sahf 
                jbe printResult
                
            rootRangeLoop: ;checks if the root is close enough to the epsilon ranges
                ; putting the next in r14 and r15 by the given formula 
                ; keep on itertating until you get closer enough to the root
                finit
                jmp subtraction
            afterSubtraction:
                mov rdi, r14
                mov rsi, r15
                mov r8, qword [originalDegree]
                push r13
                call calcEval
                pop r13
                push r15
                pop r15
                add rcx,16
                add rdx,16
                push r9
                mov r9, r13
                mov r8, qword [originalDegree]
                dec r8
                mov rdi, r14
                mov rsi, r15
                twoPush r13, r15
                call calcEval
                threePop r15, r13, r9
                sub rdx, 16
                sub rcx, 16
                pushAll rdi, rsi, rdx, rcx, r8, r9, r13
                mov rdi,rdx
                add rdx, 8
                mov rsi, rdx
                add rdx, 8
                add rcx, 16
                mov r8, rcx
                add r8, 8
                mov r9, rcx
                add r9,16
                call division
                popAll r13, r9, r8, rcx, rdx, rsi, rdi
                finit
                fld qword [rdx]
                fst st1
                fmul
                fst qword [DivTemp]
                fld qword [rcx]
                fst st1
                fmul
                fld qword [DivTemp]
                fadd
                fst qword [nirmul]
                fld qword [nirmul]
                fcomp qword [epsilonTwice]
                fstsw ax
                sahf 
                jbe printResult
                jmp rootRangeLoop
                
    multiplication:
        initialize
        ;number1_real * number2_real
        fld qword [rdi]	
        fst st1			
	fld qword [rdx]		
	fmul
	;number1_img * number2_img
	fst qword [r8]
	fld qword [rsi]	
        fst st1			
	fld qword [rcx]		
	fmul
	;subtract the numbers
	fst st1
	fld qword [r8]
	fsub st1
	;number1_img * number2_real
	fst qword [r8]
	fld qword [rsi]	
        fst st1			
	fld qword [rdx]		
	fmul
	;number1_real * number2_img
	fst qword [r9]
	fld qword [rdi]	
        fst st1			
	fld qword [rcx]		
	fmul
	; add between the results
	fst st1
	fld qword [r9]
	fadd
	fst qword [r9]
	terminate
	
    division:
        initialize
        mov qword [r8], 0
        mov qword [r9], 0
        ;number1_real * number2_real + number1_img * number2_img 
        fld qword [rdi]	
        fst st1			
	fld qword [rdx]		
	fmul
	fst qword [r8]
	fld qword [rsi]	
        fst st1			
	fld qword [rcx]		
	fmul
	fst st1
	fld qword [r8]
	fadd
	;number2_real^2 + number2_img^2
	fst qword [r8]
	fld qword [rdx]
	fst st1
	fmul
	fst qword [DivTemp]
	fld qword [rcx]
	fst st1
	fmul
	fst st1
	fld qword [DivTemp]
	fadd
	;division between the results
	fst st1
	fld qword [r8]
	fdiv st1
	;number1_img*number2_real - number1_real*number2_img
	fst qword [r8]	
	fld qword [rsi]
	fst st1
	fld qword [rdx]
	fmul
	fst qword [DivTemp]
	fld qword [rdi]
	fst st1
	fld qword [rcx]
	fmul
	fst st1
	fld qword [DivTemp]
	fsub st1
	fst qword [r9]
	fld qword [rdx]
	fst st1
	fmul
	fst qword [DivTemp]
	fld qword [rcx]
	fst st1
	fmul
	fst st1
	fld qword [DivTemp]
	fadd
	fst st1
	fld qword [r9]
	;division between the results
	fdiv st1
	fst qword [r9]
	terminate
	
    addition:
        initialize
        mov qword [r8], 0
        mov qword [r9], 0
        ;number1_real + number2_real
	fld qword [rdi]		
	fst st1			
	fld qword [rdx]		
	fadd			
	;number1_img + number2_img
	fst qword [r8]		
	fld qword [rsi]		
	fst st1			
	fld qword [rcx]		
	fadd			
	fst qword [r9]	
	terminate
	
    subtraction:
        fld qword [rdx+32]
        fst st1
        fld qword [r14]
        fsub st0, st1
        fst qword [r14]
        fld qword [rdx+40]
        fst st1
        fld qword [r15]
        fsub st0, st1
        fst qword [r15]
	jmp afterSubtraction

	
    calcDeriv: ;calculate the derivative of the polynom
        initialize
        mov r14,0 
        mov r15, rdx ; Index of the for loop
        twoPush rdx, rdi
        mov rax, 2 
        mul rdx
        mov rdi, 8
        mul rdi
        mov rdi, rax 
        mov r12, rdi
        call malloc
        mov qword [DerivTemp], rax
        mov rsi, qword [DerivTemp]
        twoPop rdi, rdx
        firstLoop: ; 
            cmp r15, 0
            je endOfLoop
            mov r13, rdx ;length
            sub r13, r14
            mov [nirmulTemp], r13
            fild qword [nirmulTemp]
            fst st1
            fld qword [rdi]
            fmul
            fst qword [rsi]
            add rsi,8
            add rdi,8
            fild qword [nirmulTemp]
            fst st1
            fld qword [rdi]
            fmul
            fst qword [rsi]
            add rdi, 8
            add rsi, 8
            inc r14
            dec r15
            jmp firstLoop
            
    endOfLoop:
            sub rsi, r12
            sub rdi, r12
            terminate
            
    calcEval:    ; evaluating the given numbers in the polynom
            initialize
            mov r14, rdi
            mov r15, rsi 
            mov r12, r8 
            mov r11, r8 
            mov r13, r9 
            mov r10, r9 
            twoPush rdi, rsi
            mov rdi, rdx 
            mov rsi, rcx 
            twoPop rcx, rdx
            add r13,8
            fld qword [r9] 
            fst qword [rdi]
            fld qword [r9+8] 
            fst qword [rsi]
            mov r8, qword [mallocation]
            add r8, 8
            mov r9, r8
            sub r8, 8
            
                secondLoop: ; looping on the each order of the polynom
                twoPush rdi, rsi
                cmp r12,0
                je endOfEval
                call multiplication 
                mov rdi, r8 
                mov rsi, r9
                add r13, 8
                mov rdx, r13 
                add r13, 8
                mov rcx, r13
                twoPop r9, r8
                call addition 
                mov rdx, r14 
                mov rcx, r15
                twoPush r8, r9
                mov r8, rdi 
                mov r9, rsi
                twoPop rsi, rdi
                dec r12
                jmp secondLoop
            
            endOfEval: 
                twoPush rdi, rsi
                mov rdi, r14
                mov rsi, r15
                twoPop rcx, rdx
                mov r8, r11
                mov r9, r10
                terminate

                
            printResult: ; printing the real and imaginary numbers of the result
                movsd xmm0, qword[r14]
                movsd xmm1, qword[r15]
                mov rdi, result
                mov rax, 2
                call printf 
                terminate
                
            orderZero: ; if the order is zero, prepare for printing the initial
                mov r14, guessReal
                mov r15, guessImg
                jmp printResult
	
	
	
	
