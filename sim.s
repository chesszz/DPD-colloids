	.file	"sim.c"
	.section	.text.unlikely,"x"
LCOLDB2:
	.text
LHOTB2:
	.p2align 4,,15
	.globl	_refold_positions
	.def	_refold_positions;	.scl	2;	.type	32;	.endef
_refold_positions:
LFB11:
	.cfi_startproc
	pushl	%ebx
	.cfi_def_cfa_offset 8
	.cfi_offset 3, -8
	subl	$64, %esp
	.cfi_def_cfa_offset 72
	movl	84(%esp), %edx
	movl	76(%esp), %eax
	movl	80(%esp), %ecx
	movl	72(%esp), %ebx
	movl	%edx, 8(%esp)
	movl	88(%esp), %edx
	movl	%eax, (%esp)
	movl	%ecx, 4(%esp)
	movl	%edx, 12(%esp)
	movl	92(%esp), %edx
	movl	%edx, 16(%esp)
	movl	96(%esp), %edx
	movl	%edx, 20(%esp)
	movl	100(%esp), %edx
	movl	%edx, 24(%esp)
	movl	104(%esp), %edx
	movl	%edx, 28(%esp)
	movl	108(%esp), %edx
	movl	%edx, 32(%esp)
	movl	112(%esp), %edx
	movl	%edx, 36(%esp)
	movl	116(%esp), %edx
	movl	%edx, 40(%esp)
	movl	120(%esp), %edx
	movl	%edx, 44(%esp)
	movl	124(%esp), %edx
	movl	%edx, 48(%esp)
	movl	128(%esp), %edx
	movl	%edx, 52(%esp)
	movl	132(%esp), %edx
	movl	%edx, 56(%esp)
	movl	136(%esp), %edx
	fldl	24(%esp)
	movl	%edx, 60(%esp)
	leal	(%eax,%eax,2), %edx
	testl	%edx, %edx
	jle	L10
	movl	(%ebx), %eax
	fldz
	leal	(%eax,%edx,8), %edx
	jmp	L9
	.p2align 4,,10
L26:
	fsub	%st(2), %st
	fstpl	(%eax)
	jmp	L7
	.p2align 4,,10
L29:
	fstp	%st(0)
L7:
	addl	$8, %eax
	cmpl	%edx, %eax
	je	L28
L9:
	fldl	(%eax)
	fucomi	%st(2), %st
	ja	L26
	fld	%st(1)
	fucomip	%st(1), %st
	jbe	L29
	fadd	%st(2), %st
	addl	$8, %eax
	fstpl	-8(%eax)
	cmpl	%edx, %eax
	jne	L9
	fstp	%st(0)
	jmp	L10
	.p2align 4,,10
L28:
	fstp	%st(0)
L10:
	leal	(%ecx,%ecx,2), %edx
	testl	%edx, %edx
	jle	L30
	movl	12(%ebx), %eax
	fldz
	leal	(%eax,%edx,8), %edx
	jmp	L15
	.p2align 4,,10
L27:
	fsub	%st(2), %st
	fstpl	(%eax)
	jmp	L13
	.p2align 4,,10
L32:
	fstp	%st(0)
L13:
	addl	$8, %eax
	cmpl	%eax, %edx
	je	L31
L15:
	fldl	(%eax)
	fucomi	%st(2), %st
	ja	L27
	fld	%st(1)
	fucomip	%st(1), %st
	jbe	L32
	fadd	%st(2), %st
	addl	$8, %eax
	fstpl	-8(%eax)
	cmpl	%eax, %edx
	jne	L15
	fstp	%st(0)
	fstp	%st(0)
	jmp	L1
L30:
	fstp	%st(0)
	jmp	L1
	.p2align 4,,10
L31:
	fstp	%st(0)
	fstp	%st(0)
L1:
	addl	$64, %esp
	.cfi_def_cfa_offset 8
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 4
	ret
	.cfi_endproc
LFE11:
	.section	.text.unlikely,"x"
LCOLDE2:
	.text
LHOTE2:
	.section .rdata,"dr"
	.align 4
LC3:
	.ascii "neigh_counter < in.N_WATER * (in.N_WATER) / 2\0"
LC4:
	.ascii "sim.c\0"
LC5:
	.ascii "j != -255\0"
	.align 4
LC6:
	.ascii "j == -1 || (j >= 0 && j < in.N_WATER)\0"
LC9:
	.ascii "dist_weight > 0\0"
LC10:
	.ascii "dist_weight <= 1\0"
	.section	.text.unlikely,"x"
LCOLDB13:
	.text
LHOTB13:
	.p2align 4,,15
	.globl	_calculate_acc
	.def	_calculate_acc;	.scl	2;	.type	32;	.endef
_calculate_acc:
LFB12:
	.cfi_startproc
	pushl	%ebp
	.cfi_def_cfa_offset 8
	.cfi_offset 5, -8
	pushl	%edi
	.cfi_def_cfa_offset 12
	.cfi_offset 7, -12
	pushl	%esi
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushl	%ebx
	.cfi_def_cfa_offset 20
	.cfi_offset 3, -20
	subl	$220, %esp
	.cfi_def_cfa_offset 240
	movl	244(%esp), %eax
	movl	248(%esp), %ebx
	movl	%eax, %edi
	movl	%eax, 28(%esp)
	movl	%eax, 144(%esp)
	movl	252(%esp), %eax
	movl	%ebx, 148(%esp)
	movl	%eax, 152(%esp)
	movl	256(%esp), %eax
	movl	%eax, 156(%esp)
	movl	260(%esp), %eax
	movl	%eax, 160(%esp)
	movl	264(%esp), %eax
	movl	%eax, 164(%esp)
	movl	268(%esp), %eax
	movl	%eax, 168(%esp)
	movl	272(%esp), %eax
	movl	%eax, 172(%esp)
	movl	276(%esp), %eax
	movl	%eax, 176(%esp)
	movl	280(%esp), %eax
	movl	%eax, 180(%esp)
	movl	284(%esp), %eax
	movl	%eax, 184(%esp)
	movl	288(%esp), %eax
	movl	%eax, 188(%esp)
	movl	292(%esp), %eax
	movl	%eax, 192(%esp)
	movl	296(%esp), %eax
	fldl	160(%esp)
	movl	%eax, 196(%esp)
	movl	300(%esp), %eax
	fstpl	32(%esp)
	fldl	184(%esp)
	movl	%eax, 200(%esp)
	movl	304(%esp), %eax
	fldl	192(%esp)
	movl	%eax, 204(%esp)
	leal	(%edi,%edi,2), %eax
	fstpl	120(%esp)
	testl	%eax, %eax
	jle	L37
	movl	240(%esp), %edi
	fstpl	40(%esp)
	sall	$3, %eax
	movl	8(%edi), %edx
	movl	%eax, 8(%esp)
	movl	$0, 4(%esp)
	movl	%edx, (%esp)
	call	_memset
	fldl	40(%esp)
L37:
	leal	(%ebx,%ebx,2), %eax
	testl	%eax, %eax
	jle	L36
	movl	240(%esp), %edi
	fstpl	40(%esp)
	sall	$3, %eax
	movl	20(%edi), %edx
	movl	%eax, 8(%esp)
	movl	$0, 4(%esp)
	movl	%edx, (%esp)
	call	_memset
	fldl	40(%esp)
L36:
	fldl	32(%esp)
	fsqrt
	fstl	112(%esp)
	fucomip	%st(0), %st
	jp	L78
L38:
	fld	%st(0)
	fadd	%st(1), %st
	fld	%st(0)
	fsqrt
	fucomi	%st(0), %st
	jp	L79
	fstp	%st(1)
L40:
	movl	28(%esp), %eax
	testl	%eax, %eax
	jle	L85
	movl	308(%esp), %eax
	movl	28(%esp), %edi
	movl	(%eax), %edx
	movl	%edi, %eax
	imull	%edi, %eax
	sarl	%eax
	testl	%eax, %eax
	movl	%eax, 56(%esp)
	jle	L86
	cmpl	$-255, %edx
	je	L87
	fxch	%st(1)
	fchs
	xorl	%ebp, %ebp
	xorl	%edi, %edi
	fstpl	128(%esp)
	fstpl	136(%esp)
	jmp	L48
	.p2align 4,,10
L81:
	cmpl	28(%esp), %edx
	jge	L69
	movl	%edx, %eax
	shrl	$31, %eax
	testb	%al, %al
	jne	L69
	movl	240(%esp), %eax
	fldl	168(%esp)
	leal	(%edi,%edi,2), %esi
	leal	(%edx,%edx,2), %ebx
	movl	(%eax), %eax
	leal	0(,%ebx,8), %ecx
	movl	%ecx, 32(%esp)
	fldl	(%eax,%esi,8)
	fsubl	(%eax,%ebx,8)
	flds	LC7
	fld	%st(2)
	fmul	%st(1), %st
	fxch	%st(2)
	fucomi	%st(2), %st
	ja	L53
	fld	%st(3)
	fchs
	fmulp	%st, %st(2)
	fxch	%st(1)
	fucomip	%st(1), %st
	jbe	L54
	fadd	%st(2), %st
L54:
	leal	1(%esi), %ecx
	leal	0(,%ecx,8), %edx
	movl	%ecx, 48(%esp)
	movl	%edx, 32(%esp)
	leal	1(%ebx), %edx
	leal	0(,%edx,8), %ecx
	movl	%ecx, 40(%esp)
	movl	48(%esp), %ecx
	fldl	(%eax,%ecx,8)
	fsubl	(%eax,%edx,8)
	fucomi	%st(2), %st
	ja	L56
	fld	%st(3)
	fchs
	fmuls	LC7
	fucomip	%st(1), %st
	jbe	L57
	fadd	%st(3), %st
L57:
	movl	32(%esp), %ecx
	fldl	8(%eax,%ecx)
	movl	40(%esp), %ecx
	fsubl	8(%eax,%ecx)
	fucomi	%st(3), %st
	fstp	%st(3)
	ja	L59
	fld	%st(3)
	fchs
	fmuls	LC7
	fucomip	%st(3), %st
	jbe	L88
	fxch	%st(2)
	faddp	%st, %st(3)
	jmp	L60
	.p2align 4,,10
L88:
	fstp	%st(3)
	fxch	%st(1)
	fxch	%st(2)
	fxch	%st(1)
L60:
	fld	%st(1)
	fmul	%st(2), %st
	fld	%st(1)
	fmul	%st(2), %st
	faddp	%st, %st(1)
	fld	%st(3)
	fmul	%st(4), %st
	faddp	%st, %st(1)
	fld1
	fucomip	%st(1), %st
	ja	L80
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
L52:
	addl	$1, %ebp
	cmpl	28(%esp), %edi
	jge	L33
L68:
	cmpl	56(%esp), %ebp
	movl	308(%esp), %eax
	movl	(%eax,%ebp,4), %edx
	je	L46
	cmpl	$-255, %edx
	je	L47
L48:
	cmpl	$-1, %edx
	jne	L81
	addl	$1, %edi
	addl	$1, %ebp
	cmpl	28(%esp), %edi
	jl	L68
	jmp	L33
L85:
	fstp	%st(0)
	fstp	%st(0)
L33:
	addl	$220, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 20
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 16
	popl	%esi
	.cfi_restore 6
	.cfi_def_cfa_offset 12
	popl	%edi
	.cfi_restore 7
	.cfi_def_cfa_offset 8
	popl	%ebp
	.cfi_restore 5
	.cfi_def_cfa_offset 4
	ret
	.p2align 4,,10
L59:
	.cfi_restore_state
	fxch	%st(2)
	fsubp	%st, %st(3)
	jmp	L60
	.p2align 4,,10
L56:
	fsub	%st(3), %st
	jmp	L57
	.p2align 4,,10
L53:
	fstp	%st(1)
	fsub	%st(2), %st
	jmp	L54
	.p2align 4,,10
L80:
	fld	%st(0)
	fsqrt
	fucomi	%st(0), %st
	jp	L82
	fstp	%st(1)
L64:
	fld1
	fld	%st(0)
	fsub	%st(2), %st
	fldz
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	L83
	fxch	%st(1)
	fucomip	%st(1), %st
	jb	L84
	fxch	%st(2)
	fdiv	%st(1), %st
	fxch	%st(3)
	leal	0(,%esi,8), %ecx
	movl	240(%esp), %eax
	leal	8(%ecx), %edx
	movl	4(%eax), %eax
	movl	%edx, 40(%esp)
	leal	0(,%ebx,8), %edx
	movl	%eax, 32(%esp)
	leal	8(%edx), %eax
	movl	%eax, 48(%esp)
	leal	16(%ecx), %eax
	movl	%eax, 60(%esp)
	leal	16(%edx), %eax
	fdiv	%st(1), %st
	fxch	%st(4)
	movl	%eax, 72(%esp)
	movl	32(%esp), %eax
	fdivp	%st, %st(1)
	fldl	120(%esp)
	fmul	%st(2), %st
	fstpl	64(%esp)
	fldl	(%eax,%esi,8)
	movl	%ecx, %esi
	fsubl	(%eax,%ebx,8)
	movl	%edx, %ebx
	fmul	%st(3), %st
	fxch	%st(3)
	fstpl	104(%esp)
	fldl	8(%eax,%ecx)
	fsubl	8(%eax,%edx)
	fmul	%st(4), %st
	fxch	%st(4)
	fstpl	96(%esp)
	fxch	%st(3)
	faddp	%st, %st(2)
	fldl	16(%eax,%ecx)
	fsubl	16(%eax,%edx)
	fmul	%st(3), %st
	fxch	%st(3)
	fstpl	88(%esp)
	fxch	%st(1)
	faddp	%st, %st(2)
	fldl	128(%esp)
	fmul	%st(1), %st
	fmul	%st(1), %st
	fxch	%st(1)
	fstpl	80(%esp)
	fmulp	%st, %st(1)
	fstpl	32(%esp)
	call	_rand
	movl	%eax, 76(%esp)
	fildl	76(%esp)
	movl	240(%esp), %ecx
	fdivl	LC11
	movl	8(%ecx), %eax
	movl	60(%esp), %ecx
	leal	(%eax,%esi), %edx
	movl	48(%esp), %esi
	fsubs	LC7
	fldl	80(%esp)
	fmull	136(%esp)
	fmull	LC12
	fmulp	%st, %st(1)
	fdivl	112(%esp)
	fldl	64(%esp)
	faddl	32(%esp)
	faddp	%st, %st(1)
	fldl	104(%esp)
	fmul	%st(1), %st
	fldl	(%edx)
	fadd	%st(1), %st
	fstpl	(%edx)
	movl	40(%esp), %edx
	fldl	96(%esp)
	addl	%eax, %edx
	fmul	%st(2), %st
	fldl	(%edx)
	fadd	%st(1), %st
	fstpl	(%edx)
	leal	(%eax,%ecx), %edx
	fldl	88(%esp)
	fmulp	%st, %st(3)
	fldl	(%edx)
	fadd	%st(3), %st
	fstpl	(%edx)
	fxch	%st(1)
	movl	%ebx, %edx
	addl	%eax, %edx
	fsubrl	(%edx)
	fstpl	(%edx)
	leal	(%eax,%esi), %edx
	addl	72(%esp), %eax
	fsubrl	(%edx)
	fstpl	(%edx)
	fsubrl	(%eax)
	fstpl	(%eax)
	jmp	L52
L69:
	movl	$LC6, 12(%esp)
	movl	$___func__.3692, 8(%esp)
	movl	$238, 4(%esp)
	movl	$LC4, (%esp)
	call	___assert_func
L82:
	fstp	%st(0)
	fxch	%st(1)
	fstpl	48(%esp)
	fxch	%st(1)
	fstpl	40(%esp)
	fxch	%st(1)
	fstpl	32(%esp)
	fstpl	(%esp)
	call	_sqrt
	fldl	48(%esp)
	fldl	40(%esp)
	fldl	32(%esp)
	fxch	%st(3)
	fxch	%st(1)
	fxch	%st(2)
	fxch	%st(1)
	jmp	L64
L86:
	fstp	%st(0)
	fstp	%st(0)
L46:
	movl	$LC3, 12(%esp)
	movl	$___func__.3692, 8(%esp)
	movl	$233, 4(%esp)
	movl	$LC4, (%esp)
	call	___assert_func
L87:
	fstp	%st(0)
	fstp	%st(0)
L47:
	movl	$LC5, 12(%esp)
	movl	$___func__.3692, 8(%esp)
	movl	$237, 4(%esp)
	movl	$LC4, (%esp)
	call	___assert_func
L79:
	fstp	%st(0)
	fxch	%st(1)
	fstpl	32(%esp)
	fstpl	(%esp)
	call	_sqrt
	fldl	32(%esp)
	fxch	%st(1)
	jmp	L40
L78:
	fstpl	40(%esp)
	fldl	32(%esp)
	fstpl	(%esp)
	call	_sqrt
	fstpl	112(%esp)
	fldl	40(%esp)
	jmp	L38
L83:
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	movl	$LC9, 12(%esp)
	movl	$___func__.3692, 8(%esp)
	movl	$271, 4(%esp)
	movl	$LC4, (%esp)
	call	___assert_func
L84:
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	movl	$LC10, 12(%esp)
	movl	$___func__.3692, 8(%esp)
	movl	$272, 4(%esp)
	movl	$LC4, (%esp)
	call	___assert_func
	.cfi_endproc
LFE12:
	.section	.text.unlikely,"x"
LCOLDE13:
	.text
LHOTE13:
	.section	.text.unlikely,"x"
LCOLDB14:
	.text
LHOTB14:
	.p2align 4,,15
	.globl	_get_rel_vector
	.def	_get_rel_vector;	.scl	2;	.type	32;	.endef
_get_rel_vector:
LFB13:
	.cfi_startproc
	pushl	%ebx
	.cfi_def_cfa_offset 8
	.cfi_offset 3, -8
	subl	$64, %esp
	.cfi_def_cfa_offset 72
	movl	76(%esp), %eax
	movl	140(%esp), %ecx
	movl	144(%esp), %edx
	movl	148(%esp), %ebx
	movl	%eax, (%esp)
	movl	80(%esp), %eax
	leal	(%ecx,%ecx,2), %ecx
	leal	(%edx,%edx,2), %edx
	movl	%eax, 4(%esp)
	movl	84(%esp), %eax
	movl	%eax, 8(%esp)
	movl	88(%esp), %eax
	movl	%eax, 12(%esp)
	movl	92(%esp), %eax
	movl	%eax, 16(%esp)
	movl	96(%esp), %eax
	movl	%eax, 20(%esp)
	movl	100(%esp), %eax
	movl	%eax, 24(%esp)
	movl	104(%esp), %eax
	movl	%eax, 28(%esp)
	movl	108(%esp), %eax
	movl	%eax, 32(%esp)
	movl	112(%esp), %eax
	movl	%eax, 36(%esp)
	movl	116(%esp), %eax
	movl	%eax, 40(%esp)
	movl	120(%esp), %eax
	movl	%eax, 44(%esp)
	movl	124(%esp), %eax
	movl	%eax, 48(%esp)
	movl	128(%esp), %eax
	movl	%eax, 52(%esp)
	movl	132(%esp), %eax
	fldl	24(%esp)
	movl	%eax, 56(%esp)
	movl	136(%esp), %eax
	flds	LC7
	movl	%eax, 60(%esp)
	fld	%st(1)
	movl	72(%esp), %eax
	fmul	%st(1), %st
	movl	(%eax), %eax
	fldl	(%eax,%ecx,8)
	fsubl	(%eax,%edx,8)
	fucomi	%st(1), %st
	ja	L90
	fld	%st(3)
	fchs
	fmulp	%st, %st(3)
	fxch	%st(2)
	fucomip	%st(2), %st
	ja	L106
	fxch	%st(1)
L103:
	fstpl	(%ebx)
L92:
	fldl	8(%eax,%ecx,8)
	fsubl	8(%eax,%edx,8)
	fucomi	%st(1), %st
	ja	L93
	fld	%st(2)
	fchs
	fmuls	LC7
	fucomip	%st(1), %st
	jbe	L104
	fadd	%st(2), %st
	fstpl	8(%ebx)
	fldl	16(%eax,%ecx,8)
	fsubl	16(%eax,%edx,8)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	L101
L107:
	fsubp	%st, %st(1)
	fstpl	16(%ebx)
	addl	$64, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 4
	ret
	.p2align 4,,10
L93:
	.cfi_restore_state
	fsub	%st(2), %st
L104:
	fstpl	8(%ebx)
	fldl	16(%eax,%ecx,8)
	fsubl	16(%eax,%edx,8)
	fucomi	%st(1), %st
	fstp	%st(1)
	ja	L107
L101:
	fld	%st(1)
	fchs
	fmuls	LC7
	fucomip	%st(1), %st
	jbe	L108
	faddp	%st, %st(1)
	jmp	L105
	.p2align 4,,10
L108:
	fstp	%st(1)
L105:
	fstpl	16(%ebx)
	addl	$64, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 4
	ret
	.p2align 4,,10
L106:
	.cfi_restore_state
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	(%ebx)
	jmp	L92
	.p2align 4,,10
L90:
	fstp	%st(2)
	fxch	%st(1)
	fsub	%st(2), %st
	jmp	L103
	.cfi_endproc
LFE13:
	.section	.text.unlikely,"x"
LCOLDE14:
	.text
LHOTE14:
	.section .rdata,"dr"
	.align 4
LC15:
	.ascii "list_pointer < neigh_arr_length\0"
	.section	.text.unlikely,"x"
LCOLDB18:
	.text
LHOTB18:
	.p2align 4,,15
	.globl	_update_neigh_list
	.def	_update_neigh_list;	.scl	2;	.type	32;	.endef
_update_neigh_list:
LFB14:
	.cfi_startproc
	pushl	%ebp
	.cfi_def_cfa_offset 8
	.cfi_offset 5, -8
	pushl	%edi
	.cfi_def_cfa_offset 12
	.cfi_offset 7, -12
	pushl	%esi
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushl	%ebx
	.cfi_def_cfa_offset 20
	.cfi_offset 3, -20
	subl	$124, %esp
	.cfi_def_cfa_offset 144
	movl	144(%esp), %eax
	movl	212(%esp), %esi
	movl	220(%esp), %ecx
	movl	%eax, 40(%esp)
	movl	216(%esp), %eax
	movl	%eax, 44(%esp)
	movl	148(%esp), %eax
	movl	%eax, 48(%esp)
	movl	152(%esp), %eax
	movl	%eax, 52(%esp)
	movl	156(%esp), %eax
	movl	%eax, 56(%esp)
	movl	160(%esp), %eax
	movl	%eax, 60(%esp)
	movl	164(%esp), %eax
	movl	%eax, 64(%esp)
	movl	168(%esp), %eax
	movl	%eax, 68(%esp)
	movl	172(%esp), %eax
	movl	%eax, 72(%esp)
	movl	176(%esp), %eax
	movl	%eax, 76(%esp)
	movl	180(%esp), %eax
	movl	%eax, 80(%esp)
	movl	184(%esp), %eax
	movl	%eax, 84(%esp)
	movl	188(%esp), %eax
	movl	%eax, 88(%esp)
	movl	192(%esp), %eax
	movl	%eax, 92(%esp)
	movl	196(%esp), %eax
	movl	%eax, 96(%esp)
	movl	200(%esp), %eax
	movl	%eax, 100(%esp)
	movl	204(%esp), %eax
	movl	%eax, 104(%esp)
	movl	208(%esp), %eax
	movl	%eax, 108(%esp)
	leal	1(%ecx), %eax
	imull	%ecx, %eax
	movl	%eax, %ebx
	shrl	$31, %ebx
	addl	%ebx, %eax
	sarl	%eax
	movl	%eax, %ebx
	leal	(%esi,%eax,4), %edx
	movl	%esi, %eax
	testl	%ebx, %ebx
	jle	L114
	.p2align 4,,10
L137:
	movl	$-255, (%eax)
	addl	$4, %eax
	cmpl	%eax, %edx
	jne	L137
L114:
	testl	%ecx, %ecx
	jle	L111
	movl	$0, 36(%esp)
	addl	$1, 36(%esp)
	xorl	%edi, %edi
	movl	36(%esp), %eax
	xorl	%edx, %edx
	cmpl	%eax, %ecx
	je	L115
	.p2align 4,,10
L149:
	cmpl	%edx, %ebx
	jle	L117
	movl	40(%esp), %eax
	movl	(%eax), %ebp
	leal	0(%ebp,%edi), %eax
	fldl	(%eax)
	fldl	8(%ebp,%edi)
	fstpl	16(%esp)
	fldl	72(%esp)
	fldl	16(%ebp,%edi)
	movl	36(%esp), %ebp
	fstpl	24(%esp)
	fld	%st(0)
	fmuls	LC7
	fld	%st(1)
	fchs
	jmp	L118
	.p2align 4,,10
L146:
	fld	%st(1)
	fmuls	LC7
	fucomip	%st(1), %st
	jbe	L120
	fadd	%st(3), %st
L120:
	fldl	16(%esp)
	fsubl	32(%eax)
	fucomi	%st(3), %st
	ja	L122
L147:
	fld	%st(2)
	fmuls	LC7
	fucomip	%st(1), %st
	jbe	L123
	fadd	%st(4), %st
L123:
	fldl	24(%esp)
	fsubl	40(%eax)
	fucomi	%st(4), %st
	ja	L125
L148:
	fld	%st(3)
	fmuls	LC7
	fucomip	%st(1), %st
	jbe	L150
	fadd	%st(5), %st
	fxch	%st(2)
	jmp	L126
	.p2align 4,,10
L150:
	fxch	%st(2)
L126:
	fmul	%st(0), %st
	fxch	%st(1)
	fmul	%st(0), %st
	faddp	%st, %st(1)
	fxch	%st(1)
	fmul	%st(0), %st
	faddp	%st, %st(1)
	fldl	LC17
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	L128
	movl	%ebp, (%esi,%edx,4)
	addl	$1, %edx
L128:
	addl	$1, %ebp
	cmpl	%ebp, %ecx
	je	L145
	addl	$24, %eax
	cmpl	%ebx, %edx
	jge	L151
L118:
	fldl	24(%eax)
	fsubr	%st(4), %st
	fucomi	%st(2), %st
	jbe	L146
	fsub	%st(3), %st
	fldl	16(%esp)
	fsubl	32(%eax)
	fucomi	%st(3), %st
	jbe	L147
L122:
	fsub	%st(4), %st
	fldl	24(%esp)
	fsubl	40(%eax)
	fucomi	%st(4), %st
	jbe	L148
L125:
	fsub	%st(5), %st
	fxch	%st(2)
	jmp	L126
	.p2align 4,,10
L145:
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	addl	$1, 36(%esp)
	movl	$-1, (%esi,%edx,4)
	addl	$24, %edi
	movl	36(%esp), %eax
	addl	$1, %edx
	cmpl	%eax, %ecx
	jne	L149
L115:
	movl	$-1, (%esi,%edx,4)
L111:
	leal	(%ecx,%ecx,2), %eax
	testl	%eax, %eax
	jle	L140
	sall	$3, %eax
	movl	$0, 148(%esp)
	movl	%eax, 152(%esp)
	movl	44(%esp), %eax
	movl	%eax, 144(%esp)
	addl	$124, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 20
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 16
	popl	%esi
	.cfi_restore 6
	.cfi_def_cfa_offset 12
	popl	%edi
	.cfi_restore 7
	.cfi_def_cfa_offset 8
	popl	%ebp
	.cfi_restore 5
	.cfi_def_cfa_offset 4
	jmp	_memset
L140:
	.cfi_restore_state
	addl	$124, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 20
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 16
	popl	%esi
	.cfi_restore 6
	.cfi_def_cfa_offset 12
	popl	%edi
	.cfi_restore 7
	.cfi_def_cfa_offset 8
	popl	%ebp
	.cfi_restore 5
	.cfi_def_cfa_offset 4
	ret
L151:
	.cfi_restore_state
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
L117:
	movl	$LC15, 12(%esp)
	movl	$___func__.3734, 8(%esp)
	movl	$370, 4(%esp)
	movl	$LC4, (%esp)
	call	___assert_func
	.cfi_endproc
LFE14:
	.section	.text.unlikely,"x"
LCOLDE18:
	.text
LHOTE18:
	.section .rdata,"dr"
LC19:
	.ascii "w\0"
LC20:
	.ascii "dyn.out\0"
LC21:
	.ascii "temp.out\0"
LC22:
	.ascii "Error in allocating memory.\12\0"
	.align 4
LC23:
	.ascii "Time Step, Particle Number, Particle Type, Position, Velocity, Acceleration\12\0"
	.align 4
LC24:
	.ascii "%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f\12\0"
LC26:
	.ascii "%d, %f\12\0"
LC27:
	.ascii "Time step %d / %d complete.\12\0"
	.section	.text.unlikely,"x"
LCOLDB28:
	.text
LHOTB28:
	.p2align 4,,15
	.globl	_evolve_system
	.def	_evolve_system;	.scl	2;	.type	32;	.endef
_evolve_system:
LFB10:
	.cfi_startproc
	pushl	%ebp
	.cfi_def_cfa_offset 8
	.cfi_offset 5, -8
	pushl	%edi
	.cfi_def_cfa_offset 12
	.cfi_offset 7, -12
	movl	$1717986919, %edx
	pushl	%esi
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushl	%ebx
	.cfi_def_cfa_offset 20
	.cfi_offset 3, -20
	subl	$252, %esp
	.cfi_def_cfa_offset 272
	movl	276(%esp), %eax
	movl	272(%esp), %esi
	movl	%eax, %ebx
	movl	%eax, 108(%esp)
	movl	%eax, 176(%esp)
	movl	280(%esp), %eax
	leal	(%ebx,%ebx,2), %edi
	movl	%eax, %ebp
	movl	%eax, 112(%esp)
	movl	%eax, 180(%esp)
	movl	284(%esp), %eax
	movl	%eax, %ecx
	movl	%eax, 184(%esp)
	movl	288(%esp), %eax
	movl	%eax, 188(%esp)
	movl	292(%esp), %eax
	movl	%eax, 192(%esp)
	movl	296(%esp), %eax
	movl	%eax, 196(%esp)
	movl	300(%esp), %eax
	movl	%eax, 200(%esp)
	movl	304(%esp), %eax
	movl	%eax, 204(%esp)
	movl	308(%esp), %eax
	movl	%eax, 208(%esp)
	movl	312(%esp), %eax
	movl	%eax, 212(%esp)
	movl	316(%esp), %eax
	movl	%eax, 216(%esp)
	movl	320(%esp), %eax
	movl	%eax, 220(%esp)
	movl	324(%esp), %eax
	fldl	192(%esp)
	movl	$LC19, 4(%esp)
	movl	$LC20, (%esp)
	movl	%ecx, 124(%esp)
	movl	%eax, 224(%esp)
	movl	328(%esp), %eax
	fstpl	96(%esp)
	movl	%eax, 228(%esp)
	movl	332(%esp), %eax
	movl	%eax, 232(%esp)
	movl	336(%esp), %eax
	movl	%eax, 236(%esp)
	movl	%ecx, %eax
	imull	%edx
	movl	%ecx, %eax
	sarl	$31, %eax
	sarl	$2, %edx
	subl	%eax, %edx
	movl	%edx, 160(%esp)
	call	_fopen
	movl	$LC19, 4(%esp)
	movl	$LC21, (%esp)
	movl	%eax, 104(%esp)
	call	_fopen
	movl	$8, 4(%esp)
	movl	%edi, (%esp)
	movl	%eax, 148(%esp)
	call	_calloc
	leal	0(%ebp,%ebp,2), %ecx
	movl	$8, 4(%esp)
	movl	%eax, 128(%esp)
	movl	%ecx, (%esp)
	movl	%ecx, 116(%esp)
	call	_calloc
	movl	%eax, 144(%esp)
	leal	1(%ebx), %eax
	movl	$4, 4(%esp)
	imull	%ebx, %eax
	movl	%eax, %edx
	shrl	$31, %edx
	addl	%edx, %eax
	sarl	%eax
	movl	%eax, (%esp)
	call	_calloc
	movl	%eax, 140(%esp)
	movl	%eax, %ebx
	leal	1(%ebp), %eax
	movl	$4, 4(%esp)
	imull	%ebp, %eax
	movl	%eax, %edx
	shrl	$31, %edx
	addl	%edx, %eax
	sarl	%eax
	movl	%eax, (%esp)
	call	_calloc
	movl	128(%esp), %edx
	movl	%eax, 172(%esp)
	testl	%edx, %edx
	je	L153
	movl	144(%esp), %ebp
	testl	%ebp, %ebp
	je	L153
	testl	%ebx, %ebx
	je	L153
	testl	%eax, %eax
	je	L153
	movl	104(%esp), %eax
	movl	$76, 8(%esp)
	movl	$1, 4(%esp)
	movl	$LC23, (%esp)
	movl	%eax, 12(%esp)
	call	_fwrite
	movl	124(%esp), %ebx
	testl	%ebx, %ebx
	jle	L157
	leal	0(,%edi,8), %eax
	fildl	108(%esp)
	movl	$0, 120(%esp)
	movl	$1, 132(%esp)
	movl	$1, 136(%esp)
	movl	%eax, 164(%esp)
	movl	116(%esp), %eax
	fstpl	152(%esp)
	sall	$3, %eax
	movl	%eax, 168(%esp)
	.p2align 4,,10
L198:
	testl	%edi, %edi
	jle	L199
	fldl	96(%esp)
	movl	4(%esi), %eax
	xorl	%edx, %edx
	movl	8(%esi), %ecx
	fmuls	LC7
	.p2align 4,,10
L158:
	fldl	(%ecx,%edx,8)
	fmul	%st(1), %st
	faddl	(%eax,%edx,8)
	fstpl	(%eax,%edx,8)
	addl	$1, %edx
	cmpl	%edx, %edi
	jne	L158
	movl	116(%esp), %ecx
	testl	%ecx, %ecx
	jle	L258
L159:
	movl	16(%esi), %eax
	movl	20(%esi), %ecx
	xorl	%edx, %edx
	movl	116(%esp), %ebx
	.p2align 4,,10
L162:
	fldl	(%ecx,%edx,8)
	fmul	%st(1), %st
	faddl	(%eax,%edx,8)
	fstpl	(%eax,%edx,8)
	addl	$1, %edx
	cmpl	%edx, %ebx
	jg	L162
	fstp	%st(0)
	testl	%edi, %edi
	jle	L234
	movl	4(%esi), %eax
	jmp	L160
L258:
	fstp	%st(0)
	.p2align 4,,10
L160:
	movl	(%esi), %ebx
	movl	164(%esp), %ebp
	movl	128(%esp), %edx
	movl	%ebx, %ecx
	addl	%eax, %ebp
	.p2align 4,,10
L164:
	fldl	96(%esp)
	addl	$8, %eax
	addl	$8, %ecx
	addl	$8, %edx
	fmull	-8(%eax)
	fldl	-8(%ecx)
	fadd	%st(1), %st
	fstpl	-8(%ecx)
	faddl	-8(%edx)
	fstpl	-8(%edx)
	cmpl	%ebp, %eax
	jne	L164
	movl	116(%esp), %eax
	testl	%eax, %eax
	jle	L243
	movl	16(%esi), %eax
L208:
	movl	12(%esi), %ebx
	movl	168(%esp), %ebp
	movl	144(%esp), %edx
	movl	%ebx, %ecx
	addl	%eax, %ebp
	.p2align 4,,10
L166:
	fldl	96(%esp)
	addl	$8, %eax
	addl	$8, %ecx
	addl	$8, %edx
	fmull	-8(%eax)
	fldl	-8(%ecx)
	fadd	%st(1), %st
	fstpl	-8(%ecx)
	faddl	-8(%edx)
	fstpl	-8(%edx)
	cmpl	%ebp, %eax
	jne	L166
	movl	108(%esp), %eax
	fldl	96(%esp)
	testl	%edi, %edi
	fstpl	192(%esp)
	movl	%eax, 176(%esp)
	movl	112(%esp), %eax
	fldl	200(%esp)
	movl	%eax, 180(%esp)
	movl	124(%esp), %eax
	movl	%eax, 184(%esp)
	jle	L202
	movl	(%esi), %ebx
L201:
	xorl	%eax, %eax
	fldz
	jmp	L172
	.p2align 4,,10
L253:
	fsub	%st(2), %st
	fstpl	(%ebx,%eax,8)
	jmp	L170
	.p2align 4,,10
L260:
	fstp	%st(0)
L170:
	addl	$1, %eax
	cmpl	%eax, %edi
	jle	L259
L172:
	fldl	(%ebx,%eax,8)
	fucomi	%st(2), %st
	ja	L253
	fld	%st(1)
	fucomip	%st(1), %st
	jbe	L260
	fadd	%st(2), %st
	fstpl	(%ebx,%eax,8)
	addl	$1, %eax
	cmpl	%eax, %edi
	jg	L172
	fstp	%st(0)
	jmp	L252
	.p2align 4,,10
L259:
	fstp	%st(0)
L252:
	movl	116(%esp), %edx
	testl	%edx, %edx
	jle	L261
	movl	12(%esi), %ebx
L202:
	xorl	%eax, %eax
	fldz
	movl	116(%esp), %edx
	jmp	L180
	.p2align 4,,10
L254:
	fsub	%st(2), %st
	fstpl	(%ebx,%eax,8)
	jmp	L178
	.p2align 4,,10
L263:
	fstp	%st(0)
L178:
	addl	$1, %eax
	cmpl	%eax, %edx
	jle	L262
L180:
	fldl	(%ebx,%eax,8)
	fucomi	%st(2), %st
	ja	L254
	fld	%st(1)
	fucomip	%st(1), %st
	jbe	L263
	fadd	%st(2), %st
	fstpl	(%ebx,%eax,8)
	addl	$1, %eax
	cmpl	%eax, %edx
	jg	L180
	fstp	%st(0)
	fstp	%st(0)
	jmp	L181
L261:
	fstp	%st(0)
	jmp	L181
	.p2align 4,,10
L262:
	fstp	%st(0)
	fstp	%st(0)
L181:
	movl	136(%esp), %eax
	testl	%eax, %eax
	jne	L255
	movl	132(%esp), %eax
	testl	%eax, %eax
	jne	L256
L182:
	movl	140(%esp), %eax
	movl	%eax, 68(%esp)
	movl	176(%esp), %eax
	movl	%eax, 4(%esp)
	movl	180(%esp), %eax
	movl	%eax, 8(%esp)
	movl	184(%esp), %eax
	movl	%eax, 12(%esp)
	movl	188(%esp), %eax
	movl	%eax, 16(%esp)
	movl	192(%esp), %eax
	movl	%eax, 20(%esp)
	movl	196(%esp), %eax
	movl	%eax, 24(%esp)
	movl	200(%esp), %eax
	movl	%eax, 28(%esp)
	movl	204(%esp), %eax
	movl	%eax, 32(%esp)
	movl	208(%esp), %eax
	movl	%eax, 36(%esp)
	movl	212(%esp), %eax
	movl	%eax, 40(%esp)
	movl	216(%esp), %eax
	movl	%eax, 44(%esp)
	movl	220(%esp), %eax
	movl	%eax, 48(%esp)
	movl	224(%esp), %eax
	movl	%eax, 52(%esp)
	movl	228(%esp), %eax
	movl	%eax, 56(%esp)
	movl	232(%esp), %eax
	movl	%eax, 60(%esp)
	movl	236(%esp), %eax
	movl	%esi, (%esp)
	movl	%eax, 64(%esp)
	call	_calculate_acc
	testl	%edi, %edi
	jle	L183
	fldl	96(%esp)
	movl	4(%esi), %edx
	xorl	%eax, %eax
	movl	8(%esi), %ecx
	fmuls	LC7
	.p2align 4,,10
L184:
	fldl	(%ecx,%eax,8)
	fmul	%st(1), %st
	faddl	(%edx,%eax,8)
	fstpl	(%edx,%eax,8)
	addl	$1, %eax
	cmpl	%eax, %edi
	jne	L184
	movl	116(%esp), %ebp
	testl	%ebp, %ebp
	jle	L264
L205:
	movl	16(%esi), %edx
	movl	20(%esi), %ecx
	xorl	%eax, %eax
	movl	116(%esp), %ebx
	.p2align 4,,10
L186:
	fldl	(%ecx,%eax,8)
	fmul	%st(1), %st
	faddl	(%edx,%eax,8)
	fstpl	(%edx,%eax,8)
	addl	$1, %eax
	cmpl	%eax, %ebx
	jg	L186
	fstp	%st(0)
	testl	%edi, %edi
	jle	L204
	jmp	L187
L264:
	fstp	%st(0)
	.p2align 4,,10
L187:
	fldz
	movl	128(%esp), %edx
	xorl	%eax, %eax
	fld	%st(0)
	jmp	L189
	.p2align 4,,10
L265:
	fxch	%st(1)
L189:
	fldl	(%edx,%eax,8)
	fucomi	%st(2), %st
	fld	%st(2)
	fcmovnbe	%st(1), %st
	fstp	%st(1)
	fxch	%st(1)
	fcmovnbe	%st(2), %st
	fstp	%st(2)
	addl	$1, %eax
	cmpl	%eax, %edi
	jg	L265
	faddp	%st, %st(1)
	xorl	%eax, %eax
	fld1
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	seta	%al
	movl	%eax, 136(%esp)
L188:
	movl	108(%esp), %eax
	testl	%eax, %eax
	jle	L190
	movl	%edi, 132(%esp)
	movl	120(%esp), %edi
	xorl	%ebx, %ebx
	xorl	%ebp, %ebp
	.p2align 4,,10
L191:
	movl	8(%esi), %ecx
	movl	(%esi), %eax
	movl	4(%esi), %edx
	fldl	16(%ecx,%ebx)
	fstpl	84(%esp)
	fldl	8(%ecx,%ebx)
	fstpl	76(%esp)
	fldl	(%ecx,%ebx)
	fstpl	68(%esp)
	fldl	16(%edx,%ebx)
	fstpl	60(%esp)
	fldl	8(%edx,%ebx)
	fstpl	52(%esp)
	fldl	(%edx,%ebx)
	fstpl	44(%esp)
	fldl	16(%eax,%ebx)
	fstpl	36(%esp)
	fldl	8(%eax,%ebx)
	fstpl	28(%esp)
	fldl	(%eax,%ebx)
	addl	$24, %ebx
	movl	104(%esp), %eax
	movl	%ebp, 12(%esp)
	addl	$1, %ebp
	movl	$0, 16(%esp)
	movl	%edi, 8(%esp)
	fstpl	20(%esp)
	movl	$LC24, 4(%esp)
	movl	%eax, (%esp)
	call	_fprintf
	cmpl	108(%esp), %ebp
	jne	L191
	movl	112(%esp), %ebx
	movl	132(%esp), %edi
	testl	%ebx, %ebx
	jle	L194
L206:
	movl	%edi, 132(%esp)
	movl	120(%esp), %edi
	xorl	%ebx, %ebx
	xorl	%ebp, %ebp
	.p2align 4,,10
L193:
	movl	20(%esi), %ecx
	movl	12(%esi), %eax
	movl	16(%esi), %edx
	fldl	16(%ecx,%ebx)
	fstpl	84(%esp)
	fldl	8(%ecx,%ebx)
	fstpl	76(%esp)
	fldl	(%ecx,%ebx)
	fstpl	68(%esp)
	fldl	16(%edx,%ebx)
	fstpl	60(%esp)
	fldl	8(%edx,%ebx)
	fstpl	52(%esp)
	fldl	(%edx,%ebx)
	fstpl	44(%esp)
	fldl	16(%eax,%ebx)
	fstpl	36(%esp)
	fldl	8(%eax,%ebx)
	fstpl	28(%esp)
	fldl	(%eax,%ebx)
	addl	$24, %ebx
	movl	104(%esp), %eax
	movl	%ebp, 12(%esp)
	addl	$1, %ebp
	movl	$1, 16(%esp)
	movl	%edi, 8(%esp)
	fstpl	20(%esp)
	movl	$LC24, 4(%esp)
	movl	%eax, (%esp)
	call	_fprintf
	cmpl	112(%esp), %ebp
	jl	L193
	movl	108(%esp), %eax
	movl	132(%esp), %edi
	testl	%eax, %eax
	jle	L207
L194:
	fldz
	movl	4(%esi), %eax
	movl	108(%esp), %ecx
	xorl	%edx, %edx
	flds	LC7
	.p2align 4,,10
L196:
	fldl	(%eax)
	addl	$1, %edx
	addl	$24, %eax
	fldl	-16(%eax)
	fldl	-8(%eax)
	fxch	%st(2)
	cmpl	%ecx, %edx
	fmul	%st(0), %st
	fxch	%st(1)
	fmul	%st(0), %st
	faddp	%st, %st(1)
	fxch	%st(1)
	fmul	%st(0), %st
	faddp	%st, %st(1)
	fmul	%st(1), %st
	faddp	%st, %st(2)
	jl	L196
	fstp	%st(0)
L195:
	fldl	152(%esp)
	movl	120(%esp), %ebx
	movl	148(%esp), %eax
	movl	$LC26, 4(%esp)
	fmuls	LC25
	movl	%ebx, 8(%esp)
	movl	%eax, (%esp)
	fdivrp	%st, %st(1)
	fstpl	12(%esp)
	call	_fprintf
	movl	%ebx, %eax
	cltd
	idivl	160(%esp)
	testl	%edx, %edx
	je	L257
L197:
	addl	$1, 120(%esp)
	movl	$0, 132(%esp)
	movl	120(%esp), %eax
	cmpl	124(%esp), %eax
	jne	L198
L157:
	movl	104(%esp), %eax
	movl	%eax, (%esp)
	call	_fclose
	movl	148(%esp), %eax
	movl	%eax, (%esp)
	call	_fclose
	movl	128(%esp), %eax
	movl	%eax, (%esp)
	call	_free
	movl	144(%esp), %eax
	movl	%eax, (%esp)
	call	_free
	movl	140(%esp), %eax
	movl	%eax, (%esp)
	call	_free
	movl	172(%esp), %eax
	movl	%eax, 272(%esp)
	addl	$252, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 20
	popl	%ebx
	.cfi_restore 3
	.cfi_def_cfa_offset 16
	popl	%esi
	.cfi_restore 6
	.cfi_def_cfa_offset 12
	popl	%edi
	.cfi_restore 7
	.cfi_def_cfa_offset 8
	popl	%ebp
	.cfi_restore 5
	.cfi_def_cfa_offset 4
	jmp	_free
	.p2align 4,,10
L257:
	.cfi_restore_state
	call	___getreent
	movl	124(%esp), %ecx
	movl	$LC27, 4(%esp)
	movl	%ecx, 12(%esp)
	movl	120(%esp), %ecx
	movl	%ecx, 8(%esp)
	movl	12(%eax), %eax
	movl	%eax, (%esp)
	call	_fprintf
	jmp	L197
	.p2align 4,,10
L255:
	movl	108(%esp), %eax
	movl	%eax, 76(%esp)
	movl	128(%esp), %eax
	movl	%eax, 72(%esp)
	movl	140(%esp), %eax
	movl	%eax, 68(%esp)
	movl	176(%esp), %eax
	movl	%eax, 4(%esp)
	movl	180(%esp), %eax
	movl	%eax, 8(%esp)
	movl	184(%esp), %eax
	movl	%eax, 12(%esp)
	movl	188(%esp), %eax
	movl	%eax, 16(%esp)
	movl	192(%esp), %eax
	movl	%eax, 20(%esp)
	movl	196(%esp), %eax
	movl	%eax, 24(%esp)
	movl	200(%esp), %eax
	movl	%eax, 28(%esp)
	movl	204(%esp), %eax
	movl	%eax, 32(%esp)
	movl	208(%esp), %eax
	movl	%eax, 36(%esp)
	movl	212(%esp), %eax
	movl	%eax, 40(%esp)
	movl	216(%esp), %eax
	movl	%eax, 44(%esp)
	movl	220(%esp), %eax
	movl	%eax, 48(%esp)
	movl	224(%esp), %eax
	movl	%eax, 52(%esp)
	movl	228(%esp), %eax
	movl	%esi, (%esp)
	movl	%eax, 56(%esp)
	movl	232(%esp), %eax
	movl	%eax, 60(%esp)
	movl	236(%esp), %eax
	movl	%eax, 64(%esp)
	call	_update_neigh_list
	movl	132(%esp), %eax
	testl	%eax, %eax
	je	L182
L256:
	movl	112(%esp), %eax
	movl	%eax, 76(%esp)
	movl	144(%esp), %eax
	movl	%eax, 72(%esp)
	movl	172(%esp), %eax
	movl	%eax, 68(%esp)
	movl	176(%esp), %eax
	movl	%eax, 4(%esp)
	movl	180(%esp), %eax
	movl	%eax, 8(%esp)
	movl	184(%esp), %eax
	movl	%eax, 12(%esp)
	movl	188(%esp), %eax
	movl	%eax, 16(%esp)
	movl	192(%esp), %eax
	movl	%eax, 20(%esp)
	movl	196(%esp), %eax
	movl	%eax, 24(%esp)
	movl	200(%esp), %eax
	movl	%eax, 28(%esp)
	movl	204(%esp), %eax
	movl	%eax, 32(%esp)
	movl	208(%esp), %eax
	movl	%eax, 36(%esp)
	movl	212(%esp), %eax
	movl	%eax, 40(%esp)
	movl	216(%esp), %eax
	movl	%eax, 44(%esp)
	movl	220(%esp), %eax
	movl	%eax, 48(%esp)
	movl	224(%esp), %eax
	movl	%eax, 52(%esp)
	movl	228(%esp), %eax
	movl	%esi, (%esp)
	movl	%eax, 56(%esp)
	movl	232(%esp), %eax
	movl	%eax, 60(%esp)
	movl	236(%esp), %eax
	movl	%eax, 64(%esp)
	call	_update_neigh_list
	jmp	L182
L183:
	movl	116(%esp), %eax
	testl	%eax, %eax
	jle	L204
	fldl	96(%esp)
	fmuls	LC7
	jmp	L205
L190:
	movl	112(%esp), %ecx
	testl	%ecx, %ecx
	jg	L206
L207:
	fldz
	jmp	L195
L204:
	movl	$0, 136(%esp)
	jmp	L188
L199:
	movl	116(%esp), %eax
	testl	%eax, %eax
	jle	L242
	fldl	96(%esp)
	fmuls	LC7
	jmp	L159
L243:
	movl	108(%esp), %eax
	fldl	96(%esp)
	fstpl	192(%esp)
	movl	%eax, 176(%esp)
	movl	112(%esp), %eax
	fldl	200(%esp)
	movl	%eax, 180(%esp)
	movl	124(%esp), %eax
	movl	%eax, 184(%esp)
	jmp	L201
L234:
	movl	116(%esp), %edx
	testl	%edx, %edx
	jg	L208
L242:
	movl	108(%esp), %eax
	fldl	96(%esp)
	fstpl	192(%esp)
	movl	%eax, 176(%esp)
	movl	112(%esp), %eax
	movl	%eax, 180(%esp)
	movl	124(%esp), %eax
	movl	%eax, 184(%esp)
	jmp	L181
L153:
	call	___getreent
	movl	12(%eax), %eax
	movl	$28, 8(%esp)
	movl	$1, 4(%esp)
	movl	$LC22, (%esp)
	movl	%eax, 12(%esp)
	call	_fwrite
	movl	$1, (%esp)
	call	_exit
	.cfi_endproc
LFE10:
	.section	.text.unlikely,"x"
LCOLDE28:
	.text
LHOTE28:
	.section	.text.unlikely,"x"
LCOLDB29:
	.text
LHOTB29:
	.p2align 4,,15
	.globl	_check_update_req
	.def	_check_update_req;	.scl	2;	.type	32;	.endef
_check_update_req:
LFB15:
	.cfi_startproc
	movl	8(%esp), %eax
	leal	(%eax,%eax,2), %edx
	testl	%edx, %edx
	jle	L269
	fldz
	movl	4(%esp), %eax
	fld	%st(0)
	leal	(%eax,%edx,8), %edx
	jmp	L268
	.p2align 4,,10
L271:
	fxch	%st(1)
L268:
	fldl	(%eax)
	fucomi	%st(2), %st
	fld	%st(2)
	fcmovnbe	%st(1), %st
	fstp	%st(1)
	fxch	%st(1)
	fcmovnbe	%st(2), %st
	fstp	%st(2)
	addl	$8, %eax
	cmpl	%eax, %edx
	jne	L271
	faddp	%st, %st(1)
	xorl	%eax, %eax
	fld1
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	seta	%al
	ret
L269:
	xorl	%eax, %eax
	ret
	.cfi_endproc
LFE15:
	.section	.text.unlikely,"x"
LCOLDE29:
	.text
LHOTE29:
	.section	.text.unlikely,"x"
LCOLDB30:
	.text
LHOTB30:
	.p2align 4,,15
	.globl	_calc_temp
	.def	_calc_temp;	.scl	2;	.type	32;	.endef
_calc_temp:
LFB16:
	.cfi_startproc
	subl	$76, %esp
	.cfi_def_cfa_offset 80
	movl	88(%esp), %eax
	movl	84(%esp), %ecx
	movl	%eax, 12(%esp)
	movl	92(%esp), %eax
	testl	%ecx, %ecx
	movl	%ecx, 8(%esp)
	movl	%eax, 16(%esp)
	movl	96(%esp), %eax
	movl	%eax, 20(%esp)
	movl	100(%esp), %eax
	movl	%eax, 24(%esp)
	movl	104(%esp), %eax
	movl	%eax, 28(%esp)
	movl	108(%esp), %eax
	movl	%eax, 32(%esp)
	movl	112(%esp), %eax
	movl	%eax, 36(%esp)
	movl	116(%esp), %eax
	movl	%eax, 40(%esp)
	movl	120(%esp), %eax
	movl	%eax, 44(%esp)
	movl	124(%esp), %eax
	movl	%eax, 48(%esp)
	movl	128(%esp), %eax
	movl	%eax, 52(%esp)
	movl	132(%esp), %eax
	movl	%eax, 56(%esp)
	movl	136(%esp), %eax
	movl	%eax, 60(%esp)
	movl	140(%esp), %eax
	movl	%eax, 64(%esp)
	movl	144(%esp), %eax
	movl	%eax, 68(%esp)
	jle	L275
	movl	80(%esp), %eax
	leal	(%ecx,%ecx,2), %edx
	fldz
	flds	LC7
	movl	4(%eax), %eax
	leal	(%eax,%edx,8), %edx
	.p2align 4,,10
L274:
	fldl	(%eax)
	addl	$24, %eax
	fldl	-16(%eax)
	fldl	-8(%eax)
	fxch	%st(2)
	cmpl	%eax, %edx
	fmul	%st(0), %st
	fxch	%st(1)
	fmul	%st(0), %st
	faddp	%st, %st(1)
	fxch	%st(1)
	fmul	%st(0), %st
	faddp	%st, %st(1)
	fmul	%st(1), %st
	faddp	%st, %st(2)
	jne	L274
	fstp	%st(0)
L273:
	movl	%ecx, 4(%esp)
	fildl	4(%esp)
	fmuls	LC25
	addl	$76, %esp
	.cfi_remember_state
	.cfi_def_cfa_offset 4
	fdivrp	%st, %st(1)
	ret
L275:
	.cfi_restore_state
	fldz
	jmp	L273
	.cfi_endproc
LFE16:
	.section	.text.unlikely,"x"
LCOLDE30:
	.text
LHOTE30:
	.section .rdata,"dr"
	.align 4
___func__.3734:
	.ascii "update_neigh_list\0"
	.align 4
___func__.3692:
	.ascii "calculate_acc\0"
	.align 4
LC7:
	.long	1056964608
	.align 8
LC11:
	.long	-4194304
	.long	1105199103
	.align 8
LC12:
	.long	-1580547965
	.long	1074509381
	.align 8
LC17:
	.long	0
	.long	1074790400
	.align 4
LC25:
	.long	1069547520
	.ident	"GCC: (GNU) 5.4.0"
	.def	_memset;	.scl	2;	.type	32;	.endef
	.def	_rand;	.scl	2;	.type	32;	.endef
	.def	___assert_func;	.scl	2;	.type	32;	.endef
	.def	_sqrt;	.scl	2;	.type	32;	.endef
	.def	_fopen;	.scl	2;	.type	32;	.endef
	.def	_calloc;	.scl	2;	.type	32;	.endef
	.def	_fwrite;	.scl	2;	.type	32;	.endef
	.def	_fprintf;	.scl	2;	.type	32;	.endef
	.def	_fclose;	.scl	2;	.type	32;	.endef
	.def	_free;	.scl	2;	.type	32;	.endef
	.def	___getreent;	.scl	2;	.type	32;	.endef
	.def	_exit;	.scl	2;	.type	32;	.endef
