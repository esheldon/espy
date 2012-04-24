	.file	"checksum-compare-khash.c"
	.section	.rodata
	.align 8
	.type	__ac_HASH_UPPER, @object
	.size	__ac_HASH_UPPER, 8
__ac_HASH_UPPER:
	.long	171798692
	.long	1072210903
.LC0:
	.string	"r"
.LC1:
	.string	"Could not open file: %s\n"
	.text
.globl file_open
	.type	file_open, @function
file_open:
.LFB14:
	pushq	%rbp
.LCFI0:
	movq	%rsp, %rbp
.LCFI1:
	subq	$32, %rsp
.LCFI2:
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rdi
	movl	$.LC0, %esi
	call	fopen
	movq	%rax, -8(%rbp)
	cmpq	$0, -8(%rbp)
	jne	.L2
	movq	stderr(%rip), %rdi
	movq	-24(%rbp), %rdx
	movl	$.LC1, %esi
	movl	$0, %eax
	call	fprintf
	movl	$45, %edi
	call	exit
.L2:
	movq	-8(%rbp), %rax
	leave
	ret
.LFE14:
	.size	file_open, .-file_open
	.section	.rodata
.LC2:
	.string	"%s %s"
	.text
.globl load_master
	.type	load_master, @function
load_master:
.LFB15:
	pushq	%rbp
.LCFI3:
	movq	%rsp, %rbp
.LCFI4:
	pushq	%rbx
.LCFI5:
	subq	$376, %rsp
.LCFI6:
	movq	%rdi, -376(%rbp)
	movl	$2, -40(%rbp)
	movl	$0, -356(%rbp)
	movl	$0, %eax
	call	kh_init_md5
	movq	%rax, -32(%rbp)
	movq	-376(%rbp), %rdi
	call	file_open
	movq	%rax, -24(%rbp)
	jmp	.L5
.L6:
	leaq	-304(%rbp), %rax
	leaq	2(%rax), %rdi
	call	strdup
	movq	%rax, %rsi
	leaq	-356(%rbp), %rdx
	movq	-32(%rbp), %rdi
	call	kh_put_md5
	movl	%eax, -36(%rbp)
	movq	-32(%rbp), %rax
	movq	32(%rax), %rdx
	mov	-36(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rbx
	leaq	-352(%rbp), %rdi
	call	strdup
	movq	%rax, (%rbx)
.L5:
	leaq	-304(%rbp), %rcx
	leaq	-352(%rbp), %rdx
	movq	-24(%rbp), %rdi
	movl	$.LC2, %esi
	movl	$0, %eax
	call	fscanf
	cmpl	$2, %eax
	je	.L6
	movq	-24(%rbp), %rdi
	call	fclose
	movq	-32(%rbp), %rax
	addq	$376, %rsp
	popq	%rbx
	leave
	ret
.LFE15:
	.size	load_master, .-load_master
	.type	kh_init_md5, @function
kh_init_md5:
.LFB7:
	pushq	%rbp
.LCFI7:
	movq	%rsp, %rbp
.LCFI8:
	movl	$40, %esi
	movl	$1, %edi
	call	calloc
	leave
	ret
.LFE7:
	.size	kh_init_md5, .-kh_init_md5
	.type	kh_put_md5, @function
kh_put_md5:
.LFB12:
	pushq	%rbp
.LCFI9:
	movq	%rsp, %rbp
.LCFI10:
	subq	$64, %rsp
.LCFI11:
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	movq	-40(%rbp), %rax
	movl	8(%rax), %edx
	movq	-40(%rbp), %rax
	movl	12(%rax), %eax
	cmpl	%eax, %edx
	jb	.L11
	movq	-40(%rbp), %rax
	movl	(%rax), %edx
	movq	-40(%rbp), %rax
	movl	4(%rax), %eax
	addl	%eax, %eax
	cmpl	%eax, %edx
	jbe	.L12
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	leal	-1(%rax), %esi
	movq	-40(%rbp), %rdi
	call	kh_resize_md5
	jmp	.L11
.L12:
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	leal	1(%rax), %esi
	movq	-40(%rbp), %rdi
	call	kh_resize_md5
.L11:
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	subl	$1, %eax
	movl	%eax, -4(%rbp)
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -12(%rbp)
	movl	-12(%rbp), %eax
	movl	%eax, -28(%rbp)
	movq	-48(%rbp), %rdi
	call	__ac_X31_hash_string
	movl	%eax, -20(%rbp)
	movl	-4(%rbp), %edx
	movl	-20(%rbp), %eax
	andl	%edx, %eax
	movl	%eax, -16(%rbp)
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-16(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-16(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$2, %eax
	testl	%eax, %eax
	je	.L13
	movl	-16(%rbp), %eax
	movl	%eax, -28(%rbp)
	jmp	.L14
.L13:
	movl	-20(%rbp), %eax
	movl	%eax, %edx
	shrl	$3, %edx
	movl	-20(%rbp), %eax
	sall	$3, %eax
	xorl	%edx, %eax
	orl	$1, %eax
	andl	-4(%rbp), %eax
	movl	%eax, -24(%rbp)
	movl	-16(%rbp), %eax
	movl	%eax, -8(%rbp)
	jmp	.L15
.L18:
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-16(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-16(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$1, %eax
	testb	%al, %al
	je	.L16
	movl	-16(%rbp), %eax
	movl	%eax, -12(%rbp)
.L16:
	movl	-24(%rbp), %edx
	movl	-16(%rbp), %eax
	addl	%edx, %eax
	andl	-4(%rbp), %eax
	movl	%eax, -16(%rbp)
	movl	-16(%rbp), %eax
	cmpl	-8(%rbp), %eax
	jne	.L15
	movl	-12(%rbp), %eax
	movl	%eax, -28(%rbp)
	jmp	.L17
.L15:
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-16(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-16(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$2, %eax
	testl	%eax, %eax
	jne	.L17
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-16(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-16(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$1, %eax
	testb	%al, %al
	jne	.L18
	movq	-40(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-16(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rdi
	movq	-48(%rbp), %rsi
	call	strcmp
	testl	%eax, %eax
	jne	.L18
.L17:
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-28(%rbp), %eax
	jne	.L14
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-16(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-16(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$2, %eax
	testl	%eax, %eax
	je	.L19
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-12(%rbp), %eax
	je	.L19
	movl	-12(%rbp), %eax
	movl	%eax, -28(%rbp)
	jmp	.L14
.L19:
	movl	-16(%rbp), %eax
	movl	%eax, -28(%rbp)
.L14:
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-28(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-28(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$2, %eax
	testl	%eax, %eax
	je	.L20
	movq	-40(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-28(%rbp), %eax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-48(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-28(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rsi
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-28(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-28(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	$3, %eax
	salq	%cl, %rax
	notl	%eax
	andl	%edx, %eax
	movl	%eax, (%rsi)
	movq	-40(%rbp), %rax
	movl	4(%rax), %eax
	leal	1(%rax), %edx
	movq	-40(%rbp), %rax
	movl	%edx, 4(%rax)
	movq	-40(%rbp), %rax
	movl	8(%rax), %eax
	leal	1(%rax), %edx
	movq	-40(%rbp), %rax
	movl	%edx, 8(%rax)
	movq	-56(%rbp), %rax
	movl	$1, (%rax)
	jmp	.L21
.L20:
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-28(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-28(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$1, %eax
	testb	%al, %al
	je	.L22
	movq	-40(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-28(%rbp), %eax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-48(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-28(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rsi
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-28(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-28(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	$3, %eax
	salq	%cl, %rax
	notl	%eax
	andl	%edx, %eax
	movl	%eax, (%rsi)
	movq	-40(%rbp), %rax
	movl	4(%rax), %eax
	leal	1(%rax), %edx
	movq	-40(%rbp), %rax
	movl	%edx, 4(%rax)
	movq	-56(%rbp), %rax
	movl	$2, (%rax)
	jmp	.L21
.L22:
	movq	-56(%rbp), %rax
	movl	$0, (%rax)
.L21:
	movl	-28(%rbp), %eax
	leave
	ret
.LFE12:
	.size	kh_put_md5, .-kh_put_md5
	.type	kh_resize_md5, @function
kh_resize_md5:
.LFB11:
	pushq	%rbp
.LCFI12:
	movq	%rsp, %rbp
.LCFI13:
	subq	$144, %rsp
.LCFI14:
	movq	%rdi, -72(%rbp)
	movl	%esi, -76(%rbp)
	movq	$0, -64(%rbp)
	movl	$1, -52(%rbp)
	subl	$1, -76(%rbp)
	movl	-76(%rbp), %eax
	shrl	%eax
	orl	%eax, -76(%rbp)
	movl	-76(%rbp), %eax
	shrl	$2, %eax
	orl	%eax, -76(%rbp)
	movl	-76(%rbp), %eax
	shrl	$4, %eax
	orl	%eax, -76(%rbp)
	movl	-76(%rbp), %eax
	shrl	$8, %eax
	orl	%eax, -76(%rbp)
	movl	-76(%rbp), %eax
	shrl	$16, %eax
	orl	%eax, -76(%rbp)
	addl	$1, -76(%rbp)
	cmpl	$3, -76(%rbp)
	ja	.L25
	movl	$4, -76(%rbp)
.L25:
	movq	-72(%rbp), %rax
	movl	4(%rax), %eax
	movl	%eax, -116(%rbp)
	mov	-76(%rbp), %eax
	movq	%rax, -128(%rbp)
	cmpq	$0, -128(%rbp)
	js	.L26
	cvtsi2sdq	-128(%rbp), %xmm0
	movsd	%xmm0, -112(%rbp)
	jmp	.L27
.L26:
	movq	-128(%rbp), %rax
	shrq	%rax
	movq	-128(%rbp), %rdx
	andl	$1, %edx
	orq	%rdx, %rax
	cvtsi2sdq	%rax, %xmm0
	movapd	%xmm0, %xmm1
	addsd	%xmm0, %xmm1
	movsd	%xmm1, -112(%rbp)
.L27:
	movsd	__ac_HASH_UPPER(%rip), %xmm0
	movsd	-112(%rbp), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	.LC3(%rip), %xmm0
	addsd	%xmm1, %xmm0
	cvttsd2siq	%xmm0, %rax
	cmpl	%eax, -116(%rbp)
	jb	.L28
	movl	$0, -52(%rbp)
	jmp	.L29
.L28:
	cmpl	$15, -76(%rbp)
	jbe	.L30
	movl	-76(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	movq	%rax, -104(%rbp)
	jmp	.L31
.L30:
	movq	$4, -104(%rbp)
.L31:
	movq	-104(%rbp), %rdi
	call	malloc
	movq	%rax, -64(%rbp)
	cmpl	$15, -76(%rbp)
	jbe	.L32
	movl	-76(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	movq	%rax, -96(%rbp)
	jmp	.L33
.L32:
	movq	$4, -96(%rbp)
.L33:
	movq	-64(%rbp), %rdi
	movq	-96(%rbp), %rdx
	movl	$170, %esi
	call	memset
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-76(%rbp), %eax
	jae	.L29
	mov	-76(%rbp), %eax
	leaq	0(,%rax,8), %rsi
	movq	-72(%rbp), %rax
	movq	24(%rax), %rdi
	call	realloc
	movq	%rax, %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, 24(%rax)
	mov	-76(%rbp), %eax
	leaq	0(,%rax,8), %rsi
	movq	-72(%rbp), %rax
	movq	32(%rax), %rdi
	call	realloc
	movq	%rax, %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, 32(%rax)
.L29:
	cmpl	$0, -52(%rbp)
	je	.L45
	movl	$0, -52(%rbp)
	jmp	.L35
.L41:
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-52(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-52(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$3, %eax
	testl	%eax, %eax
	jne	.L36
	movq	-72(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-52(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rax
	movq	%rax, -48(%rbp)
	movl	-76(%rbp), %eax
	subl	$1, %eax
	movl	%eax, -32(%rbp)
	movq	-72(%rbp), %rax
	movq	32(%rax), %rdx
	mov	-52(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rax
	movq	%rax, -40(%rbp)
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-52(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rsi
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-52(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-52(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	$1, %eax
	salq	%cl, %rax
	orl	%edx, %eax
	movl	%eax, (%rsi)
.L40:
	movq	-48(%rbp), %rdi
	call	__ac_X31_hash_string
	movl	%eax, -24(%rbp)
	movl	-32(%rbp), %edx
	movl	-24(%rbp), %eax
	andl	%edx, %eax
	movl	%eax, -20(%rbp)
	movl	-24(%rbp), %eax
	movl	%eax, %edx
	shrl	$3, %edx
	movl	-24(%rbp), %eax
	sall	$3, %eax
	xorl	%edx, %eax
	orl	$1, %eax
	andl	-32(%rbp), %eax
	movl	%eax, -28(%rbp)
	jmp	.L37
.L38:
	movl	-28(%rbp), %edx
	movl	-20(%rbp), %eax
	addl	%edx, %eax
	andl	-32(%rbp), %eax
	movl	%eax, -20(%rbp)
.L37:
	movl	-20(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %edx
	movl	-20(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$2, %eax
	testl	%eax, %eax
	je	.L38
	movl	-20(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	movq	%rax, %rsi
	addq	-64(%rbp), %rsi
	movl	-20(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %edx
	movl	-20(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	$2, %eax
	salq	%cl, %rax
	notl	%eax
	andl	%edx, %eax
	movl	%eax, (%rsi)
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-20(%rbp), %eax
	jbe	.L39
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-20(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-20(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$3, %eax
	testl	%eax, %eax
	jne	.L39
	movq	-72(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rax
	movq	%rax, -16(%rbp)
	movq	-72(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-48(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-16(%rbp), %rax
	movq	%rax, -48(%rbp)
	movq	-72(%rbp), %rax
	movq	32(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rax
	movq	%rax, -8(%rbp)
	movq	-72(%rbp), %rax
	movq	32(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-40(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-8(%rbp), %rax
	movq	%rax, -40(%rbp)
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-20(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rsi
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-20(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-20(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	$1, %eax
	salq	%cl, %rax
	orl	%edx, %eax
	movl	%eax, (%rsi)
	jmp	.L40
.L39:
	movq	-72(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-48(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-72(%rbp), %rax
	movq	32(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-40(%rbp), %rax
	movq	%rax, (%rdx)
.L36:
	addl	$1, -52(%rbp)
.L35:
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-52(%rbp), %eax
	jne	.L41
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-76(%rbp), %eax
	jbe	.L42
	mov	-76(%rbp), %eax
	leaq	0(,%rax,8), %rsi
	movq	-72(%rbp), %rax
	movq	24(%rax), %rdi
	call	realloc
	movq	%rax, %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, 24(%rax)
	mov	-76(%rbp), %eax
	leaq	0(,%rax,8), %rsi
	movq	-72(%rbp), %rax
	movq	32(%rax), %rdi
	call	realloc
	movq	%rax, %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, 32(%rax)
.L42:
	movq	-72(%rbp), %rax
	movq	16(%rax), %rdi
	call	free
	movq	-72(%rbp), %rdx
	movq	-64(%rbp), %rax
	movq	%rax, 16(%rdx)
	movq	-72(%rbp), %rdx
	movl	-76(%rbp), %eax
	movl	%eax, (%rdx)
	movq	-72(%rbp), %rax
	movl	4(%rax), %edx
	movq	-72(%rbp), %rax
	movl	%edx, 8(%rax)
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	mov	%eax, %eax
	movq	%rax, -136(%rbp)
	cmpq	$0, -136(%rbp)
	js	.L43
	cvtsi2sdq	-136(%rbp), %xmm0
	movsd	%xmm0, -88(%rbp)
	jmp	.L44
.L43:
	movq	-136(%rbp), %rax
	shrq	%rax
	movq	-136(%rbp), %rdx
	andl	$1, %edx
	orq	%rdx, %rax
	cvtsi2sdq	%rax, %xmm0
	movapd	%xmm0, %xmm1
	addsd	%xmm0, %xmm1
	movsd	%xmm1, -88(%rbp)
.L44:
	movsd	__ac_HASH_UPPER(%rip), %xmm0
	movsd	-88(%rbp), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	.LC3(%rip), %xmm0
	addsd	%xmm1, %xmm0
	cvttsd2siq	%xmm0, %rax
	movl	%eax, %edx
	movq	-72(%rbp), %rax
	movl	%edx, 12(%rax)
.L45:
	leave
	ret
.LFE11:
	.size	kh_resize_md5, .-kh_resize_md5
	.type	__ac_X31_hash_string, @function
__ac_X31_hash_string:
.LFB5:
	pushq	%rbp
.LCFI15:
	movq	%rsp, %rbp
.LCFI16:
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movzbl	(%rax), %eax
	movsbl	%al,%eax
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	je	.L47
	addq	$1, -24(%rbp)
	jmp	.L48
.L49:
	movl	-4(%rbp), %eax
	sall	$5, %eax
	movl	%eax, %edx
	subl	-4(%rbp), %edx
	movq	-24(%rbp), %rax
	movzbl	(%rax), %eax
	movsbl	%al,%eax
	leal	(%rdx,%rax), %eax
	movl	%eax, -4(%rbp)
	addq	$1, -24(%rbp)
.L48:
	movq	-24(%rbp), %rax
	movzbl	(%rax), %eax
	testb	%al, %al
	jne	.L49
.L47:
	movl	-4(%rbp), %eax
	leave
	ret
.LFE5:
	.size	__ac_X31_hash_string, .-__ac_X31_hash_string
	.section	.rodata
.LC4:
	.string	"%s not found in master list\n"
	.align 8
.LC5:
	.string	"%s has md5sum %s instead of %s\n"
	.text
.globl compare_md5sums
	.type	compare_md5sums, @function
compare_md5sums:
.LFB16:
	pushq	%rbp
.LCFI17:
	movq	%rsp, %rbp
.LCFI18:
	subq	$352, %rsp
.LCFI19:
	movq	%rdi, -344(%rbp)
	movq	%rsi, -352(%rbp)
	movq	-352(%rbp), %rdi
	call	file_open
	movq	%rax, -16(%rbp)
	jmp	.L52
.L54:
	leaq	-288(%rbp), %rax
	leaq	27(%rax), %rsi
	movq	-344(%rbp), %rdi
	call	kh_get_md5
	movl	%eax, -20(%rbp)
	movq	-344(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-20(%rbp), %eax
	jne	.L53
	movq	stderr(%rip), %rdi
	leaq	-288(%rbp), %rdx
	movl	$.LC4, %esi
	movl	$0, %eax
	call	fprintf
	jmp	.L52
.L53:
	movq	-344(%rbp), %rax
	movq	32(%rax), %rdx
	mov	-20(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rax
	movq	%rax, -8(%rbp)
	movq	-8(%rbp), %rsi
	leaq	-336(%rbp), %rdi
	call	strcmp
	testl	%eax, %eax
	je	.L52
	movq	stderr(%rip), %rdi
	movq	-8(%rbp), %rax
	leaq	-336(%rbp), %rcx
	leaq	-288(%rbp), %rdx
	movq	%rax, %r8
	movl	$.LC5, %esi
	movl	$0, %eax
	call	fprintf
.L52:
	leaq	-288(%rbp), %rcx
	leaq	-336(%rbp), %rdx
	movq	-16(%rbp), %rdi
	movl	$.LC2, %esi
	movl	$0, %eax
	call	fscanf
	cmpl	$2, %eax
	je	.L54
	leave
	ret
.LFE16:
	.size	compare_md5sums, .-compare_md5sums
	.type	kh_get_md5, @function
kh_get_md5:
.LFB10:
	pushq	%rbp
.LCFI20:
	movq	%rsp, %rbp
.LCFI21:
	subq	$64, %rsp
.LCFI22:
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	testl	%eax, %eax
	je	.L57
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	subl	$1, %eax
	movl	%eax, -4(%rbp)
	movq	-48(%rbp), %rdi
	call	__ac_X31_hash_string
	movl	%eax, -16(%rbp)
	movl	-4(%rbp), %edx
	movl	-16(%rbp), %eax
	andl	%edx, %eax
	movl	%eax, -12(%rbp)
	movl	-16(%rbp), %eax
	movl	%eax, %edx
	shrl	$3, %edx
	movl	-16(%rbp), %eax
	sall	$3, %eax
	xorl	%edx, %eax
	orl	$1, %eax
	andl	-4(%rbp), %eax
	movl	%eax, -20(%rbp)
	movl	-12(%rbp), %eax
	movl	%eax, -8(%rbp)
	jmp	.L58
.L61:
	movl	-20(%rbp), %edx
	movl	-12(%rbp), %eax
	addl	%edx, %eax
	andl	-4(%rbp), %eax
	movl	%eax, -12(%rbp)
	movl	-12(%rbp), %eax
	cmpl	-8(%rbp), %eax
	jne	.L58
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -56(%rbp)
	jmp	.L59
.L58:
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-12(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-12(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$2, %eax
	testl	%eax, %eax
	jne	.L60
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-12(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-12(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$1, %eax
	testb	%al, %al
	jne	.L61
	movq	-40(%rbp), %rax
	movq	24(%rax), %rdx
	mov	-12(%rbp), %eax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rdi
	movq	-48(%rbp), %rsi
	call	strcmp
	testl	%eax, %eax
	jne	.L61
.L60:
	movq	-40(%rbp), %rax
	movq	16(%rax), %rdx
	movl	-12(%rbp), %eax
	shrl	$4, %eax
	mov	%eax, %eax
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movl	(%rax), %edx
	movl	-12(%rbp), %eax
	andl	$15, %eax
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%edx, %eax
	shrl	%cl, %eax
	andl	$3, %eax
	testl	%eax, %eax
	je	.L62
	movq	-40(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -52(%rbp)
	jmp	.L63
.L62:
	movl	-12(%rbp), %eax
	movl	%eax, -52(%rbp)
.L63:
	movl	-52(%rbp), %eax
	movl	%eax, -56(%rbp)
	jmp	.L59
.L57:
	movl	$0, -56(%rbp)
.L59:
	movl	-56(%rbp), %eax
	leave
	ret
.LFE10:
	.size	kh_get_md5, .-kh_get_md5
	.section	.rodata
	.align 8
.LC6:
	.string	"usage: checksum-compare master_list test_list\n"
	.text
.globl main
	.type	main, @function
main:
.LFB17:
	pushq	%rbp
.LCFI23:
	movq	%rsp, %rbp
.LCFI24:
	subq	$32, %rsp
.LCFI25:
	movl	%edi, -20(%rbp)
	movq	%rsi, -32(%rbp)
	movq	$0, -8(%rbp)
	cmpl	$2, -20(%rbp)
	jg	.L66
	movl	$.LC6, %edi
	call	puts
	movl	$45, %edi
	call	exit
.L66:
	movq	-32(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rdi
	call	load_master
	movq	%rax, -8(%rbp)
	movq	-32(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rsi
	movq	-8(%rbp), %rdi
	call	compare_md5sums
	movl	$0, %eax
	leave
	ret
.LFE17:
	.size	main, .-main
	.section	.rodata
	.align 8
.LC3:
	.long	0
	.long	1071644672
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x1
	.byte	0x3
	.byte	0xc
	.uleb128 0x7
	.uleb128 0x8
	.byte	0x90
	.uleb128 0x1
	.align 8
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB14
	.long	.LFE14-.LFB14
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI0-.LFB14
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI1-.LCFI0
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE1:
.LSFDE3:
	.long	.LEFDE3-.LASFDE3
.LASFDE3:
	.long	.LASFDE3-.Lframe1
	.long	.LFB15
	.long	.LFE15-.LFB15
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI3-.LFB15
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI4-.LCFI3
	.byte	0xd
	.uleb128 0x6
	.byte	0x4
	.long	.LCFI6-.LCFI4
	.byte	0x83
	.uleb128 0x3
	.align 8
.LEFDE3:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB7
	.long	.LFE7-.LFB7
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI7-.LFB7
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI8-.LCFI7
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE5:
.LSFDE7:
	.long	.LEFDE7-.LASFDE7
.LASFDE7:
	.long	.LASFDE7-.Lframe1
	.long	.LFB12
	.long	.LFE12-.LFB12
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI9-.LFB12
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI10-.LCFI9
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE7:
.LSFDE9:
	.long	.LEFDE9-.LASFDE9
.LASFDE9:
	.long	.LASFDE9-.Lframe1
	.long	.LFB11
	.long	.LFE11-.LFB11
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI12-.LFB11
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI13-.LCFI12
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE9:
.LSFDE11:
	.long	.LEFDE11-.LASFDE11
.LASFDE11:
	.long	.LASFDE11-.Lframe1
	.long	.LFB5
	.long	.LFE5-.LFB5
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI15-.LFB5
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI16-.LCFI15
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE11:
.LSFDE13:
	.long	.LEFDE13-.LASFDE13
.LASFDE13:
	.long	.LASFDE13-.Lframe1
	.long	.LFB16
	.long	.LFE16-.LFB16
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI17-.LFB16
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI18-.LCFI17
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE13:
.LSFDE15:
	.long	.LEFDE15-.LASFDE15
.LASFDE15:
	.long	.LASFDE15-.Lframe1
	.long	.LFB10
	.long	.LFE10-.LFB10
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI20-.LFB10
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI21-.LCFI20
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE15:
.LSFDE17:
	.long	.LEFDE17-.LASFDE17
.LASFDE17:
	.long	.LASFDE17-.Lframe1
	.long	.LFB17
	.long	.LFE17-.LFB17
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI23-.LFB17
	.byte	0xe
	.uleb128 0x10
	.byte	0x86
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI24-.LCFI23
	.byte	0xd
	.uleb128 0x6
	.align 8
.LEFDE17:
	.ident	"GCC: (GNU) 4.3.2 20081007 (Red Hat 4.3.2-7)"
	.section	.note.GNU-stack,"",@progbits
