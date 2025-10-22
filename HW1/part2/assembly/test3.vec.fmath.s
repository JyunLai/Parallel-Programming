	.text
	.file	"test3.c"
	.globl	test3                           # -- Begin function test3
	.p2align	4, 0x90
	.type	test3,@function
test3:                                  # @test3
	.cfi_startproc
# %bb.0:
	pxor	%xmm0, %xmm0
	xorl	%eax, %eax
	.p2align	4, 0x90
.LBB0_1:                                # =>This Loop Header: Depth=1
                                        #     Child Loop BB0_2 Depth 2
	movq	%xmm0, %xmm0                    # xmm0 = xmm0[0],zero
	xorpd	%xmm1, %xmm1
	xorl	%ecx, %ecx
	.p2align	4, 0x90
.LBB0_2:                                #   Parent Loop BB0_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	addpd	(%rdi,%rcx,8), %xmm0
	addpd	16(%rdi,%rcx,8), %xmm1
	addpd	32(%rdi,%rcx,8), %xmm0
	addpd	48(%rdi,%rcx,8), %xmm1
	addpd	64(%rdi,%rcx,8), %xmm0
	addpd	80(%rdi,%rcx,8), %xmm1
	addpd	96(%rdi,%rcx,8), %xmm0
	addpd	112(%rdi,%rcx,8), %xmm1
	addq	$16, %rcx
	cmpq	$1024, %rcx                     # imm = 0x400
	jne	.LBB0_2
# %bb.3:                                #   in Loop: Header=BB0_1 Depth=1
	addpd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm0                    # xmm0 = xmm0[1],xmm1[1]
	addsd	%xmm1, %xmm0
	addl	$1, %eax
	cmpl	$20000000, %eax                 # imm = 0x1312D00
	jne	.LBB0_1
# %bb.4:
	retq
.Lfunc_end0:
	.size	test3, .Lfunc_end0-test3
	.cfi_endproc
                                        # -- End function
	.ident	"clang version 11.1.0"
	.section	".note.GNU-stack","",@progbits
	.addrsig
