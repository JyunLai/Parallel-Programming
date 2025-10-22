	.text
	.file	"test3.c"
	.globl	test3                           # -- Begin function test3
	.p2align	4, 0x90
	.type	test3,@function
test3:                                  # @test3
	.cfi_startproc
# %bb.0:
	xorpd	%xmm0, %xmm0
	xorl	%eax, %eax
	.p2align	4, 0x90
.LBB0_1:                                # =>This Loop Header: Depth=1
                                        #     Child Loop BB0_2 Depth 2
	xorl	%ecx, %ecx
	.p2align	4, 0x90
.LBB0_2:                                #   Parent Loop BB0_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	addsd	(%rdi,%rcx,8), %xmm0
	addsd	8(%rdi,%rcx,8), %xmm0
	addsd	16(%rdi,%rcx,8), %xmm0
	addsd	24(%rdi,%rcx,8), %xmm0
	addsd	32(%rdi,%rcx,8), %xmm0
	addsd	40(%rdi,%rcx,8), %xmm0
	addsd	48(%rdi,%rcx,8), %xmm0
	addsd	56(%rdi,%rcx,8), %xmm0
	addq	$8, %rcx
	cmpq	$1024, %rcx                     # imm = 0x400
	jne	.LBB0_2
# %bb.3:                                #   in Loop: Header=BB0_1 Depth=1
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
