### Status of Optimizations

i. [Current Method]
ii. [Status of Current Method]

## I. Current Method

First, a kernel (kernel 1) is run with `(n+1) * n / 2` threads to calculate each particles interaction with each other particle AND itself, and load this into a matrix:

```
  a         b         c
a F(a on a) F(b on a) F(c on a)
b F(a on b) F(b on b) F(c on b)
c F(a on c) F(b on c) F(c on c)
```

Because `F(x on y) = -F(y on x)`, one thread updates both of those elements of the grid.

Then, a BLAS operation is used to sum each row, getting the total force on each particle.

Finally, a kernel (kernel 2) is run with `n` threads to move and accelerate each particle.

## II. Status of Current Method

NVVP reports that kernel 2 is not a performance concerned, although dealing with edge cases requires significant branching.

Kernel 1 probably needs the most optimization, the biggest issue is probably memory access. Accessing `F(x on y)` and `F(y on x)` is non-sequential, which can be slower. It might be possible to only update the upper triangle, then generate the rest of the matrix via a BLAS operation, which might be optimized for this.
