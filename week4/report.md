In general, I find this homework to be manageable. One of the reasons is that I'm much familiar with Codon's tools and rules now, making debugging and conversion much simpler. Another reason is that because we have to implement the algorithms ourself, which I find much easier than understanding other people's thousands of lines of code.

The hardest part of this assignment was backtracking for affine alignment. Even with the help of AI, I had to draw out the simple alignment matrices (upper, lower, and middle) for visualization. I learned that each matrix has different cell scores, hence it's needed to jump from one matrix to another to retrieve the corresponding letter/alignment. 

Another problem I ran into was insufficient RAM allocation on WSL, causing the system to crash. I fixed this problem by increasing the memory size to 32GB.

Last but not least, I've never realized how efficient Codon is until now. For example, the affine alignment between mt_human and mt_orang is 30 times faster using Codon compared to Python.
