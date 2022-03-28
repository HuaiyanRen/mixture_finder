# Mixture Finder

**This is an algorithm to find the optimal number of classes in mixture model.**

1. Start with class 1, find the best (lowest AIC or BIC score) model.
2. Add a class, find the best model in the new class, which give the lowest score to full model.
   
   - If the new score is lower than previous, keep the new model.

     - Fix the new model, go back and check whether there will be a lower score by changing the model in previous classes. Keep model with the lowest score.
     - Go back to step 2.
   
   - If the new score is higher than previous, keep the previous model, go to step 3.

3. Conclude current model, this is the optimal mixture model.
   
#### Mar 4th, 2020

**1_multiple_optimisations**

In this folder, I show my consideration about variation of BIC score. It is too large, so my algorithm in first script cannot garantee it can always get optimal number of classes.

**algorithm_ver2**

I set repeatitions in my algorithm, to reduce the influence of large variation.

#### Mar 8th, 2020

**algorithm_ver3**

JC, HKY and GTR models (all GTR sub models can be included in the algorithm) will be considered in each step. The one with lowest score will be kept in class of the full model.

#### Mar 18th, 2020

**algorithm2_ver1**

In this algorithm, at first, estimate the number of conponents of 1:10, and find the lowest score, in case the score is not monotone decreasing and increasing.
After got the lowest score, fix the number of coponents and attempt to replace GTR by other classes (JC, HKY) in each components, seeking possible lower score.

This idea is not good. Abandon.

#### Mar 28th, 2020

**algorithm_ver4**

In this version, after successfully adding a class, all previous classes will be check: whether there will be a better sub model in previous classes.


