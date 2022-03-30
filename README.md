# Mixture Finder

**This is an algorithm to find the optimal number of classes in mixture model.**

### The algotithm is written as a function in python script. ###

*To use the code, just assign the parameters (at bottom of the python file) and run the file.*

**Version:**

Now, it is version 4:  *mixture_finder/1_multiple_optimisations/scripts/algorithm_ver4.py*

### **LINUX:**

This is version 5. It's the same algorithm as version 4, but good to be run in server.

In this version, the output is a txt file: _results/result.txt_. The parameter *file_position* and *iqtree_position* don't need to change.

**Name:** 

optimize_classes()

**Parameters of function:**

*file_position*: cmd statement that open (should include 'cd' comment) the folder including sequence file. Use ' \ '.

*iqtree_position*: the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use '/'. (I am windows user.)

*file_name*: the name of sequence file.

*score_type*: '(AIC)', '(AICc)', '(BIC)'

*repeats*: the times of fitting the model in each step.

**Return:** 

*The optimal numbers of classes* and *the full model*.

### function structure ###

1. Start with class 1.

   - try all submodels of GTR model, choose the one with the lowest score.
   
2. Add a new class.
   
   - try all submodels of GTR model in the new class, choose the one with the new lowest score.
   
      - if the new score is lower than the last one:

         - keep the newly added class fixed and try to change the class 1 by all submodels of GTR model.

         - if there is a lower score, change the submodel in class 1 to corresponding sbumodel.

            - if the class 1 is changed, then do these changing to class 2,3..., until there is no change or we've changed all previous classes.

   - repeat step 2 until the score increases by adding a new class.

3. Conclude the results and return.

**NOTICE**

The parameters that including commandline should be carefully checked.

In each repeatition, the lowest score will be choose. Repeatition aims to reducing random error of running IQ-TREE.

The python module *os* and *numpy* are used in this function.

This function includes some sub functions.

## Logs ##
   
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

#### Mar 28th, 2020

**algorithm_ver5**

The linux version of ver4.

**algorithm_ver4**

Fixed some bugs.
