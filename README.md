# Mixture Finder

**This is an algorithm to find the optimal number of classes in mixture model.**

it is write in python script. The required software is IQ-TREE version2. 

**Version:**

The newest version is *mixfinder.py*, it is used in commandline.

**How to use it?**

The most commonly use is: confirm *mixfinder.py* is in the folder, then input:

```
python3 mixfinder.py -n (filename) -s (score type)
```

If you input ```python3 mixfinder.py -h```, it will explain all parameters:

   ```--file_name``` or ```-n```: the name of sequence file.
   
   ```--score_type``` or ```-s```: select one of AIC, AICc, BIC.
   
   ```--repeats``` or ```-r```: the times of fitting the model in each fitting.
   
   ```--pool_num``` or ```-p```: the maximum of multiprocessing pool.
   
There are two other parameters that commondline users should not care, they were only used them in windows:

   ```--file_position``` or ```-fp```: cmd statement that open (should include 'cd' comment) the folder including sequence file. Use '\ '. Used in windows.
  
   ```--iqtree_position``` or ```-tp```: the position of the folder including iqtree files (iqtree2, iqtree2-click, libiomp5md.dll). Use '/'. Used in windows.


**Output** 

The output *result.txt* and iqtree logs will be put in a folder named *file_name.results*.

### algorithm structure ###

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



## Some logs ##
   
#### Mar 4th, 2020

**1_multiple_optimisations**

In this folder, I show my consideration about variation of BIC score. It is too large, so my algorithm in first script cannot garantee it can always get optimal number of classes.

#### Mar 18th, 2020

**algorithm2_ver1**

In this algorithm, at first, estimate the number of conponents of 1:10, and find the lowest score, in case the score is not monotone decreasing and increasing.
After got the lowest score, fix the number of coponents and attempt to replace GTR by other classes (JC, HKY) in each components, seeking possible lower score.

This idea is not good. Abandon.
