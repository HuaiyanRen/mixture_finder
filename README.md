# Mixture Finder

**This is an algorithm to find the optimal number of classes in mixture model.**

1. Start with class 1, find the best (lowest AIC or BIC score) model.
2. Add a class, find the best model in the new class, which give the lowest score to full model.
   
   - If the new score is lower than previous, keep the new model.

     - Fix the new model, go back and check whether there will be a lower score by changing the model in previous classes. Keep model with the lowest score.
     - Go back to step 2.
   
   - If the new score is higher than previous, keep the previous model, go to step 3.

3. Conclude current model, this is the optimal mixture model.
   
