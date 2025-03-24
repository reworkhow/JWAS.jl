# User-defined relationship between middle layer (intermediate traits) and output layer (phenotype)

* Firstly, a function should be pre-defined by user, where each function argument represents the input from a intermediate trait. For example, with two omics, we can have `pig_growth(omics1,omics2) = sqrt(omics1^2 / (omics1^2 + omics2^2))`)
* Then, put the name of the user-define function in the `nonlinear_function` arugument when build the model.

### example(b): user-defined relationship between middle layer (intermediate traits) and output layer (phenotype)
- number of nodes in the middle layer: 2
- nonlinear function (to define relationship between middle nodes and phenotype): y = sqrt(x1^2 / (x1^2 + x2^2))
- sample the missing omics in the middle layer: Matropolis-Hastings

All the other code are the same as before, except Step 3: Build Model Equations. Note that user-defined relationship is supported in both fully-connected and partial-connected neural network. 
```julia
# Step 3: Build Model Equations
pig_growth(x1,x2) = sqrt(x1^2 / (x1^2 + x2^2))
model = build_model(model_equation,
                    num_hidden_nodes=2,
                    nonlinear_function=pig_growth);
```
