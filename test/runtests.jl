using JWAS
using Base.Test

# write your own tests here
model_equations = "BW = intercept + age + sex;
                   CW = intercept + age + sex";
R               = [6.72   24.84
                   24.84  708.41]
models          = build_model(model_equations,R);
@test typeof(models)==JWAS.MME
