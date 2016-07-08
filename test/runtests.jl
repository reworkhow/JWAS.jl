using JWAS
using Base.Test

# write your own tests here
@test model_equations = "BW = intercept + age + sex;
                         CW = intercept + age + sex";
@test R               = [6.72   24.84
                         24.84  708.41]
@test models          = JWAS.MT.build_model(model_equations,R);
